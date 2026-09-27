from dataclasses import replace
from math import isfinite

import pytest
import torch

import oncobio_agentfm.clinical_baseline as baseline
from oncobio_agentfm.clinical import NSCLCCohort, StageGroup
from oncobio_agentfm.clinical_baseline import (
    FEATURE_NAMES,
    ClinicalEncoder,
    fit_clinical_cox,
)
from oncobio_agentfm.splitting import CohortSplit, split_patient_ids


@pytest.fixture
def study(make_record):
    stage_groups = (StageGroup.ONE, StageGroup.TWO, StageGroup.THREE, StageGroup.FOUR)
    records = (
        *(
            make_record(
                f"T{position:02d}",
                age_years=55 + 10 * (position % 3),
                stage=stage_groups[position // 3],
                ecog=position % 3,
                event=position not in (2, 8),
                follow_up_days=150 - (position // 3) * 35 + (position % 3) * 4,
            )
            for position in range(12)
        ),
        make_record("V00", age_years=60, stage=None, ecog=None, event=False),
        make_record("V01", age_years=72, stage=StageGroup.FOUR, ecog=2),
        make_record("E00", age_years=63, stage=StageGroup.TWO, ecog=1),
    )
    cohort = NSCLCCohort(records=records, ajcc_version=9)
    split = CohortSplit(
        train=tuple(f"T{position:02d}" for position in range(12)),
        validation=("V00", "V01"),
        test=("E00",),
    )
    return cohort, split


def test_encoder_uses_training_age_and_distinguishes_missing_from_zero(study) -> None:
    cohort, split = study
    train = tuple(record.baseline for record in cohort.records if record.patient_id in split.train)
    encoder = ClinicalEncoder.fit(train, ajcc_version=cohort.ajcc_version)
    missing = next(record.baseline for record in cohort.records if record.patient_id == "V00")
    observed = next(record.baseline for record in cohort.records if record.patient_id == "T00")
    missing_features, observed_features = encoder.transform((missing, observed))

    assert encoder.training_age_mean == pytest.approx(65.0)
    assert missing_features[FEATURE_NAMES.index("stage_missing")].item() == 1.0
    assert missing_features[FEATURE_NAMES.index("ecog_missing")].item() == 1.0
    assert missing_features[FEATURE_NAMES.index("ecog_0")].item() == 0.0
    assert observed_features[FEATURE_NAMES.index("stage_I")].item() == 1.0
    assert observed_features[FEATURE_NAMES.index("ecog_0")].item() == 1.0
    assert observed_features[FEATURE_NAMES.index("ecog_missing")].item() == 0.0
    assert encoder.transform(()).shape == (0, len(FEATURE_NAMES))


def test_encoder_rejects_mixed_training_stage_versions(study) -> None:
    cohort, split = study
    train = tuple(record.baseline for record in cohort.records if record.patient_id in split.train)
    with pytest.raises(ValueError, match="stage version"):
        ClinicalEncoder.fit(train, ajcc_version=8)
    with pytest.raises(ValueError, match="finite"):
        ClinicalEncoder(training_age_mean=float("nan"), ajcc_version=9)


def test_full_training_risk_set_includes_censored_patients(study, monkeypatch) -> None:
    cohort, split = study
    original_loss = baseline.efron_negative_partial_log_likelihood
    observed_inputs: list[tuple[int, int, int]] = []

    def observe(log_risk, durations, events):
        observed_inputs.append((len(log_risk), len(durations), int(events.sum())))
        return original_loss(log_risk, durations, events)

    monkeypatch.setattr(baseline, "efron_negative_partial_log_likelihood", observe)
    fitted = fit_clinical_cox(cohort, split)

    assert observed_inputs
    assert set(observed_inputs) == {(12, 12, 10)}
    assert (fitted.train_patients, fitted.train_deaths) == (12, 10)


def test_held_out_covariates_and_outcomes_cannot_change_fit(study, make_record) -> None:
    cohort, split = study
    first = fit_clinical_cox(cohort, split)
    changed = tuple(
        make_record(
            record.patient_id,
            age_years=120,
            stage=StageGroup.ONE,
            ecog=4,
            event=False,
            follow_up_days=1,
        )
        if record.patient_id not in split.train
        else record
        for record in cohort.records
    )
    second = fit_clinical_cox(NSCLCCohort(records=changed, ajcc_version=9), split)

    assert first == second
    assert first.encoder.training_age_mean == 65.0


def test_fit_is_order_invariant_and_does_not_change_global_torch_rng(study) -> None:
    cohort, split = study
    before = torch.random.get_rng_state().clone()
    first = fit_clinical_cox(cohort, split)
    after = torch.random.get_rng_state()
    reversed_split = CohortSplit(
        train=tuple(reversed(split.train)),
        validation=tuple(reversed(split.validation)),
        test=split.test,
    )
    second = fit_clinical_cox(
        NSCLCCohort(records=tuple(reversed(cohort.records)), ajcc_version=9), reversed_split
    )

    assert torch.equal(before, after)
    assert first == second


def test_existing_patient_split_integrates_with_clinical_baseline(study) -> None:
    cohort, _ = study
    split = split_patient_ids(cohort.patient_ids, seed="synthetic-cohort-check")
    fitted = fit_clinical_cox(cohort, split)

    assert split.all_patient_ids == cohort.patient_ids
    assert fitted.train_patients == len(split.train)
    assert fitted.train_deaths > 0


def test_fitted_scores_reflect_a_controlled_stage_difference(study, make_record) -> None:
    cohort, split = study
    fitted = fit_clinical_cox(cohort, split)
    lower = make_record("P-low", age_years=65, stage=StageGroup.ONE, ecog=0).baseline
    higher = make_record("P-high", age_years=65, stage=StageGroup.FOUR, ecog=0).baseline

    assert len(fitted.coefficient_table()) == len(FEATURE_NAMES)
    assert all(isfinite(value) for _, value in fitted.coefficient_table())
    assert fitted.predict_log_risk(higher) > fitted.predict_log_risk(lower)
    assert isfinite(fitted.training_objective)
    assert fitted.training_gradient_max < 1e-5


def test_training_improves_over_an_uninformative_zero_risk_baseline(study) -> None:
    cohort, split = study
    fitted = fit_clinical_cox(cohort, split)
    training = tuple(record for record in cohort.records if record.patient_id in split.train)
    durations = torch.tensor(
        [record.outcome.duration_days for record in training], dtype=torch.float64
    )
    events = torch.tensor([record.outcome.event_observed for record in training], dtype=torch.int64)
    null_loss = baseline.efron_negative_partial_log_likelihood(
        torch.zeros(len(training), dtype=torch.float64), durations, events
    )

    assert fitted.training_objective < null_loss.item()


def test_fit_rejects_unmatched_split_and_zero_event_training(study) -> None:
    cohort, split = study
    incomplete = CohortSplit(train=split.train, validation=split.validation, test=())
    with pytest.raises(ValueError, match="exactly"):
        fit_clinical_cox(cohort, incomplete)

    without_deaths = NSCLCCohort(
        records=tuple(
            replace(
                record,
                outcome=replace(
                    record.outcome,
                    event_date=None,
                    censor_date=record.outcome.end_date,
                ),
            )
            if record.patient_id in split.train and record.outcome.event_observed
            else record
            for record in cohort.records
        ),
        ajcc_version=9,
    )
    with pytest.raises(ValueError, match="observed death"):
        fit_clinical_cox(without_deaths, split)


def test_fit_rejects_uninformative_risk_set_and_nonconvergence(study, make_record) -> None:
    same_day = NSCLCCohort(
        records=(
            make_record("A", event=True, follow_up_days=30),
            make_record("B", event=False, follow_up_days=30),
        ),
        ajcc_version=9,
    )
    with pytest.raises(ValueError, match="no comparable"):
        fit_clinical_cox(same_day, CohortSplit(train=("A", "B"), validation=(), test=()))

    with pytest.raises(RuntimeError, match="did not converge"):
        fit_clinical_cox(*study, max_iterations=1)


@pytest.mark.parametrize("penalty", [0.0, -1.0, True, float("nan"), float("inf")])
def test_fit_rejects_invalid_regularization(study, penalty) -> None:
    with pytest.raises(ValueError, match="ridge_strength"):
        fit_clinical_cox(*study, ridge_strength=penalty)


def test_fitted_encoder_rejects_incompatible_stage_version(study, make_record) -> None:
    cohort, split = study
    fitted = fit_clinical_cox(cohort, split)
    incompatible = make_record("P-old", stage_version=8).baseline
    with pytest.raises(ValueError, match="stage version"):
        fitted.predict_log_risk(incompatible)


def test_prediction_refuses_categories_not_observed_in_training(study, make_record) -> None:
    cohort, split = study
    fitted = fit_clinical_cox(cohort, split)
    with pytest.raises(ValueError, match="stage category"):
        fitted.predict_log_risk(make_record("P-no-stage", stage=None, ecog=0).baseline)
    with pytest.raises(ValueError, match="ECOG category"):
        fitted.predict_log_risk(make_record("P-ecog-four", stage=StageGroup.ONE, ecog=4).baseline)
    with pytest.raises(ValueError, match="training range"):
        fitted.predict_log_risk(
            make_record("P-older", age_years=120, stage=StageGroup.ONE, ecog=0).baseline
        )

    supported = next(record.baseline for record in cohort.records if record.patient_id == "V01")
    assert isfinite(fitted.predict_log_risk(supported))


def test_fitted_parameters_reject_uninterpretable_shapes(study) -> None:
    cohort, split = study
    fitted = fit_clinical_cox(cohort, split)
    with pytest.raises(ValueError, match="coefficients"):
        replace(fitted, coefficients=(0.0,))
    with pytest.raises(ValueError, match="event count"):
        replace(fitted, train_deaths=0)
