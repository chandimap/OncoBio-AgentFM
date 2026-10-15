from dataclasses import replace
from datetime import date, datetime, timedelta

import pytest

from oncobio_agentfm.contracts import SurvivalRecord
from oncobio_agentfm.evaluation_support import (
    assess_horizon_support,
    assess_metric_support,
    training_censoring_survival,
)
from oncobio_agentfm.prognosis import ProtocolViolation
from oncobio_agentfm.research_protocol import SurvivalMetric


def outcome(patient_id: str, days: int, *, death: bool = True) -> SurvivalRecord:
    index = date(2025, 1, 1)
    end = index + timedelta(days=days)
    return SurvivalRecord(patient_id, index, end if death else None, None if death else end)


def supported_samples():
    training = (outcome("train-death", 300), outcome("train-late", 1400, death=False))
    evaluation = (
        outcome("eval-death", 400),
        *(outcome(f"eval-survivor-{number}", 1200, death=False) for number in range(10)),
    )
    return training, evaluation


def test_supported_horizon_reports_inputs_without_claiming_performance(protocol):
    training, evaluation = supported_samples()
    support = assess_horizon_support(training, evaluation, horizon_days=730, protocol=protocol)
    assert support.status == "SUPPORTED_FOR_METRIC_IMPLEMENTATION"
    assert support.training_censoring_survival == 1.0
    assert support.evaluation_deaths_by_horizon == 1
    assert support.evaluation_at_risk_before_horizon == 10
    assert support.evaluation_early_censorings == 0


def test_reverse_km_ties_remove_deaths_before_censoring_and_match_hand_calculation():
    training = (
        outcome("tied-death", 100),
        outcome("tied-censoring", 100, death=False),
        outcome("late", 1000, death=False),
    )
    # Reverse KM at day 100: 1 - 1 / (3 - 1) = 0.5, not 1 - 1 / 3.
    assert training_censoring_survival(training, 99) == 1.0
    assert training_censoring_survival(training, 100) == 0.5
    assert training_censoring_survival(tuple(reversed(training)), 730) == 0.5


def test_held_out_censoring_never_changes_training_estimated_weights(protocol):
    training, evaluation = supported_samples()
    first = assess_horizon_support(training, evaluation, horizon_days=730, protocol=protocol)
    changed = (evaluation[0], outcome("early-censoring", 100, death=False), *evaluation[2:])
    second = assess_horizon_support(training, changed, horizon_days=730, protocol=protocol)
    assert first.training_censoring_survival == second.training_censoring_survival
    assert second.evaluation_early_censorings == 1
    assert second.evaluation_at_risk_before_horizon == 9
    assert "evaluation_at_risk_below_threshold" in second.reasons


def test_extreme_training_censoring_withholds_horizon_instead_of_clipping_weights(protocol):
    training = (
        outcome("train-death", 50),
        *(outcome(f"censor-{number}", 100, death=False) for number in range(19)),
        outcome("late", 1400, death=False),
    )
    _, evaluation = supported_samples()
    support = assess_horizon_support(training, evaluation, horizon_days=730, protocol=protocol)
    assert support.training_censoring_survival == pytest.approx(0.05)
    assert "training_censoring_survival_below_threshold" in support.reasons
    assert support.status == "NOT_ESTIMABLE"


def test_empty_or_exhausted_training_support_does_not_extrapolate(protocol):
    _, evaluation = supported_samples()
    for training, reason in (
        ((), "training_support_empty"),
        ((outcome("train", 730),), "training_follow_up_support_exhausted"),
    ):
        support = assess_horizon_support(training, evaluation, horizon_days=730, protocol=protocol)
        assert support.training_censoring_survival is None
        assert reason in support.reasons


def test_patients_cannot_overlap_or_be_duplicated_in_metric_inputs(protocol):
    training, evaluation = supported_samples()
    with pytest.raises(ProtocolViolation, match="patient_overlap"):
        assess_horizon_support(training, (training[0],), horizon_days=730, protocol=protocol)
    with pytest.raises(ProtocolViolation, match="duplicate_patient"):
        assess_horizon_support(
            training, (evaluation[0], evaluation[0]), horizon_days=730, protocol=protocol
        )


def test_horizon_cannot_be_changed_to_make_support_look_better(protocol):
    training, evaluation = supported_samples()
    with pytest.raises(ProtocolViolation, match="unprespecified_evaluation_horizon"):
        assess_horizon_support(training, evaluation, horizon_days=500, protocol=protocol)
    with pytest.raises(ValueError, match="support threshold"):
        replace(protocol, minimum_evaluation_at_risk=1)


def test_empty_evaluation_and_zero_duration_are_explicit_failures(protocol):
    training, _ = supported_samples()
    empty = assess_horizon_support(training, (), horizon_days=730, protocol=protocol)
    assert empty.status == "NOT_ESTIMABLE"
    assert "evaluation_has_no_deaths_by_horizon" in empty.reasons
    with pytest.raises(ProtocolViolation, match="nonpositive_follow_up"):
        assess_horizon_support(training, (outcome("zero", 0),), horizon_days=730, protocol=protocol)


def test_death_exactly_at_tau_does_not_supply_a_primary_concordance_event(protocol):
    training, evaluation = supported_samples()
    boundary = (outcome("eval-boundary-death", 730), *evaluation[1:])
    support = assess_horizon_support(training, boundary, horizon_days=730, protocol=protocol)
    assert support.evaluation_deaths_by_horizon == 1
    assert support.evaluation_deaths_before_horizon == 0
    assert support.evaluation_at_risk_before_horizon == 11
    assert "evaluation_has_no_deaths_before_horizon" in support.reasons
    assert support.status == "NOT_ESTIMABLE"


def metric_support(training, evaluation, protocol, metric, **kwargs):

    return assess_metric_support(
        training,
        evaluation,
        metric=metric,
        horizon_days=None if metric is SurvivalMetric.HARRELL else 730,
        protocol=protocol,
        **kwargs,
    )


def test_one_event_and_many_pairs_do_not_establish_precision(protocol):

    training, evaluation = supported_samples()
    result = metric_support(training, evaluation, protocol, SurvivalMetric.UNO)
    assert result["status"] == "OUTCOME_PREREQUISITES_MET"
    assert result["details"]["comparable_pairs"] == 10
    assert result["details"]["distinct_comparable_event_patients"] == 1
    assert result["precision_status"] == "NOT_ESTABLISHED"
    assert result["prediction_outputs_validated"] is False
    assert result["model_evaluation"] == "NOT_RUN"


def test_auc_includes_cutoff_day_case_while_uno_excludes_its_event_contribution(protocol):

    training, evaluation = supported_samples()
    boundary = (outcome("boundary-event", 730), *evaluation[1:])
    uno = metric_support(training, boundary, protocol, SurvivalMetric.UNO)
    auc = metric_support(training, boundary, protocol, SurvivalMetric.AUC)
    assert uno["status"] == "NOT_ESTIMABLE"
    assert "no_comparable_event_pairs" in uno["fatal_reasons"]
    assert auc["status"] == "OUTCOME_PREREQUISITES_MET"
    assert auc["details"]["deaths_by_horizon"] == 1
    assert auc["details"]["observed_beyond_horizon"] == 10


def test_brier_and_ibs_can_have_no_early_deaths_but_auc_cannot(protocol):

    training, evaluation = supported_samples()
    survivors = evaluation[1:]
    for metric in (SurvivalMetric.BRIER, SurvivalMetric.IBS):
        result = metric_support(training, survivors, protocol, metric)
        assert result["status"] == "OUTCOME_PREREQUISITES_MET"
        assert not result["prediction_outputs_validated"]
        if metric is SurvivalMetric.IBS:
            assert result["details"]["grid"]["points"] == 730
            assert result["details"]["grid"]["first_day_without_known_contributions"] is None
    auc = metric_support(training, survivors, protocol, SurvivalMetric.AUC)
    assert "auc_has_no_cases_by_horizon" in auc["fatal_reasons"]


def test_ibs_does_not_crop_the_fixed_grid_to_avoid_late_censoring(protocol):

    training, _ = supported_samples()
    censored = tuple(outcome(f"censored-{i}", 100, death=False) for i in range(12))
    result = metric_support(training, censored, protocol, SurvivalMetric.IBS)
    assert result["status"] == "NOT_ESTIMABLE"
    assert result["details"]["grid"]["stop_day"] == 730
    assert result["details"]["grid"]["first_day_without_known_contributions"] == 100
    assert result["details"]["early_censorings"] == 12


def test_mathematical_support_and_prespecified_stability_are_different(protocol):

    training, evaluation = supported_samples()
    result = metric_support(training, evaluation[:2], protocol, SurvivalMetric.UNO)
    assert result["fatal_reasons"] == []
    assert result["status"] == "WITHHELD_BY_PROTOCOL"
    assert "evaluation_at_risk_below_threshold" in result["stability_reasons"]
    assert result["precision_status"] == "NOT_ESTABLISHED"


def test_low_positive_G_is_withheld_without_clipping_or_claiming_mathematical_impossibility(
    protocol,
):

    training = (
        outcome("train-death", 50),
        *(outcome(f"censor-{i}", 100, death=False) for i in range(19)),
        outcome("late", 1400, death=False),
    )
    _, evaluation = supported_samples()
    result = metric_support(training, evaluation, protocol, SurvivalMetric.AUC)
    assert result["details"]["training_G_at_horizon"] == pytest.approx(0.05)
    assert result["status"] == "WITHHELD_BY_PROTOCOL"
    assert result["fatal_reasons"] == []


def test_pooled_horizon_support_does_not_certify_calibration_tertiles(protocol):

    training, evaluation = supported_samples()
    missing = metric_support(training, evaluation, protocol, SurvivalMetric.CALIBRATION)
    assert missing["status"] == "NOT_ASSESSED"
    sparse = metric_support(
        training,
        evaluation,
        protocol,
        SurvivalMetric.CALIBRATION,
        training_risk=(0.0, 1.0),
        evaluation_risk=(0.1,) * len(evaluation),
    )
    assert sparse["details"]["at_risk_before_horizon"] == 10
    assert sparse["status"] == "NOT_ESTIMABLE"
    assert "calibration_group_2_empty" in sparse["fatal_reasons"]


def test_calibration_group_boundaries_are_fitted_from_training_risk_only(protocol):

    training, evaluation = supported_samples()
    first = metric_support(
        training,
        evaluation,
        protocol,
        SurvivalMetric.CALIBRATION,
        training_risk=(0.0, 1.0),
        evaluation_risk=(0.1,) * len(evaluation),
    )
    second = metric_support(
        training,
        evaluation,
        protocol,
        SurvivalMetric.CALIBRATION,
        training_risk=(0.0, 1.0),
        evaluation_risk=(0.9,) * len(evaluation),
    )
    assert (
        first["details"]["training_only_cutpoints"] == second["details"]["training_only_cutpoints"]
    )
    assert (
        first["details"]["fixed_training_risk_groups"]
        != second["details"]["fixed_training_risk_groups"]
    )


def test_tied_training_tertiles_are_not_replaced_by_test_selected_groups(protocol):

    training, evaluation = supported_samples()
    result = metric_support(
        training,
        evaluation,
        protocol,
        SurvivalMetric.CALIBRATION,
        training_risk=(0.0, 0.0),
        evaluation_risk=(0.0,) * len(evaluation),
    )
    assert "calibration_training_tertiles_tied" in result["fatal_reasons"]
    assert result["status"] == "NOT_ESTIMABLE"


def test_positive_calibration_group_support_still_has_no_precision_or_probability_validation(
    protocol,
):

    training = (
        outcome("train-event", 100),
        outcome("train-survivor", 1400, death=False),
        outcome("train-survivor-2", 1500, death=False),
    )
    records = []
    risks = []
    for group in range(3):
        records.append(outcome(f"group-{group}-death", 400))
        records.extend(outcome(f"group-{group}-survivor-{i}", 1200, death=False) for i in range(10))
        risks.extend([float(group)] * 11)
    result = metric_support(
        training,
        tuple(records),
        protocol,
        SurvivalMetric.CALIBRATION,
        training_risk=(0.0, 1.0, 2.0),
        evaluation_risk=tuple(risks),
    )
    assert result["status"] == "OUTCOME_PREREQUISITES_MET"
    assert result["precision_status"] == "NOT_ESTABLISHED"
    assert result["prediction_outputs_validated"] is False
    assert [
        g["at_risk_before_horizon"] for g in result["details"]["fixed_training_risk_groups"]
    ] == [10, 10, 10]


def test_harrell_descriptive_support_does_not_inherit_unrelated_ipcw_horizon_requirements(protocol):

    training = (outcome("train-death", 100), outcome("train-censor", 200, death=False))
    evaluation = (outcome("eval-event", 300), outcome("eval-censor", 900, death=False))
    result = metric_support(training, evaluation, protocol, SurvivalMetric.HARRELL)
    assert result["status"] == "OUTCOME_PREREQUISITES_MET"
    assert result["horizon_days"] is None
    assert result["precision_status"] == "NOT_ESTABLISHED"


def test_metric_support_rejects_mixed_date_precision_and_boolean_risks(protocol):

    training, evaluation = supported_samples()
    timestamp = SurvivalRecord("timestamp", datetime(2025, 1, 1), datetime(2025, 4, 1))
    with pytest.raises(ProtocolViolation, match="date_precision_invalid"):
        metric_support(training, (timestamp,), protocol, SurvivalMetric.UNO)
    with pytest.raises(ProtocolViolation, match="risk_scores_invalid"):
        metric_support(
            training,
            evaluation,
            protocol,
            SurvivalMetric.CALIBRATION,
            training_risk=(True, 1.0),
            evaluation_risk=(0.1,) * len(evaluation),
        )
