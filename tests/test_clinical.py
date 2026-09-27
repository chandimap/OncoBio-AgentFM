from dataclasses import replace
from datetime import date

import pytest

from oncobio_agentfm.clinical import (
    ClinicalStage,
    ECOGStatus,
    NSCLCBaseline,
    NSCLCCohort,
    NSCLCOverallSurvivalRecord,
    StageGroup,
)
from oncobio_agentfm.contracts import EvidenceRecord, MissingnessReason


def test_audit_counts_deaths_censoring_and_missingness(make_record) -> None:
    cohort = NSCLCCohort(
        records=(
            make_record("patient-a-001", stage=StageGroup.ONE, ecog=0),
            make_record("patient-b-002", stage=None, ecog=2, event=False),
            make_record("patient-c-003", stage=StageGroup.FOUR, ecog=None),
        ),
        ajcc_version=9,
    )

    audit = cohort.audit()
    assert (audit.patients, audit.deaths, audit.censored) == (3, 2, 1)
    assert audit.stage_counts == (("I", 1), ("IV", 1), ("missing", 1))
    assert audit.ecog_counts == (("0", 1), ("2", 1), ("missing", 1))
    assert audit.stage_missing_reasons == (("not_recorded", 1),)
    assert audit.ecog_missing_reasons == (("not_recorded", 1),)
    assert all(
        patient_id not in repr(audit)
        for patient_id in ("patient-a-001", "patient-b-002", "patient-c-003")
    )


@pytest.mark.parametrize("age", [True, 17, 121, 65.0])
def test_age_at_index_rejects_invalid_values(make_record, age) -> None:
    with pytest.raises(ValueError, match="age_years"):
        make_record(age_years=age)


@pytest.mark.parametrize("score", [True, -1, 5, 1.0])
def test_ecog_rejects_invalid_or_postmortem_scores(make_record, score) -> None:
    source = make_record().baseline.ecog.evidence
    with pytest.raises(ValueError, match="0 through 4"):
        ECOGStatus(evidence=source, score=score)


def test_stage_and_ecog_reject_present_value_with_missing_evidence(make_record) -> None:
    missing = make_record(stage=None).baseline.stage.evidence
    with pytest.raises(ValueError, match="presence must agree"):
        ClinicalStage(evidence=missing, group=StageGroup.TWO, ajcc_version=9)
    with pytest.raises(ValueError, match="presence must agree"):
        ECOGStatus(evidence=missing, score=1)


def test_stage_requires_explicit_version_and_enum(make_record) -> None:
    present = make_record().baseline.stage.evidence
    with pytest.raises(ValueError, match="AJCC version"):
        ClinicalStage(evidence=present, group=StageGroup.TWO, ajcc_version=None)
    with pytest.raises(TypeError, match="StageGroup"):
        ClinicalStage(evidence=present, group="II", ajcc_version=9)


def test_missingness_reason_must_have_declared_vocabulary(make_record) -> None:
    original = make_record(stage=None).baseline.stage.evidence
    incorrect = replace(original, missing_reason="not_recorded")
    with pytest.raises(TypeError, match="MissingnessReason"):
        ClinicalStage(evidence=incorrect, group=None, ajcc_version=None)


def test_index_locks_acquisition_and_availability(make_record) -> None:
    with pytest.raises(ValueError, match="acquired after"):
        make_record(stage_acquired_offset=1, stage_available_offset=1)
    with pytest.raises(ValueError, match="information leakage"):
        make_record(stage_acquired_offset=-10, stage_available_offset=1)
    with pytest.raises(ValueError, match="information leakage"):
        make_record(ecog_available_offset=1)


def test_all_evidence_must_match_baseline_patient(make_record) -> None:
    record = make_record()
    foreign = replace(record.baseline.stage.evidence, patient_id="OTHER")
    with pytest.raises(ValueError, match="same patient"):
        replace(record.baseline, stage=replace(record.baseline.stage, evidence=foreign))


def test_diagnosis_and_age_require_present_pre_index_sources(make_record) -> None:
    record = make_record()
    absent = EvidenceRecord(
        patient_id=record.patient_id,
        source_id=None,
        acquired_on=None,
        available_on=None,
        missing_reason=MissingnessReason.UNKNOWN,
    )
    with pytest.raises(ValueError, match="diagnosis and age"):
        replace(record.baseline, diagnosis=absent)
    with pytest.raises(ValueError, match="diagnosis and age"):
        replace(record.baseline, age_evidence=absent)
    late = replace(record.baseline.diagnosis, available_on=date(2025, 1, 2))
    with pytest.raises(ValueError, match="information leakage"):
        replace(record.baseline, diagnosis=late)


def test_baseline_exists_without_outcome_and_join_requires_same_index(
    make_baseline, make_record
) -> None:
    baseline: NSCLCBaseline = make_baseline("A")
    assert baseline.patient_id == "A"
    assert not hasattr(baseline, "outcome")

    outcome = make_record("A").outcome
    with pytest.raises(ValueError, match="same patient"):
        NSCLCOverallSurvivalRecord(baseline=baseline, outcome=replace(outcome, patient_id="B"))
    with pytest.raises(ValueError, match="same index date"):
        NSCLCOverallSurvivalRecord(
            baseline=baseline, outcome=replace(outcome, index_date=date(2025, 1, 2))
        )


def test_cohort_rejects_duplicates_and_mixed_stage_versions(make_record) -> None:
    first = make_record("A")
    with pytest.raises(ValueError, match="unique"):
        NSCLCCohort(records=(first, first), ajcc_version=9)
    with pytest.raises(ValueError, match="AJCC version"):
        NSCLCCohort(records=(first, make_record("B", stage_version=8)), ajcc_version=9)
    with pytest.raises(TypeError, match="tuple"):
        NSCLCCohort(records=[first], ajcc_version=9)


def test_valid_records_remain_immutable(make_record) -> None:
    original: NSCLCOverallSurvivalRecord = make_record()
    with pytest.raises(AttributeError):
        original.baseline.age_years = 70
