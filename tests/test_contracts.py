from datetime import date

import pytest

from oncobio_agentfm.contracts import (
    EvidenceRecord,
    MissingnessReason,
    SurvivalRecord,
    validate_patient_bundle,
)


def test_survival_record_encodes_right_censoring() -> None:
    record = SurvivalRecord(
        patient_id="P001",
        index_date=date(2025, 1, 1),
        censor_date=date(2025, 4, 1),
    )

    assert record.event_observed is False
    assert record.duration_days == 90


def test_survival_record_requires_exactly_one_endpoint() -> None:
    with pytest.raises(ValueError, match="exactly one"):
        SurvivalRecord(patient_id="P001", index_date=date(2025, 1, 1))


def test_present_evidence_rejects_information_available_after_index() -> None:
    evidence = EvidenceRecord(
        patient_id="P001",
        source_id="SRC-1",
        acquired_on=date(2024, 12, 20),
        available_on=date(2025, 1, 2),
    )

    with pytest.raises(ValueError, match="information leakage"):
        evidence.validate_for_index(date(2025, 1, 1))


def test_missing_evidence_requires_explicit_reason() -> None:
    evidence = EvidenceRecord(
        patient_id="P001",
        source_id=None,
        acquired_on=None,
        available_on=None,
        missing_reason=MissingnessReason.NOT_ACQUIRED,
    )

    evidence.validate_for_index(date(2025, 1, 1))
    assert evidence.is_present is False


def test_patient_bundle_rejects_cross_patient_evidence() -> None:
    outcome = SurvivalRecord(
        patient_id="P001",
        index_date=date(2025, 1, 1),
        event_date=date(2025, 6, 1),
    )
    evidence = EvidenceRecord(
        patient_id="P002",
        source_id="SRC-2",
        acquired_on=date(2024, 12, 1),
        available_on=date(2024, 12, 2),
    )

    with pytest.raises(ValueError, match="same patient"):
        validate_patient_bundle(outcome, [evidence])


def test_survival_record_rejects_endpoint_before_index() -> None:
    with pytest.raises(ValueError, match="cannot precede"):
        SurvivalRecord(
            patient_id="P001",
            index_date=date(2025, 1, 2),
            event_date=date(2025, 1, 1),
        )


def test_present_evidence_rejects_acquisition_after_index() -> None:
    evidence = EvidenceRecord(
        patient_id="P001",
        source_id="SRC-3",
        acquired_on=date(2025, 1, 2),
        available_on=date(2025, 1, 2),
    )

    with pytest.raises(ValueError, match="acquired after"):
        evidence.validate_for_index(date(2025, 1, 1))


def test_patient_id_rejects_surrounding_whitespace() -> None:
    with pytest.raises(ValueError, match="surrounding whitespace"):
        SurvivalRecord(
            patient_id=" P001",
            index_date=date(2025, 1, 1),
            censor_date=date(2025, 2, 1),
        )
