from __future__ import annotations

from datetime import date, timedelta

import pytest

from oncobio_agentfm.clinical import (
    ClinicalStage,
    ECOGStatus,
    NSCLCBaseline,
    NSCLCOverallSurvivalRecord,
    StageGroup,
)
from oncobio_agentfm.contracts import EvidenceRecord, MissingnessReason, SurvivalRecord

INDEX_DATE = date(2025, 1, 1)


@pytest.fixture
def make_baseline():
    def create(
        patient_id: str = "P001",
        *,
        age_years: int = 65,
        stage: StageGroup | None = StageGroup.TWO,
        stage_version: int | None = 9,
        ecog: int | None = 1,
        stage_acquired_offset: int = -10,
        stage_available_offset: int = -5,
        ecog_available_offset: int = -1,
    ) -> NSCLCBaseline:
        def present(name: str, acquired_offset: int, available_offset: int) -> EvidenceRecord:
            return EvidenceRecord(
                patient_id=patient_id,
                source_id=f"{name}:{patient_id}",
                acquired_on=INDEX_DATE + timedelta(days=acquired_offset),
                available_on=INDEX_DATE + timedelta(days=available_offset),
            )

        def missing() -> EvidenceRecord:
            return EvidenceRecord(
                patient_id=patient_id,
                source_id=None,
                acquired_on=None,
                available_on=None,
                missing_reason=MissingnessReason.NOT_RECORDED,
            )

        return NSCLCBaseline(
            patient_id=patient_id,
            index_date=INDEX_DATE,
            diagnosis=present("histology", -20, -15),
            age_years=age_years,
            age_evidence=present("age", -2, -2),
            stage=ClinicalStage(
                evidence=present("stage", stage_acquired_offset, stage_available_offset)
                if stage is not None
                else missing(),
                group=stage,
                ajcc_version=stage_version if stage is not None else None,
            ),
            ecog=ECOGStatus(
                evidence=present("ecog", -2, ecog_available_offset)
                if ecog is not None
                else missing(),
                score=ecog,
            ),
        )

    return create


@pytest.fixture
def make_record(make_baseline):
    def create(
        patient_id: str = "P001",
        *,
        event: bool = True,
        follow_up_days: int = 100,
        **baseline_kwargs,
    ) -> NSCLCOverallSurvivalRecord:
        baseline = make_baseline(patient_id, **baseline_kwargs)
        end_date = baseline.index_date + timedelta(days=follow_up_days)
        outcome = SurvivalRecord(
            patient_id=patient_id,
            index_date=baseline.index_date,
            event_date=end_date if event else None,
            censor_date=None if event else end_date,
        )
        return NSCLCOverallSurvivalRecord(baseline=baseline, outcome=outcome)

    return create
