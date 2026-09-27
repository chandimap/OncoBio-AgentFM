"""Index-date NSCLC Records and Aggregate Cohort Checks."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from datetime import date
from enum import StrEnum

from oncobio_agentfm.contracts import (
    EvidenceRecord,
    MissingnessReason,
    SurvivalRecord,
)


class StageGroup(StrEnum):
    """Clinical NSCLC stage group, without an inferred substage."""

    ONE = "I"
    TWO = "II"
    THREE = "III"
    FOUR = "IV"


@dataclass(frozen=True, slots=True)
class ClinicalStage:
    """An observed clinical stage or an explicitly missing stage."""

    evidence: EvidenceRecord
    group: StageGroup | None
    ajcc_version: int | None

    def __post_init__(self) -> None:
        if not isinstance(self.evidence, EvidenceRecord):
            raise TypeError("stage evidence must be an EvidenceRecord")
        if self.evidence.is_present != (self.group is not None):
            raise ValueError("stage value and evidence presence must agree")
        if self.group is None:
            if self.ajcc_version is not None:
                raise ValueError("missing stage cannot declare a staging version")
            if not isinstance(self.evidence.missing_reason, MissingnessReason):
                raise TypeError("missing stage requires a MissingnessReason")
            return
        if not isinstance(self.group, StageGroup):
            raise TypeError("group must be a StageGroup")
        if type(self.ajcc_version) is not int or self.ajcc_version not in (8, 9):
            raise ValueError("observed stage requires AJCC version 8 or 9")


@dataclass(frozen=True, slots=True)
class ECOGStatus:
    """An observed ECOG performance status or an explicitly missing value."""

    evidence: EvidenceRecord
    score: int | None

    def __post_init__(self) -> None:
        if not isinstance(self.evidence, EvidenceRecord):
            raise TypeError("ECOG evidence must be an EvidenceRecord")
        if self.evidence.is_present != (self.score is not None):
            raise ValueError("ECOG value and evidence presence must agree")
        if self.score is None and not isinstance(self.evidence.missing_reason, MissingnessReason):
            raise TypeError("missing ECOG requires a MissingnessReason")
        if self.score is not None and (type(self.score) is not int or self.score not in range(5)):
            raise ValueError("baseline ECOG score must be an integer from 0 through 4")


@dataclass(frozen=True, slots=True)
class NSCLCBaseline:
    """Index-date NSCLC evidence that can exist before any outcome is known."""

    patient_id: str
    index_date: date
    diagnosis: EvidenceRecord
    age_years: int
    age_evidence: EvidenceRecord
    stage: ClinicalStage
    ecog: ECOGStatus

    def __post_init__(self) -> None:
        if not isinstance(self.patient_id, str) or not self.patient_id.strip():
            raise ValueError("baseline patient_id must be non-empty")
        if self.patient_id != self.patient_id.strip():
            raise ValueError("baseline patient_id must not contain surrounding whitespace")
        if type(self.index_date) is not date:
            raise TypeError("index_date must be a calendar date")
        if not isinstance(self.diagnosis, EvidenceRecord) or not isinstance(
            self.age_evidence, EvidenceRecord
        ):
            raise TypeError("diagnosis and age evidence must be EvidenceRecords")
        if not isinstance(self.stage, ClinicalStage) or not isinstance(self.ecog, ECOGStatus):
            raise TypeError("stage and ECOG must be validated clinical observations")
        if type(self.age_years) is not int or not 18 <= self.age_years <= 120:
            raise ValueError("age_years must be an integer from 18 through 120 at the index date")
        if not self.diagnosis.is_present or not self.age_evidence.is_present:
            raise ValueError("NSCLC diagnosis and age require present source evidence")
        for evidence in (
            self.diagnosis,
            self.age_evidence,
            self.stage.evidence,
            self.ecog.evidence,
        ):
            if evidence.patient_id != self.patient_id:
                raise ValueError("all baseline evidence must belong to the same patient")
            evidence.validate_for_index(self.index_date)


@dataclass(frozen=True, slots=True)
class NSCLCOverallSurvivalRecord:
    """Joining a separately validated NSCLC baseline to all-cause death follow-up."""

    baseline: NSCLCBaseline
    outcome: SurvivalRecord

    def __post_init__(self) -> None:
        if not isinstance(self.baseline, NSCLCBaseline) or not isinstance(
            self.outcome, SurvivalRecord
        ):
            raise TypeError("baseline and outcome must be validated records")
        if self.baseline.patient_id != self.outcome.patient_id:
            raise ValueError("baseline and outcome must belong to the same patient")
        if self.baseline.index_date != self.outcome.index_date:
            raise ValueError("baseline and outcome must share the same index date")

    @property
    def patient_id(self) -> str:
        """Returning the identity shared by the outcome and all baseline evidence."""

        return self.baseline.patient_id


@dataclass(frozen=True, slots=True)
class ClinicalAudit:
    """Aggregating endpoint and evidence counts without patient identifiers."""

    patients: int
    deaths: int
    censored: int
    ajcc_version: int
    stage_counts: tuple[tuple[str, int], ...]
    ecog_counts: tuple[tuple[str, int], ...]
    stage_missing_reasons: tuple[tuple[str, int], ...]
    ecog_missing_reasons: tuple[tuple[str, int], ...]


def _missing_reason_value(evidence: EvidenceRecord) -> str:
    """Reading the declared reason for a validated, missing clinical observation."""

    reason = evidence.missing_reason
    if not isinstance(reason, MissingnessReason):
        raise RuntimeError("missing clinical evidence requires a validated missingness reason")
    return reason.value


@dataclass(frozen=True, slots=True)
class NSCLCCohort:
    """Unique, index-valid records sharing one clinical staging version."""

    records: tuple[NSCLCOverallSurvivalRecord, ...]
    ajcc_version: int

    def __post_init__(self) -> None:
        if type(self.ajcc_version) is not int or self.ajcc_version not in (8, 9):
            raise ValueError("cohort AJCC version must be 8 or 9")
        if type(self.records) is not tuple or any(
            not isinstance(record, NSCLCOverallSurvivalRecord) for record in self.records
        ):
            raise TypeError("cohort records must be a tuple of validated NSCLC records")
        if not self.records:
            raise ValueError("cohort must contain at least one patient")
        patient_ids = [record.patient_id for record in self.records]
        if len(set(patient_ids)) != len(patient_ids):
            raise ValueError("cohort patient identifiers must be unique")
        if any(
            record.baseline.stage.ajcc_version != self.ajcc_version
            for record in self.records
            if record.baseline.stage.group is not None
        ):
            raise ValueError("observed clinical stages must use the cohort AJCC version")

    @property
    def patient_ids(self) -> frozenset[str]:
        """Returning exactly the patients represented by this cohort."""

        return frozenset(record.patient_id for record in self.records)

    def audit(self) -> ClinicalAudit:
        """Counting endpoints, observed categories, and explicit missingness reasons."""

        stages = Counter(
            record.baseline.stage.group.value
            if record.baseline.stage.group is not None
            else "missing"
            for record in self.records
        )
        ecog = Counter(
            str(record.baseline.ecog.score) if record.baseline.ecog.score is not None else "missing"
            for record in self.records
        )
        stage_reasons = Counter(
            _missing_reason_value(record.baseline.stage.evidence)
            for record in self.records
            if record.baseline.stage.group is None
        )
        ecog_reasons = Counter(
            _missing_reason_value(record.baseline.ecog.evidence)
            for record in self.records
            if record.baseline.ecog.score is None
        )
        deaths = sum(record.outcome.event_observed for record in self.records)
        return ClinicalAudit(
            patients=len(self.records),
            deaths=deaths,
            censored=len(self.records) - deaths,
            ajcc_version=self.ajcc_version,
            stage_counts=tuple(sorted(stages.items())),
            ecog_counts=tuple(sorted(ecog.items())),
            stage_missing_reasons=tuple(sorted(stage_reasons.items())),
            ecog_missing_reasons=tuple(sorted(ecog_reasons.items())),
        )
