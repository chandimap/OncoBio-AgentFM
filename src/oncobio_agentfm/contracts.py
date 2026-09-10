"""Patient-level temporal and survival contracts for pretreatment research."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from datetime import date
from enum import StrEnum


class MissingnessReason(StrEnum):
    """Explicit reasons why a pretreatment evidence item is absent."""

    NOT_ACQUIRED = "not_acquired"
    NOT_RECORDED = "not_recorded"
    UNAVAILABLE = "unavailable"
    UNUSABLE = "unusable"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class SurvivalRecord:
    """Right-censored survival outcome anchored to a patient-specific index date."""

    patient_id: str
    index_date: date
    event_date: date | None = None
    censor_date: date | None = None

    def __post_init__(self) -> None:
        if not self.patient_id.strip():
            raise ValueError("patient_id must be non-empty")
        if self.patient_id != self.patient_id.strip():
            raise ValueError("patient_id must not contain surrounding whitespace")
        if (self.event_date is None) == (self.censor_date is None):
            raise ValueError("exactly one of event_date or censor_date must be provided")
        if self.end_date < self.index_date:
            raise ValueError("survival end date cannot precede the index date")

    @property
    def event_observed(self) -> bool:
        """Returning whether the survival endpoint was observed rather than censored."""

        return self.event_date is not None

    @property
    def end_date(self) -> date:
        """Returning the observed event date or right-censoring date."""

        endpoint = self.event_date if self.event_date is not None else self.censor_date
        if endpoint is None:  # Defensive guard; constructor validation makes this unreachable.
            raise RuntimeError("survival endpoint is missing")
        return endpoint

    @property
    def duration_days(self) -> int:
        """Returning non-negative follow-up duration in whole days from the index date."""

        return (self.end_date - self.index_date).days


@dataclass(frozen=True, slots=True)
class EvidenceRecord:
    """One patient-linked evidence item with acquisition and availability timing."""

    patient_id: str
    source_id: str | None
    acquired_on: date | None
    available_on: date | None
    missing_reason: MissingnessReason | None = None

    def __post_init__(self) -> None:
        if not self.patient_id.strip():
            raise ValueError("patient_id must be non-empty")
        if self.patient_id != self.patient_id.strip():
            raise ValueError("patient_id must not contain surrounding whitespace")

        present = self.source_id is not None
        if present:
            if not self.source_id.strip():
                raise ValueError("source_id must be non-empty when evidence is present")
            if self.source_id != self.source_id.strip():
                raise ValueError("source_id must not contain surrounding whitespace")
            if self.acquired_on is None or self.available_on is None:
                raise ValueError("present evidence requires acquisition and availability dates")
            if self.missing_reason is not None:
                raise ValueError("present evidence cannot also have a missingness reason")
            if self.available_on < self.acquired_on:
                raise ValueError("availability date cannot precede acquisition date")
        else:
            if self.acquired_on is not None or self.available_on is not None:
                raise ValueError("missing evidence cannot carry acquisition or availability dates")
            if self.missing_reason is None:
                raise ValueError("missing evidence requires an explicit missingness reason")

    @property
    def is_present(self) -> bool:
        """Returning whether this evidence item is available as a concrete source."""

        return self.source_id is not None

    def validate_for_index(self, index_date: date) -> None:
        """Rejecting evidence that would leak post-index information into pretreatment inputs."""

        if not self.is_present:
            return
        if self.acquired_on is None or self.available_on is None:
            raise RuntimeError("present evidence is missing validated timing metadata")
        if self.acquired_on > index_date:
            raise ValueError("evidence acquired after the index date is not pretreatment evidence")
        if self.available_on > index_date:
            raise ValueError("evidence unavailable at the index date creates information leakage")


def validate_patient_bundle(outcome: SurvivalRecord, evidence: Iterable[EvidenceRecord]) -> None:
    """Validating patient identity and pretreatment timing across a patient evidence bundle."""

    for item in evidence:
        if item.patient_id != outcome.patient_id:
            raise ValueError("all evidence must belong to the same patient as the outcome")
        item.validate_for_index(outcome.index_date)
