"""Deterministic patient-level cohort partitioning."""

from __future__ import annotations

import hashlib
from collections.abc import Iterable
from dataclasses import dataclass
from math import isclose


@dataclass(frozen=True, slots=True)
class CohortSplit:
    """Disjointing patient identifiers assigned to train, validation, and test partitions."""

    train: tuple[str, ...]
    validation: tuple[str, ...]
    test: tuple[str, ...]

    def __post_init__(self) -> None:
        train = set(self.train)
        validation = set(self.validation)
        test = set(self.test)
        duplicate_within_partition = (
            len(train) != len(self.train)
            or len(validation) != len(self.validation)
            or len(test) != len(self.test)
        )
        if duplicate_within_partition:
            raise ValueError("patient identifiers must be unique within each partition")
        if train & validation or train & test or validation & test:
            raise ValueError("patient partitions must be mutually disjoint")

    @property
    def all_patient_ids(self) -> frozenset[str]:
        """Return the complete set of assigned patient identifiers."""

        return frozenset((*self.train, *self.validation, *self.test))


def _stable_unit_interval(patient_id: str, seed: str) -> float:
    digest = hashlib.sha256(f"{seed}\x1f{patient_id}".encode()).digest()
    return int.from_bytes(digest, byteorder="big") / (1 << (8 * len(digest)))


def split_patient_ids(
    patient_ids: Iterable[str],
    *,
    seed: str = "oncobio-agentfm-v1",
    train_fraction: float = 0.70,
    validation_fraction: float = 0.15,
    test_fraction: float = 0.15,
) -> CohortSplit:
    """Assign unique patients deterministically without row-order dependence."""

    ids = tuple(patient_ids)
    if not ids:
        raise ValueError("patient_ids must not be empty")
    if any(not patient_id.strip() for patient_id in ids):
        raise ValueError("patient identifiers must be non-empty")
    if any(patient_id != patient_id.strip() for patient_id in ids):
        raise ValueError("patient identifiers must not contain surrounding whitespace")
    if len(ids) != len(set(ids)):
        raise ValueError("patient_ids must be unique before splitting")

    fractions = (train_fraction, validation_fraction, test_fraction)
    if any(fraction <= 0.0 or fraction >= 1.0 for fraction in fractions):
        raise ValueError("each split fraction must be strictly between 0 and 1")
    if not isclose(sum(fractions), 1.0, rel_tol=0.0, abs_tol=1e-12):
        raise ValueError("split fractions must sum to 1")

    train_cutoff = train_fraction
    validation_cutoff = train_fraction + validation_fraction
    train: list[str] = []
    validation: list[str] = []
    test: list[str] = []

    for patient_id in sorted(ids):
        value = _stable_unit_interval(patient_id, seed)
        if value < train_cutoff:
            train.append(patient_id)
        elif value < validation_cutoff:
            validation.append(patient_id)
        else:
            test.append(patient_id)

    return CohortSplit(tuple(train), tuple(validation), tuple(test))
