"""Full-risk-set, Regularized Cox Baseline for index-date Clinical Evidence."""

from __future__ import annotations

from dataclasses import dataclass
from math import fsum, isfinite

import torch
from torch import Tensor

from oncobio_agentfm.clinical import NSCLCBaseline, NSCLCCohort, StageGroup
from oncobio_agentfm.splitting import CohortSplit
from oncobio_agentfm.survival import CoxPHLinear, efron_negative_partial_log_likelihood

FEATURE_NAMES = (
    "age_per_decade",
    "stage_I",
    "stage_II",
    "stage_III",
    "stage_IV",
    "stage_missing",
    "ecog_0",
    "ecog_1",
    "ecog_2",
    "ecog_3",
    "ecog_4",
    "ecog_missing",
)


@dataclass(frozen=True, slots=True)
class ClinicalEncoder:
    """Fixed clinical categories with age centered using training patients only."""

    training_age_mean: float
    ajcc_version: int

    def __post_init__(self) -> None:
        if not isfinite(self.training_age_mean):
            raise ValueError("training_age_mean must be finite")
        if type(self.ajcc_version) is not int or self.ajcc_version not in (8, 9):
            raise ValueError("encoder AJCC version must be 8 or 9")

    @classmethod
    def fit(cls, records: tuple[NSCLCBaseline, ...], *, ajcc_version: int) -> ClinicalEncoder:
        """Estimating the sole data-dependent feature statistic from training records."""

        if not records:
            raise ValueError("encoder requires at least one training patient")
        if type(ajcc_version) is not int or ajcc_version not in (8, 9):
            raise ValueError("encoder AJCC version must be 8 or 9")
        if any(
            record.stage.ajcc_version != ajcc_version
            for record in records
            if record.stage.group is not None
        ):
            raise ValueError("training records contain a different clinical stage version")
        return cls(
            training_age_mean=fsum(sorted(record.age_years for record in records)) / len(records),
            ajcc_version=ajcc_version,
        )

    def transform(self, records: tuple[NSCLCBaseline, ...]) -> Tensor:
        """Producing float64 rows in the declared, stable feature order."""

        rows: list[list[float]] = []
        for record in records:
            if record.stage.group is not None and record.stage.ajcc_version != self.ajcc_version:
                raise ValueError("clinical stage version differs from the fitted encoder")
            stage = record.stage.group
            ecog = record.ecog.score
            rows.append(
                [
                    (record.age_years - self.training_age_mean) / 10.0,
                    *(float(stage == group) for group in StageGroup),
                    float(stage is None),
                    *(float(ecog == score) for score in range(5)),
                    float(ecog is None),
                ]
            )
        if not rows:
            return torch.empty((0, len(FEATURE_NAMES)), dtype=torch.float64)
        return torch.tensor(rows, dtype=torch.float64)


@dataclass(frozen=True, slots=True)
class FittedClinicalCox:
    """Frozen coefficients and the training-only encoder used to obtain them."""

    encoder: ClinicalEncoder
    coefficients: tuple[float, ...]
    train_patients: int
    train_deaths: int
    ridge_strength: float
    training_objective: float
    training_gradient_max: float
    training_stage_levels: tuple[str, ...]
    training_ecog_levels: tuple[str, ...]
    training_age_range: tuple[int, int]

    def __post_init__(self) -> None:
        if len(self.coefficients) != len(FEATURE_NAMES) or not all(
            isfinite(value) for value in self.coefficients
        ):
            raise ValueError("fitted coefficients must be finite and match the clinical design")
        if self.train_patients < 2 or not 0 < self.train_deaths <= self.train_patients:
            raise ValueError("fitted model requires a valid training event count")
        if not isfinite(self.ridge_strength) or self.ridge_strength <= 0:
            raise ValueError("fitted ridge_strength must be finite and positive")
        if not isfinite(self.training_objective):
            raise ValueError("fitted training objective must be finite")
        if not isfinite(self.training_gradient_max) or self.training_gradient_max > 1e-5:
            raise ValueError("fitted model requires a converged training gradient")
        if not self.training_stage_levels or not self.training_ecog_levels:
            raise ValueError("fitted model requires observed training category support")
        if not 18 <= self.training_age_range[0] <= self.training_age_range[1] <= 120:
            raise ValueError("fitted model requires a valid training age range")

    def coefficient_table(self) -> tuple[tuple[str, float], ...]:
        """Returning labeled log-hazard coefficients in the fixed design order."""

        return tuple(zip(FEATURE_NAMES, self.coefficients, strict=True))

    def predict_log_risk(self, record: NSCLCBaseline) -> float:
        """Scoring baseline evidence without requiring any known outcome."""

        features = self.encoder.transform((record,))[0].tolist()
        stage_level = record.stage.group.value if record.stage.group is not None else "missing"
        ecog_level = str(record.ecog.score) if record.ecog.score is not None else "missing"
        if stage_level not in self.training_stage_levels:
            raise ValueError("stage category was not observed in training")
        if ecog_level not in self.training_ecog_levels:
            raise ValueError("ECOG category was not observed in training")
        if not self.training_age_range[0] <= record.age_years <= self.training_age_range[1]:
            raise ValueError("age is outside the observed training range")
        return fsum(
            weight * value for weight, value in zip(self.coefficients, features, strict=True)
        )


def fit_clinical_cox(
    cohort: NSCLCCohort,
    split: CohortSplit,
    *,
    ridge_strength: float = 0.01,
    max_iterations: int = 128,
) -> FittedClinicalCox:
    """Fitting Efron Cox loss on all training patients and their complete risk sets."""

    if split.all_patient_ids != cohort.patient_ids:
        raise ValueError("split must cover exactly the cohort patient identifiers")
    if len(split.train) < 2:
        raise ValueError("training partition requires at least two patients")
    if isinstance(ridge_strength, bool) or not isinstance(ridge_strength, (int, float)):
        raise ValueError("ridge_strength must be a finite positive number")
    if not isfinite(ridge_strength) or ridge_strength <= 0:
        raise ValueError("ridge_strength must be a finite positive number")
    if type(max_iterations) is not int or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")

    by_id = {record.patient_id: record for record in cohort.records}
    train_records = tuple(by_id[patient_id] for patient_id in sorted(split.train))
    train_baselines = tuple(record.baseline for record in train_records)
    train_deaths = sum(record.outcome.event_observed for record in train_records)
    if train_deaths == 0:
        raise ValueError("training partition requires at least one observed death")
    longest_follow_up = max(record.outcome.duration_days for record in train_records)
    if not any(
        record.outcome.event_observed and record.outcome.duration_days < longest_follow_up
        for record in train_records
    ):
        raise ValueError("training partition has no comparable event and later follow-up")

    encoder = ClinicalEncoder.fit(train_baselines, ajcc_version=cohort.ajcc_version)
    features = encoder.transform(train_baselines)
    durations = torch.tensor(
        [record.outcome.duration_days for record in train_records], dtype=torch.float64
    )
    events = torch.tensor(
        [record.outcome.event_observed for record in train_records], dtype=torch.int64
    )

    with torch.random.fork_rng(devices=[]):
        model = CoxPHLinear(in_features=len(FEATURE_NAMES)).to(dtype=torch.float64)
        with torch.no_grad():
            model.linear.weight.zero_()

    optimizer = torch.optim.LBFGS(
        model.parameters(),
        lr=1.0,
        max_iter=max_iterations,
        line_search_fn="strong_wolfe",
        tolerance_grad=1e-10,
        tolerance_change=1e-12,
    )

    def objective() -> Tensor:
        loss = efron_negative_partial_log_likelihood(model(features), durations, events)
        return loss + 0.5 * ridge_strength * model.linear.weight.square().sum()

    def closure() -> Tensor:
        optimizer.zero_grad()
        loss = objective()
        if not torch.isfinite(loss):
            raise FloatingPointError("non-finite clinical Cox training objective")
        loss.backward()
        return loss

    optimizer.step(closure)
    final_loss = objective()
    final_gradient_max = torch.autograd.grad(final_loss, model.linear.weight)[0].abs().max().item()
    with torch.no_grad():
        final_objective = final_loss.item()
        coefficients = tuple(model.linear.weight.squeeze(0).tolist())
    if not isfinite(final_objective) or not all(isfinite(value) for value in coefficients):
        raise FloatingPointError("clinical Cox fitting produced non-finite parameters")
    if not isfinite(final_gradient_max) or final_gradient_max > 1e-5:
        raise RuntimeError("clinical Cox optimization did not converge")

    return FittedClinicalCox(
        encoder=encoder,
        coefficients=coefficients,
        train_patients=len(train_records),
        train_deaths=train_deaths,
        ridge_strength=float(ridge_strength),
        training_objective=final_objective,
        training_gradient_max=final_gradient_max,
        training_stage_levels=tuple(
            sorted(
                {
                    record.stage.group.value if record.stage.group is not None else "missing"
                    for record in train_baselines
                }
            )
        ),
        training_ecog_levels=tuple(
            sorted(
                {
                    str(record.ecog.score) if record.ecog.score is not None else "missing"
                    for record in train_baselines
                }
            )
        ),
        training_age_range=(
            min(record.age_years for record in train_baselines),
            max(record.age_years for record in train_baselines),
        ),
    )
