"""Training-only censoring support and prespecified horizon feasibility, not performance."""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
from datetime import date
from math import isfinite, prod
from statistics import quantiles

from oncobio_agentfm.contracts import SurvivalRecord
from oncobio_agentfm.prognosis import ProtocolViolation
from oncobio_agentfm.research_protocol import ResearchProtocol, SurvivalMetric


def _records(records: tuple[SurvivalRecord, ...]) -> None:
    """Require unique, positive-duration patient outcomes at one selected time origin."""

    if type(records) is not tuple or any(not isinstance(row, SurvivalRecord) for row in records):
        raise ProtocolViolation("support_records_invalid")
    identifiers = [row.patient_id for row in records]
    if len(identifiers) != len(set(identifiers)):
        raise ProtocolViolation("support_duplicate_patient")
    if any(type(row.index_date) is not date or type(row.end_date) is not date for row in records):
        raise ProtocolViolation("support_date_precision_invalid")
    if any(row.duration_days <= 0 for row in records):
        raise ProtocolViolation("support_nonpositive_follow_up")


def training_censoring_survival(training: tuple[SurvivalRecord, ...], horizon_days: int) -> float:
    """Reverse Kaplan-Meier G(u); remove tied deaths before censoring jumps."""

    _records(training)
    if not training:
        raise ProtocolViolation("training_support_empty")
    if type(horizon_days) is not int or horizon_days <= 0:
        raise ProtocolViolation("support_horizon_invalid")
    if horizon_days >= max(row.duration_days for row in training):
        raise ProtocolViolation("training_follow_up_support_exhausted")
    factors = []
    for time in sorted(
        {row.duration_days for row in training if row.duration_days <= horizon_days}
    ):
        risk = sum(row.duration_days >= time for row in training)
        deaths = sum(row.duration_days == time and row.event_observed for row in training)
        censorings = sum(row.duration_days == time and not row.event_observed for row in training)
        if censorings:
            denominator = risk - deaths
            factors.append(1.0 - censorings / denominator)
    return float(prod(factors))


@dataclass(frozen=True, slots=True)
class HorizonSupport:
    """Inspectable numerical feasibility; passing does not establish power or validity."""

    horizon_days: int
    training_patients: int
    evaluation_patients: int
    training_censoring_survival: float | None
    evaluation_deaths_by_horizon: int
    evaluation_deaths_before_horizon: int
    evaluation_at_risk_before_horizon: int
    evaluation_observed_beyond_horizon: int
    evaluation_early_censorings: int
    reasons: tuple[str, ...]

    @property
    def status(self) -> str:
        return "SUPPORTED_FOR_METRIC_IMPLEMENTATION" if not self.reasons else "NOT_ESTIMABLE"


def assess_horizon_support(
    training: tuple[SurvivalRecord, ...],
    evaluation: tuple[SurvivalRecord, ...],
    *,
    horizon_days: int,
    protocol: ResearchProtocol,
) -> HorizonSupport:
    """Never estimate censoring from held-out outcomes or turn early censoring into survival."""

    if not isinstance(protocol, ResearchProtocol):
        raise ProtocolViolation("support_protocol_invalid")
    if type(horizon_days) is not int or horizon_days not in (
        protocol.primary_horizon_days,
        *protocol.secondary_horizons_days,
    ):
        raise ProtocolViolation("unprespecified_evaluation_horizon")
    _records(training)
    _records(evaluation)
    if {row.patient_id for row in training} & {row.patient_id for row in evaluation}:
        raise ProtocolViolation("training_evaluation_patient_overlap")
    reasons = []
    if not any(row.event_observed for row in training):
        reasons.append("training_has_no_deaths")
    censoring_survival = None
    try:
        censoring_survival = training_censoring_survival(training, horizon_days)
    except ProtocolViolation as error:
        reasons.append(error.code)
    if censoring_survival is not None and (
        censoring_survival < protocol.minimum_training_censoring_survival
    ):
        reasons.append("training_censoring_survival_below_threshold")
    deaths = sum(row.event_observed and row.duration_days <= horizon_days for row in evaluation)
    deaths_before = sum(
        row.event_observed and row.duration_days < horizon_days for row in evaluation
    )
    at_risk = sum(row.duration_days >= horizon_days for row in evaluation)
    survivors = sum(row.duration_days > horizon_days for row in evaluation)
    early_censorings = sum(
        not row.event_observed and row.duration_days <= horizon_days for row in evaluation
    )
    if deaths == 0:
        reasons.append("evaluation_has_no_deaths_by_horizon")
    if deaths_before == 0:
        reasons.append("evaluation_has_no_deaths_before_horizon")
    if survivors == 0:
        reasons.append("evaluation_has_no_observed_survivors_beyond_horizon")
    if at_risk < protocol.minimum_evaluation_at_risk:
        reasons.append("evaluation_at_risk_below_threshold")
    return HorizonSupport(
        horizon_days=horizon_days,
        training_patients=len(training),
        evaluation_patients=len(evaluation),
        training_censoring_survival=censoring_survival,
        evaluation_deaths_by_horizon=deaths,
        evaluation_deaths_before_horizon=deaths_before,
        evaluation_at_risk_before_horizon=at_risk,
        evaluation_observed_beyond_horizon=survivors,
        evaluation_early_censorings=early_censorings,
        reasons=tuple(reasons),
    )


def _counts(records: tuple[SurvivalRecord, ...], day: int) -> dict[str, int]:
    """Separate event cases, observed survivors and unknown early-censored outcomes."""

    return {
        "patients": len(records),
        "deaths_by_horizon": sum(r.event_observed and r.duration_days <= day for r in records),
        "deaths_before_horizon": sum(r.event_observed and r.duration_days < day for r in records),
        "at_risk_before_horizon": sum(r.duration_days >= day for r in records),
        "observed_beyond_horizon": sum(r.duration_days > day for r in records),
        "early_censorings": sum(not r.event_observed and r.duration_days <= day for r in records),
    }


def _risk_scores(values: tuple[float, ...], count: int) -> None:
    """Require finite, aligned risk scores for a calibration grouping audit only."""

    if (
        type(values) is not tuple
        or len(values) != count
        or any(type(value) not in (int, float) or not isfinite(value) for value in values)
    ):
        raise ProtocolViolation("calibration_risk_scores_invalid")


def assess_metric_support(
    training: tuple[SurvivalRecord, ...],
    evaluation: tuple[SurvivalRecord, ...],
    *,
    metric: SurvivalMetric,
    horizon_days: int | None,
    protocol: ResearchProtocol,
    training_risk: tuple[float, ...] | None = None,
    evaluation_risk: tuple[float, ...] | None = None,
) -> dict[str, object]:
    """Audit metric-specific outcome prerequisites; do not calculate or certify scores.

    Records must already be anchored to one clinical setting. The cohort audit
    establishes that scope. Passing these checks cannot establish precision,
    censoring assumptions, model convergence or validity of prediction outputs.
    """

    if not isinstance(protocol, ResearchProtocol) or not isinstance(metric, SurvivalMetric):
        raise ProtocolViolation("metric_support_configuration_invalid")
    _records(training)
    _records(evaluation)
    if {r.patient_id for r in training} & {r.patient_id for r in evaluation}:
        raise ProtocolViolation("training_evaluation_patient_overlap")
    if metric is SurvivalMetric.HARRELL:
        if horizon_days is not None:
            raise ProtocolViolation("harrell_is_prespecified_without_truncation")
    elif type(horizon_days) is not int or horizon_days not in (
        protocol.primary_horizon_days,
        *protocol.secondary_horizons_days,
    ):
        raise ProtocolViolation("unprespecified_evaluation_horizon")
    if metric in (SurvivalMetric.UNO, SurvivalMetric.IBS) and (
        horizon_days != protocol.primary_horizon_days
    ):
        raise ProtocolViolation("unprespecified_metric_horizon")
    if metric is not SurvivalMetric.CALIBRATION and (
        training_risk is not None or evaluation_risk is not None
    ):
        raise ProtocolViolation("calibration_scores_supplied_to_wrong_metric")
    day = horizon_days if type(horizon_days) is int else 0
    fatal = []
    stability = []
    missing = []
    details = {}
    if len(training) < 2:
        fatal.append("training_has_fewer_than_two_patients")
    if not any(r.event_observed for r in training):
        fatal.append("training_has_no_deaths_for_prespecified_cox")
    if not evaluation:
        fatal.append("evaluation_empty")
    if metric in (SurvivalMetric.UNO, SurvivalMetric.HARRELL):
        anchors = [
            r
            for r in evaluation
            if r.event_observed and (metric is SurvivalMetric.HARRELL or r.duration_days < day)
        ]
        pairs = [
            sum(other.duration_days > r.duration_days for other in evaluation) for r in anchors
        ]
        details["comparable_pairs"] = sum(pairs)
        details["distinct_comparable_event_patients"] = sum(count > 0 for count in pairs)
        if not sum(pairs):
            fatal.append("no_comparable_event_pairs")
    if metric is not SurvivalMetric.HARRELL:
        if horizon_days is None:
            raise ProtocolViolation("support_horizon_invalid")
        counts = _counts(evaluation, horizon_days)
        details.update(counts)
        try:
            g = training_censoring_survival(training, horizon_days)
            details["training_G_at_horizon"] = g
            if g <= 0:
                fatal.append("training_censoring_survival_zero")
            elif g < protocol.minimum_training_censoring_survival:
                stability.append("training_censoring_survival_below_threshold")
        except ProtocolViolation as error:
            fatal.append(error.code)
        if counts["at_risk_before_horizon"] < protocol.minimum_evaluation_at_risk:
            stability.append("evaluation_at_risk_below_threshold")
        if metric is SurvivalMetric.AUC:
            if not counts["deaths_by_horizon"]:
                fatal.append("auc_has_no_cases_by_horizon")
            if not counts["observed_beyond_horizon"]:
                fatal.append("auc_has_no_controls_beyond_horizon")
        if metric in (SurvivalMetric.BRIER, SurvivalMetric.IBS):
            if not counts["deaths_by_horizon"] + counts["observed_beyond_horizon"]:
                fatal.append("brier_has_no_known_outcome_contributions")
        if metric is SurvivalMetric.IBS:
            # G and risk-set size are nonincreasing: their tail values bound every grid day.
            grid = range(1, protocol.primary_horizon_days + 1)
            unsupported = next(
                (
                    day
                    for day in grid
                    if not (
                        _counts(evaluation, day)["deaths_by_horizon"]
                        + _counts(evaluation, day)["observed_beyond_horizon"]
                    )
                ),
                None,
            )
            details["grid"] = {
                "start_day": 1,
                "stop_day": horizon_days,
                "step_days": 1,
                "points": len(grid),
                "first_day_without_known_contributions": unsupported,
                "G_and_risk_checked_by_monotone_tail_bounds": True,
            }
            if unsupported is not None:
                fatal.append("integrated_brier_grid_has_unsupported_day")
        if metric is SurvivalMetric.CALIBRATION:
            if training_risk is None or evaluation_risk is None:
                missing.append("frozen_training_and_evaluation_risk_scores_for_group_support")
            else:
                _risk_scores(training_risk, len(training))
                _risk_scores(evaluation_risk, len(evaluation))
                if len(training_risk) < 2:
                    fatal.append("calibration_training_risk_insufficient")
                else:
                    cuts = quantiles(training_risk, n=3, method="inclusive")
                    details["training_only_cutpoints"] = cuts
                    if cuts[0] == cuts[1]:
                        fatal.append("calibration_training_tertiles_tied")
                    groups = []
                    for group in range(3):
                        members = tuple(
                            r
                            for r, score in zip(evaluation, evaluation_risk, strict=True)
                            if bisect_right(cuts, score) == group
                        )
                        group_counts = _counts(members, horizon_days)
                        groups.append(group_counts)
                        if not members:
                            fatal.append(f"calibration_group_{group + 1}_empty")
                        if (
                            group_counts["at_risk_before_horizon"]
                            < protocol.minimum_evaluation_at_risk
                        ):
                            stability.append(
                                f"calibration_group_{group + 1}_at_risk_below_threshold"
                            )
                    details["fixed_training_risk_groups"] = groups
    status = (
        "NOT_ESTIMABLE"
        if fatal
        else "NOT_ASSESSED"
        if missing
        else "WITHHELD_BY_PROTOCOL"
        if stability
        else "OUTCOME_PREREQUISITES_MET"
    )
    probability_metric = metric in (
        SurvivalMetric.BRIER,
        SurvivalMetric.IBS,
        SurvivalMetric.CALIBRATION,
    )
    return {
        "metric": metric.value,
        "horizon_days": horizon_days,
        "status": status,
        "fatal_reasons": sorted(set(fatal)),
        "stability_reasons": sorted(set(stability)),
        "missing_inputs": missing,
        "details": details,
        "scope": "outcome_prerequisites_only",
        "precision_status": "NOT_ESTABLISHED",
        "prediction_outputs_validated": False,
        "additional_prediction_requirements": (
            "Aligned finite survival probabilities, training-only baseline hazard, "
            "and the metric's complete prespecified time grid."
            if probability_metric
            else "Aligned finite risk predictions from both frozen fitted comparison arms."
        ),
        "model_evaluation": "NOT_RUN",
        "interpretation": (
            "Numerical support, pair counts and stability floors do not establish "
            "statistical precision or power."
        ),
    }
