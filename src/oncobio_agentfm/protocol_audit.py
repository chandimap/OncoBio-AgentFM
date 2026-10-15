"""Aggregate eligibility evidence and a gated bridge to the existing clinical core."""

from __future__ import annotations

import json
from collections import Counter
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path

from oncobio_agentfm.clinical import NSCLCCohort, NSCLCOverallSurvivalRecord
from oncobio_agentfm.contracts import SurvivalRecord
from oncobio_agentfm.evaluation_support import assess_horizon_support, assess_metric_support
from oncobio_agentfm.prognosis import (
    OverallSurvivalFollowUp,
    PredictionContext,
    PredictorFeature,
    ProtocolViolation,
    validate_predictors,
)
from oncobio_agentfm.research_protocol import (
    FeatureKind,
    PredictionSetting,
    ResearchProtocol,
    SurvivalMetric,
)
from oncobio_agentfm.splitting import CohortSplit, split_patient_ids


@dataclass(frozen=True, slots=True)
class AdjudicatedCase:
    """Local source-adjudicated data; never serialize this object for publication."""

    context: PredictionContext
    primary_features: tuple[PredictorFeature, ...]
    landmark_features: tuple[PredictorFeature, ...]
    follow_up: OverallSurvivalFollowUp

    def __post_init__(self) -> None:
        if not isinstance(self.context, PredictionContext) or not isinstance(
            self.follow_up, OverallSurvivalFollowUp
        ):
            raise ProtocolViolation("adjudicated_case_type_invalid")
        for features in (self.primary_features, self.landmark_features):
            if type(features) is not tuple or any(
                not isinstance(feature, PredictorFeature) for feature in features
            ):
                raise ProtocolViolation("features_must_be_validated_tuple")
        death = self.follow_up.death_on
        if death is not None:
            treatment = self.context.first_treatment_on
            surgery = self.context.surgery
            if treatment is not None and treatment > death:
                raise ProtocolViolation("documented_treatment_after_death")
            if surgery is not None and surgery.is_present:
                # Validate precision before comparing raw surgical and outcome clocks.
                surgery_date = surgery.acquired_on
                if type(surgery_date) is not date:
                    raise ProtocolViolation("surgery_date_invalid")
                if surgery_date > death:
                    raise ProtocolViolation("documented_surgery_after_death")

    def features_for(self, setting: PredictionSetting) -> tuple[PredictorFeature, ...]:
        """Keep baseline and landmark feature sets distinct."""

        if not isinstance(setting, PredictionSetting):
            raise TypeError("setting must be a PredictionSetting")
        if setting is PredictionSetting.POSTOPERATIVE:
            for kind in (FeatureKind.CT, FeatureKind.PET_CT):
                primary = tuple(item for item in self.primary_features if item.kind is kind)
                landmark = tuple(item for item in self.landmark_features if item.kind is kind)
                if len(primary) != 1:
                    raise ProtocolViolation("frozen_baseline_imaging_not_uniquely_declared")
                # Missing/duplicate roles retain the main gate's existing diagnostics.
                if len(primary) == len(landmark) == 1 and (
                    primary[0].evidence != landmark[0].evidence
                    or primary[0].availability_basis != landmark[0].availability_basis
                ):
                    raise ProtocolViolation("baseline_imaging_not_frozen_across_settings")
        return (
            self.primary_features
            if setting is PredictionSetting.PRETREATMENT
            else self.landmark_features
        )


def _unique_cases(cases: tuple[AdjudicatedCase, ...]) -> None:
    """Reject duplicate patients instead of treating rows as independent people."""

    if (
        type(cases) is not tuple
        or not cases
        or any(not isinstance(case, AdjudicatedCase) for case in cases)
    ):
        raise ValueError("cases must be a non-empty tuple of AdjudicatedCase records")
    ids = [case.context.baseline.patient_id for case in cases]
    if len(set(ids)) != len(ids):
        raise ProtocolViolation("duplicate_patient")
    cutoff_policies = {case.follow_up.administrative_cutoff_source for case in cases}
    if len(cutoff_policies) != 1:
        raise ProtocolViolation("cohort_administrative_cutoff_policy_mismatch")


def _comparison_support(
    outcomes: dict[str, SurvivalRecord], hypotheses: tuple[str, ...], protocol: ResearchProtocol
) -> dict[str, object]:
    """Use exactly one prespecified source-available cohort in both arms and every split."""

    split = split_patient_ids(outcomes) if outcomes else CohortSplit((), (), ())
    partitions = {"train": split.train, "validation": split.validation, "test": split.test}
    training = tuple(outcomes[identifier] for identifier in split.train)
    counts = {}
    support = {}
    for name, ids in partitions.items():
        records = tuple(outcomes[identifier] for identifier in ids)
        deaths = sum(record.event_observed for record in records)
        counts[name] = {
            "patients": len(records),
            "deaths": deaths,
            "censored": len(records) - deaths,
        }
        if name == "train":
            continue
        reports = []
        for metric in SurvivalMetric:
            horizons = (
                (None,)
                if metric is SurvivalMetric.HARRELL
                else (
                    (protocol.primary_horizon_days,)
                    if metric in (SurvivalMetric.UNO, SurvivalMetric.IBS)
                    else tuple(
                        sorted((protocol.primary_horizon_days, *protocol.secondary_horizons_days))
                    )
                )
            )
            for horizon in horizons:
                reports.append(
                    assess_metric_support(
                        training, records, metric=metric, horizon_days=horizon, protocol=protocol
                    )
                )
        support[name] = reports
    return {
        "hypotheses": list(hypotheses),
        "partition_counts": counts,
        "metric_support": support,
        "comparison_policy": (
            "same eligible patients in both arms across training, validation and test"
        ),
        "precision_status": "NOT_ESTABLISHED",
        "model_evaluation": "NOT_RUN",
    }


def audit_eligibility(
    cases: tuple[AdjudicatedCase, ...],
    *,
    setting: PredictionSetting,
    protocol: ResearchProtocol,
    synthetic: bool = False,
) -> dict[str, object]:
    """Count first-failure reasons and outcome joins without exporting patient data."""

    _unique_cases(cases)
    if type(synthetic) is not bool:
        raise ProtocolViolation("evidence_scope_flag_invalid")
    protocol.required_features(setting)
    predictor_rejections: Counter[str] = Counter()
    outcome_rejections: Counter[str] = Counter()
    predictor_eligible = 0
    accepted = {}
    present_features = {}
    for case in cases:
        try:
            predictors = validate_predictors(
                case.context, case.features_for(setting), setting=setting, protocol=protocol
            )
        except ProtocolViolation as exc:
            predictor_rejections[exc.code] += 1
            continue
        predictor_eligible += 1
        try:
            outcome = case.follow_up.anchor(predictors)
        except ProtocolViolation as exc:
            outcome_rejections[exc.code] += 1
            continue
        accepted[predictors.patient_id] = outcome
        present_features[predictors.patient_id] = {
            item.kind for item in predictors.features if item.evidence.is_present
        }
    partitions = {}
    support_reports = {}
    if accepted:
        split = split_patient_ids(accepted)
        for name, ids in (
            ("train", split.train),
            ("validation", split.validation),
            ("test", split.test),
        ):
            deaths = sum(accepted[patient_id].event_observed for patient_id in ids)
            partitions[name] = {
                "patients": len(ids),
                "deaths": deaths,
                "censored": len(ids) - deaths,
                "at_risk_at_730_days": sum(
                    accepted[patient_id].duration_days >= protocol.primary_horizon_days
                    for patient_id in ids
                ),
            }
            if name != "train":
                reports = []
                for horizon in sorted(
                    (protocol.primary_horizon_days, *protocol.secondary_horizons_days)
                ):
                    support = assess_horizon_support(
                        tuple(accepted[patient_id] for patient_id in split.train),
                        tuple(accepted[patient_id] for patient_id in ids),
                        horizon_days=horizon,
                        protocol=protocol,
                    )
                    reports.append(
                        {
                            **asdict(support),
                            "reasons": list(support.reasons),
                            "status": support.status,
                        }
                    )
                support_reports[name] = reports
    comparisons = {}
    if setting is PredictionSetting.PRETREATMENT:
        comparisons["clinical_ct"] = _comparison_support(
            accepted, ("H1", "H3_hide_all_imaging"), protocol
        )
        pet = {
            identifier: outcome
            for identifier, outcome in accepted.items()
            if FeatureKind.PET_CT in present_features[identifier]
        }
        comparisons["pet_available"] = _comparison_support(pet, ("H2", "H3_hide_pet"), protocol)
    else:
        molecular = {
            identifier: outcome
            for identifier, outcome in accepted.items()
            if FeatureKind.RESECTION_MOLECULAR in present_features[identifier]
        }
        comparisons["molecular_available"] = _comparison_support(molecular, ("H4",), protocol)
    return {
        "schema_version": "1.2",
        "protocol_id": protocol.protocol_id,
        "protocol_revision": protocol.revision,
        "protocol_sha256": protocol.sha256,
        "evidence_type": "synthetic_software_validation"
        if synthetic
        else "local_adjudication_audit",
        "setting": setting.value,
        "input_patients": len(cases),
        "predictor_eligible": predictor_eligible,
        "endpoint_eligible": len(accepted),
        "predictor_rejections": dict(sorted(predictor_rejections.items())),
        "outcome_rejections": dict(sorted(outcome_rejections.items())),
        "partition_counts": partitions,
        "horizon_support": support_reports,
        "horizon_support_scope": (
            "retained conservative screen; individual metrics use comparison_feasibility"
        ),
        "comparison_feasibility": comparisons,
        "precision_status": "NOT_ESTABLISHED",
        "reason_policy": (
            "first_failure_per_gate; malformed source records must be audited separately"
        ),
        "performance_evaluation": "NOT_RUN",
        "clinical_validity": "NOT_ESTABLISHED_BY_SOFTWARE_CHECKS",
    }


def build_primary_clinical_cohort(
    cases: tuple[AdjudicatedCase, ...], *, protocol: ResearchProtocol, ajcc_version: int
) -> NSCLCCohort:
    """Gate every primary record before reusing the existing Cox cohort contract."""

    _unique_cases(cases)
    records = []
    for case in cases:
        predictors = validate_predictors(
            case.context,
            case.primary_features,
            setting=PredictionSetting.PRETREATMENT,
            protocol=protocol,
        )
        records.append(
            NSCLCOverallSurvivalRecord(predictors.baseline, case.follow_up.anchor(predictors))
        )
    return NSCLCCohort(tuple(records), ajcc_version)


def write_report(path: Path, report: dict[str, object]) -> None:
    """Create a deterministic report exclusively; never overwrite recorded evidence."""

    serialized = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(serialized)
