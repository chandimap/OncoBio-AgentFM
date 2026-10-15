"""Versioned, outcome-independent NSCLC prognostic research specification."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path


class PredictionSetting(StrEnum):
    """Separate populations and time origins; never pool their performance."""

    PRETREATMENT = "pretreatment"
    POSTOPERATIVE = "postoperative_day_30"


class FeatureKind(StrEnum):
    """Typed input roles, including prohibited roles for adversarial validation."""

    CT = "ct"
    PET_CT = "pet_ct"
    RESECTION_MOLECULAR = "resection_molecular"
    PATHOLOGICAL_STAGE = "pathological_stage"
    RECURRENCE = "recurrence"
    ADJUVANT_TREATMENT = "adjuvant_treatment"
    SURVIVAL_OUTCOME = "survival_outcome"


class SurvivalMetric(StrEnum):
    """Different metric denominators require different feasibility checks."""

    UNO = "uno_concordance"
    HARRELL = "harrell_concordance"
    AUC = "cumulative_dynamic_auc"
    BRIER = "ipcw_brier"
    IBS = "integrated_brier"
    CALIBRATION = "calibration_tertiles"


@dataclass(frozen=True, slots=True)
class ResearchProtocol:
    """The executable choices extracted from the reviewed JSON specification."""

    protocol_id: str
    sha256: str
    landmark_days: int
    primary_horizon_days: int
    secondary_horizons_days: tuple[int, ...]
    minimum_training_censoring_survival: float = 0.10
    minimum_evaluation_at_risk: int = 10
    revision: str = "2026-10-15"

    def __post_init__(self) -> None:
        if self.protocol_id != "nsclc-os-v1":
            raise ValueError("unsupported protocol identity")
        if self.revision != "2026-10-15":
            raise ValueError("unsupported protocol revision")
        if (
            not isinstance(self.sha256, str)
            or len(self.sha256) != 64
            or any(character not in "0123456789abcdef" for character in self.sha256)
        ):
            raise ValueError("protocol requires an exact SHA-256 fingerprint")
        if type(self.landmark_days) is not int or self.landmark_days != 30:
            raise ValueError("unsupported landmark")
        if type(self.primary_horizon_days) is not int or self.primary_horizon_days != 730:
            raise ValueError("unsupported primary horizon")
        if type(self.secondary_horizons_days) is not tuple or self.secondary_horizons_days != (
            365,
            1095,
        ):
            raise ValueError("unsupported secondary horizons")
        if any(type(value) is not int for value in self.secondary_horizons_days):
            raise ValueError("secondary horizons require integer days")
        if type(self.minimum_training_censoring_survival) is not float or (
            self.minimum_training_censoring_survival != 0.10
        ):
            raise ValueError("unsupported censoring support threshold")
        if (
            type(self.minimum_evaluation_at_risk) is not int
            or self.minimum_evaluation_at_risk != 10
        ):
            raise ValueError("unsupported at-risk support threshold")

    def required_features(self, setting: PredictionSetting) -> frozenset[FeatureKind]:
        """Require explicit missingness rather than silently omit a modality."""

        if not isinstance(setting, PredictionSetting):
            raise TypeError("setting must be a PredictionSetting")
        roles = {FeatureKind.CT, FeatureKind.PET_CT}
        if setting is PredictionSetting.POSTOPERATIVE:
            roles.add(FeatureKind.RESECTION_MOLECULAR)
        return frozenset(roles)


def _unique_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    """Reject duplicate JSON keys that could disguise a protocol change."""

    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("protocol JSON contains duplicate keys")
        result[key] = value
    return result


def _same_typed_value(actual: object, expected: object) -> bool:
    """Reject bool/int and nested policy substitutions, not only unequal values."""

    if type(actual) is not type(expected):
        return False
    if isinstance(actual, dict) and isinstance(expected, dict):
        return actual.keys() == expected.keys() and all(
            _same_typed_value(actual[key], value) for key, value in expected.items()
        )
    if isinstance(actual, list) and isinstance(expected, list):
        return len(actual) == len(expected) and all(
            _same_typed_value(a, b) for a, b in zip(actual, expected, strict=True)
        )
    return actual == expected


def load_research_protocol(path: Path) -> ResearchProtocol:
    """Fail closed on unsupported scientific choices; hash exact reviewed bytes."""

    raw = path.read_bytes()
    spec = json.loads(raw, object_pairs_hook=_unique_object)
    expected = {
        "schema_version": "1.2",
        "protocol_id": "nsclc-os-v1",
        "revision": "2026-10-15",
        "endpoint": "all_cause_death",
        "primary_index": "earliest_qualifying_pretreatment_assessment_end_of_day",
        "secondary_index": "surgery_plus_fixed_days",
        "landmark_days": 30,
        "primary_horizon_days": 730,
        "secondary_horizons_days": [365, 1095],
        "primary_metric": "uno_ipcw_concordance_at_730_days",
        "split_seed": "oncobio-agentfm-v1",
        "split_fractions": [0.70, 0.15, 0.15],
        "minimum_training_censoring_survival": 0.10,
        "minimum_evaluation_at_risk": 10,
        "bootstrap_replicates": 2000,
        "bootstrap_seed": 20261014,
        "chronology_policy": {
            "frozen_imaging": "unique_original_roles_and_same_selected_evidence_at_landmark",
            "administrative_cutoff": "one_reviewed_policy_reference_on_verified_patient_axes",
            "source_consistency": "no_documented_cancer_treatment_or_surgery_after_death",
        },
        "feasibility_policy": {
            "scope": "outcome_prerequisites_only",
            "precision_status": "NOT_ESTABLISHED",
            "ibs_grid_start_day": 1,
            "ibs_grid_stop_day": 730,
            "ibs_grid_step_days": 1,
            "calibration_groups": 3,
            "calibration_cutpoints": "training_only_inclusive_linear_tertiles",
            "calibration_stability": "existing_at_risk_floor_applied_per_group",
            "comparison_cohorts": "same_modality_eligible_patients_in_both_arms_and_all_splits",
        },
    }
    narrative_keys = {"population", "estimands", "hypotheses", "evaluation_rules", "claim_limits"}
    if not isinstance(spec, dict) or set(spec) != set(expected) | narrative_keys:
        raise ValueError("protocol fields differ from the supported schema")
    for key, value in expected.items():
        if not _same_typed_value(spec[key], value):
            raise ValueError(f"unsupported protocol choice: {key}")
    if not isinstance(spec["population"], str) or not spec["population"].strip():
        raise ValueError("protocol requires a population definition")
    estimands = spec["estimands"]
    if not isinstance(estimands, dict) or set(estimands) != {"primary", "secondary"}:
        raise ValueError("protocol requires distinct primary and secondary estimands")
    if any(not isinstance(value, str) or not value.strip() for value in estimands.values()):
        raise ValueError("estimand definitions must be non-empty strings")
    hypotheses = spec["hypotheses"]
    if not isinstance(hypotheses, list) or len(hypotheses) != 4:
        raise ValueError("protocol requires the four prespecified hypotheses")
    for hypothesis, identifier in zip(hypotheses, ("H1", "H2", "H3", "H4"), strict=True):
        if not isinstance(hypothesis, dict) or set(hypothesis) != {"id", "setting", "contrast"}:
            raise ValueError("hypothesis fields differ from the supported schema")
        setting = "postoperative_day_30" if identifier == "H4" else "pretreatment"
        if hypothesis["id"] != identifier or hypothesis["setting"] != setting:
            raise ValueError("hypothesis identity or prediction setting changed")
        if not isinstance(hypothesis["contrast"], str) or not hypothesis["contrast"].strip():
            raise ValueError("hypothesis requires an explicit contrast")
    for key in ("evaluation_rules", "claim_limits"):
        values = spec[key]
        if not isinstance(values, list) or not values:
            raise ValueError(f"protocol requires {key}")
        if any(not isinstance(value, str) or not value.strip() for value in values):
            raise ValueError(f"{key} must contain non-empty statements")
    return ResearchProtocol(
        protocol_id=spec["protocol_id"],
        sha256=hashlib.sha256(raw).hexdigest(),
        landmark_days=spec["landmark_days"],
        primary_horizon_days=spec["primary_horizon_days"],
        secondary_horizons_days=tuple(spec["secondary_horizons_days"]),
        revision=spec["revision"],
    )
