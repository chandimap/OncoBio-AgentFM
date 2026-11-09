"""Descriptive selection contrasts; no predictor fitting or outcome-based cohort choice."""

from __future__ import annotations

from collections import Counter
from math import isfinite, sqrt
from statistics import mean, variance

from oncobio_agentfm.cohort_sources import SourceInventory, require
from oncobio_agentfm.prognosis import ProtocolViolation
from oncobio_agentfm.source_audit import SOURCE_MISSING_TOKENS

CATEGORIES = {
    "Gender": ("Female", "Male"),
    "Patient affiliation": ("Stanford", "VA"),
    "Smoking status": ("Current", "Former", "Nonsmoker"),
}


def _age(value: str) -> float | None:
    """Describe source age at histological diagnosis; do not relabel it as index age."""

    if value in SOURCE_MISSING_TOKENS:
        return None
    try:
        result = float(value)
    except ValueError:
        raise ProtocolViolation("selection_age_invalid") from None
    require(isfinite(result) and 0 <= result <= 120, "selection_age_invalid")
    return result


def _continuous(first: list[float], second: list[float]) -> dict[str, object]:
    """Use unweighted group variances; leave undefined contrasts explicitly null."""

    def summary(values: list[float]) -> dict[str, object]:
        return {
            "observed_n": len(values),
            "mean": mean(values) if values else None,
            "sd": sqrt(variance(values)) if len(values) > 1 else None,
        }

    smd = None
    if len(first) > 1 and len(second) > 1:
        denominator = sqrt((variance(first) + variance(second)) / 2)
        if denominator > 0:
            smd = (mean(first) - mean(second)) / denominator
    return {"complete": summary(first), "incomplete": summary(second), "smd": smd}


def _compare_records(records, complete: frozenset[str]) -> dict[str, object]:
    """Compare complete-source versus incomplete patients inside the clinical denominator."""

    require(complete <= records.keys(), "selection_denominator_mismatch")
    incomplete = frozenset(records) - complete
    groups = [[records[patient] for patient in sorted(ids)] for ids in (complete, incomplete)]
    ages = [[_age(row["Age at Histological Diagnosis"]) for row in group] for group in groups]
    age = _continuous(*[[value for value in values if value is not None] for values in ages])
    age["complete"]["missing_n"] = sum(value is None for value in ages[0])
    age["incomplete"]["missing_n"] = sum(value is None for value in ages[1])
    categorical = {}
    for column, levels in CATEGORIES.items():
        counts = [
            Counter(
                row[column] if row[column] in levels else "missing_or_unrecognized" for row in group
            )
            for group in groups
        ]
        denominators = [sum(counter[level] for level in levels) for counter in counts]
        contrasts = []
        for level in levels:
            proportions = [
                counter[level] / total if total else None
                for counter, total in zip(counts, denominators, strict=True)
            ]
            first, second = proportions
            difference = smd = None
            if first is not None and second is not None:
                difference = first - second
                denominator = sqrt((first * (1 - first) + second * (1 - second)) / 2)
                if denominator > 0:
                    smd = difference / denominator
            contrasts.append(
                {
                    "level": level,
                    "complete_n": counts[0][level],
                    "incomplete_n": counts[1][level],
                    "proportions": proportions,
                    "proportion_difference": difference,
                    "smd": smd,
                }
            )
        categorical[column] = {
            "observed_denominators": denominators,
            "missing_or_unrecognized": [c["missing_or_unrecognized"] for c in counts],
            "levels": contrasts,
        }
    return {
        "scope": "source_complete_vs_incomplete_within_clinical_records",
        "complete_n": len(complete),
        "incomplete_n": len(incomplete),
        "age_at_histological_diagnosis": age,
        "categorical": categorical,
        "direction": "complete_minus_incomplete",
        "smd_definition": "mean_or_proportion_difference_over_sqrt_mean_of_group_variances",
        "inference": "DESCRIPTIVE_ONLY_NO_CAUSAL_OR_MISSINGNESS_MECHANISM_CLAIM",
        "limitations": [
            "Source completeness does not establish clinical eligibility.",
            "Observed differences cannot measure unobserved selection or remove bias.",
            "Source age and recorded gender are not adjudicated index age or gender identity.",
            "No p-value, imputation, balancing or outcome-guided selection is performed.",
        ],
    }


def compare_selection(inventory: SourceInventory, complete: frozenset[str]) -> dict[str, object]:
    """Retain the pooled description and distinguish the source recruitment cohorts."""

    result = _compare_records(inventory.clinical, complete)
    strata = {}
    for name, prefix in (("R01", "R01-"), ("AMC", "AMC-"), ("other", None)):
        records = {
            patient: row
            for patient, row in inventory.clinical.items()
            if (
                patient.startswith(prefix)
                if prefix is not None
                else not patient.startswith(("R01-", "AMC-"))
            )
        }
        if records:
            summary = _compare_records(records, complete & records.keys())
            summary["comparison_support"] = (
                "BOTH_GROUPS_OBSERVED"
                if summary["complete_n"] and summary["incomplete_n"]
                else "NO_WITHIN_STRATUM_COMPARISON"
            )
            strata[name] = summary
    result["source_cohort_strata"] = strata
    result["stratification_scope"] = (
        "source_R01_AMC_identifier_prefixes_only; "
        "descriptive_not_adjustment_or_causal_identification"
    )
    return result
