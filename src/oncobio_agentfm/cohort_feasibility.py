"""Aggregate endpoint diagnostics; never construct outcomes from unadjudicated rows."""

from __future__ import annotations

import re
from collections import Counter
from collections.abc import Mapping
from datetime import datetime

from oncobio_agentfm.cohort_sources import SourceInventory
from oncobio_agentfm.source_audit import SOURCE_MISSING_TOKENS

DATE_FIELDS = ("CT Date", "PET Date", "Date of Death", "Date of Last Known Alive")
DURATION_FIELDS = ("Days between CT and surgery", "Time to Death (days)")


def _date_ordinal(value: str) -> tuple[int | None, str]:
    """Use the reviewed CSV format only; keep missing and invalid values distinct."""

    if value in SOURCE_MISSING_TOKENS:
        return None, "source_missing_token"
    if re.fullmatch(r"[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}", value) is not None:
        try:
            return datetime.strptime(value, "%m/%d/%Y").date().toordinal(), "parsed"
        except ValueError:
            pass
    return None, "invalid"


def _integer_days(value: str) -> tuple[int | None, str]:
    """Do not round fractional days, strip tokens, or coerce unknown durations to zero."""

    if value in SOURCE_MISSING_TOKENS:
        return None, "source_missing_token"
    if re.fullmatch(r"[+-]?[0-9]{1,10}", value) is None:
        return None, "invalid"
    return int(value), "parsed"


def _profile(records: Mapping[str, Mapping[str, str]]) -> dict[str, object]:
    """Count syntax and diagnostic contradictions without selecting analysis participants."""

    states = {field: Counter() for field in (*DATE_FIELDS, *DURATION_FIELDS)}
    labels = Counter()
    counts = Counter(
        dict.fromkeys(
            (
                "alive_with_death_date",
                "dead_without_parseable_death_date",
                "last_alive_after_death",
                "alive_with_numeric_time_to_death",
                "death_CT_pairs",
                "death_on_or_before_CT",
                "alive_last_contact_CT_pairs",
                "alive_last_contact_on_or_before_CT",
                "time_to_death_CT_comparison_pairs",
                "time_to_death_matches_death_minus_CT",
                "candidate_surgery_records",
                "candidate_surgery_recorded_deaths",
                "candidate_death_on_or_before_surgery",
                "candidate_surgery_recorded_alive",
                "candidate_alive_contact_on_or_before_surgery",
                "time_to_death_surgery_comparison_pairs",
                "time_to_death_matches_death_minus_candidate_surgery",
                "time_to_death_disagrees_with_death_minus_candidate_surgery",
            ),
            0,
        )
    )
    for row in records.values():
        dates, durations = {}, {}
        for field in DATE_FIELDS:
            dates[field], state = _date_ordinal(row[field])
            states[field][state] += 1
        for field in DURATION_FIELDS:
            durations[field], state = _integer_days(row[field])
            states[field][state] += 1
            if durations[field] is not None:
                states[field]["negative"] += durations[field] < 0
                states[field]["zero"] += durations[field] == 0
        status = row["Survival Status"]
        labels[status if status in {"Dead", "Alive"} else "unrecognized"] += 1
        ct, death, last_alive = (
            dates[field] for field in ("CT Date", "Date of Death", "Date of Last Known Alive")
        )
        gap, ttd = (durations[field] for field in DURATION_FIELDS)
        counts["alive_with_death_date"] += status == "Alive" and death is not None
        counts["dead_without_parseable_death_date"] += status == "Dead" and death is None
        counts["last_alive_after_death"] += (
            death is not None and last_alive is not None and last_alive > death
        )
        counts["alive_with_numeric_time_to_death"] += status == "Alive" and ttd is not None
        if status == "Dead" and death is not None and ct is not None:
            counts["death_CT_pairs"] += 1
            counts["death_on_or_before_CT"] += death <= ct
            if ttd is not None:
                counts["time_to_death_CT_comparison_pairs"] += 1
                counts["time_to_death_matches_death_minus_CT"] += ttd == death - ct
        if status == "Alive" and last_alive is not None and ct is not None:
            counts["alive_last_contact_CT_pairs"] += 1
            counts["alive_last_contact_on_or_before_CT"] += last_alive <= ct
        if ct is None or gap is None or gap < 0:
            continue
        # An ordinal sum is a diagnostic candidate, never an adjudicated surgery date.
        # Integer arithmetic also avoids overflowing date/timedelta on malformed intervals.
        candidate_surgery = ct + gap
        counts["candidate_surgery_records"] += 1
        if status == "Dead" and death is not None:
            counts["candidate_surgery_recorded_deaths"] += 1
            counts["candidate_death_on_or_before_surgery"] += death <= candidate_surgery
            if ttd is not None:
                agrees = ttd == death - candidate_surgery
                counts["time_to_death_surgery_comparison_pairs"] += 1
                counts["time_to_death_matches_death_minus_candidate_surgery"] += agrees
                counts["time_to_death_disagrees_with_death_minus_candidate_surgery"] += not agrees
        if status == "Alive" and last_alive is not None:
            counts["candidate_surgery_recorded_alive"] += 1
            counts["candidate_alive_contact_on_or_before_surgery"] += (
                last_alive <= candidate_surgery
            )
    return {
        "clinical_rows": len(records),
        "recorded_dead": labels["Dead"],
        "recorded_alive": labels["Alive"],
        "unrecognized_status": labels["unrecognized"],
        "field_states": {
            field: {
                key: states[field][key]
                for key in (
                    ("parsed", "source_missing_token", "invalid", "negative", "zero")
                    if field in DURATION_FIELDS
                    else ("parsed", "source_missing_token", "invalid")
                )
            }
            for field in states
        },
        "arithmetic_diagnostics": dict(counts),
    }


def endpoint_feasibility(inventory: SourceInventory) -> dict[str, object]:
    """Expose non-identifiability of the existing estimands from the locked source bundle."""

    ct = {s.patient_id for s in inventory.series if s.modality == "CT"}
    pet = {s.patient_id for s in inventory.series if s.modality == "PT"}
    rna = set(inventory.rna_samples.values()) & {
        patient for patient, count in inventory.expression_observed_per_patient.items() if count > 0
    }
    complete = ct & pet & rna & inventory.clinical.keys()
    return {
        "evidence_type": "source_endpoint_feasibility_not_clinical_outcome_adjudication",
        "scope_decision": "SOURCE_COHORT_AND_ENDPOINT_FEASIBILITY_ONLY",
        "clinical_date_semantics": {
            "reference": "https://doi.org/10.1038/sdata.2018.202",
            "locator": "Methods/Subject Demographics and Clinical Data",
            "established": "clinical_dates_shifted_with_relative_chronological_order_preserved",
            "not_established": "clinical_to_NBIA_patient_axis_alignment_or_clinical_release_times",
        },
        "surgery_diagnostic_definition": (
            "candidate_ordinal = parsed_CT_Date + nonnegative_integer_gap"
        ),
        "surgery_diagnostic_limit": (
            "Arithmetic hypothesis check only: neither origin of Time to Death nor patient-level "
            "first-treatment surgery is adjudicated. Agreement is not clinical validation; "
            "disagreement identifies a query, not which source field is wrong."
        ),
        "source_populations": {
            "clinical_release": _profile(inventory.clinical),
            "all_four_source_types": _profile(
                {patient: inventory.clinical[patient] for patient in sorted(complete)}
            ),
        },
        "unresolved_required_facts": {
            "pretreatment_assessment": "no_documented_earliest_qualifying_assessment_or_t0",
            "diagnosis_and_baseline_access": "no_patient_level_result_release_chronology",
            "untreated_status": "no_affirmative_untreated_at_t0_or_first_treatment_chronology",
            "selected_images": "no_adjudicated_series_acquisition_access_and_pair_components",
            "calendar_alignment": "clinical_relative_order_does_not_verify_cross_source_axes",
            "first_treatment_resection": "CT_gap_does_not_adjudicate_first_treatment_resection",
            "molecular_access": "no_patient_level_clinical_assay_release_by_surgery_plus_30_days",
            "OS_cutoff": "no_reviewed_administrative_cutoff_policy_required_by_existing_protocol",
            "survival_origin": "no_source_defined_origin_adjudication_in_locked_bundle",
        },
        "analysis_populations": {
            setting: {
                "status": "NOT_IDENTIFIABLE_FROM_LOCKED_SOURCES",
                "eligible_patients": None,
                "survival_events": None,
            }
            for setting in ("pretreatment", "postoperative_day_30")
        },
        "interpretation": (
            "Unknown eligible/event counts are not zero. Recorded deaths are not validated "
            "events. No date imputation, origin selection, cutoff selection or outcome creation "
            "is performed. A CT- or surgery-indexed research association would change the "
            "estimand and still requires justified endpoint semantics."
        ),
        "privacy": "aggregate_only_no_source_identifiers_patient_dates_or_raw_tokens",
    }


def private_endpoint_queries(inventory: SourceInventory) -> dict[str, object]:
    """Trace diagnostic queries to immutable CSV coordinates; keep this artifact local."""

    query_counters = (
        "alive_with_death_date",
        "dead_without_parseable_death_date",
        "last_alive_after_death",
        "alive_with_numeric_time_to_death",
        "death_on_or_before_CT",
        "alive_last_contact_on_or_before_CT",
        "candidate_death_on_or_before_surgery",
        "candidate_alive_contact_on_or_before_surgery",
        "time_to_death_disagrees_with_death_minus_candidate_surgery",
    )
    queries = []
    for patient, row in sorted(inventory.clinical.items()):
        profile = _profile({patient: row})
        codes = [key for key in query_counters if profile["arithmetic_diagnostics"][key]]
        if profile["unrecognized_status"]:
            codes.append("unrecognized_status")
        for field, states in profile["field_states"].items():
            if states["invalid"]:
                codes.append(f"invalid_field:{field}")
            if states.get("negative", 0):
                codes.append(f"negative_duration:{field}")
        if codes:
            queries.append(
                {
                    "patient_id": patient,
                    "locator": f"csv/data-row/{inventory.clinical_source_rows[patient]}",
                    "query_codes": codes,
                }
            )
    return {
        "schema_version": "1.0",
        "privacy": "PRIVATE_DO_NOT_COMMIT",
        "clinical_source_sha256": next(
            receipt.sha256 for receipt in inventory.receipts if receipt.kind == "clinical"
        ),
        "scope": "queries_only_not_clinical_exclusions_or_validated_outcomes",
        "records": queries,
    }
