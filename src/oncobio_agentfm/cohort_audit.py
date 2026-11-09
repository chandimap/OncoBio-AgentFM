"""Patient-level modality overlap, selection evidence and private source lineage."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from dataclasses import asdict
from itertools import product
from pathlib import Path

from oncobio_agentfm.cohort_feasibility import endpoint_feasibility, private_endpoint_queries
from oncobio_agentfm.cohort_selection import compare_selection
from oncobio_agentfm.cohort_sources import (
    SOURCE_URLS,
    SourceInventory,
    implementation_fingerprint,
    load_inventory,
    require,
)
from oncobio_agentfm.contracts import MissingnessReason
from oncobio_agentfm.prognosis import ProtocolViolation
from oncobio_agentfm.protocol_audit import write_report
from oncobio_agentfm.source_audit import SOURCE_MISSING_TOKENS
from oncobio_agentfm.splitting import split_patient_ids

MODALITIES = ("ct", "pet", "clinical", "rna")
MISSINGNESS_FIELDS = (
    "CT Date",
    "PET Date",
    "Date of Death",
    "Date of Last Known Alive",
    "EGFR mutation status",
    "KRAS mutation status",
    "ALK translocation status",
    "Pathological T stage",
    "Pathological N stage",
    "Pathological M stage",
)


def source_sets(inventory: SourceInventory) -> dict[str, frozenset[str]]:
    """PT is PET evidence; paired same-study CT is audited independently."""

    return {
        "ct": frozenset(item.patient_id for item in inventory.series if item.modality == "CT"),
        "pet": frozenset(item.patient_id for item in inventory.series if item.modality == "PT"),
        "clinical": frozenset(inventory.clinical),
        "rna": frozenset(inventory.rna_samples.values())
        & frozenset(
            patient
            for patient, observed in inventory.expression_observed_per_patient.items()
            if observed > 0
        ),
    }


def _label_counts(inventory: SourceInventory, ids: frozenset[str]) -> dict[str, int]:
    """Count release vital-status labels, without manufacturing survival endpoints."""

    counts = Counter(
        inventory.clinical[patient]["Survival Status"]
        for patient in sorted(ids & inventory.clinical.keys())
    )
    return {
        "patients": len(ids),
        "clinical_rows": len(ids & inventory.clinical.keys()),
        "recorded_dead": counts["Dead"],
        "recorded_alive": counts["Alive"],
        "unrecognized_status": sum(counts.values()) - counts["Dead"] - counts["Alive"],
    }


def availability_report(inventory: SourceInventory) -> dict[str, object]:
    """Project validated source inventory to aggregates only; never serialize raw rows."""

    sets = source_sets(inventory)
    patients = inventory.patient_ids
    complete = frozenset.intersection(*(sets[modality] for modality in MODALITIES))
    patterns = Counter(
        tuple(patient in sets[modality] for modality in MODALITIES) for patient in patients
    )
    studies: dict[str, set[str]] = {}
    for item in inventory.series:
        studies.setdefault(item.study_uid, set()).add(item.modality)
    paired_pet = frozenset(
        item.patient_id
        for item in inventory.series
        if item.modality == "PT" and "CT" in studies[item.study_uid]
    )
    populations = {
        "ct_clinical": sets["ct"] & sets["clinical"],
        "ct_pet_clinical": sets["ct"] & sets["pet"] & sets["clinical"],
        "imaging_rna": (sets["ct"] | sets["pet"]) & sets["rna"],
        "all_four": complete,
        "same_study_pet_ct_clinical": paired_pet & sets["clinical"],
        "same_study_pet_ct_clinical_rna": paired_pet & sets["clinical"] & sets["rna"],
    }
    raw_counts = {}
    for name, ids in populations.items():
        split = split_patient_ids(ids) if ids else None
        raw_counts[name] = {
            **_label_counts(inventory, ids),
            "partitions": {
                partition: _label_counts(inventory, frozenset(getattr(split, partition)))
                if split
                else _label_counts(inventory, frozenset())
                for partition in ("train", "validation", "test")
            },
        }
    multi = {}
    for modality in ("CT", "PT"):
        counts = Counter(item.patient_id for item in inventory.series if item.modality == modality)
        multi[modality] = sum(count > 1 for count in counts.values())
    field_tokens = {}
    for column in MISSINGNESS_FIELDS:
        counts = Counter(row[column] for row in inventory.clinical.values())
        field_tokens[column] = {
            "blank": counts[""],
            "source_not_collected": counts["Not Collected"] + counts["Not collected"],
            "source_unknown": counts["Unknown"],
            "source_NA": counts["N/A"],
            "other_nonmissing_or_unreviewed": sum(counts.values())
            - sum(counts[token] for token in SOURCE_MISSING_TOKENS),
        }
    return {
        "schema_version": "1.2",
        "evidence_type": "observed_public_source_inventory",
        "implementation_sha256": implementation_fingerprint(),
        "source_receipts": [
            {**asdict(receipt), "reference": SOURCE_URLS[receipt.kind]}
            for receipt in inventory.receipts
        ],
        "denominator": "union_of_exact_patient_identifiers_across_locked_sources",
        "patients": len(patients),
        "modality_definitions": {
            "ct": "at_least_one_nonempty_NBIA_CT_series_metadata_record",
            "pet": "at_least_one_nonempty_NBIA_PT_series_metadata_record",
            "clinical": "one_unique_row_in_reviewed_TCIA_clinical_CSV",
            "rna": "GSE103584_sample_title_and_expression_column_with_observed_numeric_cells",
        },
        "availability_patterns": [
            {"present": dict(zip(MODALITIES, pattern, strict=True)), "patients": patterns[pattern]}
            for pattern in product((False, True), repeat=len(MODALITIES))
        ],
        "modality_counts": {modality: len(ids) for modality, ids in sets.items()},
        "missingness": {
            modality: {
                "absent_from_observed_inventory": len(patients - ids),
                "clinical_reason": "unknown",
                "source_nonascertainment_is_not": "not_acquired",
            }
            for modality, ids in sets.items()
        },
        "source_field_tokens": field_tokens,
        "source_field_token_scope": (
            "recorded_field_states_only_not_patient_modality_reasons; "
            "death_date_NA_in_living_patients_is_not_a_missing_OS_endpoint"
        ),
        "date_field_availability_discordance": {
            "clinical_PET_date_missing_but_PT_in_inventory": sum(
                row["PET Date"] in SOURCE_MISSING_TOKENS and patient in sets["pet"]
                for patient, row in inventory.clinical.items()
            ),
            "date_alignment": "NOT_ASSESSED_SOURCE_CALENDARS_NOT_ASSUMED_COMPARABLE",
        },
        "linkage": {
            "rule": "exact_identifier_equality_no_normalization_or_fuzzy_matching",
            "imaging_without_clinical": len((sets["ct"] | sets["pet"]) - sets["clinical"]),
            "rna_without_clinical": len(sets["rna"] - sets["clinical"]),
            "rna_samples": len(inventory.rna_samples),
            "rna_unique_patients": len(set(inventory.rna_samples.values())),
            "unresolved_cross_source_aliases": "NOT_ASSESSED_NO_IDENTITY_REPAIR_PERFORMED",
        },
        "imaging_reconciliation": {
            "canonical_series": len(inventory.series),
            "studies": len(studies),
            "exact_duplicate_rows_collapsed": inventory.duplicate_series_rows,
            "modality_series_counts": dict(
                sorted(Counter(s.modality for s in inventory.series).items())
            ),
            "patients_with_multiple_candidate_series": multi,
            "patients_with_PT_and_CT_in_same_study": len(paired_pet),
            "series_missing_series_date": sum(s.series_date is None for s in inventory.series),
            "acquisition_dates": "NOT_ESTABLISHED_BY_SERIESDATE_0008_0021",
            "pixel_duplicates": "NOT_ASSESSED_METADATA_ONLY",
            "diagnostic_CT_selection": "NOT_ADJUDICATED_NO_SERIES_SELECTED",
            "clinical_result_availability": "NOT_ESTABLISHED_BY_DICOM_OR_REPOSITORY_DATES",
        },
        "rna_matrix_integrity": {
            "gene_rows": inventory.expression_rows,
            "patient_columns": len(inventory.rna_expression_patients),
            "missing_expression_cells": inventory.expression_missing_cells,
            "total_expression_cells": inventory.expression_rows
            * len(inventory.rna_expression_patients),
            "columns_with_no_observed_numeric_cells": sum(
                count == 0 for count in inventory.expression_observed_per_patient.values()
            ),
            "clinical_result_availability": "NOT_ESTABLISHED",
            "preprocessing_and_biological_quality": "NOT_VALIDATED_BY_MATRIX_PARSE",
        },
        "source_population_label_audit": raw_counts,
        "label_audit_scope": "release_status_labels_not_adjudicated_OS_events_or_horizon_events",
        "selection": compare_selection(inventory, complete),
        "endpoint_feasibility": endpoint_feasibility(inventory),
        "analysis_populations": {
            setting: {
                "status": "BLOCKED_PENDING_SOURCE_ADJUDICATION",
                "eligible_patients": None,
                "survival_events": None,
            }
            for setting in ("pretreatment", "postoperative_day_30")
        },
        "clinical_eligibility": "NOT_ESTABLISHED",
        "model_evaluation": "NOT_RUN",
        "precision_status": "NOT_ESTABLISHED",
        "privacy_scope": "aggregate_only_no_patient_identifiers_UIDs_patient_dates_or_free_text",
    }


def private_manifest(inventory: SourceInventory) -> dict[str, object]:
    """Emit explicit source coordinates locally; this artifact must never be committed."""

    sets = source_sets(inventory)
    hashes = {receipt.kind: receipt.sha256 for receipt in inventory.receipts}
    split = split_patient_ids(inventory.patient_ids)
    partitions = {
        patient: name
        for name in ("train", "validation", "test")
        for patient in getattr(split, name)
    }
    clinical_offsets = inventory.clinical_source_rows
    expression_offsets = inventory.expression_source_columns

    def ref(kind: str, locator: str) -> dict[str, str]:
        return {"kind": kind, "source_sha256": hashes[kind], "locator": locator}

    records = []
    for patient in sorted(inventory.patient_ids):
        links = []
        if patient in clinical_offsets:
            links.append(
                {
                    "role": "clinical",
                    "ref": ref("clinical", f"csv/data-row/{clinical_offsets[patient]}"),
                }
            )
        for series in inventory.series:
            if series.patient_id == patient:
                links.append(
                    {
                        "role": series.modality,
                        "ref": ref("imaging", f"/{series.source_row}"),
                        "study_uid": series.study_uid,
                        "series_uid": series.series_uid,
                        "series_date_on_source_axis": series.series_date,
                        "acquired_on_source_axis": None,
                        "clinical_available_on": None,
                    }
                )
        for accession, owner in sorted(inventory.rna_samples.items()):
            if owner == patient:
                links.append(
                    {
                        "role": "rna_metadata",
                        "ref": ref("rna_metadata", f"sample/{accession}"),
                        "clinical_available_on": None,
                    }
                )
        if patient in expression_offsets:
            links.append(
                {
                    "role": "rna_expression",
                    "ref": ref("rna_expression", f"tsv/column/{expression_offsets[patient]}"),
                    "clinical_available_on": None,
                }
            )
        records.append(
            {
                "patient_id": patient,
                "partition": partitions[patient],
                "observed_presence": {kind: patient in ids for kind, ids in sets.items()},
                "missingness": {
                    kind: MissingnessReason.UNKNOWN.value
                    for kind, ids in sets.items()
                    if patient not in ids
                },
                "source_links": links,
                "clinical_eligibility": "NOT_ADJUDICATED",
            }
        )
    return {
        "schema_version": "1.1",
        "privacy": "PRIVATE_DO_NOT_COMMIT",
        "implementation_sha256": implementation_fingerprint(),
        "source_receipts": [asdict(receipt) for receipt in inventory.receipts],
        "split_seed": "oncobio-agentfm-v1",
        "records": records,
    }


def flow_markdown(report: dict[str, object]) -> str:
    """Draw an aggregate branch with mutually exclusive completeness groups."""

    total = report["patients"]
    complete = report["source_population_label_audit"]["all_four"]["patients"]
    return (
        "# Observed source inventory flow\n\n"
        "Counts describe repository evidence, not a clinically eligible prognostic cohort.\n\n"
        "```mermaid\nflowchart TD\n"
        f'    A["Exact linked source union: {total}"] --> B["All four source types: {complete}"]\n'
        f'    A --> C["At least one source type absent: {total - complete}"]\n'
        '    B --> D["Clinical eligibility and OS events: not established"]\n'
        "    C --> D\n```\n\n"
        "No participant is excluded from the inventory for missing optional evidence. "
        "Clinical inclusion/exclusion counts remain unassessed pending source adjudication.\n"
    )


def main() -> int:
    """Write local audits; exit 2 explicitly distinguishes inventory success from readiness."""

    parser = argparse.ArgumentParser(description=__doc__)
    for kind in SOURCE_URLS:
        parser.add_argument(f"--{kind.replace('_', '-')}", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--local-manifest", action="store_true")
    args = parser.parse_args()
    destination = args.output_dir.resolve()
    require(
        destination.is_relative_to((Path.cwd() / "outputs").resolve()),
        "audit_output_must_be_under_local_outputs",
    )
    inventory = load_inventory({kind: getattr(args, kind) for kind in SOURCE_URLS}, args.lock)
    report = availability_report(inventory)
    destination.mkdir(parents=True, exist_ok=False)
    write_report(destination / "source_availability.json", report)
    with (destination / "source_flow.md").open("x", encoding="utf-8", newline="\n") as stream:
        stream.write(flow_markdown(report))
    if args.local_manifest:
        write_report(destination / "PRIVATE_cohort_manifest.json", private_manifest(inventory))
        write_report(
            destination / "PRIVATE_endpoint_queries.json", private_endpoint_queries(inventory)
        )
    print(json.dumps({"inventory": "AUDITED", "clinical_eligibility": "NOT_ESTABLISHED"}))
    return 2


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ProtocolViolation as error:
        raise SystemExit(error.code) from None
    except (OSError, EOFError, UnicodeError, csv.Error, json.JSONDecodeError):
        raise SystemExit("audit_source_unreadable") from None
