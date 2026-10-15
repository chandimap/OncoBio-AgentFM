"""Offline source-schema and release-budget inspection; no clinical auto-import."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
from collections import Counter
from pathlib import Path

from oncobio_agentfm.protocol_audit import write_report
from oncobio_agentfm.splitting import split_patient_ids

SOURCE_REFERENCE = (
    "https://www.cancerimagingarchive.net/wp-content/uploads/"
    "NSCLCR01Radiogenomic_DATA_LABELS_2018-05-22_1500-shifted.csv"
)
REVIEWED_SOURCE_COLUMNS = (
    "Case ID",
    "Patient affiliation",
    "Age at Histological Diagnosis",
    "Weight (lbs)",
    "Gender",
    "Ethnicity",
    "Smoking status",
    "Pack Years",
    "Quit Smoking Year",
    "%GG",
    "Tumor Location (choice=RUL)",
    "Tumor Location (choice=RML)",
    "Tumor Location (choice=RLL)",
    "Tumor Location (choice=LUL)",
    "Tumor Location (choice=LLL)",
    "Tumor Location (choice=L Lingula)",
    "Tumor Location (choice=Unknown)",
    "Histology ",
    "Pathological T stage",
    "Pathological N stage",
    "Pathological M stage",
    "Histopathological Grade",
    "Lymphovascular invasion",
    "Pleural invasion (elastic, visceral, or parietal)",
    "EGFR mutation status",
    "KRAS mutation status",
    "ALK translocation status",
    "Adjuvant Treatment",
    "Chemotherapy",
    "Radiation",
    "Recurrence",
    "Recurrence Location",
    "Date of Recurrence",
    "Date of Last Known Alive",
    "Survival Status",
    "Date of Death",
    "Time to Death (days)",
    "CT Date",
    "Days between CT and surgery",
    "PET Date",
)
PATHOLOGICAL_COLUMNS = ("Pathological T stage", "Pathological N stage", "Pathological M stage")
STAGE_TOKENS = {
    "Pathological T stage": re.compile(r"T(?:is|[0-4][abc]?)\Z"),
    "Pathological N stage": re.compile(r"N[0-3][abc]?\Z"),
    "Pathological M stage": re.compile(r"M[01][abc]?\Z"),
}
SOURCE_MISSING_TOKENS = frozenset({"", "Not Collected", "Not collected", "N/A", "Unknown"})


def audit_source_schema(path: Path) -> dict[str, object]:
    """Report source metadata and aggregate release counts without patient dates, IDs or rows."""

    raw = path.read_bytes()
    reader = csv.DictReader(io.StringIO(raw.decode("utf-8-sig")), strict=True)
    headers = reader.fieldnames
    if not headers or len(set(headers)) != len(headers):
        raise ValueError("clinical CSV requires unique headers")
    if tuple(headers) != REVIEWED_SOURCE_COLUMNS:
        raise ValueError("clinical CSV does not match the reviewed TCIA header schema")
    rows = list(reader)
    if not rows or any(None in row or any(value is None for value in row.values()) for row in rows):
        raise ValueError("clinical CSV is empty or has inconsistent row widths")
    ids = [row["Case ID"] for row in rows]
    if any(not identifier.strip() or identifier != identifier.strip() for identifier in ids):
        raise ValueError("clinical CSV contains invalid patient identifiers")
    if len(set(ids)) != len(ids):
        raise ValueError("clinical CSV contains duplicate patient identifiers")
    stage_vocabulary = {
        column: sorted(
            {
                value
                if value in SOURCE_MISSING_TOKENS or STAGE_TOKENS[column].fullmatch(value)
                else "<unrecognized-token-withheld>"
                for value in (row[column] for row in rows)
            }
        )
        for column in PATHOLOGICAL_COLUMNS
        if column in headers
    }
    # These are administrative release counts, never adjudicated survival outcomes.
    status_by_id = {
        row["Case ID"]: {"Dead": "dead", "Alive": "alive"}.get(
            row["Survival Status"], "unrecognized"
        )
        for row in rows
    }
    split = split_patient_ids(ids)
    raw_partitions = {}
    for name, identifiers in (
        ("train", split.train),
        ("validation", split.validation),
        ("test", split.test),
    ):
        statuses = Counter(status_by_id[identifier] for identifier in identifiers)
        raw_partitions[name] = {
            "subjects": len(identifiers),
            "recorded_dead": statuses["dead"],
            "recorded_alive": statuses["alive"],
            "unrecognized_status": statuses["unrecognized"],
        }
    statuses = Counter(status_by_id.values())
    return {
        "schema_version": "1.2",
        "evidence_type": "source_schema_inspection_only",
        "source_reference": SOURCE_REFERENCE,
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "source_bytes": len(raw),
        "rows": len(rows),
        "columns": headers,
        "raw_pathological_stage_vocabulary": stage_vocabulary,
        "raw_source_budget": {
            "subjects": len(rows),
            "recorded_dead": statuses["dead"],
            "recorded_alive": statuses["alive"],
            "unrecognized_status": statuses["unrecognized"],
            "nominal_15_percent_subjects": 0.15 * len(rows),
            "nominal_15_percent_recorded_deaths": 0.15 * statuses["dead"],
            "fixed_hash_partition_counts": raw_partitions,
            "split_seed": "oncobio-agentfm-v1",
            "precision_status": "NOT_ESTABLISHED",
            "scope": (
                "Full-release label counts only; not eligible cohort, landmark survivors or "
                "730-day events. Do not alter split, horizon or population in response."
            ),
        },
        "real_data_temporal_eligibility": "NOT_ESTABLISHED",
        "primary_prediction_readiness": "BLOCKED_PENDING_SOURCE_ADJUDICATION",
        "secondary_prediction_readiness": "BLOCKED_PENDING_SOURCE_ADJUDICATION",
        "blockers": [
            "No verified pretreatment assessment and diagnosis-result availability chronology.",
            "CT/PET acquisition dates are not verified result-availability dates.",
            "Pathological TNM cannot supply pretreatment clinical stage; edition is not declared.",
            "Age at histological diagnosis is not automatically age at the prediction date.",
            "No supplied ECOG field; represent absence explicitly, never synthesize a score.",
            "Molecular result availability and patient-level surgery linkage require adjudication.",
            "Shifted image, treatment and follow-up dates need a verified "
            "common patient time axis.",
            "Time to Death (days) cannot be reused without a verified original time origin.",
            "No source-adjudicated administrative-cutoff policy on verified patient axes.",
        ],
        "scope": (
            "Schema inspection cannot verify dates, source meanings or endpoint consistency. "
            "The reference URL identifies the source family, not authenticated provenance "
            "of arbitrary local bytes. Header changes require renewed source review."
        ),
    }


def main() -> int:
    """Write inspectable source limitations; exit 2 means eligibility remains blocked."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit_source_schema(args.input)
    write_report(args.output, report)
    print(json.dumps({"readiness": report["primary_prediction_readiness"], "rows": report["rows"]}))
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
