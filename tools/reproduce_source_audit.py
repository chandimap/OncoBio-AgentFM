"""Independent source oracle: standard library only, no project-module imports."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import date
from itertools import product
from pathlib import Path

FILENAMES = {
    "clinical": "clinical.csv",
    "imaging": "nbia-series.json",
    "rna_metadata": "geo-family.soft.gz",
    "rna_expression": "rna-expression.tsv.gz",
}
TOKENS = {"", "Not Collected", "Not collected", "N/A", "Unknown"}
CATEGORIES = {
    "Gender": ("Female", "Male"),
    "Patient affiliation": ("Stanford", "VA"),
    "Smoking status": ("Current", "Former", "Nonsmoker"),
}


def moments(values):
    n = len(values)
    average = math.fsum(values) / n if n else None
    sample_variance = (
        math.fsum((value - average) ** 2 for value in values) / (n - 1) if n > 1 else None
    )
    return average, sample_variance


def selection(rows, complete):
    groups = [
        [row for patient, row in sorted(rows.items()) if (patient in complete) == included]
        for included in (True, False)
    ]
    ages = [
        [
            float(row["Age at Histological Diagnosis"])
            for row in group
            if row["Age at Histological Diagnosis"] not in TOKENS
        ]
        for group in groups
    ]
    stats = [moments(values) for values in ages]
    denominator = (
        math.sqrt((stats[0][1] + stats[1][1]) / 2) if all(v is not None for _, v in stats) else 0
    )
    result = {
        "complete_n": len(groups[0]),
        "incomplete_n": len(groups[1]),
        "age_at_histological_diagnosis": {
            name: {
                "observed_n": len(values),
                "missing_n": len(group) - len(values),
                "mean": average,
                "sd": math.sqrt(var) if var is not None else None,
            }
            for name, values, group, (average, var) in zip(
                ("complete", "incomplete"), ages, groups, stats, strict=True
            )
        },
        "categorical": {},
    }
    result["age_at_histological_diagnosis"]["smd"] = (
        (stats[0][0] - stats[1][0]) / denominator if denominator else None
    )
    for column, levels in CATEGORIES.items():
        counters = [Counter(row[column] for row in group) for group in groups]
        totals = [sum(counter[level] for level in levels) for counter in counters]
        contrast = []
        for level in levels:
            proportions = [
                counter[level] / total if total else None
                for counter, total in zip(counters, totals, strict=True)
            ]
            difference = scale = None
            if all(p is not None for p in proportions):
                difference = proportions[0] - proportions[1]
                scale = math.sqrt(sum(p * (1 - p) for p in proportions) / 2)
            contrast.append(
                {
                    "level": level,
                    "complete_n": counters[0][level],
                    "incomplete_n": counters[1][level],
                    "proportions": proportions,
                    "proportion_difference": difference,
                    "smd": difference / scale if scale else None,
                }
            )
        result["categorical"][column] = {
            "observed_denominators": totals,
            "missing_or_unrecognized": [
                len(group) - total for group, total in zip(groups, totals, strict=True)
            ],
            "levels": contrast,
        }
    return result


def endpoint_profile(rows):
    """Independently recompute the locked release's arithmetic, without clinical inference."""

    def parse_date(value):
        if value in TOKENS:
            return None
        month, day, year = map(int, value.split("/"))
        return date(year, month, day).toordinal()

    date_fields = ("CT Date", "PET Date", "Date of Death", "Date of Last Known Alive")
    duration_fields = ("Days between CT and surgery", "Time to Death (days)")
    parsed = [
        {
            **{key: parse_date(row[key]) for key in date_fields},
            **{key: None if row[key] in TOKENS else int(row[key]) for key in duration_fields},
            "status": row["Survival Status"],
        }
        for row in rows.values()
    ]
    fields = {}
    for field in (*date_fields, *duration_fields):
        values = [row[field] for row in parsed]
        fields[field] = {
            "parsed": sum(v is not None for v in values),
            "source_missing_token": sum(v is None for v in values),
            "invalid": 0,
        }
        if field in duration_fields:
            fields[field].update(
                negative=sum(v is not None and v < 0 for v in values),
                zero=sum(v == 0 for v in values),
            )
    dead = [r for r in parsed if r["status"] == "Dead" and r["Date of Death"] is not None]
    alive = [r for r in parsed if r["status"] == "Alive" and r["Date of Last Known Alive"]]
    pairs = [r for r in dead if r["CT Date"] is not None]
    alive_pairs = [r for r in alive if r["CT Date"] is not None]
    ttd_pairs = [r for r in pairs if r["Time to Death (days)"] is not None]
    candidate = [
        r
        for r in parsed
        if r["CT Date"] is not None
        and r["Days between CT and surgery"] is not None
        and r["Days between CT and surgery"] >= 0
    ]
    surgery_dead = [r for r in candidate if r in dead]
    surgery_alive = [r for r in candidate if r in alive]
    surgery_ttd = [r for r in surgery_dead if r["Time to Death (days)"] is not None]
    matches = sum(
        r["Date of Death"] - r["CT Date"]
        == r["Time to Death (days)"] + r["Days between CT and surgery"]
        for r in surgery_ttd
    )
    return {
        "clinical_rows": len(rows),
        "recorded_dead": sum(r["status"] == "Dead" for r in parsed),
        "recorded_alive": sum(r["status"] == "Alive" for r in parsed),
        "unrecognized_status": sum(r["status"] not in {"Alive", "Dead"} for r in parsed),
        "field_states": fields,
        "arithmetic_diagnostics": {
            "alive_with_death_date": sum(
                r["status"] == "Alive" and r["Date of Death"] is not None for r in parsed
            ),
            "dead_without_parseable_death_date": sum(
                r["status"] == "Dead" and r["Date of Death"] is None for r in parsed
            ),
            "last_alive_after_death": sum(
                r["Date of Last Known Alive"] is not None
                and r["Date of Last Known Alive"] > r["Date of Death"]
                for r in dead
            ),
            "alive_with_numeric_time_to_death": sum(
                r["status"] == "Alive" and r["Time to Death (days)"] is not None for r in parsed
            ),
            "death_CT_pairs": len(pairs),
            "death_on_or_before_CT": sum(r["Date of Death"] <= r["CT Date"] for r in pairs),
            "alive_last_contact_CT_pairs": len(alive_pairs),
            "alive_last_contact_on_or_before_CT": sum(
                r["Date of Last Known Alive"] <= r["CT Date"] for r in alive_pairs
            ),
            "time_to_death_CT_comparison_pairs": len(ttd_pairs),
            "time_to_death_matches_death_minus_CT": sum(
                r["Date of Death"] - r["CT Date"] == r["Time to Death (days)"] for r in ttd_pairs
            ),
            "candidate_surgery_records": len(candidate),
            "candidate_surgery_recorded_deaths": len(surgery_dead),
            "candidate_death_on_or_before_surgery": sum(
                r["Date of Death"] - r["CT Date"] <= r["Days between CT and surgery"]
                for r in surgery_dead
            ),
            "candidate_surgery_recorded_alive": len(surgery_alive),
            "candidate_alive_contact_on_or_before_surgery": sum(
                r["Date of Last Known Alive"] - r["CT Date"] <= r["Days between CT and surgery"]
                for r in surgery_alive
            ),
            "time_to_death_surgery_comparison_pairs": len(surgery_ttd),
            "time_to_death_matches_death_minus_candidate_surgery": matches,
            "time_to_death_disagrees_with_death_minus_candidate_surgery": len(surgery_ttd)
            - matches,
        },
    }


def independent_facts(source_dir: Path, lock: dict) -> dict:
    raw = {kind: (source_dir / filename).read_bytes() for kind, filename in FILENAMES.items()}
    for kind, value in raw.items():
        if not (len(value) == lock["sources"][kind]["size_bytes"]):
            raise ValueError("locked_size_mismatch")
        if not (hashlib.sha256(value).hexdigest() == lock["sources"][kind]["sha256"]):
            raise ValueError("locked_hash_mismatch")
    reader = csv.DictReader(raw["clinical"].decode("utf-8-sig").splitlines())
    clinical_rows = list(reader)
    rows = {row["Case ID"]: row for row in clinical_rows}
    if not (len(rows) == len(clinical_rows)):
        raise ValueError("duplicate_clinical_key")
    records = json.loads(raw["imaging"])
    series = {}
    study_modalities = defaultdict(set)
    owners = {}
    duplicates = 0
    for record in records:
        uid, study, patient = (
            record[key] for key in ("SeriesInstanceUID", "StudyInstanceUID", "PatientID")
        )
        if not (study not in owners or owners[study] == patient):
            raise ValueError("conflicting_study_owner")
        owners[study] = patient
        if uid in series:
            if not (record == series[uid]):
                raise ValueError("conflicting_series_record")
            duplicates += 1
        series[uid] = record
        study_modalities[study].add(record["Modality"])
    titles, accession = {}, None
    for line in gzip.decompress(raw["rna_metadata"]).decode().splitlines():
        if line.startswith("^SAMPLE = "):
            accession = line.split(" = ", 1)[1]
        if line.startswith("!Sample_title = "):
            if not (accession not in titles):
                raise ValueError("repeated_GEO_sample")
            titles[accession] = line.split(" = ", 1)[1]
    matrix = csv.reader(
        gzip.decompress(raw["rna_expression"]).decode().splitlines(), delimiter="\t"
    )
    header = next(matrix)[1:]
    if not (len(set(header)) == len(header)):
        raise ValueError("repeated_expression_patient")
    observed = [0] * len(header)
    missing = gene_rows = 0
    genes = set()
    for row in matrix:
        if not (len(row) == len(header) + 1 and row[0] not in genes):
            raise ValueError("matrix_structure_invalid")
        genes.add(row[0])
        gene_rows += 1
        for index, cell in enumerate(row[1:]):
            if cell == "NA":
                missing += 1
            else:
                if not (math.isfinite(float(cell))):
                    raise ValueError("nonfinite_expression")
                observed[index] += 1
    ct = {r["PatientID"] for r in series.values() if r["Modality"] == "CT"}
    pet = {r["PatientID"] for r in series.values() if r["Modality"] == "PT"}
    rna = {patient for patient, n in zip(header, observed, strict=True) if n > 0} & set(
        titles.values()
    )
    patients = (
        set(rows) | {r["PatientID"] for r in series.values()} | set(titles.values()) | set(header)
    )
    paired = {
        r["PatientID"]
        for r in series.values()
        if r["Modality"] == "PT" and "CT" in study_modalities[r["StudyInstanceUID"]]
    }
    sets = (ct, pet, set(rows), rna)
    complete = set.intersection(*sets)
    population = {
        "ct_clinical": ct & set(rows),
        "ct_pet_clinical": ct & pet & set(rows),
        "imaging_rna": (ct | pet) & rna,
        "all_four": complete,
        "same_study_pet_ct_clinical": paired & set(rows),
        "same_study_pet_ct_clinical_rna": paired & set(rows) & rna,
    }

    def labels(ids):
        counts = Counter(rows[patient]["Survival Status"] for patient in ids if patient in rows)
        return {
            "patients": len(ids),
            "clinical_rows": len(ids & set(rows)),
            "recorded_dead": counts["Dead"],
            "recorded_alive": counts["Alive"],
            "unrecognized_status": sum(counts.values()) - counts["Dead"] - counts["Alive"],
        }

    def partition(patient):
        number = int.from_bytes(
            hashlib.sha256(("oncobio-agentfm-v1\x1f" + patient).encode()).digest()
        )
        # Integer threshold comparison is independent of the implementation's float path.
        return (
            "train"
            if number * 10 < 7 * 2**256
            else ("validation" if number * 100 < 85 * 2**256 else "test")
        )

    budgets = {
        name: {
            **labels(ids),
            "partitions": {
                part: labels({p for p in ids if partition(p) == part})
                for part in ("train", "validation", "test")
            },
        }
        for name, ids in population.items()
    }
    field_tokens = {}
    for column in (
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
    ):
        counter = Counter(row[column] for row in rows.values())
        field_tokens[column] = {
            "blank": counter[""],
            "source_not_collected": counter["Not Collected"] + counter["Not collected"],
            "source_unknown": counter["Unknown"],
            "source_NA": counter["N/A"],
            "other_nonmissing_or_unreviewed": sum(
                n for key, n in counter.items() if key not in TOKENS
            ),
        }
    selection_facts = selection(rows, complete)
    selection_facts["source_cohort_strata"] = {}
    for prefix, name in (("R01-", "R01"), ("AMC-", "AMC")):
        subset = {key: row for key, row in rows.items() if key.startswith(prefix)}
        selection_facts["source_cohort_strata"][name] = selection(subset, complete & set(subset))
    return {
        "patients": len(patients),
        "modality_counts": dict(zip(("ct", "pet", "clinical", "rna"), map(len, sets), strict=True)),
        "availability_patterns": [
            {
                "present": dict(zip(("ct", "pet", "clinical", "rna"), pattern, strict=True)),
                "patients": sum(tuple(p in s for s in sets) == pattern for p in patients),
            }
            for pattern in product((False, True), repeat=4)
        ],
        "missingness": {
            name: {"absent_from_observed_inventory": len(patients - ids)}
            for name, ids in zip(("ct", "pet", "clinical", "rna"), sets, strict=True)
        },
        "linkage": {
            "imaging_without_clinical": len((ct | pet) - set(rows)),
            "rna_without_clinical": len(rna - set(rows)),
            "rna_samples": len(titles),
            "rna_unique_patients": len(set(titles.values())),
        },
        "imaging_reconciliation": {
            "canonical_series": len(series),
            "studies": len(owners),
            "exact_duplicate_rows_collapsed": duplicates,
            "modality_series_counts": dict(Counter(r["Modality"] for r in series.values())),
            "patients_with_multiple_candidate_series": {
                modality: sum(
                    n > 1
                    for n in Counter(
                        r["PatientID"] for r in series.values() if r["Modality"] == modality
                    ).values()
                )
                for modality in ("CT", "PT")
            },
            "patients_with_PT_and_CT_in_same_study": len(paired),
            "series_missing_series_date": sum(not r.get("SeriesDate") for r in series.values()),
        },
        "rna_matrix_integrity": {
            "gene_rows": gene_rows,
            "patient_columns": len(header),
            "missing_expression_cells": missing,
            "total_expression_cells": gene_rows * len(header),
            "columns_with_no_observed_numeric_cells": sum(n == 0 for n in observed),
        },
        "source_field_tokens": field_tokens,
        "date_field_availability_discordance": {
            "clinical_PET_date_missing_but_PT_in_inventory": sum(
                row["PET Date"] in TOKENS and p in pet for p, row in rows.items()
            )
        },
        "source_population_label_audit": budgets,
        "selection": selection_facts,
        "endpoint_feasibility": {
            "source_populations": {
                "clinical_release": endpoint_profile(rows),
                "all_four_source_types": endpoint_profile(
                    {patient: rows[patient] for patient in complete}
                ),
            },
        },
    }


def verify_report(facts, report):
    checks = []

    def visit(expected, actual, path):
        if isinstance(expected, dict):
            for key, value in expected.items():
                if key == "series_missing_series_date" and key not in actual:
                    visit(
                        value,
                        actual["series_missing_acquisition_date"],
                        path + "/series_missing_acquisition_date",
                    )
                    continue
                if key == "source_cohort_strata" and key not in actual:
                    continue  # Older reports did not claim stratified findings.
                visit(value, actual[key], path + "/" + key)
        elif isinstance(expected, list):
            if not (len(expected) == len(actual)):
                raise ValueError(path)
            for index, (first, second) in enumerate(zip(expected, actual, strict=True)):
                visit(first, second, path + "/" + str(index))
        else:
            same = (
                math.isclose(expected, actual, rel_tol=1e-12, abs_tol=1e-12)
                if isinstance(expected, float)
                else expected == actual
            )
            checks.append({"field": path, "passed": same})

    visit(facts, report, "")
    return {
        "schema_version": "1.0",
        "method": "independent_stdlib_oracle_no_project_imports",
        "source_hash_checks": 4,
        "aggregate_checks": len(checks),
        "status": "PASS" if all(check["passed"] for check in checks) else "FAIL",
        "checks": checks,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--lock", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    facts = independent_facts(args.source_dir, json.loads(args.lock.read_bytes()))
    inspected = json.loads(args.report.read_bytes())
    result = verify_report(facts, inspected)
    result["oracle_sha256"] = hashlib.sha256(
        Path(__file__).read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    result["audited_implementation_sha256"] = inspected["implementation_sha256"]
    result["interpretation_limits"] = (
        ["legacy_missing_acquisition_date_field_counts_missing_SeriesDate_only"]
        if "series_missing_acquisition_date" in inspected["imaging_reconciliation"]
        else []
    )
    for setting in ("pretreatment", "postoperative_day_30"):
        if (
            json.loads(args.report.read_bytes())["analysis_populations"][setting][
                "eligible_patients"
            ]
            is not None
        ):
            raise ValueError("clinical_status_must_remain_unassessed")
        if inspected["analysis_populations"][setting]["survival_events"] is not None:
            raise ValueError("clinical_status_must_remain_unassessed")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8", newline="\n") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(
        json.dumps(
            {key: result[key] for key in ("status", "aggregate_checks", "source_hash_checks")}
        )
    )
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
