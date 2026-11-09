"""Deterministic, explicitly synthetic evidence for cohort construction safeguards."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

from oncobio_agentfm.cohort_audit import availability_report, private_manifest, source_sets
from oncobio_agentfm.cohort_linkage import (
    audit_linked_cases,
    case_payload,
    case_references,
    verify_case_review,
    verify_documents,
)
from oncobio_agentfm.cohort_sources import (
    SourceInventory,
    implementation_fingerprint,
    load_inventory,
    read_clinical,
    read_expression,
    read_rna_metadata,
    read_series,
    unique_object,
)
from oncobio_agentfm.contracts import EvidenceRecord, MissingnessReason
from oncobio_agentfm.prognosis import (
    AvailabilityBasis,
    PredictorFeature,
    ProtocolViolation,
    VitalStatus,
    validate_predictors,
)
from oncobio_agentfm.protocol_audit import AdjudicatedCase, write_report
from oncobio_agentfm.protocol_validation import make_synthetic_case
from oncobio_agentfm.research_protocol import (
    FeatureKind,
    PredictionSetting,
    ResearchProtocol,
    load_research_protocol,
)
from oncobio_agentfm.source_audit import REVIEWED_SOURCE_COLUMNS
from oncobio_agentfm.splitting import split_patient_ids


def synthetic_bundle(directory: Path, count: int = 12) -> tuple[dict[str, Path], Path]:
    """Invent source files with known overlap; no real patient material is used."""

    identifiers = [f"SYN-cohort-{number}" for number in range(count)]
    clinical = io.StringIO(newline="")
    writer = csv.DictWriter(clinical, fieldnames=REVIEWED_SOURCE_COLUMNS)
    writer.writeheader()
    for number, patient in enumerate(identifiers):
        writer.writerow(
            {
                "Case ID": patient,
                "Age at Histological Diagnosis": str(50 + number),
                "Gender": "Female" if number % 2 else "Male",
                "Patient affiliation": "Stanford" if number % 3 else "VA",
                "Smoking status": "Former",
                "Survival Status": "Dead" if number % 2 else "Alive",
            }
        )
    series = []
    for number, patient in enumerate(identifiers):
        for modality, prefix in (("CT", 10), ("PT", 30)):
            if modality == "PT" and number % 2:
                continue
            series.append(
                {
                    "PatientID": patient,
                    "StudyInstanceUID": f"1.2.840.20.{number}",
                    "SeriesInstanceUID": f"1.2.840.{prefix}.{number}",
                    "Modality": modality,
                    "ImageCount": "32",
                    "Collection": "NSCLC Radiogenomics",
                    "SeriesDate": "2024-12-22 00:00:00.0",
                    "StudyDate": "2024-12-22 00:00:00.0",
                }
            )
    molecular = [(number, patient) for number, patient in enumerate(identifiers) if number % 3 == 0]
    soft = "^SERIES = GSE103584\n!Series_geo_accession = GSE103584\n"
    for number, patient in molecular:
        soft += (
            f"^SAMPLE = GSM{900000 + number}\n!Sample_title = {patient}\n"
            "!Sample_type = SRA\n!Sample_library_strategy = RNA-Seq\n"
        )
    expression = "\t" + "\t".join(patient for _, patient in molecular) + "\n"
    expression += "SYN-GENE-1\t" + "\t".join("1.5" for _ in molecular) + "\n"
    expression += "SYN-GENE-2\t" + "\t".join("NA" for _ in molecular) + "\n"
    raw = {
        "clinical": clinical.getvalue().encode(),
        "imaging": json.dumps(series).encode(),
        "rna_metadata": gzip.compress(soft.encode(), mtime=0),
        "rna_expression": gzip.compress(expression.encode(), mtime=0),
    }
    directory.mkdir(parents=True, exist_ok=True)
    paths = {}
    lock = {"schema_version": "1.0", "sources": {}}
    for kind, value in raw.items():
        paths[kind] = directory / f"SYNTHETIC-{kind}.bin"
        paths[kind].write_bytes(value)
        lock["sources"][kind] = {
            "version": "synthetic-v1",
            "retrieved_on": "2025-01-01",
            "sha256": hashlib.sha256(value).hexdigest(),
            "size_bytes": len(value),
        }
    lock_path = directory / "SYNTHETIC-lock.json"
    lock_path.write_text(json.dumps(lock), encoding="utf-8")
    return paths, lock_path


def synthetic_linked_cases(
    directory: Path,
    inventory: SourceInventory,
) -> tuple[tuple[AdjudicatedCase, ...], dict[str, tuple[Path, str]], dict[str, Path]]:
    """Invent reviewed clinical assertions and bind image/RNA IDs to synthetic sources."""

    cases = []
    for patient in sorted(inventory.patient_ids):
        case = make_synthetic_case(patient)
        if inventory.clinical[patient]["Survival Status"] == "Alive":
            case = replace(
                case, follow_up=replace(case.follow_up, status=VitalStatus.ALIVE, death_on=None)
            )
        selected_ct = next(
            s for s in inventory.series if s.patient_id == patient and s.modality == "CT"
        )
        ct = replace(
            case.primary_features[0],
            evidence=replace(case.primary_features[0].evidence, source_id=selected_ct.series_uid),
        )
        pet = case.primary_features[1]
        selected_pet = next(
            (s for s in inventory.series if s.patient_id == patient and s.modality == "PT"), None
        )
        if selected_pet is not None:
            pet = PredictorFeature(
                FeatureKind.PET_CT,
                replace(ct.evidence, source_id=selected_pet.series_uid),
                "synthetic-pet-access",
                AvailabilityBasis.CLINICAL_RECORD,
            )
        molecular = case.landmark_features[-1]
        accession = next(
            (sample for sample, owner in inventory.rna_samples.items() if owner == patient), None
        )
        if accession is None:
            molecular = PredictorFeature(
                FeatureKind.RESECTION_MOLECULAR,
                EvidenceRecord(patient, None, None, None, MissingnessReason.NOT_ACQUIRED),
                None,
            )
        else:
            molecular = replace(
                molecular, evidence=replace(molecular.evidence, source_id=accession)
            )
        case = replace(case, primary_features=(ct, pet), landmark_features=(ct, pet, molecular))
        cases.append(case)
    return tuple(cases), *write_synthetic_reviews(directory, inventory, tuple(cases))


def write_synthetic_reviews(directory, inventory, cases):
    """Write invented assertion bundles; never use this helper to adjudicate real sources."""

    bundles = {}
    for case in cases:
        patient = case.context.baseline.patient_id
        ct, pet = case.primary_features
        pairings = {}
        if pet.evidence.is_present:
            pairings[pet.evidence.source_id] = {
                "ct_evidence": case_payload(case)["primary_features"][0]["evidence"],
                "timing_source": ct.timing_source,
                "pairing_source": f"synthetic-pair-review:{patient}",
                "availability_basis": "clinical_record",
            }
        molecular = {}
        for feature in case.landmark_features:
            item = feature.evidence
            if feature.kind is FeatureKind.RESECTION_MOLECULAR and item.is_present:
                surgery = case.context.surgery
                molecular[item.source_id] = {
                    "expression_source_sha256": next(
                        r.sha256 for r in inventory.receipts if r.kind == "rna_expression"
                    ),
                    "expression_locator": (
                        f"tsv/column/{inventory.expression_source_columns[patient]}"
                    ),
                    "resection_source_id": surgery.source_id if surgery else None,
                    "specimen_collected_on": item.acquired_on.isoformat(),
                    "clinical_available_on": item.available_on.isoformat(),
                    "availability_basis": feature.availability_basis.value
                    if feature.availability_basis
                    else None,
                    "timing_source": feature.timing_source,
                }
        bundles[patient] = {
            "case": case_payload(case),
            "pet_ct_pairings": pairings,
            "molecular_provenance": molecular,
        }
    document = directory / "SYNTHETIC-adjudication.json"
    document.write_text(
        json.dumps(
            {
                "schema_version": "1.0",
                "evidence_type": "source_adjudication_assertions",
                "cases": bundles,
            },
            sort_keys=True,
        ),
        encoding="utf-8",
    )
    digest = hashlib.sha256(document.read_bytes()).hexdigest()
    reviews = {}
    for case in cases:
        patient = case.context.baseline.patient_id
        bundle = bundles[patient]
        refs = set(case_references(case))
        for pairing in bundle["pet_ct_pairings"].values():
            refs.update(
                (
                    pairing["ct_evidence"]["source_id"],
                    pairing["timing_source"],
                    pairing["pairing_source"],
                )
            )
        references = {
            ref: {"source_sha256": digest, "locator": f"/cases/{patient}"} for ref in refs
        }
        mappings = {}
        for feature in (*case.primary_features, *case.landmark_features):
            item = feature.evidence
            if not item.is_present:
                continue
            if feature.kind in (FeatureKind.CT, FeatureKind.PET_CT):
                source = next(s for s in inventory.series if s.series_uid == item.source_id)
                references[item.source_id] = {
                    "source_sha256": next(
                        r.sha256 for r in inventory.receipts if r.kind == "imaging"
                    ),
                    "locator": f"/{source.source_row}",
                }
                mappings[item.source_id] = {
                    "source_axis_date": source.series_date,
                    "source_date_kind": "SeriesDate_0008_0021",
                    "acquired_on": item.acquired_on.isoformat(),
                    "acquisition_reference": feature.timing_source,
                    "alignment_reference": case.context.calendar_alignment_source,
                }
            if feature.kind is FeatureKind.RESECTION_MOLECULAR:
                references[item.source_id] = {
                    "source_sha256": next(
                        r.sha256 for r in inventory.receipts if r.kind == "rna_metadata"
                    ),
                    "locator": f"sample/{item.source_id}",
                }
        review = {
            "schema_version": "1.1",
            **bundle,
            "references": references,
            "image_date_mappings": mappings,
        }
        path = directory / f"{patient}-PRIVATE-review.json"
        path.write_text(json.dumps(review, sort_keys=True), encoding="utf-8")
        reviews[patient] = (path, hashlib.sha256(path.read_bytes()).hexdigest())
    return reviews, {digest: document}


def validation_report(protocol: ResearchProtocol) -> dict[str, object]:
    """Run prespecified positive and adversarial probes with inspectable fixed outcomes."""

    checks = []

    def check(name, function, expected=None):
        observed = None
        try:
            result = function()
            if result is False:
                observed = "unexpected_result"
        except ProtocolViolation as exc:
            observed = exc.code
        checks.append(
            {
                "name": name,
                "expected_rejection": expected,
                "observed_rejection": observed,
                "passed": observed == expected,
            }
        )

    with TemporaryDirectory(prefix="oncobio-synthetic-") as temporary:
        root = Path(temporary)
        paths, lock = synthetic_bundle(root)
        inventory = load_inventory(paths, lock)
        audit = availability_report(inventory)
        rows = json.loads(paths["imaging"].read_bytes())
        check(
            "known_patient_overlap",
            lambda: audit["modality_counts"] == {"ct": 12, "pet": 6, "clinical": 12, "rna": 4},
        )
        check(
            "known_all_four_overlap",
            lambda: audit["source_population_label_audit"]["all_four"]["patients"] == 2,
        )
        check(
            "exact_duplicate_series_does_not_add_patients",
            lambda: len(read_series(json.dumps([*rows, rows[0]]).encode())[0]) == len(rows),
        )
        changed = {**rows[0], "PatientID": "SYN-other-owner"}
        check(
            "study_cannot_belong_to_two_patients",
            lambda: read_series(json.dumps([*rows, changed]).encode()),
            "imaging_study_owner_conflict",
        )
        changed_series = {**rows[0], "ImageCount": "99"}
        check(
            "same_series_conflicting_metadata_is_blocked",
            lambda: read_series(json.dumps([*rows, changed_series]).encode()),
            "imaging_series_metadata_conflict",
        )
        check(
            "whitespace_identity_is_not_repaired",
            lambda: read_clinical(
                paths["clinical"].read_bytes().replace(b"SYN-cohort-0,", b" SYN-cohort-0,")
            ),
            "source_patient_identifier_invalid",
        )
        check(
            "duplicate_json_keys_are_blocked",
            lambda: json.loads(
                '{"PatientID":"SYN-A","PatientID":"SYN-B"}', object_pairs_hook=unique_object
            ),
            "source_duplicate_json_key",
        )
        check("NA_cells_remain_missing", lambda: inventory.expression_missing_cells == 4)
        check(
            "nonfinite_expression_is_blocked",
            lambda: read_expression(gzip.compress(b"\tSYN-A\nGENE\tinf\n", mtime=0)),
            "rna_expression_nonfinite",
        )
        check(
            "duplicate_RNA_patient_columns_are_blocked",
            lambda: read_expression(gzip.compress(b"\tSYN-A\tSYN-A\nGENE\t1\t2\n", mtime=0)),
            "rna_duplicate_patient_column",
        )
        check(
            "wrong_GEO_series_is_blocked",
            lambda: read_rna_metadata(gzip.compress(b"!Series_geo_accession = GSE0\n", mtime=0)),
            "rna_series_mismatch",
        )
        manifest = private_manifest(inventory)
        check(
            "all_modalities_keep_the_patient_partition",
            lambda: (
                len(manifest["records"]) == 12
                and len({record["patient_id"] for record in manifest["records"]}) == 12
            ),
        )
        split = split_patient_ids(inventory.patient_ids)
        subset = split_patient_ids(source_sets(inventory)["rna"])
        check(
            "subset_never_reassigns_patient_partition",
            lambda: all(
                set(getattr(subset, name)) <= set(getattr(split, name))
                for name in ("train", "validation", "test")
            ),
        )
        cases, reviews, documents = synthetic_linked_cases(root, inventory)
        linked = audit_linked_cases(
            cases, inventory=inventory, reviews=reviews, documents=documents, protocol=protocol
        )
        check(
            "linked_primary_preserves_existing_gate",
            lambda: linked["settings"]["pretreatment"]["endpoint_eligible"] == 12,
        )
        check(
            "linked_all_four_uses_landmark_population",
            lambda: linked["all_four_landmark_patients"] == 2,
        )
        case = cases[0]
        late = replace(
            case.primary_features[0],
            evidence=replace(
                case.primary_features[0].evidence,
                available_on=case.context.baseline.index_date.replace(day=2),
            ),
        )
        check(
            "source_presence_does_not_admit_late_CT",
            lambda: validate_predictors(
                case.context,
                (late, case.primary_features[1]),
                setting=PredictionSetting.PRETREATMENT,
                protocol=protocol,
            ),
            "post_index_evidence",
        )
        path, fingerprint = reviews[case.context.baseline.patient_id]
        changed_case = replace(
            case,
            context=replace(case.context, baseline=replace(case.context.baseline, age_years=90)),
        )
        check(
            "reviewed_values_cannot_be_silently_changed",
            lambda: verify_case_review(
                changed_case,
                path,
                expected_sha256=fingerprint,
                inventory=inventory,
                verified_documents=verify_documents(documents),
            ),
            "case_differs_from_reviewed_source",
        )
        check(
            "unverified_document_cannot_supply_lineage",
            lambda: verify_case_review(
                case,
                path,
                expected_sha256=fingerprint,
                inventory=inventory,
                verified_documents={},
            ),
            "case_lineage_source_unverified",
        )
        serialized = json.dumps(audit, allow_nan=False)
        check(
            "public_projection_contains_no_patient_or_series_ids",
            lambda: (
                all(patient not in serialized for patient in inventory.patient_ids)
                and all(series.series_uid not in serialized for series in inventory.series)
            ),
        )
        check(
            "mutable_series_inventory_is_rejected",
            lambda: replace(inventory, series=list(inventory.series)),
            "inventory_series_must_be_immutable",
        )
        reordered = replace(
            inventory,
            clinical=dict(reversed(tuple(inventory.clinical.items()))),
            expression_observed_per_patient=dict(
                reversed(tuple(inventory.expression_observed_per_patient.items()))
            ),
        )
        check(
            "source_coordinates_survive_mapping_order",
            lambda: private_manifest(inventory) == private_manifest(reordered),
        )

        def probe_review(name, mutate, expected, update_assertion=False):
            review = json.loads(path.read_bytes())
            mutate(review)
            supporting = documents
            if update_assertion:
                original_hash, original_doc = next(iter(documents.items()))
                document = json.loads(original_doc.read_bytes())
                document["cases"][case.context.baseline.patient_id] = {
                    key: review[key] for key in ("case", "pet_ct_pairings", "molecular_provenance")
                }
                support_path = root / f"{name}-assertions.json"
                support_path.write_text(json.dumps(document, sort_keys=True), encoding="utf-8")
                digest = hashlib.sha256(support_path.read_bytes()).hexdigest()
                supporting = {digest: support_path}
                for link in review["references"].values():
                    if link["source_sha256"] == original_hash:
                        link["source_sha256"] = digest
            probe_path = root / f"{name}-PRIVATE-review.json"
            probe_path.write_text(json.dumps(review, sort_keys=True), encoding="utf-8")
            fingerprint = hashlib.sha256(probe_path.read_bytes()).hexdigest()
            check(
                name,
                lambda: verify_case_review(
                    case,
                    probe_path,
                    expected_sha256=fingerprint,
                    inventory=inventory,
                    verified_documents=verify_documents(supporting),
                ),
                expected,
            )

        probe_review(
            "nonexistent_clinical_locator_is_rejected",
            lambda review: review["references"][case.context.baseline.diagnosis.source_id].update(
                locator="/not-present"
            ),
            "clinical_locator_unresolved",
        )
        probe_review(
            "PET_CT_pairing_must_be_explicit",
            lambda review: review.update(pet_ct_pairings={}),
            "pet_ct_pairing_undeclared",
        )
        probe_review(
            "paired_CT_available_after_index_is_rejected",
            lambda review: next(iter(review["pet_ct_pairings"].values()))["ct_evidence"].update(
                available_on="2025-01-02"
            ),
            "paired_ct_not_available_at_index",
        )
        probe_review(
            "wrong_expression_column_is_rejected",
            lambda review: next(iter(review["molecular_provenance"].values())).update(
                expression_locator="tsv/column/999"
            ),
            "molecular_provenance_mismatch",
            update_assertion=True,
        )
    return {
        "schema_version": "1.0",
        "evidence_type": "synthetic_software_validation",
        "implementation_sha256": implementation_fingerprint(),
        "protocol_sha256": protocol.sha256,
        "checks": checks,
        "status": "PASS" if all(item["passed"] for item in checks) else "FAIL",
        "real_clinical_eligibility": "NOT_ESTABLISHED",
        "model_evaluation": "NOT_RUN",
    }


def main() -> int:
    """Recompute the synthetic snapshot offline; no source dataset is needed in CI."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=Path("configs/research/nsclc_os_v1.json"))
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-report", type=Path)
    args = parser.parse_args()
    report = validation_report(load_research_protocol(args.protocol))
    if args.check_report:
        if json.loads(args.check_report.read_bytes()) != report:
            return 1
    if args.output:
        write_report(args.output, report)
    print(json.dumps({"status": report["status"], "checks": len(report["checks"])}))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
