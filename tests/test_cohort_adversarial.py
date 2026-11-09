"""Independent fixtures and adversarial clinical/provenance oracles for Commit 04."""

import csv
import gzip
import hashlib
import io
import json
import os
import subprocess
import sys
from dataclasses import asdict, replace
from datetime import date, timedelta
from pathlib import Path
from types import SimpleNamespace

import pytest

from oncobio_agentfm.clinical import ClinicalStage, ECOGStatus, NSCLCBaseline
from oncobio_agentfm.cohort_audit import availability_report, private_manifest
from oncobio_agentfm.cohort_linkage import audit_linked_cases, verify_case_review, verify_documents
from oncobio_agentfm.cohort_selection import compare_selection
from oncobio_agentfm.cohort_sources import load_inventory
from oncobio_agentfm.contracts import EvidenceRecord, MissingnessReason
from oncobio_agentfm.prognosis import (
    AvailabilityBasis,
    OverallSurvivalFollowUp,
    PredictionContext,
    PredictorFeature,
    ProtocolViolation,
    VitalStatus,
)
from oncobio_agentfm.protocol_audit import AdjudicatedCase
from oncobio_agentfm.research_protocol import FeatureKind
from oncobio_agentfm.source_audit import REVIEWED_SOURCE_COLUMNS

T0 = date(2025, 1, 1)


def plain(value):
    return json.loads(json.dumps(value, default=lambda item: item.isoformat()))


@pytest.fixture
def scenario(tmp_path):
    """Three invented patients; separate CT/PT study IDs with explicitly reviewed pairing."""

    patients = ("SYN-B", "SYN-C", "SYN-A")
    clinical = io.StringIO(newline="")
    writer = csv.DictWriter(clinical, fieldnames=REVIEWED_SOURCE_COLUMNS)
    writer.writeheader()
    for patient in patients:
        writer.writerow(
            {
                "Case ID": patient,
                "Age at Histological Diagnosis": "60",
                "Gender": "Female",
                "Smoking status": "Former",
                "Patient affiliation": "Stanford",
                "Survival Status": "Dead" if patient == "SYN-A" else "Alive",
            }
        )
    series = []
    samples = "^SERIES = GSE103584\n!Series_geo_accession = GSE103584\n"
    cases = []
    reviews = {}
    for number, patient in enumerate(patients, 1):

        def evidence(ref, acquired, available, patient=patient):
            return EvidenceRecord(
                patient, ref, T0 + timedelta(days=acquired), T0 + timedelta(days=available)
            )

        absent = EvidenceRecord(patient, None, None, None, MissingnessReason.NOT_RECORDED)
        baseline = NSCLCBaseline(
            patient,
            T0,
            evidence("diagnosis:" + patient, -20, -15),
            60,
            evidence("age:" + patient, 0, 0),
            ClinicalStage(absent, None, None),
            ECOGStatus(absent, None),
        )
        context = PredictionContext(
            baseline,
            evidence("assessment:" + patient, 0, 0),
            evidence("referral:" + patient, -2, -1),
            evidence("untreated:" + patient, 0, 0),
            "aligned:" + patient,
            "clinical-review:" + patient,
            T0,
            "earliest:" + patient,
            T0 + timedelta(days=10),
            evidence("resection:" + patient, 10, 11),
        )
        ct_uid, pt_uid, paired_uid = (f"1.2.{prefix}.{number}" for prefix in (10, 20, 30))
        for uid, modality, offset, study in (
            (ct_uid, "CT", -10, f"1.2.100.{number}"),
            (pt_uid, "PT", -2, f"1.2.200.{number}"),
            (paired_uid, "CT", -2, f"1.2.300.{number}"),
        ):
            series.append(
                {
                    "PatientID": patient,
                    "SeriesInstanceUID": uid,
                    "StudyInstanceUID": study,
                    "Modality": modality,
                    "ImageCount": 20,
                    "Collection": "NSCLC Radiogenomics",
                    "SeriesDate": (T0 + timedelta(days=offset)).isoformat() + " 00:00:00.0",
                    "StudyDate": (T0 + timedelta(days=offset)).isoformat() + " 00:00:00.0",
                }
            )
        ct = PredictorFeature(
            FeatureKind.CT,
            evidence(ct_uid, -10, -8),
            "ct-access:" + patient,
            AvailabilityBasis.CLINICAL_RECORD,
        )
        pet = PredictorFeature(FeatureKind.PET_CT, absent, None)
        if patient != "SYN-B":
            pet = PredictorFeature(
                FeatureKind.PET_CT,
                evidence(pt_uid, -2, -1),
                "pet-access:" + patient,
                AvailabilityBasis.CLINICAL_RECORD,
            )
        accession = f"GSM{800000 + number}"
        samples += (
            f"^SAMPLE = {accession}\n!Sample_title = {patient}\n!Sample_type = SRA\n"
            "!Sample_library_strategy = RNA-Seq\n!Sample_series_id = GSE103584\n"
        )
        molecular = PredictorFeature(FeatureKind.RESECTION_MOLECULAR, absent, None)
        if patient != "SYN-C":
            molecular = PredictorFeature(
                FeatureKind.RESECTION_MOLECULAR,
                evidence(accession, 10, 25),
                "molecular-release:" + patient,
                AvailabilityBasis.CLINICAL_RECORD,
            )
        followup = OverallSurvivalFollowUp(
            patient,
            VitalStatus.DEAD if patient == "SYN-A" else VitalStatus.ALIVE,
            T0 + timedelta(days=100) if patient == "SYN-A" else None,
            T0 + timedelta(days=80 if patient == "SYN-A" else 400),
            T0 + timedelta(days=200),
            "endpoint:" + patient,
            "aligned:" + patient,
            "shared-reviewed-cutoff",
        )
        case = AdjudicatedCase(context, (ct, pet), (ct, pet, molecular), followup)
        cases.append(case)
        pairs = {}
        if pet.evidence.is_present:
            pairs[pt_uid] = {
                "ct_evidence": plain(asdict(evidence(paired_uid, -2, -1))),
                "timing_source": "pair-access:" + patient,
                "pairing_source": "pair-review:" + patient,
                "availability_basis": "clinical_record",
            }
        reviews[patient] = {
            "schema_version": "1.1",
            "case": plain(asdict(case)),
            "pet_ct_pairings": pairs,
            "molecular_provenance": {},
            "references": {},
            "image_date_mappings": {},
        }
    expression_patients = ("SYN-C", "SYN-A", "SYN-B")
    raw = {
        "clinical": clinical.getvalue().encode(),
        "imaging": json.dumps(series).encode(),
        "rna_metadata": gzip.compress(samples.encode(), mtime=0),
        "rna_expression": gzip.compress(
            (
                "\t" + "\t".join(expression_patients) + "\nGENE-A\t1\t0\tNA\nGENE-B\t2\t3\t4\n"
            ).encode(),
            mtime=0,
        ),
    }
    paths, receipts = {}, {}
    for kind, value in raw.items():
        paths[kind] = tmp_path / (kind + ".bin")
        paths[kind].write_bytes(value)
        receipts[kind] = {
            "version": "independent-invented-v1",
            "retrieved_on": "2025-01-01",
            "size_bytes": len(value),
            "sha256": hashlib.sha256(value).hexdigest(),
        }
    lock = tmp_path / "lock.json"
    lock.write_text(json.dumps({"schema_version": "1.0", "sources": receipts}))
    inventory = load_inventory(paths, lock)
    for case in cases:
        patient = case.context.baseline.patient_id
        review = reviews[patient]
        clinical_refs = {
            "diagnosis:" + patient,
            "age:" + patient,
            "assessment:" + patient,
            "referral:" + patient,
            "untreated:" + patient,
            "aligned:" + patient,
            "clinical-review:" + patient,
            "earliest:" + patient,
            "resection:" + patient,
            "endpoint:" + patient,
            "shared-reviewed-cutoff",
            "ct-access:" + patient,
        }
        image_items = [
            (feature.evidence, "CT" if feature.kind is FeatureKind.CT else "PT")
            for feature in case.primary_features
            if feature.evidence.is_present
        ]
        if case.primary_features[1].evidence.is_present:
            clinical_refs.update(
                ("pet-access:" + patient, "pair-access:" + patient, "pair-review:" + patient)
            )
            pair = next(iter(review["pet_ct_pairings"].values()))
            paired = pair["ct_evidence"]
            image_items.append(
                (
                    EvidenceRecord(
                        patient,
                        paired["source_id"],
                        date.fromisoformat(paired["acquired_on"]),
                        date.fromisoformat(paired["available_on"]),
                    ),
                    "CT",
                )
            )
        for item, _ in image_items:
            index = next(
                i for i, row in enumerate(series) if row["SeriesInstanceUID"] == item.source_id
            )
            review["references"][item.source_id] = {
                "source_sha256": receipts["imaging"]["sha256"],
                "locator": f"/{index}",
            }
            review["image_date_mappings"][item.source_id] = {
                "source_axis_date": series[index]["SeriesDate"][:10],
                "source_date_kind": "SeriesDate_0008_0021",
                "acquired_on": item.acquired_on.isoformat(),
                "acquisition_reference": (
                    "ct-access:" + patient
                    if item.source_id.startswith("1.2.10.")
                    else (
                        "pet-access:" + patient
                        if item.source_id.startswith("1.2.20.")
                        else "pair-access:" + patient
                    )
                ),
                "alignment_reference": "aligned:" + patient,
            }
        feature = case.landmark_features[-1]
        if feature.evidence.is_present:
            item = feature.evidence
            clinical_refs.add("molecular-release:" + patient)
            review["references"][item.source_id] = {
                "source_sha256": receipts["rna_metadata"]["sha256"],
                "locator": "sample/" + item.source_id,
            }
            review["molecular_provenance"][item.source_id] = {
                "expression_source_sha256": receipts["rna_expression"]["sha256"],
                "expression_locator": f"tsv/column/{expression_patients.index(patient) + 1}",
                "resection_source_id": "resection:" + patient,
                "specimen_collected_on": (T0 + timedelta(days=10)).isoformat(),
                "clinical_available_on": (T0 + timedelta(days=25)).isoformat(),
                "availability_basis": "clinical_record",
                "timing_source": "molecular-release:" + patient,
            }
        for ref in clinical_refs:
            review["references"][ref] = {"source_sha256": "PENDING", "locator": "/cases/" + patient}
    state = {
        "inventory": inventory,
        "cases": tuple(cases),
        "payloads": reviews,
        "documents": {},
        "clinical_order": patients,
        "expression_order": expression_patients,
    }

    def resign():
        bundles = {
            patient: {
                key: review[key] for key in ("case", "pet_ct_pairings", "molecular_provenance")
            }
            for patient, review in reviews.items()
        }
        document = tmp_path / "PRIVATE-assertions.json"
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
        fingerprint = hashlib.sha256(document.read_bytes()).hexdigest()
        old = set(state["documents"]) | {"PENDING"}
        state["documents"] = {fingerprint: document}
        state["reviews"] = {}
        for patient, review in reviews.items():
            for ref in review["references"].values():
                if ref["source_sha256"] in old:
                    ref["source_sha256"] = fingerprint
            path = tmp_path / (patient + "-PRIVATE-review.json")
            path.write_text(json.dumps(review, sort_keys=True), encoding="utf-8")
            state["reviews"][patient] = (path, hashlib.sha256(path.read_bytes()).hexdigest())

    state["resign"] = resign
    resign()
    return state


def run(state, protocol):
    return audit_linked_cases(
        state["cases"],
        inventory=state["inventory"],
        reviews=state["reviews"],
        documents=state["documents"],
        protocol=protocol,
    )


def test_independent_event_and_landmark_oracle(scenario, protocol):
    report = run(scenario, protocol)
    assert report["validated_event_budget"]["pretreatment"] == {
        "eligible_patients": 3,
        "deaths": 1,
        "censored": 2,
    }
    assert report["validated_event_budget"]["postoperative_day_30"] == {
        "eligible_patients": 3,
        "deaths": 1,
        "censored": 2,
    }
    assert report["all_four_landmark_patients"] == report["all_four_landmark_deaths"] == 1


def test_reviewed_cross_study_pair_does_not_impose_an_unprespecified_selection_filter(
    scenario, protocol
):
    assert (
        availability_report(scenario["inventory"])["imaging_reconciliation"][
            "patients_with_PT_and_CT_in_same_study"
        ]
        == 0
    )
    assert run(scenario, protocol)["settings"]["pretreatment"]["endpoint_eligible"] == 3


def test_source_coordinates_survive_mapping_reordering(scenario):
    inventory = scenario["inventory"]
    reordered = replace(
        inventory,
        clinical=dict(reversed(tuple(inventory.clinical.items()))),
        expression_observed_per_patient=dict(
            reversed(tuple(inventory.expression_observed_per_patient.items()))
        ),
    )
    assert private_manifest(inventory) == private_manifest(reordered)
    for record in private_manifest(reordered)["records"]:
        patient = record["patient_id"]
        refs = {link["role"]: link["ref"]["locator"] for link in record["source_links"]}
        assert refs["clinical"] == f"csv/data-row/{scenario['clinical_order'].index(patient)}"
        assert (
            refs["rna_expression"]
            == f"tsv/column/{scenario['expression_order'].index(patient) + 1}"
        )


def test_mutable_series_and_invalid_coordinates_are_rejected(scenario):
    inventory = scenario["inventory"]
    with pytest.raises(ProtocolViolation, match="inventory_series_must_be_immutable"):
        replace(inventory, series=list(inventory.series))
    with pytest.raises(ProtocolViolation, match="inventory_source_coordinates_invalid"):
        replace(inventory, clinical_source_rows={patient: 0 for patient in inventory.clinical})


@pytest.mark.parametrize(
    "pointer,code",
    [
        ("/nonexistent", "clinical_locator_unresolved"),
        ("/cases/SYN-B", "clinical_assertion_case_mismatch"),
    ],
)
def test_clinical_citations_resolve_to_the_correct_patient(scenario, protocol, pointer, code):
    scenario["payloads"]["SYN-A"]["references"]["diagnosis:SYN-A"]["locator"] = pointer
    scenario["resign"]()
    with pytest.raises(ProtocolViolation, match=code):
        run(scenario, protocol)


def test_imaging_metadata_cannot_prove_diagnosis_or_clinical_access(scenario, protocol):
    review = scenario["payloads"]["SYN-A"]
    review["references"]["diagnosis:SYN-A"] = dict(review["references"]["1.2.10.3"])
    scenario["resign"]()
    with pytest.raises(ProtocolViolation, match="clinical_assertion_source_required"):
        run(scenario, protocol)


def test_pair_component_availability_is_gated_at_the_original_index(scenario, protocol):
    pair = next(iter(scenario["payloads"]["SYN-A"]["pet_ct_pairings"].values()))
    pair["ct_evidence"]["available_on"] = (T0 + timedelta(days=1)).isoformat()
    scenario["resign"]()
    with pytest.raises(ProtocolViolation, match="paired_ct_not_available_at_index"):
        run(scenario, protocol)


def test_pair_relationship_cannot_be_left_implicit(scenario, protocol):
    scenario["payloads"]["SYN-A"]["pet_ct_pairings"] = {}
    scenario["resign"]()
    with pytest.raises(ProtocolViolation, match="pet_ct_pairing_undeclared"):
        run(scenario, protocol)


@pytest.mark.parametrize(
    "field,value",
    [
        ("expression_locator", "tsv/column/1"),
        ("expression_source_sha256", "0" * 64),
        ("resection_source_id", "other-operation"),
        ("clinical_available_on", "2010-01-01"),
    ],
)
def test_molecular_source_column_and_specimen_release_binding(scenario, protocol, field, value):
    proof = next(iter(scenario["payloads"]["SYN-A"]["molecular_provenance"].values()))
    proof[field] = value
    scenario["resign"]()
    with pytest.raises(ProtocolViolation, match="molecular_provenance_mismatch"):
        run(scenario, protocol)


def test_ambiguous_aliquots_cannot_be_assigned_to_one_patient_column(scenario, protocol):
    inventory = scenario["inventory"]
    scenario["inventory"] = replace(
        inventory, rna_samples={**inventory.rna_samples, "GSM999999": "SYN-A"}
    )
    with pytest.raises(ProtocolViolation, match="rna_aliquot_column_link_ambiguous"):
        run(scenario, protocol)


def test_bare_or_forged_document_fingerprints_are_not_proof(scenario):
    case = scenario["cases"][0]
    path, fingerprint = scenario["reviews"][case.context.baseline.patient_id]
    for supplied, code in (
        (frozenset(scenario["documents"]), "case_document_bytes_required"),
        ({next(iter(scenario["documents"])): b"altered"}, "adjudication_document_hash_mismatch"),
    ):
        with pytest.raises(ProtocolViolation, match=code):
            verify_case_review(
                case,
                path,
                expected_sha256=fingerprint,
                inventory=scenario["inventory"],
                verified_documents=supplied,
            )
    verified = verify_documents(scenario["documents"])
    with pytest.raises(TypeError):
        verified["new"] = b"unverified"


@pytest.mark.parametrize(
    "death_offset,events,post_eligible", [(40, 1, 2), (200, 1, 3), (201, 0, 3)]
)
def test_death_landmark_and_cutoff_boundaries(
    scenario, protocol, death_offset, events, post_eligible
):
    cases = list(scenario["cases"])
    index = next(i for i, case in enumerate(cases) if case.context.baseline.patient_id == "SYN-A")
    old = cases[index]
    cases[index] = replace(
        old,
        follow_up=replace(
            old.follow_up,
            death_on=T0 + timedelta(days=death_offset),
            last_known_alive_on=T0 + timedelta(days=20),
        ),
    )
    scenario["cases"] = tuple(cases)
    scenario["payloads"]["SYN-A"]["case"] = plain(asdict(cases[index]))
    scenario["resign"]()
    report = run(scenario, protocol)
    assert report["validated_event_budget"]["pretreatment"]["deaths"] == events
    assert (
        report["validated_event_budget"]["postoperative_day_30"]["eligible_patients"]
        == post_eligible
    )


@pytest.mark.parametrize("basis", [AvailabilityBasis.RESEARCH_ONLY, AvailabilityBasis.UNKNOWN])
def test_research_assay_existence_never_substitutes_for_clinical_availability(
    scenario, protocol, basis
):
    cases = list(scenario["cases"])
    index = next(i for i, case in enumerate(cases) if case.context.baseline.patient_id == "SYN-A")
    old = cases[index]
    feature = replace(old.landmark_features[-1], availability_basis=basis)
    cases[index] = replace(old, landmark_features=(*old.landmark_features[:-1], feature))
    scenario["cases"] = tuple(cases)
    scenario["payloads"]["SYN-A"]["case"] = plain(asdict(cases[index]))
    next(iter(scenario["payloads"]["SYN-A"]["molecular_provenance"].values()))[
        "availability_basis"
    ] = basis.value
    scenario["resign"]()
    report = run(scenario, protocol)
    assert report["settings"]["pretreatment"]["endpoint_eligible"] == 3
    assert report["settings"]["postoperative_day_30"]["predictor_rejections"] == {
        "feature_clinical_availability_unverified": 1
    }


def test_wrong_patient_image_cannot_pass_even_with_a_reviewed_case(scenario, protocol):
    cases = list(scenario["cases"])
    index = next(i for i, case in enumerate(cases) if case.context.baseline.patient_id == "SYN-A")
    old = cases[index]
    wrong = replace(
        old.primary_features[0],
        evidence=replace(old.primary_features[0].evidence, source_id="1.2.10.1"),
    )
    cases[index] = replace(
        old,
        primary_features=(wrong, old.primary_features[1]),
        landmark_features=(wrong, *old.landmark_features[1:]),
    )
    scenario["cases"] = tuple(cases)
    review = scenario["payloads"]["SYN-A"]
    review["case"] = plain(asdict(cases[index]))
    del review["references"]["1.2.10.3"]
    del review["image_date_mappings"]["1.2.10.3"]
    review["references"]["1.2.10.1"] = dict(scenario["payloads"]["SYN-B"]["references"]["1.2.10.1"])
    review["image_date_mappings"]["1.2.10.1"] = {
        **scenario["payloads"]["SYN-B"]["image_date_mappings"]["1.2.10.1"],
        "alignment_reference": "aligned:SYN-A",
    }
    scenario["resign"]()
    with pytest.raises(ProtocolViolation, match="selected_image_patient_or_source_mismatch"):
        run(scenario, protocol)


@pytest.mark.parametrize(
    "field,offset,code",
    [
        ("available_on", 41, "post_index_evidence"),
        ("acquired_on", 11, "molecular_specimen_not_linked_to_resection"),
    ],
)
def test_molecular_timing_keeps_original_landmark_and_specimen_contracts(
    scenario, protocol, field, offset, code
):
    cases = list(scenario["cases"])
    index = next(i for i, case in enumerate(cases) if case.context.baseline.patient_id == "SYN-A")
    old = cases[index]
    molecular = old.landmark_features[-1]
    molecular = replace(
        molecular, evidence=replace(molecular.evidence, **{field: T0 + timedelta(days=offset)})
    )
    cases[index] = replace(old, landmark_features=(*old.landmark_features[:-1], molecular))
    scenario["cases"] = tuple(cases)
    review = scenario["payloads"]["SYN-A"]
    review["case"] = plain(asdict(cases[index]))
    proof = next(iter(review["molecular_provenance"].values()))
    proof["clinical_available_on" if field == "available_on" else "specimen_collected_on"] = (
        T0 + timedelta(days=offset)
    ).isoformat()
    scenario["resign"]()
    report = run(scenario, protocol)
    assert report["settings"]["pretreatment"]["endpoint_eligible"] == 3
    assert report["settings"]["postoperative_day_30"]["predictor_rejections"] == {code: 1}


def test_source_cohort_composition_is_not_a_within_cohort_missingness_effect():
    # Independent arithmetic oracle; no source/clinical eligibility claim is made here.
    rows = {}
    for group, genders in (
        ("R01", ("Female", "Male", "Female", "Male")),
        ("AMC", ("Female", "Female", "Female", "Female")),
    ):
        for index, gender in enumerate(genders):
            rows[f"{group}-SYN-{index}"] = {
                "Age at Histological Diagnosis": "60",
                "Gender": gender,
                "Patient affiliation": "Stanford",
                "Smoking status": "Former",
            }
    result = compare_selection(
        SimpleNamespace(clinical=rows), frozenset(("R01-SYN-0", "R01-SYN-1"))
    )
    assert result["categorical"]["Gender"]["levels"][0]["proportion_difference"] == pytest.approx(
        -1 / 3
    )
    assert (
        result["source_cohort_strata"]["R01"]["categorical"]["Gender"]["levels"][0][
            "proportion_difference"
        ]
        == 0
    )
    assert (
        result["source_cohort_strata"]["AMC"]["comparison_support"]
        == "NO_WITHIN_STRATUM_COMPARISON"
    )


def test_fingerprint_covers_inherited_policy_and_normalizes_windows_newlines(tmp_path, monkeypatch):
    from oncobio_agentfm import cohort_sources

    package = Path(cohort_sources.__file__).parent
    baseline = cohort_sources.implementation_fingerprint()
    for path in package.glob("*.py"):
        normalized = path.read_bytes().replace(b"\r\n", b"\n")
        (tmp_path / path.name).write_bytes(normalized.replace(b"\n", b"\r\n"))
    monkeypatch.setattr(cohort_sources, "__file__", str(tmp_path / "cohort_sources.py"))
    assert cohort_sources.implementation_fingerprint() == baseline
    (tmp_path / "splitting.py").write_bytes(
        (tmp_path / "splitting.py").read_bytes() + b"\r\n# changed policy source\r\n"
    )
    assert cohort_sources.implementation_fingerprint() != baseline


def test_invalid_source_age_does_not_leak_raw_text_through_CLI_traceback(scenario, tmp_path):
    path = tmp_path / "clinical.bin"
    path.write_bytes(path.read_bytes().replace(b",60,", b",PRIVATE-INVALID-AGE-TEXT,"))
    lock_path = tmp_path / "lock.json"
    lock = json.loads(lock_path.read_bytes())
    lock["sources"]["clinical"].update(
        sha256=hashlib.sha256(path.read_bytes()).hexdigest(), size_bytes=path.stat().st_size
    )
    lock_path.write_text(json.dumps(lock))
    arguments = [
        sys.executable,
        "-m",
        "oncobio_agentfm.cohort_audit",
        "--lock",
        str(lock_path),
        "--output-dir",
        "outputs/invalid-age",
    ]
    for kind in ("clinical", "imaging", "rna_metadata", "rna_expression"):
        arguments.extend(("--" + kind.replace("_", "-"), str(tmp_path / (kind + ".bin"))))
    environment = {**os.environ, "PYTHONPATH": str(Path(__file__).resolve().parents[1] / "src")}
    process = subprocess.run(  # noqa: S603 -- fixed interpreter and generated local fixture paths
        arguments, cwd=tmp_path, env=environment, capture_output=True, text=True, check=False
    )
    assert process.returncode == 1
    assert "selection_age_invalid" in process.stderr
    assert "PRIVATE-INVALID-AGE-TEXT" not in process.stderr
    assert "Traceback" not in process.stderr


def test_series_date_is_metadata_and_never_an_inferred_acquisition_date(scenario):
    inventory = scenario["inventory"]
    report = availability_report(inventory)
    assert "series_missing_acquisition_date" not in report["imaging_reconciliation"]
    assert report["imaging_reconciliation"]["series_missing_series_date"] == 0
    for record in private_manifest(inventory)["records"]:
        for link in record["source_links"]:
            if link["role"] in {"CT", "PT", "SEG"}:
                assert link["series_date_on_source_axis"] is not None
                assert link["acquired_on_source_axis"] is None


def test_separately_adjudicated_acquisition_can_exist_without_SeriesDate(scenario, protocol):
    inventory = scenario["inventory"]
    scenario["inventory"] = replace(
        inventory,
        series=tuple(
            replace(series, series_date=None) if series.series_uid == "1.2.10.3" else series
            for series in inventory.series
        ),
    )
    scenario["payloads"]["SYN-A"]["image_date_mappings"]["1.2.10.3"]["source_axis_date"] = None
    scenario["resign"]()
    assert run(scenario, protocol)["settings"]["pretreatment"]["endpoint_eligible"] == 3


def test_metadata_date_tag_cannot_be_relabelled_as_acquisition(scenario, protocol):
    scenario["payloads"]["SYN-A"]["image_date_mappings"]["1.2.10.3"]["source_date_kind"] = (
        "AcquisitionDate_0008_0022"
    )
    scenario["resign"]()
    with pytest.raises(ProtocolViolation, match="selected_image_calendar_mapping_mismatch"):
        run(scenario, protocol)


def test_source_evidence_is_bound_to_current_implementation_and_independent_oracle():
    from oncobio_agentfm.cohort_sources import implementation_fingerprint

    root = Path(__file__).resolve().parents[1]
    source = json.loads((root / "docs/evidence/source_availability.json").read_bytes())
    reproduction = json.loads((root / "docs/evidence/source_reproduction.json").read_bytes())
    assert source["implementation_sha256"] == implementation_fingerprint()
    assert reproduction["audited_implementation_sha256"] == source["implementation_sha256"]
    assert (
        reproduction["oracle_sha256"]
        == hashlib.sha256(
            (root / "tools/reproduce_source_audit.py").read_bytes().replace(b"\r\n", b"\n")
        ).hexdigest()
    )
    assert reproduction["status"] == "PASS"
