"""Adversarial cohort tests: identity, lineage, missingness, attrition and timing."""

import gzip
import hashlib
import json
from dataclasses import replace

import pytest

from oncobio_agentfm.cohort_audit import (
    availability_report,
    flow_markdown,
    private_manifest,
    source_sets,
)
from oncobio_agentfm.cohort_linkage import (
    SourceExclusion,
    audit_linked_cases,
    verify_case_review,
    verify_documents,
)
from oncobio_agentfm.cohort_selection import compare_selection
from oncobio_agentfm.cohort_sources import (
    load_inventory,
    read_clinical,
    read_expression,
    read_rna_metadata,
    read_series,
)
from oncobio_agentfm.cohort_validation import (
    synthetic_bundle,
    synthetic_linked_cases,
    validation_report,
    write_synthetic_reviews,
)
from oncobio_agentfm.prognosis import ProtocolViolation


@pytest.fixture
def inputs(tmp_path):
    return synthetic_bundle(tmp_path)


@pytest.fixture
def inventory(inputs):
    return load_inventory(*inputs)


@pytest.fixture
def linked(tmp_path, inventory):
    return synthetic_linked_cases(tmp_path, inventory)


def replace_synthetic_source(inputs, kind, raw):
    """Deliberately review new synthetic bytes; never used with real source locks."""

    paths, lock_path = inputs
    paths[kind].write_bytes(raw)
    lock = json.loads(lock_path.read_bytes())
    lock["sources"][kind]["sha256"] = hashlib.sha256(raw).hexdigest()
    lock["sources"][kind]["size_bytes"] = len(raw)
    lock_path.write_text(json.dumps(lock), encoding="utf-8")


def rewrite_review(path, review):
    path.write_text(json.dumps(review), encoding="utf-8")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_source_union_and_all_four_overlap_have_known_independent_denominators(inventory):
    report = availability_report(inventory)
    assert report["patients"] == 12
    assert report["modality_counts"] == {"clinical": 12, "ct": 12, "pet": 6, "rna": 4}
    assert report["source_population_label_audit"]["all_four"]["patients"] == 2
    assert sum(row["patients"] for row in report["availability_patterns"]) == 12
    assert report["analysis_populations"]["pretreatment"]["eligible_patients"] is None
    assert report["analysis_populations"]["postoperative_day_30"]["survival_events"] is None


def test_input_bytes_cannot_change_under_an_existing_source_lock(inputs):
    paths, lock = inputs
    raw = paths["clinical"].read_bytes()
    paths["clinical"].write_bytes(raw.replace(b"Female", b"FEMALE"))
    with pytest.raises(ProtocolViolation, match="source_hash_mismatch"):
        load_inventory(paths, lock)


def test_truncated_source_is_not_silently_accepted(inputs):
    paths, lock = inputs
    paths["clinical"].write_bytes(paths["clinical"].read_bytes()[:-1])
    with pytest.raises(ProtocolViolation, match="source_size_mismatch"):
        load_inventory(paths, lock)


def test_missing_source_is_not_relabelled_as_missing_patients(inputs):
    paths, lock = inputs
    with pytest.raises(ProtocolViolation, match="source_set_incomplete"):
        load_inventory({key: value for key, value in paths.items() if key != "imaging"}, lock)


def test_boolean_cannot_replace_a_source_size(inputs):
    paths, lock_path = inputs
    lock = json.loads(lock_path.read_bytes())
    lock["sources"]["clinical"]["size_bytes"] = True
    lock_path.write_text(json.dumps(lock), encoding="utf-8")
    with pytest.raises(ProtocolViolation, match="source_size_invalid"):
        load_inventory(paths, lock_path)


def test_inventory_cannot_mutate_after_source_validation(inventory):
    with pytest.raises(TypeError):
        inventory.clinical["SYN-cohort-0"]["Gender"] = "changed"
    with pytest.raises(TypeError):
        inventory.rna_samples["GSM900000"] = "different-owner"


def test_duplicate_clinical_patient_cannot_inflate_sample_size(inputs):
    raw = inputs[0]["clinical"].read_bytes()
    duplicated = raw + raw.splitlines(keepends=True)[1]
    with pytest.raises(ProtocolViolation, match="clinical_duplicate_patient"):
        read_clinical(duplicated)


@pytest.mark.parametrize("change", [b" SYN-cohort-0", b"SYN-cohort-0 ", b"SYN cohort 0"])
def test_identity_is_exact_and_never_repaired(inputs, change):
    raw = inputs[0]["clinical"].read_bytes().replace(b"SYN-cohort-0", change)
    with pytest.raises(ProtocolViolation, match="source_patient_identifier_invalid"):
        read_clinical(raw)


def test_orphan_imaging_identity_remains_in_union_and_outside_clinical_population(inputs):
    paths, _ = inputs
    rows = json.loads(paths["imaging"].read_bytes())
    rows.append(
        {
            **rows[0],
            "PatientID": "SYN-orphan",
            "SeriesInstanceUID": "1.2.999.1",
            "StudyInstanceUID": "1.2.999.2",
        }
    )
    replace_synthetic_source(inputs, "imaging", json.dumps(rows).encode())
    report = availability_report(load_inventory(*inputs))
    assert report["patients"] == 13
    assert report["linkage"]["imaging_without_clinical"] == 1
    assert report["source_population_label_audit"]["ct_clinical"]["patients"] == 12


def test_exact_duplicate_series_is_collapsed_once(inputs):
    rows = json.loads(inputs[0]["imaging"].read_bytes())
    series, duplicates = read_series(json.dumps([*rows, rows[0]]).encode())
    assert len(series) == len(rows)
    assert duplicates == 1


@pytest.mark.parametrize(
    "field,value,code",
    [
        ("ImageCount", "99", "imaging_series_metadata_conflict"),
        ("Modality", "PT", "imaging_series_metadata_conflict"),
        ("PatientID", "SYN-wrong-owner", "imaging_study_owner_conflict"),
    ],
)
def test_same_series_cannot_conflict_with_a_second_source_row(inputs, field, value, code):
    rows = json.loads(inputs[0]["imaging"].read_bytes())
    with pytest.raises(ProtocolViolation, match=code):
        read_series(json.dumps([*rows, {**rows[0], field: value}]).encode())


def test_duplicate_series_cannot_switch_study_even_with_same_patient(inputs):
    rows = json.loads(inputs[0]["imaging"].read_bytes())
    with pytest.raises(ProtocolViolation, match="imaging_series_metadata_conflict"):
        read_series(json.dumps([*rows, {**rows[0], "StudyInstanceUID": "1.2.999"}]).encode())


def test_multiple_distinct_series_are_candidates_not_duplicates_or_selected_images(inputs):
    rows = json.loads(inputs[0]["imaging"].read_bytes())
    rows.append({**rows[0], "SeriesInstanceUID": "1.2.999.3"})
    replace_synthetic_source(inputs, "imaging", json.dumps(rows).encode())
    report = availability_report(load_inventory(*inputs))
    assert report["imaging_reconciliation"]["patients_with_multiple_candidate_series"]["CT"] == 1
    assert report["modality_counts"]["ct"] == 12
    assert report["imaging_reconciliation"]["diagnostic_CT_selection"].startswith("NOT_ADJUDICATED")


@pytest.mark.parametrize("count", [True, 0, -1, "1.5", "", None])
def test_invalid_image_counts_do_not_create_available_modalities(inputs, count):
    rows = json.loads(inputs[0]["imaging"].read_bytes())
    rows[0]["ImageCount"] = count
    with pytest.raises(ProtocolViolation, match=r"imaging_(count_invalid|empty_series)"):
        read_series(json.dumps(rows).encode())


def test_integer_image_counts_from_real_nbia_schema_are_accepted(inputs):
    rows = json.loads(inputs[0]["imaging"].read_bytes())
    rows[0]["ImageCount"] = 32
    assert len(read_series(json.dumps(rows).encode())[0]) == len(rows)


def test_missing_series_date_is_not_filled_from_study_date(inputs):
    rows = json.loads(inputs[0]["imaging"].read_bytes())
    del rows[0]["SeriesDate"]
    series, _ = read_series(json.dumps(rows).encode())
    item = next(s for s in series if s.series_uid == rows[0]["SeriesInstanceUID"])
    assert item.series_date is None
    assert item.study_date == "2024-12-22"


def test_absence_from_repository_is_unknown_clinical_missingness(inventory):
    report = availability_report(inventory)
    assert report["missingness"]["pet"]["absent_from_observed_inventory"] == 6
    assert report["missingness"]["pet"]["clinical_reason"] == "unknown"
    assert report["rna_matrix_integrity"]["missing_expression_cells"] == 4
    assert report["missingness"]["rna"]["absent_from_observed_inventory"] == 8


def test_recorded_field_absence_does_not_override_observed_PET_inventory(inventory):
    report = availability_report(inventory)
    assert report["source_field_tokens"]["PET Date"]["blank"] == 12
    assert report["modality_counts"]["pet"] == 6
    assert (
        report["date_field_availability_discordance"][
            "clinical_PET_date_missing_but_PT_in_inventory"
        ]
        == 6
    )


@pytest.mark.parametrize("value", ["inf", "-inf", "nan", "NaN", "bad", ""])
def test_invalid_expression_cells_cannot_become_finite_evidence(value):
    raw = gzip.compress(f"\tSYN-A\nGENE\t{value}\n".encode(), mtime=0)
    with pytest.raises(ProtocolViolation, match=r"rna_expression_(nonfinite|nonnumeric)"):
        read_expression(raw)


def test_RNA_NA_and_true_zero_have_distinct_meanings():
    raw = gzip.compress(b"\tSYN-A\tSYN-B\nGENE\tNA\t0\n", mtime=0)
    observed, genes, missing = read_expression(raw)
    assert observed == {"SYN-A": 0, "SYN-B": 1}
    assert genes == missing == 1


def test_empty_expression_column_is_not_usable_source_presence(inputs):
    raw = gzip.compress(
        b"\tSYN-cohort-0\tSYN-cohort-3\tSYN-cohort-6\tSYN-cohort-9\nGENE\tNA\t1\t1\t1\n", mtime=0
    )
    replace_synthetic_source(inputs, "rna_expression", raw)
    report = availability_report(load_inventory(*inputs))
    assert report["rna_matrix_integrity"]["columns_with_no_observed_numeric_cells"] == 1
    assert report["modality_counts"]["rna"] == 3


def test_RNA_matrix_and_GEO_patient_sets_must_agree(inputs):
    raw = gzip.compress(b"\tSYN-other\nGENE\t1\n", mtime=0)
    replace_synthetic_source(inputs, "rna_expression", raw)
    with pytest.raises(ProtocolViolation, match="rna_matrix_metadata_patient_set_mismatch"):
        load_inventory(*inputs)


def test_repeated_GEO_aliquots_retain_one_patient_identity(inputs):
    raw = gzip.decompress(inputs[0]["rna_metadata"].read_bytes())
    raw += (
        b"^SAMPLE = GSM999999\n!Sample_title = SYN-cohort-0\n!Sample_type = SRA\n"
        b"!Sample_library_strategy = RNA-Seq\n"
    )
    samples = read_rna_metadata(gzip.compress(raw, mtime=0))
    assert len(samples) == 5
    assert len(set(samples.values())) == 4


def test_gene_duplicates_and_broken_width_are_blocked():
    for raw, code in [
        (b"\tSYN-A\nGENE\t1\nGENE\t2\n", "rna_gene_identifier_invalid_or_duplicate"),
        (b"\tSYN-A\nGENE\t1\t2\n", "rna_matrix_row_width_invalid"),
    ]:
        with pytest.raises(ProtocolViolation, match=code):
            read_expression(gzip.compress(raw, mtime=0))


def test_selection_description_has_explicit_observed_and_missing_denominators(inventory):
    complete = frozenset.intersection(*source_sets(inventory).values())
    selection = compare_selection(inventory, complete)
    assert selection["complete_n"] == 2
    assert selection["incomplete_n"] == 10
    assert selection["age_at_histological_diagnosis"]["complete"]["observed_n"] == 2
    assert selection["age_at_histological_diagnosis"]["complete"]["missing_n"] == 0
    assert selection["direction"] == "complete_minus_incomplete"
    assert selection["categorical"]["Gender"]["levels"][0]["proportion_difference"] == -0.6


def test_empty_comparison_group_withholds_standardized_difference(inventory):
    report = compare_selection(inventory, frozenset())
    assert report["age_at_histological_diagnosis"]["smd"] is None
    assert all(level["smd"] is None for level in report["categorical"]["Gender"]["levels"])
    json.dumps(report, allow_nan=False)


def test_public_projection_withholds_unrecognized_source_text(inputs):
    raw = inputs[0]["clinical"].read_bytes().replace(b"Female", b"SYN-PRIVATE-FREE-TEXT")
    replace_synthetic_source(inputs, "clinical", raw)
    serialized = json.dumps(availability_report(load_inventory(*inputs)))
    assert "SYN-PRIVATE-FREE-TEXT" not in serialized


def test_public_audit_is_order_independent_and_patient_free(inventory):
    report = availability_report(inventory)
    reversed_inventory = replace(
        inventory,
        series=tuple(reversed(inventory.series)),
        clinical=dict(reversed(tuple(inventory.clinical.items()))),
    )
    assert report == availability_report(reversed_inventory)
    serialized = json.dumps(report, allow_nan=False)
    for patient in inventory.patient_ids:
        assert patient not in serialized
    for series in inventory.series:
        assert series.series_uid not in serialized and series.study_uid not in serialized
    assert "2024-12-22" not in serialized
    assert "clinical_available_on" not in serialized


def test_private_manifest_traces_sources_and_keeps_modality_partition_together(inventory):
    manifest = private_manifest(inventory)
    hashes = {receipt.sha256 for receipt in inventory.receipts}
    assert manifest["privacy"] == "PRIVATE_DO_NOT_COMMIT"
    assert len(manifest["records"]) == 12
    for record in manifest["records"]:
        assert record["partition"] in {"train", "validation", "test"}
        assert all(link["ref"]["source_sha256"] in hashes for link in record["source_links"])
        assert all(link.get("clinical_available_on") is None for link in record["source_links"])


def test_aggregate_flow_uses_exclusive_complete_and_incomplete_groups(inventory):
    flow = flow_markdown(availability_report(inventory))
    assert "All four source types: 2" in flow
    assert "At least one source type absent: 10" in flow
    assert "not established" in flow
    assert all(patient not in flow for patient in inventory.patient_ids)


def test_CLI_creates_a_local_audit_and_preserves_the_blocked_readiness_status(
    tmp_path,
    inputs,
    monkeypatch,
):
    from oncobio_agentfm.cohort_audit import main

    paths, lock = inputs
    monkeypatch.chdir(tmp_path)
    argv = [
        "cohort-audit",
        "--lock",
        str(lock),
        "--output-dir",
        "outputs/review",
        "--local-manifest",
    ]
    for kind, path in paths.items():
        argv.extend((f"--{kind.replace('_', '-')}", str(path)))
    monkeypatch.setattr("sys.argv", argv)
    assert main() == 2
    public = tmp_path / "outputs/review/source_availability.json"
    private = tmp_path / "outputs/review/PRIVATE_cohort_manifest.json"
    assert "SYN-cohort-0" not in public.read_text()
    assert "SYN-cohort-0" in private.read_text()
    before = public.read_bytes()
    with pytest.raises(FileExistsError):
        main()
    assert public.read_bytes() == before


def test_CLI_refuses_patient_manifest_outside_local_outputs(tmp_path, monkeypatch):
    from oncobio_agentfm.cohort_audit import main

    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(
        "sys.argv",
        [
            "cohort-audit",
            "--clinical",
            "missing",
            "--imaging",
            "missing",
            "--rna-metadata",
            "missing",
            "--rna-expression",
            "missing",
            "--lock",
            "missing",
            "--output-dir",
            "docs/evidence",
        ],
    )
    with pytest.raises(ProtocolViolation, match="audit_output_must_be_under_local_outputs"):
        main()
    assert not (tmp_path / "docs").exists()


def test_linked_cohort_reuses_survival_gates_and_reports_modality_event_budget(
    inventory, linked, protocol
):
    cases, reviews, documents = linked
    report = audit_linked_cases(
        cases, inventory=inventory, reviews=reviews, documents=documents, protocol=protocol
    )
    assert report["input_patients"] == report["lineage_verified_patients"] == 12
    assert report["settings"]["pretreatment"]["endpoint_eligible"] == 12
    primary_events = sum(
        part["deaths"] for part in report["settings"]["pretreatment"]["partition_counts"].values()
    )
    assert primary_events == 6
    assert report["all_four_landmark_patients"] == 2
    assert (
        sum(part["deaths"] for part in report["all_four_landmark"]["partition_counts"].values())
        == 0
    )
    assert "SYN-cohort" not in json.dumps(report)


def test_missing_case_review_cannot_silently_drop_a_patient(inventory, linked, protocol):
    cases, reviews, documents = linked
    with pytest.raises(ProtocolViolation, match="case_review_coverage_mismatch"):
        audit_linked_cases(
            cases, inventory=inventory, reviews={}, documents=documents, protocol=protocol
        )
    with pytest.raises(ProtocolViolation, match="case_inventory_coverage_incomplete"):
        audit_linked_cases(
            cases[1:],
            inventory=inventory,
            reviews={
                case.context.baseline.patient_id: reviews[case.context.baseline.patient_id]
                for case in cases[1:]
            },
            documents=documents,
            protocol=protocol,
        )


def test_duplicate_case_cannot_inflate_clinical_events(inventory, linked, protocol):
    cases, reviews, documents = linked
    with pytest.raises(ProtocolViolation, match="linked_duplicate_patient"):
        audit_linked_cases(
            (*cases, cases[0]),
            inventory=inventory,
            reviews=reviews,
            documents=documents,
            protocol=protocol,
        )


def test_review_and_supporting_document_tampering_are_distinct_failures(inventory, linked):
    cases, reviews, documents = linked
    path, fingerprint = reviews[cases[0].context.baseline.patient_id]
    with pytest.raises(ProtocolViolation, match="case_review_hash_mismatch"):
        verify_case_review(
            cases[0],
            path,
            expected_sha256="0" * 64,
            inventory=inventory,
            verified_documents=verify_documents(documents),
        )
    with pytest.raises(ProtocolViolation, match="case_lineage_source_unverified"):
        verify_case_review(
            cases[0],
            path,
            expected_sha256=fingerprint,
            inventory=inventory,
            verified_documents={},
        )


def test_selected_image_source_locator_must_resolve_the_correct_metadata_row(inventory, linked):
    cases, reviews, documents = linked
    case = cases[0]
    path, _ = reviews[case.context.baseline.patient_id]
    review = json.loads(path.read_bytes())
    review["references"][case.primary_features[0].evidence.source_id]["locator"] = "/999999"
    fingerprint = rewrite_review(path, review)
    with pytest.raises(ProtocolViolation, match="selected_image_locator_mismatch"):
        verify_case_review(
            case,
            path,
            expected_sha256=fingerprint,
            inventory=inventory,
            verified_documents=verify_documents(documents),
        )


def test_calendar_mapping_must_preserve_source_and_aligned_dates(inventory, linked):
    cases, reviews, documents = linked
    case = cases[0]
    path, _ = reviews[case.context.baseline.patient_id]
    review = json.loads(path.read_bytes())
    review["image_date_mappings"][case.primary_features[0].evidence.source_id][
        "source_axis_date"
    ] = "2000-01-01"
    fingerprint = rewrite_review(path, review)
    with pytest.raises(ProtocolViolation, match="selected_image_calendar_mapping_mismatch"):
        verify_case_review(
            case,
            path,
            expected_sha256=fingerprint,
            inventory=inventory,
            verified_documents=verify_documents(documents),
        )


def test_reviewed_case_values_cannot_change_without_source_review(inventory, linked):
    cases, reviews, documents = linked
    case = cases[0]
    changed = replace(
        case, context=replace(case.context, baseline=replace(case.context.baseline, age_years=90))
    )
    path, fingerprint = reviews[case.context.baseline.patient_id]
    with pytest.raises(ProtocolViolation, match="case_differs_from_reviewed_source"):
        verify_case_review(
            changed,
            path,
            expected_sha256=fingerprint,
            inventory=inventory,
            verified_documents=verify_documents(documents),
        )


def test_preconstruction_exclusion_requires_a_verified_patient_specific_record(
    tmp_path, inventory, linked, protocol
):
    cases, reviews, documents = linked
    excluded = cases[0].context.baseline.patient_id
    document = tmp_path / "SYN-exclusion.json"
    document.write_text(
        json.dumps({"entry": {"patient_id": excluded, "reason": "source_adjudication_unresolved"}}),
        encoding="utf-8",
    )
    fingerprint = hashlib.sha256(document.read_bytes()).hexdigest()
    documents = {**documents, fingerprint: document}
    report = audit_linked_cases(
        cases[1:],
        inventory=inventory,
        reviews={key: value for key, value in reviews.items() if key != excluded},
        documents=documents,
        protocol=protocol,
        exclusions=(
            SourceExclusion(excluded, "source_adjudication_unresolved", fingerprint, "entry"),
        ),
    )
    assert report["source_exclusions"] == {"source_adjudication_unresolved": 1}
    assert report["input_patients"] == 12
    assert report["lineage_verified_patients"] == 11


def test_all_patients_can_be_accounted_for_by_documented_source_exclusions(
    tmp_path,
    inventory,
    protocol,
):
    document = tmp_path / "SYN-all-exclusions.json"
    entries = {
        f"entry-{number}": {"patient_id": patient, "reason": "source_adjudication_unresolved"}
        for number, patient in enumerate(sorted(inventory.patient_ids))
    }
    document.write_text(json.dumps(entries), encoding="utf-8")
    fingerprint = hashlib.sha256(document.read_bytes()).hexdigest()
    exclusions = tuple(
        SourceExclusion(value["patient_id"], value["reason"], fingerprint, key)
        for key, value in entries.items()
    )
    report = audit_linked_cases(
        (),
        inventory=inventory,
        reviews={},
        documents={fingerprint: document},
        protocol=protocol,
        exclusions=exclusions,
    )
    assert report["source_exclusions"] == {"source_adjudication_unresolved": 12}
    assert report["validated_event_budget"]["pretreatment"]["eligible_patients"] == 0
    assert report["lineage_verified_patients"] == 0


def test_exclusion_cannot_be_both_a_case_or_an_unsupported_patient(
    inventory,
    linked,
    protocol,
):
    cases, reviews, documents = linked
    first = cases[0].context.baseline.patient_id
    digest = next(iter(documents))
    exclusion = SourceExclusion(first, "source_adjudication_unresolved", digest, "entry")
    with pytest.raises(ProtocolViolation, match="source_exclusion_identity_overlap"):
        audit_linked_cases(
            cases,
            inventory=inventory,
            reviews=reviews,
            documents=documents,
            protocol=protocol,
            exclusions=(exclusion,),
        )


def test_outcome_identity_mismatch_is_counted_after_predictor_eligibility(
    inventory, linked, protocol
):
    cases, reviews, documents = linked
    original = cases[0]
    changed = replace(
        original, follow_up=replace(original.follow_up, patient_id="SYN-wrong-outcome")
    )
    from pathlib import Path

    directory = next(iter(documents.values())).parent
    reviews, documents = write_synthetic_reviews(Path(directory), inventory, (changed, *cases[1:]))
    result = audit_linked_cases(
        (changed, *cases[1:]),
        inventory=inventory,
        reviews=reviews,
        documents=documents,
        protocol=protocol,
    )["settings"]["pretreatment"]
    assert result["predictor_eligible"] == 12
    assert result["endpoint_eligible"] == 11
    assert result["outcome_rejections"] == {"outcome_patient_mismatch": 1}


def test_synthetic_snapshot_is_recomputed_offline_and_separate_from_real_inventory(protocol):
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    observed = validation_report(protocol)
    assert observed["status"] == "PASS"
    assert all(check["passed"] for check in observed["checks"])
    assert observed == json.loads((root / "docs/evidence/cohort_validation.json").read_bytes())
    real = json.loads((root / "docs/evidence/source_availability.json").read_bytes())
    assert real["evidence_type"] == "observed_public_source_inventory"
    assert real["patients"] == 211
    assert real["modality_counts"] == {"ct": 211, "pet": 201, "clinical": 211, "rna": 130}
    assert real["clinical_eligibility"] == "NOT_ESTABLISHED"
    assert real["analysis_populations"]["pretreatment"]["survival_events"] is None
    lock = json.loads((root / "configs/cohort/nsclc_source_lock_v1.json").read_bytes())
    assert {
        receipt["kind"]: {
            key: value for key, value in receipt.items() if key not in {"kind", "reference"}
        }
        for receipt in real["source_receipts"]
    } == lock["sources"]
