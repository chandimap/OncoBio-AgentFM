import csv
import json
from dataclasses import replace
from pathlib import Path

import pytest

from oncobio_agentfm.prognosis import ProtocolViolation, VitalStatus
from oncobio_agentfm.protocol_audit import (
    audit_eligibility,
    build_primary_clinical_cohort,
    write_report,
)
from oncobio_agentfm.protocol_validation import (
    make_synthetic_case,
    source_integrity_checks,
    validation_report,
)
from oncobio_agentfm.research_protocol import PredictionSetting
from oncobio_agentfm.source_audit import REVIEWED_SOURCE_COLUMNS, audit_source_schema
from oncobio_agentfm.splitting import split_patient_ids

ROOT = Path(__file__).resolve().parents[1]


def test_committed_adversarial_report_matches_recomputed_evidence(protocol) -> None:
    observed = validation_report(protocol)
    recorded = json.loads((ROOT / "docs/evidence/temporal_validation.json").read_text())
    assert observed == recorded
    assert observed["status"] == "PASS"
    assert observed["real_data_temporal_eligibility"] == "NOT_ESTABLISHED"
    assert observed["model_performance"] == "NOT_EVALUATED"


def test_duplicate_patients_cannot_inflate_eligibility_counts(protocol) -> None:
    case = make_synthetic_case()
    with pytest.raises(ProtocolViolation, match="duplicate_patient"):
        audit_eligibility((case, case), setting=PredictionSetting.PRETREATMENT, protocol=protocol)


def test_audit_is_order_independent_and_does_not_export_identifiers_sources_or_dates(protocol):
    cases = tuple(make_synthetic_case(f"SYN-private-id-{number}") for number in range(8))
    first = audit_eligibility(cases, setting=PredictionSetting.PRETREATMENT, protocol=protocol)
    second = audit_eligibility(
        tuple(reversed(cases)), setting=PredictionSetting.PRETREATMENT, protocol=protocol
    )
    assert first == second
    serialized = json.dumps(first)
    for case in cases:
        assert case.context.baseline.patient_id not in serialized
        assert case.follow_up.endpoint_source not in serialized
    assert "2025-01-01" not in serialized
    assert first["endpoint_eligible"] == 8


def test_predictor_eligibility_is_unchanged_by_endpoint_status(protocol) -> None:
    case = make_synthetic_case()
    censored = replace(
        case, follow_up=replace(case.follow_up, status=VitalStatus.ALIVE, death_on=None)
    )
    first = audit_eligibility((case,), setting=PredictionSetting.PRETREATMENT, protocol=protocol)
    second = audit_eligibility(
        (censored,), setting=PredictionSetting.PRETREATMENT, protocol=protocol
    )
    assert first["predictor_eligible"] == second["predictor_eligible"] == 1
    assert first["predictor_rejections"] == second["predictor_rejections"] == {}


def test_gated_cohort_reuses_existing_contract_and_split_without_seed_search(protocol) -> None:
    cases = tuple(make_synthetic_case(f"SYN-{number}") for number in range(12))
    cohort = build_primary_clinical_cohort(cases, protocol=protocol, ajcc_version=9)
    assert cohort.patient_ids == {case.context.baseline.patient_id for case in cases}
    assert cohort.audit().deaths == 12
    split = split_patient_ids(cohort.patient_ids)
    subset = split_patient_ids(tuple(sorted(cohort.patient_ids))[:8])
    for name in ("train", "validation", "test"):
        assert set(getattr(subset, name)) <= set(getattr(split, name))
    invalid = replace(cases[0], context=replace(cases[0].context, calendar_alignment_source=None))
    with pytest.raises(ProtocolViolation, match="calendar_alignment_unverified"):
        build_primary_clinical_cohort((invalid,), protocol=protocol, ajcc_version=9)


def test_report_creation_is_deterministic_and_refuses_overwrite(tmp_path, protocol) -> None:
    report = validation_report(protocol)
    path = tmp_path / "nested/report.json"
    write_report(path, report)
    assert json.loads(path.read_text()) == report
    with pytest.raises(FileExistsError):
        write_report(path, {"status": "different"})
    assert json.loads(path.read_text()) == report


def write_source_fixture(path, identifiers=("SYN-A", "SYN-B")) -> None:
    headers = REVIEWED_SOURCE_COLUMNS
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=headers)
        writer.writeheader()
        for identifier in identifiers:
            writer.writerow(
                {
                    "Case ID": identifier,
                    "CT Date": "SYNTHETIC-DATE",
                    "Survival Status": "SYNTHETIC-STATUS",
                    "Pathological T stage": "T1a",
                }
            )


def test_schema_audit_never_claims_real_temporal_validation_or_exports_row_values(tmp_path) -> None:
    path = tmp_path / "source.csv"
    write_source_fixture(path)
    report = audit_source_schema(path)
    serialized = json.dumps(report)
    assert "SYN-A" not in serialized and "SYNTHETIC-DATE" not in serialized
    assert "SYNTHETIC-STATUS" not in serialized
    assert report["real_data_temporal_eligibility"] == "NOT_ESTABLISHED"
    assert report["primary_prediction_readiness"] == "BLOCKED_PENDING_SOURCE_ADJUDICATION"
    assert report["raw_pathological_stage_vocabulary"]["Pathological T stage"] == ["T1a"]


def test_source_audit_rejects_duplicate_patients_and_broken_schema(tmp_path) -> None:
    path = tmp_path / "source.csv"
    write_source_fixture(path, identifiers=("SYN-A", "SYN-A"))
    with pytest.raises(ValueError, match="duplicate patient"):
        audit_source_schema(path)
    path.write_text("Case ID,Case ID\nSYN-A,SYN-A\n")
    with pytest.raises(ValueError, match="unique headers"):
        audit_source_schema(path)


def test_unrecognized_stage_text_is_withheld_from_aggregate_reports(tmp_path) -> None:
    path = tmp_path / "source.csv"
    write_source_fixture(path)
    path.write_text(path.read_text().replace("T1a", "SYN-private-source-value"))
    report = audit_source_schema(path)
    assert "SYN-private-source-value" not in json.dumps(report)
    assert report["raw_pathological_stage_vocabulary"]["Pathological T stage"] == [
        "<unrecognized-token-withheld>"
    ]


def test_source_snapshot_records_schema_limits_without_patient_rows() -> None:
    report = json.loads((ROOT / "docs/evidence/source_schema_audit.json").read_text())
    assert report["rows"] == 211
    assert (
        report["source_sha256"]
        == "0a1f7cdf02cef58450db17f116a3fb0ecfd4d737ca6843bb1939dbe421071476"
    )
    assert "Pathological T stage" in report["columns"]
    assert not any("AJCC" in column or "ECOG" in column for column in report["columns"])
    assert "patients" not in report and "records" not in report


def test_unknown_header_requires_review_without_exporting_it(tmp_path):
    path = tmp_path / "source.csv"
    write_source_fixture(path)
    path.write_text(path.read_text().replace("CT Date", "SYN-private-header"))
    with pytest.raises(ValueError, match="reviewed TCIA header schema") as error:
        audit_source_schema(path)
    assert "SYN-private-header" not in str(error.value)


def test_invalid_report_does_not_leave_an_empty_evidence_file(tmp_path):
    path = tmp_path / "invalid.json"
    with pytest.raises(ValueError, match="Out of range float"):
        write_report(path, {"not_valid": float("nan")})
    assert not path.exists()


def test_paired_modality_comparisons_audit_their_own_identical_split_populations(protocol):
    from oncobio_agentfm.contracts import EvidenceRecord, MissingnessReason
    from oncobio_agentfm.prognosis import AvailabilityBasis, PredictorFeature
    from oncobio_agentfm.research_protocol import FeatureKind

    cases = []
    pet_ids = set()
    molecular_ids = set()
    for number in range(24):
        case = make_synthetic_case(f"SYN-paired-{number}")
        ct, pet = case.primary_features
        molecular = case.landmark_features[-1]
        if number % 2 == 0:
            pet = PredictorFeature(
                FeatureKind.PET_CT,
                ct.evidence,
                "synthetic-pet-timing",
                AvailabilityBasis.CLINICAL_RECORD,
            )
            pet_ids.add(case.context.baseline.patient_id)
        if number % 3 != 0:
            absent = EvidenceRecord(
                case.context.baseline.patient_id, None, None, None, MissingnessReason.NOT_ACQUIRED
            )
            molecular = PredictorFeature(FeatureKind.RESECTION_MOLECULAR, absent, None)
        else:
            molecular_ids.add(case.context.baseline.patient_id)
        cases.append(
            replace(case, primary_features=(ct, pet), landmark_features=(ct, pet, molecular))
        )
    primary = audit_eligibility(
        tuple(cases), setting=PredictionSetting.PRETREATMENT, protocol=protocol
    )
    landmark = audit_eligibility(
        tuple(cases), setting=PredictionSetting.POSTOPERATIVE, protocol=protocol
    )
    assert primary["endpoint_eligible"] == landmark["endpoint_eligible"] == 24
    for audit, key, ids in (
        (primary, "pet_available", pet_ids),
        (landmark, "molecular_available", molecular_ids),
    ):
        comparison = audit["comparison_feasibility"][key]
        split = split_patient_ids(ids)
        for partition in ("train", "validation", "test"):
            assert comparison["partition_counts"][partition]["patients"] == len(
                getattr(split, partition)
            )
        assert comparison["precision_status"] == "NOT_ESTABLISHED"
        assert comparison["model_evaluation"] == "NOT_RUN"
    serialized = json.dumps(primary) + json.dumps(landmark)
    assert "SYN-paired-" not in serialized


def test_empty_pet_subset_is_reported_without_replacing_test_patients_or_split_seed(protocol):
    report = audit_eligibility(
        (make_synthetic_case(),), setting=PredictionSetting.PRETREATMENT, protocol=protocol
    )
    pet = report["comparison_feasibility"]["pet_available"]
    assert sum(part["patients"] for part in pet["partition_counts"].values()) == 0
    assert all(metric["status"] == "NOT_ESTIMABLE" for metric in pet["metric_support"]["test"])


def test_truthy_text_cannot_mislabel_real_audit_as_synthetic(protocol):
    with pytest.raises(ProtocolViolation, match="evidence_scope_flag_invalid"):
        audit_eligibility(
            (make_synthetic_case(),),
            setting=PredictionSetting.PRETREATMENT,
            protocol=protocol,
            synthetic="false",
        )


def test_raw_release_budget_is_not_eligible_or_730_day_event_evidence(tmp_path):
    path = tmp_path / "source.csv"
    write_source_fixture(path)
    report = audit_source_schema(path)
    budget = report["raw_source_budget"]
    assert budget["subjects"] == 2
    assert budget["unrecognized_status"] == 2
    assert budget["recorded_dead"] == budget["recorded_alive"] == 0
    assert budget["precision_status"] == "NOT_ESTABLISHED"
    assert report["real_data_temporal_eligibility"] == "NOT_ESTABLISHED"
    assert "SYNTHETIC-STATUS" not in json.dumps(budget)


@pytest.mark.parametrize("kind_name", ["CT", "PET_CT"])
def test_landmark_cannot_swap_frozen_image_source_or_introduce_baseline_pet(kind_name, protocol):
    from oncobio_agentfm.research_protocol import FeatureKind

    kind = FeatureKind[kind_name]
    case = make_synthetic_case()
    template = case.primary_features[0]
    other = replace(
        template,
        kind=kind,
        evidence=replace(template.evidence, source_id="synthetic-different-old-scan"),
    )
    changed = replace(
        case,
        landmark_features=tuple(
            other if item.kind is kind else item for item in case.landmark_features
        ),
    )
    primary = audit_eligibility(
        (changed,), setting=PredictionSetting.PRETREATMENT, protocol=protocol
    )
    landmark = audit_eligibility(
        (changed,), setting=PredictionSetting.POSTOPERATIVE, protocol=protocol
    )
    assert primary["endpoint_eligible"] == 1
    assert landmark["predictor_rejections"] == {"baseline_imaging_not_frozen_across_settings": 1}
    assert landmark["endpoint_eligible"] == 0


def test_cohort_refuses_mixed_cutoff_policies_before_outcome_analysis(protocol):
    first = make_synthetic_case("SYN-policy-first")
    other = make_synthetic_case("SYN-policy-other")
    other = replace(
        other,
        follow_up=replace(other.follow_up, administrative_cutoff_source="synthetic-other-policy"),
    )
    with pytest.raises(ProtocolViolation, match="cohort_administrative_cutoff_policy_mismatch"):
        audit_eligibility((first, other), setting=PredictionSetting.PRETREATMENT, protocol=protocol)
    with pytest.raises(ProtocolViolation, match="cohort_administrative_cutoff_policy_mismatch"):
        build_primary_clinical_cohort((first, other), protocol=protocol, ajcc_version=9)


def test_distinct_patient_cutoff_dates_can_follow_one_reviewed_policy(protocol):
    first = make_synthetic_case("SYN-cutoff-first")
    other = make_synthetic_case("SYN-cutoff-other")
    from datetime import date

    other = replace(
        other, follow_up=replace(other.follow_up, administrative_cutoff=date(2027, 1, 1))
    )
    result = audit_eligibility(
        (first, other), setting=PredictionSetting.PRETREATMENT, protocol=protocol
    )
    assert result["endpoint_eligible"] == 2
    assert result["clinical_validity"] == "NOT_ESTABLISHED_BY_SOFTWARE_CHECKS"


def test_source_consistency_probes_are_inspectable_and_pass(protocol):
    checks = source_integrity_checks(protocol)
    assert len(checks) == 5
    assert all(check["passed"] for check in checks)
    assert any(
        check["observed_rejection"] == "documented_treatment_after_death" for check in checks
    )
    assert any(check["observed_rejection"] == "documented_surgery_after_death" for check in checks)
