"""Deterministic synthetic adversarial evidence for the reviewed research protocol."""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass, replace
from datetime import date, timedelta
from pathlib import Path

from oncobio_agentfm.clinical import ClinicalStage, ECOGStatus, NSCLCBaseline, StageGroup
from oncobio_agentfm.contracts import EvidenceRecord, MissingnessReason, SurvivalRecord
from oncobio_agentfm.evaluation_support import assess_metric_support
from oncobio_agentfm.prognosis import (
    AvailabilityBasis,
    OverallSurvivalFollowUp,
    PredictionContext,
    PredictorFeature,
    ProtocolViolation,
    StageBasis,
    StageProvenance,
    VitalStatus,
    validate_predictors,
)
from oncobio_agentfm.protocol_audit import AdjudicatedCase, audit_eligibility, write_report
from oncobio_agentfm.research_protocol import (
    FeatureKind,
    PredictionSetting,
    ResearchProtocol,
    SurvivalMetric,
    load_research_protocol,
)


def make_synthetic_case(patient_id: str = "SYNTHETIC-ONLY") -> AdjudicatedCase:
    """Invent an explicitly synthetic patient; this is not a TCIA importer."""

    index = date(2025, 1, 1)

    def present(role: str, acquired: int, available: int) -> EvidenceRecord:
        return EvidenceRecord(
            patient_id,
            f"synthetic:{role}:{patient_id}",
            index + timedelta(days=acquired),
            index + timedelta(days=available),
        )

    missing = EvidenceRecord(patient_id, None, None, None, MissingnessReason.NOT_ACQUIRED)
    stage_evidence = present("clinical-stage", -5, -3)
    baseline = NSCLCBaseline(
        patient_id=patient_id,
        index_date=index,
        diagnosis=present("biopsy-diagnosis", -20, -15),
        age_years=65,
        age_evidence=present("age-at-assessment", 0, 0),
        stage=ClinicalStage(stage_evidence, StageGroup.TWO, 9),
        ecog=ECOGStatus(missing, None),
    )
    context = PredictionContext(
        baseline=baseline,
        assessment=present("pretreatment-assessment", 0, 0),
        surgical_referral=present("referral", -2, -1),
        untreated_at_assessment=present("no-prior-treatment", 0, 0),
        calendar_alignment_source="synthetic-shared-calendar-v1",
        clinical_timing_source="synthetic-diagnosis-age-and-clinical-timing-review",
        age_reference_date=index,
        index_selection_source="synthetic-earliest-qualifying-assessment-review",
        first_treatment_on=index + timedelta(days=10),
        surgery=present("first-treatment-resection", 10, 11),
        stage_provenance=StageProvenance(
            stage_evidence, "II", StageBasis.CLINICAL, "AJCC", 9, StageGroup.TWO, "synthetic-stage"
        ),
    )
    ct = PredictorFeature(
        FeatureKind.CT,
        present("ct", -10, -8),
        "synthetic-ct-access",
        AvailabilityBasis.CLINICAL_RECORD,
    )
    pet = PredictorFeature(FeatureKind.PET_CT, missing, None)
    molecular = PredictorFeature(
        FeatureKind.RESECTION_MOLECULAR,
        present("resection-specimen", 10, 25),
        "synthetic-result",
        AvailabilityBasis.CLINICAL_RECORD,
    )
    follow_up = OverallSurvivalFollowUp(
        patient_id=patient_id,
        status=VitalStatus.DEAD,
        death_on=index + timedelta(days=900),
        last_known_alive_on=index + timedelta(days=850),
        administrative_cutoff=index + timedelta(days=1000),
        endpoint_source="synthetic-all-cause-follow-up",
        calendar_alignment_source="synthetic-shared-calendar-v1",
        administrative_cutoff_source="synthetic-prespecified-cutoff-policy-v1",
    )
    return AdjudicatedCase(context, (ct, pet), (ct, pet, molecular), follow_up)


@dataclass(frozen=True, slots=True)
class ValidationProbe:
    """A synthetic case and prespecified first-failure codes for both settings."""

    name: str
    case: AdjudicatedCase
    expected_primary: str | None
    expected_secondary: str | None


def synthetic_probes() -> tuple[ValidationProbe, ...]:
    """Challenge structural leakage as well as date comparisons."""

    probes = []

    def add(name: str, case: AdjudicatedCase, primary: str | None, secondary: str | None) -> None:
        probes.append(ValidationProbe(name, case, primary, secondary))

    positive = make_synthetic_case("SYN-positive")
    add("valid_timed_evidence", positive, None, None)
    late_ct_case = make_synthetic_case("SYN-late-ct")
    ct, pet = late_ct_case.primary_features
    late_ct = replace(ct, evidence=replace(ct.evidence, available_on=date(2025, 1, 2)))
    add(
        "ct_result_available_after_index",
        replace(
            late_ct_case,
            primary_features=(late_ct, pet),
            landmark_features=(late_ct, pet, late_ct_case.landmark_features[-1]),
        ),
        "post_index_evidence",
        "post_index_evidence",
    )
    for name, basis, edition, code in (
        ("pathological_stage_at_baseline", StageBasis.PATHOLOGICAL, 9, "nonclinical_stage"),
        (
            "historical_stage_not_relabelled",
            StageBasis.CLINICAL,
            7,
            "unsupported_or_unknown_stage_edition",
        ),
        (
            "unknown_stage_edition",
            StageBasis.CLINICAL,
            None,
            "unsupported_or_unknown_stage_edition",
        ),
    ):
        case = make_synthetic_case(f"SYN-{name}")
        original_provenance = case.context.stage_provenance
        if original_provenance is None:
            raise RuntimeError("synthetic stage provenance is missing")
        provenance = replace(original_provenance, basis=basis, edition=edition)
        add(
            name,
            replace(case, context=replace(case.context, stage_provenance=provenance)),
            code,
            code,
        )
    case = make_synthetic_case("SYN-forbidden-molecular")
    add(
        "resection_molecular_forbidden_at_pretreatment",
        replace(case, primary_features=case.landmark_features),
        "forbidden_predictor_role",
        None,
    )
    case = make_synthetic_case("SYN-late-molecular")
    ct, pet, molecular = case.landmark_features
    late_molecular = replace(
        molecular, evidence=replace(molecular.evidence, available_on=date(2025, 2, 11))
    )
    add(
        "molecular_result_after_fixed_landmark",
        replace(case, landmark_features=(ct, pet, late_molecular)),
        None,
        "post_index_evidence",
    )
    case = make_synthetic_case("SYN-death-before-landmark")
    add(
        "death_before_landmark_not_immortal_time",
        replace(
            case,
            follow_up=replace(
                case.follow_up, death_on=date(2025, 1, 21), last_known_alive_on=date(2025, 1, 20)
            ),
        ),
        None,
        "no_positive_follow_up_from_prediction_time",
    )
    case = make_synthetic_case("SYN-calendar-unknown")
    add(
        "unverified_shifted_calendar",
        replace(case, context=replace(case.context, calendar_alignment_source=None)),
        "calendar_alignment_unverified",
        "calendar_alignment_unverified",
    )
    case = make_synthetic_case("SYN-outcome-axis")
    add(
        "follow_up_from_different_calendar_axis",
        replace(
            case, follow_up=replace(case.follow_up, calendar_alignment_source="different-axis")
        ),
        "outcome_calendar_alignment_mismatch",
        "outcome_calendar_alignment_mismatch",
    )
    for name, field, value, code in (
        (
            "unreviewed_clinical_timing",
            "clinical_timing_source",
            None,
            "clinical_timing_unverified",
        ),
        (
            "age_from_another_event",
            "age_reference_date",
            date(2024, 12, 1),
            "age_not_referenced_to_index",
        ),
        (
            "same_day_treatment",
            "first_treatment_on",
            date(2025, 1, 1),
            "treatment_before_or_on_index",
        ),
    ):
        case = make_synthetic_case(f"SYN-{name}")
        add(name, replace(case, context=replace(case.context, **{field: value})), code, code)
    case = make_synthetic_case("SYN-research-assay")
    ct, pet, molecular = case.landmark_features
    research_only = replace(molecular, availability_basis=AvailabilityBasis.RESEARCH_ONLY)
    add(
        "research_only_assay_not_clinical_result",
        replace(case, landmark_features=(ct, pet, research_only)),
        None,
        "feature_clinical_availability_unverified",
    )
    case = make_synthetic_case("SYN-death-on-landmark")
    add(
        "death_on_landmark_day",
        replace(
            case,
            follow_up=replace(case.follow_up, death_on=date(2025, 2, 10), last_known_alive_on=None),
        ),
        None,
        "no_positive_follow_up_from_prediction_time",
    )
    case = make_synthetic_case("SYN-follow-up-before-landmark")
    add(
        "last_known_alive_before_landmark",
        replace(
            case,
            follow_up=replace(
                case.follow_up,
                status=VitalStatus.ALIVE,
                death_on=None,
                last_known_alive_on=date(2025, 1, 20),
            ),
        ),
        None,
        "no_positive_follow_up_from_prediction_time",
    )
    case = make_synthetic_case("SYN-outcome-as-feature")
    ct, pet = case.primary_features
    outcome_feature = replace(ct, kind=FeatureKind.SURVIVAL_OUTCOME)
    add(
        "outcome_predictor_rejected_despite_early_date",
        replace(
            case,
            primary_features=(ct, pet, outcome_feature),
            landmark_features=(*case.landmark_features, outcome_feature),
        ),
        "forbidden_predictor_role",
        "forbidden_predictor_role",
    )
    for kind in (FeatureKind.CT, FeatureKind.PET_CT):
        for acquired, available, label in (
            (1, 2, "acquired_after_assessment"),
            (-5, 2, "access_after_assessment"),
        ):
            case = make_synthetic_case(f"SYN-frozen-{kind.value}-{label}")
            late = replace(
                case.primary_features[0],
                kind=kind,
                evidence=EvidenceRecord(
                    case.context.baseline.patient_id,
                    "synthetic-pre-surgery-image",
                    date(2025, 1, 1) + timedelta(days=acquired),
                    date(2025, 1, 1) + timedelta(days=available),
                ),
            )
            primary = tuple(late if item.kind is kind else item for item in case.primary_features)
            landmark = tuple(late if item.kind is kind else item for item in case.landmark_features)
            add(
                f"frozen_{kind.value}_{label}",
                replace(case, primary_features=primary, landmark_features=landmark),
                "post_index_evidence",
                "post_index_evidence",
            )
    for kind in (FeatureKind.CT, FeatureKind.PET_CT):
        case = make_synthetic_case(f"SYN-changed-baseline-{kind.value}")
        template = case.primary_features[0]
        other = replace(
            template,
            kind=kind,
            evidence=replace(template.evidence, source_id=f"synthetic-other-baseline-{kind.value}"),
        )
        changed = tuple(other if item.kind is kind else item for item in case.landmark_features)
        add(
            f"changed_frozen_baseline_{kind.value}",
            replace(case, landmark_features=changed),
            None,
            "baseline_imaging_not_frozen_across_settings",
        )
    for label, field, value, code in (
        (
            "stale_untreated_status",
            "untreated_at_assessment",
            replace(
                positive.context.untreated_at_assessment,
                acquired_on=date(2024, 12, 1),
                available_on=date(2024, 12, 1),
            ),
            "untreated_status_not_verified_at_index",
        ),
        (
            "first_treatment_after_surgery",
            "first_treatment_on",
            date(2025, 1, 20),
            "first_treatment_after_documented_surgery",
        ),
    ):
        case = make_synthetic_case(f"SYN-{label}")
        if isinstance(value, EvidenceRecord):
            value = replace(
                value,
                patient_id=case.context.baseline.patient_id,
                source_id=case.context.untreated_at_assessment.source_id,
            )
        add(label, replace(case, context=replace(case.context, **{field: value})), code, code)
    for kind in (FeatureKind.CT, FeatureKind.PET_CT):
        for duplicate in (False, True):
            label = "duplicate" if duplicate else "undeclared"
            case = make_synthetic_case(f"SYN-original-{kind.value}-{label}")
            selected = next(item for item in case.primary_features if item.kind is kind)
            original = (
                (*case.primary_features, selected)
                if duplicate
                else tuple(item for item in case.primary_features if item.kind is not kind)
            )
            add(
                f"original_{kind.value}_{label}",
                replace(case, primary_features=original),
                "duplicate_feature_role" if duplicate else "modality_missingness_undeclared",
                "frozen_baseline_imaging_not_uniquely_declared",
            )
    return tuple(probes)


def source_integrity_checks(protocol: ResearchProtocol) -> list[dict[str, object]]:
    """Challenge malformed source bundles before outcome analysis, without real records."""

    checks = []
    for name, expected in (
        ("missing_cutoff_policy", "administrative_cutoff_source_missing"),
        ("mixed_cutoff_policies", "cohort_administrative_cutoff_policy_mismatch"),
        ("treatment_after_death", "documented_treatment_after_death"),
        ("surgery_after_death", "documented_surgery_after_death"),
        ("same_policy_with_different_adjudicated_cutoffs", None),
    ):
        case = make_synthetic_case("SYN-source-integrity")
        other = make_synthetic_case("SYN-source-integrity-other")
        observed = None
        try:
            if name == "missing_cutoff_policy":
                replace(case.follow_up, administrative_cutoff_source="")
            elif name == "mixed_cutoff_policies":
                other = replace(
                    other,
                    follow_up=replace(
                        other.follow_up, administrative_cutoff_source="synthetic-other-policy"
                    ),
                )
                audit_eligibility(
                    (case, other),
                    setting=PredictionSetting.PRETREATMENT,
                    protocol=protocol,
                    synthetic=True,
                )
            elif name == "treatment_after_death":
                replace(
                    case,
                    context=replace(case.context, first_treatment_on=date(2028, 1, 1)),
                )
            elif name == "surgery_after_death":
                surgery = case.context.surgery
                if surgery is None:
                    raise RuntimeError("synthetic surgery is missing")
                replace(
                    case,
                    context=replace(
                        case.context,
                        surgery=replace(
                            surgery, acquired_on=date(2028, 1, 1), available_on=date(2028, 1, 2)
                        ),
                    ),
                )
            else:
                other = replace(
                    other,
                    follow_up=replace(other.follow_up, administrative_cutoff=date(2027, 1, 1)),
                )
                audit_eligibility(
                    (case, other),
                    setting=PredictionSetting.PRETREATMENT,
                    protocol=protocol,
                    synthetic=True,
                )
        except ProtocolViolation as error:
            observed = error.code
        checks.append(
            {
                "probe": name,
                "expected_rejection": expected,
                "observed_rejection": observed,
                "passed": observed == expected,
            }
        )
    return checks


def metric_validation_checks(protocol: ResearchProtocol) -> list[dict[str, object]]:
    """Invent independent denominator/precision probes; no model performance is scored."""

    index = date(2025, 1, 1)

    def record(identifier: str, days: int, death: bool = True) -> SurvivalRecord:
        end = index + timedelta(days=days)
        return SurvivalRecord(identifier, index, end if death else None, None if death else end)

    training = (record("SYN-train-death", 300), record("SYN-train-censor", 1400, False))
    survivors = tuple(record(f"SYN-control-{i}", 1200, False) for i in range(10))
    ordinary = (record("SYN-event", 400), *survivors)
    boundary = (record("SYN-boundary-event", 730), *survivors)
    lost = tuple(record(f"SYN-lost-{i}", 100, False) for i in range(12))
    scenarios = (
        (
            "one_event_is_not_precision_evidence",
            SurvivalMetric.UNO,
            ordinary,
            "OUTCOME_PREREQUISITES_MET",
            {},
        ),
        ("uno_excludes_event_at_tau", SurvivalMetric.UNO, boundary, "NOT_ESTIMABLE", {}),
        (
            "auc_includes_event_at_horizon",
            SurvivalMetric.AUC,
            boundary,
            "OUTCOME_PREREQUISITES_MET",
            {},
        ),
        ("auc_without_cases", SurvivalMetric.AUC, survivors, "NOT_ESTIMABLE", {}),
        (
            "brier_without_early_deaths",
            SurvivalMetric.BRIER,
            survivors,
            "OUTCOME_PREREQUISITES_MET",
            {},
        ),
        (
            "ibs_complete_grid_without_early_deaths",
            SurvivalMetric.IBS,
            survivors,
            "OUTCOME_PREREQUISITES_MET",
            {},
        ),
        ("ibs_cannot_crop_unsupported_tail", SurvivalMetric.IBS, lost, "NOT_ESTIMABLE", {}),
        (
            "calibration_groups_not_assessed_from_pooled_counts",
            SurvivalMetric.CALIBRATION,
            ordinary,
            "NOT_ASSESSED",
            {},
        ),
        (
            "calibration_empty_tertiles",
            SurvivalMetric.CALIBRATION,
            ordinary,
            "NOT_ESTIMABLE",
            {"training_risk": (0.0, 1.0), "evaluation_risk": (0.1,) * len(ordinary)},
        ),
        (
            "small_comparable_set_withheld_by_protocol",
            SurvivalMetric.UNO,
            ordinary[:2],
            "WITHHELD_BY_PROTOCOL",
            {},
        ),
    )
    checks = []
    for name, metric, records, expected, kwargs in scenarios:
        observed = assess_metric_support(
            training, records, metric=metric, horizon_days=730, protocol=protocol, **kwargs
        )
        checks.append(
            {
                "probe": name,
                "expected_status": expected,
                "observed": observed,
                "passed": observed["status"] == expected
                and observed["precision_status"] == "NOT_ESTABLISHED"
                and observed["prediction_outputs_validated"] is False,
            }
        )
    return checks


def validation_report(protocol: ResearchProtocol) -> dict[str, object]:
    """Record expected versus observed rejection and aggregate eligibility."""

    probes = synthetic_probes()
    checks = []
    for probe in probes:
        for setting, expected in (
            (PredictionSetting.PRETREATMENT, probe.expected_primary),
            (PredictionSetting.POSTOPERATIVE, probe.expected_secondary),
        ):
            observed = None
            try:
                predictors = validate_predictors(
                    probe.case.context,
                    probe.case.features_for(setting),
                    setting=setting,
                    protocol=protocol,
                )
                probe.case.follow_up.anchor(predictors)
            except ProtocolViolation as exc:
                observed = exc.code
            checks.append(
                {
                    "probe": probe.name,
                    "setting": setting.value,
                    "expected_rejection": expected,
                    "observed_rejection": observed,
                    "passed": observed == expected,
                }
            )
    cases = tuple(probe.case for probe in probes)
    metric_checks = metric_validation_checks(protocol)
    integrity_checks = source_integrity_checks(protocol)
    return {
        "schema_version": "1.2",
        "protocol_id": protocol.protocol_id,
        "protocol_revision": protocol.revision,
        "protocol_sha256": protocol.sha256,
        "evidence_type": "synthetic_adversarial_validation_only",
        "status": "PASS"
        if all(check["passed"] for check in (*checks, *metric_checks, *integrity_checks))
        else "FAIL",
        "checks": checks,
        "metric_feasibility_checks": metric_checks,
        "source_integrity_checks": integrity_checks,
        "cohort_audits": [
            audit_eligibility(cases, setting=setting, protocol=protocol, synthetic=True)
            for setting in PredictionSetting
        ],
        "real_data_temporal_eligibility": "NOT_ESTABLISHED",
        "model_performance": "NOT_EVALUATED",
    }


def main() -> int:
    """Generate evidence or compare it byte-independently with a committed snapshot."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol", type=Path, default=Path("configs/research/nsclc_os_v1.json"))
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--output", type=Path)
    group.add_argument("--check-report", type=Path)
    args = parser.parse_args()
    report = validation_report(load_research_protocol(args.protocol))
    if args.output:
        write_report(args.output, report)
    if args.check_report and json.loads(args.check_report.read_text(encoding="utf-8")) != report:
        print("Protocol evidence differs from the committed snapshot.")
        return 1
    print(
        json.dumps(
            {
                "status": report["status"],
                "checks": len(report["checks"]),
                "metric_checks": len(report["metric_feasibility_checks"]),
                "integrity_checks": len(report["source_integrity_checks"]),
            }
        )
    )
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
