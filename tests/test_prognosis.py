from dataclasses import replace
from datetime import date, datetime

import pytest

from oncobio_agentfm.clinical import ClinicalStage, StageGroup
from oncobio_agentfm.contracts import EvidenceRecord, MissingnessReason
from oncobio_agentfm.prognosis import (
    AvailabilityBasis,
    EligiblePredictors,
    PredictorFeature,
    ProtocolViolation,
    StageBasis,
    VitalStatus,
    validate_predictors,
)
from oncobio_agentfm.protocol_validation import make_synthetic_case, synthetic_probes
from oncobio_agentfm.research_protocol import FeatureKind, PredictionSetting


@pytest.fixture
def case():
    return make_synthetic_case()


def predictors(case, protocol, setting=PredictionSetting.PRETREATMENT):
    return validate_predictors(
        case.context, case.features_for(setting), setting=setting, protocol=protocol
    )


@pytest.mark.parametrize("probe", synthetic_probes(), ids=lambda probe: probe.name)
def test_adversarial_temporal_and_structural_failures(probe, protocol) -> None:
    for setting, code in (
        (PredictionSetting.PRETREATMENT, probe.expected_primary),
        (PredictionSetting.POSTOPERATIVE, probe.expected_secondary),
    ):
        if code is None:
            probe.case.follow_up.anchor(predictors(probe.case, protocol, setting))
        else:
            with pytest.raises(ProtocolViolation, match=f"^{code}$"):
                probe.case.follow_up.anchor(predictors(probe.case, protocol, setting))


def test_prediction_validation_needs_no_outcome_and_does_not_expose_later_treatment(case, protocol):
    validated = predictors(case, protocol)
    assert validated.prediction_date == date(2025, 1, 1)
    assert not hasattr(validated, "follow_up")
    assert not hasattr(validated, "surgery")
    assert not hasattr(validated, "first_treatment_on")
    assert validated.baseline is case.context.baseline


def test_eligible_container_cannot_be_constructed_without_running_the_gate(case, protocol):
    with pytest.raises(ProtocolViolation, match="predictor_gate_required"):
        EligiblePredictors(
            case.context.baseline,
            date(2025, 1, 1),
            PredictionSetting.PRETREATMENT,
            case.primary_features,
            "synthetic-shared-calendar-v1",
            protocol.sha256,
            None,
        )


def test_validated_prediction_date_cannot_be_relabelled_by_dataclass_replace(case, protocol):
    validated = predictors(case, protocol)
    with pytest.raises(TypeError, match="InitVar"):
        replace(validated, prediction_date=date(2025, 2, 10))


def test_landmark_restarts_clock_and_retains_original_baseline_information(case, protocol):
    primary = case.follow_up.anchor(predictors(case, protocol))
    landmark_inputs = predictors(case, protocol, PredictionSetting.POSTOPERATIVE)
    secondary = case.follow_up.anchor(landmark_inputs)
    assert secondary.index_date == date(2025, 2, 10)
    assert secondary.duration_days == primary.duration_days - 40
    assert landmark_inputs.baseline.index_date == primary.index_date


def test_index_must_identify_the_actual_sourced_assessment(case, protocol):
    assessment = replace(case.context.assessment, acquired_on=date(2024, 12, 31))
    changed = replace(case, context=replace(case.context, assessment=assessment))
    with pytest.raises(ProtocolViolation, match="index_not_documented_assessment"):
        predictors(changed, protocol)


def test_declared_but_absent_ct_cannot_establish_primary_imaging_eligibility(case, protocol):
    ct, pet = case.primary_features
    absent = EvidenceRecord(
        case.context.baseline.patient_id, None, None, None, MissingnessReason.UNAVAILABLE
    )
    with pytest.raises(ProtocolViolation, match="required_ct_missing"):
        predictors(
            replace(
                case,
                primary_features=(
                    replace(ct, evidence=absent, timing_source=None, availability_basis=None),
                    pet,
                ),
            ),
            protocol,
        )


@pytest.mark.parametrize("last_alive", [date(2025, 2, 9), date(2025, 2, 10)])
def test_censoring_before_or_on_landmark_does_not_establish_landmark_follow_up(
    case, protocol, last_alive
):
    follow_up = replace(
        case.follow_up, status=VitalStatus.ALIVE, death_on=None, last_known_alive_on=last_alive
    )
    assert follow_up.anchor(predictors(case, protocol)).duration_days > 0
    with pytest.raises(ProtocolViolation, match="no_positive_follow_up"):
        follow_up.anchor(predictors(case, protocol, PredictionSetting.POSTOPERATIVE))


@pytest.mark.parametrize(
    "kind",
    [
        FeatureKind.PATHOLOGICAL_STAGE,
        FeatureKind.RECURRENCE,
        FeatureKind.ADJUVANT_TREATMENT,
        FeatureKind.SURVIVAL_OUTCOME,
    ],
)
def test_semantically_forbidden_roles_rejected_even_with_early_forged_dates(case, protocol, kind):
    ct, pet = case.primary_features
    forbidden = PredictorFeature(kind, ct.evidence, "synthetic-forged-timing")
    changed = replace(case, primary_features=(ct, pet, forbidden))
    with pytest.raises(ProtocolViolation, match="forbidden_predictor_role"):
        predictors(changed, protocol)


def test_missing_modality_requires_explicit_record_and_duplicate_roles_are_rejected(case, protocol):
    ct, pet = case.primary_features
    with pytest.raises(ProtocolViolation, match="modality_missingness_undeclared"):
        predictors(replace(case, primary_features=(ct,)), protocol)
    with pytest.raises(ProtocolViolation, match="duplicate_feature_role"):
        predictors(replace(case, primary_features=(ct, pet, pet)), protocol)


def test_missing_molecular_evidence_keeps_landmark_survivor_in_coverage(case, protocol):
    ct, pet, molecular = case.landmark_features
    absent = EvidenceRecord(
        case.context.baseline.patient_id, None, None, None, MissingnessReason.UNKNOWN
    )
    missing = PredictorFeature(molecular.kind, absent, None)
    validated = predictors(
        replace(case, landmark_features=(ct, pet, missing)),
        protocol,
        PredictionSetting.POSTOPERATIVE,
    )
    assert not next(
        item for item in validated.features if item.kind is molecular.kind
    ).evidence.is_present


def test_acquisition_without_result_availability_cannot_be_imported(case):
    ct = case.primary_features[0]
    with pytest.raises(ValueError, match="availability dates"):
        replace(ct.evidence, available_on=None)
    with pytest.raises(ProtocolViolation, match="feature_timing_unverified"):
        replace(ct, timing_source=None)


@pytest.mark.parametrize("offset_date", [date(2025, 1, 1), date(2024, 12, 31)])
def test_treatment_before_or_same_calendar_day_is_excluded(case, protocol, offset_date):
    changed = replace(case, context=replace(case.context, first_treatment_on=offset_date))
    with pytest.raises(ProtocolViolation, match="treatment_before_or_on_index"):
        predictors(changed, protocol)


def test_primary_does_not_require_future_resection_but_landmark_does(case, protocol):
    changed = replace(case, context=replace(case.context, surgery=None, first_treatment_on=None))
    assert predictors(changed, protocol).prediction_date == date(2025, 1, 1)
    with pytest.raises(ProtocolViolation, match="landmark_surgery_missing"):
        predictors(changed, protocol, PredictionSetting.POSTOPERATIVE)


def test_neoadjuvant_first_treatment_disqualifies_resection_molecular_setting(case, protocol):
    changed = replace(case, context=replace(case.context, first_treatment_on=date(2025, 1, 5)))
    predictors(changed, protocol)
    with pytest.raises(ProtocolViolation, match="resection_not_verified_first_treatment"):
        predictors(changed, protocol, PredictionSetting.POSTOPERATIVE)


def test_molecular_specimen_must_match_documented_resection(case, protocol):
    ct, pet, molecular = case.landmark_features
    unrelated = replace(
        molecular, evidence=replace(molecular.evidence, acquired_on=date(2025, 1, 12))
    )
    with pytest.raises(ProtocolViolation, match="molecular_specimen_not_linked_to_resection"):
        predictors(
            replace(case, landmark_features=(ct, pet, unrelated)),
            protocol,
            PredictionSetting.POSTOPERATIVE,
        )


def test_foreign_patient_cannot_supply_ct_or_follow_up(case, protocol):
    ct, pet = case.primary_features
    foreign = replace(ct, evidence=replace(ct.evidence, patient_id="OTHER"))
    with pytest.raises(ProtocolViolation, match="patient_identity_mismatch"):
        predictors(replace(case, primary_features=(foreign, pet)), protocol)
    with pytest.raises(ProtocolViolation, match="outcome_patient_mismatch"):
        replace(case.follow_up, patient_id="OTHER").anchor(predictors(case, protocol))


def test_quarantined_historical_stage_is_retained_without_fabricating_a_new_edition(case, protocol):
    original = case.context.stage_provenance
    historical = replace(original, edition=7, basis=StageBasis.PATHOLOGICAL, adjudicated_group=None)
    missing = EvidenceRecord(
        case.context.baseline.patient_id, None, None, None, MissingnessReason.UNUSABLE
    )
    baseline = replace(case.context.baseline, stage=ClinicalStage(missing, None, None))
    context = replace(case.context, baseline=baseline, stage_provenance=historical)
    assert (
        validate_predictors(
            context,
            case.primary_features,
            setting=PredictionSetting.PRETREATMENT,
            protocol=protocol,
        ).baseline.stage.group
        is None
    )
    assert context.stage_provenance.edition == 7
    assert context.stage_provenance.raw_label == "II"
    with pytest.raises(ProtocolViolation, match="nonclinical_stage"):
        historical.as_clinical_stage(baseline.index_date)


def test_stage_provenance_must_match_encoded_stage(case, protocol):
    changed = replace(case.context.stage_provenance, adjudicated_group=StageGroup.FOUR)
    with pytest.raises(ProtocolViolation, match="stage_provenance_mismatch"):
        predictors(replace(case, context=replace(case.context, stage_provenance=changed)), protocol)


@pytest.mark.parametrize(
    ("changes", "code"),
    [
        ({"status": VitalStatus.UNKNOWN}, "vital_status_unknown"),
        ({"death_on": None}, "death_date_missing"),
        ({"status": VitalStatus.ALIVE}, "alive_with_death_date"),
        ({"last_known_alive_on": date(2028, 1, 1)}, "alive_date_after_death"),
        ({"endpoint_kind": "recurrence"}, "wrong_survival_endpoint"),
        ({"administrative_cutoff": datetime(2026, 1, 1)}, "cutoff_date_invalid"),
    ],
)
def test_endpoint_semantics_and_date_precision_are_not_guessed(case, changes, code):
    with pytest.raises(ProtocolViolation, match=code):
        replace(case.follow_up, **changes)


def test_death_after_administrative_cutoff_is_censored_at_cutoff(case, protocol):
    follow_up = replace(case.follow_up, administrative_cutoff=date(2027, 1, 1))
    outcome = follow_up.anchor(predictors(case, protocol))
    assert not outcome.event_observed
    assert outcome.end_date == date(2027, 1, 1)


def test_alive_censoring_uses_last_known_alive_not_later_registry_export(case, protocol):
    follow_up = replace(
        case.follow_up,
        status=VitalStatus.ALIVE,
        death_on=None,
        last_known_alive_on=date(2025, 3, 1),
    )
    outcome = follow_up.anchor(predictors(case, protocol))
    assert outcome.censor_date == date(2025, 3, 1)
    assert outcome.duration_days == 59


@pytest.mark.parametrize("end", [date(2024, 12, 31), date(2025, 1, 1)])
def test_endpoint_before_or_on_index_is_not_positive_follow_up(case, protocol, end):
    follow_up = replace(case.follow_up, death_on=end, last_known_alive_on=None)
    with pytest.raises(ProtocolViolation, match="no_positive_follow_up"):
        follow_up.anchor(predictors(case, protocol))


def test_prediction_day_evidence_is_allowed_at_end_of_day_but_landmark_day_death_is_excluded(
    case, protocol
):
    ct, pet = case.primary_features
    same_day_ct = replace(
        ct,
        evidence=replace(ct.evidence, acquired_on=date(2025, 1, 1), available_on=date(2025, 1, 1)),
    )
    predictors(replace(case, primary_features=(same_day_ct, pet)), protocol)
    follow_up = replace(case.follow_up, death_on=date(2025, 2, 10), last_known_alive_on=None)
    with pytest.raises(ProtocolViolation, match="no_positive_follow_up"):
        follow_up.anchor(predictors(case, protocol, PredictionSetting.POSTOPERATIVE))


def test_age_at_another_event_cannot_be_silently_treated_as_index_age(case, protocol):
    context = replace(case.context, age_reference_date=date(2024, 12, 1))
    with pytest.raises(ProtocolViolation, match="age_not_referenced_to_index"):
        predictors(replace(case, context=context), protocol)


@pytest.mark.parametrize(
    "field", ["clinical_timing_source", "age_reference_date", "index_selection_source"]
)
def test_baseline_clinical_timing_needs_source_review_and_index_reference(case, protocol, field):
    context = replace(case.context, **{field: None})
    with pytest.raises(
        ProtocolViolation,
        match=r"clinical_timing_unverified|age_reference_date_missing|index_selection_unverified",
    ):
        predictors(replace(case, context=context), protocol)


@pytest.mark.parametrize(
    "basis", [None, AvailabilityBasis.UNKNOWN, AvailabilityBasis.RESEARCH_ONLY]
)
def test_early_research_assay_date_does_not_establish_clinical_availability(case, protocol, basis):
    ct, pet, molecular = case.landmark_features
    unsupported = replace(molecular, availability_basis=basis)
    with pytest.raises(ProtocolViolation, match="feature_clinical_availability_unverified"):
        predictors(
            replace(case, landmark_features=(ct, pet, unsupported)),
            protocol,
            PredictionSetting.POSTOPERATIVE,
        )


def test_known_usable_stage_cannot_be_hidden_as_unusable(case, protocol):
    missing = EvidenceRecord(
        case.context.baseline.patient_id, None, None, None, MissingnessReason.UNUSABLE
    )
    baseline = replace(case.context.baseline, stage=ClinicalStage(missing, None, None))
    context = replace(case.context, baseline=baseline)
    with pytest.raises(ProtocolViolation, match="usable_stage_cannot_be_quarantined"):
        predictors(replace(case, context=context), protocol)


def test_missing_feature_does_not_carry_a_false_availability_declaration(case):
    pet = case.primary_features[1]
    with pytest.raises(ProtocolViolation, match="missing_feature_has_availability_claim"):
        replace(pet, timing_source="invented-result-release")


@pytest.mark.parametrize("kind", [FeatureKind.CT, FeatureKind.PET_CT])
@pytest.mark.parametrize("setting", list(PredictionSetting))
@pytest.mark.parametrize(
    "acquired,available",
    [(date(2025, 1, 2), date(2025, 1, 3)), (date(2024, 12, 28), date(2025, 1, 3))],
)
def test_pre_surgery_imaging_cannot_expand_frozen_pretreatment_information(
    case, protocol, kind, setting, acquired, available
):
    template = case.primary_features[0]
    evidence = EvidenceRecord(
        case.context.baseline.patient_id, "synthetic-late-image", acquired, available
    )
    replacement = replace(template, kind=kind, evidence=evidence)
    features = tuple(
        replacement if item.kind is kind else item for item in case.features_for(setting)
    )
    assert available < case.context.surgery.acquired_on
    with pytest.raises(ProtocolViolation, match="post_index_evidence"):
        validate_predictors(case.context, features, setting=setting, protocol=protocol)


@pytest.mark.parametrize("setting", list(PredictionSetting))
def test_old_untreated_status_does_not_prove_no_treatment_at_assessment(case, protocol, setting):
    old_status = replace(
        case.context.untreated_at_assessment,
        acquired_on=date(2024, 12, 1),
        available_on=date(2024, 12, 1),
    )
    changed = replace(case, context=replace(case.context, untreated_at_assessment=old_status))
    with pytest.raises(ProtocolViolation, match="untreated_status_not_verified_at_index"):
        predictors(changed, protocol, setting)


@pytest.mark.parametrize("setting", list(PredictionSetting))
def test_first_treatment_cannot_be_later_than_documented_cancer_surgery(case, protocol, setting):
    changed = replace(case, context=replace(case.context, first_treatment_on=date(2025, 1, 20)))
    with pytest.raises(ProtocolViolation, match="first_treatment_after_documented_surgery"):
        predictors(changed, protocol, setting)


@pytest.mark.parametrize("field", ["diagnosis", "age_evidence", "stage", "ecog"])
def test_original_clinical_contract_rejects_late_baseline_inputs_before_landmark_scoring(
    case, field
):
    baseline = case.context.baseline
    late = EvidenceRecord(
        baseline.patient_id, "synthetic-post-index-clinical", date(2025, 1, 2), date(2025, 1, 3)
    )
    if field == "stage":
        value = replace(baseline.stage, evidence=late)
    elif field == "ecog":
        value = replace(baseline.ecog, evidence=late, score=1)
    else:
        value = late
    with pytest.raises(ValueError, match=r"after the index|information leakage"):
        replace(baseline, **{field: value})


@pytest.mark.parametrize("source", [None, "", " unreviewed "])
def test_administrative_cutoff_requires_its_own_source_policy(case, source):
    with pytest.raises(ProtocolViolation, match="administrative_cutoff_source_missing"):
        replace(case.follow_up, administrative_cutoff_source=source)
