"""Clinical time-origin, stage-provenance and all-cause survival eligibility gates."""

from __future__ import annotations

from dataclasses import InitVar, dataclass
from datetime import date, timedelta
from enum import StrEnum

from oncobio_agentfm.clinical import ClinicalStage, NSCLCBaseline, StageGroup
from oncobio_agentfm.contracts import EvidenceRecord, MissingnessReason, SurvivalRecord
from oncobio_agentfm.research_protocol import FeatureKind, PredictionSetting, ResearchProtocol

_VALIDATED_PREDICTOR_TOKEN = object()


class ProtocolViolation(ValueError):
    """Stable reason codes, with no patient identifiers or source values in messages."""

    def __init__(self, code: str) -> None:
        self.code = code
        super().__init__(code)


def _source(value: str | None, code: str) -> str:
    """Require a concrete local adjudication reference; a date alone is not proof."""

    if not isinstance(value, str) or not value.strip() or value != value.strip():
        raise ProtocolViolation(code)
    return value


def _date(value: date | None, code: str) -> date:
    """Exclude timestamps, strings and booleans from a calendar-day contract."""

    if type(value) is not date:
        raise ProtocolViolation(code)
    return value


def _evidence(item: EvidenceRecord, patient_id: str, index: date) -> None:
    """Reuse the original evidence contract for identity and both timing clocks."""

    if not isinstance(item, EvidenceRecord):
        raise ProtocolViolation("evidence_type_invalid")
    if item.patient_id != patient_id:
        raise ProtocolViolation("patient_identity_mismatch")
    if item.is_present:
        _date(item.acquired_on, "acquisition_date_invalid")
        _date(item.available_on, "availability_date_invalid")
    elif not isinstance(item.missing_reason, MissingnessReason):
        raise ProtocolViolation("missingness_reason_invalid")
    try:
        item.validate_for_index(index)
    except ValueError as exc:
        raise ProtocolViolation("post_index_evidence") from exc


class StageBasis(StrEnum):
    """Pathological and unknown staging cannot substitute for clinical staging."""

    CLINICAL = "clinical"
    PATHOLOGICAL = "pathological"
    UNKNOWN = "unknown"


class AvailabilityBasis(StrEnum):
    """Clinical access and retrospective research generation are different claims."""

    CLINICAL_RECORD = "clinical_record"
    RESEARCH_ONLY = "research_only"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class StageProvenance:
    """Preserve the original label and edition; do not implement cross-edition mapping."""

    evidence: EvidenceRecord
    raw_label: str
    basis: StageBasis
    staging_system: str
    edition: int | None
    adjudicated_group: StageGroup | None
    adjudication_source: str

    def __post_init__(self) -> None:
        if not isinstance(self.evidence, EvidenceRecord) or not self.evidence.is_present:
            raise ProtocolViolation("raw_stage_source_missing")
        _source(self.raw_label, "raw_stage_label_missing")
        _source(self.staging_system, "staging_system_missing")
        _source(self.adjudication_source, "stage_adjudication_missing")
        if not isinstance(self.basis, StageBasis):
            raise ProtocolViolation("stage_basis_invalid")
        if self.edition is not None and (type(self.edition) is not int or self.edition < 1):
            raise ProtocolViolation("stage_edition_invalid")
        if self.adjudicated_group is not None and not isinstance(
            self.adjudicated_group, StageGroup
        ):
            raise ProtocolViolation("stage_group_invalid")

    def as_clinical_stage(self, index: date) -> ClinicalStage:
        """Admit only sourced clinical AJCC 8/9; retain historical values separately."""

        if self.basis is not StageBasis.CLINICAL:
            raise ProtocolViolation("nonclinical_stage")
        if self.staging_system != "AJCC" or self.edition not in (8, 9):
            raise ProtocolViolation("unsupported_or_unknown_stage_edition")
        if self.adjudicated_group is None:
            raise ProtocolViolation("stage_group_unadjudicated")
        _evidence(self.evidence, self.evidence.patient_id, index)
        return ClinicalStage(self.evidence, self.adjudicated_group, self.edition)


@dataclass(frozen=True, slots=True)
class PredictionContext:
    """Source-adjudicated clinical context; contains no survival outcome."""

    baseline: NSCLCBaseline
    assessment: EvidenceRecord
    surgical_referral: EvidenceRecord
    untreated_at_assessment: EvidenceRecord
    calendar_alignment_source: str | None
    clinical_timing_source: str | None
    age_reference_date: date | None
    index_selection_source: str | None
    first_treatment_on: date | None = None
    surgery: EvidenceRecord | None = None
    stage_provenance: StageProvenance | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.baseline, NSCLCBaseline):
            raise ProtocolViolation("baseline_type_invalid")
        if self.first_treatment_on is not None:
            _date(self.first_treatment_on, "first_treatment_date_invalid")
        if self.surgery is not None and not isinstance(self.surgery, EvidenceRecord):
            raise ProtocolViolation("surgery_type_invalid")
        if self.stage_provenance is not None and not isinstance(
            self.stage_provenance, StageProvenance
        ):
            raise ProtocolViolation("stage_provenance_type_invalid")


@dataclass(frozen=True, slots=True)
class PredictorFeature:
    """A source evidence item; molecular acquisition means specimen collection."""

    kind: FeatureKind
    evidence: EvidenceRecord
    timing_source: str | None
    availability_basis: AvailabilityBasis | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.kind, FeatureKind) or not isinstance(self.evidence, EvidenceRecord):
            raise ProtocolViolation("feature_type_invalid")
        if self.evidence.is_present:
            _source(self.timing_source, "feature_timing_unverified")
        elif self.timing_source is not None or self.availability_basis is not None:
            raise ProtocolViolation("missing_feature_has_availability_claim")
        if self.availability_basis is not None and not isinstance(
            self.availability_basis, AvailabilityBasis
        ):
            raise ProtocolViolation("availability_basis_invalid")


@dataclass(frozen=True, slots=True)
class EligiblePredictors:
    """Only eligible evidence is exposed; no outcomes or later treatment variables."""

    baseline: NSCLCBaseline
    prediction_date: date
    setting: PredictionSetting
    features: tuple[PredictorFeature, ...]
    calendar_alignment_source: str
    protocol_sha256: str
    _gate_token: InitVar[object]

    def __post_init__(self, _gate_token: object) -> None:
        if _gate_token is not _VALIDATED_PREDICTOR_TOKEN:
            raise ProtocolViolation("predictor_gate_required")

    @property
    def patient_id(self) -> str:
        """Identity is inherited from the already validated baseline."""

        return self.baseline.patient_id


def validate_predictors(
    context: PredictionContext,
    features: tuple[PredictorFeature, ...],
    *,
    setting: PredictionSetting,
    protocol: ResearchProtocol,
) -> EligiblePredictors:
    """Apply an outcome-independent gate before joining any follow-up."""

    if not isinstance(context, PredictionContext) or not isinstance(protocol, ResearchProtocol):
        raise ProtocolViolation("context_or_protocol_type_invalid")
    required = protocol.required_features(setting)
    baseline = context.baseline
    index = baseline.index_date
    patient_id = baseline.patient_id
    alignment = _source(context.calendar_alignment_source, "calendar_alignment_unverified")
    _source(context.index_selection_source, "index_selection_unverified")
    _source(context.clinical_timing_source, "clinical_timing_unverified")
    if _date(context.age_reference_date, "age_reference_date_missing") != index:
        raise ProtocolViolation("age_not_referenced_to_index")
    for item in (
        baseline.diagnosis,
        baseline.age_evidence,
        baseline.stage.evidence,
        baseline.ecog.evidence,
    ):
        _evidence(item, patient_id, index)
    for item in (context.assessment, context.surgical_referral, context.untreated_at_assessment):
        _evidence(item, patient_id, index)
        if not item.is_present:
            raise ProtocolViolation("clinical_context_missing")
    if context.assessment.acquired_on != index:
        raise ProtocolViolation("index_not_documented_assessment")
    if context.untreated_at_assessment.acquired_on != index:
        raise ProtocolViolation("untreated_status_not_verified_at_index")
    if context.first_treatment_on is not None and context.first_treatment_on <= index:
        raise ProtocolViolation("treatment_before_or_on_index")
    if context.surgery is not None:
        if context.surgery.patient_id != patient_id:
            raise ProtocolViolation("patient_identity_mismatch")
        if context.surgery.is_present:
            surgery_date = _date(context.surgery.acquired_on, "surgery_date_invalid")
            if surgery_date <= index:
                raise ProtocolViolation("surgery_before_or_on_index")
            if context.first_treatment_on is not None and context.first_treatment_on > surgery_date:
                raise ProtocolViolation("first_treatment_after_documented_surgery")
    stage = baseline.stage
    provenance = context.stage_provenance
    if stage.group is not None:
        if provenance is None:
            raise ProtocolViolation("stage_provenance_missing")
        if provenance.as_clinical_stage(index) != stage:
            raise ProtocolViolation("stage_provenance_mismatch")
    elif provenance is not None and stage.evidence.missing_reason is not MissingnessReason.UNUSABLE:
        raise ProtocolViolation("quarantined_stage_requires_unusable_reason")
    if provenance is not None and stage.group is None:
        usable_at_index = (
            provenance.basis is StageBasis.CLINICAL
            and provenance.staging_system == "AJCC"
            and provenance.edition in (8, 9)
            and provenance.adjudicated_group is not None
            and provenance.evidence.available_on is not None
            and provenance.evidence.available_on <= index
        )
        if usable_at_index:
            raise ProtocolViolation("usable_stage_cannot_be_quarantined")
    if provenance is not None and provenance.evidence.patient_id != patient_id:
        raise ProtocolViolation("patient_identity_mismatch")
    prediction_date = index
    resection_date = None
    if setting is PredictionSetting.POSTOPERATIVE:
        surgery = context.surgery
        if surgery is None or not surgery.is_present:
            raise ProtocolViolation("landmark_surgery_missing")
        resection_date = _date(surgery.acquired_on, "surgery_date_invalid")
        if context.first_treatment_on != resection_date:
            raise ProtocolViolation("resection_not_verified_first_treatment")
        prediction_date = resection_date + timedelta(days=protocol.landmark_days)
        _evidence(surgery, patient_id, prediction_date)
    if type(features) is not tuple or any(
        not isinstance(item, PredictorFeature) for item in features
    ):
        raise ProtocolViolation("features_must_be_validated_tuple")
    roles = [item.kind for item in features]
    if len(set(roles)) != len(roles):
        raise ProtocolViolation("duplicate_feature_role")
    if set(roles) - required:
        raise ProtocolViolation("forbidden_predictor_role")
    if set(roles) != required:
        raise ProtocolViolation("modality_missingness_undeclared")
    for item in features:
        _evidence(item.evidence, patient_id, prediction_date)
        if (
            item.evidence.is_present
            and item.availability_basis is not AvailabilityBasis.CLINICAL_RECORD
        ):
            raise ProtocolViolation("feature_clinical_availability_unverified")
        if item.kind is FeatureKind.CT and not item.evidence.is_present:
            raise ProtocolViolation("required_ct_missing")
        if item.kind in (FeatureKind.CT, FeatureKind.PET_CT):
            _evidence(item.evidence, patient_id, index)
        if item.kind is FeatureKind.RESECTION_MOLECULAR and item.evidence.is_present:
            if item.evidence.acquired_on != resection_date:
                raise ProtocolViolation("molecular_specimen_not_linked_to_resection")
    return EligiblePredictors(
        baseline=baseline,
        prediction_date=prediction_date,
        setting=setting,
        features=tuple(sorted(features, key=lambda item: item.kind.value)),
        calendar_alignment_source=alignment,
        protocol_sha256=protocol.sha256,
        _gate_token=_VALIDATED_PREDICTOR_TOKEN,
    )


class VitalStatus(StrEnum):
    """Unknown vital status cannot be manufactured into right censoring."""

    ALIVE = "alive"
    DEAD = "dead"
    UNKNOWN = "unknown"


@dataclass(frozen=True, slots=True)
class OverallSurvivalFollowUp:
    """Adjudicated all-cause death or last-known-alive follow-up, outside predictor data."""

    patient_id: str
    status: VitalStatus
    death_on: date | None
    last_known_alive_on: date | None
    administrative_cutoff: date
    endpoint_source: str
    calendar_alignment_source: str
    administrative_cutoff_source: str
    endpoint_kind: str = "all_cause_death"

    def __post_init__(self) -> None:
        _source(self.patient_id, "patient_id_invalid")
        _source(self.endpoint_source, "endpoint_source_missing")
        _source(self.calendar_alignment_source, "endpoint_calendar_alignment_missing")
        _source(self.administrative_cutoff_source, "administrative_cutoff_source_missing")
        if self.death_on is not None:
            _date(self.death_on, "death_date_invalid")
        if self.last_known_alive_on is not None:
            _date(self.last_known_alive_on, "last_alive_date_invalid")
        _date(self.administrative_cutoff, "cutoff_date_invalid")
        if not isinstance(self.status, VitalStatus):
            raise ProtocolViolation("vital_status_invalid")
        if self.endpoint_kind != "all_cause_death":
            raise ProtocolViolation("wrong_survival_endpoint")
        if self.status is VitalStatus.UNKNOWN:
            raise ProtocolViolation("vital_status_unknown")
        if self.status is VitalStatus.DEAD and self.death_on is None:
            raise ProtocolViolation("death_date_missing")
        if self.status is VitalStatus.ALIVE and self.death_on is not None:
            raise ProtocolViolation("alive_with_death_date")
        if self.status is VitalStatus.ALIVE and self.last_known_alive_on is None:
            raise ProtocolViolation("last_alive_date_missing")
        if self.death_on is not None and self.last_known_alive_on is not None:
            if self.last_known_alive_on > self.death_on:
                raise ProtocolViolation("alive_date_after_death")

    def anchor(self, predictors: EligiblePredictors) -> SurvivalRecord:
        """Restart survival at the selected prediction time without immortal person-time."""

        if not isinstance(predictors, EligiblePredictors):
            raise ProtocolViolation("predictors_not_validated")
        if self.patient_id != predictors.patient_id:
            raise ProtocolViolation("outcome_patient_mismatch")
        if self.calendar_alignment_source != predictors.calendar_alignment_source:
            raise ProtocolViolation("outcome_calendar_alignment_mismatch")
        index = predictors.prediction_date
        if self.death_on is not None and self.death_on <= self.administrative_cutoff:
            event = True
            endpoint = self.death_on
        elif self.death_on is not None:
            event = False
            endpoint = self.administrative_cutoff
        else:
            event = False
            if self.last_known_alive_on is None:
                raise ProtocolViolation("last_alive_date_missing")
            endpoint = min(self.last_known_alive_on, self.administrative_cutoff)
        if endpoint <= index:
            raise ProtocolViolation("no_positive_follow_up_from_prediction_time")
        return SurvivalRecord(
            patient_id=self.patient_id,
            index_date=index,
            event_date=endpoint if event else None,
            censor_date=None if event else endpoint,
        )
