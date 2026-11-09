"""Bind adjudicated cases to verified bytes before reusing the existing protocol gates."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from datetime import date
from pathlib import Path
from types import MappingProxyType

from oncobio_agentfm.cohort_sources import (
    SourceInventory,
    implementation_fingerprint,
    require,
    unique_object,
)
from oncobio_agentfm.contracts import EvidenceRecord
from oncobio_agentfm.prognosis import ProtocolViolation, validate_predictors
from oncobio_agentfm.protocol_audit import AdjudicatedCase, audit_eligibility
from oncobio_agentfm.research_protocol import FeatureKind, PredictionSetting, ResearchProtocol


def case_payload(case: AdjudicatedCase) -> dict[str, object]:
    """Canonical local review representation; never publish this patient-level payload."""

    def calendar(value: object) -> str:
        require(type(value) is date, "review_payload_type_invalid")
        return value.isoformat()

    return json.loads(json.dumps(asdict(case), default=calendar, allow_nan=False))


def case_references(case: AdjudicatedCase) -> frozenset[str]:
    """Collect all concrete evidence and review references, including endpoint policy."""

    context = case.context
    baseline = context.baseline
    evidence = [
        baseline.diagnosis,
        baseline.age_evidence,
        baseline.stage.evidence,
        baseline.ecog.evidence,
        context.assessment,
        context.surgical_referral,
        context.untreated_at_assessment,
    ]
    if context.surgery is not None:
        evidence.append(context.surgery)
    if context.stage_provenance is not None:
        evidence.append(context.stage_provenance.evidence)
    features = (*case.primary_features, *case.landmark_features)
    evidence.extend(feature.evidence for feature in features)
    sources = [item.source_id for item in evidence]
    sources.extend(feature.timing_source for feature in features)
    sources.extend(
        (
            context.calendar_alignment_source,
            context.clinical_timing_source,
            context.index_selection_source,
            case.follow_up.endpoint_source,
            case.follow_up.calendar_alignment_source,
            case.follow_up.administrative_cutoff_source,
        )
    )
    if context.stage_provenance is not None:
        sources.append(context.stage_provenance.adjudication_source)
    return frozenset(source for source in sources if source is not None)


def verify_documents(paths: dict[str, Path]) -> Mapping[str, bytes]:
    """Verify external adjudication-document fingerprints without exposing their paths."""

    verified = {}
    for expected, path in paths.items():
        raw = path.read_bytes()
        require(
            hashlib.sha256(raw).hexdigest() == expected,
            "adjudication_document_hash_mismatch",
        )
        verified[expected] = raw
    return MappingProxyType(verified)


def clinical_references(case: AdjudicatedCase) -> frozenset[str]:
    """Clinical facts and access assertions cannot be supported by an image/GEO row."""

    context = case.context
    baseline = context.baseline
    evidence = [
        baseline.diagnosis,
        baseline.age_evidence,
        baseline.stage.evidence,
        baseline.ecog.evidence,
        context.assessment,
        context.surgical_referral,
        context.untreated_at_assessment,
    ]
    if context.surgery is not None:
        evidence.append(context.surgery)
    if context.stage_provenance is not None:
        evidence.append(context.stage_provenance.evidence)
    refs = [item.source_id for item in evidence]
    refs.extend(
        feature.timing_source for feature in (*case.primary_features, *case.landmark_features)
    )
    refs.extend(
        (
            context.calendar_alignment_source,
            context.clinical_timing_source,
            context.index_selection_source,
            case.follow_up.endpoint_source,
            case.follow_up.calendar_alignment_source,
            case.follow_up.administrative_cutoff_source,
        )
    )
    if context.stage_provenance is not None:
        refs.append(context.stage_provenance.adjudication_source)
    return frozenset(ref for ref in refs if ref is not None)


def _resolve_pointer(document: object, pointer: str) -> object:
    """Resolve an RFC 6901 pointer; a merely nonempty location is not provenance."""

    require(pointer.startswith("/") and len(pointer) <= 4096, "clinical_locator_unresolved")
    current = document
    for part in pointer[1:].split("/"):
        require("~" not in part.replace("~0", "").replace("~1", ""), "clinical_locator_unresolved")
        token = part.replace("~1", "/").replace("~0", "~")
        require(type(current) is dict and token in current, "clinical_locator_unresolved")
        current = current[token]
    return current


def _canonical(value: object) -> str:
    try:
        return json.dumps(value, sort_keys=True, allow_nan=False)
    except (TypeError, ValueError):
        raise ProtocolViolation("case_review_value_invalid") from None


def _paired_evidence(value: object) -> EvidenceRecord:
    require(
        type(value) is dict
        and set(value)
        == {"patient_id", "source_id", "acquired_on", "available_on", "missing_reason"},
        "paired_ct_evidence_invalid",
    )
    require(value["missing_reason"] is None, "paired_ct_evidence_invalid")
    require(
        type(value["patient_id"]) is str
        and type(value["source_id"]) is str
        and bool(value["patient_id"].strip())
        and bool(value["source_id"].strip()),
        "paired_ct_evidence_invalid",
    )
    try:
        acquired = date.fromisoformat(value["acquired_on"])
        available = date.fromisoformat(value["available_on"])
        require(
            acquired.isoformat() == value["acquired_on"]
            and available.isoformat() == value["available_on"],
            "paired_ct_evidence_invalid",
        )
        return EvidenceRecord(value["patient_id"], value["source_id"], acquired, available)
    except (TypeError, ValueError):
        raise ProtocolViolation("paired_ct_evidence_invalid") from None


def verify_case_review(
    case: AdjudicatedCase,
    review_path: Path,
    *,
    expected_sha256: str,
    inventory: SourceInventory,
    verified_documents: Mapping[str, bytes],
) -> None:
    """Require exact reviewed case values and a source locator for every assertion."""

    raw = review_path.read_bytes()
    require(hashlib.sha256(raw).hexdigest() == expected_sha256, "case_review_hash_mismatch")
    review = json.loads(raw, object_pairs_hook=unique_object)
    require(
        type(review) is dict
        and set(review)
        == {
            "schema_version",
            "case",
            "references",
            "image_date_mappings",
            "pet_ct_pairings",
            "molecular_provenance",
        },
        "case_review_schema_invalid",
    )
    require(review["schema_version"] == "1.1", "case_review_version_invalid")
    expected = json.dumps(case_payload(case), sort_keys=True, allow_nan=False)
    require(
        json.dumps(review["case"], sort_keys=True, allow_nan=False) == expected,
        "case_differs_from_reviewed_source",
    )
    references = review["references"]
    require(
        isinstance(verified_documents, Mapping)
        and all(type(raw) is bytes for raw in verified_documents.values()),
        "case_document_bytes_required",
    )
    require(
        all(
            hashlib.sha256(raw).hexdigest() == fingerprint
            for fingerprint, raw in verified_documents.items()
        ),
        "adjudication_document_hash_mismatch",
    )
    features = (*case.primary_features, *case.landmark_features)
    paired = review["pet_ct_pairings"]
    pet_ids = {
        feature.evidence.source_id
        for feature in features
        if feature.kind is FeatureKind.PET_CT and feature.evidence.is_present
    }
    require(type(paired) is dict and set(paired) == pet_ids, "pet_ct_pairing_undeclared")
    pair_features = []
    pair_timing = {}
    extra_refs = set()
    for pair in paired.values():
        require(
            type(pair) is dict
            and set(pair)
            == {"ct_evidence", "timing_source", "pairing_source", "availability_basis"},
            "pet_ct_pairing_invalid",
        )
        item = _paired_evidence(pair["ct_evidence"])
        require(
            item.patient_id == case.context.baseline.patient_id and item.is_present,
            "paired_ct_patient_mismatch",
        )
        require(
            pair["availability_basis"] == "clinical_record", "paired_ct_availability_unverified"
        )
        try:
            item.validate_for_index(case.context.baseline.index_date)
        except ValueError:
            raise ProtocolViolation("paired_ct_not_available_at_index") from None
        for key in ("timing_source", "pairing_source"):
            require(
                type(pair[key]) is str
                and bool(pair[key].strip())
                and pair[key] == pair[key].strip(),
                "pet_ct_pairing_reference_missing",
            )
            extra_refs.add(pair[key])
        pair_features.append(item)
        pair_timing[item.source_id] = pair["timing_source"]
    required_refs = case_references(case) | extra_refs | {item.source_id for item in pair_features}
    require(
        type(references) is dict and set(references) == required_refs, "case_lineage_incomplete"
    )
    fingerprints = set(verified_documents) | {receipt.sha256 for receipt in inventory.receipts}
    for link in references.values():
        require(
            type(link) is dict and set(link) == {"source_sha256", "locator"},
            "case_lineage_shape_invalid",
        )
        require(
            type(link["source_sha256"]) is str and link["source_sha256"] in fingerprints,
            "case_lineage_source_unverified",
        )
        require(
            type(link["locator"]) is str
            and bool(link["locator"].strip())
            and link["locator"] == link["locator"].strip(),
            "case_lineage_locator_missing",
        )
    clinical_refs = clinical_references(case) | extra_refs
    feature_ids = {
        feature.evidence.source_id for feature in features if feature.evidence.is_present
    }
    feature_ids.update(item.source_id for item in pair_features)
    require(not clinical_refs & feature_ids, "source_reference_role_collision")
    bundle = {
        "case": case_payload(case),
        "pet_ct_pairings": paired,
        "molecular_provenance": review["molecular_provenance"],
    }
    for ref in clinical_refs:
        link = references[ref]
        require(link["source_sha256"] in verified_documents, "clinical_assertion_source_required")
        try:
            document = json.loads(
                verified_documents[link["source_sha256"]], object_pairs_hook=unique_object
            )
        except (UnicodeError, json.JSONDecodeError):
            raise ProtocolViolation("clinical_assertion_document_invalid") from None
        require(
            type(document) is dict
            and document.get("schema_version") == "1.0"
            and document.get("evidence_type") == "source_adjudication_assertions",
            "clinical_assertion_document_invalid",
        )
        require(
            _canonical(_resolve_pointer(document, link["locator"])) == _canonical(bundle),
            "clinical_assertion_case_mismatch",
        )
    patient = case.context.baseline.patient_id
    require(patient in inventory.clinical, "case_not_in_clinical_source")
    series = {item.series_uid: item for item in inventory.series}
    selected_ids = {
        feature.evidence.source_id
        for feature in (*case.primary_features, *case.landmark_features)
        if feature.kind in (FeatureKind.CT, FeatureKind.PET_CT) and feature.evidence.is_present
    }
    selected_ids.update(item.source_id for item in pair_features)
    mappings = review["image_date_mappings"]
    require(
        type(mappings) is dict and set(mappings) == selected_ids,
        "image_calendar_mapping_incomplete",
    )
    selected_evidence = [
        (feature.kind, feature.evidence, feature.timing_source) for feature in features
    ]
    selected_evidence.extend(
        (FeatureKind.CT, item, pair_timing[item.source_id]) for item in pair_features
    )
    unique_evidence = {}
    for kind, item, timing_source in selected_evidence:
        if not item.is_present:
            continue
        require(
            item.source_id not in unique_evidence or unique_evidence[item.source_id] == item,
            "selected_source_evidence_conflict",
        )
        unique_evidence[item.source_id] = item
        if kind in (FeatureKind.CT, FeatureKind.PET_CT):
            selected = series.get(item.source_id)
            require(
                selected is not None and selected.patient_id == patient,
                "selected_image_patient_or_source_mismatch",
            )
            expected_modality = "CT" if kind is FeatureKind.CT else "PT"
            require(selected.modality == expected_modality, "selected_image_modality_mismatch")
            require(
                references[item.source_id]
                == {
                    "source_sha256": next(
                        r.sha256 for r in inventory.receipts if r.kind == "imaging"
                    ),
                    "locator": f"/{selected.source_row}",
                },
                "selected_image_locator_mismatch",
            )
            require(
                item.acquired_on is not None,
                "selected_image_acquisition_missing",
            )
            require(
                mappings[item.source_id]
                == {
                    "source_axis_date": selected.series_date,
                    "source_date_kind": "SeriesDate_0008_0021",
                    "acquired_on": item.acquired_on.isoformat(),
                    "acquisition_reference": timing_source,
                    "alignment_reference": case.context.calendar_alignment_source,
                },
                "selected_image_calendar_mapping_mismatch",
            )
        if kind is FeatureKind.RESECTION_MOLECULAR:
            require(
                inventory.rna_samples.get(item.source_id) == patient
                and inventory.expression_observed_per_patient.get(patient, 0) > 0,
                "selected_rna_patient_or_source_mismatch",
            )
            require(
                references[item.source_id]
                == {
                    "source_sha256": next(
                        r.sha256 for r in inventory.receipts if r.kind == "rna_metadata"
                    ),
                    "locator": f"sample/{item.source_id}",
                },
                "selected_rna_locator_mismatch",
            )
            require(
                sum(owner == patient for owner in inventory.rna_samples.values()) == 1,
                "rna_aliquot_column_link_ambiguous",
            )
    molecular = review["molecular_provenance"]
    molecular_features = [
        feature
        for feature in features
        if feature.kind is FeatureKind.RESECTION_MOLECULAR and feature.evidence.is_present
    ]
    require(
        type(molecular) is dict
        and set(molecular) == {feature.evidence.source_id for feature in molecular_features},
        "molecular_provenance_incomplete",
    )
    for feature in molecular_features:
        item = feature.evidence
        surgery = case.context.surgery
        expected_provenance = {
            "expression_source_sha256": next(
                r.sha256 for r in inventory.receipts if r.kind == "rna_expression"
            ),
            "expression_locator": f"tsv/column/{inventory.expression_source_columns[patient]}",
            "resection_source_id": surgery.source_id if surgery is not None else None,
            "specimen_collected_on": item.acquired_on.isoformat(),
            "clinical_available_on": item.available_on.isoformat(),
            "availability_basis": feature.availability_basis.value
            if feature.availability_basis
            else None,
            "timing_source": feature.timing_source,
        }
        require(
            _canonical(molecular[item.source_id]) == _canonical(expected_provenance),
            "molecular_provenance_mismatch",
        )


@dataclass(frozen=True, slots=True)
class SourceExclusion:
    """A documented first-failure exclusion before validated case construction."""

    patient_id: str
    reason: str
    document_sha256: str
    locator: str

    def __post_init__(self) -> None:
        require(
            self.reason
            in {
                "source_adjudication_unresolved",
                "invalid_identity",
                "contradictory_survival",
                "mandatory_context_unavailable",
                "calendar_alignment_unverified",
            },
            "source_exclusion_reason_invalid",
        )
        require(
            type(self.locator) is str and bool(self.locator.strip()),
            "source_exclusion_locator_missing",
        )


def audit_linked_cases(
    cases: tuple[AdjudicatedCase, ...],
    *,
    inventory: SourceInventory,
    reviews: dict[str, tuple[Path, str]],
    documents: dict[str, Path],
    protocol: ResearchProtocol,
    exclusions: tuple[SourceExclusion, ...] = (),
) -> dict[str, object]:
    """Reject incomplete lineage before reporting any clinically eligible population."""

    require(
        type(cases) is tuple and all(isinstance(case, AdjudicatedCase) for case in cases),
        "linked_cases_invalid",
    )
    patients = [case.context.baseline.patient_id for case in cases]
    require(len(set(patients)) == len(patients), "linked_duplicate_patient")
    require(set(patients) == set(reviews), "case_review_coverage_mismatch")
    require(
        type(exclusions) is tuple and all(isinstance(item, SourceExclusion) for item in exclusions),
        "source_exclusions_invalid",
    )
    excluded = [item.patient_id for item in exclusions]
    require(
        len(set(excluded)) == len(excluded) and not set(excluded) & set(patients),
        "source_exclusion_identity_overlap",
    )
    require(
        frozenset((*patients, *excluded)) == inventory.patient_ids,
        "case_inventory_coverage_incomplete",
    )
    verified = verify_documents(documents)
    for item in exclusions:
        require(item.document_sha256 in verified, "source_exclusion_document_unverified")
        document = json.loads(verified[item.document_sha256], object_pairs_hook=unique_object)
        require(
            type(document) is dict
            and document.get(item.locator)
            == {"patient_id": item.patient_id, "reason": item.reason},
            "source_exclusion_review_mismatch",
        )
    for case in cases:
        path, fingerprint = reviews[case.context.baseline.patient_id]
        verify_case_review(
            case,
            path,
            expected_sha256=fingerprint,
            inventory=inventory,
            verified_documents=verified,
        )
    reports = {
        setting.value: (
            audit_eligibility(cases, setting=setting, protocol=protocol)
            if cases
            else {
                "status": "NO_CASES_AFTER_DOCUMENTED_SOURCE_EXCLUSIONS",
                "input_patients": 0,
                "predictor_eligible": 0,
                "endpoint_eligible": 0,
                "partition_counts": {},
            }
        )
        for setting in PredictionSetting
    }
    # An additional complete-modality population retains the existing landmark origin.
    complete = []
    for case in cases:
        try:
            predictors = validate_predictors(
                case.context,
                case.features_for(PredictionSetting.POSTOPERATIVE),
                setting=PredictionSetting.POSTOPERATIVE,
                protocol=protocol,
            )
            case.follow_up.anchor(predictors)
        except ProtocolViolation:
            continue
        if all(feature.evidence.is_present for feature in predictors.features):
            complete.append(case)
    complete_report = (
        audit_eligibility(
            tuple(complete), setting=PredictionSetting.POSTOPERATIVE, protocol=protocol
        )
        if complete
        else None
    )
    return {
        "schema_version": "1.0",
        "evidence_type": "linked_local_adjudication_audit",
        "implementation_sha256": implementation_fingerprint(),
        "input_patients": len(inventory.patient_ids),
        "lineage_verified_patients": len(cases),
        "source_exclusions": dict(sorted(Counter(item.reason for item in exclusions).items())),
        "settings": reports,
        "validated_event_budget": {
            setting: {
                "eligible_patients": report["endpoint_eligible"],
                "deaths": sum(part["deaths"] for part in report["partition_counts"].values()),
                "censored": sum(part["censored"] for part in report["partition_counts"].values()),
            }
            for setting, report in reports.items()
        },
        "all_four_landmark": complete_report,
        "all_four_landmark_patients": len(complete),
        "all_four_landmark_deaths": sum(
            part["deaths"] for part in complete_report["partition_counts"].values()
        )
        if complete_report
        else 0,
        "lineage_scope": "byte_identity_and_reviewed_assertions_not_source_authentication",
        "clinical_validity": "NOT_ESTABLISHED_BY_SOFTWARE_CHECKS",
    }
