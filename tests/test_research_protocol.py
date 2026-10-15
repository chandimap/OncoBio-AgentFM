import json
from dataclasses import replace
from pathlib import Path

import pytest

from oncobio_agentfm.research_protocol import (
    FeatureKind,
    PredictionSetting,
    load_research_protocol,
)

SPEC_PATH = Path(__file__).resolve().parents[1] / "configs/research/nsclc_os_v1.json"


def test_reviewed_protocol_has_distinct_settings_and_explicit_modality_absence(protocol) -> None:
    assert protocol.required_features(PredictionSetting.PRETREATMENT) == {
        FeatureKind.CT,
        FeatureKind.PET_CT,
    }
    assert protocol.required_features(PredictionSetting.POSTOPERATIVE) == {
        FeatureKind.CT,
        FeatureKind.PET_CT,
        FeatureKind.RESECTION_MOLECULAR,
    }
    assert protocol.landmark_days == 30
    assert protocol.primary_horizon_days == 730


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("endpoint", "recurrence"),
        ("landmark_days", 60),
        ("primary_horizon_days", True),
        ("secondary_horizons_days", [365.0, 1095]),
    ],
)
def test_unsupported_estimand_changes_fail_closed(tmp_path, key, value) -> None:
    spec = json.loads(SPEC_PATH.read_text())
    spec[key] = value
    path = tmp_path / "changed.json"
    path.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match="unsupported protocol choice"):
        load_research_protocol(path)


def test_protocol_rejects_unknown_and_duplicate_keys(tmp_path) -> None:
    path = tmp_path / "changed.json"
    spec = json.loads(SPEC_PATH.read_text())
    spec["ignore_late_results"] = True
    path.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match="schema"):
        load_research_protocol(path)
    path.write_text('{"schema_version":"1.0","schema_version":"2.0"}')
    with pytest.raises(ValueError, match="duplicate"):
        load_research_protocol(path)


def test_protocol_hash_tracks_exact_reviewed_bytes(tmp_path, protocol) -> None:
    path = tmp_path / "reformatted.json"
    path.write_bytes(SPEC_PATH.read_bytes() + b"\n")
    updated = load_research_protocol(path)
    assert updated.sha256 != protocol.sha256
    assert updated.primary_horizon_days == protocol.primary_horizon_days


def test_direct_protocol_construction_cannot_change_landmark_or_remove_fingerprint(
    protocol,
) -> None:
    with pytest.raises(ValueError, match="landmark"):
        replace(protocol, landmark_days=60)
    with pytest.raises(ValueError, match="SHA-256"):
        replace(protocol, sha256="")


def test_setting_strings_are_not_silently_interpreted(protocol) -> None:
    with pytest.raises(TypeError, match="PredictionSetting"):
        protocol.required_features("pretreatment")


def test_nested_feasibility_policy_cannot_replace_integer_grid_step_with_boolean(tmp_path):
    spec = json.loads(SPEC_PATH.read_text())
    spec["feasibility_policy"]["ibs_grid_step_days"] = True
    path = tmp_path / "changed.json"
    path.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match="feasibility_policy"):
        load_research_protocol(path)


def test_nested_policy_and_revision_changes_fail_closed(tmp_path):
    for key, value in (("revision", "unreviewed"), ("feasibility_policy", {})):
        spec = json.loads(SPEC_PATH.read_text())
        spec[key] = value
        path = tmp_path / "changed.json"
        path.write_text(json.dumps(spec))
        with pytest.raises(ValueError, match="unsupported protocol choice"):
            load_research_protocol(path)


def test_administrative_and_frozen_information_policy_cannot_be_silently_relaxed(tmp_path):
    spec = json.loads(SPEC_PATH.read_text())
    spec["chronology_policy"]["administrative_cutoff"] = "choose_per_patient_after_results"
    path = tmp_path / "changed.json"
    path.write_text(json.dumps(spec))
    with pytest.raises(ValueError, match="chronology_policy"):
        load_research_protocol(path)
