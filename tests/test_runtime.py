import pytest

import oncobio_agentfm.runtime as runtime
from oncobio_agentfm.runtime import EXPECTED_PACKAGES, EXPECTED_PYTHON, assert_expected_runtime


def test_canonical_runtime_versions() -> None:
    observed = assert_expected_runtime()

    assert observed["python"] == EXPECTED_PYTHON
    for package, expected_version in EXPECTED_PACKAGES.items():
        assert observed[package] in (
            (expected_version, f"{expected_version}+cpu")
            if package == "torch"
            else (expected_version,)
        )


def test_cpu_build_keeps_exact_public_torch_version(monkeypatch) -> None:
    expected = {"python": EXPECTED_PYTHON, **EXPECTED_PACKAGES}
    monkeypatch.setattr(runtime, "observed_runtime", lambda: {**expected, "torch": "2.14.0+cpu"})
    assert runtime.assert_expected_runtime()["torch"] == "2.14.0+cpu"

    monkeypatch.setattr(runtime, "observed_runtime", lambda: {**expected, "torch": "2.13.0+cpu"})
    with pytest.raises(RuntimeError, match="torch"):
        runtime.assert_expected_runtime()
