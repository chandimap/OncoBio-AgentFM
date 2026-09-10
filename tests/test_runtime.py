from oncobio_agentfm.runtime import EXPECTED_PACKAGES, EXPECTED_PYTHON, assert_expected_runtime


def test_canonical_runtime_versions() -> None:
    observed = assert_expected_runtime()

    assert observed["python"] == EXPECTED_PYTHON
    for package, expected_version in EXPECTED_PACKAGES.items():
        assert observed[package] == expected_version
