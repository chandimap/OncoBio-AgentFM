from __future__ import annotations

import platform
from importlib.metadata import PackageNotFoundError, version

EXPECTED_PYTHON = "3.13.15"
EXPECTED_PACKAGES = {
    "torch": "2.14.0",
    "monai": "1.6.0",
    "numpy": "2.5.3",
}


def observed_runtime() -> dict[str, str]:
    """Collecting the interpreter and direct scientific dependency versions."""

    observed = {"python": platform.python_version()}
    for package in EXPECTED_PACKAGES:
        try:
            observed[package] = version(package)
        except PackageNotFoundError:
            observed[package] = "NOT_INSTALLED"
    return observed


def assert_expected_runtime() -> dict[str, str]:
    """Raising when the active environment differs from the canonical direct-version contract."""

    observed = observed_runtime()
    expected = {"python": EXPECTED_PYTHON, **EXPECTED_PACKAGES}
    mismatches = {
        name: (expected[name], observed[name])
        for name in expected
        if observed[name] != expected[name]
    }
    if mismatches:
        details = ", ".join(
            f"{name}: expected {wanted}, observed {actual}"
            for name, (wanted, actual) in mismatches.items()
        )
        raise RuntimeError(f"runtime verification failed ({details})")
    return observed


def main() -> None:
    """Verifying and printing the canonical runtime."""

    observed = assert_expected_runtime()
    print("OncoBio-AgentFM runtime verified:")
    for name, observed_version in observed.items():
        print(f"  {name}={observed_version}")


if __name__ == "__main__":
    main()
