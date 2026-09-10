import pytest

from oncobio_agentfm.splitting import split_patient_ids


def test_split_is_disjoint_complete_and_order_invariant() -> None:
    patient_ids = [f"P{i:04d}" for i in range(200)]
    first = split_patient_ids(patient_ids, seed="fixed-study-seed")
    second = split_patient_ids(reversed(patient_ids), seed="fixed-study-seed")

    assert first == second
    assert first.all_patient_ids == frozenset(patient_ids)
    assert set(first.train).isdisjoint(first.validation)
    assert set(first.train).isdisjoint(first.test)
    assert set(first.validation).isdisjoint(first.test)


def test_split_rejects_duplicate_patient_ids() -> None:
    with pytest.raises(ValueError, match="unique"):
        split_patient_ids(["P001", "P001", "P002"])


def test_split_seed_changes_at_least_one_assignment() -> None:
    patient_ids = [f"P{i:04d}" for i in range(200)]
    first = split_patient_ids(patient_ids, seed="seed-a")
    second = split_patient_ids(patient_ids, seed="seed-b")

    assert first != second


def test_split_rejects_invalid_fraction_sum() -> None:
    with pytest.raises(ValueError, match="sum to 1"):
        split_patient_ids(
            ["P001", "P002", "P003"],
            train_fraction=0.60,
            validation_fraction=0.20,
            test_fraction=0.10,
        )


def test_split_rejects_noncanonical_patient_identifier() -> None:
    with pytest.raises(ValueError, match="surrounding whitespace"):
        split_patient_ids(["P001", " P002", "P003"])
