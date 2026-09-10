import math

import pytest
import torch

from oncobio_agentfm.survival import CoxPHLinear, efron_negative_partial_log_likelihood


def test_cox_linear_returns_one_log_risk_per_patient() -> None:
    model = CoxPHLinear(in_features=3)
    features = torch.ones((4, 3), dtype=torch.float64)
    model = model.to(dtype=torch.float64)

    assert model(features).shape == (4,)


def test_efron_loss_matches_two_patient_no_tie_case() -> None:
    log_risk = torch.tensor([math.log(2.0), 0.0], dtype=torch.float64)
    durations = torch.tensor([1.0, 2.0], dtype=torch.float64)
    events = torch.tensor([1, 1], dtype=torch.int64)

    loss = efron_negative_partial_log_likelihood(log_risk, durations, events)
    expected = math.log(1.5) / 2.0

    assert loss.item() == pytest.approx(expected, rel=1e-12, abs=1e-12)


def test_efron_loss_handles_tied_events() -> None:
    log_risk = torch.zeros(3, dtype=torch.float64)
    durations = torch.tensor([1.0, 1.0, 2.0], dtype=torch.float64)
    events = torch.tensor([1, 1, 0], dtype=torch.int64)

    loss = efron_negative_partial_log_likelihood(log_risk, durations, events)

    assert loss.item() == pytest.approx(math.log(6.0) / 2.0, rel=1e-12, abs=1e-12)


def test_efron_loss_is_row_order_invariant_and_differentiable() -> None:
    log_risk = torch.tensor([0.4, -0.2, 0.1, 0.8], dtype=torch.float64, requires_grad=True)
    durations = torch.tensor([2.0, 1.0, 2.0, 3.0], dtype=torch.float64)
    events = torch.tensor([1, 1, 1, 0], dtype=torch.int64)
    order = torch.tensor([3, 1, 0, 2])

    first = efron_negative_partial_log_likelihood(log_risk, durations, events)
    second = efron_negative_partial_log_likelihood(log_risk[order], durations[order], events[order])
    first.backward()

    assert second.item() == pytest.approx(first.item(), rel=1e-12, abs=1e-12)
    assert log_risk.grad is not None
    assert torch.isfinite(log_risk.grad).all()


def test_efron_loss_rejects_zero_event_cohort() -> None:
    log_risk = torch.zeros(3, dtype=torch.float64)
    durations = torch.tensor([1.0, 2.0, 3.0], dtype=torch.float64)
    events = torch.zeros(3, dtype=torch.int64)

    with pytest.raises(ValueError, match="observed event"):
        efron_negative_partial_log_likelihood(log_risk, durations, events)


def test_efron_loss_rejects_nonfinite_inputs() -> None:
    log_risk = torch.tensor([0.0, float("nan")], dtype=torch.float64)
    durations = torch.tensor([1.0, 2.0], dtype=torch.float64)
    events = torch.tensor([1, 0], dtype=torch.int64)

    with pytest.raises(ValueError, match="finite"):
        efron_negative_partial_log_likelihood(log_risk, durations, events)
