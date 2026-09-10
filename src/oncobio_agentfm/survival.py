"""Minimal PyTorch Cox proportional-hazards baseline with Efron tie handling."""

from __future__ import annotations

from typing import Literal

import torch
from torch import Tensor, nn


class CoxPHLinear(nn.Module):
    """Linear log-risk model without an unidentifiable Cox intercept."""

    def __init__(self, in_features: int) -> None:
        super().__init__()
        if in_features <= 0:
            raise ValueError("in_features must be positive")
        self.linear = nn.Linear(in_features, 1, bias=False)

    def forward(self, features: Tensor) -> Tensor:
        """Returning one log-risk score per patient."""

        if features.ndim != 2:
            raise ValueError("features must have shape [patients, features]")
        if features.shape[1] != self.linear.in_features:
            raise ValueError("feature dimension does not match model input size")
        return self.linear(features).squeeze(-1)


def efron_negative_partial_log_likelihood(
    log_risk: Tensor,
    durations: Tensor,
    events: Tensor,
    *,
    reduction: Literal["mean", "sum"] = "mean",
) -> Tensor:
    """Compute the Cox negative partial log-likelihood using Efron's tie correction."""

    log_risk = log_risk.squeeze(-1) if log_risk.ndim == 2 and log_risk.shape[1] == 1 else log_risk
    if log_risk.ndim != 1 or durations.ndim != 1 or events.ndim != 1:
        raise ValueError("log_risk, durations, and events must be one-dimensional")
    if not (log_risk.shape == durations.shape == events.shape):
        raise ValueError("log_risk, durations, and events must have identical shapes")
    if not (log_risk.device == durations.device == events.device):
        raise ValueError("log_risk, durations, and events must be on the same device")
    if log_risk.numel() == 0:
        raise ValueError("survival tensors must not be empty")
    if not log_risk.dtype.is_floating_point:
        raise TypeError("log_risk must use a floating-point dtype")
    if not torch.isfinite(log_risk).all() or not torch.isfinite(durations).all():
        raise ValueError("log_risk and durations must be finite")
    if torch.any(durations < 0):
        raise ValueError("durations must be non-negative")
    if not torch.all((events == 0) | (events == 1)):
        raise ValueError("events must contain only 0/1 values")
    if reduction not in {"mean", "sum"}:
        raise ValueError("reduction must be 'mean' or 'sum'")

    event_mask = events.to(dtype=torch.bool)
    event_count = int(event_mask.sum().item())
    if event_count == 0:
        raise ValueError("at least one observed event is required")

    event_times = torch.unique(durations[event_mask], sorted=True)
    partial_log_likelihood = log_risk.new_zeros(())

    for event_time in event_times:
        tied_mask = event_mask & (durations == event_time)
        risk_mask = durations >= event_time
        tied_scores = log_risk[tied_mask]
        risk_scores = log_risk[risk_mask]
        tied_count = tied_scores.numel()

        stabilizer = risk_scores.max()
        risk_sum = torch.exp(risk_scores - stabilizer).sum()
        tied_sum = torch.exp(tied_scores - stabilizer).sum()
        fractions = torch.arange(
            tied_count,
            dtype=log_risk.dtype,
            device=log_risk.device,
        ) / tied_count
        denominators = risk_sum - fractions * tied_sum
        tiny = torch.finfo(log_risk.dtype).tiny
        log_denominators = stabilizer + torch.log(torch.clamp(denominators, min=tiny))

        partial_log_likelihood = partial_log_likelihood + tied_scores.sum() - log_denominators.sum()

    negative_log_likelihood = -partial_log_likelihood
    if reduction == "mean":
        return negative_log_likelihood / event_count
    return negative_log_likelihood
