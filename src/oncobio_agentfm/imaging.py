"""Narrow in-memory CT intensity preprocessing built from MONAI transforms."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
from monai.transforms import Compose, EnsureChannelFirst, ScaleIntensityRange
from torch import Tensor


@dataclass(frozen=True, slots=True)
class CTIntensityWindow:
    """Hounsfield-unit interval mapped linearly into [0, 1]."""

    lower_hu: float = -1000.0
    upper_hu: float = 400.0

    def __post_init__(self) -> None:
        if not np.isfinite(self.lower_hu) or not np.isfinite(self.upper_hu):
            raise ValueError("CT window bounds must be finite")
        if self.lower_hu >= self.upper_hu:
            raise ValueError("lower_hu must be smaller than upper_hu")


def preprocess_ct_volume(
    volume_hu: np.ndarray | Tensor,
    *,
    window: CTIntensityWindow | None = None,
) -> Tensor:
    """Validating a 3D HU array, add a channel axis, clip, and scale into [0, 1]."""

    window = CTIntensityWindow() if window is None else window
    source = torch.as_tensor(volume_hu)
    if source.ndim != 3:
        raise ValueError("volume_hu must be a single 3D CT volume without a channel axis")
    if source.numel() == 0:
        raise ValueError("volume_hu must not be empty")
    if not torch.isfinite(source).all():
        raise ValueError("volume_hu must contain only finite intensities")

    transform = Compose(
        [
            EnsureChannelFirst(channel_dim="no_channel"),
            ScaleIntensityRange(
                a_min=window.lower_hu,
                a_max=window.upper_hu,
                b_min=0.0,
                b_max=1.0,
                clip=True,
                dtype=np.float32,
            ),
        ]
    )
    transformed = transform(source)
    return torch.as_tensor(transformed, dtype=torch.float32)
