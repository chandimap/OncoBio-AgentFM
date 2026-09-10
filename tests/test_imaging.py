import numpy as np
import pytest
import torch

from oncobio_agentfm.imaging import CTIntensityWindow, preprocess_ct_volume


def test_ct_preprocessing_adds_channel_clips_and_scales() -> None:
    volume = np.array(
        [
            [[-1200.0, -1000.0], [-300.0, 400.0]],
            [[900.0, -300.0], [-1000.0, 400.0]],
        ],
        dtype=np.float32,
    )

    transformed = preprocess_ct_volume(volume)

    assert transformed.shape == (1, 2, 2, 2)
    assert transformed.dtype == torch.float32
    assert transformed.min().item() == pytest.approx(0.0)
    assert transformed.max().item() == pytest.approx(1.0)
    assert transformed[0, 0, 1, 0].item() == pytest.approx(0.5)


def test_ct_preprocessing_rejects_non_3d_input() -> None:
    with pytest.raises(ValueError, match="3D CT volume"):
        preprocess_ct_volume(np.zeros((1, 8, 8, 8), dtype=np.float32))


def test_ct_window_rejects_invalid_bounds() -> None:
    with pytest.raises(ValueError, match="smaller"):
        CTIntensityWindow(lower_hu=100.0, upper_hu=100.0)


def test_ct_preprocessing_respects_custom_window() -> None:
    volume = np.array([[[0.0, 100.0]]], dtype=np.float32)
    window = CTIntensityWindow(lower_hu=0.0, upper_hu=200.0)

    transformed = preprocess_ct_volume(volume, window=window)

    assert transformed[0, 0, 0, 0].item() == pytest.approx(0.0)
    assert transformed[0, 0, 0, 1].item() == pytest.approx(0.5)


def test_ct_preprocessing_rejects_nonfinite_values() -> None:
    volume = np.array([[[0.0, np.nan]]], dtype=np.float32)

    with pytest.raises(ValueError, match="finite intensities"):
        preprocess_ct_volume(volume)
