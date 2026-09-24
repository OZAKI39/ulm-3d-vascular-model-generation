"""Dependency-free robust feature scaling with a serializable state."""

from __future__ import annotations

import numpy as np

from .sampling_types import ScalerState


def robust_scale(
    matrix: np.ndarray,
    feature_names: tuple[str, ...],
) -> tuple[np.ndarray, ScalerState]:
    values = np.asarray(matrix, dtype=float)
    if values.ndim != 2 or values.shape[1] != len(feature_names):
        raise ValueError("Feature matrix shape does not match feature names")
    median = np.median(values, axis=0)
    lower = np.quantile(values, 0.25, axis=0)
    upper = np.quantile(values, 0.75, axis=0)
    iqr = upper - lower
    safe_iqr = np.where(iqr > 1.0e-12, iqr, 1.0)
    return (values - median) / safe_iqr, ScalerState(feature_names, median, safe_iqr)


def apply_group_weights(
    scaled: np.ndarray,
    groups: tuple[str, ...],
    *,
    radius_weight: float,
    structure_weight: float,
) -> np.ndarray:
    weights = np.asarray(
        [radius_weight if group == "radius" else structure_weight for group in groups],
        dtype=float,
    )
    return np.asarray(scaled, dtype=float) * np.sqrt(weights)[None, :]

