"""Physical geometry calculations for coarse vascular branches."""

from __future__ import annotations

import numpy as np
from scipy.signal import savgol_filter


def indices_to_lps(
    indices_xyz: np.ndarray,
    origin_lps_um: tuple[float, float, float],
    spacing_um: tuple[float, float, float],
) -> np.ndarray:
    return np.asarray(origin_lps_um, dtype=float) + np.asarray(indices_xyz, dtype=float) * np.asarray(
        spacing_um, dtype=float
    )


def cumulative_length(points: np.ndarray) -> np.ndarray:
    if len(points) == 0:
        return np.empty(0, dtype=float)
    if len(points) == 1:
        return np.zeros(1, dtype=float)
    steps = np.linalg.norm(np.diff(points, axis=0), axis=1)
    return np.concatenate(([0.0], np.cumsum(steps)))


def _resample_sequence(
    points: np.ndarray,
    values: np.ndarray,
    step_um: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    source_s = cumulative_length(points)
    total = float(source_s[-1])
    if len(points) <= 2 or total <= 0:
        return points.copy(), values.copy(), source_s
    sample_count = max(2, int(np.ceil(total / step_um)) + 1)
    target_s = np.linspace(0.0, total, sample_count)
    target_points = np.column_stack(
        [np.interp(target_s, source_s, points[:, axis]) for axis in range(3)]
    )
    target_values = np.interp(target_s, source_s, values)
    return target_points, target_values, target_s


def _odd_window(requested: int, count: int) -> int:
    window = min(requested, count if count % 2 else count - 1)
    return window if window >= 3 else 0


def smooth_and_measure(
    points_raw: np.ndarray,
    radius_raw: np.ndarray,
    *,
    resample_step_um: float,
    smoothing_enabled: bool,
    smoothing_window_points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float]:
    points, radii, _ = _resample_sequence(points_raw, radius_raw, resample_step_um)
    unsmoothed_resampled = points.copy()
    window = _odd_window(smoothing_window_points, len(points))
    if smoothing_enabled and window:
        polynomial_order = min(2, window - 1)
        points = savgol_filter(points, window, polynomial_order, axis=0, mode="interp")
        radii = savgol_filter(radii, window, polynomial_order, mode="interp")
        points[0] = points_raw[0]
        points[-1] = points_raw[-1]
        radii[0] = radius_raw[0]
        radii[-1] = radius_raw[-1]
        radii = np.maximum(radii, np.finfo(float).eps)
    arc_length = cumulative_length(points)
    deviation = float(np.max(np.linalg.norm(points - unsmoothed_resampled, axis=1)))
    if len(points) < 2 or arc_length[-1] <= 0:
        directions = np.zeros_like(points)
        curvature = np.zeros(len(points), dtype=float)
        return points, arc_length, radii, directions, curvature, deviation
    directions = np.gradient(points, arc_length, axis=0, edge_order=1)
    norms = np.linalg.norm(directions, axis=1)
    valid = norms > np.finfo(float).eps
    directions[valid] /= norms[valid, None]
    directions[~valid] = 0.0
    curvature_vectors = np.gradient(directions, arc_length, axis=0, edge_order=1)
    curvature = np.linalg.norm(curvature_vectors, axis=1)
    curvature[~np.isfinite(curvature)] = 0.0
    return points, arc_length, radii, directions, curvature, deviation


def tortuosity(length_um: float, straight_distance_um: float) -> float:
    if straight_distance_um <= np.finfo(float).eps:
        return float("nan")
    return max(1.0, float(length_um / straight_distance_um))
