"""Direction-preserving branch geometry calculations."""

from __future__ import annotations

import numpy as np


def cumulative_length(points: np.ndarray) -> np.ndarray:
    values = np.asarray(points, dtype=float)
    if len(values) == 0:
        return np.empty(0, dtype=float)
    return np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(values, axis=0), axis=1))))


def _remove_consecutive_duplicates(
    points: np.ndarray, radii: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    if len(points) <= 1:
        return points.copy(), radii.copy()
    keep = np.concatenate(([True], np.linalg.norm(np.diff(points, axis=0), axis=1) > 1.0e-12))
    return points[keep], radii[keep]


def fill_invalid_radius(radius: np.ndarray, arc: np.ndarray) -> np.ndarray:
    """Return a derived positive radius array without changing raw SWC values."""

    values = np.asarray(radius, dtype=float).copy()
    valid = np.isfinite(values) & (values > 0)
    if np.all(valid):
        return values
    if not np.any(valid):
        return np.full_like(values, np.nan)
    values[~valid] = np.interp(arc[~valid], arc[valid], values[valid])
    return values


def derive_branch_geometry(
    points_um: np.ndarray,
    radius_um: np.ndarray,
    *,
    smoothing_enabled: bool,
    smoothing_window_points: int,
    resample_step_um: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float, float, float, float]:
    """Create a derived sequence while retaining parent-to-child point order."""

    raw_points = np.asarray(points_um, dtype=float)
    raw_radius = np.asarray(radius_um, dtype=float)
    raw_arc = cumulative_length(raw_points)
    raw_length = float(raw_arc[-1]) if len(raw_arc) else 0.0
    points, radius = _remove_consecutive_duplicates(raw_points, raw_radius)
    arc = cumulative_length(points)
    radius = fill_invalid_radius(radius, arc)

    if len(points) >= 2 and arc[-1] > 0:
        count = max(2, int(np.ceil(arc[-1] / resample_step_um)) + 1)
        sample_arc = np.linspace(0.0, float(arc[-1]), count)
        derived = np.column_stack(
            [np.interp(sample_arc, arc, points[:, axis]) for axis in range(3)]
        )
        derived_radius = (
            np.interp(sample_arc, arc, radius)
            if np.all(np.isfinite(radius))
            else np.full(count, np.nan)
        )
    else:
        derived = points.copy()
        derived_radius = radius.copy()

    if smoothing_enabled and len(derived) >= smoothing_window_points:
        half = smoothing_window_points // 2
        kernel = np.full(smoothing_window_points, 1.0 / smoothing_window_points)
        padded = np.pad(derived, ((half, half), (0, 0)), mode="edge")
        smoothed = np.column_stack(
            [np.convolve(padded[:, axis], kernel, mode="valid") for axis in range(3)]
        )
        smoothed[0] = derived[0]
        smoothed[-1] = derived[-1]
        derived = smoothed

    derived_arc = cumulative_length(derived)
    derived_length = float(derived_arc[-1]) if len(derived_arc) else 0.0
    directions = np.zeros_like(derived)
    if len(derived) >= 2:
        tangent = np.gradient(derived, axis=0)
        magnitude = np.linalg.norm(tangent, axis=1)
        nonzero = magnitude > 1.0e-12
        directions[nonzero] = tangent[nonzero] / magnitude[nonzero, None]
    chord = float(np.linalg.norm(derived[-1] - derived[0])) if len(derived) >= 2 else 0.0
    tortuosity = derived_length / chord if chord > 1.0e-12 else float("nan")
    finite_radius = derived_radius[np.isfinite(derived_radius) & (derived_radius > 0)]
    mean_radius = float(np.mean(finite_radius)) if len(finite_radius) else float("nan")
    return (
        derived,
        derived_radius,
        directions,
        raw_length,
        derived_length,
        tortuosity,
        mean_radius,
    )
