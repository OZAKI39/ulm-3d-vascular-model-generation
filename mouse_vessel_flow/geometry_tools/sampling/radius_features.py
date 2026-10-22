"""Arc-length-weighted radius distributions for clipped ROI edges."""

from __future__ import annotations

import numpy as np

from .sampling_types import ROIRecord


def arc_length_radius_samples(roi: ROIRecord) -> tuple[np.ndarray, np.ndarray]:
    """Return segment-midpoint radii and physical arc-length weights."""

    points = np.asarray(roi.local_edge_points_um, dtype=float)
    radii = np.asarray(roi.local_edge_radius_um, dtype=float)
    lengths = np.linalg.norm(points[:, 1] - points[:, 0], axis=1)
    midpoint_radius = np.mean(radii, axis=1)
    valid = (
        np.isfinite(lengths)
        & (lengths > 1.0e-12)
        & np.isfinite(midpoint_radius)
        & (midpoint_radius > 0)
    )
    return midpoint_radius[valid], lengths[valid]


def weighted_quantile(
    values: np.ndarray,
    weights: np.ndarray,
    quantiles: np.ndarray | tuple[float, ...],
) -> np.ndarray:
    """Inverse of the weighted empirical CDF defined in the sampling plan."""

    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    q = np.asarray(quantiles, dtype=float)
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(valid):
        return np.full(q.shape, np.nan, dtype=float)
    order = np.argsort(values[valid], kind="stable")
    sorted_values = values[valid][order]
    cumulative = np.cumsum(weights[valid][order])
    targets = np.clip(q, 0.0, 1.0) * cumulative[-1]
    indices = np.searchsorted(cumulative, targets, side="left")
    return sorted_values[np.clip(indices, 0, len(sorted_values) - 1)]


def compute_radius_features(roi: ROIRecord) -> dict[str, float]:
    values, weights = arc_length_radius_samples(roi)
    if not len(values):
        return {name: float("nan") for name in (
            "r05", "r10", "r25", "r50", "r75", "r90", "r95",
            "radius_min_um", "radius_max_um", "radius_mean_um", "radius_std_um",
        )}
    quantile_values = weighted_quantile(
        values,
        weights,
        np.asarray((0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95)),
    )
    mean = float(np.average(values, weights=weights))
    variance = float(np.average((values - mean) ** 2, weights=weights))
    return {
        "r05": float(quantile_values[0]),
        "r10": float(quantile_values[1]),
        "r25": float(quantile_values[2]),
        "r50": float(quantile_values[3]),
        "r75": float(quantile_values[4]),
        "r90": float(quantile_values[5]),
        "r95": float(quantile_values[6]),
        "radius_min_um": float(np.min(values)),
        "radius_max_um": float(np.max(values)),
        "radius_mean_um": mean,
        "radius_std_um": float(np.sqrt(max(variance, 0.0))),
    }


def aggregate_radius_distribution(rois: list[ROIRecord] | tuple[ROIRecord, ...]) -> tuple[np.ndarray, np.ndarray]:
    values: list[np.ndarray] = []
    weights: list[np.ndarray] = []
    for roi in rois:
        local_values, local_weights = arc_length_radius_samples(roi)
        if len(local_values):
            values.append(local_values)
            weights.append(local_weights)
    if not values:
        return np.empty(0, dtype=float), np.empty(0, dtype=float)
    return np.concatenate(values), np.concatenate(weights)

