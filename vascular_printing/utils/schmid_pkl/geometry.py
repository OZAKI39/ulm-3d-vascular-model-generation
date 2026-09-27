"""Step 2 equivalent: lossless centerline assembly and derived geometry measurements."""

from __future__ import annotations

import numpy as np

from ..graph.geometry import cumulative_length, smooth_and_measure, tortuosity


def concatenate_edge_geometry(
    ordered: list[tuple[np.ndarray, np.ndarray]],
) -> tuple[np.ndarray, np.ndarray]:
    """Join edge paths while retaining mismatched junction samples for transparent QC."""

    point_parts: list[np.ndarray] = []
    radius_parts: list[np.ndarray] = []
    for index, (points, diameters) in enumerate(ordered):
        radii = np.asarray(diameters, dtype=float) / 2.0
        points = np.asarray(points, dtype=float)
        if index and len(point_parts) and np.allclose(
            point_parts[-1][-1], points[0], rtol=0.0, atol=1.0e-9
        ):
            points = points[1:]
            radii = radii[1:]
        if len(points):
            point_parts.append(points)
            radius_parts.append(radii)
    if not point_parts:
        return np.empty((0, 3), dtype=float), np.empty(0, dtype=float)
    return np.concatenate(point_parts), np.concatenate(radius_parts)


def measure_branch_geometry(
    points_raw: np.ndarray,
    radius_raw: np.ndarray,
    *,
    resample_step_um: float,
    smoothing_enabled: bool,
    smoothing_window_points: int,
) -> dict[str, np.ndarray | float | None]:
    arc_raw = cumulative_length(points_raw)
    geometric_length = float(arc_raw[-1]) if len(arc_raw) else 0.0
    straight_distance = (
        float(np.linalg.norm(points_raw[-1] - points_raw[0])) if len(points_raw) >= 2 else 0.0
    )
    if geometric_length <= np.finfo(float).eps or len(points_raw) < 2:
        points_smooth = points_raw.copy()
        radius_smooth = radius_raw.copy()
        arc_smooth = arc_raw.copy()
        directions = np.zeros_like(points_raw)
        curvature = np.zeros(len(points_raw), dtype=float)
        deviation = 0.0
    else:
        (
            points_smooth,
            arc_smooth,
            radius_smooth,
            directions,
            curvature,
            deviation,
        ) = smooth_and_measure(
            points_raw,
            radius_raw,
            resample_step_um=resample_step_um,
            smoothing_enabled=smoothing_enabled,
            smoothing_window_points=smoothing_window_points,
        )
    raw_tortuosity = tortuosity(geometric_length, straight_distance)
    smooth_length = float(arc_smooth[-1]) if len(arc_smooth) else 0.0
    smooth_tortuosity = tortuosity(smooth_length, straight_distance)
    return {
        "arc_length_raw_um": arc_raw,
        "points_smoothed_um": points_smooth,
        "radius_smoothed_um": radius_smooth,
        "arc_length_smoothed_um": arc_smooth,
        "local_direction_smoothed": directions,
        "curvature_smoothed_per_um": curvature,
        "geometric_length_um": geometric_length,
        "straight_distance_um": straight_distance,
        "tortuosity_raw": float(raw_tortuosity) if np.isfinite(raw_tortuosity) else None,
        "tortuosity_smoothed": float(smooth_tortuosity) if np.isfinite(smooth_tortuosity) else None,
        "smoothing_max_deviation_um": float(deviation),
    }
