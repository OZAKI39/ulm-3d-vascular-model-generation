"""Exact line-segment intersection with an axis-aligned ROI boundary."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True, slots=True)
class ClippedSegment:
    start_um: np.ndarray
    end_um: np.ndarray
    start_t: float
    end_t: float
    start_face: str | None
    end_face: str | None


def boundary_face(
    point_um: np.ndarray,
    bbox_min_um: np.ndarray,
    bbox_max_um: np.ndarray,
    *,
    tolerance: float = 1.0e-7,
) -> str:
    """Return a deterministic boundary-face label for a point on a box."""

    point = np.asarray(point_um, dtype=float)
    lower = np.asarray(bbox_min_um, dtype=float)
    upper = np.asarray(bbox_max_um, dtype=float)
    labels = (("x_min", 0, lower[0]), ("x_max", 0, upper[0]),
              ("y_min", 1, lower[1]), ("y_max", 1, upper[1]),
              ("z_min", 2, lower[2]), ("z_max", 2, upper[2]))
    for label, axis, value in labels:
        if abs(float(point[axis] - value)) <= tolerance:
            return label
    return "component_boundary"


def point_in_box(
    point_um: np.ndarray,
    bbox_min_um: np.ndarray,
    bbox_max_um: np.ndarray,
    *,
    tolerance: float = 1.0e-9,
) -> bool:
    point = np.asarray(point_um, dtype=float)
    return bool(
        np.all(point >= np.asarray(bbox_min_um, dtype=float) - tolerance)
        and np.all(point <= np.asarray(bbox_max_um, dtype=float) + tolerance)
    )


def clip_segment_to_box(
    start_um: np.ndarray,
    end_um: np.ndarray,
    bbox_min_um: np.ndarray,
    bbox_max_um: np.ndarray,
    *,
    tolerance: float = 1.0e-12,
) -> ClippedSegment | None:
    """Clip a segment with the slab method and retain exact interpolation parameters."""

    start = np.asarray(start_um, dtype=float)
    end = np.asarray(end_um, dtype=float)
    lower = np.asarray(bbox_min_um, dtype=float)
    upper = np.asarray(bbox_max_um, dtype=float)
    direction = end - start
    enter = 0.0
    leave = 1.0
    for axis in range(3):
        if abs(float(direction[axis])) <= tolerance:
            if start[axis] < lower[axis] - tolerance or start[axis] > upper[axis] + tolerance:
                return None
            continue
        first = float((lower[axis] - start[axis]) / direction[axis])
        second = float((upper[axis] - start[axis]) / direction[axis])
        axis_enter, axis_leave = min(first, second), max(first, second)
        enter = max(enter, axis_enter)
        leave = min(leave, axis_leave)
        if enter > leave + tolerance:
            return None
    enter = float(np.clip(enter, 0.0, 1.0))
    leave = float(np.clip(leave, 0.0, 1.0))
    clipped_start = start + enter * direction
    clipped_end = start + leave * direction
    if float(np.linalg.norm(clipped_end - clipped_start)) <= tolerance:
        return None
    return ClippedSegment(
        clipped_start,
        clipped_end,
        enter,
        leave,
        boundary_face(clipped_start, lower, upper) if enter > tolerance else None,
        boundary_face(clipped_end, lower, upper) if leave < 1.0 - tolerance else None,
    )

