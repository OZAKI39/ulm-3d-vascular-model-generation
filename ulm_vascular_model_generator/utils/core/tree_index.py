"""Three-dimensional spatial index and segment-distance utilities."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .models import Vessel

try:  # pragma: no cover
    from scipy.spatial import cKDTree
except Exception:  # pragma: no cover
    cKDTree = None


@dataclass
class SegmentIndex:
    ids: np.ndarray
    starts: np.ndarray
    ends: np.ndarray
    midpoints: np.ndarray
    half_lengths: np.ndarray
    radii: np.ndarray
    main_trunk: np.ndarray
    branchable: np.ndarray
    midpoint_tree: object | None

    @property
    def max_half_length(self) -> float:
        return float(np.max(self.half_lengths)) if self.half_lengths.size else 0.0

    @property
    def max_radius(self) -> float:
        return float(np.max(self.radii)) if self.radii.size else 0.0


def build_segment_index(vessels: list[Vessel]) -> SegmentIndex:
    """Build a full XYZ midpoint index for the current vessel segments."""
    if not vessels:
        empty_i = np.empty(0, dtype=int)
        empty_xyz = np.empty((0, 3), dtype=float)
        empty_f = np.empty(0, dtype=float)
        empty_b = np.empty(0, dtype=bool)
        return SegmentIndex(
            empty_i,
            empty_xyz,
            empty_xyz,
            empty_xyz,
            empty_f,
            empty_f,
            empty_b,
            empty_b,
            None,
        )

    ids = np.fromiter((v.vid for v in vessels), dtype=int, count=len(vessels))
    starts = np.asarray([v.x_p for v in vessels], dtype=float)
    ends = np.asarray([v.x_d for v in vessels], dtype=float)
    midpoints = 0.5 * (starts + ends)
    half_lengths = 0.5 * np.linalg.norm(ends - starts, axis=1)
    radii = np.asarray([v.radius for v in vessels], dtype=float)
    main_trunk = np.asarray([v.is_main_trunk for v in vessels], dtype=bool)
    branchable = np.asarray(
        [v.branching_mode != "non_branching" for v in vessels], dtype=bool
    )
    midpoint_tree = (
        cKDTree(midpoints) if cKDTree is not None and len(vessels) else None
    )
    return SegmentIndex(
        ids,
        starts,
        ends,
        midpoints,
        half_lengths,
        radii,
        main_trunk,
        branchable,
        midpoint_tree,
    )


def point_segment_distances_3d(
    p: np.ndarray, starts: np.ndarray, ends: np.ndarray
) -> np.ndarray:
    """Return the Euclidean distance from one XYZ point to many XYZ segments."""
    if starts.size == 0:
        return np.empty(0, dtype=float)
    p = np.asarray(p, dtype=float)
    starts = np.asarray(starts, dtype=float)
    ends = np.asarray(ends, dtype=float)
    ab = ends - starts
    denom = np.einsum("ij,ij->i", ab, ab)
    t = np.zeros(len(starts), dtype=float)
    safe = denom > 1.0e-30
    t[safe] = np.einsum(
        "ij,ij->i", p - starts[safe], ab[safe]
    ) / denom[safe]
    t = np.clip(t, 0.0, 1.0)
    closest = starts + t[:, None] * ab
    return np.linalg.norm(p - closest, axis=1)


def point_segment_distances_3d_many_points(
    points: np.ndarray, a: np.ndarray, b: np.ndarray
) -> np.ndarray:
    """Return distances from many XYZ points to one XYZ segment."""
    if points.size == 0:
        return np.empty(0, dtype=float)
    points = np.asarray(points, dtype=float)
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    ab = b - a
    denom = float(np.dot(ab, ab))
    if denom < 1.0e-30:
        return np.linalg.norm(points - a, axis=1)
    t = np.einsum("ij,j->i", points - a, ab) / denom
    t = np.clip(t, 0.0, 1.0)
    closest = a + t[:, None] * ab
    return np.linalg.norm(points - closest, axis=1)


def _segment_segment_distance(
    p0: np.ndarray, p1: np.ndarray, q0: np.ndarray, q1: np.ndarray
) -> float:
    """Shortest distance between two finite segments in any Euclidean dimension."""
    u = p1 - p0
    v = q1 - q0
    w = p0 - q0
    a = float(np.dot(u, u))
    b = float(np.dot(u, v))
    c = float(np.dot(v, v))
    d = float(np.dot(u, w))
    e = float(np.dot(v, w))
    eps = 1.0e-30

    if a <= eps and c <= eps:
        return float(np.linalg.norm(p0 - q0))
    if a <= eps:
        t = np.clip(e / c, 0.0, 1.0)
        return float(np.linalg.norm(p0 - (q0 + t * v)))
    if c <= eps:
        s = np.clip(-d / a, 0.0, 1.0)
        return float(np.linalg.norm((p0 + s * u) - q0))

    denominator = a * c - b * b
    if denominator > eps:
        s = np.clip((b * e - c * d) / denominator, 0.0, 1.0)
    else:
        s = 0.0
    t = (b * s + e) / c
    if t < 0.0:
        t = 0.0
        s = np.clip(-d / a, 0.0, 1.0)
    elif t > 1.0:
        t = 1.0
        s = np.clip((b - d) / a, 0.0, 1.0)
    return float(np.linalg.norm((p0 + s * u) - (q0 + t * v)))


def segment_distances_3d(
    a1: np.ndarray, a2: np.ndarray, starts: np.ndarray, ends: np.ndarray
) -> np.ndarray:
    """Return shortest distances from one XYZ segment to many XYZ segments."""
    if starts.size == 0:
        return np.empty(0, dtype=float)
    return np.fromiter(
        (
            _segment_segment_distance(
                np.asarray(a1, dtype=float),
                np.asarray(a2, dtype=float),
                start,
                end,
            )
            for start, end in zip(starts, ends, strict=True)
        ),
        dtype=float,
        count=len(starts),
    )


def segment_intersections_3d(
    a1: np.ndarray,
    a2: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
    eps: float = 1.0e-7,
) -> np.ndarray:
    """Treat two centreline segments as intersecting when their distance is tiny."""
    return segment_distances_3d(a1, a2, starts, ends) <= float(eps)

