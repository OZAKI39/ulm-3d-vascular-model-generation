"""Dimensionless ray/triangle geometry and research-stage contracts.

The numerical kernel mirrors the final Seeder/TreElm implementation: translate
the triangle to the ray origin, scale every coordinate by the lattice-link
component scale, and evaluate Moller--Trumbore in that local O(1) frame.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np


DIMENSIONLESS_TOLERANCE_FACTOR = 64.0
NO_INTERSECTION_Q = -1.0


@dataclass(frozen=True, slots=True)
class RayTriangleResult:
    hit: bool
    q: float | None
    barycentric_u: float | None
    barycentric_v: float | None
    classification: str


def dimensionless_ray_triangle(
    origin_world: np.ndarray,
    link_world: np.ndarray,
    triangle_world: np.ndarray,
    *,
    tolerance_factor: float = DIMENSIONLESS_TOLERANCE_FACTOR,
) -> RayTriangleResult:
    """Intersect one lattice link ray with a triangle in local coordinates.

    ``q`` is the fraction along the supplied, non-normalized lattice link.  A
    diagonal ``link_world`` therefore keeps exactly the Seeder/TreElm qVal
    contract: ``origin + q * link_world`` is the intersection point.
    """

    origin = np.asarray(origin_world, dtype=np.float64)
    link = np.asarray(link_world, dtype=np.float64)
    triangle = np.asarray(triangle_world, dtype=np.float64)
    if origin.shape != (3,) or link.shape != (3,) or triangle.shape != (3, 3):
        raise ValueError("origin/link/triangle shapes must be (3,), (3,), (3,3)")
    scale = float(np.max(np.abs(link)))
    if not np.isfinite(scale) or scale <= 0.0:
        return RayTriangleResult(False, None, None, None, "ZERO_LENGTH_LINK")

    direction = link / scale
    vertices = (triangle - origin) / scale
    edge1 = vertices[1] - vertices[0]
    edge2 = vertices[2] - vertices[0]
    edge_product = float(np.linalg.norm(edge1) * np.linalg.norm(edge2))
    direction_norm = float(np.linalg.norm(direction))
    if edge_product <= np.finfo(np.float64).tiny:
        return RayTriangleResult(False, None, None, None, "DEGENERATE_TRIANGLE")

    pvec = np.cross(direction, edge2)
    determinant = float(np.dot(edge1, pvec))
    tolerance = float(tolerance_factor) * np.finfo(np.float64).eps
    if abs(determinant) <= tolerance * edge_product * direction_norm:
        return RayTriangleResult(False, None, None, None, "PARALLEL_OR_GRAZING")

    inverse = 1.0 / determinant
    tvec = -vertices[0]
    u = float(np.dot(tvec, pvec) * inverse)
    if u < -tolerance or u > 1.0 + tolerance:
        return RayTriangleResult(False, None, u, None, "OUTSIDE_TRIANGLE")
    qvec = np.cross(tvec, edge1)
    v = float(np.dot(direction, qvec) * inverse)
    if v < -tolerance or u + v > 1.0 + tolerance:
        return RayTriangleResult(False, None, u, v, "OUTSIDE_TRIANGLE")
    q = float(np.dot(edge2, qvec) * inverse)
    if q < -tolerance:
        return RayTriangleResult(False, q, u, v, "BEHIND_RAY")
    q = max(q, 0.0)
    barycentric = np.asarray((1.0 - u - v, u, v), dtype=np.float64)
    near_count = int(np.count_nonzero(np.abs(barycentric) <= tolerance))
    if near_count >= 2:
        classification = "HIT_VERTEX"
    elif near_count == 1:
        classification = "HIT_EDGE"
    else:
        classification = "HIT_INTERIOR"
    return RayTriangleResult(True, q, u, v, classification)


def compare_qvalue_support(
    seeder_q: np.ndarray, oracle_q: np.ndarray
) -> dict[str, Any]:
    """Separate intersection classification mismatches from numeric q error."""

    seeder = np.asarray(seeder_q, dtype=np.float64)
    oracle = np.asarray(oracle_q, dtype=np.float64)
    if seeder.shape != oracle.shape:
        raise ValueError("q arrays must have identical shapes")
    seeder_hit = np.isfinite(seeder) & (seeder >= 0.0)
    oracle_hit = np.isfinite(oracle) & (oracle >= 0.0)
    both = seeder_hit & oracle_hit
    neither = ~seeder_hit & ~oracle_hit
    seeder_only = seeder_hit & ~oracle_hit
    oracle_only = ~seeder_hit & oracle_hit
    error = seeder[both] - oracle[both]
    absolute = np.abs(error)
    return {
        "TRUE_INTERSECTION_BOTH": int(np.count_nonzero(both)),
        "NO_INTERSECTION_BOTH": int(np.count_nonzero(neither)),
        "SEEDER_ONLY_INTERSECTION": int(np.count_nonzero(seeder_only)),
        "ORACLE_ONLY_INTERSECTION": int(np.count_nonzero(oracle_only)),
        "numeric_q_error": {
            "count": int(len(error)),
            "mean_bias": float(np.mean(error)) if len(error) else None,
            "median_absolute": float(np.median(absolute)) if len(error) else None,
            "rms": float(np.sqrt(np.mean(error * error))) if len(error) else None,
            "p95": float(np.percentile(absolute, 95.0)) if len(error) else None,
            "max": float(np.max(absolute)) if len(error) else None,
        },
    }


@dataclass(slots=True)
class OperationalRetryBudget:
    stage_limit: int = 5
    total_limit: int = 20
    total_attempts: int = 0
    stage_attempts: dict[str, int] = field(default_factory=dict)

    def consume(self, stage: str) -> tuple[int, int]:
        stage_count = self.stage_attempts.get(stage, 0) + 1
        total = self.total_attempts + 1
        if stage_count > self.stage_limit or total > self.total_limit:
            raise RuntimeError("CFD_FLOW_OPERATIONAL_INFRASTRUCTURE_BLOCKED")
        self.stage_attempts[stage] = stage_count
        self.total_attempts = total
        return stage_count, total


def completed_stage_reusable(
    checkpoint: Mapping[str, Any], current_input_hashes: Mapping[str, str]
) -> bool:
    return bool(
        checkpoint.get("status") == "PASS"
        and dict(checkpoint.get("input_hashes", {})) == dict(current_input_hashes)
    )


def semantic_files_success(
    root: Path, required_relative_paths: tuple[str, ...]
) -> dict[str, Any]:
    files = {
        name: {
            "exists": (Path(root) / name).is_file(),
            "size_bytes": (Path(root) / name).stat().st_size
            if (Path(root) / name).is_file()
            else 0,
        }
        for name in required_relative_paths
    }
    success = all(row["exists"] and row["size_bytes"] > 0 for row in files.values())
    return {"semantic_success": bool(success), "files": files}
