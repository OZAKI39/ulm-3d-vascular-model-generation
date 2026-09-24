"""Population contract for the research-only adaptive-flux Musubi patch.

The target is the mean of one mass-flow sample per *active* local boundary
element.  Musubi's ``remove_solid_in_bc`` may make that population smaller
than the mesh-time ``globBC%nElems_total`` value, so the denominator is the
MPI sum of the same local active counts used for sampling and pressure_eq.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import Any

import numpy as np


class AdaptiveTargetPopulationError(ValueError):
    """Raised when a rank-local adaptive target population is inconsistent."""


def audit_adaptive_target_population(
    samples_by_rank: Sequence[Sequence[float]],
    point_indices_by_rank: Sequence[Sequence[int]],
    active_links_by_rank: Sequence[Sequence[int]],
    *,
    stale_mesh_population: int | None = None,
) -> dict[str, Any]:
    """Validate and reduce the source-defined active sampling population.

    Duplicate positive point indices are legal: a spacetime variable may map
    multiple active boundary elements to the same registered evaluation point.
    Missing/non-positive indices, non-finite samples, and active-list elements
    without a pressure boundary link fail rather than silently biasing a mean.
    """

    rank_count = len(samples_by_rank)
    if len(point_indices_by_rank) != rank_count or len(active_links_by_rank) != rank_count:
        raise AdaptiveTargetPopulationError("rank population arrays differ in length")

    per_rank: list[dict[str, Any]] = []
    all_samples: list[float] = []
    for rank, (raw_samples, raw_indices, raw_links) in enumerate(
        zip(samples_by_rank, point_indices_by_rank, active_links_by_rank, strict=True)
    ):
        samples = np.asarray(raw_samples, dtype=np.float64)
        indices = np.asarray(raw_indices, dtype=np.int64)
        links = np.asarray(raw_links, dtype=np.int64)
        if samples.ndim != 1 or indices.ndim != 1 or links.ndim != 1:
            raise AdaptiveTargetPopulationError(f"rank {rank}: populations must be one-dimensional")
        if not (samples.size == indices.size == links.size):
            raise AdaptiveTargetPopulationError(f"rank {rank}: sample population mismatch")
        if np.any(indices <= 0):
            raise AdaptiveTargetPopulationError(f"rank {rank}: invalid point index")
        if not np.all(np.isfinite(samples)):
            raise AdaptiveTargetPopulationError(f"rank {rank}: non-finite mass-flow sample")
        if np.any(links <= 0):
            raise AdaptiveTargetPopulationError(f"rank {rank}: inactive element in active list")
        all_samples.extend(float(value) for value in samples)
        unique = int(np.unique(indices).size)
        per_rank.append(
            {
                "rank": rank,
                "active_count": int(samples.size),
                "point_count": int(indices.size),
                "point_unique_count": unique,
                "point_duplicate_count": int(indices.size) - unique,
                "active_link_count": int(np.sum(links, dtype=np.int64)),
                "mass_local": math.fsum(float(value) for value in samples),
            }
        )

    active_count = sum(int(item["active_count"]) for item in per_rank)
    if active_count <= 0:
        raise AdaptiveTargetPopulationError("adaptive target has no active boundary elements")
    mass_global = math.fsum(all_samples)
    target_mass_flow = mass_global / active_count
    return {
        "status": "PASS",
        "per_rank": per_rank,
        "active_count_global": active_count,
        "sample_count_global": len(all_samples),
        "stale_mesh_population": stale_mesh_population,
        "mass_global": mass_global,
        "target_mass_flow": target_mass_flow,
    }
