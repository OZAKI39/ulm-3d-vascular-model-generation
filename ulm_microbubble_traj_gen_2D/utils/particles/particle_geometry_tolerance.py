"""Shared floating-point tolerance for authoritative particle geometry."""

from __future__ import annotations

import numpy as np


PARTICLE_GEOMETRY_ROUNDOFF_FACTOR = 8192.0


def particle_geometry_roundoff_tolerance_um(*values: object) -> float:
    """Return one conservative coordinate-scale roundoff tolerance in micrometres.

    The tolerance may decide whether a tiny negative gap is eligible for an
    inward projection. It never licenses storing or consuming a negative
    physical gap.
    """

    scale = 1.0
    for value in values:
        array = np.asarray(value, dtype=np.float64)
        if array.size == 0:
            continue
        if np.any(~np.isfinite(array)):
            raise ValueError("Particle geometry tolerance inputs must be finite.")
        scale = max(scale, float(np.max(np.abs(array))))
    return (
        PARTICLE_GEOMETRY_ROUNDOFF_FACTOR
        * np.finfo(np.float64).eps
        * scale
    )


def project_roundoff_negative_wall_gap_xz_um(
    position_xz_um: object,
    radius_um: float,
    boundary_geometry: object,
    *,
    maximum_iterations: int = 8,
) -> tuple[np.ndarray, float, float]:
    """Repair only a coordinate-roundoff negative solid-wall gap.

    The returned point is accepted as repaired only when an authoritative wall
    re-query reports a non-negative gap.  A materially penetrating point is
    returned unchanged with its negative gap, so callers cannot hide a real
    geometry or integration failure.
    """

    point = np.asarray(position_xz_um, dtype=np.float64).copy()
    if point.shape != (2,) or np.any(~np.isfinite(point)):
        raise ValueError("position_xz_um must contain two finite coordinates.")
    radius = float(radius_um)
    if not np.isfinite(radius) or radius <= 0.0:
        raise ValueError("radius_um must be finite and positive.")
    iteration_limit = int(maximum_iterations)
    if iteration_limit < 1:
        raise ValueError("maximum_iterations must be positive.")

    spacing_um = float(boundary_geometry.domain.spacing_um)
    tolerance_um = particle_geometry_roundoff_tolerance_um(
        point,
        radius,
        spacing_um,
    )

    state = boundary_geometry.exact_solid_wall_state_xz_um(point)
    gap_um = float(np.asarray(state.distance_um, dtype=np.float64))
    gap_um -= radius
    if not np.isfinite(gap_um):
        raise ValueError("The authoritative solid-wall gap must be finite.")
    if gap_um >= 0.0 or gap_um < -tolerance_um:
        return point, gap_um, tolerance_um

    for _ in range(iteration_limit):
        normal = np.asarray(state.inward_normal_xz, dtype=np.float64)
        normal_norm_squared = float(np.dot(normal, normal))
        if (
            normal.shape != (2,)
            or np.any(~np.isfinite(normal))
            or not np.isfinite(normal_norm_squared)
            or normal_norm_squared <= 0.0
        ):
            break
        point += (
            (-gap_um + tolerance_um)
            * normal
            / normal_norm_squared
        )
        state = boundary_geometry.exact_solid_wall_state_xz_um(point)
        gap_um = float(np.asarray(state.distance_um, dtype=np.float64))
        gap_um -= radius
        if not np.isfinite(gap_um):
            raise ValueError("The authoritative solid-wall gap must be finite.")
        if gap_um >= 0.0:
            return point, gap_um, tolerance_um
        if gap_um < -tolerance_um:
            break

    return point, gap_um, tolerance_um
