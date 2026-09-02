"""Fused Eulerian-field sampling for the particle RHS hot path.

Every continuous particle stage needs the same grid cell and four bilinear
weights for several CFD fields.  Computing those indices independently for
each field creates many short NumPy arrays and repeated clipping operations.
This module computes the shared interpolation geometry once per particle and
reads all required fields in one compiled loop.
"""

from __future__ import annotations

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - the production environment includes Numba
    njit = None


def sample_bilinear(
    field: np.ndarray,
    positions_grid: np.ndarray,
) -> np.ndarray:
    """Sample scalar, vector, or tensor grid values at fractional X-Z positions."""

    values = np.asarray(field)
    positions = np.asarray(positions_grid, dtype=np.float64)
    nx, nz = values.shape[:2]
    gx = np.clip(positions[:, 0], 0.0, nx - 1.0)
    gz = np.clip(positions[:, 1], 0.0, nz - 1.0)
    i0 = np.clip(np.floor(gx).astype(np.int64), 0, nx - 2)
    j0 = np.clip(np.floor(gz).astype(np.int64), 0, nz - 2)
    tx = gx - i0
    tz = gz - j0
    trailing_axes = (1,) * (values.ndim - 2)
    w00 = ((1.0 - tx) * (1.0 - tz)).reshape((-1, *trailing_axes))
    w10 = (tx * (1.0 - tz)).reshape((-1, *trailing_axes))
    w01 = ((1.0 - tx) * tz).reshape((-1, *trailing_axes))
    w11 = (tx * tz).reshape((-1, *trailing_axes))
    return (
        w00 * values[i0, j0]
        + w10 * values[i0 + 1, j0]
        + w01 * values[i0, j0 + 1]
        + w11 * values[i0 + 1, j0 + 1]
    )


def sample_regular_grid_fields_batch(
    positions_grid: np.ndarray,
    active: np.ndarray,
    velocity_xz_um_s: np.ndarray,
    velocity_gradient_s_inv: np.ndarray,
    dynamic_viscosity_pa_s: np.ndarray,
    local_vessel_radius_um: np.ndarray,
    wall_shear_stress_pa: np.ndarray,
    *,
    use_numba: bool,
) -> tuple[
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
    np.ndarray,
]:
    """Sample the far-field cache and grid-carried material fields."""

    kernel = (
        _sample_regular_grid_fields_numba
        if use_numba and njit is not None
        else _sample_regular_grid_fields_kernel
    )
    return kernel(
        np.ascontiguousarray(positions_grid, dtype=np.float64),
        np.ascontiguousarray(active, dtype=np.bool_),
        velocity_xz_um_s,
        velocity_gradient_s_inv,
        dynamic_viscosity_pa_s,
        local_vessel_radius_um,
        wall_shear_stress_pa,
    )


def _sample_regular_grid_fields_kernel(
    positions,
    active,
    velocity,
    velocity_gradient,
    viscosity,
    vessel_radius,
    wall_shear,
):
    count = positions.shape[0]
    nx = viscosity.shape[0]
    nz = viscosity.shape[1]
    sampled_velocity = np.zeros((count, 2), dtype=np.float64)
    sampled_gradient = np.zeros((count, 2, 2), dtype=np.float64)
    sampled_viscosity = np.zeros(count, dtype=np.float64)
    sampled_radius = np.zeros(count, dtype=np.float64)
    sampled_shear = np.zeros(count, dtype=np.float64)

    for lane in range(count):
        if not active[lane]:
            continue

        gx = positions[lane, 0]
        gz = positions[lane, 1]
        if gx < 0.0:
            gx = 0.0
        elif gx > nx - 1.0:
            gx = nx - 1.0
        if gz < 0.0:
            gz = 0.0
        elif gz > nz - 1.0:
            gz = nz - 1.0

        i0 = int(np.floor(gx))
        j0 = int(np.floor(gz))
        if i0 > nx - 2:
            i0 = nx - 2
        if j0 > nz - 2:
            j0 = nz - 2
        i1 = i0 + 1
        j1 = j0 + 1
        tx = gx - i0
        tz = gz - j0
        w00 = (1.0 - tx) * (1.0 - tz)
        w10 = tx * (1.0 - tz)
        w01 = (1.0 - tx) * tz
        w11 = tx * tz

        for component in range(2):
            sampled_velocity[lane, component] = (
                w00 * velocity[i0, j0, component]
                + w10 * velocity[i1, j0, component]
                + w01 * velocity[i0, j1, component]
                + w11 * velocity[i1, j1, component]
            )
            for derivative in range(2):
                sampled_gradient[lane, component, derivative] = (
                    w00 * velocity_gradient[i0, j0, component, derivative]
                    + w10 * velocity_gradient[i1, j0, component, derivative]
                    + w01 * velocity_gradient[i0, j1, component, derivative]
                    + w11 * velocity_gradient[i1, j1, component, derivative]
                )

        sampled_viscosity[lane] = (
            w00 * viscosity[i0, j0]
            + w10 * viscosity[i1, j0]
            + w01 * viscosity[i0, j1]
            + w11 * viscosity[i1, j1]
        )
        sampled_radius[lane] = (
            w00 * vessel_radius[i0, j0]
            + w10 * vessel_radius[i1, j0]
            + w01 * vessel_radius[i0, j1]
            + w11 * vessel_radius[i1, j1]
        )
        sampled_shear[lane] = (
            w00 * wall_shear[i0, j0]
            + w10 * wall_shear[i1, j0]
            + w01 * wall_shear[i0, j1]
            + w11 * wall_shear[i1, j1]
        )

    return (
        sampled_velocity,
        sampled_gradient,
        sampled_viscosity,
        sampled_radius,
        sampled_shear,
    )


if njit is not None:
    _sample_regular_grid_fields_numba = njit(cache=True)(
        _sample_regular_grid_fields_kernel
    )
else:  # pragma: no cover - retained for deliberately minimal environments
    _sample_regular_grid_fields_numba = _sample_regular_grid_fields_kernel
