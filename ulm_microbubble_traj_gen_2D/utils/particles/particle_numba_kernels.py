"""Allocation-conscious batch kernels for the particle mobility hot path.

The public transport code keeps lifecycle, event timing, and named diagnostic
objects in Python.  This module only fuses repeated numeric work that used to
cross the Python/Numba boundary once per bubble.  Every permanent-ID and
finite-radius rule remains owned by the existing transport modules.
"""

from __future__ import annotations

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - the production environment includes Numba
    njit = None

from .particle_mobility import (
    FREE_ROTATIONAL_SCALED_MOBILITY,
    apply_scaled_mobility_scalar,
    background_hydrodynamic_velocity_scalar,
    bulk_angular_velocity_scalar,
    effective_dimensionless_mobility_entries_scalar,
    translational_stokes_mobility_scalar,
    wall_to_global_components_scalar,
)


def evaluate_hydrodynamic_mobility_batch(
    fluid_velocity_xz_um_s: np.ndarray,
    velocity_gradient_s_inv: np.ndarray,
    dynamic_viscosity_pa_s: np.ndarray,
    exact_wall_distance_um: np.ndarray,
    sampled_wall_normal_xz: np.ndarray,
    local_vessel_radius_um: np.ndarray,
    bubble_radii_um: np.ndarray,
    active: np.ndarray,
    *,
    near_wall_enabled: bool,
    xi_min: float,
    xi_near: float,
    xi_far: float,
    two_wall_warning_gap_ratio: float,
) -> tuple[np.ndarray, ...]:
    """Run the independent per-bubble mobility preparation in one JIT call."""

    result = _hydrodynamic_mobility_numba(
        np.ascontiguousarray(fluid_velocity_xz_um_s, dtype=np.float64),
        np.ascontiguousarray(velocity_gradient_s_inv, dtype=np.float64),
        np.ascontiguousarray(dynamic_viscosity_pa_s, dtype=np.float64),
        np.ascontiguousarray(exact_wall_distance_um, dtype=np.float64),
        np.ascontiguousarray(sampled_wall_normal_xz, dtype=np.float64),
        np.ascontiguousarray(local_vessel_radius_um, dtype=np.float64),
        np.ascontiguousarray(bubble_radii_um, dtype=np.float64),
        np.ascontiguousarray(active, dtype=np.bool_),
        bool(near_wall_enabled),
        float(xi_min),
        float(xi_near),
        float(xi_far),
        float(two_wall_warning_gap_ratio),
    )
    if bool(result[-1]):
        raise ValueError("Sampled particle viscosity must be finite and positive.")
    return result[:-1]


def apply_particle_loads_batch(
    particle_velocity_xz_um_s: np.ndarray,
    angular_velocity_rad_s: np.ndarray,
    dynamic_viscosity_pa_s: np.ndarray,
    bubble_radii_um: np.ndarray,
    wall_normal_xz: np.ndarray,
    mobility_entries: np.ndarray,
    collision_force_xz_pn: np.ndarray,
    active: np.ndarray,
    bond_force_tangent_pn: np.ndarray | None = None,
    bond_force_normal_pn: np.ndarray | None = None,
    bond_torque_pn_um: np.ndarray | None = None,
) -> np.ndarray:
    """Apply collision and optional molecular loads without per-lane dispatch."""

    count = int(np.asarray(active).size)
    zero = np.zeros(count, dtype=np.float64)
    return _apply_particle_loads_numba(
        np.ascontiguousarray(particle_velocity_xz_um_s, dtype=np.float64),
        np.ascontiguousarray(angular_velocity_rad_s, dtype=np.float64),
        np.ascontiguousarray(dynamic_viscosity_pa_s, dtype=np.float64),
        np.ascontiguousarray(bubble_radii_um, dtype=np.float64),
        np.ascontiguousarray(wall_normal_xz, dtype=np.float64),
        np.ascontiguousarray(mobility_entries, dtype=np.float64),
        np.ascontiguousarray(collision_force_xz_pn, dtype=np.float64),
        np.ascontiguousarray(active, dtype=np.bool_),
        zero if bond_force_tangent_pn is None else np.ascontiguousarray(bond_force_tangent_pn, dtype=np.float64),
        zero if bond_force_normal_pn is None else np.ascontiguousarray(bond_force_normal_pn, dtype=np.float64),
        zero if bond_torque_pn_um is None else np.ascontiguousarray(bond_torque_pn_um, dtype=np.float64),
    )


def build_generalized_mobility_batch(
    viscosity_pa_s: np.ndarray,
    radius_um: np.ndarray,
    inward_normal_xz: np.ndarray,
    dimensionless_entries: np.ndarray,
) -> np.ndarray:
    """Build contact mobility matrices without broadcast temporaries."""

    matrices, invalid = _generalized_mobility_numba(
        np.ascontiguousarray(viscosity_pa_s, dtype=np.float64),
        np.ascontiguousarray(radius_um, dtype=np.float64),
        np.ascontiguousarray(inward_normal_xz, dtype=np.float64),
        np.ascontiguousarray(dimensionless_entries, dtype=np.float64),
    )
    if invalid:
        raise ValueError(
            "Generalized mobility inputs must be finite with positive viscosity, "
            "radius, and wall-normal length."
        )
    return matrices


def _hydrodynamic_mobility_kernel(
    fluid_velocity,
    gradients,
    viscosity,
    wall_distance,
    sampled_normals,
    local_vessel_radius,
    radii,
    active,
    near_wall_enabled,
    xi_min,
    xi_near,
    xi_far,
    two_wall_warning_gap_ratio,
):
    count = active.size
    particle_velocity = np.zeros((count, 2), dtype=np.float64)
    angular_velocity = np.zeros(count, dtype=np.float64)
    normalized_normals = np.zeros((count, 2), dtype=np.float64)
    wall_gap = np.full(count, np.nan, dtype=np.float64)
    gap_ratio = np.full(count, np.nan, dtype=np.float64)
    wall_weight = np.zeros(count, dtype=np.float64)
    two_wall_warning = np.zeros(count, dtype=np.bool_)
    mobility_entries = np.zeros((count, 5), dtype=np.float64)
    translation_mobility_global = np.zeros((count, 2, 2), dtype=np.float64)
    reciprocity_error = np.zeros(count, dtype=np.float64)
    degenerate_normal = np.zeros(count, dtype=np.uint8)
    invalid_viscosity = False

    for lane in range(count):
        if not active[lane]:
            continue
        radius_um = radii[lane]
        mu_pa_s = viscosity[lane]
        if not np.isfinite(mu_pa_s) or mu_pa_s <= 0.0:
            invalid_viscosity = True
            continue

        gap_um = wall_distance[lane] - radius_um
        wall_gap[lane] = gap_um
        gap_ratio[lane] = gap_um / radius_um
        normal_x = sampled_normals[lane, 0]
        normal_z = sampled_normals[lane, 1]
        normal_norm = np.sqrt(normal_x * normal_x + normal_z * normal_z)
        scaled = (gap_ratio[lane] - xi_near) / (xi_far - xi_near)
        scaled = min(1.0, max(0.0, scaled))
        potential_wall_weight = 1.0 - 3.0 * scaled * scaled + 2.0 * scaled * scaled * scaled

        if normal_norm <= 1.0e-12 or not np.isfinite(normal_norm):
            normal_x = 0.0
            normal_z = 1.0
            if near_wall_enabled and potential_wall_weight > 0.0:
                degenerate_normal[lane] = 1
        else:
            normal_x /= normal_norm
            normal_z /= normal_norm
        normalized_normals[lane, 0] = normal_x
        normalized_normals[lane, 1] = normal_z

        if near_wall_enabled:
            vx, vz, omega, weight, _ = background_hydrodynamic_velocity_scalar(
                fluid_velocity[lane, 0],
                fluid_velocity[lane, 1],
                gradients[lane, 0, 0],
                gradients[lane, 0, 1],
                gradients[lane, 1, 0],
                gradients[lane, 1, 1],
                normal_x,
                normal_z,
                radius_um,
                max(gap_um, 0.0),
                mu_pa_s,
                xi_min,
                xi_near,
                xi_far,
            )
            entries = effective_dimensionless_mobility_entries_scalar(
                max(gap_um, 0.0) / radius_um,
                xi_min,
                xi_near,
                xi_far,
            )
        else:
            vx = fluid_velocity[lane, 0]
            vz = fluid_velocity[lane, 1]
            omega = bulk_angular_velocity_scalar(
                gradients[lane, 0, 1],
                gradients[lane, 1, 0],
            )
            weight = 0.0
            entries = (1.0, 1.0, 0.0, 0.0, FREE_ROTATIONAL_SCALED_MOBILITY)

        particle_velocity[lane, 0] = vx
        particle_velocity[lane, 1] = vz
        angular_velocity[lane] = omega
        wall_weight[lane] = weight
        for entry_index in range(5):
            mobility_entries[lane, entry_index] = entries[entry_index]

        denominator = max(abs(entries[2]), abs(entries[3]), 1.0e-30)
        reciprocity_error[lane] = abs(entries[2] - entries[3]) / denominator
        tangent_x = -normal_z
        tangent_z = normal_x
        base = translational_stokes_mobility_scalar(mu_pa_s, radius_um)
        m_tt = entries[0]
        m_nn = entries[1]
        translation_mobility_global[lane, 0, 0] = base * (
            m_tt * tangent_x * tangent_x + m_nn * normal_x * normal_x
        )
        translation_mobility_global[lane, 0, 1] = base * (
            m_tt * tangent_x * tangent_z + m_nn * normal_x * normal_z
        )
        translation_mobility_global[lane, 1, 0] = translation_mobility_global[lane, 0, 1]
        translation_mobility_global[lane, 1, 1] = base * (
            m_tt * tangent_z * tangent_z + m_nn * normal_z * normal_z
        )
        opposite_gap_um = 2.0 * local_vessel_radius[lane] - wall_distance[lane] - radius_um
        opposite_gap_ratio = opposite_gap_um / radius_um
        two_wall_warning[lane] = bool(
            near_wall_enabled
            and weight > 0.0
            and opposite_gap_ratio <= two_wall_warning_gap_ratio
        )

    return (
        particle_velocity,
        angular_velocity,
        normalized_normals,
        wall_gap,
        gap_ratio,
        wall_weight,
        two_wall_warning,
        mobility_entries,
        translation_mobility_global,
        reciprocity_error,
        degenerate_normal,
        invalid_viscosity,
    )


def _apply_particle_loads_kernel(
    particle_velocity,
    angular_velocity,
    viscosity,
    radii,
    normalized_normals,
    mobility_entries,
    collision_force,
    active,
    bond_force_tangent,
    bond_force_normal,
    bond_torque,
):
    count = active.size
    collision_speed = np.zeros(count, dtype=np.float64)
    for lane in range(count):
        if not active[lane]:
            continue
        normal_x = normalized_normals[lane, 0]
        normal_z = normalized_normals[lane, 1]
        tangent_x = -normal_z
        tangent_z = normal_x
        force_x = collision_force[lane, 0]
        force_z = collision_force[lane, 1]
        force_t = tangent_x * force_x + tangent_z * force_z
        force_n = normal_x * force_x + normal_z * force_z
        entries = mobility_entries[lane]
        delta_t, delta_n, delta_omega = apply_scaled_mobility_scalar(
            viscosity[lane],
            radii[lane],
            force_t + bond_force_tangent[lane],
            force_n + bond_force_normal[lane],
            bond_torque[lane],
            entries[0],
            entries[1],
            entries[2],
            entries[3],
            entries[4],
        )
        delta_x, delta_z = wall_to_global_components_scalar(
            delta_t,
            delta_n,
            normal_x,
            normal_z,
        )
        particle_velocity[lane, 0] += delta_x
        particle_velocity[lane, 1] += delta_z
        angular_velocity[lane] += delta_omega

        collision_delta_t, collision_delta_n, _ = apply_scaled_mobility_scalar(
            viscosity[lane],
            radii[lane],
            force_t,
            force_n,
            0.0,
            entries[0],
            entries[1],
            entries[2],
            entries[3],
            entries[4],
        )
        collision_speed[lane] = np.sqrt(
            collision_delta_t * collision_delta_t
            + collision_delta_n * collision_delta_n
        )
    return collision_speed


def _generalized_mobility_kernel(viscosity, radius, normal, entries):
    count = viscosity.size
    matrices = np.zeros((count, 3, 3), dtype=np.float64)
    invalid = False
    for lane in range(count):
        mu = viscosity[lane]
        particle_radius = radius[lane]
        normal_x = normal[lane, 0]
        normal_z = normal[lane, 1]
        normal_length = np.sqrt(normal_x * normal_x + normal_z * normal_z)
        if (
            not np.isfinite(mu)
            or not np.isfinite(particle_radius)
            or not np.isfinite(normal_length)
            or mu <= 0.0
            or particle_radius <= 0.0
            or normal_length <= 0.0
        ):
            invalid = True
            continue
        normal_x /= normal_length
        normal_z /= normal_length
        tangent_x = -normal_z
        tangent_z = normal_x
        base = 1.0 / (6.0 * np.pi * mu * particle_radius)
        m_tt = entries[lane, 0]
        m_nn = entries[lane, 1]
        m_t_r = entries[lane, 2]
        m_r_t = entries[lane, 3]
        m_r_r = entries[lane, 4]
        for value_index in range(5):
            if not np.isfinite(entries[lane, value_index]):
                invalid = True

        matrices[lane, 0, 0] = base * (
            m_tt * tangent_x * tangent_x + m_nn * normal_x * normal_x
        )
        matrices[lane, 0, 1] = base * (
            m_tt * tangent_x * tangent_z + m_nn * normal_x * normal_z
        )
        matrices[lane, 1, 0] = matrices[lane, 0, 1]
        matrices[lane, 1, 1] = base * (
            m_tt * tangent_z * tangent_z + m_nn * normal_z * normal_z
        )
        matrices[lane, 0, 2] = base * m_t_r * tangent_x / particle_radius
        matrices[lane, 1, 2] = base * m_t_r * tangent_z / particle_radius
        matrices[lane, 2, 0] = base * m_r_t * tangent_x / particle_radius
        matrices[lane, 2, 1] = base * m_r_t * tangent_z / particle_radius
        matrices[lane, 2, 2] = base * m_r_r / (particle_radius * particle_radius)
    return matrices, invalid


if njit is not None:
    _hydrodynamic_mobility_numba = njit(cache=True)(_hydrodynamic_mobility_kernel)
    _apply_particle_loads_numba = njit(cache=True)(_apply_particle_loads_kernel)
    _generalized_mobility_numba = njit(cache=True)(_generalized_mobility_kernel)
else:  # pragma: no cover - retained for deliberately minimal environments
    _hydrodynamic_mobility_numba = _hydrodynamic_mobility_kernel
    _apply_particle_loads_numba = _apply_particle_loads_kernel
    _generalized_mobility_numba = _generalized_mobility_kernel
