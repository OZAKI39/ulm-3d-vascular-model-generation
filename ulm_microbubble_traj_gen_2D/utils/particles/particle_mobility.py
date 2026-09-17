"""Near-wall mobility formulas for two-dimensional microbubble transport.

The formulas follow ``references/Revised v5.md`` and the supplementary
material for Borden et al. (2018), with the contact shear-force coefficient
restored to the Goldman et al. (1967) Part II value.  A bubble is treated as a
three-dimensional rigid sphere whose centre moves in the X-Z plane and whose
retained rotation is about +Y.  For an inward unit normal ``n = (n_x, n_z)``,
the tangent is always ``t = (-n_z, n_x)`` so that ``t cross n`` points along
+Y.

Lengths are expressed in micrometers, velocities in micrometers per second,
dynamic viscosity in pascal-seconds, forces in piconewtons, and torque in
piconewton-micrometers.  The scalar functions are allocation-free and are
JIT-compiled when Numba is available.
"""

from __future__ import annotations

import math
from typing import Callable

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - exercised only in minimal environments
    def njit(*_args: object, **_kwargs: object) -> Callable:
        def decorate(function: Callable) -> Callable:
            return function

        return decorate


DEFAULT_XI_MIN = 1.0e-3
DEFAULT_XI_NEAR = 0.1
DEFAULT_XI_FAR = 1.0
# Goldman, Cox, and Brenner (1967), Part II, contact limit H/R = 1.
# The downstream Borden supplementary material transcribes F_s as 1.0075;
# Goldman's original extrapolated table value is F_s = 1.7005.
SHEAR_FORCE_COEFFICIENT = 1.7005
SHEAR_TORQUE_COEFFICIENT = 0.9440
FREE_ROTATIONAL_SCALED_MOBILITY = 0.75


@njit(cache=True, inline="always")
def tangent_from_inward_normal_scalar(normal_x: float, normal_z: float) -> tuple[float, float]:
    """Return ``t=(-n_z,n_x)``, for which ``t cross n`` is along +Y."""

    return -normal_z, normal_x


@njit(cache=True, inline="always")
def global_to_wall_components_scalar(
    vector_x: float,
    vector_z: float,
    normal_x: float,
    normal_z: float,
) -> tuple[float, float]:
    """Project an X-Z vector onto the local tangent and inward normal."""

    tangent_x, tangent_z = tangent_from_inward_normal_scalar(normal_x, normal_z)
    return (
        tangent_x * vector_x + tangent_z * vector_z,
        normal_x * vector_x + normal_z * vector_z,
    )


@njit(cache=True, inline="always")
def wall_to_global_components_scalar(
    tangent_component: float,
    normal_component: float,
    normal_x: float,
    normal_z: float,
) -> tuple[float, float]:
    """Transform local tangent-normal components back to the X-Z plane."""

    tangent_x, tangent_z = tangent_from_inward_normal_scalar(normal_x, normal_z)
    return (
        tangent_component * tangent_x + normal_component * normal_x,
        tangent_component * tangent_z + normal_component * normal_z,
    )


@njit(cache=True, inline="always")
def bulk_angular_velocity_scalar(du_x_dz: float, du_z_dx: float) -> float:
    """Return the +Y angular velocity, one half of the signed fluid vorticity."""

    return 0.5 * (du_x_dz - du_z_dx)


@njit(cache=True, inline="always")
def local_shear_rate_scalar(
    du_x_dx: float,
    du_x_dz: float,
    du_z_dx: float,
    du_z_dz: float,
    normal_x: float,
    normal_z: float,
) -> float:
    """Evaluate the signed near-wall shear rate ``t.T @ grad(u) @ n``."""

    grad_u_n_x = du_x_dx * normal_x + du_x_dz * normal_z
    grad_u_n_z = du_z_dx * normal_x + du_z_dz * normal_z
    tangent_x, tangent_z = tangent_from_inward_normal_scalar(normal_x, normal_z)
    return tangent_x * grad_u_n_x + tangent_z * grad_u_n_z


@njit(cache=True, inline="always")
def wall_blend_weight_scalar(
    gap_ratio: float,
    xi_near: float = DEFAULT_XI_NEAR,
    xi_far: float = DEFAULT_XI_FAR,
) -> float:
    """Return the cubic smoothstep weight of the near-wall model."""

    s = (gap_ratio - xi_near) / (xi_far - xi_near)
    if s <= 0.0:
        return 1.0
    if s >= 1.0:
        return 0.0
    return 1.0 - 3.0 * s * s + 2.0 * s * s * s


@njit(cache=True, inline="always")
def wall_resistance_coefficients_scalar(
    model_gap_ratio: float,
) -> tuple[float, float, float, float, float, float]:
    """Return ``C_Ft, C_Tt, C_Fr, C_Tr, Delta, Lambda`` at positive xi."""

    log_xi = math.log(model_gap_ratio)
    c_force_translation = (8.0 / 15.0) * log_xi - 0.9588
    c_torque_translation = (1.0 / 10.0) * log_xi - 0.1895
    c_force_rotation = (2.0 / 15.0) * log_xi - 0.2526
    c_torque_rotation = (2.0 / 5.0) * log_xi - 0.3187
    determinant = (
        c_torque_translation * c_force_rotation
        - c_force_translation * c_torque_rotation
    )
    normal_resistance = 1.0 + 1.0 / model_gap_ratio
    return (
        c_force_translation,
        c_torque_translation,
        c_force_rotation,
        c_torque_rotation,
        determinant,
        normal_resistance,
    )


@njit(cache=True, inline="always")
def wall_dimensionless_mobility_entries_scalar(
    gap_ratio: float,
    xi_min: float = DEFAULT_XI_MIN,
) -> tuple[float, float, float, float, float]:
    """Return the five non-zero entries of the scaled wall mobility matrix.

    The tuple order is ``(M_tt, M_nn, M_tR, M_Rt, M_RR)`` for the scaled
    state ``(V_t, V_n, R*Omega)`` and load ``(F_t, F_n, T/R)``.  The two
    translation-rotation entries retain their distinct published formulas and
    are deliberately not symmetrized.
    """

    model_gap_ratio = gap_ratio
    if model_gap_ratio < xi_min:
        model_gap_ratio = xi_min
    c_ft, c_tt, c_fr, c_tr, determinant, normal_resistance = (
        wall_resistance_coefficients_scalar(model_gap_ratio)
    )
    return (
        c_tr / determinant,
        1.0 / normal_resistance,
        3.0 * c_fr / (4.0 * determinant),
        c_tt / determinant,
        3.0 * c_ft / (4.0 * determinant),
    )


@njit(cache=True, inline="always")
def effective_dimensionless_mobility_entries_scalar(
    gap_ratio: float,
    xi_min: float = DEFAULT_XI_MIN,
    xi_near: float = DEFAULT_XI_NEAR,
    xi_far: float = DEFAULT_XI_FAR,
) -> tuple[float, float, float, float, float]:
    """Blend the frozen near-wall asymptote into free-space mobility."""

    model_gap_ratio = gap_ratio
    if model_gap_ratio < xi_min:
        model_gap_ratio = xi_min
    if model_gap_ratio > xi_near:
        model_gap_ratio = xi_near
    m_tt_wall, m_nn_wall, m_t_r_wall, m_r_t_wall, m_r_r_wall = (
        wall_dimensionless_mobility_entries_scalar(model_gap_ratio, xi_min)
    )
    weight = wall_blend_weight_scalar(gap_ratio, xi_near, xi_far)
    one_minus_weight = 1.0 - weight
    return (
        one_minus_weight + weight * m_tt_wall,
        one_minus_weight + weight * m_nn_wall,
        weight * m_t_r_wall,
        weight * m_r_t_wall,
        one_minus_weight * FREE_ROTATIONAL_SCALED_MOBILITY + weight * m_r_r_wall,
    )


@njit(cache=True, inline="always")
def translational_stokes_mobility_scalar(
    viscosity_pa_s: float,
    radius_um: float,
) -> float:
    """Return free translation mobility in ``(um/s)/pN``."""

    return 1.0 / (6.0 * math.pi * viscosity_pa_s * radius_um)


@njit(cache=True, inline="always")
def wall_shear_load_scalar(
    viscosity_pa_s: float,
    radius_um: float,
    gap_um: float,
    shear_rate_per_s: float,
    force_coefficient: float = SHEAR_FORCE_COEFFICIENT,
    torque_coefficient: float = SHEAR_TORQUE_COEFFICIENT,
) -> tuple[float, float, float]:
    """Return tangential shear force, torque, and sphere-centre wall height.

    ``gap_um`` is the surface-to-wall clearance used by the mobility formula.
    The load formula instead uses the sphere-centre height
    ``H_center = radius + max(gap, 0)``.  Keeping these quantities separate
    prevents the shear force from incorrectly vanishing at wall contact.
    """

    model_gap_um = gap_um
    if model_gap_um < 0.0:
        model_gap_um = 0.0
    center_height_um = radius_um + model_gap_um
    force_t_pn = (
        force_coefficient
        * 6.0
        * math.pi
        * viscosity_pa_s
        * radius_um
        * center_height_um
        * shear_rate_per_s
    )
    torque_y_pn_um = (
        torque_coefficient
        * 4.0
        * math.pi
        * viscosity_pa_s
        * radius_um
        * radius_um
        * radius_um
        * shear_rate_per_s
    )
    return force_t_pn, torque_y_pn_um, center_height_um


@njit(cache=True, inline="always")
def apply_scaled_mobility_scalar(
    viscosity_pa_s: float,
    radius_um: float,
    force_t_pn: float,
    force_n_pn: float,
    torque_y_pn_um: float,
    m_tt: float,
    m_nn: float,
    m_t_r: float,
    m_r_t: float,
    m_r_r: float,
) -> tuple[float, float, float]:
    """Map a local generalized load to ``V_t, V_n, Omega_y``."""

    base_mobility = translational_stokes_mobility_scalar(viscosity_pa_s, radius_um)
    scaled_torque_load_pn = torque_y_pn_um / radius_um
    velocity_t_um_s = base_mobility * (
        m_tt * force_t_pn + m_t_r * scaled_torque_load_pn
    )
    velocity_n_um_s = base_mobility * m_nn * force_n_pn
    radius_angular_velocity_um_s = base_mobility * (
        m_r_t * force_t_pn + m_r_r * scaled_torque_load_pn
    )
    return (
        velocity_t_um_s,
        velocity_n_um_s,
        radius_angular_velocity_um_s / radius_um,
    )


def generalized_mobility_matrix_xz(
    viscosity_pa_s: np.ndarray,
    radius_um: np.ndarray,
    inward_normal_xz: np.ndarray,
    dimensionless_entries: np.ndarray,
) -> np.ndarray:
    """Build physical ``(Vx,Vz,Omega) <- (Fx,Fz,T)`` mobility matrices.

    The published near-wall entries use the scaled local variables
    ``(V_t, V_n, R*Omega)`` and ``(F_t, F_n, T/R)``.  Contact must act through
    the same mobility relation as collisions and molecular loads, so this
    helper performs the radius scaling and rotates the local tangent/normal
    matrix into the global X-Z basis once per RHS evaluation.
    """

    viscosity = np.asarray(viscosity_pa_s, dtype=np.float64)
    radius = np.asarray(radius_um, dtype=np.float64)
    normal = np.asarray(inward_normal_xz, dtype=np.float64)
    entries = np.asarray(dimensionless_entries, dtype=np.float64)
    if viscosity.ndim != 1 or radius.shape != viscosity.shape:
        raise ValueError("viscosity_pa_s and radius_um must be matching one-dimensional arrays.")
    if normal.shape != (viscosity.size, 2):
        raise ValueError("inward_normal_xz must have shape (particle_count, 2).")
    if entries.shape != (viscosity.size, 5):
        raise ValueError("dimensionless_entries must have shape (particle_count, 5).")
    if (
        not np.all(np.isfinite(viscosity))
        or not np.all(np.isfinite(radius))
        or not np.all(np.isfinite(normal))
        or not np.all(np.isfinite(entries))
        or np.any(viscosity <= 0.0)
        or np.any(radius <= 0.0)
    ):
        raise ValueError("Generalized mobility inputs must be finite with positive viscosity and radius.")

    count = int(viscosity.size)
    matrices = np.zeros((count, 3, 3), dtype=np.float64)
    if count == 0:
        return matrices

    normal_norm = np.linalg.norm(normal, axis=1)
    if np.any(~np.isfinite(normal_norm)) or np.any(normal_norm <= 0.0):
        raise ValueError("Every inward wall normal must have a finite, non-zero length.")
    unit_normal = normal / normal_norm[:, None]
    tangent = np.column_stack((-unit_normal[:, 1], unit_normal[:, 0]))
    base = 1.0 / (6.0 * np.pi * viscosity * radius)

    m_tt = entries[:, 0]
    m_nn = entries[:, 1]
    m_t_r = entries[:, 2]
    m_r_t = entries[:, 3]
    m_r_r = entries[:, 4]

    matrices[:, :2, :2] = base[:, None, None] * (
        m_tt[:, None, None] * tangent[:, :, None] * tangent[:, None, :]
        + m_nn[:, None, None] * unit_normal[:, :, None] * unit_normal[:, None, :]
    )
    matrices[:, :2, 2] = (
        base * m_t_r / radius
    )[:, None] * tangent
    matrices[:, 2, :2] = (
        base * m_r_t / radius
    )[:, None] * tangent
    matrices[:, 2, 2] = base * m_r_r / (radius * radius)
    return matrices


@njit(cache=True, inline="always")
def bulk_background_velocity_local_scalar(
    velocity_x_um_s: float,
    velocity_z_um_s: float,
    du_x_dz_per_s: float,
    du_z_dx_per_s: float,
    normal_x: float,
    normal_z: float,
    radius_um: float,
) -> tuple[float, float, float]:
    """Return the scaled local bulk state ``(V_t,V_n,R*Omega_y)``."""

    velocity_t, velocity_n = global_to_wall_components_scalar(
        velocity_x_um_s,
        velocity_z_um_s,
        normal_x,
        normal_z,
    )
    angular_velocity = bulk_angular_velocity_scalar(du_x_dz_per_s, du_z_dx_per_s)
    return velocity_t, velocity_n, radius_um * angular_velocity


@njit(cache=True, inline="always")
def wall_shear_velocity_local_scalar(
    viscosity_pa_s: float,
    radius_um: float,
    gap_um: float,
    shear_rate_per_s: float,
    xi_min: float = DEFAULT_XI_MIN,
    xi_near: float = DEFAULT_XI_NEAR,
) -> tuple[float, float, float]:
    """Return the scaled wall-shear state ``(V_t,V_n,R*Omega_y)``."""

    gap_ratio = gap_um / radius_um
    model_gap_ratio = gap_ratio
    if model_gap_ratio < xi_min:
        model_gap_ratio = xi_min
    if model_gap_ratio > xi_near:
        model_gap_ratio = xi_near
    entries = wall_dimensionless_mobility_entries_scalar(model_gap_ratio, xi_min)
    force_t_pn, torque_y_pn_um, _ = wall_shear_load_scalar(
        viscosity_pa_s,
        radius_um,
        gap_um,
        shear_rate_per_s,
    )
    velocity_t, velocity_n, angular_velocity = apply_scaled_mobility_scalar(
        viscosity_pa_s,
        radius_um,
        force_t_pn,
        0.0,
        torque_y_pn_um,
        entries[0],
        entries[1],
        entries[2],
        entries[3],
        entries[4],
    )
    return velocity_t, velocity_n, radius_um * angular_velocity


@njit(cache=True, inline="always")
def background_hydrodynamic_velocity_scalar(
    velocity_x_um_s: float,
    velocity_z_um_s: float,
    du_x_dx_per_s: float,
    du_x_dz_per_s: float,
    du_z_dx_per_s: float,
    du_z_dz_per_s: float,
    normal_x: float,
    normal_z: float,
    radius_um: float,
    gap_um: float,
    viscosity_pa_s: float,
    xi_min: float = DEFAULT_XI_MIN,
    xi_near: float = DEFAULT_XI_NEAR,
    xi_far: float = DEFAULT_XI_FAR,
) -> tuple[float, float, float, float, float]:
    """Return global background velocity, rotation, wall weight, and shear.

    The first three values are ``V_x``, ``V_z`` and ``Omega_y``.  The final
    two diagnostic values are the near-wall blend weight and signed local
    shear rate.  The inward normal is expected to be a unit vector.
    """

    gap_ratio = gap_um / radius_um
    weight = wall_blend_weight_scalar(gap_ratio, xi_near, xi_far)
    bulk_omega = bulk_angular_velocity_scalar(du_x_dz_per_s, du_z_dx_per_s)
    if weight <= 0.0:
        shear_rate = local_shear_rate_scalar(
            du_x_dx_per_s,
            du_x_dz_per_s,
            du_z_dx_per_s,
            du_z_dz_per_s,
            normal_x,
            normal_z,
        )
        return velocity_x_um_s, velocity_z_um_s, bulk_omega, 0.0, shear_rate
    bulk_t, bulk_n, bulk_r_omega = bulk_background_velocity_local_scalar(
        velocity_x_um_s,
        velocity_z_um_s,
        du_x_dz_per_s,
        du_z_dx_per_s,
        normal_x,
        normal_z,
        radius_um,
    )
    shear_rate = local_shear_rate_scalar(
        du_x_dx_per_s,
        du_x_dz_per_s,
        du_z_dx_per_s,
        du_z_dz_per_s,
        normal_x,
        normal_z,
    )
    wall_t, wall_n, wall_r_omega = wall_shear_velocity_local_scalar(
        viscosity_pa_s,
        radius_um,
        gap_um,
        shear_rate,
        xi_min,
        xi_near,
    )
    one_minus_weight = 1.0 - weight
    velocity_t = one_minus_weight * bulk_t + weight * wall_t
    velocity_n = one_minus_weight * bulk_n + weight * wall_n
    radius_angular_velocity = (
        one_minus_weight * bulk_r_omega + weight * wall_r_omega
    )
    velocity_x, velocity_z = wall_to_global_components_scalar(
        velocity_t,
        velocity_n,
        normal_x,
        normal_z,
    )
    return (
        velocity_x,
        velocity_z,
        radius_angular_velocity / radius_um,
        weight,
        shear_rate,
    )
