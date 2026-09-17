"""Deterministic mean-field molecular bonds for targeted microbubbles.

This module implements the numerical core specified by the final section of
``references/Revised v7.md``.  It deliberately knows nothing about vessel
topology, particle registries, or the mobility matrix.  A caller supplies the
target-positive reaction area and the instantaneous wall-frame slip velocity;
the returned local force, torque, and state derivatives can then be coupled to
the existing particle transport code.

External units follow the rest of the particle solver:

* lengths, gaps, and tangential extensions are in micrometers;
* velocities are in micrometers per second;
* force is in piconewtons and torque in piconewton-micrometers;
* surface densities are in molecules per square meter;
* the effective two-dimensional association rate is in
  square-meters per molecule per second;
* the Bell sensitivity length is in nanometers; and
* temperature is in kelvin.

The square-micrometer-to-square-meter and piconewton-nanometer-to-joule
conversions are centralized below.  The scalar and batch kernels never mutate
their inputs.  Numba is optional; the Python/NumPy fallback follows exactly the
same scalar implementation so it remains a useful correctness oracle.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable

import numpy as np

try:
    from numba import njit
    from numba.extending import register_jitable
except ImportError:  # pragma: no cover - exercised only in minimal environments
    njit = None

    def register_jitable(function: Callable) -> Callable:
        return function


BOLTZMANN_CONSTANT_J_PER_K = 1.380649e-23
SQUARE_METERS_PER_SQUARE_MICROMETER = 1.0e-12
JOULES_PER_PICONEWTON_NANOMETER = 1.0e-21
MAX_FLOAT64 = float(np.finfo(np.float64).max)
LOG_MAX_FLOAT64 = math.log(MAX_FLOAT64)


@dataclass(frozen=True)
class MolecularBindingParameters:
    """Physical parameters for one effective ligand-target molecular pair."""

    ligand_density_molecules_m2: float
    target_density_molecules_m2: float
    association_rate_m2_per_molecule_s: float
    zero_force_dissociation_rate_s: float
    rest_length_um: float
    spring_stiffness_pn_per_um: float
    force_sensitivity_length_nm: float
    temperature_k: float
    bell_exponent_limit: float = 80.0

    def __post_init__(self) -> None:
        nonnegative = {
            "ligand_density_molecules_m2": self.ligand_density_molecules_m2,
            "target_density_molecules_m2": self.target_density_molecules_m2,
            "association_rate_m2_per_molecule_s": (
                self.association_rate_m2_per_molecule_s
            ),
            "zero_force_dissociation_rate_s": (
                self.zero_force_dissociation_rate_s
            ),
            "rest_length_um": self.rest_length_um,
            "spring_stiffness_pn_per_um": self.spring_stiffness_pn_per_um,
            "force_sensitivity_length_nm": self.force_sensitivity_length_nm,
        }
        for name, value in nonnegative.items():
            if not math.isfinite(float(value)) or float(value) < 0.0:
                raise ValueError(f"{name} must be finite and non-negative.")
        if not math.isfinite(float(self.temperature_k)) or self.temperature_k <= 0.0:
            raise ValueError("temperature_k must be finite and greater than zero.")
        if (
            not math.isfinite(float(self.bell_exponent_limit))
            or self.bell_exponent_limit <= 0.0
        ):
            raise ValueError("bell_exponent_limit must be finite and greater than zero.")


@dataclass(frozen=True)
class MolecularBindingEvaluation:
    """Immutable scalar or broadcast-array evaluation of one bond RHS stage."""

    expected_bond_count: np.ndarray
    total_tangential_extension_um: np.ndarray
    mean_tangential_extension_um: np.ndarray
    reaction_area_um2: np.ndarray
    ligand_count: np.ndarray
    target_count: np.ndarray
    formation_capacity: np.ndarray
    formation_rate_bonds_s: np.ndarray
    dissociation_rate_s: np.ndarray
    slip_velocity_um_s: np.ndarray
    extension_source_um_s: np.ndarray
    single_bond_tension_pn: np.ndarray
    force_t_pn: np.ndarray
    force_n_pn: np.ndarray
    torque_y_pn_um: np.ndarray
    bell_exponent: np.ndarray
    bell_rate_saturated: np.ndarray
    formation_rate_saturated: np.ndarray


@dataclass(frozen=True)
class BondStateUpdate:
    """A newly allocated, canonical bond state from an accepted-step helper."""

    expected_bond_count: np.ndarray
    total_tangential_extension_um: np.ndarray
    capacity_limited: np.ndarray


def reaction_disk_radius_um(
    radius_um: float | np.ndarray,
    gap_um: float | np.ndarray,
    capture_distance_um: float | np.ndarray,
) -> np.ndarray:
    """Return ``sqrt(2 R delta - delta**2)`` for the capture projection disk.

    A round-off-sized negative gap at an accepted wall contact is interpreted
    as zero.  A materially negative gap is rejected so an invalid penetrated
    state cannot silently acquire a molecular reaction area.  The transport
    geometry remains responsible for preventing wall penetration. ``capture_distance``
    may not exceed ``2*radius`` because the local spherical-cap formula would
    otherwise leave its physical branch.
    """

    radius, gap, capture = np.broadcast_arrays(
        np.asarray(radius_um, dtype=np.float64),
        np.asarray(gap_um, dtype=np.float64),
        np.asarray(capture_distance_um, dtype=np.float64),
    )
    if np.any(~np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("radius_um must contain only finite positive values.")
    if np.any(~np.isfinite(gap)):
        raise ValueError("gap_um must contain only finite values.")
    if np.any(~np.isfinite(capture)) or np.any(capture < 0.0):
        raise ValueError(
            "capture_distance_um must contain only finite non-negative values."
        )
    if np.any(capture > 2.0 * radius):
        raise ValueError("capture_distance_um must not exceed twice radius_um.")

    roundoff_tolerance = 256.0 * np.finfo(np.float64).eps * np.maximum.reduce(
        (np.ones_like(radius), radius, capture)
    )
    if np.any(gap < -roundoff_tolerance):
        raise ValueError(
            "gap_um contains a materially negative true wall gap; molecular "
            "capture cannot repair an invalid particle position."
        )
    physical_gap = np.maximum(gap, 0.0)
    delta = np.maximum(capture - physical_gap, 0.0)
    radicand = np.maximum(2.0 * radius * delta - delta * delta, 0.0)
    return np.sqrt(radicand)


def target_positive_reaction_area_um2(
    reaction_radius_um: float,
    target_interval_starts_um: np.ndarray,
    target_interval_ends_um: np.ndarray,
) -> float:
    r"""Evaluate the Revised-v7 target-positive line integral exactly.

    The supplied intervals are target-positive portions of the local wall,
    expressed as signed arc-length offsets ``eta`` from the closest wall point.
    For an indicator function that is one on their union, this function returns

    .. math::

       2\int_{-a}^{a}\chi_T(s+\eta)\sqrt{a^2-\eta^2}\,d\eta.

    Overlapping intervals are merged on local copies, so they cannot double
    count area and the caller's arrays are never modified.  The analytic disk
    antiderivative makes a fully positive interval exactly ``pi*a**2`` and
    avoids quadrature noise as a bubble crosses a target boundary.
    """

    radius = float(reaction_radius_um)
    if not math.isfinite(radius) or radius < 0.0:
        raise ValueError("reaction_radius_um must be finite and non-negative.")
    starts = np.asarray(target_interval_starts_um, dtype=np.float64)
    ends = np.asarray(target_interval_ends_um, dtype=np.float64)
    if starts.ndim != 1 or ends.ndim != 1 or starts.shape != ends.shape:
        raise ValueError("Target interval starts and ends must be matching 1D arrays.")
    if np.any(~np.isfinite(starts)) or np.any(~np.isfinite(ends)):
        raise ValueError("Target intervals must contain only finite values.")
    if np.any(ends < starts):
        raise ValueError("Every target interval end must be >= its start.")
    if radius == 0.0 or starts.size == 0:
        return 0.0

    clipped_starts = np.maximum(starts, -radius)
    clipped_ends = np.minimum(ends, radius)
    intersects = clipped_ends > clipped_starts
    if not np.any(intersects):
        return 0.0
    clipped_starts = clipped_starts[intersects]
    clipped_ends = clipped_ends[intersects]
    order = np.argsort(clipped_starts, kind="mergesort")
    sorted_starts = clipped_starts[order]
    sorted_ends = clipped_ends[order]

    area = 0.0
    merged_start = float(sorted_starts[0])
    merged_end = float(sorted_ends[0])
    for index in range(1, sorted_starts.size):
        start = float(sorted_starts[index])
        end = float(sorted_ends[index])
        if start <= merged_end:
            merged_end = max(merged_end, end)
            continue
        area += _target_interval_area_scalar(radius, merged_start, merged_end)
        merged_start = start
        merged_end = end
    area += _target_interval_area_scalar(radius, merged_start, merged_end)
    return min(max(float(area), 0.0), math.pi * radius * radius)


def surface_slip_velocity_um_s(
    tangential_velocity_um_s: float | np.ndarray,
    radius_um: float | np.ndarray,
    angular_velocity_y_s: float | np.ndarray,
) -> np.ndarray:
    """Return ``V_t + R*Omega_y`` under the project's +Y sign convention."""

    tangential_velocity, radius, angular_velocity = np.broadcast_arrays(
        np.asarray(tangential_velocity_um_s, dtype=np.float64),
        np.asarray(radius_um, dtype=np.float64),
        np.asarray(angular_velocity_y_s, dtype=np.float64),
    )
    if np.any(~np.isfinite(tangential_velocity)):
        raise ValueError("tangential_velocity_um_s must contain only finite values.")
    if np.any(~np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("radius_um must contain only finite positive values.")
    if np.any(~np.isfinite(angular_velocity)):
        raise ValueError("angular_velocity_y_s must contain only finite values.")
    return tangential_velocity + radius * angular_velocity


def evaluate_mean_field_bonds(
    expected_bond_count: float | np.ndarray,
    total_tangential_extension_um: float | np.ndarray,
    gap_um: float | np.ndarray,
    radius_um: float | np.ndarray,
    slip_velocity_um_s: float | np.ndarray,
    reaction_area_um2: float | np.ndarray,
    parameters: MolecularBindingParameters,
    *,
    use_numba: bool = True,
) -> MolecularBindingEvaluation:
    """Evaluate one immutable deterministic bond RHS stage.

    Existing bonds are not clipped when the reaction area shrinks.  If the
    area is zero, only formation stops; the returned Bell dissociation terms
    continue to decay ``n`` and ``m``.  When an inconsistent input has
    ``n == 0`` but nonzero ``m``, the evaluation uses the canonical Revised-v7
    state ``n=m=q=0`` without modifying the caller's storage.
    """

    if not isinstance(parameters, MolecularBindingParameters):
        raise TypeError("parameters must be a MolecularBindingParameters instance.")
    arrays = np.broadcast_arrays(
        np.asarray(expected_bond_count, dtype=np.float64),
        np.asarray(total_tangential_extension_um, dtype=np.float64),
        np.asarray(gap_um, dtype=np.float64),
        np.asarray(radius_um, dtype=np.float64),
        np.asarray(slip_velocity_um_s, dtype=np.float64),
        np.asarray(reaction_area_um2, dtype=np.float64),
    )
    n, m, gap, radius, slip, area = arrays
    _validate_state_inputs(n, m, gap, radius, slip, area)
    shape = n.shape
    flat = tuple(np.ascontiguousarray(value.reshape(-1)) for value in arrays)
    kernel = (
        _evaluate_mean_field_batch_numba
        if use_numba and njit is not None
        else _evaluate_mean_field_batch_kernel
    )
    values = kernel(
        *flat,
        float(parameters.ligand_density_molecules_m2),
        float(parameters.target_density_molecules_m2),
        float(parameters.association_rate_m2_per_molecule_s),
        float(parameters.zero_force_dissociation_rate_s),
        float(parameters.rest_length_um),
        float(parameters.spring_stiffness_pn_per_um),
        float(parameters.force_sensitivity_length_nm),
        float(parameters.temperature_k),
        float(parameters.bell_exponent_limit),
    )
    reshaped = tuple(np.asarray(value).reshape(shape) for value in values)
    return MolecularBindingEvaluation(*reshaped)


def predict_bond_state_exponential_euler(
    expected_bond_count: float | np.ndarray,
    total_tangential_extension_um: float | np.ndarray,
    rhs: MolecularBindingEvaluation,
    dt_s: float,
    *,
    extension_source_um_s: float | np.ndarray | None = None,
    use_numba: bool = True,
) -> BondStateUpdate:
    """Build a positive first-stage predictor without committing any state.

    Formation and extension sources and the Bell off-rate are frozen at the
    current RHS stage.  The linear dissociation term is integrated exactly,
    which remains non-negative even when ``k_off*dt`` is very large.  A purely
    numerical formation overshoot is limited to ``max(n_old, capacity)``;
    consequently a shrinking contact area never deletes pre-existing bonds.
    """

    extension_source = (
        rhs.extension_source_um_s
        if extension_source_um_s is None
        else np.asarray(extension_source_um_s, dtype=np.float64)
    )
    return _advance_bond_state(
        expected_bond_count,
        total_tangential_extension_um,
        rhs.formation_rate_bonds_s,
        rhs.dissociation_rate_s,
        extension_source,
        rhs.formation_capacity,
        dt_s,
        use_numba=use_numba,
    )


def accept_bond_state_exponential_heun(
    expected_bond_count: float | np.ndarray,
    total_tangential_extension_um: float | np.ndarray,
    first_rhs: MolecularBindingEvaluation,
    predictor_rhs: MolecularBindingEvaluation,
    dt_s: float,
    *,
    extension_source_um_s: float | np.ndarray | None = None,
    use_numba: bool = True,
) -> BondStateUpdate:
    """Return a positive accepted state from two Heun-style RHS stages.

    The source and off-rate are trapezoid-averaged, after which the averaged
    linear sink is integrated exponentially.  The second RHS must have been
    evaluated at the complete predicted particle state (including predicted
    position, velocity, ``n``, and ``m``); this function never evaluates or
    mutates that state itself.
    """

    formation_rate = _finite_half_sum(
        first_rhs.formation_rate_bonds_s,
        predictor_rhs.formation_rate_bonds_s,
    )
    off_rate = _finite_half_sum(
        first_rhs.dissociation_rate_s,
        predictor_rhs.dissociation_rate_s,
    )
    extension_source = (
        _finite_half_sum(
            first_rhs.extension_source_um_s,
            predictor_rhs.extension_source_um_s,
        )
        if extension_source_um_s is None
        else np.asarray(extension_source_um_s, dtype=np.float64)
    )
    capacity = np.maximum(first_rhs.formation_capacity, predictor_rhs.formation_capacity)
    return _advance_bond_state(
        expected_bond_count,
        total_tangential_extension_um,
        formation_rate,
        off_rate,
        extension_source,
        capacity,
        dt_s,
        use_numba=use_numba,
    )


def _validate_state_inputs(
    n: np.ndarray,
    m: np.ndarray,
    gap: np.ndarray,
    radius: np.ndarray,
    slip: np.ndarray,
    area: np.ndarray,
) -> None:
    if np.any(~np.isfinite(n)) or np.any(n < 0.0):
        raise ValueError("expected_bond_count must contain finite non-negative values.")
    if np.any(~np.isfinite(m)):
        raise ValueError("total_tangential_extension_um must contain only finite values.")
    if np.any(~np.isfinite(gap)):
        raise ValueError("gap_um must contain only finite values.")
    if np.any(~np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("radius_um must contain only finite positive values.")
    if np.any(~np.isfinite(slip)):
        raise ValueError("slip_velocity_um_s must contain only finite values.")
    if np.any(~np.isfinite(area)) or np.any(area < 0.0):
        raise ValueError("reaction_area_um2 must contain finite non-negative values.")


@register_jitable
def _target_disk_antiderivative_scalar(radius: float, eta: float) -> float:
    clipped = min(max(eta, -radius), radius)
    root = math.sqrt(max(radius * radius - clipped * clipped, 0.0))
    return 0.5 * (
        clipped * root + radius * radius * math.asin(clipped / radius)
    )


@register_jitable
def _target_interval_area_scalar(radius: float, start: float, end: float) -> float:
    return 2.0 * (
        _target_disk_antiderivative_scalar(radius, end)
        - _target_disk_antiderivative_scalar(radius, start)
    )


@register_jitable
def _safe_nonnegative_product(a: float, b: float) -> float:
    if a <= 0.0 or b <= 0.0:
        return 0.0
    # If either factor is <= 1 the product cannot exceed the other finite
    # factor.  This guard also avoids an overflowing MAX_FLOAT64/tiny division.
    if a > 1.0 and b > 1.0 and a > MAX_FLOAT64 / b:
        return MAX_FLOAT64
    return a * b


@register_jitable
def _safe_signed_product(a: float, b: float) -> float:
    if a == 0.0 or b == 0.0:
        return 0.0
    magnitude = _safe_nonnegative_product(abs(a), abs(b))
    return math.copysign(magnitude, a * math.copysign(1.0, b))


@register_jitable
def _safe_signed_sum(a: float, b: float) -> float:
    if a > 0.0 and b > 0.0 and a > MAX_FLOAT64 - b:
        return MAX_FLOAT64
    if a < 0.0 and b < 0.0 and a < -MAX_FLOAT64 - b:
        return -MAX_FLOAT64
    return a + b


@register_jitable
def _safe_count(density_molecules_m2: float, area_m2: float) -> float:
    return _safe_nonnegative_product(density_molecules_m2, area_m2)


@register_jitable
def _formation_rate_scalar(
    association_rate: float,
    area_m2: float,
    available_ligands: float,
    available_targets: float,
) -> tuple[float, bool]:
    if (
        association_rate <= 0.0
        or area_m2 <= 0.0
        or available_ligands <= 0.0
        or available_targets <= 0.0
    ):
        return 0.0, False
    log_rate = (
        math.log(association_rate)
        + math.log(available_ligands)
        + math.log(available_targets)
        - math.log(area_m2)
    )
    if log_rate >= LOG_MAX_FLOAT64:
        return MAX_FLOAT64, True
    return math.exp(log_rate), False


@register_jitable
def _bell_rate_scalar(
    tension_pn: float,
    zero_force_rate: float,
    sensitivity_nm: float,
    temperature_k: float,
    exponent_limit: float,
) -> tuple[float, float, bool]:
    if tension_pn <= 0.0 or sensitivity_nm <= 0.0:
        return zero_force_rate, 0.0, False

    log_exponent = (
        math.log(tension_pn)
        + math.log(sensitivity_nm)
        + math.log(JOULES_PER_PICONEWTON_NANOMETER)
        - math.log(BOLTZMANN_CONSTANT_J_PER_K)
        - math.log(temperature_k)
    )
    if log_exponent >= math.log(exponent_limit):
        exponent = exponent_limit
        exponent_was_limited = True
    else:
        exponent = math.exp(log_exponent)
        exponent_was_limited = False
    if zero_force_rate <= 0.0:
        return 0.0, exponent, exponent_was_limited
    maximum_exponent = LOG_MAX_FLOAT64 - math.log(zero_force_rate)
    if exponent >= maximum_exponent:
        return MAX_FLOAT64, exponent, True
    return (
        zero_force_rate * math.exp(exponent),
        exponent,
        exponent_was_limited,
    )


@register_jitable
def _evaluate_mean_field_scalar(
    n_input: float,
    m_input: float,
    gap_input: float,
    radius: float,
    slip: float,
    area_um2: float,
    ligand_density: float,
    target_density: float,
    association_rate: float,
    zero_force_off_rate: float,
    rest_length: float,
    spring_stiffness: float,
    sensitivity_nm: float,
    temperature_k: float,
    bell_exponent_limit: float,
) -> tuple:
    n = n_input
    m = m_input
    if n <= 0.0:
        n = 0.0
        m = 0.0
        q = 0.0
    else:
        q = m / n
    gap = max(gap_input, 0.0)
    area_m2 = _safe_nonnegative_product(
        area_um2, SQUARE_METERS_PER_SQUARE_MICROMETER
    )
    ligand_count = _safe_count(ligand_density, area_m2)
    target_count = _safe_count(target_density, area_m2)
    capacity = min(ligand_count, target_count)
    available_ligands = max(ligand_count - n, 0.0)
    available_targets = max(target_count - n, 0.0)
    formation_rate, formation_saturated = _formation_rate_scalar(
        association_rate,
        area_m2,
        available_ligands,
        available_targets,
    )

    bond_length = math.hypot(gap, q)
    bond_extension = (
        0.0
        if n <= 0.0 and formation_rate <= 0.0
        else max(bond_length - rest_length, 0.0)
    )
    tension = _safe_nonnegative_product(spring_stiffness, bond_extension)
    off_rate, bell_exponent, bell_saturated = _bell_rate_scalar(
        tension,
        zero_force_off_rate,
        sensitivity_nm,
        temperature_k,
        bell_exponent_limit,
    )

    extension_source = _safe_signed_product(n, slip)

    if n <= 0.0 or tension <= 0.0 or bond_length <= 0.0:
        force_t = 0.0
        force_n = 0.0
    else:
        total_tension = _safe_nonnegative_product(n, tension)
        force_t = -_safe_signed_product(total_tension, q / bond_length)
        force_n = -_safe_nonnegative_product(total_tension, gap / bond_length)
    torque_y = _safe_signed_product(radius, force_t)
    return (
        n,
        m,
        q,
        area_um2,
        ligand_count,
        target_count,
        capacity,
        formation_rate,
        off_rate,
        slip,
        extension_source,
        tension,
        force_t,
        force_n,
        torque_y,
        bell_exponent,
        bell_saturated,
        formation_saturated,
    )


def _evaluate_mean_field_batch_kernel(
    n,
    m,
    gap,
    radius,
    slip,
    area,
    ligand_density,
    target_density,
    association_rate,
    zero_force_off_rate,
    rest_length,
    spring_stiffness,
    sensitivity_nm,
    temperature_k,
    bell_exponent_limit,
):
    size = n.size
    n_out = np.empty(size, dtype=np.float64)
    m_out = np.empty(size, dtype=np.float64)
    q_out = np.empty(size, dtype=np.float64)
    area_out = np.empty(size, dtype=np.float64)
    ligand_count = np.empty(size, dtype=np.float64)
    target_count = np.empty(size, dtype=np.float64)
    capacity = np.empty(size, dtype=np.float64)
    formation_rate = np.empty(size, dtype=np.float64)
    off_rate = np.empty(size, dtype=np.float64)
    slip_out = np.empty(size, dtype=np.float64)
    extension_source = np.empty(size, dtype=np.float64)
    tension = np.empty(size, dtype=np.float64)
    force_t = np.empty(size, dtype=np.float64)
    force_n = np.empty(size, dtype=np.float64)
    torque_y = np.empty(size, dtype=np.float64)
    bell_exponent = np.empty(size, dtype=np.float64)
    bell_saturated = np.empty(size, dtype=np.bool_)
    formation_saturated = np.empty(size, dtype=np.bool_)
    for index in range(size):
        result = _evaluate_mean_field_scalar(
            n[index],
            m[index],
            gap[index],
            radius[index],
            slip[index],
            area[index],
            ligand_density,
            target_density,
            association_rate,
            zero_force_off_rate,
            rest_length,
            spring_stiffness,
            sensitivity_nm,
            temperature_k,
            bell_exponent_limit,
        )
        n_out[index] = result[0]
        m_out[index] = result[1]
        q_out[index] = result[2]
        area_out[index] = result[3]
        ligand_count[index] = result[4]
        target_count[index] = result[5]
        capacity[index] = result[6]
        formation_rate[index] = result[7]
        off_rate[index] = result[8]
        slip_out[index] = result[9]
        extension_source[index] = result[10]
        tension[index] = result[11]
        force_t[index] = result[12]
        force_n[index] = result[13]
        torque_y[index] = result[14]
        bell_exponent[index] = result[15]
        bell_saturated[index] = result[16]
        formation_saturated[index] = result[17]
    return (
        n_out,
        m_out,
        q_out,
        area_out,
        ligand_count,
        target_count,
        capacity,
        formation_rate,
        off_rate,
        slip_out,
        extension_source,
        tension,
        force_t,
        force_n,
        torque_y,
        bell_exponent,
        bell_saturated,
        formation_saturated,
    )


@register_jitable
def _decay_survival_scalar(off_rate: float, dt: float) -> float:
    if off_rate <= 0.0 or dt <= 0.0:
        return 1.0
    product = _safe_nonnegative_product(off_rate, dt)
    if product >= -math.log(float(np.nextafter(0.0, 1.0))):
        return 0.0
    return math.exp(-product)


@register_jitable
def _integrated_source_scalar(
    source: float,
    off_rate: float,
    dt: float,
    survival: float,
) -> float:
    if source == 0.0 or dt <= 0.0:
        return 0.0
    if off_rate <= 0.0:
        return _safe_signed_product(source, dt)
    factor = 1.0 - survival
    if factor <= 0.0:
        return 0.0
    log_magnitude = math.log(abs(source)) + math.log(factor) - math.log(off_rate)
    if log_magnitude >= LOG_MAX_FLOAT64:
        return math.copysign(MAX_FLOAT64, source)
    return math.copysign(math.exp(log_magnitude), source)


@register_jitable
def _advance_state_scalar(
    n_old: float,
    m_old: float,
    formation_rate: float,
    off_rate: float,
    extension_source: float,
    capacity: float,
    dt: float,
) -> tuple[float, float, bool]:
    if n_old <= 0.0:
        n_old = 0.0
        m_old = 0.0
    survival = _decay_survival_scalar(off_rate, dt)
    retained_n = _safe_nonnegative_product(n_old, survival)
    added_n = _integrated_source_scalar(formation_rate, off_rate, dt, survival)
    raw_n = _safe_signed_sum(retained_n, added_n)
    retained_m = _safe_signed_product(m_old, survival)
    added_m = _integrated_source_scalar(extension_source, off_rate, dt, survival)
    raw_m = _safe_signed_sum(retained_m, added_m)

    upper_bound = max(n_old, capacity)
    capacity_limited = raw_n > upper_bound
    if capacity_limited:
        raw_n = upper_bound
    n_new = max(raw_n, 0.0)
    m_new = raw_m if n_new > 0.0 else 0.0
    return n_new, m_new, capacity_limited


def _advance_state_batch_kernel(
    n,
    m,
    formation_rate,
    off_rate,
    extension_source,
    capacity,
    dt,
):
    size = n.size
    n_new = np.empty(size, dtype=np.float64)
    m_new = np.empty(size, dtype=np.float64)
    limited = np.empty(size, dtype=np.bool_)
    for index in range(size):
        n_new[index], m_new[index], limited[index] = _advance_state_scalar(
            n[index],
            m[index],
            formation_rate[index],
            off_rate[index],
            extension_source[index],
            capacity[index],
            dt,
        )
    return n_new, m_new, limited


def _advance_bond_state(
    expected_bond_count,
    total_tangential_extension_um,
    formation_rate,
    off_rate,
    extension_source,
    capacity,
    dt_s: float,
    *,
    use_numba: bool,
) -> BondStateUpdate:
    dt = float(dt_s)
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("dt_s must be finite and greater than zero.")
    arrays = np.broadcast_arrays(
        np.asarray(expected_bond_count, dtype=np.float64),
        np.asarray(total_tangential_extension_um, dtype=np.float64),
        np.asarray(formation_rate, dtype=np.float64),
        np.asarray(off_rate, dtype=np.float64),
        np.asarray(extension_source, dtype=np.float64),
        np.asarray(capacity, dtype=np.float64),
    )
    n, m, formation, dissociation, extension, cap = arrays
    if np.any(~np.isfinite(n)) or np.any(n < 0.0):
        raise ValueError("expected_bond_count must contain finite non-negative values.")
    if np.any(~np.isfinite(m)):
        raise ValueError("total_tangential_extension_um must contain finite values.")
    for name, value in (
        ("formation_rate", formation),
        ("off_rate", dissociation),
        ("capacity", cap),
    ):
        if np.any(~np.isfinite(value)) or np.any(value < 0.0):
            raise ValueError(f"{name} must contain finite non-negative values.")
    if np.any(~np.isfinite(extension)):
        raise ValueError("extension_source must contain finite values.")
    shape = n.shape
    flat = tuple(np.ascontiguousarray(value.reshape(-1)) for value in arrays)
    kernel = (
        _advance_state_batch_numba
        if use_numba and njit is not None
        else _advance_state_batch_kernel
    )
    updated_n, updated_m, limited = kernel(*flat, dt)
    return BondStateUpdate(
        expected_bond_count=np.asarray(updated_n).reshape(shape),
        total_tangential_extension_um=np.asarray(updated_m).reshape(shape),
        capacity_limited=np.asarray(limited).reshape(shape),
    )


def _finite_half_sum(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    first_array, second_array = np.broadcast_arrays(
        np.asarray(first, dtype=np.float64), np.asarray(second, dtype=np.float64)
    )
    # Halving before addition cannot overflow when both inputs are finite.
    return 0.5 * first_array + 0.5 * second_array


if njit is not None:
    _evaluate_mean_field_batch_numba = njit(cache=True)(
        _evaluate_mean_field_batch_kernel
    )
    _advance_state_batch_numba = njit(cache=True)(_advance_state_batch_kernel)
else:  # pragma: no cover - the project environment includes Numba
    _evaluate_mean_field_batch_numba = _evaluate_mean_field_batch_kernel
    _advance_state_batch_numba = _advance_state_batch_kernel
