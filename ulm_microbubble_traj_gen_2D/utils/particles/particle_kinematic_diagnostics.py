"""Step-level kinematic checks for Revised-v15 particle integration.

These helpers never change a trajectory.  They compare the displacement that
was actually committed with the displacement implied by the generalized
velocity accepted for the same positive physical interval.  Keeping this
calculation outside the integrator prevents diagnostics from becoming another
position-correction rule.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np

try:
    from numba import njit
except ImportError:  # pragma: no cover - used only in minimal environments
    njit = None


@dataclass(frozen=True, slots=True)
class StepKinematicMetrics:
    """Per-lane distances and finite-step wall-gap consistency residuals."""

    position_path_um: np.ndarray
    velocity_path_um: np.ndarray
    position_to_velocity_path_ratio: np.ndarray
    gap_change_um: np.ndarray
    velocity_predicted_gap_change_um: np.ndarray
    gap_kinematic_residual_um: np.ndarray
    nonzero_velocity_zero_progress: np.ndarray


@dataclass(frozen=True, slots=True)
class ReducedStepKinematicDiagnostics:
    """Scalar diagnostics accumulated without allocating per-lane outputs."""

    contact_position_path_sum_um: float
    contact_velocity_path_sum_um: float
    contact_nonzero_velocity_zero_progress_count: int
    minimum_contact_path_ratio: float
    maximum_contact_path_ratio: float
    maximum_free_gap_residual_um: float
    maximum_live_displacement_um: float


@dataclass(frozen=True, slots=True)
class ReducedPredictiveStepDiagnostics:
    """All scalar v15 transport audits reduced in one lane pass."""

    live_count: int
    active_contact_count: int
    maximum_contact_reaction_force_pn: float
    minimum_path_gap_um: float
    maximum_complementarity_residual_pn_um: float
    contact_interval_count: int
    contact_position_path_sum_um: float
    contact_velocity_path_sum_um: float
    contact_nonzero_velocity_zero_progress_count: int
    minimum_contact_path_ratio: float
    maximum_contact_path_ratio: float
    maximum_free_gap_residual_um: float
    maximum_live_displacement_um: float


def evaluate_step_kinematics(
    start_positions_grid: np.ndarray,
    end_positions_grid: np.ndarray,
    constrained_velocity_xz_um_s: np.ndarray,
    start_true_gap_um: np.ndarray,
    end_true_gap_um: np.ndarray,
    start_inward_normal_xz: np.ndarray,
    accepted: np.ndarray,
    *,
    spacing_um: float,
    duration_s: float,
) -> StepKinematicMetrics:
    """Compare one accepted physical step with its accepted velocity.

    The gap prediction uses the first-order identity ``dg/dt = n dot V``.
    It is not used to accept or reject a state.  Its residual should decrease
    during a physical-time refinement study when the geometry and velocity are
    sampled consistently.
    """

    start = np.asarray(start_positions_grid, dtype=np.float64)
    end = np.asarray(end_positions_grid, dtype=np.float64)
    velocity = np.asarray(constrained_velocity_xz_um_s, dtype=np.float64)
    gap_start = np.asarray(start_true_gap_um, dtype=np.float64)
    gap_end = np.asarray(end_true_gap_um, dtype=np.float64)
    normal = np.asarray(start_inward_normal_xz, dtype=np.float64)
    active = np.asarray(accepted, dtype=bool)

    count = start.shape[0]
    if start.shape != (count, 2) or end.shape != (count, 2):
        raise ValueError("Step positions must have shape (n, 2).")
    if velocity.shape != (count, 2) or normal.shape != (count, 2):
        raise ValueError("Step velocities and wall normals must have shape (n, 2).")
    if gap_start.shape != (count,) or gap_end.shape != (count,):
        raise ValueError("Step wall gaps must have shape (n,).")
    if active.shape != (count,):
        raise ValueError("accepted must have shape (n,).")

    spacing = float(spacing_um)
    duration = float(duration_s)
    if not np.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("spacing_um must be finite and positive.")
    if not np.isfinite(duration) or duration < 0.0:
        raise ValueError("duration_s must be finite and non-negative.")

    displacement_um = (end - start) * spacing
    position_path = np.linalg.norm(displacement_um, axis=1)
    velocity_path = duration * np.linalg.norm(velocity, axis=1)
    ratio = np.full(count, np.nan, dtype=np.float64)
    np.divide(
        position_path,
        velocity_path,
        out=ratio,
        where=active & np.isfinite(velocity_path) & (velocity_path > 0.0),
    )

    gap_change = gap_end - gap_start
    predicted_gap_change = duration * np.sum(normal * velocity, axis=1)
    gap_residual = gap_change - predicted_gap_change

    scale = np.maximum.reduce(
        (
            np.ones(count, dtype=np.float64),
            position_path,
            velocity_path,
            np.abs(gap_start),
            np.abs(gap_end),
        )
    )
    numerical_zero = 128.0 * np.finfo(np.float64).eps * scale
    zero_progress = active & (velocity_path > numerical_zero) & (
        position_path <= numerical_zero
    )

    for values in (
        position_path,
        velocity_path,
        ratio,
        gap_change,
        predicted_gap_change,
        gap_residual,
    ):
        values[~active] = np.nan

    return StepKinematicMetrics(
        position_path_um=position_path,
        velocity_path_um=velocity_path,
        position_to_velocity_path_ratio=ratio,
        gap_change_um=gap_change,
        velocity_predicted_gap_change_um=predicted_gap_change,
        gap_kinematic_residual_um=gap_residual,
        nonzero_velocity_zero_progress=zero_progress,
    )


def reduce_step_kinematic_diagnostics(
    start_positions_grid: np.ndarray,
    end_positions_grid: np.ndarray,
    constrained_velocity_xz_um_s: np.ndarray,
    start_true_gap_um: np.ndarray,
    end_true_gap_um: np.ndarray,
    start_inward_normal_xz: np.ndarray,
    contact_lanes: np.ndarray,
    free_gap_lanes: np.ndarray,
    live_lanes: np.ndarray,
    *,
    spacing_um: float,
    duration_s: float,
    use_numba: bool,
    inputs_prevalidated: bool = False,
) -> ReducedStepKinematicDiagnostics:
    """Reduce the two transport kinematic audits in one allocation-free pass."""

    start = np.asarray(start_positions_grid, dtype=np.float64)
    end = np.asarray(end_positions_grid, dtype=np.float64)
    velocity = np.asarray(constrained_velocity_xz_um_s, dtype=np.float64)
    gap_start = np.asarray(start_true_gap_um, dtype=np.float64)
    gap_end = np.asarray(end_true_gap_um, dtype=np.float64)
    normal = np.asarray(start_inward_normal_xz, dtype=np.float64)
    contact = np.asarray(contact_lanes, dtype=np.bool_)
    free = np.asarray(free_gap_lanes, dtype=np.bool_)
    live = np.asarray(live_lanes, dtype=np.bool_)
    count = int(start.shape[0]) if start.ndim == 2 else -1
    spacing = float(spacing_um)
    duration = float(duration_s)
    if not inputs_prevalidated:
        if start.shape != (count, 2) or end.shape != (count, 2):
            raise ValueError("Step positions must have shape (n, 2).")
        if velocity.shape != (count, 2) or normal.shape != (count, 2):
            raise ValueError(
                "Step velocities and wall normals must have shape (n, 2)."
            )
        if gap_start.shape != (count,) or gap_end.shape != (count,):
            raise ValueError("Step wall gaps must have shape (n,).")
        if contact.shape != (count,) or free.shape != (count,) or live.shape != (count,):
            raise ValueError("Kinematic diagnostic masks must have shape (n,).")
        if not np.isfinite(spacing) or spacing <= 0.0:
            raise ValueError("spacing_um must be finite and positive.")
        if not np.isfinite(duration) or duration < 0.0:
            raise ValueError("duration_s must be finite and non-negative.")

    kernel = (
        _reduce_step_kinematic_diagnostics_numba
        if use_numba and njit is not None
        else _reduce_step_kinematic_diagnostics_kernel
    )
    values = kernel(
        start,
        end,
        velocity,
        gap_start,
        gap_end,
        normal,
        contact,
        free,
        live,
        spacing,
        duration,
    )
    return ReducedStepKinematicDiagnostics(
        contact_position_path_sum_um=float(values[0]),
        contact_velocity_path_sum_um=float(values[1]),
        contact_nonzero_velocity_zero_progress_count=int(values[2]),
        minimum_contact_path_ratio=float(values[3]),
        maximum_contact_path_ratio=float(values[4]),
        maximum_free_gap_residual_um=float(values[5]),
        maximum_live_displacement_um=float(values[6]),
    )


def reduce_predictive_step_diagnostics(
    start_positions_grid: np.ndarray,
    end_positions_grid: np.ndarray,
    constrained_velocity_xz_um_s: np.ndarray,
    start_true_gap_um: np.ndarray,
    end_true_gap_um: np.ndarray,
    start_inward_normal_xz: np.ndarray,
    alive: np.ndarray,
    active_contact: np.ndarray,
    reaction_force_pn: np.ndarray,
    minimum_path_gap_um: np.ndarray,
    complementarity_residual_pn_um: np.ndarray,
    *,
    spacing_um: float,
    duration_s: float,
    use_numba: bool,
    inputs_prevalidated: bool = False,
) -> ReducedPredictiveStepDiagnostics:
    """Reduce every per-step v15 diagnostic without temporary masks.

    The production transport loop invokes this function once per accepted
    internal step.  Keeping lane classification and all scalar reductions in
    one serial kernel avoids many short NumPy allocations while preserving the
    deterministic lane order used by the reference implementation.
    """

    start = np.asarray(start_positions_grid, dtype=np.float64)
    end = np.asarray(end_positions_grid, dtype=np.float64)
    velocity = np.asarray(constrained_velocity_xz_um_s, dtype=np.float64)
    gap_start = np.asarray(start_true_gap_um, dtype=np.float64)
    gap_end = np.asarray(end_true_gap_um, dtype=np.float64)
    normal = np.asarray(start_inward_normal_xz, dtype=np.float64)
    live = np.asarray(alive, dtype=np.bool_)
    contact = np.asarray(active_contact, dtype=np.bool_)
    reaction = np.asarray(reaction_force_pn, dtype=np.float64)
    path_gap = np.asarray(minimum_path_gap_um, dtype=np.float64)
    complementarity = np.asarray(
        complementarity_residual_pn_um, dtype=np.float64
    )
    count = int(start.shape[0]) if start.ndim == 2 else -1
    spacing = float(spacing_um)
    duration = float(duration_s)
    if not inputs_prevalidated:
        if start.shape != (count, 2) or end.shape != (count, 2):
            raise ValueError("Step positions must have shape (n, 2).")
        if velocity.shape != (count, 2) or normal.shape != (count, 2):
            raise ValueError(
                "Step velocities and wall normals must have shape (n, 2)."
            )
        one_dimensional = (
            gap_start,
            gap_end,
            live,
            contact,
            reaction,
            path_gap,
            complementarity,
        )
        if any(values.shape != (count,) for values in one_dimensional):
            raise ValueError("V15 diagnostic arrays must have shape (n,).")
        if not np.isfinite(spacing) or spacing <= 0.0:
            raise ValueError("spacing_um must be finite and positive.")
        if not np.isfinite(duration) or duration < 0.0:
            raise ValueError("duration_s must be finite and non-negative.")

    kernel = (
        _reduce_predictive_step_diagnostics_numba
        if use_numba and njit is not None
        else _reduce_predictive_step_diagnostics_kernel
    )
    values = kernel(
        start,
        end,
        velocity,
        gap_start,
        gap_end,
        normal,
        live,
        contact,
        reaction,
        path_gap,
        complementarity,
        spacing,
        duration,
    )
    return ReducedPredictiveStepDiagnostics(
        live_count=int(values[0]),
        active_contact_count=int(values[1]),
        maximum_contact_reaction_force_pn=float(values[2]),
        minimum_path_gap_um=float(values[3]),
        maximum_complementarity_residual_pn_um=float(values[4]),
        contact_interval_count=int(values[5]),
        contact_position_path_sum_um=float(values[6]),
        contact_velocity_path_sum_um=float(values[7]),
        contact_nonzero_velocity_zero_progress_count=int(values[8]),
        minimum_contact_path_ratio=float(values[9]),
        maximum_contact_path_ratio=float(values[10]),
        maximum_free_gap_residual_um=float(values[11]),
        maximum_live_displacement_um=float(values[12]),
    )


def _reduce_step_kinematic_diagnostics_kernel(
    start,
    end,
    velocity,
    gap_start,
    gap_end,
    normal,
    contact,
    free,
    live,
    spacing,
    duration,
):
    position_sum = 0.0
    velocity_sum = 0.0
    zero_progress_count = 0
    minimum_ratio = math.inf
    maximum_ratio = -math.inf
    maximum_free_residual = -math.inf
    maximum_displacement = 0.0
    epsilon = 128.0 * np.finfo(np.float64).eps

    for lane in range(start.shape[0]):
        dx = (end[lane, 0] - start[lane, 0]) * spacing
        dz = (end[lane, 1] - start[lane, 1]) * spacing
        position_path = math.hypot(dx, dz)
        velocity_path = duration * math.hypot(
            velocity[lane, 0], velocity[lane, 1]
        )

        if live[lane] and math.isfinite(position_path):
            maximum_displacement = max(maximum_displacement, position_path)

        if contact[lane]:
            if math.isfinite(position_path):
                position_sum += position_path
            if math.isfinite(velocity_path):
                velocity_sum += velocity_path
            if math.isfinite(velocity_path) and velocity_path > 0.0:
                ratio = position_path / velocity_path
                if math.isfinite(ratio):
                    minimum_ratio = min(minimum_ratio, ratio)
                    maximum_ratio = max(maximum_ratio, ratio)
            scale = max(
                1.0,
                position_path,
                velocity_path,
                abs(gap_start[lane]),
                abs(gap_end[lane]),
            )
            numerical_zero = epsilon * scale
            if velocity_path > numerical_zero and position_path <= numerical_zero:
                zero_progress_count += 1

        if free[lane]:
            residual = (
                gap_end[lane]
                - gap_start[lane]
                - duration
                * (
                    normal[lane, 0] * velocity[lane, 0]
                    + normal[lane, 1] * velocity[lane, 1]
                )
            )
            if math.isfinite(residual):
                maximum_free_residual = max(
                    maximum_free_residual, abs(residual)
                )

    if minimum_ratio == math.inf:
        minimum_ratio = math.nan
    if maximum_ratio == -math.inf:
        maximum_ratio = math.nan
    if maximum_free_residual == -math.inf:
        maximum_free_residual = math.nan
    return (
        position_sum,
        velocity_sum,
        zero_progress_count,
        minimum_ratio,
        maximum_ratio,
        maximum_free_residual,
        maximum_displacement,
    )


def _reduce_predictive_step_diagnostics_kernel(
    start,
    end,
    velocity,
    gap_start,
    gap_end,
    normal,
    live,
    active_contact,
    reaction,
    path_gap,
    complementarity,
    spacing,
    duration,
):
    live_count = 0
    active_contact_count = 0
    maximum_reaction = -math.inf
    minimum_gap = math.inf
    maximum_complementarity = -math.inf
    contact_interval_count = 0
    position_sum = 0.0
    velocity_sum = 0.0
    zero_progress_count = 0
    minimum_ratio = math.inf
    maximum_ratio = -math.inf
    maximum_free_residual = -math.inf
    maximum_displacement = 0.0
    machine_epsilon = np.finfo(np.float64).eps

    for lane in range(start.shape[0]):
        if not live[lane]:
            continue
        live_count += 1
        if active_contact[lane]:
            active_contact_count += 1
        if math.isfinite(reaction[lane]):
            maximum_reaction = max(maximum_reaction, reaction[lane])
        if math.isfinite(path_gap[lane]):
            minimum_gap = min(minimum_gap, path_gap[lane])
        if math.isfinite(complementarity[lane]):
            maximum_complementarity = max(
                maximum_complementarity, abs(complementarity[lane])
            )
        dx = (end[lane, 0] - start[lane, 0]) * spacing
        dz = (end[lane, 1] - start[lane, 1]) * spacing
        position_path = math.hypot(dx, dz)
        velocity_path = duration * math.hypot(
            velocity[lane, 0], velocity[lane, 1]
        )
        if math.isfinite(position_path):
            maximum_displacement = max(maximum_displacement, position_path)

        metric_lane = math.isfinite(gap_start[lane]) and math.isfinite(
            gap_end[lane]
        )
        if not metric_lane:
            continue
        gap_scale = max(1.0, abs(gap_start[lane]), abs(gap_end[lane]))
        contact_roundoff = 256.0 * machine_epsilon * gap_scale
        contact_lane = (
            active_contact[lane]
            or gap_start[lane] <= contact_roundoff
            or gap_end[lane] <= contact_roundoff
        )
        if contact_lane:
            contact_interval_count += 1
            if math.isfinite(position_path):
                position_sum += position_path
            if math.isfinite(velocity_path):
                velocity_sum += velocity_path
            if math.isfinite(velocity_path) and velocity_path > 0.0:
                ratio = position_path / velocity_path
                if math.isfinite(ratio):
                    minimum_ratio = min(minimum_ratio, ratio)
                    maximum_ratio = max(maximum_ratio, ratio)
            zero_scale = max(
                1.0,
                position_path,
                velocity_path,
                abs(gap_start[lane]),
                abs(gap_end[lane]),
            )
            numerical_zero = 128.0 * machine_epsilon * zero_scale
            if velocity_path > numerical_zero and position_path <= numerical_zero:
                zero_progress_count += 1

        if not active_contact[lane]:
            residual = (
                gap_end[lane]
                - gap_start[lane]
                - duration
                * (
                    normal[lane, 0] * velocity[lane, 0]
                    + normal[lane, 1] * velocity[lane, 1]
                )
            )
            if math.isfinite(residual):
                maximum_free_residual = max(
                    maximum_free_residual, abs(residual)
                )

    if maximum_reaction == -math.inf:
        maximum_reaction = math.nan
    if minimum_gap == math.inf:
        minimum_gap = math.nan
    if maximum_complementarity == -math.inf:
        maximum_complementarity = math.nan
    if minimum_ratio == math.inf:
        minimum_ratio = math.nan
    if maximum_ratio == -math.inf:
        maximum_ratio = math.nan
    if maximum_free_residual == -math.inf:
        maximum_free_residual = math.nan
    return (
        live_count,
        active_contact_count,
        maximum_reaction,
        minimum_gap,
        maximum_complementarity,
        contact_interval_count,
        position_sum,
        velocity_sum,
        zero_progress_count,
        minimum_ratio,
        maximum_ratio,
        maximum_free_residual,
        maximum_displacement,
    )


if njit is not None:
    _reduce_step_kinematic_diagnostics_numba = njit(cache=True, nogil=True)(
        _reduce_step_kinematic_diagnostics_kernel
    )
    _reduce_predictive_step_diagnostics_numba = njit(cache=True, nogil=True)(
        _reduce_predictive_step_diagnostics_kernel
    )
else:  # pragma: no cover - exercised only without Numba
    _reduce_step_kinematic_diagnostics_numba = (
        _reduce_step_kinematic_diagnostics_kernel
    )
    _reduce_predictive_step_diagnostics_numba = (
        _reduce_predictive_step_diagnostics_kernel
    )
