"""Time-grid planning and resolution diagnostics for particle transport.

The saved trajectory cadence and the numerical particle step serve different
purposes.  A saved frame is an observation for downstream analysis, whereas an
internal substep controls how often position-dependent mobility, collisions,
and geometry constraints are recomputed.  This module keeps that distinction
explicit without depending on a particular particle-force implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral


@dataclass(frozen=True)
class ParticleTimeStepPlan:
    """Exact relationship between saved frames and internal integration steps."""

    output_intervals: int
    stored_frames: int
    output_dt_s: float
    integration_substeps: int
    internal_dt_s: float
    total_internal_steps: int
    expected_rhs_evaluations: int


@dataclass(frozen=True)
class ParticleStepResolution:
    """Dimensionless diagnostics for the largest accepted integration proposal."""

    maximum_internal_step_displacement_um: float
    grid_displacement_ratio: float
    radius_displacement_ratio: float
    field_size_displacement_ratio: float
    conservative_collision_layer_ratio: float
    collision_dt_over_relaxation_time: float


def build_particle_time_step_plan(
    output_intervals: int,
    output_dt_s: float,
    integration_substeps: int,
    time_integrator: str,
) -> ParticleTimeStepPlan:
    """Validate and construct an exactly aligned fixed-substep schedule.

    Dividing the output interval by an integer number of substeps guarantees
    that every internal trajectory lands exactly on the next saved-frame time;
    no shortened remainder step or cumulative clock drift is required.
    """

    intervals = int(output_intervals)
    output_dt = float(output_dt_s)
    if isinstance(integration_substeps, bool) or not isinstance(
        integration_substeps, Integral
    ):
        raise ValueError("Particle integration_substeps must be an integer.")
    substeps = int(integration_substeps)
    integrator = str(time_integrator).lower()
    if intervals < 0:
        raise ValueError("Particle output intervals must be non-negative.")
    if not math.isfinite(output_dt) or output_dt <= 0.0:
        raise ValueError("Particle output dt must be finite and positive.")
    if substeps < 1:
        raise ValueError("Particle integration_substeps must be at least one.")
    if integrator not in {"euler", "heun"}:
        raise ValueError("Particle time integrator must be 'euler' or 'heun'.")

    total_internal_steps = intervals * substeps
    stage_evaluations = 1 if integrator == "euler" else 2
    return ParticleTimeStepPlan(
        output_intervals=intervals,
        stored_frames=intervals + 1,
        output_dt_s=output_dt,
        integration_substeps=substeps,
        internal_dt_s=output_dt / substeps,
        total_internal_steps=total_internal_steps,
        expected_rhs_evaluations=total_internal_steps * stage_evaluations + 1,
    )


def build_particle_step_resolution(
    maximum_internal_step_displacement_um: float,
    grid_spacing_um: float,
    minimum_radius_um: float,
    collision_layer_um: float,
    collisions_enabled: bool,
    internal_dt_s: float,
    collision_relaxation_time_s: float,
) -> ParticleStepResolution:
    """Separate transport resolution from the conservative collision-layer ratio.

    Absolute motion relative to the grid and particle radius measures whether
    spatially varying fields and finite-size motion are time-resolved.  Absolute
    advection relative to a collision layer is retained as a conservative legacy
    diagnostic, but it is not used as the sole convergence warning because two
    bubbles may translate together without changing their relative separation.
    """

    displacement = max(float(maximum_internal_step_displacement_um), 0.0)
    spacing = float(grid_spacing_um)
    radius = float(minimum_radius_um)
    if not math.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("grid_spacing_um must be finite and positive.")
    if not math.isfinite(radius) or radius <= 0.0:
        raise ValueError("minimum_radius_um must be finite and positive.")

    grid_ratio = displacement / spacing
    radius_ratio = displacement / radius
    conservative_ratio = 0.0
    if bool(collisions_enabled) and float(collision_layer_um) > 0.0:
        conservative_ratio = displacement / float(collision_layer_um)

    relaxation_ratio = 0.0
    if bool(collisions_enabled):
        relaxation = float(collision_relaxation_time_s)
        if not math.isfinite(relaxation) or relaxation <= 0.0:
            raise ValueError("collision_relaxation_time_s must be finite and positive.")
        relaxation_ratio = float(internal_dt_s) / relaxation

    return ParticleStepResolution(
        maximum_internal_step_displacement_um=displacement,
        grid_displacement_ratio=grid_ratio,
        radius_displacement_ratio=radius_ratio,
        field_size_displacement_ratio=max(grid_ratio, radius_ratio),
        conservative_collision_layer_ratio=conservative_ratio,
        collision_dt_over_relaxation_time=relaxation_ratio,
    )
