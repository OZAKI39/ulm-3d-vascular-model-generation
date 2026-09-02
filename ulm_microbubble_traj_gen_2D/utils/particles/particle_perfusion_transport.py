"""Event-driven continuous perfusion around the established mobility RHS."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from time import perf_counter

import numpy as np

from ulm_vascular_model_generator.utils.core.models import Vessel

from ..cardiac.cardiac_pulsatility import CardiacPulsatility
from ..core.config import MolecularBindingConfig, ParticleConfig, ParticleDynamicsConfig
from ..runtime.console_output import print_key_values, print_section, print_warning
from ..molecular.molecular_binding import (
    MolecularBindingParameters,
    accept_bond_state_exponential_heun,
    predict_bond_state_exponential_euler,
    reaction_disk_radius_um,
)
from ..molecular.molecular_target_field import MolecularTargetField
from .particle_hydrodynamic_fields import ParticleHydrodynamicFields
from .particle_inlet_flux import build_inlet_flux_model
from .particle_geometry_tolerance import (
    particle_geometry_roundoff_tolerance_um,
    project_roundoff_negative_wall_gap_xz_um,
)
from .particle_perfusion_schedule import PerfusionSchedule, build_perfusion_schedule
from .particle_topological_ownership import (
    TopologicalCommitmentCatalog,
    TopologicalCrossingBatch,
    build_topological_commitment_catalog,
    inspect_topological_crossings,
)
from ..runtime.progress import create_particle_progress_bar
from .particle_backend import (
    configure_particle_numba_worker_threads,
    particle_backend_details,
    resolve_particle_backend,
)
from .particle_trajectory_kinematics import realized_center_velocities_um_s
from .particle_time_stepping import build_particle_step_resolution, build_particle_time_step_plan
from .particle_kinematic_diagnostics import (
    evaluate_step_kinematics,
    reduce_predictive_step_diagnostics,
)
from .red_blood_cell_transport import (
    CFL_DIAMETER_MAX_UM,
    CFL_DIAMETER_MIN_UM,
    CFL_PLATEAU_DIAMETER_UM,
    CFL_PLATEAU_THICKNESS_UM,
    RBC_MARGINATION_STRAIN,
    RBC_MAJOR_DIAMETER_UM,
    REFERENCE_SHEAR_RATE_S_INV,
    REFERENCE_TUBE_HEMATOCRIT,
    REFERENCE_TRANSVERSE_DIFFUSIVITY_UM2_S,
    SHEAR_RATE_QUANTITATIVE_MAX_S_INV,
    SHEAR_RATE_QUANTITATIVE_MIN_S_INV,
    TUBE_HEMATOCRIT_QUANTITATIVE_MAX,
    TUBE_HEMATOCRIT_QUANTITATIVE_MIN,
    RedBloodCellNetworkState,
    evaluate_red_blood_cell_transport,
)
from .particle_counter_rng import COUNTER_RNG_ALGORITHM, counter_normal_batch
from .particle_rbc_stochastic import (
    RbcStochasticSweepResult,
    sweep_rbc_stochastic_displacement,
)
from .particle_internal_target_exposure import integrate_internal_target_exposure
from .particle_predictive_contact import (
    PredictiveBoundaryLifecycleError,
    PredictiveContactStep,
    PredictiveWallContactFailure,
    PredictiveWallContactGeometryError,
    inspect_predictive_wall_contact_failure,
    solve_predictive_contact_step,
)
from .particle_constrained_step import (
    BatchLocalState,
    PhysicalTimeInterval,
    PhysicalTimeRefinementError,
    split_physical_time_interval,
)
from .particle_field_sampling import sample_bilinear
from ..core.types import FlowField, GridDomain, ParticleTrajectories, RasterizedVessels


_TERMINATION_OUTLET = np.uint8(1)


@dataclass(frozen=True, slots=True)
class ParticleWallContactFailure:
    """A v15 failure enriched with permanent lifecycle identity and time."""

    permanent_microbubble_id: int
    physical_time_s: float
    integration_stage: str
    contact: PredictiveWallContactFailure
    interval_start_position_xz_um: tuple[float, float] | None = None
    rejected_trial_displacement_um: float | None = None


class ParticleWallContactGeometryError(RuntimeError):
    """Wall-contact failure whose message is complete enough for a traceback."""

    def __init__(
        self,
        message: str,
        failures: tuple[ParticleWallContactFailure, ...],
    ) -> None:
        self.failures = failures
        super().__init__(_format_wall_contact_failure_message(message, failures))


@dataclass(frozen=True, slots=True)
class ParticleBoundaryLifecycleFailure:
    """An unclaimed continuous-boundary exit with permanent identity and time."""

    permanent_microbubble_id: int
    physical_time_s: float
    integration_stage: str
    local_lane: int
    start_position_xz_um: tuple[float, float]
    end_position_xz_um: tuple[float, float]
    bubble_radius_um: float
    nearest_open_section_index: int
    nearest_open_section_kind: int
    nearest_open_section_label: int
    section_signed_start_um: float
    section_signed_end_um: float
    section_lateral_end_um: float
    section_half_width_um: float


class ParticleBoundaryLifecycleError(RuntimeError):
    """An anatomical opening crossing was not classified by the lifecycle."""

    def __init__(
        self, failures: tuple[ParticleBoundaryLifecycleFailure, ...]
    ) -> None:
        self.failures = failures
        lines = [
            "A v16 particle path ended outside the continuous lumen without "
            "crossing an authoritative outlet.",
            "Boundary-lifecycle failure diagnostics:",
        ]
        for index, failure in enumerate(failures):
            kind = (
                "outlet"
                if failure.nearest_open_section_kind > 0
                else "inlet"
            )
            lines.extend(
                (
                    f"failure[{index}]:",
                    "  permanent_microbubble_id="
                    f"{failure.permanent_microbubble_id}",
                    f"  physical_time_s={_format_float(failure.physical_time_s)}",
                    f"  integration_stage={failure.integration_stage}",
                    f"  local_lane={failure.local_lane}",
                    "  start_position_xz_um="
                    f"{_format_pair(failure.start_position_xz_um)}",
                    "  end_position_xz_um="
                    f"{_format_pair(failure.end_position_xz_um)}",
                    f"  bubble_radius_um={_format_float(failure.bubble_radius_um)}",
                    "  nearest_open_section_index="
                    f"{failure.nearest_open_section_index}",
                    f"  nearest_open_section_kind={kind}",
                    "  nearest_open_section_label="
                    f"{failure.nearest_open_section_label}",
                    "  section_signed_start_um="
                    f"{_format_float(failure.section_signed_start_um)}",
                    "  section_signed_end_um="
                    f"{_format_float(failure.section_signed_end_um)}",
                    "  section_lateral_end_um="
                    f"{_format_float(failure.section_lateral_end_um)}",
                    "  section_half_width_um="
                    f"{_format_float(failure.section_half_width_um)}",
                )
            )
        super().__init__("\n".join(lines))


@dataclass(slots=True)
class _PerfusionState:
    position_grid: np.ndarray
    active: np.ndarray
    rotation_angle_rad: np.ndarray
    bond_count_expected: np.ndarray
    bond_total_tangential_extension_um: np.ndarray
    admission_time_s: np.ndarray
    exit_time_s: np.ndarray
    termination_reason: np.ndarray
    vessel_id: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int32)
    )
    last_generalized_velocity: np.ndarray = field(
        default_factory=lambda: np.empty((0, 3), dtype=np.float64)
    )
    last_contact_reaction_force_pn: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.float64)
    )
    last_contact_active: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=bool)
    )
    last_free_normal_velocity_um_s: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.float64)
    )
    last_constrained_normal_velocity_um_s: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.float64)
    )
    last_step_valid: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=bool)
    )
    active_id_buffer: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int64)
    )
    active_count: int = 0
    waiting_ids: list[int] = field(default_factory=list)
    admission_event_ids: list[int] = field(default_factory=list)
    next_event: int = 0
    target_exposure_time_s: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.float64)
    )
    target_reaction_area_time_um2_s: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.float64)
    )
    target_exposure_event_count: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int64)
    )
    target_exposure_event_end_count: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int64)
    )
    target_exposure_open: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=bool)
    )
    target_exposure_applicable_time_s: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.float64)
    )


@dataclass(slots=True)
class _Diagnostics:
    maximum_physical_overlap_um: float = 0.0
    maximum_collision_compression_um: float = 0.0
    maximum_interacting_pairs: int = 0
    total_interacting_pair_evaluations: int = 0
    all_pairs_rhs_evaluations: int = 0
    cell_list_rhs_evaluations: int = 0
    maximum_reciprocity_relative_error: float = 0.0
    degenerate_near_wall_normal_evaluations: int = 0
    maximum_collision_speed_um_s: float = 0.0
    maximum_step_displacement_um: float = 0.0
    two_wall_observations: int = 0
    inlet_wait_events: int = 0
    maximum_inlet_wait_s: float = 0.0
    maximum_expected_bond_count: float = 0.0
    maximum_single_bond_tension_pn: float = 0.0
    maximum_total_bond_force_pn: float = 0.0
    maximum_bond_torque_pn_um: float = 0.0
    molecular_bell_saturation_count: int = 0
    molecular_formation_saturation_count: int = 0
    molecular_capacity_limited_accepted_step_count: int = 0
    contact_constraint_evaluations: int = 0
    active_contact_constraint_evaluations: int = 0
    maximum_contact_reaction_force_pn: float = 0.0
    contact_time_refinement_count: int = 0
    maximum_contact_time_refinement_depth: int = 0
    minimum_accepted_internal_wall_gap_um: float = float("inf")
    contact_nonzero_velocity_zero_progress_count: int = 0
    contact_residual_projection_count: int = 0
    maximum_contact_residual_projection_um: float = 0.0
    maximum_contact_complementarity_residual_pn_um: float = 0.0
    contact_kinematic_interval_evaluations: int = 0
    contact_cumulative_position_path_um: float = 0.0
    contact_cumulative_velocity_path_um: float = 0.0
    minimum_contact_interval_position_to_velocity_path_ratio: float = float("inf")
    maximum_contact_interval_position_to_velocity_path_ratio: float = float("-inf")
    maximum_free_gap_kinematic_residual_um: float = 0.0
    outlet_event_count: int = 0
    topological_transition_count: int = 0
    topological_event_bubble_id: list[int] = field(default_factory=list)
    topological_event_time_s: list[float] = field(default_factory=list)
    topological_event_from_vessel_id: list[int] = field(default_factory=list)
    topological_event_to_vessel_id: list[int] = field(default_factory=list)
    topological_event_section_index: list[int] = field(default_factory=list)
    topological_event_position_xz_um: list[tuple[float, float]] = field(default_factory=list)
    active_outside_lumen_violations: int = 0
    active_outside_accessible_domain_violations: int = 0
    discrete_accessibility_disagreement_records: int = 0
    maximum_red_blood_cell_speed_um_s: float = 0.0
    red_blood_cell_rhs_observations: int = 0
    red_blood_cell_nonunique_wall_suppressions: int = 0
    red_blood_cell_hematocrit_out_of_range_observations: int = 0
    red_blood_cell_shear_out_of_range_observations: int = 0
    red_blood_cell_quantitative_applicability_observations: int = 0
    red_blood_cell_random_wall_reflection_count: int = 0
    maximum_red_blood_cell_random_displacement_um: float = 0.0
    red_blood_cell_random_displacement_squared_sum_um2: float = 0.0
    red_blood_cell_random_displacement_observations: int = 0
    red_blood_cell_diffusion_enabled_observations: int = 0
    red_blood_cell_stochastic_observations: int = 0
    red_blood_cell_invalid_transverse_space_observations: int = 0


@dataclass(slots=True)
class _FrameRecords:
    offsets: list[int] = field(default_factory=lambda: [0])
    position_grid: list[np.ndarray] = field(default_factory=list)
    particle_velocity_xz: list[np.ndarray] = field(default_factory=list)
    wall_shear: list[np.ndarray] = field(default_factory=list)
    vessel_id: list[np.ndarray] = field(default_factory=list)
    bubble_id: list[np.ndarray] = field(default_factory=list)
    wall_gap: list[np.ndarray] = field(default_factory=list)
    wall_contact: list[np.ndarray] = field(default_factory=list)
    wall_normal: list[np.ndarray] = field(default_factory=list)
    fluid_velocity_xz: list[np.ndarray] = field(default_factory=list)
    angular_velocity: list[np.ndarray] = field(default_factory=list)
    rotation_angle: list[np.ndarray] = field(default_factory=list)
    collision_force: list[np.ndarray] = field(default_factory=list)
    collision_neighbors: list[np.ndarray] = field(default_factory=list)
    gap_ratio: list[np.ndarray] = field(default_factory=list)
    near_wall_weight: list[np.ndarray] = field(default_factory=list)
    two_wall_warning: list[np.ndarray] = field(default_factory=list)
    cardiac_multiplier: list[np.ndarray] = field(default_factory=list)
    bond_count_expected: list[np.ndarray] = field(default_factory=list)
    bond_total_tangential_extension_um: list[np.ndarray] = field(default_factory=list)
    bond_mean_tangential_extension_um: list[np.ndarray] = field(default_factory=list)
    bond_force_xz_pn: list[np.ndarray] = field(default_factory=list)
    bond_force_tangent_pn: list[np.ndarray] = field(default_factory=list)
    bond_force_normal_pn: list[np.ndarray] = field(default_factory=list)
    bond_torque_pn_um: list[np.ndarray] = field(default_factory=list)
    single_bond_tension_pn: list[np.ndarray] = field(default_factory=list)
    bond_formation_rate_bonds_s: list[np.ndarray] = field(default_factory=list)
    bond_dissociation_rate_s_inv: list[np.ndarray] = field(default_factory=list)
    target_reaction_area_um2: list[np.ndarray] = field(default_factory=list)
    available_ligand_count: list[np.ndarray] = field(default_factory=list)
    available_target_count: list[np.ndarray] = field(default_factory=list)
    target_overlap_fraction: list[np.ndarray] = field(default_factory=list)
    contact_constraint_active: list[np.ndarray] = field(default_factory=list)
    contact_reaction_force_pn: list[np.ndarray] = field(default_factory=list)
    contact_free_normal_velocity_um_s: list[np.ndarray] = field(default_factory=list)
    contact_constrained_normal_velocity_um_s: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_velocity_xz_um_s: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_drift_velocity_xz_um_s: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_fick_velocity_xz_um_s: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_local_vessel_diameter_um: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_discharge_hematocrit: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_tube_hematocrit: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_shear_rate_s_inv: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_cfl_width_um: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_target_gap_um: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_transverse_diffusivity_um2_s: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_margination_length_um: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_margination_time_s: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_scale_activation: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_nearest_wall_unique: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_hematocrit_in_quantitative_range: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_shear_rate_in_quantitative_range: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_quantitative_applicability: list[np.ndarray] = field(default_factory=list)
    red_blood_cell_transverse_space_valid: list[np.ndarray] = field(default_factory=list)


@dataclass(frozen=True, slots=True)
class _BatchAdvanceState:
    """Accepted state for one synchronous active-particle time interval."""

    particles: BatchLocalState
    alive: np.ndarray
    termination_reason: np.ndarray
    exit_time_s: np.ndarray
    vessel_id: np.ndarray = field(
        default_factory=lambda: np.empty(0, dtype=np.int32)
    )
    last_generalized_velocity: np.ndarray | None = None
    last_contact_reaction_force_pn: np.ndarray | None = None
    last_contact_active: np.ndarray | None = None
    last_free_normal_velocity_um_s: np.ndarray | None = None
    last_constrained_normal_velocity_um_s: np.ndarray | None = None
    last_step_valid: np.ndarray | None = None


@dataclass(frozen=True, slots=True)
class _BatchAdvanceResult:
    """Accepted local state plus diagnostics that are safe to commit."""

    state: _BatchAdvanceState
    diagnostics: _Diagnostics
    stochastic_sweep: RbcStochasticSweepResult | None = None
    red_blood_cell_quantitative_applicability: np.ndarray | None = None
    accepted_path_position_grid: tuple[np.ndarray, ...] = ()
    accepted_path_active: tuple[np.ndarray, ...] = ()


@dataclass(frozen=True, slots=True)
class _BatchAttemptResult:
    """One accepted chronological piece of a requested physical interval."""

    state: _BatchAdvanceState
    diagnostics: _Diagnostics
    accepted_end_time_s: float


@dataclass(frozen=True, slots=True)
class _GlobalRbcStochasticResult:
    """Committed stochastic outcome for one nominal internal step."""

    terminated_count: int
    permanent_ids: np.ndarray
    sweep: RbcStochasticSweepResult | None
    quantitative_applicability: np.ndarray | None


def _accepted_batch_state(
    start_state: _BatchAdvanceState,
    positions: np.ndarray,
    angles: np.ndarray,
    bond_count: np.ndarray,
    bond_extension: np.ndarray,
    *,
    alive: np.ndarray | None = None,
    vessel_id: np.ndarray | None = None,
    termination_reason: np.ndarray | None = None,
    exit_time_s: np.ndarray | None = None,
    last_generalized_velocity: np.ndarray | None = None,
    last_contact_reaction_force_pn: np.ndarray | None = None,
    last_contact_active: np.ndarray | None = None,
    last_free_normal_velocity_um_s: np.ndarray | None = None,
    last_constrained_normal_velocity_um_s: np.ndarray | None = None,
    last_step_valid: np.ndarray | None = None,
) -> _BatchAdvanceState:
    """Package newly accepted particle arrays with their lifecycle state."""

    return _BatchAdvanceState(
        particles=BatchLocalState(
            position_grid=positions,
            rotation_angle_rad=angles,
            bond_count_expected=bond_count,
            bond_total_tangential_extension_um=bond_extension,
        ),
        alive=start_state.alive if alive is None else alive,
        vessel_id=start_state.vessel_id if vessel_id is None else vessel_id,
        termination_reason=(
            start_state.termination_reason
            if termination_reason is None
            else termination_reason
        ),
        exit_time_s=start_state.exit_time_s if exit_time_s is None else exit_time_s,
        last_generalized_velocity=(
            start_state.last_generalized_velocity
            if last_generalized_velocity is None
            else last_generalized_velocity
        ),
        last_contact_reaction_force_pn=(
            start_state.last_contact_reaction_force_pn
            if last_contact_reaction_force_pn is None
            else last_contact_reaction_force_pn
        ),
        last_contact_active=(
            start_state.last_contact_active
            if last_contact_active is None
            else last_contact_active
        ),
        last_free_normal_velocity_um_s=(
            start_state.last_free_normal_velocity_um_s
            if last_free_normal_velocity_um_s is None
            else last_free_normal_velocity_um_s
        ),
        last_constrained_normal_velocity_um_s=(
            start_state.last_constrained_normal_velocity_um_s
            if last_constrained_normal_velocity_um_s is None
            else last_constrained_normal_velocity_um_s
        ),
        last_step_valid=(
            start_state.last_step_valid
            if last_step_valid is None
            else last_step_valid
        ),
    )
class _RefineContactTimeStep(RuntimeError):
    """Signal that the whole synchronous batch needs two physical half steps.

    The failed trial is intentionally discarded, so any diagnostic that
    describes that rejected trial must travel with this exception.  Explicit
    fields keep that bookkeeping independent of human-readable error text.
    """

    def __init__(
        self,
        message: str,
        *,
        failed_lanes: np.ndarray | None = None,
        failure_codes: np.ndarray | None = None,
        failure_positions_grid: np.ndarray | None = None,
        failure_event_fractions: np.ndarray | None = None,
        integration_stage: str = "predictor",
    ) -> None:
        super().__init__(message)
        self.failed_lanes = (
            np.empty(0, dtype=np.int64)
            if failed_lanes is None
            else np.asarray(failed_lanes, dtype=np.int64).copy()
        )
        self.failure_codes = (
            np.empty(0, dtype=np.int64)
            if failure_codes is None
            else np.asarray(failure_codes, dtype=np.int64).copy()
        )
        self.failure_positions_grid = (
            np.empty((0, 2), dtype=np.float64)
            if failure_positions_grid is None
            else np.asarray(failure_positions_grid, dtype=np.float64).reshape(-1, 2).copy()
        )
        self.failure_event_fractions = (
            np.empty(0, dtype=np.float64)
            if failure_event_fractions is None
            else np.asarray(failure_event_fractions, dtype=np.float64).reshape(-1).copy()
        )
        self.integration_stage = str(integration_stage)


class _SplitAtBoundaryEvent(RuntimeError):
    """Retry an interval after cutting it at a directed boundary event.

    This is a chronological lifecycle split, not a failed-contact refinement.
    No speculative state is committed.  Re-running the prefix ensures every
    particle and every molecular state reaches one shared physical event time.
    """

    def __init__(self, split_time_s: float) -> None:
        super().__init__("A particle interval must be recomputed at its boundary event time.")
        self.split_time_s = float(split_time_s)


def advect_particles_with_continuous_perfusion(
    domain: GridDomain,
    raster: RasterizedVessels,
    flow_field: FlowField,
    vessels: list[Vessel],
    particle_cfg: ParticleConfig,
    dynamics_cfg: ParticleDynamicsConfig,
    hydrodynamic_fields: ParticleHydrodynamicFields,
    *,
    effective_thickness_um: float,
    boundary_depth_cells: float,
    molecular_binding_cfg: MolecularBindingConfig,
    random_seed: int = 42,
    cardiac: CardiacPulsatility | None = None,
    molecular_target_field: MolecularTargetField | None = None,
    red_blood_cell_network: RedBloodCellNetworkState | None = None,
    topological_ownership: TopologicalCommitmentCatalog | None = None,
) -> ParticleTrajectories:
    """
    Record continuous perfusion with deterministic mechanics and optional RBC
    drift--diffusion from an empty lumen at time zero.
    """
    from .particle_mobility_transport import _EvaluationContext, _evaluate_rhs

    # ============================================================================================
    # ======= Build the time plan, inlet flux, and perfusion schedule.
    # ============================================================================================
    # Turn the requested frame interval into a clear time plan.
    # The plan also divides each saved-frame interval into smaller movement steps, because bubbles may move too far if they are updated only once per frame.
    time_plan = build_particle_time_step_plan(particle_cfg.n_steps, particle_cfg.dt_s, dynamics_cfg.integration_substeps, dynamics_cfg.time_integrator)

    # Calculate the physical duration of the saved trajectory.
    # Start a timer too, so the final report can show the real computing time.
    record_duration = float(time_plan.output_intervals * time_plan.output_dt_s)
    started         = perf_counter()

    # Build a physical description of the root inlet.
    # It combines flow, bubble size, concentration, and usable opening size to calculate how many finite-size bubbles should enter each second.
    inlet = build_inlet_flux_model(domain, flow_field, vessels, hydrodynamic_fields.boundary_geometry, particle_cfg, effective_thickness_um=effective_thickness_um,
                                   boundary_depth_cells=boundary_depth_cells)

    # Create every planned entry event for this recording period.
    # When heartbeat pulsation is enabled, entry times follow the changing flow rather than one fixed time interval.
    schedule = build_perfusion_schedule(inlet, record_duration, cardiac)

    # Respect the optional safety limit on the total number of unique bubbles.
    # Each event needs its own permanent ID, so stop before allocating memory if the schedule is larger than the user allowed.
    configured_capacity = int(particle_cfg.max_unique_bubbles)
    if configured_capacity > 0 and schedule.count > configured_capacity:
        raise ValueError(
            f"The deterministic schedule requires {schedule.count} permanent IDs, exceeding "
            f"particles.max_unique_bubbles={configured_capacity}."
        )

    # Choose the calculation engine. Numba can speed up repeated number work, while the Python path keeps the same model as a fallback.
    backend     = resolve_particle_backend(particle_cfg.acceleration_backend)
    backend_details = particle_backend_details(backend)
    backend_details = dict(backend_details)
    backend_details["particle_numba_worker_threads"] = int(
        configure_particle_numba_worker_threads(backend)
    )

    # Reserve one permanent storage position for every planned bubble ID.
    # A finished bubble's position is never reused for a new ID, which prevents two different bubble life histories from being mixed together.
    capacity    = schedule.count

    # Store the fixed radius that belongs to each permanent bubble ID.
    radii       = np.ascontiguousarray(schedule.radius_um, dtype=np.float64)
    boundary_geometry = hydrodynamic_fields.boundary_geometry
    maximum_radius = float(np.max(radii, initial=0.5 * particle_cfg.bubble_diameter_max_um))
    capture_distance = (
        float(molecular_binding_cfg.capture_distance_um)
        if molecular_binding_cfg.enabled
        else 0.0
    )
    required_finite_element_distance = max(
        maximum_radius
        + capture_distance
        + float(particle_cfg.contact_geometry_tolerance_um),
        math.sqrt(2.0) * float(domain.spacing_um)
        + float(particle_cfg.contact_geometry_tolerance_um),
    )
    configured_finite_element_distance = float(
        hydrodynamic_fields.hybrid_velocity.finite_element_distance_um
    )
    if configured_finite_element_distance < required_finite_element_distance:
        raise ValueError(
            "field.hybrid_finite_element_distance_um is too small for the "
            "largest bubble, molecular capture range, and Cartesian "
            "interpolation footprint: configured="
            f"{configured_finite_element_distance:.6g} um, required at least "
            f"{required_finite_element_distance:.6g} um."
        )
    if topological_ownership is None:
        topological_ownership = build_topological_commitment_catalog(
            vessels,
            boundary_geometry,
        )
    scheduled_owner = np.asarray(schedule.initial_vessel_id, dtype=np.int32)
    root_owner = set(np.asarray(topological_ownership.root_vessel_id, dtype=np.int32).tolist())
    if scheduled_owner.shape != (schedule.count,) or any(
        int(value) not in root_owner for value in scheduled_owner
    ):
        raise RuntimeError(
            "Every perfusion event must be born with a valid Revised-v20 root vessel ID."
        )
    backend_details = dict(backend_details)
    backend_details["particle_continuous_wall_numba_bin_width_cells"] = int(
        round(
            float(boundary_geometry._exact_bin_size_um)
            / float(domain.spacing_um)
        )
    )
    scheduled_world_xz = boundary_geometry.grid_to_world_xz(
        schedule.position_grid
    )
    scheduled_gaps = np.asarray(
        boundary_geometry.true_gap_at_xz_um(scheduled_world_xz, radii),
        dtype=np.float64,
    )
    scheduled_inside = np.asarray(
        boundary_geometry.contains_xz_um(scheduled_world_xz), dtype=bool
    )
    scheduled_accessible = np.asarray(
        boundary_geometry.is_accessible_grid(schedule.position_grid, radii),
        dtype=bool,
    )
    scheduled_gap_roundoff_um = particle_geometry_roundoff_tolerance_um(
        scheduled_world_xz,
        radii,
        float(domain.spacing_um),
    )
    if (
        np.any(~np.isfinite(scheduled_gaps))
        or np.any(scheduled_gaps < -scheduled_gap_roundoff_um)
        or not np.all(scheduled_inside)
        or not np.all(scheduled_accessible)
    ):
        raise RuntimeError(
            "The inlet schedule contains a finite-size bubble outside the "
            "canonical inlet-connected feasible domain."
        )

    # Molecular forces are optional. 
    molecular_parameters = None
    if molecular_binding_cfg.enabled:
        # Put the checked bond settings into one small object used by every time step.
        # This avoids passing a long list of separate numbers again and again.
        molecular_parameters = MolecularBindingParameters(
            ligand_density_molecules_m2=float(molecular_binding_cfg.ligand_density_molecules_per_m2),
            target_density_molecules_m2=float(molecular_target_field.target_density_molecules_per_m2),
            association_rate_m2_per_molecule_s=float(molecular_binding_cfg.association_rate_m2_per_molecule_s),
            zero_force_dissociation_rate_s=float(molecular_binding_cfg.zero_force_dissociation_rate_s),
            rest_length_um=float(molecular_binding_cfg.rest_length_um),
            spring_stiffness_pn_per_um=float(molecular_binding_cfg.bond_stiffness_pn_per_um),
            force_sensitivity_length_nm=float(molecular_binding_cfg.reactive_compliance_nm),
            temperature_k=float(molecular_binding_cfg.temperature_k),
            bell_exponent_limit=float(molecular_binding_cfg.bell_exponent_limit))

    # Gather all read-only fields needed by the movement-speed calculator.
    # "Context" simply means one box of shared information. 
    # Reusing this box is faster and clearer than rebuilding these large arrays at every small step.
    context = _EvaluationContext(
        velocity_xz_um_s=np.ascontiguousarray(flow_field.velocity_xz_um_s),
        wall_shear_stress_pa=np.ascontiguousarray(flow_field.wall_shear_stress_pa),
        vessel_id=np.ascontiguousarray(raster.vessel_id, dtype=np.int32),
        local_vessel_radius_um=np.ascontiguousarray(raster.radius_um),
        velocity_gradient_s_inv=np.ascontiguousarray(hydrodynamic_fields.velocity_gradient_s_inv),
        dynamic_viscosity_pa_s=np.ascontiguousarray(hydrodynamic_fields.dynamic_viscosity_pa_s),
        all_bubble_radii_um=radii,
        spacing_um=float(domain.spacing_um),
        dynamics=dynamics_cfg,
        use_numba=backend == "numba_cpu",
        random_seed=int(random_seed),
        contact_geometry_tolerance_um=float(
            particle_cfg.contact_geometry_tolerance_um
        ),
        cardiac=cardiac,
        origin_xz_um=(float(domain.origin_um[0]), float(domain.origin_um[2])),
        molecular_target_field=molecular_target_field,
        molecular_binding_parameters=molecular_parameters,
        molecular_capture_distance_um=(
            float(molecular_binding_cfg.capture_distance_um)
            if molecular_target_field is not None
            else 0.0
        ),
        molecular_mean_field_warning_count=(float(molecular_binding_cfg.mean_field_warning_count) if molecular_parameters is not None else 10.0),
        boundary_geometry=hydrodynamic_fields.boundary_geometry,
        hybrid_velocity=hydrodynamic_fields.hybrid_velocity,
        red_blood_cell_network=(
            red_blood_cell_network
            if red_blood_cell_network is not None and red_blood_cell_network.enabled
            else None
        ),
        topological_ownership=topological_ownership,
    )

    # Create the live state for every permanent bubble ID.
    # All bubbles begin outside the vessel: 
    # active flags are false, bond values are zero, and entry or exit times use NaN, a special marker meaning "not known yet".
    state = _PerfusionState(
        position_grid=np.zeros((capacity, 2), dtype=np.float64),
        active=np.zeros(capacity, dtype=bool),
        vessel_id=np.zeros(capacity, dtype=np.int32),
        rotation_angle_rad=np.zeros(capacity, dtype=np.float64),
        bond_count_expected=np.zeros(capacity, dtype=np.float64),
        bond_total_tangential_extension_um=np.zeros(capacity, dtype=np.float64),
        admission_time_s=np.full(capacity, np.nan, dtype=np.float64),
        exit_time_s=np.full(capacity, np.nan, dtype=np.float64),
        termination_reason=np.zeros(capacity, dtype=np.uint8),
        last_generalized_velocity=np.zeros((capacity, 3), dtype=np.float64),
        last_contact_reaction_force_pn=np.zeros(capacity, dtype=np.float64),
        last_contact_active=np.zeros(capacity, dtype=bool),
        last_free_normal_velocity_um_s=np.zeros(capacity, dtype=np.float64),
        last_constrained_normal_velocity_um_s=np.zeros(capacity, dtype=np.float64),
        last_step_valid=np.zeros(capacity, dtype=bool),
        active_id_buffer=np.empty(capacity, dtype=np.int64),
        target_exposure_time_s=np.zeros(capacity, dtype=np.float64),
        target_reaction_area_time_um2_s=np.zeros(capacity, dtype=np.float64),
        target_exposure_event_count=np.zeros(capacity, dtype=np.int64),
        target_exposure_event_end_count=np.zeros(capacity, dtype=np.int64),
        target_exposure_open=np.zeros(capacity, dtype=bool),
        target_exposure_applicable_time_s=np.zeros(capacity, dtype=np.float64),
    )

    # Start empty diagnostic counters. 
    # They do not change the movement, they tell us later whether the time step or a model assumption may need closer review.
    diagnostics = _Diagnostics()

    # Print the plan before expensive movement begins so the user can check the
    # concentration, entry rate, duration, frame count, and permanent ID count.
    print_section("Continuous perfusion plan")
    perfusion_rows: list[tuple[str, object]] = [
        ("Inlet concentration", f"{inlet.number_concentration_mb_per_ml:.6g} MB/mL"),
        ("Number concentration", f"{inlet.number_concentration_mb_per_m3:.6g} bubble/m^3"),
        ("Finite-size injection rate", f"{inlet.injection_rate_per_s:.6g} bubble/s"),
        ("Mean planned injection interval", f"{inlet.mean_injection_interval_s:.6g} s"),
        ("Stored frames", time_plan.stored_frames),
        ("Formal recording duration", f"{record_duration:.6g} s"),
        ("Scheduled permanent IDs", capacity),
        ("Initial condition", "empty lumen at formal frame 0"),
        ("Accepted inlet-section flow", f"{inlet.raw_section_flow_um2_s:.6g} um^2/s"),
    ]

    # Show the inlet-flow consistency check only when it produced a real number.
    if math.isfinite(inlet.section_flux_relative_error):
        perfusion_rows.append(("Inlet/reference flux difference", f"{inlet.section_flux_relative_error:.6g}"))

    print_key_values(perfusion_rows)
    if cardiac is not None:
        # Report the weakest and strongest heartbeat multipliers so the user can
        # see the real speed and entry-rate range before bubbles are advanced.
        minimum_multiplier = float(np.min(cardiac.waveform.multiplier))
        maximum_multiplier = float(np.max(cardiac.waveform.multiplier))
        print_section("Cardiac pulsatility plan")
        print_key_values(
            [
                ("Waveform", cardiac.waveform_name),
                ("Heart rate", f"{cardiac.waveform.bpm:.6g} BPM"),
                ("Cardiac period", f"{cardiac.waveform.period_s:.6g} s"),
                ("Modulation strength", f"{cardiac.modulation_strength:.6g}"),
                ("Formal cardiac cycles", f"{record_duration / cardiac.waveform.period_s:.6g}"),
                ("Inlet multiplier range", f"{minimum_multiplier:.6g}..{maximum_multiplier:.6g}"),
                ("Maximum root-path delay", f"{float(np.max(cardiac.delay_s)):.6g} s"),
                (
                    "Instantaneous injection-rate range",
                    f"{inlet.injection_rate_per_s * minimum_multiplier:.6g}.."
                    f"{inlet.injection_rate_per_s * maximum_multiplier:.6g} bubble/s",
                ),
            ]
        )

    # Collect only bubbles that truly exist in each saved frame.
    records             = _FrameRecords()

    # Keep simple per-frame counts for summaries and later plots.
    active_counts       = np.zeros(time_plan.stored_frames, dtype=np.int32)
    injected_counts     = np.zeros(time_plan.stored_frames, dtype=np.int32)
    terminated_counts   = np.zeros(time_plan.stored_frames, dtype=np.int32)

    # Save frame zero before any entry event is accepted.
    # This deliberately records the requested empty-vessel starting condition.
    _record_current_frame(state, context, particle_cfg, records, diagnostics, _evaluate_rhs, 0.0,)

    # Show progress in small internal movement steps.
    # absolute_time is physical simulation time, not real computing time.
    formal_progress = create_particle_progress_bar(time_plan.total_internal_steps)
    absolute_time   = 0.0
    global_internal_step = 0
    batched_saved_frame_count = 0
    cumulative_terminated_count = 0

    # ============================================================================================
    # ======= Build every saved frame after the empty starting frame.
    # ============================================================================================
    for frame in range(1, time_plan.stored_frames):
        # Reset counts that belong only to this saved-frame interval.
        frame_admitted   = 0
        frame_terminated = 0

        # Take several small movement steps before saving the next frame.
        # Small steps make walls, collisions, and quickly changing forces easier to follow than one large jump between saved frames.
        frame_end_time = absolute_time + time_plan.output_dt_s
        if (
            context.use_numba
            and time_plan.integration_substeps > 1
            and _can_batch_injection_free_frame(
                frame_end_time, state, schedule
            )
        ):
            frame_terminated, absolute_time = _advance_injection_free_substeps(
                absolute_time,
                int(time_plan.integration_substeps),
                float(time_plan.internal_dt_s),
                state,
                context,
                particle_cfg,
                dynamics_cfg,
                diagnostics,
                _evaluate_rhs,
                first_global_internal_step=global_internal_step,
            )
            global_internal_step += int(time_plan.integration_substeps)
            batched_saved_frame_count += 1
            formal_progress.update(time_plan.integration_substeps)
        else:
            for _ in range(time_plan.integration_substeps):
                next_time = absolute_time + time_plan.internal_dt_s

                target_exposure_enabled = (
                    context.molecular_target_field is not None
                    and bool(context.molecular_target_field.enabled)
                    and float(context.molecular_capture_distance_um) > 0.0
                )
                track_internal_path = (
                    context.red_blood_cell_network is not None
                    or target_exposure_enabled
                )
                step_start_ids = np.empty(0, dtype=np.int64)
                step_start_position_grid = np.empty((0, 2), dtype=np.float64)
                step_start_vessel_id = np.empty(0, dtype=np.int32)
                admission_log_start = len(state.admission_event_ids)
                accepted_target_path_batches = [] if target_exposure_enabled else None
                if track_internal_path:
                    step_start_ids = _active_ids_view(state).copy()
                    step_start_position_grid = np.asarray(
                        state.position_grid[step_start_ids], dtype=np.float64
                    ).copy()
                    step_start_vessel_id = np.asarray(
                        state.vessel_id[step_start_ids], dtype=np.int32
                    ).copy()

                # Advance exactly to the next small-step time. 
                # The helper may split this interval at an entry event, admit waiting bubbles, move active bubbles, 
                # apply wall rules, and remove bubbles that reached an exit.
                admitted, terminated = _advance_interval(
                    absolute_time, next_time, state, schedule, context, particle_cfg, dynamics_cfg,
                    diagnostics,
                    _evaluate_rhs,
                    accepted_target_path_batches=accepted_target_path_batches,
                )
                if track_internal_path:
                    newly_admitted_ids = np.asarray(
                        state.admission_event_ids[admission_log_start:], dtype=np.int64
                    )
                    if newly_admitted_ids.size:
                        step_start_ids = np.concatenate(
                            (step_start_ids, newly_admitted_ids)
                        )
                        step_start_position_grid = np.vstack(
                            (
                                step_start_position_grid,
                                np.asarray(
                                    schedule.position_grid[newly_admitted_ids],
                                    dtype=np.float64,
                                ),
                            )
                        )
                        step_start_vessel_id = np.concatenate(
                            (
                                step_start_vessel_id,
                                np.asarray(
                                    schedule.initial_vessel_id[newly_admitted_ids],
                                    dtype=np.int32,
                                ),
                            )
                        )
                    deterministic_end_position_grid = np.asarray(
                        state.position_grid[step_start_ids], dtype=np.float64
                    ).copy()
                else:
                    deterministic_end_position_grid = np.empty(
                        (0, 2), dtype=np.float64
                    )
                stochastic_result = _apply_rbc_stochastic_global_state(
                    state,
                    context,
                    diagnostics,
                    step_start_time_s=absolute_time,
                    step_end_time_s=next_time,
                    global_internal_step=global_internal_step,
                    step_start_ids=step_start_ids,
                    step_start_position_grid=step_start_position_grid,
                    step_start_vessel_id=step_start_vessel_id,
                    schedule=schedule,
                )
                terminated += stochastic_result.terminated_count
                if target_exposure_enabled and step_start_ids.size:
                    admission_time = np.asarray(
                        state.admission_time_s[step_start_ids], dtype=np.float64
                    )
                    exit_time = np.asarray(
                        state.exit_time_s[step_start_ids], dtype=np.float64
                    )
                    physical_end = np.where(
                        np.isfinite(exit_time),
                        np.minimum(exit_time, float(next_time)),
                        float(next_time),
                    )
                    active_duration = np.maximum(
                        physical_end
                        - np.maximum(admission_time, float(absolute_time)),
                        0.0,
                    )
                    _accumulate_internal_target_exposure(
                        state,
                        context,
                        permanent_ids=step_start_ids,
                        start_position_grid=step_start_position_grid,
                        start_vessel_id=step_start_vessel_id,
                        deterministic_end_position_grid=(
                            deterministic_end_position_grid
                        ),
                        final_position_grid=np.asarray(
                            state.position_grid[step_start_ids], dtype=np.float64
                        ),
                        active_duration_s=active_duration,
                        terminated_at_end=~np.asarray(
                            state.active[step_start_ids], dtype=bool
                        ),
                        stochastic_sweep=stochastic_result.sweep,
                        stochastic_permanent_ids=stochastic_result.permanent_ids,
                        deterministic_path_batches=(
                            accepted_target_path_batches
                        ),
                    )

                # Add this small step's changes to the current frame totals.
                frame_admitted      += admitted
                frame_terminated    += terminated

                # Move the clock forward only after the step has been accepted.
                absolute_time       = next_time
                global_internal_step += 1
                formal_progress.update(1)

        # Save one snapshot after all internal steps for this frame are complete.
        _record_current_frame(state, context, particle_cfg, records, diagnostics, _evaluate_rhs, absolute_time)

        # Store the number of active, newly entered, and finished bubbles.
        active_counts[frame]        = _active_ids_view(state).size
        injected_counts[frame]      = frame_admitted
        terminated_counts[frame]    = frame_terminated
        cumulative_terminated_count += frame_terminated

        # Recheck memory use after adding this variable-size frame.
        _check_record_limit(records, particle_cfg)

        # Update only the progress display; this does not change the simulation.
        formal_progress.set_postfix(
            time=f"{absolute_time:.3f}/{record_duration:.3f}s",
            frame=f"{frame}/{time_plan.output_intervals}",
            created=f"{len(state.admission_event_ids)}/{capacity}",
            active=int(active_counts[frame]),
            exited=int(cumulative_terminated_count),
            waiting=len(state.waiting_ids),
            refresh=False,
        )

    # Close the progress display cleanly after the final step.
    formal_progress.close()

    # Keep only entry events reached during the formal recording duration.
    # The registry is the master list with one row per permanent bubble ID.
    registry_count          = state.next_event
    registry_ids            = np.arange(registry_count, dtype=np.int64)

    # Convert each ID's fixed radius to the diameter used in public outputs.
    registry_diameter_um    = np.asarray(2.0 * radii[:registry_count], dtype=np.float32)

    # Copy real entry and exit times before converting them to frame numbers.
    relative_admission      = state.admission_time_s[:registry_count].copy()
    relative_exit           = state.exit_time_s[:registry_count].copy()

    # Begin with -1, which means that no observed birth frame exists yet.
    birth_frame             = np.full(registry_count, -1, dtype=np.int32)

    # A finite admission time means that the bubble really entered the vessel.
    admitted                = np.isfinite(relative_admission)

    # An entry at or before formal time zero belongs to frame zero.
    birth_frame[admitted & (relative_admission <= 0.0)] = 0

    # Later entries belong to the first saved frame at or after their entry time.
    # The tiny subtraction stops harmless number rounding from moving an event
    # exactly on a frame boundary into the next frame.
    after_start     = admitted & (relative_admission > 0.0)
    birth_frame[after_start] = np.ceil(relative_admission[after_start] / time_plan.output_dt_s - 1.0e-12).astype(np.int32)

    # Apply the same frame rule to real exit times. A bubble still inside the
    # vessel at the end keeps -1 because it has no death frame in this recording.
    death_frame     = np.full(registry_count, -1, dtype=np.int32)
    formal_death    = np.isfinite(relative_exit) & (relative_exit >= 0.0)
    death_frame[formal_death] = np.ceil(relative_exit[formal_death] / time_plan.output_dt_s - 1.0e-12).astype(np.int32)

    # Join the variable-size frame lists into flat arrays for compact saving.
    # Frame offsets stored elsewhere mark where each frame starts and ends, so no
    # information is lost even though different frames contain different counts.
    flat_ids                = _concat(records.bubble_id, (0,), np.int64)
    flat_position_grid      = _concat(records.position_grid, (0, 2), np.float64)
    flat_frame_offsets      = np.asarray(records.offsets, dtype=np.int64)

    # Change internal grid positions into real micrometre positions for users.
    positions_world         = np.empty((flat_position_grid.shape[0], 3), dtype=np.float32)
    positions_world[:, 0]   = domain.origin_um[0] + flat_position_grid[:, 0] * domain.spacing_um
    positions_world[:, 1]   = domain.fixed_y_um
    positions_world[:, 2]   = domain.origin_um[2] + flat_position_grid[:, 1] * domain.spacing_um

    # Add the fixed third coordinate expected by the public three-number output.
    # The physical motion is still the two-dimensional X-Z motion calculated above.
    particle_velocity_world = _embed_xz(_concat(records.particle_velocity_xz, (0, 2), np.float64))

    # Measure what the geometry-accepted bubble centres actually did between
    # stored frames. This remains separate from the instantaneous mobility RHS.
    realized_velocity_world = realized_center_velocities_um_s(
        flat_frame_offsets,
        flat_ids,
        positions_world,
        float(time_plan.output_dt_s),
    ).astype(np.float32)

    # Decide whether optional detailed arrays and molecular arrays should be saved.
    full_diagnostics        = bool(dynamics_cfg.store_full_diagnostics)
    molecular_enabled       = molecular_parameters is not None
    red_blood_cell_enabled  = context.red_blood_cell_network is not None
    target_exposure_enabled = (
        context.molecular_target_field is not None
        and bool(context.molecular_target_field.enabled)
        and float(context.molecular_capture_distance_um) > 0.0
    )
    _ensure_target_exposure_storage(state)
    registry_target_exposure_time_s = np.asarray(
        state.target_exposure_time_s[:registry_count], dtype=np.float64
    ).copy()
    registry_target_reaction_area_time_um2_s = np.asarray(
        state.target_reaction_area_time_um2_s[:registry_count], dtype=np.float64
    ).copy()
    registry_target_exposure_event_count = np.asarray(
        state.target_exposure_event_count[:registry_count], dtype=np.int64
    ).copy()
    registry_target_exposure_right_censored = np.asarray(
        state.target_exposure_open[:registry_count]
        & state.active[:registry_count],
        dtype=bool,
    )
    registry_target_exposure_applicability_fraction = np.divide(
        np.asarray(
            state.target_exposure_applicable_time_s[:registry_count],
            dtype=np.float64,
        ),
        registry_target_exposure_time_s,
        out=np.zeros(registry_count, dtype=np.float64),
        where=registry_target_exposure_time_s > 0.0,
    )
    target_exposure_total_time_s = float(
        np.sum(registry_target_exposure_time_s)
    )
    target_reaction_area_time_um2_s = float(
        np.sum(registry_target_reaction_area_time_um2_s)
    )
    target_exposure_event_count = int(
        np.sum(registry_target_exposure_event_count)
    )
    target_exposure_encountered_bubble_count = int(
        np.count_nonzero(registry_target_exposure_event_count > 0)
    )
    target_exposure_right_censored_count = int(
        np.count_nonzero(registry_target_exposure_right_censored)
    )
    red_blood_cell_exposure_weighted_applicability_fraction = (
        float(
            np.sum(
                state.target_exposure_applicable_time_s[:registry_count]
            )
            / target_exposure_total_time_s
        )
        if target_exposure_total_time_s > 0.0
        else 0.0
    )

    # Keep the heartbeat multiplier seen by each saved bubble observation only
    # when pulsation was active.
    cardiac_multiplier = (
        _concat(records.cardiac_multiplier, (0,), np.float32)
        if cardiac is not None
        else None
    )

    # Fluid velocity is useful for comparing carrier flow with actual bubble speed,
    # but it is omitted when detailed diagnostics were not requested.
    fluid_velocity_world    = (_embed_xz(_concat(records.fluid_velocity_xz, (0, 2), np.float64)) if full_diagnostics else None)

    # Every row in the flat records describes a bubble that was active in that frame.
    flat_active             = np.ones(flat_ids.size, dtype=bool)

    # Gather wall gaps and their smallest observed value for geometry checks.
    flat_gap                = _concat(records.wall_gap, (0,), np.float64)
    minimum_gap             = float(np.min(flat_gap)) if flat_gap.size else float("nan")
    flat_red_blood_cell_velocity_xz_um_s = (
        _concat(records.red_blood_cell_velocity_xz_um_s, (0, 2), np.float64)
        if red_blood_cell_enabled
        else None
    )
    flat_red_blood_cell_target_gap_um = (
        _concat(records.red_blood_cell_target_gap_um, (0,), np.float64)
        if red_blood_cell_enabled
        else None
    )
    flat_red_blood_cell_diffusivity_um2_s = (
        _concat(
            records.red_blood_cell_transverse_diffusivity_um2_s,
            (0,),
            np.float64,
        )
        if red_blood_cell_enabled
        else None
    )
    if red_blood_cell_enabled and flat_ids.size:
        reached_cfl = flat_gap <= flat_red_blood_cell_target_gap_um
        red_blood_cell_cfl_record_fraction = float(np.mean(reached_cfl))
        red_blood_cell_cfl_reaching_bubble_fraction = float(
            np.unique(flat_ids[reached_cfl]).size / np.unique(flat_ids).size
        )
        red_blood_cell_quantitative_record_fraction = float(
            np.mean(
                _concat(
                    records.red_blood_cell_quantitative_applicability,
                    (0,),
                    bool,
                )
            )
        )
        red_blood_cell_diffusion_enabled_record_fraction = float(
            np.mean(flat_red_blood_cell_diffusivity_um2_s > 0.0)
        )
    else:
        red_blood_cell_cfl_record_fraction = 0.0
        red_blood_cell_cfl_reaching_bubble_fraction = 0.0
        red_blood_cell_quantitative_record_fraction = 0.0
        red_blood_cell_diffusion_enabled_record_fraction = 0.0
    red_blood_cell_rms_random_displacement_um = (
        math.sqrt(
            diagnostics.red_blood_cell_random_displacement_squared_sum_um2
            / diagnostics.red_blood_cell_random_displacement_observations
        )
        if diagnostics.red_blood_cell_random_displacement_observations > 0
        else 0.0
    )
    flat_contact_constraint_active = _concat(
        records.contact_constraint_active, (0,), bool
    )
    flat_contact_reaction_force_pn = _concat(
        records.contact_reaction_force_pn, (0,), np.float64
    )
    flat_contact_free_normal_velocity_um_s = _concat(
        records.contact_free_normal_velocity_um_s, (0,), np.float64
    )
    flat_contact_constrained_normal_velocity_um_s = _concat(
        records.contact_constrained_normal_velocity_um_s, (0,), np.float64
    )

    # Use the smallest actual radius as the strictest bubble-size scale.
    # If no bubble entered, fall back to the inlet model's minimum possible radius.
    minimum_radius          = float(np.min(radii[:registry_count])) if registry_count else inlet.radius_min_um

    # Compare the largest internal movement with grid size, bubble radius, and
    # collision scales. These ratios help the user judge whether smaller time
    # steps are needed; they do not alter an already completed trajectory.
    step_resolution = build_particle_step_resolution(
        diagnostics.maximum_step_displacement_um, float(domain.spacing_um), minimum_radius,
        float(dynamics_cfg.collision_layer_um), bool(dynamics_cfg.collisions_enabled), float(time_plan.internal_dt_s),
        float(dynamics_cfg.collision_relaxation_time_s),
    )

    # Finish the real computing-time measurement and record the diameter range
    # that actually appeared, which may differ slightly from the allowed range.
    transport_seconds       = perf_counter() - started
    sampled_diameter_min    = (float(np.min(registry_diameter_um)) if registry_diameter_um.size else float("nan"))
    sampled_diameter_max    = (float(np.max(registry_diameter_um)) if registry_diameter_um.size else float("nan"))

    # Start molecular quality checks at zero. They remain zero when bonds are off.
    molecular_mean_field_low_ligand_observations        = 0
    molecular_mean_field_low_target_observations        = 0
    molecular_maximum_reaction_radius_curvature_proxy   = 0.0
    molecular_capture_distance_to_minimum_radius_ratio  = 0.0
    if molecular_enabled:
        # Gather the wall area close enough to each saved bubble to take part in bonding.
        area_records_um2        = _concat(records.target_reaction_area_um2, (0,), np.float64)
        positive_reaction_area  = area_records_um2 > 0.0

        # Convert density times area into an effective number of available molecules.
        # The factor 1e-12 changes square micrometres into square metres.
        ligand_counts = float(molecular_parameters.ligand_density_molecules_m2) * area_records_um2 * 1.0e-12
        target_counts = float(molecular_parameters.target_density_molecules_m2) * area_records_um2 * 1.0e-12
        warning_count = float(context.molecular_mean_field_warning_count)

        # Count observations where the average-population bond model may be weak
        # because too few ligand or target molecules are present in the contact area.
        molecular_mean_field_low_ligand_observations = int(np.count_nonzero(positive_reaction_area & (ligand_counts < warning_count)))
        molecular_mean_field_low_target_observations = int(np.count_nonzero(positive_reaction_area & (target_counts < warning_count)))

        if registry_count:
            # Compare capture distance with the smallest bubble radius. A large
            # value warns that the small-contact-area picture may not be suitable.
            molecular_capture_distance_to_minimum_radius_ratio = float(context.molecular_capture_distance_um / minimum_radius)

        if flat_ids.size:
            # Calculate the largest reaction-circle size relative to local vessel
            # radius. This is another check of the locally flat wall assumption.
            record_radii_um     = radii[flat_ids]
            reaction_radii_um   = reaction_disk_radius_um(record_radii_um, flat_gap, float(context.molecular_capture_distance_um))
            local_radius_um     = sample_bilinear(context.local_vessel_radius_um, flat_position_grid)

            molecular_maximum_reaction_radius_curvature_proxy = float(
                np.max(
                    np.divide(reaction_radii_um, local_radius_um, out=np.zeros_like(reaction_radii_um), where=local_radius_um > 0.0)
                )
            )
            
    # Build a plain text-friendly record of how this trajectory was produced.
    # Metadata lets a future reader understand and repeat the run without guessing
    # which optional physics or storage rules were active.
    metadata = {
        # Name the file layout and the two main movement choices first.
        "trajectory_schema": (
            "continuous_perfusion_revised_v20_topological_molecular_rbc_drift_diffusion_records_v14"
            if molecular_enabled and red_blood_cell_enabled
            else (
                "continuous_perfusion_revised_v20_topological_rbc_drift_diffusion_records_v14"
                if red_blood_cell_enabled
                else (
                    "continuous_perfusion_revised_v20_topological_molecular_records_v13"
                    if molecular_enabled
                    else "continuous_perfusion_revised_v20_topological_records_v13"
                )
            )
        ),
        "motion_model": "mobility",
        **topological_ownership.to_metadata(),
        "topological_transition_count": diagnostics.topological_transition_count,
        "wall_contact_integrator": (
            "revised_v16_continuous_geometry_predictive_mobility_unilateral_single_wall"
        ),
        "maximum_simultaneous_wall_constraints": 1,
        "perfusion_model": (
            "deterministic_pulsatile_flux_halton_empty_start_v7"
            if cardiac is not None
            else "deterministic_equal_flux_halton_empty_start_v6"
        ),
        # Save both the public frame spacing and the smaller internal step plan.
        "time_integrator": str(dynamics_cfg.time_integrator),
        "dt_s": float(time_plan.output_dt_s),
        "output_dt_s": float(time_plan.output_dt_s),
        "integration_substeps": int(time_plan.integration_substeps),
        "internal_dt_s": float(time_plan.internal_dt_s),
        "total_internal_integration_steps": int(time_plan.total_internal_steps),
        "expected_rhs_evaluations": int(time_plan.expected_rhs_evaluations),
        "n_steps": int(particle_cfg.n_steps),
        "stored_frames": int(time_plan.stored_frames),
        # Record how concentration and usable inlet flow became an entry rate.
        "inlet_number_concentration_mb_per_ml": inlet.number_concentration_mb_per_ml,
        "inlet_number_concentration_mb_per_m3": inlet.number_concentration_mb_per_m3,
        "inlet_number_concentration_bubbles_per_m3": inlet.number_concentration_mb_per_m3,
        "injection_rate_bubbles_per_s": inlet.injection_rate_per_s,
        "mean_planned_injection_interval_s": inlet.mean_injection_interval_s,
        "effective_thickness_um": inlet.effective_thickness_um,
        "inlet_section_flow_um2_s": inlet.raw_section_flow_um2_s,
        "size_accessible_inlet_flow_um2_s": inlet.size_accessible_section_flow_um2_s,
        "reference_inlet_flow_um2_s": inlet.reference_inlet_flow_um2_s,
        "inlet_section_flux_relative_error": inlet.section_flux_relative_error,
        # State clearly that the recording begins before the first bubble enters.
        "initial_condition": "empty_lumen_at_formal_time_zero",
        "active_bubbles_at_formal_time_zero": 0,
        # Summarize the changing population and any delayed entry events.
        "minimum_active_bubbles": int(np.min(active_counts)),
        "maximum_active_bubbles": int(np.max(active_counts)),
        "mean_active_bubbles": float(np.mean(active_counts)),
        "unique_bubbles_created": registry_count,
        "formal_admissions": int(np.sum(injected_counts)),
        "formal_terminations": int(np.sum(terminated_counts)),
        "inlet_wait_events": diagnostics.inlet_wait_events,
        "maximum_inlet_wait_s": diagnostics.maximum_inlet_wait_s,
        "waiting_bubbles_at_end": len(state.waiting_ids),
        # Confirm that one storage position always belonged to one permanent ID.
        "state_storage_slots_reused": False,
        # Save both the allowed size range and the range that actually appeared.
        "bubble_diameter_min_um": 2.0 * inlet.radius_min_um,
        "bubble_diameter_max_um": 2.0 * inlet.radius_max_um,
        "bubble_diameter_sample_min_um": sampled_diameter_min,
        "bubble_diameter_sample_max_um": sampled_diameter_max,
        # Record the calculation engine and measured computing time.
        "acceleration_backend": f"{backend}_synchronous_mobility",
        "particle_transport_seconds": transport_seconds,
        "particle_internal_steps_per_wall_second": (
            float(time_plan.total_internal_steps) / max(transport_seconds, 1.0e-30)
        ),
        "particle_numba_batched_saved_frame_count": int(
            batched_saved_frame_count
        ),
        # The unregularized gap is the only geometry/capture gap.  xi_min
        # changes hydrodynamic coefficients only and never moves the wall.
        "particle_true_gap_definition": (
            "g_R=distance_to_continuous_closed_vessel_boundary-radius_um"
        ),
        "particle_hydrodynamic_gap_definition": (
            "max(g_R,radius_um*xi_min)_mobility_coefficients_only"
        ),
        "near_wall_xi_min": float(dynamics_cfg.xi_min),
        "near_wall_xi_near": float(dynamics_cfg.xi_near),
        "near_wall_xi_far": float(dynamics_cfg.xi_far),
        "minimum_hydrodynamic_regularization_gap_um": (
            float(inlet.radius_min_um) * float(dynamics_cfg.xi_min)
        ),
        "maximum_hydrodynamic_regularization_gap_um": (
            float(inlet.radius_max_um) * float(dynamics_cfg.xi_min)
        ),
        # Store the single feasible-side geometry tolerance and the diagnostics
        # produced by the mobility-based unilateral contact constraint.
        "contact_geometry_tolerance_um": float(
            particle_cfg.contact_geometry_tolerance_um
        ),
        "contact_max_time_refinements": int(
            particle_cfg.contact_max_time_refinements
        ),
        "wall_contact_threshold_um": float(particle_cfg.wall_contact_threshold_um),
        "minimum_wall_gap_um": minimum_gap,
        "contact_observations": int(
            np.count_nonzero(_concat(records.wall_contact, (0,), bool))
        ),
        "contact_constraint_evaluations": diagnostics.contact_constraint_evaluations,
        "active_contact_constraint_evaluations": (
            diagnostics.active_contact_constraint_evaluations
        ),
        "maximum_contact_reaction_force_pn": (
            diagnostics.maximum_contact_reaction_force_pn
        ),
        "contact_time_refinement_count": diagnostics.contact_time_refinement_count,
        "maximum_contact_time_refinement_depth": (
            diagnostics.maximum_contact_time_refinement_depth
        ),
        "contact_residual_projection_count": diagnostics.contact_residual_projection_count,
        "maximum_contact_residual_projection_um": diagnostics.maximum_contact_residual_projection_um,
        "maximum_contact_complementarity_residual_pn_um": diagnostics.maximum_contact_complementarity_residual_pn_um,
        "minimum_accepted_internal_wall_gap_um": (
            diagnostics.minimum_accepted_internal_wall_gap_um
            if math.isfinite(diagnostics.minimum_accepted_internal_wall_gap_um)
            else float("nan")
        ),
        "accepted_negative_gap_count": 0,
        "contact_nonzero_velocity_zero_progress_count": (
            diagnostics.contact_nonzero_velocity_zero_progress_count
        ),
        "contact_kinematic_interval_evaluations": (
            diagnostics.contact_kinematic_interval_evaluations
        ),
        "contact_cumulative_position_path_um": (
            diagnostics.contact_cumulative_position_path_um
        ),
        "contact_cumulative_velocity_path_um": (
            diagnostics.contact_cumulative_velocity_path_um
        ),
        "contact_position_to_velocity_path_ratio": (
            diagnostics.contact_cumulative_position_path_um
            / diagnostics.contact_cumulative_velocity_path_um
            if diagnostics.contact_cumulative_velocity_path_um > 0.0
            else float("nan")
        ),
        "minimum_contact_interval_position_to_velocity_path_ratio": (
            diagnostics.minimum_contact_interval_position_to_velocity_path_ratio
            if math.isfinite(
                diagnostics.minimum_contact_interval_position_to_velocity_path_ratio
            )
            else float("nan")
        ),
        "maximum_contact_interval_position_to_velocity_path_ratio": (
            diagnostics.maximum_contact_interval_position_to_velocity_path_ratio
            if math.isfinite(
                diagnostics.maximum_contact_interval_position_to_velocity_path_ratio
            )
            else float("nan")
        ),
        "maximum_free_gap_kinematic_residual_um": (
            diagnostics.maximum_free_gap_kinematic_residual_um
        ),
        "directed_outlet_event_count": diagnostics.outlet_event_count,
        "active_outside_lumen_violations": diagnostics.active_outside_lumen_violations,
        "active_outside_accessible_domain_violations": diagnostics.active_outside_accessible_domain_violations,
        # The widest-path field is a nearest-node admission model.  A safe
        # continuous position may disagree with it at a half-grid boundary; such
        # records are audited here without rejecting a certified swept path.
        "discrete_accessibility_disagreement_records": (
            diagnostics.discrete_accessibility_disagreement_records
        ),
        # Store the switches and scales used for near-wall and collision movement.
        "near_wall_enabled": bool(dynamics_cfg.near_wall_enabled),
        "collisions_enabled": bool(dynamics_cfg.collisions_enabled),
        "collision_layer_um": float(dynamics_cfg.collision_layer_um),
        "collision_relaxation_time_s": float(dynamics_cfg.collision_relaxation_time_s),
        "collision_neighbor_search_configured": str(dynamics_cfg.neighbor_search),
        "collision_all_pairs_rhs_evaluations": diagnostics.all_pairs_rhs_evaluations,
        "collision_cell_list_rhs_evaluations": diagnostics.cell_list_rhs_evaluations,
        # Keep maximum overlap, collision, wall-normal, and opposite-wall
        # hydrodynamic-validity observations.  This last value is only a warning
        # about the one-planar-wall mobility approximation; true simultaneous
        # contact with distinct walls is detected separately and stops the run.
        # These are quality checks and do not add an extra force after the run.
        "maximum_physical_overlap_um": diagnostics.maximum_physical_overlap_um,
        "maximum_collision_compression_um": diagnostics.maximum_collision_compression_um,
        "maximum_interacting_pairs": diagnostics.maximum_interacting_pairs,
        "total_interacting_pair_evaluations": diagnostics.total_interacting_pair_evaluations,
        "maximum_mobility_reciprocity_relative_error": (
            diagnostics.maximum_reciprocity_relative_error
        ),
        "degenerate_near_wall_normal_evaluations": (
            diagnostics.degenerate_near_wall_normal_evaluations
        ),
        "maximum_collision_speed_um_s": diagnostics.maximum_collision_speed_um_s,
        "opposite_wall_hydrodynamic_validity_warning_observations": (
            diagnostics.two_wall_observations
        ),
        # Save step-size ratios used to decide whether a smaller internal time step
        # should be tested in a follow-up convergence study.
        "maximum_internal_step_displacement_um": diagnostics.maximum_step_displacement_um,
        "maximum_internal_step_grid_ratio": step_resolution.grid_displacement_ratio,
        "maximum_internal_step_radius_ratio": step_resolution.radius_displacement_ratio,
        "maximum_internal_step_field_size_ratio": step_resolution.field_size_displacement_ratio,
        "conservative_absolute_collision_layer_ratio": (
            step_resolution.conservative_collision_layer_ratio
        ),
        "collision_dt_over_relaxation_time": step_resolution.collision_dt_over_relaxation_time,
        # Explain what the saved particle velocity contains so it is not mistaken
        # for plain CFD fluid velocity or for a simple frame-to-frame difference.
        "particle_velocity_semantics": (
            (
                (
                    "last_accepted_internal_v16_continuous_geometry_predictive_velocity_from_retarded_phase_cardiac_background_plus_mobility_collision_molecular_bond_and_unilateral_contact"
                    if cardiac is not None
                    else "last_accepted_internal_v16_continuous_geometry_predictive_velocity_from_background_hydrodynamics_plus_mobility_collision_molecular_bond_and_unilateral_contact"
                )
                if molecular_enabled
                else (
                    "last_accepted_internal_v16_continuous_geometry_predictive_velocity_from_retarded_phase_cardiac_background_plus_mobility_collision_and_unilateral_contact"
                    if cardiac is not None
                    else "last_accepted_internal_v16_continuous_geometry_predictive_velocity_from_background_hydrodynamics_plus_mobility_collision_and_unilateral_contact"
                )
            )
            + (
                "_plus_rbc_drift_and_fick_translation"
                if red_blood_cell_enabled
                else ""
            )
        ),
        "realized_particle_velocity_semantics": (
            "same_permanent_id_geometry_accepted_saved_frame_center_displacement_including_rbc_stochastic_geometry_when_enabled"
        ),
        # Explain the time origin and which optional output groups were active.
        "saved_frame_time_semantics": "formal time begins with an empty lumen at t=0",
        "scheduled_injection_time_semantics": "time relative to empty-lumen formal t=0",
        "store_full_diagnostics": full_diagnostics,
        "cardiac_pulsatility_enabled": cardiac is not None,
        "molecular_binding_enabled": molecular_enabled,
        "red_blood_cell_transport_enabled": red_blood_cell_enabled,
        "red_blood_cell_transport_model": (
            "RBC-induced drift–diffusion reduced-order transport model"
            if red_blood_cell_enabled
            else "disabled"
        ),
        "red_blood_cell_root_discharge_hematocrit": (
            float(context.red_blood_cell_network.root_discharge_hematocrit)
            if red_blood_cell_enabled
            else 0.0
        ),
        "red_blood_cell_effective_diameter_um": (
            float(context.red_blood_cell_network.effective_rbc_diameter_um)
            if red_blood_cell_enabled
            else RBC_MAJOR_DIAMETER_UM
        ),
        "red_blood_cell_major_diameter_um": RBC_MAJOR_DIAMETER_UM,
        "red_blood_cell_cfl_model": (
            "thg_3pef_cortical_mouse_linear_to_plateau_v1"
        ),
        "red_blood_cell_cfl_plateau_diameter_um": CFL_PLATEAU_DIAMETER_UM,
        "red_blood_cell_cfl_plateau_thickness_um": CFL_PLATEAU_THICKNESS_UM,
        "red_blood_cell_cfl_measured_diameter_range_um": (
            f"{CFL_DIAMETER_MIN_UM:g}..{CFL_DIAMETER_MAX_UM:g}"
        ),
        "red_blood_cell_coupling_scope": (
            "one_way_particle_translation_only_without_cfd_viscosity_mobility_torque_or_molecular_parameter_changes"
        ),
        "red_blood_cell_shear_rate_definition": "sqrt(2*E:E)_from_cycle_mean_CFD_velocity_gradient_without_cardiac_spatial_modulation",
        "red_blood_cell_compatibility_margination_length_semantics": (
            "g_center_minus_g_star_diagnostic_only"
        ),
        "red_blood_cell_compatibility_margination_time_semantics": (
            "Gamma_m_over_eta_gamma_diagnostic_only_never_used_to_compute_drift_or_diffusion"
        ),
        "red_blood_cell_quantitative_tube_hematocrit_range": "0.15..0.30",
        "red_blood_cell_quantitative_shear_rate_range_s_inv": "100..10000",
        "red_blood_cell_reference_cumulative_shear_strain": RBC_MARGINATION_STRAIN,
        "red_blood_cell_reference_transverse_diffusivity_um2_s": (
            REFERENCE_TRANSVERSE_DIFFUSIVITY_UM2_S
        ),
        "red_blood_cell_reference_tube_hematocrit": REFERENCE_TUBE_HEMATOCRIT,
        "red_blood_cell_reference_shear_rate_s_inv": (
            REFERENCE_SHEAR_RATE_S_INV
        ),
        "red_blood_cell_rng_seed": int(context.random_seed),
        "red_blood_cell_rng_algorithm_version": COUNTER_RNG_ALGORITHM,
        "red_blood_cell_random_wall_reflection_count": (
            diagnostics.red_blood_cell_random_wall_reflection_count
        ),
        "red_blood_cell_maximum_random_displacement_um": (
            diagnostics.maximum_red_blood_cell_random_displacement_um
        ),
        "red_blood_cell_rms_random_displacement_um": (
            red_blood_cell_rms_random_displacement_um
        ),
        "red_blood_cell_random_displacement_observation_count": (
            diagnostics.red_blood_cell_random_displacement_observations
        ),
        "red_blood_cell_diffusion_enabled_record_fraction": (
            red_blood_cell_diffusion_enabled_record_fraction
        ),
        "red_blood_cell_diffusion_enabled_internal_fraction": (
            diagnostics.red_blood_cell_diffusion_enabled_observations
            / diagnostics.red_blood_cell_stochastic_observations
            if diagnostics.red_blood_cell_stochastic_observations > 0
            else 0.0
        ),
        "red_blood_cell_invalid_transverse_space_observations": (
            diagnostics.red_blood_cell_invalid_transverse_space_observations
        ),
        "red_blood_cell_transverse_space_valid_semantics": (
            "g_center_greater_than_g_star_and_g_center_greater_than_zero"
        ),
        "red_blood_cell_saved_frame_applicability_semantics": (
            "quantitative_range_and_transverse_space_valid_are_separate_flags"
        ),
        "red_blood_cell_maximum_speed_um_s": (
            diagnostics.maximum_red_blood_cell_speed_um_s
        ),
        "red_blood_cell_rhs_observations": diagnostics.red_blood_cell_rhs_observations,
        "red_blood_cell_nonunique_wall_suppressions": (
            diagnostics.red_blood_cell_nonunique_wall_suppressions
        ),
        "red_blood_cell_hematocrit_out_of_range_rhs_observations": (
            diagnostics.red_blood_cell_hematocrit_out_of_range_observations
        ),
        "red_blood_cell_shear_out_of_range_rhs_observations": (
            diagnostics.red_blood_cell_shear_out_of_range_observations
        ),
        "red_blood_cell_quantitative_applicability_rhs_observations": (
            diagnostics.red_blood_cell_quantitative_applicability_observations
        ),
        "red_blood_cell_quantitative_applicability_record_fraction": (
            red_blood_cell_quantitative_record_fraction
        ),
        "red_blood_cell_exposure_weighted_quantitative_applicability_fraction": (
            red_blood_cell_exposure_weighted_applicability_fraction
        ),
        "red_blood_cell_exposure_applicability_semantics": (
            "target_exposure_time_weighted_diameter_H_t_and_shear_range_flag_at_nominal_internal_step_start_using_continuous_topological_owner"
        ),
        "red_blood_cell_cfl_record_fraction": red_blood_cell_cfl_record_fraction,
        "red_blood_cell_cfl_reaching_bubble_fraction": (
            red_blood_cell_cfl_reaching_bubble_fraction
        ),
        # Describe the target field even when bonding is off, because a geometry-only
        # target can still be useful for contact studies.
        "molecular_target_region_mode": (
            str(molecular_target_field.region_mode)
            if molecular_target_field is not None
            else "disabled"
        ),
        "molecular_target_density_molecules_per_m2": (
            float(molecular_target_field.target_density_molecules_per_m2)
            if molecular_target_field is not None
            else 0.0
        ),
        "target_internal_exposure_enabled": bool(target_exposure_enabled),
        "target_internal_exposure_integration": (
            "accepted_deterministic_transaction_endpoints_plus_stochastic_geometry_event_vertices_and_final_endpoint_batched_bisection_v2_plus_eight_point_gauss_legendre_over_exposed_subintervals_with_one_shared_physical_dt"
            if target_exposure_enabled
            else "disabled"
        ),
        "target_exposure_N_enc": target_exposure_encountered_bubble_count,
        "target_exposure_encountered_bubble_count": (
            target_exposure_encountered_bubble_count
        ),
        "target_exposure_event_count": target_exposure_event_count,
        "target_exposure_total_time_s": target_exposure_total_time_s,
        "target_reaction_area_time_um2_s": target_reaction_area_time_um2_s,
        "target_exposure_right_censored_event_count": (
            target_exposure_right_censored_count
        ),
        # State that bond count is a continuous expected amount, not a yes/no flag
        # for one individually tracked molecular bond.
        "molecular_state_variables": "expected_bond_count,total_tangential_extension_um",
        "molecular_dynamics_semantics": (
            "deterministic_mean_field_bonds_without_binary_bound_state"
            if molecular_enabled
            else "disabled"
        ),
        # Save the largest molecular loads and counts seen during the run.
        "molecular_maximum_expected_bond_count": diagnostics.maximum_expected_bond_count,
        "molecular_maximum_single_bond_tension_pn": (
            diagnostics.maximum_single_bond_tension_pn
        ),
        "molecular_maximum_total_bond_force_pn": diagnostics.maximum_total_bond_force_pn,
        "molecular_maximum_bond_torque_pn_um": diagnostics.maximum_bond_torque_pn_um,
        "molecular_bell_saturation_rhs_observations": (
            diagnostics.molecular_bell_saturation_count
        ),
        "molecular_formation_rate_float_overflow_rhs_observations": (
            diagnostics.molecular_formation_saturation_count
        ),
        "molecular_capacity_limited_accepted_step_observations": (
            diagnostics.molecular_capacity_limited_accepted_step_count
        ),
        # Save geometry ratios and low-molecule warnings used to judge whether the
        # simplified local bond model was used outside its comfortable range.
        "molecular_capture_distance_to_minimum_radius_ratio": (
            molecular_capture_distance_to_minimum_radius_ratio
        ),
        "molecular_maximum_a_times_curvature_proxy": (
            molecular_maximum_reaction_radius_curvature_proxy
        ),
        "molecular_curvature_proxy_semantics": (
            "a/local_vessel_radius; conservative diagnostic, not an exact wall-curvature fit"
        ),
        "molecular_mean_field_warning_count": float(
            context.molecular_mean_field_warning_count
        ),
        "molecular_mean_field_warning_is_physical_threshold": False,
        "molecular_low_ligand_count_record_observations": (
            molecular_mean_field_low_ligand_observations
        ),
        "molecular_low_target_count_record_observations": (
            molecular_mean_field_low_target_observations
        ),
        # List deliberate engineering additions that are not part of the basic
        # single-bubble near-wall reference model.
        "engineering_extensions": (
            "whole-vessel wall/free-space blending; pair collision relaxation; "
            "event-split synchronous integration; deterministic upstream inlet waiting; "
            "position-level mobility-based unilateral wall contact; directed outlet planes; "
            "inlet-connected finite-radius accessibility; physical-time refinement"
        ),
    }
    metadata.update(backend_details)

    if molecular_enabled:
        # Add the exact bond settings only when molecular calculations were active.
        # Keeping disabled fields out of this section avoids making zeros look like
        # physical parameters that were actually used.
        metadata.update(
            {
                "molecular_binding_model": str(molecular_binding_cfg.model),
                "molecular_ligand_density_molecules_per_m2": float(
                    molecular_binding_cfg.ligand_density_molecules_per_m2
                ),
                "molecular_capture_distance_um": float(
                    molecular_binding_cfg.capture_distance_um
                ),
                "molecular_rest_length_um": float(molecular_binding_cfg.rest_length_um),
                "molecular_association_rate_m2_per_molecule_s": float(
                    molecular_binding_cfg.association_rate_m2_per_molecule_s
                ),
                "molecular_zero_force_dissociation_rate_s": float(
                    molecular_binding_cfg.zero_force_dissociation_rate_s
                ),
                "molecular_bond_stiffness_pn_per_um": float(
                    molecular_binding_cfg.bond_stiffness_pn_per_um
                ),
                "molecular_reactive_compliance_nm": float(
                    molecular_binding_cfg.reactive_compliance_nm
                ),
                "molecular_temperature_k": float(molecular_binding_cfg.temperature_k),
                "molecular_reaction_area_integration": (
                    "analytic_target_roi_intersected_with_eligible_closed_wall_intervals"
                ),
                "molecular_bell_exponent_limit": float(
                    molecular_binding_cfg.bell_exponent_limit
                ),
            }
        )

    if cardiac is not None:
        # Add the complete heartbeat description only when pulsation was active.
        # These values allow the local speed changes and entry timing to be repeated.
        metadata.update(
            {
                "cardiac_waveform": cardiac.waveform_name,
                "cardiac_waveform_semantics": "positive ECG-shaped surrogate for carrier-flow modulation",
                "cardiac_bpm": float(cardiac.waveform.bpm),
                "cardiac_period_s": float(cardiac.waveform.period_s),
                "cardiac_formal_cycle_count": float(record_duration / cardiac.waveform.period_s),
                "cardiac_phase_offset_s": float(cardiac.phase_offset_s),
                "cardiac_waveform_samples_per_cycle": int(cardiac.waveform.multiplier.size),
                "cardiac_cycle_mean_multiplier": float(cardiac.waveform.cycle_mean),
                "cardiac_minimum_multiplier": float(np.min(cardiac.waveform.multiplier)),
                "cardiac_maximum_multiplier": float(np.max(cardiac.waveform.multiplier)),
                "cardiac_minimum_injection_rate_bubbles_per_s": float(
                    inlet.injection_rate_per_s * np.min(cardiac.waveform.multiplier)
                ),
                "cardiac_maximum_injection_rate_bubbles_per_s": float(
                    inlet.injection_rate_per_s * np.max(cardiac.waveform.multiplier)
                ),
                "cardiac_pulse_propagation_velocity_um_s": float(
                    cardiac.pulse_propagation_velocity_um_s
                ),
                "cardiac_maximum_root_path_um": float(np.max(cardiac.path_distance_um)),
                "cardiac_maximum_propagation_delay_s": float(np.max(cardiac.delay_s)),
                "cardiac_preserve_cycle_mean_flow": bool(cardiac.preserve_cycle_mean_flow),
                "cardiac_modulation_strength": float(cardiac.modulation_strength),
                "background_flow_semantics": "accepted steady CFD field interpreted as cycle-mean reference",
                "cardiac_model_scope": "retarded-phase quasi-steady kinematic modulation; not compliant transient CFD",
            }
        )
        
    # Save the physical centre and inward direction of every inlet section.
    # Multiple roots are written separately so their geometry remains unambiguous.
    for section_index, section in enumerate(inlet.sections):
        metadata[f"inlet_section_{section_index}_vessel_id"] = int(section.vessel_id)
        metadata[f"inlet_section_{section_index}_center_x_um"] = float(section.center_xz_um[0])
        metadata[f"inlet_section_{section_index}_center_z_um"] = float(section.center_xz_um[1])
        metadata[f"inlet_section_{section_index}_normal_x"] = float(section.inward_normal_xz[0])
        metadata[f"inlet_section_{section_index}_normal_z"] = float(section.inward_normal_xz[1])

    # Print a short human-readable summary after all movement has finished.
    # These rows make population, waiting, time-step size, and computing cost easy
    # to inspect without opening the saved numerical file.
    print_section("Continuous perfusion transport summary")
    transport_rows: list[tuple[str, object]] = [
        ("Unique scheduled bubbles", registry_count),
        ("Active population range", f"{int(np.min(active_counts))}..{int(np.max(active_counts))}"),
        ("Mean active population", f"{float(np.mean(active_counts)):.6g}"),
        ("Inlet waiting events", diagnostics.inlet_wait_events),
        ("Maximum inlet wait", f"{diagnostics.maximum_inlet_wait_s:.6g} s"),
        ("Trajectory records", int(records.offsets[-1])),
        ("Particle backend", backend),
        ("Numeric kernels", backend_details["particle_numeric_kernel_family"]),
        (
            "Batched saved-frame intervals",
            f"{batched_saved_frame_count}/{time_plan.output_intervals}",
        ),
        ("Collision search", str(dynamics_cfg.neighbor_search)),
        (
            "Collision RHS paths",
            (
                f"all-pairs={diagnostics.all_pairs_rhs_evaluations}, "
                f"cell-list={diagnostics.cell_list_rhs_evaluations}"
            ),
        ),
        ("Maximum internal-step displacement", f"{diagnostics.maximum_step_displacement_um:.6g} um"),
        ("Internal-step / grid ratio", f"{step_resolution.grid_displacement_ratio:.6g}"),
        ("Internal-step / radius ratio", f"{step_resolution.radius_displacement_ratio:.6g}"),
        (
            "Active contact-constraint evaluations",
            diagnostics.active_contact_constraint_evaluations,
        ),
        (
            "Maximum contact reaction",
            f"{diagnostics.maximum_contact_reaction_force_pn:.6g} pN",
        ),
        (
            "Contact physical-time refinements",
            diagnostics.contact_time_refinement_count,
        ),
        (
            "Maximum contact-refinement depth",
            diagnostics.maximum_contact_time_refinement_depth,
        ),
    ]
    if backend == "numba_cpu":
        transport_rows.append(
            (
                "Numba particle workers",
                backend_details["particle_numba_worker_threads"],
            )
        )
    transport_rows.append(
        (
            "Continuous-wall bin width",
            f"{backend_details['particle_continuous_wall_numba_bin_width_cells']} grid cells",
        )
    )

    # The collision-layer comparison has meaning only when that model is active.
    if dynamics_cfg.collisions_enabled and dynamics_cfg.collision_layer_um > 0.0:
        transport_rows.append(
            (
                "Displacement / collision-layer ratio",
                f"{step_resolution.conservative_collision_layer_ratio:.6g}",
            )
        )

    # Always report whether the internal time step is small compared with the
    # chosen collision relaxation time.
    transport_rows.append(
        (
            "Internal dt / collision relaxation",
            f"{step_resolution.collision_dt_over_relaxation_time:.6g}",
        )
    )

    # Estimate how many frames the record safety limit could hold at the observed
    # peak population. This is a memory guide, not a physical simulation limit.
    record_limit = int(particle_cfg.max_particle_frame_records)
    if record_limit > 0:
        peak = max(int(np.max(active_counts)), 1)
        transport_rows.append(("Record-limit frame capacity", f">= {record_limit // peak}"))
    transport_rows.append((f"Particle transport time ({backend})", f"{transport_seconds:.3f} s"))
    transport_rows.append(
        ("Topological vessel transitions", diagnostics.topological_transition_count)
    )
    transport_rows.append(
        (
            "Internal steps / wall second",
            f"{float(time_plan.total_internal_steps) / max(transport_seconds, 1.0e-30):.6g}",
        )
    )
    if red_blood_cell_enabled:
        transport_rows.extend(
            [
                (
                    "Maximum RBC deterministic drift/Fick speed",
                    f"{diagnostics.maximum_red_blood_cell_speed_um_s:.6g} um/s",
                ),
                (
                    "Bubbles reaching the mean CFL",
                    f"{red_blood_cell_cfl_reaching_bubble_fraction:.6g}",
                ),
                (
                    "RBC quantitative-applicability records",
                    f"{red_blood_cell_quantitative_record_fraction:.6g}",
                ),
            ]
        )
    print_key_values(transport_rows)

    # Warn when planned bubbles had to wait outside because the inlet was occupied.
    if diagnostics.inlet_wait_events > 0:
        print_warning(
            "One or more deterministic inlet events had to wait upstream for "
            "a non-overlapping admission state."
        )

    # Warn where one-flat-wall hydrodynamics may be weak because a bubble is also
    # close to the opposite side of a narrow vessel.  This is not the true
    # simultaneous-contact error, which is rejected by the v15 geometry step.
    if diagnostics.two_wall_observations > 0:
        print_warning(
            "Some stored states are close to an opposite wall; the referenced "
            "single-planar-wall hydrodynamic mobility remains an approximation "
            "there. This warning is distinct from simultaneous wall contact."
        )

    if red_blood_cell_enabled and (
        diagnostics.red_blood_cell_hematocrit_out_of_range_observations > 0
        or diagnostics.red_blood_cell_shear_out_of_range_observations > 0
    ):
        print_warning(
            "The RBC-induced drift–diffusion closure was evaluated outside its quantitative "
            "tube-haematocrit or shear-rate range for some observations. Values "
            "were recorded without clipping or coefficient tuning."
        )

    # Warn when a wall direction could not be read cleanly and a repeatable backup
    # direction had to be used instead.
    if diagnostics.degenerate_near_wall_normal_evaluations > 0:
        print_warning(
            "Near-wall mobility encountered degenerate distance-field normals; "
            "a deterministic fallback basis was used and counted in metadata."
        )

    # Warn when one internal move is large compared with the grid or bubble size.
    # The user should then repeat the run with more internal steps and compare results.
    if step_resolution.field_size_displacement_ratio >= 0.2:
        print_warning(
            "An internal mobility step reached at least 0.2 of the grid/radius "
            "scale. Increase particle_dynamics.integration_substeps and repeat the "
            "substep-convergence study."
        )

    # Collision forces can also change too quickly when the internal time step is
    # not small compared with their relaxation time.
    if step_resolution.collision_dt_over_relaxation_time >= 0.2:
        print_warning(
            "The internal particle step reached at least 0.2 of the collision "
            "relaxation time. Increase particle_dynamics.integration_substeps while "
            "keeping collision_relaxation_time_s fixed."
        )

    # Molecular summaries and warnings are shown only when bond physics was active.
    if molecular_enabled:
        print_section("Molecular binding summary")
        print_key_values(
            [
                ("Maximum expected bond population", f"{diagnostics.maximum_expected_bond_count:.6g}"),
                ("Maximum single-bond tension", f"{diagnostics.maximum_single_bond_tension_pn:.6g} pN"),
                ("Maximum total bond torque", f"{diagnostics.maximum_bond_torque_pn_um:.6g} pN um"),
                (
                    "Capture distance / minimum radius",
                    f"{molecular_capture_distance_to_minimum_radius_ratio:.6g}",
                ),
                (
                    "Reaction radius / curvature proxy",
                    f"{molecular_maximum_reaction_radius_curvature_proxy:.6g}",
                ),
            ]
        )

        # Check whether the capture layer is small compared with bubble size.
        if molecular_capture_distance_to_minimum_radius_ratio >= 0.1:
            print_warning(
                "Molecular capture distance is not much smaller than the minimum "
                "bubble radius; review the local spherical-cap approximation."
            )

        # Check whether the reaction circle is small compared with the local vessel.
        if molecular_maximum_reaction_radius_curvature_proxy >= 0.1:
            print_warning(
                "A reaction disk reached at least 0.1 of the local vessel-radius "
                "curvature proxy; review the locally planar target-wall approximation."
            )

        # Report contact areas with too few molecules for a smooth average-count model.
        if (
            molecular_mean_field_low_ligand_observations > 0
            or molecular_mean_field_low_target_observations > 0
        ):
            print_warning(
                "Some positive reaction areas contained fewer effective ligand or "
                "target molecules than molecular_binding.mean_field_warning_count. This is "
                "an engineering validity warning, not a physical binding threshold."
            )

        # Report when the force-sensitive breaking rule produced a value that grew
        # too large and therefore reached its numerical safety cap.
        if diagnostics.molecular_bell_saturation_count > 0:
            print_warning(
                "The Bell exponent reached molecular_binding.bell_exponent_limit; "
                "the configured overflow guard was applied and counted in metadata."
            )

        # Report when an accepted bond step had to be kept within available molecule counts.
        if diagnostics.molecular_capacity_limited_accepted_step_count > 0:
            print_warning(
                "A molecular accepted step required the ligand/target capacity "
                "limiter. Increase particle_dynamics.integration_substeps and repeat the "
                "bond-state convergence study."
            )

    # Pack all public results into one trajectory object for the file writer.
    # The first group describes every active observation in the flat frame records.
    return ParticleTrajectories(
        frame_offsets=flat_frame_offsets,
        bubble_id=flat_ids,
        positions_um=positions_world,
        velocities_um_s=particle_velocity_world,
        realized_velocities_um_s=realized_velocity_world,
        wall_shear_stress_pa=_concat(records.wall_shear, (0,), np.float32),
        vessel_id=_concat(records.vessel_id, (0,), np.int32),
        active=flat_active,
        diameter_um=registry_diameter_um[flat_ids] if flat_ids.size else np.empty(0, dtype=np.float32),
        wall_gap_um=np.asarray(flat_gap, dtype=np.float32),
        wall_contact=_concat(records.wall_contact, (0,), bool),
        wall_normal_xz=_concat(records.wall_normal, (0, 2), np.float32),
        # This group is the permanent-ID master list and each bubble's lifetime.
        registry_bubble_id=registry_ids,
        registry_diameter_um=registry_diameter_um,
        birth_frame=birth_frame,
        death_frame=death_frame,
        termination_reason=np.asarray(state.termination_reason[:registry_count], dtype=np.uint8),
        # These arrays give a quick population summary for each saved frame.
        active_count_per_frame=active_counts,
        injected_count_per_frame=injected_counts,
        terminated_count_per_frame=terminated_counts,
        metadata=metadata,
        # Detailed flow, rotation, collision, and near-wall arrays are included
        # only when full diagnostics were requested; otherwise they remain None.
        fluid_velocities_um_s=fluid_velocity_world,
        angular_velocity_rad_s=(
            _concat(records.angular_velocity, (0,), np.float32) if full_diagnostics else None
        ),
        rotation_angle_rad=(
            _concat(records.rotation_angle, (0,), np.float32) if full_diagnostics else None
        ),
        collision_force_xz_pn=(
            _concat(records.collision_force, (0, 2), np.float32)
            if full_diagnostics
            else None
        ),
        collision_neighbor_count=(
            _concat(records.collision_neighbors, (0,), np.int32)
            if full_diagnostics
            else None
        ),
        gap_ratio=_concat(records.gap_ratio, (0,), np.float32) if full_diagnostics else None,
        near_wall_weight=(
            _concat(records.near_wall_weight, (0,), np.float32)
            if full_diagnostics
            else None
        ),
        two_wall_warning=(
            _concat(records.two_wall_warning, (0,), bool) if full_diagnostics else None
        ),
        contact_constraint_active=flat_contact_constraint_active,
        contact_reaction_force_pn=np.asarray(
            flat_contact_reaction_force_pn, dtype=np.float32
        ),
        contact_free_normal_velocity_um_s=np.asarray(
            flat_contact_free_normal_velocity_um_s, dtype=np.float32
        ),
        contact_constrained_normal_velocity_um_s=np.asarray(
            flat_contact_constrained_normal_velocity_um_s, dtype=np.float32
        ),
        # Save planned and real entry/exit times so inlet waiting can be measured
        # without guessing from the nearest stored frame.
        registry_scheduled_injection_time_s=np.asarray(
            schedule.planned_time_s[:registry_count], dtype=np.float64
        ),
        registry_admission_time_s=np.asarray(relative_admission, dtype=np.float64),
        registry_exit_time_s=np.asarray(relative_exit, dtype=np.float64),
        registry_inlet_wait_time_s=np.asarray(
            state.admission_time_s[:registry_count] - schedule.planned_time_s[:registry_count],
            dtype=np.float64,
        ),
        # Save the heartbeat curve and spatial delay maps only when pulsation was used.
        cardiac_multiplier=cardiac_multiplier,
        cardiac_waveform_time_s=(
            np.asarray(cardiac.waveform.sample_time_s, dtype=np.float64)
            if cardiac is not None
            else None
        ),
        cardiac_waveform_multiplier=(
            np.asarray(cardiac.waveform.multiplier, dtype=np.float64)
            if cardiac is not None
            else None
        ),
        cardiac_path_distance_um=(
            np.asarray(cardiac.path_distance_um, dtype=np.float32)
            if cardiac is not None
            else None
        ),
        cardiac_delay_s=(
            np.asarray(cardiac.delay_s, dtype=np.float32) if cardiac is not None else None
        ),
        # Save expected bond population and average bond stretching when bonding was active.
        bond_count_expected=(
            _concat(records.bond_count_expected, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        bond_total_tangential_extension_um=(
            _concat(records.bond_total_tangential_extension_um, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        bond_mean_tangential_extension_um=(
            _concat(records.bond_mean_tangential_extension_um, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        # Save molecular force, torque, and single-bond tension observations.
        bond_force_xz_pn=(
            _concat(records.bond_force_xz_pn, (0, 2), np.float64)
            if molecular_enabled
            else None
        ),
        bond_force_tangent_pn=(
            _concat(records.bond_force_tangent_pn, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        bond_force_normal_pn=(
            _concat(records.bond_force_normal_pn, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        bond_torque_pn_um=(
            _concat(records.bond_torque_pn_um, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        single_bond_tension_pn=(
            _concat(records.single_bond_tension_pn, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        # Save bond formation and breaking rates at every active observation.
        bond_formation_rate_bonds_s=(
            _concat(records.bond_formation_rate_bonds_s, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        bond_dissociation_rate_s_inv=(
            _concat(records.bond_dissociation_rate_s_inv, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        # Save molecular contact geometry and the available molecule counts.
        target_reaction_area_um2=(
            _concat(records.target_reaction_area_um2, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        available_ligand_count=(
            _concat(records.available_ligand_count, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        available_target_count=(
            _concat(records.available_target_count, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        target_overlap_fraction=(
            _concat(records.target_overlap_fraction, (0,), np.float64)
            if molecular_enabled
            else None
        ),
        # Keep the final bond state for every permanent ID, including bubbles that
        # are no longer present in the last saved frame.
        registry_final_bond_count_expected=(
            np.asarray(state.bond_count_expected[:registry_count], dtype=np.float64)
            if molecular_enabled
            else None
        ),
        registry_final_bond_total_tangential_extension_um=(
            np.asarray(
                state.bond_total_tangential_extension_um[:registry_count], dtype=np.float64
            )
            if molecular_enabled
            else None
        ),
        registry_target_exposure_time_s=(
            registry_target_exposure_time_s if target_exposure_enabled else None
        ),
        registry_target_exposure_event_count=(
            registry_target_exposure_event_count
            if target_exposure_enabled
            else None
        ),
        registry_target_reaction_area_time_um2_s=(
            registry_target_reaction_area_time_um2_s
            if target_exposure_enabled
            else None
        ),
        registry_target_exposure_right_censored=(
            registry_target_exposure_right_censored
            if target_exposure_enabled
            else None
        ),
        registry_target_exposure_quantitative_applicability_fraction=(
            registry_target_exposure_applicability_fraction
            if target_exposure_enabled
            else None
        ),
        red_blood_cell_velocity_xz_um_s=(
            np.asarray(flat_red_blood_cell_velocity_xz_um_s, dtype=np.float32)
            if red_blood_cell_enabled
            else None
        ),
        red_blood_cell_drift_velocity_xz_um_s=(
            _concat(
                records.red_blood_cell_drift_velocity_xz_um_s,
                (0, 2),
                np.float32,
            )
            if red_blood_cell_enabled
            else None
        ),
        red_blood_cell_fick_velocity_xz_um_s=(
            _concat(
                records.red_blood_cell_fick_velocity_xz_um_s,
                (0, 2),
                np.float32,
            )
            if red_blood_cell_enabled
            else None
        ),
        red_blood_cell_local_vessel_diameter_um=(
            _concat(records.red_blood_cell_local_vessel_diameter_um, (0,), np.float32)
            if red_blood_cell_enabled and full_diagnostics
            else None
        ),
        red_blood_cell_discharge_hematocrit=(
            _concat(records.red_blood_cell_discharge_hematocrit, (0,), np.float32)
            if red_blood_cell_enabled and full_diagnostics
            else None
        ),
        red_blood_cell_tube_hematocrit=(
            _concat(records.red_blood_cell_tube_hematocrit, (0,), np.float32)
            if red_blood_cell_enabled and full_diagnostics
            else None
        ),
        red_blood_cell_shear_rate_s_inv=(
            _concat(records.red_blood_cell_shear_rate_s_inv, (0,), np.float32)
            if red_blood_cell_enabled and full_diagnostics
            else None
        ),
        red_blood_cell_cfl_width_um=(
            _concat(records.red_blood_cell_cfl_width_um, (0,), np.float32)
            if red_blood_cell_enabled and full_diagnostics
            else None
        ),
        red_blood_cell_target_gap_um=(
            np.asarray(flat_red_blood_cell_target_gap_um, dtype=np.float32)
            if red_blood_cell_enabled and full_diagnostics
            else None
        ),
        red_blood_cell_transverse_diffusivity_um2_s=(
            _concat(
                records.red_blood_cell_transverse_diffusivity_um2_s,
                (0,),
                np.float32,
            )
            if red_blood_cell_enabled
            else None
        ),
        red_blood_cell_margination_length_um=(
            _concat(records.red_blood_cell_margination_length_um, (0,), np.float32)
            if red_blood_cell_enabled and full_diagnostics
            else None
        ),
        red_blood_cell_margination_time_s=(
            _concat(records.red_blood_cell_margination_time_s, (0,), np.float32)
            if red_blood_cell_enabled and full_diagnostics
            else None
        ),
        red_blood_cell_scale_activation=(
            _concat(records.red_blood_cell_scale_activation, (0,), np.float32)
            if red_blood_cell_enabled and full_diagnostics
            else None
        ),
        red_blood_cell_nearest_wall_unique=(
            _concat(records.red_blood_cell_nearest_wall_unique, (0,), bool)
            if red_blood_cell_enabled and full_diagnostics
            else None
        ),
        red_blood_cell_hematocrit_in_quantitative_range=(
            _concat(
                records.red_blood_cell_hematocrit_in_quantitative_range,
                (0,),
                bool,
            )
            if red_blood_cell_enabled and full_diagnostics
            else None
        ),
        red_blood_cell_shear_rate_in_quantitative_range=(
            _concat(
                records.red_blood_cell_shear_rate_in_quantitative_range,
                (0,),
                bool,
            )
            if red_blood_cell_enabled and full_diagnostics
            else None
        ),
        red_blood_cell_quantitative_applicability=(
            _concat(records.red_blood_cell_quantitative_applicability, (0,), bool)
            if red_blood_cell_enabled
            else None
        ),
        red_blood_cell_transverse_space_valid=(
            _concat(records.red_blood_cell_transverse_space_valid, (0,), bool)
            if red_blood_cell_enabled
            else None
        ),
        registry_final_vessel_id=np.asarray(
            state.vessel_id[:registry_count], dtype=np.int32
        ),
        topological_commitment_parent_vessel_id=np.asarray(
            topological_ownership.parent_vessel_id, dtype=np.int32
        ),
        topological_commitment_child_vessel_id=np.asarray(
            topological_ownership.child_vessel_id, dtype=np.int32
        ),
        topological_commitment_point_xz_um=np.asarray(
            topological_ownership.point_xz_um, dtype=np.float64
        ),
        topological_commitment_downstream_normal_xz=np.asarray(
            topological_ownership.downstream_normal_xz, dtype=np.float64
        ),
        topological_commitment_tangent_xz=np.asarray(
            topological_ownership.tangent_xz, dtype=np.float64
        ),
        topological_commitment_half_width_um=np.asarray(
            topological_ownership.half_width_um, dtype=np.float64
        ),
        topological_commitment_transition_end_distance_um=np.asarray(
            topological_ownership.transition_end_distance_um, dtype=np.float64
        ),
        topological_commitment_distance_um=np.asarray(
            topological_ownership.commitment_distance_um, dtype=np.float64
        ),
        topological_event_bubble_id=np.asarray(
            diagnostics.topological_event_bubble_id, dtype=np.int64
        ),
        topological_event_time_s=np.asarray(
            diagnostics.topological_event_time_s, dtype=np.float64
        ),
        topological_event_from_vessel_id=np.asarray(
            diagnostics.topological_event_from_vessel_id, dtype=np.int32
        ),
        topological_event_to_vessel_id=np.asarray(
            diagnostics.topological_event_to_vessel_id, dtype=np.int32
        ),
        topological_event_section_index=np.asarray(
            diagnostics.topological_event_section_index, dtype=np.int32
        ),
        topological_event_position_xz_um=np.asarray(
            diagnostics.topological_event_position_xz_um, dtype=np.float64
        ).reshape(-1, 2),
    )


def _advance_interval(
    start_time: float,
    end_time: float,
    state: _PerfusionState,
    schedule: PerfusionSchedule,
    context: object,
    particle_cfg: ParticleConfig,
    dynamics_cfg: ParticleDynamicsConfig,
    diagnostics: _Diagnostics,
    evaluate_rhs: object,
    accepted_target_path_batches: list[tuple[
        np.ndarray, tuple[np.ndarray, ...], tuple[np.ndarray, ...]
    ]] | None = None,
) -> tuple[int, int]:
    # ==================================================================================
    # ===== Prepare time counter and population change counters for this interval.
    # ==================================================================================
    # Count how many bubbles enter or finish during this one small time interval.
    admitted_count      = 0
    terminated_count    = 0

    # cursor is a time bookmark. It shows how far this function has already moved the system.
    cursor              = float(start_time)

    # Computer decimal numbers can differ by a tiny rounding amount even when they
    # should represent the same time. epsilon is a very small allowance that stops
    # such harmless rounding from creating an extra movement step or missing an event.
    epsilon             = 32.0 * np.finfo(float).eps * max(abs(end_time), 1.0)

    # ==================================================================================
    # ====== Read and handle all planned events that fall inside this interval.
    # ==================================================================================
    # Read planned entry events in permanent-ID order until the next event lies beyond this interval. 
    while state.next_event < schedule.count:

        # Look up the planned arrival time of the next bubble that has not yet been added to the inlet waiting list.
        event_time = float(schedule.planned_time_s[state.next_event])
        if event_time > end_time + epsilon:                             # Stop when the next event belongs to a future interval. 
            break
        event_time = max(event_time, cursor)                            # Never move the clock backwards.

        if event_time > cursor + epsilon:
            terminated_count    += _advance_segment(
                cursor,
                event_time,
                state,
                context,
                particle_cfg,
                dynamics_cfg,
                diagnostics,
                evaluate_rhs,
                accepted_target_path_batches=accepted_target_path_batches,
            )

            # Existing bubbles may have moved away from the inlet during this segment, 
            # so give older waiting bubbles another chance to enter.
            admitted_count      += _admit_waiting(event_time, state, schedule, context, diagnostics)
            cursor              = event_time            # Move the time bookmark to the event that has just been reached.

        # Give the newly reached schedule event its permanent bubble ID.
        new_event_id        = state.next_event

        # Put every new arrival into the same waiting list first, even when the
        # inlet is empty. This makes new and older arrivals follow one safety rule.
        state.waiting_ids.append(new_event_id)

        # Point to the following schedule event so this event cannot be added twice.
        state.next_event    += 1

        # Try to admit all currently waiting bubbles without placing any bubble
        # where it would overlap another finite-size bubble at the inlet.
        admitted_count      += _admit_waiting(cursor, state, schedule, context, diagnostics)

        # If the new ID is still in the list, its first entry attempt was blocked.
        # Record that fact for diagnostics; the bubble remains available for later retries.
        if not state.active[new_event_id]:
            diagnostics.inlet_wait_events += 1

    # After all events inside this interval are handled, 
    # move active bubbles through any time that remains between the last event and end_time.
    if end_time > cursor + epsilon:
        terminated_count    += _advance_segment(
            cursor,
            end_time,
            state,
            context,
            particle_cfg,
            dynamics_cfg,
            diagnostics,
            evaluate_rhs,
            accepted_target_path_batches=accepted_target_path_batches,
        )

        # Movement may have cleared the inlet again, 
        # so make one final entry attempt at the interval end instead of making waiting bubbles wait longer.
        admitted_count      += _admit_waiting(end_time, state, schedule, context, diagnostics)

    return admitted_count, terminated_count


def _ensure_last_step_storage(state: _PerfusionState) -> None:
    """Create accepted-step observation arrays for legacy in-memory fixtures."""

    capacity = int(state.position_grid.shape[0])
    if state.last_generalized_velocity.shape != (capacity, 3):
        state.last_generalized_velocity = np.zeros(
            (capacity, 3), dtype=np.float64
        )
    for name, dtype in (
        ("last_contact_reaction_force_pn", np.float64),
        ("last_contact_active", bool),
        ("last_free_normal_velocity_um_s", np.float64),
        ("last_constrained_normal_velocity_um_s", np.float64),
        ("last_step_valid", bool),
    ):
        values = getattr(state, name)
        if values.shape != (capacity,):
            setattr(state, name, np.zeros(capacity, dtype=dtype))


def _ensure_target_exposure_storage(state: _PerfusionState) -> None:
    """Create permanent-ID exposure accumulators for legacy in-memory fixtures."""

    capacity = int(state.position_grid.shape[0])
    for name, dtype in (
        ("target_exposure_time_s", np.float64),
        ("target_reaction_area_time_um2_s", np.float64),
        ("target_exposure_event_count", np.int64),
        ("target_exposure_event_end_count", np.int64),
        ("target_exposure_open", bool),
        ("target_exposure_applicable_time_s", np.float64),
    ):
        value = np.asarray(getattr(state, name, np.empty(0)), dtype=dtype)
        if value.shape != (capacity,):
            setattr(state, name, np.zeros(capacity, dtype=dtype))


def _ensure_vessel_id_storage(state: _PerfusionState) -> None:
    """Create explicit root ownership for legacy single-vessel test fixtures."""

    capacity = int(state.position_grid.shape[0])
    if state.vessel_id.shape == (capacity,):
        return
    if state.vessel_id.size != 0:
        raise RuntimeError("Persistent topological vessel state has an invalid shape.")
    state.vessel_id = np.zeros(capacity, dtype=np.int32)
    state.vessel_id[np.asarray(state.active, dtype=bool)] = 1


def _active_ids_view(state: _PerfusionState) -> np.ndarray:
    """Return the stable active-ID prefix without rebuilding a Python list."""

    capacity = int(state.active.shape[0])
    buffer = np.asarray(state.active_id_buffer, dtype=np.int64)
    count = int(state.active_count)
    needs_rebuild = (
        buffer.shape != (capacity,)
        or count < 0
        or count > capacity
        or (count == 0 and bool(np.any(state.active)))
    )
    if needs_rebuild:
        active = np.flatnonzero(state.active).astype(np.int64)
        buffer = np.empty(capacity, dtype=np.int64)
        buffer[: active.size] = active
        state.active_id_buffer = buffer
        state.active_count = int(active.size)
        count = int(active.size)
    return buffer[:count]


def _append_active_id(state: _PerfusionState, bubble_id: int) -> None:
    """Append one admitted permanent ID to the preallocated active prefix."""

    _active_ids_view(state)
    count = int(state.active_count)
    if count >= state.active_id_buffer.size:
        raise RuntimeError("The active-particle ID buffer exceeded its capacity.")
    state.active_id_buffer[count] = int(bubble_id)
    state.active_count = count + 1


def _replace_active_ids(state: _PerfusionState, bubble_ids: np.ndarray) -> None:
    """Stable-compact the active prefix after one or more lifecycle events."""

    ids = np.asarray(bubble_ids, dtype=np.int64)
    _active_ids_view(state)
    if ids.size > state.active_id_buffer.size:
        raise RuntimeError("The active-particle ID buffer exceeded its capacity.")
    state.active_id_buffer[: ids.size] = ids
    state.active_count = int(ids.size)


def _apply_rbc_stochastic_batch_state(
    batch_state: _BatchAdvanceState,
    permanent_ids: np.ndarray,
    active_dt_s: np.ndarray,
    context: object,
    *,
    global_internal_step: int,
    physical_time_s: float,
    coefficient_position_grid: np.ndarray | None = None,
    coefficient_vessel_id: np.ndarray | None = None,
    coefficient_evaluation: object | None = None,
) -> _BatchAdvanceResult:
    """Apply one counter-based random geometry transaction after deterministic motion.

    The diffusion coefficient belongs to the beginning of the nominal physical
    substep.  Only the transverse direction is re-read after the deterministic
    transaction, from the continuously committed vessel owner.  This keeps the
    operator order independent of contact/time refinement and topology splits.
    """

    network = getattr(context, "red_blood_cell_network", None)
    if network is None:
        return _BatchAdvanceResult(batch_state, _Diagnostics())
    alive = np.asarray(batch_state.alive, dtype=bool)
    if not np.any(alive):
        return _BatchAdvanceResult(batch_state, _Diagnostics())

    ids = np.ascontiguousarray(permanent_ids, dtype=np.int64)
    positions = np.ascontiguousarray(
        batch_state.particles.position_grid, dtype=np.float64
    )
    vessel_id = np.ascontiguousarray(batch_state.vessel_id, dtype=np.int32)
    coefficient_positions = np.ascontiguousarray(
        positions if coefficient_position_grid is None else coefficient_position_grid,
        dtype=np.float64,
    )
    coefficient_owner = np.ascontiguousarray(
        vessel_id if coefficient_vessel_id is None else coefficient_vessel_id,
        dtype=np.int32,
    )
    if coefficient_positions.shape != positions.shape or coefficient_owner.shape != alive.shape:
        raise ValueError("RBC stochastic start-coefficient arrays have inconsistent shapes.")
    dt_by_lane = np.ascontiguousarray(active_dt_s, dtype=np.float64)
    if dt_by_lane.shape != alive.shape or np.any(dt_by_lane < 0.0):
        raise ValueError("RBC stochastic active durations must be non-negative.")
    radii = np.asarray(context.all_bubble_radii_um[ids], dtype=np.float64)

    cached_diffusivity = getattr(
        coefficient_evaluation,
        "red_blood_cell_transverse_diffusivity_um2_s",
        None,
    )
    cached_enabled = getattr(
        coefficient_evaluation, "red_blood_cell_diffusion_enabled", None
    )
    cached_applicability = getattr(
        coefficient_evaluation,
        "red_blood_cell_quantitative_applicability",
        None,
    )
    cached_transverse_valid = getattr(
        coefficient_evaluation,
        "red_blood_cell_transverse_space_valid",
        None,
    )
    use_cached_coefficients = all(
        value is not None
        for value in (
            cached_diffusivity,
            cached_enabled,
            cached_applicability,
            cached_transverse_valid,
        )
    )
    if use_cached_coefficients:
        transverse_diffusivity = np.asarray(
            cached_diffusivity, dtype=np.float64
        ).reshape(-1)
        diffusion_enabled = np.asarray(cached_enabled, dtype=bool).reshape(-1)
        quantitative_applicability = np.asarray(
            cached_applicability, dtype=bool
        ).reshape(-1)
        transverse_space_valid = np.asarray(
            cached_transverse_valid, dtype=bool
        ).reshape(-1)
        if not (
            transverse_diffusivity.shape
            == diffusion_enabled.shape
            == quantitative_applicability.shape
            == transverse_space_valid.shape
            == alive.shape
        ):
            raise ValueError(
                "Cached RBC stochastic coefficient arrays have inconsistent shapes."
            )
    else:
        gradients = np.ascontiguousarray(
            sample_bilinear(context.velocity_gradient_s_inv, coefficient_positions),
            dtype=np.float64,
        )
        wall_gap = np.full(ids.size, np.nan, dtype=np.float64)
        wall_normal = np.zeros((ids.size, 2), dtype=np.float64)
        wall_unique = np.zeros(ids.size, dtype=bool)
        live_lanes = np.flatnonzero(alive)
        live_world = context.boundary_geometry.grid_to_world_xz(
            coefficient_positions[live_lanes]
        )
        exact = context.boundary_geometry.exact_solid_wall_state_xz_um_accelerated(
            live_world
        )
        wall_gap[live_lanes] = (
            np.asarray(exact.distance_um, dtype=np.float64).reshape(-1)
            - radii[live_lanes]
        )
        wall_normal[live_lanes] = np.asarray(
            exact.inward_normal_xz, dtype=np.float64
        ).reshape(-1, 2)
        wall_unique[live_lanes] = np.asarray(
            getattr(
                exact,
                "unique_nearest_wall",
                np.zeros(live_lanes.size, dtype=bool),
            ),
            dtype=bool,
        ).reshape(-1)

        evaluation = evaluate_red_blood_cell_transport(
            coefficient_owner,
            gradients,
            wall_gap,
            wall_normal,
            wall_unique,
            radii,
            alive,
            network,
            use_numba=bool(context.use_numba),
        )
        transverse_diffusivity = np.asarray(
            evaluation.transverse_diffusivity_um2_s, dtype=np.float64
        )
        diffusion_enabled = np.asarray(evaluation.diffusion_enabled, dtype=bool)
        quantitative_applicability = np.asarray(
            evaluation.quantitative_applicability, dtype=bool
        )
        transverse_space_valid = np.asarray(
            evaluation.transverse_space_valid, dtype=bool
        )
    displacement = np.zeros((ids.size, 2), dtype=np.float64)
    stochastic_lanes = np.flatnonzero(
        alive
        & (dt_by_lane > 0.0)
        & diffusion_enabled
    )
    if stochastic_lanes.size:
        normal = counter_normal_batch(
            int(context.random_seed),
            ids[stochastic_lanes],
            int(global_internal_step),
            0,
            use_numba=bool(context.use_numba),
        )
        amplitude = np.sqrt(
            2.0
            * transverse_diffusivity[stochastic_lanes]
            * dt_by_lane[stochastic_lanes]
        )
        final_tangent = np.asarray(
            network.dense_downstream_tangent_xz_by_vessel_id[
                vessel_id[stochastic_lanes].astype(np.int64) - 1
            ],
            dtype=np.float64,
        )
        final_transverse = np.column_stack(
            (-final_tangent[:, 1], final_tangent[:, 0])
        )
        displacement[stochastic_lanes] = (
            amplitude[:, None]
            * normal[:, None]
            * final_transverse
        )

    local = _Diagnostics()
    local.red_blood_cell_diffusion_enabled_observations = int(
        np.count_nonzero(alive & diffusion_enabled)
    )
    local.red_blood_cell_stochastic_observations = int(np.count_nonzero(alive))
    local.red_blood_cell_invalid_transverse_space_observations = int(
        np.count_nonzero(alive & ~transverse_space_valid)
    )
    if stochastic_lanes.size == 0:
        return _BatchAdvanceResult(
            batch_state,
            local,
            stochastic_sweep=None,
            red_blood_cell_quantitative_applicability=np.asarray(
                quantitative_applicability, dtype=bool
            ).copy(),
        )

    sweep = sweep_rbc_stochastic_displacement(
        positions,
        displacement,
        ids,
        alive,
        vessel_id,
        radii,
        batch_state.termination_reason,
        batch_state.exit_time_s,
        physical_time_s=float(physical_time_s),
        boundary_geometry=context.boundary_geometry,
        topological_ownership=context.topological_ownership,
        use_numba=bool(context.use_numba),
    )
    accepted = _accepted_batch_state(
        batch_state,
        sweep.position_grid,
        np.asarray(batch_state.particles.rotation_angle_rad, dtype=np.float64),
        np.asarray(batch_state.particles.bond_count_expected, dtype=np.float64),
        np.asarray(
            batch_state.particles.bond_total_tangential_extension_um,
            dtype=np.float64,
        ),
        alive=sweep.alive,
        vessel_id=sweep.vessel_id,
        termination_reason=sweep.termination_reason,
        exit_time_s=sweep.exit_time_s,
    )
    requested = sweep.requested_displacement_um[stochastic_lanes]
    local.red_blood_cell_random_wall_reflection_count = int(
        np.sum(sweep.reflection_count)
    )
    local.maximum_red_blood_cell_random_displacement_um = (
        float(np.max(requested)) if requested.size else 0.0
    )
    local.red_blood_cell_random_displacement_squared_sum_um2 = float(
        np.sum(requested * requested)
    )
    local.red_blood_cell_random_displacement_observations = int(requested.size)
    local.outlet_event_count = int(np.count_nonzero(alive & ~sweep.alive))
    if sweep.topology_event_bubble_id.size:
        local.topological_transition_count = int(
            sweep.topology_event_bubble_id.size
        )
        local.topological_event_bubble_id.extend(
            sweep.topology_event_bubble_id.astype(int).tolist()
        )
        local.topological_event_time_s.extend(
            [float(physical_time_s)] * sweep.topology_event_bubble_id.size
        )
        local.topological_event_from_vessel_id.extend(
            sweep.topology_event_from_vessel_id.astype(int).tolist()
        )
        local.topological_event_to_vessel_id.extend(
            sweep.topology_event_to_vessel_id.astype(int).tolist()
        )
        local.topological_event_section_index.extend(
            sweep.topology_event_section_index.astype(int).tolist()
        )
        local.topological_event_position_xz_um.extend(
            (float(row[0]), float(row[1]))
            for row in sweep.topology_event_position_xz_um
        )
    return _BatchAdvanceResult(
        accepted,
        local,
        stochastic_sweep=sweep,
        red_blood_cell_quantitative_applicability=np.asarray(
            quantitative_applicability, dtype=bool
        ).copy(),
    )


def _apply_rbc_stochastic_global_state(
    state: _PerfusionState,
    context: object,
    diagnostics: _Diagnostics,
    *,
    step_start_time_s: float,
    step_end_time_s: float,
    global_internal_step: int,
    step_start_ids: np.ndarray,
    step_start_position_grid: np.ndarray,
    step_start_vessel_id: np.ndarray,
    schedule: PerfusionSchedule,
) -> _GlobalRbcStochasticResult:
    """Gather, apply, and commit the stochastic half of one nominal step."""

    if getattr(context, "red_blood_cell_network", None) is None:
        return _GlobalRbcStochasticResult(
            0, np.empty(0, dtype=np.int64), None, None
        )
    ids = _active_ids_view(state).copy()
    if ids.size == 0:
        return _GlobalRbcStochasticResult(0, ids, None, None)
    _ensure_last_step_storage(state)
    _ensure_vessel_id_storage(state)
    batch = _BatchAdvanceState(
        particles=BatchLocalState(
            position_grid=np.asarray(state.position_grid[ids], dtype=np.float64),
            rotation_angle_rad=np.asarray(state.rotation_angle_rad[ids], dtype=np.float64),
            bond_count_expected=np.asarray(state.bond_count_expected[ids], dtype=np.float64),
            bond_total_tangential_extension_um=np.asarray(
                state.bond_total_tangential_extension_um[ids], dtype=np.float64
            ),
        ),
        alive=np.ones(ids.size, dtype=bool),
        vessel_id=np.asarray(state.vessel_id[ids], dtype=np.int32).copy(),
        termination_reason=np.asarray(state.termination_reason[ids], dtype=np.uint8).copy(),
        exit_time_s=np.asarray(state.exit_time_s[ids], dtype=np.float64).copy(),
        last_generalized_velocity=np.asarray(state.last_generalized_velocity[ids], dtype=np.float64).copy(),
        last_contact_reaction_force_pn=np.asarray(state.last_contact_reaction_force_pn[ids], dtype=np.float64).copy(),
        last_contact_active=np.asarray(state.last_contact_active[ids], dtype=bool).copy(),
        last_free_normal_velocity_um_s=np.asarray(state.last_free_normal_velocity_um_s[ids], dtype=np.float64).copy(),
        last_constrained_normal_velocity_um_s=np.asarray(state.last_constrained_normal_velocity_um_s[ids], dtype=np.float64).copy(),
        last_step_valid=np.asarray(state.last_step_valid[ids], dtype=bool).copy(),
    )
    admission = np.asarray(state.admission_time_s[ids], dtype=np.float64)
    active_dt = np.maximum(
        float(step_end_time_s) - np.maximum(admission, float(step_start_time_s)),
        0.0,
    )
    start_ids = np.asarray(step_start_ids, dtype=np.int64)
    start_positions = np.asarray(step_start_position_grid, dtype=np.float64)
    start_owner = np.asarray(step_start_vessel_id, dtype=np.int32)
    if start_positions.shape != (start_ids.size, 2) or start_owner.shape != start_ids.shape:
        raise ValueError("RBC nominal-step start snapshots have inconsistent shapes.")
    coefficient_positions = np.empty((ids.size, 2), dtype=np.float64)
    coefficient_owner = np.empty(ids.size, dtype=np.int32)
    if start_ids.size:
        order = np.argsort(start_ids, kind="stable")
        sorted_ids = start_ids[order]
        locations = np.searchsorted(sorted_ids, ids)
        matched = locations < sorted_ids.size
        matched[matched] &= sorted_ids[locations[matched]] == ids[matched]
        if np.any(matched):
            source = order[locations[matched]]
            coefficient_positions[matched] = start_positions[source]
            coefficient_owner[matched] = start_owner[source]
    else:
        matched = np.zeros(ids.size, dtype=bool)
    newly_admitted = ~matched
    if np.any(newly_admitted):
        coefficient_positions[newly_admitted] = np.asarray(
            schedule.position_grid[ids[newly_admitted]], dtype=np.float64
        )
        coefficient_owner[newly_admitted] = np.asarray(
            schedule.initial_vessel_id[ids[newly_admitted]], dtype=np.int32
        )
    result = _apply_rbc_stochastic_batch_state(
        batch,
        ids,
        active_dt,
        context,
        global_internal_step=int(global_internal_step),
        physical_time_s=float(step_end_time_s),
        coefficient_position_grid=coefficient_positions,
        coefficient_vessel_id=coefficient_owner,
    )
    state.position_grid[ids] = result.state.particles.position_grid
    state.active[ids] = result.state.alive
    state.vessel_id[ids] = result.state.vessel_id
    state.termination_reason[ids] = result.state.termination_reason
    state.exit_time_s[ids] = result.state.exit_time_s
    if np.any(~result.state.alive):
        _replace_active_ids(state, ids[result.state.alive])
    _merge_diagnostics(diagnostics, result.diagnostics)
    return _GlobalRbcStochasticResult(
        int(np.count_nonzero(~result.state.alive)),
        ids,
        result.stochastic_sweep,
        result.red_blood_cell_quantitative_applicability,
    )


def _rbc_quantitative_applicability_at_step_start(
    context: object,
    positions_grid: np.ndarray,
    vessel_id: np.ndarray,
) -> np.ndarray:
    """Return the existing H_t/shear applicability flag without raster ownership."""

    positions = np.asarray(positions_grid, dtype=np.float64).reshape(-1, 2)
    owner = np.asarray(vessel_id, dtype=np.int32).reshape(-1)
    if owner.shape != (positions.shape[0],):
        raise ValueError("RBC exposure-applicability arrays have inconsistent shapes.")
    network = getattr(context, "red_blood_cell_network", None)
    if network is None or positions.shape[0] == 0:
        return np.zeros(positions.shape[0], dtype=bool)
    row = owner.astype(np.int64) - 1
    dense_tube = np.asarray(
        network.dense_tube_hematocrit_by_vessel_id, dtype=np.float64
    )
    dense_diameter = np.asarray(
        network.dense_diameter_um_by_vessel_id, dtype=np.float64
    )
    if (
        np.any(row < 0)
        or np.any(row >= dense_tube.size)
        or np.any(row >= dense_diameter.size)
    ):
        raise ValueError(
            "An internal exposure path has no continuous topological vessel owner."
        )
    gradients = np.asarray(
        sample_bilinear(context.velocity_gradient_s_inv, positions),
        dtype=np.float64,
    )
    symmetric_01 = 0.5 * (gradients[:, 0, 1] + gradients[:, 1, 0])
    shear = np.sqrt(
        np.maximum(
            2.0
            * (
                gradients[:, 0, 0] ** 2
                + gradients[:, 1, 1] ** 2
                + 2.0 * symmetric_01**2
            ),
            0.0,
        )
    )
    tube = dense_tube[row]
    diameter = dense_diameter[row]
    return (
        np.isfinite(tube)
        & np.isfinite(shear)
        & np.isfinite(diameter)
        & (diameter >= CFL_DIAMETER_MIN_UM)
        & (diameter <= CFL_DIAMETER_MAX_UM)
        & (tube >= TUBE_HEMATOCRIT_QUANTITATIVE_MIN)
        & (tube <= TUBE_HEMATOCRIT_QUANTITATIVE_MAX)
        & (shear >= SHEAR_RATE_QUANTITATIVE_MIN_S_INV)
        & (shear <= SHEAR_RATE_QUANTITATIVE_MAX_S_INV)
    )


def _accumulate_internal_target_exposure(
    state: _PerfusionState,
    context: object,
    *,
    permanent_ids: np.ndarray,
    start_position_grid: np.ndarray,
    start_vessel_id: np.ndarray,
    deterministic_end_position_grid: np.ndarray,
    final_position_grid: np.ndarray,
    active_duration_s: np.ndarray,
    terminated_at_end: np.ndarray,
    quantitatively_applicable: np.ndarray | None = None,
    stochastic_sweep: RbcStochasticSweepResult | None = None,
    stochastic_permanent_ids: np.ndarray | None = None,
    deterministic_path_batches: list[tuple[
        np.ndarray, tuple[np.ndarray, ...], tuple[np.ndarray, ...]
    ]] | None = None,
) -> None:
    """Accumulate one nominal step without adding saved-frame records.

    The CSR path preserves every accepted deterministic transaction endpoint,
    followed by every stochastic topology/wall/outlet event and the final
    endpoint.  All pieces share the single physical duration supplied for that
    permanent ID; the random geometry sweep consumes no second time interval.
    """

    target = getattr(context, "molecular_target_field", None)
    capture = float(getattr(context, "molecular_capture_distance_um", 0.0))
    if target is None or not bool(getattr(target, "enabled", False)) or capture <= 0.0:
        return
    ids = np.asarray(permanent_ids, dtype=np.int64).reshape(-1)
    if ids.size == 0:
        return
    start_grid = np.asarray(start_position_grid, dtype=np.float64).reshape(-1, 2)
    start_owner = np.asarray(start_vessel_id, dtype=np.int32).reshape(-1)
    deterministic_grid = np.asarray(
        deterministic_end_position_grid, dtype=np.float64
    ).reshape(-1, 2)
    final_grid = np.asarray(final_position_grid, dtype=np.float64).reshape(-1, 2)
    duration = np.asarray(active_duration_s, dtype=np.float64).reshape(-1)
    terminated = np.asarray(terminated_at_end, dtype=bool).reshape(-1)
    lane_count = int(ids.size)
    if not (
        start_grid.shape == deterministic_grid.shape == final_grid.shape == (lane_count, 2)
        and start_owner.shape == duration.shape == terminated.shape == (lane_count,)
    ):
        raise ValueError("Internal target-exposure step arrays have inconsistent shapes.")
    supplied_applicability = None
    if quantitatively_applicable is not None:
        supplied_applicability = np.asarray(
            quantitatively_applicable, dtype=bool
        ).reshape(-1)
        if supplied_applicability.shape != (lane_count,):
            raise ValueError(
                "Internal target-exposure applicability must match the particle lanes."
            )
    positive_time = duration > 0.0
    if not np.any(positive_time):
        return
    if not np.all(positive_time):
        selected = np.flatnonzero(positive_time)
        ids = ids[selected]
        start_grid = start_grid[selected]
        start_owner = start_owner[selected]
        deterministic_grid = deterministic_grid[selected]
        final_grid = final_grid[selected]
        duration = duration[selected]
        terminated = terminated[selected]
        if supplied_applicability is not None:
            supplied_applicability = supplied_applicability[selected]
        lane_count = int(ids.size)

    order = np.argsort(ids, kind="stable")
    sorted_ids = ids[order]

    deterministic_count = np.zeros(lane_count, dtype=np.int64)
    deterministic_step_lane_parts: list[np.ndarray] = []
    deterministic_ordinal_parts: list[np.ndarray] = []
    deterministic_position_grid_parts: list[np.ndarray] = []
    if deterministic_path_batches is None:
        deterministic_step_lane_parts.append(
            np.arange(lane_count, dtype=np.int64)
        )
        deterministic_ordinal_parts.append(
            np.zeros(lane_count, dtype=np.int64)
        )
        deterministic_position_grid_parts.append(deterministic_grid)
        deterministic_count[:] = 1
    else:
        next_ordinal = np.zeros(lane_count, dtype=np.int64)
        for batch_ids_raw, path_positions, path_active in deterministic_path_batches:
            batch_ids = np.asarray(batch_ids_raw, dtype=np.int64).reshape(-1)
            if len(path_positions) != len(path_active):
                raise ValueError(
                    "Accepted deterministic path positions and masks do not align."
                )
            for position_raw, active_raw in zip(
                path_positions, path_active, strict=True
            ):
                position = np.asarray(position_raw, dtype=np.float64)
                active_mask = np.asarray(active_raw, dtype=bool).reshape(-1)
                if (
                    position.shape != (batch_ids.size, 2)
                    or active_mask.shape != (batch_ids.size,)
                ):
                    raise ValueError(
                        "Accepted deterministic path batch arrays are inconsistent."
                    )
                active_batch_lane = np.flatnonzero(active_mask)
                if active_batch_lane.size == 0:
                    continue
                active_ids = batch_ids[active_batch_lane]
                locations = np.searchsorted(sorted_ids, active_ids)
                found = locations < sorted_ids.size
                found[found] &= sorted_ids[locations[found]] == active_ids[found]
                if not np.any(found):
                    continue
                mapped = order[locations[found]]
                ordinal = next_ordinal[mapped].copy()
                next_ordinal[mapped] += 1
                deterministic_step_lane_parts.append(mapped)
                deterministic_ordinal_parts.append(ordinal)
                deterministic_position_grid_parts.append(
                    position[active_batch_lane[found]]
                )
        deterministic_count[:] = next_ordinal

    deterministic_step_lane = (
        np.concatenate(deterministic_step_lane_parts)
        if deterministic_step_lane_parts
        else np.empty(0, dtype=np.int64)
    )
    deterministic_ordinal = (
        np.concatenate(deterministic_ordinal_parts)
        if deterministic_ordinal_parts
        else np.empty(0, dtype=np.int64)
    )
    deterministic_position_grid = (
        np.vstack(deterministic_position_grid_parts)
        if deterministic_position_grid_parts
        else np.empty((0, 2), dtype=np.float64)
    )

    stochastic_event_count = np.zeros(lane_count, dtype=np.int64)
    event_step_lane = np.empty(0, dtype=np.int64)
    event_ordinal = np.empty(0, dtype=np.int64)
    event_position = np.empty((0, 2), dtype=np.float64)
    if stochastic_sweep is not None:
        sweep_ids = np.asarray(stochastic_permanent_ids, dtype=np.int64).reshape(-1)
        if sweep_ids.shape != np.asarray(stochastic_sweep.event_count).shape:
            raise ValueError("Stochastic exposure IDs do not match the sweep lanes.")
        locations = np.searchsorted(sorted_ids, sweep_ids)
        found = locations < sorted_ids.size
        found[found] &= sorted_ids[locations[found]] == sweep_ids[found]
        sweep_to_step = np.full(sweep_ids.size, -1, dtype=np.int64)
        sweep_to_step[found] = order[locations[found]]
        if np.any(found):
            stochastic_event_count[sweep_to_step[found]] = np.asarray(
                stochastic_sweep.event_count, dtype=np.int64
            )[found]
        event_mask = np.asarray(stochastic_sweep.event_code, dtype=np.uint8) != 0
        event_ordinal_matrix = np.cumsum(event_mask, axis=0, dtype=np.int64) - 1
        event_round, event_sweep_lane = np.nonzero(event_mask)
        mapped = sweep_to_step[event_sweep_lane]
        keep = mapped >= 0
        event_step_lane = mapped[keep]
        event_ordinal = event_ordinal_matrix[
            event_round[keep], event_sweep_lane[keep]
        ]
        event_position = np.asarray(
            stochastic_sweep.event_position_xz_um, dtype=np.float64
        )[event_round[keep], event_sweep_lane[keep]]

    point_count = deterministic_count + stochastic_event_count + 2
    offsets = np.concatenate(
        (np.zeros(1, dtype=np.int64), np.cumsum(point_count, dtype=np.int64))
    )
    points = np.empty((int(offsets[-1]), 2), dtype=np.float64)
    starts = offsets[:-1]
    points[starts] = context.boundary_geometry.grid_to_world_xz(start_grid)
    if deterministic_step_lane.size:
        points[
            starts[deterministic_step_lane] + 1 + deterministic_ordinal
        ] = context.boundary_geometry.grid_to_world_xz(
            deterministic_position_grid
        )
    if event_step_lane.size:
        points[
            starts[event_step_lane]
            + 1
            + deterministic_count[event_step_lane]
            + event_ordinal
        ] = event_position
    points[offsets[1:] - 1] = context.boundary_geometry.grid_to_world_xz(final_grid)

    _ensure_target_exposure_storage(state)
    radii = np.asarray(context.all_bubble_radii_um[ids], dtype=np.float64)
    candidate_query = getattr(target, "polyline_target_candidate_mask", None)
    if callable(candidate_query):
        candidate = np.asarray(
            candidate_query(
                points,
                offsets,
                radii,
                capture,
            ),
            dtype=bool,
        ).reshape(-1)
        if candidate.shape != (lane_count,):
            raise ValueError(
                "The molecular target broad phase returned an incompatible lane mask."
            )
        # A previously open event must still be evaluated so a departure can
        # close it.  All other rejected lanes are rigorously too far from every
        # target-positive wall element along their already accepted path.
        candidate |= np.asarray(state.target_exposure_open[ids], dtype=bool)
        if not np.all(candidate):
            selected = np.flatnonzero(candidate)
            if selected.size == 0:
                return
            point_counts = np.diff(offsets)
            points = points[np.repeat(candidate, point_counts)]
            offsets = np.concatenate(
                (
                    np.zeros(1, dtype=np.int64),
                    np.cumsum(point_counts[selected], dtype=np.int64),
                )
            )
            ids = ids[selected]
            start_grid = start_grid[selected]
            start_owner = start_owner[selected]
            duration = duration[selected]
            terminated = terminated[selected]
            radii = radii[selected]
            if supplied_applicability is not None:
                supplied_applicability = supplied_applicability[selected]
            lane_count = int(ids.size)
    applicable = (
        _rbc_quantitative_applicability_at_step_start(
            context, start_grid, start_owner
        )
        if supplied_applicability is None
        else supplied_applicability
    )
    exposure = integrate_internal_target_exposure(
        permanent_ids=ids,
        path_points_xz_um=points,
        path_point_offsets=offsets,
        dt_s=duration,
        radius_um=radii,
        capture_distance_um=capture,
        boundary_geometry=context.boundary_geometry,
        molecular_target_field=target,
        initially_exposed=state.target_exposure_open[ids],
        quantitatively_applicable=applicable,
        terminated_at_end=terminated,
    )
    state.target_exposure_time_s[ids] += exposure.exposure_time_s
    state.target_reaction_area_time_um2_s[ids] += (
        exposure.reaction_area_time_um2_s
    )
    state.target_exposure_event_count[ids] += exposure.event_start_count
    state.target_exposure_event_end_count[ids] += exposure.event_end_count
    state.target_exposure_open[ids] = exposure.exposure_open_at_end
    state.target_exposure_applicable_time_s[ids] += (
        exposure.quantitatively_applicable_exposure_time_s
    )


def _advance_segment(
    start_time: float,
    end_time: float,
    state: _PerfusionState,
    context: object,
    particle_cfg: ParticleConfig,
    dynamics_cfg: ParticleDynamicsConfig,
    diagnostics: _Diagnostics,
    evaluate_rhs: object,
    *,
    accepted_target_path_batches: list[tuple[
        np.ndarray, tuple[np.ndarray, ...], tuple[np.ndarray, ...]
    ]] | None = None,
) -> int:
    """Advance one injection-free interval and commit it as one transaction.

    Every speculative Euler/Heun state lives in detached arrays.  The permanent
    simulator state is updated only after the complete physical interval has
    been accepted, so a rejected position-level trial cannot leak position,
    rotation, bond, lifecycle, or ordinary RHS diagnostics into the run.
    """

    ids = _active_ids_view(state)
    if ids.size == 0:
        return 0

    _ensure_last_step_storage(state)
    _ensure_vessel_id_storage(state)

    initial_alive = np.ones(ids.size, dtype=bool)
    local_state = _BatchAdvanceState(
        particles=BatchLocalState(
            position_grid=np.asarray(state.position_grid[ids], dtype=np.float64),
            rotation_angle_rad=np.asarray(
                state.rotation_angle_rad[ids], dtype=np.float64
            ),
            bond_count_expected=np.asarray(
                state.bond_count_expected[ids], dtype=np.float64
            ),
            bond_total_tangential_extension_um=np.asarray(
                state.bond_total_tangential_extension_um[ids], dtype=np.float64
            ),
        ),
        alive=initial_alive,
        vessel_id=np.asarray(state.vessel_id[ids], dtype=np.int32).copy(),
        termination_reason=np.asarray(
            state.termination_reason[ids], dtype=np.uint8
        ).copy(),
        exit_time_s=np.asarray(state.exit_time_s[ids], dtype=np.float64).copy(),
        last_generalized_velocity=np.asarray(
            state.last_generalized_velocity[ids], dtype=np.float64
        ).copy(),
        last_contact_reaction_force_pn=np.asarray(
            state.last_contact_reaction_force_pn[ids], dtype=np.float64
        ).copy(),
        last_contact_active=np.asarray(
            state.last_contact_active[ids], dtype=bool
        ).copy(),
        last_free_normal_velocity_um_s=np.asarray(
            state.last_free_normal_velocity_um_s[ids], dtype=np.float64
        ).copy(),
        last_constrained_normal_velocity_um_s=np.asarray(
            state.last_constrained_normal_velocity_um_s[ids], dtype=np.float64
        ).copy(),
        last_step_valid=np.asarray(
            state.last_step_valid[ids], dtype=bool
        ).copy(),
    )
    interval = PhysicalTimeInterval(float(start_time), float(end_time), 0)
    result = _advance_batch_interval(
        interval,
        local_state,
        ids,
        context,
        particle_cfg,
        dynamics_cfg,
        evaluate_rhs,
    )
    if (
        accepted_target_path_batches is not None
        and result.accepted_path_position_grid
    ):
        accepted_target_path_batches.append(
            (
                np.asarray(ids, dtype=np.int64).copy(),
                result.accepted_path_position_grid,
                result.accepted_path_active,
            )
        )

    # Commit all coupled particle state only after both refined halves (if any)
    # have covered the complete requested physical interval successfully.
    state.position_grid[ids] = result.state.particles.position_grid
    state.rotation_angle_rad[ids] = result.state.particles.rotation_angle_rad
    if result.state.particles.bond_count_expected is not None:
        state.bond_count_expected[ids] = (
            result.state.particles.bond_count_expected
        )
    if result.state.particles.bond_total_tangential_extension_um is not None:
        state.bond_total_tangential_extension_um[ids] = (
            result.state.particles.bond_total_tangential_extension_um
        )
    state.active[ids] = result.state.alive
    state.vessel_id[ids] = result.state.vessel_id
    state.termination_reason[ids] = result.state.termination_reason
    state.exit_time_s[ids] = result.state.exit_time_s
    if result.state.last_generalized_velocity is not None:
        state.last_generalized_velocity[ids] = result.state.last_generalized_velocity
    if result.state.last_contact_reaction_force_pn is not None:
        state.last_contact_reaction_force_pn[ids] = (
            result.state.last_contact_reaction_force_pn
        )
    if result.state.last_contact_active is not None:
        state.last_contact_active[ids] = result.state.last_contact_active
    if result.state.last_free_normal_velocity_um_s is not None:
        state.last_free_normal_velocity_um_s[ids] = (
            result.state.last_free_normal_velocity_um_s
        )
    if result.state.last_constrained_normal_velocity_um_s is not None:
        state.last_constrained_normal_velocity_um_s[ids] = (
            result.state.last_constrained_normal_velocity_um_s
        )
    if result.state.last_step_valid is not None:
        state.last_step_valid[ids] = result.state.last_step_valid
    if np.any(~result.state.alive):
        _replace_active_ids(state, ids[result.state.alive])
    _merge_diagnostics(diagnostics, result.diagnostics)
    return int(np.count_nonzero(initial_alive & ~result.state.alive))


def _can_batch_injection_free_frame(
    frame_end_time_s: float,
    state: _PerfusionState,
    schedule: PerfusionSchedule,
) -> bool:
    """Return whether a complete saved-frame interval has no admission work."""

    if state.waiting_ids:
        return False
    if state.next_event >= schedule.count:
        return True
    epsilon = 32.0 * np.finfo(float).eps * max(abs(frame_end_time_s), 1.0)
    return bool(
        float(schedule.planned_time_s[state.next_event])
        > float(frame_end_time_s) + epsilon
    )


def _advance_injection_free_substeps(
    start_time_s: float,
    substep_count: int,
    internal_dt_s: float,
    state: _PerfusionState,
    context: object,
    particle_cfg: ParticleConfig,
    dynamics_cfg: ParticleDynamicsConfig,
    diagnostics: _Diagnostics,
    evaluate_rhs: object,
    *,
    first_global_internal_step: int,
) -> tuple[int, float]:
    """Advance several event-free internal intervals with one gather and commit."""

    count = int(substep_count)
    if count < 1:
        raise ValueError("An injection-free frame batch must contain at least one step.")
    dt = float(internal_dt_s)
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("An injection-free frame batch requires a positive finite dt.")
    cursor = float(start_time_s)
    ids = _active_ids_view(state)
    if ids.size == 0:
        for _ in range(count):
            cursor += dt
        return 0, cursor

    _ensure_last_step_storage(state)
    _ensure_vessel_id_storage(state)
    initial_alive = np.ones(ids.size, dtype=bool)
    local_state = _BatchAdvanceState(
        particles=BatchLocalState(
            position_grid=np.asarray(state.position_grid[ids], dtype=np.float64),
            rotation_angle_rad=np.asarray(
                state.rotation_angle_rad[ids], dtype=np.float64
            ),
            bond_count_expected=np.asarray(
                state.bond_count_expected[ids], dtype=np.float64
            ),
            bond_total_tangential_extension_um=np.asarray(
                state.bond_total_tangential_extension_um[ids], dtype=np.float64
            ),
        ),
        alive=initial_alive,
        vessel_id=np.asarray(state.vessel_id[ids], dtype=np.int32).copy(),
        termination_reason=np.asarray(
            state.termination_reason[ids], dtype=np.uint8
        ).copy(),
        exit_time_s=np.asarray(state.exit_time_s[ids], dtype=np.float64).copy(),
        last_generalized_velocity=np.asarray(
            state.last_generalized_velocity[ids], dtype=np.float64
        ).copy(),
        last_contact_reaction_force_pn=np.asarray(
            state.last_contact_reaction_force_pn[ids], dtype=np.float64
        ).copy(),
        last_contact_active=np.asarray(
            state.last_contact_active[ids], dtype=bool
        ).copy(),
        last_free_normal_velocity_um_s=np.asarray(
            state.last_free_normal_velocity_um_s[ids], dtype=np.float64
        ).copy(),
        last_constrained_normal_velocity_um_s=np.asarray(
            state.last_constrained_normal_velocity_um_s[ids], dtype=np.float64
        ).copy(),
        last_step_valid=np.asarray(
            state.last_step_valid[ids], dtype=bool
        ).copy(),
    )
    batch_diagnostics = _Diagnostics()
    target_exposure_enabled = (
        context.molecular_target_field is not None
        and bool(context.molecular_target_field.enabled)
        and float(context.molecular_capture_distance_um) > 0.0
    )
    for substep_index in range(count):
        next_time = cursor + dt
        if np.any(local_state.alive):
            substep_start_alive = np.asarray(local_state.alive, dtype=bool).copy()
            rbc_start_positions = None
            rbc_start_owner = None
            start_evaluation_cache: list[object] = []
            if context.red_blood_cell_network is not None or target_exposure_enabled:
                rbc_start_positions = np.asarray(
                    local_state.particles.position_grid, dtype=np.float64
                ).copy()
                rbc_start_owner = np.asarray(
                    local_state.vessel_id, dtype=np.int32
                ).copy()
            result = _advance_batch_interval(
                PhysicalTimeInterval(cursor, next_time, 0),
                local_state,
                ids,
                context,
                particle_cfg,
                dynamics_cfg,
                evaluate_rhs,
                start_evaluation_cache=start_evaluation_cache,
            )
            local_state = result.state
            _merge_diagnostics(batch_diagnostics, result.diagnostics)
            deterministic_end_positions = np.asarray(
                local_state.particles.position_grid, dtype=np.float64
            ).copy()
            stochastic = _apply_rbc_stochastic_batch_state(
                local_state,
                ids,
                np.where(local_state.alive, dt, 0.0),
                context,
                global_internal_step=(
                    int(first_global_internal_step) + substep_index
                ),
                physical_time_s=next_time,
                coefficient_position_grid=rbc_start_positions,
                coefficient_vessel_id=rbc_start_owner,
                coefficient_evaluation=(
                    start_evaluation_cache[0] if start_evaluation_cache else None
                ),
            )
            local_state = stochastic.state
            _merge_diagnostics(batch_diagnostics, stochastic.diagnostics)
            if target_exposure_enabled:
                lanes = np.flatnonzero(substep_start_alive)
                exit_time = np.asarray(
                    local_state.exit_time_s[lanes], dtype=np.float64
                )
                physical_end = np.where(
                    np.isfinite(exit_time),
                    np.minimum(exit_time, float(next_time)),
                    float(next_time),
                )
                _accumulate_internal_target_exposure(
                    state,
                    context,
                    permanent_ids=ids[lanes],
                    start_position_grid=np.asarray(
                        rbc_start_positions[lanes], dtype=np.float64
                    ),
                    start_vessel_id=np.asarray(
                        rbc_start_owner[lanes], dtype=np.int32
                    ),
                    deterministic_end_position_grid=np.asarray(
                        deterministic_end_positions[lanes], dtype=np.float64
                    ),
                    final_position_grid=np.asarray(
                        local_state.particles.position_grid[lanes], dtype=np.float64
                    ),
                    active_duration_s=np.maximum(
                        physical_end - float(cursor), 0.0
                    ),
                    terminated_at_end=~np.asarray(
                        local_state.alive[lanes], dtype=bool
                    ),
                    quantitatively_applicable=(
                        None
                        if stochastic.red_blood_cell_quantitative_applicability
                        is None
                        else np.asarray(
                            stochastic.red_blood_cell_quantitative_applicability,
                            dtype=bool,
                        )[lanes]
                    ),
                    stochastic_sweep=stochastic.stochastic_sweep,
                    stochastic_permanent_ids=ids,
                    deterministic_path_batches=[
                        (
                            np.asarray(ids, dtype=np.int64),
                            result.accepted_path_position_grid,
                            result.accepted_path_active,
                        )
                    ],
                )
        cursor = next_time

    state.position_grid[ids] = local_state.particles.position_grid
    state.rotation_angle_rad[ids] = local_state.particles.rotation_angle_rad
    if local_state.particles.bond_count_expected is not None:
        state.bond_count_expected[ids] = local_state.particles.bond_count_expected
    if local_state.particles.bond_total_tangential_extension_um is not None:
        state.bond_total_tangential_extension_um[ids] = (
            local_state.particles.bond_total_tangential_extension_um
        )
    state.active[ids] = local_state.alive
    state.vessel_id[ids] = local_state.vessel_id
    state.termination_reason[ids] = local_state.termination_reason
    state.exit_time_s[ids] = local_state.exit_time_s
    if local_state.last_generalized_velocity is not None:
        state.last_generalized_velocity[ids] = local_state.last_generalized_velocity
    if local_state.last_contact_reaction_force_pn is not None:
        state.last_contact_reaction_force_pn[ids] = (
            local_state.last_contact_reaction_force_pn
        )
    if local_state.last_contact_active is not None:
        state.last_contact_active[ids] = local_state.last_contact_active
    if local_state.last_free_normal_velocity_um_s is not None:
        state.last_free_normal_velocity_um_s[ids] = (
            local_state.last_free_normal_velocity_um_s
        )
    if local_state.last_constrained_normal_velocity_um_s is not None:
        state.last_constrained_normal_velocity_um_s[ids] = (
            local_state.last_constrained_normal_velocity_um_s
        )
    if local_state.last_step_valid is not None:
        state.last_step_valid[ids] = local_state.last_step_valid
    if np.any(~local_state.alive):
        _replace_active_ids(state, ids[local_state.alive])
    _merge_diagnostics(diagnostics, batch_diagnostics)
    return int(np.count_nonzero(initial_alive & ~local_state.alive)), cursor


def _advance_batch_interval(
    interval: PhysicalTimeInterval,
    start_state: _BatchAdvanceState,
    ids: np.ndarray,
    context: object,
    particle_cfg: ParticleConfig,
    dynamics_cfg: ParticleDynamicsConfig,
    evaluate_rhs: object,
    *,
    start_evaluation_cache: list[object] | None = None,
) -> _BatchAdvanceResult:
    """Accept one physical interval, recursively bisecting unsafe trials."""

    if not np.any(start_state.alive):
        return _BatchAdvanceResult(start_state, _Diagnostics())

    try:
        attempt = _attempt_batch_step(
            interval,
            start_state,
            ids,
            context,
            particle_cfg,
            dynamics_cfg,
            evaluate_rhs,
            start_evaluation_cache=start_evaluation_cache,
        )
    except _SplitAtBoundaryEvent as event_split:
        split_time = float(event_split.split_time_s)
        if not (float(interval.start_time_s) < split_time < float(interval.end_time_s)):
            raise RuntimeError(
                "A boundary event split did not consume a representable positive "
                "physical duration."
            ) from event_split
        first = PhysicalTimeInterval(
            float(interval.start_time_s),
            split_time,
            int(interval.refinement_depth),
        )
        second = PhysicalTimeInterval(
            split_time,
            float(interval.end_time_s),
            int(interval.refinement_depth),
        )
        first_result = _advance_batch_interval(
            first,
            start_state,
            ids,
            context,
            particle_cfg,
            dynamics_cfg,
            evaluate_rhs,
            start_evaluation_cache=start_evaluation_cache,
        )
        second_result = _advance_batch_interval(
            second,
            first_result.state,
            ids,
            context,
            particle_cfg,
            dynamics_cfg,
            evaluate_rhs,
            start_evaluation_cache=start_evaluation_cache,
        )
        merged = _Diagnostics()
        _merge_diagnostics(merged, first_result.diagnostics)
        _merge_diagnostics(merged, second_result.diagnostics)
        return _BatchAdvanceResult(
            second_result.state,
            merged,
            accepted_path_position_grid=(
                first_result.accepted_path_position_grid
                + second_result.accepted_path_position_grid
            ),
            accepted_path_active=(
                first_result.accepted_path_active
                + second_result.accepted_path_active
            ),
        )
    except _RefineContactTimeStep as error:
        try:
            first, second = split_physical_time_interval(
                interval,
                int(particle_cfg.contact_max_time_refinements),
            )
        except PhysicalTimeRefinementError as refinement_error:
            maximum_refinement_depth = int(particle_cfg.contact_max_time_refinements)
            refinement_stop_reason = str(getattr(refinement_error, "reason", "unknown"))
            rejected_duration_s = getattr(refinement_error, "duration_s", None)
            if rejected_duration_s is None:
                rejected_duration_s = float(interval.duration_s)
            local_time_ulp_s = getattr(refinement_error, "local_time_ulp_s", None)
            if local_time_ulp_s is None:
                local_time_ulp_s = max(
                    math.ulp(float(interval.start_time_s)),
                    math.ulp(float(interval.end_time_s)),
                )
            failed_lanes = np.asarray(error.failed_lanes, dtype=np.int64)
            failed_ids = (
                np.asarray(ids, dtype=np.int64)[failed_lanes]
                if failed_lanes.size
                else np.empty(0, dtype=np.int64)
            )
            records: list[ParticleWallContactFailure] = []
            for compact, lane_value in enumerate(failed_lanes):
                lane = int(lane_value)
                bubble_id = int(failed_ids[compact])
                position_grid = (
                    error.failure_positions_grid[compact]
                    if compact < error.failure_positions_grid.shape[0]
                    and np.all(np.isfinite(error.failure_positions_grid[compact]))
                    else start_state.particles.position_grid[lane]
                )
                event_fraction = (
                    float(error.failure_event_fractions[compact])
                    if compact < error.failure_event_fractions.size
                    and np.isfinite(error.failure_event_fractions[compact])
                    else 0.0
                )
                code = (
                    int(error.failure_codes[compact])
                    if compact < error.failure_codes.size
                    else -1
                )
                reason = {
                    2: "refinement_exhausted_invalid_contact_mobility",
                    3: "refinement_exhausted_endpoint_penetration",
                    4: "refinement_exhausted_swept_path_penetration",
                }.get(code, f"refinement_exhausted_failure_code_{code}")
                world = np.asarray(
                    context.boundary_geometry.grid_to_world_xz(position_grid),
                    dtype=np.float64,
                ).reshape(2)
                contact = inspect_predictive_wall_contact_failure(
                    local_lane=lane,
                    event_fraction=event_fraction,
                    position_xz_um=world,
                    bubble_radius_um=float(context.all_bubble_radii_um[bubble_id]),
                    boundary_geometry=context.boundary_geometry,
                    reason=reason,
                )
                start_world = np.asarray(
                    context.boundary_geometry.grid_to_world_xz(
                        start_state.particles.position_grid[lane]
                    ),
                    dtype=np.float64,
                ).reshape(2)
                records.append(
                    ParticleWallContactFailure(
                        permanent_microbubble_id=bubble_id,
                        physical_time_s=float(interval.start_time_s)
                        + event_fraction * float(interval.duration_s),
                        integration_stage=error.integration_stage,
                        contact=contact,
                        interval_start_position_xz_um=(
                            float(start_world[0]),
                            float(start_world[1]),
                        ),
                        rejected_trial_displacement_um=float(
                            np.linalg.norm(world - start_world)
                        ),
                    )
                )
            message = (
                "Wall-contact integration could not produce a strictly feasible "
                "state before physical-time refinement had to stop. "
                "No failed trial was committed and no particle was silently held. "
                f"Rejected interval=[{interval.start_time_s:.17g}, "
                f"{interval.end_time_s:.17g}] s, "
                f"refinement_stop_reason={refinement_stop_reason}, "
                f"duration_s={float(rejected_duration_s):.17g}, "
                f"local_time_ulp_s={float(local_time_ulp_s):.17g}, "
                f"depth={interval.refinement_depth}, "
                f"maximum_refinement_depth={maximum_refinement_depth}, "
                f"failure_codes={error.failure_codes.tolist()}; reason={error}."
            )
            raise ParticleWallContactGeometryError(
                message, tuple(records)
            ) from refinement_error

        first_result = _advance_batch_interval(
            first,
            start_state,
            ids,
            context,
            particle_cfg,
            dynamics_cfg,
            evaluate_rhs,
            start_evaluation_cache=start_evaluation_cache,
        )
        second_result = _advance_batch_interval(
            second,
            first_result.state,
            ids,
            context,
            particle_cfg,
            dynamics_cfg,
            evaluate_rhs,
            start_evaluation_cache=start_evaluation_cache,
        )
        merged = _Diagnostics()
        _merge_diagnostics(merged, first_result.diagnostics)
        _merge_diagnostics(merged, second_result.diagnostics)
        merged.contact_time_refinement_count += 1
        merged.maximum_contact_time_refinement_depth = max(
            merged.maximum_contact_time_refinement_depth,
            first.refinement_depth,
        )
        return _BatchAdvanceResult(
            second_result.state,
            merged,
            accepted_path_position_grid=(
                first_result.accepted_path_position_grid
                + second_result.accepted_path_position_grid
            ),
            accepted_path_active=(
                first_result.accepted_path_active
                + second_result.accepted_path_active
            ),
        )

    accepted_end = float(attempt.accepted_end_time_s)
    time_epsilon = _time_epsilon(interval.start_time_s, interval.end_time_s)
    if accepted_end >= interval.end_time_s - time_epsilon:
        path_position, path_active = _accepted_target_path_piece(
            start_state, attempt.state, context
        )
        return _BatchAdvanceResult(
            attempt.state,
            attempt.diagnostics,
            accepted_path_position_grid=path_position,
            accepted_path_active=path_active,
        )

    # A geometric event may finish only a prefix of this interval.  Re-evaluate
    # every still-alive lane at that event state for the untouched remainder.
    if accepted_end <= interval.start_time_s + time_epsilon:
        lifecycle_changed = not np.array_equal(
            attempt.state.alive,
            start_state.alive,
        )
        topology_changed = not np.array_equal(
            attempt.state.vessel_id,
            start_state.vessel_id,
        )
        if not lifecycle_changed and not topology_changed:
            raise RuntimeError(
                "A boundary event consumed zero physical time without changing "
                "lifecycle state or topological state. The step was not silently accepted."
            )
        accepted_end = float(interval.start_time_s)

    remaining = PhysicalTimeInterval(
        accepted_end,
        float(interval.end_time_s),
        int(interval.refinement_depth),
    )
    remainder = _advance_batch_interval(
        remaining,
        attempt.state,
        ids,
        context,
        particle_cfg,
        dynamics_cfg,
        evaluate_rhs,
        start_evaluation_cache=start_evaluation_cache,
    )
    merged = _Diagnostics()
    _merge_diagnostics(merged, attempt.diagnostics)
    _merge_diagnostics(merged, remainder.diagnostics)
    path_position, path_active = _accepted_target_path_piece(
        start_state, attempt.state, context
    )
    return _BatchAdvanceResult(
        remainder.state,
        merged,
        accepted_path_position_grid=(
            path_position + remainder.accepted_path_position_grid
        ),
        accepted_path_active=(path_active + remainder.accepted_path_active),
    )


def _accepted_target_path_piece(
    start_state: _BatchAdvanceState,
    accepted_state: _BatchAdvanceState,
    context: object,
) -> tuple[tuple[np.ndarray, ...], tuple[np.ndarray, ...]]:
    """Retain one accepted deterministic vertex only for internal exposure."""

    target = getattr(context, "molecular_target_field", None)
    if (
        target is None
        or not bool(getattr(target, "enabled", False))
        or float(getattr(context, "molecular_capture_distance_um", 0.0)) <= 0.0
    ):
        return (), ()
    position = np.asarray(
        accepted_state.particles.position_grid, dtype=np.float64
    ).copy()
    active = np.asarray(start_state.alive, dtype=bool).copy()
    if position.shape != (active.size, 2):
        raise ValueError("Accepted target-path state arrays are inconsistent.")
    return (position,), (active,)


def _attempt_batch_step(
    interval: PhysicalTimeInterval,
    start_state: _BatchAdvanceState,
    ids: np.ndarray,
    context: object,
    particle_cfg: ParticleConfig,
    dynamics_cfg: ParticleDynamicsConfig,
    evaluate_rhs: object,
    *,
    start_evaluation_cache: list[object] | None = None,
) -> _BatchAttemptResult:
    """Try one complete Revised-v15 predictive Euler/Heun interval.

    Production lifecycle queries come exclusively from
    ``context.boundary_geometry``.
    """

    alive = np.asarray(start_state.alive, dtype=bool)
    # The batch state is already detached from permanent simulator storage by
    # the advanced-index gather in ``_advance_segment``.  Trials below treat
    # these arrays as immutable and build new accepted arrays, so a second full
    # copy per internal step only adds memory traffic.
    positions = np.asarray(start_state.particles.position_grid, dtype=np.float64)
    angles = np.asarray(start_state.particles.rotation_angle_rad, dtype=np.float64)
    bond_count = np.asarray(
        start_state.particles.bond_count_expected, dtype=np.float64
    )
    bond_extension = np.asarray(
        start_state.particles.bond_total_tangential_extension_um,
        dtype=np.float64,
    )
    vessel_id = np.asarray(start_state.vessel_id, dtype=np.int32)
    if vessel_id.size == 0:
        vessel_id = np.ones(alive.shape, dtype=np.int32)
    if vessel_id.shape != alive.shape or np.any(alive & (vessel_id <= 0)):
        raise RuntimeError(
            "Every active particle must carry a valid Revised-v20 topological vessel ID."
        )
    dt = float(interval.duration_s)
    spacing = float(context.spacing_um)
    radii = np.asarray(context.all_bubble_radii_um[ids], dtype=np.float64)
    tolerance = float(particle_cfg.contact_geometry_tolerance_um)
    geometry = context.boundary_geometry

    start_evaluation = _call_particle_rhs(
        evaluate_rhs,
        positions,
        ids,
        alive,
        context,
        float(interval.start_time_s),
        bond_count,
        bond_extension,
        vessel_id,
    )
    if start_evaluation_cache is not None and not start_evaluation_cache:
        start_evaluation_cache.append(start_evaluation)
    local_diagnostics = _Diagnostics()
    _accumulate_diagnostics(start_evaluation, local_diagnostics)
    start_generalized = _generalized_velocity(start_evaluation)
    cached_start_distance = None
    cached_start_normal = None
    start_normal_valid = getattr(start_evaluation, "wall_normal_valid", None)
    if start_normal_valid is not None and np.all(
        ~alive | np.asarray(start_normal_valid, dtype=bool)
    ):
        cached_start_distance = (
            np.asarray(start_evaluation.wall_gap_um, dtype=np.float64) + radii
        )
        cached_start_normal = np.asarray(
            start_evaluation.wall_normal_xz, dtype=np.float64
        )
    predictor = _solve_v15_trial_with_failure_context(
        positions,
        alive,
        dt,
        start_generalized,
        start_evaluation.generalized_mobility,
        radii,
        context,
        tolerance,
        permanent_ids=ids,
        interval=interval,
        integration_stage="predictor",
        start_wall_distance_um=cached_start_distance,
        start_wall_normal_xz=cached_start_normal,
    )
    _raise_for_v15_refinement(predictor, integration_stage="predictor")
    predictor_topology = _inspect_topological_step(
        positions,
        predictor,
        vessel_id,
        alive,
        context,
    )
    _split_for_internal_events(interval, predictor, predictor_topology)

    immediate = _commit_immediate_outlets(
        interval,
        start_state,
        predictor,
        positions,
        angles,
        bond_count,
        bond_extension,
        alive,
    )
    if immediate is not None:
        local_diagnostics.outlet_event_count += int(
            np.count_nonzero(alive & ~immediate.state.alive)
        )
        return _BatchAttemptResult(
            immediate.state,
            local_diagnostics,
            immediate.accepted_end_time_s,
        )
    immediate_topology = _commit_immediate_topology(
        interval,
        start_state,
        predictor_topology,
        positions,
        angles,
        bond_count,
        bond_extension,
        alive,
        ids,
    )
    if immediate_topology is not None:
        _merge_diagnostics(local_diagnostics, immediate_topology.diagnostics)
        return _BatchAttemptResult(
            immediate_topology.state,
            local_diagnostics,
            immediate_topology.accepted_end_time_s,
        )

    predictor_positions = predictor.accepted_position_grid
    predictor_angles = angles + predictor.angle_increment_rad

    if str(dynamics_cfg.time_integrator) != "heun":
        new_count, new_extension, limited = _accept_euler_bonds(
            bond_count,
            bond_extension,
            start_evaluation,
            positions,
            predictor_positions,
            angles,
            predictor_angles,
            radii,
            alive,
            spacing,
            dt,
            bool(context.use_numba),
        )
        local_diagnostics.molecular_capacity_limited_accepted_step_count += int(
            limited
        )
        _accumulate_v15_position_diagnostics(
            positions,
            predictor_positions,
            start_evaluation,
            predictor,
            alive,
            dt,
            spacing,
            local_diagnostics,
            use_numba=bool(context.use_numba),
        )
        accepted_alive, termination, exit_time = _lifecycle_after_endpoint(
            interval,
            start_state,
            predictor,
        )
        accepted_vessel_id = _vessel_ids_after_endpoint(
            interval,
            start_state,
            predictor_topology,
            alive,
            ids,
            local_diagnostics,
        )
        local_diagnostics.outlet_event_count += int(
            np.count_nonzero(alive & ~accepted_alive)
        )
        return _BatchAttemptResult(
            _accepted_batch_state(
                start_state,
                predictor_positions,
                predictor_angles,
                new_count,
                new_extension,
                alive=accepted_alive,
                vessel_id=accepted_vessel_id,
                termination_reason=termination,
                exit_time_s=exit_time,
                **_step_observation_kwargs(
                    start_generalized,
                    predictor.contact_normal_xz,
                    predictor,
                    alive,
                ),
            ),
            local_diagnostics,
            float(interval.end_time_s),
        )

    # Heun's second RHS belongs to the position-level predictor.  Predictor
    # bond values are speculative and are discarded after building the final
    # corrected bond state below.
    predictor_count, predictor_extension, _ = _accept_euler_bonds(
        bond_count,
        bond_extension,
        start_evaluation,
        positions,
        predictor_positions,
        angles,
        predictor_angles,
        radii,
        alive,
        spacing,
        dt,
        bool(context.use_numba),
    )
    predictor_vessel_id = _vessel_ids_after_endpoint(
        interval,
        start_state,
        predictor_topology,
        alive,
        ids,
        None,
    )
    predicted_evaluation = _call_particle_rhs(
        evaluate_rhs,
        predictor_positions,
        ids,
        alive,
        context,
        float(interval.end_time_s),
        predictor_count,
        predictor_extension,
        predictor_vessel_id,
    )
    _accumulate_diagnostics(predicted_evaluation, local_diagnostics)
    predicted_generalized = _generalized_velocity(predicted_evaluation)
    corrected_generalized = 0.5 * (start_generalized + predicted_generalized)
    corrected_mobility = 0.5 * (
        np.asarray(start_evaluation.generalized_mobility, dtype=np.float64)
        + np.asarray(predicted_evaluation.generalized_mobility, dtype=np.float64)
    )
    corrected = _solve_v15_trial_with_failure_context(
        positions,
        alive,
        dt,
        corrected_generalized,
        corrected_mobility,
        radii,
        context,
        tolerance,
        permanent_ids=ids,
        interval=interval,
        integration_stage="corrector",
        start_wall_distance_um=cached_start_distance,
        start_wall_normal_xz=cached_start_normal,
    )
    _raise_for_v15_refinement(corrected, integration_stage="corrector")
    corrected_topology = _inspect_topological_step(
        positions,
        corrected,
        vessel_id,
        alive,
        context,
    )
    _split_for_internal_events(interval, corrected, corrected_topology)
    immediate = _commit_immediate_outlets(
        interval,
        start_state,
        corrected,
        positions,
        angles,
        bond_count,
        bond_extension,
        alive,
    )
    if immediate is not None:
        local_diagnostics.outlet_event_count += int(
            np.count_nonzero(alive & ~immediate.state.alive)
        )
        return _BatchAttemptResult(
            immediate.state,
            local_diagnostics,
            immediate.accepted_end_time_s,
        )
    immediate_topology = _commit_immediate_topology(
        interval,
        start_state,
        corrected_topology,
        positions,
        angles,
        bond_count,
        bond_extension,
        alive,
        ids,
    )
    if immediate_topology is not None:
        _merge_diagnostics(local_diagnostics, immediate_topology.diagnostics)
        return _BatchAttemptResult(
            immediate_topology.state,
            local_diagnostics,
            immediate_topology.accepted_end_time_s,
        )

    corrected_positions = corrected.accepted_position_grid
    corrected_angles = angles + corrected.angle_increment_rad
    new_count, new_extension, limited = _accept_heun_bonds(
        bond_count,
        bond_extension,
        predictor_count,
        start_evaluation,
        predicted_evaluation,
        positions,
        corrected_positions,
        angles,
        corrected_angles,
        radii,
        alive,
        spacing,
        dt,
        bool(context.use_numba),
    )
    local_diagnostics.molecular_capacity_limited_accepted_step_count += int(limited)
    _accumulate_v15_position_diagnostics(
        positions,
        corrected_positions,
        start_evaluation,
        corrected,
        alive,
        dt,
        spacing,
        local_diagnostics,
        use_numba=bool(context.use_numba),
    )
    accepted_alive, termination, exit_time = _lifecycle_after_endpoint(
        interval,
        start_state,
        corrected,
    )
    accepted_vessel_id = _vessel_ids_after_endpoint(
        interval,
        start_state,
        corrected_topology,
        alive,
        ids,
        local_diagnostics,
    )
    local_diagnostics.outlet_event_count += int(
        np.count_nonzero(alive & ~accepted_alive)
    )
    return _BatchAttemptResult(
        _accepted_batch_state(
            start_state,
            corrected_positions,
            corrected_angles,
            new_count,
            new_extension,
            alive=accepted_alive,
            vessel_id=accepted_vessel_id,
            termination_reason=termination,
            exit_time_s=exit_time,
            **_step_observation_kwargs(
                corrected_generalized,
                corrected.contact_normal_xz,
                corrected,
                alive,
            ),
        ),
        local_diagnostics,
        float(interval.end_time_s),
    )


def _generalized_velocity(evaluation: object) -> np.ndarray:
    """Join the RHS translation and rotation into ``(Vx, Vz, Omega)``."""

    return np.column_stack(
        (
            np.asarray(evaluation.particle_velocity_xz_um_s, dtype=np.float64),
            np.asarray(evaluation.angular_velocity_rad_s, dtype=np.float64),
        )
    )


def _call_particle_rhs(
    evaluate_rhs: object,
    positions_grid: np.ndarray,
    bubble_ids: np.ndarray,
    active: np.ndarray,
    context: object,
    time_s: float,
    bond_count: np.ndarray,
    bond_extension: np.ndarray,
    vessel_id: np.ndarray,
) -> object:
    """Pass Revised-v20 ownership while retaining narrow legacy test adapters."""

    try:
        return evaluate_rhs(
            positions_grid,
            bubble_ids,
            active,
            context,
            time_s,
            bond_count,
            bond_extension,
            topological_vessel_id=vessel_id,
        )
    except TypeError as error:
        message = str(error)
        if "topological_vessel_id" not in message and "unexpected keyword" not in message:
            raise
        return evaluate_rhs(
            positions_grid,
            bubble_ids,
            active,
            context,
            time_s,
            bond_count,
            bond_extension,
        )


def _step_observation_kwargs(
    free_generalized_velocity: np.ndarray,
    inward_normals_xz: np.ndarray,
    step: PredictiveContactStep,
    alive: np.ndarray,
) -> dict[str, np.ndarray]:
    """Package the last accepted internal step for aligned saved-frame output."""

    live = np.asarray(alive, dtype=bool)
    free = np.asarray(free_generalized_velocity, dtype=np.float64)
    normal = np.asarray(inward_normals_xz, dtype=np.float64)
    constrained = np.asarray(
        step.constrained_generalized_velocity, dtype=np.float64
    )
    free_normal = np.sum(normal * free[:, :2], axis=1)
    constrained_normal = np.sum(normal * constrained[:, :2], axis=1)
    valid = live & np.all(np.isfinite(constrained), axis=1)
    return {
        "last_generalized_velocity": constrained.copy(),
        "last_contact_reaction_force_pn": np.where(
            valid, step.reaction_force_pn, 0.0
        ),
        "last_contact_active": valid & step.active_contact,
        "last_free_normal_velocity_um_s": np.where(valid, free_normal, 0.0),
        "last_constrained_normal_velocity_um_s": np.where(
            valid, constrained_normal, 0.0
        ),
        "last_step_valid": valid,
    }


def _solve_v15_trial_with_failure_context(
    positions: np.ndarray,
    alive: np.ndarray,
    dt_s: float,
    generalized_velocity: np.ndarray,
    mobility: np.ndarray,
    radii: np.ndarray,
    context: object,
    tolerance_um: float,
    *,
    permanent_ids: np.ndarray,
    interval: PhysicalTimeInterval,
    integration_stage: str,
    start_wall_distance_um: np.ndarray | None = None,
    start_wall_normal_xz: np.ndarray | None = None,
) -> PredictiveContactStep:
    """Attach permanent identity and physical time to a geometry failure."""

    try:
        cached_geometry_kwargs = (
            {
                "start_wall_distance_um": start_wall_distance_um,
                "start_wall_normal_xz": start_wall_normal_xz,
            }
            if _solve_v15_trial is _NATIVE_SOLVE_V15_TRIAL
            else {}
        )
        return _solve_v15_trial(
            positions,
            alive,
            dt_s,
            generalized_velocity,
            mobility,
            radii,
            context,
            tolerance_um,
            **cached_geometry_kwargs,
        )
    except PredictiveWallContactGeometryError as error:
        ids = np.asarray(permanent_ids, dtype=np.int64).reshape(-1)
        enriched: list[ParticleWallContactFailure] = []
        for failure in error.failures:
            lane = int(failure.local_lane)
            if lane < 0 or lane >= ids.size:
                raise RuntimeError(
                    "A v15 wall-contact diagnostic contains an invalid local "
                    f"lane {lane} for a batch of size {ids.size}."
                ) from error
            physical_time_s = float(interval.start_time_s) + (
                float(failure.event_fraction) * float(interval.duration_s)
            )
            enriched.append(
                ParticleWallContactFailure(
                    permanent_microbubble_id=int(ids[lane]),
                    physical_time_s=physical_time_s,
                    integration_stage=str(integration_stage),
                    contact=failure,
                )
            )
        raise ParticleWallContactGeometryError(
            str(error), tuple(enriched)
        ) from error
    except PredictiveBoundaryLifecycleError as error:
        ids = np.asarray(permanent_ids, dtype=np.int64).reshape(-1)
        enriched_boundary: list[ParticleBoundaryLifecycleFailure] = []
        for failure in error.failures:
            lane = int(failure.local_lane)
            if lane < 0 or lane >= ids.size:
                raise RuntimeError(
                    "A v16 boundary-lifecycle diagnostic contains an invalid "
                    f"local lane {lane} for a batch of size {ids.size}."
                ) from error
            enriched_boundary.append(
                ParticleBoundaryLifecycleFailure(
                    permanent_microbubble_id=int(ids[lane]),
                    physical_time_s=float(interval.end_time_s),
                    integration_stage=str(integration_stage),
                    local_lane=lane,
                    start_position_xz_um=failure.start_position_xz_um,
                    end_position_xz_um=failure.end_position_xz_um,
                    bubble_radius_um=float(failure.bubble_radius_um),
                    nearest_open_section_index=int(
                        failure.nearest_open_section_index
                    ),
                    nearest_open_section_kind=int(failure.nearest_open_section_kind),
                    nearest_open_section_label=int(
                        failure.nearest_open_section_label
                    ),
                    section_signed_start_um=float(failure.section_signed_start_um),
                    section_signed_end_um=float(failure.section_signed_end_um),
                    section_lateral_end_um=float(failure.section_lateral_end_um),
                    section_half_width_um=float(failure.section_half_width_um),
                )
            )
        raise ParticleBoundaryLifecycleError(tuple(enriched_boundary)) from error


def _solve_v15_trial(
    positions: np.ndarray,
    alive: np.ndarray,
    dt_s: float,
    generalized_velocity: np.ndarray,
    mobility: np.ndarray,
    radii: np.ndarray,
    context: object,
    tolerance_um: float,
    *,
    start_wall_distance_um: np.ndarray | None = None,
    start_wall_normal_xz: np.ndarray | None = None,
) -> PredictiveContactStep:
    return solve_predictive_contact_step(
        positions,
        alive,
        dt_s,
        generalized_velocity,
        mobility,
        radii,
        grid_spacing_um=float(context.spacing_um),
        boundary_geometry=context.boundary_geometry,
        geometry_tolerance_um=float(tolerance_um),
        use_numba=bool(context.use_numba),
        inputs_prevalidated=True,
        start_wall_distance_um=start_wall_distance_um,
        start_wall_normal_xz=start_wall_normal_xz,
    )


_NATIVE_SOLVE_V15_TRIAL = _solve_v15_trial


def _format_wall_contact_failure_message(
    message: str,
    failures: tuple[ParticleWallContactFailure, ...],
) -> str:
    """Format diagnostics as stable key/value lines in the final traceback."""

    lines = [str(message), "Wall-contact failure diagnostics:"]
    for record_index, record in enumerate(failures):
        failure = record.contact
        lines.extend(
            (
                f"failure[{record_index}]:",
                "  permanent_microbubble_id="
                f"{record.permanent_microbubble_id}",
                f"  physical_time_s={_format_float(record.physical_time_s)}",
                f"  integration_stage={record.integration_stage}",
                "  interval_start_position_xz_um="
                f"{_format_optional_pair(record.interval_start_position_xz_um)}",
                "  rejected_trial_displacement_um="
                f"{_format_optional_float(record.rejected_trial_displacement_um)}",
                f"  local_lane={failure.local_lane}",
                "  event_fraction="
                f"{_format_float(failure.event_fraction)}",
                "  position_grid="
                f"{_format_pair(failure.position_grid)}",
                "  position_xz_um="
                f"{_format_pair(failure.position_xz_um)}",
                f"  gap_um={_format_float(failure.gap_um)}",
                "  bubble_radius_um="
                f"{_format_float(failure.bubble_radius_um)}",
                "  multiple_wall_contact="
                f"{str(failure.multiple_wall_contact).lower()}",
                f"  reason={failure.reason}",
                "  candidate_wall_source="
                f"{failure.candidate_wall_source}",
                "  candidate_wall_count="
                f"{len(failure.candidate_walls)}",
                "  candidate_walls:",
            )
        )
        if not failure.candidate_walls:
            lines.append("    []")
            continue
        for wall in failure.candidate_walls:
            lines.append(
                "    - solid_face_index="
                f"{wall.solid_face_index}, "
                "center_xz_um="
                f"{_format_optional_pair(wall.center_xz_um)}, "
                f"length_um={_format_optional_float(wall.length_um)}, "
                "inward_normal_xz="
                f"{_format_optional_pair(wall.inward_normal_xz)}, "
                "exact_distance_um="
                f"{_format_float(wall.exact_distance_um)}, "
                f"exact_gap_um={_format_float(wall.exact_gap_um)}"
            )
    return "\n".join(lines)


def _format_float(value: float) -> str:
    return format(float(value), ".17g")


def _format_pair(values: tuple[float, float]) -> str:
    return f"[{_format_float(values[0])}, {_format_float(values[1])}]"


def _format_optional_pair(
    values: tuple[int, int] | tuple[float, float] | None,
) -> str:
    if values is None:
        return "unavailable"
    if all(isinstance(value, (int, np.integer)) for value in values):
        return f"[{int(values[0])}, {int(values[1])}]"
    return f"[{_format_float(values[0])}, {_format_float(values[1])}]"


def _format_optional_float(value: float | None) -> str:
    return "unavailable" if value is None else _format_float(value)


def _raise_for_v15_refinement(
    step: PredictiveContactStep, *, integration_stage: str
) -> None:
    failed = np.flatnonzero(step.need_time_refinement)
    if failed.size:
        raise _RefineContactTimeStep(
            "The predictive mobility constraint or complete trial chord could "
            "not be certified; retry two physical half intervals.",
            failed_lanes=failed,
            failure_codes=step.failure_codes[failed],
            failure_positions_grid=step.failure_position_grid[failed],
            failure_event_fractions=step.failure_event_fraction[failed],
            integration_stage=integration_stage,
        )


def _inspect_topological_step(
    start_positions_grid: np.ndarray,
    step: PredictiveContactStep,
    current_vessel_id: np.ndarray,
    alive: np.ndarray,
    context: object,
) -> TopologicalCrossingBatch:
    """Inspect the accepted centre chord without consulting a raster label."""

    catalog = getattr(context, "topological_ownership", None)
    if catalog is None:
        count = int(np.asarray(current_vessel_id).size)
        return TopologicalCrossingBatch(
            fraction=np.full(count, np.nan, dtype=np.float64),
            new_vessel_id=np.asarray(current_vessel_id, dtype=np.int32).copy(),
            section_index=np.full(count, -1, dtype=np.int32),
            position_xz_um=np.full((count, 2), np.nan, dtype=np.float64),
        )
    geometry = context.boundary_geometry
    start_world = np.asarray(
        geometry.grid_to_world_xz(start_positions_grid), dtype=np.float64
    )
    end_world = np.asarray(
        geometry.grid_to_world_xz(step.accepted_position_grid), dtype=np.float64
    )
    crossing = inspect_topological_crossings(
        start_world,
        end_world,
        current_vessel_id,
        alive,
        catalog,
        use_numba=bool(context.use_numba),
    )
    chord_fraction = np.asarray(crossing.fraction, dtype=np.float64)
    outlet_fraction = np.asarray(step.outlet_fraction, dtype=np.float64)
    endpoint_time_fraction = np.where(
        np.isfinite(outlet_fraction), outlet_fraction, 1.0
    )
    physical_fraction = chord_fraction * endpoint_time_fraction

    # The outlet is authoritative when a topology plane meets the terminal
    # opening at the same physical instant.  Such a lane terminates without a
    # redundant ownership transition at the zero-measure outlet rim.
    tie_scale = np.maximum(np.abs(endpoint_time_fraction), 1.0)
    outlet_tie = (
        np.isfinite(physical_fraction)
        & np.isfinite(outlet_fraction)
        & (
            np.abs(physical_fraction - outlet_fraction)
            <= 128.0 * np.finfo(np.float64).eps * tie_scale
        )
    )
    if np.any(outlet_tie):
        physical_fraction = physical_fraction.copy()
        physical_fraction[outlet_tie] = np.nan
        new_owner = np.asarray(crossing.new_vessel_id, dtype=np.int32).copy()
        new_owner[outlet_tie] = np.asarray(current_vessel_id, dtype=np.int32)[
            outlet_tie
        ]
        section_index = np.asarray(crossing.section_index, dtype=np.int32).copy()
        section_index[outlet_tie] = -1
    else:
        new_owner = crossing.new_vessel_id
        section_index = crossing.section_index
    return TopologicalCrossingBatch(
        fraction=np.asarray(physical_fraction, dtype=np.float64),
        new_vessel_id=np.asarray(new_owner, dtype=np.int32),
        section_index=np.asarray(section_index, dtype=np.int32),
        position_xz_um=np.asarray(crossing.position_xz_um, dtype=np.float64),
    )


def _split_for_internal_events(
    interval: PhysicalTimeInterval,
    step: PredictiveContactStep,
    topology: TopologicalCrossingBatch,
) -> None:
    """Split at the earliest positive topology or outlet event in the batch."""

    outlet_fraction = step.earliest_outlet_fraction
    topology_fraction = topology.earliest_fraction
    if outlet_fraction is None:
        earliest = topology_fraction
    elif topology_fraction is None:
        earliest = outlet_fraction
    else:
        tolerance = 128.0 * np.finfo(np.float64).eps * max(
            abs(float(outlet_fraction)), abs(float(topology_fraction)), 1.0
        )
        earliest = (
            float(outlet_fraction)
            if float(outlet_fraction) <= float(topology_fraction) + tolerance
            else float(topology_fraction)
        )
    if earliest is None:
        return
    time_tolerance = _time_epsilon(interval.start_time_s, interval.end_time_s)
    event_duration = float(earliest) * float(interval.duration_s)
    if (
        event_duration > time_tolerance
        and event_duration < float(interval.duration_s) - time_tolerance
    ):
        raise _SplitAtBoundaryEvent(float(interval.start_time_s) + event_duration)


def _commit_immediate_topology(
    interval: PhysicalTimeInterval,
    start_state: _BatchAdvanceState,
    topology: TopologicalCrossingBatch,
    positions: np.ndarray,
    angles: np.ndarray,
    bond_count: np.ndarray,
    bond_extension: np.ndarray,
    alive: np.ndarray,
    permanent_ids: np.ndarray,
) -> _BatchAttemptResult | None:
    """Commit zero-time directed crossings before re-evaluating the remainder."""

    epsilon = _time_epsilon(interval.start_time_s, interval.end_time_s)
    fractions = np.asarray(topology.fraction, dtype=np.float64)
    immediate = np.asarray(alive, dtype=bool) & np.isfinite(fractions) & (
        fractions * float(interval.duration_s) <= epsilon
    )
    if not np.any(immediate):
        return None
    vessel_id = np.asarray(start_state.vessel_id, dtype=np.int32).copy()
    vessel_id[immediate] = topology.new_vessel_id[immediate]
    diagnostics = _Diagnostics()
    _append_topological_events(
        diagnostics,
        immediate,
        permanent_ids,
        start_state.vessel_id,
        vessel_id,
        topology,
        np.full(fractions.shape, float(interval.start_time_s), dtype=np.float64),
    )
    state = _accepted_batch_state(
        start_state,
        positions.copy(),
        angles.copy(),
        bond_count.copy(),
        bond_extension.copy(),
        vessel_id=vessel_id,
    )
    return _BatchAttemptResult(state, diagnostics, float(interval.start_time_s))


def _vessel_ids_after_endpoint(
    interval: PhysicalTimeInterval,
    start_state: _BatchAdvanceState,
    topology: TopologicalCrossingBatch,
    alive: np.ndarray,
    permanent_ids: np.ndarray,
    diagnostics: _Diagnostics | None,
) -> np.ndarray:
    """Apply only crossings that coincide with the accepted interval endpoint."""

    fractions = np.asarray(topology.fraction, dtype=np.float64)
    epsilon = _time_epsilon(interval.start_time_s, interval.end_time_s)
    endpoint = np.asarray(alive, dtype=bool) & np.isfinite(fractions) & (
        (1.0 - fractions) * float(interval.duration_s) <= epsilon
    )
    if not np.any(endpoint):
        return np.asarray(start_state.vessel_id, dtype=np.int32)
    vessel_id = np.asarray(start_state.vessel_id, dtype=np.int32).copy()
    vessel_id[endpoint] = topology.new_vessel_id[endpoint]
    if diagnostics is not None:
        event_time = (
            float(interval.start_time_s)
            + fractions * float(interval.duration_s)
        )
        _append_topological_events(
            diagnostics,
            endpoint,
            permanent_ids,
            start_state.vessel_id,
            vessel_id,
            topology,
            event_time,
        )
    return vessel_id


def _append_topological_events(
    diagnostics: _Diagnostics,
    event_mask: np.ndarray,
    permanent_ids: np.ndarray,
    previous_vessel_id: np.ndarray,
    new_vessel_id: np.ndarray,
    topology: TopologicalCrossingBatch,
    event_time_s: np.ndarray,
) -> None:
    lanes = np.flatnonzero(event_mask)
    diagnostics.topological_transition_count += int(lanes.size)
    for lane_value in lanes:
        lane = int(lane_value)
        diagnostics.topological_event_bubble_id.append(int(permanent_ids[lane]))
        diagnostics.topological_event_time_s.append(float(event_time_s[lane]))
        diagnostics.topological_event_from_vessel_id.append(
            int(previous_vessel_id[lane])
        )
        diagnostics.topological_event_to_vessel_id.append(int(new_vessel_id[lane]))
        diagnostics.topological_event_section_index.append(
            int(topology.section_index[lane])
        )
        position = topology.position_xz_um[lane]
        diagnostics.topological_event_position_xz_um.append(
            (float(position[0]), float(position[1]))
        )


def _split_for_internal_outlet(
    interval: PhysicalTimeInterval,
    step: PredictiveContactStep,
) -> None:
    fraction = step.earliest_outlet_fraction
    if fraction is None:
        return
    time_tolerance = _time_epsilon(interval.start_time_s, interval.end_time_s)
    event_duration = float(fraction) * float(interval.duration_s)
    if (
        event_duration > time_tolerance
        and event_duration < float(interval.duration_s) - time_tolerance
    ):
        raise _SplitAtBoundaryEvent(
            float(interval.start_time_s) + event_duration
        )


def _commit_immediate_outlets(
    interval: PhysicalTimeInterval,
    start_state: _BatchAdvanceState,
    step: PredictiveContactStep,
    positions: np.ndarray,
    angles: np.ndarray,
    bond_count: np.ndarray,
    bond_extension: np.ndarray,
    alive: np.ndarray,
) -> _BatchAttemptResult | None:
    epsilon = _time_epsilon(interval.start_time_s, interval.end_time_s)
    fractions = np.asarray(step.outlet_fraction, dtype=np.float64)
    immediate = alive & np.isfinite(fractions) & (
        fractions * float(interval.duration_s) <= epsilon
    )
    if not np.any(immediate):
        return None
    accepted_positions = positions.copy()
    accepted_positions[immediate] = step.outlet_position_grid[immediate]
    accepted_alive = alive.copy()
    accepted_alive[immediate] = False
    termination = np.asarray(
        start_state.termination_reason, dtype=np.uint8
    ).copy()
    termination[immediate] = _TERMINATION_OUTLET
    exit_time = np.asarray(start_state.exit_time_s, dtype=np.float64).copy()
    exit_time[immediate] = float(interval.start_time_s)
    state = _accepted_batch_state(
        start_state,
        accepted_positions,
        angles.copy(),
        bond_count.copy(),
        bond_extension.copy(),
        alive=accepted_alive,
        termination_reason=termination,
        exit_time_s=exit_time,
    )
    return _BatchAttemptResult(
        state,
        _Diagnostics(),
        float(interval.start_time_s),
    )


def _lifecycle_after_endpoint(
    interval: PhysicalTimeInterval,
    start_state: _BatchAdvanceState,
    step: PredictiveContactStep,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    alive = np.asarray(start_state.alive, dtype=bool)
    fractions = np.asarray(step.outlet_fraction, dtype=np.float64)
    exiting = alive & np.isfinite(fractions)
    if not np.any(exiting):
        # The overwhelmingly common path has no lifecycle event.  These arrays
        # are immutable batch state, so reuse them and copy only when an outlet
        # event actually needs a write.
        return (
            alive,
            np.asarray(start_state.termination_reason, dtype=np.uint8),
            np.asarray(start_state.exit_time_s, dtype=np.float64),
        )
    alive = alive.copy()
    termination = np.asarray(
        start_state.termination_reason, dtype=np.uint8
    ).copy()
    exit_time = np.asarray(start_state.exit_time_s, dtype=np.float64).copy()
    alive[exiting] = False
    termination[exiting] = _TERMINATION_OUTLET
    exit_time[exiting] = (
        float(interval.start_time_s)
        + fractions[exiting] * float(interval.duration_s)
    )
    return alive, termination, exit_time


def _accumulate_v15_position_diagnostics(
    start_positions: np.ndarray,
    end_positions: np.ndarray,
    evaluation: object,
    step: PredictiveContactStep,
    alive: np.ndarray,
    dt_s: float,
    spacing_um: float,
    diagnostics: _Diagnostics,
    *,
    use_numba: bool,
) -> None:
    live = np.asarray(alive, dtype=bool)
    projected = np.asarray(step.residual_correction_um, dtype=np.float64)[live]
    positive_projection = projected[np.isfinite(projected) & (projected > 0.0)]
    diagnostics.contact_residual_projection_count += int(
        positive_projection.size
    )
    if positive_projection.size:
        diagnostics.maximum_contact_residual_projection_um = max(
            diagnostics.maximum_contact_residual_projection_um,
            float(np.max(positive_projection)),
        )
    if use_numba:
        reduced = reduce_predictive_step_diagnostics(
            start_positions,
            end_positions,
            step.constrained_generalized_velocity[:, :2],
            np.asarray(evaluation.wall_gap_um, dtype=np.float64),
            step.endpoint_gap_um,
            np.asarray(step.contact_normal_xz, dtype=np.float64),
            live,
            step.active_contact,
            step.reaction_force_pn,
            step.minimum_path_gap_um,
            step.complementarity_residual_pn_um,
            spacing_um=float(spacing_um),
            duration_s=float(dt_s),
            use_numba=True,
            inputs_prevalidated=True,
        )
        diagnostics.contact_constraint_evaluations += reduced.live_count
        diagnostics.active_contact_constraint_evaluations += (
            reduced.active_contact_count
        )
        if np.isfinite(reduced.maximum_contact_reaction_force_pn):
            diagnostics.maximum_contact_reaction_force_pn = max(
                diagnostics.maximum_contact_reaction_force_pn,
                reduced.maximum_contact_reaction_force_pn,
            )
        if np.isfinite(reduced.minimum_path_gap_um):
            diagnostics.minimum_accepted_internal_wall_gap_um = min(
                diagnostics.minimum_accepted_internal_wall_gap_um,
                reduced.minimum_path_gap_um,
            )
        if np.isfinite(reduced.maximum_complementarity_residual_pn_um):
            diagnostics.maximum_contact_complementarity_residual_pn_um = max(
                diagnostics.maximum_contact_complementarity_residual_pn_um,
                reduced.maximum_complementarity_residual_pn_um,
            )
        diagnostics.contact_kinematic_interval_evaluations += (
            reduced.contact_interval_count
        )
        diagnostics.contact_cumulative_position_path_um += (
            reduced.contact_position_path_sum_um
        )
        diagnostics.contact_cumulative_velocity_path_um += (
            reduced.contact_velocity_path_sum_um
        )
        diagnostics.contact_nonzero_velocity_zero_progress_count += (
            reduced.contact_nonzero_velocity_zero_progress_count
        )
        if np.isfinite(reduced.minimum_contact_path_ratio):
            diagnostics.minimum_contact_interval_position_to_velocity_path_ratio = min(
                diagnostics.minimum_contact_interval_position_to_velocity_path_ratio,
                reduced.minimum_contact_path_ratio,
            )
            diagnostics.maximum_contact_interval_position_to_velocity_path_ratio = max(
                diagnostics.maximum_contact_interval_position_to_velocity_path_ratio,
                reduced.maximum_contact_path_ratio,
            )
        if np.isfinite(reduced.maximum_free_gap_residual_um):
            diagnostics.maximum_free_gap_kinematic_residual_um = max(
                diagnostics.maximum_free_gap_kinematic_residual_um,
                reduced.maximum_free_gap_residual_um,
            )
        if reduced.live_count:
            diagnostics.maximum_step_displacement_um = max(
                diagnostics.maximum_step_displacement_um,
                reduced.maximum_live_displacement_um,
            )
        return

    diagnostics.contact_constraint_evaluations += int(np.count_nonzero(live))
    diagnostics.active_contact_constraint_evaluations += int(
        np.count_nonzero(live & step.active_contact)
    )
    finite_reaction = step.reaction_force_pn[
        live & np.isfinite(step.reaction_force_pn)
    ]
    if finite_reaction.size:
        diagnostics.maximum_contact_reaction_force_pn = max(
            diagnostics.maximum_contact_reaction_force_pn,
            float(np.max(finite_reaction)),
        )
    finite_gap = step.minimum_path_gap_um[
        live & np.isfinite(step.minimum_path_gap_um)
    ]
    if finite_gap.size:
        diagnostics.minimum_accepted_internal_wall_gap_um = min(
            diagnostics.minimum_accepted_internal_wall_gap_um,
            float(np.min(finite_gap)),
        )
    finite_complementarity = np.abs(
        step.complementarity_residual_pn_um[
            live & np.isfinite(step.complementarity_residual_pn_um)
        ]
    )
    if finite_complementarity.size:
        diagnostics.maximum_contact_complementarity_residual_pn_um = max(
            diagnostics.maximum_contact_complementarity_residual_pn_um,
            float(np.max(finite_complementarity)),
        )
    start_gap = np.asarray(evaluation.wall_gap_um, dtype=np.float64)
    end_gap = np.asarray(step.endpoint_gap_um, dtype=np.float64).copy()
    metric_lanes = live & np.isfinite(start_gap) & np.isfinite(end_gap)
    gap_scale = np.maximum.reduce(
        (
            np.ones(start_gap.shape, dtype=np.float64),
            np.abs(start_gap),
            np.abs(end_gap),
        )
    )
    contact_roundoff = 256.0 * np.finfo(np.float64).eps * gap_scale
    contact_lanes = metric_lanes & (
        np.asarray(step.active_contact, dtype=bool)
        | (start_gap <= contact_roundoff)
        | (end_gap <= contact_roundoff)
    )
    contact_lane_count = int(np.count_nonzero(contact_lanes))
    diagnostics.contact_kinematic_interval_evaluations += contact_lane_count
    free_gap_lanes = metric_lanes & ~np.asarray(step.active_contact, dtype=bool)
    if contact_lane_count:
        contact_metrics = evaluate_step_kinematics(
            start_positions,
            end_positions,
            step.constrained_generalized_velocity[:, :2],
            start_gap,
            end_gap,
            np.asarray(step.contact_normal_xz, dtype=np.float64),
            contact_lanes,
            spacing_um=float(spacing_um),
            duration_s=float(dt_s),
        )
        diagnostics.contact_cumulative_position_path_um += float(
            np.nansum(contact_metrics.position_path_um)
        )
        diagnostics.contact_cumulative_velocity_path_um += float(
            np.nansum(contact_metrics.velocity_path_um)
        )
        diagnostics.contact_nonzero_velocity_zero_progress_count += int(
            np.count_nonzero(contact_metrics.nonzero_velocity_zero_progress)
        )
        finite_contact_ratio = contact_metrics.position_to_velocity_path_ratio[
            np.isfinite(contact_metrics.position_to_velocity_path_ratio)
        ]
        if finite_contact_ratio.size:
            diagnostics.minimum_contact_interval_position_to_velocity_path_ratio = min(
                diagnostics.minimum_contact_interval_position_to_velocity_path_ratio,
                float(np.min(finite_contact_ratio)),
            )
            diagnostics.maximum_contact_interval_position_to_velocity_path_ratio = max(
                diagnostics.maximum_contact_interval_position_to_velocity_path_ratio,
                float(np.max(finite_contact_ratio)),
            )

    free_gap_metrics = evaluate_step_kinematics(
        start_positions,
        end_positions,
        step.constrained_generalized_velocity[:, :2],
        start_gap,
        end_gap,
        np.asarray(evaluation.wall_normal_xz, dtype=np.float64),
        free_gap_lanes,
        spacing_um=float(spacing_um),
        duration_s=float(dt_s),
    )
    residual = np.abs(free_gap_metrics.gap_kinematic_residual_um)
    if np.any(np.isfinite(residual)):
        diagnostics.maximum_free_gap_kinematic_residual_um = max(
            diagnostics.maximum_free_gap_kinematic_residual_um,
            float(np.nanmax(residual)),
        )
    displacement = np.linalg.norm(
        (np.asarray(end_positions) - np.asarray(start_positions))
        * float(spacing_um),
        axis=1,
    )
    if np.any(live):
        diagnostics.maximum_step_displacement_um = max(
            diagnostics.maximum_step_displacement_um,
            float(np.max(displacement[live])),
        )


def _accept_euler_bonds(
    bond_count: np.ndarray,
    bond_extension: np.ndarray,
    evaluation: object,
    start_positions: np.ndarray,
    accepted_positions: np.ndarray,
    start_angles: np.ndarray,
    accepted_angles: np.ndarray,
    radii_um: np.ndarray,
    moving: np.ndarray,
    spacing_um: float,
    dt_s: float,
    use_numba: bool,
) -> tuple[np.ndarray, np.ndarray, int]:
    molecular_rhs = evaluation.molecular_binding_evaluation
    if molecular_rhs is None or dt_s <= 0.0:
        return bond_count, bond_extension, 0
    new_count = np.asarray(bond_count, dtype=np.float64).copy()
    new_extension = np.asarray(bond_extension, dtype=np.float64).copy()
    source = _accepted_surface_extension_source_um_s(
        start_positions,
        accepted_positions,
        start_angles,
        accepted_angles,
        radii_um,
        bond_count,
        evaluation.wall_normal_xz,
        None,
        spacing_um,
        dt_s,
    )
    candidate = predict_bond_state_exponential_euler(
        bond_count,
        bond_extension,
        molecular_rhs,
        dt_s,
        extension_source_um_s=source,
        use_numba=use_numba,
    )
    new_count[moving] = candidate.expected_bond_count[moving]
    new_extension[moving] = (
        candidate.total_tangential_extension_um[moving]
    )
    limited = int(np.count_nonzero(candidate.capacity_limited[moving]))
    return new_count, new_extension, limited


def _accept_heun_bonds(
    bond_count: np.ndarray,
    bond_extension: np.ndarray,
    predictor_bond_count: np.ndarray,
    start_evaluation: object,
    predicted_evaluation: object,
    start_positions: np.ndarray,
    accepted_positions: np.ndarray,
    start_angles: np.ndarray,
    accepted_angles: np.ndarray,
    radii_um: np.ndarray,
    moving: np.ndarray,
    spacing_um: float,
    dt_s: float,
    use_numba: bool,
) -> tuple[np.ndarray, np.ndarray, int]:
    start_rhs = start_evaluation.molecular_binding_evaluation
    predicted_rhs = predicted_evaluation.molecular_binding_evaluation
    if start_rhs is None:
        return bond_count, bond_extension, 0
    if predicted_rhs is None:
        raise RuntimeError("The molecular Heun predictor did not produce a bond RHS.")
    new_count = np.asarray(bond_count, dtype=np.float64).copy()
    new_extension = np.asarray(bond_extension, dtype=np.float64).copy()
    source = _accepted_surface_extension_source_um_s(
        start_positions,
        accepted_positions,
        start_angles,
        accepted_angles,
        radii_um,
        0.5 * (bond_count + predictor_bond_count),
        start_evaluation.wall_normal_xz,
        predicted_evaluation.wall_normal_xz,
        spacing_um,
        dt_s,
    )
    candidate = accept_bond_state_exponential_heun(
        bond_count,
        bond_extension,
        start_rhs,
        predicted_rhs,
        dt_s,
        extension_source_um_s=source,
        use_numba=use_numba,
    )
    new_count[moving] = candidate.expected_bond_count[moving]
    new_extension[moving] = (
        candidate.total_tangential_extension_um[moving]
    )
    limited = int(np.count_nonzero(candidate.capacity_limited[moving]))
    return new_count, new_extension, limited


def _merge_diagnostics(target: _Diagnostics, source: _Diagnostics) -> None:
    """Merge accepted diagnostics; rejected speculative values never reach here."""

    # This function runs once per accepted internal particle step.  Explicit
    # assignments avoid repeatedly building a dictionary, inspecting every
    # field name, and doing dynamic attribute lookups in this hot loop.
    target.maximum_physical_overlap_um = max(target.maximum_physical_overlap_um, source.maximum_physical_overlap_um)
    target.maximum_collision_compression_um = max(target.maximum_collision_compression_um, source.maximum_collision_compression_um)
    target.maximum_interacting_pairs = max(target.maximum_interacting_pairs, source.maximum_interacting_pairs)
    target.total_interacting_pair_evaluations += source.total_interacting_pair_evaluations
    target.all_pairs_rhs_evaluations += source.all_pairs_rhs_evaluations
    target.cell_list_rhs_evaluations += source.cell_list_rhs_evaluations
    target.maximum_reciprocity_relative_error = max(target.maximum_reciprocity_relative_error, source.maximum_reciprocity_relative_error)
    target.degenerate_near_wall_normal_evaluations += source.degenerate_near_wall_normal_evaluations
    target.maximum_collision_speed_um_s = max(target.maximum_collision_speed_um_s, source.maximum_collision_speed_um_s)
    target.maximum_step_displacement_um = max(target.maximum_step_displacement_um, source.maximum_step_displacement_um)
    target.two_wall_observations += source.two_wall_observations
    target.inlet_wait_events += source.inlet_wait_events
    target.maximum_inlet_wait_s = max(target.maximum_inlet_wait_s, source.maximum_inlet_wait_s)
    target.maximum_expected_bond_count = max(target.maximum_expected_bond_count, source.maximum_expected_bond_count)
    target.maximum_single_bond_tension_pn = max(target.maximum_single_bond_tension_pn, source.maximum_single_bond_tension_pn)
    target.maximum_total_bond_force_pn = max(target.maximum_total_bond_force_pn, source.maximum_total_bond_force_pn)
    target.maximum_bond_torque_pn_um = max(target.maximum_bond_torque_pn_um, source.maximum_bond_torque_pn_um)
    target.molecular_bell_saturation_count += source.molecular_bell_saturation_count
    target.molecular_formation_saturation_count += source.molecular_formation_saturation_count
    target.molecular_capacity_limited_accepted_step_count += source.molecular_capacity_limited_accepted_step_count
    target.contact_constraint_evaluations += source.contact_constraint_evaluations
    target.active_contact_constraint_evaluations += source.active_contact_constraint_evaluations
    target.maximum_contact_reaction_force_pn = max(target.maximum_contact_reaction_force_pn, source.maximum_contact_reaction_force_pn)
    target.contact_time_refinement_count += source.contact_time_refinement_count
    target.maximum_contact_time_refinement_depth = max(target.maximum_contact_time_refinement_depth, source.maximum_contact_time_refinement_depth)
    target.minimum_accepted_internal_wall_gap_um = min(target.minimum_accepted_internal_wall_gap_um, source.minimum_accepted_internal_wall_gap_um)
    target.contact_nonzero_velocity_zero_progress_count += source.contact_nonzero_velocity_zero_progress_count
    target.contact_residual_projection_count += source.contact_residual_projection_count
    target.maximum_contact_residual_projection_um = max(
        target.maximum_contact_residual_projection_um,
        source.maximum_contact_residual_projection_um,
    )
    target.maximum_contact_complementarity_residual_pn_um = max(target.maximum_contact_complementarity_residual_pn_um, source.maximum_contact_complementarity_residual_pn_um)
    target.contact_kinematic_interval_evaluations += source.contact_kinematic_interval_evaluations
    target.contact_cumulative_position_path_um += source.contact_cumulative_position_path_um
    target.contact_cumulative_velocity_path_um += source.contact_cumulative_velocity_path_um
    target.minimum_contact_interval_position_to_velocity_path_ratio = min(
        target.minimum_contact_interval_position_to_velocity_path_ratio,
        source.minimum_contact_interval_position_to_velocity_path_ratio,
    )
    target.maximum_contact_interval_position_to_velocity_path_ratio = max(
        target.maximum_contact_interval_position_to_velocity_path_ratio,
        source.maximum_contact_interval_position_to_velocity_path_ratio,
    )
    target.maximum_free_gap_kinematic_residual_um = max(target.maximum_free_gap_kinematic_residual_um, source.maximum_free_gap_kinematic_residual_um)
    target.outlet_event_count += source.outlet_event_count
    target.topological_transition_count += source.topological_transition_count
    target.topological_event_bubble_id.extend(source.topological_event_bubble_id)
    target.topological_event_time_s.extend(source.topological_event_time_s)
    target.topological_event_from_vessel_id.extend(
        source.topological_event_from_vessel_id
    )
    target.topological_event_to_vessel_id.extend(
        source.topological_event_to_vessel_id
    )
    target.topological_event_section_index.extend(
        source.topological_event_section_index
    )
    target.topological_event_position_xz_um.extend(
        source.topological_event_position_xz_um
    )
    target.active_outside_lumen_violations += source.active_outside_lumen_violations
    target.active_outside_accessible_domain_violations += source.active_outside_accessible_domain_violations
    target.discrete_accessibility_disagreement_records += (
        source.discrete_accessibility_disagreement_records
    )
    target.maximum_red_blood_cell_speed_um_s = max(
        target.maximum_red_blood_cell_speed_um_s,
        source.maximum_red_blood_cell_speed_um_s,
    )
    target.red_blood_cell_rhs_observations += source.red_blood_cell_rhs_observations
    target.red_blood_cell_nonunique_wall_suppressions += (
        source.red_blood_cell_nonunique_wall_suppressions
    )
    target.red_blood_cell_hematocrit_out_of_range_observations += (
        source.red_blood_cell_hematocrit_out_of_range_observations
    )
    target.red_blood_cell_shear_out_of_range_observations += (
        source.red_blood_cell_shear_out_of_range_observations
    )
    target.red_blood_cell_quantitative_applicability_observations += (
        source.red_blood_cell_quantitative_applicability_observations
    )
    target.red_blood_cell_random_wall_reflection_count += (
        source.red_blood_cell_random_wall_reflection_count
    )
    target.maximum_red_blood_cell_random_displacement_um = max(
        target.maximum_red_blood_cell_random_displacement_um,
        source.maximum_red_blood_cell_random_displacement_um,
    )
    target.red_blood_cell_random_displacement_squared_sum_um2 += (
        source.red_blood_cell_random_displacement_squared_sum_um2
    )
    target.red_blood_cell_random_displacement_observations += (
        source.red_blood_cell_random_displacement_observations
    )
    target.red_blood_cell_diffusion_enabled_observations += (
        source.red_blood_cell_diffusion_enabled_observations
    )
    target.red_blood_cell_stochastic_observations += (
        source.red_blood_cell_stochastic_observations
    )
    target.red_blood_cell_invalid_transverse_space_observations += (
        source.red_blood_cell_invalid_transverse_space_observations
    )


def _time_epsilon(start_time_s: float, end_time_s: float) -> float:
    scale = max(abs(float(start_time_s)), abs(float(end_time_s)), 1.0)
    return 32.0 * np.finfo(float).eps * scale


def _accepted_surface_extension_source_um_s(
    start_positions_grid: np.ndarray,
    accepted_positions_grid: np.ndarray,
    start_angles_rad: np.ndarray,
    accepted_angles_rad: np.ndarray,
    radii_um: np.ndarray,
    representative_bond_count: np.ndarray,
    start_wall_normal_xz: np.ndarray,
    end_wall_normal_xz: np.ndarray | None,
    spacing_um: float,
    dt_s: float,
) -> np.ndarray:
    """Convert the accepted position-level motion into the bond-extension source.

    The mobility RHS proposes a velocity, but finite-radius contact may hold or
    redirect the centre before the step is accepted.  Molecular extension must
    therefore use the realized tangential centre displacement plus the accepted
    rotational surface displacement, rather than the unconstrained proposal.
    """

    normals = np.asarray(start_wall_normal_xz, dtype=np.float64)
    if end_wall_normal_xz is not None:
        normals = normals + np.asarray(end_wall_normal_xz, dtype=np.float64)
    norm = np.linalg.norm(normals, axis=1)
    valid = np.isfinite(norm) & (norm > 1.0e-12)
    unit_normals = np.zeros_like(normals)
    unit_normals[valid] = normals[valid] / norm[valid, None]
    unit_normals[~valid, 1] = 1.0
    tangents = np.column_stack((-unit_normals[:, 1], unit_normals[:, 0]))
    displacement_um = (
        np.asarray(accepted_positions_grid, dtype=np.float64)
        - np.asarray(start_positions_grid, dtype=np.float64)
    ) * float(spacing_um)
    tangential_speed_um_s = np.sum(displacement_um * tangents, axis=1) / float(dt_s)
    angular_speed_rad_s = (
        np.asarray(accepted_angles_rad, dtype=np.float64)
        - np.asarray(start_angles_rad, dtype=np.float64)
    ) / float(dt_s)
    mean_slip_um_s = tangential_speed_um_s + np.asarray(
        radii_um, dtype=np.float64
    ) * angular_speed_rad_s
    return np.asarray(representative_bond_count, dtype=np.float64) * mean_slip_um_s


def _admit_waiting(
    time_s: float,
    state: _PerfusionState,
    schedule: PerfusionSchedule,
    context: object,
    diagnostics: _Diagnostics,
) -> int:
    if not state.waiting_ids:
        return 0
    _ensure_last_step_storage(state)
    _ensure_vessel_id_storage(state)
    admitted = 0
    still_waiting: list[int] = []
    for bubble_id in state.waiting_ids:
        candidate = np.asarray(schedule.position_grid[bubble_id], dtype=np.float64).copy()
        radius_um = float(context.all_bubble_radii_um[bubble_id])
        candidate_world_xz, candidate_gap_um, _ = (
            project_roundoff_negative_wall_gap_xz_um(
                context.boundary_geometry.grid_to_world_xz(candidate),
                radius_um,
                context.boundary_geometry,
            )
        )
        candidate = np.asarray(
            context.boundary_geometry.world_xz_to_grid(candidate_world_xz),
            dtype=np.float64,
        )
        active_ids = _active_ids_view(state)
        can_enter = True
        if active_ids.size:
            distances_um = np.linalg.norm(state.position_grid[active_ids] - candidate, axis=1) * float(
                context.spacing_um
            )
            required = context.all_bubble_radii_um[active_ids] + context.all_bubble_radii_um[bubble_id]
            can_enter = bool(np.all(distances_um >= required - 1.0e-12))
        if not can_enter:
            still_waiting.append(bubble_id)
            continue
        if not bool(context.boundary_geometry.contains_grid(candidate)):
            raise RuntimeError(
                f"Scheduled inlet position for bubble {bubble_id} is outside the canonical lumen."
            )
        if not bool(
            context.boundary_geometry.is_accessible_grid(candidate, radius_um)
        ):
            raise RuntimeError(
                f"Scheduled inlet position for bubble {bubble_id} is not in Omega_R^in."
            )
        if not np.isfinite(candidate_gap_um) or candidate_gap_um < 0.0:
            raise RuntimeError(
                f"Scheduled inlet position for bubble {bubble_id} has a negative true wall gap."
            )
        state.position_grid[bubble_id] = candidate
        state.vessel_id[bubble_id] = int(schedule.initial_vessel_id[bubble_id])
        state.rotation_angle_rad[bubble_id] = 0.0
        state.bond_count_expected[bubble_id] = 0.0
        state.bond_total_tangential_extension_um[bubble_id] = 0.0
        state.last_generalized_velocity[bubble_id] = 0.0
        state.last_contact_reaction_force_pn[bubble_id] = 0.0
        state.last_contact_active[bubble_id] = False
        state.last_free_normal_velocity_um_s[bubble_id] = 0.0
        state.last_constrained_normal_velocity_um_s[bubble_id] = 0.0
        state.last_step_valid[bubble_id] = False
        _append_active_id(state, int(bubble_id))
        state.active[bubble_id] = True
        state.admission_time_s[bubble_id] = time_s
        state.admission_event_ids.append(int(bubble_id))
        wait = max(float(time_s - schedule.planned_time_s[bubble_id]), 0.0)
        diagnostics.maximum_inlet_wait_s = max(diagnostics.maximum_inlet_wait_s, wait)
        admitted += 1
    state.waiting_ids = still_waiting
    return admitted


def _record_current_frame(
    state: _PerfusionState,
    context: object,
    particle_cfg: ParticleConfig,
    records: _FrameRecords,
    diagnostics: _Diagnostics,
    evaluate_rhs: object,
    time_s: float,
) -> None:
    ids = _active_ids_view(state)
    count = int(ids.size)
    if count == 0:
        records.offsets.append(records.offsets[-1])
        return
    _ensure_vessel_id_storage(state)
    positions = np.asarray(state.position_grid[ids], dtype=np.float64)
    active = np.ones(count, dtype=bool)
    evaluation = _call_particle_rhs(
        evaluate_rhs,
        positions,
        ids,
        active,
        context,
        time_s,
        state.bond_count_expected[ids],
        state.bond_total_tangential_extension_um[ids],
        state.vessel_id[ids],
    )
    _accumulate_diagnostics(evaluation, diagnostics)
    geometry = context.boundary_geometry
    inside = np.asarray(geometry.contains_grid(positions), dtype=bool)
    diagnostics.active_outside_lumen_violations += int(
        np.count_nonzero(~inside)
    )
    radii = np.asarray(context.all_bubble_radii_um[ids], dtype=np.float64)
    accessible = np.asarray(
        geometry.is_accessible_grid(positions, radii), dtype=bool
    )
    continuous_safe = (
        np.isfinite(evaluation.wall_gap_um)
        & (
            evaluation.wall_gap_um
            >= -float(particle_cfg.contact_geometry_tolerance_um)
        )
    )
    diagnostics.discrete_accessibility_disagreement_records += int(
        np.count_nonzero(inside & continuous_safe & ~accessible)
    )
    _ensure_last_step_storage(state)
    saved_valid = np.asarray(state.last_step_valid[ids], dtype=bool)
    saved_generalized = np.asarray(
        state.last_generalized_velocity[ids], dtype=np.float64
    )
    displayed_velocity = np.asarray(
        evaluation.particle_velocity_xz_um_s, dtype=np.float64
    ).copy()
    displayed_angular_velocity = np.asarray(
        evaluation.angular_velocity_rad_s, dtype=np.float64
    ).copy()
    displayed_velocity[saved_valid] = saved_generalized[saved_valid, :2]
    displayed_angular_velocity[saved_valid] = saved_generalized[saved_valid, 2]
    records.position_grid.append(positions.copy())
    records.particle_velocity_xz.append(displayed_velocity)
    records.wall_shear.append(evaluation.sampled_wall_shear_stress_pa.copy())
    records.vessel_id.append(np.asarray(state.vessel_id[ids], dtype=np.int32).copy())
    records.bubble_id.append(ids.copy())
    records.wall_gap.append(evaluation.wall_gap_um.copy())
    records.wall_contact.append(
        evaluation.wall_gap_um <= particle_cfg.wall_contact_threshold_um
    )
    records.wall_normal.append(evaluation.wall_normal_xz.copy())
    if context.dynamics.store_full_diagnostics:
        records.fluid_velocity_xz.append(evaluation.fluid_velocity_xz_um_s.copy())
        records.angular_velocity.append(displayed_angular_velocity)
        records.rotation_angle.append(state.rotation_angle_rad[ids].copy())
        records.collision_force.append(evaluation.collision_force_xz_pn.copy())
        records.collision_neighbors.append(evaluation.collision_neighbor_count.copy())
        records.gap_ratio.append(evaluation.gap_ratio.copy())
        records.near_wall_weight.append(evaluation.near_wall_weight.copy())
        records.two_wall_warning.append(evaluation.two_wall_warning.copy())
    if context.cardiac is not None:
        records.cardiac_multiplier.append(evaluation.cardiac_multiplier.copy())
    records.contact_constraint_active.append(
        np.asarray(state.last_contact_active[ids], dtype=bool).copy()
    )
    records.contact_reaction_force_pn.append(
        np.asarray(
            state.last_contact_reaction_force_pn[ids], dtype=np.float64
        ).copy()
    )
    records.contact_free_normal_velocity_um_s.append(
        np.asarray(
            state.last_free_normal_velocity_um_s[ids], dtype=np.float64
        ).copy()
    )
    records.contact_constrained_normal_velocity_um_s.append(
        np.asarray(
            state.last_constrained_normal_velocity_um_s[ids], dtype=np.float64
        ).copy()
    )
    if evaluation.bond_count_expected is not None:
        records.bond_count_expected.append(evaluation.bond_count_expected.copy())
        records.bond_total_tangential_extension_um.append(
            evaluation.bond_total_tangential_extension_um.copy()
        )
        records.bond_mean_tangential_extension_um.append(
            evaluation.bond_mean_tangential_extension_um.copy()
        )
        records.bond_force_xz_pn.append(evaluation.bond_force_xz_pn.copy())
        records.bond_force_tangent_pn.append(evaluation.bond_force_tangent_pn.copy())
        records.bond_force_normal_pn.append(evaluation.bond_force_normal_pn.copy())
        records.bond_torque_pn_um.append(evaluation.bond_torque_pn_um.copy())
        records.single_bond_tension_pn.append(evaluation.single_bond_tension_pn.copy())
        records.bond_formation_rate_bonds_s.append(
            evaluation.bond_formation_rate_bonds_s.copy()
        )
        records.bond_dissociation_rate_s_inv.append(
            evaluation.bond_dissociation_rate_s_inv.copy()
        )
        records.target_reaction_area_um2.append(evaluation.target_reaction_area_um2.copy())
        records.available_ligand_count.append(evaluation.available_ligand_count.copy())
        records.available_target_count.append(evaluation.available_target_count.copy())
        records.target_overlap_fraction.append(evaluation.target_overlap_fraction.copy())
    if evaluation.red_blood_cell_velocity_xz_um_s is not None:
        records.red_blood_cell_velocity_xz_um_s.append(
            evaluation.red_blood_cell_velocity_xz_um_s.copy()
        )
        records.red_blood_cell_drift_velocity_xz_um_s.append(
            evaluation.red_blood_cell_drift_velocity_xz_um_s.copy()
        )
        records.red_blood_cell_fick_velocity_xz_um_s.append(
            evaluation.red_blood_cell_fick_velocity_xz_um_s.copy()
        )
        # The target gap remains an internal saved-frame statistic used for the
        # run-level CFL encounter fractions.  It is exposed in the NPZ only
        # when the existing detailed-diagnostics switch is enabled.
        records.red_blood_cell_target_gap_um.append(
            evaluation.red_blood_cell_target_gap_um.copy()
        )
        records.red_blood_cell_transverse_diffusivity_um2_s.append(
            evaluation.red_blood_cell_transverse_diffusivity_um2_s.copy()
        )
        records.red_blood_cell_quantitative_applicability.append(
            evaluation.red_blood_cell_quantitative_applicability.copy()
        )
        records.red_blood_cell_transverse_space_valid.append(
            evaluation.red_blood_cell_transverse_space_valid.copy()
        )
        if context.dynamics.store_full_diagnostics:
            records.red_blood_cell_local_vessel_diameter_um.append(
                evaluation.red_blood_cell_local_vessel_diameter_um.copy()
            )
            records.red_blood_cell_discharge_hematocrit.append(
                evaluation.red_blood_cell_discharge_hematocrit.copy()
            )
            records.red_blood_cell_tube_hematocrit.append(
                evaluation.red_blood_cell_tube_hematocrit.copy()
            )
            records.red_blood_cell_shear_rate_s_inv.append(
                evaluation.red_blood_cell_shear_rate_s_inv.copy()
            )
            records.red_blood_cell_cfl_width_um.append(
                evaluation.red_blood_cell_cfl_width_um.copy()
            )
            records.red_blood_cell_margination_length_um.append(
                evaluation.red_blood_cell_margination_length_um.copy()
            )
            records.red_blood_cell_margination_time_s.append(
                evaluation.red_blood_cell_margination_time_s.copy()
            )
            records.red_blood_cell_scale_activation.append(
                evaluation.red_blood_cell_scale_activation.copy()
            )
            records.red_blood_cell_nearest_wall_unique.append(
                evaluation.red_blood_cell_nearest_wall_unique.copy()
            )
            records.red_blood_cell_hematocrit_in_quantitative_range.append(
                evaluation.red_blood_cell_hematocrit_in_quantitative_range.copy()
            )
            records.red_blood_cell_shear_rate_in_quantitative_range.append(
                evaluation.red_blood_cell_shear_rate_in_quantitative_range.copy()
            )
    records.offsets.append(records.offsets[-1] + count)


def _accumulate_diagnostics(evaluation: object, diagnostics: _Diagnostics) -> None:
    diagnostics.maximum_physical_overlap_um = max(
        diagnostics.maximum_physical_overlap_um,
        float(evaluation.maximum_physical_overlap_um),
    )
    diagnostics.maximum_collision_compression_um = max(
        diagnostics.maximum_collision_compression_um,
        float(evaluation.maximum_collision_compression_um),
    )
    diagnostics.maximum_interacting_pairs = max(
        diagnostics.maximum_interacting_pairs, int(evaluation.interacting_pair_count)
    )
    diagnostics.total_interacting_pair_evaluations += int(evaluation.interacting_pair_count)
    if evaluation.collision_search_strategy == "all_pairs":
        diagnostics.all_pairs_rhs_evaluations += 1
    elif evaluation.collision_search_strategy == "cell_list":
        diagnostics.cell_list_rhs_evaluations += 1
    diagnostics.maximum_reciprocity_relative_error = max(
        diagnostics.maximum_reciprocity_relative_error,
        float(evaluation.maximum_reciprocity_relative_error),
    )
    diagnostics.degenerate_near_wall_normal_evaluations += int(
        evaluation.degenerate_near_wall_normal_count
    )
    diagnostics.maximum_collision_speed_um_s = max(
        diagnostics.maximum_collision_speed_um_s,
        float(evaluation.maximum_collision_speed_um_s),
    )
    diagnostics.two_wall_observations += int(np.count_nonzero(evaluation.two_wall_warning))
    red_blood_cell_velocity = getattr(
        evaluation, "red_blood_cell_velocity_xz_um_s", None
    )
    if red_blood_cell_velocity is not None:
        valid = np.isfinite(evaluation.red_blood_cell_local_vessel_diameter_um)
        diagnostics.red_blood_cell_rhs_observations += int(np.count_nonzero(valid))
        if np.any(valid):
            diagnostics.maximum_red_blood_cell_speed_um_s = max(
                diagnostics.maximum_red_blood_cell_speed_um_s,
                float(
                    np.max(
                        np.linalg.norm(
                            red_blood_cell_velocity[valid], axis=1
                        )
                    )
                ),
            )
        diagnostics.red_blood_cell_nonunique_wall_suppressions += int(
            np.count_nonzero(valid & ~evaluation.red_blood_cell_nearest_wall_unique)
        )
        diagnostics.red_blood_cell_hematocrit_out_of_range_observations += int(
            np.count_nonzero(
                valid
                & ~evaluation.red_blood_cell_hematocrit_in_quantitative_range
            )
        )
        diagnostics.red_blood_cell_shear_out_of_range_observations += int(
            np.count_nonzero(
                valid & ~evaluation.red_blood_cell_shear_rate_in_quantitative_range
            )
        )
        diagnostics.red_blood_cell_quantitative_applicability_observations += int(
            np.count_nonzero(
                valid & evaluation.red_blood_cell_quantitative_applicability
            )
        )
    if evaluation.bond_count_expected is not None:
        if evaluation.bond_count_expected.size:
            diagnostics.maximum_expected_bond_count = max(
                diagnostics.maximum_expected_bond_count,
                float(np.max(evaluation.bond_count_expected)),
            )
            diagnostics.maximum_single_bond_tension_pn = max(
                diagnostics.maximum_single_bond_tension_pn,
                float(np.max(evaluation.single_bond_tension_pn)),
            )
            diagnostics.maximum_total_bond_force_pn = max(
                diagnostics.maximum_total_bond_force_pn,
                float(np.max(np.linalg.norm(evaluation.bond_force_xz_pn, axis=1))),
            )
            diagnostics.maximum_bond_torque_pn_um = max(
                diagnostics.maximum_bond_torque_pn_um,
                float(np.max(np.abs(evaluation.bond_torque_pn_um))),
            )
        diagnostics.molecular_bell_saturation_count += int(
            np.count_nonzero(evaluation.molecular_bell_rate_saturated)
        )
        diagnostics.molecular_formation_saturation_count += int(
            np.count_nonzero(evaluation.molecular_formation_rate_saturated)
        )




def _check_record_limit(records: _FrameRecords, cfg: ParticleConfig) -> None:
    limit = int(cfg.max_particle_frame_records)
    if limit > 0 and records.offsets[-1] > limit:
        raise RuntimeError(
            f"Continuous perfusion produced {records.offsets[-1]} particle-frame records, "
            f"exceeding particles.max_particle_frame_records={limit}."
        )


def _concat(arrays: list[np.ndarray], empty_shape: tuple[int, ...], dtype: object) -> np.ndarray:
    if not arrays:
        return np.empty(empty_shape, dtype=dtype)
    return np.asarray(np.concatenate(arrays, axis=0), dtype=dtype)


def _embed_xz(values_xz: np.ndarray) -> np.ndarray:
    values = np.asarray(values_xz)
    result = np.zeros((values.shape[0], 3), dtype=np.float32)
    result[:, 0] = values[:, 0]
    result[:, 2] = values[:, 1]
    return result
