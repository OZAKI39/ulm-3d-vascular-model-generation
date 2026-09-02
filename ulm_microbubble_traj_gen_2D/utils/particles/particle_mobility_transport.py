"""Finite-size mobility right-hand-side evaluation.

The continuous-perfusion driver owns particle lifecycle, recording, and time
integration.  This module only evaluates the coupled hydrodynamic/collision/
molecular right-hand side.  Contact handling is owned by the dedicated
unilateral-contact stepping modules.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..cardiac.cardiac_pulsatility import CardiacPulsatility
from ..core.config import ParticleDynamicsConfig
from ..geometry.continuous_vessel_geometry import ContinuousVesselGeometry
from ..core.types import HybridVelocityField
from ..flow.hybrid_velocity import (
    blend_velocity_and_gradient,
    sample_finite_element_velocity,
)
from ..molecular.molecular_binding import (
    evaluate_mean_field_bonds,
    reaction_disk_radius_um,
    surface_slip_velocity_um_s,
)
from .particle_collisions import (
    compute_collision_forces,
    compute_collision_forces_prevalidated_numba,
)
from .particle_field_sampling import (
    sample_bilinear,
    sample_regular_grid_fields_batch,
)
from .particle_mobility import (
    FREE_ROTATIONAL_SCALED_MOBILITY,
    apply_scaled_mobility_scalar,
    background_hydrodynamic_velocity_scalar,
    bulk_angular_velocity_scalar,
    effective_dimensionless_mobility_entries_scalar,
    generalized_mobility_matrix_xz,
    translational_stokes_mobility_scalar,
    wall_to_global_components_scalar,
)
from .particle_numba_kernels import (
    apply_particle_loads_batch,
    build_generalized_mobility_batch,
    evaluate_hydrodynamic_mobility_batch,
)
from .particle_geometry_tolerance import (
    particle_geometry_roundoff_tolerance_um,
)
from .red_blood_cell_transport import evaluate_red_blood_cell_transport


@dataclass(frozen=True, slots=True)
class ParticleVesselOwnershipFailure:
    """One active bubble whose persistent Revised-v20 owner is invalid."""

    permanent_microbubble_id: int
    physical_time_s: float
    local_lane: int
    position_grid: tuple[float, float]
    position_xz_um: tuple[float, float]
    persistent_topological_vessel_id: int


class ParticleVesselOwnershipError(ValueError):
    """Persistent topological-owner failure with identity, time, and position."""

    def __init__(self, failures: tuple[ParticleVesselOwnershipFailure, ...]) -> None:
        self.failures = failures
        super().__init__(_format_vessel_ownership_failures(failures))


@dataclass(frozen=True, slots=True)
class ParticleWallGapFailure:
    """One active bubble whose authoritative centre state penetrates a wall."""

    permanent_microbubble_id: int
    physical_time_s: float
    local_lane: int
    position_grid: tuple[float, float]
    position_xz_um: tuple[float, float]
    bubble_radius_um: float
    wall_gap_um: float
    geometry_roundoff_tolerance_um: float


class ParticleWallGapInvariantError(ValueError):
    """An RHS start state has a wall penetration larger than roundoff."""

    def __init__(self, failures: tuple[ParticleWallGapFailure, ...]) -> None:
        self.failures = failures
        super().__init__(_format_wall_gap_failures(failures))


@dataclass(frozen=True)
class _EvaluationContext:
    """Read-only inputs shared by all right-hand-side evaluations.

    Keeping these large Eulerian arrays in one immutable context makes it clear
    that an RHS evaluation only *reads* the solved flow field.  It also avoids
    repeatedly converting field dtypes or memory layouts inside the frame loop.
    """

    velocity_xz_um_s: np.ndarray
    wall_shear_stress_pa: np.ndarray
    vessel_id: np.ndarray
    local_vessel_radius_um: np.ndarray
    velocity_gradient_s_inv: np.ndarray
    dynamic_viscosity_pa_s: np.ndarray
    all_bubble_radii_um: np.ndarray
    spacing_um: float
    dynamics: ParticleDynamicsConfig
    contact_geometry_tolerance_um: float
    use_numba: bool
    boundary_geometry: ContinuousVesselGeometry
    hybrid_velocity: HybridVelocityField
    random_seed: int = 42
    cardiac: CardiacPulsatility | None = None
    origin_xz_um: tuple[float, float] = (0.0, 0.0)
    molecular_target_field: object | None = None
    molecular_binding_parameters: object | None = None
    molecular_capture_distance_um: float = 0.0
    molecular_mean_field_warning_count: float = 10.0
    red_blood_cell_network: object | None = None
    topological_ownership: object | None = None


@dataclass(frozen=True)
class _RhsEvaluation:
    """Velocities and diagnostics calculated for one particle snapshot.

    The result contains both the particle velocity used by the integrator and
    the sampled carrier-fluid velocity.  They need not be equal: near-wall
    mobility and bubble-bubble forces can change translation and rotation.
    """

    particle_velocity_xz_um_s: np.ndarray
    fluid_velocity_xz_um_s: np.ndarray
    angular_velocity_rad_s: np.ndarray
    generalized_mobility: np.ndarray
    collision_force_xz_pn: np.ndarray
    collision_neighbor_count: np.ndarray
    wall_gap_um: np.ndarray
    hydrodynamic_gap_um: np.ndarray
    wall_normal_xz: np.ndarray
    gap_ratio: np.ndarray
    near_wall_weight: np.ndarray
    two_wall_warning: np.ndarray
    sampled_wall_shear_stress_pa: np.ndarray
    topological_vessel_id: np.ndarray
    cardiac_multiplier: np.ndarray
    maximum_physical_overlap_um: float
    maximum_collision_compression_um: float
    interacting_pair_count: int
    collision_search_strategy: str
    maximum_reciprocity_relative_error: float
    degenerate_near_wall_normal_count: int
    maximum_collision_speed_um_s: float
    bond_count_expected: np.ndarray | None = None
    bond_total_tangential_extension_um: np.ndarray | None = None
    bond_mean_tangential_extension_um: np.ndarray | None = None
    bond_force_xz_pn: np.ndarray | None = None
    bond_force_tangent_pn: np.ndarray | None = None
    bond_force_normal_pn: np.ndarray | None = None
    bond_torque_pn_um: np.ndarray | None = None
    single_bond_tension_pn: np.ndarray | None = None
    bond_formation_rate_bonds_s: np.ndarray | None = None
    bond_dissociation_rate_s_inv: np.ndarray | None = None
    target_reaction_area_um2: np.ndarray | None = None
    available_ligand_count: np.ndarray | None = None
    available_target_count: np.ndarray | None = None
    target_overlap_fraction: np.ndarray | None = None
    molecular_bell_rate_saturated: np.ndarray | None = None
    molecular_formation_rate_saturated: np.ndarray | None = None
    molecular_binding_evaluation: object | None = None
    wall_normal_valid: np.ndarray | None = None
    red_blood_cell_velocity_xz_um_s: np.ndarray | None = None
    red_blood_cell_drift_velocity_xz_um_s: np.ndarray | None = None
    red_blood_cell_fick_velocity_xz_um_s: np.ndarray | None = None
    red_blood_cell_transverse_direction_xz: np.ndarray | None = None
    red_blood_cell_local_vessel_diameter_um: np.ndarray | None = None
    red_blood_cell_discharge_hematocrit: np.ndarray | None = None
    red_blood_cell_tube_hematocrit: np.ndarray | None = None
    red_blood_cell_shear_rate_s_inv: np.ndarray | None = None
    red_blood_cell_cfl_width_um: np.ndarray | None = None
    red_blood_cell_target_gap_um: np.ndarray | None = None
    red_blood_cell_center_gap_um: np.ndarray | None = None
    red_blood_cell_core_transverse_diffusivity_um2_s: np.ndarray | None = None
    red_blood_cell_transverse_diffusivity_um2_s: np.ndarray | None = None
    red_blood_cell_diffusion_taper_coordinate: np.ndarray | None = None
    red_blood_cell_diffusion_taper: np.ndarray | None = None
    red_blood_cell_margination_length_um: np.ndarray | None = None
    red_blood_cell_margination_time_s: np.ndarray | None = None
    red_blood_cell_scale_activation: np.ndarray | None = None
    red_blood_cell_nearest_wall_unique: np.ndarray | None = None
    red_blood_cell_hematocrit_in_quantitative_range: np.ndarray | None = None
    red_blood_cell_shear_rate_in_quantitative_range: np.ndarray | None = None
    red_blood_cell_quantitative_applicability: np.ndarray | None = None
    red_blood_cell_transverse_space_valid: np.ndarray | None = None
    red_blood_cell_diffusion_enabled: np.ndarray | None = None

def _evaluate_rhs(
    positions_grid: np.ndarray,
    bubble_ids: np.ndarray,
    active: np.ndarray,
    context: _EvaluationContext,
    time_s: float = 0.0,
    bond_count_expected: np.ndarray | None = None,
    bond_total_tangential_extension_um: np.ndarray | None = None,
    topological_vessel_id: np.ndarray | None = None,
) -> _RhsEvaluation:
    """
    Evaluate background motion, effective mobility, and pair forces.
    "RHS" means the time derivative used by Euler or Heun: translational and angular particle velocity.
    """

    # =================================================================================
    # ===== Read the current local flow field and wall geometry at every bubble centre. 
    # =================================================================================
    positions   = np.ascontiguousarray(positions_grid, dtype=np.float64)
    ids         = np.ascontiguousarray(bubble_ids, dtype=np.int64)
    active_mask = np.ascontiguousarray(active, dtype=bool)
    count       = int(ids.size)
    all_active  = bool(count == 0 or np.all(active_mask))
    radii       = np.asarray(context.all_bubble_radii_um[ids], dtype=np.float64)
    sampled_vessel_id = (
        np.full(count, -1, dtype=np.int32)
        if topological_vessel_id is None
        else (
            np.ascontiguousarray(topological_vessel_id, dtype=np.int32)
            if all_active
            else np.ascontiguousarray(topological_vessel_id, dtype=np.int32).copy()
        )
    )
    if sampled_vessel_id.shape != (count,):
        raise ValueError("Persistent topological vessel IDs must match the active particle lanes.")

    bond_count      = (np.zeros(count, dtype=np.float64) if bond_count_expected is None 
                       else np.ascontiguousarray(bond_count_expected, dtype=np.float64))
    bond_extension  = (np.zeros(count, dtype=np.float64) if bond_total_tangential_extension_um is None 
                       else np.ascontiguousarray(bond_total_tangential_extension_um, dtype=np.float64))

    if bond_count.shape != (count,) or bond_extension.shape != (count,):
        raise ValueError("Molecular bond state arrays must match the active particle lanes.")
    
    # The continuous wall is queried first because its exact distance chooses
    # the velocity representation. The Cartesian mask never makes this choice.
    nearest_wall_unique = np.zeros(count, dtype=bool)
    wall_distance = np.zeros(count, dtype=np.float64)
    sampled_normals = np.zeros((count, 2), dtype=np.float64)
    world_positions = np.column_stack(
        (
            float(context.origin_xz_um[0])
            + positions[:, 0] * float(context.spacing_um),
            float(context.origin_xz_um[1])
            + positions[:, 1] * float(context.spacing_um),
        )
    )
    if np.any(active_mask):
        active_lanes = None if all_active else np.flatnonzero(active_mask)
        active_world = world_positions if all_active else world_positions[active_lanes]
        exact_wall = (
            context.boundary_geometry.exact_solid_wall_state_xz_um_accelerated(
                active_world
            )
            if context.use_numba
            else context.boundary_geometry.exact_solid_wall_state_xz_um(
                active_world
            )
        )
        exact_distance = np.asarray(
            exact_wall.distance_um, dtype=np.float64
        ).reshape(-1)
        exact_normals = np.asarray(
            exact_wall.inward_normal_xz, dtype=np.float64
        ).reshape(-1, 2)
        exact_unique = np.asarray(
            exact_wall.unique_nearest_wall,
            dtype=bool,
        ).reshape(-1)
        if all_active:
            wall_distance = exact_distance
            sampled_normals = exact_normals
            nearest_wall_unique = exact_unique
        else:
            wall_distance[active_lanes] = exact_distance
            sampled_normals[active_lanes] = exact_normals
            nearest_wall_unique[active_lanes] = exact_unique
    sampled_normal_norm = np.linalg.norm(sampled_normals, axis=1)
    wall_normal_valid = (
        np.isfinite(sampled_normal_norm)
        & (sampled_normal_norm > 1.0e-12)
    )
    geometry_roundoff = particle_geometry_roundoff_tolerance_um(
        world_positions[active_mask],
        radii[active_mask],
        context.spacing_um,
    )
    authoritative_gap = wall_distance - radii
    materially_negative = active_mask & (
        authoritative_gap < -geometry_roundoff
    )
    if np.any(materially_negative):
        raise ParticleWallGapInvariantError(
            tuple(
                ParticleWallGapFailure(
                    permanent_microbubble_id=int(ids[lane]),
                    physical_time_s=float(time_s),
                    local_lane=int(lane),
                    position_grid=(
                        float(positions[lane, 0]),
                        float(positions[lane, 1]),
                    ),
                    position_xz_um=(
                        float(world_positions[lane, 0]),
                        float(world_positions[lane, 1]),
                    ),
                    bubble_radius_um=float(radii[lane]),
                    wall_gap_um=float(authoritative_gap[lane]),
                    geometry_roundoff_tolerance_um=float(geometry_roundoff),
                )
                for lane in np.flatnonzero(materially_negative)
            )
        )
    roundoff_negative = active_mask & (authoritative_gap < 0.0)
    if np.any(roundoff_negative):
        # Geometry transactions project these residuals at their source. This
        # guard also makes old/reused states obey one downstream gap contract.
        wall_distance = wall_distance.copy()
        wall_distance[roundoff_negative] = radii[roundoff_negative]

    # Auxiliary material fields always use the Cartesian cache. Velocity and
    # gradient use that cache only in the far field; every point with non-zero
    # finite-element weight is also evaluated in its DOLFINx triangle.
    if context.use_numba:
        (
            grid_velocity,
            grid_gradient,
            viscosity,
            local_vessel_radius,
            base_sampled_shear,
        ) = sample_regular_grid_fields_batch(
            positions,
            active_mask,
            context.velocity_xz_um_s,
            context.velocity_gradient_s_inv,
            context.dynamic_viscosity_pa_s,
            context.local_vessel_radius_um,
            context.wall_shear_stress_pa,
            use_numba=True,
        )
    else:
        grid_velocity = sample_bilinear(context.velocity_xz_um_s, positions)
        grid_gradient = sample_bilinear(context.velocity_gradient_s_inv, positions)
        viscosity = sample_bilinear(context.dynamic_viscosity_pa_s, positions)
        local_vessel_radius = sample_bilinear(
            context.local_vessel_radius_um, positions
        )
        base_sampled_shear = sample_bilinear(
            context.wall_shear_stress_pa, positions
        )

    finite_element_required = active_mask & (
        wall_distance < float(context.hybrid_velocity.regular_grid_distance_um)
    )
    finite_element_velocity, finite_element_gradient, _ = (
        sample_finite_element_velocity(
            context.hybrid_velocity.finite_element,
            world_positions,
            finite_element_required,
            use_numba=context.use_numba,
        )
    )
    base_fluid_velocity, base_gradients, _ = blend_velocity_and_gradient(
        grid_velocity,
        grid_gradient,
        finite_element_velocity,
        finite_element_gradient,
        wall_distance,
        sampled_normals,
        context.hybrid_velocity,
    )
    base_fluid_velocity[~active_mask] = 0.0
    base_gradients[~active_mask] = 0.0

    if context.cardiac is None:
        cardiac_multiplier  = np.ones(count, dtype=np.float64)
        fluid_velocity      = base_fluid_velocity
        gradients           = base_gradients
        sampled_shear       = base_sampled_shear
    else:
        # With heartbeat modulation, ask how strong the pulse is at this place and time. 
        # The pulse may arrive later downstream because it travels through the vessel.
        cardiac_sample      = context.cardiac.sample(
            positions, float(time_s), use_numba=context.use_numba
        )
        cardiac_multiplier  = np.asarray(cardiac_sample.multiplier, dtype=np.float64)
        fluid_velocity      = base_fluid_velocity * cardiac_multiplier[:, None]          # Apply the local pulse strength to the background fluid velocity.

        # A velocity gradient says how velocity changes over a small distance.
        gradients = (base_gradients * cardiac_multiplier[:, None, None] + base_fluid_velocity[:, :, None] * cardiac_sample.gradient_per_um[:, None, :])

        # The steady wall-shear value is adjusted by the same local pulse strength.
        sampled_shear = base_sampled_shear * cardiac_multiplier

    maximum_reciprocity_error       = 0.0
    degenerate_near_wall_normals    = 0

    if context.use_numba:
        (
            particle_velocity,
            angular_velocity,
            normalized_normals,
            wall_gap,
            gap_ratio,
            wall_weight,
            two_wall_warning,
            mobility_entries,
            translation_mobility_global,
            reciprocity_error_by_lane,
            degenerate_normal_by_lane,
        ) = evaluate_hydrodynamic_mobility_batch(
            fluid_velocity,
            gradients,
            viscosity,
            wall_distance,
            sampled_normals,
            local_vessel_radius,
            radii,
            active_mask,
            near_wall_enabled=bool(context.dynamics.near_wall_enabled),
            xi_min=float(context.dynamics.xi_min),
            xi_near=float(context.dynamics.xi_near),
            xi_far=float(context.dynamics.xi_far),
            two_wall_warning_gap_ratio=float(
                context.dynamics.two_wall_warning_gap_ratio
            ),
        )
        maximum_reciprocity_error = (
            float(np.max(reciprocity_error_by_lane))
            if reciprocity_error_by_lane.size
            else 0.0
        )
        degenerate_near_wall_normals = int(
            np.count_nonzero(degenerate_normal_by_lane)
        )
        if not all_active:
            sampled_vessel_id[~active_mask] = -1
    else:
        # The Python reference path fills these arrays one lane at a time. The
        # compiled path above returns its own arrays and does not allocate them twice.
        particle_velocity = np.zeros((count, 2), dtype=np.float64)
        angular_velocity = np.zeros(count, dtype=np.float64)
        normalized_normals = np.zeros((count, 2), dtype=np.float64)
        wall_gap = np.full(count, np.nan, dtype=np.float64)
        gap_ratio = np.full(count, np.nan, dtype=np.float64)
        wall_weight = np.zeros(count, dtype=np.float64)
        two_wall_warning = np.zeros(count, dtype=bool)
        mobility_entries = np.zeros((count, 5), dtype=np.float64)
        translation_mobility_global = np.zeros((count, 2, 2), dtype=np.float64)

    # =================================================================================
    # Inspect one bubble at a time to build its background movement and its ability to respond to later forces. 
    # =================================================================================
    
    # Compute the effect on microbubbles of the local flow field, wall proximity, and near-wall mobility.
    for lane in range(0 if context.use_numba else count):
        # An inactive lane is just an empty seat. Mark its vessel label as invalid and skip physical calculations.
        if not active_mask[lane]:
            sampled_vessel_id[lane] = -1
            continue
        radius_um   = float(radii[lane])
        mu_pa_s     = float(viscosity[lane])

        # Viscosity tells us how strongly the liquid resists motion.
        if not np.isfinite(mu_pa_s) or mu_pa_s <= 0.0:
            raise ValueError("Sampled particle viscosity must be finite and positive.")
        
        # The oriented distance field is positive inside the lumen and zero at the true wall.
        # Subtracting the bubble radius gives the unregularized surface gap shared by geometry and capture; mobility regularizes only its coefficients later.
        gap_um          = float(wall_distance[lane] - radius_um)
        wall_gap[lane]  = gap_um
        gap_ratio[lane] = gap_um / radius_um

        # Read the direction pointing from the nearest wall into the blood.
        normal_x    = float(sampled_normals[lane, 0])
        normal_z    = float(sampled_normals[lane, 1])
        normal_norm = float(np.hypot(normal_x, normal_z))
        potential_wall_weight = _smooth_wall_weight(gap_ratio[lane], float(context.dynamics.xi_near), float(context.dynamics.xi_far))

        # At a vessel centre line, two walls can be equally close and the nearest-wall direction may be unclear. 
        # Use a fixed safe direction in this rare case.
        if normal_norm <= 1.0e-12 or not np.isfinite(normal_norm):
            normal_x = 0.0
            normal_z = 1.0
            if context.dynamics.near_wall_enabled and potential_wall_weight > 0.0:
                degenerate_near_wall_normals += 1
        else:
            normal_x /= normal_norm
            normal_z /= normal_norm

        normalized_normals[lane, 0] = normal_x
        normalized_normals[lane, 1] = normal_z

        if context.dynamics.near_wall_enabled:
            # A thin liquid layer between a bubble and a wall makes sliding and rolling harder. 
            # This helper includes that effect and smoothly returns to ordinary free-flow behaviour as the bubble moves farther away.
            # Near a wall, translation and rotation differ from passive fluid advection.  
            vx, vz, omega, weight, _ = background_hydrodynamic_velocity_scalar(
                float(fluid_velocity[lane, 0]), float(fluid_velocity[lane, 1]),
                float(gradients[lane, 0, 0]), float(gradients[lane, 0, 1]), float(gradients[lane, 1, 0]), float(gradients[lane, 1, 1]),
                normal_x, normal_z,
                radius_um,
                max(gap_um, 0.0),
                mu_pa_s,
                float(context.dynamics.xi_min), float(context.dynamics.xi_near), float(context.dynamics.xi_far),
            )
            entries = effective_dimensionless_mobility_entries_scalar(
                max(gap_um, 0.0) / radius_um,
                float(context.dynamics.xi_min), float(context.dynamics.xi_near), float(context.dynamics.xi_far),
            )
        else:
            # If the near-wall model is switched off, use the simple bulk-flow rule: 
            # follow the local liquid and rotate with its local swirl.
            vx      = float(fluid_velocity[lane, 0])
            vz      = float(fluid_velocity[lane, 1])
            omega   = bulk_angular_velocity_scalar(float(gradients[lane, 0, 1]), float(gradients[lane, 1, 0]))
            weight  = 0.0
            entries = (1.0, 1.0, 0.0, 0.0, FREE_ROTATIONAL_SCALED_MOBILITY)

        # Save the background translation, rotation, and wall influence for this bubble. 
        particle_velocity[lane, 0]  = vx
        particle_velocity[lane, 1]  = vz
        angular_velocity[lane]      = omega
        wall_weight[lane]           = weight
        mobility_entries[lane]      = entries

        # Two coupling values are expected to match because of a physical symmetry:
        # "a push causing rotation" and "a twist causing movement" must agree in
        # the matching reversed experiment. 
        denominator                 = max(abs(float(entries[2])), abs(float(entries[3])), 1.0e-30)
        maximum_reciprocity_error   = max(maximum_reciprocity_error, abs(float(entries[2]) - float(entries[3])) / denominator)

        # The wall formulas use two easy local directions: along the wall and away from the wall. 
        tangent_x = -normal_z
        tangent_z = normal_x

        base = translational_stokes_mobility_scalar(mu_pa_s, radius_um)
        m_tt = float(entries[0])
        m_nn = float(entries[1])

        translation_mobility_global[lane, 0, 0] = base * (m_tt * tangent_x * tangent_x + m_nn * normal_x * normal_x)
        translation_mobility_global[lane, 0, 1] = base * (m_tt * tangent_x * tangent_z + m_nn * normal_x * normal_z)
        translation_mobility_global[lane, 1, 0] = translation_mobility_global[lane, 0, 1]
        translation_mobility_global[lane, 1, 1] = base * (m_tt * tangent_z * tangent_z + m_nn * normal_z * normal_z)

        # The selected near-wall formula assumes that one wall is the main wall.
        # In a narrow vessel, the wall on the opposite side can also be close.
        opposite_gap_um         = 2.0 * float(local_vessel_radius[lane]) - float(wall_distance[lane]) - radius_um
        opposite_gap_ratio      = opposite_gap_um / radius_um
        two_wall_warning[lane]  = bool(context.dynamics.near_wall_enabled and weight > 0.0 and opposite_gap_ratio <= context.dynamics.two_wall_warning_gap_ratio)

    # Keep the geometric gap separate from the regularized gap used only by the
    # asymptotic mobility coefficients.  The latter prevents logarithmic and
    # lubrication singularities; it must never become a minimum accepted wall
    # distance or a molecular capture distance.
    hydrodynamic_gap = np.where(
        active_mask,
        np.maximum(wall_gap, radii * float(context.dynamics.xi_min)),
        np.nan,
    )

    # =================================================================================
    # ===== Evaluate bubble-bubble collisions and optional molecular binding.
    # =================================================================================
    if context.dynamics.collisions_enabled:
        # Collision distances must use micrometres, not grid-cell units.
        positions_um = positions * float(context.spacing_um)
        # Find nearby pairs and create gentle repulsive forces when bubbles become too close.
        if context.use_numba:
            (
                collision_force,
                collision_neighbors,
                maximum_physical_overlap,
                maximum_collision_compression,
                pair_count,
                collision_search_strategy,
            ) = compute_collision_forces_prevalidated_numba(
                positions_um,
                radii,
                translation_mobility_global,
                active_mask,
                ids,
                collision_layer_um=float(context.dynamics.collision_layer_um),
                relaxation_time_s=float(
                    context.dynamics.collision_relaxation_time_s
                ),
                strategy=str(context.dynamics.neighbor_search),
            )
        else:
            collision = compute_collision_forces(
                positions_um, radii, translation_mobility_global, active_mask, ids,
                collision_layer_um=float(context.dynamics.collision_layer_um),
                relaxation_time_s=float(context.dynamics.collision_relaxation_time_s),
                strategy=str(context.dynamics.neighbor_search), use_numba=False,
            )
            collision_force                 = collision.force_xz_pn
            collision_neighbors             = collision.neighbor_count
            maximum_physical_overlap        = collision.maximum_physical_overlap_um
            maximum_collision_compression   = collision.maximum_collision_compression_um
            pair_count                      = collision.interacting_pair_count
            collision_search_strategy       = collision.search_strategy
    else:
        # If collisions are disabled, create correctly shaped zero arrays. 
        # Zero force naturally means that collisions add no extra movement.
        collision_force                 = np.zeros((count, 2), dtype=np.float64)
        collision_neighbors             = np.zeros(count, dtype=np.int32)
        maximum_physical_overlap        = 0.0
        maximum_collision_compression   = 0.0
        pair_count                      = 0
        collision_search_strategy       = "disabled"

    # Start with empty molecular results.
    molecular_evaluation        = None
    target_reaction_area_um2    = None
    target_overlap_fraction     = None
    bond_force_xz_pn            = None
    if context.molecular_binding_parameters is not None:
        # Estimate the small contact/reaction disk that each finite-size bubble can present to the nearby wall. 
        # A closer or larger bubble can expose a different reaction area, which changes how many bonds may form.
        reaction_radius_um = reaction_disk_radius_um(
            radii, np.where(active_mask, wall_gap, context.molecular_capture_distance_um),
            float(context.molecular_capture_distance_um),
        )

        # Inactive lanes must not form bonds. 
        molecular_gap_um            = np.where(active_mask, wall_gap, 0.0)
        molecular_bond_count        = np.where(active_mask, bond_count, 0.0)
        molecular_bond_extension    = np.where(active_mask, bond_extension, 0.0)

        # Convert grid positions into real X-Z coordinates. 
        centres_xz_um       = np.empty((count, 2), dtype=np.float64)
        centres_xz_um[:, 0] = float(context.origin_xz_um[0]) + positions[:, 0] * float(context.spacing_um)
        centres_xz_um[:, 1] = float(context.origin_xz_um[1]) + positions[:, 1] * float(context.spacing_um)

        # Keep the mobility tangent for force/rotation bookkeeping. The target
        # field receives the particle centre itself and projects it onto exact
        # canonical solid faces; reconstructing a wall point from a bilinear
        # distance and gradient can select the wrong face at a raster corner.
        tangents_xz         = np.column_stack((-normalized_normals[:, 1], normalized_normals[:, 0]))

        # Ask the target map how much of each reaction disk actually overlaps a
        # target-positive wall. Bubbles outside the capture layer have a zero
        # reaction radius, so exclude them before the wall-tree query and interval
        # integration instead of repeatedly sending the full active population.
        target_reaction_area_um2 = np.zeros(count, dtype=np.float64)
        reaction_lanes = np.flatnonzero(active_mask & (reaction_radius_um > 0.0))
        if reaction_lanes.size:
            target_reaction_area_um2[reaction_lanes] = np.asarray(
                context.molecular_target_field.reaction_area_um2(
                    centres_xz_um[reaction_lanes],
                    tangents_xz[reaction_lanes],
                    reaction_radius_um[reaction_lanes],
                ),
                dtype=np.float64,
            )

        # Compare target-positive area with the entire possible reaction disk.
        # The result ranges from no target overlap to full target overlap and is
        # useful both for bond formation and for saved diagnostic output.
        full_reaction_area_um2 = np.pi * reaction_radius_um * reaction_radius_um
        target_overlap_fraction = np.divide(
            target_reaction_area_um2,
            full_reaction_area_um2,
            out=np.zeros(count, dtype=np.float64),
            where=full_reaction_area_um2 > 0.0,
        )

        # First estimate how fast the bubble surface slides across the wall using
        # background motion only. Translation and rotation both contribute to this
        # surface slip, just as a rolling wheel's rim speed depends on both.
        preliminary_tangential_velocity = np.sum(particle_velocity * tangents_xz, axis=1)
        preliminary_slip = surface_slip_velocity_um_s(
            preliminary_tangential_velocity,
            radii,
            angular_velocity,
        )

        # Use current bonds, wall gap, available target area, and estimated slip to
        # calculate bond formation, breakage, force, torque, and state-change rates.
        # "Mean field" means many individual bonds are represented by smooth
        # expected totals instead of tracking every molecule one by one.
        molecular_evaluation = evaluate_mean_field_bonds(
            molecular_bond_count,
            molecular_bond_extension,
            molecular_gap_um,
            radii,
            preliminary_slip,
            target_reaction_area_um2,
            context.molecular_binding_parameters,
            use_numba=context.use_numba,
        )

        # Bond forces are first reported along and away from the wall. Rebuild their
        # global X-Z vector so it can be combined with collision forces consistently.
        bond_force_xz_pn = (
            molecular_evaluation.force_t_pn[:, None] * tangents_xz
            + molecular_evaluation.force_n_pn[:, None] * normalized_normals
        )

        # Remove any placeholder-lane force so an inactive row can never look like
        # a physical bond measurement in saved results.
        bond_force_xz_pn[~active_mask] = 0.0

    # A force is a push or pull, not a velocity. Mobility tells us how much motion
    # that push or pull creates for this bubble in this liquid and near this wall.
    # Use the same local mobility for collision and bond loads so wall resistance
    # affects every source of bubble motion consistently.
    # A repulsive pair force does not directly prescribe a velocity.  Pass it
    # through the *same local wall mobility* as the background response, because
    # a nearby wall also resists collision-driven translation and can couple it
    # to rotation.
    maximum_collision_speed = 0.0
    if context.use_numba:
        collision_speed_by_lane = apply_particle_loads_batch(
            particle_velocity,
            angular_velocity,
            viscosity,
            radii,
            normalized_normals,
            mobility_entries,
            collision_force,
            active_mask,
            (
                molecular_evaluation.force_t_pn
                if molecular_evaluation is not None
                else None
            ),
            (
                molecular_evaluation.force_n_pn
                if molecular_evaluation is not None
                else None
            ),
            (
                molecular_evaluation.torque_y_pn_um
                if molecular_evaluation is not None
                else None
            ),
        )
        maximum_collision_speed = (
            float(np.max(collision_speed_by_lane))
            if collision_speed_by_lane.size
            else 0.0
        )
    for lane in range(0 if context.use_numba else count):
        # Empty lanes cannot receive forces or move.
        if not active_mask[lane]:
            continue

        # Rebuild the two easy wall directions for this bubble: the normal points
        # into the blood and the tangent runs along the wall.
        normal_x = float(normalized_normals[lane, 0])
        normal_z = float(normalized_normals[lane, 1])
        tangent_x = -normal_z
        tangent_z = normal_x

        # Read the collision force in global X-Z coordinates, then split it into
        # the along-wall part and the away-from-wall part. Near a wall, these two
        # directions face different amounts of liquid resistance.
        force_x = float(collision_force[lane, 0])
        force_z = float(collision_force[lane, 1])
        force_t = tangent_x * force_x + tangent_z * force_z
        force_n = normal_x * force_x + normal_z * force_z

        # Read the matching bond force and bond turning effect when binding is on.
        # If no molecular calculation ran, zeros mean that bonds add no load.
        bond_force_t = (
            float(molecular_evaluation.force_t_pn[lane])
            if molecular_evaluation is not None
            else 0.0
        )
        bond_force_n = (
            float(molecular_evaluation.force_n_pn[lane])
            if molecular_evaluation is not None
            else 0.0
        )
        bond_torque = (
            float(molecular_evaluation.torque_y_pn_um[lane])
            if molecular_evaluation is not None
            else 0.0
        )

        # Retrieve the five local mobility values prepared earlier. They describe
        # how along-wall force, normal force, and turning force affect this bubble.
        entries = mobility_entries[lane]

        # Combine collision and bond loads, then convert them into extra along-wall
        # speed, normal speed, and rotation speed. This is the central job of the
        # mobility model: changing a physical push or twist into motion.
        delta_t, delta_n, delta_omega = apply_scaled_mobility_scalar(
            float(viscosity[lane]),
            float(radii[lane]),
            force_t + bond_force_t,
            force_n + bond_force_n,
            bond_torque,
            float(entries[0]),
            float(entries[1]),
            float(entries[2]),
            float(entries[3]),
            float(entries[4]),
        )

        # Change the two local speed components back into global X-Z components.
        delta_x, delta_z = wall_to_global_components_scalar(
            delta_t,
            delta_n,
            normal_x,
            normal_z,
        )

        # Add the force-created movement to the background flow movement calculated
        # earlier. The result is the bubble's full instantaneous motion estimate.
        particle_velocity[lane, 0] += delta_x
        particle_velocity[lane, 1] += delta_z
        angular_velocity[lane] += delta_omega

        # Repeat the conversion with collision force alone only for a diagnostic.
        # This lets the caller see the largest speed caused specifically by
        # collision handling, without changing the already accepted total motion.
        collision_delta_t, collision_delta_n, _ = apply_scaled_mobility_scalar(
            float(viscosity[lane]),
            float(radii[lane]),
            force_t,
            force_n,
            0.0,
            float(entries[0]),
            float(entries[1]),
            float(entries[2]),
            float(entries[3]),
            float(entries[4]),
        )

        # Keep only the largest collision-only speed seen in this snapshot. A very
        # large value can warn that the time step may need closer checking.
        maximum_collision_speed = max(
            maximum_collision_speed,
            float(np.hypot(collision_delta_t, collision_delta_n)),
        )

    if molecular_evaluation is not None:
        # Bond force may have slowed translation or changed rotation, so the bubble
        # surface can now slide at a different speed than the first estimate.
        # Recalculate slip using the completed post-mobility motion.
        # The extension equation uses the actual post-mobility surface slip.
        # Re-evaluating the pure bond RHS changes dm/dt but not the current-state
        # load, and is therefore compatible with both Euler and Heun stages.
        final_tangential_velocity = np.sum(particle_velocity * tangents_xz, axis=1)
        final_slip = surface_slip_velocity_um_s(
            final_tangential_velocity,
            radii,
            angular_velocity,
        )

        # Re-evaluate only the bond-state change rates with the corrected slip.
        # The present bond load still belongs to the same frozen-time snapshot;
        # this second call makes future bond stretching consistent with final speed.
        molecular_evaluation = evaluate_mean_field_bonds(
            molecular_bond_count,
            molecular_bond_extension,
            molecular_gap_um,
            radii,
            final_slip,
            target_reaction_area_um2,
            context.molecular_binding_parameters,
            use_numba=context.use_numba,
        )

    red_blood_cell_evaluation = None
    if context.red_blood_cell_network is not None:
        _validate_red_blood_cell_vessel_owners(
            sampled_vessel_id,
            active_mask,
            ids,
            positions,
            float(time_s),
            context,
        )
        # RBC transport is closed against the accepted cycle-mean CFD field.
        # Cardiac carrier-flow modulation must not introduce a second spatial
        # propagation model into either drift or diffusion.
        red_blood_cell_gradients = (
            gradients
            if context.cardiac is None
            else sample_bilinear(context.velocity_gradient_s_inv, positions)
        )
        red_blood_cell_evaluation = evaluate_red_blood_cell_transport(
            sampled_vessel_id,
            red_blood_cell_gradients,
            wall_gap,
            normalized_normals,
            nearest_wall_unique & wall_normal_valid,
            radii,
            active_mask,
            context.red_blood_cell_network,
            use_numba=context.use_numba,
        )
        particle_velocity += red_blood_cell_evaluation.velocity_xz_um_s

    # Inactive lanes are retained to keep a fixed frame width.  Zeroing their
    # sampled outputs prevents stale or extrapolated values from being mistaken
    # for physical measurements.
    if not all_active:
        particle_velocity[~active_mask] = 0.0
        fluid_velocity[~active_mask] = 0.0
        angular_velocity[~active_mask] = 0.0
        sampled_shear[~active_mask] = 0.0
        cardiac_multiplier[~active_mask] = 0.0

    # These availability arrays exist only when molecular binding was evaluated.
    available_ligand_count = None
    available_target_count = None
    if molecular_evaluation is not None:
        # A ligand belongs to the bubble and a target belongs to the wall. Subtract
        # bonds already in use to find how many partners remain available for new
        # bonds. Clamp at zero because an available count cannot be negative.
        available_ligand_count = np.maximum(
            molecular_evaluation.ligand_count - molecular_evaluation.expected_bond_count,
            0.0,
        )
        available_target_count = np.maximum(
            molecular_evaluation.target_count - molecular_evaluation.expected_bond_count,
            0.0,
        )

    # Package everything into one named result object. The time integrator mainly
    # uses translation and rotation speeds, while the remaining fields explain why
    # that motion occurred and help detect questionable numerical or physical cases.
    generalized_mobility = np.zeros((count, 3, 3), dtype=np.float64)
    # Revised v15 predicts rigid-wall contact before committing a move. A bubble may begin a
    # step away from the wall and still reach it before that physical interval
    # ends, so the contact solver must be able to read M for every live lane.
    # Building M only for lanes already touching the wall would leave
    # a first-contact event with a zero matrix and no physically meaningful
    # reaction direction.
    all_mobility_valid = bool(all_active and np.all(np.isfinite(wall_gap)))
    mobility_lanes = (
        None
        if all_mobility_valid
        else np.flatnonzero(active_mask & np.isfinite(wall_gap))
    )
    if all_mobility_valid:
        mobility_builder = (
            build_generalized_mobility_batch
            if context.use_numba
            else generalized_mobility_matrix_xz
        )
        generalized_mobility = mobility_builder(
            viscosity, radii, normalized_normals, mobility_entries
        )
    elif mobility_lanes.size:
        mobility_builder = (
            build_generalized_mobility_batch
            if context.use_numba
            else generalized_mobility_matrix_xz
        )
        generalized_mobility[mobility_lanes] = mobility_builder(
            viscosity[mobility_lanes], radii[mobility_lanes],
            normalized_normals[mobility_lanes], mobility_entries[mobility_lanes],
        )

    return _RhsEvaluation(
        # These are the main answers used to move and rotate bubbles.
        particle_velocity_xz_um_s=particle_velocity,
        fluid_velocity_xz_um_s=fluid_velocity,
        angular_velocity_rad_s=angular_velocity,
        generalized_mobility=generalized_mobility,

        # These values describe bubble collisions and local wall geometry.
        collision_force_xz_pn=collision_force,
        collision_neighbor_count=collision_neighbors,
        wall_gap_um=wall_gap,
        hydrodynamic_gap_um=hydrodynamic_gap,
        wall_normal_xz=normalized_normals,
        gap_ratio=gap_ratio,
        near_wall_weight=wall_weight,
        two_wall_warning=two_wall_warning,

        # These values record the local flow sample seen by each bubble.
        sampled_wall_shear_stress_pa=sampled_shear,
        topological_vessel_id=sampled_vessel_id,
        cardiac_multiplier=cardiac_multiplier,

        # These summary numbers help judge whether collisions, wall assumptions,
        # or mobility calculations became unusually severe or inconsistent.
        maximum_physical_overlap_um=float(maximum_physical_overlap),
        maximum_collision_compression_um=float(maximum_collision_compression),
        interacting_pair_count=int(pair_count),
        collision_search_strategy=str(collision_search_strategy),
        maximum_reciprocity_relative_error=float(maximum_reciprocity_error),
        degenerate_near_wall_normal_count=int(degenerate_near_wall_normals),
        maximum_collision_speed_um_s=float(maximum_collision_speed),

        # These fields describe the current expected bond population and stretch.
        # They are None when molecular binding is disabled.
        bond_count_expected=(
            molecular_evaluation.expected_bond_count if molecular_evaluation is not None else None
        ),
        bond_total_tangential_extension_um=(
            molecular_evaluation.total_tangential_extension_um
            if molecular_evaluation is not None
            else None
        ),
        bond_mean_tangential_extension_um=(
            molecular_evaluation.mean_tangential_extension_um
            if molecular_evaluation is not None
            else None
        ),

        # These fields report the total force and turning effect made by bonds.
        bond_force_xz_pn=bond_force_xz_pn,
        bond_force_tangent_pn=(
            molecular_evaluation.force_t_pn if molecular_evaluation is not None else None
        ),
        bond_force_normal_pn=(
            molecular_evaluation.force_n_pn if molecular_evaluation is not None else None
        ),
        bond_torque_pn_um=(
            molecular_evaluation.torque_y_pn_um if molecular_evaluation is not None else None
        ),

        # These fields explain the stress on one bond and how quickly bonds are
        # being created or broken at this instant.
        single_bond_tension_pn=(
            molecular_evaluation.single_bond_tension_pn
            if molecular_evaluation is not None
            else None
        ),
        bond_formation_rate_bonds_s=(
            molecular_evaluation.formation_rate_bonds_s
            if molecular_evaluation is not None
            else None
        ),
        bond_dissociation_rate_s_inv=(
            molecular_evaluation.dissociation_rate_s
            if molecular_evaluation is not None
            else None
        ),

        # These fields describe how much target wall is reachable and how many
        # unused binding partners remain on the bubble and wall.
        target_reaction_area_um2=target_reaction_area_um2,
        available_ligand_count=available_ligand_count,
        available_target_count=available_target_count,
        target_overlap_fraction=target_overlap_fraction,

        # A saturated flag means a safety limit in the molecular formula was
        # reached. Saving it makes such cases visible during later review.
        molecular_bell_rate_saturated=(
            molecular_evaluation.bell_rate_saturated
            if molecular_evaluation is not None
            else None
        ),
        molecular_formation_rate_saturated=(
            molecular_evaluation.formation_rate_saturated
            if molecular_evaluation is not None
            else None
        ),

        # Keep the full molecular result as well, so later code can use details
        # without repeating the expensive bond calculation.
        molecular_binding_evaluation=molecular_evaluation,
        wall_normal_valid=wall_normal_valid,
        red_blood_cell_velocity_xz_um_s=(
            red_blood_cell_evaluation.velocity_xz_um_s
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_drift_velocity_xz_um_s=(
            red_blood_cell_evaluation.drift_velocity_xz_um_s
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_fick_velocity_xz_um_s=(
            red_blood_cell_evaluation.fick_velocity_xz_um_s
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_transverse_direction_xz=(
            red_blood_cell_evaluation.transverse_direction_xz
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_local_vessel_diameter_um=(
            red_blood_cell_evaluation.local_vessel_diameter_um
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_discharge_hematocrit=(
            red_blood_cell_evaluation.discharge_hematocrit
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_tube_hematocrit=(
            red_blood_cell_evaluation.tube_hematocrit
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_shear_rate_s_inv=(
            red_blood_cell_evaluation.shear_rate_s_inv
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_cfl_width_um=(
            red_blood_cell_evaluation.cfl_width_um
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_target_gap_um=(
            red_blood_cell_evaluation.target_gap_um
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_center_gap_um=(
            red_blood_cell_evaluation.center_gap_um
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_core_transverse_diffusivity_um2_s=(
            red_blood_cell_evaluation.core_transverse_diffusivity_um2_s
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_transverse_diffusivity_um2_s=(
            red_blood_cell_evaluation.transverse_diffusivity_um2_s
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_diffusion_taper_coordinate=(
            red_blood_cell_evaluation.diffusion_taper_coordinate
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_diffusion_taper=(
            red_blood_cell_evaluation.diffusion_taper
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_margination_length_um=(
            red_blood_cell_evaluation.margination_length_um
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_margination_time_s=(
            red_blood_cell_evaluation.margination_time_s
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_scale_activation=(
            red_blood_cell_evaluation.scale_activation
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_nearest_wall_unique=(
            red_blood_cell_evaluation.nearest_wall_unique
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_hematocrit_in_quantitative_range=(
            red_blood_cell_evaluation.hematocrit_in_quantitative_range
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_shear_rate_in_quantitative_range=(
            red_blood_cell_evaluation.shear_rate_in_quantitative_range
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_quantitative_applicability=(
            red_blood_cell_evaluation.quantitative_applicability
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_transverse_space_valid=(
            red_blood_cell_evaluation.transverse_space_valid
            if red_blood_cell_evaluation is not None
            else None
        ),
        red_blood_cell_diffusion_enabled=(
            red_blood_cell_evaluation.diffusion_enabled
            if red_blood_cell_evaluation is not None
            else None
        ),
    )


def _validate_red_blood_cell_vessel_owners(
    sampled_trajectory_vessel_id: np.ndarray,
    active: np.ndarray,
    permanent_microbubble_ids: np.ndarray,
    positions_grid: np.ndarray,
    physical_time_s: float,
    context: _EvaluationContext,
) -> None:
    """Raise a complete diagnostic before the RBC kernel sees an invalid owner."""

    owner = np.asarray(sampled_trajectory_vessel_id, dtype=np.int64) - 1
    owner_count = int(
        context.red_blood_cell_network.dense_diameter_um_by_vessel_id.size
    )
    invalid_lanes = np.flatnonzero(
        np.asarray(active, dtype=bool) & ((owner < 0) | (owner >= owner_count))
    )
    if invalid_lanes.size == 0:
        return

    positions = np.asarray(positions_grid, dtype=np.float64)
    origin = np.asarray(context.origin_xz_um, dtype=np.float64)
    positions_xz_um = origin[None, :] + positions * float(context.spacing_um)
    failures = tuple(
        ParticleVesselOwnershipFailure(
            permanent_microbubble_id=int(permanent_microbubble_ids[lane]),
            physical_time_s=float(physical_time_s),
            local_lane=int(lane),
            position_grid=(float(positions[lane, 0]), float(positions[lane, 1])),
            position_xz_um=(
                float(positions_xz_um[lane, 0]),
                float(positions_xz_um[lane, 1]),
            ),
            persistent_topological_vessel_id=int(
                sampled_trajectory_vessel_id[lane]
            ),
        )
        for lane in invalid_lanes
    )
    raise ParticleVesselOwnershipError(failures)


def _format_vessel_ownership_failures(
    failures: tuple[ParticleVesselOwnershipFailure, ...],
) -> str:
    """Format stable key/value diagnostics for the command-line error log."""

    lines = [
        "An active particle has no valid persistent topological vessel owner.",
        "Revised-v20 vessel-owner failure diagnostics:",
    ]
    for record_index, failure in enumerate(failures):
        lines.extend(
            (
                f"failure[{record_index}]:",
                "  permanent_microbubble_id="
                f"{failure.permanent_microbubble_id}",
                f"  physical_time_s={failure.physical_time_s:.17g}",
                f"  local_lane={failure.local_lane}",
                "  position_grid=["
                f"{failure.position_grid[0]:.17g}, "
                f"{failure.position_grid[1]:.17g}]",
                "  position_xz_um=["
                f"{failure.position_xz_um[0]:.17g}, "
                f"{failure.position_xz_um[1]:.17g}]",
                "  persistent_topological_vessel_id="
                f"{failure.persistent_topological_vessel_id}",
            )
        )
    return "\n".join(lines)


def _format_wall_gap_failures(
    failures: tuple[ParticleWallGapFailure, ...],
) -> str:
    """Format complete diagnostics for a true accepted-state penetration."""

    lines = [
        "An active particle has a materially negative authoritative wall gap.",
        "Particle wall-gap invariant failure diagnostics:",
    ]
    for record_index, failure in enumerate(failures):
        lines.extend(
            (
                f"failure[{record_index}]:",
                "  permanent_microbubble_id="
                f"{failure.permanent_microbubble_id}",
                f"  physical_time_s={failure.physical_time_s:.17g}",
                f"  local_lane={failure.local_lane}",
                "  position_grid=["
                f"{failure.position_grid[0]:.17g}, "
                f"{failure.position_grid[1]:.17g}]",
                "  position_xz_um=["
                f"{failure.position_xz_um[0]:.17g}, "
                f"{failure.position_xz_um[1]:.17g}]",
                f"  bubble_radius_um={failure.bubble_radius_um:.17g}",
                f"  wall_gap_um={failure.wall_gap_um:.17g}",
                "  geometry_roundoff_tolerance_um="
                f"{failure.geometry_roundoff_tolerance_um:.17g}",
            )
        )
    return "\n".join(lines)


def _nearest_indices(
    positions_grid: np.ndarray,
    shape: tuple[int, int],
) -> tuple[np.ndarray, np.ndarray]:
    """Return clipped nearest-cell indices for categorical grid data."""

    positions = np.asarray(positions_grid, dtype=np.float64)
    ii = np.clip(np.floor(positions[:, 0] + 0.5).astype(np.int64), 0, shape[0] - 1)
    jj = np.clip(np.floor(positions[:, 1] + 0.5).astype(np.int64), 0, shape[1] - 1)
    return ii, jj


def _smooth_wall_weight(gap_ratio: float, xi_near: float, xi_far: float) -> float:
    """Blend wall and free-space models with a smooth zero-slope transition.

    The cubic smoothstep avoids a sudden velocity or derivative jump when a
    bubble crosses the near/far thresholds.  A value of one means fully
    wall-corrected behaviour; zero means free-space behaviour.
    """

    scaled = (float(gap_ratio) - xi_near) / (xi_far - xi_near)
    scaled = min(1.0, max(0.0, scaled))
    return 1.0 - 3.0 * scaled * scaled + 2.0 * scaled * scaled * scaled
