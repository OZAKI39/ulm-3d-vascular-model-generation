"""Shared data containers for field-based microbubble transport."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from ulm_vascular_model_generator.utils.core.models import Vessel


@dataclass(frozen=True)
class PhysicsInput:
    """Resolved input files and loaded vessel objects."""

    swc_path: Path
    vessel_data_path: Path
    vessels: list[Vessel]
    vessel_metadata: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class GridDomain:
    """Two-dimensional X-Z grid that represents the planar vascular model."""

    origin_um: np.ndarray
    spacing_um: float
    shape: tuple[int, int]
    fixed_y_um: float
    x_coordinates_um: np.ndarray
    z_coordinates_um: np.ndarray


@dataclass(frozen=True)
class RasterizedVessels:
    """Grid fields obtained by converting cylindrical vessels to pixels."""

    lumen_mask: np.ndarray
    wall_mask: np.ndarray
    vessel_id: np.ndarray
    radius_um: np.ndarray
    flow_rate_um3_s: np.ndarray
    q2d_flow_um2_s: np.ndarray
    viscosity_mpas: np.ndarray
    direction_xz: np.ndarray
    distance_to_centerline_um: np.ndarray
    distance_to_wall_um: np.ndarray
    wall_normal_xz: np.ndarray
    lumen_fraction: np.ndarray
    junction_core_mask: np.ndarray | None = None


@dataclass(frozen=True)
class FiniteElementVelocityField:
    """Cell-local polynomial representation of the solved DOLFINx velocity."""

    degree: int
    cell_vertices_xz_um: np.ndarray
    polynomial_exponents: np.ndarray
    velocity_coefficients_um_s: np.ndarray
    bin_origin_xz_um: np.ndarray
    bin_size_um: float
    bin_shape: tuple[int, int]
    bin_offsets: np.ndarray
    bin_cell_indices: np.ndarray


@dataclass(frozen=True)
class HybridVelocityField:
    """Near-wall finite-element field and its two switching distances."""

    finite_element: FiniteElementVelocityField
    finite_element_distance_um: float
    regular_grid_distance_um: float


@dataclass(frozen=True)
class FlowField:
    """Hybrid velocity field plus Cartesian diagnostic and far-field caches."""

    velocity_xz_um_s: np.ndarray
    speed_um_s: np.ndarray
    wall_shear_stress_pa: np.ndarray
    hybrid_velocity: HybridVelocityField
    local_shear_stress_pa: np.ndarray | None = None
    initial_velocity_xz_um_s: np.ndarray | None = None
    initial_speed_um_s: np.ndarray | None = None
    divergence_s_inv: np.ndarray | None = None
    wall_penetration_um_s: np.ndarray | None = None
    pressure: np.ndarray | None = None
    inlet_label: np.ndarray | None = None
    outlet_label: np.ndarray | None = None
    boundary_velocity_xz_um_s: np.ndarray | None = None
    boundary_normal_xz: np.ndarray | None = None
    boundary_weight: np.ndarray | None = None
    boundary_edge_length_um: np.ndarray | None = None
    open_boundary_flux_um2_s: np.ndarray | None = None
    face_flux_x_um2_s: np.ndarray | None = None
    face_flux_z_um2_s: np.ndarray | None = None
    inlet_target_by_label_um2_s: np.ndarray | None = None
    outlet_target_by_label_um2_s: np.ndarray | None = None
    inlet_actual_by_label_um2_s: np.ndarray | None = None
    outlet_actual_by_label_um2_s: np.ndarray | None = None
    open_face_cell_ij: np.ndarray | None = None
    open_face_index_ij: np.ndarray | None = None
    open_face_axis: np.ndarray | None = None
    open_face_normal_xz: np.ndarray | None = None
    open_face_center_xz_um: np.ndarray | None = None
    open_face_length_um: np.ndarray | None = None
    open_face_label: np.ndarray | None = None
    open_face_kind: np.ndarray | None = None
    open_section_point_xz_um: np.ndarray | None = None
    open_section_outward_normal_xz: np.ndarray | None = None
    open_section_tangent_xz: np.ndarray | None = None
    open_section_half_width_um: np.ndarray | None = None
    open_section_label: np.ndarray | None = None
    open_section_kind: np.ndarray | None = None
    solver_metadata: dict[str, object] = field(default_factory=dict)


@dataclass(frozen=True)
class ParticleTrajectories:
    """Append-only per-frame observations and permanent bubble lifecycles."""

    frame_offsets: np.ndarray
    bubble_id: np.ndarray
    positions_um: np.ndarray
    velocities_um_s: np.ndarray
    wall_shear_stress_pa: np.ndarray
    vessel_id: np.ndarray
    active: np.ndarray
    diameter_um: np.ndarray
    wall_gap_um: np.ndarray
    wall_contact: np.ndarray
    wall_normal_xz: np.ndarray
    registry_bubble_id: np.ndarray
    registry_diameter_um: np.ndarray
    birth_frame: np.ndarray
    death_frame: np.ndarray
    termination_reason: np.ndarray
    active_count_per_frame: np.ndarray
    injected_count_per_frame: np.ndarray
    terminated_count_per_frame: np.ndarray
    metadata: dict[str, object]
    realized_velocities_um_s: np.ndarray | None = None
    fluid_velocities_um_s: np.ndarray | None = None
    angular_velocity_rad_s: np.ndarray | None = None
    rotation_angle_rad: np.ndarray | None = None
    collision_force_xz_pn: np.ndarray | None = None
    collision_neighbor_count: np.ndarray | None = None
    gap_ratio: np.ndarray | None = None
    near_wall_weight: np.ndarray | None = None
    two_wall_warning: np.ndarray | None = None
    contact_constraint_active: np.ndarray | None = None
    contact_reaction_force_pn: np.ndarray | None = None
    contact_free_normal_velocity_um_s: np.ndarray | None = None
    contact_constrained_normal_velocity_um_s: np.ndarray | None = None
    registry_scheduled_injection_time_s: np.ndarray | None = None
    registry_admission_time_s: np.ndarray | None = None
    registry_exit_time_s: np.ndarray | None = None
    registry_inlet_wait_time_s: np.ndarray | None = None
    cardiac_multiplier: np.ndarray | None = None
    cardiac_waveform_time_s: np.ndarray | None = None
    cardiac_waveform_multiplier: np.ndarray | None = None
    cardiac_path_distance_um: np.ndarray | None = None
    cardiac_delay_s: np.ndarray | None = None
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
    registry_final_bond_count_expected: np.ndarray | None = None
    registry_final_bond_total_tangential_extension_um: np.ndarray | None = None
    registry_target_exposure_time_s: np.ndarray | None = None
    registry_target_exposure_event_count: np.ndarray | None = None
    registry_target_reaction_area_time_um2_s: np.ndarray | None = None
    registry_target_exposure_right_censored: np.ndarray | None = None
    registry_target_exposure_quantitative_applicability_fraction: np.ndarray | None = None
    red_blood_cell_velocity_xz_um_s: np.ndarray | None = None
    red_blood_cell_drift_velocity_xz_um_s: np.ndarray | None = None
    red_blood_cell_fick_velocity_xz_um_s: np.ndarray | None = None
    red_blood_cell_local_vessel_diameter_um: np.ndarray | None = None
    red_blood_cell_discharge_hematocrit: np.ndarray | None = None
    red_blood_cell_tube_hematocrit: np.ndarray | None = None
    red_blood_cell_shear_rate_s_inv: np.ndarray | None = None
    red_blood_cell_cfl_width_um: np.ndarray | None = None
    red_blood_cell_target_gap_um: np.ndarray | None = None
    red_blood_cell_transverse_diffusivity_um2_s: np.ndarray | None = None
    red_blood_cell_margination_length_um: np.ndarray | None = None
    red_blood_cell_margination_time_s: np.ndarray | None = None
    red_blood_cell_scale_activation: np.ndarray | None = None
    red_blood_cell_nearest_wall_unique: np.ndarray | None = None
    red_blood_cell_hematocrit_in_quantitative_range: np.ndarray | None = None
    red_blood_cell_shear_rate_in_quantitative_range: np.ndarray | None = None
    red_blood_cell_quantitative_applicability: np.ndarray | None = None
    red_blood_cell_transverse_space_valid: np.ndarray | None = None
    registry_final_vessel_id: np.ndarray | None = None
    topological_commitment_parent_vessel_id: np.ndarray | None = None
    topological_commitment_child_vessel_id: np.ndarray | None = None
    topological_commitment_point_xz_um: np.ndarray | None = None
    topological_commitment_downstream_normal_xz: np.ndarray | None = None
    topological_commitment_tangent_xz: np.ndarray | None = None
    topological_commitment_half_width_um: np.ndarray | None = None
    topological_commitment_transition_end_distance_um: np.ndarray | None = None
    topological_commitment_distance_um: np.ndarray | None = None
    topological_event_bubble_id: np.ndarray | None = None
    topological_event_time_s: np.ndarray | None = None
    topological_event_from_vessel_id: np.ndarray | None = None
    topological_event_to_vessel_id: np.ndarray | None = None
    topological_event_section_index: np.ndarray | None = None
    topological_event_position_xz_um: np.ndarray | None = None
