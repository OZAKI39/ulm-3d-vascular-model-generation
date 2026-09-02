"""
Boundary-fitted Gmsh / DOLFINx Stokes backend.

The authoritative lumen remains ``ContinuousVesselGeometry``.  Gmsh meshes its
piecewise-linear exterior ring and preserves separate physical tags for solid
wall, root inlets, and terminal outlets.  The exported vascular flow is imposed
as a parabolic velocity profile at every inlet and terminal outlet, while one
pressure degree of freedom is fixed only to remove the constant pressure
nullspace. DOLFINx solves a Taylor--Hood mixed velocity-pressure problem on
that conforming mesh. The solved velocity is
exported both as cell-local finite-element polynomials and as a Cartesian
far-field cache. Particles select between them from their exact continuous-wall
distance.

Gmsh, DOLFINx, PETSc, and MPI are imported lazily so importing the rest of the
trajectory package does not initialize the finite-element stack.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable

import numpy as np
from ulm_vascular_model_generator.utils.core.models import Vessel

from ..core.config import FieldConfig
from ..core.types import (
    FlowField,
    GridDomain,
    HybridVelocityField,
    RasterizedVessels,
)
from .hybrid_velocity import (
    build_finite_element_velocity_field,
    interpolation_nodes,
)
from .face_flux_projection import (
    cell_velocity_to_face_flux,
    divergence_from_face_flux,
    normalized_divergence_error_from_flux,
    solid_wall_penetration_field,
)
from .flow_boundaries import (
    apply_open_boundary_fluxes,
    build_flux_boundaries,
    realized_open_boundary_flux_field,
)
from .initial_velocity import build_initial_velocity_field

if TYPE_CHECKING:
    from ..geometry.continuous_vessel_geometry import ContinuousVesselGeometry


FLUID_PHYSICAL_TAG = 1
WALL_PHYSICAL_TAG = 10
INLET_PHYSICAL_TAG_OFFSET = 1_000
OUTLET_PHYSICAL_TAG_OFFSET = 2_000
MMHG_TO_PA = 133.322387415
SQUARE_METRES_TO_SQUARE_MICROMETRES = 1.0e12


@dataclass(frozen=True)
class BoundaryMeshBlueprint:
    """Pure-NumPy description used to construct and test the Gmsh model."""

    points_xz_um: np.ndarray
    segment_physical_tag: np.ndarray
    physical_tag_to_section_index: dict[int, int]

    @property
    def segment_count(self) -> int:
        return int(self.points_xz_um.shape[0])


@dataclass(frozen=True)
class _DolfinxMeshBundle:
    mesh: object
    cell_tags: object
    facet_tags: object
    blueprint: BoundaryMeshBlueprint
    gmsh_version: str
    bulk_size_um: float
    wall_size_um: float
    refinement_distance_um: float


def section_physical_tag(kind: int, label: int) -> int:
    """Return a stable Gmsh facet tag for one inlet or outlet section."""

    if int(label) < 1:
        raise ValueError("Open-section labels must be positive.")
    if int(kind) < 0:
        return INLET_PHYSICAL_TAG_OFFSET + int(label)
    if int(kind) > 0:
        return OUTLET_PHYSICAL_TAG_OFFSET + int(label)
    raise ValueError("Open-section kind must be -1 for inlet or +1 for outlet.")


def build_boundary_mesh_blueprint(
    continuous_geometry: "ContinuousVesselGeometry",
) -> BoundaryMeshBlueprint:
    """
    Classify every continuous boundary segment as wall, inlet, or outlet.

    ``ContinuousVesselGeometry`` already records the indices of all solid wall
    rows.  The complementary rows are flat anatomical openings and are matched
    to their exact section by plane and transverse-coordinate distance.
    """

    starts = np.asarray(continuous_geometry.full_boundary_start_xz_um, dtype=np.float64)
    ends = np.asarray(continuous_geometry.full_boundary_end_xz_um, dtype=np.float64)
    if starts.ndim != 2 or starts.shape[1] != 2 or starts.shape != ends.shape:
        raise ValueError("Continuous boundary rows must have shape (N, 2).")
    if starts.shape[0] < 3:
        raise ValueError("A Gmsh lumen boundary requires at least three segments.")

    scale = max(float(np.max(np.abs(starts), initial=0.0)), 1.0)
    continuity_tolerance = 2_048.0 * np.finfo(np.float64).eps * scale
    next_starts = np.roll(starts, -1, axis=0)
    if not np.allclose(ends, next_starts, rtol=0.0, atol=continuity_tolerance):
        raise ValueError(
            "Continuous boundary segments do not form one ordered closed ring."
        )

    segment_count = int(starts.shape[0])
    tags = np.full(segment_count, WALL_PHYSICAL_TAG, dtype=np.int32)
    solid_ring_indices = np.asarray(
        continuous_geometry.solid_face_ring_index, dtype=np.int64
    )
    if np.any((solid_ring_indices < 0) | (solid_ring_indices >= segment_count)):
        raise ValueError("Solid-wall ring indices fall outside the full boundary.")
    is_solid = np.zeros(segment_count, dtype=bool)
    is_solid[solid_ring_indices] = True
    open_indices = np.flatnonzero(~is_solid)
    if open_indices.size == 0:
        raise ValueError("The boundary-fitted lumen contains no open inlet or outlet.")

    section_points = np.asarray(
        continuous_geometry.open_section_point_xz_um, dtype=np.float64
    )
    section_normals = np.asarray(
        continuous_geometry.open_section_outward_normal_xz, dtype=np.float64
    )
    section_tangents = np.asarray(
        continuous_geometry.open_section_tangent_xz, dtype=np.float64
    )
    section_half_widths = np.asarray(
        continuous_geometry.open_section_half_width_um, dtype=np.float64
    )
    section_labels = np.asarray(continuous_geometry.open_section_label, dtype=np.int32)
    section_kinds = np.asarray(continuous_geometry.open_section_kind, dtype=np.int8)
    if (
        section_points.shape != section_normals.shape
        or section_points.shape != section_tangents.shape
    ):
        raise ValueError(
            "Continuous open-section point, normal, and tangent arrays disagree."
        )
    if section_points.shape[0] == 0:
        raise ValueError(
            "The boundary-fitted lumen requires at least one open section."
        )

    tag_to_section: dict[int, int] = {}
    for section_index, (kind, label) in enumerate(
        zip(section_kinds, section_labels, strict=True)
    ):
        tag = section_physical_tag(int(kind), int(label))
        if tag in tag_to_section:
            raise ValueError(f"Duplicate open-section physical tag {tag}.")
        tag_to_section[tag] = int(section_index)

    midpoints = 0.5 * (starts[open_indices] + ends[open_indices])
    segment_lengths = np.linalg.norm(ends[open_indices] - starts[open_indices], axis=1)
    for row, (ring_index, midpoint) in enumerate(
        zip(open_indices, midpoints, strict=True)
    ):
        relative = midpoint[None, :] - section_points
        plane_distance = np.abs(np.einsum("ij,ij->i", relative, section_normals))
        transverse = np.abs(np.einsum("ij,ij->i", relative, section_tangents))
        transverse_excess = np.maximum(transverse - section_half_widths, 0.0)
        score = plane_distance + transverse_excess
        section_index = int(np.argmin(score))
        tolerance = max(
            continuity_tolerance,
            1.0e-7,
            0.51 * float(segment_lengths[row]),
        )
        if (
            float(plane_distance[section_index]) > tolerance
            or float(transverse_excess[section_index]) > tolerance
        ):
            raise ValueError(
                "Could not associate continuous opening segment "
                f"{int(ring_index)} with an inlet or outlet section."
            )
        tags[ring_index] = section_physical_tag(
            int(section_kinds[section_index]),
            int(section_labels[section_index]),
        )

    observed_open_tags = set(
        int(value) for value in np.unique(tags) if int(value) != WALL_PHYSICAL_TAG
    )
    missing = set(tag_to_section) - observed_open_tags
    if missing:
        raise ValueError(f"Open sections {sorted(missing)} have no boundary segments.")
    return BoundaryMeshBlueprint(
        points_xz_um=np.ascontiguousarray(starts),
        segment_physical_tag=np.ascontiguousarray(tags),
        physical_tag_to_section_index=tag_to_section,
    )


def solve_dolfinx_stokes_gmsh_2d(
    domain: GridDomain,
    raster: RasterizedVessels,
    cfg: FieldConfig,
    vessels: list[Vessel] | tuple[Vessel, ...],
    continuous_geometry: "ContinuousVesselGeometry",
    *,
    vessel_metadata: dict[str, str] | None = None,
    progress_callback: Callable[[str], None] | None = None,
    subprogress_callback: Callable[[str, int, int], None] | None = None,
) -> FlowField:
    """Generate a boundary-fitted mesh, solve Stokes flow, and sample to the grid.

    ``progress_callback`` receives solver module names as those modules start.
    ``subprogress_callback`` receives the active operation and quantitative
    checkpoints within long Gmsh and Stokes modules. Both are optional so
    library callers and numerical tests do not need a progress UI.
    """

    if continuous_geometry is None:
        raise ValueError(
            "The Gmsh/DOLFINx backend requires ContinuousVesselGeometry; "
            "pass the same geometry used by the particle solver."
        )
    if not vessels:
        raise ValueError("The Gmsh/DOLFINx backend requires the exported Vessel list.")
    # Kept in the public signature so existing callers remain compatible. The
    # fixed-flow formulation deliberately does not consume pressure metadata.
    del vessel_metadata
    _report_progress(progress_callback, "fixed-flow boundaries")

    _report_progress(progress_callback, "DOLFINx stack")
    optional = _import_dolfinx_stack()
    _report_progress(progress_callback, "Gmsh mesh")
    mesh_bundle = _build_dolfinx_mesh(
        continuous_geometry=continuous_geometry,
        domain=domain,
        cfg=cfg,
        optional=optional,
        subprogress_callback=subprogress_callback,
    )
    _report_progress(progress_callback, "Stokes linear solve")
    solution = _solve_taylor_hood_stokes(
        mesh_bundle=mesh_bundle,
        continuous_geometry=continuous_geometry,
        vessels=vessels,
        cfg=cfg,
        optional=optional,
        subprogress_callback=subprogress_callback,
    )
    _report_progress(progress_callback, "finite-element cache")
    finite_element_velocity = _export_finite_element_velocity(
        solution["velocity"],
        mesh_bundle.mesh,
        degree=int(cfg.dolfinx_velocity_degree),
        preferred_bin_size_um=min(
            float(mesh_bundle.bulk_size_um),
            float(domain.spacing_um),
        ),
    )
    hybrid_velocity = HybridVelocityField(
        finite_element=finite_element_velocity,
        finite_element_distance_um=float(
            cfg.hybrid_finite_element_distance_um
        ),
        regular_grid_distance_um=float(
            cfg.hybrid_finite_element_distance_um
            + cfg.hybrid_transition_width_um
        ),
    )

    _report_progress(progress_callback, "Cartesian sampling")
    sample_points, sample_indices = _cartesian_lumen_points(domain, raster.lumen_mask)
    velocity_samples, located_velocity = _sample_dolfinx_function(
        solution["velocity"], mesh_bundle.mesh, sample_points, optional
    )
    pressure_samples, located_pressure = _sample_dolfinx_function(
        solution["pressure"], mesh_bundle.mesh, sample_points, optional
    )
    wall_state = continuous_geometry.exact_solid_wall_state_xz_um_accelerated(
        sample_points[:, :2]
    )
    wall_sample_offset_um = max(1.0e-7, 1.0e-5 * mesh_bundle.wall_size_um)
    wall_sample_points = np.zeros_like(sample_points)
    wall_sample_points[:, :2] = np.asarray(
        wall_state.nearest_point_xz_um, dtype=np.float64
    ).reshape(-1, 2) + wall_sample_offset_um * np.asarray(
        wall_state.inward_normal_xz, dtype=np.float64
    ).reshape(-1, 2)
    local_gradient_samples, located_local_gradient = _sample_dolfinx_function(
        solution["velocity_gradient"],
        mesh_bundle.mesh,
        sample_points,
        optional,
    )
    wall_gradient_samples, located_wall_gradient = _sample_dolfinx_function(
        solution["velocity_gradient"],
        mesh_bundle.mesh,
        wall_sample_points,
        optional,
    )
    located = (
        located_velocity
        & located_pressure
        & located_local_gradient
        & located_wall_gradient
    )
    missing_count = int(np.count_nonzero(~located))
    if missing_count:
        diagnostics, cause = _cartesian_sampling_failure_diagnostics(
            domain=domain,
            continuous_geometry=continuous_geometry,
            mesh=mesh_bundle.mesh,
            sample_points=sample_points,
            sample_indices=sample_indices,
            wall_state=wall_state,
            wall_sample_points=wall_sample_points,
            wall_sample_offset_um=wall_sample_offset_um,
            located_velocity=located_velocity,
            located_pressure=located_pressure,
            located_local_gradient=located_local_gradient,
            located_wall_gradient=located_wall_gradient,
        )
        print(diagnostics, flush=True)
        raise ValueError(
            f"DOLFINx Cartesian sampling failed for {missing_count} of "
            f"{sample_points.shape[0]} accepted lumen samples. {cause} "
            "Detailed point diagnostics were printed immediately above."
        )

    velocity = np.zeros((*domain.shape, 2), dtype=np.float64)
    pressure = np.zeros(domain.shape, dtype=np.float64)
    local_gradient = np.zeros((*domain.shape, 2, 2), dtype=np.float64)
    wall_gradient = np.zeros((*domain.shape, 2, 2), dtype=np.float64)
    velocity[sample_indices[:, 0], sample_indices[:, 1]] = np.asarray(
        velocity_samples, dtype=np.float64
    ).reshape(-1, 2)
    pressure_kinematic_to_mmHg = (
        float(cfg.blood_density_kg_m3)
        / SQUARE_METRES_TO_SQUARE_MICROMETRES
        / MMHG_TO_PA
    )
    pressure[sample_indices[:, 0], sample_indices[:, 1]] = (
        np.asarray(pressure_samples, dtype=np.float64).reshape(-1)
        * pressure_kinematic_to_mmHg
    )
    local_gradient[sample_indices[:, 0], sample_indices[:, 1]] = np.asarray(
        local_gradient_samples, dtype=np.float64
    ).reshape(-1, 2, 2)
    wall_gradient[sample_indices[:, 0], sample_indices[:, 1]] = np.asarray(
        wall_gradient_samples, dtype=np.float64
    ).reshape(-1, 2, 2)
    velocity[~raster.lumen_mask] = 0.0
    pressure[~raster.lumen_mask] = 0.0

    _report_progress(progress_callback, "wall shear")
    speed = np.linalg.norm(velocity, axis=-1)
    wall_shear = _continuous_wall_shear_on_grid(
        gradient=wall_gradient,
        raster=raster,
        domain=domain,
        continuous_geometry=continuous_geometry,
    )
    local_shear = _local_shear_stress_on_grid(
        gradient=local_gradient,
        raster=raster,
        domain=domain,
    )

    _report_progress(progress_callback, "flux diagnostics")
    boundaries = build_flux_boundaries(
        domain,
        raster,
        vessels,
        effective_thickness_um=cfg.effective_thickness_um,
        depth_cells=cfg.boundary_depth_cells,
    )
    inlet_actual = np.asarray(
        solution["inlet_actual_by_label_um2_s"], dtype=np.float64
    )
    outlet_actual = np.asarray(
        solution["outlet_actual_by_label_um2_s"], dtype=np.float64
    )
    face_flux_x, face_flux_z = cell_velocity_to_face_flux(
        velocity, raster.lumen_mask, domain.spacing_um
    )
    face_flux_x, face_flux_z = apply_open_boundary_fluxes(
        face_flux_x,
        face_flux_z,
        boundaries,
    )
    divergence = divergence_from_face_flux(
        face_flux_x, face_flux_z, domain, raster.lumen_mask
    )
    wall_penetration = solid_wall_penetration_field(
        face_flux_x,
        face_flux_z,
        raster.lumen_mask,
        boundaries,
        domain.spacing_um,
    )
    normalized_grid_divergence = normalized_divergence_error_from_flux(
        velocity,
        divergence,
        domain,
        raster.lumen_mask,
    )

    _report_progress(progress_callback, "physical acceptance")
    initial_velocity = build_initial_velocity_field(
        raster, direction_smoothing_cells=1.0
    )
    initial_velocity[~raster.lumen_mask] = 0.0
    actual_inlet_total = float(np.sum(inlet_actual))
    actual_outlet_total = float(np.sum(outlet_actual))
    actual_flux_scale = max(
        abs(actual_inlet_total),
        abs(actual_outlet_total),
        np.finfo(np.float64).eps,
    )
    actual_flux_balance_error = (
        abs(actual_inlet_total - actual_outlet_total) / actual_flux_scale
    )
    inlet_target_error = abs(
        actual_inlet_total - float(solution["total_inlet_q2d_um2_s"])
    ) / max(
        abs(float(solution["total_inlet_q2d_um2_s"])),
        np.finfo(np.float64).eps,
    )
    outlet_target_error = abs(
        actual_outlet_total - float(solution["total_outlet_q2d_um2_s"])
    ) / max(
        abs(float(solution["total_outlet_q2d_um2_s"])),
        np.finfo(np.float64).eps,
    )
    inlet_max_target_error = _maximum_labeled_flux_relative_error(
        inlet_actual,
        np.asarray(solution["inlet_target_by_label_um2_s"], dtype=np.float64),
    )
    outlet_max_target_error = _maximum_labeled_flux_relative_error(
        outlet_actual,
        np.asarray(solution["outlet_target_by_label_um2_s"], dtype=np.float64),
    )
    physical_converged = bool(
        solution["linear_solver_converged"]
        and actual_flux_balance_error <= float(cfg.flux_tolerance)
        and inlet_target_error <= float(cfg.flux_tolerance)
        and outlet_target_error <= float(cfg.flux_tolerance)
        and inlet_max_target_error <= float(cfg.flux_tolerance)
        and outlet_max_target_error <= float(cfg.flux_tolerance)
    )

    metadata = {
        "solver_mode": "dolfinx_stokes_gmsh_2d",
        "vascular_geometry_mode": "planar_2d",
        "solver_dimension": 2,
        "flow_conversion": (
            "Q3D_over_inlet_circular_area_equivalent_planar_depth"
        ),
        "physical_model": "steady_incompressible_stokes_boundary_fitted_fem_2d",
        "finite_element_pair": (
            f"P{int(cfg.dolfinx_velocity_degree)}-"
            f"P{int(cfg.dolfinx_pressure_degree)}_taylor_hood"
        ),
        "velocity_gradient_projection": (
            f"discontinuous_P{max(0, int(cfg.dolfinx_velocity_degree) - 1)}"
        ),
        "geometry_source": "continuous_vessel_geometry_exterior_ring",
        "geometry_hash_sha256": str(continuous_geometry.geometry_hash_sha256),
        "gmsh_version": mesh_bundle.gmsh_version,
        "gmsh_element_order": int(cfg.gmsh_element_order),
        "gmsh_bulk_mesh_size_um": float(mesh_bundle.bulk_size_um),
        "gmsh_wall_mesh_size_um": float(mesh_bundle.wall_size_um),
        "gmsh_wall_refinement_distance_um": float(mesh_bundle.refinement_distance_um),
        "linear_solver_backend": str(solution["linear_solver_backend"]),
        "linear_solver_converged_reason": int(
            solution["linear_solver_converged_reason"]
        ),
        "linear_solver_iterations": int(solution["linear_solver_iterations"]),
        "linear_solver_relative_residual": float(
            solution["linear_solver_relative_residual"]
        ),
        "linear_solver_residual_tolerance": float(
            solution.get(
                "linear_solver_residual_tolerance",
                cfg.dolfinx_ksp_rtol,
            )
        ),
        "dolfinx_ksp_converged_reason": int(solution.get("ksp_converged_reason", 0)),
        "dolfinx_ksp_iterations": int(solution.get("ksp_iterations", 0)),
        "dolfinx_ksp_rtol": float(cfg.dolfinx_ksp_rtol),
        "physical_acceptance_schema": str(solution["physical_acceptance_schema"]),
        "kinematic_viscosity_um2_s": float(cfg.kinematic_viscosity_um2_s),
        "blood_density_kg_m3": float(cfg.blood_density_kg_m3),
        "effective_thickness_um": float(cfg.effective_thickness_um),
        "target_inlet_q2d_um2_s": float(solution["total_inlet_q2d_um2_s"]),
        "target_outlet_q2d_um2_s": float(solution["total_outlet_q2d_um2_s"]),
        "prescribed_flux_balance_relative_error": float(
            solution["flux_balance_relative_error"]
        ),
        "actual_inlet_q2d_um2_s": actual_inlet_total,
        "actual_outlet_q2d_um2_s": actual_outlet_total,
        "actual_flux_balance_relative_error": actual_flux_balance_error,
        "inlet_flux_relative_error": inlet_target_error,
        "outlet_flux_relative_error": outlet_target_error,
        "maximum_inlet_flux_relative_error": inlet_max_target_error,
        "maximum_outlet_flux_relative_error": outlet_max_target_error,
        "boundary_condition_type": "inlet_and_terminal_outlet_fixed_flow_dirichlet",
        "wall_boundary_condition": "no_slip_dirichlet",
        "inlet_boundary_condition": "parabolic_velocity_dirichlet_from_exported_flow",
        "outlet_boundary_condition": "parabolic_velocity_dirichlet_from_exported_flow",
        "vascular_flow_model": "fixed_total_inlet_equal_terminal_shares",
        "pressure_gauge": "single_pressure_degree_of_freedom_zero",
        "pressure_boundary_source": "none_flow_prescribed_at_all_openings",
        "fem_pressure_unit": "um2_s-2_kinematic_gauge",
        "sampled_pressure_unit": "mmHg",
        "sampled_pressure_semantics": "gauge_relative_to_pinned_pressure_dof",
        "wall_shear_definition": "fem_velocity_gradient_with_continuous_wall_normal",
        "wall_shear_gradient_sampling": "nearest_solid_wall_point_offset_inward",
        "wall_shear_sampling_offset_um": float(wall_sample_offset_um),
        "local_shear_definition": "mu*sqrt(2*D:D)_newtonian_viscous_stress_magnitude",
        "local_shear_gradient_sampling": "fem_velocity_gradient_at_cartesian_lumen_centres",
        "particle_velocity_sampling": (
            "continuous_wall_distance_hybrid_fem_transition_cartesian"
        ),
        "hybrid_finite_element_distance_um": float(
            hybrid_velocity.finite_element_distance_um
        ),
        "hybrid_regular_grid_distance_um": float(
            hybrid_velocity.regular_grid_distance_um
        ),
        "hybrid_transition_width_um": float(cfg.hybrid_transition_width_um),
        "finite_element_velocity_representation": (
            "cell_local_complete_reference_triangle_polynomial"
        ),
        "finite_element_velocity_cell_count": int(
            finite_element_velocity.cell_vertices_xz_um.shape[0]
        ),
        "sampling_stage": (
            "dolfinx_solution_to_hybrid_fem_and_cartesian_far_field_cache"
        ),
        "sampled_lumen_point_count": int(sample_indices.shape[0]),
        "sampled_grid_normalized_divergence_error": float(normalized_grid_divergence),
        "sampled_grid_divergence_is_acceptance_metric": False,
        "physical_converged": physical_converged,
        "converged": physical_converged,
        "mesh_cell_count": int(
            mesh_bundle.mesh.topology.index_map(
                mesh_bundle.mesh.topology.dim
            ).size_global
        ),
        "mesh_vertex_count": int(mesh_bundle.mesh.topology.index_map(0).size_global),
    }
    if not bool(solution["linear_solver_converged"]):
        raise ValueError(
            "DOLFINx linear solve did not converge: backend="
            f"{solution['linear_solver_backend']}, reason="
            f"{int(solution['linear_solver_converged_reason'])}, relative_residual="
            f"{float(solution['linear_solver_relative_residual']):.6g}."
        )
    if not physical_converged:
        raise ValueError(
            "DOLFINx fixed inlet/outlet-flow solution failed flux acceptance: "
            f"inlet_total_error={inlet_target_error:.6g}, "
            f"outlet_total_error={outlet_target_error:.6g}, "
            f"inlet_max_error={inlet_max_target_error:.6g}, "
            f"outlet_max_error={outlet_max_target_error:.6g}, "
            f"balance_error={actual_flux_balance_error:.6g}, "
            f"limit={float(cfg.flux_tolerance):.6g}."
        )
    return FlowField(
        velocity_xz_um_s=velocity.astype(np.float32),
        speed_um_s=speed.astype(np.float32),
        wall_shear_stress_pa=wall_shear.astype(np.float32),
        local_shear_stress_pa=local_shear.astype(np.float32),
        hybrid_velocity=hybrid_velocity,
        initial_velocity_xz_um_s=initial_velocity.astype(np.float32),
        initial_speed_um_s=np.linalg.norm(initial_velocity, axis=-1).astype(np.float32),
        divergence_s_inv=divergence.astype(np.float32),
        wall_penetration_um_s=wall_penetration.astype(np.float32),
        pressure=pressure.astype(np.float32),
        inlet_label=boundaries.inlet_label,
        outlet_label=boundaries.outlet_label,
        boundary_velocity_xz_um_s=boundaries.boundary_velocity_xz_um_s,
        boundary_normal_xz=boundaries.boundary_normal_xz,
        boundary_weight=boundaries.boundary_weight,
        boundary_edge_length_um=boundaries.boundary_edge_length_um,
        open_boundary_flux_um2_s=realized_open_boundary_flux_field(
            boundaries
        ),
        face_flux_x_um2_s=face_flux_x.astype(np.float32),
        face_flux_z_um2_s=face_flux_z.astype(np.float32),
        inlet_target_by_label_um2_s=boundaries.inlet_target_by_label_um2_s,
        outlet_target_by_label_um2_s=boundaries.outlet_target_by_label_um2_s,
        inlet_actual_by_label_um2_s=inlet_actual,
        outlet_actual_by_label_um2_s=outlet_actual,
        open_face_cell_ij=boundaries.open_face_cell_ij,
        open_face_index_ij=boundaries.open_face_index_ij,
        open_face_axis=boundaries.open_face_axis,
        open_face_normal_xz=boundaries.open_face_normal_xz,
        open_face_center_xz_um=boundaries.open_face_center_xz_um,
        open_face_length_um=boundaries.open_face_length_um,
        open_face_label=boundaries.open_face_label,
        open_face_kind=boundaries.open_face_kind,
        open_section_point_xz_um=boundaries.open_section_point_xz_um,
        open_section_outward_normal_xz=boundaries.open_section_outward_normal_xz,
        open_section_tangent_xz=boundaries.open_section_tangent_xz,
        open_section_half_width_um=boundaries.open_section_half_width_um,
        open_section_label=boundaries.open_section_label,
        open_section_kind=boundaries.open_section_kind,
        solver_metadata=metadata,
    )


def _report_progress(
    callback: Callable[[str], None] | None,
    completed_module: str,
) -> None:
    """Report a solver module without coupling the numerical backend to tqdm."""

    if callback is not None:
        callback(str(completed_module))


def _report_subprogress(
    callback: Callable[[str, int, int], None] | None,
    submodule: str,
    completed: int,
    total: int,
) -> None:
    """Report a checkpoint within one long-running solver module."""

    if callback is not None:
        callback(str(submodule), int(completed), int(total))


def _import_dolfinx_stack() -> dict[str, object]:
    import os

    if os.name == "nt":
        # Intel MPI may otherwise try to initialize an unavailable OFI network
        # provider when an IDE launches python.exe without `conda activate`.
        # This backend is serial on Windows, so shared-memory transport is enough.
        os.environ.setdefault("I_MPI_FABRICS", "shm")
    try:
        import basix.ufl
        import gmsh
        import ufl
        from dolfinx import fem, geometry
        from dolfinx.io import gmsh as gmshio
        from mpi4py import MPI
    except ImportError as exc:
        raise ImportError(
            "The boundary-fitted solver requires Gmsh and DOLFINx. "
            "Install the conda environment described by "
            "`ulm_microbubble_traj_gen/environment-dolfinx.yml`."
        ) from exc
    try:
        from dolfinx.fem.petsc import LinearProblem
        from petsc4py import PETSc
    except ImportError:
        LinearProblem = None
        PETSc = None
    try:
        from scipy.sparse import linalg as scipy_sparse_linalg
    except ImportError:
        scipy_sparse_linalg = None
    return {
        "basix_ufl": basix.ufl,
        "fem": fem,
        "geometry": geometry,
        "gmshio": gmshio,
        "LinearProblem": LinearProblem,
        "gmsh": gmsh,
        "MPI": MPI,
        "PETSc": PETSc,
        "scipy_sparse_linalg": scipy_sparse_linalg,
        "ufl": ufl,
    }


def _build_dolfinx_mesh(
    *,
    continuous_geometry: "ContinuousVesselGeometry",
    domain: GridDomain,
    cfg: FieldConfig,
    optional: dict[str, object],
    subprogress_callback: Callable[[str, int, int], None] | None = None,
) -> _DolfinxMeshBundle:
    progress_total = 8
    gmsh = optional["gmsh"]
    gmshio = optional["gmshio"]
    MPI = optional["MPI"]
    _report_subprogress(
        subprogress_callback,
        "build boundary blueprint",
        0,
        progress_total,
    )
    blueprint = build_boundary_mesh_blueprint(continuous_geometry)
    _report_subprogress(
        subprogress_callback,
        "resolve mesh-size controls",
        1,
        progress_total,
    )
    bulk_size = (
        float(cfg.gmsh_bulk_mesh_size_um)
        if float(cfg.gmsh_bulk_mesh_size_um) > 0.0
        else float(domain.spacing_um)
    )
    wall_size = (
        float(cfg.gmsh_wall_mesh_size_um)
        if float(cfg.gmsh_wall_mesh_size_um) > 0.0
        else 0.5 * bulk_size
    )
    refinement_distance = (
        float(cfg.gmsh_wall_refinement_distance_um)
        if float(cfg.gmsh_wall_refinement_distance_um) > 0.0
        else 3.0 * bulk_size
    )
    if wall_size > bulk_size:
        raise ValueError("field.gmsh_wall_mesh_size_um cannot exceed the bulk size.")

    _report_subprogress(
        subprogress_callback,
        "initialize Gmsh model",
        2,
        progress_total,
    )
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("continuous_vessel_lumen")
        _report_subprogress(
            subprogress_callback,
            "create boundary geometry",
            3,
            progress_total,
        )
        point_scale = max(
            float(np.max(np.abs(blueprint.points_xz_um), initial=0.0)), 1.0
        )
        point_merge_tolerance = (
            2_048.0 * np.finfo(np.float64).eps * point_scale
        )
        point_tags = []
        point_tag_by_bin: dict[tuple[int, int], list[tuple[np.ndarray, int]]] = {}
        for point in blueprint.points_xz_um:
            point_bin = np.rint(point / point_merge_tolerance).astype(np.int64)
            point_tag = None
            for offset_x in (-1, 0, 1):
                for offset_z in (-1, 0, 1):
                    nearby = point_tag_by_bin.get(
                        (
                            int(point_bin[0]) + offset_x,
                            int(point_bin[1]) + offset_z,
                        ),
                        (),
                    )
                    for existing_point, existing_tag in nearby:
                        if np.linalg.norm(point - existing_point) <= point_merge_tolerance:
                            point_tag = existing_tag
                            break
                    if point_tag is not None:
                        break
                if point_tag is not None:
                    break
            if point_tag is None:
                point_tag = gmsh.model.geo.addPoint(
                    float(point[0]), float(point[1]), 0.0, bulk_size
                )
                point_tag_by_bin.setdefault(
                    (int(point_bin[0]), int(point_bin[1])), []
                ).append((point.copy(), int(point_tag)))
            point_tags.append(int(point_tag))
        line_tags = [
            gmsh.model.geo.addLine(
                point_tags[index],
                point_tags[(index + 1) % blueprint.segment_count],
            )
            for index in range(blueprint.segment_count)
        ]
        loop_tag = gmsh.model.geo.addCurveLoop(line_tags)
        surface_tag = gmsh.model.geo.addPlaneSurface([loop_tag])
        gmsh.model.geo.synchronize()

        _report_subprogress(
            subprogress_callback,
            "assign physical boundary groups",
            4,
            progress_total,
        )
        gmsh.model.addPhysicalGroup(2, [surface_tag], FLUID_PHYSICAL_TAG)
        gmsh.model.setPhysicalName(2, FLUID_PHYSICAL_TAG, "fluid")
        tag_to_lines: dict[int, list[int]] = {}
        for physical_tag, line_tag in zip(
            blueprint.segment_physical_tag, line_tags, strict=True
        ):
            tag_to_lines.setdefault(int(physical_tag), []).append(int(line_tag))
        for physical_tag, tagged_lines in tag_to_lines.items():
            gmsh.model.addPhysicalGroup(1, tagged_lines, physical_tag)
            if physical_tag == WALL_PHYSICAL_TAG:
                name = "wall"
            elif physical_tag >= OUTLET_PHYSICAL_TAG_OFFSET:
                name = f"outlet_{physical_tag - OUTLET_PHYSICAL_TAG_OFFSET}"
            else:
                name = f"inlet_{physical_tag - INLET_PHYSICAL_TAG_OFFSET}"
            gmsh.model.setPhysicalName(1, physical_tag, name)

        wall_lines = tag_to_lines.get(WALL_PHYSICAL_TAG, [])
        _report_subprogress(
            subprogress_callback,
            "configure wall refinement",
            5,
            progress_total,
        )
        if wall_size < bulk_size and wall_lines:
            distance_field = gmsh.model.mesh.field.add("Distance")
            try:
                gmsh.model.mesh.field.setNumbers(
                    distance_field, "CurvesList", wall_lines
                )
            except Exception:
                gmsh.model.mesh.field.setNumbers(
                    distance_field, "EdgesList", wall_lines
                )
            threshold_field = gmsh.model.mesh.field.add("Threshold")
            gmsh.model.mesh.field.setNumber(threshold_field, "InField", distance_field)
            gmsh.model.mesh.field.setNumber(threshold_field, "SizeMin", wall_size)
            gmsh.model.mesh.field.setNumber(threshold_field, "SizeMax", bulk_size)
            gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", 0.0)
            gmsh.model.mesh.field.setNumber(
                threshold_field, "DistMax", refinement_distance
            )
            gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)

        gmsh.option.setNumber("Mesh.MeshSizeMin", wall_size)
        gmsh.option.setNumber("Mesh.MeshSizeMax", bulk_size)
        gmsh.option.setNumber("Mesh.Algorithm", 6)
        gmsh.option.setNumber("Mesh.Smoothing", 10)
        _report_subprogress(
            subprogress_callback,
            "generate and smooth 2D mesh",
            6,
            progress_total,
        )
        gmsh.model.mesh.generate(2)
        if int(cfg.gmsh_element_order) > 1:
            gmsh.model.mesh.setOrder(int(cfg.gmsh_element_order))
        gmsh_version = str(getattr(gmsh, "__version__", "unknown"))
        _report_subprogress(
            subprogress_callback,
            "convert mesh to DOLFINx",
            7,
            progress_total,
        )
        mesh_data = gmshio.model_to_mesh(gmsh.model, MPI.COMM_SELF, 0, gdim=2)
        if hasattr(mesh_data, "mesh"):
            mesh = mesh_data.mesh
            cell_tags = mesh_data.cell_tags
            facet_tags = mesh_data.facet_tags
        else:
            mesh, cell_tags, facet_tags = mesh_data[:3]
    finally:
        gmsh.finalize()

    if facet_tags is None:
        raise ValueError("DOLFINx received no Gmsh facet markers.")
    return _DolfinxMeshBundle(
        mesh=mesh,
        cell_tags=cell_tags,
        facet_tags=facet_tags,
        blueprint=blueprint,
        gmsh_version=gmsh_version,
        bulk_size_um=bulk_size,
        wall_size_um=wall_size,
        refinement_distance_um=refinement_distance,
    )


def _solve_taylor_hood_stokes(
    *,
    mesh_bundle: _DolfinxMeshBundle,
    continuous_geometry: "ContinuousVesselGeometry",
    vessels: list[Vessel] | tuple[Vessel, ...],
    cfg: FieldConfig,
    optional: dict[str, object],
    subprogress_callback: Callable[[str, int, int], None] | None = None,
) -> dict[str, object]:
    progress_total = 12
    basix_ufl = optional["basix_ufl"]
    fem = optional["fem"]
    LinearProblem = optional["LinearProblem"]
    PETSc = optional["PETSc"]
    scipy_sparse_linalg = optional["scipy_sparse_linalg"]
    ufl = optional["ufl"]
    mesh = mesh_bundle.mesh
    facet_tags = mesh_bundle.facet_tags
    gdim = int(mesh.geometry.dim)
    if gdim != 2:
        raise ValueError(f"Expected a two-dimensional DOLFINx mesh, got gdim={gdim}.")

    _report_subprogress(
        subprogress_callback,
        "create mixed finite-element spaces",
        0,
        progress_total,
    )
    velocity_element = basix_ufl.element(
        "Lagrange",
        mesh.basix_cell(),
        int(cfg.dolfinx_velocity_degree),
        shape=(gdim,),
    )
    pressure_element = basix_ufl.element(
        "Lagrange",
        mesh.basix_cell(),
        int(cfg.dolfinx_pressure_degree),
    )
    mixed_element = basix_ufl.mixed_element([velocity_element, pressure_element])
    W = fem.functionspace(mesh, mixed_element)
    V, _ = W.sub(0).collapse()
    Q, _ = W.sub(1).collapse()
    fdim = int(mesh.topology.dim - 1)

    _report_subprogress(
        subprogress_callback,
        "apply no-slip wall conditions",
        1,
        progress_total,
    )
    bcs = []
    wall_facets = facet_tags.find(WALL_PHYSICAL_TAG)
    if wall_facets.size == 0:
        raise ValueError("The Gmsh mesh contains no tagged solid-wall facets.")
    wall_velocity = fem.Function(V)
    wall_velocity.x.array[:] = 0.0
    wall_dofs = fem.locate_dofs_topological((W.sub(0), V), fdim, wall_facets)
    bcs.append(fem.dirichletbc(wall_velocity, wall_dofs, W.sub(0)))

    vessel_flow = {
        int(vessel.vid): max(float(vessel.flow_rate), 0.0) for vessel in vessels
    }
    section_points = np.asarray(
        continuous_geometry.open_section_point_xz_um, dtype=np.float64
    )
    section_normals = np.asarray(
        continuous_geometry.open_section_outward_normal_xz, dtype=np.float64
    )
    section_tangents = np.asarray(
        continuous_geometry.open_section_tangent_xz, dtype=np.float64
    )
    section_half_widths = np.asarray(
        continuous_geometry.open_section_half_width_um, dtype=np.float64
    )
    section_kinds = np.asarray(continuous_geometry.open_section_kind, dtype=np.int8)
    section_labels = np.asarray(
        continuous_geometry.open_section_label, dtype=np.int32
    )
    section_vessel_ids = np.asarray(
        continuous_geometry.open_section_vessel_id, dtype=np.int32
    )
    section_q2d = np.asarray(
        [
            vessel_flow[int(vessel_id)]
            / max(float(cfg.effective_thickness_um), np.finfo(float).eps)
            for vessel_id in section_vessel_ids
        ],
        dtype=np.float64,
    )
    total_inlet_q2d = float(np.sum(section_q2d[section_kinds < 0]))
    total_outlet_q2d = float(np.sum(section_q2d[section_kinds > 0]))
    flux_scale = max(total_inlet_q2d, total_outlet_q2d, np.finfo(float).eps)
    flux_imbalance = abs(total_inlet_q2d - total_outlet_q2d) / flux_scale
    if flux_imbalance > float(cfg.flux_tolerance):
        raise ValueError(
            "The exported vascular reference flows are incompatible: total "
            "inlet and terminal 2D fluxes "
            f"differ by {flux_imbalance:.6g}, above field.flux_tolerance="
            f"{float(cfg.flux_tolerance):.6g}."
        )
    inlet_target_by_label = np.bincount(
        section_labels[section_kinds < 0],
        weights=section_q2d[section_kinds < 0],
        minlength=int(np.max(section_labels[section_kinds < 0], initial=0)) + 1,
    ).astype(np.float64, copy=False)
    outlet_target_by_label = np.bincount(
        section_labels[section_kinds > 0],
        weights=section_q2d[section_kinds > 0],
        minlength=int(np.max(section_labels[section_kinds > 0], initial=0)) + 1,
    ).astype(np.float64, copy=False)
    _report_subprogress(
        subprogress_callback,
        "apply inlet and outlet conditions",
        2,
        progress_total,
    )
    for (
        physical_tag,
        section_index,
    ) in mesh_bundle.blueprint.physical_tag_to_section_index.items():
        facets = facet_tags.find(int(physical_tag))
        if facets.size == 0:
            raise ValueError(f"Gmsh physical section {physical_tag} has no facets.")
        point = section_points[section_index].copy()
        normal = section_normals[section_index].copy()
        tangent = section_tangents[section_index].copy()
        half_width = float(section_half_widths[section_index])
        kind = int(section_kinds[section_index])
        q2d = float(section_q2d[section_index])
        maximum_speed = 3.0 * q2d / max(4.0 * half_width, np.finfo(float).eps)
        # The stored normal points out of the lumen. Inlets flow opposite that
        # normal; terminal outlets flow along it.
        flow_direction = -normal if kind < 0 else normal

        def parabolic_profile(
            x,
            *,
            point=point,
            tangent=tangent,
            half_width=half_width,
            maximum_speed=maximum_speed,
            flow_direction=flow_direction,
        ):
            transverse = (x[0] - point[0]) * tangent[0] + (x[1] - point[1]) * tangent[1]
            shape = np.maximum(
                1.0 - (transverse / max(half_width, np.finfo(float).eps)) ** 2,
                0.0,
            )
            return np.vstack(
                (
                    maximum_speed * shape * flow_direction[0],
                    maximum_speed * shape * flow_direction[1],
                )
            )

        section_velocity = fem.Function(V)
        section_velocity.interpolate(parabolic_profile)
        section_dofs = fem.locate_dofs_topological((W.sub(0), V), fdim, facets)
        bcs.append(fem.dirichletbc(section_velocity, section_dofs, W.sub(0)))

    # Prescribing velocity on every boundary leaves pressure determined only up
    # to an additive constant. Pin one pressure vertex to zero as a numerical
    # gauge; it does not impose an anatomical inlet or outlet pressure.
    pressure_reference = fem.Function(Q)
    pressure_reference.x.array[:] = 0.0
    mesh.topology.create_connectivity(0, mesh.topology.dim)
    pressure_reference_dofs = fem.locate_dofs_topological(
        (W.sub(1), Q),
        0,
        np.asarray([0], dtype=np.int32),
    )
    if np.asarray(pressure_reference_dofs).size == 0:
        raise ValueError("Could not locate a pressure degree of freedom for the gauge.")
    bcs.append(
        fem.dirichletbc(pressure_reference, pressure_reference_dofs, W.sub(1))
    )

    _report_subprogress(
        subprogress_callback,
        "construct Stokes variational forms",
        3,
        progress_total,
    )
    trial_velocity, trial_pressure = ufl.TrialFunctions(W)
    test_velocity, test_pressure = ufl.TestFunctions(W)
    dx = ufl.Measure("dx", domain=mesh)
    scalar_type = PETSc.ScalarType if PETSc is not None else np.float64
    viscosity = fem.Constant(mesh, scalar_type(float(cfg.kinematic_viscosity_um2_s)))
    zero_force = fem.Constant(mesh, np.zeros(gdim, dtype=scalar_type))
    a = (
        2.0
        * viscosity
        * ufl.inner(ufl.sym(ufl.grad(trial_velocity)), ufl.sym(ufl.grad(test_velocity)))
        * dx
        - trial_pressure * ufl.div(test_velocity) * dx
        + test_pressure * ufl.div(trial_velocity) * dx
    )
    L = ufl.inner(zero_force, test_velocity) * dx
    if LinearProblem is not None and PETSc is not None:
        solution, linear_diagnostics = _solve_with_petsc(
            a=a,
            L=L,
            bcs=bcs,
            W=W,
            cfg=cfg,
            LinearProblem=LinearProblem,
            subprogress_callback=subprogress_callback,
            progress_total=progress_total,
        )
    else:
        solution, linear_diagnostics = _solve_with_scipy(
            a=a,
            L=L,
            bcs=bcs,
            W=W,
            cfg=cfg,
            fem=fem,
            scipy_sparse_linalg=scipy_sparse_linalg,
            subprogress_callback=subprogress_callback,
            progress_total=progress_total,
        )
    _report_subprogress(
        subprogress_callback,
        "extract velocity and pressure fields",
        9,
        progress_total,
    )
    solution.x.scatter_forward()
    velocity = solution.sub(0).collapse()
    pressure = solution.sub(1).collapse()
    _report_subprogress(
        subprogress_callback,
        "integrate open-boundary fluxes",
        10,
        progress_total,
    )
    inlet_actual, outlet_actual = _integrate_open_section_fluxes(
        velocity=velocity,
        mesh=mesh,
        facet_tags=facet_tags,
        blueprint=mesh_bundle.blueprint,
        section_kinds=section_kinds,
        section_labels=section_labels,
        fem=fem,
        ufl=ufl,
    )

    gradient_degree = max(0, int(cfg.dolfinx_velocity_degree) - 1)
    _report_subprogress(
        subprogress_callback,
        "project the velocity gradient",
        11,
        progress_total,
    )
    gradient_element = basix_ufl.element(
        "DG", mesh.basix_cell(), gradient_degree, shape=(gdim, gdim)
    )
    gradient_space = fem.functionspace(mesh, gradient_element)
    interpolation_points = gradient_space.element.interpolation_points
    if callable(interpolation_points):
        interpolation_points = interpolation_points()
    gradient_expression = fem.Expression(ufl.grad(velocity), interpolation_points)
    velocity_gradient = fem.Function(gradient_space)
    velocity_gradient.interpolate(gradient_expression)
    velocity_gradient.x.scatter_forward()
    return {
        "velocity": velocity,
        "pressure": pressure,
        "velocity_gradient": velocity_gradient,
        "total_inlet_q2d_um2_s": total_inlet_q2d,
        "total_outlet_q2d_um2_s": total_outlet_q2d,
        "flux_balance_relative_error": flux_imbalance,
        "inlet_target_by_label_um2_s": inlet_target_by_label,
        "outlet_target_by_label_um2_s": outlet_target_by_label,
        "inlet_actual_by_label_um2_s": inlet_actual,
        "outlet_actual_by_label_um2_s": outlet_actual,
        **linear_diagnostics,
    }


def _integrate_open_section_fluxes(
    *,
    velocity,
    mesh,
    facet_tags,
    blueprint: BoundaryMeshBlueprint,
    section_kinds: np.ndarray,
    section_labels: np.ndarray,
    fem,
    ufl,
) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the solved FEM velocity on every tagged inlet and outlet."""

    maximum_inlet_label = int(
        np.max(section_labels[section_kinds < 0], initial=0)
    )
    maximum_outlet_label = int(
        np.max(section_labels[section_kinds > 0], initial=0)
    )
    inlet = np.zeros(maximum_inlet_label + 1, dtype=np.float64)
    outlet = np.zeros(maximum_outlet_label + 1, dtype=np.float64)
    ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags)
    outward_normal = ufl.FacetNormal(mesh)
    normal_flux = ufl.dot(velocity, outward_normal)

    for physical_tag, section_index in blueprint.physical_tag_to_section_index.items():
        local_flux = float(
            fem.assemble_scalar(
                fem.form(normal_flux * ds(int(physical_tag)))
            )
        )
        communicator = getattr(mesh, "comm", None)
        signed_outward_flux = (
            float(communicator.allreduce(local_flux))
            if communicator is not None
            else local_flux
        )
        kind = int(section_kinds[section_index])
        label = int(section_labels[section_index])
        if kind < 0:
            inlet[label] += -signed_outward_flux
        else:
            outlet[label] += signed_outward_flux
    return inlet, outlet


def _maximum_labeled_flux_relative_error(
    actual_by_label: np.ndarray,
    target_by_label: np.ndarray,
) -> float:
    """Return the worst relative error across positive boundary labels."""

    actual = np.asarray(actual_by_label, dtype=np.float64).reshape(-1)
    target = np.asarray(target_by_label, dtype=np.float64).reshape(-1)
    size = max(actual.size, target.size)
    if size <= 1:
        return float("inf")
    actual_padded = np.zeros(size, dtype=np.float64)
    target_padded = np.zeros(size, dtype=np.float64)
    actual_padded[: actual.size] = actual
    target_padded[: target.size] = target
    positive = target_padded[1:] > np.finfo(np.float64).eps
    if not np.any(positive):
        return float("inf")
    errors = np.abs(actual_padded[1:][positive] - target_padded[1:][positive])
    errors /= target_padded[1:][positive]
    return float(np.max(errors, initial=0.0))


def _solve_with_petsc(
    *,
    a,
    L,
    bcs,
    W,
    cfg,
    LinearProblem,
    subprogress_callback: Callable[[str, int, int], None] | None = None,
    progress_total: int = 12,
):
    """Solve the mixed system with PETSc when petsc4py is available."""

    petsc_options = {
        "ksp_type": "preonly",
        "pc_type": "lu",
        "pc_factor_mat_solver_type": "mumps",
        "ksp_error_if_not_converged": True,
    }
    _report_subprogress(
        subprogress_callback,
        "compile and assemble PETSc system",
        4,
        progress_total,
    )
    try:
        problem = LinearProblem(
            a,
            L,
            bcs=bcs,
            petsc_options_prefix="vessel_stokes_",
            petsc_options=petsc_options,
        )
    except TypeError:
        problem = LinearProblem(a, L, bcs=bcs, petsc_options=petsc_options)
    _report_subprogress(
        subprogress_callback,
        "solve PETSc linear system",
        8,
        progress_total,
    )
    solution = problem.solve()
    solver = problem.solver
    converged_reason = int(solver.getConvergedReason())
    return solution, {
        "physical_acceptance_schema": "dolfinx_fixed_flow_petsc_stokes_v3",
        "linear_solver_backend": "petsc_preonly_mumps_lu",
        "linear_solver_converged": converged_reason > 0,
        "linear_solver_converged_reason": converged_reason,
        "linear_solver_iterations": int(solver.getIterationNumber()),
        "linear_solver_relative_residual": float("nan"),
        "linear_solver_residual_tolerance": float(cfg.dolfinx_ksp_rtol),
        "linear_solver_residual_norm": float(solver.getResidualNorm()),
        "ksp_converged_reason": converged_reason,
        "ksp_iterations": int(solver.getIterationNumber()),
    }


def _solve_with_scipy(
    *,
    a,
    L,
    bcs,
    W,
    cfg,
    fem,
    scipy_sparse_linalg,
    subprogress_callback: Callable[[str, int, int], None] | None = None,
    progress_total: int = 12,
):
    """
    Assemble with DOLFINx and solve serially through SciPy.

    This is the supported fallback for Windows DOLFINx builds, for which
    PETSc/petsc4py is not distributed.  Gmsh meshing, finite-element spaces,
    variational forms, boundary conditions, and assembly remain in DOLFINx.
    """

    if scipy_sparse_linalg is None:
        raise ImportError(
            "This DOLFINx build has no petsc4py. Install scipy (or pyamg, which "
            "depends on scipy) to enable the serial Windows solver."
        )
    if int(W.mesh.comm.size) != 1:
        raise ValueError(
            "The SciPy DOLFINx fallback is serial-only. Install petsc4py or run "
            "with one MPI rank."
        )
    _report_subprogress(
        subprogress_callback,
        "JIT-compile the bilinear form",
        4,
        progress_total,
    )
    a_form = fem.form(a)
    _report_subprogress(
        subprogress_callback,
        "JIT-compile the linear form",
        5,
        progress_total,
    )
    L_form = fem.form(L)
    _report_subprogress(
        subprogress_callback,
        "assemble the sparse Stokes matrix",
        6,
        progress_total,
    )
    assembled_matrix = fem.assemble_matrix(a_form, bcs=bcs)
    matrix = assembled_matrix.to_scipy().tocsr()
    nonfinite_entries = np.flatnonzero(~np.isfinite(matrix.data))
    if nonfinite_entries.size:
        raise ValueError(
            "The assembled mixed Stokes matrix contains non-finite entries; "
            f"count={nonfinite_entries.size}, first stored indices="
            f"{nonfinite_entries[:10].tolist()}."
        )
    zero_rows = np.flatnonzero(np.diff(matrix.indptr) == 0)
    if zero_rows.size:
        raise ValueError(
            "The assembled mixed Stokes matrix contains unconstrained zero "
            f"rows; first rows={zero_rows[:10].tolist()}."
        )
    _report_subprogress(
        subprogress_callback,
        "assemble and constrain the right-hand side",
        7,
        progress_total,
    )
    rhs = fem.assemble_vector(L_form)
    fem.apply_lifting(rhs.array, [a_form], bcs=[bcs])
    for bc in bcs:
        bc.set(rhs.array)

    _report_subprogress(
        subprogress_callback,
        "solve the sparse linear system",
        8,
        progress_total,
    )
    solution_values = np.asarray(
        scipy_sparse_linalg.spsolve(matrix, rhs.array), dtype=np.float64
    )
    solution = fem.Function(W)
    solution.x.array[:] = solution_values
    residual = np.asarray(matrix @ solution_values - rhs.array, dtype=np.float64)
    rhs_norm = max(float(np.linalg.norm(rhs.array)), np.finfo(np.float64).eps)
    relative_residual = float(np.linalg.norm(residual) / rhs_norm)
    residual_limit = max(100.0 * float(cfg.dolfinx_ksp_rtol), 1.0e-10)
    converged = bool(
        np.all(np.isfinite(solution_values))
        and np.isfinite(relative_residual)
        and relative_residual <= residual_limit
    )
    return solution, {
        "physical_acceptance_schema": "dolfinx_fixed_flow_scipy_stokes_v3",
        "linear_solver_backend": "scipy_sparse_direct",
        "linear_solver_converged": converged,
        "linear_solver_converged_reason": 1 if converged else -1,
        "linear_solver_iterations": 1,
        "linear_solver_relative_residual": relative_residual,
        "linear_solver_residual_tolerance": residual_limit,
    }


def _cartesian_lumen_points(
    domain: GridDomain, lumen_mask: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    indices = np.argwhere(np.asarray(lumen_mask, dtype=bool)).astype(
        np.int32, copy=False
    )
    points = np.zeros((indices.shape[0], 3), dtype=np.float64)
    points[:, 0] = np.asarray(domain.x_coordinates_um)[indices[:, 0]]
    points[:, 1] = np.asarray(domain.z_coordinates_um)[indices[:, 1]]
    return points, indices


def _cartesian_sampling_failure_diagnostics(
    *,
    domain: GridDomain,
    continuous_geometry: "ContinuousVesselGeometry",
    mesh,
    sample_points: np.ndarray,
    sample_indices: np.ndarray,
    wall_state,
    wall_sample_points: np.ndarray,
    wall_sample_offset_um: float,
    located_velocity: np.ndarray,
    located_pressure: np.ndarray,
    located_local_gradient: np.ndarray,
    located_wall_gradient: np.ndarray,
    maximum_point_rows: int = 12,
) -> tuple[str, str]:
    """Describe exactly which Cartesian or near-wall mesh queries failed."""

    centre_located = (
        np.asarray(located_velocity, dtype=bool)
        & np.asarray(located_pressure, dtype=bool)
        & np.asarray(located_local_gradient, dtype=bool)
    )
    wall_located = np.asarray(located_wall_gradient, dtype=bool)
    located = centre_located & wall_located
    missing_rows = np.flatnonzero(~located)

    centre_xz = np.asarray(sample_points, dtype=np.float64)[:, :2]
    continuous_inside = np.asarray(
        continuous_geometry.contains_xz_um(centre_xz), dtype=bool
    ).reshape(-1)
    outside_count = int(np.count_nonzero(~continuous_inside))

    mesh_coordinates = np.asarray(mesh.geometry.x, dtype=np.float64)
    if mesh_coordinates.size:
        mesh_minimum = np.min(mesh_coordinates[:, :2], axis=0)
        mesh_maximum = np.max(mesh_coordinates[:, :2], axis=0)
    else:
        mesh_minimum = np.full(2, np.nan, dtype=np.float64)
        mesh_maximum = np.full(2, np.nan, dtype=np.float64)

    if outside_count:
        cause = (
            f"raster.lumen_mask contains {outside_count} grid centre(s) outside "
            "the authoritative continuous lumen; this is a raster/geometry "
            "consistency error, not evidence that the Gmsh mesh is too coarse."
        )
    elif np.any(~centre_located):
        cause = (
            "At least one continuous-lumen centre was not located in the Gmsh "
            "mesh; inspect boundary coincidence and mesh collision tolerances."
        )
    else:
        cause = (
            "All lumen centres were located, but at least one inward-offset "
            "near-wall gradient point was not; inspect its wall projection, "
            "normal, and sampling offset."
        )

    def _count_missing(values: np.ndarray) -> int:
        return int(np.count_nonzero(~np.asarray(values, dtype=bool)))

    lines = [
        "",
        "[DOLFINx Cartesian sampling diagnostics]",
        f"  accepted lumen samples: {sample_points.shape[0]}",
        (
            "  location failures: "
            f"velocity={_count_missing(located_velocity)}, "
            f"pressure={_count_missing(located_pressure)}, "
            f"centre_gradient={_count_missing(located_local_gradient)}, "
            f"wall_gradient={_count_missing(located_wall_gradient)}, "
            f"union={missing_rows.size}"
        ),
        (
            "  accepted grid centres outside continuous lumen: "
            f"{outside_count}"
        ),
        (
            "  domain: "
            f"shape={tuple(int(value) for value in domain.shape)}, "
            f"spacing_um={float(domain.spacing_um):.17g}, "
            f"origin_xz_um=({float(domain.origin_um[0]):.17g}, "
            f"{float(domain.origin_um[2]):.17g})"
        ),
        (
            "  Gmsh mesh x-z bounds um: "
            f"[({mesh_minimum[0]:.17g}, {mesh_minimum[1]:.17g}), "
            f"({mesh_maximum[0]:.17g}, {mesh_maximum[1]:.17g})]"
        ),
        f"  wall sampling offset um: {float(wall_sample_offset_um):.17g}",
        f"  classification: {cause}",
    ]

    nearest_points = np.asarray(
        wall_state.nearest_point_xz_um, dtype=np.float64
    ).reshape(-1, 2)
    inward_normals = np.asarray(
        wall_state.inward_normal_xz, dtype=np.float64
    ).reshape(-1, 2)
    wall_distances = np.asarray(
        wall_state.distance_um, dtype=np.float64
    ).reshape(-1)
    unique_nearest = np.asarray(
        wall_state.unique_nearest_wall, dtype=bool
    ).reshape(-1)
    displayed_rows = missing_rows[:maximum_point_rows]
    signed_distances = np.asarray(
        continuous_geometry.signed_lumen_distance_at_xz_um(
            centre_xz[displayed_rows]
        ),
        dtype=np.float64,
    ).reshape(-1)

    for display_index, (row, signed_distance) in enumerate(
        zip(displayed_rows, signed_distances, strict=True),
        start=1,
    ):
        grid_index = np.asarray(sample_indices[row], dtype=np.int64)
        centre = centre_xz[row]
        nearest = nearest_points[row]
        normal = inward_normals[row]
        wall_sample = np.asarray(wall_sample_points[row, :2], dtype=np.float64)
        lines.extend(
            [
                f"  failed sample {display_index}:",
                (
                    f"    flat_row={int(row)}, "
                    f"grid_index=({int(grid_index[0])}, {int(grid_index[1])})"
                ),
                (
                    f"    centre_xz_um=({centre[0]:.17g}, {centre[1]:.17g}), "
                    f"continuous_inside={bool(continuous_inside[row])}, "
                    f"signed_lumen_distance_um={signed_distance:.17g}"
                ),
                (
                    f"    nearest_wall_xz_um=({nearest[0]:.17g}, "
                    f"{nearest[1]:.17g}), "
                    f"wall_distance_um={wall_distances[row]:.17g}, "
                    f"unique_nearest_wall={bool(unique_nearest[row])}"
                ),
                (
                    f"    inward_normal_xz=({normal[0]:.17g}, "
                    f"{normal[1]:.17g}), "
                    f"wall_sample_xz_um=({wall_sample[0]:.17g}, "
                    f"{wall_sample[1]:.17g})"
                ),
                (
                    "    located: "
                    f"velocity={bool(located_velocity[row])}, "
                    f"pressure={bool(located_pressure[row])}, "
                    f"centre_gradient={bool(located_local_gradient[row])}, "
                    f"wall_gradient={bool(located_wall_gradient[row])}"
                ),
            ]
        )
    if missing_rows.size > maximum_point_rows:
        lines.append(
            f"  ... {missing_rows.size - maximum_point_rows} additional "
            "failed samples omitted"
        )
    return "\n".join(lines), cause


def _export_finite_element_velocity(
    velocity_function,
    mesh,
    *,
    degree: int,
    preferred_bin_size_um: float,
):
    """Reconstruct the solved velocity as one polynomial per affine triangle."""

    geometry_dofmap = np.asarray(mesh.geometry.dofmap, dtype=np.int64)
    if geometry_dofmap.ndim != 2 or geometry_dofmap.shape[1] != 3:
        raise ValueError(
            "Hybrid velocity export requires affine three-node triangle "
            "geometry. Keep field.gmsh_element_order equal to 1."
        )
    geometry_points = np.asarray(mesh.geometry.x, dtype=np.float64)
    cell_vertices = np.ascontiguousarray(
        geometry_points[geometry_dofmap, :2],
        dtype=np.float64,
    )
    reference_nodes = interpolation_nodes(int(degree))
    cell_count = int(cell_vertices.shape[0])
    node_count = int(reference_nodes.shape[0])
    sampled = np.empty((cell_count, node_count, 2), dtype=np.float64)
    chunk_cells = max(1, 250_000 // node_count)
    for start in range(0, cell_count, chunk_cells):
        stop = min(start + chunk_cells, cell_count)
        local_vertices = cell_vertices[start:stop]
        physical_xz = (
            local_vertices[:, None, 0, :]
            + reference_nodes[None, :, 0, None]
            * (
                local_vertices[:, None, 1, :]
                - local_vertices[:, None, 0, :]
            )
            + reference_nodes[None, :, 1, None]
            * (
                local_vertices[:, None, 2, :]
                - local_vertices[:, None, 0, :]
            )
        )
        points = np.zeros((physical_xz.shape[0] * node_count, 3), dtype=np.float64)
        points[:, :2] = physical_xz.reshape(-1, 2)
        cells = np.repeat(
            np.arange(start, stop, dtype=np.int32),
            node_count,
        )
        values = velocity_function.eval(points, cells)
        sampled[start:stop] = np.asarray(values, dtype=np.float64).reshape(
            stop - start, node_count, 2
        )
    return build_finite_element_velocity_field(
        cell_vertices,
        sampled,
        int(degree),
        preferred_bin_size_um=float(preferred_bin_size_um),
    )


def _sample_dolfinx_function(function, mesh, points, optional):
    geometry = optional["geometry"]
    points = np.ascontiguousarray(points, dtype=np.float64)
    tree = geometry.bb_tree(mesh, mesh.topology.dim)
    candidates = geometry.compute_collisions_points(tree, points)
    colliding = geometry.compute_colliding_cells(mesh, candidates, points)
    cells = np.full(points.shape[0], -1, dtype=np.int32)
    for index in range(points.shape[0]):
        links = colliding.links(index)
        if len(links):
            cells[index] = int(links[0])
    located = cells >= 0
    value_shape = tuple(int(value) for value in function.ufl_shape)
    value_size = int(np.prod(value_shape, dtype=np.int64)) if value_shape else 1
    values = np.zeros((points.shape[0], value_size), dtype=np.float64)
    if np.any(located):
        evaluated = function.eval(points[located], cells[located])
        values[located] = np.asarray(evaluated).reshape(-1, value_size)
    return values, located


def _continuous_wall_shear_on_grid(
    *,
    gradient: np.ndarray,
    raster: RasterizedVessels,
    domain: GridDomain,
    continuous_geometry: "ContinuousVesselGeometry",
) -> np.ndarray:
    indices = np.argwhere(raster.lumen_mask)
    result = np.zeros(domain.shape, dtype=np.float64)
    if indices.size == 0:
        return result
    points = np.column_stack(
        (
            np.asarray(domain.x_coordinates_um)[indices[:, 0]],
            np.asarray(domain.z_coordinates_um)[indices[:, 1]],
        )
    )
    wall_state = continuous_geometry.exact_solid_wall_state_xz_um_accelerated(points)
    normals = np.asarray(wall_state.inward_normal_xz, dtype=np.float64).reshape(-1, 2)
    tangents = np.column_stack((-normals[:, 1], normals[:, 0]))
    local_gradient = gradient[indices[:, 0], indices[:, 1]]
    strain_twice = local_gradient + np.swapaxes(local_gradient, 1, 2)
    strain_on_normal = np.einsum("nij,nj->ni", strain_twice, normals)
    shear_rate = np.abs(np.einsum("ni,ni->n", tangents, strain_on_normal))
    viscosity_mpas = np.asarray(raster.viscosity_mpas, dtype=np.float64)[
        indices[:, 0], indices[:, 1]
    ]
    positive = viscosity_mpas[np.isfinite(viscosity_mpas) & (viscosity_mpas > 0.0)]
    default_viscosity = float(np.median(positive)) if positive.size else 3.0
    viscosity_pa_s = (
        np.where(viscosity_mpas > 0.0, viscosity_mpas, default_viscosity) * 1.0e-3
    )
    result[indices[:, 0], indices[:, 1]] = viscosity_pa_s * shear_rate
    return result


def _local_shear_stress_on_grid(
    *,
    gradient: np.ndarray,
    raster: RasterizedVessels,
    domain: GridDomain,
) -> np.ndarray:
    """Return Newtonian local viscous shear magnitude throughout the lumen.

    The scalar is ``mu * sqrt(2 * D:D)``, where
    ``D = 0.5 * (grad(u) + grad(u).T)``.  It equals ``mu * |du/dy|`` for
    ordinary simple shear, is zero for rigid-body rotation, and is distinct
    from wall shear stress, which is defined only on a solid wall.
    """

    gradient_array = np.asarray(gradient, dtype=np.float64)
    expected_shape = (*domain.shape, 2, 2)
    if gradient_array.shape != expected_shape:
        raise ValueError(
            f"gradient must have shape {expected_shape}, got {gradient_array.shape}."
        )

    lumen = np.asarray(raster.lumen_mask, dtype=bool)
    viscosity_mpas = np.asarray(raster.viscosity_mpas, dtype=np.float64)
    if lumen.shape != domain.shape or viscosity_mpas.shape != domain.shape:
        raise ValueError("Raster lumen and viscosity arrays must match the grid domain.")

    result = np.zeros(domain.shape, dtype=np.float64)
    if not np.any(lumen):
        return result

    local_gradient = gradient_array[lumen]
    strain_rate = 0.5 * (
        local_gradient + np.swapaxes(local_gradient, 1, 2)
    )
    strain_invariant = np.einsum("nij,nij->n", strain_rate, strain_rate)
    equivalent_shear_rate = np.sqrt(np.maximum(2.0 * strain_invariant, 0.0))

    local_viscosity_mpas = viscosity_mpas[lumen]
    positive = local_viscosity_mpas[
        np.isfinite(local_viscosity_mpas) & (local_viscosity_mpas > 0.0)
    ]
    default_viscosity_mpas = float(np.median(positive)) if positive.size else 3.0
    viscosity_pa_s = (
        np.where(
            np.isfinite(local_viscosity_mpas) & (local_viscosity_mpas > 0.0),
            local_viscosity_mpas,
            default_viscosity_mpas,
        )
        * 1.0e-3
    )
    result[lumen] = viscosity_pa_s * equivalent_shear_rate
    if not np.isfinite(result[lumen]).all():
        raise ValueError("Local shear stress contains non-finite lumen values.")
    return result
