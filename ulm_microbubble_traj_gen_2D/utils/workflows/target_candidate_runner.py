"""Flow-only preparation workflow for manual or automatic molecular targets."""

from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path

from ..core.config import PhysicsConfig
from ..flow.connectivity import validate_fluid_connectivity
from ..flow.dolfinx_gmsh_solver import solve_dolfinx_stokes_gmsh_2d
from ..flow.flow_diagnostics import write_flow_diagnostics
from ..geometry.continuous_vessel_geometry import build_continuous_vessel_geometry
from ..geometry.grid_domain import build_domain_from_vessels
from ..geometry.vessel_rasterizer import rasterize_vessels
from ..io.field_io import (
    particle_boundary_metadata,
    save_domain_metadata,
    save_field_npz,
    save_run_config,
)
from ..io.vascular_io import load_physics_input, validate_physics_input_geometry
from ..molecular.molecular_target_auto_selection import (
    select_automatic_influence_anchor,
)
from ..molecular.molecular_target_candidate_io import (
    save_candidate_catalog,
    save_spatially_heterogeneous_target_mask,
    save_spatially_heterogeneous_target_report,
)
from ..molecular.molecular_target_candidates import build_molecular_target_candidates
from ..molecular.molecular_target_spatial_heterogeneity import (
    build_spatially_heterogeneous_target,
)
from ..particles.particle_hydrodynamic_fields import build_particle_hydrodynamic_fields
from ..particles.particle_inlet_flux import build_inlet_flux_model
from ..runtime.console_output import print_key_values, print_stage, print_warning
from ..visualization.vtk.pyvista_flow import validate_cfd_flow_dependencies
from ..visualization.vtk.pyvista_wall_shear import render_wall_shear_visualization
from ..visualization.vtk.vtk_flow_grid import build_vtk_stage_grid


@dataclass(frozen=True)
class TargetCandidatePreparationResult:
    """Candidate files and any configured automatic synthetic target."""

    output_dir: Path
    field_npz_path: Path
    candidate_npz_path: Path
    candidate_json_path: Path
    final_flow_vti_path: Path
    candidate_count: int
    automatic_target_mask_path: Path | None
    automatic_selection_report_path: Path | None


def run_target_candidate_preparation(
    cfg: PhysicsConfig,
) -> TargetCandidatePreparationResult:
    """Solve the accepted steady field and build selectable vessel-bed candidates."""

    total_stages = 6
    print_stage(1, total_stages, "Initialize target-candidate preparation")
    if not cfg.save_npz:
        raise ValueError(
            "Target-candidate preparation requires output.save_npz=true so the selector can load its catalog."
        )
    validate_cfd_flow_dependencies()
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    if cfg.save_run_config:
        resolved_run_config = deepcopy(cfg.raw)
        resolved_run_config["_resolved_source_config_path"] = str(cfg.source_path)
        resolved_run_config["_workflow"] = "molecular_target_candidate_preparation"
        save_run_config(cfg.output_dir / "run_config.yaml", resolved_run_config)

    print_stage(
        2,
        total_stages,
        "Load exported vascular model",
        detail=str(cfg.model_dir),
    )
    physics_input = load_physics_input(
        cfg.model_dir,
        planar_extrusion_depth_um=cfg.field.effective_thickness_um,
    )
    validate_physics_input_geometry(
        cfg.model_dir,
        physics_input,
        expected_mode="planar_2d",
    )

    print_stage(3, total_stages, "Build and validate the X-Z lumen domain")
    domain = build_domain_from_vessels(physics_input.vessels, cfg.domain)
    print_key_values(
        [
            ("Grid shape", f"{domain.shape[0]} x {domain.shape[1]}"),
            ("Grid spacing", f"{domain.spacing_um:.3f} um"),
        ]
    )
    continuous_geometry = build_continuous_vessel_geometry(
        physics_input.vessels,
        domain,
        maximum_boundary_element_length_um=(
            cfg.domain.continuous_boundary_maximum_element_length_um
        ),
    )
    raster = rasterize_vessels(
        physics_input.vessels,
        domain,
        cfg.domain,
        effective_thickness_um=cfg.field.effective_thickness_um,
        continuous_geometry=continuous_geometry,
        dynamic_viscosity_mpas=(
            cfg.field.kinematic_viscosity_um2_s
            * cfg.field.blood_density_kg_m3
            * 1.0e-9
        ),
    )
    validate_fluid_connectivity(
        raster,
        domain,
        continuous_geometry=continuous_geometry,
    )

    print_stage(4, total_stages, "Solve and diagnose the velocity field")
    flow = solve_dolfinx_stokes_gmsh_2d(
        domain,
        raster,
        cfg.field,
        physics_input.vessels,
        continuous_geometry,
        vessel_metadata=physics_input.vessel_metadata,
    )
    write_flow_diagnostics(
        cfg.output_dir,
        domain,
        raster,
        flow,
        physics_input.vessels,
    )

    print_stage(5, total_stages, "Build candidate vessel-bed catalog")
    hydrodynamic_fields = build_particle_hydrodynamic_fields(
        domain,
        raster,
        flow,
        continuous_geometry=continuous_geometry,
    )
    boundary_geometry = hydrodynamic_fields.boundary_geometry
    inlet_flux = build_inlet_flux_model(
        domain,
        flow,
        physics_input.vessels,
        boundary_geometry,
        cfg.particles,
        effective_thickness_um=cfg.field.effective_thickness_um,
        boundary_depth_cells=cfg.field.boundary_depth_cells,
    )
    observation_time_s = cfg.particles.n_steps * cfg.particles.dt_s
    catalog = build_molecular_target_candidates(
        domain,
        raster,
        flow,
        hydrodynamic_fields,
        physics_input.vessels,
        injection_rate_per_s=inlet_flux.injection_rate_per_s,
        observation_time_s=observation_time_s,
    )
    field_path = cfg.output_dir / "velocity_and_wall_shear_field.npz"
    candidate_npz_path = cfg.output_dir / "molecular_target_candidates.npz"
    candidate_json_path = cfg.output_dir / "molecular_target_candidates.json"
    automatic_target_mask_path = None
    automatic_selection_report_path = None
    save_field_npz(
        field_path,
        domain,
        raster,
        flow,
        continuous_geometry=continuous_geometry,
    )
    save_candidate_catalog(candidate_npz_path, candidate_json_path, catalog)
    if cfg.molecular_target_selection.default_mode == "automatic":
        influence_fraction = (
            cfg.molecular_target_selection
            .influence_region_endothelial_wall_area_fraction
        )
        positive_fraction = (
            cfg.molecular_target_selection
            .target_positive_wall_fraction_within_influence
        )
        correlation_length_um = (
            cfg.molecular_target_selection.target_correlation_length_um
        )
        automatic_anchor = select_automatic_influence_anchor(
            catalog,
            influence_fraction,
        )
        automatic_target = build_spatially_heterogeneous_target(
            catalog,
            automatic_anchor,
            influence_wall_area_fraction=influence_fraction,
            positive_wall_fraction_within_influence=positive_fraction,
            correlation_length_um=correlation_length_um,
            random_seed=cfg.molecular_target_selection.random_seed,
            random_field_modes=cfg.molecular_target_selection.random_field_modes,
        )
        automatic_target_mask_path = (
            cfg.output_dir / "selected_molecular_target_mask.npz"
        )
        automatic_selection_report_path = (
            cfg.output_dir / "automatic_molecular_target_selection.json"
        )
        save_spatially_heterogeneous_target_mask(
            automatic_target_mask_path,
            catalog,
            automatic_target,
        )
        save_spatially_heterogeneous_target_report(
            automatic_selection_report_path,
            catalog,
            automatic_anchor,
            automatic_target,
        )
        print_key_values(
            [
                ("Influence anchor", automatic_anchor.anchor_candidate_id),
                ("Influence radius", f"{automatic_target.influence_radius_um:.6g} um"),
                (
                    "Positive wall fraction",
                    f"{automatic_target.achieved_positive_wall_fraction_within_influence:.6g}",
                ),
                ("Spatial target patches", automatic_target.patch_count),
                ("Random realization seed", automatic_target.random_seed),
            ]
        )
        if automatic_target.correlation_length_grid_cells < 3.0:
            print_warning(
                "The target correlation length spans fewer than three grid cells. "
                "Repeat at finer grid resolution before interpreting patch morphology."
            )
        if catalog.unmapped_endothelial_wall_area_um2 > max(
            1.0e-9,
            1.0e-10 * catalog.network_endothelial_wall_area_um2,
        ):
            print_warning(
                "Some theoretical endothelial area has no solid-wall sample on this grid. "
                "The audit report records that area; refine the grid for quantitative "
                "wall-area fractions."
            )
    save_domain_metadata(
        cfg.output_dir / "domain_metadata.yaml",
        domain,
        {
            "workflow": "molecular_target_candidate_preparation",
            "swc_path": str(physics_input.swc_path),
            "vessel_data_path": str(physics_input.vessel_data_path),
            "vascular_geometry_mode": "planar_2d",
            "flow_solver_dimension": 2,
            "flow_solver_paradigm": "boundary_fitted_dolfinx_stokes_xz_2d",
            "planar_effective_thickness_mode": (
                "configured_extrusion_depth_for_true_planar_flux"
                if physics_input.vessel_metadata.get("source_flow_quantity")
                == "planar_flux_per_unit_depth"
                else "legacy_inlet_circular_area_over_planar_width"
            ),
            "planar_effective_thickness_um": float(
                cfg.field.effective_thickness_um
            ),
            "source_flow_quantity": physics_input.vessel_metadata.get(
                "source_flow_quantity",
                physics_input.vessel_metadata.get("flow_quantity", "legacy_volume_flow"),
            ),
            "n_vessels": len(physics_input.vessels),
            "candidate_count": len(catalog.candidates),
            "vessel_bed_unit_count": len(catalog.topology.units),
            "unresolved_junction_cells_filled_from_nearest_flow_resolved_basin": int(
                catalog.unresolved_junction_cells
            ),
            "candidate_npz_path": str(candidate_npz_path),
            "candidate_json_path": str(candidate_json_path),
            "candidate_injection_rate_per_s": float(inlet_flux.injection_rate_per_s),
            "candidate_observation_time_s": float(observation_time_s),
            **particle_boundary_metadata(boundary_geometry),
            "automatic_target_mask_path": (
                str(automatic_target_mask_path)
                if automatic_target_mask_path is not None
                else None
            ),
            "automatic_selection_report_path": (
                str(automatic_selection_report_path)
                if automatic_selection_report_path is not None
                else None
            ),
            "field_solver_metadata": dict(flow.solver_metadata),
        },
    )

    print_stage(6, total_stages, "Render candidate-selection background fields")
    render_wall_shear_visualization(
        cfg.output_dir,
        domain=domain,
        raster=raster,
        flow=flow,
    )
    final_flow_vti_path = cfg.output_dir / "final_flow_field.vti"
    final_stage_grid = build_vtk_stage_grid(
        domain,
        raster,
        flow,
        stage="final",
        include_lic=False,
    )
    final_stage_grid.image_grid.save(final_flow_vti_path)
    return TargetCandidatePreparationResult(
        output_dir=cfg.output_dir,
        field_npz_path=field_path,
        candidate_npz_path=candidate_npz_path,
        candidate_json_path=candidate_json_path,
        final_flow_vti_path=final_flow_vti_path,
        candidate_count=len(catalog.candidates),
        automatic_target_mask_path=automatic_target_mask_path,
        automatic_selection_report_path=automatic_selection_report_path,
    )
