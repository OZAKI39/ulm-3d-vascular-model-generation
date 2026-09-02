"""
Run the DOLFINx flow and accelerated microbubble pipeline.
"""

from copy import deepcopy
import numpy as np
from ..flow.connectivity import validate_fluid_connectivity
from ..flow.dolfinx_gmsh_solver import solve_dolfinx_stokes_gmsh_2d
from ..flow.flow_diagnostics import write_flow_diagnostics
from ..geometry.continuous_vessel_geometry import build_continuous_vessel_geometry
from ..geometry.grid_domain import build_domain_from_vessels
from ..geometry.vessel_rasterizer import rasterize_vessels
from ..io.field_io import (
    build_field_reuse_config_contract,
    load_reusable_flow_field,
    particle_boundary_metadata,
    save_domain_metadata,
    save_field_npz,
    save_molecular_target_npz,
    save_red_blood_cell_transport_npz,
    save_run_config,
    save_trajectories_npz,
)
from ..io.vascular_io import load_physics_input, validate_physics_input_geometry
from ..molecular.molecular_target_field import build_molecular_target_field
from ..particles.particle_hydrodynamic_fields import build_particle_hydrodynamic_fields
from ..particles.particle_perfusion_transport import advect_particles_with_continuous_perfusion
from ..particles.particle_topological_ownership import build_topological_commitment_catalog
from ..particles.red_blood_cell_transport import build_red_blood_cell_network
from ..runtime.console_output import print_key_values, print_section, print_stage, print_warning
from ..runtime.progress import create_stage_progress_bar
from ..visualization.vtk.pyvista_flow import render_cfd_flow_fields, validate_cfd_flow_dependencies
from ..visualization.vtk.pyvista_wall_shear import render_wall_shear_visualization
from .molecular_pilot_runner import run_configured_contact_pilot


def _begin_progress_module(progress, module: str) -> None:
    """Show which module is currently running before it performs long work."""

    progress.set_postfix(module=str(module), refresh=True)


def run_generation(cfg, *, render_artifacts=True, reuse_field_from=None):
    """
    Run vascular loading, field generation, accelerated particle advection, and saving.

    ``render_artifacts=False`` keeps every numerical stage enabled while avoiding
    the optional PyVista/VTK dependency check and the final rendering stage. This
    is intended for parameter searches where only the saved numerical artifacts
    are required.

    ``reuse_field_from`` skips the DOLFINx solve after the saved field has passed
    strict configuration, grid, raster, continuous-geometry, convergence, and
    array-contract validation.  Omitting it preserves the original solve path.
    """
    total_stages = 8

    # ============================================================================
    # ======= Initialize output directory and validate visualization dependencies
    # ============================================================================
    if render_artifacts:
        print_stage(1, total_stages, "Initialize output and validate visualization dependencies")
    else:
        print_stage(1, total_stages, "Initialize numerical output", detail="Visualization rendering disabled")

    stage_progress = create_stage_progress_bar(1, 3)
    if render_artifacts:
        _begin_progress_module(stage_progress, "visualization dependencies")
        validate_cfd_flow_dependencies()
        stage_progress.update()
    else:
        _begin_progress_module(stage_progress, "rendering disabled")
        stage_progress.update()

    # Create the output directory and save the resolved run configuration for reproducibility.
    _begin_progress_module(stage_progress, "output directory")
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    stage_progress.update()
    _begin_progress_module(
        stage_progress,
        "run configuration" if cfg.save_run_config else "configuration saving skipped",
    )
    if cfg.save_run_config:
        resolved_run_config                                     = deepcopy(cfg.raw)
        resolved_run_config["random_seed"]                      = cfg.random_seed
        resolved_run_config["_resolved_source_config_path"]     = str(cfg.source_path)
        resolved_run_config["_resolved_field_reuse_contract"]   = build_field_reuse_config_contract(cfg)

        if cfg.molecular_target.mask_npz_path is not None:
            resolved_run_config.setdefault("molecular_target", {})["mask_npz_path"] = str(cfg.molecular_target.mask_npz_path.resolve())

        save_run_config(cfg.output_dir / "run_config.yaml", resolved_run_config)
    stage_progress.update()
    stage_progress.close()

    # ============================================================================
    # ======= Load the existing vessel model and build the X-Z grid domain
    # ============================================================================
    print_stage(2, total_stages, "Load exported vascular model", detail=str(cfg.model_dir))
    stage_progress = create_stage_progress_bar(2, 2)
    _begin_progress_module(stage_progress, "vascular model")
    physics_input = load_physics_input(
        cfg.model_dir,
        planar_extrusion_depth_um=cfg.field.effective_thickness_um,
    )
    validate_physics_input_geometry(
        cfg.model_dir,
        physics_input,
        expected_mode="planar_2d",
    )
    stage_progress.update()

    # compute the red blood cell local haematocrit for each vessel segment
    _begin_progress_module(stage_progress, "RBC transport network")
    red_blood_cell_network  = build_red_blood_cell_network(physics_input.vessels, cfg.red_blood_cell_root_discharge_hematocrit)
    stage_progress.update()
    stage_progress.close()

    if red_blood_cell_network.enabled:
        print_section("RBC-induced drift–diffusion reduced-order transport model")
        print_key_values(
            [
                ("Root discharge haematocrit", f"{red_blood_cell_network.root_discharge_hematocrit:.6g}"),
                ("Effective RBC diameter", f"{red_blood_cell_network.effective_rbc_diameter_um:.6g} um"),
                ("Local discharge-haematocrit range",
                 f"{float(red_blood_cell_network.discharge_hematocrit.min()):.6g}.."
                 f"{float(red_blood_cell_network.discharge_hematocrit.max()):.6g}",
                ),
                ("Local tube-haematocrit range",
                 f"{float(red_blood_cell_network.tube_hematocrit.min()):.6g}.."
                 f"{float(red_blood_cell_network.tube_hematocrit.max()):.6g}",
                ),
                ("Maximum RBC-flux conservation error", f"{red_blood_cell_network.maximum_rbc_flux_conservation_relative_error:.6g}"),
            ]
        )

        # warn if any vessel tube-haematocrit values are outside the RBC model's normal quantitative range
        if not red_blood_cell_network.tube_hematocrit_in_quantitative_range.all():
            print_warning(
                "Some vessel tube-haematocrit values are outside the RBC model's "
                "quantitative range 0.15..0.30. They are retained and flagged "
                "without clipping."
            )

        # warn if any vessel diameters are outside the paper's approximate 3~50 um range
        if not red_blood_cell_network.cfl_diameter_in_range.all():
            print_warning(
                "Some vessel diameters are outside the approximately 3..50 um "
                "range measured in the cortical THG/3PEF CFL study. Values remain "
                "recorded without clipping; RBC "
                "drift and diffusion are nevertheless exactly disabled where "
                "D/d_R <= 1.25."
            )

    # ============================================================================
    # ======= Build the grid domain to calculate how many grids are needed to cover the entire vessel model,
    # ======= and how many grids are needed in the X and Z directions.
    # ============================================================================
    print_stage(3, total_stages, "Build and validate the X-Z lumen domain")
    stage_progress = create_stage_progress_bar(3, 5)

    # Build grid domain for the X-Z plane, which is used for the CFD solver and particle advection.
    _begin_progress_module(stage_progress, "Cartesian domain")
    domain = build_domain_from_vessels(physics_input.vessels, cfg.domain)
    stage_progress.update()

    # Combine and smooth the vessel geometry into a continuous representation.
    # Build precise query indexing of the continuous geometry for each particles in vessel segment
    _begin_progress_module(stage_progress, "continuous geometry")
    continuous_geometry = build_continuous_vessel_geometry(
        physics_input.vessels,
        domain,
        maximum_boundary_element_length_um=(
            cfg.domain.continuous_boundary_maximum_element_length_um
        ),
    )
    stage_progress.update()

    # Build the virtual gate, when particle hasn't cross the gate, it is still in the vessel segment, 
    # when it cross the gate, it is in the next vessel segment
    _begin_progress_module(stage_progress, "topological ownership")
    topological_ownership   = build_topological_commitment_catalog(physics_input.vessels, continuous_geometry)
    stage_progress.update()

    # Build Cartesian far-field-cache support and diagnostic material arrays.
    # This grid never defines particle contact or the physical vessel wall.
    _begin_progress_module(stage_progress, "vessel rasterization")
    raster = rasterize_vessels(physics_input.vessels, domain, cfg.domain,
                               effective_thickness_um=cfg.field.effective_thickness_um,
                               continuous_geometry=continuous_geometry,
                               dynamic_viscosity_mpas=(
                                   cfg.field.kinematic_viscosity_um2_s
                                   * cfg.field.blood_density_kg_m3
                                   * 1.0e-9
                               ))
    stage_progress.update()

    _begin_progress_module(stage_progress, "fluid connectivity")
    validate_fluid_connectivity(
        raster,
        domain,
        continuous_geometry=continuous_geometry,
    )
    stage_progress.update()
    stage_progress.close()
    print_key_values(
        [("Grid shape", f"{domain.shape[0]} x {domain.shape[1]}"),
         ("Grid spacing", f"{domain.spacing_um:.3f} um"),
        ]
    )

    # ============================================================================================
    # ======= Solve or reuse the velocity field and wall-shear proxy
    # ============================================================================================
    field_reused                        = reuse_field_from is not None
    field_reuse_source_npz_path         = None
    field_reuse_source_run_config_path  = None

    if reuse_field_from is None:
        print_stage(4, total_stages, "Solve and diagnose the velocity field")
        stage_progress = create_stage_progress_bar(4, 10)
        completed_solver_modules = 0
        started_solver_modules = 0

        def report_solver_progress(module: str) -> None:
            nonlocal completed_solver_modules, started_solver_modules
            if started_solver_modules:
                stage_progress.update()
                completed_solver_modules += 1
            _begin_progress_module(stage_progress, module)
            started_solver_modules += 1

        def report_solver_subprogress(
            submodule: str,
            completed: int,
            total: int,
        ) -> None:
            stage_progress.set_submodule_progress(
                completed=completed,
                total=total,
                submodule=submodule,
            )

        flow = solve_dolfinx_stokes_gmsh_2d(
            domain, raster, cfg.field, physics_input.vessels, continuous_geometry,
            vessel_metadata=physics_input.vessel_metadata,
            progress_callback=report_solver_progress,
            subprogress_callback=report_solver_subprogress,
        )
        # Test doubles and third-party wrappers may not emit the optional
        # milestones.  A successful return still completes the solver portion.
        if started_solver_modules:
            stage_progress.update()
            completed_solver_modules += 1
        if completed_solver_modules < 9:
            stage_progress.update(9 - completed_solver_modules)
    else:
        print_stage(4, total_stages, "Reuse and diagnose the validated velocity field", detail=str(reuse_field_from))
        stage_progress = create_stage_progress_bar(4, 2)
        _begin_progress_module(stage_progress, "validated field reuse")
        reused  = load_reusable_flow_field(reuse_field_from, cfg=cfg, domain=domain, raster=raster, continuous_geometry=continuous_geometry)
        flow    = reused.flow
        field_reuse_source_npz_path         = reused.field_npz_path
        field_reuse_source_run_config_path  = reused.run_config_path
        stage_progress.update()

    _begin_progress_module(stage_progress, "flow diagnostics")
    write_flow_diagnostics(cfg.output_dir, domain, raster, flow, physics_input.vessels)
    stage_progress.update()
    stage_progress.close()

    # ============================================================================================
    # ======= Prepare the hydrodynamic fields for particle advection and the optional molecular target
    # ============================================================================================
    print_stage(5, total_stages, "Prepare particle hydrodynamics and molecular target")
    stage_progress = create_stage_progress_bar(5, 2)
    _begin_progress_module(stage_progress, "particle hydrodynamics")
    hydrodynamic_fields = build_particle_hydrodynamic_fields(domain, raster, flow, continuous_geometry=continuous_geometry)
    boundary_geometry   = hydrodynamic_fields.boundary_geometry
    stage_progress.update()

    # Check whether this simulation needs a molecular target on the vessel wall.
    _begin_progress_module(
        stage_progress,
        "molecular target" if cfg.molecular_target.enabled else "molecular target disabled",
    )
    molecular_target_field  = (
        build_molecular_target_field(domain, hydrodynamic_fields, cfg.molecular_target)
        if cfg.molecular_target.enabled else None
    )
    stage_progress.update()
    stage_progress.close()

    # ============================================================================================
    # ======= Advect concentration-driven microbubbles through the accepted velocity field.
    # ============================================================================================
    print_stage(6, total_stages, "Advance finite-size microbubbles", detail="Deterministic continuous perfusion with permanent, never-reused IDs")
    trajectories = advect_particles_with_continuous_perfusion(
        domain, raster, flow, physics_input.vessels, cfg.particles, cfg.particle_dynamics,
        hydrodynamic_fields,
        effective_thickness_um=cfg.field.effective_thickness_um,
        boundary_depth_cells=cfg.field.boundary_depth_cells,
        random_seed=cfg.random_seed,
        molecular_binding_cfg=cfg.molecular_binding,
        molecular_target_field=molecular_target_field,
        red_blood_cell_network=red_blood_cell_network,
        topological_ownership=topological_ownership,
    )
    trajectory_metadata = trajectories.metadata

    # ============================================================================================
    # ======= Run the optional molecular contact pilot and save the numerical results.
    # ============================================================================================
    print_stage(7, total_stages, "Run optional contact pilot and save numerical results")
    stage_7_modules = (
        2
        + (1 if molecular_target_field is not None and cfg.binding_scenario_sweep.enabled else 0)
        + (
            2
            + (1 if molecular_target_field is not None else 0)
            + (1 if red_blood_cell_network.enabled else 0)
            if cfg.save_npz
            else 0
        )
    )
    stage_progress = create_stage_progress_bar(7, stage_7_modules)
    if molecular_target_field is not None and cfg.binding_scenario_sweep.enabled:
        _begin_progress_module(stage_progress, "molecular contact pilot")
    molecular_contact_pilot_paths = (
        run_configured_contact_pilot(
            cfg.output_dir, trajectories, molecular_target_field,
            cfg.molecular_binding, cfg.binding_scenario_sweep,
        )
        if molecular_target_field is not None and cfg.binding_scenario_sweep.enabled
        else ()
    )
    if molecular_target_field is not None and cfg.binding_scenario_sweep.enabled:
        stage_progress.update()

    # ============================================================================================
    # ======= Save the accepted field, particle trajectories, and domain metadata.
    # ============================================================================================
    field_path                      = cfg.output_dir / "velocity_and_wall_shear_field.npz"
    trajectories_path               = cfg.output_dir / "microbubble_field_trajectories.npz"
    metadata_path                   = cfg.output_dir / "domain_metadata.yaml"
    molecular_target_path           = cfg.output_dir / "molecular_target_field.npz"
    red_blood_cell_transport_path   = cfg.output_dir / "red_blood_cell_transport.npz"

    if cfg.save_npz:
        _begin_progress_module(stage_progress, "flow field NPZ")
        save_field_npz(field_path, domain, raster, flow, continuous_geometry=continuous_geometry)
        stage_progress.update()
        _begin_progress_module(stage_progress, "trajectories NPZ")
        save_trajectories_npz(trajectories_path, trajectories)
        stage_progress.update()

        if molecular_target_field is not None:
            _begin_progress_module(stage_progress, "molecular target NPZ")
            save_molecular_target_npz(molecular_target_path, molecular_target_field)
            stage_progress.update()
        if red_blood_cell_network.enabled:
            _begin_progress_module(stage_progress, "RBC transport NPZ")
            save_red_blood_cell_transport_npz(red_blood_cell_transport_path, red_blood_cell_network)
            stage_progress.update()
        
    _begin_progress_module(stage_progress, "domain metadata")
    save_domain_metadata(
        metadata_path,
        domain,
        {
            "swc_path":             str(physics_input.swc_path),
            "vessel_data_path":     str(physics_input.vessel_data_path),
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
            "n_vessels":            len(physics_input.vessels),
            "grid_lumen_rasterization": "shapely_exact_cell_intersection_area_v1",
            "grid_lumen_boolean_threshold": 0.5,
            "grid_lumen_fractional_area_um2": float(np.sum(raster.lumen_fraction, dtype=np.float64) * domain.spacing_um**2),
            "continuous_lumen_area_um2": float(continuous_geometry.lumen_polygon.area),
            "inlet_number_concentration_mb_per_ml": cfg.particles.inlet_number_concentration_mb_per_ml,
            "n_steps":              cfg.particles.n_steps,
            "dt_s":                 cfg.particles.dt_s,
            "stored_frames":        cfg.particles.n_steps + 1,
            "particle_output_dt_s": float(trajectory_metadata["output_dt_s"]),
            "particle_integration_substeps": int(trajectory_metadata["integration_substeps"]),
            "particle_internal_dt_s": float(trajectory_metadata["internal_dt_s"]),
            "particle_total_internal_steps": int(trajectory_metadata["total_internal_integration_steps"]),
            "particle_expected_rhs_evaluations": int(trajectory_metadata["expected_rhs_evaluations"]),
            "unique_bubbles_created": int(trajectory_metadata["unique_bubbles_created"]),
            "formally_recorded_bubbles_injected": int(trajectory_metadata["formal_admissions"]),
            "minimum_active_bubbles":   int(trajectory_metadata["minimum_active_bubbles"]),
            "maximum_active_bubbles":   int(trajectory_metadata["maximum_active_bubbles"]),
            "mean_active_bubbles":      float(trajectory_metadata["mean_active_bubbles"]),
            "perfusion_model":          str(trajectory_metadata["perfusion_model"]),
            "injection_rate_bubbles_per_s": float(trajectory_metadata["injection_rate_bubbles_per_s"]),
            "particle_initial_condition": str(trajectory_metadata["initial_condition"]),
            "inlet_wait_events": int(trajectory_metadata["inlet_wait_events"]),
            "maximum_inlet_wait_s": float(trajectory_metadata["maximum_inlet_wait_s"]),
            "particle_motion_model":    "mobility",
            "particle_vessel_ownership_model": topological_ownership.schema,
            "particle_vessel_id_convention": "one_based_persistent_topological_state",
            "particle_vessel_id_raster_authority": False,
            "topological_commitment_section_count": topological_ownership.section_count,
            "topological_commitment_geometry_hash_sha256": topological_ownership.geometry_hash_sha256,
            **particle_boundary_metadata(boundary_geometry),
            "molecular_target_enabled": molecular_target_field is not None,
            "molecular_target_region_mode": (
                str(molecular_target_field.region_mode)
                if molecular_target_field is not None
                else "disabled"
            ),
            "molecular_target_wall_sites": (
                int(molecular_target_field.target_wall_mask.sum())
                if molecular_target_field is not None
                else 0
            ),
            "molecular_target_positive_solid_faces": (
                int(molecular_target_field.wall_target_positive.sum())
                if molecular_target_field is not None
                else 0
            ),
            "molecular_binding_enabled":    cfg.molecular_binding.enabled,
            "molecular_binding_model":      cfg.molecular_binding.model,
            "molecular_target_density_molecules_per_m2": (
                cfg.molecular_target.target_density_molecules_per_m2
            ),
            "molecular_target_npz_path": (
                str(molecular_target_path)
                if molecular_target_field is not None and cfg.save_npz
                else None
            ),
            "molecular_target_source_mask_npz_path": (
                str(molecular_target_field.source_mask_npz_path)
                if molecular_target_field is not None
                and molecular_target_field.source_mask_npz_path is not None
                else None
            ),
            "molecular_contact_pilot_report_count": len(molecular_contact_pilot_paths),
            "molecular_contact_pilot_paths": [
                str(path) for path in molecular_contact_pilot_paths
            ],
            "red_blood_cell_transport_enabled": bool(
                red_blood_cell_network.enabled
            ),
            "red_blood_cell_root_discharge_hematocrit": float(
                red_blood_cell_network.root_discharge_hematocrit
            ),
            "red_blood_cell_effective_diameter_um": float(
                red_blood_cell_network.effective_rbc_diameter_um
            ),
            "red_blood_cell_major_diameter_um": float(
                trajectory_metadata.get(
                    "red_blood_cell_major_diameter_um", 8.0
                )
            ),
            "red_blood_cell_minimum_discharge_hematocrit": float(
                red_blood_cell_network.discharge_hematocrit.min()
            ),
            "red_blood_cell_maximum_discharge_hematocrit": float(
                red_blood_cell_network.discharge_hematocrit.max()
            ),
            "red_blood_cell_minimum_tube_hematocrit": float(
                red_blood_cell_network.tube_hematocrit.min()
            ),
            "red_blood_cell_maximum_tube_hematocrit": float(
                red_blood_cell_network.tube_hematocrit.max()
            ),
            "red_blood_cell_vessels_outside_quantitative_hematocrit_range": int(
                (~red_blood_cell_network.tube_hematocrit_in_quantitative_range).sum()
            ),
            "red_blood_cell_maximum_flux_conservation_relative_error": float(
                red_blood_cell_network.maximum_rbc_flux_conservation_relative_error
            ),
            "red_blood_cell_maximum_speed_um_s": float(
                trajectory_metadata.get("red_blood_cell_maximum_speed_um_s", 0.0)
            ),
            "red_blood_cell_rng_seed": int(
                trajectory_metadata.get(
                    "red_blood_cell_rng_seed", cfg.random_seed
                )
            ),
            "red_blood_cell_rng_algorithm_version": str(
                trajectory_metadata.get(
                    "red_blood_cell_rng_algorithm_version", "disabled"
                )
            ),
            "red_blood_cell_reference_cumulative_shear_strain": float(
                trajectory_metadata.get(
                    "red_blood_cell_reference_cumulative_shear_strain", 450.0
                )
            ),
            "red_blood_cell_reference_transverse_diffusivity_um2_s": float(
                trajectory_metadata.get(
                    "red_blood_cell_reference_transverse_diffusivity_um2_s", 45.0
                )
            ),
            "red_blood_cell_random_wall_reflection_count": int(
                trajectory_metadata.get(
                    "red_blood_cell_random_wall_reflection_count", 0
                )
            ),
            "red_blood_cell_maximum_random_displacement_um": float(
                trajectory_metadata.get(
                    "red_blood_cell_maximum_random_displacement_um", 0.0
                )
            ),
            "red_blood_cell_rms_random_displacement_um": float(
                trajectory_metadata.get(
                    "red_blood_cell_rms_random_displacement_um", 0.0
                )
            ),
            "red_blood_cell_diffusion_enabled_record_fraction": float(
                trajectory_metadata.get(
                    "red_blood_cell_diffusion_enabled_record_fraction", 0.0
                )
            ),
            "red_blood_cell_cfl_record_fraction": float(
                trajectory_metadata.get("red_blood_cell_cfl_record_fraction", 0.0)
            ),
            "red_blood_cell_cfl_reaching_bubble_fraction": float(
                trajectory_metadata.get(
                    "red_blood_cell_cfl_reaching_bubble_fraction", 0.0
                )
            ),
            "red_blood_cell_quantitative_applicability_record_fraction": float(
                trajectory_metadata.get(
                    "red_blood_cell_quantitative_applicability_record_fraction",
                    0.0,
                )
            ),
            "red_blood_cell_exposure_weighted_quantitative_applicability_fraction": float(
                trajectory_metadata.get(
                    "red_blood_cell_exposure_weighted_quantitative_applicability_fraction",
                    0.0,
                )
            ),
            "target_exposure_N_enc": int(
                trajectory_metadata.get("target_exposure_N_enc", 0)
            ),
            "target_exposure_event_count": int(
                trajectory_metadata.get("target_exposure_event_count", 0)
            ),
            "target_exposure_total_time_s": float(
                trajectory_metadata.get("target_exposure_total_time_s", 0.0)
            ),
            "target_reaction_area_time_um2_s": float(
                trajectory_metadata.get(
                    "target_reaction_area_time_um2_s", 0.0
                )
            ),
            "target_exposure_right_censored_event_count": int(
                trajectory_metadata.get(
                    "target_exposure_right_censored_event_count", 0
                )
            ),
            "red_blood_cell_transport_npz_path": (
                str(red_blood_cell_transport_path)
                if red_blood_cell_network.enabled and cfg.save_npz
                else None
            ),
            "particle_time_integrator": str(
                trajectory_metadata["time_integrator"]
            ),
            "particle_trajectory_schema": str(trajectory_metadata["trajectory_schema"]),
            "wall_contact_integrator": str(
                trajectory_metadata["wall_contact_integrator"]
            ),
            "maximum_simultaneous_wall_constraints": int(
                trajectory_metadata["maximum_simultaneous_wall_constraints"]
            ),
            "particle_numeric_kernel_family": str(
                trajectory_metadata["particle_numeric_kernel_family"]
            ),
            "particle_numba_thread_capacity": int(
                trajectory_metadata.get("particle_numba_thread_capacity", 1)
            ),
            "particle_numba_worker_threads": int(
                trajectory_metadata.get("particle_numba_worker_threads", 1)
            ),
            "particle_numba_parallel_swept_path_queries": bool(
                trajectory_metadata.get(
                    "particle_numba_parallel_swept_path_queries", False
                )
            ),
            "particle_numba_parallel_exact_wall_state_queries": bool(
                trajectory_metadata.get(
                    "particle_numba_parallel_exact_wall_state_queries", False
                )
            ),
            "particle_numba_exact_wall_state_parallel_min_batch": int(
                trajectory_metadata.get(
                    "particle_numba_exact_wall_state_parallel_min_batch", 0
                )
            ),
            "particle_numba_directed_inlet_crossing_guard": bool(
                trajectory_metadata.get(
                    "particle_numba_directed_inlet_crossing_guard", False
                )
            ),
            "particle_continuous_wall_inward_normals_cached": bool(
                trajectory_metadata.get(
                    "particle_continuous_wall_inward_normals_cached", False
                )
            ),
            "particle_numba_outlet_spatial_index": str(
                trajectory_metadata.get(
                    "particle_numba_outlet_spatial_index", "not_used"
                )
            ),
            "particle_numba_exact_solid_face_queries": bool(
                trajectory_metadata.get(
                    "particle_numba_exact_solid_face_queries", False
                )
            ),
            "particle_numba_swept_disc_audit": bool(
                trajectory_metadata.get(
                    "particle_numba_swept_disc_audit", False
                )
            ),
            "particle_numba_scalar_diagnostic_reduction": bool(
                trajectory_metadata.get(
                    "particle_numba_scalar_diagnostic_reduction", False
                )
            ),
            "particle_numba_cached_predictive_wall_endpoints": bool(
                trajectory_metadata.get(
                    "particle_numba_cached_predictive_wall_endpoints", False
                )
            ),
            "particle_numba_frame_transaction_batching": bool(
                trajectory_metadata.get(
                    "particle_numba_frame_transaction_batching", False
                )
            ),
            "particle_numba_batched_saved_frame_count": int(
                trajectory_metadata.get(
                    "particle_numba_batched_saved_frame_count", 0
                )
            ),
            "particle_continuous_wall_numba_bin_width_cells": int(
                trajectory_metadata.get(
                    "particle_continuous_wall_numba_bin_width_cells", 0
                )
            ),
            "particle_transport_seconds": float(
                trajectory_metadata.get("particle_transport_seconds", 0.0)
            ),
            "particle_internal_steps_per_wall_second": float(
                trajectory_metadata.get(
                    "particle_internal_steps_per_wall_second", 0.0
                )
            ),
            "particle_velocity_semantics": str(
                trajectory_metadata["particle_velocity_semantics"]
            ),
            "realized_particle_velocity_semantics": str(
                trajectory_metadata["realized_particle_velocity_semantics"]
            ),
            "bubble_diameter_min_um":           float(trajectory_metadata["bubble_diameter_min_um"]),
            "bubble_diameter_max_um":           float(trajectory_metadata["bubble_diameter_max_um"]),
            "bubble_diameter_sample_min_um":    float(trajectory_metadata["bubble_diameter_sample_min_um"]),
            "bubble_diameter_sample_max_um":    float(trajectory_metadata["bubble_diameter_sample_max_um"]),
            "near_wall_xi_min": float(trajectory_metadata["near_wall_xi_min"]),
            "near_wall_xi_near": float(trajectory_metadata["near_wall_xi_near"]),
            "near_wall_xi_far": float(trajectory_metadata["near_wall_xi_far"]),
            "minimum_hydrodynamic_regularization_gap_um": float(
                trajectory_metadata["minimum_hydrodynamic_regularization_gap_um"]
            ),
            "maximum_hydrodynamic_regularization_gap_um": float(
                trajectory_metadata["maximum_hydrodynamic_regularization_gap_um"]
            ),
            "contact_geometry_tolerance_um":   float(
                trajectory_metadata["contact_geometry_tolerance_um"]
            ),
            "contact_max_time_refinements": int(
                trajectory_metadata["contact_max_time_refinements"]
            ),
            "wall_contact_threshold_um":        float(trajectory_metadata["wall_contact_threshold_um"]),
            "minimum_wall_gap_um":              float(trajectory_metadata["minimum_wall_gap_um"]),
            "contact_observations":             int(trajectory_metadata["contact_observations"]),
            "contact_constraint_evaluations": int(
                trajectory_metadata["contact_constraint_evaluations"]
            ),
            "active_contact_constraint_evaluations": int(
                trajectory_metadata["active_contact_constraint_evaluations"]
            ),
            "maximum_contact_reaction_force_pn": float(
                trajectory_metadata["maximum_contact_reaction_force_pn"]
            ),
            "contact_time_refinement_count": int(
                trajectory_metadata["contact_time_refinement_count"]
            ),
            "maximum_contact_time_refinement_depth": int(
                trajectory_metadata["maximum_contact_time_refinement_depth"]
            ),
            "contact_residual_projection_count": int(
                trajectory_metadata["contact_residual_projection_count"]
            ),
            "maximum_contact_residual_projection_um": float(
                trajectory_metadata["maximum_contact_residual_projection_um"]
            ),
            "maximum_contact_complementarity_residual_pn_um": float(
                trajectory_metadata["maximum_contact_complementarity_residual_pn_um"]
            ),
            "minimum_accepted_internal_wall_gap_um": float(
                trajectory_metadata["minimum_accepted_internal_wall_gap_um"]
            ),
            "accepted_negative_gap_count": int(
                trajectory_metadata["accepted_negative_gap_count"]
            ),
            "contact_nonzero_velocity_zero_progress_count": int(
                trajectory_metadata["contact_nonzero_velocity_zero_progress_count"]
            ),
            "contact_kinematic_interval_evaluations": int(
                trajectory_metadata["contact_kinematic_interval_evaluations"]
            ),
            "contact_cumulative_position_path_um": float(
                trajectory_metadata["contact_cumulative_position_path_um"]
            ),
            "contact_cumulative_velocity_path_um": float(
                trajectory_metadata["contact_cumulative_velocity_path_um"]
            ),
            "contact_position_to_velocity_path_ratio": float(
                trajectory_metadata["contact_position_to_velocity_path_ratio"]
            ),
            "minimum_contact_interval_position_to_velocity_path_ratio": float(
                trajectory_metadata[
                    "minimum_contact_interval_position_to_velocity_path_ratio"
                ]
            ),
            "maximum_contact_interval_position_to_velocity_path_ratio": float(
                trajectory_metadata[
                    "maximum_contact_interval_position_to_velocity_path_ratio"
                ]
            ),
            "maximum_free_gap_kinematic_residual_um": float(
                trajectory_metadata["maximum_free_gap_kinematic_residual_um"]
            ),
            "directed_outlet_event_count": int(
                trajectory_metadata["directed_outlet_event_count"]
            ),
            "active_outside_lumen_violations": int(
                trajectory_metadata["active_outside_lumen_violations"]
            ),
            "active_outside_accessible_domain_violations": int(
                trajectory_metadata[
                    "active_outside_accessible_domain_violations"
                ]
            ),
            "discrete_accessibility_disagreement_records": int(
                trajectory_metadata[
                    "discrete_accessibility_disagreement_records"
                ]
            ),
            "field_solver_metadata":    dict(flow.solver_metadata),
            "field_reused":             field_reused,
            "field_reuse_source_npz_path": (
                str(field_reuse_source_npz_path)
                if field_reuse_source_npz_path is not None
                else None
            ),
            "field_reuse_source_run_config_path": (
                str(field_reuse_source_run_config_path)
                if field_reuse_source_run_config_path is not None
                else None
            ),
            "field_reuse_validation": (
                "strict_input_domain_field_grid_raster_continuous_geometry_and_convergence"
                if field_reused
                else "not_requested"
            ),
        },
    )
    stage_progress.update()
    _begin_progress_module(stage_progress, "numerical results complete")
    stage_progress.update()
    stage_progress.close()

    # ===========================================================================================
    # ======= Render the final wall-shear and CFD flow visualization artifacts.
    # ===========================================================================================
    if render_artifacts:
        print_stage(8, total_stages, "Render wall-shear and CFD visualization artifacts")
        stage_progress = create_stage_progress_bar(8, 2)
        _begin_progress_module(stage_progress, "wall-shear visualization")
        wall_shear_artifacts = render_wall_shear_visualization(
            cfg.output_dir,
            domain=domain,
            raster=raster,
            flow=flow,
        )
        stage_progress.update()

        _begin_progress_module(stage_progress, "CFD visualizations")
        initial_flow_html_path, final_flow_html_path = render_cfd_flow_fields(
            domain,
            raster,
            flow,
            cfg.output_dir,
            vessels=physics_input.vessels,
            effective_thickness_um=cfg.field.effective_thickness_um,
            boundary_depth_cells=cfg.field.boundary_depth_cells,
        )
        stage_progress.update()
        stage_progress.close()
        final_wall_shear_html_path = wall_shear_artifacts.html_path
    else:
        print_stage(
            8,
            total_stages,
            "Skip visualization rendering",
            detail="Numerical artifacts are complete",
        )
        stage_progress = create_stage_progress_bar(8, 1)
        _begin_progress_module(stage_progress, "rendering skipped")
        stage_progress.update()
        stage_progress.close()
        initial_flow_html_path = None
        final_flow_html_path = None
        final_wall_shear_html_path = None

    return {
        "output_dir": cfg.output_dir,
        "field_npz_path": field_path,
        "trajectories_npz_path": trajectories_path,
        "metadata_path": metadata_path,
        "molecular_target_npz_path": (
            molecular_target_path
            if molecular_target_field is not None and cfg.save_npz
            else None
        ),
        "molecular_contact_pilot_paths": molecular_contact_pilot_paths,
        "initial_flow_html_path": initial_flow_html_path,
        "final_flow_html_path": final_flow_html_path,
        "final_wall_shear_html_path": final_wall_shear_html_path,
        "red_blood_cell_transport_npz_path": (
            red_blood_cell_transport_path
            if red_blood_cell_network.enabled and cfg.save_npz
            else None
        ),
    }
