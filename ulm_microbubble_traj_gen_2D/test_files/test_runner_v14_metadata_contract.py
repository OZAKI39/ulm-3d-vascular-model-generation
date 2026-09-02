from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import numpy as np
from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.test_files.particle_fixtures import (
    advect_test_particles,
    particle_config,
    straight_channel_case,
)
from ulm_microbubble_traj_gen_2D.utils.core.config import (
    ParticleDynamicsConfig,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.particles import particle_perfusion_transport
from ulm_microbubble_traj_gen_2D.utils.workflows import runner


class RunnerV15MetadataContractTests(unittest.TestCase):
    def test_transport_v15_contract_reaches_domain_metadata(self) -> None:
        domain, raster, flow = straight_channel_case()
        root = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([0.0, 0.0, 3.0]),
            x_d=np.asarray([19.0, 0.0, 3.0]),
            radius=3.0,
            flow_rate=50.0,
        )
        particle_cfg = particle_config(
            n_steps=1,
            dt_s=0.01,
            acceleration_backend="python",
        )
        dynamics_cfg = ParticleDynamicsConfig(
            time_integrator="euler",
            integration_substeps=1,
            near_wall_enabled=False,
            collisions_enabled=False,
            store_full_diagnostics=True,
        )

        with mock.patch.object(
            particle_perfusion_transport,
            "create_particle_progress_bar",
            return_value=mock.Mock(),
        ):
            trajectories = advect_test_particles(
                domain,
                raster,
                flow,
                [root],
                particle_cfg,
                dynamics_cfg,
            )

        continuous_geometry = build_continuous_vessel_geometry([root], domain)
        hydrodynamic_fields = SimpleNamespace(
            boundary_geometry=continuous_geometry
        )

        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)
            cfg = SimpleNamespace(
                source_path=output_dir / "config.yaml",
                raw={},
                random_seed=42,
                model_dir=output_dir,
                output_dir=output_dir,
                save_run_config=False,
                save_npz=False,
                domain=SimpleNamespace(
                    continuous_boundary_maximum_element_length_um=4.0,
                ),
                field=SimpleNamespace(
                    effective_thickness_um=1.0,
                    boundary_depth_cells=1.5,
                    kinematic_viscosity_um2_s=3.0e6,
                    blood_density_kg_m3=1060.0,
                ),
                particles=particle_cfg,
                particle_dynamics=dynamics_cfg,
                red_blood_cell_root_discharge_hematocrit=0.0,
                molecular_target=SimpleNamespace(
                    enabled=False,
                    target_density_molecules_per_m2=0.0,
                ),
                molecular_binding=SimpleNamespace(
                    enabled=False,
                    model="deterministic_mean_field",
                ),
                binding_scenario_sweep=SimpleNamespace(enabled=False),
            )
            physics_input = SimpleNamespace(
                vessels=[root],
                swc_path=output_dir / "input.swc",
                vessel_data_path=output_dir / "input.vessels.npz",
                vessel_metadata={},
            )
            wall_render = SimpleNamespace(
                html_path=output_dir / "wall_shear.html"
            )

            with (
                mock.patch.object(
                    runner, "validate_cfd_flow_dependencies"
                ) as validate_dependencies,
                mock.patch.object(
                    runner, "load_physics_input", return_value=physics_input
                ),
                mock.patch.object(runner, "validate_physics_input_geometry"),
                mock.patch.object(
                    runner, "build_domain_from_vessels", return_value=domain
                ),
                mock.patch.object(
                    runner, "rasterize_vessels", return_value=raster
                ),
                mock.patch.object(runner, "validate_fluid_connectivity"),
                mock.patch.object(
                    runner, "solve_dolfinx_stokes_gmsh_2d", return_value=flow
                ) as solve_velocity,
                mock.patch.object(
                    runner,
                    "load_reusable_flow_field",
                    return_value=SimpleNamespace(
                        flow=flow,
                        field_npz_path=(output_dir / "source" / "velocity_and_wall_shear_field.npz").resolve(),
                        run_config_path=(output_dir / "source" / "run_config.yaml").resolve(),
                    ),
                ) as load_reusable,
                mock.patch.object(runner, "write_flow_diagnostics"),
                mock.patch.object(
                    runner,
                    "build_particle_hydrodynamic_fields",
                    return_value=hydrodynamic_fields,
                ),
                mock.patch.object(
                    runner,
                    "advect_particles_with_continuous_perfusion",
                    return_value=trajectories,
                ) as advect_particles,
                mock.patch.object(runner, "save_domain_metadata") as save_metadata,
                mock.patch.object(
                    runner,
                    "render_wall_shear_visualization",
                    return_value=wall_render,
                ) as render_wall_shear,
                mock.patch.object(
                    runner,
                    "render_cfd_flow_fields",
                    return_value=(
                        output_dir / "initial.html",
                        output_dir / "final.html",
                    ),
                ) as render_flow,
            ):
                rendered_result = runner.run_generation(cfg)
                numerical_result = runner.run_generation(
                    cfg,
                    render_artifacts=False,
                )
                reused_result = runner.run_generation(
                    cfg,
                    render_artifacts=False,
                    reuse_field_from=output_dir / "source",
                )

            validate_dependencies.assert_called_once_with()
            self.assertEqual(solve_velocity.call_count, 2)
            self.assertEqual(advect_particles.call_count, 3)
            for call in advect_particles.call_args_list:
                self.assertNotIn("cardiac", call.kwargs)
            load_reusable.assert_called_once_with(
                output_dir / "source",
                cfg=cfg,
                domain=domain,
                raster=raster,
                continuous_geometry=mock.ANY,
            )
            render_wall_shear.assert_called_once()
            render_flow.assert_called_once()
            self.assertEqual(
                rendered_result["initial_flow_html_path"],
                output_dir / "initial.html",
            )
            self.assertEqual(
                rendered_result["final_flow_html_path"],
                output_dir / "final.html",
            )
            self.assertEqual(
                rendered_result["final_wall_shear_html_path"],
                output_dir / "wall_shear.html",
            )
            self.assertIsNone(numerical_result["initial_flow_html_path"])
            self.assertIsNone(numerical_result["final_flow_html_path"])
            self.assertIsNone(numerical_result["final_wall_shear_html_path"])
            self.assertIsNone(reused_result["initial_flow_html_path"])
            self.assertIsNone(reused_result["final_flow_html_path"])
            self.assertIsNone(reused_result["final_wall_shear_html_path"])

        saved_metadata = save_metadata.call_args.args[2]
        trajectory_metadata = trajectories.metadata

        self.assertTrue(saved_metadata["field_reused"])
        self.assertEqual(
            saved_metadata["grid_lumen_rasterization"],
            "shapely_exact_cell_intersection_area_v1",
        )
        self.assertEqual(
            saved_metadata["grid_lumen_boolean_threshold"],
            0.5,
        )
        self.assertEqual(
            saved_metadata["field_reuse_source_npz_path"],
            str((output_dir / "source" / "velocity_and_wall_shear_field.npz").resolve()),
        )
        self.assertEqual(
            saved_metadata["field_reuse_source_run_config_path"],
            str((output_dir / "source" / "run_config.yaml").resolve()),
        )

        self.assertEqual(
            saved_metadata["particle_trajectory_schema"],
            "continuous_perfusion_revised_v20_topological_records_v13",
        )
        self.assertEqual(
            saved_metadata["wall_contact_integrator"],
            "revised_v16_continuous_geometry_predictive_mobility_unilateral_single_wall",
        )
        self.assertEqual(saved_metadata["maximum_simultaneous_wall_constraints"], 1)
        self.assertEqual(
            saved_metadata["particle_boundary_geometry_schema"],
            "v16_continuous_swept_vessel_boundary",
        )
        self.assertNotIn(
            "particle_boundary_legacy_label_fallback",
            saved_metadata,
        )
        self.assertEqual(
            saved_metadata["particle_true_gap_definition"],
            "g_R=distance_to_continuous_closed_vessel_boundary-radius_um",
        )
        self.assertEqual(
            saved_metadata["particle_true_gap_discretization"],
            "pre_raster_continuous_boundary_brep_with_recorded_curve_tessellation",
        )
        self.assertEqual(
            saved_metadata["particle_hydrodynamic_gap_definition"],
            "max(g_R,radius_um*xi_min)_mobility_coefficients_only",
        )
        self.assertEqual(
            saved_metadata["particle_molecular_capture_gap_definition"],
            "same_unregularized_continuous_g_R_used_by_position_contact",
        )
        self.assertEqual(
            saved_metadata["molecular_target_surface_support"],
            "subset_of_continuous_closed_wall_excluding_anatomical_open_sections",
        )

        for diagnostic_name in (
            "contact_residual_projection_count",
            "maximum_contact_residual_projection_um",
            "maximum_contact_complementarity_residual_pn_um",
            "minimum_accepted_internal_wall_gap_um",
            "accepted_negative_gap_count",
            "contact_nonzero_velocity_zero_progress_count",
            "contact_kinematic_interval_evaluations",
            "contact_cumulative_position_path_um",
            "contact_cumulative_velocity_path_um",
            "contact_position_to_velocity_path_ratio",
            "minimum_contact_interval_position_to_velocity_path_ratio",
            "maximum_contact_interval_position_to_velocity_path_ratio",
            "maximum_free_gap_kinematic_residual_um",
            "directed_outlet_event_count",
            "active_outside_lumen_violations",
            "active_outside_accessible_domain_violations",
            "discrete_accessibility_disagreement_records",
        ):
            expected = trajectory_metadata[diagnostic_name]
            actual = saved_metadata[diagnostic_name]
            if isinstance(expected, float) and np.isnan(expected):
                self.assertTrue(np.isnan(actual), diagnostic_name)
            else:
                self.assertEqual(actual, expected, diagnostic_name)

        for removed_v13_name in (
            "contact_manifold_projection_count",
            "maximum_contact_projection_iterations",
            "contact_projection_failure_attempts",
        ):
            self.assertNotIn(removed_v13_name, saved_metadata)


if __name__ == "__main__":
    unittest.main()
