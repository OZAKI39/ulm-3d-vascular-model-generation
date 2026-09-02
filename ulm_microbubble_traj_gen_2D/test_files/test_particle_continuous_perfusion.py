from __future__ import annotations

from dataclasses import replace
import tempfile
import unittest
from unittest import mock
from pathlib import Path

import numpy as np

from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.utils.core.config import (
    CardiacPulsatilityConfig,
    ParticleDynamicsConfig,
)
from ulm_microbubble_traj_gen_2D.utils.io.field_io import save_trajectories_npz
from ulm_microbubble_traj_gen_2D.utils.particles.particle_perfusion_schedule import (
    radical_inverse,
)
from ulm_microbubble_traj_gen_2D.utils.particles import particle_perfusion_transport
from ulm_microbubble_traj_gen_2D.test_files.particle_fixtures import (
    advect_test_particles as advect_particles,
    particle_config as _particle_config,
    straight_channel_case as _straight_channel_case,
)
from ulm_microbubble_traj_gen_2D.utils.visualization.results.result_loader import _trajectory_arrays


class ContinuousPerfusionMathematicsTests(unittest.TestCase):
    def test_two_dimensional_halton_coordinates_are_exact_and_seed_free(self) -> None:
        expected = (
            (0.5, 1.0 / 3.0),
            (0.25, 2.0 / 3.0),
            (0.75, 1.0 / 9.0),
            (0.125, 4.0 / 9.0),
        )
        actual = tuple(
            (radical_inverse(index, 2), radical_inverse(index, 3))
            for index in range(1, 5)
        )
        np.testing.assert_allclose(actual, expected, rtol=0.0, atol=1.0e-15)


class ContinuousPerfusionTransportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.domain, cls.raster, cls.flow = _straight_channel_case()
        cls.root = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([0.0, 0.0, 3.0]),
            x_d=np.asarray([19.0, 0.0, 3.0]),
            radius=3.0,
            flow_rate=50.0,
        )

    def test_formal_simulation_starts_empty_and_builds_variable_population(self) -> None:
        particle_cfg = replace(
            _particle_config(n_steps=10, dt_s=0.1),
            inlet_number_concentration_mb_per_ml=5.0e10,
            acceleration_backend="python",
        )
        dynamics_cfg = ParticleDynamicsConfig(
            time_integrator="euler",
            integration_substeps=1,
            near_wall_enabled=False,
            collisions_enabled=False,
            store_full_diagnostics=True,
        )
        formal_progress = mock.Mock()

        with mock.patch.object(
            particle_perfusion_transport,
            "create_particle_progress_bar",
            return_value=formal_progress,
        ):
            trajectories = advect_particles(
                self.domain,
                self.raster,
                self.flow,
                [self.root],
                particle_cfg,
                dynamics_cfg,
            )

        self.assertEqual(
            trajectories.metadata["perfusion_model"],
            "deterministic_equal_flux_halton_empty_start_v6",
        )
        self.assertAlmostEqual(
            trajectories.metadata["injection_rate_bubbles_per_s"], 2.375
        )
        self.assertEqual(
            trajectories.metadata["initial_condition"],
            "empty_lumen_at_formal_time_zero",
        )
        self.assertEqual(trajectories.frame_offsets.size, particle_cfg.n_steps + 2)
        np.testing.assert_array_equal(
            np.diff(trajectories.frame_offsets), trajectories.active_count_per_frame
        )
        np.testing.assert_array_equal(
            trajectories.registry_bubble_id,
            np.arange(trajectories.registry_bubble_id.size),
        )
        self.assertEqual(int(trajectories.active_count_per_frame[0]), 0)
        self.assertEqual(int(trajectories.frame_offsets[1]), 0)
        self.assertTrue(np.all(trajectories.registry_scheduled_injection_time_s > 0.0))
        self.assertTrue(np.all(trajectories.registry_admission_time_s > 0.0))
        self.assertFalse(np.any(trajectories.registry_admission_time_s < 0.0))
        self.assertTrue(np.all(trajectories.active))
        self.assertEqual(int(trajectories.metadata["inlet_wait_events"]), 0)
        self.assertFalse(bool(trajectories.metadata["state_storage_slots_reused"]))
        self.assertEqual(
            trajectories.metadata["trajectory_schema"],
            "continuous_perfusion_revised_v20_topological_records_v13",
        )
        self.assertEqual(
            trajectories.metadata["wall_contact_integrator"],
            "revised_v16_continuous_geometry_predictive_mobility_unilateral_single_wall",
        )
        self.assertEqual(
            trajectories.metadata["maximum_simultaneous_wall_constraints"],
            1,
        )
        self.assertEqual(
            trajectories.metadata["particle_velocity_semantics"],
            "last_accepted_internal_v16_continuous_geometry_predictive_velocity_from_background_hydrodynamics_plus_mobility_collision_and_unilateral_contact",
        )
        self.assertIn("contact_constraint_evaluations", trajectories.metadata)
        for diagnostic_name in (
            "contact_residual_projection_count",
            "maximum_contact_residual_projection_um",
            "maximum_contact_complementarity_residual_pn_um",
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
        ):
            self.assertIn(diagnostic_name, trajectories.metadata)
        self.assertIn(
            "minimum_accepted_internal_wall_gap_um",
            trajectories.metadata,
        )
        self.assertEqual(trajectories.metadata["accepted_negative_gap_count"], 0)
        self.assertEqual(
            trajectories.metadata["scheduled_injection_time_semantics"],
            "time relative to empty-lumen formal t=0",
        )
        formal_progress.update.assert_called_with(dynamics_cfg.integration_substeps)
        self.assertEqual(
            formal_progress.update.call_count,
            particle_cfg.n_steps,
        )

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "continuous_perfusion.npz"
            save_trajectories_npz(path, trajectories)
            with np.load(path) as saved:
                self.assertIn("record_realized_velocities_um_s", saved.files)
                self.assertNotIn("record_wall_slide_failure_duration_s", saved.files)
                self.assertIn("record_contact_constraint_active", saved.files)
                self.assertIn("record_contact_reaction_force_pn", saved.files)
                self.assertIn(
                    "record_contact_free_normal_velocity_um_s",
                    saved.files,
                )
                self.assertIn(
                    "record_contact_constrained_normal_velocity_um_s",
                    saved.files,
                )
                self.assertIn("registry_scheduled_injection_time_s", saved.files)
                self.assertIn("registry_admission_time_s", saved.files)
                self.assertIn("registry_exit_time_s", saved.files)
                self.assertIn("registry_inlet_wait_time_s", saved.files)
                display = _trajectory_arrays(saved)

        for bubble in np.unique(display["bubble_id"][display["bubble_id"] >= 0]):
            display_lanes = np.argwhere(display["bubble_id"] == bubble)[:, 1]
            self.assertEqual(np.unique(display_lanes).size, 1)

        for frame in range(particle_cfg.n_steps + 1):
            start, end = trajectories.frame_offsets[frame : frame + 2]
            ids = trajectories.bubble_id[int(start) : int(end)]
            positions = trajectories.positions_um[int(start) : int(end)][:, (0, 2)]
            radii = 0.5 * trajectories.registry_diameter_um[ids]
            for left in range(ids.size):
                for right in range(left + 1, ids.size):
                    distance = float(np.linalg.norm(positions[left] - positions[right]))
                    self.assertGreaterEqual(distance + 1.0e-9, radii[left] + radii[right])

    def test_progress_advances_each_accepted_internal_substep(self) -> None:
        particle_cfg = replace(
            _particle_config(n_steps=2, dt_s=0.1),
            inlet_number_concentration_mb_per_ml=5.0e10,
            acceleration_backend="python",
        )
        dynamics_cfg = ParticleDynamicsConfig(
            time_integrator="euler",
            integration_substeps=3,
            near_wall_enabled=False,
            collisions_enabled=False,
            store_full_diagnostics=False,
        )
        formal_progress = mock.Mock()

        with mock.patch.object(
            particle_perfusion_transport,
            "create_particle_progress_bar",
            return_value=formal_progress,
        ):
            advect_particles(
                self.domain,
                self.raster,
                self.flow,
                [self.root],
                particle_cfg,
                dynamics_cfg,
            )

        self.assertEqual(
            formal_progress.update.call_args_list,
            [mock.call(1)] * 6,
        )
        final_postfix = formal_progress.set_postfix.call_args.kwargs
        self.assertEqual(final_postfix["time"], "0.200/0.200s")
        self.assertEqual(final_postfix["frame"], "2/2")
        self.assertIn("created", final_postfix)
        self.assertIn("active", final_postfix)
        self.assertIn("exited", final_postfix)
        self.assertIn("waiting", final_postfix)
        self.assertFalse(final_postfix["refresh"])

    def test_short_horizon_with_no_injection_event_writes_valid_empty_frames(self) -> None:
        particle_cfg = replace(
            _particle_config(n_steps=2, dt_s=0.1),
            inlet_number_concentration_mb_per_ml=5.0e10,
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
            trajectories = advect_particles(
                self.domain,
                self.raster,
                self.flow,
                [self.root],
                particle_cfg,
                dynamics_cfg,
            )

        np.testing.assert_array_equal(trajectories.frame_offsets, np.zeros(4, dtype=np.int64))
        np.testing.assert_array_equal(
            trajectories.active_count_per_frame, np.zeros(3, dtype=np.int32)
        )
        self.assertEqual(trajectories.registry_bubble_id.size, 0)
        self.assertEqual(trajectories.positions_um.shape, (0, 3))
        self.assertTrue(np.isnan(trajectories.metadata["bubble_diameter_sample_min_um"]))

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "empty_continuous_perfusion.npz"
            save_trajectories_npz(path, trajectories)
            with np.load(path) as saved:
                display = _trajectory_arrays(saved)
        self.assertEqual(display["positions_um"].shape, (3, 0, 3))

    def test_enabled_cardiac_flow_modulates_transport_and_saved_records(self) -> None:
        particle_cfg = replace(
            _particle_config(n_steps=10, dt_s=0.02),
            inlet_number_concentration_mb_per_ml=1.0e11,
            acceleration_backend="python",
        )
        dynamics_cfg = ParticleDynamicsConfig(
            time_integrator="heun",
            integration_substeps=2,
            near_wall_enabled=False,
            collisions_enabled=False,
            store_full_diagnostics=True,
        )
        cardiac_cfg = CardiacPulsatilityConfig(
            enabled=True,
            bpm=300.0,
            pulse_propagation_velocity_um_s=25_000.0,
            waveform_samples_per_cycle=1024,
            preserve_cycle_mean_flow=True,
        )

        with mock.patch.object(
            particle_perfusion_transport,
            "create_particle_progress_bar",
            return_value=mock.Mock(),
        ):
            trajectories = advect_particles(
                self.domain,
                self.raster,
                self.flow,
                [self.root],
                particle_cfg,
                dynamics_cfg,
                cardiac_cfg=cardiac_cfg,
            )

        self.assertTrue(trajectories.metadata["cardiac_pulsatility_enabled"])
        self.assertEqual(
            trajectories.metadata["perfusion_model"],
            "deterministic_pulsatile_flux_halton_empty_start_v7",
        )
        self.assertAlmostEqual(trajectories.metadata["cardiac_period_s"], 0.2)
        self.assertAlmostEqual(trajectories.metadata["cardiac_cycle_mean_multiplier"], 1.0)
        self.assertAlmostEqual(trajectories.metadata["cardiac_modulation_strength"], 1.0)
        self.assertIsNotNone(trajectories.cardiac_multiplier)
        self.assertIsNotNone(trajectories.cardiac_waveform_multiplier)
        self.assertIsNotNone(trajectories.cardiac_path_distance_um)
        self.assertEqual(trajectories.cardiac_multiplier.size, trajectories.bubble_id.size)
        self.assertGreater(float(np.ptp(trajectories.cardiac_multiplier)), 0.0)
        self.assertEqual(int(trajectories.active_count_per_frame[0]), 0)
        np.testing.assert_array_equal(
            trajectories.registry_bubble_id,
            np.arange(trajectories.registry_bubble_id.size),
        )

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cardiac_continuous_perfusion.npz"
            save_trajectories_npz(path, trajectories)
            with np.load(path) as saved:
                self.assertIn("record_cardiac_multiplier", saved.files)
                self.assertIn("cardiac_waveform_time_s", saved.files)
                self.assertIn("cardiac_waveform_multiplier", saved.files)
                self.assertIn("cardiac_path_distance_um", saved.files)
                self.assertIn("cardiac_delay_s", saved.files)

    def test_numba_transport_matches_python_reference_lifecycle(self) -> None:
        base_particle_cfg = _particle_config(
            n_steps=12,
            dt_s=0.02,
            inlet_number_concentration_mb_per_ml=5.0e11,
            bubble_diameter_min_um=0.8,
            bubble_diameter_max_um=1.0,
        )
        dynamics_cfg = ParticleDynamicsConfig(
            time_integrator="euler",
            integration_substeps=4,
            near_wall_enabled=True,
            collisions_enabled=True,
            collision_layer_um=0.05,
            collision_relaxation_time_s=0.01,
            neighbor_search="all_pairs",
            store_full_diagnostics=True,
        )

        def run(backend: str):
            with mock.patch.object(
                particle_perfusion_transport,
                "create_particle_progress_bar",
                return_value=mock.Mock(),
            ):
                return advect_particles(
                    self.domain,
                    self.raster,
                    self.flow,
                    [self.root],
                    replace(base_particle_cfg, acceleration_backend=backend),
                    dynamics_cfg,
                )

        reference = run("python")
        accelerated = run("numba_cpu")

        for field_name in (
            "frame_offsets",
            "bubble_id",
            "active_count_per_frame",
            "injected_count_per_frame",
            "terminated_count_per_frame",
            "registry_bubble_id",
            "birth_frame",
            "death_frame",
            "termination_reason",
        ):
            np.testing.assert_array_equal(
                getattr(accelerated, field_name),
                getattr(reference, field_name),
            )
        for field_name in (
            "positions_um",
            "velocities_um_s",
            "realized_velocities_um_s",
            "wall_gap_um",
            "contact_reaction_force_pn",
            "contact_free_normal_velocity_um_s",
            "contact_constrained_normal_velocity_um_s",
            "angular_velocity_rad_s",
            "rotation_angle_rad",
            "collision_force_xz_pn",
            "registry_admission_time_s",
            "registry_exit_time_s",
        ):
            np.testing.assert_allclose(
                getattr(accelerated, field_name),
                getattr(reference, field_name),
                rtol=5.0e-12,
                atol=5.0e-12,
                equal_nan=True,
            )
        np.testing.assert_array_equal(
            accelerated.contact_constraint_active,
            reference.contact_constraint_active,
        )
        self.assertEqual(
            accelerated.metadata["particle_numeric_kernel_family"],
            "numba_batched_component_kernels_v18",
        )
        self.assertEqual(
            accelerated.metadata["particle_numba_outlet_spatial_index"],
            "uniform_grid_csr",
        )
        self.assertTrue(
            bool(
                accelerated.metadata[
                    "particle_numba_scalar_diagnostic_reduction"
                ]
            )
        )
        self.assertFalse(bool(accelerated.metadata["state_storage_slots_reused"]))


if __name__ == "__main__":
    unittest.main()
