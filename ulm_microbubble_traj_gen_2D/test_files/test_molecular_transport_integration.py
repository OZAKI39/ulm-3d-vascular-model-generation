from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.utils.particles import particle_perfusion_transport
from ulm_microbubble_traj_gen_2D.utils.particles.particle_perfusion_transport import (
    _accepted_surface_extension_source_um_s,
)
from ulm_microbubble_traj_gen_2D.utils.core.config import (
    MolecularBindingConfig,
    MolecularTargetConfig,
    ParticleDynamicsConfig,
)
from ulm_microbubble_traj_gen_2D.utils.io.field_io import save_trajectories_npz
from ulm_microbubble_traj_gen_2D.test_files.particle_fixtures import (
    advect_test_particles as advect_particles,
    particle_config as _particle_config,
    straight_channel_case as _straight_channel_case,
)


class MolecularTransportIntegrationTests(unittest.TestCase):
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

    def test_heun_perfusion_forms_and_saves_continuous_mean_field_bonds(self) -> None:
        particle_cfg = replace(
            _particle_config(
                n_steps=10,
                dt_s=0.1,
                bubble_diameter_min_um=2.0,
                bubble_diameter_max_um=2.0,
            ),
            inlet_number_concentration_mb_per_ml=1.0e11,
            acceleration_backend="auto",
        )
        dynamics_cfg = ParticleDynamicsConfig(
            time_integrator="heun",
            integration_substeps=2,
            near_wall_enabled=True,
            collisions_enabled=False,
            store_full_diagnostics=True,
        )
        target_directory = tempfile.TemporaryDirectory()
        self.addCleanup(target_directory.cleanup)
        target_path = Path(target_directory.name) / "selected_target.npz"
        np.savez(
            target_path,
            x_um=self.domain.x_coordinates_um,
            z_um=self.domain.z_coordinates_um,
            target_mask=np.ones(self.domain.shape, dtype=bool),
        )
        target_cfg = MolecularTargetConfig(
            enabled=True,
            region_mode="mask_npz",
            mask_npz_path=target_path,
            target_density_molecules_per_m2=1.0e16,
        )
        binding_cfg = MolecularBindingConfig(
            enabled=True,
            ligand_density_molecules_per_m2=1.0e16,
            capture_distance_um=1.2,
            rest_length_um=0.05,
            association_rate_m2_per_molecule_s=1.0e-12,
            zero_force_dissociation_rate_s=1.0,
            bond_stiffness_pn_per_um=1.0e-12,
            reactive_compliance_nm=0.01,
            temperature_k=310.0,
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
                molecular_target_cfg=target_cfg,
                molecular_binding_cfg=binding_cfg,
            )

        self.assertEqual(
            trajectories.metadata["trajectory_schema"],
            "continuous_perfusion_revised_v20_topological_molecular_records_v13",
        )
        self.assertTrue(trajectories.metadata["molecular_binding_enabled"])
        self.assertEqual(
            trajectories.metadata["target_internal_exposure_integration"],
            "accepted_deterministic_transaction_endpoints_plus_stochastic_geometry_event_vertices_and_final_endpoint_batched_bisection_v2_plus_eight_point_gauss_legendre_over_exposed_subintervals_with_one_shared_physical_dt",
        )
        self.assertIn(
            trajectories.metadata["acceleration_backend"],
            {"numba_cpu_synchronous_mobility", "python_synchronous_mobility"},
        )
        self.assertEqual(
            trajectories.metadata["molecular_association_rate_m2_per_molecule_s"],
            binding_cfg.association_rate_m2_per_molecule_s,
        )
        self.assertEqual(
            trajectories.metadata["molecular_capture_distance_um"],
            binding_cfg.capture_distance_um,
        )
        self.assertIsNotNone(trajectories.bond_count_expected)
        self.assertEqual(trajectories.bond_count_expected.dtype, np.float64)
        self.assertEqual(trajectories.bond_count_expected.size, trajectories.bubble_id.size)
        self.assertGreater(float(np.max(trajectories.target_reaction_area_um2)), 0.0)
        self.assertGreater(float(np.max(trajectories.bond_count_expected)), 0.0)
        self.assertTrue(np.all(trajectories.bond_count_expected >= 0.0))
        self.assertTrue(np.all(np.isfinite(trajectories.bond_force_xz_pn)))
        self.assertGreater(
            float(np.max(np.linalg.norm(trajectories.bond_force_xz_pn, axis=1))),
            0.0,
        )
        self.assertGreater(
            trajectories.metadata["molecular_maximum_bond_torque_pn_um"],
            0.0,
        )
        self.assertGreater(
            trajectories.metadata[
                "molecular_capacity_limited_accepted_step_observations"
            ],
            0,
        )
        self.assertEqual(
            trajectories.registry_final_bond_count_expected.size,
            trajectories.registry_bubble_id.size,
        )
        self.assertIsNotNone(trajectories.registry_target_exposure_time_s)
        self.assertIsNotNone(trajectories.registry_target_exposure_event_count)
        self.assertIsNotNone(
            trajectories.registry_target_reaction_area_time_um2_s
        )
        self.assertGreater(
            float(np.sum(trajectories.registry_target_exposure_time_s)), 0.0
        )
        self.assertGreater(
            float(
                np.sum(
                    trajectories.registry_target_reaction_area_time_um2_s
                )
            ),
            0.0,
        )
        self.assertGreater(
            int(np.sum(trajectories.registry_target_exposure_event_count)), 0
        )
        self.assertEqual(
            trajectories.metadata["target_exposure_N_enc"], 1
        )

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "molecular_trajectories.npz"
            save_trajectories_npz(path, trajectories)
            with np.load(path, allow_pickle=False) as saved:
                self.assertIn("record_bond_count_expected", saved.files)
                self.assertIn("record_bond_formation_rate_bonds_s", saved.files)
                self.assertIn("record_target_reaction_area_um2", saved.files)
                self.assertIn("record_bond_force_xz_pn", saved.files)
                self.assertIn("registry_final_bond_count_expected", saved.files)
                self.assertIn("registry_target_exposure_time_s", saved.files)
                self.assertIn("registry_target_exposure_event_count", saved.files)
                self.assertIn(
                    "registry_target_reaction_area_time_um2_s", saved.files
                )

    def test_projected_hold_uses_realized_motion_for_extension_source(self) -> None:
        source = _accepted_surface_extension_source_um_s(
            start_positions_grid=np.asarray([[1.0, 2.0]]),
            accepted_positions_grid=np.asarray([[1.0, 2.0]]),
            start_angles_rad=np.asarray([0.2]),
            accepted_angles_rad=np.asarray([0.2]),
            radii_um=np.asarray([1.0]),
            representative_bond_count=np.asarray([3.0]),
            start_wall_normal_xz=np.asarray([[0.0, 1.0]]),
            end_wall_normal_xz=None,
            spacing_um=2.0,
            dt_s=0.1,
        )
        np.testing.assert_array_equal(source, np.zeros(1))


if __name__ == "__main__":
    unittest.main()
