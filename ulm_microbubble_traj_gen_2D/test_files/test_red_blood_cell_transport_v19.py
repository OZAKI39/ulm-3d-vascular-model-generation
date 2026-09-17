from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np

from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.utils.core.config import DEFAULT_CONFIG_PATH, load_config
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_contact_pilot import (
    analyze_molecular_contact_pilot,
    contact_pilot_report_mapping,
)
from ulm_microbubble_traj_gen_2D.utils.particles.red_blood_cell_transport import (
    CFL_PLATEAU_DIAMETER_UM,
    CFL_PLATEAU_THICKNESS_UM,
    RBC_MAJOR_DIAMETER_UM,
    RBC_MARGINATION_STRAIN,
    REFERENCE_SHEAR_RATE_S_INV,
    REFERENCE_TRANSVERSE_DIFFUSIVITY_UM2_S,
    REFERENCE_TUBE_HEMATOCRIT,
    build_red_blood_cell_network,
    cfl_width_um,
    evaluate_red_blood_cell_transport,
    fahraeus_tube_hematocrit,
    phase_separation_fraction,
    scale_activation,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_mobility_transport import (
    ParticleVesselOwnershipError,
    _evaluate_rhs,
)
from ulm_microbubble_traj_gen_2D.utils.particles import particle_perfusion_transport
from ulm_microbubble_traj_gen_2D.utils.core.types import GridDomain
from ulm_microbubble_traj_gen_2D.test_files.particle_fixtures import (
    advect_test_particles,
    particle_config,
    straight_channel_case,
)
from ulm_microbubble_traj_gen_2D.test_files.hybrid_velocity_fixture import (
    rectangular_hybrid_velocity,
)
from ulm_microbubble_traj_gen_2D.test_files.test_particle_mobility_transport import (
    _dynamics_config,
    _evaluation_context,
)


class RedBloodCellNetworkV19Tests(unittest.TestCase):
    def test_default_physics_configuration_enables_formal_rbc_transport(self) -> None:
        for quick_test in (False, True):
            with self.subTest(quick_test=quick_test):
                cfg = load_config(DEFAULT_CONFIG_PATH, quick_test=quick_test)

                self.assertEqual(
                    cfg.red_blood_cell_root_discharge_hematocrit,
                    0.35,
                )
                self.assertEqual(
                    set(cfg.raw["red_blood_cell_transport"]),
                    {"root_discharge_hematocrit"},
                )

    def test_formal_configuration_uses_only_hd0_035_and_fixed_8_um_rbc(self) -> None:
        config_path = (
            Path(__file__).resolve().parents[1]
            / "configs"
            / "molecular_contact_fixed_target_20s.yaml"
        )
        cfg = load_config(config_path)

        self.assertEqual(
            cfg.red_blood_cell_root_discharge_hematocrit, 0.35
        )
        self.assertEqual(RBC_MAJOR_DIAMETER_UM, 8.0)
        self.assertEqual(
            set(cfg.raw["red_blood_cell_transport"]),
            {"root_discharge_hematocrit"},
        )

    def test_symmetric_split_preserves_hd_and_rbc_flux(self) -> None:
        vessels = [
            _vessel(0, -1, [1, 2], 100.0, 15.0),
            _vessel(1, 0, [], 50.0, 10.0),
            _vessel(2, 0, [], 50.0, 10.0),
        ]
        network = build_red_blood_cell_network(vessels, 0.35)

        np.testing.assert_allclose(network.discharge_hematocrit, 0.35)
        self.assertAlmostEqual(network.phase_flow_fraction_child_1[0], 0.5)
        self.assertAlmostEqual(network.phase_rbc_flux_fraction_child_1[0], 0.5)
        self.assertLessEqual(
            network.maximum_rbc_flux_conservation_relative_error, 1.0e-15
        )

    def test_multilevel_split_uses_each_local_parent_hematocrit(self) -> None:
        vessels = [
            _vessel(0, -1, [1, 2], 100.0, 15.0),
            _vessel(1, 0, [3, 4], 60.0, 10.0),
            _vessel(2, 0, [], 40.0, 9.0),
            _vessel(3, 1, [], 35.0, 8.0),
            _vessel(4, 1, [], 25.0, 7.0),
        ]
        network = build_red_blood_cell_network(vessels, 0.35)
        row = {int(vessel_id): index for index, vessel_id in enumerate(network.vessel_id)}

        self.assertNotAlmostEqual(network.discharge_hematocrit[row[1]], 0.35)
        q_parent = network.flow_rate_um3_s[row[1]]
        expected_flux = q_parent * network.discharge_hematocrit[row[1]]
        child_flux = sum(
            network.flow_rate_um3_s[row[child]]
            * network.discharge_hematocrit[row[child]]
            for child in (3, 4)
        )
        self.assertAlmostEqual(expected_flux, child_flux, places=12)

    def test_child_swap_keeps_physical_branch_values(self) -> None:
        first = build_red_blood_cell_network(
            [
                _vessel(0, -1, [1, 2], 100.0, 15.0),
                _vessel(1, 0, [], 70.0, 11.0),
                _vessel(2, 0, [], 30.0, 8.0),
            ],
            0.35,
        )
        swapped = build_red_blood_cell_network(
            [
                _vessel(0, -1, [1, 2], 100.0, 15.0),
                _vessel(1, 0, [], 30.0, 8.0),
                _vessel(2, 0, [], 70.0, 11.0),
            ],
            0.35,
        )

        self.assertAlmostEqual(
            first.discharge_hematocrit[1], swapped.discharge_hematocrit[2]
        )
        self.assertAlmostEqual(
            first.discharge_hematocrit[2], swapped.discharge_hematocrit[1]
        )

    def test_phase_separation_uses_pries_1990_parameters(self) -> None:
        q_parent = 100.0
        q_1 = 30.0
        hd = 0.35
        d_parent = 20.0
        d_1 = 10.0
        d_2 = 15.0

        fq, fe, coefficient_a, coefficient_b, cutoff_x0 = (
            phase_separation_fraction(
                q_parent,
                q_1,
                hd,
                d_parent,
                d_1,
                d_2,
            )
        )

        expected_a = -6.96 * np.log(d_1 / d_2) / d_parent
        expected_b = 1.0 + 6.98 * (1.0 - hd) / d_parent
        expected_x0 = 0.4 / d_parent
        expected_z = (fq - expected_x0) / (1.0 - 2.0 * expected_x0)
        expected_fe = 1.0 / (
            1.0
            + np.exp(
                -expected_a
                - expected_b * np.log(expected_z / (1.0 - expected_z))
            )
        )

        self.assertAlmostEqual(fq, q_1 / q_parent)
        self.assertAlmostEqual(coefficient_a, expected_a)
        self.assertAlmostEqual(coefficient_b, expected_b)
        self.assertAlmostEqual(cutoff_x0, expected_x0)
        self.assertAlmostEqual(fe, expected_fe)

    def test_fahraeus_cfl_and_scale_relations_match_declared_values(self) -> None:
        tube = fahraeus_tube_hematocrit(0.35, 20.0)
        self.assertGreater(float(tube), 0.0)
        self.assertLess(float(tube), 0.35)
        np.testing.assert_allclose(
            cfl_width_um(np.asarray([4.0, 8.0, 10.0, 50.0])),
            np.asarray([0.5, 1.0, 1.0, 1.0]),
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            scale_activation(np.asarray([10.0, 11.0, 12.0])),
            np.asarray([0.0, 0.5, 1.0]),
        )

    def test_cfl_rejects_invalid_vessel_diameters(self) -> None:
        for diameter in (0.0, -1.0, np.nan, np.inf):
            with self.subTest(diameter=diameter):
                with self.assertRaises(ValueError):
                    cfl_width_um(diameter)

    def test_cfl_and_scale_activation_are_precomputed_per_vessel(self) -> None:
        network = build_red_blood_cell_network(
            [_vessel(0, -1, [], 100.0, 5.5)], 0.35
        )
        np.testing.assert_array_equal(network.cfl_width_um, cfl_width_um([11.0]))
        np.testing.assert_array_equal(
            network.scale_activation, scale_activation([11.0])
        )
        self.assertEqual(
            network.dense_cfl_width_um_by_vessel_id[0], network.cfl_width_um[0]
        )
        self.assertEqual(
            network.dense_scale_activation_by_vessel_id[0],
            network.scale_activation[0],
        )
        payload = network.to_npz_payload()
        np.testing.assert_array_equal(payload["cfl_width_um"], network.cfl_width_um)
        np.testing.assert_array_equal(
            payload["scale_activation"], network.scale_activation
        )
        self.assertEqual(payload["reference_tube_hematocrit"][0], 0.24)
        self.assertEqual(payload["reference_shear_rate_s_inv"][0], 1000.0)
        self.assertEqual(
            payload["cfl_model"][0],
            "thg_3pef_cortical_mouse_linear_to_plateau_v1",
        )
        self.assertEqual(
            payload["cfl_plateau_diameter_um"][0],
            CFL_PLATEAU_DIAMETER_UM,
        )
        self.assertEqual(
            payload["cfl_plateau_thickness_um"][0],
            CFL_PLATEAU_THICKNESS_UM,
        )

    def test_scale_activation_has_exact_endpoints_and_zero_endpoint_slopes(self) -> None:
        threshold_low_um = 1.25 * RBC_MAJOR_DIAMETER_UM
        threshold_high_um = 1.5 * RBC_MAJOR_DIAMETER_UM
        self.assertEqual(float(scale_activation(threshold_low_um)), 0.0)
        self.assertEqual(float(scale_activation(threshold_high_um)), 1.0)

        step_um = 1.0e-6
        low_one_sided_slope = (
            float(scale_activation(threshold_low_um + step_um))
            - float(scale_activation(threshold_low_um))
        ) / step_um
        high_one_sided_slope = (
            float(scale_activation(threshold_high_um))
            - float(scale_activation(threshold_high_um - step_um))
        ) / step_um
        self.assertAlmostEqual(low_one_sided_slope, 0.0, delta=1.0e-5)
        self.assertAlmostEqual(high_one_sided_slope, 0.0, delta=1.0e-5)

    def test_zero_hd0_marks_network_disabled(self) -> None:
        network = build_red_blood_cell_network(
            [_vessel(0, -1, [], 100.0, 3.0)], 0.0
        )

        self.assertFalse(network.enabled)
        np.testing.assert_array_equal(network.discharge_hematocrit, 0.0)
        np.testing.assert_array_equal(network.tube_hematocrit, 0.0)


class RedBloodCellParticleV19Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.network = build_red_blood_cell_network(
            [_vessel(0, -1, [], 100.0, 10.0)], 0.35
        )

    def test_reference_state_has_tapered_diffusivity_drift_and_fick_velocity(self) -> None:
        evaluation = _evaluate_particle(
            self.network,
            gap_um=5.0,
            radius_um=1.0,
            unique=True,
            use_numba=False,
        )

        self.assertAlmostEqual(evaluation.shear_rate_s_inv[0], 1000.0)
        expected_core_diffusivity = (
            REFERENCE_TRANSVERSE_DIFFUSIVITY_UM2_S
            * evaluation.tube_hematocrit[0]
            / REFERENCE_TUBE_HEMATOCRIT
            * evaluation.shear_rate_s_inv[0]
            / REFERENCE_SHEAR_RATE_S_INV
            * evaluation.scale_activation[0]
        )
        expected_r = (
            5.0 - evaluation.target_gap_um[0]
        ) / (
            evaluation.center_gap_um[0] - evaluation.target_gap_um[0]
        )
        expected_q = 3.0 * expected_r**2 - 2.0 * expected_r**3
        self.assertAlmostEqual(
            evaluation.core_transverse_diffusivity_um2_s[0],
            expected_core_diffusivity,
        )
        self.assertAlmostEqual(
            evaluation.transverse_diffusivity_um2_s[0],
            expected_core_diffusivity * expected_q,
        )
        expected_drift_z = -(
            evaluation.scale_activation[0]
            * evaluation.shear_rate_s_inv[0]
            / RBC_MARGINATION_STRAIN
            * (5.0 - evaluation.target_gap_um[0])
        )
        expected_fick_z = (
            expected_core_diffusivity
            * 6.0
            * expected_r
            * (1.0 - expected_r)
            / (evaluation.center_gap_um[0] - evaluation.target_gap_um[0])
        )
        np.testing.assert_allclose(
            evaluation.drift_velocity_xz_um_s[0], [0.0, expected_drift_z]
        )
        np.testing.assert_allclose(
            evaluation.fick_velocity_xz_um_s[0], [0.0, expected_fick_z]
        )
        np.testing.assert_allclose(
            evaluation.velocity_xz_um_s[0], [0.0, expected_drift_z + expected_fick_z]
        )
        self.assertTrue(evaluation.quantitative_applicability[0])
        self.assertTrue(evaluation.diffusion_enabled[0])

    def test_conditional_mean_drift_converges_to_declared_exponential(self) -> None:
        initial_gap_um = 5.0
        gap_um = initial_gap_um
        elapsed_s = 0.1
        step_count = 5000
        dt_s = elapsed_s / step_count
        first = _evaluate_particle(
            self.network,
            gap_um=gap_um,
            radius_um=1.0,
            unique=True,
            use_numba=False,
        )
        target_gap_um = float(first.target_gap_um[0])
        decay_rate_s_inv = float(
            first.scale_activation[0]
            * first.shear_rate_s_inv[0]
            / RBC_MARGINATION_STRAIN
        )
        for _ in range(step_count):
            evaluation = _evaluate_particle(
                self.network,
                gap_um=gap_um,
                radius_um=1.0,
                unique=True,
                use_numba=False,
            )
            gap_um += dt_s * float(
                np.dot(
                    evaluation.drift_velocity_xz_um_s[0],
                    np.asarray([0.0, 1.0]),
                )
            )

        expected_gap_um = target_gap_um + (
            initial_gap_um - target_gap_um
        ) * np.exp(-decay_rate_s_inv * elapsed_s)
        self.assertAlmostEqual(gap_um, expected_gap_um, delta=3.0e-5)

    def test_target_gap_suppresses_all_transport_and_nonunique_wall_only_deterministic_terms(self) -> None:
        target_gap = float(cfl_width_um(20.0) - 1.0)
        at_target = _evaluate_particle(
            self.network,
            gap_um=target_gap,
            radius_um=0.5,
            unique=True,
            use_numba=False,
        )
        nonunique = _evaluate_particle(
            self.network,
            gap_um=5.0,
            radius_um=1.0,
            unique=False,
            use_numba=False,
        )

        np.testing.assert_array_equal(at_target.velocity_xz_um_s, 0.0)
        np.testing.assert_array_equal(at_target.drift_velocity_xz_um_s, 0.0)
        np.testing.assert_array_equal(at_target.fick_velocity_xz_um_s, 0.0)
        np.testing.assert_array_equal(at_target.transverse_diffusivity_um2_s, 0.0)
        np.testing.assert_array_equal(nonunique.velocity_xz_um_s, 0.0)
        np.testing.assert_array_equal(nonunique.drift_velocity_xz_um_s, 0.0)
        np.testing.assert_array_equal(nonunique.fick_velocity_xz_um_s, 0.0)
        self.assertGreater(nonunique.transverse_diffusivity_um2_s[0], 0.0)

    def test_capillary_activation_closes_drift_diffusion_at_chi_1_25(self) -> None:
        capillary_network = build_red_blood_cell_network(
            [_vessel(0, -1, [], 100.0, 4.0)], 0.35
        )
        evaluation = _evaluate_particle(
            capillary_network,
            gap_um=2.0,
            radius_um=0.5,
            unique=True,
            use_numba=False,
        )

        self.assertEqual(evaluation.scale_activation[0], 0.0)
        self.assertTrue(capillary_network.cfl_diameter_in_range[0])
        self.assertTrue(evaluation.quantitative_applicability[0])
        np.testing.assert_array_equal(evaluation.drift_velocity_xz_um_s, 0.0)
        np.testing.assert_array_equal(evaluation.fick_velocity_xz_um_s, 0.0)
        np.testing.assert_array_equal(evaluation.transverse_diffusivity_um2_s, 0.0)
        self.assertFalse(evaluation.diffusion_enabled[0])

    def test_diffusion_taper_endpoints_and_fick_centreline_are_exact(self) -> None:
        target_gap = max(float(cfl_width_um(20.0)) - 2.0, 0.0)
        center_gap = 10.0 - 1.0
        edge = _evaluate_particle(
            self.network,
            gap_um=target_gap,
            radius_um=1.0,
            unique=True,
            use_numba=False,
        )
        centre = _evaluate_particle(
            self.network,
            gap_um=center_gap,
            radius_um=1.0,
            unique=False,
            use_numba=False,
        )

        self.assertEqual(edge.diffusion_taper_coordinate[0], 0.0)
        self.assertEqual(edge.diffusion_taper[0], 0.0)
        self.assertEqual(edge.transverse_diffusivity_um2_s[0], 0.0)
        self.assertEqual(centre.diffusion_taper_coordinate[0], 1.0)
        self.assertEqual(centre.diffusion_taper[0], 1.0)
        self.assertEqual(
            centre.transverse_diffusivity_um2_s[0],
            centre.core_transverse_diffusivity_um2_s[0],
        )
        np.testing.assert_array_equal(centre.fick_velocity_xz_um_s, 0.0)

        step = 1.0e-7
        q_near_zero = 3.0 * step**2 - 2.0 * step**3
        q_near_one = 3.0 * (1.0 - step) ** 2 - 2.0 * (1.0 - step) ** 3
        self.assertAlmostEqual(q_near_zero / step, 0.0, delta=1.0e-6)
        self.assertAlmostEqual((1.0 - q_near_one) / step, 0.0, delta=1.0e-6)

    def test_invalid_transverse_space_is_diagnosed_and_inert(self) -> None:
        evaluation = _evaluate_particle(
            self.network,
            gap_um=0.0,
            radius_um=10.0,
            unique=True,
            use_numba=False,
        )

        self.assertFalse(evaluation.transverse_space_valid[0])
        np.testing.assert_array_equal(evaluation.velocity_xz_um_s, 0.0)
        np.testing.assert_array_equal(evaluation.transverse_diffusivity_um2_s, 0.0)

    def test_numba_and_python_particle_models_match(self) -> None:
        reference = _evaluate_particle(
            self.network,
            gap_um=5.0,
            radius_um=1.0,
            unique=True,
            use_numba=False,
        )
        accelerated = _evaluate_particle(
            self.network,
            gap_um=5.0,
            radius_um=1.0,
            unique=True,
            use_numba=True,
        )

        for name in reference.__dataclass_fields__:
            np.testing.assert_allclose(
                getattr(accelerated, name),
                getattr(reference, name),
                rtol=0.0,
                atol=0.0,
            )

    def test_continuous_geometry_flags_centreline_as_nonunique(self) -> None:
        domain = GridDomain(
            origin_um=np.asarray([-20.0, 0.0, -20.0]),
            spacing_um=1.0,
            shape=(141, 41),
            fixed_y_um=0.0,
            x_coordinates_um=np.arange(-20.0, 121.0),
            z_coordinates_um=np.arange(-20.0, 21.0),
        )
        vessel = _vessel(0, -1, [], 100.0, 10.0)
        geometry = build_continuous_vessel_geometry([vessel], domain)

        reference = geometry.exact_solid_wall_state_xz_um(
            np.asarray([[50.0, 0.0], [50.0, 8.0]])
        )
        accelerated = geometry.exact_solid_wall_state_xz_um_accelerated(
            np.asarray([[50.0, 0.0], [50.0, 8.0]])
        )

        np.testing.assert_array_equal(
            reference.unique_nearest_wall, np.asarray([False, True])
        )
        np.testing.assert_array_equal(
            accelerated.unique_nearest_wall, reference.unique_nearest_wall
        )

    def test_primary_target_area_time_exposure_is_reported(self) -> None:
        pilot = analyze_molecular_contact_pilot(
            frame_offsets=np.asarray([0, 1, 2], dtype=np.int64),
            bubble_id=np.asarray([7, 7], dtype=np.int64),
            reaction_area_um2=np.asarray([2.0, 3.0]),
            tangential_slip_um_s=np.asarray([0.0, 0.0]),
            output_dt_s=0.1,
        )
        report = contact_pilot_report_mapping(pilot)

        self.assertAlmostEqual(pilot.target_area_time_exposure_um2_s, 0.5)
        self.assertTrue(report["pilot"]["molecular_comparison_ready"])

    def test_rhs_adds_only_rbc_translation(self) -> None:
        domain, raster, flow = straight_channel_case()
        context = _evaluation_context(
            domain,
            raster,
            flow,
            _dynamics_config(near_wall_enabled=False, collisions_enabled=False),
            radii_um=np.asarray([1.0]),
        )
        gradient = np.zeros_like(context.velocity_gradient_s_inv)
        gradient[..., 0, 1] = 1000.0
        baseline_context = replace(context, velocity_gradient_s_inv=gradient)
        rbc_context = replace(
            baseline_context, red_blood_cell_network=self.network
        )
        arguments = (
            np.asarray([[5.0, 4.0]]),
            np.asarray([0], dtype=np.int64),
            np.asarray([True]),
        )

        owner = np.asarray([1], dtype=np.int32)
        baseline = _evaluate_rhs(
            *arguments, baseline_context, topological_vessel_id=owner
        )
        corrected = _evaluate_rhs(
            *arguments, rbc_context, topological_vessel_id=owner
        )

        np.testing.assert_allclose(
            corrected.particle_velocity_xz_um_s
            - baseline.particle_velocity_xz_um_s,
            corrected.red_blood_cell_velocity_xz_um_s,
        )
        np.testing.assert_array_equal(
            corrected.angular_velocity_rad_s, baseline.angular_velocity_rad_s
        )
        np.testing.assert_array_equal(
            corrected.generalized_mobility, baseline.generalized_mobility
        )

    def test_persistent_owner_is_independent_of_missing_raster_label(self) -> None:
        domain, raster, flow = straight_channel_case()
        context = _evaluation_context(
            domain,
            raster,
            flow,
            _dynamics_config(near_wall_enabled=False, collisions_enabled=False),
            radii_um=np.ones(1, dtype=np.float64),
        )
        ownerless_context = replace(
            context,
            vessel_id=np.full_like(context.vessel_id, -1),
            red_blood_cell_network=self.network,
        )

        evaluation = _evaluate_rhs(
            np.asarray([[5.25, 4.0]], dtype=np.float64),
            np.asarray([0], dtype=np.int64),
            np.asarray([True]),
            ownerless_context,
            time_s=6.7683,
            topological_vessel_id=np.asarray([1], dtype=np.int32),
        )

        np.testing.assert_array_equal(
            evaluation.topological_vessel_id, np.asarray([1], dtype=np.int32)
        )
        self.assertIsNotNone(evaluation.red_blood_cell_velocity_xz_um_s)

    def test_invalid_persistent_owner_reports_particle_id_time_and_position(self) -> None:
        domain, raster, flow = straight_channel_case()
        context = _evaluation_context(
            domain,
            raster,
            flow,
            _dynamics_config(near_wall_enabled=False, collisions_enabled=False),
            radii_um=np.ones(48, dtype=np.float64),
        )
        ownerless_context = replace(
            context,
            vessel_id=np.full_like(context.vessel_id, -1),
            red_blood_cell_network=self.network,
        )

        with self.assertRaises(ParticleVesselOwnershipError) as caught:
            _evaluate_rhs(
                np.asarray([[5.25, 20.0]], dtype=np.float64),
                np.asarray([47], dtype=np.int64),
                np.asarray([True]),
                ownerless_context,
                time_s=6.7683,
                topological_vessel_id=np.asarray([0], dtype=np.int32),
            )

        message = str(caught.exception)
        self.assertIn("permanent_microbubble_id=47", message)
        self.assertIn("physical_time_s=6.7683", message)
        self.assertIn("position_grid=[5.25, 20]", message)
        self.assertIn("position_xz_um=[5.25, 20]", message)
        self.assertEqual(caught.exception.failures[0].permanent_microbubble_id, 47)

    def test_zero_hd0_is_bitwise_trajectory_baseline(self) -> None:
        domain, raster, flow = straight_channel_case()
        vessel = _channel_vessel(radius_um=3.0)
        disabled = build_red_blood_cell_network([vessel], 0.0)
        particle_cfg = particle_config(
            n_steps=5,
            dt_s=0.02,
            inlet_number_concentration_mb_per_ml=5.0e11,
            acceleration_backend="python",
        )
        dynamics = _dynamics_config(
            integration_substeps=2,
            near_wall_enabled=False,
            collisions_enabled=False,
        )
        with mock.patch.object(
            particle_perfusion_transport,
            "create_particle_progress_bar",
            return_value=mock.Mock(),
        ):
            baseline = advect_test_particles(
                domain, raster, flow, [vessel], particle_cfg, dynamics
            )
            zero_hd = advect_test_particles(
                domain,
                raster,
                flow,
                [vessel],
                particle_cfg,
                dynamics,
                red_blood_cell_network=disabled,
            )

        for name in (
            "frame_offsets",
            "bubble_id",
            "positions_um",
            "velocities_um_s",
            "wall_gap_um",
            "wall_normal_xz",
            "birth_frame",
            "death_frame",
            "termination_reason",
        ):
            np.testing.assert_array_equal(getattr(zero_hd, name), getattr(baseline, name))
        self.assertIsNone(zero_hd.red_blood_cell_velocity_xz_um_s)

    def test_enabled_perfusion_saves_rbc_record_diagnostics(self) -> None:
        domain, raster, flow = straight_channel_case()
        vessel = _channel_vessel(radius_um=10.0)
        network = build_red_blood_cell_network([vessel], 0.35)
        particle_cfg = particle_config(
            n_steps=5,
            dt_s=0.02,
            inlet_number_concentration_mb_per_ml=5.0e11,
            acceleration_backend="python",
        )
        dynamics = _dynamics_config(
            integration_substeps=2,
            near_wall_enabled=False,
            collisions_enabled=False,
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
                [vessel],
                particle_cfg,
                dynamics,
                red_blood_cell_network=network,
            )

        count = trajectories.bubble_id.size
        self.assertGreater(count, 0)
        self.assertEqual(
            trajectories.red_blood_cell_velocity_xz_um_s.shape, (count, 2)
        )
        self.assertEqual(
            trajectories.red_blood_cell_drift_velocity_xz_um_s.shape, (count, 2)
        )
        self.assertEqual(
            trajectories.red_blood_cell_fick_velocity_xz_um_s.shape, (count, 2)
        )
        self.assertEqual(
            trajectories.red_blood_cell_quantitative_applicability.shape,
            (count,),
        )
        self.assertTrue(trajectories.metadata["red_blood_cell_transport_enabled"])
        self.assertEqual(
            trajectories.metadata["trajectory_schema"],
            "continuous_perfusion_revised_v20_topological_rbc_drift_diffusion_records_v14",
        )
        self.assertEqual(
            trajectories.metadata["red_blood_cell_transport_model"],
            "RBC-induced drift–diffusion reduced-order transport model",
        )

    def test_default_compact_rbc_records_keep_only_transport_essentials(self) -> None:
        domain, raster, flow = straight_channel_case()
        vessel = _channel_vessel(radius_um=10.0)
        network = build_red_blood_cell_network([vessel], 0.35)
        particle_cfg = particle_config(
            n_steps=3,
            dt_s=0.02,
            inlet_number_concentration_mb_per_ml=5.0e11,
            acceleration_backend="python",
        )
        dynamics = _dynamics_config(
            integration_substeps=2,
            near_wall_enabled=False,
            collisions_enabled=False,
            store_full_diagnostics=False,
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
                [vessel],
                particle_cfg,
                dynamics,
                red_blood_cell_network=network,
            )

        count = trajectories.bubble_id.size
        for name, shape in (
            ("red_blood_cell_velocity_xz_um_s", (count, 2)),
            ("red_blood_cell_drift_velocity_xz_um_s", (count, 2)),
            ("red_blood_cell_fick_velocity_xz_um_s", (count, 2)),
            ("red_blood_cell_transverse_diffusivity_um2_s", (count,)),
            ("red_blood_cell_quantitative_applicability", (count,)),
            ("red_blood_cell_transverse_space_valid", (count,)),
        ):
            self.assertEqual(getattr(trajectories, name).shape, shape)
        for name in (
            "red_blood_cell_local_vessel_diameter_um",
            "red_blood_cell_discharge_hematocrit",
            "red_blood_cell_tube_hematocrit",
            "red_blood_cell_shear_rate_s_inv",
            "red_blood_cell_cfl_width_um",
            "red_blood_cell_target_gap_um",
            "red_blood_cell_margination_length_um",
            "red_blood_cell_margination_time_s",
            "red_blood_cell_scale_activation",
            "red_blood_cell_nearest_wall_unique",
            "red_blood_cell_hematocrit_in_quantitative_range",
            "red_blood_cell_shear_rate_in_quantitative_range",
        ):
            self.assertIsNone(getattr(trajectories, name))
        self.assertEqual(trajectories.metadata["red_blood_cell_rng_seed"], 42)
        self.assertEqual(
            trajectories.metadata["red_blood_cell_rng_algorithm_version"],
            "splitmix64_box_muller_v1",
        )

    def test_production_advect_replays_same_seed_and_changes_one_id_for_new_seed(self) -> None:
        domain, raster, flow = straight_channel_case()
        velocity = np.asarray(flow.velocity_xz_um_s, dtype=np.float32).copy()
        velocity[..., 0] = (
            10.0
            + 2.0
            * np.asarray(domain.z_coordinates_um, dtype=np.float32)[None, :]
        )
        velocity[..., 1] = 0.0
        flow = replace(
            flow,
            velocity_xz_um_s=velocity,
            speed_um_s=np.linalg.norm(velocity, axis=-1),
            hybrid_velocity=rectangular_hybrid_velocity(
                domain,
                velocity_xz=(10.0, 0.0),
                bounds_xz=(-1.0, 20.0, -8.0, 14.0),
            ),
        )
        vessel = _channel_vessel(radius_um=10.0)
        network = build_red_blood_cell_network([vessel], 0.35)
        particle_cfg = particle_config(
            n_steps=6,
            dt_s=0.02,
            inlet_number_concentration_mb_per_ml=1.0e11,
            acceleration_backend="python",
        )
        dynamics = _dynamics_config(
            integration_substeps=4,
            near_wall_enabled=False,
            collisions_enabled=False,
        )

        def run(seed: int):
            return advect_test_particles(
                domain,
                raster,
                flow,
                [vessel],
                particle_cfg,
                dynamics,
                red_blood_cell_network=network,
                random_seed=seed,
            )

        with mock.patch.object(
            particle_perfusion_transport,
            "create_particle_progress_bar",
            return_value=mock.Mock(),
        ):
            first = run(731)
            replay = run(731)
            changed = run(732)

        for name in first.__dataclass_fields__:
            first_value = getattr(first, name)
            replay_value = getattr(replay, name)
            if isinstance(first_value, np.ndarray):
                np.testing.assert_array_equal(
                    replay_value,
                    first_value,
                    err_msg=f"same-seed replay changed array {name}",
                )

        common_ids = np.intersect1d(
            np.unique(first.bubble_id), np.unique(changed.bubble_id)
        )
        self.assertGreater(common_ids.size, 0)
        one_id_changed = False
        for permanent_id in common_ids:
            first_path = first.positions_um[first.bubble_id == permanent_id]
            changed_path = changed.positions_um[changed.bubble_id == permanent_id]
            if first_path.shape != changed_path.shape or not np.array_equal(
                first_path, changed_path
            ):
                one_id_changed = True
                break
        self.assertTrue(one_id_changed)


def _vessel(
    vessel_id: int,
    parent_id: int,
    children: list[int],
    flow_rate_um3_s: float,
    radius_um: float,
) -> Vessel:
    start = float(vessel_id) * 100.0
    return Vessel(
        vid=vessel_id,
        parent_id=parent_id,
        children=children,
        x_p=np.asarray([start, 0.0, 0.0]),
        x_d=np.asarray([start + 100.0, 0.0, 0.0]),
        radius=radius_um,
        flow_rate=flow_rate_um3_s,
    )


def _channel_vessel(*, radius_um: float) -> Vessel:
    return Vessel(
        vid=0,
        parent_id=-1,
        children=[],
        x_p=np.asarray([0.0, 0.0, 3.0]),
        x_d=np.asarray([19.0, 0.0, 3.0]),
        radius=radius_um,
        flow_rate=50.0,
    )


def _evaluate_particle(
    network,
    *,
    gap_um: float,
    radius_um: float,
    unique: bool,
    use_numba: bool,
):
    return evaluate_red_blood_cell_transport(
        sampled_trajectory_vessel_id=np.asarray([1], dtype=np.int32),
        velocity_gradient_s_inv=np.asarray(
            [[[0.0, 1000.0], [0.0, 0.0]]], dtype=np.float64
        ),
        wall_gap_um=np.asarray([gap_um], dtype=np.float64),
        inward_wall_normal_xz=np.asarray([[0.0, 1.0]], dtype=np.float64),
        nearest_wall_unique=np.asarray([unique], dtype=bool),
        bubble_radius_um=np.asarray([radius_um], dtype=np.float64),
        active=np.asarray([True], dtype=bool),
        network=network,
        use_numba=use_numba,
    )


if __name__ == "__main__":
    unittest.main()
