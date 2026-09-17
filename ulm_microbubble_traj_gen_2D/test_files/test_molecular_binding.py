from __future__ import annotations

import math
import unittest

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_binding import (
    BOLTZMANN_CONSTANT_J_PER_K,
    JOULES_PER_PICONEWTON_NANOMETER,
    MolecularBindingParameters,
    accept_bond_state_exponential_heun,
    evaluate_mean_field_bonds,
    predict_bond_state_exponential_euler,
    reaction_disk_radius_um,
    surface_slip_velocity_um_s,
    target_positive_reaction_area_um2,
)


def _parameters(**overrides: float) -> MolecularBindingParameters:
    values = {
        "ligand_density_molecules_m2": 2.0e12,
        "target_density_molecules_m2": 3.0e12,
        "association_rate_m2_per_molecule_s": 1.0e-12,
        "zero_force_dissociation_rate_s": 0.25,
        "rest_length_um": 0.1,
        "spring_stiffness_pn_per_um": 10.0,
        "force_sensitivity_length_nm": 0.1,
        "temperature_k": 310.0,
        "bell_exponent_limit": 80.0,
    }
    values.update(overrides)
    return MolecularBindingParameters(**values)


class ReactionGeometryTests(unittest.TestCase):
    def test_capture_projection_uses_physical_gap_and_spherical_cap(self) -> None:
        radius = 1.0
        capture = 0.2
        result = reaction_disk_radius_um(
            radius,
            np.asarray([-1.0e-15, 0.1, 0.2, 0.3]),
            capture,
        )
        expected_delta = np.asarray([0.2, 0.1, 0.0, 0.0])
        expected = np.sqrt(2.0 * radius * expected_delta - expected_delta**2)
        np.testing.assert_allclose(result, expected, rtol=0.0, atol=0.0)

    def test_capture_projection_rejects_nonphysical_spherical_cap(self) -> None:
        with self.assertRaisesRegex(ValueError, "twice radius"):
            reaction_disk_radius_um(1.0, 0.0, 2.01)

    def test_capture_projection_rejects_material_wall_penetration(self) -> None:
        with self.assertRaisesRegex(ValueError, "materially negative"):
            reaction_disk_radius_um(1.0, -1.0e-6, 0.2)

    def test_full_half_and_absent_target_intervals_have_exact_disk_areas(self) -> None:
        radius = 0.7
        full = target_positive_reaction_area_um2(
            radius, np.asarray([-radius]), np.asarray([radius])
        )
        half = target_positive_reaction_area_um2(
            radius, np.asarray([0.0]), np.asarray([radius])
        )
        absent = target_positive_reaction_area_um2(
            radius, np.asarray([radius + 1.0]), np.asarray([radius + 2.0])
        )
        self.assertAlmostEqual(full, math.pi * radius**2, places=15)
        self.assertAlmostEqual(half, 0.5 * math.pi * radius**2, places=15)
        self.assertEqual(absent, 0.0)

    def test_overlapping_target_intervals_are_merged_without_mutating_inputs(self) -> None:
        starts = np.asarray([-1.0, -0.2, 0.7])
        ends = np.asarray([0.1, 0.8, 2.0])
        starts_before = starts.copy()
        ends_before = ends.copy()
        area = target_positive_reaction_area_um2(1.0, starts, ends)

        # The clipped union spans the complete projection disk.
        self.assertAlmostEqual(area, math.pi, places=15)
        np.testing.assert_array_equal(starts, starts_before)
        np.testing.assert_array_equal(ends, ends_before)

    def test_target_area_vanishes_continuously_at_patch_edge(self) -> None:
        radius = 1.0
        areas = [
            target_positive_reaction_area_um2(
                radius, np.asarray([edge]), np.asarray([2.0])
            )
            for edge in (0.0, 0.5, 0.9, 0.99, 1.0)
        ]
        self.assertTrue(np.all(np.diff(areas) < 0.0))
        self.assertEqual(areas[-1], 0.0)


class MeanFieldKineticsTests(unittest.TestCase):
    def test_area_units_capacity_and_formation_rate_follow_revised_v7(self) -> None:
        result = evaluate_mean_field_bonds(
            expected_bond_count=1.0,
            total_tangential_extension_um=0.0,
            gap_um=0.05,
            radius_um=1.0,
            slip_velocity_um_s=4.0,
            reaction_area_um2=1.0,
            parameters=_parameters(),
            use_numba=False,
        )

        # 1 um^2 = 1e-12 m^2, hence N_L=2 and N_T=3.  With n=1,
        # r+ = (1e-12/1e-12) * (2-1) * (3-1) = 2 bonds/s.
        self.assertEqual(float(result.ligand_count), 2.0)
        self.assertEqual(float(result.target_count), 3.0)
        self.assertEqual(float(result.formation_capacity), 2.0)
        self.assertAlmostEqual(float(result.formation_rate_bonds_s), 2.0)
        self.assertAlmostEqual(float(result.dissociation_rate_s), 0.25)
        self.assertAlmostEqual(
            float(result.formation_rate_bonds_s - result.dissociation_rate_s),
            1.75,
        )
        self.assertAlmostEqual(float(result.extension_source_um_s), 4.0)

    def test_no_reaction_area_stops_formation_but_preserves_natural_decay(self) -> None:
        result = evaluate_mean_field_bonds(
            expected_bond_count=2.0,
            total_tangential_extension_um=0.0,
            gap_um=0.05,
            radius_um=1.0,
            slip_velocity_um_s=0.0,
            reaction_area_um2=0.0,
            parameters=_parameters(),
            use_numba=False,
        )
        self.assertEqual(float(result.formation_rate_bonds_s), 0.0)
        self.assertEqual(float(result.expected_bond_count), 2.0)
        self.assertAlmostEqual(
            float(-result.dissociation_rate_s * result.expected_bond_count),
            -0.5,
        )

    def test_existing_bonds_above_a_shrunken_capacity_are_not_clipped(self) -> None:
        result = evaluate_mean_field_bonds(
            expected_bond_count=5.0,
            total_tangential_extension_um=1.0,
            gap_um=0.05,
            radius_um=1.0,
            slip_velocity_um_s=0.0,
            reaction_area_um2=1.0,
            parameters=_parameters(),
            use_numba=False,
        )
        self.assertEqual(float(result.formation_capacity), 2.0)
        self.assertEqual(float(result.expected_bond_count), 5.0)
        self.assertEqual(float(result.formation_rate_bonds_s), 0.0)
        self.assertGreater(float(result.dissociation_rate_s), 0.0)

    def test_zero_bonds_canonicalize_inconsistent_total_extension_locally(self) -> None:
        original_m = np.asarray([7.0])
        result = evaluate_mean_field_bonds(
            expected_bond_count=np.asarray([0.0]),
            total_tangential_extension_um=original_m,
            gap_um=0.2,
            radius_um=1.0,
            slip_velocity_um_s=5.0,
            reaction_area_um2=0.0,
            parameters=_parameters(),
            use_numba=False,
        )
        self.assertEqual(float(result.total_tangential_extension_um[0]), 0.0)
        self.assertEqual(float(result.mean_tangential_extension_um[0]), 0.0)
        self.assertEqual(float(result.extension_source_um_s[0]), 0.0)
        self.assertEqual(float(result.force_t_pn[0]), 0.0)
        np.testing.assert_array_equal(original_m, np.asarray([7.0]))

    def test_zero_bonds_without_formation_do_not_create_a_bell_load(self) -> None:
        results = []
        for use_numba in (False, True):
            with self.subTest(use_numba=use_numba):
                result = evaluate_mean_field_bonds(
                    expected_bond_count=0.0,
                    total_tangential_extension_um=0.0,
                    gap_um=2.5,
                    radius_um=1.0,
                    slip_velocity_um_s=0.0,
                    reaction_area_um2=0.0,
                    parameters=_parameters(
                        spring_stiffness_pn_per_um=1.0e4,
                        force_sensitivity_length_nm=0.039,
                    ),
                    use_numba=use_numba,
                )
                self.assertEqual(float(result.formation_rate_bonds_s), 0.0)
                self.assertEqual(float(result.single_bond_tension_pn), 0.0)
                self.assertEqual(float(result.bell_exponent), 0.0)
                self.assertFalse(bool(result.bell_rate_saturated))
                self.assertEqual(float(result.dissociation_rate_s), 0.25)
                self.assertEqual(float(result.force_t_pn), 0.0)
                self.assertEqual(float(result.force_n_pn), 0.0)
                self.assertEqual(float(result.torque_y_pn_um), 0.0)
                results.append(result)
        for field_name in MolecularBindingEvaluationFields:
            np.testing.assert_allclose(
                getattr(results[1], field_name),
                getattr(results[0], field_name),
                rtol=0.0,
                atol=0.0,
            )

    def test_existing_bonds_without_formation_keep_their_bell_load(self) -> None:
        result = evaluate_mean_field_bonds(
            expected_bond_count=2.0,
            total_tangential_extension_um=0.0,
            gap_um=1.0,
            radius_um=1.0,
            slip_velocity_um_s=0.0,
            reaction_area_um2=0.0,
            parameters=_parameters(
                spring_stiffness_pn_per_um=1.0e4,
                force_sensitivity_length_nm=0.039,
            ),
            use_numba=False,
        )
        self.assertEqual(float(result.formation_rate_bonds_s), 0.0)
        self.assertGreater(float(result.single_bond_tension_pn), 0.0)
        self.assertTrue(bool(result.bell_rate_saturated))
        self.assertGreater(float(result.dissociation_rate_s), 0.25)
        self.assertLess(float(result.force_n_pn), 0.0)

    def test_available_formation_keeps_the_new_bond_bell_load(self) -> None:
        result = evaluate_mean_field_bonds(
            expected_bond_count=0.0,
            total_tangential_extension_um=0.0,
            gap_um=1.0,
            radius_um=1.0,
            slip_velocity_um_s=0.0,
            reaction_area_um2=1.0,
            parameters=_parameters(
                spring_stiffness_pn_per_um=1.0e4,
                force_sensitivity_length_nm=0.039,
            ),
            use_numba=False,
        )
        self.assertGreater(float(result.formation_rate_bonds_s), 0.0)
        self.assertGreater(float(result.single_bond_tension_pn), 0.0)
        self.assertTrue(bool(result.bell_rate_saturated))
        self.assertGreater(float(result.dissociation_rate_s), 0.25)


class BondMechanicsTests(unittest.TestCase):
    def test_tension_only_spring_does_not_resist_compression(self) -> None:
        result = evaluate_mean_field_bonds(
            2.0,
            0.0,
            gap_um=0.05,
            radius_um=1.0,
            slip_velocity_um_s=0.0,
            reaction_area_um2=0.0,
            parameters=_parameters(rest_length_um=0.1),
            use_numba=False,
        )
        self.assertEqual(float(result.single_bond_tension_pn), 0.0)
        self.assertEqual(float(result.force_t_pn), 0.0)
        self.assertEqual(float(result.force_n_pn), 0.0)

    def test_force_and_torque_follow_local_tangent_normal_signs(self) -> None:
        positive = evaluate_mean_field_bonds(
            2.0,
            1.0,
            gap_um=0.3,
            radius_um=1.2,
            slip_velocity_um_s=0.0,
            reaction_area_um2=0.0,
            parameters=_parameters(),
            use_numba=False,
        )
        negative = evaluate_mean_field_bonds(
            2.0,
            -1.0,
            gap_um=0.3,
            radius_um=1.2,
            slip_velocity_um_s=0.0,
            reaction_area_um2=0.0,
            parameters=_parameters(),
            use_numba=False,
        )
        self.assertLess(float(positive.force_t_pn), 0.0)
        self.assertGreater(float(negative.force_t_pn), 0.0)
        self.assertLess(float(positive.force_n_pn), 0.0)
        self.assertLess(float(negative.force_n_pn), 0.0)
        self.assertAlmostEqual(
            float(positive.torque_y_pn_um),
            1.2 * float(positive.force_t_pn),
        )
        self.assertAlmostEqual(
            float(negative.torque_y_pn_um),
            1.2 * float(negative.force_t_pn),
        )

    def test_bell_exponent_uses_piconewton_nanometer_energy_conversion(self) -> None:
        result = evaluate_mean_field_bonds(
            1.0,
            0.0,
            gap_um=0.5,
            radius_um=1.0,
            slip_velocity_um_s=0.0,
            reaction_area_um2=0.0,
            parameters=_parameters(
                rest_length_um=0.1,
                spring_stiffness_pn_per_um=10.0,
                force_sensitivity_length_nm=0.5,
                zero_force_dissociation_rate_s=2.0,
            ),
            use_numba=False,
        )
        tension = 10.0 * (0.5 - 0.1)
        expected_exponent = (
            tension
            * 0.5
            * JOULES_PER_PICONEWTON_NANOMETER
            / (BOLTZMANN_CONSTANT_J_PER_K * 310.0)
        )
        self.assertAlmostEqual(float(result.bell_exponent), expected_exponent)
        self.assertAlmostEqual(
            float(result.dissociation_rate_s), 2.0 * math.exp(expected_exponent)
        )

    def test_bell_exponent_limit_keeps_extreme_loads_finite(self) -> None:
        result = evaluate_mean_field_bonds(
            1.0,
            1.0,
            gap_um=1.0,
            radius_um=1.0,
            slip_velocity_um_s=1.0e300,
            reaction_area_um2=1.0,
            parameters=_parameters(
                spring_stiffness_pn_per_um=1.0e300,
                force_sensitivity_length_nm=1.0e300,
                bell_exponent_limit=40.0,
            ),
            use_numba=False,
        )
        self.assertEqual(float(result.bell_exponent), 40.0)
        self.assertTrue(bool(result.bell_rate_saturated))
        for value in (
            result.dissociation_rate_s,
            result.single_bond_tension_pn,
            result.force_t_pn,
            result.force_n_pn,
            result.torque_y_pn_um,
        ):
            self.assertTrue(np.all(np.isfinite(value)))

    def test_surface_slip_uses_positive_y_rotation_convention(self) -> None:
        slip = surface_slip_velocity_um_s(
            np.asarray([3.0, 3.0]),
            np.asarray([2.0, 2.0]),
            np.asarray([-1.5, 0.5]),
        )
        np.testing.assert_allclose(slip, np.asarray([0.0, 4.0]))


class BatchAndTimeIntegrationTests(unittest.TestCase):
    def test_numba_and_numpy_batch_paths_match_and_do_not_mutate_inputs(self) -> None:
        n = np.asarray([0.0, 0.5, 3.0])
        m = np.asarray([0.0, -0.1, 0.6])
        gaps = np.asarray([-1.0e-9, 0.04, 0.3])
        areas = np.asarray([0.0, 0.8, 1.2])
        originals = tuple(value.copy() for value in (n, m, gaps, areas))
        kwargs = dict(
            expected_bond_count=n,
            total_tangential_extension_um=m,
            gap_um=gaps,
            radius_um=1.0,
            slip_velocity_um_s=np.asarray([0.0, -2.0, 4.0]),
            reaction_area_um2=areas,
            parameters=_parameters(),
        )
        numpy_result = evaluate_mean_field_bonds(**kwargs, use_numba=False)
        numba_result = evaluate_mean_field_bonds(**kwargs, use_numba=True)
        for field_name in MolecularBindingEvaluationFields:
            np.testing.assert_allclose(
                getattr(numba_result, field_name),
                getattr(numpy_result, field_name),
                rtol=1.0e-13,
                atol=1.0e-13,
            )
        for value, original in zip((n, m, gaps, areas), originals):
            np.testing.assert_array_equal(value, original)

    def test_exponential_predictor_integrates_stiff_pure_decay_positively(self) -> None:
        parameters = _parameters(
            association_rate_m2_per_molecule_s=0.0,
            zero_force_dissociation_rate_s=1.0e6,
            rest_length_um=10.0,
        )
        rhs = evaluate_mean_field_bonds(
            4.0, 2.0, 0.0, 1.0, 0.0, 0.0, parameters, use_numba=False
        )
        updated = predict_bond_state_exponential_euler(
            4.0, 2.0, rhs, dt_s=1.0e-3, use_numba=False
        )
        self.assertGreaterEqual(float(updated.expected_bond_count), 0.0)
        self.assertEqual(float(updated.expected_bond_count), 0.0)
        self.assertEqual(float(updated.total_tangential_extension_um), 0.0)

    def test_capacity_limiter_does_not_artificially_scale_total_extension(self) -> None:
        parameters = _parameters(
            association_rate_m2_per_molecule_s=1.0,
            zero_force_dissociation_rate_s=0.0,
            rest_length_um=10.0,
        )
        rhs = evaluate_mean_field_bonds(
            1.0,
            0.5,
            gap_um=0.0,
            radius_um=1.0,
            slip_velocity_um_s=0.0,
            reaction_area_um2=1.0,
            parameters=parameters,
            use_numba=False,
        )
        updated = predict_bond_state_exponential_euler(
            1.0, 0.5, rhs, dt_s=1.0, use_numba=False
        )
        self.assertTrue(bool(updated.capacity_limited))
        self.assertEqual(float(updated.expected_bond_count), 2.0)
        # New bonds enter with zero extension; limiting their numerical count
        # must not delete extension carried by the pre-existing bond group.
        self.assertEqual(float(updated.total_tangential_extension_um), 0.5)

    def test_projected_motion_can_override_the_unconstrained_extension_source(self) -> None:
        parameters = _parameters(
            association_rate_m2_per_molecule_s=0.0,
            zero_force_dissociation_rate_s=0.0,
            rest_length_um=10.0,
        )
        rhs = evaluate_mean_field_bonds(
            2.0, 0.0, 0.0, 1.0, 10.0, 0.0, parameters, use_numba=False
        )
        unconstrained = predict_bond_state_exponential_euler(
            2.0, 0.0, rhs, dt_s=1.0, use_numba=False
        )
        projected_hold = predict_bond_state_exponential_euler(
            2.0,
            0.0,
            rhs,
            dt_s=1.0,
            extension_source_um_s=0.0,
            use_numba=False,
        )
        self.assertEqual(float(unconstrained.total_tangential_extension_um), 20.0)
        self.assertEqual(float(projected_hold.total_tangential_extension_um), 0.0)

    def test_exponential_heun_acceptance_uses_both_rhs_stages(self) -> None:
        parameters = _parameters(
            association_rate_m2_per_molecule_s=0.0,
            zero_force_dissociation_rate_s=2.0,
            rest_length_um=10.0,
        )
        first = evaluate_mean_field_bonds(
            3.0, 1.5, 0.0, 1.0, 0.0, 0.0, parameters, use_numba=False
        )
        predictor = predict_bond_state_exponential_euler(
            3.0, 1.5, first, dt_s=0.2, use_numba=False
        )
        second = evaluate_mean_field_bonds(
            predictor.expected_bond_count,
            predictor.total_tangential_extension_um,
            0.0,
            1.0,
            0.0,
            0.0,
            parameters,
            use_numba=False,
        )
        accepted = accept_bond_state_exponential_heun(
            3.0, 1.5, first, second, dt_s=0.2, use_numba=False
        )
        survival = math.exp(-2.0 * 0.2)
        self.assertAlmostEqual(float(accepted.expected_bond_count), 3.0 * survival)
        self.assertAlmostEqual(
            float(accepted.total_tangential_extension_um), 1.5 * survival
        )


MolecularBindingEvaluationFields = (
    "expected_bond_count",
    "total_tangential_extension_um",
    "mean_tangential_extension_um",
    "reaction_area_um2",
    "ligand_count",
    "target_count",
    "formation_capacity",
    "formation_rate_bonds_s",
    "dissociation_rate_s",
    "slip_velocity_um_s",
    "extension_source_um_s",
    "single_bond_tension_pn",
    "force_t_pn",
    "force_n_pn",
    "torque_y_pn_um",
    "bell_exponent",
    "bell_rate_saturated",
    "formation_rate_saturated",
)


class ValidationTests(unittest.TestCase):
    def test_parameters_reject_nonphysical_values(self) -> None:
        for name, value in (
            ("ligand_density_molecules_m2", -1.0),
            ("temperature_k", 0.0),
            ("bell_exponent_limit", math.inf),
        ):
            with self.subTest(name=name):
                with self.assertRaises(ValueError):
                    _parameters(**{name: value})

    def test_rhs_rejects_negative_bond_count_and_reaction_area(self) -> None:
        with self.assertRaisesRegex(ValueError, "expected_bond_count"):
            evaluate_mean_field_bonds(
                -1.0, 0.0, 0.0, 1.0, 0.0, 1.0, _parameters()
            )
        with self.assertRaisesRegex(ValueError, "reaction_area_um2"):
            evaluate_mean_field_bonds(
                0.0, 0.0, 0.0, 1.0, 0.0, -1.0, _parameters()
            )


if __name__ == "__main__":
    unittest.main()
