from __future__ import annotations

import math
import unittest

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.particles.particle_mobility import (
    DEFAULT_XI_FAR,
    DEFAULT_XI_MIN,
    DEFAULT_XI_NEAR,
    FREE_ROTATIONAL_SCALED_MOBILITY,
    apply_scaled_mobility_scalar as apply_scaled_mobility,
    background_hydrodynamic_velocity_scalar,
    bulk_angular_velocity_scalar as bulk_angular_velocity,
    effective_dimensionless_mobility_entries_scalar as effective_dimensionless_mobility_entries,
    global_to_wall_components_scalar as global_to_wall_components,
    generalized_mobility_matrix_xz,
    local_shear_rate_scalar as local_shear_rate,
    tangent_from_inward_normal_scalar as tangent_from_inward_normal,
    translational_stokes_mobility_scalar as translational_stokes_mobility,
    wall_blend_weight_scalar as wall_blend_weight,
    wall_dimensionless_mobility_entries_scalar as wall_dimensionless_mobility_entries,
    wall_resistance_coefficients_scalar as wall_resistance_coefficients,
    wall_shear_load_scalar as wall_shear_load,
    wall_to_global_components_scalar as wall_to_global_components,
)


class ParticleUnitTests(unittest.TestCase):
    def test_micrometer_piconewton_stokes_mobility_matches_si(self) -> None:
        viscosity_pa_s = 0.0035
        radius_um = 1.0
        force_pn = 2.0
        velocity_um_s = translational_stokes_mobility(viscosity_pa_s, radius_um) * force_pn
        velocity_si_um_s = (
            force_pn * 1.0e-12
            / (6.0 * math.pi * viscosity_pa_s * radius_um * 1.0e-6)
            * 1.0e6
        )
        self.assertAlmostEqual(float(velocity_um_s), velocity_si_um_s, places=12)

    def test_generalized_matrix_matches_the_existing_scaled_mobility_map(self) -> None:
        viscosity = np.asarray([0.004], dtype=np.float64)
        radius = np.asarray([1.25], dtype=np.float64)
        normal = np.asarray([[0.6, 0.8]], dtype=np.float64)
        entries = np.asarray([[0.7, 0.2, -0.15, -0.1, 0.5]], dtype=np.float64)
        force_global = np.asarray([2.0, -3.0], dtype=np.float64)
        torque = 0.4

        matrix = generalized_mobility_matrix_xz(
            viscosity,
            radius,
            normal,
            entries,
        )[0]
        matrix_result = matrix @ np.asarray(
            [force_global[0], force_global[1], torque],
            dtype=np.float64,
        )

        tangent = np.asarray([-normal[0, 1], normal[0, 0]])
        force_t = float(np.dot(tangent, force_global))
        force_n = float(np.dot(normal[0], force_global))
        velocity_t, velocity_n, omega = apply_scaled_mobility(
            float(viscosity[0]),
            float(radius[0]),
            force_t,
            force_n,
            torque,
            *entries[0],
        )
        expected = np.asarray(
            [
                tangent[0] * velocity_t + normal[0, 0] * velocity_n,
                tangent[1] * velocity_t + normal[0, 1] * velocity_n,
                omega,
            ]
        )
        np.testing.assert_allclose(matrix_result, expected, rtol=1.0e-13, atol=1.0e-13)


class WallCoordinateTests(unittest.TestCase):
    def test_tangent_cross_normal_points_along_positive_y(self) -> None:
        normal_x = 0.6
        normal_z = 0.8
        tangent_x, tangent_z = tangent_from_inward_normal(normal_x, normal_z)
        cross_y = tangent_z * normal_x - tangent_x * normal_z
        self.assertAlmostEqual(cross_y, 1.0)

    def test_local_global_transform_is_reversible(self) -> None:
        normal_x = np.array([1.0, 0.6])
        normal_z = np.array([0.0, 0.8])
        vector_x = np.array([3.0, -2.0])
        vector_z = np.array([-4.0, 5.0])
        local = [
            global_to_wall_components(vector_x[i], vector_z[i], normal_x[i], normal_z[i])
            for i in range(2)
        ]
        recovered = [
            wall_to_global_components(local[i][0], local[i][1], normal_x[i], normal_z[i])
            for i in range(2)
        ]
        recovered_x, recovered_z = np.asarray(recovered).T
        np.testing.assert_allclose(recovered_x, vector_x)
        np.testing.assert_allclose(recovered_z, vector_z)

    def test_signed_shear_and_angular_velocity_follow_xz_convention(self) -> None:
        shear = local_shear_rate(0.0, 0.0, 4.0, 0.0, 1.0, 0.0)
        omega = bulk_angular_velocity(0.0, 4.0)
        self.assertEqual(float(shear), 4.0)
        self.assertEqual(float(omega), -2.0)


class NearWallMobilityTests(unittest.TestCase):
    def test_published_resistance_coefficients(self) -> None:
        xi = 0.01
        c_ft, c_tt, c_fr, c_tr, determinant, normal_resistance = (
            wall_resistance_coefficients(xi)
        )
        log_xi = math.log(xi)
        self.assertAlmostEqual(float(c_ft), (8.0 / 15.0) * log_xi - 0.9588)
        self.assertAlmostEqual(float(c_tt), 0.1 * log_xi - 0.1895)
        self.assertAlmostEqual(float(c_fr), (2.0 / 15.0) * log_xi - 0.2526)
        self.assertAlmostEqual(float(c_tr), 0.4 * log_xi - 0.3187)
        self.assertAlmostEqual(float(determinant), c_tt * c_fr - c_ft * c_tr)
        self.assertAlmostEqual(float(normal_resistance), 1.0 + 1.0 / xi)

    def test_wall_matrix_retains_distinct_published_cross_entries(self) -> None:
        entries = wall_dimensionless_mobility_entries(0.01)
        self.assertAlmostEqual(float(entries[2]), 0.09536646220140073)
        self.assertAlmostEqual(float(entries[3]), 0.09537379845535458)
        self.assertNotEqual(float(entries[2]), float(entries[3]))

    def test_normal_mobility_tends_to_zero_with_gap(self) -> None:
        gap_ratios = np.array([1.0e-1, 1.0e-2, 1.0e-3])
        m_nn = np.asarray(
            [wall_dimensionless_mobility_entries(float(gap))[1] for gap in gap_ratios]
        )
        np.testing.assert_allclose(m_nn, gap_ratios / (1.0 + gap_ratios))
        self.assertTrue(np.all(np.diff(m_nn) < 0.0))

    def test_effective_matrix_has_nonnegative_symmetric_quadratic_form(self) -> None:
        loads = np.array(
            [
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 1.0],
                [2.0, -3.0, 4.0],
                [-5.0, 1.5, -0.25],
            ]
        )
        for xi in np.geomspace(DEFAULT_XI_MIN, DEFAULT_XI_FAR, 40):
            m_tt, m_nn, m_t_r, m_r_t, m_r_r = (
                effective_dimensionless_mobility_entries(xi)
            )
            matrix = np.array(
                [
                    [m_tt, 0.0, m_t_r],
                    [0.0, m_nn, 0.0],
                    [m_r_t, 0.0, m_r_r],
                ],
                dtype=float,
            )
            dissipation = np.einsum("ni,ij,nj->n", loads, matrix, loads)
            self.assertTrue(np.all(dissipation >= -1.0e-12))

    def test_smoothstep_weight_and_derivatives_at_transition_boundaries(self) -> None:
        self.assertEqual(float(wall_blend_weight(DEFAULT_XI_NEAR)), 1.0)
        self.assertEqual(float(wall_blend_weight(DEFAULT_XI_FAR)), 0.0)
        midpoint = 0.5 * (DEFAULT_XI_NEAR + DEFAULT_XI_FAR)
        self.assertAlmostEqual(float(wall_blend_weight(midpoint)), 0.5)
        epsilon = 1.0e-6
        near_slope = (
            wall_blend_weight(DEFAULT_XI_NEAR + epsilon)
            - wall_blend_weight(DEFAULT_XI_NEAR - epsilon)
        ) / (2.0 * epsilon)
        far_slope = (
            wall_blend_weight(DEFAULT_XI_FAR + epsilon)
            - wall_blend_weight(DEFAULT_XI_FAR - epsilon)
        ) / (2.0 * epsilon)
        self.assertAlmostEqual(float(near_slope), 0.0, places=4)
        self.assertAlmostEqual(float(far_slope), 0.0, places=4)

    def test_effective_matrix_freezes_wall_asymptote_and_reaches_free_space(self) -> None:
        at_near = effective_dimensionless_mobility_entries(DEFAULT_XI_NEAR)
        wall_at_near = wall_dimensionless_mobility_entries(DEFAULT_XI_NEAR)
        np.testing.assert_allclose(at_near, wall_at_near)

        at_far = effective_dimensionless_mobility_entries(DEFAULT_XI_FAR)
        np.testing.assert_allclose(at_far, (1.0, 1.0, 0.0, 0.0, 0.75))

        transition_xi = 0.55
        transition = effective_dimensionless_mobility_entries(transition_xi)
        weight = wall_blend_weight(transition_xi)
        expected = (
            (1.0 - weight) + weight * wall_at_near[0],
            (1.0 - weight) + weight * wall_at_near[1],
            weight * wall_at_near[2],
            weight * wall_at_near[3],
            (1.0 - weight) * FREE_ROTATIONAL_SCALED_MOBILITY
            + weight * wall_at_near[4],
        )
        np.testing.assert_allclose(transition, expected)

        below_cutoff = effective_dimensionless_mobility_entries(-0.5)
        at_cutoff = effective_dimensionless_mobility_entries(DEFAULT_XI_MIN)
        np.testing.assert_allclose(below_cutoff, at_cutoff)

    def test_free_space_scaled_matrix_recovers_translational_and_rotational_stokes_laws(self) -> None:
        viscosity_pa_s = 0.004
        radius_um = 1.2
        force_t_pn = 3.0
        force_n_pn = -2.0
        torque_pn_um = 5.0
        velocity_t, velocity_n, omega = apply_scaled_mobility(
            viscosity_pa_s,
            radius_um,
            force_t_pn,
            force_n_pn,
            torque_pn_um,
            1.0,
            1.0,
            0.0,
            0.0,
            0.75,
        )
        self.assertAlmostEqual(
            float(velocity_t),
            force_t_pn / (6.0 * math.pi * viscosity_pa_s * radius_um),
        )
        self.assertAlmostEqual(
            float(velocity_n),
            force_n_pn / (6.0 * math.pi * viscosity_pa_s * radius_um),
        )
        self.assertAlmostEqual(
            float(omega),
            torque_pn_um / (8.0 * math.pi * viscosity_pa_s * radius_um**3),
        )


class WallShearBackgroundTests(unittest.TestCase):
    def test_contact_shear_load_uses_goldman_coefficients_and_center_height(self) -> None:
        viscosity_pa_s = 0.0035
        radius_um = 1.0
        shear_rate_per_s = 100.0
        force_t, torque_y, center_height = wall_shear_load(
            viscosity_pa_s,
            radius_um,
            0.0,
            shear_rate_per_s,
        )
        self.assertEqual(float(center_height), radius_um)
        self.assertNotEqual(float(force_t), 0.0)
        self.assertNotEqual(float(torque_y), 0.0)
        expected_force = (
            1.7005
            * 6.0
            * math.pi
            * viscosity_pa_s
            * radius_um
            * radius_um
            * shear_rate_per_s
        )
        expected_torque = (
            0.9440
            * 4.0
            * math.pi
            * viscosity_pa_s
            * radius_um**3
            * shear_rate_per_s
        )
        self.assertAlmostEqual(float(force_t), expected_force)
        self.assertAlmostEqual(float(torque_y), expected_torque)

    def test_background_velocity_is_exactly_bulk_at_far_boundary(self) -> None:
        result = background_hydrodynamic_velocity_scalar(
            7.0,
            -2.0,
            0.0,
            3.0,
            -1.0,
            0.0,
            1.0,
            0.0,
            1.0,
            DEFAULT_XI_FAR,
            0.004,
        )
        self.assertAlmostEqual(result[0], 7.0)
        self.assertAlmostEqual(result[1], -2.0)
        self.assertAlmostEqual(result[2], 2.0)
        self.assertEqual(result[3], 0.0)
        self.assertEqual(result[4], -1.0)

    def test_far_bulk_velocity_does_not_depend_on_a_degenerate_wall_normal(self) -> None:
        result = background_hydrodynamic_velocity_scalar(
            7.0,
            -2.0,
            0.0,
            3.0,
            -1.0,
            0.0,
            0.0,
            0.0,
            1.0,
            4.0,
            0.004,
        )
        self.assertEqual(result[:4], (7.0, -2.0, 2.0, 0.0))

    def test_contact_wall_shear_velocity_is_finite_and_signed(self) -> None:
        result = background_hydrodynamic_velocity_scalar(
            0.0,
            0.0,
            0.0,
            0.0,
            10.0,
            0.0,
            1.0,
            0.0,
            1.0,
            0.0,
            0.004,
        )
        self.assertTrue(np.all(np.isfinite(result)))
        self.assertGreater(result[1], 0.0)
        self.assertEqual(result[3], 1.0)
        self.assertEqual(result[4], 10.0)


if __name__ == "__main__":
    unittest.main()
