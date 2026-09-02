from __future__ import annotations

import unittest
from unittest import mock
from types import SimpleNamespace

import numpy as np
from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.utils.core.types import GridDomain
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_predictive_contact import (
    PredictiveWallContactGeometryError,
    solve_predictive_contact_step,
)
from ulm_microbubble_traj_gen_2D.utils.particles import particle_perfusion_transport as transport
from ulm_microbubble_traj_gen_2D.utils.particles.particle_constrained_step import (
    PhysicalTimeInterval,
)
class ParticlePredictiveContactV15Tests(unittest.TestCase):
    def setUp(self) -> None:
        self.geometry = _straight_channel_geometry()

    def solve(
        self,
        position: tuple[float, float],
        velocity: tuple[float, float, float],
        *,
        dt_s: float = 1.0,
        radius_um: float = 0.5,
        mobility: np.ndarray | None = None,
        use_numba: bool = False,
        geometry_tolerance_um: float = 1.0e-8,
    ):
        return solve_predictive_contact_step(
            np.asarray([position], dtype=np.float64),
            np.asarray([True]),
            dt_s,
            np.asarray([velocity], dtype=np.float64),
            np.asarray(
                [np.eye(3) if mobility is None else mobility],
                dtype=np.float64,
            ),
            np.asarray([radius_um], dtype=np.float64),
            grid_spacing_um=1.0,
            boundary_geometry=self.geometry,
            geometry_tolerance_um=geometry_tolerance_um,
            use_numba=use_numba,
        )

    def test_free_step_consumes_the_complete_physical_interval(self) -> None:
        result = self.solve((2.0, 2.0), (0.25, 0.1, 0.3), dt_s=0.4)

        np.testing.assert_allclose(
            result.accepted_position_grid[0], [2.1, 2.04], atol=1.0e-14
        )
        np.testing.assert_allclose(
            result.constrained_generalized_velocity[0], [0.25, 0.1, 0.3]
        )
        self.assertAlmostEqual(result.angle_increment_rad[0], 0.12)
        self.assertEqual(result.reaction_force_pn[0], 0.0)
        self.assertFalse(result.active_contact[0])
        self.assertFalse(result.need_time_refinement[0])
        self.assertGreaterEqual(result.minimum_path_gap_um[0], 0.0)

    def test_closed_form_reaction_preserves_full_mobility_coupling(self) -> None:
        mobility = np.asarray(
            [
                [1.0, 0.5, 0.0],
                [0.5, 1.5, 0.0],
                [0.0, -0.4, 1.0],
            ],
            dtype=np.float64,
        )
        result = self.solve(
            (3.0, 1.2), (0.0, -1.0, 0.5), mobility=mobility
        )

        expected_reaction = 0.8 / 1.5
        expected_velocity = np.asarray([0.5, 1.5, -0.4]) * expected_reaction
        expected_velocity += np.asarray([0.0, -1.0, 0.5])
        self.assertAlmostEqual(
            result.reaction_force_pn[0], expected_reaction, places=12
        )
        np.testing.assert_allclose(
            result.constrained_generalized_velocity[0],
            expected_velocity,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            result.accepted_position_grid[0],
            np.asarray([3.0, 1.2]) + expected_velocity[:2],
            atol=1.0e-12,
        )
        self.assertTrue(result.active_contact[0])
        self.assertGreaterEqual(result.endpoint_gap_um[0], 0.0)
        self.assertGreaterEqual(result.minimum_path_gap_um[0], 0.0)
        self.assertLessEqual(
            result.complementarity_residual_pn_um[0], 1.0e-12
        )

    def test_numba_and_reference_transactions_agree(self) -> None:
        reference = self.solve((3.0, 1.2), (0.4, -1.0, 0.2))
        accelerated = self.solve(
            (3.0, 1.2), (0.4, -1.0, 0.2), use_numba=True
        )

        for field in (
            "accepted_position_grid",
            "constrained_generalized_velocity",
            "reaction_force_pn",
            "endpoint_gap_um",
            "minimum_path_gap_um",
            "contact_normal_xz",
        ):
            np.testing.assert_allclose(
                getattr(accelerated, field), getattr(reference, field), atol=1.0e-13
            )

    def test_numba_and_reference_direct_outlet_prefixes_agree(self) -> None:
        reference = self.solve((5.0, 2.0), (2.0, 0.0, 0.0))
        accelerated = self.solve(
            (5.0, 2.0),
            (2.0, 0.0, 0.0),
            use_numba=True,
        )

        for field in (
            "accepted_position_grid",
            "outlet_fraction",
            "outlet_position_grid",
            "outlet_label",
            "endpoint_gap_um",
            "minimum_path_gap_um",
        ):
            np.testing.assert_allclose(
                getattr(accelerated, field),
                getattr(reference, field),
                atol=1.0e-13,
            )
        self.assertAlmostEqual(float(reference.outlet_fraction[0]), 0.5)
        np.testing.assert_allclose(
            reference.accepted_position_grid[0],
            [6.0, 2.0],
        )

    def test_wall_tangent_outlet_corner_is_claimed_as_an_outlet(self) -> None:
        # Regression for permanent bubble 259 at t=0.571775 s.  Its disc was
        # tangent to the wall beside an outlet, so the wall event preceded the
        # centre-plane crossing even though the full swept disc never
        # penetrated the wall.  That legal sliding exit must not become an
        # unclaimed boundary-lifecycle failure.
        for use_numba in (False, True):
            with self.subTest(use_numba=use_numba):
                result = self.solve(
                    (5.999, 1.0),
                    (0.002, 0.0, 0.0),
                    radius_um=0.5,
                    use_numba=use_numba,
                )

                self.assertAlmostEqual(float(result.outlet_fraction[0]), 0.5)
                self.assertEqual(int(result.outlet_label[0]), 1)
                np.testing.assert_allclose(
                    result.accepted_position_grid[0],
                    [6.0, 1.0],
                    atol=1.0e-13,
                )
                self.assertFalse(result.need_time_refinement[0])
                self.assertGreaterEqual(result.minimum_path_gap_um[0], 0.0)

    def test_outlet_crossing_does_not_hide_prior_wall_penetration(self) -> None:
        # The centre line still intersects the outlet cap, but the finite disc
        # reaches and attempts to penetrate the lower wall first.  A singular
        # normal mobility prevents a valid contact correction, so the trial
        # must request refinement instead of being mislabeled as an outlet.
        mobility = np.diag([1.0, 0.0, 1.0])
        for use_numba in (False, True):
            with self.subTest(use_numba=use_numba):
                result = self.solve(
                    (5.0, 1.1),
                    (2.0, -0.3, 0.0),
                    radius_um=0.5,
                    mobility=mobility,
                    use_numba=use_numba,
                )

                self.assertTrue(result.need_time_refinement[0])
                self.assertFalse(np.isfinite(result.outlet_fraction[0]))

    def test_true_simultaneous_distinct_walls_is_a_hard_model_error(self) -> None:
        with self.assertRaises(PredictiveWallContactGeometryError) as captured:
            self.solve(
                (3.0, 2.0),
                (0.0, -0.1, 0.0),
                radius_um=1.5,
                use_numba=True,
            )

        failure = captured.exception.failures[0]
        self.assertTrue(failure.multiple_wall_contact)
        self.assertEqual(
            failure.reason, "simultaneous_distinct_solid_wall_contact"
        )
        self.assertEqual(failure.local_lane, 0)

    def test_nearby_second_wall_inside_residual_tolerance_is_not_contact(self) -> None:
        # The 0.0015 um configuration tolerance limits endpoint roundoff
        # projection.  It must not turn a nearby, still separated wall into a
        # simultaneous physical contact.  In this 3 um channel the disc touches
        # the lower wall while retaining a 0.001 um gap to the upper wall.
        result = self.solve(
            (3.0, 1.9995),
            (0.0, -0.1, 0.0),
            radius_um=1.4995,
            use_numba=True,
            geometry_tolerance_um=0.0015,
        )

        self.assertTrue(result.active_contact[0])
        self.assertFalse(result.need_time_refinement[0])
        self.assertGreaterEqual(result.endpoint_gap_um[0], 0.0)
        np.testing.assert_allclose(
            result.contact_normal_xz[0], [0.0, 1.0], atol=1.0e-13
        )

    def test_roundoff_negative_endpoint_is_projected_strictly_feasible(self) -> None:
        for use_numba in (False, True):
            with self.subTest(use_numba=use_numba):
                result = self.solve(
                    (3.0, 1.0001),
                    (0.0, -0.000100000001, 0.0),
                    radius_um=0.5,
                    use_numba=use_numba,
                )
                accepted_world = self.geometry.grid_to_world_xz(
                    result.accepted_position_grid[0]
                )
                exact_gap = float(
                    self.geometry.exact_true_gap_at_xz_um(
                        accepted_world, 0.5
                    )
                )

                self.assertGreaterEqual(exact_gap, 0.0)
                self.assertGreaterEqual(result.endpoint_gap_um[0], 0.0)
                self.assertFalse(result.need_time_refinement[0])

    def test_reported_id47_location_is_exactly_penetrating_one_wall(self) -> None:
        # Regression for the original v14 invalid_gap_bracket traceback.  Its
        # bilinear raster gap was approximately zero, while the authoritative
        # finite face already gave a negative physical gap.
        position = np.asarray(
            [1276.6890241357419, 1829.358070443591], dtype=np.float64
        )
        center = np.asarray(
            [1276.5760868277539, 1830.1080443518688], dtype=np.float64
        )
        radius = 0.77336554823843295
        face_start = center + np.asarray([-0.75, 0.0])
        face_end = center + np.asarray([0.75, 0.0])

        delta = face_end - face_start
        fraction = np.clip(
            ((position - face_start) @ delta) / (delta @ delta), 0.0, 1.0
        )
        exact_distance = np.linalg.norm(position - (face_start + fraction * delta))
        exact_gap = exact_distance - radius

        self.assertAlmostEqual(exact_distance, 0.74997390827775234, places=12)
        self.assertAlmostEqual(exact_gap, -0.02339163996068061, places=12)
        self.assertLess(exact_gap, 0.0)

    def test_transport_failure_adds_permanent_id_time_and_wall_candidates(self) -> None:
        context = SimpleNamespace(
            spacing_um=1.0,
            use_numba=True,
            boundary_geometry=self.geometry,
        )
        with self.assertRaises(
            transport.ParticleWallContactGeometryError
        ) as captured:
            transport._solve_v15_trial_with_failure_context(
                np.asarray([[3.0, 2.0]], dtype=np.float64),
                np.asarray([True]),
                0.5,
                np.asarray([[0.0, -0.1, 0.0]], dtype=np.float64),
                np.eye(3, dtype=np.float64)[None, :, :],
                np.asarray([1.5], dtype=np.float64),
                context,
                1.0e-8,
                permanent_ids=np.asarray([47], dtype=np.int64),
                interval=PhysicalTimeInterval(6.0, 6.5),
                integration_stage="predictor",
            )

        message = str(captured.exception)
        self.assertIn("permanent_microbubble_id=47", message)
        self.assertIn("physical_time_s=6", message)
        self.assertIn("multiple_wall_contact=true", message)
        self.assertIn("candidate_wall_count=", message)
        self.assertIn("position_xz_um=", message)
        self.assertIn("gap_um=", message)

    def test_refinement_exhaustion_keeps_complete_wall_diagnostics(self) -> None:
        context = SimpleNamespace(
            all_bubble_radii_um=np.asarray([0.5]),
            boundary_geometry=self.geometry,
        )
        start = transport._BatchAdvanceState(
            particles=transport.BatchLocalState(
                position_grid=np.asarray([[3.0, 1.0]], dtype=np.float64),
                rotation_angle_rad=np.asarray([0.0]),
                bond_count_expected=np.asarray([0.0]),
                bond_total_tangential_extension_um=np.asarray([0.0]),
            ),
            alive=np.asarray([True]),
            termination_reason=np.asarray([0], dtype=np.uint8),
            exit_time_s=np.asarray([np.nan]),
        )
        rejected = transport._RefineContactTimeStep(
            "rejected swept chord",
            failed_lanes=np.asarray([0]),
            failure_codes=np.asarray([4]),
            failure_positions_grid=np.asarray([[3.0, 0.9]]),
            failure_event_fractions=np.asarray([1.0]),
            integration_stage="predictor",
        )

        with (
            mock.patch.object(transport, "_attempt_batch_step", side_effect=rejected),
            self.assertRaises(transport.ParticleWallContactGeometryError) as captured,
        ):
            transport._advance_batch_interval(
                PhysicalTimeInterval(6.0, 6.5),
                start,
                np.asarray([0], dtype=np.int64),
                context,
                SimpleNamespace(contact_max_time_refinements=0),
                SimpleNamespace(time_integrator="euler"),
                object(),
            )

        message = str(captured.exception)
        self.assertIn("permanent_microbubble_id=0", message)
        self.assertIn("physical_time_s=6.5", message)
        self.assertIn("integration_stage=predictor", message)
        self.assertIn("position_grid=[3, 0.90000000000000002]", message)
        self.assertIn("candidate_wall_count=", message)
        self.assertIn("multiple_wall_contact=", message)
        self.assertIn("refinement_exhausted_swept_path_penetration", message)
        self.assertIn("refinement_stop_reason=depth_limit", message)
        self.assertIn("duration_s=0.5", message)
        self.assertIn("depth=0", message)
        self.assertIn("maximum_refinement_depth=0", message)

    def test_unrepresentable_refinement_reports_time_floor_not_depth_limit(
        self,
    ) -> None:
        context = SimpleNamespace(
            all_bubble_radii_um=np.asarray([0.5]),
            boundary_geometry=self.geometry,
        )
        start = transport._BatchAdvanceState(
            particles=transport.BatchLocalState(
                position_grid=np.asarray([[3.0, 1.0]], dtype=np.float64),
                rotation_angle_rad=np.asarray([0.0]),
                bond_count_expected=np.asarray([0.0]),
                bond_total_tangential_extension_um=np.asarray([0.0]),
            ),
            alive=np.asarray([True]),
            termination_reason=np.asarray([0], dtype=np.uint8),
            exit_time_s=np.asarray([np.nan]),
        )
        rejected = transport._RefineContactTimeStep(
            "rejected swept chord",
            failed_lanes=np.asarray([0]),
            failure_codes=np.asarray([4]),
            failure_positions_grid=np.asarray([[3.0, 0.9]]),
            failure_event_fractions=np.asarray([1.0]),
            integration_stage="predictor",
        )
        interval_start = 0.46333585277642331
        interval_end = np.nextafter(interval_start, np.inf)

        with (
            mock.patch.object(transport, "_attempt_batch_step", side_effect=rejected),
            self.assertRaises(transport.ParticleWallContactGeometryError) as captured,
        ):
            transport._advance_batch_interval(
                PhysicalTimeInterval(interval_start, interval_end, refinement_depth=7),
                start,
                np.asarray([0], dtype=np.int64),
                context,
                SimpleNamespace(contact_max_time_refinements=12),
                SimpleNamespace(time_integrator="euler"),
                object(),
            )

        message = str(captured.exception)
        self.assertIn("refinement_stop_reason=unrepresentable_midpoint", message)
        self.assertIn("duration_s=5.5511151231257827e-17", message)
        self.assertIn("local_time_ulp_s=5.5511151231257827e-17", message)
        self.assertIn("depth=7", message)
        self.assertIn("maximum_refinement_depth=12", message)
        self.assertNotIn("contact_max_time_refinements was reached", message)


def _straight_channel_geometry():
    shape = (7, 5)
    domain = GridDomain(
        origin_um=np.asarray([0.0, 0.0, 0.0]),
        spacing_um=1.0,
        shape=shape,
        fixed_y_um=0.0,
        x_coordinates_um=np.arange(shape[0], dtype=np.float64),
        z_coordinates_um=np.arange(shape[1], dtype=np.float64),
    )
    vessel = Vessel(vid=0, parent_id=-1, children=[])
    vessel.x_p = np.asarray([0.0, 0.0, 2.0])
    vessel.x_d = np.asarray([6.0, 0.0, 2.0])
    vessel.radius = 1.5
    return build_continuous_vessel_geometry([vessel], domain)


if __name__ == "__main__":
    unittest.main()
