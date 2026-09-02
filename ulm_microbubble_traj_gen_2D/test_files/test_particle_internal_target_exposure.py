from __future__ import annotations

from dataclasses import dataclass
from types import SimpleNamespace
import unittest

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.particles.particle_internal_target_exposure import (
    integrate_internal_target_exposure,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_perfusion_transport import (
    _accumulate_internal_target_exposure,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_rbc_stochastic import (
    RbcStochasticSweepResult,
)


@dataclass(frozen=True)
class _WallState:
    distance_um: np.ndarray
    inward_normal_xz: np.ndarray


class _PlanarWallGeometry:
    """Continuous wall z=0 with lumen on the positive-z side."""

    def __init__(self) -> None:
        self.query_count = 0
        self.queried_point_count = 0
        self.query_sizes: list[int] = []

    def exact_solid_wall_state_xz_um_accelerated(
        self, points_xz_um: np.ndarray
    ) -> _WallState:
        self.query_count += 1
        points = np.asarray(points_xz_um, dtype=np.float64)
        self.queried_point_count += int(points.shape[0])
        self.query_sizes.append(int(points.shape[0]))
        normals = np.zeros_like(points)
        normals[:, 1] = 1.0
        return _WallState(distance_um=points[:, 1], inward_normal_xz=normals)

    @staticmethod
    def grid_to_world_xz(points_grid: np.ndarray) -> np.ndarray:
        return np.asarray(points_grid, dtype=np.float64)


class _FullyTargetedWall:
    enabled = True

    def __init__(self) -> None:
        self.query_count = 0
        self.queried_point_count = 0
        self.query_sizes: list[int] = []

    def reaction_area_um2(
        self,
        points_xz_um: np.ndarray,
        tangents_xz: np.ndarray,
        reaction_radius_um: np.ndarray,
    ) -> np.ndarray:
        self.query_count += 1
        point_count = int(len(points_xz_um))
        self.queried_point_count += point_count
        self.query_sizes.append(point_count)
        np.testing.assert_allclose(
            tangents_xz,
            np.tile(np.asarray([-1.0, 0.0]), (len(points_xz_um), 1)),
            atol=0.0,
            rtol=0.0,
        )
        radius = np.asarray(reaction_radius_um, dtype=np.float64)
        return np.pi * radius * radius


class InternalTargetExposureTests(unittest.TestCase):
    def test_stationary_exposure_uses_full_dt_and_marks_right_censoring(self) -> None:
        geometry = _PlanarWallGeometry()
        target = _FullyTargetedWall()
        result = integrate_internal_target_exposure(
            permanent_ids=np.asarray([19]),
            path_points_xz_um=np.asarray([[2.0, 1.1]]),
            path_point_offsets=np.asarray([0, 1]),
            dt_s=2.0,
            radius_um=1.0,
            capture_distance_um=0.2,
            boundary_geometry=geometry,
            molecular_target_field=target,
            quantitatively_applicable=False,
            observation_ends_at_step_end=True,
        )

        self.assertAlmostEqual(float(result.exposure_time_s[0]), 2.0, places=14)
        self.assertAlmostEqual(
            float(result.reaction_area_time_um2_s[0]),
            2.0 * np.pi * 0.19,
            places=13,
        )
        self.assertEqual(float(result.quantitatively_applicable_exposure_time_s[0]), 0.0)
        self.assertEqual(int(result.event_start_count[0]), 1)
        self.assertEqual(int(result.event_end_count[0]), 0)
        self.assertTrue(bool(result.exposure_open_at_end[0]))
        self.assertTrue(bool(result.right_censored_at_end[0]))
        self.assertEqual(geometry.query_count, 1)
        self.assertEqual(target.query_count, 1)

    def test_transient_enter_and_exit_is_seen_when_saved_endpoints_are_outside(self) -> None:
        geometry = _PlanarWallGeometry()
        target = _FullyTargetedWall()
        result = integrate_internal_target_exposure(
            permanent_ids=np.asarray([7]),
            path_points_xz_um=np.asarray(
                [[0.0, 1.3], [0.0, 1.05], [0.0, 1.3]], dtype=np.float64
            ),
            path_point_offsets=np.asarray([0, 3]),
            dt_s=1.0,
            radius_um=1.0,
            capture_distance_um=0.2,
            boundary_geometry=geometry,
            molecular_target_field=target,
        )

        self.assertEqual(int(result.event_start_count[0]), 1)
        self.assertEqual(int(result.event_end_count[0]), 1)
        self.assertFalse(bool(result.exposure_open_at_end[0]))
        self.assertGreater(float(result.exposure_time_s[0]), 0.5)
        self.assertLess(float(result.exposure_time_s[0]), 0.7)
        self.assertGreater(float(result.reaction_area_time_um2_s[0]), 0.0)
        self.assertGreater(geometry.query_count, 1)
        self.assertGreaterEqual(geometry.query_count, target.query_count)
        self.assertGreater(target.query_count, 0)

    def test_all_zero_initial_area_closes_only_initially_open_events(self) -> None:
        geometry = _PlanarWallGeometry()
        target = _FullyTargetedWall()
        result = integrate_internal_target_exposure(
            permanent_ids=np.asarray([31, 9]),
            path_points_xz_um=np.asarray(
                [[0.0, 1.4], [1.0, 1.4], [2.0, 1.4]], dtype=np.float64
            ),
            path_point_offsets=np.asarray([0, 2, 3]),
            dt_s=np.asarray([0.25, 0.5]),
            radius_um=np.asarray([1.0, 1.0]),
            capture_distance_um=0.2,
            boundary_geometry=geometry,
            molecular_target_field=target,
            initially_exposed=np.asarray([True, False]),
            quantitatively_applicable=np.asarray([True, True]),
            terminated_at_end=np.asarray([True, False]),
            observation_ends_at_step_end=np.asarray([True, True]),
        )

        np.testing.assert_array_equal(result.permanent_ids, np.asarray([31, 9]))
        np.testing.assert_array_equal(result.exposure_time_s, np.zeros(2))
        np.testing.assert_array_equal(result.reaction_area_time_um2_s, np.zeros(2))
        np.testing.assert_array_equal(
            result.quantitatively_applicable_exposure_time_s, np.zeros(2)
        )
        np.testing.assert_array_equal(
            result.event_start_count, np.zeros(2, dtype=np.int64)
        )
        np.testing.assert_array_equal(
            result.event_end_count, np.asarray([1, 0], dtype=np.int64)
        )
        np.testing.assert_array_equal(
            result.event_started, np.asarray([False, False])
        )
        np.testing.assert_array_equal(result.event_ended, np.asarray([True, False]))
        np.testing.assert_array_equal(
            result.exposure_open_at_end, np.asarray([False, False])
        )
        np.testing.assert_array_equal(
            result.right_censored_at_end, np.asarray([False, False])
        )
        self.assertEqual(geometry.query_count, 1)
        self.assertGreater(geometry.queried_point_count, 3)
        self.assertEqual(target.query_count, 0)
        self.assertEqual(target.queried_point_count, 0)

    def test_target_query_receives_only_positive_reaction_radius_points(self) -> None:
        geometry = _PlanarWallGeometry()
        target = _FullyTargetedWall()
        result = integrate_internal_target_exposure(
            permanent_ids=np.asarray([4, 8]),
            path_points_xz_um=np.asarray([[0.0, 1.1], [0.0, 1.4]]),
            path_point_offsets=np.asarray([0, 1, 2]),
            dt_s=np.asarray([0.2, 0.3]),
            radius_um=np.asarray([1.0, 1.0]),
            capture_distance_um=0.2,
            boundary_geometry=geometry,
            molecular_target_field=target,
        )

        np.testing.assert_allclose(result.exposure_time_s, np.asarray([0.2, 0.0]))
        np.testing.assert_allclose(
            result.reaction_area_time_um2_s,
            np.asarray([0.2 * np.pi * 0.19, 0.0]),
        )
        self.assertEqual(geometry.query_sizes, [2])
        self.assertEqual(target.query_sizes, [1])
        self.assertEqual(geometry.queried_point_count, 2)
        self.assertEqual(target.queried_point_count, 1)

    def test_extremely_short_vertex_enter_and_exit_is_root_resolved(self) -> None:
        geometry = _PlanarWallGeometry()
        target = _FullyTargetedWall()
        penetration_um = 1.0e-7
        outside_height_um = 1.3
        inside_height_um = 1.2 - penetration_um
        result = integrate_internal_target_exposure(
            permanent_ids=np.asarray([73]),
            path_points_xz_um=np.asarray(
                [
                    [0.0, outside_height_um],
                    [0.0, inside_height_um],
                    [0.0, outside_height_um],
                ],
                dtype=np.float64,
            ),
            path_point_offsets=np.asarray([0, 3]),
            dt_s=1.0,
            radius_um=1.0,
            capture_distance_um=0.2,
            boundary_geometry=geometry,
            molecular_target_field=target,
        )

        # Each symmetric segment consumes half of the step.  The two exposed
        # tips together therefore occupy penetration / segment-height-change
        # of the full dt: about one part per million, far inside the blind spot
        # of an eight-node whole-segment Gauss rule.
        height_change_um = outside_height_um - inside_height_um
        expected_exposure_s = penetration_um / height_change_um
        expected_mean_area_um2 = np.pi * (
            penetration_um - penetration_um**2 / 3.0
        )
        expected_area_time_um2_s = (
            expected_exposure_s * expected_mean_area_um2
        )

        self.assertEqual(int(result.event_start_count[0]), 1)
        self.assertEqual(int(result.event_end_count[0]), 1)
        self.assertFalse(bool(result.exposure_open_at_end[0]))
        self.assertAlmostEqual(
            float(result.exposure_time_s[0]),
            expected_exposure_s,
            delta=2.0e-13,
        )
        self.assertAlmostEqual(
            float(result.reaction_area_time_um2_s[0]),
            expected_area_time_um2_s,
            delta=2.0e-19,
        )
        self.assertGreater(float(result.exposure_time_s[0]), 0.0)
        self.assertGreater(float(result.reaction_area_time_um2_s[0]), 0.0)

    def test_polyline_segments_share_one_physical_time(self) -> None:
        result = integrate_internal_target_exposure(
            permanent_ids=np.asarray([3]),
            path_points_xz_um=np.asarray(
                [[0.0, 1.1], [1.0, 1.1], [2.0, 1.1], [4.0, 1.1]]
            ),
            path_point_offsets=np.asarray([0, 4]),
            dt_s=0.75,
            radius_um=1.0,
            capture_distance_um=0.2,
            boundary_geometry=_PlanarWallGeometry(),
            molecular_target_field=_FullyTargetedWall(),
        )

        self.assertAlmostEqual(float(result.exposure_time_s[0]), 0.75, places=14)
        self.assertAlmostEqual(
            float(result.reaction_area_time_um2_s[0]),
            0.75 * np.pi * 0.19,
            places=13,
        )
        self.assertEqual(int(result.event_start_count[0]), 1)
        self.assertEqual(int(result.event_end_count[0]), 0)

    def test_termination_closes_an_open_exposure_event(self) -> None:
        result = integrate_internal_target_exposure(
            permanent_ids=np.asarray([11]),
            path_points_xz_um=np.asarray([[0.0, 1.1], [1.0, 1.1]]),
            path_point_offsets=np.asarray([0, 2]),
            dt_s=0.1,
            radius_um=1.0,
            capture_distance_um=0.2,
            boundary_geometry=_PlanarWallGeometry(),
            molecular_target_field=_FullyTargetedWall(),
            initially_exposed=True,
            terminated_at_end=True,
            observation_ends_at_step_end=True,
        )

        self.assertEqual(int(result.event_start_count[0]), 0)
        self.assertEqual(int(result.event_end_count[0]), 1)
        self.assertFalse(bool(result.exposure_open_at_end[0]))
        self.assertFalse(bool(result.right_censored_at_end[0]))

    def test_batch_preserves_permanent_id_mapping_and_applicability_weight(self) -> None:
        result = integrate_internal_target_exposure(
            permanent_ids=np.asarray([42, 5]),
            path_points_xz_um=np.asarray([[0.0, 1.1], [0.0, 1.4]]),
            path_point_offsets=np.asarray([0, 1, 2]),
            dt_s=np.asarray([0.2, 0.3]),
            radius_um=np.asarray([1.0, 1.0]),
            capture_distance_um=0.2,
            boundary_geometry=_PlanarWallGeometry(),
            molecular_target_field=_FullyTargetedWall(),
            quantitatively_applicable=np.asarray([True, False]),
        )

        np.testing.assert_array_equal(result.permanent_ids, np.asarray([42, 5]))
        np.testing.assert_allclose(result.exposure_time_s, np.asarray([0.2, 0.0]))
        np.testing.assert_allclose(
            result.quantitatively_applicable_exposure_time_s,
            np.asarray([0.2, 0.0]),
        )

    def test_material_wall_penetration_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "materially negative true wall gap"):
            integrate_internal_target_exposure(
                permanent_ids=np.asarray([1]),
                path_points_xz_um=np.asarray([[0.0, 0.9]]),
                path_point_offsets=np.asarray([0, 1]),
                dt_s=0.1,
                radius_um=1.0,
                capture_distance_um=0.2,
                boundary_geometry=_PlanarWallGeometry(),
                molecular_target_field=_FullyTargetedWall(),
            )

    def test_coordinate_roundoff_negative_gap_is_treated_as_contact(self) -> None:
        result = integrate_internal_target_exposure(
            permanent_ids=np.asarray([1]),
            path_points_xz_um=np.asarray([[3000.0, 1.0 - 1.0e-10]]),
            path_point_offsets=np.asarray([0, 1]),
            dt_s=0.1,
            radius_um=1.0,
            capture_distance_um=0.2,
            boundary_geometry=_PlanarWallGeometry(),
            molecular_target_field=_FullyTargetedWall(),
        )

        self.assertAlmostEqual(float(result.exposure_time_s[0]), 0.1)
        self.assertGreater(float(result.reaction_area_time_um2_s[0]), 0.0)

    def test_retry_rounds_do_not_leave_holes_in_stochastic_exposure_polyline(self) -> None:
        event_code = np.zeros((4, 1), dtype=np.uint8)
        event_code[1:3, 0] = 3
        event_position = np.full((4, 1, 2), np.nan, dtype=np.float64)
        event_position[1, 0] = [0.0, 1.1]
        event_position[2, 0] = [1.0, 1.1]
        sweep = RbcStochasticSweepResult(
            position_grid=np.asarray([[0.0, 1.3]]),
            alive=np.asarray([True]),
            vessel_id=np.asarray([1], dtype=np.int32),
            termination_reason=np.zeros(1, dtype=np.uint8),
            exit_time_s=np.asarray([np.nan]),
            reflection_count=np.asarray([2], dtype=np.int32),
            requested_displacement_um=np.asarray([1.0]),
            event_count=np.asarray([2], dtype=np.int32),
            event_code=event_code,
            event_fraction=np.full((4, 1), np.nan, dtype=np.float64),
            event_position_xz_um=event_position,
            topology_event_bubble_id=np.empty(0, dtype=np.int64),
            topology_event_from_vessel_id=np.empty(0, dtype=np.int32),
            topology_event_to_vessel_id=np.empty(0, dtype=np.int32),
            topology_event_section_index=np.empty(0, dtype=np.int32),
            topology_event_position_xz_um=np.empty((0, 2), dtype=np.float64),
        )
        state = SimpleNamespace(position_grid=np.zeros((8, 2), dtype=np.float64))
        context = SimpleNamespace(
            boundary_geometry=_PlanarWallGeometry(),
            molecular_target_field=_FullyTargetedWall(),
            molecular_capture_distance_um=0.2,
            all_bubble_radii_um=np.ones(8, dtype=np.float64),
            red_blood_cell_network=None,
        )

        _accumulate_internal_target_exposure(
            state,
            context,
            permanent_ids=np.asarray([7]),
            start_position_grid=np.asarray([[0.0, 1.3]]),
            start_vessel_id=np.asarray([1], dtype=np.int32),
            deterministic_end_position_grid=np.asarray([[0.0, 1.3]]),
            final_position_grid=np.asarray([[0.0, 1.3]]),
            active_duration_s=np.asarray([1.0]),
            terminated_at_end=np.asarray([False]),
            stochastic_sweep=sweep,
            stochastic_permanent_ids=np.asarray([7]),
        )

        self.assertGreater(float(state.target_exposure_time_s[7]), 0.0)
        self.assertEqual(int(state.target_exposure_event_count[7]), 1)
        self.assertEqual(int(state.target_exposure_event_end_count[7]), 1)
        self.assertFalse(bool(state.target_exposure_open[7]))

    def test_deterministic_transaction_vertex_is_kept_in_exposure_polyline(self) -> None:
        state = SimpleNamespace(position_grid=np.zeros((8, 2), dtype=np.float64))
        context = SimpleNamespace(
            boundary_geometry=_PlanarWallGeometry(),
            molecular_target_field=_FullyTargetedWall(),
            molecular_capture_distance_um=0.2,
            all_bubble_radii_um=np.ones(8, dtype=np.float64),
            red_blood_cell_network=None,
        )
        deterministic_batches = [
            (
                np.asarray([7], dtype=np.int64),
                (
                    np.asarray([[0.0, 1.2 - 1.0e-7]], dtype=np.float64),
                    np.asarray([[0.0, 1.3]], dtype=np.float64),
                ),
                (
                    np.asarray([True]),
                    np.asarray([True]),
                ),
            )
        ]

        _accumulate_internal_target_exposure(
            state,
            context,
            permanent_ids=np.asarray([7]),
            start_position_grid=np.asarray([[0.0, 1.3]]),
            start_vessel_id=np.asarray([1], dtype=np.int32),
            deterministic_end_position_grid=np.asarray([[0.0, 1.3]]),
            final_position_grid=np.asarray([[0.0, 1.3]]),
            active_duration_s=np.asarray([1.0]),
            terminated_at_end=np.asarray([False]),
            deterministic_path_batches=deterministic_batches,
        )

        self.assertGreater(float(state.target_exposure_time_s[7]), 0.0)
        self.assertEqual(int(state.target_exposure_event_count[7]), 1)
        self.assertEqual(int(state.target_exposure_event_end_count[7]), 1)
        self.assertFalse(bool(state.target_exposure_open[7]))


if __name__ == "__main__":
    unittest.main()
