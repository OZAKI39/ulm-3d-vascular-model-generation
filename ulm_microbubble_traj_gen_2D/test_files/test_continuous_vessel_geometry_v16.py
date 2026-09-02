from __future__ import annotations

import tempfile
import unittest
import warnings
from pathlib import Path
from unittest import mock

import numpy as np
from scipy.spatial import cKDTree

from ulm_microbubble_traj_gen_2D.utils.geometry import continuous_vessel_numba
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    _nearest_segment_rows,
    _segment_segment_distance as _reference_segment_segment_distance,
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_numba import (
    _segment_segment_distance as _accelerated_segment_segment_distance,
)
from ulm_microbubble_traj_gen_2D.utils.core.types import GridDomain
from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_target_field import (
    build_molecular_target_field,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_hydrodynamic_fields import (
    ParticleHydrodynamicFields,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_boundary_support import (
    sample_continuous_boundary_support_to_grid,
)
from ulm_vascular_model_generator.utils.core.models import Vessel
from ulm_microbubble_traj_gen_2D.test_files.hybrid_velocity_fixture import (
    rectangular_hybrid_velocity,
)


class ContinuousVesselGeometryV16Tests(unittest.TestCase):
    def test_straight_channel_has_analytic_distance_normal_and_open_caps(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (10.0, 0.0), 2.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(1.0))

        state = geometry.exact_solid_wall_state_xz_um(
            np.asarray([[5.0, 0.0], [5.0, 1.5]], dtype=np.float64)
        )
        np.testing.assert_allclose(state.distance_um, [2.0, 0.5], atol=1.0e-12)
        np.testing.assert_allclose(
            state.inward_normal_xz,
            [[0.0, -1.0], [0.0, -1.0]],
            atol=1.0e-12,
        )
        self.assertEqual(int(np.sum(geometry.open_section_kind < 0)), 1)
        self.assertEqual(int(np.sum(geometry.open_section_kind > 0)), 1)
        starts_x = geometry.solid_face_start_xz_um[:, 0]
        ends_x = geometry.solid_face_end_xz_um[:, 0]
        self.assertFalse(
            np.any((np.abs(starts_x) < 1.0e-12) & (np.abs(ends_x) < 1.0e-12))
        )
        self.assertFalse(
            np.any(
                (np.abs(starts_x - 10.0) < 1.0e-12) & (np.abs(ends_x - 10.0) < 1.0e-12)
            )
        )

    def test_accelerated_inlet_guard_detects_only_directed_cap_exits(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (10.0, 0.0), 2.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(1.0))
        starts_world = np.asarray(
            [[0.5, 0.0], [0.5, 0.0], [0.5, 3.0]], dtype=np.float64
        )
        ends_world = np.asarray(
            [[-0.5, 0.0], [1.5, 0.0], [-0.5, 3.0]], dtype=np.float64
        )

        crossed = geometry.directed_inlet_crossing_mask_grid_accelerated(
            geometry.world_xz_to_grid(starts_world),
            geometry.world_xz_to_grid(ends_world),
        )

        np.testing.assert_array_equal(crossed, [True, False, False])

    def test_inclined_wall_normal_is_grid_independent(self) -> None:
        direction = np.asarray([3.0, 4.0], dtype=np.float64)
        direction /= np.linalg.norm(direction)
        inward = np.asarray([direction[1], -direction[0]])
        midpoint = np.asarray([6.0, 8.0])
        probe = midpoint + 1.25 * inward
        vessel = _vessel(0, -1, [], (0.0, 0.0), (12.0, 16.0), 2.0)

        fine = build_continuous_vessel_geometry([vessel], _domain(0.5))
        coarse = build_continuous_vessel_geometry([vessel], _domain(2.0))
        fine_state = fine.exact_solid_wall_state_xz_um(probe)
        coarse_state = coarse.exact_solid_wall_state_xz_um(probe)

        self.assertEqual(fine.geometry_hash_sha256, coarse.geometry_hash_sha256)
        self.assertAlmostEqual(float(fine_state.distance_um), 0.75, places=11)
        self.assertAlmostEqual(float(coarse_state.distance_um), 0.75, places=11)
        np.testing.assert_allclose(fine_state.inward_normal_xz, -inward, atol=1.0e-12)
        np.testing.assert_allclose(coarse_state.inward_normal_xz, -inward, atol=1.0e-12)

    def test_whole_chord_sweep_detects_first_wall_contact(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (10.0, 0.0), 2.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(1.0))

        contact = geometry.inspect_swept_solid_wall_path_xz_um(
            np.asarray([5.0, 0.0]),
            np.asarray([5.0, 3.0]),
            0.5,
        )

        self.assertAlmostEqual(float(contact.minimum_gap_um), -0.5, places=12)
        self.assertAlmostEqual(float(contact.first_contact_fraction), 0.5, places=12)
        np.testing.assert_allclose(contact.contact_position_xz_um, [5.0, 1.5])
        self.assertFalse(contact.multiple_wall_contact)

    def test_stationary_sweep_far_from_wall_is_feasible_in_both_backends(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (10.0, 0.0), 2.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(1.0))
        point = np.asarray([5.0, 0.0], dtype=np.float64)
        radius = 0.5

        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            reference = geometry.inspect_swept_solid_wall_path_xz_um(
                point, point, radius
            )
            accelerated = geometry.inspect_swept_solid_wall_paths_xz_um_accelerated(
                point.reshape(1, 2),
                point.reshape(1, 2),
                np.asarray([radius], dtype=np.float64),
            )

        self.assertAlmostEqual(reference.minimum_gap_um, 1.5, places=12)
        self.assertIsNone(reference.first_contact_fraction)
        self.assertIsNone(reference.contact_position_xz_um)
        self.assertEqual(reference.primary_face_index, -1)
        self.assertFalse(reference.multiple_wall_contact)
        gaps, fractions, faces, multiple = accelerated
        self.assertAlmostEqual(float(gaps[0]), 1.5, places=12)
        self.assertTrue(np.isnan(fractions[0]))
        self.assertEqual(int(faces[0]), -1)
        self.assertFalse(bool(multiple[0]))

    def test_stationary_sweep_touching_wall_is_contact_in_both_backends(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (10.0, 0.0), 2.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(1.0))
        point = np.asarray([5.0, 1.5], dtype=np.float64)
        radius = 0.5

        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            reference = geometry.inspect_swept_solid_wall_path_xz_um(
                point, point, radius
            )
            accelerated = geometry.inspect_swept_solid_wall_paths_xz_um_accelerated(
                point.reshape(1, 2),
                point.reshape(1, 2),
                np.asarray([radius], dtype=np.float64),
            )

        self.assertAlmostEqual(reference.minimum_gap_um, 0.0, places=12)
        self.assertEqual(reference.first_contact_fraction, 0.0)
        np.testing.assert_allclose(reference.contact_position_xz_um, point, atol=0.0)
        self.assertGreaterEqual(reference.primary_face_index, 0)
        self.assertFalse(reference.multiple_wall_contact)
        gaps, fractions, faces, multiple = accelerated
        self.assertAlmostEqual(float(gaps[0]), 0.0, places=12)
        self.assertEqual(float(fractions[0]), 0.0)
        self.assertGreaterEqual(int(faces[0]), 0)
        self.assertFalse(bool(multiple[0]))

    def test_segment_distance_handles_zero_and_microscopic_first_segment(self) -> None:
        point = np.asarray([0.5, 2.0], dtype=np.float64)
        wall_start = np.asarray([0.0, 0.0], dtype=np.float64)
        wall_end = np.asarray([1.0, 0.0], dtype=np.float64)
        microscopic_end = np.asarray(
            [np.nextafter(point[0], np.inf), point[1]], dtype=np.float64
        )

        for label, path_end in (
            ("zero", point.copy()),
            ("microscopic", microscopic_end),
        ):
            with self.subTest(length=label), warnings.catch_warnings():
                warnings.simplefilter("error", RuntimeWarning)
                reference = _reference_segment_segment_distance(
                    point, path_end, wall_start, wall_end
                )
                accelerated = _accelerated_segment_segment_distance(
                    float(point[0]),
                    float(point[1]),
                    float(path_end[0]),
                    float(path_end[1]),
                    float(wall_start[0]),
                    float(wall_start[1]),
                    float(wall_end[0]),
                    float(wall_end[1]),
                )
                reverse_reference = _reference_segment_segment_distance(
                    wall_start, wall_end, point, path_end
                )
                reverse_accelerated = _accelerated_segment_segment_distance(
                    float(wall_start[0]),
                    float(wall_start[1]),
                    float(wall_end[0]),
                    float(wall_end[1]),
                    float(point[0]),
                    float(point[1]),
                    float(path_end[0]),
                    float(path_end[1]),
                )

                self.assertAlmostEqual(reference, 2.0, places=15)
                self.assertAlmostEqual(float(accelerated), 2.0, places=15)
                self.assertAlmostEqual(reverse_reference, 2.0, places=15)
                self.assertAlmostEqual(float(reverse_accelerated), 2.0, places=15)

        other_point = np.asarray([0.5, -3.0], dtype=np.float64)
        self.assertAlmostEqual(
            _reference_segment_segment_distance(point, point, other_point, other_point),
            5.0,
            places=15,
        )
        self.assertAlmostEqual(
            _accelerated_segment_segment_distance(
                float(point[0]),
                float(point[1]),
                float(point[0]),
                float(point[1]),
                float(other_point[0]),
                float(other_point[1]),
                float(other_point[0]),
                float(other_point[1]),
            ),
            5.0,
            places=15,
        )

    def test_reference_nearest_query_treats_zero_length_segment_as_point(self) -> None:
        starts = np.asarray([[0.0, 0.0], [10.0, 0.0]], dtype=np.float64)
        ends = np.asarray([[0.0, 0.0], [12.0, 0.0]], dtype=np.float64)
        inward = np.asarray([[0.0, 1.0], [0.0, 1.0]], dtype=np.float64)
        centres = 0.5 * (starts + ends)

        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            distance, normal, primary, projected = _nearest_segment_rows(
                np.asarray([[1.0, 0.0]], dtype=np.float64),
                starts,
                ends,
                inward,
                cKDTree(centres),
            )

        np.testing.assert_allclose(distance, [1.0], atol=0.0)
        np.testing.assert_allclose(normal, [[0.0, 1.0]], atol=0.0)
        np.testing.assert_array_equal(primary, [0])
        np.testing.assert_allclose(projected, [[0.0, 0.0]], atol=0.0)

    def test_microscopic_chord_at_model_scale_does_not_intersect_a_remote_wall(
        self,
    ) -> None:
        point = np.asarray([1751.5290799475301, 1422.9795810707387], dtype=np.float64)
        path_end = np.asarray(
            [np.nextafter(point[0], np.inf), point[1]], dtype=np.float64
        )
        wall_start = np.asarray(
            [1750.8214714601572, 1405.1399054157594], dtype=np.float64
        )
        wall_end = np.asarray(
            [1751.8196869308701, 1405.1653539885551], dtype=np.float64
        )
        expected_distance = 17.81659728888684

        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            reference = _reference_segment_segment_distance(
                point, path_end, wall_start, wall_end
            )
            accelerated = _accelerated_segment_segment_distance(
                float(point[0]),
                float(point[1]),
                float(path_end[0]),
                float(path_end[1]),
                float(wall_start[0]),
                float(wall_start[1]),
                float(wall_end[0]),
                float(wall_end[1]),
            )

        self.assertAlmostEqual(reference, expected_distance, places=12)
        self.assertAlmostEqual(float(accelerated), expected_distance, places=12)

    def test_outlet_crossing_uses_directed_continuous_section(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (10.0, 0.0), 2.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(1.0))

        crossing = geometry.first_outlet_crossing_xz_um(
            np.asarray([9.0, 0.0]), np.asarray([11.0, 0.0])
        )
        reverse = geometry.first_outlet_crossing_xz_um(
            np.asarray([11.0, 0.0]), np.asarray([9.0, 0.0])
        )

        self.assertIsNotNone(crossing)
        self.assertAlmostEqual(crossing.fraction, 0.5, places=12)
        np.testing.assert_allclose(crossing.position_xz_um, [10.0, 0.0])
        self.assertIsNone(reverse)

    def test_outlet_event_accepts_only_roundoff_distance_from_section(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (10.0, 0.0), 2.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(1.0))
        start = np.asarray([9.9, 0.0])
        roundoff_inside = np.asarray([10.0 - 5.0e-14, 0.0])

        reference = geometry.first_outlet_crossing_xz_um(
            start, roundoff_inside
        )
        fractions, indices, positions, _, _ = (
            geometry.first_outlet_crossing_arrays_grid_accelerated(
                geometry.world_xz_to_grid(start).reshape(1, 2),
                geometry.world_xz_to_grid(roundoff_inside).reshape(1, 2),
            )
        )

        self.assertIsNotNone(reference)
        self.assertAlmostEqual(reference.fraction, 1.0, places=12)
        self.assertAlmostEqual(float(fractions[0]), 1.0, places=12)
        self.assertGreaterEqual(int(indices[0]), 0)
        np.testing.assert_allclose(
            geometry.grid_to_world_xz(positions[0]), [10.0 - 5.0e-14, 0.0],
            atol=1.0e-12,
        )
        self.assertIsNone(
            geometry.first_outlet_crossing_xz_um(
                start, np.asarray([10.0 - 1.0e-6, 0.0])
            )
        )

    def test_accessible_domain_shrinks_monotonically_with_radius(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (10.0, 0.0), 2.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(1.0))
        probes = np.asarray([[5.0, 0.0], [5.0, 1.0], [5.0, 1.6]])

        small = geometry.is_accessible_xz_um(probes, 0.25)
        large = geometry.is_accessible_xz_um(probes, 1.25)

        self.assertTrue(np.all(~large | small))
        np.testing.assert_array_equal(small, [True, True, True])
        np.testing.assert_array_equal(large, [True, False, False])

    def test_adjacent_curve_elements_do_not_create_false_multiwall_contact(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (10.0, 0.0), 2.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(1.0))
        vertex = geometry.solid_face_end_xz_um[0]
        result = geometry.multiple_wall_contact_mask_xz_um(vertex, 1.0e-8)
        self.assertFalse(bool(result))

    def test_numba_point_and_sweep_queries_match_reference(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (12.0, 16.0), 2.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(1.0))
        rng = np.random.default_rng(1601)
        points = rng.uniform([-1.0, -2.0], [13.0, 18.0], size=(256, 2))

        reference = geometry.exact_solid_wall_state_xz_um(points)
        accelerated = geometry.exact_solid_wall_state_xz_um_accelerated(points)
        np.testing.assert_allclose(accelerated.distance_um, reference.distance_um, atol=1.0e-12)
        np.testing.assert_allclose(
            accelerated.inward_normal_xz, reference.inward_normal_xz, atol=1.0e-12
        )
        np.testing.assert_array_equal(
            accelerated.primary_face_index, reference.primary_face_index
        )

        starts = rng.uniform([1.0, 1.0], [11.0, 15.0], size=(128, 2))
        ends = starts + rng.normal(0.0, 1.0, size=(128, 2))
        radii = rng.uniform(0.1, 1.0, size=128)
        gaps, fractions, faces, multiple = (
            geometry.inspect_swept_solid_wall_paths_xz_um_accelerated(
                starts, ends, radii
            )
        )
        for lane in range(starts.shape[0]):
            expected = geometry.inspect_swept_solid_wall_path_xz_um(
                starts[lane], ends[lane], float(radii[lane])
            )
            self.assertAlmostEqual(float(gaps[lane]), expected.minimum_gap_um, places=12)
            expected_fraction = expected.first_contact_fraction
            if expected_fraction is None:
                self.assertTrue(np.isnan(fractions[lane]))
            else:
                self.assertAlmostEqual(float(fractions[lane]), expected_fraction, places=12)
            self.assertEqual(int(faces[lane]), expected.primary_face_index)
            self.assertEqual(bool(multiple[lane]), expected.multiple_wall_contact)

    @unittest.skipUnless(
        continuous_vessel_numba.NUMBA_AVAILABLE,
        "Numba is required to test serial/parallel exact-state dispatch.",
    )
    def test_exact_state_parallel_threshold_is_bitwise_equivalent_and_routed(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (12.0, 16.0), 2.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(1.0))
        threshold = continuous_vessel_numba._EXACT_STATE_PARALLEL_MIN_BATCH
        self.assertGreater(threshold, 1)
        rng = np.random.default_rng(1621)
        points = np.ascontiguousarray(
            rng.uniform(
                [-1.0, -2.0],
                [13.0, 18.0],
                size=(threshold + 1, 2),
            )
        )
        common = (
            geometry.solid_face_start_xz_um,
            geometry.solid_face_end_xz_um,
            geometry.solid_face_inward_normal_xz,
            geometry.solid_face_ring_index,
            int(geometry.full_boundary_segment_count),
            geometry._exact_bin_edge_origin_xz_um,
            float(geometry._exact_bin_size_um),
            geometry._exact_bin_shape,
            geometry._exact_bin_offsets,
            geometry._exact_bin_segment_indices,
        )

        serial = continuous_vessel_numba._exact_continuous_wall_states_numba(
            points, *common
        )
        parallel = (
            continuous_vessel_numba._exact_continuous_wall_states_parallel_numba(
                points, *common
            )
        )
        for serial_field, parallel_field in zip(serial, parallel):
            np.testing.assert_array_equal(parallel_field, serial_field)

        serial_dispatcher = (
            continuous_vessel_numba._exact_continuous_wall_states_numba
        )
        parallel_dispatcher = (
            continuous_vessel_numba._exact_continuous_wall_states_parallel_numba
        )
        calls: list[tuple[str, int]] = []

        def tracked_serial(*args):
            calls.append(("serial", int(args[0].shape[0])))
            return serial_dispatcher(*args)

        def tracked_parallel(*args):
            calls.append(("parallel", int(args[0].shape[0])))
            return parallel_dispatcher(*args)

        with (
            mock.patch.object(
                continuous_vessel_numba,
                "_exact_continuous_wall_states_numba",
                side_effect=tracked_serial,
            ),
            mock.patch.object(
                continuous_vessel_numba,
                "_exact_continuous_wall_states_parallel_numba",
                side_effect=tracked_parallel,
            ),
        ):
            below = continuous_vessel_numba.exact_continuous_wall_states(
                points[: threshold - 1], *common, use_numba=True
            )
            at_threshold = continuous_vessel_numba.exact_continuous_wall_states(
                points[:threshold], *common, use_numba=True
            )

        self.assertEqual(
            calls,
            [("serial", threshold - 1), ("parallel", threshold)],
        )
        for actual, expected in zip(below, serial):
            np.testing.assert_array_equal(actual, expected[: threshold - 1])
        for actual, expected in zip(at_threshold, serial):
            np.testing.assert_array_equal(actual, expected[:threshold])

    def test_fused_endpoint_distance_and_parallel_sweep_are_bitwise_equivalent(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (12.0, 16.0), 2.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(1.0))
        rng = np.random.default_rng(1618)
        starts = rng.uniform([1.0, 1.0], [11.0, 15.0], size=(128, 2))
        ends = starts + rng.normal(0.0, 0.25, size=(128, 2))
        radii = rng.uniform(0.1, 1.0, size=128)
        start_distance = np.asarray(
            geometry.exact_solid_wall_state_xz_um_accelerated(starts).distance_um
        )
        end_distance = np.asarray(
            geometry.exact_solid_wall_state_xz_um_accelerated(ends).distance_um
        )
        expected_sweep = geometry.inspect_swept_solid_wall_paths_xz_um_precomputed(
            starts, ends, radii, start_distance, end_distance
        )
        fused = geometry.inspect_swept_solid_wall_paths_with_end_distance_xz_um_precomputed(
            starts, ends, radii, start_distance
        )

        np.testing.assert_array_equal(fused[0], end_distance)
        for actual, expected in zip(fused[1:], expected_sweep):
            np.testing.assert_array_equal(actual, expected)

        serial_chunks = [
            geometry.inspect_swept_solid_wall_paths_with_end_distance_xz_um_precomputed(
                starts[first : first + 32],
                ends[first : first + 32],
                radii[first : first + 32],
                start_distance[first : first + 32],
            )
            for first in range(0, starts.shape[0], 32)
        ]
        for field_index, parallel_value in enumerate(fused):
            serial_value = np.concatenate(
                [chunk[field_index] for chunk in serial_chunks], axis=0
            )
            np.testing.assert_array_equal(parallel_value, serial_value)

    def test_one_smooth_wall_cluster_is_not_reported_as_multiwall(self) -> None:
        vessels = [
            _vessel(0, -1, [1, 2], (0.0, 0.0), (10.0, 0.0), 3.0),
            _vessel(1, 0, [], (10.0, 0.0), (20.0, 8.0), 1.5),
            _vessel(2, 0, [], (10.0, 0.0), (20.0, -8.0), 2.0),
        ]
        geometry = build_continuous_vessel_geometry(vessels, _domain(0.5))
        face = int(np.argmin(np.linalg.norm(
            geometry.solid_face_center_xz_um - np.asarray([10.0, 3.0]), axis=1
        )))
        radius = 0.8
        point = (
            geometry.solid_face_center_xz_um[face]
            - radius * geometry.solid_face_outward_normal_xz[face]
        )
        self.assertFalse(
            bool(geometry.multiple_wall_contact_mask_xz_um(point, radius, tolerance_um=1.0e-8))
        )
        gaps, fractions, _, multiple = (
            geometry.inspect_swept_solid_wall_paths_xz_um_accelerated(
                point.reshape(1, 2),
                (point + np.asarray([0.01, 0.0])).reshape(1, 2),
                np.asarray([radius]),
            )
        )
        self.assertGreaterEqual(float(gaps[0]), -1.0e-8)
        self.assertEqual(float(fractions[0]), 0.0)
        self.assertFalse(bool(multiple[0]))

    def test_genuine_opposite_wall_contact_remains_multiwall(self) -> None:
        vessel = _vessel(0, -1, [], (0.0, 0.0), (10.0, 0.0), 1.0)
        geometry = build_continuous_vessel_geometry([vessel], _domain(0.5))
        self.assertTrue(
            bool(
                geometry.multiple_wall_contact_mask_xz_um(
                    np.asarray([5.0, 0.0]), 1.0, tolerance_um=1.0e-10
                )
            )
        )

    def test_continuous_target_npz_is_hash_guarded_and_uses_wall_arclength(self) -> None:
        domain = _domain(1.0)
        geometry = build_continuous_vessel_geometry(
            [_vessel(0, -1, [], (0.0, 0.0), (10.0, 0.0), 2.0)], domain
        )
        fields = _hydrodynamic_fields(domain, geometry)
        with tempfile.TemporaryDirectory() as directory:
            target_path = Path(directory) / "continuous_target.npz"
            np.savez_compressed(
                target_path,
                target_geometry_schema=np.asarray(
                    "v16_continuous_wall_arclength_target"
                ),
                continuous_geometry_hash_sha256=np.asarray(
                    geometry.geometry_hash_sha256
                ),
                wall_ring_index=geometry.solid_face_ring_index,
                wall_arclength_start_um=geometry.solid_face_arclength_start_um,
                wall_arclength_end_um=geometry.solid_face_arclength_end_um,
                wall_target_positive=np.ones(
                    geometry.solid_face_length_um.size, dtype=bool
                ),
            )
            target = build_molecular_target_field(
                domain,
                fields,
                {
                    "enabled": True,
                    "region_mode": "continuous_wall_npz",
                    "mask_npz_path": target_path,
                    "target_density_molecules_per_m2": 1.0e15,
                },
            )
            area = target.reaction_area_um2(
                np.asarray([[5.0, 1.5]]),
                np.asarray([[1.0, 0.0]]),
                np.asarray([0.4]),
            )
            self.assertAlmostEqual(float(area[0]), np.pi * 0.4**2, places=12)
            self.assertTrue(np.all(target.wall_target_positive))

            with np.load(target_path, allow_pickle=False) as saved:
                payload = {name: saved[name] for name in saved.files}
            payload["continuous_geometry_hash_sha256"] = np.asarray("BAD")
            np.savez_compressed(target_path, **payload)
            with self.assertRaisesRegex(ValueError, "geometry hash"):
                build_molecular_target_field(
                    domain,
                    fields,
                    {
                        "enabled": True,
                        "region_mode": "continuous_wall_npz",
                        "mask_npz_path": target_path,
                        "target_density_molecules_per_m2": 1.0e15,
                    },
                )


def _domain(spacing_um: float) -> GridDomain:
    x = np.arange(-4.0, 15.0 + 0.5 * spacing_um, spacing_um)
    z = np.arange(-5.0, 21.0 + 0.5 * spacing_um, spacing_um)
    return GridDomain(
        origin_um=np.asarray([x[0], 0.0, z[0]], dtype=np.float64),
        spacing_um=float(spacing_um),
        shape=(x.size, z.size),
        fixed_y_um=0.0,
        x_coordinates_um=x,
        z_coordinates_um=z,
    )


def _vessel(
    vid: int,
    parent_id: int,
    children: list[int],
    start_xz: tuple[float, float],
    end_xz: tuple[float, float],
    radius_um: float,
) -> Vessel:
    vessel = Vessel(vid=vid, parent_id=parent_id, children=list(children))
    vessel.x_p = np.asarray([start_xz[0], 0.0, start_xz[1]], dtype=np.float64)
    vessel.x_d = np.asarray([end_xz[0], 0.0, end_xz[1]], dtype=np.float64)
    vessel.radius = float(radius_um)
    return vessel


def _hydrodynamic_fields(
    domain: GridDomain, geometry
) -> ParticleHydrodynamicFields:
    shape = domain.shape
    solid_sites, open_boundary = sample_continuous_boundary_support_to_grid(
        domain, geometry
    )
    return ParticleHydrodynamicFields(
        velocity_gradient_s_inv=np.zeros((*shape, 2, 2), dtype=np.float32),
        dynamic_viscosity_pa_s=np.ones(shape, dtype=np.float32),
        solid_site_mask=solid_sites,
        open_boundary_mask=open_boundary,
        boundary_geometry=geometry,
        hybrid_velocity=rectangular_hybrid_velocity(domain),
    )


if __name__ == "__main__":
    unittest.main()
