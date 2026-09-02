from __future__ import annotations

import unittest
import numpy as np

from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.test_files.particle_fixtures import (
    straight_channel_case,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.core.types import GridDomain
from ulm_microbubble_traj_gen_2D.utils.particles.particle_counter_rng import (
    COUNTER_RNG_ALGORITHM,
    counter_hashes,
    counter_normal,
    counter_normal_batch,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_rbc_stochastic import (
    RbcStochasticGeometryError,
    sweep_rbc_stochastic_displacement,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_topological_ownership import (
    build_topological_commitment_catalog,
)
from ulm_microbubble_traj_gen_2D.test_files.test_particle_topological_ownership_v20 import (
    _branching_tree,
)


class RbcCounterRandomNumberTests(unittest.TestCase):
    def test_fixed_key_has_versioned_known_answer(self) -> None:
        self.assertEqual(COUNTER_RNG_ALGORITHM, "splitmix64_box_muller_v1")
        self.assertEqual(
            counter_hashes(42, 7, 19, 0),
            (15964014261161886031, 8341215932289832289),
        )
        self.assertAlmostEqual(
            counter_normal(42, 7, 19, 0),
            -0.5135917105267578,
            places=15,
        )

    def test_permanent_id_sequence_is_independent_of_lane_order(self) -> None:
        ids = np.asarray([11, 2, 99, 7, 34], dtype=np.int64)
        permutation = np.asarray([3, 0, 4, 1, 2], dtype=np.int64)
        reference = counter_normal_batch(
            731, ids, 51, component=0, use_numba=False
        )
        reordered = counter_normal_batch(
            731, ids[permutation], 51, component=0, use_numba=False
        )

        restored = np.empty_like(reordered)
        restored[permutation] = reordered
        np.testing.assert_array_equal(restored, reference)

    def test_different_seed_changes_particle_sequence(self) -> None:
        ids = np.arange(64, dtype=np.int64)
        first = counter_normal_batch(731, ids, 51, use_numba=False)
        second = counter_normal_batch(732, ids, 51, use_numba=False)

        self.assertFalse(np.array_equal(first, second))
        self.assertTrue(np.all(first != second))

    def test_compiled_and_reference_counter_normals_are_bitwise_equal(self) -> None:
        ids = np.asarray([0, 1, 2, 7, 99, (1 << 31) - 1], dtype=np.int64)
        reference = counter_normal_batch(
            731, ids, 51, component=0, use_numba=False
        )
        compiled = counter_normal_batch(
            731, ids, 51, component=0, use_numba=True
        )

        np.testing.assert_array_equal(compiled, reference)

    def test_constant_diffusivity_increment_has_expected_mean_and_variance(self) -> None:
        particle_ids = np.arange(100_000, dtype=np.int64)
        normal = counter_normal_batch(
            9281, particle_ids, 3, component=0, use_numba=False
        )
        diffusivity_um2_s = 45.0
        elapsed_s = 1.0e-3
        displacement_um = np.sqrt(
            2.0 * diffusivity_um2_s * elapsed_s
        ) * normal

        self.assertAlmostEqual(float(np.mean(displacement_um)), 0.0, delta=2.0e-3)
        self.assertAlmostEqual(
            float(np.var(displacement_um)),
            2.0 * diffusivity_um2_s * elapsed_s,
            delta=2.0e-3,
        )


class RbcStochasticGeometryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        domain, _, _ = straight_channel_case()
        vessel = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([0.0, 0.0, 3.0]),
            x_d=np.asarray([19.0, 0.0, 3.0]),
            radius=3.0,
            flow_rate=50.0,
        )
        continuous = build_continuous_vessel_geometry([vessel], domain)
        cls.geometry = continuous
        cls.topology = build_topological_commitment_catalog([vessel], continuous)

    def test_flat_wall_reflects_remaining_normal_chord_without_freezing(self) -> None:
        result = self._sweep(
            start_grid=np.asarray([[5.0, 3.0]]),
            displacement_um=np.asarray([[2.0, 4.0]]),
        )

        self.assertTrue(result.alive[0])
        self.assertEqual(result.vessel_id[0], 1)
        self.assertEqual(result.reflection_count[0], 1)
        self.assertEqual(result.event_count[0], 1)
        self.assertEqual(result.event_code[0, 0], 3)
        # The X component is tangential to this wall and is unchanged.  The
        # uncompleted Z component is mirrored, so the particle finishes at Z=2.
        np.testing.assert_allclose(result.position_grid[0], [7.0, 3.0], atol=2.0e-12)
        self.assertNotEqual(float(result.position_grid[0, 1]), 4.5)

        final_world = self.geometry.grid_to_world_xz(result.position_grid)
        exact = self.geometry.exact_solid_wall_state_xz_um(final_world)
        final_gap_um = float(exact.distance_um[0]) - 1.0
        self.assertGreaterEqual(final_gap_um, -1.0e-12)
        self.assertFalse(hasattr(result, "mechanical_contact_force"))

    def test_contact_start_moving_inward_has_no_zero_time_reflection(self) -> None:
        result = self._sweep(
            start_grid=np.asarray([[5.0, 4.5]]),
            displacement_um=np.asarray([[2.0, -1.0]]),
        )

        self.assertTrue(result.alive[0])
        self.assertEqual(result.reflection_count[0], 0)
        self.assertEqual(result.event_count[0], 0)
        np.testing.assert_allclose(result.position_grid[0], [7.0, 3.5], atol=2.0e-12)

    def test_roundoff_negative_final_gap_is_projected_strictly_inside(self) -> None:
        result = self._sweep(
            start_grid=np.asarray([[5.0, 5.0 + 1.0e-12]]),
            displacement_um=np.asarray([[0.0, 0.0]]),
        )

        final_world = self.geometry.grid_to_world_xz(result.position_grid)
        exact = self.geometry.exact_solid_wall_state_xz_um(final_world)
        final_gap_um = float(exact.distance_um[0]) - 1.0
        self.assertGreaterEqual(final_gap_um, 0.0)
        self.assertEqual(result.reflection_count[0], 0)

    def test_subresolution_random_remainder_does_not_exhaust_event_limit(self) -> None:
        # Regression for bubble 259: a 1e-23 um reflected tail remained
        # bitwise nonzero and used all 32 geometry rounds. Both the requested
        # motion and the final negative gap are below coordinate resolution.
        result = self._sweep(
            start_grid=np.asarray([[5.0, 5.0 + 3.7e-13]]),
            displacement_um=np.asarray([[0.0, 3.9e-13]]),
        )

        final_world = self.geometry.grid_to_world_xz(result.position_grid)
        exact = self.geometry.exact_solid_wall_state_xz_um(final_world)
        final_gap_um = float(exact.distance_um[0]) - 1.0
        self.assertGreaterEqual(final_gap_um, 0.0)
        self.assertTrue(result.alive[0])
        self.assertEqual(result.event_count[0], 0)
        self.assertEqual(result.reflection_count[0], 0)

    def test_near_wall_roundoff_stress_batch_is_strictly_feasible(self) -> None:
        count = 96
        magnitudes = np.geomspace(1.0e-24, 1.0e-7, count)
        signs = np.where(np.arange(count) % 2, 1.0, -1.0)
        offsets = np.where(
            np.arange(count) % 3 == 0,
            3.0e-13,
            np.where(np.arange(count) % 3 == 1, -3.0e-13, 0.0),
        )
        positions = np.column_stack(
            (np.full(count, 5.0), 5.0 + offsets)
        )
        displacement = np.column_stack(
            (
                0.25
                * magnitudes
                * np.where(np.arange(count) % 4 < 2, 1.0, -1.0),
                signs * magnitudes,
            )
        )

        for use_numba in (False, True):
            with self.subTest(use_numba=use_numba):
                result = sweep_rbc_stochastic_displacement(
                    positions_grid=positions,
                    displacement_xz_um=displacement,
                    permanent_ids=np.arange(count, dtype=np.int64),
                    alive=np.ones(count, dtype=bool),
                    vessel_id=np.ones(count, dtype=np.int32),
                    radii_um=np.ones(count, dtype=np.float64),
                    termination_reason=np.zeros(count, dtype=np.uint8),
                    exit_time_s=np.full(count, np.nan),
                    physical_time_s=0.125,
                    boundary_geometry=self.geometry,
                    topological_ownership=self.topology,
                    use_numba=use_numba,
                )

                final_world = self.geometry.grid_to_world_xz(
                    result.position_grid
                )
                exact = self.geometry.exact_solid_wall_state_xz_um(final_world)
                final_gap_um = (
                    np.asarray(exact.distance_um, dtype=np.float64) - 1.0
                )
                self.assertTrue(np.all(result.alive))
                self.assertTrue(np.all(final_gap_um >= 0.0))
                self.assertLessEqual(int(np.max(result.event_count)), 1)

    def test_large_world_coordinates_reproduce_bubble_259_scale_safely(self) -> None:
        domain = GridDomain(
            origin_um=np.asarray([2830.0, 0.0, 1270.0]),
            spacing_um=1.0,
            shape=(20, 7),
            fixed_y_um=0.0,
            x_coordinates_um=np.arange(2830.0, 2850.0),
            z_coordinates_um=np.arange(1270.0, 1277.0),
        )
        vessel = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([2830.0, 0.0, 1273.0]),
            x_d=np.asarray([2849.0, 0.0, 1273.0]),
            radius=3.0,
            flow_rate=50.0,
        )
        geometry = build_continuous_vessel_geometry([vessel], domain)
        topology = build_topological_commitment_catalog([vessel], geometry)
        result = sweep_rbc_stochastic_displacement(
            positions_grid=np.asarray([[5.0, 5.0 + 3.7e-13]]),
            displacement_xz_um=np.asarray(
                [[-3.437150835301031e-13, -1.7518758943204562e-13]]
            ),
            permanent_ids=np.asarray([259], dtype=np.int64),
            alive=np.asarray([True]),
            vessel_id=np.asarray([1], dtype=np.int32),
            radii_um=np.asarray([1.0]),
            termination_reason=np.asarray([0], dtype=np.uint8),
            exit_time_s=np.asarray([np.nan]),
            physical_time_s=0.441,
            boundary_geometry=geometry,
            topological_ownership=topology,
            use_numba=True,
        )

        final_world = geometry.grid_to_world_xz(result.position_grid)
        final_gap_um = (
            float(geometry.exact_solid_wall_state_xz_um(final_world).distance_um[0])
            - 1.0
        )
        self.assertTrue(result.alive[0])
        self.assertGreaterEqual(final_gap_um, 0.0)
        self.assertEqual(result.event_count[0], 0)

    def test_material_negative_final_gap_remains_a_geometry_error(self) -> None:
        with self.assertRaises(RbcStochasticGeometryError) as raised:
            self._sweep(
                start_grid=np.asarray([[5.0, 5.0 + 1.0e-6]]),
                displacement_um=np.asarray([[0.0, 0.0]]),
            )

        self.assertEqual(
            raised.exception.failures[0].permanent_microbubble_id, 71
        )
        self.assertLess(raised.exception.failures[0].wall_gap_um, 0.0)

    def test_earlier_wall_event_precedes_later_outlet_event(self) -> None:
        result = self._sweep(
            start_grid=np.asarray([[15.0, 3.0]]),
            displacement_um=np.asarray([[6.0, 4.0]]),
        )

        self.assertFalse(result.alive[0])
        self.assertEqual(result.reflection_count[0], 1)
        self.assertEqual(result.event_code[0, 0], 3)
        self.assertEqual(result.event_code[1, 0], 1)

    def test_reverse_inlet_escape_is_rejected(self) -> None:
        with self.assertRaises(RbcStochasticGeometryError) as raised:
            self._sweep(
                start_grid=np.asarray([[1.0, 3.0]]),
                displacement_um=np.asarray([[-3.0, 0.0]]),
            )

        self.assertEqual(raised.exception.failures[0].permanent_microbubble_id, 71)
        self.assertEqual(raised.exception.failures[0].event_codes, ())
        self.assertEqual(
            raised.exception.failures[0].requested_displacement_xz_um,
            (-3.0, 0.0),
        )

    def test_wall_before_reverse_inlet_is_processed_chronologically(self) -> None:
        with self.assertRaises(RbcStochasticGeometryError) as raised:
            self._sweep(
                start_grid=np.asarray([[4.0, 3.0]]),
                displacement_um=np.asarray([[-6.0, 4.0]]),
            )

        # The unreﬂected chord intersects the inlet only after the upper wall.
        # The reflected path is processed first and the failure diagnostic
        # therefore retains the actual wall event.
        self.assertEqual(raised.exception.failures[0].event_codes, (3,))

    def test_outward_open_outlet_crossing_terminates_without_reflection(self) -> None:
        result = self._sweep(
            start_grid=np.asarray([[18.0, 3.0]]),
            displacement_um=np.asarray([[3.0, 0.0]]),
        )

        self.assertFalse(result.alive[0])
        self.assertEqual(result.termination_reason[0], 1)
        self.assertEqual(result.exit_time_s[0], 0.125)
        self.assertEqual(result.reflection_count[0], 0)
        self.assertEqual(result.event_count[0], 1)
        self.assertEqual(result.event_code[0, 0], 1)
        self.assertEqual(result.vessel_id[0], 1)
        np.testing.assert_allclose(result.position_grid[0], [19.0, 3.0])

    def test_random_path_changes_owner_only_at_continuous_commitment_section(self) -> None:
        vessels = _branching_tree()
        domain = GridDomain(
            origin_um=np.asarray([-5.0, 0.0, -20.0]),
            spacing_um=1.0,
            shape=(56, 41),
            fixed_y_um=0.0,
            x_coordinates_um=np.arange(-5.0, 51.0),
            z_coordinates_um=np.arange(-20.0, 21.0),
        )
        continuous = build_continuous_vessel_geometry(vessels, domain)
        topology = build_topological_commitment_catalog(vessels, continuous)
        point = topology.point_xz_um[0]
        normal = topology.downstream_normal_xz[0]
        start = continuous.world_xz_to_grid((point - normal)[None, :])
        end = continuous.world_xz_to_grid((point + normal)[None, :])

        result = sweep_rbc_stochastic_displacement(
            positions_grid=start,
            displacement_xz_um=end - start,
            permanent_ids=np.asarray([55], dtype=np.int64),
            alive=np.asarray([True]),
            vessel_id=np.asarray([1], dtype=np.int32),
            radii_um=np.asarray([0.1]),
            termination_reason=np.asarray([0], dtype=np.uint8),
            exit_time_s=np.asarray([np.nan]),
            physical_time_s=0.25,
            boundary_geometry=continuous,
            topological_ownership=topology,
            use_numba=False,
        )

        self.assertTrue(result.alive[0])
        self.assertEqual(result.vessel_id[0], 2)
        self.assertEqual(result.event_code[0, 0], 2)
        np.testing.assert_allclose(result.position_grid[0], end[0])
        np.testing.assert_array_equal(result.topology_event_bubble_id, [55])
        np.testing.assert_array_equal(result.topology_event_from_vessel_id, [1])
        np.testing.assert_array_equal(result.topology_event_to_vessel_id, [2])

    def _sweep(
        self,
        *,
        start_grid: np.ndarray,
        displacement_um: np.ndarray,
    ):
        return sweep_rbc_stochastic_displacement(
            positions_grid=start_grid,
            displacement_xz_um=displacement_um,
            permanent_ids=np.asarray([71], dtype=np.int64),
            alive=np.asarray([True]),
            vessel_id=np.asarray([1], dtype=np.int32),
            radii_um=np.asarray([1.0]),
            termination_reason=np.asarray([0], dtype=np.uint8),
            exit_time_s=np.asarray([np.nan]),
            physical_time_s=0.125,
            boundary_geometry=self.geometry,
            topological_ownership=self.topology,
            use_numba=False,
        )
if __name__ == "__main__":
    unittest.main()
