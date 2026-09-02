from __future__ import annotations

from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.particles import particle_perfusion_transport as transport
from ulm_microbubble_traj_gen_2D.utils.particles.particle_predictive_contact import (
    PredictiveContactStep,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_constrained_step import (
    BatchLocalState,
    PhysicalTimeInterval,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_topological_ownership import (
    TopologicalCommitmentCatalog,
)


class _BoundaryGeometry:
    def __init__(self, outlet_x: float | None = None, maximum_x: float = 100.0):
        self.outlet_x = outlet_x
        self.maximum_x = maximum_x

    @staticmethod
    def world_xz_to_grid(position):
        return np.asarray(position, dtype=np.float64)

    @staticmethod
    def grid_to_world_xz(position):
        return np.asarray(position, dtype=np.float64)

    @staticmethod
    def exact_solid_wall_state_xz_um(positions):
        points = np.asarray(positions, dtype=np.float64)
        scalar = points.shape == (2,)
        points = points.reshape(-1, 2)
        lower = points[:, 1]
        upper = 5.0 - points[:, 1]
        use_lower = lower <= upper
        distance = np.where(use_lower, lower, upper)
        normal = np.zeros_like(points)
        normal[:, 1] = np.where(use_lower, 1.0, -1.0)
        primary = np.where(use_lower, 0, 1)
        if scalar:
            return SimpleNamespace(
                distance_um=np.asarray(distance[0]),
                inward_normal_xz=normal[0],
                primary_face_index=np.asarray(primary[0]),
            )
        return SimpleNamespace(
            distance_um=distance,
            inward_normal_xz=normal,
            primary_face_index=primary,
        )

    def exact_true_gap_at_xz_um(self, positions, radius_um):
        return np.asarray(
            self.exact_solid_wall_state_xz_um(positions).distance_um
        ) - np.asarray(radius_um)

    def inspect_swept_solid_wall_path_xz_um(
        self, start, end, radius_um, *, tolerance_um=0.0
    ):
        del tolerance_um
        start = np.asarray(start, dtype=np.float64)
        end = np.asarray(end, dtype=np.float64)
        start_gap = float(self.exact_true_gap_at_xz_um(start, radius_um))
        end_gap = float(self.exact_true_gap_at_xz_um(end, radius_um))
        fraction = None
        if start_gap <= 0.0:
            fraction = 0.0
        elif end_gap <= 0.0:
            fraction = start_gap / (start_gap - end_gap)
        return SimpleNamespace(
            minimum_gap_um=min(start_gap, end_gap),
            first_contact_fraction=fraction,
            contact_position_xz_um=(
                None if fraction is None else start + fraction * (end - start)
            ),
            multiple_wall_contact=False,
        )

    @staticmethod
    def multiple_wall_contact_mask_xz_um(
        positions, radius_um, *, tolerance_um=1.0e-9
    ):
        del radius_um, tolerance_um
        points = np.asarray(positions).reshape(-1, 2)
        result = np.zeros(points.shape[0], dtype=bool)
        return bool(result[0]) if np.asarray(positions).shape == (2,) else result

    def first_outlet_crossing_grid(self, start, end):
        if self.outlet_x is None:
            return None
        start = np.asarray(start, dtype=np.float64)
        end = np.asarray(end, dtype=np.float64)
        delta = float(end[0] - start[0])
        if delta <= 0.0 or start[0] > self.outlet_x or end[0] < self.outlet_x:
            return None
        fraction = float((self.outlet_x - start[0]) / delta)
        return SimpleNamespace(
            fraction=fraction,
            label=1,
            position_xz_um=start + fraction * (end - start),
        )

    def first_outlet_crossings_grid(self, starts, ends):
        return tuple(
            self.first_outlet_crossing_grid(start, end)
            for start, end in zip(
                np.asarray(starts, dtype=np.float64),
                np.asarray(ends, dtype=np.float64),
                strict=True,
            )
        )

    def contains_grid(self, positions, *, tolerance_um=0.0):
        points = np.asarray(positions, dtype=np.float64)
        scalar = points.shape == (2,)
        points = points.reshape(-1, 2)
        result = (
            (points[:, 0] >= -tolerance_um)
            & (points[:, 0] <= self.maximum_x + tolerance_um)
            & (points[:, 1] >= -tolerance_um)
            & (points[:, 1] <= 5.0 + tolerance_um)
        )
        return bool(result[0]) if scalar else result

    @staticmethod
    def is_accessible_grid(positions, radius_um):
        points = np.asarray(positions)
        result = np.ones(points.shape[:-1] or (), dtype=bool)
        return bool(result) if result.shape == () else result


def _context(
    particle_count: int,
    shape: tuple[int, int] = (12, 6),
    *,
    outlet_x: float | None = None,
    maximum_x: float = 100.0,
) -> SimpleNamespace:
    return SimpleNamespace(
        all_bubble_radii_um=np.ones(particle_count, dtype=np.float64),
        spacing_um=1.0,
        use_numba=False,
        boundary_geometry=_BoundaryGeometry(outlet_x, maximum_x),
    )


def _particle_config(maximum_refinements: int = 8) -> SimpleNamespace:
    return SimpleNamespace(
        contact_geometry_tolerance_um=1.0e-6,
        contact_max_time_refinements=maximum_refinements,
        wall_contact_threshold_um=0.05,
    )


def _state(positions: np.ndarray) -> transport._PerfusionState:
    particle_count = int(positions.shape[0])
    return transport._PerfusionState(
        position_grid=np.asarray(positions, dtype=np.float64).copy(),
        active=np.ones(particle_count, dtype=bool),
        rotation_angle_rad=np.zeros(particle_count, dtype=np.float64),
        bond_count_expected=np.zeros(particle_count, dtype=np.float64),
        bond_total_tangential_extension_um=np.zeros(
            particle_count, dtype=np.float64
        ),
        admission_time_s=np.zeros(particle_count, dtype=np.float64),
        exit_time_s=np.full(particle_count, np.nan, dtype=np.float64),
        termination_reason=np.zeros(particle_count, dtype=np.uint8),
    )


def _rhs(
    velocities_um_s: np.ndarray,
    sampled_times: list[float] | None = None,
    sampled_positions: list[np.ndarray] | None = None,
    sampled_bond_counts: list[np.ndarray] | None = None,
):
    velocities_um_s = np.asarray(velocities_um_s, dtype=np.float64)

    def evaluate(
        positions,
        bubble_ids,
        active,
        context,
        time_s,
        bond_count,
        bond_extension,
    ):
        del bond_extension
        if sampled_times is not None:
            sampled_times.append(float(time_s))
        if sampled_positions is not None:
            sampled_positions.append(np.asarray(positions).copy())
        if sampled_bond_counts is not None:
            sampled_bond_counts.append(np.asarray(bond_count).copy())
        count = len(bubble_ids)
        velocity = velocities_um_s[np.asarray(bubble_ids, dtype=np.int64)].copy()
        velocity[~active] = 0.0
        gaps = np.asarray(
            [
                min(float(position[1]), 5.0 - float(position[1]))
                - context.all_bubble_radii_um[bubble_ids[lane]]
                if active[lane]
                else np.nan
                for lane, position in enumerate(positions)
            ],
            dtype=np.float64,
        )
        normals = np.tile(np.asarray([0.0, 1.0]), (count, 1))
        return SimpleNamespace(
            particle_velocity_xz_um_s=velocity,
            angular_velocity_rad_s=np.zeros(count),
            generalized_mobility=np.tile(np.eye(3), (count, 1, 1)),
            wall_gap_um=gaps,
            wall_normal_xz=normals,
            maximum_physical_overlap_um=0.0,
            maximum_collision_compression_um=0.0,
            interacting_pair_count=0,
            collision_search_strategy="all_pairs",
            maximum_reciprocity_relative_error=0.0,
            degenerate_near_wall_normal_count=0,
            maximum_collision_speed_um_s=0.0,
            two_wall_warning=np.zeros(count, dtype=bool),
            molecular_binding_evaluation=None,
            bond_count_expected=None,
            sampled_wall_shear_stress_pa=np.zeros(count),
            sampled_vessel_id=np.ones(count, dtype=np.int32),
            fluid_velocity_xz_um_s=velocity.copy(),
            collision_force_xz_pn=np.zeros((count, 2)),
            collision_neighbor_count=np.zeros(count, dtype=np.int32),
            gap_ratio=np.zeros(count),
            near_wall_weight=np.zeros(count),
            cardiac_multiplier=np.ones(count),
        )

    return evaluate


class PerfusionPhysicalTimeTransactionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.euler = SimpleNamespace(time_integrator="euler")

    def test_rbc_stochastic_step_reuses_cached_start_coefficients(self) -> None:
        context = SimpleNamespace(
            red_blood_cell_network=SimpleNamespace(),
            all_bubble_radii_um=np.asarray([1.0]),
            use_numba=False,
            random_seed=17,
        )
        state = transport._BatchAdvanceState(
            particles=BatchLocalState(
                position_grid=np.asarray([[2.0, 2.0]]),
                rotation_angle_rad=np.asarray([0.0]),
                bond_count_expected=np.asarray([0.0]),
                bond_total_tangential_extension_um=np.asarray([0.0]),
            ),
            alive=np.asarray([True]),
            termination_reason=np.asarray([0], dtype=np.uint8),
            exit_time_s=np.asarray([np.nan]),
            vessel_id=np.asarray([1], dtype=np.int32),
        )
        cached = SimpleNamespace(
            red_blood_cell_transverse_diffusivity_um2_s=np.asarray([0.0]),
            red_blood_cell_diffusion_enabled=np.asarray([False]),
            red_blood_cell_quantitative_applicability=np.asarray([True]),
            red_blood_cell_transverse_space_valid=np.asarray([True]),
        )

        with mock.patch.object(
            transport,
            "evaluate_red_blood_cell_transport",
            side_effect=AssertionError("cached coefficients must avoid reevaluation"),
        ) as evaluate:
            result = transport._apply_rbc_stochastic_batch_state(
                state,
                np.asarray([0], dtype=np.int64),
                np.asarray([0.1]),
                context,
                global_internal_step=4,
                physical_time_s=0.4,
                coefficient_evaluation=cached,
            )

        evaluate.assert_not_called()
        np.testing.assert_array_equal(
            result.red_blood_cell_quantitative_applicability,
            np.asarray([True]),
        )
        self.assertIsNone(result.stochastic_sweep)

    def test_position_contact_consumes_one_complete_shared_interval(self) -> None:
        context = _context(2)
        state = _state(np.asarray([[2.0, 2.0], [4.0, 2.0]]))
        sampled_times: list[float] = []
        sampled_positions: list[np.ndarray] = []
        sampled_bond_counts: list[np.ndarray] = []
        accepted_bond_durations: list[float] = []

        def accept_bonds(
            bond_count,
            bond_extension,
            evaluation,
            start_positions,
            accepted_positions,
            start_angles,
            accepted_angles,
            radii_um,
            moving,
            spacing_um,
            dt_s,
            use_numba,
        ):
            del (
                evaluation,
                start_positions,
                accepted_positions,
                start_angles,
                accepted_angles,
                radii_um,
                spacing_um,
                use_numba,
            )
            accepted_bond_durations.append(float(dt_s))
            count = np.asarray(bond_count, dtype=np.float64).copy()
            extension = np.asarray(bond_extension, dtype=np.float64).copy()
            count[moving] += float(dt_s)
            extension[moving] += 2.0 * float(dt_s)
            return count, extension, 0

        with mock.patch.object(
            transport, "_accept_euler_bonds", side_effect=accept_bonds
        ):
            transport._advance_segment(
                0.0,
                1.0,
                state,
                context,
                _particle_config(),
                self.euler,
                transport._Diagnostics(),
                _rhs(
                    np.asarray([[0.0, -2.0], [2.0, 0.0]]),
                    sampled_times,
                    sampled_positions,
                    sampled_bond_counts,
                ),
            )

        # Revised v15 solves the predictive contact equation for the complete
        # physical interval. A wall touch is no longer a zero-time lifecycle
        # event and does not split the other synchronous lanes.
        self.assertEqual(len(sampled_times), 1)
        self.assertAlmostEqual(sampled_times[0], 0.0)
        np.testing.assert_allclose(sampled_bond_counts[0], [0.0, 0.0])
        np.testing.assert_allclose(state.position_grid[0], [2.0, 1.0], atol=1.0e-6)
        np.testing.assert_allclose(state.position_grid[1], [6.0, 2.0], atol=1.0e-6)

        self.assertEqual(len(accepted_bond_durations), 1)
        self.assertAlmostEqual(sum(accepted_bond_durations), 1.0, places=12)
        np.testing.assert_allclose(state.bond_count_expected, [1.0, 1.0])
        np.testing.assert_allclose(
            state.bond_total_tangential_extension_um, [2.0, 2.0]
        )

    def test_heun_contact_direction_is_owned_by_geometry_not_rhs_stages(self) -> None:
        context = _context(1)
        ids = np.asarray([0], dtype=np.int64)
        start = transport._BatchAdvanceState(
            particles=BatchLocalState(
                position_grid=np.asarray([[2.0, 2.0]]),
                rotation_angle_rad=np.asarray([0.0]),
                bond_count_expected=np.asarray([0.0]),
                bond_total_tangential_extension_um=np.asarray([0.0]),
            ),
            alive=np.asarray([True]),
            termination_reason=np.asarray([0], dtype=np.uint8),
            exit_time_s=np.asarray([np.nan]),
        )
        base_rhs = _rhs(np.asarray([[0.1, 0.0]]))
        rhs_calls = 0

        def staged_rhs(*args, **kwargs):
            nonlocal rhs_calls
            evaluation = base_rhs(*args, **kwargs)
            evaluation.wall_normal_xz = np.asarray(
                [[1.0, 0.0]] if rhs_calls == 0 else [[0.0, 1.0]],
                dtype=np.float64,
            )
            rhs_calls += 1
            return evaluation

        trial_calls = 0

        def solve_trial(
            positions,
            alive,
            dt_s,
            generalized_velocity,
            mobility,
            radii,
            context,
            tolerance_um,
        ):
            nonlocal trial_calls
            del mobility, radii, context, tolerance_um
            trial_calls += 1
            count = int(len(alive))
            accepted = np.asarray(positions, dtype=np.float64) + (
                float(dt_s)
                * np.asarray(generalized_velocity, dtype=np.float64)[:, :2]
            )
            return PredictiveContactStep(
                accepted_position_grid=accepted,
                angle_increment_rad=(
                    float(dt_s)
                    * np.asarray(generalized_velocity, dtype=np.float64)[:, 2]
                ),
                constrained_generalized_velocity=np.asarray(
                    generalized_velocity, dtype=np.float64
                ).copy(),
                reaction_force_pn=np.zeros(count),
                endpoint_gap_um=np.ones(count),
                minimum_path_gap_um=np.ones(count),
                predicted_normal_gap_um=np.ones(count),
                complementarity_residual_pn_um=np.zeros(count),
                residual_correction_um=np.zeros(count),
                contact_normal_xz=np.tile([0.0, 1.0], (count, 1)),
                active_contact=np.zeros(count, dtype=bool),
                outlet_fraction=np.full(count, np.nan),
                outlet_position_grid=np.full((count, 2), np.nan),
                outlet_label=np.full(count, -1, dtype=np.int32),
                need_time_refinement=np.zeros(count, dtype=bool),
                failure_codes=np.zeros(count, dtype=np.int8),
                failure_position_grid=np.full((count, 2), np.nan),
                failure_event_fraction=np.full(count, np.nan),
            )

        with mock.patch.object(
            transport, "_solve_v15_trial", side_effect=solve_trial
        ):
            transport._attempt_batch_step(
                PhysicalTimeInterval(0.0, 1.0),
                start,
                ids,
                context,
                _particle_config(),
                SimpleNamespace(time_integrator="heun"),
                staged_rhs,
            )

        self.assertEqual(trial_calls, 2)

    def test_failed_trial_is_discarded_before_two_refined_halves_are_committed(self) -> None:
        context = _context(1)
        ids = np.asarray([0], dtype=np.int64)
        start = transport._BatchAdvanceState(
            particles=BatchLocalState(
                position_grid=np.asarray([[2.0, 2.0]]),
                rotation_angle_rad=np.asarray([0.0]),
                bond_count_expected=np.asarray([4.0]),
                bond_total_tangential_extension_um=np.asarray([0.5]),
            ),
            alive=np.asarray([True]),
            termination_reason=np.asarray([0], dtype=np.uint8),
            exit_time_s=np.asarray([np.nan]),
        )
        calls: list[tuple[float, float, float]] = []

        def attempt(interval, start_state, *args, **kwargs):
            del args, kwargs
            calls.append(
                (
                    float(interval.start_time_s),
                    float(interval.end_time_s),
                    float(start_state.particles.bond_count_expected[0]),
                )
            )
            trial = transport._accepted_batch_state(
                start_state,
                start_state.particles.position_grid.copy(),
                start_state.particles.rotation_angle_rad.copy(),
                start_state.particles.bond_count_expected.copy(),
                start_state.particles.bond_total_tangential_extension_um.copy(),
            )
            if interval.refinement_depth == 0:
                # This deliberately imitates a failed speculative trial that
                # had already changed its private position and bond arrays.
                trial.particles.position_grid[0, 0] = 999.0
                trial.particles.bond_count_expected[0] = 999.0
                raise transport._RefineContactTimeStep(
                    "position-level chord failed",
                )

            duration = float(interval.duration_s)
            trial.particles.position_grid[0, 0] += duration
            trial.particles.bond_count_expected[0] += duration
            trial.particles.bond_total_tangential_extension_um[0] += duration
            diagnostics = transport._Diagnostics()
            diagnostics.maximum_step_displacement_um = duration
            return transport._BatchAttemptResult(
                trial,
                diagnostics,
                float(interval.end_time_s),
            )

        with mock.patch.object(transport, "_attempt_batch_step", side_effect=attempt):
            result = transport._advance_batch_interval(
                PhysicalTimeInterval(2.0, 3.0),
                start,
                ids,
                context,
                _particle_config(maximum_refinements=4),
                self.euler,
                object(),
            )

        self.assertEqual(
            [(start_time, end_time) for start_time, end_time, _ in calls],
            [(2.0, 3.0), (2.0, 2.5), (2.5, 3.0)],
        )
        self.assertEqual([bond for _, _, bond in calls], [4.0, 4.0, 4.5])
        self.assertAlmostEqual(result.state.particles.position_grid[0, 0], 3.0)
        self.assertAlmostEqual(result.state.particles.bond_count_expected[0], 5.0)
        self.assertAlmostEqual(
            result.state.particles.bond_total_tangential_extension_um[0], 1.5
        )
        self.assertEqual(result.diagnostics.contact_time_refinement_count, 1)

    def test_zero_time_lifecycle_termination_is_a_legal_state_change(self) -> None:
        context = _context(1)
        ids = np.asarray([0], dtype=np.int64)
        start = transport._BatchAdvanceState(
            particles=BatchLocalState(
                position_grid=np.asarray([[2.0, 2.0]]),
                rotation_angle_rad=np.asarray([0.0]),
                bond_count_expected=np.asarray([0.0]),
                bond_total_tangential_extension_um=np.asarray([0.0]),
            ),
            alive=np.asarray([True]),
            termination_reason=np.asarray([0], dtype=np.uint8),
            exit_time_s=np.asarray([np.nan]),
        )
        calls = 0

        def attempt(interval, start_state, *args, **kwargs):
            nonlocal calls
            del args, kwargs
            calls += 1
            terminated = transport._accepted_batch_state(
                start_state,
                start_state.particles.position_grid.copy(),
                start_state.particles.rotation_angle_rad.copy(),
                start_state.particles.bond_count_expected.copy(),
                start_state.particles.bond_total_tangential_extension_um.copy(),
                alive=np.asarray([False]),
                termination_reason=np.asarray(
                    [transport._TERMINATION_OUTLET], dtype=np.uint8
                ),
                exit_time_s=np.asarray([float(interval.start_time_s)]),
            )
            return transport._BatchAttemptResult(
                terminated,
                transport._Diagnostics(),
                float(interval.start_time_s),
            )

        with mock.patch.object(transport, "_attempt_batch_step", side_effect=attempt):
            result = transport._advance_batch_interval(
                PhysicalTimeInterval(0.0, 1.0),
                start,
                ids,
                context,
                _particle_config(maximum_refinements=4),
                self.euler,
                object(),
            )

        # A zero-duration prefix is legal only because the particle really
        # leaves the active population.  No persistent wall-contact state is
        # created or compared in Revised v15.
        self.assertEqual(calls, 1)
        self.assertFalse(result.state.alive[0])
        self.assertEqual(
            result.state.termination_reason[0], transport._TERMINATION_OUTLET
        )
        self.assertEqual(result.state.exit_time_s[0], 0.0)
        self.assertEqual(result.diagnostics.contact_time_refinement_count, 0)

    def test_zero_time_event_without_lifecycle_change_still_raises(self) -> None:
        context = _context(1)
        ids = np.asarray([0], dtype=np.int64)
        start = transport._BatchAdvanceState(
            particles=BatchLocalState(
                position_grid=np.asarray([[2.0, 1.0]]),
                rotation_angle_rad=np.asarray([0.0]),
                bond_count_expected=np.asarray([0.0]),
                bond_total_tangential_extension_um=np.asarray([0.0]),
            ),
            alive=np.asarray([True]),
            termination_reason=np.asarray([0], dtype=np.uint8),
            exit_time_s=np.asarray([np.nan]),
        )

        with mock.patch.object(
            transport,
            "_attempt_batch_step",
            return_value=transport._BatchAttemptResult(
                start,
                transport._Diagnostics(),
                0.0,
            ),
        ):
            with self.assertRaisesRegex(RuntimeError, "changing lifecycle state"):
                transport._advance_batch_interval(
                    PhysicalTimeInterval(0.0, 1.0),
                    start,
                    ids,
                    context,
                    _particle_config(),
                    self.euler,
                    object(),
                )

    def test_exit_time_uses_the_exact_event_fraction_not_the_step_end(self) -> None:
        context = _context(1, outlet_x=11.0, maximum_x=11.0)
        state = _state(np.asarray([[10.0, 2.0]]))

        transport._advance_segment(
            7.0,
            9.0,
            state,
            context,
            _particle_config(),
            self.euler,
            transport._Diagnostics(),
            _rhs(np.asarray([[2.0, 0.0]])),
        )

        # The outlet is at x=11.  Starting from x=10 at 2 um/s reaches it
        # after 0.5 s, even though the requested parent step ends at t=9 s.
        self.assertFalse(state.active[0])
        self.assertEqual(state.termination_reason[0], transport._TERMINATION_OUTLET)
        self.assertAlmostEqual(state.exit_time_s[0], 7.5, places=12)
        np.testing.assert_allclose(state.position_grid[0], [11.0, 2.0])

    def test_internal_commitment_crossing_splits_and_commits_one_owner_change(self) -> None:
        context = _context(1)
        context.topological_ownership = _single_commitment_catalog()
        context.molecular_target_field = SimpleNamespace(enabled=True)
        context.molecular_capture_distance_um = 0.2
        state = _state(np.asarray([[2.0, 2.0]]))
        state.vessel_id = np.asarray([1], dtype=np.int32)
        diagnostics = transport._Diagnostics()
        accepted_path_batches = []

        transport._advance_segment(
            0.0,
            1.0,
            state,
            context,
            _particle_config(),
            self.euler,
            diagnostics,
            _rhs(np.asarray([[2.0, 0.0]])),
            accepted_target_path_batches=accepted_path_batches,
        )

        np.testing.assert_allclose(state.position_grid[0], [4.0, 2.0])
        self.assertEqual(int(state.vessel_id[0]), 2)
        self.assertEqual(diagnostics.topological_transition_count, 1)
        self.assertEqual(diagnostics.topological_event_bubble_id, [0])
        self.assertEqual(diagnostics.topological_event_from_vessel_id, [1])
        self.assertEqual(diagnostics.topological_event_to_vessel_id, [2])
        self.assertAlmostEqual(diagnostics.topological_event_time_s[0], 0.5)
        self.assertEqual(len(accepted_path_batches), 1)
        batch_ids, path_positions, path_active = accepted_path_batches[0]
        np.testing.assert_array_equal(batch_ids, [0])
        accepted_vertices = np.vstack(path_positions)
        self.assertTrue(
            np.any(np.all(np.isclose(accepted_vertices, [3.0, 2.0]), axis=1))
        )
        np.testing.assert_allclose(accepted_vertices[-1], [4.0, 2.0])
        self.assertTrue(all(bool(mask[0]) for mask in path_active))

    def test_heun_recomputes_the_remainder_with_the_committed_child_owner(self) -> None:
        context = _context(1)
        context.topological_ownership = _single_commitment_catalog()
        state = _state(np.asarray([[2.0, 2.0]]))
        state.vessel_id = np.asarray([1], dtype=np.int32)
        diagnostics = transport._Diagnostics()

        transport._advance_segment(
            0.0,
            1.0,
            state,
            context,
            _particle_config(),
            SimpleNamespace(time_integrator="heun"),
            diagnostics,
            _rhs(np.asarray([[2.0, 0.0]])),
        )

        np.testing.assert_allclose(state.position_grid[0], [4.0, 2.0])
        self.assertEqual(int(state.vessel_id[0]), 2)
        self.assertEqual(diagnostics.topological_transition_count, 1)
        self.assertAlmostEqual(diagnostics.topological_event_time_s[0], 0.5)


def _single_commitment_catalog() -> TopologicalCommitmentCatalog:
    return TopologicalCommitmentCatalog(
        schema="revised_v20_continuous_commitment_sections_v1",
        root_vessel_id=np.asarray([1], dtype=np.int32),
        terminal_vessel_id=np.asarray([2], dtype=np.int32),
        parent_vessel_id=np.asarray([1], dtype=np.int32),
        child_vessel_id=np.asarray([2], dtype=np.int32),
        point_xz_um=np.asarray([[3.0, 2.5]], dtype=np.float64),
        downstream_normal_xz=np.asarray([[1.0, 0.0]], dtype=np.float64),
        tangent_xz=np.asarray([[0.0, 1.0]], dtype=np.float64),
        half_width_um=np.asarray([2.5], dtype=np.float64),
        transition_end_distance_um=np.asarray([1.0], dtype=np.float64),
        commitment_distance_um=np.asarray([2.0], dtype=np.float64),
        child_section_by_vessel_id=np.asarray([-1, -1, 0], dtype=np.int32),
        child_section_offsets_by_parent_id=np.asarray([0, 0, 1, 1], dtype=np.int64),
        child_section_indices_by_parent_id=np.asarray([0], dtype=np.int32),
        geometry_hash_sha256="test",
    )


if __name__ == "__main__":
    unittest.main()
