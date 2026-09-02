from __future__ import annotations

import unittest

from ulm_microbubble_traj_gen_2D.utils.particles.particle_time_stepping import (
    build_particle_step_resolution,
    build_particle_time_step_plan,
)


class ParticleTimeSteppingTests(unittest.TestCase):
    def test_fixed_substeps_land_exactly_on_saved_frame_times(self) -> None:
        plan = build_particle_time_step_plan(1000, 0.001, 40, "euler")

        self.assertEqual(plan.output_intervals, 1000)
        self.assertEqual(plan.stored_frames, 1001)
        self.assertEqual(plan.integration_substeps, 40)
        self.assertAlmostEqual(plan.internal_dt_s, 2.5e-5)
        self.assertEqual(plan.total_internal_steps, 40_000)
        self.assertEqual(plan.expected_rhs_evaluations, 40_001)
        self.assertAlmostEqual(
            plan.total_internal_steps * plan.internal_dt_s,
            plan.output_intervals * plan.output_dt_s,
        )

    def test_heun_counts_two_rhs_stages_per_internal_step(self) -> None:
        plan = build_particle_time_step_plan(3, 0.01, 4, "heun")

        self.assertEqual(plan.total_internal_steps, 12)
        self.assertEqual(plan.expected_rhs_evaluations, 25)

    def test_invalid_substep_counts_are_rejected(self) -> None:
        for substeps in (0, -1):
            with self.subTest(substeps=substeps):
                with self.assertRaisesRegex(ValueError, "at least one"):
                    build_particle_time_step_plan(3, 0.01, substeps, "euler")
        with self.assertRaisesRegex(ValueError, "integer"):
            build_particle_time_step_plan(3, 0.01, 1.5, "euler")

    def test_resolution_separates_transport_and_collision_diagnostics(self) -> None:
        resolution = build_particle_step_resolution(
            maximum_internal_step_displacement_um=0.1,
            grid_spacing_um=1.5,
            minimum_radius_um=0.75,
            collision_layer_um=0.05,
            collisions_enabled=True,
            internal_dt_s=2.5e-5,
            collision_relaxation_time_s=0.01,
        )

        self.assertAlmostEqual(resolution.grid_displacement_ratio, 1.0 / 15.0)
        self.assertAlmostEqual(resolution.radius_displacement_ratio, 2.0 / 15.0)
        self.assertAlmostEqual(resolution.field_size_displacement_ratio, 2.0 / 15.0)
        self.assertAlmostEqual(resolution.conservative_collision_layer_ratio, 2.0)
        self.assertAlmostEqual(resolution.collision_dt_over_relaxation_time, 0.0025)


if __name__ == "__main__":
    unittest.main()
