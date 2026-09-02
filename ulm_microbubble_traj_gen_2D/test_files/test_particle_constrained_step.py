from __future__ import annotations

import math
import unittest

from ulm_microbubble_traj_gen_2D.utils.particles.particle_constrained_step import (
    PhysicalTimeInterval,
    PhysicalTimeRefinementError,
    split_physical_time_interval,
    validate_physical_time_partition,
)


class PhysicalTimeRefinementTests(unittest.TestCase):
    def test_bisection_exactly_tiles_parent_time_interval(self) -> None:
        parent = PhysicalTimeInterval(0.125, 0.126, refinement_depth=2)
        first, second = split_physical_time_interval(parent, 8)

        self.assertEqual(first.start_time_s, parent.start_time_s)
        self.assertEqual(first.end_time_s, second.start_time_s)
        self.assertEqual(second.end_time_s, parent.end_time_s)
        self.assertEqual(first.refinement_depth, 3)
        self.assertEqual(second.refinement_depth, 3)
        validate_physical_time_partition(parent, first, second)

    def test_refinement_limit_fails_instead_of_returning_a_hold_step(self) -> None:
        interval = PhysicalTimeInterval(1.0, 1.1, refinement_depth=4)
        with self.assertRaisesRegex(
            PhysicalTimeRefinementError,
            "not accepted.*not silently held",
        ) as captured:
            split_physical_time_interval(interval, 4)

        error = captured.exception
        self.assertEqual(error.reason, "depth_limit")
        self.assertEqual(error.refinement_depth, 4)
        self.assertEqual(error.maximum_refinement_depth, 4)
        self.assertAlmostEqual(error.duration_s, interval.duration_s)
        self.assertEqual(
            error.local_time_ulp_s,
            max(math.ulp(interval.start_time_s), math.ulp(interval.end_time_s)),
        )

    def test_unrepresentably_small_interval_fails_fast(self) -> None:
        start = 1.0
        end = math.nextafter(start, math.inf)
        interval = PhysicalTimeInterval(start, end, refinement_depth=7)
        with self.assertRaisesRegex(
            PhysicalTimeRefinementError, "too small"
        ) as captured:
            split_physical_time_interval(interval, 10)

        error = captured.exception
        self.assertEqual(error.reason, "unrepresentable_midpoint")
        self.assertEqual(error.refinement_depth, 7)
        self.assertEqual(error.maximum_refinement_depth, 10)
        self.assertEqual(error.duration_s, math.ulp(start))
        self.assertEqual(
            error.local_time_ulp_s,
            max(math.ulp(start), math.ulp(end)),
        )

    def test_partition_validator_rejects_a_gap(self) -> None:
        parent = PhysicalTimeInterval(0.0, 1.0)
        first = PhysicalTimeInterval(0.0, 0.4, 1)
        second = PhysicalTimeInterval(0.5, 1.0, 1)
        with self.assertRaisesRegex(ValueError, "gap, overlap"):
            validate_physical_time_partition(parent, first, second)

    def test_invalid_intervals_and_depth_limits_are_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "positive duration"):
            PhysicalTimeInterval(2.0, 2.0)
        with self.assertRaisesRegex(ValueError, "non-negative integer"):
            PhysicalTimeInterval(0.0, 1.0, -1)
        with self.assertRaisesRegex(ValueError, "non-negative integer"):
            split_physical_time_interval(PhysicalTimeInterval(0.0, 1.0), -1)


if __name__ == "__main__":
    unittest.main()
