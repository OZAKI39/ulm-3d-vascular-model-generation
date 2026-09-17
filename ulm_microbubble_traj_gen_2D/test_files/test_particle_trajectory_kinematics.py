from __future__ import annotations

import unittest

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.particles.particle_trajectory_kinematics import (
    realized_center_velocities_um_s,
)


class ParticleTrajectoryKinematicsTests(unittest.TestCase):
    def test_realized_velocity_uses_only_adjacent_records_of_the_same_id(self) -> None:
        offsets = np.asarray([0, 2, 4, 5], dtype=np.int64)
        bubble_id = np.asarray([0, 1, 0, 2, 2], dtype=np.int64)
        positions = np.asarray(
            [
                [0.0, 0.0, 0.0],
                [10.0, 0.0, 0.0],
                [2.0, 0.0, 0.0],
                [100.0, 0.0, 0.0],
                [103.0, 0.0, 0.0],
            ]
        )

        realized = realized_center_velocities_um_s(
            offsets,
            bubble_id,
            positions,
            0.5,
        )

        np.testing.assert_allclose(realized[0], [4.0, 0.0, 0.0])
        np.testing.assert_allclose(realized[2], [4.0, 0.0, 0.0])
        np.testing.assert_allclose(realized[1], [0.0, 0.0, 0.0])
        np.testing.assert_allclose(realized[3], [6.0, 0.0, 0.0])
        np.testing.assert_allclose(realized[4], [6.0, 0.0, 0.0])


if __name__ == "__main__":
    unittest.main()
