from __future__ import annotations

import unittest

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.core.types import HybridVelocityField
from ulm_microbubble_traj_gen_2D.utils.flow.hybrid_velocity import (
    blend_velocity_and_gradient,
    build_finite_element_velocity_field,
    hybrid_region_map,
    interpolation_nodes,
    sample_finite_element_velocity,
)


class HybridVelocityTests(unittest.TestCase):
    def test_quadratic_triangle_velocity_and_gradient_are_reconstructed(self):
        vertices = np.asarray([[[0.0, 0.0], [2.0, 0.0], [0.0, 3.0]]])
        nodes = interpolation_nodes(2)
        world = (
            vertices[:, None, 0]
            + nodes[None, :, 0, None]
            * (vertices[:, None, 1] - vertices[:, None, 0])
            + nodes[None, :, 1, None]
            * (vertices[:, None, 2] - vertices[:, None, 0])
        )
        x = world[..., 0]
        z = world[..., 1]
        values = np.stack(
            (
                1.0 + 2.0 * x + 3.0 * z + 0.5 * x * x + x * z,
                -2.0 + x - 4.0 * z + 0.25 * z * z,
            ),
            axis=-1,
        )
        field = build_finite_element_velocity_field(
            vertices,
            values,
            2,
            preferred_bin_size_um=1.0,
        )

        points = np.asarray([[0.4, 0.6], [1.0, 0.3]])
        velocity, gradient, cells = sample_finite_element_velocity(
            field,
            points,
            np.ones(points.shape[0], dtype=bool),
            use_numba=False,
        )
        px = points[:, 0]
        pz = points[:, 1]
        expected_velocity = np.stack(
            (
                1.0 + 2.0 * px + 3.0 * pz + 0.5 * px * px + px * pz,
                -2.0 + px - 4.0 * pz + 0.25 * pz * pz,
            ),
            axis=-1,
        )
        expected_gradient = np.empty((points.shape[0], 2, 2))
        expected_gradient[:, 0, 0] = 2.0 + px + pz
        expected_gradient[:, 0, 1] = 3.0 + px
        expected_gradient[:, 1, 0] = 1.0
        expected_gradient[:, 1, 1] = -4.0 + 0.5 * pz
        np.testing.assert_allclose(velocity, expected_velocity, atol=1.0e-12)
        np.testing.assert_allclose(gradient, expected_gradient, atol=1.0e-12)
        np.testing.assert_array_equal(cells, np.zeros(2, dtype=np.int32))

    def test_region_map_and_blend_are_continuous(self):
        vertices = np.asarray([[[0.0, 0.0], [2.0, 0.0], [0.0, 2.0]]])
        values = np.zeros((1, interpolation_nodes(1).shape[0], 2))
        finite_element = build_finite_element_velocity_field(
            vertices,
            values,
            1,
            preferred_bin_size_um=1.0,
        )
        hybrid = HybridVelocityField(
            finite_element=finite_element,
            finite_element_distance_um=1.0,
            regular_grid_distance_um=3.0,
        )
        distance = np.asarray([[0.5, 1.0, 2.0, 3.0, 4.0]])
        lumen = np.ones(distance.shape, dtype=bool)
        region, weight = hybrid_region_map(distance, lumen, hybrid)
        np.testing.assert_array_equal(region, [[1, 1, 2, 3, 3]])
        np.testing.assert_allclose(weight, [[1.0, 1.0, 0.5, 0.0, 0.0]])

        grid_velocity = np.zeros((3, 2))
        finite_element_velocity = np.ones((3, 2))
        zero_gradient = np.zeros((3, 2, 2))
        normals = np.tile([1.0, 0.0], (3, 1))
        velocity, gradient, sampled_weight = blend_velocity_and_gradient(
            grid_velocity,
            zero_gradient,
            finite_element_velocity,
            zero_gradient,
            np.asarray([1.0, 2.0, 3.0]),
            normals,
            hybrid,
        )
        np.testing.assert_allclose(sampled_weight, [1.0, 0.5, 0.0])
        np.testing.assert_allclose(velocity[:, 0], [1.0, 0.5, 0.0])
        np.testing.assert_allclose(gradient[[0, 2]], 0.0)
        np.testing.assert_allclose(gradient[1, :, 0], -0.75)


if __name__ == "__main__":
    unittest.main()
