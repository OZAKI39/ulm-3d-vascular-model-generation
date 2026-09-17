from __future__ import annotations

import unittest

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.particles.particle_collisions import (
    compute_collision_forces,
    resolve_collision_search_strategy,
)


class ParticleCollisionTests(unittest.TestCase):
    def test_cell_list_matches_all_pairs_for_variable_radii(self) -> None:
        rng = np.random.default_rng(20260717)
        positions = rng.uniform(-5.0, 15.0, size=(96, 2))
        radii = rng.uniform(0.75, 1.25, size=96)
        mobility = np.asarray(
            [
                [[scale, 0.02 * scale], [0.02 * scale, 1.2 * scale]]
                for scale in rng.uniform(0.5, 1.5, size=96)
            ],
            dtype=np.float64,
        )
        active = rng.random(96) > 0.1
        ids = rng.permutation(np.arange(96, dtype=np.int64)) + 1000
        common = dict(
            positions_xz_um=positions,
            radii_um=radii,
            translational_mobility_xz_um_pn_s=mobility,
            active=active,
            bubble_ids=ids,
            collision_layer_um=0.15,
            relaxation_time_s=0.02,
            use_numba=True,
        )

        reference = compute_collision_forces(strategy="all_pairs", **common)
        sparse = compute_collision_forces(strategy="cell_list", **common)

        np.testing.assert_allclose(
            sparse.force_xz_pn,
            reference.force_xz_pn,
            rtol=2.0e-13,
            atol=2.0e-13,
        )
        np.testing.assert_array_equal(sparse.neighbor_count, reference.neighbor_count)
        self.assertEqual(sparse.interacting_pair_count, reference.interacting_pair_count)
        self.assertAlmostEqual(
            sparse.maximum_physical_overlap_um,
            reference.maximum_physical_overlap_um,
        )
        self.assertAlmostEqual(
            sparse.maximum_collision_compression_um,
            reference.maximum_collision_compression_um,
        )
        self.assertEqual(sparse.search_strategy, "cell_list")

    def test_auto_selects_sparse_search_for_large_population(self) -> None:
        self.assertEqual(resolve_collision_search_strategy("auto", 128), "all_pairs")
        self.assertEqual(resolve_collision_search_strategy("auto", 1024), "cell_list")

    def test_interacting_pair_receives_equal_and_opposite_force(self) -> None:
        result = compute_collision_forces(
            positions_xz_um=np.asarray([[0.0, 0.0], [1.5, 0.0]]),
            radii_um=np.asarray([1.0, 1.0]),
            translational_mobility_xz_um_pn_s=_isotropic_mobility([1.0, 1.0]),
            active=np.asarray([True, True]),
            bubble_ids=np.asarray([4, 9]),
            collision_layer_um=0.0,
            relaxation_time_s=0.5,
            use_numba=False,
        )

        np.testing.assert_allclose(result.force_xz_pn[0], -result.force_xz_pn[1], rtol=0.0, atol=0.0)
        np.testing.assert_allclose(result.force_xz_pn.sum(axis=0), 0.0, rtol=0.0, atol=0.0)
        np.testing.assert_array_equal(result.neighbor_count, np.asarray([1, 1]))
        self.assertEqual(result.interacting_pair_count, 1)
        self.assertAlmostEqual(result.maximum_physical_overlap_um, 0.5)
        self.assertAlmostEqual(result.maximum_collision_compression_um, 0.5)

    def test_coincident_centres_use_a_repeatable_permanent_id_direction(self) -> None:
        ids = np.asarray([101, 7])
        kwargs = dict(
            positions_xz_um=np.zeros((2, 2), dtype=float),
            radii_um=np.asarray([1.0, 1.0]),
            translational_mobility_xz_um_pn_s=_isotropic_mobility([1.0, 1.0]),
            active=np.asarray([True, True]),
            bubble_ids=ids,
            collision_layer_um=0.0,
            relaxation_time_s=1.0,
            use_numba=False,
        )

        first = compute_collision_forces(**kwargs)
        second = compute_collision_forces(**kwargs)

        np.testing.assert_array_equal(first.force_xz_pn, second.force_xz_pn)
        self.assertTrue(np.all(np.isfinite(first.force_xz_pn)))
        self.assertGreater(float(np.linalg.norm(first.force_xz_pn[0])), 0.0)
        np.testing.assert_allclose(first.force_xz_pn[0], -first.force_xz_pn[1], rtol=0.0, atol=0.0)

        low_id = int(np.min(ids))
        high_id = int(np.max(ids))
        code = (low_id * 73_856_093 + high_id * 19_349_663) % 104_729
        expected_low_id_direction = np.asarray(
            [
                np.cos(2.0 * np.pi * (code + 0.5) / 104_729.0),
                np.sin(2.0 * np.pi * (code + 0.5) / 104_729.0),
            ]
        )
        low_id_force = first.force_xz_pn[int(np.flatnonzero(ids == low_id)[0])]
        np.testing.assert_allclose(
            low_id_force / np.linalg.norm(low_id_force),
            expected_low_id_direction,
            rtol=0.0,
            atol=1.0e-14,
        )

    def test_storage_order_permutation_is_physically_equivalent_by_id(self) -> None:
        positions = np.asarray([[0.0, 0.0], [1.1, 0.2], [0.45, 1.0]])
        radii = np.asarray([0.9, 0.8, 0.75])
        mobilities = np.asarray(
            [
                [[1.0, 0.1], [0.1, 1.4]],
                [[0.8, -0.05], [-0.05, 1.2]],
                [[1.5, 0.2], [0.2, 0.9]],
            ]
        )
        ids = np.asarray([21, 5, 12])
        active = np.ones(3, dtype=bool)
        original = compute_collision_forces(
            positions,
            radii,
            mobilities,
            active,
            ids,
            collision_layer_um=0.15,
            relaxation_time_s=0.4,
            use_numba=False,
        )
        permutation = np.asarray([2, 0, 1])
        permuted = compute_collision_forces(
            positions[permutation],
            radii[permutation],
            mobilities[permutation],
            active[permutation],
            ids[permutation],
            collision_layer_um=0.15,
            relaxation_time_s=0.4,
            use_numba=False,
        )

        original_force_by_id = {int(bubble_id): original.force_xz_pn[index] for index, bubble_id in enumerate(ids)}
        permuted_force_by_id = {
            int(bubble_id): permuted.force_xz_pn[index] for index, bubble_id in enumerate(ids[permutation])
        }
        original_neighbors_by_id = {
            int(bubble_id): int(original.neighbor_count[index]) for index, bubble_id in enumerate(ids)
        }
        permuted_neighbors_by_id = {
            int(bubble_id): int(permuted.neighbor_count[index]) for index, bubble_id in enumerate(ids[permutation])
        }
        self.assertEqual(original_force_by_id.keys(), permuted_force_by_id.keys())
        for bubble_id in original_force_by_id:
            np.testing.assert_allclose(
                original_force_by_id[bubble_id],
                permuted_force_by_id[bubble_id],
                rtol=1.0e-14,
                atol=1.0e-14,
            )
        self.assertEqual(original_neighbors_by_id, permuted_neighbors_by_id)
        self.assertEqual(original.interacting_pair_count, permuted.interacting_pair_count)
        self.assertAlmostEqual(original.maximum_physical_overlap_um, permuted.maximum_physical_overlap_um)
        self.assertAlmostEqual(original.maximum_collision_compression_um, permuted.maximum_collision_compression_um)

    def test_particles_beyond_collision_layer_have_zero_force(self) -> None:
        result = compute_collision_forces(
            positions_xz_um=np.asarray([[0.0, 0.0], [2.26, 0.0]]),
            radii_um=np.asarray([1.0, 1.0]),
            translational_mobility_xz_um_pn_s=_isotropic_mobility([1.0, 1.0]),
            active=np.asarray([True, True]),
            bubble_ids=np.asarray([0, 1]),
            collision_layer_um=0.25,
            relaxation_time_s=0.5,
            use_numba=False,
        )

        np.testing.assert_array_equal(result.force_xz_pn, np.zeros((2, 2)))
        np.testing.assert_array_equal(result.neighbor_count, np.zeros(2, dtype=np.int32))
        self.assertEqual(result.interacting_pair_count, 0)
        self.assertEqual(result.maximum_physical_overlap_um, 0.0)
        self.assertEqual(result.maximum_collision_compression_um, 0.0)

    def test_force_magnitude_uses_both_particles_different_mobilities(self) -> None:
        mobilities = np.asarray(
            [
                [[1.0, 0.0], [0.0, 2.0]],
                [[3.0, 0.0], [0.0, 4.0]],
            ]
        )
        result = compute_collision_forces(
            positions_xz_um=np.asarray([[0.0, 0.0], [1.5, 0.0]]),
            radii_um=np.asarray([1.0, 1.0]),
            translational_mobility_xz_um_pn_s=mobilities,
            active=np.asarray([True, True]),
            bubble_ids=np.asarray([0, 1]),
            collision_layer_um=0.0,
            relaxation_time_s=0.25,
            use_numba=False,
        )

        # Compression is 0.5 um and relative mobility along X is 1 + 3 = 4 um/(pN s).
        # Therefore |F| = 0.5 / (4 * 0.25) = 0.5 pN.
        np.testing.assert_allclose(result.force_xz_pn, [[-0.5, 0.0], [0.5, 0.0]], rtol=0.0, atol=0.0)

    def test_nonpositive_or_nonfinite_relative_mobility_is_rejected(self) -> None:
        invalid_mobilities = {
            "zero": np.zeros((2, 2, 2), dtype=float),
            "negative": _isotropic_mobility([-1.0, -1.0]),
            "nan": np.full((2, 2, 2), np.nan),
            "infinite": np.full((2, 2, 2), np.inf),
        }
        for label, mobility in invalid_mobilities.items():
            with self.subTest(label=label):
                with np.errstate(all="ignore"):
                    with self.assertRaisesRegex(ValueError, "finite and strictly positive"):
                        compute_collision_forces(
                            positions_xz_um=np.asarray([[0.0, 0.0], [1.5, 0.0]]),
                            radii_um=np.asarray([1.0, 1.0]),
                            translational_mobility_xz_um_pn_s=mobility,
                            active=np.asarray([True, True]),
                            bubble_ids=np.asarray([0, 1]),
                            collision_layer_um=0.0,
                            relaxation_time_s=0.5,
                            use_numba=False,
                        )


def _isotropic_mobility(scales: list[float]) -> np.ndarray:
    return np.asarray([np.eye(2, dtype=float) * scale for scale in scales])


if __name__ == "__main__":
    unittest.main()
