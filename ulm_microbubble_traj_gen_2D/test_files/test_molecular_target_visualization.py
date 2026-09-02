from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.visualization.results.molecular_target_overlay import (
    load_molecular_target_overlay,
    target_display_masks,
)


class MolecularTargetVisualizationTests(unittest.TestCase):
    def test_missing_target_file_returns_invisible_empty_overlay(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            overlay = load_molecular_target_overlay(
                Path(directory),
                np.asarray([0.0, 1.0]),
                np.asarray([0.0, 1.0, 2.0]),
            )

        self.assertFalse(overlay.visible)
        self.assertEqual(overlay.target_wall_mask.shape, (2, 3))
        self.assertEqual(overlay.target_wall_site_count, 0)

    def test_target_file_loads_mask_density_and_visibility(self) -> None:
        x_um = np.asarray([0.0, 1.0, 2.0])
        z_um = np.asarray([10.0, 11.0])
        mask = np.zeros((3, 2), dtype=bool)
        mask[1, 1] = True
        density_field = np.zeros((3, 2), dtype=float)
        density_field[mask] = 2.7e14

        with tempfile.TemporaryDirectory() as directory:
            target_path = Path(directory) / "molecular_target_field.npz"
            np.savez_compressed(
                target_path,
                enabled=np.asarray([True]),
                x_coordinates_um=x_um,
                z_coordinates_um=z_um,
                target_wall_mask=mask,
                target_density_field_molecules_per_m2=density_field,
                target_density_molecules_per_m2=np.asarray([2.7e14]),
            )
            overlay = load_molecular_target_overlay(Path(directory), x_um, z_um)

        self.assertTrue(overlay.visible)
        self.assertEqual(overlay.target_wall_site_count, 1)
        self.assertEqual(overlay.target_density_molecules_per_m2, 2.7e14)
        np.testing.assert_array_equal(overlay.target_wall_mask, mask)
        np.testing.assert_array_equal(
            overlay.target_density_field_molecules_per_m2,
            density_field,
        )

    def test_mismatched_target_coordinates_are_rejected(self) -> None:
        x_um = np.asarray([0.0, 1.0])
        z_um = np.asarray([0.0, 1.0])
        with tempfile.TemporaryDirectory() as directory:
            np.savez_compressed(
                Path(directory) / "molecular_target_field.npz",
                enabled=np.asarray([True]),
                x_coordinates_um=np.asarray([0.0, 2.0]),
                z_coordinates_um=z_um,
                target_wall_mask=np.zeros((2, 2), dtype=bool),
                target_density_field_molecules_per_m2=np.zeros((2, 2), dtype=float),
                target_density_molecules_per_m2=np.asarray([0.0]),
            )
            with self.assertRaisesRegex(ValueError, "coordinates do not match"):
                load_molecular_target_overlay(Path(directory), x_um, z_um)

    def test_display_halo_never_changes_the_exact_target_mask(self) -> None:
        mask = np.zeros((7, 7), dtype=bool)
        mask[3, 3] = True

        exact, halo = target_display_masks(mask, halo_cells=2)

        np.testing.assert_array_equal(exact, mask)
        self.assertFalse(np.any(exact & halo))
        self.assertGreater(int(np.count_nonzero(halo)), 0)


if __name__ == "__main__":
    unittest.main()
