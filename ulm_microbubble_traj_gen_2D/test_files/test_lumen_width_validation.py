"""Validation tests for local lumen diameter mask quality metrics."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from ulm_microbubble_traj_gen_2D.utils.geometry.lumen_width import classify_narrow_lumen_cells


class LocalLumenDiameterTests(unittest.TestCase):
    def test_graph_diameter_marks_underresolved_segment(self) -> None:
        lumen = np.zeros((48, 48), dtype=bool)
        lumen[8:40, 12:36] = True
        radius_um = np.full(lumen.shape, np.nan, dtype=float)
        radius_um[lumen] = 3.0

        result = classify_narrow_lumen_cells(lumen, spacing_um=1.0, radius_um=radius_um, min_diameter_px=8.0)

        self.assertTrue(np.all(result.narrow_mask[lumen]))
        self.assertAlmostEqual(float(np.nanmedian(result.graph_diameter_px[lumen])), 6.0)
        self.assertLess(float(np.nanmax(result.effective_diameter_px[lumen])), 8.0)

    def test_mask_diameter_marks_local_channel_width_not_wall_distance(self) -> None:
        lumen = np.zeros((90, 60), dtype=bool)
        lumen[8:58, 12:28] = True
        lumen[58:80, 18:21] = True
        radius_um = np.full(lumen.shape, np.nan, dtype=float)
        radius_um[lumen] = 20.0

        result = classify_narrow_lumen_cells(lumen, spacing_um=1.0, radius_um=radius_um, min_diameter_px=8.0)

        central_wide = lumen[18:48, 15:25]
        near_wall_wide = lumen[18:48, 12:14]
        narrow_channel = lumen[62:76, 18:21]
        self.assertLess(np.count_nonzero(result.narrow_mask[18:48, 15:25]), 0.05 * np.count_nonzero(central_wide))
        self.assertLess(np.count_nonzero(result.narrow_mask[18:48, 12:14]), 0.10 * np.count_nonzero(near_wall_wide))
        self.assertGreater(np.count_nonzero(result.narrow_mask[62:76, 18:21]), 0.80 * np.count_nonzero(narrow_channel))
        self.assertGreater(float(np.nanmedian(result.mask_diameter_px[18:48, 15:25])), 8.0)
        self.assertLess(float(np.nanmedian(result.mask_diameter_px[62:76, 18:21])), 8.0)

    def test_junction_cells_use_mask_diameter_for_narrow_classification(self) -> None:
        lumen = np.zeros((80, 80), dtype=bool)
        lumen[20:60, 32:48] = True
        lumen[32:48, 20:60] = True
        radius_um = np.full(lumen.shape, np.nan, dtype=float)
        radius_um[lumen] = 3.0
        junction = np.zeros(lumen.shape, dtype=bool)
        junction[26:54, 26:54] = lumen[26:54, 26:54]

        without_junction = classify_narrow_lumen_cells(lumen, spacing_um=1.0, radius_um=radius_um, min_diameter_px=8.0)
        with_junction = classify_narrow_lumen_cells(
            lumen,
            spacing_um=1.0,
            radius_um=radius_um,
            min_diameter_px=8.0,
            junction_mask=junction,
        )

        self.assertGreater(np.count_nonzero(without_junction.narrow_mask & junction), 0)
        self.assertEqual(int(np.count_nonzero(with_junction.narrow_mask & junction)), 0)

    def test_open_boundary_cells_use_graph_diameter_not_cap_distance(self) -> None:
        lumen = np.zeros((80, 40), dtype=bool)
        lumen[10:70, 16:24] = True
        radius_um = np.full(lumen.shape, np.nan, dtype=float)
        radius_um[lumen] = 4.1
        open_boundary = np.zeros(lumen.shape, dtype=bool)
        open_boundary[10:13, 16:24] = True
        open_boundary[67:70, 16:24] = True

        without_open_boundary = classify_narrow_lumen_cells(lumen, spacing_um=1.0, radius_um=radius_um, min_diameter_px=8.0)
        with_open_boundary = classify_narrow_lumen_cells(
            lumen,
            spacing_um=1.0,
            radius_um=radius_um,
            min_diameter_px=8.0,
            open_boundary_mask=open_boundary,
        )

        self.assertGreater(np.count_nonzero(without_open_boundary.narrow_mask & open_boundary), 0)
        self.assertEqual(int(np.count_nonzero(with_open_boundary.narrow_mask & open_boundary)), 0)


if __name__ == "__main__":
    unittest.main()
