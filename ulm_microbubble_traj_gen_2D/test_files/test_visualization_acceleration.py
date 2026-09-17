from __future__ import annotations

from pathlib import Path
import tempfile
import unittest

import numpy as np
from PIL import Image

from ulm_microbubble_traj_gen_2D.utils.visualization.results.plotting import (
    _PillowGifWriter,
    _wall_shear_color_limits,
)
from ulm_microbubble_traj_gen_2D.utils.visualization.results.result_loader import (
    _assign_display_lanes,
    _assign_display_lanes_impl,
)


class VisualizationAccelerationTests(unittest.TestCase):
    def test_wall_shear_scale_is_shared_and_robust_to_one_extreme_cell(self) -> None:
        class _Data:
            lumen_mask = np.ones((20, 20), dtype=bool)
            wall_shear_grid_pa = np.ones((20, 20), dtype=float)

        data = _Data()
        data.wall_shear_grid_pa[-1, -1] = 1.0e6
        lower, upper = _wall_shear_color_limits(data)

        self.assertEqual(lower, 0.0)
        self.assertGreaterEqual(upper, 1.0)
        self.assertLess(upper, 1.0e6)

    def test_accelerated_lane_assignment_preserves_recycling_order(self) -> None:
        offsets = np.asarray([0, 2, 4, 5, 7], dtype=np.int64)
        compact_ids = np.asarray([0, 1, 0, 1, 1, 1, 2], dtype=np.int64)
        expected = np.asarray([0, 1, 0, 1, 1, 1, 0], dtype=np.int64)

        np.testing.assert_array_equal(
            _assign_display_lanes_impl(offsets, compact_ids, 2),
            expected,
        )
        np.testing.assert_array_equal(
            _assign_display_lanes(offsets, compact_ids, 2),
            expected,
        )

    def test_streaming_gif_preserves_frame_count_and_pixels(self) -> None:
        first = np.full((24, 32, 3), 12, dtype=np.uint8)
        second = first.copy()
        first[4:10, 3:9] = np.asarray([240, 40, 20], dtype=np.uint8)
        second[12:18, 20:26] = np.asarray([240, 40, 20], dtype=np.uint8)

        with tempfile.TemporaryDirectory() as directory:
            output_path = Path(directory) / "streamed.gif"
            with _PillowGifWriter(output_path, fps=20) as writer:
                writer.append(first)
                writer.append(second)

            with Image.open(output_path) as animation:
                self.assertEqual(animation.n_frames, 2)
                self.assertEqual(animation.size, (32, 24))
                animation.seek(0)
                decoded_first = np.asarray(animation.convert("RGB"))
                animation.seek(1)
                decoded_second = np.asarray(animation.convert("RGB"))

        np.testing.assert_array_equal(decoded_first, first)
        np.testing.assert_array_equal(decoded_second, second)


if __name__ == "__main__":
    unittest.main()
