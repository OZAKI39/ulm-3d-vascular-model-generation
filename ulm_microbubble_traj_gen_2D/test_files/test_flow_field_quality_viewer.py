from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from ulm_microbubble_traj_gen_2D.figure_draw.quality_check.flow_field_quality_viewer import (
    FlowFieldQualityViewer,
    load_flow_field_quality_data,
)


class FlowFieldQualityViewerTests(unittest.TestCase):
    def _archive(
        self,
        directory: Path,
        *,
        include_optional_fields: bool = True,
    ) -> Path:
        shape = (8, 7)
        spacing = 0.5
        x = np.arange(shape[0], dtype=np.float64) * spacing
        z = np.arange(shape[1], dtype=np.float64) * spacing
        lumen = np.zeros(shape, dtype=bool)
        lumen[1:7, 1:6] = True

        velocity = np.zeros((*shape, 2), dtype=np.float32)
        velocity[..., 0] = np.arange(shape[0], dtype=np.float32)[:, None]
        velocity[..., 1] = -np.arange(shape[1], dtype=np.float32)[None, :]
        velocity[~lumen] = 0.0
        speed = np.linalg.norm(velocity, axis=2).astype(np.float32)
        wall_shear = np.where(lumen, 0.25 + speed, 0.0).astype(np.float32)

        data: dict[str, np.ndarray] = {
            "x_coordinates_um": x,
            "z_coordinates_um": z,
            "spacing_um": np.asarray([spacing]),
            "lumen_mask": lumen,
            "velocity_xz_um_s": velocity,
            "speed_um_s": speed,
            "wall_shear_stress_pa": wall_shear,
        }
        if include_optional_fields:
            initial_velocity = 0.8 * velocity
            boundary_weight = np.zeros(shape, dtype=np.float32)
            boundary_weight[1, 2:5] = 1.0
            boundary_velocity = np.zeros_like(velocity)
            boundary_velocity[1, 2:5, 0] = 4.0
            open_flux = np.zeros(shape, dtype=np.float32)
            open_flux[1, 2:5] = -2.0
            data.update(
                {
                    "pressure": np.where(lumen, 35.0 + x[:, None], 0.0).astype(
                        np.float32
                    ),
                    "divergence_s_inv": np.where(
                        lumen, x[:, None] - z[None, :], 0.0
                    ).astype(np.float32),
                    "wall_penetration_um_s": np.zeros(shape, dtype=np.float32),
                    "initial_velocity_xz_um_s": initial_velocity,
                    "initial_speed_um_s": np.linalg.norm(
                        initial_velocity, axis=2
                    ).astype(np.float32),
                    "boundary_velocity_xz_um_s": boundary_velocity,
                    "boundary_weight": boundary_weight,
                    "open_boundary_flux_um2_s": open_flux,
                    "continuous_wall_start_xz_um": np.asarray(
                        [[0.5, 0.5], [3.0, 0.5]], dtype=np.float64
                    ),
                    "continuous_wall_end_xz_um": np.asarray(
                        [[3.0, 0.5], [3.0, 2.5]], dtype=np.float64
                    ),
                    "continuous_open_section_point_xz_um": np.asarray(
                        [[0.5, 1.5], [3.0, 1.5]], dtype=np.float64
                    ),
                    "continuous_open_section_tangent_xz": np.asarray(
                        [[0.0, 1.0], [0.0, 1.0]], dtype=np.float64
                    ),
                    "continuous_open_section_half_width_um": np.asarray(
                        [0.5, 0.5], dtype=np.float64
                    ),
                    "continuous_open_section_kind": np.asarray(
                        [-1, 1], dtype=np.int8
                    ),
                    "solver_sampled_pressure_unit": np.asarray(["mmHg"]),
                    "solver_sampled_pressure_semantics": np.asarray(
                        ["gauge_relative_to_pinned_pressure_dof"]
                    ),
                }
            )

        path = directory / "velocity_and_wall_shear_field.npz"
        np.savez_compressed(path, **data)
        return path

    def test_loads_solver_fields_and_builds_interactive_layers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            archive = self._archive(Path(temporary))
            data = load_flow_field_quality_data(archive)

            self.assertEqual(data.shape, (8, 7))
            self.assertEqual(data.pressure_unit, "mmHg")
            self.assertIn("final_velocity", data.vector_fields)
            self.assertIn("initial_velocity", data.vector_fields)
            self.assertIn("speed_change_um_s", data.scalar_fields)
            self.assertEqual(int(np.count_nonzero(data.boundary_mask)), 3)

            viewer = FlowFieldQualityViewer(data)
            keys = {layer.key for layer in viewer.layers}
            expected = {
                "speed_um_s",
                "velocity_x_um_s",
                "velocity_z_um_s",
                "pressure",
                "wall_shear_stress_pa",
                "divergence_s_inv",
                "wall_penetration_um_s",
                "speed_change_um_s",
                "initial_speed_um_s",
                "initial_velocity_x_um_s",
                "initial_velocity_z_um_s",
                "boundary_speed_um_s",
                "open_boundary_flux_um2_s",
            }
            self.assertEqual(keys, expected)

            viewer.select_layer("divergence_s_inv")
            lower, upper = viewer.image.get_clim()
            self.assertAlmostEqual(lower, -upper)
            self.assertTrue(viewer.velocity_quiver.get_visible())

            viewer.select_layer("boundary_speed_um_s")
            self.assertEqual(viewer.current_layer.vector_key, "boundary_velocity")
            viewer._show_cell_grid = True
            viewer.axes.set_xlim(0.25, 3.25)
            viewer.axes.set_ylim(0.25, 2.75)
            viewer._refresh_dynamic_overlays()
            self.assertGreater(len(viewer.grid_collection.get_segments()), 0)

            output = Path(temporary) / "flow_snapshot.png"
            viewer.save_snapshot(output, dpi=72)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 0)
            plt.close(viewer.figure)

    def test_accepts_minimum_legacy_flow_archive(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            archive = self._archive(
                Path(temporary), include_optional_fields=False
            )
            data = load_flow_field_quality_data(archive)
            viewer = FlowFieldQualityViewer(data)
            self.assertEqual(
                [layer.key for layer in viewer.layers],
                [
                    "speed_um_s",
                    "velocity_x_um_s",
                    "velocity_z_um_s",
                    "wall_shear_stress_pa",
                ],
            )
            self.assertEqual(
                data.continuous_wall_start_xz_um.shape,
                (0, 2),
            )
            self.assertFalse(viewer.wall_collection.get_visible())
            self.assertFalse(viewer.open_section_collection.get_visible())
            plt.close(viewer.figure)

    def test_rejects_inconsistent_velocity_shape(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            path = directory / "velocity_and_wall_shear_field.npz"
            shape = (3, 2)
            np.savez(
                path,
                x_coordinates_um=np.arange(shape[0], dtype=float),
                z_coordinates_um=np.arange(shape[1], dtype=float),
                spacing_um=np.asarray([1.0]),
                lumen_mask=np.ones(shape, dtype=bool),
                velocity_xz_um_s=np.zeros(shape, dtype=np.float32),
                speed_um_s=np.ones(shape, dtype=np.float32),
                wall_shear_stress_pa=np.ones(shape, dtype=np.float32),
            )
            with self.assertRaisesRegex(ValueError, "velocity_xz_um_s"):
                load_flow_field_quality_data(path)


if __name__ == "__main__":
    unittest.main()
