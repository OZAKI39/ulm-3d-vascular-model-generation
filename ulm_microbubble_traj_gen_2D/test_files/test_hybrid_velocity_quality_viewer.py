from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path

import matplotlib

matplotlib.use("Agg", force=True)
import matplotlib.pyplot as plt
import numpy as np

from ulm_microbubble_traj_gen_2D.figure_draw.quality_check.hybrid_velocity_quality_viewer import (
    HybridVelocityQualityViewer,
    _finite_element_mesh_fields,
    hybrid_velocity_quality_data_from_objects,
    load_hybrid_velocity_quality_data,
)
from ulm_microbubble_traj_gen_2D.test_files.hybrid_velocity_fixture import (
    rectangular_hybrid_velocity,
)


class HybridVelocityQualityViewerTests(unittest.TestCase):
    def _archive(self, directory: Path) -> Path:
        shape = (7, 6)
        x = np.arange(shape[0], dtype=np.float64) * 0.5
        z = np.arange(shape[1], dtype=np.float64) * 0.5
        lumen = np.zeros(shape, dtype=bool)
        lumen[1:6, 1:5] = True
        wall = lumen.copy()
        wall[2:5, 2:4] = False
        vessel_id = np.full(shape, -1, dtype=np.int32)
        vessel_id[lumen] = 3
        direction = np.zeros((*shape, 2), dtype=np.float32)
        direction[..., 0] = 1.0
        normal = np.zeros((*shape, 2), dtype=np.float32)
        normal[..., 1] = 1.0
        region = np.zeros(shape, dtype=np.uint8)
        region[lumen] = 2
        region[2:5, 2:4] = 3
        weight = np.full(shape, np.nan, dtype=np.float32)
        weight[lumen] = 0.5
        triangles = np.asarray(
            [
                [[0.5, 0.5], [1.0, 0.5], [0.5, 1.0]],
                [[1.0, 1.0], [0.5, 1.0], [1.0, 0.5]],
            ],
            dtype=np.float64,
        )
        path = directory / "velocity_and_wall_shear_field.npz"
        lumen_fraction = lumen.astype(np.float32)
        lumen_fraction[0, 2] = 0.25
        lumen_fraction[1, 2] = 0.75
        np.savez_compressed(
            path,
            x_coordinates_um=x,
            z_coordinates_um=z,
            spacing_um=np.asarray([0.5]),
            lumen_mask=lumen,
            lumen_fraction=lumen_fraction,
            wall_mask=wall,
            junction_core_mask=np.zeros(shape, dtype=bool),
            vessel_id=vessel_id,
            radius_um=np.where(lumen, 2.0, 0.0).astype(np.float32),
            flow_rate_um3_s=np.where(lumen, 30.0, 0.0).astype(np.float32),
            q2d_flow_um2_s=np.where(lumen, 3.0, 0.0).astype(np.float32),
            viscosity_mpas=np.where(lumen, 3.5, 0.0).astype(np.float32),
            direction_xz=direction,
            distance_to_centerline_um=np.where(lumen, 0.5, np.nan).astype(
                np.float32
            ),
            distance_to_wall_um=np.where(lumen, 1.0, 0.0).astype(np.float32),
            wall_normal_xz=normal,
            hybrid_velocity_region=region,
            hybrid_finite_element_weight=weight,
            hybrid_finite_element_distance_um=np.asarray(1.0),
            hybrid_regular_grid_distance_um=np.asarray(2.0),
            fem_cell_vertices_xz_um=triangles,
            continuous_wall_start_xz_um=np.asarray([[0.5, 0.5], [2.5, 0.5]]),
            continuous_wall_end_xz_um=np.asarray([[2.5, 0.5], [2.5, 2.0]]),
        )
        return path

    def test_builds_cartesian_and_dolfinx_mesh_layers(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            archive = self._archive(Path(temporary))
            data = load_hybrid_velocity_quality_data(archive)
            self.assertEqual(data.shape, (7, 6))
            self.assertIn("cartesian_cell_type", data.cartesian_fields)
            self.assertIn(
                "fem_triangle_shape_quality", data.finite_element_fields
            )
            np.testing.assert_allclose(
                data.finite_element_fields["fem_triangle_area_um2"],
                [0.125, 0.125],
            )
            np.testing.assert_allclose(
                data.finite_element_fields["fem_triangle_shape_quality"],
                [math.sqrt(3.0) / 2.0, math.sqrt(3.0) / 2.0],
            )
            cell_types = data.cartesian_fields["cartesian_cell_type"]
            self.assertEqual(int(cell_types[0, 0]), 0)
            self.assertEqual(int(cell_types[0, 2]), 1)
            self.assertEqual(int(cell_types[1, 2]), 2)
            self.assertEqual(int(cell_types[2, 2]), 3)

            viewer = HybridVelocityQualityViewer(data)
            keys = {layer.key for layer in viewer.layers}
            self.assertIn("cartesian_cell_type", keys)
            self.assertIn("cartesian_hybrid_region", keys)
            self.assertIn("cartesian_lumen_fraction", keys)
            self.assertIn("fem_triangle_region", keys)
            self.assertIn("fem_triangle_size_class", keys)
            self.assertIn("fem_triangle_area_um2", keys)
            self.assertIn("fem_triangle_minimum_edge_um", keys)
            self.assertIn("fem_triangle_maximum_edge_um", keys)
            self.assertIn("fem_triangle_shape_quality", keys)
            self.assertNotIn("vessel_id", keys)
            self.assertNotIn("direction_angle", keys)

            viewer.select_layer("cartesian_cell_type")
            self.assertTrue(viewer.image.get_visible())
            self.assertFalse(viewer.finite_element_surface.get_visible())
            viewer.axes.set_xlim(0.25, 2.75)
            viewer.axes.set_ylim(0.25, 2.25)
            viewer._refresh_dynamic_overlays()
            self.assertGreater(len(viewer.grid_collection.get_segments()), 0)
            self.assertGreater(
                len(viewer.finite_element_collection.get_segments()), 0
            )

            viewer.select_layer("fem_triangle_region")
            self.assertFalse(viewer.image.get_visible())
            self.assertTrue(viewer.finite_element_surface.get_visible())
            output = Path(temporary) / "snapshot.png"
            viewer.save_snapshot(output, dpi=72)
            self.assertTrue(output.is_file())
            self.assertGreater(output.stat().st_size, 0)
            plt.close(viewer.figure)

    def test_classifies_fem_regions_and_size_relative_to_grid(self) -> None:
        triangles = np.asarray(
            [
                [[0.0, 0.1], [0.2, 0.1], [0.0, 0.3]],
                [[0.0, 0.5], [1.0, 0.5], [0.0, 1.1]],
                [[0.0, 2.1], [2.0, 2.1], [0.0, 2.3]],
            ],
            dtype=np.float64,
        )
        fields = _finite_element_mesh_fields(
            triangles,
            np.asarray([[-10.0, 0.0]]),
            np.asarray([[10.0, 0.0]]),
            finite_element_distance_um=1.0,
            regular_grid_distance_um=2.0,
            cartesian_spacing_um=1.0,
        )
        np.testing.assert_array_equal(
            fields["fem_triangle_region"], [1, 2, 3]
        )
        np.testing.assert_array_equal(
            fields["fem_triangle_size_class"], [1, 2, 3]
        )

    def test_rejects_inconsistent_grid_array_shape(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "velocity_and_wall_shear_field.npz"
            np.savez(
                path,
                x_coordinates_um=np.asarray([0.0, 1.0]),
                z_coordinates_um=np.asarray([0.0, 1.0]),
                spacing_um=np.asarray([1.0]),
                lumen_mask=np.ones((3, 2), dtype=bool),
            )
            with self.assertRaisesRegex(ValueError, "lumen_mask"):
                load_hybrid_velocity_quality_data(path)

    def test_builds_data_directly_from_domain_and_raster_objects(self) -> None:
        from types import SimpleNamespace

        shape = (3, 2)
        domain = SimpleNamespace(
            shape=shape,
            spacing_um=1.0,
            x_coordinates_um=np.asarray([0.0, 1.0, 2.0]),
            z_coordinates_um=np.asarray([0.0, 1.0]),
        )
        raster = SimpleNamespace(
            lumen_mask=np.ones(shape, dtype=bool),
            lumen_fraction=np.ones(shape, dtype=np.float32),
            direction_xz=np.ones((*shape, 2), dtype=np.float32),
            wall_normal_xz=np.ones((*shape, 2), dtype=np.float32),
        )
        geometry = SimpleNamespace(
            solid_face_start_xz_um=np.asarray([[0.0, 0.0]]),
            solid_face_end_xz_um=np.asarray([[2.0, 0.0]]),
        )
        flow = SimpleNamespace(
            hybrid_velocity=rectangular_hybrid_velocity(domain)
        )
        raster.distance_to_wall_um = np.ones(shape, dtype=np.float32)
        data = hybrid_velocity_quality_data_from_objects(
            domain,
            raster,
            flow,
            continuous_geometry=geometry,
        )
        self.assertEqual(data.shape, shape)
        self.assertIn("cartesian_lumen_fraction", data.cartesian_fields)
        self.assertIn("cartesian_hybrid_region", data.cartesian_fields)
        self.assertIn("fem_triangle_region", data.finite_element_fields)
        self.assertEqual(data.continuous_wall_start_xz_um.shape, (1, 2))


if __name__ == "__main__":
    unittest.main()
