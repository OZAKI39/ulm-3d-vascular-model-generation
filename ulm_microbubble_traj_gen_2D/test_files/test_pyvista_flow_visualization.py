"""PyVista/VTK CFD 流场网格、连续流线和渲染产物的独立验证。"""

from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from ulm_microbubble_traj_gen_2D.utils.visualization.vtk.pyvista_flow import (
    render_cfd_flow_field,
    render_cfd_flow_fields,
    validate_cfd_flow_dependencies,
)
from ulm_microbubble_traj_gen_2D.utils.visualization.vtk.pyvista_wall_shear import render_wall_shear_visualization
from ulm_microbubble_traj_gen_2D.utils.visualization.apps.trame_flow_viewer import build_probe_plotter, sample_native_field
from ulm_microbubble_traj_gen_2D.utils.visualization.vtk.vtk_flow_grid import (
    LUMEN_ARRAY,
    PRESSURE_ARRAY,
    SPEED_ARRAY,
    VELOCITY_ARRAY,
    WALL_SHEAR_ARRAY,
    WALL_SHEAR_DISPLAY_MASK_ARRAY,
    build_vtk_stage_grid,
)
from ulm_microbubble_traj_gen_2D.utils.visualization.vtk.vtk_streamlines import (
    StreamlineContinuityError,
    _globally_observable_vessel_ids,
    _is_ordered_subsequence,
    _ordered_covered_vessels_outside_junctions,
    _outlet_vessel_map,
    _target_open_face_intersection,
    trace_root_to_outlets,
)


PYVISTA_AVAILABLE = importlib.util.find_spec("pyvista") is not None


class DependencyValidationTests(unittest.TestCase):
    """即使环境没有 PyVista，也必须给出可执行的安装提示。"""

    def test_dependency_error_contains_install_command(self) -> None:
        real_find_spec = importlib.util.find_spec

        def fake_find_spec(name: str):
            if name == "pyvista":
                return None
            return real_find_spec(name)

        with patch("importlib.util.find_spec", side_effect=fake_find_spec):
            with self.assertRaisesRegex(ImportError, r"pyvista\[jupyter\]"):
                validate_cfd_flow_dependencies()

    def test_flow_renderer_facade_keeps_two_path_return_contract(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output_dir = Path(directory)

            def fake_stage_renderer(*_args, stage: str, **_kwargs):
                return SimpleNamespace(html_path=output_dir / f"{stage}.html")

            with (
                patch(
                    "ulm_microbubble_traj_gen.utils.visualization.vtk.pyvista_flow.validate_cfd_flow_dependencies"
                ),
                patch(
                    "ulm_microbubble_traj_gen.utils.flow.flow_boundaries.build_flux_boundaries",
                    return_value=SimpleNamespace(),
                ),
                patch(
                    "ulm_microbubble_traj_gen.utils.visualization.vtk.pyvista_flow.render_cfd_flow_field",
                    side_effect=fake_stage_renderer,
                ),
            ):
                paths = render_cfd_flow_fields(
                    SimpleNamespace(),
                    SimpleNamespace(),
                    SimpleNamespace(),
                    output_dir,
                    vessels=(),
                    effective_thickness_um=10.0,
                    boundary_depth_cells=1.5,
                )
        self.assertEqual(paths, (output_dir / "initial.html", output_dir / "final.html"))


@unittest.skipUnless(PYVISTA_AVAILABLE, "PyVista is not installed")
class PyVistaFlowVisualizationTests(unittest.TestCase):
    """使用可解析的直管场检查关键物理与几何约束。"""

    def setUp(self) -> None:
        nx, nz = 64, 32
        spacing = 2.0
        x_um = 10.0 + spacing * np.arange(nx, dtype=float)
        z_um = -20.0 + spacing * np.arange(nz, dtype=float)

        lumen = np.zeros((nx, nz), dtype=bool)
        lumen[2:-2, 9:23] = True
        interior = np.zeros_like(lumen)
        interior[3:-3, 10:22] = True
        wall = lumen & ~interior

        velocity = np.zeros((nx, nz, 2), dtype=np.float32)
        velocity[..., 0] = 100.0
        velocity[~lumen] = 0.0
        pressure = np.broadcast_to(np.linspace(1.0, -1.0, nx)[:, None], (nx, nz)).copy()
        pressure[~lumen] = 0.0
        wall_shear = np.broadcast_to(np.linspace(0.05, 0.85, nz)[None, :], (nx, nz)).copy()
        wall_shear[~lumen] = 0.0
        vessel_id = np.full((nx, nz), -1, dtype=np.int32)
        vessel_id[lumen] = 0

        inlet_label = np.zeros((nx, nz), dtype=np.int32)
        outlet_label = np.zeros((nx, nz), dtype=np.int32)
        inlet_label[2, 9:23] = 1
        outlet_label[-3, 9:23] = 1
        boundary_normal = np.zeros((nx, nz, 2), dtype=np.float32)
        boundary_normal[2, 9:23, 0] = -1.0
        boundary_normal[-3, 9:23, 0] = 1.0

        self.domain = SimpleNamespace(
            shape=(nx, nz),
            spacing_um=spacing,
            x_coordinates_um=x_um,
            z_coordinates_um=z_um,
        )
        self.raster = SimpleNamespace(
            lumen_mask=lumen,
            lumen_fraction=lumen.astype(np.float32),
            wall_mask=wall,
            vessel_id=vessel_id,
            distance_to_wall_um=np.where(lumen, spacing, 0.0),
            junction_core_mask=np.zeros_like(lumen),
        )
        self.flow = SimpleNamespace(
            velocity_xz_um_s=velocity,
            initial_velocity_xz_um_s=0.8 * velocity,
            pressure=pressure,
            wall_shear_stress_pa=wall_shear.astype(np.float32),
            inlet_label=inlet_label,
            outlet_label=outlet_label,
            boundary_normal_xz=boundary_normal,
            boundary_velocity_xz_um_s=velocity.copy(),
        )
        outlet_cells = np.argwhere(outlet_label == 1).astype(np.int32)
        self.open_boundaries = SimpleNamespace(
            open_face_cell_ij=outlet_cells,
            open_face_normal_xz=np.tile(np.asarray([[1.0, 0.0]], dtype=np.float32), (outlet_cells.shape[0], 1)),
            open_face_length_um=np.full(outlet_cells.shape[0], spacing, dtype=float),
            open_face_label=np.ones(outlet_cells.shape[0], dtype=np.int32),
            open_face_kind=np.ones(outlet_cells.shape[0], dtype=np.int8),
        )

    def test_cell_grid_uses_physical_coordinates_and_nan_solid_mask(self) -> None:
        grid = build_vtk_stage_grid(self.domain, self.raster, self.flow, stage="initial", include_lic=False)
        image = grid.image_grid
        nx, nz = self.domain.shape
        spacing = self.domain.spacing_um

        self.assertEqual(image.dimensions, (nx + 1, nz + 1, 1))
        self.assertEqual(image.n_cells, nx * nz)
        np.testing.assert_allclose(
            image.origin,
            (
                self.domain.x_coordinates_um[0] - 0.5 * spacing,
                self.domain.z_coordinates_um[0] - 0.5 * spacing,
                0.0,
            ),
        )

        # VTK 的 X 索引变化最快，因此 (ix, iz) 的 cell id 是 ix + nx * iz。
        inside_ix, inside_iz = 12, 15
        inside_id = inside_ix + nx * inside_iz
        outside_id = 0
        centers = image.cell_centers().points
        np.testing.assert_allclose(
            centers[inside_id, :2],
            (self.domain.x_coordinates_um[inside_ix], self.domain.z_coordinates_um[inside_iz]),
        )
        self.assertTrue(np.isnan(np.asarray(image.cell_data[SPEED_ARRAY])[outside_id]))
        self.assertTrue(np.isfinite(np.asarray(image.cell_data[SPEED_ARRAY])[inside_id]))
        self.assertTrue(np.isnan(np.asarray(image.cell_data[PRESSURE_ARRAY])[outside_id]))
        self.assertEqual(float(np.asarray(image.cell_data[PRESSURE_ARRAY])[inside_id]), 0.0)
        self.assertEqual(str(np.asarray(image.field_data["pressure_semantics"])[0]), "zero_reference_unsolved")
        self.assertTrue(np.isnan(np.asarray(image.cell_data[WALL_SHEAR_ARRAY])[inside_id]))
        self.assertEqual(
            str(np.asarray(image.field_data["wall_shear_semantics"])[0]),
            "not_computed_for_initial_velocity",
        )
        np.testing.assert_allclose(np.asarray(image.cell_data[VELOCITY_ARRAY])[inside_id], (80.0, 0.0, 0.0))

        side_wall_id = 12 + nx * 9
        inlet_cap_id = 2 + nx * 12
        outlet_cap_id = (nx - 3) + nx * 12
        display_mask = np.asarray(image.cell_data[WALL_SHEAR_DISPLAY_MASK_ARRAY], dtype=bool)
        self.assertTrue(display_mask[side_wall_id])
        self.assertFalse(display_mask[inlet_cap_id])
        self.assertFalse(display_mask[outlet_cap_id])
        self.assertFalse(display_mask[inside_id])

        self.assertEqual(grid.fluid_grid.n_cells, int(np.count_nonzero(self.raster.lumen_mask)))
        # wall_mask 是近壁流体带，不是固体；这些单元必须仍留在 threshold 后的网格中。
        self.assertGreater(int(np.count_nonzero(self.raster.wall_mask)), 0)
        self.assertEqual(
            int(np.count_nonzero(np.asarray(grid.fluid_grid.cell_data["wall_mask"]))),
            int(np.count_nonzero(self.raster.wall_mask)),
        )

    def test_formal_streamlines_are_single_continuous_root_to_outlet_polylines(self) -> None:
        grid = build_vtk_stage_grid(self.domain, self.raster, self.flow, stage="final", include_lic=False)
        nx = self.domain.shape[0]
        side_wall_id = 12 + nx * 9
        self.assertAlmostEqual(
            float(np.asarray(grid.image_grid.cell_data[WALL_SHEAR_ARRAY])[side_wall_id]),
            float(self.flow.wall_shear_stress_pa[12, 9]),
        )
        self.assertEqual(
            str(np.asarray(grid.image_grid.field_data["wall_shear_semantics"])[0]),
            "accepted_velocity_in_plane_wall_shear_stress_magnitude_proxy",
        )
        self.assertEqual(
            str(np.asarray(grid.image_grid.field_data["wall_shear_source_stage"])[0]),
            "final_accepted_velocity",
        )
        result = trace_root_to_outlets(
            self.domain,
            self.raster,
            self.flow,
            grid,
            open_boundaries=self.open_boundaries,
        )

        self.assertEqual(result.expected_outlet_labels, (1,))
        self.assertEqual(result.reached_outlet_labels, (1,))
        self.assertEqual(result.missing_outlet_labels, ())
        self.assertEqual(result.covered_vessel_ids, (0,))
        self.assertGreater(result.formal_lines.n_lines, 0)
        self.assertTrue(np.all(result.formal_lines.cell_data["is_complete_root_to_outlet"] == 1))
        self.assertTrue(np.all(result.formal_lines.cell_data["destination_outlet_id"] == 1))

        connectivity = np.asarray(result.formal_lines.lines, dtype=np.int64)
        cursor = 0
        for _ in range(result.formal_lines.n_lines):
            count = int(connectivity[cursor])
            point_ids = connectivity[cursor + 1 : cursor + 1 + count]
            cursor += count + 1
            points = result.formal_lines.points[point_ids]
            self.assertGreaterEqual(count, 2)
            self.assertTrue(np.all(np.isfinite(points)))
            self.assertLessEqual(float(points[0, 0]), self.domain.x_coordinates_um[3])
            self.assertGreaterEqual(float(points[-1, 0]), self.domain.x_coordinates_um[-4])
            self.assertTrue(np.all(np.diff(points[:, 0]) >= -1.0e-7))

    def test_topology_helpers_preserve_order_and_use_authoritative_outlet_ids(self) -> None:
        vessel_id = self.raster.vessel_id.copy()
        for start, stop, value in ((2, 22, 0), (22, 43, 1), (43, -2, 2)):
            region = vessel_id[start:stop]
            region[self.raster.lumen_mask[start:stop]] = value
        raster = SimpleNamespace(**{**self.raster.__dict__, "vessel_id": vessel_id})
        center_z = float(self.domain.z_coordinates_um[self.domain.shape[1] // 2])
        points = np.column_stack(
            (
                self.domain.x_coordinates_um[2:-2],
                np.full(self.domain.shape[0] - 4, center_z),
            )
        )

        self.assertEqual(
            _ordered_covered_vessels_outside_junctions(points, self.domain, raster),
            (0, 1, 2),
        )
        self.assertEqual(_globally_observable_vessel_ids(raster, self.domain), frozenset((0, 1, 2)))

        centerline_distance = np.zeros_like(vessel_id, dtype=float)
        centerline_distance[vessel_id == 1] = self.domain.spacing_um
        edge_only_raster = SimpleNamespace(
            **{**raster.__dict__, "distance_to_centerline_um": centerline_distance}
        )
        self.assertEqual(
            _globally_observable_vessel_ids(edge_only_raster, self.domain),
            frozenset((0, 2)),
        )
        self.assertTrue(_is_ordered_subsequence((0, 2), (0, 1, 2)))
        self.assertFalse(_is_ordered_subsequence((0, 2, 1), (0, 1, 2)))
        self.assertFalse(_is_ordered_subsequence((0, 1, 0), (0, 1, 2)))

        authoritative = SimpleNamespace(outlet_ids=(42,))
        outlet_map = _outlet_vessel_map(self.flow.outlet_label, vessel_id, authoritative)
        self.assertEqual(outlet_map, {1: 42})

    def test_wrong_target_outlet_is_diagnostic_not_formal(self) -> None:
        outlet_label = self.flow.outlet_label.copy()
        # label 2 位于管内而不是开放出口；forward 路径会继续到 label 1，不能把它
        # 截断或改标签来伪造“已覆盖第二个出口”。
        outlet_label[self.domain.shape[0] // 2, self.domain.shape[1] // 2] = 2
        flow = SimpleNamespace(**{**self.flow.__dict__, "outlet_label": outlet_label})
        grid = build_vtk_stage_grid(self.domain, self.raster, flow, stage="final", include_lic=False)

        with self.assertRaises(StreamlineContinuityError) as raised:
            trace_root_to_outlets(
                self.domain,
                self.raster,
                flow,
                grid,
                open_boundaries=self.open_boundaries,
            )
        result = raised.exception.result
        self.assertIn(2, result.missing_outlet_labels)
        self.assertGreater(result.diagnostic_lines.n_lines, 0)
        if result.formal_lines.n_lines:
            self.assertFalse(np.any(result.formal_lines.cell_data["destination_outlet_id"] == 2))

    @unittest.skipUnless(
        importlib.util.find_spec("trame_vtk") is not None and importlib.util.find_spec("trame_vuetify") is not None,
        "PyVista HTML export dependencies are not installed",
    )
    def test_renderer_writes_complete_native_and_html_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            artifacts = render_cfd_flow_field(
                self.domain,
                self.raster,
                self.flow,
                Path(directory),
                stage="final",
                open_boundaries=self.open_boundaries,
            )
            wall_shear_artifacts = render_wall_shear_visualization(
                Path(directory),
                domain=self.domain,
                raster=self.raster,
                flow=self.flow,
            )
            expected = (
                artifacts.html_path,
                artifacts.preview_path,
                artifacts.field_vti_path,
                artifacts.formal_streamlines_vtp_path,
                artifacts.diagnostic_streamlines_vtp_path,
                artifacts.continuity_csv_path,
                wall_shear_artifacts.html_path,
                wall_shear_artifacts.preview_path,
            )
            for path in expected:
                self.assertTrue(path.is_file(), path)
                self.assertGreater(path.stat().st_size, 0, path)
            self.assertGreater(artifacts.html_path.stat().st_size, 100_000)
            html = artifacts.html_path.read_text(encoding="utf-8")
            self.assertIn("<title>Converged microvascular CFD field</title>", html)
            self.assertIn('id="cfd-global-title"', html)
            wall_shear_html = wall_shear_artifacts.html_path.read_text(encoding="utf-8")
            self.assertIn("<title>Final 2D in-plane wall-shear-stress proxy</title>", wall_shear_html)

            import pyvista as pv

            field = pv.read(artifacts.field_vti_path)
            streamlines = pv.read(artifacts.formal_streamlines_vtp_path)
            self.assertIn(LUMEN_ARRAY, field.cell_data)
            self.assertIn(VELOCITY_ARRAY, field.cell_data)
            self.assertIn(WALL_SHEAR_ARRAY, field.cell_data)
            self.assertIn(WALL_SHEAR_DISPLAY_MASK_ARRAY, field.cell_data)
            self.assertIn("boundary_normal_xz", field.cell_data)
            self.assertIn("is_complete_root_to_outlet", streamlines.cell_data)
            self.assertTrue(np.all(streamlines.cell_data["is_complete_root_to_outlet"] == 1))
            probe = sample_native_field(
                field,
                np.asarray([self.domain.x_coordinates_um[12], self.domain.z_coordinates_um[15], 0.0]),
            )
            self.assertIsNotNone(probe)
            self.assertAlmostEqual(float(probe["vx_um_s"]), 100.0)
            self.assertAlmostEqual(float(probe["vz_um_s"]), 0.0)
            self.assertAlmostEqual(
                float(probe["wall_shear_stress_pa"]),
                float(self.flow.wall_shear_stress_pa[12, 15]),
            )
            self.assertEqual(
                str(probe["wall_shear_semantics"]),
                "accepted_velocity_in_plane_wall_shear_stress_magnitude_proxy",
            )
            self.assertEqual(int(probe["vessel_id"]), 0)

            wall_shear_plotter, _ = build_probe_plotter(
                Path(directory),
                stage="final",
                view="wall-shear",
            )
            wall_shear_plotter.close()
            with self.assertRaisesRegex(ValueError, "only available for --stage final"):
                build_probe_plotter(Path(directory), stage="initial", view="wall-shear")

    def test_side_wall_near_outlet_is_not_an_open_face_hit(self) -> None:
        outlet_cell = np.argwhere(self.flow.outlet_label == 1)[len(np.argwhere(self.flow.outlet_label == 1)) // 2]
        ix, iz = (int(value) for value in outlet_cell)
        cell_center = np.asarray(
            [self.domain.x_coordinates_um[ix], self.domain.z_coordinates_um[iz]],
            dtype=float,
        )
        spacing = float(self.domain.spacing_um)

        cap_segment = np.asarray(
            [
                cell_center + np.asarray([0.00, 0.00]),
                cell_center + np.asarray([0.35 * spacing, 0.00]),
            ]
        )
        cap_hit = _target_open_face_intersection(
            cap_segment,
            np.asarray([100.0, 0.0]),
            1,
            self.flow.outlet_label,
            self.raster.lumen_mask,
            self.flow.boundary_normal_xz,
            self.domain,
            self.open_boundaries,
            maximum_extension_um=spacing,
        )
        self.assertIsNotNone(cap_hit)
        self.assertAlmostEqual(float(cap_hit[0]), float(cell_center[0] + 0.5 * spacing))

        side_segment = np.asarray(
            [
                cell_center + np.asarray([0.00, 0.00]),
                cell_center + np.asarray([0.00, 0.35 * spacing]),
            ]
        )
        side_hit = _target_open_face_intersection(
            side_segment,
            np.asarray([0.0, 100.0]),
            1,
            self.flow.outlet_label,
            self.raster.lumen_mask,
            self.flow.boundary_normal_xz,
            self.domain,
            self.open_boundaries,
            maximum_extension_um=spacing,
        )
        self.assertIsNone(side_hit)

if __name__ == "__main__":
    unittest.main()
