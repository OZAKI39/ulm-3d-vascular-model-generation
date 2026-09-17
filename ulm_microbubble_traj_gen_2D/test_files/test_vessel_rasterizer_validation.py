"""Validation tests for exact Shapely vessel-cell coverage."""

from __future__ import annotations

import math
import sys
import unittest
from dataclasses import replace
from pathlib import Path

import numpy as np
import shapely
from scipy import ndimage
from shapely.geometry import Point
from shapely.ops import unary_union

if __package__ in (None, ""):
    sys.path.append(str(Path(__file__).resolve().parents[2]))

from ulm_microbubble_traj_gen_2D.utils.core.config import DomainConfig
from ulm_microbubble_traj_gen_2D.utils.flow.connectivity import _hole_cell_count
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.grid_domain import build_domain_from_vessels
from ulm_microbubble_traj_gen_2D.utils.core.types import GridDomain
from ulm_microbubble_traj_gen_2D.utils.geometry.vessel_rasterizer import (
    _cleanup_lumen_mask,
    _exact_polygon_cell_coverage,
    rasterize_vessels,
)
from ulm_vascular_model_generator.utils.core.models import Vessel


class VesselRasterizerTests(unittest.TestCase):
    def test_cleanup_never_fills_a_geometry_rejected_enclosed_cell(self) -> None:
        lumen = np.ones((5, 5), dtype=bool)
        lumen[2, 2] = False

        cleaned = _cleanup_lumen_mask(lumen)

        self.assertFalse(cleaned[2, 2])
        self.assertTrue(np.all(cleaned <= lumen))

    def test_hole_count_ignores_enclosed_exterior_grid_centres(self) -> None:
        coordinates = np.arange(5, dtype=np.float64)
        domain = GridDomain(
            origin_um=np.asarray([0.0, 0.0, 0.0]),
            spacing_um=1.0,
            shape=(5, 5),
            fixed_y_um=0.0,
            x_coordinates_um=coordinates,
            z_coordinates_um=coordinates,
        )
        lumen = np.ones(domain.shape, dtype=bool)
        lumen[2, 2] = False

        class ContinuousGeometry:
            @staticmethod
            def contains_xz_um(points):
                return np.zeros(np.asarray(points).shape[0], dtype=bool)

        self.assertEqual(_hole_cell_count(lumen), 1)
        self.assertEqual(
            _hole_cell_count(
                lumen,
                domain=domain,
                continuous_geometry=ContinuousGeometry(),
            ),
            0,
        )

    def test_sparse_tiled_coverage_matches_full_bbox_intersections(self) -> None:
        spacing = 1.0
        x = -10.0 + spacing * np.arange(101)
        z = -10.0 + spacing * np.arange(91)
        domain = GridDomain(
            origin_um=np.asarray([x[0], 0.0, z[0]]),
            spacing_um=spacing,
            shape=(x.size, z.size),
            fixed_y_um=0.0,
            x_coordinates_um=x,
            z_coordinates_um=z,
        )
        polygon = unary_union(
            [
                Point(-3.25, -1.75).buffer(3.7, quad_segs=32),
                Point(73.5, 65.25).buffer(5.2, quad_segs=32),
            ]
        )

        actual = _exact_polygon_cell_coverage(polygon, domain)

        expected = np.zeros(domain.shape, dtype=np.float64)
        half = 0.5 * spacing
        minimum_x, minimum_z, maximum_x, maximum_z = polygon.bounds
        ix = np.flatnonzero(
            (x + half >= minimum_x) & (x - half <= maximum_x)
        )
        iz = np.flatnonzero(
            (z + half >= minimum_z) & (z - half <= maximum_z)
        )
        center_x, center_z = np.meshgrid(x[ix], z[iz], indexing="ij")
        cells = shapely.box(
            center_x - half,
            center_z - half,
            center_x + half,
            center_z + half,
        )
        expected[np.ix_(ix, iz)] = np.clip(
            np.asarray(
                shapely.area(shapely.intersection(cells, polygon)),
                dtype=np.float64,
            )
            / spacing**2,
            0.0,
            1.0,
        )

        np.testing.assert_array_equal(
            actual,
            np.ascontiguousarray(expected, dtype=np.float32),
        )

    def test_exact_cell_fractions_preserve_continuous_lumen_area(self) -> None:
        vessel = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([0.0, 0.0, 0.0]),
            x_d=np.asarray([10.0, 0.0, 7.0]),
            radius=2.3,
            flow_rate=1.0,
        )
        cfg = DomainConfig(
            grid_spacing_um=1.0,
            padding_um=4.0,
            min_lumen_radius_cells=0.5,
            min_resolved_diameter_cells=0.0,
            max_grid_cells=100_000,
        )
        domain = build_domain_from_vessels([vessel], cfg)
        geometry = build_continuous_vessel_geometry([vessel], domain)
        raster = rasterize_vessels(
            [vessel],
            domain,
            cfg,
            effective_thickness_um=1.0,
            continuous_geometry=geometry,
            dynamic_viscosity_mpas=3.18,
        )

        fraction = np.asarray(raster.lumen_fraction, dtype=np.float64)
        self.assertEqual(raster.lumen_fraction.dtype, np.float32)
        self.assertTrue(raster.lumen_fraction.flags.c_contiguous)
        self.assertTrue(np.all((fraction >= 0.0) & (fraction <= 1.0)))
        self.assertTrue(np.any((fraction > 0.0) & (fraction < 1.0)))
        rasterized_area = float(np.sum(fraction)) * domain.spacing_um**2
        self.assertAlmostEqual(
            rasterized_area,
            float(geometry.lumen_polygon.area),
            delta=2.0e-6 * float(geometry.lumen_polygon.area),
        )
        x, z = np.meshgrid(
            domain.x_coordinates_um,
            domain.z_coordinates_um,
            indexing="ij",
        )
        centre_inside = geometry.contains_xz_um(
            np.column_stack((x.ravel(), z.ravel()))
        ).reshape(domain.shape)
        np.testing.assert_array_equal(
            raster.lumen_mask,
            (fraction >= 0.5) & centre_inside,
        )
        rows = np.argwhere(raster.lumen_mask)
        centres = np.column_stack(
            (
                domain.x_coordinates_um[rows[:, 0]],
                domain.z_coordinates_um[rows[:, 1]],
            )
        )
        self.assertTrue(np.all(geometry.contains_xz_um(centres)))
        np.testing.assert_allclose(
            raster.viscosity_mpas[raster.lumen_mask],
            3.18,
            rtol=0.0,
            atol=1.0e-6,
        )

    def test_grid_wall_distance_and_normal_are_sampled_from_continuous_wall(
        self,
    ) -> None:
        vessel = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([0.0, 0.0, 0.0]),
            x_d=np.asarray([10.0, 0.0, 7.0]),
            radius=2.3,
            flow_rate=1.0,
        )
        cfg = DomainConfig(
            grid_spacing_um=1.0,
            padding_um=4.0,
            min_lumen_radius_cells=0.5,
            min_resolved_diameter_cells=0.0,
            max_grid_cells=100_000,
        )
        domain = build_domain_from_vessels([vessel], cfg)
        geometry = build_continuous_vessel_geometry([vessel], domain)
        raster = rasterize_vessels(
            [vessel],
            domain,
            cfg,
            effective_thickness_um=1.0,
            continuous_geometry=geometry,
        )

        rows = np.argwhere(raster.lumen_mask)
        points = np.column_stack(
            (
                domain.x_coordinates_um[rows[:, 0]],
                domain.z_coordinates_um[rows[:, 1]],
            )
        )
        exact = geometry.exact_solid_wall_state_xz_um_accelerated(points)
        np.testing.assert_allclose(
            raster.distance_to_wall_um[rows[:, 0], rows[:, 1]],
            exact.distance_um,
            atol=2.0e-6,
        )
        np.testing.assert_allclose(
            raster.wall_normal_xz[rows[:, 0], rows[:, 1]],
            exact.inward_normal_xz,
            atol=2.0e-6,
        )
        expected_wall = raster.lumen_mask & (
            raster.distance_to_wall_um
            <= math.sqrt(0.5) * domain.spacing_um + 1.0e-6
        )
        np.testing.assert_array_equal(raster.wall_mask, expected_wall)
        wall_normals = np.asarray(raster.wall_normal_xz[raster.wall_mask])
        self.assertTrue(
            np.any(
                (np.abs(wall_normals[:, 0]) > 0.1)
                & (np.abs(wall_normals[:, 1]) > 0.1)
            )
        )

    def test_grid_domain_rejects_unreliable_settings(self) -> None:
        vessel = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([0.0, 0.0, 0.0]),
            x_d=np.asarray([10.0, 0.0, 10.0]),
            radius=2.0,
            flow_rate=1.0,
        )
        cfg = DomainConfig(
            grid_spacing_um=1.0,
            padding_um=1.0,
            min_lumen_radius_cells=0.5,
            min_resolved_diameter_cells=8.0,
            max_grid_cells=10_000,
        )

        for changes in (
            {"grid_spacing_um": 0.0},
            {"grid_spacing_um": float("nan")},
            {"padding_um": -1.0},
            {"padding_um": float("nan")},
            {"max_grid_cells": 0},
        ):
            with self.subTest(changes=changes):
                with self.assertRaises(ValueError):
                    build_domain_from_vessels([vessel], replace(cfg, **changes))

    def test_continuous_bifurcation_has_no_single_cell_holes(self) -> None:
        vessels = [
            Vessel(
                vid=0,
                parent_id=-1,
                children=[1, 2],
                x_p=np.asarray([1063.769, 0.0, 1466.891]),
                x_d=np.asarray([867.181, 0.0, 1355.275]),
                radius=12.133,
                flow_rate=1.0,
            ),
            Vessel(
                vid=1,
                parent_id=0,
                children=[],
                x_p=np.asarray([867.181, 0.0, 1355.275]),
                x_d=np.asarray([752.454, 0.0, 1214.017]),
                radius=7.643,
                flow_rate=1.0,
            ),
            Vessel(
                vid=2,
                parent_id=0,
                children=[],
                x_p=np.asarray([867.181, 0.0, 1355.275]),
                x_d=np.asarray([566.045, 0.0, 1310.095]),
                radius=11.024,
                flow_rate=1.0,
            ),
        ]
        cfg = DomainConfig(
            grid_spacing_um=1.9,
            padding_um=25.0,
            min_lumen_radius_cells=0.5,
            min_resolved_diameter_cells=8.0,
            max_grid_cells=200000,
        )
        domain = build_domain_from_vessels(vessels, cfg)
        continuous_geometry = build_continuous_vessel_geometry(vessels, domain)
        raster = rasterize_vessels(
            vessels,
            domain,
            cfg,
            effective_thickness_um=10.0,
            continuous_geometry=continuous_geometry,
        )
        holes = ndimage.binary_fill_holes(
            raster.lumen_mask,
            structure=ndimage.generate_binary_structure(2, 1),
        ) & ~raster.lumen_mask

        self.assertEqual(int(np.count_nonzero(holes)), 0)

    def test_rejects_invalid_effective_thickness(self) -> None:
        vessel = Vessel(
            vid=0,
            parent_id=-1,
            children=[],
            x_p=np.asarray([0.0, 0.0, 0.0]),
            x_d=np.asarray([10.0, 0.0, 10.0]),
            radius=2.0,
            flow_rate=1.0,
        )
        cfg = DomainConfig(
            grid_spacing_um=1.0,
            padding_um=1.0,
            min_lumen_radius_cells=0.5,
            min_resolved_diameter_cells=8.0,
            max_grid_cells=10_000,
        )
        domain = build_domain_from_vessels([vessel], cfg)
        continuous_geometry = build_continuous_vessel_geometry([vessel], domain)

        for thickness in (0.0, -1.0, float("nan"), float("inf")):
            with self.subTest(thickness=thickness):
                with self.assertRaises(ValueError):
                    rasterize_vessels(
                        [vessel],
                        domain,
                        cfg,
                        effective_thickness_um=thickness,
                        continuous_geometry=continuous_geometry,
                    )

        for changes in (
            {"min_lumen_radius_cells": -1.0},
            {"min_lumen_radius_cells": float("nan")},
            {"min_resolved_diameter_cells": -1.0},
            {"min_resolved_diameter_cells": float("nan")},
        ):
            with self.subTest(changes=changes):
                with self.assertRaises(ValueError):
                    rasterize_vessels(
                        [vessel],
                        domain,
                        replace(cfg, **changes),
                        effective_thickness_um=1.0,
                        continuous_geometry=continuous_geometry,
                    )


if __name__ == "__main__":
    unittest.main()
