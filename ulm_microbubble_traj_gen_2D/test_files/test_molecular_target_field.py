from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from ulm_vascular_model_generator.utils.core.models import Vessel

from ulm_microbubble_traj_gen_2D.utils.core.types import (
    FlowField,
    GridDomain,
    RasterizedVessels,
)
from ulm_microbubble_traj_gen_2D.utils.flow.flow_boundaries import (
    enumerate_lumen_boundary_faces,
)
from ulm_microbubble_traj_gen_2D.utils.io.field_io import save_molecular_target_npz
from ulm_microbubble_traj_gen_2D.utils.molecular.molecular_target_field import (
    MolecularTargetField,
    build_molecular_target_field,
)
from ulm_microbubble_traj_gen_2D.utils.geometry.continuous_vessel_geometry import (
    build_continuous_vessel_geometry,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_boundary_support import (
    sample_continuous_boundary_support_to_grid,
)
from ulm_microbubble_traj_gen_2D.utils.particles.particle_hydrodynamic_fields import (
    ParticleHydrodynamicFields,
)
from ulm_microbubble_traj_gen_2D.test_files.hybrid_velocity_fixture import (
    rectangular_hybrid_velocity,
)


class MolecularTargetFieldTests(unittest.TestCase):
    def test_disabled_configuration_is_a_strict_no_op(self) -> None:
        domain, fields = _channel_case()

        target_field = build_molecular_target_field(domain, fields, None)

        self.assertFalse(target_field.enabled)
        self.assertEqual(target_field.region_mode, "disabled")
        self.assertEqual(target_field.target_density_molecules_per_m2, 0.0)
        self.assertFalse(np.any(target_field.target_wall_mask))
        self.assertFalse(np.any(target_field.target_density_field_molecules_per_m2))
        area = target_field.reaction_area_um2(
            np.asarray([[6.0, 1.5]]),
            np.asarray([[1.0, 0.0]]),
            np.asarray([0.5]),
        )
        np.testing.assert_array_equal(area, np.zeros(1))

    def test_mask_targets_only_solid_side_walls_and_excludes_open_caps(self) -> None:
        domain, fields = _channel_case()
        density = 2.7e15
        region_mask = np.zeros(domain.shape, dtype=bool)
        region_mask[5:8, :] = True
        target_field = _build_mask_target(domain, fields, region_mask, density)

        self.assertTrue(np.all(target_field.target_wall_mask <= target_field.solid_wall_mask))
        self.assertGreater(int(np.count_nonzero(target_field.target_wall_mask)), 0)
        self.assertFalse(np.any(target_field.solid_wall_mask[0, 2:7]))
        self.assertFalse(np.any(target_field.solid_wall_mask[-1, 2:7]))
        self.assertFalse(np.any(target_field.wall_axis == 0))
        np.testing.assert_array_equal(
            target_field.target_density_field_molecules_per_m2[target_field.target_wall_mask],
            np.full(np.count_nonzero(target_field.target_wall_mask), density),
        )
        self.assertFalse(
            np.any(target_field.target_density_field_molecules_per_m2[~target_field.target_wall_mask])
        )

    def test_reaction_area_is_continuous_at_mask_edge(self) -> None:
        domain, fields = _channel_case()
        region_mask = np.zeros(domain.shape, dtype=bool)
        region_mask[5:8, :] = True
        target_field = _build_mask_target(domain, fields, region_mask, 1.0e15)
        points = np.asarray(
            [
                [6.0, 1.5],
                [4.0, 1.5],
                [5.0, 1.5],
                [0.5, 4.0],
            ]
        )
        areas = target_field.reaction_area_um2(
            points,
            np.asarray([1.0, 0.0]),
            np.asarray([0.4, 0.25, 1.0, 1.0]),
        )

        self.assertAlmostEqual(float(areas[0]), np.pi * 0.4**2, places=13)
        self.assertEqual(float(areas[1]), 0.0)
        self.assertAlmostEqual(float(areas[2]), 0.5 * np.pi, places=13)
        self.assertEqual(float(areas[3]), 0.0)

        epsilon = 1.0e-6
        edge_limits = target_field.reaction_area_um2(
            np.asarray([[5.0 - epsilon, 1.5], [5.0 + epsilon, 1.5]]),
            np.asarray([1.0, 0.0]),
            np.asarray([1.0, 1.0]),
        )
        self.assertLess(abs(float(edge_limits[1] - edge_limits[0])), 5.0e-6)
        np.testing.assert_allclose(edge_limits, 0.5 * np.pi, atol=3.0e-6, rtol=0.0)

    def test_mask_npz_uses_physical_axes_and_boolean_region_values(self) -> None:
        domain, fields = _channel_case()
        region_mask = np.zeros(domain.shape, dtype=bool)
        region_mask[6:, :] = True

        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "target_region.npz"
            np.savez(
                path,
                x_coordinates_um=domain.x_coordinates_um,
                z_coordinates_um=domain.z_coordinates_um,
                target_region_mask=region_mask,
            )
            target_field = build_molecular_target_field(
                domain,
                fields,
                {
                    "enabled": True,
                    "region_mode": "mask_npz",
                    "target_density_molecules_per_m2": 3.0e14,
                    "mask_npz_path": path,
                },
            )

        self.assertEqual(target_field.region_mode, "mask_npz")
        self.assertTrue(np.all(target_field.target_wall_mask <= target_field.solid_wall_mask))
        targeted_indices = np.argwhere(target_field.target_wall_mask)
        self.assertGreater(targeted_indices.shape[0], 0)
        targeted_x = domain.x_coordinates_um[targeted_indices[:, 0]]
        self.assertTrue(np.all(targeted_x >= 6.0))
        np.testing.assert_array_equal(target_field.region_mask, region_mask)
        self.assertEqual(
            str(target_field.to_npz_payload()["source_mask_npz_path"]),
            str(path.resolve()),
        )

        epsilon = 1.0e-6
        mask_edge_areas = target_field.reaction_area_um2(
            np.asarray([[6.0 - epsilon, 1.5], [6.0, 1.5], [6.0 + epsilon, 1.5]]),
            np.asarray([1.0, 0.0]),
            np.ones(3),
        )
        self.assertAlmostEqual(float(mask_edge_areas[1]), 0.5 * np.pi, places=13)
        self.assertLess(float(np.ptp(mask_edge_areas)), 5.0e-6)

    def test_selector_output_schema_is_accepted_as_formal_wall_target_input(self) -> None:
        domain, fields = _channel_case()
        wall_only_target = np.asarray(fields.solid_site_mask, dtype=bool).copy()

        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "selected_molecular_target_mask.npz"
            np.savez_compressed(
                path,
                x_um=domain.x_coordinates_um,
                z_um=domain.z_coordinates_um,
                target_mask=wall_only_target,
                selection_mode=np.asarray("automatic_spatial_heterogeneity"),
                random_seed=np.asarray(42, dtype=np.int64),
            )
            target_field = build_molecular_target_field(
                domain,
                fields,
                {
                    "enabled": True,
                    "region_mode": "mask_npz",
                    "target_density_molecules_per_m2": 3.0e14,
                    "mask_npz_path": path,
                },
            )

        np.testing.assert_array_equal(target_field.region_mask, wall_only_target)
        np.testing.assert_array_equal(
            target_field.target_wall_mask,
            target_field.solid_wall_mask,
        )

    def test_wall_only_target_maps_both_channel_sides_symmetrically(self) -> None:
        domain, fields = _channel_case()
        target_field = _build_mask_target(
            domain,
            fields,
            np.asarray(fields.solid_site_mask, dtype=bool),
            1.0e15,
        )

        self.assertTrue(np.all(target_field.wall_target_positive))
        radius_um = 0.4
        areas = target_field.reaction_area_um2(
            np.asarray([[6.0, 1.5], [6.0, 6.5]]),
            np.asarray([[1.0, 0.0], [1.0, 0.0]]),
            np.asarray([radius_um, radius_um]),
        )
        np.testing.assert_allclose(
            areas,
            np.full(2, np.pi * radius_um**2),
            rtol=0.0,
            atol=1.0e-13,
        )

    def test_polyline_target_candidate_mask_rejects_only_distant_paths(self) -> None:
        domain, fields = _channel_case()
        region_mask = np.zeros(domain.shape, dtype=bool)
        region_mask[5:8, :] = True
        target_field = _build_mask_target(
            domain, fields, region_mask, 1.0e15
        )
        positive_faces = np.flatnonzero(target_field.wall_target_positive)
        self.assertGreater(positive_faces.size, 0)
        self.assertIsNotNone(target_field._target_positive_wall_tree)

        face_index = int(positive_faces[positive_faces.size // 2])
        target_centre = target_field.wall_coordinates_xz_um[face_index]
        outward_direction = -target_field.wall_normal_xz[face_index]
        points = np.vstack(
            (
                target_centre + 100.0 * outward_direction,
                target_centre + 101.0 * outward_direction,
                target_centre + 0.05 * outward_direction,
                target_centre + 10.05 * outward_direction,
                target_centre + 0.2 * outward_direction,
            )
        )

        candidate = target_field.polyline_target_candidate_mask(
            points,
            np.asarray([0, 2, 4, 5], dtype=np.int64),
            np.asarray([0.1, 0.1, 0.1], dtype=np.float64),
            0.1,
        )

        np.testing.assert_array_equal(
            candidate, np.asarray([False, True, True], dtype=bool)
        )

    def test_polyline_target_candidate_mask_without_index_keeps_every_lane(self) -> None:
        target_field = _two_face_corner_target_field(
            np.asarray([True, False], dtype=bool)
        )
        points = np.asarray(
            [[20.0, 20.0], [21.0, 20.0], [30.0, 30.0]],
            dtype=np.float64,
        )

        candidate = target_field.polyline_target_candidate_mask(
            points,
            np.asarray([0, 2, 3], dtype=np.int64),
            0.1,
            0.1,
        )

        np.testing.assert_array_equal(candidate, np.ones(2, dtype=bool))

    def test_polyline_target_candidate_mask_includes_corner_reaction_reach(self) -> None:
        radius_um = 1.0
        capture_um = 0.5
        target_field = _two_face_corner_target_field(
            np.asarray([True, False], dtype=bool),
            include_target_broad_phase_index=True,
        )
        target_centre = target_field.wall_coordinates_xz_um[0]
        old_straight_wall_bound = radius_um + capture_um + 0.5
        point = target_centre + np.asarray(
            [old_straight_wall_bound + 0.25, 0.0], dtype=np.float64
        )

        candidate = target_field.polyline_target_candidate_mask(
            point.reshape(1, 2),
            np.asarray([0, 1], dtype=np.int64),
            radius_um,
            capture_um,
        )

        np.testing.assert_array_equal(candidate, np.asarray([True], dtype=bool))

    def test_corner_projection_uses_one_exact_solid_face_without_normal_averaging(self) -> None:
        target_field = _two_face_corner_target_field(
            np.asarray([True, True], dtype=bool)
        )
        approximate = np.asarray([[0.2, 0.2]], dtype=np.float64)

        projected, face_indices, normals, tangents = (
            target_field.nearest_solid_face_frame_xz_um(approximate)
        )

        self.assertEqual(int(face_indices[0]), 0)
        np.testing.assert_allclose(projected[0], np.asarray([0.0, 0.2]), atol=1.0e-15)
        np.testing.assert_array_equal(normals[0], np.asarray([1.0, 0.0]))
        np.testing.assert_array_equal(tangents[0], np.asarray([0.0, 1.0]))
        area = target_field.reaction_area_um2(
            approximate,
            np.asarray([[1.0, 1.0]]),
            np.asarray([0.1]),
        )
        self.assertAlmostEqual(float(area[0]), np.pi * 0.1**2, places=14)

    def test_non_target_projected_face_does_not_create_reaction_area(self) -> None:
        target_field = _two_face_corner_target_field(
            np.asarray([False, True], dtype=bool)
        )
        area = target_field.reaction_area_um2(
            np.asarray([[0.2, 0.3]]),
            np.asarray([[0.0, 1.0]]),
            np.asarray([0.05]),
        )
        self.assertEqual(float(area[0]), 0.0)

    def test_equal_distance_corner_checks_each_face_without_averaging_normals(self) -> None:
        target_field = _two_face_corner_target_field(
            np.asarray([False, True], dtype=bool)
        )

        area = target_field.reaction_area_um2(
            np.asarray([[0.2, 0.2]]),
            np.asarray([[1.0, 1.0]]),
            np.asarray([0.1]),
        )

        self.assertAlmostEqual(float(area[0]), np.pi * 0.1**2, places=14)

    def test_mask_npz_honours_configured_array_keys(self) -> None:
        domain, fields = _channel_case()
        region_mask = np.zeros(domain.shape, dtype=bool)
        region_mask[5:8, :] = True

        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "custom_keys.npz"
            np.savez(
                path,
                physical_x=domain.x_coordinates_um,
                physical_z=domain.z_coordinates_um,
                lesion_wall_region=region_mask,
            )
            target_field = build_molecular_target_field(
                domain,
                fields,
                {
                    "enabled": True,
                    "region_mode": "mask_npz",
                    "target_density_molecules_per_m2": 4.0e14,
                    "mask_npz_path": path,
                    "x_coordinates_key": "physical_x",
                    "z_coordinates_key": "physical_z",
                    "mask_array_key": "lesion_wall_region",
                },
            )

        np.testing.assert_array_equal(target_field.region_mask, region_mask)
        self.assertGreater(int(np.count_nonzero(target_field.target_wall_mask)), 0)

    def test_zero_density_enabled_mask_supports_contact_only_pilot(self) -> None:
        domain, fields = _channel_case()
        region_mask = np.zeros(domain.shape, dtype=bool)
        region_mask[5:8, :] = True
        target_field = _build_mask_target(domain, fields, region_mask, 0.0)

        self.assertTrue(target_field.enabled)
        self.assertGreater(int(np.count_nonzero(target_field.target_wall_mask)), 0)
        self.assertFalse(np.any(target_field.target_density_field_molecules_per_m2))
        area = target_field.reaction_area_um2(
            np.asarray([[6.0, 1.5]]),
            np.asarray([[1.0, 0.0]]),
            0.4,
        )
        self.assertAlmostEqual(float(area[0]), np.pi * 0.4**2, places=13)

    def test_mask_npz_rejects_non_boolean_region_data(self) -> None:
        domain, fields = _channel_case()

        with tempfile.TemporaryDirectory() as temporary_directory:
            path = Path(temporary_directory) / "invalid_target_region.npz"
            np.savez(
                path,
                x_coordinates_um=domain.x_coordinates_um,
                z_coordinates_um=domain.z_coordinates_um,
                target_region_mask=np.ones(domain.shape, dtype=np.uint8),
            )
            with self.assertRaisesRegex(ValueError, "Boolean"):
                build_molecular_target_field(
                    domain,
                    fields,
                    {
                        "enabled": True,
                        "region_mode": "mask_npz",
                        "target_density_molecules_per_m2": 1.0e15,
                        "mask_npz_path": path,
                    },
                )

    def test_npz_payload_contains_physical_coordinates_and_wall_fields(self) -> None:
        domain, fields = _channel_case()
        region_mask = np.zeros(domain.shape, dtype=bool)
        region_mask[5:8, :] = True
        target_field = _build_mask_target(domain, fields, region_mask, 9.0e14)

        payload = target_field.to_npz_payload(prefix="molecular_target_")

        expected_keys = {
            "molecular_target_x_coordinates_um",
            "molecular_target_z_coordinates_um",
            "molecular_target_solid_wall_mask",
            "molecular_target_target_wall_mask",
            "molecular_target_open_boundary_mask",
            "molecular_target_target_density_molecules_per_m2",
            "molecular_target_region_mode",
            "molecular_target_wall_target_positive",
        }
        self.assertTrue(expected_keys <= payload.keys())
        np.testing.assert_array_equal(payload["molecular_target_x_coordinates_um"], domain.x_coordinates_um)
        np.testing.assert_array_equal(payload["molecular_target_target_wall_mask"], target_field.target_wall_mask)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "molecular_target_field.npz"
            save_molecular_target_npz(path, target_field)
            with np.load(path, allow_pickle=False) as saved:
                self.assertIn("open_boundary_mask", saved.files)
                np.testing.assert_array_equal(
                    saved["target_wall_mask"],
                    target_field.target_wall_mask,
                )


def _build_mask_target(
    domain: GridDomain,
    fields: ParticleHydrodynamicFields,
    region_mask: np.ndarray,
    density: float,
):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "selected_target.npz"
        np.savez(
            path,
            x_um=domain.x_coordinates_um,
            z_um=domain.z_coordinates_um,
            target_mask=np.asarray(region_mask, dtype=bool),
        )
        return build_molecular_target_field(
            domain,
            fields,
            {
                "enabled": True,
                "region_mode": "mask_npz",
                "target_density_molecules_per_m2": density,
                "mask_npz_path": path,
            },
        )


def _two_face_corner_target_field(
    wall_target_positive: np.ndarray,
    *,
    include_target_broad_phase_index: bool = False,
) -> MolecularTargetField:
    wall_coordinates = np.asarray([[0.0, 0.5], [0.5, 0.0]], dtype=np.float64)
    wall_normals = np.asarray([[1.0, 0.0], [0.0, 1.0]], dtype=np.float64)
    wall_axis = np.asarray([0, 1], dtype=np.int8)
    wall_length = np.ones(2, dtype=np.float64)
    wall_start = np.asarray([[0.0, 0.0], [0.0, 0.0]], dtype=np.float64)
    wall_end = np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=np.float64)
    wall_tangent = np.asarray([[0.0, 1.0], [1.0, 0.0]], dtype=np.float64)
    shape = (2, 2)
    zeros = np.zeros(shape, dtype=bool)
    target_flags = np.asarray(wall_target_positive, dtype=bool)
    positive_faces = np.flatnonzero(target_flags)
    target_tree = (
        cKDTree(wall_coordinates[positive_faces], copy_data=True)
        if include_target_broad_phase_index and positive_faces.size
        else None
    )
    return MolecularTargetField(
        enabled=True,
        region_mode="mask_npz",
        target_density_molecules_per_m2=1.0e15,
        spacing_um=1.0,
        x_coordinates_um=np.arange(2, dtype=np.float64),
        z_coordinates_um=np.arange(2, dtype=np.float64),
        solid_wall_mask=zeros.copy(),
        target_wall_mask=zeros.copy(),
        target_density_field_molecules_per_m2=np.zeros(shape, dtype=np.float64),
        wall_coordinates_xz_um=wall_coordinates,
        wall_normal_xz=wall_normals,
        wall_axis=wall_axis,
        wall_length_um=wall_length,
        wall_start_xz_um=wall_start,
        wall_end_xz_um=wall_end,
        wall_tangent_xz=wall_tangent,
        wall_ring_index=np.arange(2, dtype=np.int64),
        wall_arclength_start_um=np.full(2, np.nan, dtype=np.float64),
        wall_arclength_end_um=np.full(2, np.nan, dtype=np.float64),
        boundary_ring_length_um=0.0,
        wall_target_positive=target_flags,
        open_boundary_mask=zeros.copy(),
        _maximum_wall_half_length_um=0.5,
        _boundary_geometry=_TwoFaceCornerGeometry(),
        _wall_tree=cKDTree(wall_coordinates),
        _target_positive_wall_tree=target_tree,
        _maximum_target_face_half_length_um=(
            0.5 if target_tree is not None else float("nan")
        ),
    )


class _TwoFaceCornerGeometry:
    """Small exact-face stand-in used only by the molecular corner unit tests."""

    def nearest_solid_face_projection_xz_um(
        self,
        positions_xz_um: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, tuple[np.ndarray, ...]]:
        points = np.asarray(positions_xz_um, dtype=np.float64).reshape(-1, 2)
        projected = np.empty_like(points)
        primary = np.empty(points.shape[0], dtype=np.int64)
        tied_faces: list[np.ndarray] = []
        for lane, point in enumerate(points):
            candidates = np.asarray(
                [
                    [0.0, np.clip(point[1], 0.0, 1.0)],
                    [np.clip(point[0], 0.0, 1.0), 0.0],
                ],
                dtype=np.float64,
            )
            distances = np.linalg.norm(candidates - point[None, :], axis=1)
            minimum = float(np.min(distances))
            tied = np.flatnonzero(
                np.abs(distances - minimum)
                <= 256.0 * np.finfo(np.float64).eps
            ).astype(np.int64)
            primary[lane] = int(tied[0])
            projected[lane] = candidates[int(tied[0])]
            tied_faces.append(tied)
        return projected, primary, tuple(tied_faces)


def _channel_case() -> tuple[GridDomain, ParticleHydrodynamicFields]:
    nx, nz = 13, 9
    shape = (nx, nz)
    spacing_um = 1.0
    x_um = np.arange(nx, dtype=np.float64) * spacing_um
    z_um = np.arange(nz, dtype=np.float64) * spacing_um
    domain = GridDomain(
        origin_um=np.asarray([0.0, 0.0, 0.0]),
        spacing_um=spacing_um,
        shape=shape,
        fixed_y_um=0.0,
        x_coordinates_um=x_um,
        z_coordinates_um=z_um,
    )

    lumen = np.zeros(shape, dtype=bool)
    lumen[:, 2:7] = True
    catalog = enumerate_lumen_boundary_faces(domain, lumen)
    inlet_rows = (catalog.axis == 0) & (catalog.normal_xz[:, 0] < 0.0)
    outlet_rows = (catalog.axis == 0) & (catalog.normal_xz[:, 0] > 0.0)
    selected = inlet_rows | outlet_rows
    open_kind = np.where(outlet_rows[selected], 1, -1).astype(np.int8)

    zeros = np.zeros(shape, dtype=np.float32)
    raster = RasterizedVessels(
        lumen_mask=lumen,
        wall_mask=lumen.copy(),
        vessel_id=np.zeros(shape, dtype=np.int32),
        radius_um=np.full(shape, 2.5, dtype=np.float32),
        flow_rate_um3_s=zeros.copy(),
        q2d_flow_um2_s=zeros.copy(),
        viscosity_mpas=np.where(lumen, 3.5, 0.0).astype(np.float32),
        direction_xz=np.zeros((*shape, 2), dtype=np.float32),
        distance_to_centerline_um=zeros.copy(),
        distance_to_wall_um=zeros.copy(),
        wall_normal_xz=np.zeros((*shape, 2), dtype=np.float32),
        lumen_fraction=lumen.astype(np.float32),
    )
    flow = FlowField(
        velocity_xz_um_s=np.zeros((*shape, 2), dtype=np.float32),
        speed_um_s=zeros.copy(),
        wall_shear_stress_pa=zeros.copy(),
        hybrid_velocity=rectangular_hybrid_velocity(domain),
        open_face_cell_ij=catalog.cell_ij[selected],
        open_face_index_ij=catalog.face_index_ij[selected],
        open_face_axis=catalog.axis[selected],
        open_face_normal_xz=catalog.normal_xz[selected],
        open_face_center_xz_um=catalog.center_xz_um[selected],
        open_face_length_um=catalog.length_um[selected],
        open_face_label=np.ones(int(np.count_nonzero(selected)), dtype=np.int32),
        open_face_kind=open_kind,
        open_section_point_xz_um=np.asarray([[-0.5, 4.0], [12.5, 4.0]], dtype=np.float64),
        open_section_outward_normal_xz=np.asarray([[-1.0, 0.0], [1.0, 0.0]], dtype=np.float64),
        open_section_tangent_xz=np.asarray([[0.0, -1.0], [0.0, 1.0]], dtype=np.float64),
        open_section_half_width_um=np.asarray([2.5, 2.5], dtype=np.float64),
        open_section_label=np.asarray([1, 1], dtype=np.int32),
        open_section_kind=np.asarray([-1, 1], dtype=np.int8),
    )
    vessel = Vessel(vid=0, parent_id=-1, children=[])
    vessel.x_p = np.asarray([0.0, 0.0, 4.0])
    vessel.x_d = np.asarray([12.0, 0.0, 4.0])
    vessel.radius = 2.5
    boundary_geometry = build_continuous_vessel_geometry([vessel], domain)
    solid_sites, open_boundary = sample_continuous_boundary_support_to_grid(
        domain, boundary_geometry
    )
    fields = ParticleHydrodynamicFields(
        velocity_gradient_s_inv=np.zeros((*shape, 2, 2), dtype=np.float32),
        dynamic_viscosity_pa_s=np.full(shape, 3.5e-3, dtype=np.float32),
        solid_site_mask=solid_sites,
        open_boundary_mask=open_boundary,
        boundary_geometry=boundary_geometry,
        hybrid_velocity=rectangular_hybrid_velocity(domain),
    )
    return domain, fields


if __name__ == "__main__":
    unittest.main()
