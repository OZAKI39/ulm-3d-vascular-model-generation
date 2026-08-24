from __future__ import annotations

from pathlib import Path

import numpy as np
import trimesh

from utils.cfd_lumen.mask_assisted_refinement import INFLUENCE_CONFIGS, refine_surface_with_mask
from utils.cfd_lumen.mask_comparison import binary_overlap_metrics
from utils.cfd_lumen.mask_reference import RadiusCalibration, SWCReference, reconstruct_mask_from_swc
from utils.cfd_lumen.mask_surface import mask_surfaces, signed_distance_field


def test_binary_overlap_identity() -> None:
    mask = np.zeros((8, 9, 10), dtype=bool)
    mask[2:6, 3:7, 4:8] = True
    metrics = binary_overlap_metrics(mask, mask)
    assert metrics["dice"] == 1.0
    assert metrics["iou"] == 1.0
    assert metrics["volume_difference_um3"] == 0.0


def test_swc_mask_reconstruction_has_radius_aware_foreground() -> None:
    swc = SWCReference(
        path=Path("synthetic.swc"),
        node_ids=np.asarray((1, 2), dtype=np.int64),
        points_voxel_xyz=np.asarray(((8.0, 8.0, 4.0), (16.0, 8.0, 4.0))),
        radius_um=np.asarray((2.0, 4.0)),
        parent_ids=np.asarray((-1, 1), dtype=np.int64),
    )
    mask, metadata = reconstruct_mask_from_swc(
        swc,
        (12, 24, 24),
        calibration=RadiusCalibration(radius_scale=0.5, voxel_floor=1.0),
    )
    assert mask.shape == (12, 24, 24)
    assert mask[4, 8, 8] == 255
    assert np.count_nonzero(mask[:, :, 14:] > 103) > np.count_nonzero(mask[:, :, :4] > 103)
    assert metadata["experimental_only"] is True


def test_physical_sdf_mask_surfaces_are_closed() -> None:
    z, y, x = np.indices((13, 21, 21))
    mask = ((x - 10) ** 2 + (y - 10) ** 2 + ((z - 6) * 2.0) ** 2) <= 5.0**2
    padded = np.pad(mask, 1)
    raw, smoothed, metadata = mask_surfaces(
        padded,
        origin_xyz_um=(-1.0, -1.0, -2.0),
        smoothing_um=0.5,
    )
    assert raw.is_watertight
    assert smoothed.is_watertight
    assert raw.volume > 0
    assert metadata["smoothed_field_sigma_um"] == 0.5


def test_refinement_preserves_faces_and_fixes_explicit_vertices() -> None:
    mesh = trimesh.creation.icosphere(subdivisions=2, radius=3.0)
    face_region = np.ones(len(mesh.faces), dtype=np.uint8)
    explicit_faces = np.flatnonzero(mesh.triangles_center[:, 0] > 1.0)
    face_region[explicit_faces] = 3
    explicit_vertices = np.unique(mesh.faces[explicit_faces])
    original = np.asarray(mesh.vertices).copy()
    z, y, x = np.indices((13, 25, 25))
    mask = ((x - 12) ** 2 + (y - 12) ** 2 + ((z - 6) * 2.0) ** 2) <= 3.4**2
    sdf = signed_distance_field(mask)
    refined, report = refine_surface_with_mask(
        mesh,
        face_region,
        sdf,
        origin_xyz_um=(-12.0, -12.0, -12.0),
        junction_center_xyz_um=np.zeros(3),
        junction_radius_um=1.0,
        influence_outer_radius_um=4.0,
        config=INFLUENCE_CONFIGS[0],
    )
    assert np.array_equal(refined.faces, mesh.faces)
    assert np.array_equal(np.asarray(refined.vertices)[explicit_vertices], original[explicit_vertices])
    assert report["source_junction_connectivity_unchanged"] is True
    assert report["explicit_branch_vertex_count_moved"] == 0
