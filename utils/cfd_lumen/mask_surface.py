"""Physical-space Mask surfaces and local junction artifact helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh
from scipy import ndimage
from skimage.measure import marching_cubes

from .mask_reference import SPACING_XYZ_UM, SPACING_ZYX_UM, SWCReference, swc_edges


@dataclass(frozen=True, slots=True)
class VolumeCrop:
    array: np.ndarray
    origin_index_zyx: tuple[int, int, int]
    origin_xyz_um: tuple[float, float, float]
    bounds_xyz_um: tuple[float, float, float, float, float, float]


def signed_distance_field(mask: np.ndarray) -> np.ndarray:
    """Return a physical SDF: negative inside and positive outside."""

    binary = np.asarray(mask, dtype=bool)
    outside = ndimage.distance_transform_edt(~binary, sampling=SPACING_ZYX_UM)
    inside = ndimage.distance_transform_edt(binary, sampling=SPACING_ZYX_UM)
    return (outside - inside).astype(np.float32)


def _mesh_from_field(
    field: np.ndarray,
    *,
    origin_xyz_um: tuple[float, float, float],
) -> trimesh.Trimesh:
    vertices_zyx, faces, _, _ = marching_cubes(
        np.asarray(field, dtype=np.float32),
        level=0.0,
        spacing=SPACING_ZYX_UM,
        allow_degenerate=False,
    )
    vertices_xyz = vertices_zyx[:, ::-1] + np.asarray(origin_xyz_um, dtype=float)
    mesh = trimesh.Trimesh(vertices=vertices_xyz, faces=faces, process=False)
    mesh.remove_unreferenced_vertices()
    trimesh.repair.fix_normals(mesh, multibody=True)
    if mesh.volume < 0:
        mesh.invert()
    return mesh


def mask_surfaces(
    binary_mask: np.ndarray,
    *,
    origin_xyz_um: tuple[float, float, float],
    smoothing_um: float = 0.5,
) -> tuple[trimesh.Trimesh, trimesh.Trimesh, dict[str, Any]]:
    """Build raw and minimally SDF-smoothed Mask-only surfaces."""

    field = signed_distance_field(binary_mask)
    raw = _mesh_from_field(field, origin_xyz_um=origin_xyz_um)
    sigma_zyx = tuple(float(smoothing_um / spacing) for spacing in SPACING_ZYX_UM)
    smoothed_field = ndimage.gaussian_filter(field, sigma=sigma_zyx, mode="nearest")
    smoothed = _mesh_from_field(smoothed_field, origin_xyz_um=origin_xyz_um)
    metadata = {
        "method": "physical_signed_distance_zero_isosurface",
        "spacing_xyz_um": list(SPACING_XYZ_UM),
        "raw_field_smoothed": False,
        "smoothed_field_sigma_um": float(smoothing_um),
        "smoothed_field_sigma_zyx_voxel": list(sigma_zyx),
        "over_smoothing_forbidden": True,
        "raw_triangle_count": int(len(raw.faces)),
        "smoothed_triangle_count": int(len(smoothed.faces)),
    }
    return raw, smoothed, metadata


def save_mesh_vtp(
    mesh: trimesh.Trimesh,
    path: Path,
    *,
    cell_data: dict[str, np.ndarray] | None = None,
) -> Path:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    faces = np.column_stack(
        (np.full(len(mesh.faces), 3, dtype=np.int64), np.asarray(mesh.faces, dtype=np.int64))
    ).ravel()
    polydata = pv.PolyData(np.asarray(mesh.vertices, dtype=float), faces)
    for name, values in (cell_data or {}).items():
        if len(values) == polydata.n_cells:
            polydata.cell_data[name] = np.asarray(values)
    polydata.save(target)
    return target


def load_mesh_vtp(path: Path) -> tuple[trimesh.Trimesh, dict[str, np.ndarray]]:
    polydata = pv.read(path).extract_surface(algorithm="dataset_surface").triangulate()
    faces = np.asarray(polydata.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    mesh = trimesh.Trimesh(
        vertices=np.asarray(polydata.points, dtype=float),
        faces=faces,
        process=False,
    )
    cell_data = {name: np.asarray(values) for name, values in polydata.cell_data.items()}
    return mesh, cell_data


def physical_bounds_to_indices(
    bounds_xyz_um: tuple[float, float, float, float, float, float],
    shape_zyx: tuple[int, int, int],
    *,
    padding_voxels: int = 1,
) -> tuple[slice, slice, slice]:
    xmin, xmax, ymin, ymax, zmin, zmax = bounds_xyz_um
    lower_xyz = np.floor(
        np.asarray((xmin, ymin, zmin), dtype=float) / np.asarray(SPACING_XYZ_UM)
    ).astype(int) - int(padding_voxels)
    upper_xyz = np.ceil(
        np.asarray((xmax, ymax, zmax), dtype=float) / np.asarray(SPACING_XYZ_UM)
    ).astype(int) + int(padding_voxels) + 1
    shape_xyz = np.asarray(shape_zyx[::-1], dtype=int)
    lower_xyz = np.maximum(lower_xyz, 0)
    upper_xyz = np.minimum(upper_xyz, shape_xyz)
    return (
        slice(int(lower_xyz[2]), int(upper_xyz[2])),
        slice(int(lower_xyz[1]), int(upper_xyz[1])),
        slice(int(lower_xyz[0]), int(upper_xyz[0])),
    )


def crop_volume_physical(
    array: np.ndarray,
    bounds_xyz_um: tuple[float, float, float, float, float, float],
    *,
    padding_voxels: int = 2,
) -> VolumeCrop:
    slices = physical_bounds_to_indices(
        bounds_xyz_um,
        tuple(map(int, array.shape)),
        padding_voxels=padding_voxels,
    )
    origin_zyx = tuple(int(part.start or 0) for part in slices)
    origin_xyz = tuple(
        np.asarray(origin_zyx[::-1], dtype=float) * np.asarray(SPACING_XYZ_UM)
    )
    cropped = np.asarray(array[slices]).copy()
    return VolumeCrop(cropped, origin_zyx, origin_xyz, bounds_xyz_um)


def select_junction_component(binary: np.ndarray, center_local_zyx: np.ndarray) -> tuple[np.ndarray, dict[str, int]]:
    """Keep the component nearest the SWC junction; never infer topology from it."""

    labels, count = ndimage.label(np.asarray(binary, dtype=bool), structure=np.ones((3, 3, 3)))
    if count == 0:
        raise ValueError("Local junction Mask has no foreground component")
    center = np.rint(center_local_zyx).astype(int)
    center = np.clip(center, 0, np.asarray(binary.shape) - 1)
    selected = int(labels[tuple(center)])
    if selected == 0:
        foreground = np.argwhere(labels > 0)
        scale = np.asarray(SPACING_ZYX_UM, dtype=float)
        nearest = int(np.argmin(np.linalg.norm((foreground - center[None, :]) * scale, axis=1)))
        selected = int(labels[tuple(foreground[nearest])])
    output = labels == selected
    return output, {
        "local_mask_component_count_26": int(count),
        "selected_component_label": selected,
        "selected_component_voxel_count": int(np.count_nonzero(output)),
        "discarded_component_count_for_mask_only_surface": int(count - 1),
        "topology_source": "SWC (component choice is only for local volumetric visualization)",
    }


def local_surface(
    mesh: trimesh.Trimesh,
    bounds_xyz_um: tuple[float, float, float, float, float, float],
) -> tuple[trimesh.Trimesh, np.ndarray]:
    xmin, xmax, ymin, ymax, zmin, zmax = bounds_xyz_um
    centers = np.asarray(mesh.triangles_center, dtype=float)
    selected = (
        (centers[:, 0] >= xmin)
        & (centers[:, 0] <= xmax)
        & (centers[:, 1] >= ymin)
        & (centers[:, 1] <= ymax)
        & (centers[:, 2] >= zmin)
        & (centers[:, 2] <= zmax)
    )
    face_ids = np.flatnonzero(selected)
    if not len(face_ids):
        raise ValueError("No surface triangles intersect the local junction bounds")
    subset = mesh.submesh([face_ids], append=True, repair=False)
    return subset, face_ids


def save_local_swc_vtp(
    swc: SWCReference,
    bounds_xyz_um: tuple[float, float, float, float, float, float],
    path: Path,
) -> Path:
    """Save source SWC edges touching a physical junction box."""

    points = swc.points_um
    xmin, xmax, ymin, ymax, zmin, zmax = bounds_xyz_um
    inside = (
        (points[:, 0] >= xmin)
        & (points[:, 0] <= xmax)
        & (points[:, 1] >= ymin)
        & (points[:, 1] <= ymax)
        & (points[:, 2] >= zmin)
        & (points[:, 2] <= zmax)
    )
    edge_list = [(first, second) for first, second in swc_edges(swc) if inside[first] or inside[second]]
    if not edge_list:
        raise ValueError("No SWC edges touch the local junction box")
    used = sorted({index for edge in edge_list for index in edge})
    local_by_global = {global_index: local_index for local_index, global_index in enumerate(used)}
    lines: list[int] = []
    for first, second in edge_list:
        lines.extend((2, local_by_global[first], local_by_global[second]))
    polydata = pv.PolyData(points[np.asarray(used)], lines=np.asarray(lines, dtype=np.int64))
    polydata.point_data["radius_um"] = swc.radius_um[np.asarray(used)]
    polydata.point_data["source_swc_node_id"] = swc.node_ids[np.asarray(used)]
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    polydata.save(target)
    return target


def sample_sdf_and_gradient(
    field: np.ndarray,
    points_xyz_um: np.ndarray,
    *,
    origin_xyz_um: tuple[float, float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """Trilinearly sample a physical SDF and its physical x-y-z gradient."""

    points = np.asarray(points_xyz_um, dtype=float)
    origin = np.asarray(origin_xyz_um, dtype=float)
    coordinates_xyz = (points - origin[None, :]) / np.asarray(SPACING_XYZ_UM)
    coordinates_zyx = coordinates_xyz[:, ::-1].T
    values = ndimage.map_coordinates(
        field, coordinates_zyx, order=1, mode="nearest", prefilter=False
    )
    gradient_zyx = np.gradient(field.astype(float), *SPACING_ZYX_UM)
    gradient = np.column_stack(
        [
            ndimage.map_coordinates(
                component, coordinates_zyx, order=1, mode="nearest", prefilter=False
            )
            for component in gradient_zyx[::-1]
        ]
    )
    norm = np.linalg.norm(gradient, axis=1)
    gradient = np.divide(
        gradient,
        norm[:, None],
        out=np.zeros_like(gradient),
        where=norm[:, None] > 1.0e-12,
    )
    return values, gradient
