"""Unified volumetric, boundary, surface, and anisotropy metrics for v5-1."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import trimesh
from scipy import ndimage
from scipy.spatial import cKDTree

from .mask_reference import SPACING_XYZ_UM, SPACING_ZYX_UM
from .mask_surface import sample_sdf_and_gradient
from .mesh_defects import diagnose_mesh_defects, triangle_quality
from .surface_qc import _section_polygon


@dataclass(frozen=True, slots=True)
class SectionLocation:
    name: str
    center_xyz_um: tuple[float, float, float]
    tangent_xyz: tuple[float, float, float]
    target_radius_um: float
    branch_role: str


def binary_overlap_metrics(reference: np.ndarray, candidate: np.ndarray) -> dict[str, float | int]:
    truth = np.asarray(reference, dtype=bool)
    predicted = np.asarray(candidate, dtype=bool)
    if truth.shape != predicted.shape:
        raise ValueError(f"Mask shapes differ: {truth.shape} != {predicted.shape}")
    true_positive = int(np.count_nonzero(truth & predicted))
    false_positive = int(np.count_nonzero(~truth & predicted))
    false_negative = int(np.count_nonzero(truth & ~predicted))
    union = true_positive + false_positive + false_negative
    denominator = 2 * true_positive + false_positive + false_negative
    voxel_volume = float(np.prod(SPACING_XYZ_UM))
    truth_count = int(np.count_nonzero(truth))
    predicted_count = int(np.count_nonzero(predicted))
    return {
        "dice": float(2 * true_positive / denominator) if denominator else 1.0,
        "iou": float(true_positive / union) if union else 1.0,
        "precision": (
            float(true_positive / (true_positive + false_positive))
            if true_positive + false_positive
            else 1.0
        ),
        "recall": (
            float(true_positive / (true_positive + false_negative))
            if true_positive + false_negative
            else 1.0
        ),
        "true_positive_voxels": true_positive,
        "false_positive_voxels": false_positive,
        "false_negative_voxels": false_negative,
        "reference_voxels": truth_count,
        "candidate_voxels": predicted_count,
        "reference_volume_um3": truth_count * voxel_volume,
        "candidate_volume_um3": predicted_count * voxel_volume,
        "volume_difference_um3": (predicted_count - truth_count) * voxel_volume,
        "volume_relative_difference": (
            float((predicted_count - truth_count) / truth_count) if truth_count else 0.0
        ),
    }


def _boundary(mask: np.ndarray) -> np.ndarray:
    binary = np.asarray(mask, dtype=bool)
    return binary ^ ndimage.binary_erosion(binary, structure=np.ones((3, 3, 3)))


def boundary_distance_metrics(reference: np.ndarray, candidate: np.ndarray) -> dict[str, float | int | None]:
    truth_boundary = _boundary(reference)
    candidate_boundary = _boundary(candidate)
    if not np.any(truth_boundary) or not np.any(candidate_boundary):
        return {
            "surface_distance_mean_um": None,
            "surface_distance_median_um": None,
            "surface_distance_p95_um": None,
            "hausdorff_distance_um": None,
            "hausdorff_95_um": None,
            "reference_boundary_voxels": int(np.count_nonzero(truth_boundary)),
            "candidate_boundary_voxels": int(np.count_nonzero(candidate_boundary)),
        }
    distance_to_truth = ndimage.distance_transform_edt(~truth_boundary, sampling=SPACING_ZYX_UM)
    distance_to_candidate = ndimage.distance_transform_edt(
        ~candidate_boundary, sampling=SPACING_ZYX_UM
    )
    directed_candidate = distance_to_truth[candidate_boundary]
    directed_truth = distance_to_candidate[truth_boundary]
    symmetric = np.concatenate((directed_candidate, directed_truth))
    return {
        "surface_distance_mean_um": float(np.mean(symmetric)),
        "surface_distance_median_um": float(np.median(symmetric)),
        "surface_distance_p95_um": float(np.percentile(symmetric, 95)),
        "hausdorff_distance_um": float(np.max(symmetric)),
        "hausdorff_95_um": float(
            max(np.percentile(directed_candidate, 95), np.percentile(directed_truth, 95))
        ),
        "reference_to_candidate_mean_um": float(np.mean(directed_truth)),
        "candidate_to_reference_mean_um": float(np.mean(directed_candidate)),
        "reference_boundary_voxels": int(np.count_nonzero(truth_boundary)),
        "candidate_boundary_voxels": int(np.count_nonzero(candidate_boundary)),
    }


def compare_binary_masks(reference: np.ndarray, candidate: np.ndarray) -> dict[str, Any]:
    return {**binary_overlap_metrics(reference, candidate), **boundary_distance_metrics(reference, candidate)}


def _swc_degree(node_ids: np.ndarray, parent_ids: np.ndarray) -> np.ndarray:
    degree = np.zeros(len(node_ids), dtype=np.int32)
    index = {int(node_id): position for position, node_id in enumerate(node_ids)}
    for child, parent_id in enumerate(parent_ids):
        parent = index.get(int(parent_id))
        if parent is not None:
            degree[child] += 1
            degree[parent] += 1
    return degree


def regional_mask_metrics(
    reference: np.ndarray,
    candidate: np.ndarray,
    *,
    dense_points_voxel_xyz: np.ndarray,
    dense_radius_um: np.ndarray,
    dense_owner_indices: np.ndarray,
    node_ids: np.ndarray,
    parent_ids: np.ndarray,
) -> list[dict[str, Any]]:
    """Partition union voxels into junction/terminal/large/capillary/branch regions."""

    union_indices = np.argwhere(np.asarray(reference, dtype=bool) | np.asarray(candidate, dtype=bool))
    if not len(union_indices):
        return []
    query_xyz = union_indices[:, ::-1].astype(float)
    _, nearest = cKDTree(np.asarray(dense_points_voxel_xyz, dtype=float)).query(
        query_xyz, workers=-1
    )
    nearest = np.asarray(nearest, dtype=np.int64)
    owner = np.asarray(dense_owner_indices, dtype=np.int64)[nearest]
    radius = np.asarray(dense_radius_um, dtype=float)[nearest]
    degree = _swc_degree(np.asarray(node_ids), np.asarray(parent_ids))
    owner_degree = degree[owner]
    labels = np.full(len(union_indices), "ordinary_branch", dtype=object)
    labels[radius <= 2.0] = "capillary_or_small_vessel"
    labels[radius >= 4.0] = "large_vessel"
    labels[owner_degree == 1] = "terminal"
    labels[owner_degree >= 3] = "junction"
    truth_values = np.asarray(reference, dtype=bool)[tuple(union_indices.T)]
    candidate_values = np.asarray(candidate, dtype=bool)[tuple(union_indices.T)]
    rows: list[dict[str, Any]] = []
    for name in (
        "ordinary_branch",
        "junction",
        "terminal",
        "large_vessel",
        "capillary_or_small_vessel",
    ):
        selected = labels == name
        if not np.any(selected):
            continue
        rows.append({"region": name, **binary_overlap_metrics(truth_values[selected], candidate_values[selected])})
    return rows


def mesh_from_voxel_centers(
    mesh: trimesh.Trimesh,
    shape_zyx: tuple[int, int, int],
    *,
    origin_xyz_um: tuple[float, float, float],
) -> np.ndarray:
    """Voxelize a watertight mesh at dataset voxel centers in a local crop."""

    z, y, x = np.indices(shape_zyx, dtype=float)
    points = np.column_stack(
        (
            x.ravel() * SPACING_XYZ_UM[0] + origin_xyz_um[0],
            y.ravel() * SPACING_XYZ_UM[1] + origin_xyz_um[1],
            z.ravel() * SPACING_XYZ_UM[2] + origin_xyz_um[2],
        )
    )
    surface = _polydata(mesh)
    point_cloud = pv.PolyData(points)
    classified = point_cloud.select_interior_points(
        surface,
        method="cell_locator",
        locator_tolerance=0.0,
        check_surface=False,
    )
    inside = np.asarray(classified.point_data["selected_points"], dtype=bool)
    return inside.reshape(shape_zyx)


def surface_to_mask_metrics(
    mesh: trimesh.Trimesh,
    sdf: np.ndarray,
    *,
    origin_xyz_um: tuple[float, float, float],
    bounds_xyz_um: tuple[float, float, float, float, float, float] | None = None,
) -> dict[str, float | int | None]:
    points = np.asarray(mesh.vertices, dtype=float)
    if bounds_xyz_um is not None:
        xmin, xmax, ymin, ymax, zmin, zmax = bounds_xyz_um
        selected = (
            (points[:, 0] >= xmin)
            & (points[:, 0] <= xmax)
            & (points[:, 1] >= ymin)
            & (points[:, 1] <= ymax)
            & (points[:, 2] >= zmin)
            & (points[:, 2] <= zmax)
        )
        points = points[selected]
    if not len(points):
        return {
            "surface_sample_count": 0,
            "surface_to_mask_mean_um": None,
            "surface_to_mask_median_um": None,
            "surface_to_mask_p95_um": None,
            "surface_to_mask_max_um": None,
            "mask_fitting_energy_um2": None,
        }
    values, _ = sample_sdf_and_gradient(sdf, points, origin_xyz_um=origin_xyz_um)
    absolute = np.abs(values)
    return {
        "surface_sample_count": int(len(points)),
        "surface_to_mask_mean_um": float(np.mean(absolute)),
        "surface_to_mask_median_um": float(np.median(absolute)),
        "surface_to_mask_p95_um": float(np.percentile(absolute, 95)),
        "surface_to_mask_max_um": float(np.max(absolute)),
        "mask_fitting_energy_um2": float(np.mean(values**2)),
    }


def local_normal_roughness(
    mesh: trimesh.Trimesh,
    *,
    center_xyz_um: np.ndarray,
    influence_radius_um: float,
) -> dict[str, float | int | None]:
    center = np.asarray(center_xyz_um, dtype=float)
    angles = np.degrees(np.asarray(mesh.face_adjacency_angles, dtype=float))
    midpoints = np.asarray(mesh.vertices)[np.asarray(mesh.face_adjacency_edges)].mean(axis=1)
    selected_edges = np.linalg.norm(midpoints - center[None, :], axis=1) <= influence_radius_um
    local_angles = angles[selected_edges]
    vertices = np.asarray(mesh.vertices, dtype=float)
    selected_vertices = np.linalg.norm(vertices - center[None, :], axis=1) <= influence_radius_um
    roughness: list[float] = []
    for vertex_id in np.flatnonzero(selected_vertices):
        neighbors = np.asarray(mesh.vertex_neighbors[int(vertex_id)], dtype=np.int64)
        if not len(neighbors):
            continue
        distances = np.linalg.norm(vertices[neighbors] - vertices[vertex_id], axis=1)
        scale = float(np.mean(distances))
        if scale > 0:
            roughness.append(float(np.linalg.norm(vertices[vertex_id] - vertices[neighbors].mean(0)) / scale))
    values = np.asarray(roughness, dtype=float)
    return {
        "junction_normal_edge_count": int(len(local_angles)),
        "normal_jump_p95_deg": float(np.percentile(local_angles, 95)) if len(local_angles) else None,
        "normal_jump_p99_deg": float(np.percentile(local_angles, 99)) if len(local_angles) else None,
        "transition_roughness_mean": float(np.mean(values)) if len(values) else None,
        "transition_roughness_p95": float(np.percentile(values, 95)) if len(values) else None,
    }


def section_radius_metrics(
    mesh: trimesh.Trimesh,
    locations: Iterable[SectionLocation],
) -> tuple[list[dict[str, Any]], dict[str, float | int | None]]:
    rows: list[dict[str, Any]] = []
    for location in locations:
        section = _section_polygon(
            mesh,
            np.asarray(location.center_xyz_um, dtype=float),
            np.asarray(location.tangent_xyz, dtype=float),
        )
        area = float(section[0]) if section is not None else None
        equivalent = float(np.sqrt(area / np.pi)) if area is not None else None
        error = (
            (equivalent - location.target_radius_um) / location.target_radius_um
            if equivalent is not None
            else None
        )
        rows.append(
            {
                "location": location.name,
                "branch_role": location.branch_role,
                "center_x_um": location.center_xyz_um[0],
                "center_y_um": location.center_xyz_um[1],
                "center_z_um": location.center_xyz_um[2],
                "target_radius_um": location.target_radius_um,
                "section_area_um2": area,
                "equivalent_radius_um": equivalent,
                "radius_relative_error": error,
                "radius_absolute_relative_error": abs(error) if error is not None else None,
            }
        )
    errors = np.asarray(
        [
            row["radius_absolute_relative_error"]
            for row in rows
            if row["branch_role"] != "junction_core"
            and row["radius_absolute_relative_error"] is not None
        ],
        dtype=float,
    )
    collar_errors = np.asarray(
        [
            row["radius_absolute_relative_error"]
            for row in rows
            if row["branch_role"] != "junction_core"
            and row["radius_absolute_relative_error"] is not None
        ],
        dtype=float,
    )
    return rows, {
        "radius_sample_count": int(len(errors)),
        "radius_p95_absolute_relative_error": (
            float(np.percentile(errors, 95)) if len(errors) else None
        ),
        "collar_max_radius_absolute_relative_error": (
            float(np.max(collar_errors)) if len(collar_errors) else None
        ),
        "radius_energy_mean_squared_um2": float(
            np.mean(
                [
                    (float(row["equivalent_radius_um"]) - float(row["target_radius_um"])) ** 2
                    for row in rows
                    if row["branch_role"] != "junction_core"
                    and row["equivalent_radius_um"] is not None
                ]
            )
        ) if len(errors) else None,
    }


def unified_mesh_qc(
    mesh: trimesh.Trimesh,
    *,
    junctions: list[tuple[int, np.ndarray, float]],
    internal_caps_known: int | None = None,
) -> dict[str, Any]:
    defects, _ = diagnose_mesh_defects(
        mesh,
        junctions,
        ray_sample_limit=1024,
    )
    quality = triangle_quality(mesh)["summary"]
    internal_faces = int(defects["suspected_internal_face_count"])
    internal_caps = (
        int(internal_caps_known)
        if internal_caps_known is not None
        else (0 if mesh.is_watertight and internal_faces == 0 else None)
    )
    return {
        "self_intersections": int(defects["self_intersection_count"]),
        "internal_faces": internal_faces,
        "internal_caps": internal_caps,
        "boundary_edges": int(defects["boundary_edge_count"]),
        "non_manifold_edges": int(defects["non_manifold_edge_count"]),
        "surface_components": int(defects["surface_connected_component_count"]),
        "degenerate_triangles": int(quality["degenerate_triangle_count"]),
        "triangle_count": int(len(mesh.faces)),
        "surface_volume_um3": float(abs(mesh.volume)) if mesh.is_watertight else None,
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "topology_qc_pass": bool(
            defects["self_intersection_count"] == 0
            and internal_faces == 0
            and internal_caps == 0
            and defects["boundary_edge_count"] == 0
            and defects["non_manifold_edge_count"] == 0
            and defects["surface_connected_component_count"] == 1
            and quality["degenerate_triangle_count"] == 0
        ),
    }


def mask_anisotropy_metrics(binary_mask: np.ndarray, mesh: trimesh.Trimesh) -> dict[str, Any]:
    mask = np.asarray(binary_mask, dtype=bool)
    slice_rows: dict[str, list[float]] = {}
    for axis, name, pixel_area in (
        (2, "x", SPACING_XYZ_UM[1] * SPACING_XYZ_UM[2]),
        (1, "y", SPACING_XYZ_UM[0] * SPACING_XYZ_UM[2]),
        (0, "z", SPACING_XYZ_UM[0] * SPACING_XYZ_UM[1]),
    ):
        other = tuple(index for index in range(3) if index != axis)
        areas = np.sum(mask, axis=other) * pixel_area
        slice_rows[name] = np.asarray(areas, dtype=float).tolist()
    normals = np.abs(np.asarray(mesh.face_normals, dtype=float))
    area = np.asarray(mesh.area_faces, dtype=float)
    total_area = float(area.sum())
    axis_locked = np.max(normals, axis=1) >= 0.99
    z_locked = normals[:, 2] >= 0.99
    z_areas = np.asarray(slice_rows["z"], dtype=float)
    occupied_z = z_areas[z_areas > 0]
    return {
        "spacing_xyz_um": list(SPACING_XYZ_UM),
        "z_spacing_to_xy_spacing_ratio": 2.0,
        "axis_locked_normal_area_fraction": (
            float(area[axis_locked].sum() / total_area) if total_area else None
        ),
        "z_locked_normal_area_fraction": (
            float(area[z_locked].sum() / total_area) if total_area else None
        ),
        "occupied_z_slice_count": int(len(occupied_z)),
        "z_slice_area_successive_absolute_change_mean_um2": (
            float(np.mean(np.abs(np.diff(occupied_z)))) if len(occupied_z) > 1 else None
        ),
        "slice_area_profiles_um2": slice_rows,
        "interpretation": (
            "Axis-locked facets and 2-um z sampling quantify voxel staircase; they are "
            "sampling artifacts, not independent vascular-wall detail."
        ),
    }


def provenance_summary_figure(rows: list[dict[str, Any]], path: Path) -> Path:
    cohorts = sorted({str(row["cohort"]) for row in rows})
    classes = [
        "SWC_DERIVED_ANNOTATION_MASK",
        "UNET_PREDICTED_MASK",
        "UNKNOWN_MASK_PROVENANCE",
    ]
    colors = ("#2b8cbe", "#e34a33", "#969696")
    figure, axis = plt.subplots(figsize=(10, 5.2))
    bottom = np.zeros(len(cohorts), dtype=float)
    for name, color in zip(classes, colors):
        values = np.asarray(
            [sum(row["cohort"] == cohort and row["mask_provenance"] == name for row in rows) for cohort in cohorts],
            dtype=float,
        )
        axis.bar(cohorts, values, bottom=bottom, color=color, label=name)
        bottom += values
    axis.set_ylabel("Mask files")
    axis.set_title("Dataset Mask provenance (document + file-structure evidence)")
    axis.tick_params(axis="x", rotation=18)
    axis.legend(fontsize=8, loc="upper left")
    axis.grid(axis="y", alpha=0.2)
    figure.tight_layout()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=180)
    plt.close(figure)
    return target


def swc_mask_overlay_figure(original: np.ndarray, reconstructed: np.ndarray, path: Path) -> Path:
    truth = np.asarray(original, dtype=float) / 255.0
    candidate = np.asarray(reconstructed, dtype=float) / 255.0
    truth_projection = truth.max(axis=0)
    candidate_projection = candidate.max(axis=0)
    overlay = np.zeros((*truth_projection.shape, 3), dtype=float)
    overlay[..., 0] = truth_projection
    overlay[..., 1] = candidate_projection
    figure, axes = plt.subplots(1, 3, figsize=(13, 4.4))
    axes[0].imshow(truth_projection, cmap="gray", origin="lower")
    axes[0].set_title("Original annotation map")
    axes[1].imshow(candidate_projection, cmap="gray", origin="lower")
    axes[1].set_title("SWC-reconstructed map")
    axes[2].imshow(np.clip(overlay, 0.0, 1.0), origin="lower")
    axes[2].set_title("Overlay: original red / SWC green")
    for axis in axes:
        axis.set_xlabel("x voxel")
        axis.set_ylabel("y voxel")
    figure.tight_layout()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=180)
    plt.close(figure)
    return target


def _polydata(mesh: trimesh.Trimesh) -> pv.PolyData:
    faces = np.column_stack(
        (np.full(len(mesh.faces), 3, dtype=np.int64), np.asarray(mesh.faces, dtype=np.int64))
    ).ravel()
    return pv.PolyData(np.asarray(mesh.vertices, dtype=float), faces)


def abc_overlay_figure(
    swc_only: trimesh.Trimesh,
    mask_only: trimesh.Trimesh,
    assisted: trimesh.Trimesh,
    swc_path: Path,
    path: Path,
) -> Path:
    plotter = pv.Plotter(off_screen=True, shape=(1, 3), window_size=(1800, 620))
    centerline = pv.read(swc_path)
    panels = (
        (swc_only, "#3182bd", "A  SWC_ONLY"),
        (mask_only, "#de2d26", "B  MASK_ONLY"),
        (assisted, "#31a354", "C  SWC_MASK_ASSISTED"),
    )
    for index, (mesh, color, title) in enumerate(panels):
        plotter.subplot(0, index)
        plotter.add_mesh(_polydata(mesh), color=color, smooth_shading=True, opacity=0.88)
        if index != 1:
            plotter.add_mesh(_polydata(mask_only), color="#fb6a4a", style="wireframe", opacity=0.35)
        plotter.add_mesh(centerline, color="black", line_width=4)
        plotter.add_text(title, font_size=12)
        plotter.show_axes()
    plotter.link_views()
    plotter.view_isometric()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    plotter.screenshot(target)
    plotter.close()
    return target


def _orthogonal_basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    unit = np.asarray(normal, dtype=float)
    unit /= np.linalg.norm(unit)
    helper = np.asarray((1.0, 0.0, 0.0)) if abs(unit[0]) < 0.8 else np.asarray((0.0, 1.0, 0.0))
    first = np.cross(unit, helper)
    first /= np.linalg.norm(first)
    return first, np.cross(unit, first)


def _mesh_section_lines(
    mesh: trimesh.Trimesh,
    center: np.ndarray,
    normal: np.ndarray,
    first: np.ndarray,
    second: np.ndarray,
) -> list[np.ndarray]:
    section = mesh.section(plane_origin=center, plane_normal=normal)
    if section is None:
        return []
    output: list[np.ndarray] = []
    for line in section.discrete:
        relative = np.asarray(line, dtype=float) - center[None, :]
        output.append(np.column_stack((relative @ first, relative @ second)))
    return output


def cross_section_figure(
    probability_mask: np.ndarray,
    *,
    origin_xyz_um: tuple[float, float, float],
    locations: list[SectionLocation],
    meshes: dict[str, trimesh.Trimesh],
    path: Path,
    half_width_um: float = 5.0,
) -> Path:
    from scipy.ndimage import map_coordinates

    colors = {
        "SWC_ONLY": "#3182bd",
        "MASK_ONLY": "#de2d26",
        "SWC_MASK_ASSISTED": "#31a354",
    }
    figure, axes = plt.subplots(1, len(locations), figsize=(5.0 * len(locations), 4.8))
    axes = np.atleast_1d(axes)
    coordinate = np.linspace(-half_width_um, half_width_um, 181)
    uu, vv = np.meshgrid(coordinate, coordinate, indexing="xy")
    origin = np.asarray(origin_xyz_um, dtype=float)
    for axis, location in zip(axes, locations):
        center = np.asarray(location.center_xyz_um, dtype=float)
        normal = np.asarray(location.tangent_xyz, dtype=float)
        normal /= np.linalg.norm(normal)
        first, second = _orthogonal_basis(normal)
        points = center[None, None, :] + uu[..., None] * first + vv[..., None] * second
        coordinates_xyz = (points.reshape((-1, 3)) - origin[None, :]) / np.asarray(SPACING_XYZ_UM)
        values = map_coordinates(
            np.asarray(probability_mask, dtype=float),
            coordinates_xyz[:, ::-1].T,
            order=1,
            mode="constant",
            cval=0.0,
            prefilter=False,
        ).reshape(uu.shape)
        axis.contour(uu, vv, values, levels=[0.404 * 255.0], colors="black", linewidths=2.0)
        circle = plt.Circle(
            (0.0, 0.0),
            location.target_radius_um,
            fill=False,
            color="#756bb1",
            linestyle="--",
            linewidth=2.0,
            label="SWC target radius",
        )
        axis.add_patch(circle)
        for name, mesh in meshes.items():
            for line_number, line in enumerate(
                _mesh_section_lines(mesh, center, normal, first, second)
            ):
                axis.plot(
                    line[:, 0],
                    line[:, 1],
                    color=colors[name],
                    linewidth=1.4,
                    label=name if line_number == 0 else None,
                )
        axis.set_title(f"{location.name}\n{location.branch_role}")
        axis.set_aspect("equal")
        axis.set_xlim(-half_width_um, half_width_um)
        axis.set_ylim(-half_width_um, half_width_um)
        axis.set_xlabel("u (um)")
        axis.set_ylabel("v (um)")
        axis.grid(alpha=0.2)
    handles = [
        plt.Line2D((0,), (0,), color="black", linewidth=2, label="Dataset Mask boundary"),
        plt.Line2D((0,), (0,), color="#756bb1", linestyle="--", label="SWC target radius"),
        *[plt.Line2D((0,), (0,), color=color, label=name) for name, color in colors.items()],
    ]
    figure.legend(handles=handles, loc="lower center", ncol=5, fontsize=9)
    figure.tight_layout(rect=(0, 0.10, 1, 1))
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=190)
    plt.close(figure)
    return target


def radius_comparison_figure(rows: list[dict[str, Any]], path: Path) -> Path:
    methods = [str(row["method"]) for row in rows]
    radius = [100.0 * float(row["radius_p95_absolute_relative_error"]) for row in rows]
    collar = [100.0 * float(row["collar_max_radius_absolute_relative_error"]) for row in rows]
    x = np.arange(len(methods))
    figure, axis = plt.subplots(figsize=(10, 5.2))
    axis.bar(x - 0.2, radius, 0.4, label="branch-local radius P95")
    axis.bar(x + 0.2, collar, 0.4, label="collar max radius error")
    axis.set_xticks(x, methods, rotation=18)
    axis.set_ylabel("Absolute relative error (%)")
    axis.set_title("SWC radius fidelity under identical section evaluation")
    axis.legend()
    axis.grid(axis="y", alpha=0.2)
    figure.tight_layout()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=180)
    plt.close(figure)
    return target


def anisotropy_figure(metrics: dict[str, Any], path: Path) -> Path:
    profiles = metrics["slice_area_profiles_um2"]
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    for name, spacing in (("x", 1.0), ("y", 1.0), ("z", 2.0)):
        values = np.asarray(profiles[name], dtype=float)
        axes[0].plot(np.arange(len(values)) * spacing, values, label=f"{name} slices ({spacing:g} um)")
    axes[0].set_xlabel("Local coordinate (um)")
    axes[0].set_ylabel("Foreground slice area (um2)")
    axes[0].set_title("Mask slice quantization")
    axes[0].legend()
    labels = ("all axis-locked", "z-normal locked")
    values = (
        float(metrics["axis_locked_normal_area_fraction"]),
        float(metrics["z_locked_normal_area_fraction"]),
    )
    axes[1].bar(labels, values, color=("#756bb1", "#e6550d"))
    axes[1].set_ylim(0.0, max(0.05, 1.1 * max(values)))
    axes[1].set_ylabel("Mask-only surface area fraction")
    axes[1].set_title("Marching-cubes facet axis locking")
    axes[1].tick_params(axis="x", rotation=12)
    for axis in axes:
        axis.grid(axis="y", alpha=0.2)
    figure.tight_layout()
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(target, dpi=180)
    plt.close(figure)
    return target
