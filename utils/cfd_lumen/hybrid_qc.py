"""Hybrid-specific topology, cap, collar, and junction area-profile QC."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv
import trimesh

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .export import write_csv, write_json
from .local_implicit_junction import sample_from_junction
from .mesh_defects import diagnose_mesh_defects
from .surface_qc import _section_polygon
from .types import BranchGeometry, HybridBuildDetails
from .visualization import _equal_3d, _save, _surface_collection


def _mesh_polydata(mesh: trimesh.Trimesh) -> pv.PolyData:
    faces = np.column_stack((np.full(len(mesh.faces), 3, dtype=np.int64), mesh.faces)).ravel()
    return pv.PolyData(np.asarray(mesh.vertices), faces)


def _cap_candidate_face_ids(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
) -> tuple[np.ndarray, list[tuple[dict[str, Any], np.ndarray]]]:
    centers = np.asarray(mesh.triangles_center)
    normals = np.asarray(mesh.face_normals)
    rows: list[tuple[dict[str, Any], np.ndarray]] = []
    for cap in details.explicit_cap_planes:
        center = np.asarray(cap["center_um"], dtype=float)
        normal = np.asarray(cap["normal_away_from_junction"], dtype=float)
        radius = float(cap["radius_um"])
        relative = centers - center
        axial = np.abs(relative @ normal)
        radial = np.linalg.norm(relative - (relative @ normal)[:, None] * normal, axis=1)
        alignment = np.abs(normals @ normal)
        tolerance = max(radius * 1.0e-4, 1.0e-8)
        candidates = np.flatnonzero(
            (axial <= tolerance) & (radial <= 1.01 * radius) & (alignment >= 0.95)
        )
        rows.append((cap, candidates))
    combined = np.unique(
        np.concatenate([row[1] for row in rows]) if rows else np.empty(0, dtype=np.int64)
    )
    return combined, rows


def _confirmed_internal_caps(
    cap_candidates: list[tuple[dict[str, Any], np.ndarray]],
    internal_face_ids: np.ndarray,
) -> tuple[int, list[dict[str, Any]]]:
    internal = set(map(int, internal_face_ids))
    rows: list[dict[str, Any]] = []
    total = 0
    for cap, candidates in cap_candidates:
        confirmed = sum(int(face_id) in internal for face_id in candidates)
        total += confirmed
        rows.append(
            {
                "junction_node_id": cap["junction_node_id"],
                "branch_id": cap["branch_id"],
                "candidate_cap_plane_faces": int(len(candidates)),
                "confirmed_internal_cap_faces": int(confirmed),
            }
        )
    return total, rows


def evaluate_hybrid_surface_qc(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    roi: ROIRecord,
    config: CFDLumenConfig,
) -> tuple[dict[str, Any], dict[str, Any]]:
    junctions = [
        (
            node_id,
            np.asarray(roi.local_node_positions_um[node_id], dtype=float),
            float(roi.local_node_radius_um[node_id]),
        )
        for node_id in details.patches
    ]
    cap_face_ids, cap_candidates = _cap_candidate_face_ids(mesh, details)
    defects, artifacts = diagnose_mesh_defects(
        mesh,
        junctions,
        ray_face_ids=cap_face_ids,
        ray_sample_limit=1024,
    )
    internal_caps, cap_rows = _confirmed_internal_caps(
        cap_candidates, artifacts["suspected_internal_faces"]
    )
    triangle_quality = defects["triangle_quality"]
    checks = {
        "watertight": bool(mesh.is_watertight),
        "single_component": defects["surface_connected_component_count"] == 1,
        "zero_boundary_edges": defects["boundary_edge_count"] == 0,
        "zero_nonmanifold_edges": defects["non_manifold_edge_count"] == 0,
        "zero_self_intersections": defects["self_intersection_count"] == 0,
        "zero_internal_faces": defects["suspected_internal_face_count"] == 0,
        "zero_internal_caps": internal_caps == 0,
        "zero_degenerate_triangles": triangle_quality["degenerate_triangle_count"] == 0,
    }
    required = {
        "watertight": config.surface_qc.require_watertight,
        "single_component": config.surface_qc.require_single_component,
        "zero_boundary_edges": config.surface_qc.require_zero_boundary_edges,
        "zero_nonmanifold_edges": config.surface_qc.require_zero_nonmanifold_edges,
        "zero_self_intersections": config.surface_qc.require_zero_self_intersections,
        "zero_internal_faces": config.surface_qc.require_zero_internal_faces,
        "zero_internal_caps": config.surface_qc.require_zero_internal_caps,
        "zero_degenerate_triangles": True,
    }
    status = "PASS" if all(checks[name] or not required[name] for name in checks) else "FAIL"
    report = {
        "status": status,
        "backend": (
            "hybrid_loop_stitch"
            if details.transition_backend == "loop_stitch"
            else "hybrid_manifold_boolean_fallback"
            if details.transition_fallback_reason
            else "hybrid_local_implicit"
        ),
        "transition_backend": details.transition_backend,
        "transition_fallback_reason": details.transition_fallback_reason,
        "checks": checks,
        "self_intersection_pairs": defects["self_intersection_count"],
        "internal_face_count": defects["suspected_internal_face_count"],
        "internal_cap_face_count": internal_caps,
        "degenerate_triangle_count": triangle_quality["degenerate_triangle_count"],
        "boundary_edge_count": defects["boundary_edge_count"],
        "nonmanifold_edge_count": defects["non_manifold_edge_count"],
        "surface_component_count": defects["surface_connected_component_count"],
        "triangle_count": int(len(mesh.faces)),
        "surface_volume_um3": float(abs(mesh.volume)) if mesh.is_watertight else None,
        "cap_planes": cap_rows,
        "mesh_defects": defects,
    }
    return report, artifacts


def collar_radius_rows(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    branches: list[BranchGeometry],
) -> list[dict[str, Any]]:
    branch_by_id = {branch.branch_id: branch for branch in branches}
    rows: list[dict[str, Any]] = []
    for patch in details.patches.values():
        for collar in patch.collars:
            branch = branch_by_id[collar.branch_id]
            diameter = 2.0 * collar.collar_radius_um
            for offset_label, offset in (("collar_minus_0.5D", -0.5 * diameter), ("collar", 0.0), ("collar_plus_0.5D", 0.5 * diameter)):
                distance = float(np.clip(collar.collar_distance_um + offset, 0.0, branch.length_um))
                point, source_radius, tangent = sample_from_junction(
                    branch, collar.endpoint_index, distance
                )
                section = _section_polygon(mesh, point, tangent)
                area = float(section[0]) if section else None
                reconstructed = float(np.sqrt(area / np.pi)) if area else None
                error = (
                    (reconstructed - source_radius) / source_radius
                    if reconstructed is not None
                    else None
                )
                rows.append(
                    {
                        "junction_node_id": collar.junction_node_id,
                        "branch_id": collar.branch_id,
                        "sample_location": offset_label,
                        "distance_from_junction_um": distance,
                        "collar_distance_um": collar.collar_distance_um,
                        "source_radius_um": source_radius,
                        "cross_section_area_um2": area,
                        "reconstructed_radius_um": reconstructed,
                        "radius_relative_error": error,
                        "absolute_radius_relative_error": abs(error) if error is not None else None,
                    }
                )
    return rows


def junction_area_profile_rows(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    branches: list[BranchGeometry],
    *,
    samples_per_branch: int = 25,
) -> list[dict[str, Any]]:
    branch_by_id = {branch.branch_id: branch for branch in branches}
    rows: list[dict[str, Any]] = []
    for patch in details.patches.values():
        for collar in patch.collars:
            branch = branch_by_id[collar.branch_id]
            diameter = 2.0 * collar.collar_radius_um
            maximum = min(3.0 * diameter, 0.9 * branch.length_um)
            minimum = min(0.1 * diameter, 0.5 * maximum)
            for sample_index, distance in enumerate(np.linspace(minimum, maximum, samples_per_branch)):
                point, source_radius, tangent = sample_from_junction(
                    branch, collar.endpoint_index, float(distance)
                )
                section = _section_polygon(mesh, point, tangent)
                area = float(section[0]) if section else None
                source_area = float(np.pi * source_radius**2)
                if distance < collar.explicit_cap_distance_um:
                    region = "implicit"
                elif distance <= collar.implicit_extent_um:
                    region = "collar_overlap"
                else:
                    region = "explicit"
                rows.append(
                    {
                        "junction_node_id": collar.junction_node_id,
                        "branch_id": collar.branch_id,
                        "sample_index": sample_index,
                        "distance_from_junction_um": float(distance),
                        "distance_in_diameters": float(distance / diameter),
                        "region": region,
                        "collar_distance_um": collar.collar_distance_um,
                        "source_radius_um": source_radius,
                        "source_area_um2": source_area,
                        "cross_section_area_um2": area,
                        "area_to_source_ratio": area / source_area if area is not None else None,
                    }
                )
    return rows


def max_junction_area_ratio(rows: list[dict[str, Any]]) -> float | None:
    values = [row["area_to_source_ratio"] for row in rows if row["area_to_source_ratio"] is not None]
    return max(values, default=None)


def write_hybrid_artifacts(
    hybrid_root: Path,
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    qc: dict[str, Any],
    collar_rows: list[dict[str, Any]],
    area_rows: list[dict[str, Any]],
) -> list[Path]:
    junction_root = hybrid_root / "junctions"
    figures = hybrid_root / "figures"
    qc_root = hybrid_root / "qc"
    for folder in (junction_root, figures, qc_root):
        folder.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for node_id, patch in details.patches.items():
        root = junction_root / f"junction_{node_id}"
        root.mkdir(parents=True, exist_ok=True)
        centerline = pv.PolyData(patch.centerline_points_um)
        centerline.point_data["radius_um"] = patch.centerline_radius_um
        centerline_path = root / "centerline.vtp"
        centerline.save(centerline_path)
        raw_path = root / "implicit_raw.vtp"
        clean_path = root / "implicit_clean.vtp"
        merged_path = root / "merged.vtp"
        _mesh_polydata(patch.raw_mesh).save(raw_path)
        _mesh_polydata(patch.clean_mesh).save(clean_path)
        _mesh_polydata(details.merged_junction_meshes[node_id]).save(merged_path)
        metadata_path = write_json(root / "local_field_metadata.json", patch.metadata)
        paths.extend((centerline_path, raw_path, clean_path, merged_path, metadata_path))

        figure = plt.figure(figsize=(9, 7))
        axis = figure.add_subplot(111, projection="3d")
        points = patch.centerline_points_um
        axis.scatter(points[:, 0], points[:, 1], points[:, 2], c=patch.centerline_radius_um, s=7, cmap="viridis")
        _equal_3d(axis, points)
        axis.set_title(f"Junction {node_id}: local polyball field samples")
        paths.append(_save(figure, figures / f"junction_{node_id}_field.png"))

        for candidate, suffix, title in (
            (patch.clean_mesh, "implicit", "local implicit patch"),
            (details.merged_junction_meshes[node_id], "hybrid", "patch + explicit collars"),
        ):
            figure = plt.figure(figsize=(9, 7))
            axis = figure.add_subplot(111, projection="3d")
            _surface_collection(axis, candidate, color="#67E8F9", alpha=0.45)
            _equal_3d(axis, np.asarray(candidate.vertices))
            axis.set_title(f"Junction {node_id}: {title}")
            paths.append(_save(figure, figures / f"junction_{node_id}_{suffix}.png"))

    profile_path = write_csv(qc_root / "junction_area_profile.csv", area_rows)
    collar_path = write_csv(qc_root / "collar_radius_error.csv", collar_rows)
    qc_path = write_json(qc_root / "hybrid_surface_qc.json", qc)
    merge_path = write_json(qc_root / "merge_steps.json", details.merge_steps)
    final_path = hybrid_root / "hybrid_surface.vtp"
    _mesh_polydata(mesh).save(final_path)
    paths.extend((profile_path, collar_path, qc_path, merge_path, final_path))

    figure, axis = plt.subplots(figsize=(11, 6))
    for key in sorted({(row["junction_node_id"], row["branch_id"]) for row in area_rows}):
        selected = [row for row in area_rows if (row["junction_node_id"], row["branch_id"]) == key and row["area_to_source_ratio"] is not None]
        axis.plot(
            [row["distance_in_diameters"] for row in selected],
            [row["area_to_source_ratio"] for row in selected],
            marker=".",
            label=f"J{key[0]} / B{key[1]}",
        )
    axis.axhline(1.0, color="black", linewidth=1, linestyle="--", label="source area")
    axis.set(xlabel="distance from junction (D)", ylabel="A / source A", title="Junction area profiles")
    axis.legend(ncol=2, fontsize=8)
    paths.append(_save(figure, figures / "junction_area_profile.png"))

    for node_id in sorted({row["junction_node_id"] for row in area_rows}):
        figure, axis = plt.subplots(figsize=(11, 6))
        for branch_id in sorted(
            {row["branch_id"] for row in area_rows if row["junction_node_id"] == node_id}
        ):
            selected = [
                row
                for row in area_rows
                if row["junction_node_id"] == node_id
                and row["branch_id"] == branch_id
                and row["area_to_source_ratio"] is not None
            ]
            axis.plot(
                [row["distance_in_diameters"] for row in selected],
                [row["area_to_source_ratio"] for row in selected],
                marker=".",
                label=f"branch {branch_id}",
            )
        axis.axhline(1.0, color="black", linewidth=1, linestyle="--", label="source area")
        axis.set(
            xlabel="distance from junction (D)",
            ylabel="A / source A",
            title=f"Junction {node_id} area profile",
        )
        axis.legend(fontsize=8)
        paths.append(_save(figure, figures / f"junction_{node_id}_area_profile.png"))

    figure, axis = plt.subplots(figsize=(10, 5))
    labels = [f"J{row['junction_node_id']}/B{row['branch_id']}/{row['sample_location']}" for row in collar_rows]
    values = [100.0 * row["absolute_radius_relative_error"] if row["absolute_radius_relative_error"] is not None else np.nan for row in collar_rows]
    axis.bar(np.arange(len(labels)), values, color="#0EA5E9")
    axis.set_xticks(np.arange(len(labels)), labels, rotation=75, ha="right", fontsize=7)
    axis.set(ylabel="absolute radius error (%)", title="Hybrid collar transition")
    paths.append(_save(figure, figures / "collar_transition.png"))
    return paths
