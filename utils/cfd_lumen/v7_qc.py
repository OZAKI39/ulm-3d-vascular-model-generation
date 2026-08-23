"""Topology and former-merge-ring diagnostics for unified v7 surfaces."""

from __future__ import annotations

from typing import Any

import numpy as np
import trimesh

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .local_implicit_junction import sample_from_junction
from .mesh_defects import diagnose_mesh_defects, triangle_quality
from .types import BranchGeometry, HybridBuildDetails


def evaluate_unified_topology(
    mesh: trimesh.Trimesh,
    roi: ROIRecord,
    v6_details: HybridBuildDetails,
    config: CFDLumenConfig,
) -> dict[str, Any]:
    junctions = [
        (
            int(node_id),
            np.asarray(roi.local_node_positions_um[node_id], dtype=float),
            float(roi.local_node_radius_um[node_id]),
        )
        for node_id in v6_details.patches
    ]
    defects, _ = diagnose_mesh_defects(
        mesh, junctions, ray_sample_limit=2_048
    )
    quality = defects["triangle_quality"]
    checks = {
        "watertight": bool(mesh.is_watertight),
        "single_component": defects["surface_connected_component_count"] == 1,
        "zero_boundary_edges": defects["boundary_edge_count"] == 0,
        "zero_nonmanifold_edges": defects["non_manifold_edge_count"] == 0,
        "zero_self_intersections": defects["self_intersection_count"] == 0,
        "zero_internal_faces": defects["suspected_internal_face_count"] == 0,
        # A unified field creates no intermediate cap primitive. The only caps
        # are the labeled external CFD ports, so an internal cap count is zero
        # when the internal-face test is zero.
        "zero_internal_caps": defects["suspected_internal_face_count"] == 0,
        "zero_degenerate_triangles": quality["degenerate_triangle_count"] == 0,
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
    return {
        "status": "PASS"
        if all(checks[name] or not required[name] for name in checks)
        else "FAIL",
        "backend": "unified_polyball",
        "checks": checks,
        "self_intersection_pairs": defects["self_intersection_count"],
        "internal_face_count": defects["suspected_internal_face_count"],
        "internal_cap_face_count": 0
        if defects["suspected_internal_face_count"] == 0
        else defects["suspected_internal_face_count"],
        "boundary_edge_count": defects["boundary_edge_count"],
        "nonmanifold_edge_count": defects["non_manifold_edge_count"],
        "surface_component_count": defects["surface_connected_component_count"],
        "degenerate_triangle_count": quality["degenerate_triangle_count"],
        "mesh_defects": defects,
    }


def basic_unified_topology(mesh: trimesh.Trimesh) -> dict[str, Any]:
    _, counts = np.unique(np.asarray(mesh.edges_sorted), axis=0, return_counts=True)
    quality = triangle_quality(mesh)["summary"]
    return {
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "boundary_edge_count": int(np.count_nonzero(counts == 1)),
        "nonmanifold_edge_count": int(np.count_nonzero(counts > 2)),
        "component_count": int(len(mesh.split(only_watertight=False))),
        "degenerate_triangle_count": int(quality["degenerate_triangle_count"]),
    }


def global_wall_metrics(
    mesh: trimesh.Trimesh,
    wall_face_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    adjacency = np.asarray(mesh.face_adjacency, dtype=np.int64)
    angles = np.degrees(np.asarray(mesh.face_adjacency_angles, dtype=float))
    if wall_face_mask is not None:
        wall = np.asarray(wall_face_mask, dtype=bool)
        selected = wall[adjacency[:, 0]] & wall[adjacency[:, 1]]
        angles = angles[selected]
        face_ids = np.flatnonzero(wall)
    else:
        face_ids = None
    quality = triangle_quality(mesh, face_ids)["summary"]
    return {
        "normal_jump_mean_deg": float(np.mean(angles)) if len(angles) else None,
        "normal_jump_p95_deg": float(np.percentile(angles, 95)) if len(angles) else None,
        "normal_jump_p99_deg": float(np.percentile(angles, 99)) if len(angles) else None,
        "normal_jump_max_deg": float(np.max(angles)) if len(angles) else None,
        "triangle_quality": quality,
    }


def former_merge_ring_rows(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    branches: list[BranchGeometry],
    *,
    version: str,
) -> list[dict[str, Any]]:
    """Measure geometry at every v6 weld location without inventing a v7 seam."""

    branch_by_id = {branch.branch_id: branch for branch in branches}
    adjacency = np.asarray(mesh.face_adjacency, dtype=np.int64)
    edges = np.asarray(mesh.face_adjacency_edges, dtype=np.int64)
    edge_points = np.asarray(mesh.vertices, dtype=float)[edges]
    midpoint = edge_points.mean(axis=1)
    angle = np.degrees(np.asarray(mesh.face_adjacency_angles, dtype=float))
    quality = triangle_quality(mesh)
    rows: list[dict[str, Any]] = []
    for node_id, patch in details.patches.items():
        for collar in patch.collars:
            branch = branch_by_id[collar.branch_id]
            distance = float(
                collar.collar_distance_um + 0.5 * collar.overlap_length_um
            )
            center, radius, tangent = sample_from_junction(
                branch, collar.endpoint_index, distance
            )
            relative = midpoint - center[None, :]
            axial = relative @ tangent
            radial = np.linalg.norm(
                relative - axial[:, None] * tangent[None, :], axis=1
            )
            selected = (np.abs(axial) <= 0.35 * radius) & (
                radial <= 1.75 * radius
            )
            local_faces = np.unique(adjacency[selected]) if np.any(selected) else np.empty(0, dtype=np.int64)
            local_aspect = quality["aspect_ratio"][local_faces]
            local_angle = angle[selected]
            rows.append(
                {
                    "version": version,
                    "location": "V6_PURE_BRANCH_TO_EXPLICIT_WELD_WINDOW",
                    "junction_node_id": int(node_id),
                    "branch_id": int(collar.branch_id),
                    "reconstruction_interface_present": version == "v6",
                    "hybrid_interface_edge_count": int(np.count_nonzero(selected))
                    if version == "v6"
                    else 0,
                    "measurement_edge_count": int(np.count_nonzero(selected)),
                    "normal_jump_p95_deg": float(np.percentile(local_angle, 95))
                    if len(local_angle)
                    else None,
                    "normal_jump_p99_deg": float(np.percentile(local_angle, 99))
                    if len(local_angle)
                    else None,
                    "normal_jump_max_deg": float(np.max(local_angle))
                    if len(local_angle)
                    else None,
                    "triangle_aspect_ratio_p95": float(np.percentile(local_aspect, 95))
                    if len(local_aspect)
                    else None,
                    "triangle_aspect_ratio_max": float(np.max(local_aspect))
                    if len(local_aspect)
                    else None,
                }
            )
    return rows


def summarize_ring_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    def maximum(key: str) -> float | None:
        values = [float(row[key]) for row in rows if row.get(key) is not None]
        return max(values, default=None)

    return {
        "ring_count": len(rows),
        "hybrid_interface_edge_count": int(
            sum(int(row["hybrid_interface_edge_count"]) for row in rows)
        ),
        "worst_normal_jump_p95_deg": maximum("normal_jump_p95_deg"),
        "worst_normal_jump_p99_deg": maximum("normal_jump_p99_deg"),
        "worst_normal_jump_max_deg": maximum("normal_jump_max_deg"),
        "worst_triangle_aspect_ratio_p95": maximum("triangle_aspect_ratio_p95"),
        "worst_triangle_aspect_ratio_max": maximum("triangle_aspect_ratio_max"),
    }
