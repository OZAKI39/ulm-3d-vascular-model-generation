"""v5 port/collar normal-continuity and transition roughness measurements."""

from __future__ import annotations

from typing import Any

import numpy as np
import trimesh

from .branch_local_qc import _BranchDistanceIndex, _all_branch_distances, _ownership
from .config import CFDLumenConfig
from .local_implicit_junction import sample_from_junction
from .surface_qc import _section_polygon
from .surface_transition import (
    EXPLICIT_BRANCH_FACE,
    JUNCTION_CORE_FACE,
    TRANSITION_COLLAR_FACE,
)
from .types import BranchGeometry, HybridBuildDetails, PortGeometry


REGION_NAMES = {
    int(JUNCTION_CORE_FACE): "JUNCTION_CORE",
    int(TRANSITION_COLLAR_FACE): "TRANSITION_COLLAR",
    int(EXPLICIT_BRANCH_FACE): "EXPLICIT_BRANCH",
}


def _statistics(values: np.ndarray) -> dict[str, float | int | None]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if not len(finite):
        return {
            "sample_count": 0,
            "normal_jump_mean_deg": None,
            "normal_jump_p95_deg": None,
            "normal_jump_p99_deg": None,
            "normal_jump_max_deg": None,
        }
    return {
        "sample_count": int(len(finite)),
        "normal_jump_mean_deg": float(np.mean(finite)),
        "normal_jump_p95_deg": float(np.percentile(finite, 95)),
        "normal_jump_p99_deg": float(np.percentile(finite, 99)),
        "normal_jump_max_deg": float(np.max(finite)),
    }


def _classify_v4_faces(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    branches: list[BranchGeometry],
) -> np.ndarray:
    centers = np.asarray(mesh.triangles_center, dtype=float)
    indices = [_BranchDistanceIndex(branch) for branch in branches]
    normalized, arcs = _all_branch_distances(centers, indices)
    _, owner_column, _ = _ownership(normalized, branches)
    labels = np.full(len(mesh.faces), EXPLICIT_BRANCH_FACE, dtype=np.uint8)
    collars_by_branch: dict[int, list[Any]] = {}
    for patch in details.patches.values():
        for collar in patch.collars:
            collars_by_branch.setdefault(collar.branch_id, []).append(collar)
    for column, branch in enumerate(branches):
        selected = owner_column == column
        if not np.any(selected):
            continue
        branch_arc = arcs[selected, column]
        selected_labels = np.full(np.count_nonzero(selected), EXPLICIT_BRANCH_FACE, dtype=np.uint8)
        for collar in collars_by_branch.get(branch.branch_id, []):
            distance = (
                branch_arc
                if collar.endpoint_index == 0
                else float(branch.arc_length_um[-1]) - branch_arc
            )
            selected_labels[distance <= collar.explicit_cap_distance_um] = JUNCTION_CORE_FACE
            transition = (
                (distance > collar.explicit_cap_distance_um)
                & (distance <= collar.implicit_extent_um)
            )
            selected_labels[transition] = TRANSITION_COLLAR_FACE
        labels[np.flatnonzero(selected)] = selected_labels
    return labels


def face_region_labels(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    branches: list[BranchGeometry],
) -> np.ndarray:
    if len(details.face_region) == len(mesh.faces):
        return np.asarray(details.face_region, dtype=np.uint8)
    return _classify_v4_faces(mesh, details, branches)


def junction_normal_jump_report(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    branches: list[BranchGeometry],
) -> dict[str, Any]:
    labels = face_region_labels(mesh, details, branches)
    adjacency = np.asarray(mesh.face_adjacency, dtype=np.int64)
    angles = np.degrees(np.asarray(mesh.face_adjacency_angles, dtype=float))
    report: dict[str, Any] = {"whole_surface": _statistics(angles)}
    for code, name in REGION_NAMES.items():
        if code == int(TRANSITION_COLLAR_FACE):
            selected = (labels[adjacency[:, 0]] == code) | (
                labels[adjacency[:, 1]] == code
            )
        else:
            selected = (labels[adjacency[:, 0]] == code) & (
                labels[adjacency[:, 1]] == code
            )
        report[name] = _statistics(angles[selected])
        report[name]["triangle_count"] = int(np.count_nonzero(labels == code))
    report["face_region_codes"] = REGION_NAMES
    return report


def _port_branch(
    branches: list[BranchGeometry], port: PortGeometry
) -> tuple[BranchGeometry, int]:
    for branch in branches:
        if branch.local_node_ids[0] == port.local_node_id:
            return branch, 0
        if branch.local_node_ids[-1] == port.local_node_id:
            return branch, -1
    raise ValueError(f"Port {port.cut_port_id} has no source branch")


def _port_profile_location(
    branch: BranchGeometry,
    endpoint: int,
    port: PortGeometry,
    offset_diameters: float,
) -> tuple[np.ndarray, float, np.ndarray]:
    signed_distance = offset_diameters * 2.0 * port.radius_um
    if signed_distance >= 0:
        return (
            port.original_position_um + signed_distance * port.outward_tangent,
            port.radius_um,
            port.outward_tangent,
        )
    center, radius, tangent_inward = sample_from_junction(
        branch, endpoint, -signed_distance
    )
    return center, radius, -tangent_inward


def port_continuity_report(
    mesh: trimesh.Trimesh,
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    *,
    continuous_centerline: bool,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    adjacency_angles = np.degrees(np.asarray(mesh.face_adjacency_angles, dtype=float))
    adjacency_edge_points = np.asarray(mesh.vertices)[
        np.asarray(mesh.face_adjacency_edges, dtype=np.int64)
    ]
    adjacency_midpoints = adjacency_edge_points.mean(axis=1)
    adjacency_vectors = np.diff(adjacency_edge_points, axis=1)[:, 0, :]
    adjacency_vectors /= np.linalg.norm(adjacency_vectors, axis=1)[:, None]
    cut_angles: list[np.ndarray] = []
    for port in ports:
        branch, endpoint = _port_branch(branches, port)
        diameter = 2.0 * port.radius_um
        for offset in config.port_transition.profile_diameter_offsets:
            center, source_radius, tangent = _port_profile_location(
                branch, endpoint, port, float(offset)
            )
            section = _section_polygon(mesh, center, tangent)
            area = float(section[0]) if section else None
            source_area = float(np.pi * source_radius**2)
            relative = adjacency_midpoints - center[None, :]
            axial = relative @ tangent
            radial = np.linalg.norm(relative - axial[:, None] * tangent[None, :], axis=1)
            selected = (np.abs(axial) <= 0.25 * diameter) & (
                radial <= 1.75 * max(source_radius, port.radius_um)
            ) & (np.abs(adjacency_vectors @ tangent) <= 0.35)
            normal = _statistics(adjacency_angles[selected])
            if float(offset) == 0.0:
                cut_angles.append(adjacency_angles[selected])
            rows.append(
                {
                    "port_id": port.port_id,
                    "cut_port_id": port.cut_port_id,
                    "source_core_cut_port_id": port.source_core_cut_port_id,
                    "offset_diameters": float(offset),
                    "sample_x_um": float(center[0]),
                    "sample_y_um": float(center[1]),
                    "sample_z_um": float(center[2]),
                    "source_radius_um": source_radius,
                    "source_area_um2": source_area,
                    "cross_section_area_um2": area,
                    "area_relative_error": (
                        (area - source_area) / source_area if area is not None else None
                    ),
                    "equivalent_radius_um": (
                        float(np.sqrt(area / np.pi)) if area is not None else None
                    ),
                    **normal,
                }
            )
    pooled = np.concatenate(cut_angles) if cut_angles else np.empty(0, dtype=float)
    cut_rows = [row for row in rows if row["offset_diameters"] == 0.0]
    summary = {
        "status": "PASS" if all(row["cross_section_area_um2"] is not None for row in rows) else "FAIL",
        "port_count": len(ports),
        "profile_sample_count": len(rows),
        "maximum_cut_area_absolute_relative_error": max(
            (abs(float(row["area_relative_error"])) for row in cut_rows), default=None
        ),
        "port_normal_jump": _statistics(pooled),
        "separate_cylinder_primitive_count": (
            0 if continuous_centerline else len(ports)
        ),
        "construction": (
            "continuous centerline -> one vtkTubeFilter"
            if continuous_centerline
            else "separate branch tube + cylinder Boolean"
        ),
    }
    return rows, summary


def transition_roughness_report(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    branches: list[BranchGeometry],
) -> dict[str, Any]:
    labels = face_region_labels(mesh, details, branches)
    vertices = np.asarray(mesh.vertices, dtype=float)
    roughness = np.zeros(len(vertices), dtype=float)
    for vertex_id, neighbors in enumerate(mesh.vertex_neighbors):
        if not len(neighbors):
            roughness[vertex_id] = np.nan
            continue
        neighbor_points = vertices[np.asarray(neighbors, dtype=np.int64)]
        scale = float(np.mean(np.linalg.norm(neighbor_points - vertices[vertex_id], axis=1)))
        roughness[vertex_id] = (
            float(np.linalg.norm(vertices[vertex_id] - neighbor_points.mean(axis=0))) / scale
            if scale > 0
            else np.nan
        )
    face_roughness = np.nanmean(roughness[np.asarray(mesh.faces, dtype=np.int64)], axis=1)
    region: dict[str, Any] = {}
    for code, name in REGION_NAMES.items():
        values = face_roughness[labels == code]
        region[name] = {
            "triangle_count": int(len(values)),
            "laplacian_roughness_mean": float(np.nanmean(values)) if len(values) else None,
            "laplacian_roughness_p95": (
                float(np.nanpercentile(values, 95)) if len(values) else None
            ),
        }
    transition = region["TRANSITION_COLLAR"]["laplacian_roughness_mean"]
    explicit = region["EXPLICIT_BRANCH"]["laplacian_roughness_mean"]
    return {
        "method": "dimensionless umbrella-Laplacian magnitude",
        "regions": region,
        "transition_to_explicit_roughness_ratio": (
            transition / explicit
            if transition is not None and explicit not in {None, 0.0}
            else None
        ),
    }


def evaluate_surface_continuity(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    port_rows, port_report = port_continuity_report(
        mesh,
        branches,
        ports,
        config,
        continuous_centerline=bool(details.port_extension_rows),
    )
    report = {
        "transition_backend": details.transition_backend,
        "transition_fallback_reason": details.transition_fallback_reason,
        "port": port_report,
        "normal_jump": junction_normal_jump_report(mesh, details, branches),
        "roughness": transition_roughness_report(mesh, details, branches),
        "transition_triangle_count": sum(
            int(row.get("transition_triangle_count", 0))
            for row in details.transition_rows
        ),
        "transition_rows": details.transition_rows,
    }
    return port_rows, report
