"""Interface-local, transition-region, and silhouette QC for v6."""

from __future__ import annotations

from typing import Any

import cv2
import networkx as nx
import numpy as np
import trimesh

from .config import CFDLumenConfig
from .continuous_field_transition import PURE_BRANCH_FACE
from .local_implicit_junction import sample_from_junction
from .surface_continuity_qc import face_region_labels
from .surface_qc import _orthogonal_basis, _section_polygon
from .surface_transition import (
    EXPLICIT_BRANCH_FACE,
    JUNCTION_CORE_FACE,
    TRANSITION_COLLAR_FACE,
)
from .types import BranchGeometry, HybridBuildDetails, PortGeometry


def _normal_statistics(
    angles: np.ndarray,
    lengths: np.ndarray | None = None,
) -> dict[str, float | int | None]:
    values = np.asarray(angles, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {
            "interface_edge_count": 0,
            "interface_length_um": 0.0,
            "dihedral_mean_deg": None,
            "dihedral_p95_deg": None,
            "dihedral_p99_deg": None,
            "dihedral_max_deg": None,
        }
    return {
        "interface_edge_count": int(len(values)),
        "interface_length_um": (
            float(np.sum(lengths)) if lengths is not None else None
        ),
        "dihedral_mean_deg": float(np.mean(values)),
        "dihedral_p95_deg": float(np.percentile(values, 95)),
        "dihedral_p99_deg": float(np.percentile(values, 99)),
        "dihedral_max_deg": float(np.max(values)),
    }


def _adjacency_arrays(
    mesh: trimesh.Trimesh,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    adjacency = np.asarray(mesh.face_adjacency, dtype=np.int64)
    edges = np.asarray(mesh.face_adjacency_edges, dtype=np.int64)
    edge_points = np.asarray(mesh.vertices, dtype=float)[edges]
    lengths = np.linalg.norm(np.diff(edge_points, axis=1)[:, 0, :], axis=1)
    angles = np.degrees(np.asarray(mesh.face_adjacency_angles, dtype=float))
    return adjacency, edges, edge_points, lengths, angles


def port_interface_rows(
    mesh: trimesh.Trimesh,
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    *,
    version: str,
) -> list[dict[str, Any]]:
    """Measure the exact CORE→PORT_EXTENSION ring, not a wide port window."""

    _, edges, edge_points, lengths, angles = _adjacency_arrays(mesh)
    rows: list[dict[str, Any]] = []
    branch_by_endpoint: dict[int, BranchGeometry] = {}
    for branch in branches:
        branch_by_endpoint[int(branch.local_node_ids[0])] = branch
        branch_by_endpoint[int(branch.local_node_ids[-1])] = branch
    for port in ports:
        radius = float(port.radius_um)
        branch = branch_by_endpoint.get(int(port.local_node_id))
        # The named interface is the immutable CUT_PORT plane.  Its normal is
        # the stored port normal; the weighted source fit belongs in the C1
        # diagnostic and must not rotate the measurement plane away from the
        # actual seam on curved upstream branches.
        tangent = np.asarray(port.outward_tangent, dtype=float)
        relative = edge_points - np.asarray(port.original_position_um)[None, None, :]
        axial = relative @ tangent
        radial_vector = relative - axial[:, :, None] * tangent[None, None, :]
        radial = np.linalg.norm(radial_vector, axis=2)
        edge_vectors = np.diff(edge_points, axis=1)[:, 0, :]
        edge_vectors /= np.maximum(
            np.linalg.norm(edge_vectors, axis=1)[:, None], 1.0e-15
        )
        # A tube triangulation also contains oblique inter-ring diagonals.  The
        # earlier 0.35 cosine admitted those diagonals and fused several nearby
        # stations into one graph component, so the reported seam was displaced
        # by roughly one centerline sample.  Circumferential ring edges are much
        # closer to orthogonal even on the curved v6 profile.
        ring_candidate = (
            np.abs(edge_vectors @ tangent) <= 0.20
        ) & (np.mean(radial, axis=1) >= 0.5 * radius) & (
            np.mean(radial, axis=1) <= 1.5 * radius
        )
        candidate_ids = np.flatnonzero(ring_candidate)
        graph = nx.Graph()
        graph.add_edges_from(map(tuple, edges[candidate_ids].tolist()))
        edge_lookup = {
            tuple(sorted(map(int, edge))): int(edge_id)
            for edge_id, edge in zip(candidate_ids, edges[candidate_ids])
        }
        components: list[tuple[float, np.ndarray]] = []
        for nodes in nx.connected_components(graph):
            node_set = set(map(int, nodes))
            component_edges = np.asarray(
                [
                    edge_lookup[tuple(sorted((int(a), int(b))))]
                    for a, b in graph.subgraph(node_set).edges
                ],
                dtype=np.int64,
            )
            if not len(component_edges):
                continue
            component_center = np.asarray(mesh.vertices)[list(node_set)].mean(axis=0)
            components.append(
                (
                    float(
                        np.linalg.norm(
                            component_center
                            - np.asarray(port.original_position_um, dtype=float)
                        )
                    ),
                    component_edges,
                )
            )
        expected = int(config.geometry.tube_sides)
        near_expected = [
            item
            for item in components
            if max(3, int(0.75 * expected))
            <= len(item[1])
            <= int(1.25 * expected)
        ]
        if near_expected:
            chosen = min(near_expected, key=lambda item: item[0])[1]
        else:
            score = np.abs(np.mean(axial, axis=1))
            chosen = candidate_ids[np.argsort(score[candidate_ids])[:expected]]
        selected = np.zeros(len(edges), dtype=bool)
        selected[chosen] = True
        chosen_midpoints = edge_points[chosen].mean(axis=1)
        section = _section_polygon(
            mesh,
            np.asarray(port.original_position_um, dtype=float),
            tangent,
        )
        area = float(section[0]) if section else None
        source_area = float(np.pi * radius**2)
        rows.append(
            {
                "version": version,
                "interface": "CORE_TO_PORT_EXTENSION_INTERFACE",
                "port_id": port.port_id,
                "cut_port_id": port.cut_port_id,
                "branch_id": branch.branch_id if branch is not None else None,
                "ring_detection_method": (
                    "circumferential_connected_component"
                    if near_expected
                    else "nearest_circumferential_edges_fallback"
                ),
                "ring_mean_center_offset_um": float(
                    np.linalg.norm(
                        chosen_midpoints.mean(axis=0)
                        - np.asarray(port.original_position_um, dtype=float)
                    )
                ),
                "ring_mean_abs_axial_offset_um": float(
                    np.mean(np.abs(np.mean(axial[chosen], axis=1)))
                ),
                "source_radius_um": radius,
                "section_area_um2": area,
                "area_relative_error": (
                    (area - source_area) / source_area if area is not None else None
                ),
                "equivalent_radius_relative_error": (
                    (float(np.sqrt(area / np.pi)) - radius) / radius
                    if area is not None
                    else None
                ),
                **_normal_statistics(angles[selected], lengths[selected]),
            }
        )
    return rows


def _interface_pairs(version: str) -> tuple[tuple[str, int, int], ...]:
    if version == "v5":
        return (
            (
                "IMPLICIT_TO_LOOP_STITCH_INTERFACE",
                int(JUNCTION_CORE_FACE),
                int(TRANSITION_COLLAR_FACE),
            ),
            (
                "LOOP_STITCH_TO_EXPLICIT_INTERFACE",
                int(TRANSITION_COLLAR_FACE),
                int(EXPLICIT_BRANCH_FACE),
            ),
        )
    return (
        (
            "CORE_TO_CONTINUOUS_TRANSITION_INTERFACE",
            int(JUNCTION_CORE_FACE),
            int(TRANSITION_COLLAR_FACE),
        ),
        (
            "CONTINUOUS_TRANSITION_TO_PURE_BRANCH_INTERFACE",
            int(TRANSITION_COLLAR_FACE),
            int(PURE_BRANCH_FACE),
        ),
        (
            "PURE_BRANCH_TO_EXPLICIT_WELD_INTERFACE",
            int(PURE_BRANCH_FACE),
            int(EXPLICIT_BRANCH_FACE),
        ),
    )


def junction_interface_rows(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    branches: list[BranchGeometry],
    *,
    version: str,
) -> tuple[list[dict[str, Any]], np.ndarray]:
    """Record edge-local dihedral statistics for every named junction interface."""

    labels = face_region_labels(mesh, details, branches)
    adjacency, edges, _, lengths, angles = _adjacency_arrays(mesh)
    rows: list[dict[str, Any]] = []
    highlighted: list[np.ndarray] = []
    for name, first, second in _interface_pairs(version):
        pair = (
            ((labels[adjacency[:, 0]] == first) & (labels[adjacency[:, 1]] == second))
            | ((labels[adjacency[:, 0]] == second) & (labels[adjacency[:, 1]] == first))
        )
        if np.any(pair):
            highlighted.append(edges[pair])
        row = {
            "version": version,
            "interface": name,
            "junction_node_id": "ALL",
            "branch_id": "ALL",
            **_normal_statistics(angles[pair], lengths[pair]),
        }
        rows.append(row)
        for node_id, patch in details.patches.items():
            for collar in patch.collars:
                branch = next(
                    item for item in branches if item.branch_id == collar.branch_id
                )
                if name.startswith("IMPLICIT") or name.startswith("CORE"):
                    distance = collar.explicit_cap_distance_um
                elif "PURE_BRANCH_TO_EXPLICIT" in name:
                    distance = collar.collar_distance_um + 0.5 * collar.overlap_length_um
                elif version == "v6":
                    distance = collar.collar_distance_um
                else:
                    distance = collar.implicit_extent_um
                center, radius, tangent = sample_from_junction(
                    branch, collar.endpoint_index, distance
                )
                edge_midpoints = np.asarray(mesh.vertices)[edges].mean(axis=1)
                relative = edge_midpoints - center[None, :]
                axial = relative @ tangent
                radial = np.linalg.norm(
                    relative - axial[:, None] * tangent[None, :], axis=1
                )
                local = pair & (np.abs(axial) <= 0.35 * radius) & (
                    radial <= 1.75 * radius
                )
                rows.append(
                    {
                        "version": version,
                        "interface": name,
                        "junction_node_id": node_id,
                        "branch_id": collar.branch_id,
                        **_normal_statistics(angles[local], lengths[local]),
                    }
                )
    return rows, (
        np.unique(np.vstack(highlighted), axis=0)
        if highlighted
        else np.empty((0, 2), dtype=np.int64)
    )


def _triangle_aspect_ratio(mesh: trimesh.Trimesh) -> np.ndarray:
    triangles = np.asarray(mesh.triangles, dtype=float)
    lengths = np.linalg.norm(
        np.stack(
            (
                triangles[:, 1] - triangles[:, 0],
                triangles[:, 2] - triangles[:, 1],
                triangles[:, 0] - triangles[:, 2],
            ),
            axis=1,
        ),
        axis=2,
    )
    maximum = np.max(lengths, axis=1)
    area = np.asarray(mesh.area_faces, dtype=float)
    return maximum**2 / np.maximum(2.0 * np.sqrt(3.0) * area, 1.0e-20)


def _face_roughness(mesh: trimesh.Trimesh) -> np.ndarray:
    vertices = np.asarray(mesh.vertices, dtype=float)
    values = np.full(len(vertices), np.nan, dtype=float)
    for vertex_id, neighbors in enumerate(mesh.vertex_neighbors):
        if not len(neighbors):
            continue
        points = vertices[np.asarray(neighbors, dtype=np.int64)]
        scale = float(np.mean(np.linalg.norm(points - vertices[vertex_id], axis=1)))
        if scale > 0.0:
            values[vertex_id] = (
                float(np.linalg.norm(vertices[vertex_id] - points.mean(axis=0))) / scale
            )
    return np.nanmean(values[np.asarray(mesh.faces, dtype=np.int64)], axis=1)


def transition_region_rows(
    mesh: trimesh.Trimesh,
    details: HybridBuildDetails,
    branches: list[BranchGeometry],
    *,
    version: str,
) -> list[dict[str, Any]]:
    labels = face_region_labels(mesh, details, branches)
    adjacency = np.asarray(mesh.face_adjacency, dtype=np.int64)
    angles = np.degrees(np.asarray(mesh.face_adjacency_angles, dtype=float))
    aspect = _triangle_aspect_ratio(mesh)
    roughness = _face_roughness(mesh)
    region_codes = (
        ("CORE", int(JUNCTION_CORE_FACE)),
        ("TRANSITION", int(TRANSITION_COLLAR_FACE)),
        ("PURE_BRANCH", int(PURE_BRANCH_FACE)),
        ("EXPLICIT", int(EXPLICIT_BRANCH_FACE)),
    )
    rows: list[dict[str, Any]] = []
    for region, code in region_codes:
        selected_faces = labels == code
        selected_edges = selected_faces[adjacency[:, 0]] & selected_faces[adjacency[:, 1]]
        selected_angles = angles[selected_edges]
        rows.append(
            {
                "version": version,
                "region": region,
                "triangle_count": int(np.count_nonzero(selected_faces)),
                "normal_jump_p95_deg": (
                    float(np.percentile(selected_angles, 95))
                    if len(selected_angles)
                    else None
                ),
                "normal_jump_p99_deg": (
                    float(np.percentile(selected_angles, 99))
                    if len(selected_angles)
                    else None
                ),
                "normal_jump_max_deg": (
                    float(np.max(selected_angles)) if len(selected_angles) else None
                ),
                "mean_roughness": (
                    float(np.nanmean(roughness[selected_faces]))
                    if np.any(selected_faces)
                    else None
                ),
                "triangle_aspect_ratio_p95": (
                    float(np.percentile(aspect[selected_faces], 95))
                    if np.any(selected_faces)
                    else None
                ),
                "triangle_aspect_ratio_max": (
                    float(np.max(aspect[selected_faces]))
                    if np.any(selected_faces)
                    else None
                ),
            }
        )
    return rows


def _view_basis(direction: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    direction = np.asarray(direction, dtype=float)
    direction /= np.linalg.norm(direction)
    return _orthogonal_basis(direction)


def silhouette_rows(
    mesh: trimesh.Trimesh,
    *,
    version: str,
    large_corner_deg: float,
    comparison_meshes: tuple[trimesh.Trimesh, ...] | None = None,
    raster_size: int = 768,
) -> list[dict[str, Any]]:
    """Measure six outer silhouettes on a fixed, shared raster.

    Directly counting mesh silhouette vertices is tessellation dependent: a
    higher-resolution continuous field can look smoother while producing more
    graph vertices.  A fixed raster and shared projection bounds make v5/v6
    curvature and corner counts comparable at the same visual scale.
    """

    directions = {
        "x": np.asarray((1.0, 0.0, 0.0)),
        "y": np.asarray((0.0, 1.0, 0.0)),
        "z": np.asarray((0.0, 0.0, 1.0)),
        "diag_xyz": np.asarray((1.0, 1.0, 1.0)),
        "diag_xmy": np.asarray((1.0, -1.0, 1.0)),
        "diag_mxz": np.asarray((-1.0, 1.0, 1.0)),
    }
    vertices = np.asarray(mesh.vertices, dtype=float)
    faces = np.asarray(mesh.faces, dtype=np.int64)
    comparison = comparison_meshes or (mesh,)
    rows: list[dict[str, Any]] = []
    for name, direction in directions.items():
        direction = direction / np.linalg.norm(direction)
        first, second = _view_basis(direction)
        projected = np.column_stack((vertices @ first, vertices @ second))
        projected_all = [
            np.column_stack(
                (
                    np.asarray(candidate.vertices, dtype=float) @ first,
                    np.asarray(candidate.vertices, dtype=float) @ second,
                )
            )
            for candidate in comparison
        ]
        stacked = np.vstack(projected_all)
        minimum = stacked.min(axis=0)
        maximum = stacked.max(axis=0)
        span = np.maximum(maximum - minimum, 1.0e-9)
        padding = 0.04 * float(np.max(span))
        minimum -= padding
        maximum += padding
        shared_span = maximum - minimum
        scale = (float(raster_size) - 1.0) / float(np.max(shared_span))
        offset = 0.5 * (
            np.asarray((raster_size - 1, raster_size - 1), dtype=float)
            - shared_span * scale
        )
        pixels = np.rint((projected - minimum) * scale + offset).astype(np.int32)
        pixels[:, 1] = raster_size - 1 - pixels[:, 1]
        mask = np.zeros((raster_size, raster_size), dtype=np.uint8)
        cv2.fillPoly(mask, pixels[faces], 255)
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        contour = (
            max(contours, key=cv2.contourArea)[:, 0, :].astype(float)
            if contours
            else np.empty((0, 2), dtype=float)
        )
        # Measure over a five-pixel physical neighbourhood to suppress pixel-grid
        # 45/90 degree quantisation without smoothing the source surface itself.
        window = 5
        if len(contour) > 2 * window:
            previous = contour - np.roll(contour, window, axis=0)
            following = np.roll(contour, -window, axis=0) - contour
            previous_length = np.linalg.norm(previous, axis=1)
            following_length = np.linalg.norm(following, axis=1)
            valid = (previous_length > 0.0) & (following_length > 0.0)
            cosine = np.sum(previous * following, axis=1) / np.maximum(
                previous_length * following_length, 1.0e-15
            )
            turning_array = np.degrees(
                np.arccos(np.clip(cosine[valid], -1.0, 1.0))
            )
            physical_length = (
                0.5 * (previous_length[valid] + following_length[valid]) / scale
            )
            curvature_array = np.radians(turning_array) / np.maximum(
                physical_length, 1.0e-15
            )
        else:
            turning_array = np.empty(0, dtype=float)
            curvature_array = np.empty(0, dtype=float)
        rows.append(
            {
                "version": version,
                "view": name,
                "measurement_method": "shared_fixed_raster_outer_contour",
                "raster_size_px": int(raster_size),
                "silhouette_edge_count": int(len(contour)),
                "ordered_vertex_count": int(len(turning_array)),
                "silhouette_curvature_variation": (
                    float(np.std(curvature_array)) if len(curvature_array) else None
                ),
                "silhouette_turning_p95_deg": (
                    float(np.percentile(turning_array, 95))
                    if len(turning_array)
                    else None
                ),
                "large_corner_count": int(
                    np.count_nonzero(turning_array >= large_corner_deg)
                ),
                "large_corner_fraction": (
                    float(np.mean(turning_array >= large_corner_deg))
                    if len(turning_array)
                    else None
                ),
            }
        )
    return rows


def pooled_interface_metric(
    rows: list[dict[str, Any]],
    version: str,
    key: str,
) -> float | None:
    values = [
        float(row[key])
        for row in rows
        if row["version"] == version
        and row["junction_node_id"] == "ALL"
        and row[key] is not None
    ]
    return max(values, default=None)
