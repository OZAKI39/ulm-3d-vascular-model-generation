"""Boundary patch labeling, watertight/manifold checks, and radius fidelity."""

from __future__ import annotations

from typing import Any

import networkx as nx
import numpy as np
import pyvista as pv
import trimesh
from shapely.geometry import Point, Polygon

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .types import BranchGeometry, PatchResult, PortGeometry, RadiusFidelitySample


def identify_port_patches(
    mesh: trimesh.Trimesh,
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    *,
    plane_tolerance_um: float | None = None,
    face_normal_alignment: float | None = None,
) -> PatchResult:
    centroids = np.asarray(mesh.triangles_center, dtype=float)
    patch_id = np.zeros(len(mesh.faces), dtype=np.int32)
    patch_type = np.zeros(len(mesh.faces), dtype=np.int8)
    face_port_id = np.full(len(mesh.faces), -1, dtype=np.int32)
    port_rows: list[dict[str, Any]] = []
    all_pass = True
    plane_tolerance = float(
        config.ports.plane_tolerance_um
        if plane_tolerance_um is None
        else plane_tolerance_um
    )
    face_alignment = float(
        config.ports.minimum_normal_alignment
        if face_normal_alignment is None
        else face_normal_alignment
    )
    for port in ports:
        relative = centroids - port.cap_center_um
        signed = relative @ port.outward_tangent
        planar = relative - signed[:, None] * port.outward_tangent[None, :]
        radial = np.linalg.norm(planar, axis=1)
        mask = (
            (np.abs(signed) <= plane_tolerance)
            & (radial <= port.radius_um * (1.0 + config.ports.radial_tolerance_fraction))
            & ((mesh.face_normals @ port.outward_tangent) >= face_alignment)
            & (face_port_id < 0)
        )
        face_indices = np.flatnonzero(mask)
        patch_id[mask] = port.port_id + 1
        patch_type[mask] = 1
        face_port_id[mask] = port.port_id
        patch_area = float(mesh.area_faces[mask].sum()) if len(face_indices) else 0.0
        expected_area = float(np.pi * port.radius_um**2)
        area_error = abs(patch_area - expected_area) / expected_area
        components = _face_component_count(mesh, face_indices)
        if len(face_indices):
            weighted_normal = np.sum(
                mesh.face_normals[mask] * mesh.area_faces[mask, None], axis=0
            )
            normal_norm = float(np.linalg.norm(weighted_normal))
            alignment = (
                float(np.dot(weighted_normal / normal_norm, port.outward_tangent))
                if normal_norm > 0
                else -1.0
            )
        else:
            alignment = -1.0
        port_pass = (
            len(face_indices) > 0
            and components == 1
            and alignment >= config.ports.minimum_normal_alignment
            and (
                not config.surface_qc.check_port_area
                or area_error <= config.ports.area_relative_tolerance
            )
        )
        all_pass &= port_pass
        row = port.metadata(patch_area_um2=patch_area)
        row.update(
            {
                "face_count": int(len(face_indices)),
                "surface_component_count": components,
                "expected_area_um2": expected_area,
                "area_relative_error": area_error,
                "normal_alignment": alignment,
                "qc_status": "PASS" if port_pass else "FAIL",
            }
        )
        port_rows.append(row)
    detected = sum(row["face_count"] > 0 for row in port_rows)
    return PatchResult(
        patch_id=patch_id,
        patch_type=patch_type,
        port_id=face_port_id,
        port_rows=port_rows,
        detected_port_count=detected,
        all_ports_pass=bool(all_pass and detected == len(ports)),
    )


def _face_component_count(mesh: trimesh.Trimesh, face_indices: np.ndarray) -> int:
    if len(face_indices) == 0:
        return 0
    if len(face_indices) == 1:
        return 1
    selected = set(map(int, face_indices))
    graph = nx.Graph()
    graph.add_nodes_from(selected)
    graph.add_edges_from(
        (int(first), int(second))
        for first, second in mesh.face_adjacency
        if int(first) in selected and int(second) in selected
    )
    return nx.number_connected_components(graph)


def _pyvista_feature_edges(mesh: trimesh.Trimesh) -> tuple[int, int]:
    faces = np.column_stack((np.full(len(mesh.faces), 3), mesh.faces)).ravel()
    polydata = pv.PolyData(np.asarray(mesh.vertices), faces)
    boundary = polydata.extract_feature_edges(
        feature_edges=False,
        manifold_edges=False,
        non_manifold_edges=False,
        boundary_edges=True,
    )
    non_manifold = polydata.extract_feature_edges(
        feature_edges=False,
        manifold_edges=False,
        non_manifold_edges=True,
        boundary_edges=False,
    )
    return int(boundary.n_cells), int(non_manifold.n_cells)


def evaluate_surface_qc(
    mesh: trimesh.Trimesh,
    patch: PatchResult,
    roi: ROIRecord,
    branches: list[BranchGeometry],
    config: CFDLumenConfig,
) -> dict[str, Any]:
    _, edge_counts = np.unique(mesh.edges_sorted, axis=0, return_counts=True)
    boundary_edge_count = int(np.count_nonzero(edge_counts == 1))
    non_manifold_edge_count = int(np.count_nonzero(edge_counts > 2))
    vtk_boundary_count, vtk_non_manifold_count = _pyvista_feature_edges(mesh)
    components = mesh.split(only_watertight=False)
    connected_count = len(components)
    expected_ports = len(roi.cut_ports)
    checks = {
        "connected": connected_count == 1,
        "watertight_trimesh": bool(mesh.is_watertight),
        "watertight_vtk": vtk_boundary_count == 0,
        "manifold": non_manifold_edge_count == 0 and vtk_non_manifold_count == 0,
        "positive_volume": bool(np.isfinite(mesh.volume) and mesh.volume > 0.0),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "port_count_matches": patch.detected_port_count == expected_ports,
        "port_patches_pass": patch.all_ports_pass,
        "source_topology_preserved": len(branches) > 0,
    }
    required = ["positive_volume", "winding_consistent", "port_count_matches", "port_patches_pass"]
    if config.surface_qc.check_connected:
        required.append("connected")
    if config.surface_qc.check_watertight:
        required.extend(("watertight_trimesh", "watertight_vtk"))
    if config.surface_qc.check_manifold:
        required.append("manifold")
    return {
        "status": "PASS" if all(checks[name] for name in required) else "FAIL",
        "checks": checks,
        "required_checks": required,
        "vertex_count": int(len(mesh.vertices)),
        "triangle_count": int(len(mesh.faces)),
        "surface_component_count": connected_count,
        "boundary_edge_count": boundary_edge_count,
        "vtk_boundary_feature_edge_count": vtk_boundary_count,
        "non_manifold_edge_count": non_manifold_edge_count,
        "vtk_non_manifold_feature_edge_count": vtk_non_manifold_count,
        "surface_area_um2": float(mesh.area),
        "enclosed_volume_um3": float(mesh.volume),
        "expected_port_count": expected_ports,
        "detected_port_patch_count": patch.detected_port_count,
        "source_branch_count": len(branches),
        "source_bifurcation_count": int(
            sum(
                degree >= 3
                for degree in np.bincount(np.asarray(roi.local_edges).ravel(), minlength=roi.node_count)
            )
        ),
        "source_cut_port_count": expected_ports,
    }


def _orthogonal_basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    normal = normal / np.linalg.norm(normal)
    helper = np.asarray((1.0, 0.0, 0.0)) if abs(normal[0]) < 0.8 else np.asarray((0.0, 1.0, 0.0))
    first = np.cross(normal, helper)
    first /= np.linalg.norm(first)
    second = np.cross(normal, first)
    return first, second


def _section_polygon(
    mesh: trimesh.Trimesh,
    center: np.ndarray,
    tangent: np.ndarray,
) -> tuple[float, tuple[tuple[float, float], ...]] | None:
    section = mesh.section(plane_origin=center, plane_normal=tangent)
    if section is None:
        return None
    basis_first, basis_second = _orthogonal_basis(tangent)
    candidates: list[tuple[bool, float, float, np.ndarray]] = []
    for discrete in section.discrete:
        points = np.asarray(discrete, dtype=float)
        if len(points) < 4:
            continue
        relative = points - center
        projected = np.column_stack((relative @ basis_first, relative @ basis_second))
        polygon = Polygon(projected)
        if not polygon.is_valid:
            polygon = polygon.buffer(0)
        if polygon.is_empty or polygon.area <= 0:
            continue
        # ``buffer(0)`` can legitimately repair a self-touching section into a
        # MultiPolygon.  Treat its connected polygons as separate candidates;
        # the section belonging to this branch is selected below by containment
        # of, then distance to, the sampled centreline point.
        polygons = polygon.geoms if polygon.geom_type == "MultiPolygon" else (polygon,)
        for component in polygons:
            if component.is_empty or component.area <= 0:
                continue
            contains_center = bool(
                component.buffer(1.0e-9).contains(Point(0.0, 0.0))
            )
            distance = float(component.distance(Point(0.0, 0.0)))
            coordinates = np.asarray(component.exterior.coords, dtype=float)
            candidates.append(
                (contains_center, distance, float(component.area), coordinates)
            )
    if not candidates:
        return None
    candidates.sort(key=lambda item: (not item[0], item[1], -item[2]))
    _, _, area, coordinates = candidates[0]
    return area, tuple((float(point[0]), float(point[1])) for point in coordinates)


def _tangent(points: np.ndarray, index: int) -> np.ndarray:
    if index == 0:
        vector = points[1] - points[0]
    elif index == len(points) - 1:
        vector = points[-1] - points[-2]
    else:
        vector = points[index + 1] - points[index - 1]
    return vector / np.linalg.norm(vector)


def evaluate_radius_fidelity(
    mesh: trimesh.Trimesh,
    branches: list[BranchGeometry],
    roi: ROIRecord,
    config: CFDLumenConfig,
) -> tuple[list[RadiusFidelitySample], dict[str, Any]]:
    degree = np.bincount(np.asarray(roi.local_edges).ravel(), minlength=roi.node_count)
    cut_ids = {int(port.local_node_id) for port in roi.cut_ports}
    samples: list[RadiusFidelitySample] = []
    attempted = 0
    missed = 0
    for branch in branches:
        length = float(branch.arc_length_um[-1])
        first_local = branch.local_node_ids[0]
        last_local = branch.local_node_ids[-1]
        start_exclusion = (
            2.0 * branch.radius_um[0] * config.surface_qc.radius_fidelity_skip_diameters
            if degree[first_local] >= 3 or first_local in cut_ids
            else 0.0
        )
        end_exclusion = (
            2.0 * branch.radius_um[-1] * config.surface_qc.radius_fidelity_skip_diameters
            if degree[last_local] >= 3 or last_local in cut_ids
            else 0.0
        )
        valid = np.flatnonzero(
            (branch.arc_length_um > start_exclusion)
            & ((length - branch.arc_length_um) > end_exclusion)
        )
        if len(valid) == 0:
            continue
        count = min(config.surface_qc.radius_fidelity_samples_per_branch, len(valid))
        chosen = np.unique(valid[np.linspace(0, len(valid) - 1, count).round().astype(int)])
        for index in chosen:
            attempted += 1
            center = branch.points_um[index]
            tangent = _tangent(branch.points_um, int(index))
            section = _section_polygon(mesh, center, tangent)
            if section is None:
                missed += 1
                continue
            area, coordinates = section
            reconstructed = float(np.sqrt(area / np.pi))
            source = float(branch.radius_um[index])
            error = (reconstructed - source) / source
            samples.append(
                RadiusFidelitySample(
                    branch_id=branch.branch_id,
                    sample_index=int(index),
                    arc_length_um=float(branch.arc_length_um[index]),
                    center_um=tuple(map(float, center)),
                    tangent=tuple(map(float, tangent)),
                    source_radius_um=source,
                    reconstructed_radius_um=reconstructed,
                    relative_error=float(error),
                    section_xy_um=coordinates,
                )
            )
    absolute = np.asarray([abs(sample.relative_error) for sample in samples], dtype=float)
    if len(absolute):
        metrics = {
            "median_absolute_relative_error": float(np.median(absolute)),
            "mean_absolute_relative_error": float(np.mean(absolute)),
            "p95_absolute_relative_error": float(np.percentile(absolute, 95)),
            "max_absolute_relative_error": float(np.max(absolute)),
        }
    else:
        metrics = {
            "median_absolute_relative_error": None,
            "mean_absolute_relative_error": None,
            "p95_absolute_relative_error": None,
            "max_absolute_relative_error": None,
        }
    branch_metrics: dict[str, dict[str, float | int]] = {}
    for branch in branches:
        values = np.asarray(
            [abs(sample.relative_error) for sample in samples if sample.branch_id == branch.branch_id],
            dtype=float,
        )
        if len(values):
            branch_metrics[str(branch.branch_id)] = {
                "sample_count": int(len(values)),
                "median_absolute_relative_error": float(np.median(values)),
                "p95_absolute_relative_error": float(np.percentile(values, 95)),
                "max_absolute_relative_error": float(values.max()),
            }
    threshold = config.surface_qc.max_radius_p95_error
    threshold_pass = (
        True
        if threshold is None
        else metrics["p95_absolute_relative_error"] is not None
        and float(metrics["p95_absolute_relative_error"]) <= threshold
    )
    report = {
        "status": "PASS" if threshold_pass else "FAIL",
        "threshold_enforced": threshold is not None,
        "max_radius_p95_error": threshold,
        "attempted_sample_count": attempted,
        "successful_sample_count": len(samples),
        "missed_section_count": missed,
        **metrics,
        "per_branch": branch_metrics,
    }
    return samples, report
