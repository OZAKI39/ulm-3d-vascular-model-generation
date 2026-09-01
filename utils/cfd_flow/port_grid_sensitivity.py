"""Read-only Coarse/Base/Fine vascular-port grid-sensitivity forensics.

The module consumes the accepted Seeder meshes and continuous port patches.  It
does not launch Seeder, Musubi, or Harvester and is intentionally disconnected
from the production CFD pipeline.
"""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import trimesh
from scipy.spatial import ConvexHull

from .io import sha256_file


RESEARCH_RUN = (
    "healthy_mouse_capillary_port_grid_sensitivity_research_"
    "anchor003274_20260830"
)
GRID_RUN = "healthy_mouse_capillary_grid_convergence_anchor003274_20260829"
BASE_MESH_RUN = "axis_aligned_ideal_plane_inlet_preflight_anchor003274_20260829_120444"
AXIS_GEOMETRY_RUN = "axis_aligned_inlet_geometry_anchor003274_20260829_111451"
PORT_CLASSIFICATION_RUN = "global_to_roi_anchor003274_20260825_183628"
MODEL_RUN = "ultraliser_anchor003274_20260825_133350"

GRID_SPECS = {
    "coarse": {"dx_m": 2.6e-7, "cells": 106_350},
    "base": {"dx_m": 2.0e-7, "cells": 221_309},
    "fine": {"dx_m": 1.5384615384615385e-7, "cells": 466_023},
}
PORTS = ("inlet", "outlet_01", "outlet_02", "outlet_03")
OUTLETS = PORTS[1:]
PRESSURE_NORMALS = {
    "inlet": np.asarray((0, 0, -1), dtype=np.int64),
    "outlet_01": np.asarray((0, -1, -1), dtype=np.int64),
    "outlet_02": np.asarray((1, 0, 1), dtype=np.int64),
    "outlet_03": np.asarray((-1, 0, 0), dtype=np.int64),
}
FLOW_SPLITS = {
    "coarse": {"outlet_01": 0.0981700, "outlet_02": 0.3455414, "outlet_03": 0.5422822},
    "base": {"outlet_01": 0.0851009, "outlet_02": 0.4720324, "outlet_03": 0.4400351},
    "fine": {"outlet_01": 0.0753790, "outlet_02": 0.5525517, "outlet_03": 0.3700206},
}
PORT_ROW_MAP = {
    "inlet": "cut_000",
    "outlet_01": "cut_001",
    "outlet_02": "cut_002",
    "outlet_03": "terminal_000",
}


@dataclass(frozen=True, slots=True)
class PlaneFrame:
    origin: np.ndarray
    normal: np.ndarray
    basis_u: np.ndarray
    basis_v: np.ndarray


def _unit(vector: Iterable[float]) -> np.ndarray:
    value = np.asarray(tuple(vector), dtype=np.float64).reshape(3)
    length = float(np.linalg.norm(value))
    if not np.isfinite(length) or length <= 0.0:
        raise ValueError("A finite non-zero vector is required")
    return value / length


def orthonormal_plane_basis(normal: Iterable[float]) -> tuple[np.ndarray, np.ndarray]:
    n = _unit(normal)
    helper = np.asarray((1.0, 0.0, 0.0))
    if abs(float(np.dot(n, helper))) > 0.85:
        helper = np.asarray((0.0, 1.0, 0.0))
    u = _unit(np.cross(n, helper))
    v = _unit(np.cross(n, u))
    return u, v


def _project(points: np.ndarray, frame: PlaneFrame) -> np.ndarray:
    delta = np.asarray(points, dtype=np.float64) - frame.origin
    return np.column_stack((delta @ frame.basis_u, delta @ frame.basis_v))


def _hull_metrics(points_2d: np.ndarray) -> dict[str, float | None]:
    points = np.asarray(points_2d, dtype=np.float64).reshape((-1, 2))
    if len(points) < 3 or np.linalg.matrix_rank(points - points.mean(axis=0)) < 2:
        return {"area_m2": None, "perimeter_m": None, "compactness": None}
    hull = ConvexHull(points)
    area = float(hull.volume)
    perimeter = float(hull.area)
    compactness = 4.0 * math.pi * area / perimeter**2 if perimeter > 0 else None
    return {"area_m2": area, "perimeter_m": perimeter, "compactness": compactness}


def _mesh_boundary_perimeter(mesh: trimesh.Trimesh) -> float:
    edges = np.sort(np.asarray(mesh.edges, dtype=np.int64), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    boundary = unique[counts == 1]
    if len(boundary) == 0:
        return 0.0
    vertices = np.asarray(mesh.vertices, dtype=np.float64)
    return float(np.linalg.norm(vertices[boundary[:, 0]] - vertices[boundary[:, 1]], axis=1).sum())


def _area_weighted_geometry(path: Path, expected_normal: np.ndarray) -> tuple[PlaneFrame, dict[str, Any]]:
    # STL repeats vertices per triangle.  Processing merges only coincident
    # vertices, which is required to identify the physical outer seam rather
    # than counting every internal triangle edge as perimeter.
    loaded = trimesh.load_mesh(path, process=True)
    if not isinstance(loaded, trimesh.Trimesh):
        raise ValueError(f"Expected one triangle mesh: {path}")
    triangles = np.asarray(loaded.triangles, dtype=np.float64)
    cross = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    double_area = np.linalg.norm(cross, axis=1)
    area = float(0.5 * double_area.sum())
    if area <= 0.0:
        raise ValueError(f"Zero-area port patch: {path}")
    centroid = np.average(triangles.mean(axis=1), axis=0, weights=0.5 * double_area)
    normal = _unit(cross.sum(axis=0))
    reference = _unit(expected_normal)
    if float(np.dot(normal, reference)) < 0.0:
        normal = -normal
    u, v = orthonormal_plane_basis(normal)
    frame = PlaneFrame(centroid, normal, u, v)
    vertices = np.asarray(loaded.vertices, dtype=np.float64)
    projected = _project(vertices, frame)
    hull = _hull_metrics(projected)
    covariance = np.cov(projected.T)
    eigenvalues = np.sort(np.linalg.eigvalsh(covariance))[::-1]
    axis_ratio = float(math.sqrt(max(eigenvalues[0], 0.0) / max(eigenvalues[1], np.finfo(float).tiny)))
    perimeter = _mesh_boundary_perimeter(loaded)
    diameter_h = 4.0 * area / perimeter if perimeter > 0 else None
    signed = (vertices - centroid) @ normal
    return frame, {
        "source_path": str(path.resolve()),
        "source_sha256": sha256_file(path),
        "physical_centroid_m": centroid.tolist(),
        "physical_normal": normal.tolist(),
        "continuous_area_m2": area,
        "equivalent_radius_m": math.sqrt(area / math.pi),
        "hydraulic_diameter_m": diameter_h,
        "perimeter_m": perimeter,
        "planarity_max_abs_m": float(np.max(np.abs(signed))),
        "local_cross_section_shape": {
            "convex_hull_area_m2": hull["area_m2"],
            "convex_hull_perimeter_m": hull["perimeter_m"],
            "compactness": hull["compactness"],
            "principal_axis_ratio": axis_ratio,
        },
        "plane_basis_u": u.tolist(),
        "plane_basis_v": v.tolist(),
    }


def _read_port_rows(root: Path) -> dict[str, dict[str, str]]:
    path = root / "outputs" / "cfd_preprocess" / PORT_CLASSIFICATION_RUN / "roi" / "port_classification.csv"
    rows = list(csv.DictReader(path.read_text(encoding="utf-8-sig").splitlines()))
    result: dict[str, dict[str, str]] = {}
    for label, suffix in PORT_ROW_MAP.items():
        result[label] = next(row for row in rows if row["port_id"].endswith("__" + suffix))
    return result


def _read_transform(root: Path) -> np.ndarray:
    path = root / "outputs" / "cfd_flow" / AXIS_GEOMETRY_RUN / "transform" / "anatomical_to_cfd_transform.json"
    return np.asarray(json.loads(path.read_text(encoding="utf-8"))["rotation_matrix_3x3"], dtype=np.float64)


def _swc_distance_to_bifurcation(root: Path, position_um: np.ndarray) -> dict[str, Any]:
    path = root / "outputs" / "model_generate" / MODEL_RUN / "input" / "roi_core.swc"
    nodes: dict[int, np.ndarray] = {}
    adjacency: dict[int, list[int]] = defaultdict(list)
    for raw in path.read_text(encoding="utf-8").splitlines():
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        fields = raw.split()
        node = int(fields[0])
        nodes[node] = np.asarray(tuple(float(value) for value in fields[2:5]), dtype=np.float64)
        parent = int(fields[6])
        if parent >= 0:
            adjacency[node].append(parent)
            adjacency[parent].append(node)
    nearest = min(nodes, key=lambda node: float(np.linalg.norm(nodes[node] - position_um)))
    queue: list[tuple[int, int, float]] = [(nearest, -1, 0.0)]
    visited = {nearest}
    while queue:
        node, prior, distance = queue.pop(0)
        if node != nearest and len(adjacency[node]) >= 3:
            return {
                "method": "NEAREST_SWC_NODE_GRAPH_PATH_TO_FIRST_DEGREE_GE_3",
                "nearest_swc_node": nearest,
                "nearest_node_error_um": float(np.linalg.norm(nodes[nearest] - position_um)),
                "nearest_bifurcation_node": node,
                "distance_to_nearest_bifurcation_um": distance,
            }
        for neighbor in adjacency[node]:
            if neighbor == prior or neighbor in visited:
                continue
            visited.add(neighbor)
            segment = float(np.linalg.norm(nodes[node] - nodes[neighbor]))
            queue.append((neighbor, node, distance + segment))
    return {
        "method": "NEAREST_SWC_NODE_GRAPH_PATH_TO_FIRST_DEGREE_GE_3",
        "nearest_swc_node": nearest,
        "nearest_node_error_um": float(np.linalg.norm(nodes[nearest] - position_um)),
        "nearest_bifurcation_node": None,
        "distance_to_nearest_bifurcation_um": None,
    }


def recover_continuous_ports(project_root: Path) -> tuple[dict[str, PlaneFrame], dict[str, Any]]:
    root = Path(project_root).resolve()
    rows = _read_port_rows(root)
    rotation = _read_transform(root)
    geometry = root / "outputs" / "cfd_flow" / AXIS_GEOMETRY_RUN / "geometry" / "geometry_solver_m"
    frames: dict[str, PlaneFrame] = {}
    records: dict[str, Any] = {}
    for label in PORTS:
        row = rows[label]
        source_normal = np.asarray(
            tuple(float(row[f"outward_normal_{axis}"]) for axis in "xyz"), dtype=np.float64
        )
        expected = rotation @ source_normal
        frame, record = _area_weighted_geometry(geometry / f"{label}.stl", expected)
        tangent = rotation @ np.asarray(
            tuple(float(row[f"simulation_tangent_{axis}"]) for axis in "xyz"), dtype=np.float64
        )
        record["port_id"] = row["port_id"]
        record["boundary_origin"] = row["boundary_origin"]
        record["local_branch_direction"] = _unit(tangent).tolist()
        record["outward_normal_from_source_then_rigid_transform"] = _unit(expected).tolist()
        record["cap_normal_alignment_abs_dot"] = abs(float(np.dot(frame.normal, _unit(expected))))
        record["branch_distance"] = _swc_distance_to_bifurcation(
            root,
            np.asarray(tuple(float(row[f"{axis}_um"]) for axis in "xyz"), dtype=np.float64),
        )
        records[label] = record
        frames[label] = frame
    return frames, records
