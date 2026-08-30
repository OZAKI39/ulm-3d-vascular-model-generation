"""Read-only Coarse/Base/Fine vascular-port grid-sensitivity forensics.

The module consumes the accepted Seeder meshes and continuous port patches.  It
does not launch Seeder, Musubi, or Harvester and is intentionally disconnected
from the production CFD pipeline.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.metadata
import json
import math
import re
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import trimesh
from scipy.spatial import ConvexHull

from .adaptive_flux_steady_exact_audit import _pressure_neighbor_indices
from .apes import parse_mesh_header
from .exact_link_flux import closest_discrete_direction
from .io import sha256_file, write_json
from .musubi_boundary_mass_referee import load_mesh_contract, runtime_solid_cells
from .restart_decode import D3Q19_DIRECTIONS


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


def _header_grid_geometry(mesh: Path) -> tuple[np.ndarray, float, int]:
    text = (mesh / "header.lua").read_text(encoding="utf-8")
    block = re.search(r"boundingbox\s*=\s*\{(.*?)\n\}", text, re.DOTALL)
    if not block:
        raise ValueError(f"No bounding box in {mesh / 'header.lua'}")
    origin_match = re.search(r"origin\s*=\s*\{(.*?)\}", block.group(1), re.DOTALL)
    length_match = re.search(r"length\s*=\s*([^\s]+)", block.group(1))
    level_match = re.search(r"(?m)^\s*minLevel\s*=\s*(\d+)", text)
    if not origin_match or not length_match or not level_match:
        raise ValueError("Incomplete TreElm header geometry")
    origin = np.asarray(
        [float(value.replace("D", "E")) for value in origin_match.group(1).replace(",", " ").split()],
        dtype=np.float64,
    )
    length = float(length_match.group(1).replace("D", "E"))
    level = int(level_match.group(1))
    return origin, length / (2**level), level


def _mesh_paths(root: Path) -> dict[str, Path]:
    return {
        "coarse": root / "outputs" / "cfd_flow" / GRID_RUN / "grids" / "coarse" / "seeder" / "mesh",
        "base": root / "outputs" / "cfd_flow" / BASE_MESH_RUN / "seeder" / "mesh",
        "fine": root / "outputs" / "cfd_flow" / GRID_RUN / "grids" / "fine" / "seeder" / "mesh",
    }


def _vector_distribution(indices: np.ndarray) -> dict[str, int]:
    rows = []
    for index in np.asarray(indices, dtype=np.int64):
        if int(index) < 0:
            rows.append("ZERO")
        else:
            rows.append(",".join(str(int(value)) for value in D3Q19_DIRECTIONS[int(index)]))
    return dict(sorted(Counter(rows).items()))


def _relative(left: float | None, right: float | None) -> float | None:
    if left is None or right is None or left == 0:
        return None
    return abs(float(right) - float(left)) / abs(float(left))


def discrete_port_metrics(
    mesh_path: Path,
    *,
    dx_m: float,
    frames: Mapping[str, PlaneFrame],
    continuous: Mapping[str, Mapping[str, Any]],
    cloud_output: Path | None = None,
) -> dict[str, Any]:
    mesh = load_mesh_contract(
        mesh_path,
        expected_cells=None,
        allow_zero_normals=True,
        require_runtime_order=False,
    )
    origin, header_dx, level = _header_grid_geometry(mesh_path)
    header_summary = parse_mesh_header(mesh_path)
    if not math.isclose(header_dx, dx_m, rel_tol=1e-11, abs_tol=1e-18):
        raise ValueError(f"Header dx {header_dx} != requested dx {dx_m}")
    centers = origin + (mesh.cell_ijk.astype(np.float64) + 0.5) * dx_m
    solid = runtime_solid_cells(mesh)
    result: dict[str, Any] = {
        "mesh_path": str(mesh_path.resolve()),
        "dx_m": dx_m,
        "root_level": level,
        "fluid_cell_count": len(mesh.tree_ids),
        "runtime_solid_count_global": len(solid),
        "ports": {},
    }
    cloud_payload: dict[str, np.ndarray] = {}
    for label in PORTS:
        boundary = mesh.boundaries[label]
        normal = PRESSURE_NORMALS[label]
        valid, _, _ = _pressure_neighbor_indices(
            boundary.cell_indices, mesh.cell_ijk, mesh.lookup, normal
        )
        valid &= np.fromiter(
            (int(cell) not in solid for cell in boundary.cell_indices),
            dtype=bool,
            count=len(boundary.cell_indices),
        )
        boundary_centers = centers[boundary.cell_indices]
        valid_centers = boundary_centers[valid]
        frame = frames[label]
        projected = _project(boundary_centers, frame)
        valid_projected = _project(valid_centers, frame)
        hull = _hull_metrics(projected)
        valid_hull = _hull_metrics(valid_projected)
        centroid = boundary_centers.mean(axis=0)
        centroid_error = float(np.linalg.norm(centroid - frame.origin))
        raw_sum = boundary.raw_normals.sum(axis=0).astype(np.float64)
        if np.linalg.norm(raw_sum) > 0:
            aggregate = _unit(raw_sum)
            if np.dot(aggregate, frame.normal) < 0:
                aggregate = -aggregate
            angle = math.degrees(math.acos(float(np.clip(np.dot(aggregate, frame.normal), -1, 1))))
            aggregate_list: list[float] | None = aggregate.tolist()
        else:
            aggregate_list = None
            angle = None
        signed_phase = ((boundary_centers - frame.origin) @ frame.normal) / dx_m
        valid_count = int(np.count_nonzero(valid))
        count = len(boundary.cell_indices)
        continuous_area = float(continuous[label]["continuous_area_m2"])
        proxy = count * dx_m**2
        valid_proxy = valid_count * dx_m**2
        record = {
            "seeder_boundary_cell_count_all_lattice_sides": int(
                header_summary["boundary_cell_counts"][label]
            ),
            "d3q19_globBC_count": count,
            "boundary_cell_count": count,
            "pressure_valid_cell_count": valid_count,
            "pressure_valid_fraction": valid_count / count,
            "runtime_solid_count_on_port": int(np.count_nonzero(~valid)),
            "zero_normal_count": int(np.count_nonzero(boundary.normal_indices < 0)),
            "discrete_centroid_m": centroid.tolist(),
            "centroid_error_m": centroid_error,
            "centroid_error_over_dx": centroid_error / dx_m,
            "discrete_aggregate_normal": aggregate_list,
            "normal_angle_error_deg": angle,
            "normalInd_distribution": _vector_distribution(boundary.normal_indices),
            "projected_aperture_area_proxy_m2": hull["area_m2"],
            "projected_aperture_area_over_continuous": (
                float(hull["area_m2"]) / continuous_area if hull["area_m2"] is not None else None
            ),
            "pressure_valid_projected_aperture_proxy_m2": valid_hull["area_m2"],
            "pressure_valid_projected_area_over_continuous": (
                float(valid_hull["area_m2"]) / continuous_area
                if valid_hull["area_m2"] is not None
                else None
            ),
            "boundary_count_dx2_proxy_m2": proxy,
            "boundary_count_dx2_over_continuous": proxy / continuous_area,
            "pressure_valid_count_dx2_proxy_m2": valid_proxy,
            "pressure_valid_count_dx2_over_continuous": valid_proxy / continuous_area,
            "equivalent_discrete_diameter_m": (
                2.0 * math.sqrt(float(hull["area_m2"]) / math.pi)
                if hull["area_m2"] is not None
                else None
            ),
            "perimeter_compactness_proxy": {
                "perimeter_m": hull["perimeter_m"],
                "compactness": hull["compactness"],
            },
            "true_plane_signed_cell_center_distance_over_dx": {
                "minimum": float(np.min(signed_phase)),
                "median": float(np.median(signed_phase)),
                "maximum": float(np.max(signed_phase)),
                "maximum_absolute": float(np.max(np.abs(signed_phase))),
                "median_fractional_phase": float(np.mod(np.median(signed_phase), 1.0)),
            },
            "pressure_neighbor_discrete_direction": normal.astype(int).tolist(),
            "diagnostic_resistance_proxy_L_over_A2": (
                float(continuous[label]["branch_distance"]["distance_to_nearest_bifurcation_um"])
                * 1e-6
                / max(valid_proxy, np.finfo(float).tiny) ** 2
                if continuous[label]["branch_distance"]["distance_to_nearest_bifurcation_um"] is not None
                else None
            ),
            "diagnostic_resistance_proxy_scope": "DIAGNOSTIC_ONLY_NOT_A_BIFURCATION_SOLUTION",
        }
        result["ports"][label] = record
        cloud_payload[f"{label}_boundary_centers_m"] = boundary_centers
        cloud_payload[f"{label}_pressure_valid_centers_m"] = valid_centers
    if cloud_output is not None:
        cloud_output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cloud_output, **cloud_payload)
        result["boundary_cell_center_cloud_npz"] = str(cloud_output.resolve())
    return result


def compare_grid_ports(grids: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    fields = (
        "boundary_cell_count",
        "pressure_valid_cell_count",
        "pressure_valid_fraction",
        "centroid_error_over_dx",
        "projected_aperture_area_proxy_m2",
        "pressure_valid_projected_aperture_proxy_m2",
        "boundary_count_dx2_proxy_m2",
        "pressure_valid_count_dx2_proxy_m2",
        "equivalent_discrete_diameter_m",
    )
    changes: dict[str, Any] = {}
    sensitivity_scores: dict[str, float] = {}
    correlations: dict[str, Any] = {}
    for port in PORTS:
        changes[port] = {}
        values_by_field: dict[str, list[float]] = defaultdict(list)
        for transition, left, right in (("coarse_to_base", "coarse", "base"), ("base_to_fine", "base", "fine")):
            changes[port][transition] = {}
            for field in fields:
                a = grids[left]["ports"][port][field]
                b = grids[right]["ports"][port][field]
                changes[port][transition][field + "_relative_change"] = _relative(a, b)
        for field in fields:
            values_by_field[field] = [float(grids[g]["ports"][port][field]) for g in ("coarse", "base", "fine")]
        score_fields = (
            "pressure_valid_fraction",
            "projected_aperture_area_proxy_m2",
            "pressure_valid_projected_aperture_proxy_m2",
            "equivalent_discrete_diameter_m",
        )
        score_values = [
            changes[port][transition][field + "_relative_change"] or 0.0
            for transition in ("coarse_to_base", "base_to_fine")
            for field in score_fields
        ]
        if port in OUTLETS:
            flow_bf = _relative(FLOW_SPLITS["base"][port], FLOW_SPLITS["fine"][port])
            score_values.append(float(flow_bf or 0.0))
        sensitivity_scores[port] = max(score_values)
        if port in OUTLETS:
            flow = np.asarray([FLOW_SPLITS[g][port] for g in ("coarse", "base", "fine")])
            correlations[port] = {}
            for field in (
                "pressure_valid_count_dx2_proxy_m2",
                "pressure_valid_projected_aperture_proxy_m2",
                "diagnostic_resistance_proxy_L_over_A2",
            ):
                geometry = np.asarray([grids[g]["ports"][port][field] for g in ("coarse", "base", "fine")])
                correlations[port][field] = {
                    "pearson_r": float(np.corrcoef(geometry, flow)[0, 1]),
                    "n": 3,
                    "interpretation_limit": "DIAGNOSTIC_ONLY_THREE_GRID_POINTS",
                    "flow_and_geometry_change_same_direction_base_to_fine": bool(
                        np.sign(geometry[2] - geometry[1]) == np.sign(flow[2] - flow[1])
                    ),
                }
    most_sensitive = max(OUTLETS, key=lambda label: sensitivity_scores[label])
    same_direction = sum(
        int(correlations[port]["pressure_valid_count_dx2_proxy_m2"]["flow_and_geometry_change_same_direction_base_to_fine"])
        for port in OUTLETS
    )
    pre_classification = "GEOMETRY_AND_BC" if same_direction >= 2 else "UNRESOLVED"
    return {
        "relative_changes": changes,
        "sensitivity_scores": sensitivity_scores,
        "most_grid_sensitive_port": most_sensitive,
        "geometry_flow_split_correlations": correlations,
        "same_direction_outlet_count": same_direction,
        "root_cause_classification_before_benchmark": pre_classification,
        "classification_basis": (
            "Port aperture/pressure-stencil metrics are grid-sensitive and co-vary with at least two "
            "flow splits; pressure_eq itself is also orientation- and two-neighbor-sensitive. The "
            "three-point diagnostic cannot uniquely separate geometry from BC."
            if pre_classification == "GEOMETRY_AND_BC"
            else "Existing zero-CFD evidence does not uniquely separate port, wall, bifurcation, and BC effects."
        ),
    }


def _git_sha(path: Path) -> str:
    return subprocess.run(
        ["git", "-C", str(path), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    ).stdout.strip()


def open_source_research_manifest() -> dict[str, Any]:
    """Return the pinned clean-room source audit used by this research run."""

    musubi_sha = "4e8b277b66226277171ef93bf054d21270812793"
    candidates = {
        "pressure_eq": {
            "configuration_token": "pressure_eq",
            "source_path": "mus/source/bc/mus_bc_fluid_module.fpp:788 pressure_eq; :880 pressure_eq_velocity",
            "mathematical_idea": "prescribed pressure/density plus first-order 1.5*u1-0.5*u2 velocity extrapolation, equilibrium replacement on missing incoming links",
            "D3Q19_compatible": True,
            "pull_fetch_semantics": "globBC incoming-link writes use FETCH; PULL storage is inverse because the kernel performs missing-neighbor bounce-back",
            "curved_false_requirement": "not source-enforced; numerical plane still must supply a valid constant discrete normal",
            "neighbor_requirements": 2,
            "uses_q_value": False,
            "straight_boundary_requirement": "not declared, but two neighbors are selected along one normalInd",
            "oblique_plane_suitability": "conditional on a stable D3Q19 normalInd and two interior neighbors",
            "pressure_dirichlet": True,
            "example_test": "mus/examples/fluid/benchmark/Channel2D/C2D_BoundaryConditions/C2D_BC_MfrEq_PressEq",
            "citation": "Izquierdo, Martinez-Lera & Fueyo (2009), Eq. 6",
        },
        "pressure_nonEqExpol": {
            "configuration_token": "pressure_noneq_expol",
            "source_path": "mus/source/bc/mus_bc_fluid_nonEqExpol_module.fpp:454 pressure_nonEqExpol",
            "mathematical_idea": "prescribed pressure equilibrium plus the first interior neighbor's non-equilibrium distribution",
            "D3Q19_compatible": True,
            "pull_fetch_semantics": "updates all local directions through FETCH from one pre-collision overnext-neighbor buffer",
            "curved_false_requirement": "source header calls it a straight-boundary method and fixes qVal=0",
            "neighbor_requirements": 1,
            "uses_q_value": False,
            "straight_boundary_requirement": True,
            "oblique_plane_suitability": "only when the rasterized boundary has one coherent normal direction",
            "pressure_dirichlet": True,
            "example_test": "mus/examples/fluid/benchmark/Channel2D/C2D_Simple/C2D_Simple_rBGK",
            "citation": "Guo, Zheng & Shi (2002), DOI 10.1063/1.1471914",
        },
        "pressure_antiBounceBack": {
            "configuration_token": "pressure_antibounceback",
            "source_path": "mus/source/bc/mus_bc_fluid_module.fpp:2356 pressure_antiBounceBack",
            "mathematical_idea": "anti-bounce-back Dirichlet pressure term with non-equilibrium correction and extrapolated boundary velocity",
            "D3Q19_compatible": True,
            "pull_fetch_semantics": "writes missing incoming links with FETCH and reads one post-collision neighbor buffer",
            "curved_false_requirement": "no fractional boundary distance is read",
            "neighbor_requirements": 1,
            "uses_q_value": False,
            "straight_boundary_requirement": "assumes a consistent boundary normal/neighbor direction",
            "oblique_plane_suitability": "source-implemented but must be benchmarked for orientation sensitivity",
            "pressure_dirichlet": True,
            "example_test": "documented in mus/examples/tutorials/tut_04_boundaries.md",
            "citation": "Musubi native implementation; no copied code",
        },
        "press_neq": {
            "source_path": "mus/source/bc/mus_bc_fluid_module.fpp:440 press_neq",
            "mathematical_idea": "prescribed pressure equilibrium plus locally copied non-equilibrium component scaled by 1-omega",
            "D3Q19_compatible": True,
            "pull_fetch_semantics": "updates only incoming globBC links through FETCH",
            "curved_false_requirement": "no q-value or geometric intersection distance is used",
            "neighbor_requirements": 1,
            "uses_q_value": False,
            "straight_boundary_requirement": "not explicitly declared",
            "oblique_plane_suitability": "link-local but not selected for the bounded benchmark candidate set",
            "pressure_dirichlet": True,
            "example_test": "no dedicated pinned example found",
            "citation": "Guo & Shi non-equilibrium pressure treatment as documented in source",
        },
        "pressure_expol": {
            "source_path": "mus/source/bc/mus_bc_fluid_module.fpp:1360 pressure_expol",
            "mathematical_idea": "pressure and velocity extrapolation using two interior neighbors",
            "D3Q19_compatible": True,
            "pull_fetch_semantics": "two pre-collision neighbor buffers; qVal forcibly zero",
            "curved_false_requirement": "qVal=0 in boundary loader",
            "neighbor_requirements": 2,
            "uses_q_value": False,
            "straight_boundary_requirement": True,
            "oblique_plane_suitability": "conditional on a coherent discrete normal",
            "pressure_dirichlet": True,
            "example_test": "mus/examples/fluid_incompressible/benchmark/Pipe/PIP_Simple",
            "citation": "Musubi native extrapolation implementation",
        },
    }
    return {
        "status": "PASS",
        "clean_room_policy": "ALGORITHM_AND_GEOMETRY_CONTRACTS_ONLY_NO_EXTERNAL_CODE_COPIED",
        "repositories": [
            {
                "repository": "Musubi",
                "url": "https://github.com/apes-suite/musubi",
                "actual_commit_sha": musubi_sha,
                "embedded_mus_source_commit_sha": "81f8c4f13772f6d4af31f335e1e3f99b02726e25",
                "working_tree_state": "embedded mus source contains the already-validated isolated adaptive_flux_pressure patch; it was inspected and not modified by this research task",
                "validated_binary_sha256": "e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588",
                "license": "BSD-2-Clause",
                "source_paths": ["mus/source/bc/mus_bc_fluid_module.fpp", "mus/source/bc/mus_bc_fluid_nonEqExpol_module.fpp", "mus/source/bc/mus_bc_header_module.fpp"],
                "relevant_function": list(candidates),
                "relevant_method": "native D3Q19 pressure boundary implementations and their neighbor/q-value contracts",
                "what_we_reuse": "the already-built pinned implementation and source-proven configuration names",
                "what_we_do_not_copy": "no external source code is copied",
                "relevance": "direct solver implementation under test",
            },
            {
                "repository": "HemePure",
                "url": "https://github.com/UCL-CCS/HemePure",
                "actual_commit_sha": "0bf67b16b23b41a06507810337a445f9916d62bf",
                "license": "BSD-3-Clause",
                "source_paths": ["cases/Benchmarking", "README.md"],
                "relevant_function": "image-derived vascular case representation",
                "relevant_method": "vascular geometry and pressure/flow case reference",
                "what_we_reuse": "research context only",
                "what_we_do_not_copy": "no solver or case source copied",
                "relevance": "geometry representation reference, not a second production solver",
            },
            {
                "repository": "HemeLB",
                "url": "https://github.com/hemelb-codes/hemelb",
                "actual_commit_sha": "432d3386e32571487a1521f88a407db5dc9ca171",
                "license": "LGPL-3.0",
                "source_paths": ["Code/lb/streamers/NashZerothOrderPressure.h", "Code/lb/streamers/BouzidiFirdaousLallemand.h", "Code/lb/streamers/GuoZhengShi.h", "Code/configuration/SimConfig.cc"],
                "relevant_function": "NashZerothOrderPressureLink, BouzidiFirdaousLallemandLink, GuoZhengShiLink, DoIOForBaseInOutlet",
                "relevant_method": "ghost-site pressure equilibrium; fractional-link BFL/GZS wall treatments; explicit iolet position and normal",
                "what_we_reuse": "method-level comparison and plane/normal representation concepts",
                "what_we_do_not_copy": "LGPL implementation code is not copied or linked",
                "relevance": "distinguishes pressure-iolet semantics from fractional wall distance",
            },
            {
                "repository": "VMTK",
                "url": "https://github.com/vmtk/vmtk",
                "actual_commit_sha": "06c7fb60f8bb873718d5f8b54af0ded4841dae50",
                "license": "BSD (repository LICENSE; bundled components separately licensed)",
                "source_paths": ["vtkVmtk/ComputationalGeometry/vtkvmtkBoundaryReferenceSystems.cxx", "vtkVmtk/ComputationalGeometry/vtkvmtkCenterlineEndpointExtractor.cxx", "vmtkScripts/vmtkflowextensions.py", "vtkVmtk/ComputationalGeometry/vtkvmtkPolyDataFlowExtensionsFilter.cxx"],
                "relevant_function": "boundary barycenter/normal/radius, centerline endpoint extraction, centerline-direction or boundary-normal extensions",
                "relevant_method": "standardized port reference systems and optional CFD flow extensions",
                "what_we_reuse": "geometric concepts; current audit implements an independent VTK/PCA-equivalent plane contract",
                "what_we_do_not_copy": "no VMTK source code copied and no extension inserted",
                "relevance": "defines clean physical ports without assuming an extension is required",
            },
        ],
        "published_work": [
            {"work": "Feiger et al., Suitability of lattice Boltzmann inlet and outlet boundary conditions for image-derived vasculature", "doi": "10.1002/cnm.3198", "finding_used": "pressure BC accuracy/stability must be tested in image-derived vascular geometry"},
            {"work": "Nash et al., Choice of boundary condition for lattice-Boltzmann simulation of moderate-Re flow in complex domains", "doi": "10.1103/PhysRevE.89.023303", "finding_used": "BFL/GZS wall accuracy and orientation effects motivate a controlled rotated-pipe benchmark"},
            {"work": "Bouzidi, Firdaouss & Lallemand, Momentum transfer of a Boltzmann-lattice fluid with boundaries", "doi": "10.1063/1.1399290", "finding_used": "fractional link distance belongs to interpolated wall bounce-back"},
            {"work": "Guo, Zheng & Shi, An extrapolation method for boundary conditions in lattice Boltzmann method", "doi": "10.1063/1.1471914", "finding_used": "source-proven non-equilibrium extrapolation alternative"},
        ],
        "musubi_pressure_bc_candidates_found": candidates,
        "bounded_benchmark_set": ["pressure_eq", "pressure_nonEqExpol", "pressure_antiBounceBack"],
        "selection_reason": "CURRENT plus two pressure-Dirichlet, D3Q19-compatible, source-proven alternatives; pressure_nonEqExpol is prioritized by the task and pressure_antiBounceBack offers a distinct one-neighbor reconstruction.",
    }


def dependency_manifest() -> dict[str, Any]:
    packages = ("numpy", "scipy", "trimesh", "pyvista", "vtk", "shapely", "networkx")
    records = []
    purposes = {
        "numpy": "array geometry and mesh metrics",
        "scipy": "convex hulls and correlations",
        "trimesh": "continuous STL port recovery",
        "pyvista": "available VTK-facing inspection support",
        "vtk": "available centerline/surface inspection support",
        "shapely": "available planar geometry validation support",
        "networkx": "available graph inspection support",
    }
    for package in packages:
        try:
            version = importlib.metadata.version(package)
            state = "PRE_EXISTING_NO_INSTALL"
        except importlib.metadata.PackageNotFoundError:
            version = None
            state = "NOT_INSTALLED_NOT_REQUIRED"
        records.append({"name": package, "version": version, "install_command": state, "purpose": purposes[package]})
    return {"status": "PASS", "packages": records, "new_packages_installed": []}


def run_port_grid_sensitivity(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / RESEARCH_RUN
    qc = run_root / "qc"
    clouds = qc / "port_clouds"
    qc.mkdir(parents=True, exist_ok=True)
    clouds.mkdir(parents=True, exist_ok=True)
    frames, continuous = recover_continuous_ports(root)
    grid_records: dict[str, Any] = {}
    for label, mesh_path in _mesh_paths(root).items():
        grid_records[label] = discrete_port_metrics(
            mesh_path,
            dx_m=float(GRID_SPECS[label]["dx_m"]),
            frames=frames,
            continuous=continuous,
            cloud_output=clouds / f"{label}_port_cell_centers.npz",
        )
    comparison = compare_grid_ports(grid_records)
    result = {
        "status": "PASS",
        "method": "READ_ONLY_CONTINUOUS_PORT_PLUS_EXISTING_TREELM_D3Q19_FORENSICS",
        "seeder_calls": 0,
        "musubi_calls": 0,
        "harvester_calls": 0,
        "extra_fine_vascular_cfd_calls": 0,
        "continuous_ports": continuous,
        "grids": grid_records,
        "flow_splits": FLOW_SPLITS,
        "comparison": comparison,
    }
    write_json(qc / "open_source_boundary_research.json", open_source_research_manifest())
    write_json(qc / "python_dependency_manifest.json", dependency_manifest())
    write_json(qc / "port_grid_sensitivity_forensics.json", result)
    return result
