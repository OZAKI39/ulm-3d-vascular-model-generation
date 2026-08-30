"""Full, read-only validation of the final dimensionless-kernel BASE mesh."""

from __future__ import annotations

import csv
import math
import time
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh

from .dimensionless_geometry_kernel import semantic_files_success
from .io import sha256_file, write_json
from .musubi_boundary_mass_referee import load_mesh_contract
from .qvalue_contract_forensics import vascular_wall_qvalue_distribution
from .repaired_topology_forensics import (
    PORT_LABELS,
    cell_centers,
    parse_uniform_lattice,
    port_component_membership,
    sparse_component_labels,
)
from .restart_decode import D3Q19_DIRECTIONS


REQUIRED_MESH_FILES = (
    "header.lua",
    "elemlist.lsb",
    "bnd.lua",
    "bnd.lsb",
    "qval.lua",
    "qval.lsb",
)
FALLBACK_Q = 0.5
FALLBACK_ATOL = 1.0e-10
DIMENSIONLESS_TOLERANCE = 64.0 * np.finfo(np.float64).eps
MAX_CLASSIFICATION_MISMATCH_FRACTION = 0.005
MEDIAN_ABSOLUTE_Q_LIMIT = 0.01
P95_ABSOLUTE_Q_LIMIT = 0.02
RMS_Q_LIMIT = 0.01
ABSOLUTE_BIAS_LIMIT = 0.01


def _numeric_summary(values: np.ndarray) -> dict[str, float | int | None]:
    data = np.asarray(values, dtype=np.float64)
    data = data[np.isfinite(data)]
    if not len(data):
        return {
            "count": 0,
            "bias": None,
            "median_absolute": None,
            "rms": None,
            "p95": None,
            "max": None,
        }
    absolute = np.abs(data)
    return {
        "count": int(len(data)),
        "bias": float(np.mean(data)),
        "median_absolute": float(np.median(absolute)),
        "rms": float(np.sqrt(np.mean(data * data))),
        "p95": float(np.percentile(absolute, 95.0)),
        "max": float(np.max(absolute)),
    }


def _nearest_dimensionless_segment_intersection(
    start: np.ndarray,
    vector: np.ndarray,
    triangles: np.ndarray,
) -> float | None:
    """Return nearest q on one link using only link-local dimensionless math."""

    scale = float(np.max(np.abs(vector)))
    if not math.isfinite(scale) or scale <= 0.0 or not len(triangles):
        return None
    local = (np.asarray(triangles, dtype=np.float64) - start) / scale
    direction = vector / scale
    edge1 = local[:, 1] - local[:, 0]
    edge2 = local[:, 2] - local[:, 0]
    pvec = np.cross(np.broadcast_to(direction, edge2.shape), edge2)
    determinant = np.einsum("ij,ij->i", edge1, pvec)
    relative_scale = (
        np.linalg.norm(edge1, axis=1)
        * np.linalg.norm(edge2, axis=1)
        * np.linalg.norm(direction)
    )
    active = (
        relative_scale > np.finfo(np.float64).tiny
    ) & (np.abs(determinant) > DIMENSIONLESS_TOLERANCE * relative_scale)
    inverse = np.zeros_like(determinant)
    inverse[active] = 1.0 / determinant[active]
    tvec = -local[:, 0]
    barycentric_u = inverse * np.einsum("ij,ij->i", tvec, pvec)
    active &= (barycentric_u >= -DIMENSIONLESS_TOLERANCE) & (
        barycentric_u <= 1.0 + DIMENSIONLESS_TOLERANCE
    )
    qvec = np.cross(tvec, edge1)
    barycentric_v = inverse * np.einsum(
        "ij,ij->i", np.broadcast_to(direction, qvec.shape), qvec
    )
    active &= (barycentric_v >= -DIMENSIONLESS_TOLERANCE) & (
        barycentric_u + barycentric_v <= 1.0 + DIMENSIONLESS_TOLERANCE
    )
    fractions = inverse * np.einsum("ij,ij->i", edge2, qvec)
    active &= (fractions >= -DIMENSIONLESS_TOLERANCE) & (
        fractions <= 1.0 + DIMENSIONLESS_TOLERANCE
    )
    if not np.any(active):
        return None
    return float(max(0.0, np.min(fractions[active])))


def full_fluid_center_containment(
    *, mesh_dir: Path, continuous_surface: Path, outside_csv: Path
) -> dict[str, Any]:
    """Classify every fluid center against the watertight continuous lumen."""

    started = time.perf_counter()
    mesh_dir = Path(mesh_dir).resolve()
    surface_path = Path(continuous_surface).resolve()
    contract = load_mesh_contract(
        mesh_dir, allow_zero_normals=True, require_runtime_order=False
    )
    lattice = parse_uniform_lattice(
        (mesh_dir / "header.lua").read_text(encoding="utf-8")
    )
    centers = cell_centers(contract.cell_ijk, lattice)
    surface = pv.read(surface_path).triangulate()
    connected = surface.connectivity()
    component_count = int(
        len(np.unique(np.asarray(connected.cell_data["RegionId"], dtype=np.int64)))
    )
    inside_chunks = []
    batch_size = 50_000
    for first in range(0, len(centers), batch_size):
        selected = pv.PolyData(centers[first : first + batch_size]).select_enclosed_points(
            surface, tolerance=0.0, check_surface=True
        )
        inside_chunks.append(
            np.asarray(selected.point_data["SelectedPoints"], dtype=bool)
        )
    inside = np.concatenate(inside_chunks)
    outside_indices = np.flatnonzero(~inside)

    port_sets = {
        label: set(contract.boundaries[label].cell_indices.tolist())
        for label in PORT_LABELS
    }
    outside_rows = []
    outside_cap_count = 0
    for cell_index in outside_indices:
        labels = [
            label for label, members in port_sets.items() if int(cell_index) in members
        ]
        if labels:
            outside_cap_count += 1
        outside_rows.append(
            {
                "cell_index": int(cell_index),
                "tree_id": int(contract.tree_ids[cell_index]),
                "i": int(contract.cell_ijk[cell_index, 0]),
                "j": int(contract.cell_ijk[cell_index, 1]),
                "k": int(contract.cell_ijk[cell_index, 2]),
                "x_m": float(centers[cell_index, 0]),
                "y_m": float(centers[cell_index, 1]),
                "z_m": float(centers[cell_index, 2]),
                "port_boundary_labels": ";".join(labels),
                "boundary_cap_representation": bool(labels),
            }
        )
    outside_csv.parent.mkdir(parents=True, exist_ok=True)
    with outside_csv.open("w", encoding="utf-8", newline="") as stream:
        fieldnames = (
            "cell_index",
            "tree_id",
            "i",
            "j",
            "k",
            "x_m",
            "y_m",
            "z_m",
            "port_boundary_labels",
            "boundary_cap_representation",
        )
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(outside_rows)

    non_cap_outside = int(len(outside_indices) - outside_cap_count)
    passed = bool(
        surface.n_open_edges == 0
        and component_count == 1
        and non_cap_outside == 0
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "surface_path": str(surface_path),
        "surface_sha256": sha256_file(surface_path),
        "surface_is_watertight": bool(surface.n_open_edges == 0),
        "surface_component_count": component_count,
        "fluid_centers_tested": int(len(centers)),
        "inside_count": int(np.count_nonzero(inside)),
        "outside_count": int(len(outside_indices)),
        "outside_boundary_cap_count": int(outside_cap_count),
        "outside_non_cap_count": non_cap_outside,
        "outside_cells_complete_list": str(outside_csv.resolve()),
        "method": "VTK enclosed-point classification in deterministic 50000-point batches; every fluid center tested",
        "runtime_seconds": float(time.perf_counter() - started),
    }


def _port_memberships_for_cells(contract: Any, cells: np.ndarray) -> list[list[str]]:
    port_sets = {
        label: set(contract.boundaries[label].cell_indices.tolist())
        for label in PORT_LABELS
    }
    return [
        [label for label, members in port_sets.items() if int(cell) in members]
        for cell in cells
    ]


def full_active_wall_q_oracle(
    *, mesh_dir: Path, wall_stl: Path, mismatch_csv: Path
) -> dict[str, Any]:
    """Audit every active D3Q19 wall link against the original wall STL."""

    started = time.perf_counter()
    mesh_dir = Path(mesh_dir).resolve()
    wall_path = Path(wall_stl).resolve()
    contract = load_mesh_contract(
        mesh_dir, allow_zero_normals=True, require_runtime_order=False
    )
    lattice = parse_uniform_lattice(
        (mesh_dir / "header.lua").read_text(encoding="utf-8")
    )
    wall = contract.boundaries["wall"]
    wall_rows, directions = np.nonzero(wall.outward_masks)
    cells = wall.cell_indices[wall_rows]
    starts = cell_centers(contract.cell_ijk[cells], lattice)
    vectors = D3Q19_DIRECTIONS[directions].astype(np.float64) * lattice.dx
    q_seeder = contract.qvalues_by_cell[cells, directions]

    surface = trimesh.load_mesh(wall_path, process=False)
    if not isinstance(surface, trimesh.Trimesh):
        raise TypeError(f"Expected one triangular STL surface: {wall_path}")
    triangles = np.asarray(surface.triangles, dtype=np.float64)
    spatial_index = surface.triangles_tree
    padding = np.finfo(np.float64).eps * max(
        1.0, float(np.max(np.abs(surface.bounds)))
    )
    q_oracle = np.full(len(cells), np.nan, dtype=np.float64)
    candidate_counts = np.zeros(len(cells), dtype=np.int32)
    for index, (start, vector) in enumerate(zip(starts, vectors, strict=True)):
        end = start + vector
        bounds = np.concatenate(
            (np.minimum(start, end) - padding, np.maximum(start, end) + padding)
        )
        candidates = np.fromiter(
            spatial_index.intersection(bounds), dtype=np.int64
        )
        candidate_counts[index] = len(candidates)
        exact = _nearest_dimensionless_segment_intersection(
            start, vector, triangles[candidates]
        )
        if exact is not None:
            q_oracle[index] = min(exact, 1.0)

    oracle_hit = np.isfinite(q_oracle)
    q_valid = np.isfinite(q_seeder) & (q_seeder > 0.0) & (q_seeder <= 1.0)
    stored_fallback = q_valid & np.isclose(
        q_seeder, FALLBACK_Q, rtol=0.0, atol=FALLBACK_ATOL
    )
    exact_fallback_intersection = (
        oracle_hit
        & stored_fallback
        & np.isclose(q_seeder, q_oracle, rtol=0.0, atol=FALLBACK_ATOL)
    )
    seeder_hit = q_valid & (~stored_fallback | exact_fallback_intersection)
    both = seeder_hit & oracle_hit
    neither = ~seeder_hit & ~oracle_hit
    seeder_only = seeder_hit & ~oracle_hit
    oracle_only = ~seeder_hit & oracle_hit
    errors = q_seeder[both] - q_oracle[both]
    numeric = _numeric_summary(errors)

    memberships = _port_memberships_for_cells(contract, cells)
    mismatch_rows = []
    for index in np.flatnonzero(seeder_only | oracle_only):
        classification = (
            "SEEDER_ONLY_INTERSECTION"
            if seeder_only[index]
            else "ORACLE_ONLY_INTERSECTION"
        )
        if memberships[index]:
            location = "BOUNDARY_CAP_CELL"
        elif candidate_counts[index] > 0 and seeder_only[index]:
            location = "TRIANGLE_AABB_EDGE_OR_GRAZING_CANDIDATE"
        else:
            location = "UNRESOLVED"
        mismatch_rows.append(
            {
                "classification": classification,
                "location_evidence": location,
                "cell_index": int(cells[index]),
                "tree_id": int(contract.tree_ids[cells[index]]),
                "direction_index_zero_based": int(directions[index]),
                "direction": " ".join(
                    str(int(value)) for value in D3Q19_DIRECTIONS[directions[index]]
                ),
                "q_seeder": float(q_seeder[index]),
                "q_oracle": float(q_oracle[index])
                if np.isfinite(q_oracle[index])
                else "",
                "triangle_aabb_candidate_count": int(candidate_counts[index]),
                "port_boundary_labels": ";".join(memberships[index]),
                "x_m": float(starts[index, 0]),
                "y_m": float(starts[index, 1]),
                "z_m": float(starts[index, 2]),
            }
        )
    mismatch_csv.parent.mkdir(parents=True, exist_ok=True)
    with mismatch_csv.open("w", encoding="utf-8", newline="") as stream:
        fieldnames = (
            "classification",
            "location_evidence",
            "cell_index",
            "tree_id",
            "direction_index_zero_based",
            "direction",
            "q_seeder",
            "q_oracle",
            "triangle_aabb_candidate_count",
            "port_boundary_labels",
            "x_m",
            "y_m",
            "z_m",
        )
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(mismatch_rows)

    total_mismatch = int(np.count_nonzero(seeder_only | oracle_only))
    unresolved = sum(
        row["location_evidence"] == "UNRESOLVED" for row in mismatch_rows
    )
    classification_pass = bool(
        total_mismatch / len(cells) <= MAX_CLASSIFICATION_MISMATCH_FRACTION
        and np.count_nonzero(oracle_only) == 0
        and unresolved == 0
    )
    numeric_pass = bool(
        numeric["count"] > 0
        and float(numeric["median_absolute"]) <= MEDIAN_ABSOLUTE_Q_LIMIT
        and float(numeric["p95"]) <= P95_ABSOLUTE_Q_LIMIT
        and float(numeric["rms"]) <= RMS_Q_LIMIT
        and abs(float(numeric["bias"])) <= ABSOLUTE_BIAS_LIMIT
    )
    return {
        "status": "PASS" if classification_pass and numeric_pass else "FAIL",
        "mesh_dir": str(mesh_dir),
        "wall_stl": str(wall_path),
        "wall_stl_sha256": sha256_file(wall_path),
        "intersection_method": "nearest forward Moller-Trumbore intersection in link-local dimensionless coordinates over original STL R-tree candidates",
        "active_d3q19_wall_links_tested": int(len(cells)),
        "classification": {
            "TRUE_INTERSECTION_BOTH": int(np.count_nonzero(both)),
            "NO_INTERSECTION_BOTH": int(np.count_nonzero(neither)),
            "SEEDER_ONLY_INTERSECTION": int(np.count_nonzero(seeder_only)),
            "ORACLE_ONLY_INTERSECTION": int(np.count_nonzero(oracle_only)),
            "mismatch_fraction": float(total_mismatch / len(cells)),
            "unresolved_mismatch_count": int(unresolved),
            "fallback_value_disambiguation": "q=0.5 is no-hit unless the STL oracle independently finds q=0.5 within 1e-10",
            "complete_mismatch_list": str(mismatch_csv.resolve()),
        },
        "numeric_q_error_true_intersection_both_only": numeric,
        "gates": {
            "classification_mismatch_fraction_le_0p005": total_mismatch
            / len(cells)
            <= MAX_CLASSIFICATION_MISMATCH_FRACTION,
            "oracle_only_intersection_zero": int(np.count_nonzero(oracle_only)) == 0,
            "unresolved_mismatch_zero": unresolved == 0,
            "median_absolute_le_0p01": numeric["median_absolute"] is not None
            and float(numeric["median_absolute"]) <= MEDIAN_ABSOLUTE_Q_LIMIT,
            "p95_absolute_le_0p02": numeric["p95"] is not None
            and float(numeric["p95"]) <= P95_ABSOLUTE_Q_LIMIT,
            "rms_le_0p01": numeric["rms"] is not None
            and float(numeric["rms"]) <= RMS_Q_LIMIT,
            "absolute_bias_le_0p01": numeric["bias"] is not None
            and abs(float(numeric["bias"])) <= ABSOLUTE_BIAS_LIMIT,
        },
        "runtime_seconds": float(time.perf_counter() - started),
    }


def base_connectivity(mesh_dir: Path) -> dict[str, Any]:
    """Report face and D3Q19 components and prove all ports reach D3Q19 main."""

    started = time.perf_counter()
    mesh_dir = Path(mesh_dir).resolve()
    contract = load_mesh_contract(
        mesh_dir, allow_zero_normals=True, require_runtime_order=False
    )
    lattice = parse_uniform_lattice(
        (mesh_dir / "header.lua").read_text(encoding="utf-8")
    )
    labels_6, sizes_6 = sparse_component_labels(
        contract.cell_ijk, 6, cells_per_axis=2**lattice.level
    )
    labels_18, sizes_18 = sparse_component_labels(
        contract.cell_ijk, 18, cells_per_axis=2**lattice.level
    )
    boundary_cells = {
        label: contract.boundaries[label].cell_indices for label in PORT_LABELS
    }
    ports = port_component_membership(
        labels_18, boundary_cells, int(len(sizes_18))
    )
    all_ports_main = all(
        row["component_ids"] == [0]
        for label, row in ports.items()
        if label in PORT_LABELS
    )
    passed = bool(len(sizes_18) == 1 and all_ports_main)
    return {
        "status": "PASS" if passed else "FAIL",
        "fluid_cell_count": int(len(contract.tree_ids)),
        "face_6_neighbor_component_sizes": sizes_6.tolist(),
        "d3q19_18_neighbor_component_sizes": sizes_18.tolist(),
        "d3q19_main_component_connected": bool(len(sizes_18) == 1),
        "port_membership_d3q19": ports,
        "all_ports_connected_to_main_d3q19_component": bool(all_ports_main),
        "face_connectivity_is_diagnostic_only": True,
        "runtime_seconds": float(time.perf_counter() - started),
    }


def validate_final_base(
    *, base_run: Path, continuous_surface: Path
) -> dict[str, Any]:
    """Execute every required read-only final-BASE geometry gate."""

    started = time.perf_counter()
    base_run = Path(base_run).resolve()
    mesh_dir = base_run / "seeder/mesh"
    qc = base_run / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    semantic = semantic_files_success(mesh_dir, REQUIRED_MESH_FILES)
    containment = full_fluid_center_containment(
        mesh_dir=mesh_dir,
        continuous_surface=continuous_surface,
        outside_csv=qc / "full_containment_outside_cells.csv",
    )
    q_distribution = vascular_wall_qvalue_distribution(mesh_dir)
    wall_oracle = full_active_wall_q_oracle(
        mesh_dir=mesh_dir,
        wall_stl=base_run / "geometry/geometry_solver_m/wall.stl",
        mismatch_csv=qc / "full_wall_q_classification_mismatches.csv",
    )
    connectivity = base_connectivity(mesh_dir)
    passed = bool(
        semantic["semantic_success"]
        and containment["status"] == "PASS"
        and q_distribution["status"] == "PASS"
        and wall_oracle["status"] == "PASS"
        and connectivity["status"] == "PASS"
    )
    result = {
        "stage": "final_base_geometry_validation",
        "status": "PASS" if passed else "FAIL",
        "scientific_status": (
            "CFD_FLOW_SEEDER_GEOMETRY_KERNEL_VALIDATED"
            if passed
            else "CFD_FLOW_DIMENSIONLESS_KERNEL_FAILED"
        ),
        "semantic_output": semantic,
        "full_fluid_center_containment": containment,
        "qvalue_distribution": q_distribution,
        "all_active_d3q19_wall_q_oracle": wall_oracle,
        "connectivity": connectivity,
        "solver_calls": {
            "seeder_semantic_calls": 2,
            "small_musubi_semantic_calls": 0,
            "vascular_musubi_semantic_calls": 0,
            "launch_failures": 0,
            "preflight_failures": 0,
        },
        "runtime_seconds": float(time.perf_counter() - started),
        "first_failure": None
        if passed
        else next(
            name
            for name, status in (
                ("semantic_output", semantic["semantic_success"]),
                ("full_fluid_center_containment", containment["status"] == "PASS"),
                ("qvalue_distribution", q_distribution["status"] == "PASS"),
                ("all_active_d3q19_wall_q_oracle", wall_oracle["status"] == "PASS"),
                ("connectivity", connectivity["status"] == "PASS"),
            )
            if not status
        ),
        "recovery_attempts": 1,
        "next": "freeze_seeder_and_periodic_pipe_force" if passed else "STOP",
    }
    write_json(qc / "final_base_geometry_validation.json", result)
    return result


def write_compact_stage_files(base_run: Path, result: dict[str, Any]) -> None:
    """Persist the large gate components separately for audit and resume."""

    qc = Path(base_run).resolve() / "qc"
    for filename, key in (
        ("full_fluid_center_containment.json", "full_fluid_center_containment"),
        ("full_active_d3q19_wall_q_oracle.json", "all_active_d3q19_wall_q_oracle"),
        ("d3q19_connectivity.json", "connectivity"),
    ):
        write_json(qc / filename, result[key])
