"""Read-only localization audit for ideal-plane inlet D3Q19 diagonals."""

from __future__ import annotations

import csv
import subprocess
from collections import Counter, defaultdict
from dataclasses import replace
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .apes import parse_mesh_header
from .config import load_cfd_flow_config
from .exact_link_flux import reconstruct_boundary
from .geometry import load_frozen_surface_partition
from .ideal_inlet_plane import file_snapshot
from .io import FlowError, load_flow_inputs, read_json, write_json
from .port_flux_audit import (
    extract_boundary_property_indices,
    parse_bnd_header,
    parse_boundary_property_header,
    read_boundary_ids,
)
from .restart_decode import (
    D3Q19_DIRECTIONS,
    read_treelm_elemlist,
    tree_ids_to_ijk,
    tree_levels,
)


OUTPUT_GLOB = "axis_aligned_ideal_plane_inlet_preflight_anchor003274_*"
AUDIT_FILENAME = "inlet_rim_audit.json"
CSV_FILENAME = "inlet_rim_audit_cells.csv"

EXPECTED_TOTAL = 287
EXPECTED_CARDINAL = 213
EXPECTED_DIAGONAL = 74
TARGET_CARDINAL = np.asarray((0, 0, -1), dtype=np.int8)

INPUT_MISMATCH = "CFD_FLOW_INLET_RIM_AUDIT_INPUT_MISMATCH"
LOCALIZED = "CFD_FLOW_AXIS_ALIGNED_INLET_DIAGONALS_RIM_LOCALIZED"
NOT_LOCALIZED = "CFD_FLOW_AXIS_ALIGNED_INLET_DIAGONALS_NOT_RIM_LOCALIZED"
UNRESOLVED = "CFD_FLOW_AXIS_ALIGNED_INLET_RIM_LOCALIZATION_UNRESOLVED"

LOCALIZED_NEXT = (
    "REVIEW WHETHER 100_PERCENT_CARDINAL_GATE_IS_VALID_FOR_FINITE_CIRCULAR_INLET"
)
NOT_LOCALIZED_NEXT = "BUILD SHORT AXIS-ALIGNED STRAIGHT INLET EXTENSION PREFLIGHT"
UNRESOLVED_NEXT = "REVIEW PHYSICAL INLET SEAM RECONSTRUCTION"


def _git_value(repository: Path, *arguments: str) -> str:
    process = subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return process.stdout.strip()


def find_latest_successful_preflight(output_root: Path) -> tuple[Path, dict[str, Any]]:
    """Find the newest existing ideal-plane run whose single Seeder call succeeded."""

    for run_root in sorted(Path(output_root).glob(OUTPUT_GLOB), reverse=True):
        manifest_path = run_root / "qc" / "ideal_plane_preflight_manifest.json"
        if not manifest_path.is_file():
            continue
        manifest = read_json(manifest_path)
        mesh_dir = run_root / "seeder" / "mesh"
        required = ("header.lua", "elemlist.lsb", "bnd.lua", "bnd.lsb")
        if manifest.get("seeder_return_code") == 0 and all(
            (mesh_dir / name).is_file() for name in required
        ):
            return run_root.resolve(), manifest
    raise FileNotFoundError("No successful axis-aligned ideal-plane Seeder preflight found")


def classify_inlet_normals(normal_indices: np.ndarray) -> dict[str, Any]:
    """Preserve the existing cardinal/diagonal D3Q19 preflight classification."""

    values = np.asarray(normal_indices, dtype=np.int64).reshape(-1)
    vectors = D3Q19_DIRECTIONS[values].astype(np.int8, copy=False)
    cardinal_mask = np.all(vectors == TARGET_CARDINAL[None, :], axis=1)
    diagonal_mask = np.linalg.norm(vectors.astype(np.float64), axis=1) > 1.0
    other_mask = ~(cardinal_mask | diagonal_mask)
    return {
        "vectors": vectors,
        "cardinal_mask": cardinal_mask,
        "diagonal_mask": diagonal_mask,
        "other_mask": other_mask,
        "total_count": int(len(values)),
        "cardinal_count": int(np.count_nonzero(cardinal_mask)),
        "diagonal_count": int(np.count_nonzero(diagonal_mask)),
        "other_count": int(np.count_nonzero(other_mask)),
    }


def ordered_boundary_loop(faces: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the single ordered boundary loop and its undirected seam edges."""

    triangles = np.asarray(faces, dtype=np.int64)
    if triangles.ndim != 2 or triangles.shape[1] != 3 or len(triangles) == 0:
        raise ValueError("Inlet patch must contain triangular faces")
    edges = np.concatenate(
        (triangles[:, (0, 1)], triangles[:, (1, 2)], triangles[:, (2, 0)]), axis=0
    )
    undirected = np.sort(edges, axis=1)
    seam_edges = np.asarray(
        [edge for edge, count in Counter(map(tuple, undirected)).items() if count == 1],
        dtype=np.int64,
    )
    if seam_edges.ndim != 2 or seam_edges.shape[1] != 2 or len(seam_edges) < 3:
        raise ValueError("Physical inlet cap has no recoverable boundary loop")

    adjacency: dict[int, list[int]] = defaultdict(list)
    for first, second in seam_edges:
        adjacency[int(first)].append(int(second))
        adjacency[int(second)].append(int(first))
    if any(len(neighbors) != 2 for neighbors in adjacency.values()):
        raise ValueError("Physical inlet seam is not a degree-two closed loop")

    start = min(adjacency)
    ordered = [start]
    previous: int | None = None
    current = start
    while True:
        candidates = adjacency[current]
        following = candidates[0] if candidates[0] != previous else candidates[1]
        if following == start:
            break
        if following in ordered:
            raise ValueError("Physical inlet seam self-closes before visiting all vertices")
        ordered.append(following)
        previous, current = current, following
    if len(ordered) != len(adjacency):
        raise ValueError("Physical inlet cap contains multiple boundary loops")
    return np.asarray(ordered, dtype=np.int64), seam_edges


def seam_distances_xy(points_xy: np.ndarray, polygon_xy: np.ndarray) -> np.ndarray:
    """Minimum Euclidean distance from XY points to a closed polygonal seam."""

    points = np.asarray(points_xy, dtype=np.float64)
    polygon = np.asarray(polygon_xy, dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points_xy must have shape (n, 2)")
    if polygon.ndim != 2 or polygon.shape[1] != 2 or len(polygon) < 3:
        raise ValueError("polygon_xy must contain at least three XY vertices")
    starts = polygon
    deltas = np.roll(polygon, -1, axis=0) - starts
    squared_lengths = np.einsum("ij,ij->i", deltas, deltas)
    if np.any(squared_lengths <= 0.0):
        raise ValueError("Physical inlet seam contains a zero-length segment")
    offsets = points[:, None, :] - starts[None, :, :]
    fractions = np.einsum("nsi,si->ns", offsets, deltas) / squared_lengths[None, :]
    fractions = np.clip(fractions, 0.0, 1.0)
    closest = starts[None, :, :] + fractions[:, :, None] * deltas[None, :, :]
    return np.min(np.linalg.norm(points[:, None, :] - closest, axis=2), axis=1)


def points_in_polygon_xy(
    points_xy: np.ndarray, polygon_xy: np.ndarray, *, tolerance: float = 1.0e-15
) -> np.ndarray:
    """Ray-crossing containment; points on the physical seam count as inside."""

    points = np.asarray(points_xy, dtype=np.float64)
    polygon = np.asarray(polygon_xy, dtype=np.float64)
    distances = seam_distances_xy(points, polygon)
    on_boundary = distances <= float(tolerance)
    inside = np.zeros(len(points), dtype=bool)
    x = points[:, 0]
    y = points[:, 1]
    following = np.roll(polygon, -1, axis=0)
    for start, end in zip(polygon, following, strict=True):
        crosses_y = (start[1] > y) != (end[1] > y)
        denominator = end[1] - start[1]
        if denominator == 0.0:
            continue
        crossing_x = start[0] + (y - start[1]) * (end[0] - start[0]) / denominator
        inside ^= crosses_y & (x < crossing_x)
    inside[on_boundary] = True
    return inside


def _reconstruct_inlet(mesh_dir: Path) -> tuple[Any, np.ndarray, dict[str, Any]]:
    """Use the established TreElm and D3Q19 boundary reconstruction path."""

    mesh = parse_mesh_header(mesh_dir)
    header_path = Path(mesh["header"])
    property_header = parse_boundary_property_header(
        header_path.read_text(encoding="utf-8")
    )
    boundary_header = parse_bnd_header(
        (mesh_dir / "bnd.lua").read_text(encoding="utf-8")
    )
    tree_ids, property_bits, contract = read_treelm_elemlist(
        mesh_dir / "elemlist.lsb", n_elems=int(mesh["fluid_element_count"])
    )
    property_indices = extract_boundary_property_indices(
        property_bits, property_header.bit_position
    )
    if len(property_indices) != property_header.element_count:
        raise ValueError("Boundary property mapping count differs from the mesh header")
    boundary_ids = read_boundary_ids(
        mesh_dir / "bnd.lsb",
        element_count=property_header.element_count,
        side_count=boundary_header.side_count,
    )
    label_to_id = {
        label: index for index, label in enumerate(boundary_header.labels, start=1)
    }
    inlet = reconstruct_boundary(
        boundary_ids,
        property_indices,
        label="inlet",
        boundary_id=label_to_id["inlet"],
    )
    return inlet, tree_ids, {
        "mesh": mesh,
        "elemlist_contract": contract,
        "boundary_parser": (
            "utils.cfd_flow.port_flux_audit + "
            "utils.cfd_flow.exact_link_flux.reconstruct_boundary"
        ),
        "tree_geometry_parser": "utils.cfd_flow.restart_decode.tree_ids_to_ijk",
        "third_boundary_parser_written": False,
    }


def _frozen_paths_from_manifest(manifest: dict[str, Any]) -> tuple[Path, ...]:
    records = manifest.get("frozen_files_before", {})
    if not isinstance(records, dict) or not records:
        raise ValueError("Successful preflight has no frozen-file manifest")
    return tuple(Path(path) for path in records)


def _mesh_paths(mesh_dir: Path) -> tuple[Path, ...]:
    paths = tuple(sorted(path for path in mesh_dir.rglob("*") if path.is_file()))
    if not paths:
        raise ValueError("Successful preflight mesh directory is empty")
    return paths


def _physical_inlet_polygon(
    project_root: Path, run_root: Path
) -> tuple[np.ndarray, dict[str, Any], tuple[Path, ...]]:
    config = load_cfd_flow_config(
        project_root / "configs" / "cfd_flow.yaml", project_root=project_root
    )
    inputs = load_flow_inputs(config.paths.source_surface_run)
    provenance = read_json(run_root / "input" / "source_provenance.json")
    rotated_run = Path(provenance["previous_rotated_run"]).resolve()
    rotated_inputs = replace(
        inputs,
        run_root=rotated_run,
        tagged_surface_vtp=rotated_run
        / "geometry"
        / "cfd_surface_axis_aligned_inlet_um.vtp",
        meter_surface_stl=rotated_run
        / "geometry"
        / "cfd_surface_axis_aligned_inlet_m.stl",
    )
    solver_geometry = rotated_run / "geometry" / "geometry_solver_m"
    partition = load_frozen_surface_partition(rotated_inputs, solver_geometry)
    inlet = partition.patch("inlet")
    inlet_faces = partition.faces[inlet.face_indices]
    loop, seam_edges = ordered_boundary_loop(inlet_faces)
    loop_points_m = partition.points_um[loop] * 1.0e-6
    z_span_m = float(np.ptp(loop_points_m[:, 2]))
    signed_area = 0.5 * float(
        np.sum(
            loop_points_m[:, 0] * np.roll(loop_points_m[:, 1], -1)
            - np.roll(loop_points_m[:, 0], -1) * loop_points_m[:, 1]
        )
    )
    if abs(signed_area) <= np.finfo(float).tiny:
        raise ValueError("Physical inlet seam has zero XY polygon area")
    source_paths = (
        rotated_inputs.tagged_surface_vtp,
        rotated_inputs.meter_surface_stl,
        *(solver_geometry / f"{label}.stl" for label in (
            "wall", "inlet", "outlet_01", "outlet_02", "outlet_03"
        )),
    )
    return loop_points_m[:, :2], {
        "source_rotated_run": str(rotated_run),
        "source_tagged_geometry": str(rotated_inputs.tagged_surface_vtp),
        "source_patch": "inlet",
        "source_patch_entity_id": int(inlet.entity_id),
        "source_patch_triangle_count": int(inlet.triangle_count),
        "seam_edge_count": int(len(seam_edges)),
        "seam_vertex_count": int(len(loop)),
        "closed_single_loop": True,
        "maximum_seam_z_span_m": z_span_m,
        "xy_signed_polygon_area_m2": signed_area,
        "circle_fit_used": False,
        "cano_nd_rectangle_used": False,
    }, source_paths


def _distribution(normal_indices: np.ndarray) -> list[dict[str, int]]:
    unique, counts = np.unique(np.asarray(normal_indices, dtype=np.int64), return_counts=True)
    return [
        {
            "direction_index": int(index + 1),
            "cx": int(D3Q19_DIRECTIONS[index, 0]),
            "cy": int(D3Q19_DIRECTIONS[index, 1]),
            "cz": int(D3Q19_DIRECTIONS[index, 2]),
            "element_count": int(count),
        }
        for index, count in zip(unique, counts, strict=True)
    ]


def _write_cells_csv(
    path: Path,
    *,
    inlet: Any,
    tree_ids: np.ndarray,
    centers_m: np.ndarray,
    vectors: np.ndarray,
    classification: np.ndarray,
    distances_m: np.ndarray,
    dx_m: float,
    inside: np.ndarray,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            (
                "cell_index", "tree_id", "center_x_m", "center_y_m", "center_z_m",
                "normalInd", "normal_vector", "classification", "rim_distance_m",
                "rim_distance_over_dx", "inside_physical_inlet",
            )
        )
        for local in range(len(inlet.cell_indices)):
            vector = vectors[local]
            writer.writerow(
                (
                    int(inlet.cell_indices[local]),
                    int(tree_ids[inlet.cell_indices[local]]),
                    *(f"{value:.17g}" for value in centers_m[local]),
                    int(inlet.normal_indices[local] + 1),
                    f"({int(vector[0])},{int(vector[1])},{int(vector[2])})",
                    str(classification[local]),
                    f"{distances_m[local]:.17g}",
                    f"{distances_m[local] / dx_m:.17g}",
                    bool(inside[local]),
                )
            )


def run_inlet_rim_audit(project_root: Path) -> dict[str, Any]:
    """Audit one existing successful mesh without invoking any APES executable."""

    root = Path(project_root).resolve()
    config = load_cfd_flow_config(
        root / "configs" / "cfd_flow.yaml", project_root=root
    )
    run_root, manifest = find_latest_successful_preflight(config.paths.output_root)
    mesh_dir = run_root / "seeder" / "mesh"
    output_path = run_root / "qc" / AUDIT_FILENAME
    csv_path = run_root / "qc" / CSV_FILENAME

    frozen_paths = _frozen_paths_from_manifest(manifest)
    expected_frozen = manifest["frozen_files_before"]
    frozen_before = file_snapshot(frozen_paths)
    mesh_paths = _mesh_paths(mesh_dir)
    mesh_before = file_snapshot(mesh_paths)

    report: dict[str, Any] = {
        "status": UNRESOLVED,
        "next": UNRESOLVED_NEXT,
        "actual_head": _git_value(root, "rev-parse", "HEAD"),
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "musubi_calls": 0,
        "harvester_calls": 0,
        "grid_convergence": "NOT_RUN",
        "mesh_source_path": str(mesh_dir),
        "source_preflight_manifest": str(
            run_root / "qc" / "ideal_plane_preflight_manifest.json"
        ),
        "dx_m": float(manifest["dx_m"]),
        "rim_localized": "UNRESOLVED",
        "source_frozen_files_unchanged": False,
    }

    inlet, tree_ids, parser_contract = _reconstruct_inlet(mesh_dir)
    classes = classify_inlet_normals(inlet.normal_indices)
    report.update(
        {
            "parser_reuse": parser_contract,
            "total_inlet_d3q19_cells": classes["total_count"],
            "negative_z_cardinal_count": classes["cardinal_count"],
            "diagonal_count": classes["diagonal_count"],
            "other_count": classes["other_count"],
            "normal_ind_distribution": _distribution(inlet.normal_indices),
            "classification_contract": (
                "Existing preflight contract: exact (0,0,-1) is cardinal; all "
                "non-cardinal D3Q19 prevailing directions are diagonal; remainder is other"
            ),
        }
    )
    actual = (
        classes["total_count"], classes["cardinal_count"], classes["diagonal_count"]
    )
    expected = (EXPECTED_TOTAL, EXPECTED_CARDINAL, EXPECTED_DIAGONAL)
    if actual != expected:
        report.update(
            {
                "status": INPUT_MISMATCH,
                "next": UNRESOLVED_NEXT,
                "input_consistency": {
                    "status": "FAIL", "expected_total_cardinal_diagonal": expected,
                    "actual_total_cardinal_diagonal": actual,
                },
            }
        )
        report["source_frozen_files_unchanged"] = (
            frozen_before == expected_frozen
            and file_snapshot(frozen_paths) == frozen_before
            and file_snapshot(mesh_paths) == mesh_before
        )
        write_json(output_path, report)
        return report
    report["input_consistency"] = {
        "status": "PASS",
        "expected_total_cardinal_diagonal": list(expected),
        "actual_total_cardinal_diagonal": list(actual),
    }

    dx_m = float(manifest["dx_m"])
    levels = tree_levels(tree_ids)
    root_level = int(manifest["root_level"])
    if not np.all(levels == root_level):
        raise FlowError(UNRESOLVED, "Inlet rim audit requires the existing uniform-level mesh")
    phase_qc = read_json(
        Path(read_json(run_root / "input" / "source_provenance.json")["previous_rotated_run"])
        / "qc"
        / "inlet_grid_phase_qc.json"
    )
    origin_m = np.asarray(phase_qc["bounding_cube_origin_m"], dtype=np.float64)
    cell_tree_ids = tree_ids[inlet.cell_indices]
    cell_ijk = tree_ids_to_ijk(cell_tree_ids, levels[inlet.cell_indices])
    centers_m = origin_m[None, :] + (cell_ijk.astype(np.float64) + 0.5) * dx_m

    polygon_xy, seam_qc, physical_source_paths = _physical_inlet_polygon(root, run_root)
    physical_before = file_snapshot(physical_source_paths)
    distances_m = seam_distances_xy(centers_m[:, :2], polygon_xy)
    tolerance = max(np.finfo(float).eps * float(np.max(np.abs(polygon_xy))) * 16.0, 1.0e-15)
    inside = points_in_polygon_xy(centers_m[:, :2], polygon_xy, tolerance=tolerance)

    class_labels = np.full(len(inlet.cell_indices), "other", dtype=object)
    class_labels[classes["cardinal_mask"]] = "clean_cardinal"
    class_labels[classes["diagonal_mask"]] = "diagonal"
    diagonal = classes["diagonal_mask"]
    diagonal_over_dx = distances_m[diagonal] / dx_m
    diagonal_inside = inside[diagonal]
    count_le_1 = int(np.count_nonzero(diagonal_over_dx <= 1.0 + 1.0e-12))
    count_le_2 = int(np.count_nonzero(diagonal_over_dx <= 2.0 + 1.0e-12))
    count_gt_2 = int(np.count_nonzero(diagonal_over_dx > 2.0 + 1.0e-12))
    count_gt_3 = int(np.count_nonzero(diagonal_over_dx > 3.0 + 1.0e-12))
    diagonal_total = int(len(diagonal_over_dx))
    fraction_le_1 = count_le_1 / diagonal_total
    fraction_le_2 = count_le_2 / diagonal_total
    fraction_gt_2 = count_gt_2 / diagonal_total
    if fraction_le_2 >= 0.95 and count_gt_3 == 0:
        localized, status, next_step = "YES", LOCALIZED, LOCALIZED_NEXT
    elif fraction_gt_2 > 0.05 or count_gt_3 > 0:
        localized, status, next_step = "NO", NOT_LOCALIZED, NOT_LOCALIZED_NEXT
    else:
        localized, status, next_step = "UNRESOLVED", UNRESOLVED, UNRESOLVED_NEXT

    _write_cells_csv(
        csv_path,
        inlet=inlet,
        tree_ids=tree_ids,
        centers_m=centers_m,
        vectors=classes["vectors"],
        classification=class_labels,
        distances_m=distances_m,
        dx_m=dx_m,
        inside=inside,
    )
    frozen_after = file_snapshot(frozen_paths)
    mesh_after = file_snapshot(mesh_paths)
    physical_after = file_snapshot(physical_source_paths)
    unchanged = (
        frozen_before == expected_frozen
        and frozen_after == frozen_before
        and mesh_after == mesh_before
        and physical_after == physical_before
    )
    report.update(
        {
            "status": status,
            "next": next_step,
            "rim_localized": localized,
            "physical_inlet_seam": seam_qc,
            "cell_center_geometry": {
                "bounding_cube_origin_m": origin_m.tolist(),
                "root_level": root_level,
                "formula": "origin + (treeID_to_ijk + 0.5) * dx",
            },
            "diagonal_statistics": {
                "total": diagonal_total,
                "distance_le_1dx_count": count_le_1,
                "distance_le_1dx_fraction": fraction_le_1,
                "distance_le_2dx_count": count_le_2,
                "distance_le_2dx_fraction": fraction_le_2,
                "distance_gt_2dx_count": count_gt_2,
                "distance_gt_2dx_fraction": fraction_gt_2,
                "distance_gt_3dx_count": count_gt_3,
                "maximum_rim_distance_over_dx": float(np.max(diagonal_over_dx)),
                "inside_physical_footprint_count": int(np.count_nonzero(diagonal_inside)),
                "outside_physical_footprint_count": int(
                    diagonal_total - np.count_nonzero(diagonal_inside)
                ),
            },
            "decision_rule": {
                "localized": ">=95% diagonal <=2dx and zero diagonal >3dx",
                "not_localized": ">5% diagonal >2dx or any diagonal >3dx",
            },
            "cell_csv": str(csv_path),
            "source_frozen_files_unchanged": unchanged,
        }
    )
    write_json(output_path, report)
    return report
