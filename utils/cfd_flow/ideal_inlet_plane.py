"""Preflight one Cartesian numerical inlet plane on frozen rotated geometry.

This experiment is deliberately isolated from the production CFD pipeline.  It
reuses the accepted rigidly rotated wall and outlet patches, replaces only the
Seeder inlet boundary object with a two-triangle XY plane, launches Seeder at
most once, and never launches Musubi or Harvester.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import asdict, dataclass, replace
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from .apes import inspect_apes_environment, parse_mesh_header, run_wsl_tool
from .axis_aligned_inlet import (
    EXPECTED_DENSITY_KG_M3,
    EXPECTED_DX_M,
    EXPECTED_MASS_FLOW_KG_S,
    EXPECTED_SMOOTH_INLET_AREA_M2,
    EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
    OLD_BASELINE_STATUS,
    _mesh_boundary_preflight,
)
from .config import load_cfd_flow_config
from .geometry import BoundingCube, find_seed_point, load_frozen_surface_partition
from .io import FlowError, load_flow_inputs, read_json, sha256_file, write_json
from .restart_decode import read_treelm_elemlist, tree_ids_to_ijk, tree_levels


EXPECTED_BRANCH = "codex/cfd-flow-axis-aligned-inlet-preflight-20260829"
EXPECTED_HEAD = "478fab2aa727ea5e2e2c5d17859960a51beadd56"
PREVIOUS_ROTATED_RUN = "axis_aligned_inlet_geometry_anchor003274_20260829_111451"
FROZEN_STEADY_RUN = "musubi_project_steady_confirmation_anchor003274_20260828_225334"
FROZEN_DIRECT_FIELD_RUN = "musubi_direct_restart_field_anchor003274_20260829_000544"
OUTPUT_PREFIX = "axis_aligned_ideal_plane_inlet_preflight_anchor003274"
REVISION = "AXIS_ALIGNED_IDEAL_NUMERICAL_INLET_PLANE_SEEDER_PREFLIGHT_V1"

PASS_STATUS = "CFD_FLOW_AXIS_ALIGNED_IDEAL_PLANE_INLET_PREFLIGHT_PASS"
FAILED_STATUS = "CFD_FLOW_AXIS_ALIGNED_IDEAL_PLANE_INLET_PREFLIGHT_FAILED"
GEOMETRY_UNSAFE = "CFD_FLOW_IDEAL_INLET_PLANE_GEOMETRY_UNSAFE"
AREA_PROXY_FAILED = "CFD_FLOW_IDEAL_INLET_PLANE_AREA_PROXY_FAILED"
NORMALIND_NOT_CLEAN = "CFD_FLOW_IDEAL_INLET_PLANE_NORMALIND_NOT_CLEAN"
DOMAIN_CHANGED = "CFD_FLOW_IDEAL_INLET_PLANE_DOMAIN_CHANGED"
SEEDER_INPUT_FORMAT_FAILED = "CFD_FLOW_IDEAL_INLET_PLANE_SEEDER_INPUT_FORMAT_FAILED"
PASS_NEXT = "RUN ONE NEW AXIS_ALIGNED MUSUBI BASELINE"
FAILED_NEXT = "BUILD SHORT AXIS-ALIGNED STRAIGHT INLET EXTENSION"
CURRENT_SEEDER_FAILURE_NEXT = "FIX CURRENT SEEDER CANOND PREFLIGHT FAILURE"
TOPOLOGY_FAILURE_NEXT = "REVIEW NEED FOR SHORT AXIS-ALIGNED STRAIGHT INLET EXTENSION"

PREVIOUS_FLUID_CELL_COUNT = 221_359
RECTANGLE_MARGIN_CELLS = 2
FLUID_COUNT_MAXIMUM_RELATIVE_DIFFERENCE = 0.005
AREA_RATIO_MAXIMUM_ERROR = 0.05
VELOCITY_RATIO_MAXIMUM_ERROR = 0.05


@dataclass(frozen=True, slots=True)
class IdealPlaneLayout:
    root: Path
    geometry: Path
    geometry_solver_m: Path
    seeder: Path
    mesh: Path
    qc: Path
    input: Path


def _git_value(repository: Path, *arguments: str) -> str:
    process = subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return process.stdout.strip()


def _create_layout(output_root: Path) -> IdealPlaneLayout:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    root = Path(output_root) / f"{OUTPUT_PREFIX}_{stamp}"
    if root.exists():
        raise FlowError(FAILED_STATUS, f"Output already exists: {root}")
    geometry = root / "geometry"
    geometry_solver_m = geometry / "geometry_solver_m"
    seeder = root / "seeder"
    mesh = seeder / "mesh"
    qc = root / "qc"
    input_dir = root / "input"
    for directory in (geometry_solver_m, mesh, qc, input_dir):
        directory.mkdir(parents=True, exist_ok=False)
    return IdealPlaneLayout(root, geometry, geometry_solver_m, seeder, mesh, qc, input_dir)


def file_snapshot(paths: Iterable[Path]) -> dict[str, dict[str, Any]]:
    """Return a stable byte/SHA snapshot without changing any source file."""

    unique = sorted({Path(path).resolve() for path in paths}, key=str)
    missing = [str(path) for path in unique if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing frozen source files: {missing}")
    return {
        str(path): {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
        for path in unique
    }


def directory_manifest(directory: Path) -> dict[str, Any]:
    paths = sorted(path for path in Path(directory).rglob("*") if path.is_file())
    return {
        "root": str(Path(directory).resolve()),
        "file_count": len(paths),
        "total_bytes": int(sum(path.stat().st_size for path in paths)),
        "files": [
            {
                "relative_path": path.relative_to(directory).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
            for path in paths
        ],
    }


def build_grid_snapped_rectangle(
    physical_xy_bounds_m: np.ndarray,
    grid_origin_xy_m: np.ndarray,
    *,
    dx_m: float = EXPECTED_DX_M,
    margin_cells: int = RECTANGLE_MARGIN_CELLS,
) -> dict[str, Any]:
    """Expand by exactly ``margin_cells`` and snap every edge outward."""

    physical = np.asarray(physical_xy_bounds_m, dtype=np.float64).reshape(2, 2)
    origin = np.asarray(grid_origin_xy_m, dtype=np.float64).reshape(2)
    if np.any(physical[:, 1] <= physical[:, 0]):
        raise ValueError("Physical XY bounds must have positive extent")
    margin_m = int(margin_cells) * float(dx_m)
    expanded = physical.copy()
    expanded[:, 0] -= margin_m
    expanded[:, 1] += margin_m
    lower_phase = (expanded[:, 0] - origin) / float(dx_m)
    upper_phase = (expanded[:, 1] - origin) / float(dx_m)
    lower_index = np.floor(np.nextafter(lower_phase, np.inf)).astype(np.int64)
    upper_index = np.ceil(np.nextafter(upper_phase, -np.inf)).astype(np.int64)
    snapped = np.column_stack(
        (
            origin + lower_index * float(dx_m),
            origin + upper_index * float(dx_m),
        )
    )
    if np.any(snapped[:, 0] > expanded[:, 0]) or np.any(
        snapped[:, 1] < expanded[:, 1]
    ):
        raise AssertionError("Outward grid snapping contracted the rectangle")
    phases = (snapped - origin[:, None]) / float(dx_m)
    phase_error = float(np.max(np.abs(phases - np.rint(phases))))
    return {
        "physical_xy_bounds_m": physical.tolist(),
        "requested_expansion_margin_cells": int(margin_cells),
        "requested_expansion_margin_m": margin_m,
        "expanded_xy_bounds_before_snap_m": expanded.tolist(),
        "snapped_xy_bounds_m": snapped.tolist(),
        "snapped_grid_face_indices": np.column_stack(
            (lower_index, upper_index)
        ).tolist(),
        "maximum_grid_face_phase_error_over_dx": phase_error,
        "outward_only": True,
    }


def plane_vertices(snapped_xy_bounds_m: np.ndarray, plane_z_m: float) -> np.ndarray:
    bounds = np.asarray(snapped_xy_bounds_m, dtype=np.float64).reshape(2, 2)
    xmin, xmax = bounds[0]
    ymin, ymax = bounds[1]
    z = float(plane_z_m)
    return np.asarray(
        (
            (xmin, ymin, z),
            (xmax, ymin, z),
            (xmax, ymax, z),
            (xmin, ymax, z),
        ),
        dtype=np.float64,
    )


def write_two_triangle_plane_stl(path: Path, vertices_m: np.ndarray) -> dict[str, Any]:
    """Write an ASCII STL so the grid-aligned z value retains double precision."""

    vertices = np.asarray(vertices_m, dtype=np.float64).reshape(4, 3)
    faces = np.asarray(((0, 1, 2), (0, 2, 3)), dtype=np.int64)
    crosses = np.cross(
        vertices[faces[:, 1]] - vertices[faces[:, 0]],
        vertices[faces[:, 2]] - vertices[faces[:, 0]],
    )
    norms = np.linalg.norm(crosses, axis=1)
    if np.any(norms <= 0.0):
        raise ValueError("Numerical inlet plane has a degenerate triangle")
    normals = crosses / norms[:, None]
    if not np.all(normals[:, 2] > 0.0):
        raise ValueError("Numerical inlet triangles are not consistently +Z oriented")
    lines = ["solid numerical_inlet_plane"]
    for face in faces:
        lines.extend(("  facet normal 0 0 1", "    outer loop"))
        for index in face:
            x, y, z = vertices[index]
            lines.append(f"      vertex {x:.17g} {y:.17g} {z:.17g}")
        lines.extend(("    endloop", "  endfacet"))
    lines.append("endsolid numerical_inlet_plane")
    path.write_text("\n".join(lines) + "\n", encoding="ascii", newline="\n")
    return {
        "path": str(path.resolve()),
        "boundary_representation": "NUMERICAL_BOUNDARY_COVERAGE_PLANE",
        "boundary_type": "EXACTLY_TWO_TRIANGLE_CARTESIAN_PLANAR_STL",
        "triangle_count": 2,
        "vertex_count": 4,
        "vertices_m": vertices.tolist(),
        "plane_z_m": float(vertices[0, 2]),
        "all_vertex_z_exactly_equal": bool(np.all(vertices[:, 2] == vertices[0, 2])),
        "triangle_normals": normals.tolist(),
        "plane_normal": [0.0, 0.0, 1.0],
        "triangle_orientation_consistent": bool(np.all(normals[:, 2] > 0.0)),
        "sha256": sha256_file(path),
    }


def _triangle_plane_intersection_xy(
    triangle_m: np.ndarray, plane_z_m: float, tolerance_m: float
) -> np.ndarray:
    triangle = np.asarray(triangle_m, dtype=np.float64).reshape(3, 3)
    signed = triangle[:, 2] - float(plane_z_m)
    points: list[np.ndarray] = []
    for index, value in enumerate(signed):
        if abs(float(value)) <= tolerance_m:
            points.append(triangle[index, :2])
    for first, second in ((0, 1), (1, 2), (2, 0)):
        left = float(signed[first])
        right = float(signed[second])
        if left * right < 0.0:
            fraction = left / (left - right)
            point = triangle[first] + fraction * (triangle[second] - triangle[first])
            points.append(point[:2])
    unique: list[np.ndarray] = []
    for point in points:
        if not any(np.linalg.norm(point - prior) <= tolerance_m for prior in unique):
            unique.append(point)
    return np.asarray(unique, dtype=np.float64).reshape((-1, 2))


def _point_in_rectangle(point: np.ndarray, bounds: np.ndarray, tolerance: float) -> bool:
    return bool(
        bounds[0, 0] - tolerance <= point[0] <= bounds[0, 1] + tolerance
        and bounds[1, 0] - tolerance <= point[1] <= bounds[1, 1] + tolerance
    )


def _segment_intersects_rectangle(
    first: np.ndarray, second: np.ndarray, bounds: np.ndarray, tolerance: float
) -> bool:
    if _point_in_rectangle(first, bounds, tolerance) or _point_in_rectangle(
        second, bounds, tolerance
    ):
        return True
    delta = second - first
    for axis in (0, 1):
        other = 1 - axis
        if abs(float(delta[axis])) <= np.finfo(float).tiny:
            continue
        for edge in bounds[axis]:
            parameter = (edge - first[axis]) / delta[axis]
            if -tolerance <= parameter <= 1.0 + tolerance:
                value = first[other] + parameter * delta[other]
                if bounds[other, 0] - tolerance <= value <= bounds[other, 1] + tolerance:
                    return True
    return False


def assess_plane_coverage_safety(
    points_m: np.ndarray,
    faces: np.ndarray,
    inlet_face_indices: np.ndarray,
    *,
    plane_z_m: float,
    rectangle_xy_bounds_m: np.ndarray,
    tolerance_m: float = 1.0e-14,
) -> dict[str, Any]:
    """Reject a plane covering a cross-section disconnected from the inlet seam."""

    points = np.asarray(points_m, dtype=np.float64)
    triangles = np.asarray(faces, dtype=np.int64)
    inlet_faces = np.asarray(inlet_face_indices, dtype=np.int64).reshape(-1)
    inlet_face_set = set(int(value) for value in inlet_faces)
    inlet_vertices = set(int(value) for value in np.unique(triangles[inlet_faces]))
    bounds = np.asarray(rectangle_xy_bounds_m, dtype=np.float64).reshape(2, 2)

    intersecting: list[int] = []
    intersections: dict[int, np.ndarray] = {}
    for face_index, face in enumerate(triangles):
        if face_index in inlet_face_set:
            continue
        triangle = points[face]
        if float(np.min(triangle[:, 2])) > plane_z_m + tolerance_m:
            continue
        if float(np.max(triangle[:, 2])) < plane_z_m - tolerance_m:
            continue
        xy = _triangle_plane_intersection_xy(triangle, plane_z_m, tolerance_m)
        if len(xy):
            intersecting.append(face_index)
            intersections[face_index] = xy

    parent = {face_index: face_index for face_index in intersecting}

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = parent[value]
        return value

    def union(left: int, right: int) -> None:
        first = find(left)
        second = find(right)
        if first != second:
            parent[second] = first

    vertex_owner: dict[int, int] = {}
    for face_index in intersecting:
        for vertex in triangles[face_index]:
            value = int(vertex)
            if value in vertex_owner:
                union(face_index, vertex_owner[value])
            else:
                vertex_owner[value] = face_index

    components: dict[int, list[int]] = {}
    for face_index in intersecting:
        components.setdefault(find(face_index), []).append(face_index)
    records: list[dict[str, Any]] = []
    unsafe: list[dict[str, Any]] = []
    for component_id, component_faces in enumerate(components.values()):
        xy_values = np.vstack([intersections[index] for index in component_faces])
        touches_rectangle = False
        for face_index in component_faces:
            xy = intersections[face_index]
            if any(_point_in_rectangle(point, bounds, tolerance_m) for point in xy):
                touches_rectangle = True
                break
            if len(xy) >= 2 and _segment_intersects_rectangle(
                xy[0], xy[1], bounds, tolerance_m
            ):
                touches_rectangle = True
                break
        seam_connected = any(
            int(vertex) in inlet_vertices
            for face_index in component_faces
            for vertex in triangles[face_index]
        )
        record = {
            "component_id": component_id,
            "intersecting_face_count": len(component_faces),
            "xy_bounds_m": [
                [float(np.min(xy_values[:, 0])), float(np.max(xy_values[:, 0]))],
                [float(np.min(xy_values[:, 1])), float(np.max(xy_values[:, 1]))],
            ],
            "touches_numerical_rectangle": touches_rectangle,
            "connected_to_physical_inlet_seam": seam_connected,
        }
        records.append(record)
        if touches_rectangle and not seam_connected:
            unsafe.append(record)
    return {
        "status": "PASS" if not unsafe else GEOMETRY_UNSAFE,
        "method": "PLANE_INTERSECTION_COMPONENTS_WITH_INLET_SEAM_ADJACENCY",
        "allows_inlet_wall_seam_intersection": True,
        "plane_z_m": float(plane_z_m),
        "rectangle_xy_bounds_m": bounds.tolist(),
        "intersection_component_count": len(records),
        "components": records,
        "unsafe_component_count": len(unsafe),
        "unsafe_components": unsafe,
    }


def generate_ideal_plane_seeder_lua(
    *,
    seed_m: np.ndarray,
    cube: BoundingCube,
    rectangle_xy_bounds_m: np.ndarray,
    plane_z_m: float,
) -> str:
    def number(value: float) -> str:
        return f"{float(value):.17g}"

    def vector(values: np.ndarray) -> str:
        return "{ " + ", ".join(number(item) for item in values) + " }"

    bounds = np.asarray(rectangle_xy_bounds_m, dtype=np.float64).reshape(2, 2)
    plane_origin = np.asarray((bounds[0, 0], bounds[1, 0], plane_z_m))
    plane_x_vector = np.asarray((bounds[0, 1] - bounds[0, 0], 0.0, 0.0))
    plane_y_vector = np.asarray((0.0, bounds[1, 1] - bounds[1, 0], 0.0))
    plane_vectors = (
        "{ " + vector(plane_x_vector) + ", " + vector(plane_y_vector) + " }"
    )
    objects = [
        "  {\n"
        "    attribute = { kind = 'boundary', label = 'wall', level = minlevel, calc_dist = true },\n"
        "    geometry = { kind = 'stl', object = { filename = '../geometry/geometry_solver_m/wall.stl' } }\n"
        "  }",
        "  {\n"
        "    attribute = { kind = 'boundary', label = 'inlet', level = minlevel },\n"
        "    geometry = { kind = 'canoND', object = { origin = "
        f"{vector(plane_origin)}, vec = {plane_vectors}"
        " } }\n"
        "  }",
    ]
    for label in ("outlet_01", "outlet_02", "outlet_03"):
        objects.append(
            "  {\n"
            f"    attribute = {{ kind = 'boundary', label = '{label}', level = minlevel }},\n"
            "    geometry = { kind = 'stl', object = { filename = "
            f"'../geometry/geometry_solver_m/{label}.stl' }} }}\n"
            "  }"
        )
    objects.append(
        "  {\n"
        "    attribute = { kind = 'seed' },\n"
        f"    geometry = {{ kind = 'canoND', object = {{ origin = {vector(seed_m)} }} }}\n"
        "  }"
    )
    return (
        "-- Axis-aligned ideal numerical inlet plane preflight; Seeder only.\n"
        "folder = 'mesh/'\n"
        "comment = 'ROI003274 ideal numerical inlet plane, uniform 0.20 um lattice'\n"
        "debug = { debugMode = false, debugFiles = false, debugMesh = 'debug/' }\n"
        f"minlevel = {cube.level}\n"
        f"bounding_cube = {{ origin = {vector(cube.origin_m)}, length = {number(cube.side_m)} }}\n"
        "spatial_object = {\n"
        + ",\n".join(objects)
        + "\n}\n"
    )


def fluid_count_comparison_qc(
    new_count: int, previous_count: int = PREVIOUS_FLUID_CELL_COUNT
) -> dict[str, Any]:
    difference = abs(int(new_count) - int(previous_count)) / int(previous_count)
    return {
        "status": "PASS"
        if difference <= FLUID_COUNT_MAXIMUM_RELATIVE_DIFFERENCE
        else "FAIL",
        "previous_fluid_cell_count": int(previous_count),
        "new_fluid_cell_count": int(new_count),
        "absolute_cell_count_difference": abs(int(new_count) - int(previous_count)),
        "relative_difference": difference,
        "maximum_relative_difference": FLUID_COUNT_MAXIMUM_RELATIVE_DIFFERENCE,
    }


def connected_fluid_region_count(tree_ids: np.ndarray, level: int) -> int:
    """Count face-connected regions of one uniform TreElm level."""

    ids = np.asarray(tree_ids, dtype=np.int64).reshape(-1)
    levels = tree_levels(ids)
    if len(ids) == 0 or not np.all(levels == int(level)):
        raise ValueError("Fluid connectivity requires a non-empty uniform-level mesh")
    coordinates = tree_ids_to_ijk(ids, levels)
    cells_per_axis = 2**int(level)
    codes = (
        coordinates[:, 0].astype(np.int64)
        + cells_per_axis * coordinates[:, 1].astype(np.int64)
        + cells_per_axis**2 * coordinates[:, 2].astype(np.int64)
    )
    if len(np.unique(codes)) != len(codes):
        raise ValueError("TreElm fluid cells contain duplicate coordinates")
    lookup = {int(code): index for index, code in enumerate(codes)}
    parent = np.arange(len(codes), dtype=np.int64)
    ranks = np.zeros(len(codes), dtype=np.int8)

    def find(value: int) -> int:
        while parent[value] != value:
            parent[value] = parent[parent[value]]
            value = int(parent[value])
        return value

    def union(left: int, right: int) -> None:
        first = find(left)
        second = find(right)
        if first == second:
            return
        if ranks[first] < ranks[second]:
            first, second = second, first
        parent[second] = first
        if ranks[first] == ranks[second]:
            ranks[first] += 1

    offsets = (1, cells_per_axis, cells_per_axis**2)
    limits = (
        coordinates[:, 0] + 1 < cells_per_axis,
        coordinates[:, 1] + 1 < cells_per_axis,
        coordinates[:, 2] + 1 < cells_per_axis,
    )
    for index, code in enumerate(codes):
        for offset, valid in zip(offsets, limits, strict=True):
            if bool(valid[index]):
                neighbor = lookup.get(int(code + offset))
                if neighbor is not None:
                    union(index, neighbor)
    return len({find(index) for index in range(len(codes))})


def _frozen_source_paths(
    *,
    output_root: Path,
    previous_run: Path,
    original_surface_files: Iterable[Path],
) -> tuple[Path, ...]:
    previous_files = tuple(
        path for path in previous_run.rglob("*") if path.is_file()
    )
    restart_dir = output_root / FROZEN_STEADY_RUN / "restart"
    direct_field_dir = output_root / FROZEN_DIRECT_FIELD_RUN / "flow"
    restart_files = tuple(path for path in restart_dir.rglob("*") if path.is_file())
    direct_files = tuple(path for path in direct_field_dir.rglob("*") if path.is_file())
    if not restart_files or not direct_files:
        raise FileNotFoundError("Frozen restart or direct field is missing")
    return tuple(original_surface_files) + previous_files + restart_files + direct_files


def _copy_unchanged_solver_boundaries(source: Path, destination: Path) -> None:
    for label in ("wall", "outlet_01", "outlet_02", "outlet_03"):
        shutil.copy2(source / f"{label}.stl", destination / f"{label}.stl")


def run_ideal_inlet_plane_preflight(project_root: Path) -> dict[str, Any]:
    """Build the ideal plane, run one Seeder, and stop before all solvers."""

    root = Path(project_root).resolve()
    branch = _git_value(root, "branch", "--show-current")
    head = _git_value(root, "rev-parse", "HEAD")
    if branch != EXPECTED_BRANCH or head != EXPECTED_HEAD:
        raise FlowError(
            FAILED_STATUS,
            f"Expected {EXPECTED_BRANCH}@{EXPECTED_HEAD}, found {branch}@{head}",
        )
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = config.paths.output_root.resolve()
    previous_run = output_root / PREVIOUS_ROTATED_RUN
    previous_manifest = read_json(
        previous_run / "qc" / "axis_aligned_inlet_preflight_manifest.json"
    )
    previous_geometry_qc = read_json(
        previous_run / "qc" / "geometry_rigid_transform_qc.json"
    )
    previous_phase_qc = read_json(previous_run / "qc" / "inlet_grid_phase_qc.json")
    if previous_geometry_qc.get("status") != "PASS":
        raise FlowError(FAILED_STATUS, "Previous rigid rotation geometry QC is not PASS")
    if previous_phase_qc.get("status") != "PASS":
        raise FlowError(FAILED_STATUS, "Previous inlet grid phase QC is not PASS")

    inputs = load_flow_inputs(config.paths.source_surface_run)
    previous_solver_geometry = previous_run / "geometry" / "geometry_solver_m"
    rotated_inputs = replace(
        inputs,
        run_root=previous_run,
        tagged_surface_vtp=previous_run
        / "geometry"
        / "cfd_surface_axis_aligned_inlet_um.vtp",
        meter_surface_stl=previous_run
        / "geometry"
        / "cfd_surface_axis_aligned_inlet_m.stl",
    )
    partition = load_frozen_surface_partition(rotated_inputs, previous_solver_geometry)
    inlet = partition.patch("inlet")
    frozen_paths = _frozen_source_paths(
        output_root=output_root,
        previous_run=previous_run,
        original_surface_files=(inputs.tagged_surface_vtp, inputs.meter_surface_stl),
    )
    frozen_before = file_snapshot(frozen_paths)
    layout = _create_layout(output_root)
    manifest_path = layout.qc / "ideal_plane_preflight_manifest.json"
    summary: dict[str, Any] = {
        "status": FAILED_STATUS,
        "next": CURRENT_SEEDER_FAILURE_NEXT,
        "revision": REVISION,
        "run_root": str(layout.root),
        "branch": branch,
        "actual_head": head,
        "production_pipeline_modified": False,
        "previous_rotated_geometry_modified": False,
        "numerical_inlet_boundary_type": "OFFICIAL_CANOND_CARTESIAN_PLANE",
        "seeder_run_count": 0,
        "seeder_return_code": "NOT_RUN",
        "musubi_run_count": 0,
        "harvester_run_count": 0,
        "grid_convergence": "NOT_RUN",
        "frozen_files_before": frozen_before,
        "started_at": datetime.now().isoformat(),
    }
    write_json(manifest_path, summary)

    try:
        _copy_unchanged_solver_boundaries(
            previous_solver_geometry, layout.geometry_solver_m
        )
        copy_qc = {
            label: {
                "source": str(previous_solver_geometry / f"{label}.stl"),
                "copy": str(layout.geometry_solver_m / f"{label}.stl"),
                "source_sha256": sha256_file(previous_solver_geometry / f"{label}.stl"),
                "copy_sha256": sha256_file(layout.geometry_solver_m / f"{label}.stl"),
            }
            for label in ("wall", "outlet_01", "outlet_02", "outlet_03")
        }
        if not all(
            item["source_sha256"] == item["copy_sha256"] for item in copy_qc.values()
        ):
            raise FlowError(FAILED_STATUS, "A copied wall/outlet patch changed bytes")

        phase_origin = np.asarray(
            previous_phase_qc["bounding_cube_origin_m"], dtype=np.float64
        )
        cube = BoundingCube(
            phase_origin,
            float(previous_phase_qc["bounding_cube_side_m"]),
            int(previous_phase_qc["root_level"]),
            int(previous_phase_qc["cells_per_axis"]),
            float(previous_phase_qc["minimum_margin_cells"]),
        )
        plane_z_m = float(previous_phase_qc["inlet_plane_z_m"])
        inlet_vertices = np.unique(partition.faces[inlet.face_indices])
        inlet_points_m = partition.points_um[inlet_vertices] * 1.0e-6
        physical_bounds = np.asarray(
            (
                (float(np.min(inlet_points_m[:, 0])), float(np.max(inlet_points_m[:, 0]))),
                (float(np.min(inlet_points_m[:, 1])), float(np.max(inlet_points_m[:, 1]))),
            )
        )
        rectangle = build_grid_snapped_rectangle(
            physical_bounds,
            cube.origin_m[:2],
            dx_m=EXPECTED_DX_M,
            margin_cells=RECTANGLE_MARGIN_CELLS,
        )
        snapped_bounds = np.asarray(rectangle["snapped_xy_bounds_m"])
        vertices = plane_vertices(snapped_bounds, plane_z_m)
        plane_path = layout.geometry / "numerical_inlet_plane.stl"
        plane_qc = write_two_triangle_plane_stl(plane_path, vertices)
        phase_coordinate = (plane_z_m - cube.origin_m[2]) / EXPECTED_DX_M
        phase_error = abs(phase_coordinate - round(phase_coordinate))
        cube_upper = cube.origin_m + cube.side_m
        rectangle_inside_cube = bool(
            snapped_bounds[0, 0] >= cube.origin_m[0]
            and snapped_bounds[0, 1] <= cube_upper[0]
            and snapped_bounds[1, 0] >= cube.origin_m[1]
            and snapped_bounds[1, 1] <= cube_upper[1]
        )
        safety_qc = assess_plane_coverage_safety(
            partition.points_um * 1.0e-6,
            partition.faces,
            inlet.face_indices,
            plane_z_m=plane_z_m,
            rectangle_xy_bounds_m=snapped_bounds,
        )
        geometry_qc = {
            "status": "PASS"
            if (
                plane_qc["all_vertex_z_exactly_equal"]
                and plane_qc["triangle_orientation_consistent"]
                and phase_error == 0.0
                and rectangle_inside_cube
                and safety_qc["status"] == "PASS"
            )
            else GEOMETRY_UNSAFE,
            "representation": "NUMERICAL_BOUNDARY_COVERAGE_PLANE",
            "not_a_new_biological_lumen": True,
            "physical_smooth_inlet_area_m2": EXPECTED_SMOOTH_INLET_AREA_M2,
            "rectangle_area_not_used_as_physical_area": True,
            "plane": plane_qc,
            "rectangle": rectangle,
            "grid_phase_coordinate": phase_coordinate,
            "grid_phase_error_over_dx": phase_error,
            "bounding_cube_unchanged": True,
            "rectangle_inside_previous_bounding_cube": rectangle_inside_cube,
            "unchanged_wall_and_outlet_copies": copy_qc,
            "coverage_safety": safety_qc,
        }
        write_json(layout.qc / "ideal_inlet_plane_geometry_qc.json", geometry_qc)
        if geometry_qc["status"] != "PASS":
            raise FlowError(
                GEOMETRY_UNSAFE,
                "Ideal numerical rectangle failed plane, cube, or unrelated-lumen safety QC",
            )

        seed = find_seed_point(partition)
        seeder_lua = generate_ideal_plane_seeder_lua(
            seed_m=seed.coordinates_m,
            cube=cube,
            rectangle_xy_bounds_m=snapped_bounds,
            plane_z_m=plane_z_m,
        )
        if "geometry_solver_m/inlet.stl" in seeder_lua:
            raise FlowError(FAILED_STATUS, "Seeder still references the physical inlet STL")
        if "label = 'inlet'" not in seeder_lua or "kind = 'canoND'" not in seeder_lua:
            raise FlowError(FAILED_STATUS, "Seeder does not contain the official canoND plane")
        (layout.seeder / "seeder.lua").write_text(seeder_lua, encoding="utf-8")
        shutil.copy2(
            inputs.boundary_conditions_json,
            layout.input / "source_boundary_conditions.json",
        )
        shutil.copy2(
            inputs.boundary_manifest_csv,
            layout.input / "source_boundary_manifest.csv",
        )
        write_json(
            layout.input / "source_provenance.json",
            {
                "previous_rotated_run": str(previous_run),
                "previous_rotated_manifest_status": previous_manifest["status"],
                "previous_rotation_geometry_qc_status": previous_geometry_qc["status"],
                "previous_grid_phase_qc_status": previous_phase_qc["status"],
                "source_files_read_only": True,
                "frozen_file_count": len(frozen_before),
            },
        )

        environment = inspect_apes_environment(config.apes)
        write_json(layout.input / "apes_environment.json", asdict(environment))
        if (
            environment.status != "PASS"
            or not environment.binaries.get("seeder")
            or not environment.binaries.get("lua_compiler")
        ):
            raise FlowError(FAILED_STATUS, "Seeder or Lua compiler is unavailable")

        lua_check = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=layout.seeder,
            command=[str(environment.binaries["lua_compiler"]), "-p", "seeder.lua"],
            stdout_path=layout.seeder / "luac_stdout.log",
            stderr_path=layout.seeder / "luac_stderr.log",
            timeout_s=30,
        )
        write_json(
            layout.qc / "seeder_static_preflight_qc.json",
            {
                "status": "PASS" if lua_check.returncode == 0 else "FAIL",
                "returncode": lua_check.returncode,
            },
        )
        if lua_check.returncode != 0:
            raise FlowError(FAILED_STATUS, "Generated Seeder Lua failed syntax preflight")

        summary["seeder_run_count"] = 1
        write_json(manifest_path, summary)
        seeder_run = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=layout.seeder,
            command=[str(environment.binaries["seeder"]), "seeder.lua"],
            stdout_path=layout.seeder / "seeder_stdout.log",
            stderr_path=layout.seeder / "seeder_stderr.log",
            timeout_s=config.solver.wallclock_limit_s,
        )
        summary["seeder_return_code"] = seeder_run.returncode
        summary["seeder_wall_time_s"] = seeder_run.wall_time_s
        write_json(manifest_path, summary)
        if seeder_run.returncode != 0:
            raise FlowError(FAILED_STATUS, f"Seeder return code {seeder_run.returncode}")

        mesh_preflight = _mesh_boundary_preflight(
            layout.mesh, EXPECTED_SMOOTH_INLET_AREA_M2, EXPECTED_DX_M
        )
        area_qc = mesh_preflight["d3q19_inlet"]["area_proxy"]
        normal_qc = mesh_preflight["d3q19_inlet"]["normal_ind"]
        normal_qc["other_direction_count"] = int(
            mesh_preflight["d3q19_inlet"]["element_count"]
            - normal_qc["target_cardinal_count"]
            - normal_qc["diagonal_normal_ind_count"]
        )
        mesh_header = parse_mesh_header(layout.mesh)
        tree_ids, _, elemlist_contract = read_treelm_elemlist(
            layout.mesh / "elemlist.lsb",
            n_elems=int(mesh_header["fluid_element_count"]),
        )
        uniform_level = bool(
            mesh_header["minimum_level"] == cube.level
            and mesh_header["maximum_level"] == cube.level
        )
        region_count = (
            connected_fluid_region_count(tree_ids, cube.level)
            if uniform_level
            else None
        )
        domain_qc = fluid_count_comparison_qc(
            int(mesh_header["fluid_element_count"])
        )
        domain_qc.update(
            {
                "connected_fluid_regions": region_count,
                "single_connected_fluid_domain": region_count == 1,
                "connectivity": "6-neighbor face connectivity",
                "uniform_root_level": uniform_level,
                "root_level": cube.level,
                "elemlist_contract": elemlist_contract,
                "boundary_element_count": int(mesh_header["boundary_element_count"]),
                "qval_element_count": int(
                    mesh_header["property_element_counts"].get("has qVal", 0)
                ),
            }
        )
        if region_count != 1 or not uniform_level:
            domain_qc["status"] = "FAIL"
        write_json(layout.qc / "fluid_domain_comparison_qc.json", domain_qc)

        topology_qc = {
            "status": "PASS"
            if area_qc["status"] == "PASS"
            and normal_qc["status"] == "PASS"
            and domain_qc["status"] == "PASS"
            else "FAIL",
            "parser_reuse": {
                "boundary": "utils.cfd_flow.port_flux_audit + utils.cfd_flow.exact_link_flux.reconstruct_boundary",
                "elemlist": "utils.cfd_flow.restart_decode.read_treelm_elemlist",
                "third_parser_written": False,
            },
            "total_inlet_d3q19_elements": int(
                mesh_preflight["d3q19_inlet"]["element_count"]
            ),
            "area_proxy": area_qc,
            "normal_ind": normal_qc,
            "fluid_domain": domain_qc,
            "boundary_element_count": int(mesh_header["boundary_element_count"]),
            "qval_element_count": int(
                mesh_header["property_element_counts"].get("has qVal", 0)
            ),
        }
        write_json(layout.qc / "inlet_topology_qc.json", topology_qc)

        proxy_area = float(area_qc["mfr_eq_area_proxy_m2"])
        predicted_velocity = EXPECTED_MASS_FLOW_KG_S / (
            EXPECTED_DENSITY_KG_M3 * proxy_area
        )
        velocity_ratio = predicted_velocity / EXPECTED_SMOOTH_MEAN_VELOCITY_M_S
        prediction_qc = {
            "status": "PASS"
            if abs(velocity_ratio - 1.0) <= VELOCITY_RATIO_MAXIMUM_ERROR
            else "FAIL",
            "topology_available": True,
            "prediction_only_no_musubi_run": True,
            "target_mass_flow_kg_s": EXPECTED_MASS_FLOW_KG_S,
            "density_kg_m3": EXPECTED_DENSITY_KG_M3,
            "mfr_eq_area_proxy_m2": proxy_area,
            "predicted_u_mfr_new_m_s": predicted_velocity,
            "expected_smooth_mean_velocity_m_s": EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
            "velocity_ratio": velocity_ratio,
            "maximum_velocity_ratio_error": VELOCITY_RATIO_MAXIMUM_ERROR,
        }
        write_json(layout.qc / "mfr_eq_static_prediction_qc.json", prediction_qc)

        failures: list[dict[str, str]] = []
        if area_qc["status"] != "PASS":
            failures.append(
                {
                    "status": AREA_PROXY_FAILED,
                    "reason": f"Area proxy ratio is {area_qc['area_proxy_over_smooth_area']:.17g}",
                }
            )
        if normal_qc["status"] != "PASS":
            failures.append(
                {
                    "status": NORMALIND_NOT_CLEAN,
                    "reason": f"normalInd distribution is {normal_qc['distribution']}",
                }
            )
        if domain_qc["status"] != "PASS":
            failures.append(
                {
                    "status": DOMAIN_CHANGED,
                    "reason": "Fluid count changed by more than 0.5%, mesh is nonuniform, or fluid has multiple regions",
                }
            )
        if prediction_qc["status"] != "PASS":
            failures.append(
                {
                    "status": FAILED_STATUS,
                    "reason": f"Predicted velocity ratio is {velocity_ratio:.17g}",
                }
            )

        frozen_after = file_snapshot(frozen_paths)
        frozen_unchanged = frozen_before == frozen_after
        write_json(
            layout.qc / "frozen_source_read_only_sha_qc.json",
            {
                "status": "PASS" if frozen_unchanged else "FAIL",
                "source_frozen_files_unchanged": frozen_unchanged,
                "previous_rotated_geometry_modified": not frozen_unchanged,
                "before": frozen_before,
                "after": frozen_after,
            },
        )
        if not frozen_unchanged:
            failures.append(
                {"status": FAILED_STATUS, "reason": "A frozen source file changed"}
            )

        mesh_manifest = directory_manifest(layout.mesh)
        mesh_manifest.update(
            {
                "status": "PASS" if not failures else "PREFLIGHT_HARD_GATE_FAILED",
                "dx_m": EXPECTED_DX_M,
                "root_level": cube.level,
                "fluid_element_count": int(mesh_header["fluid_element_count"]),
                "boundary_element_count": int(mesh_header["boundary_element_count"]),
                "qval_element_count": int(
                    mesh_header["property_element_counts"].get("has qVal", 0)
                ),
                "inlet_d3q19_element_count": int(
                    mesh_preflight["d3q19_inlet"]["element_count"]
                ),
            }
        )
        write_json(layout.seeder / "mesh_manifest.json", mesh_manifest)

        final_status = PASS_STATUS if not failures else FAILED_STATUS
        final_next = PASS_NEXT if not failures else TOPOLOGY_FAILURE_NEXT
        summary.update(
            {
                "status": final_status,
                "next": final_next,
                "hard_gate_failures": failures,
                "first_failure_status": failures[0]["status"] if failures else None,
                "first_failure_reason": failures[0]["reason"] if failures else None,
                "production_pipeline_modified": False,
                "previous_rotated_geometry_modified": not frozen_unchanged,
                "numerical_inlet_plane_path": str(plane_path),
                "plane_z_m": plane_z_m,
                "plane_normal": [0.0, 0.0, 1.0],
                "rectangle_xy_bounds_m": rectangle["snapped_xy_bounds_m"],
                "expansion_margin_m": rectangle["requested_expansion_margin_m"],
                "expansion_margin_cells": RECTANGLE_MARGIN_CELLS,
                "grid_phase_error_over_dx": phase_error,
                "dx_m": EXPECTED_DX_M,
                "root_level": cube.level,
                "previous_fluid_cell_count": PREVIOUS_FLUID_CELL_COUNT,
                "new_fluid_cell_count": int(mesh_header["fluid_element_count"]),
                "fluid_count_relative_difference": domain_qc["relative_difference"],
                "connected_fluid_regions": region_count,
                "physical_smooth_inlet_area_m2": EXPECTED_SMOOTH_INLET_AREA_M2,
                "theoretical_area_over_dx_squared": EXPECTED_SMOOTH_INLET_AREA_M2
                / EXPECTED_DX_M**2,
                "new_inlet_d3q19_globbc_count": int(
                    mesh_preflight["d3q19_inlet"]["element_count"]
                ),
                "new_mfr_eq_area_proxy_m2": proxy_area,
                "area_proxy_over_smooth_area": area_qc[
                    "area_proxy_over_smooth_area"
                ],
                "normal_ind_distribution": normal_qc["distribution"],
                "negative_z_count": normal_qc["target_cardinal_count"],
                "negative_z_fraction": normal_qc[
                    "target_cardinal_normal_fraction"
                ],
                "diagonal_count": normal_qc["diagonal_normal_ind_count"],
                "other_direction_count": normal_qc["other_direction_count"],
                "predicted_u_mfr_new_m_s": predicted_velocity,
                "expected_smooth_mean_velocity_m_s": EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
                "velocity_ratio": velocity_ratio,
                "source_frozen_files_unchanged": frozen_unchanged,
                "old_oblique_baseline_status": OLD_BASELINE_STATUS,
                "previous_rotation_only_preflight_status": previous_manifest["status"],
                "grid_convergence": "NOT_RUN",
                "musubi_run_count": 0,
                "harvester_run_count": 0,
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
    except Exception as error:
        frozen_after = file_snapshot(frozen_paths)
        frozen_unchanged = frozen_before == frozen_after
        status = error.status if isinstance(error, FlowError) else FAILED_STATUS
        if status not in {GEOMETRY_UNSAFE, FAILED_STATUS}:
            status = FAILED_STATUS
        summary.update(
            {
                "status": status,
                "next": CURRENT_SEEDER_FAILURE_NEXT,
                "first_failure_reason": str(error),
                "previous_rotated_geometry_modified": not frozen_unchanged,
                "source_frozen_files_unchanged": frozen_unchanged,
                "musubi_run_count": 0,
                "harvester_run_count": 0,
                "grid_convergence": "NOT_RUN",
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(
            layout.qc / "frozen_source_read_only_sha_qc.json",
            {
                "status": "PASS" if frozen_unchanged else "FAIL",
                "source_frozen_files_unchanged": frozen_unchanged,
                "before": frozen_before,
                "after": frozen_after,
            },
        )
        write_json(manifest_path, summary)
        return summary


def finalize_failed_ideal_plane_preflight(
    project_root: Path, run_root: Path
) -> dict[str, Any]:
    """Complete truthful NOT_AVAILABLE QC after the consumed Seeder call failed."""

    root = Path(project_root).resolve()
    branch = _git_value(root, "branch", "--show-current")
    head = _git_value(root, "rev-parse", "HEAD")
    if branch != EXPECTED_BRANCH or head != EXPECTED_HEAD:
        raise FlowError(
            FAILED_STATUS,
            f"Expected {EXPECTED_BRANCH}@{EXPECTED_HEAD}, found {branch}@{head}",
        )
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    run = Path(run_root).resolve()
    if run.parent != config.paths.output_root.resolve() or not run.name.startswith(
        f"{OUTPUT_PREFIX}_"
    ):
        raise FlowError(FAILED_STATUS, f"Refusing non-ideal-plane run: {run}")
    manifest_path = run / "qc" / "ideal_plane_preflight_manifest.json"
    summary = read_json(manifest_path)
    if (
        summary.get("seeder_run_count") != 1
        or summary.get("seeder_return_code") != 134
        or summary.get("musubi_run_count") != 0
        or summary.get("harvester_run_count") != 0
    ):
        raise FlowError(
            FAILED_STATUS,
            "Run is not the consumed one-Seeder/no-solver format-failure attempt",
        )

    geometry_qc = read_json(run / "qc" / "ideal_inlet_plane_geometry_qc.json")
    previous_run = config.paths.output_root.resolve() / PREVIOUS_ROTATED_RUN
    previous_manifest = read_json(
        previous_run / "qc" / "axis_aligned_inlet_preflight_manifest.json"
    )
    previous_phase_qc = read_json(previous_run / "qc" / "inlet_grid_phase_qc.json")
    frozen_before = summary["frozen_files_before"]
    frozen_after = file_snapshot(Path(path) for path in frozen_before)
    frozen_unchanged = frozen_before == frozen_after
    unavailable_reason = (
        "Seeder 1.6.2 returned 134 while its binary-only tem_stlb reader tried "
        "to read the standards-valid ASCII two-triangle STL; no TreElm mesh was written"
    )

    source_contract = {
        "status": "SOURCE_PROVEN",
        "pinned_seeder_commit": "667109df6fafdcb39f4409e3f5d90f04d75cd33c",
        "tracked_source_clean_at_review": True,
        "failed_attempt_boundary": "ASCII_TWO_TRIANGLE_STL",
        "failed_attempt_stderr": str(run / "seeder" / "seeder_stderr.log"),
        "binary_stl_contract": {
            "source": "/home/lzy/apes-pinned/seeder_official/tem/source/tem_stlbIO_module.f90",
            "line_ranges": [[56, 87], [96, 124]],
            "finding": "Pinned STL reader is explicitly binary and reads single-precision triangle records",
        },
        "official_planar_canond_contract": {
            "source": "/home/lzy/apes-pinned/seeder_official/tem/source/shapes/tem_plane_module.fpp",
            "line_range": [133, 149],
            "finding": "Official finite plane uses canoND origin plus two vectors and is internally represented by two triangles",
        },
        "corrected_code_boundary": "OFFICIAL_CANOND_CARTESIAN_PLANE",
        "additional_seeder_call_after_correction": 0,
    }
    write_json(run / "qc" / "pinned_seeder_geometry_contract.json", source_contract)

    inlet_topology_qc = {
        "status": "NOT_AVAILABLE",
        "reason": unavailable_reason,
        "seeder_return_code": 134,
        "treelm_mesh_written": False,
        "total_inlet_d3q19_elements": None,
        "normal_ind_distribution": None,
        "negative_z_count": None,
        "negative_z_fraction": None,
        "diagonal_count": None,
        "other_direction_count": None,
        "third_parser_written": False,
    }
    domain_qc = {
        "status": "NOT_AVAILABLE",
        "reason": unavailable_reason,
        "previous_fluid_cell_count": PREVIOUS_FLUID_CELL_COUNT,
        "new_fluid_cell_count": None,
        "relative_difference": None,
        "maximum_relative_difference": FLUID_COUNT_MAXIMUM_RELATIVE_DIFFERENCE,
        "connected_fluid_regions": None,
        "single_connected_fluid_domain": None,
        "boundary_element_count": None,
        "qval_element_count": None,
    }
    prediction_qc = {
        "status": "NOT_AVAILABLE",
        "reason": unavailable_reason,
        "prediction_only_no_musubi_run": True,
        "physical_smooth_inlet_area_m2": EXPECTED_SMOOTH_INLET_AREA_M2,
        "theoretical_area_over_dx_squared": EXPECTED_SMOOTH_INLET_AREA_M2
        / EXPECTED_DX_M**2,
        "new_inlet_d3q19_globbc_count": None,
        "mfr_eq_area_proxy_m2": None,
        "area_proxy_over_smooth_area": None,
        "predicted_u_mfr_new_m_s": None,
        "expected_smooth_mean_velocity_m_s": EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
        "velocity_ratio": None,
    }
    write_json(run / "qc" / "inlet_topology_qc.json", inlet_topology_qc)
    write_json(run / "qc" / "fluid_domain_comparison_qc.json", domain_qc)
    write_json(run / "qc" / "mfr_eq_static_prediction_qc.json", prediction_qc)
    write_json(
        run / "seeder" / "mesh_manifest.json",
        {
            "status": "NOT_AVAILABLE",
            "reason": unavailable_reason,
            "file_count": 0,
            "total_bytes": 0,
            "seeder_return_code": 134,
        },
    )
    write_json(
        run / "qc" / "frozen_source_read_only_sha_qc.json",
        {
            "status": "PASS" if frozen_unchanged else "FAIL",
            "source_frozen_files_unchanged": frozen_unchanged,
            "previous_rotated_geometry_modified": not frozen_unchanged,
            "before": frozen_before,
            "after": frozen_after,
        },
    )

    rectangle = geometry_qc["rectangle"]
    plane = geometry_qc["plane"]
    summary.update(
        {
            "status": FAILED_STATUS,
            "next": FAILED_NEXT,
            "hard_gate_failures": [
                {"status": SEEDER_INPUT_FORMAT_FAILED, "reason": unavailable_reason}
            ],
            "first_failure_status": SEEDER_INPUT_FORMAT_FAILED,
            "first_failure_reason": unavailable_reason,
            "production_pipeline_modified": False,
            "previous_rotated_geometry_modified": not frozen_unchanged,
            "attempted_numerical_inlet_boundary_type": "ASCII_TWO_TRIANGLE_CARTESIAN_PLANAR_STL",
            "corrected_code_numerical_inlet_boundary_type": "OFFICIAL_CANOND_CARTESIAN_PLANE",
            "numerical_inlet_plane_path": str(
                run / "geometry" / "numerical_inlet_plane.stl"
            ),
            "plane_z_m": plane["plane_z_m"],
            "plane_normal": plane["plane_normal"],
            "rectangle_xy_bounds_m": rectangle["snapped_xy_bounds_m"],
            "expansion_margin_m": rectangle["requested_expansion_margin_m"],
            "expansion_margin_cells": RECTANGLE_MARGIN_CELLS,
            "grid_phase_error_over_dx": geometry_qc["grid_phase_error_over_dx"],
            "dx_m": EXPECTED_DX_M,
            "root_level": int(previous_phase_qc["root_level"]),
            "previous_fluid_cell_count": PREVIOUS_FLUID_CELL_COUNT,
            "new_fluid_cell_count": None,
            "fluid_count_relative_difference": None,
            "connected_fluid_regions": None,
            "physical_smooth_inlet_area_m2": EXPECTED_SMOOTH_INLET_AREA_M2,
            "theoretical_area_over_dx_squared": EXPECTED_SMOOTH_INLET_AREA_M2
            / EXPECTED_DX_M**2,
            "new_inlet_d3q19_globbc_count": None,
            "new_mfr_eq_area_proxy_m2": None,
            "area_proxy_over_smooth_area": None,
            "normal_ind_distribution": None,
            "negative_z_count": None,
            "negative_z_fraction": None,
            "diagonal_count": None,
            "other_direction_count": None,
            "predicted_u_mfr_new_m_s": None,
            "expected_smooth_mean_velocity_m_s": EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
            "velocity_ratio": None,
            "source_frozen_files_unchanged": frozen_unchanged,
            "old_oblique_baseline_status": OLD_BASELINE_STATUS,
            "previous_rotation_only_preflight_status": previous_manifest["status"],
            "grid_convergence": "NOT_RUN",
            "seeder_run_count": 1,
            "seeder_return_code": 134,
            "musubi_run_count": 0,
            "harvester_run_count": 0,
            "post_failure_additional_seeder_calls": 0,
            "completed_at": datetime.now().isoformat(),
        }
    )
    write_json(manifest_path, summary)
    return summary
