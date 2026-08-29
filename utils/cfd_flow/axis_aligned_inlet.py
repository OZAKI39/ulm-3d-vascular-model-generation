"""Build and preflight one rigidly rotated, axis-aligned inlet geometry.

This experimental path is deliberately separate from the production CFD
pipeline.  It applies one global rigid transform, runs Seeder at most once,
and never launches Musubi or harvesting.
"""

from __future__ import annotations

import math
import shutil
import subprocess
from dataclasses import asdict, dataclass, replace
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pyvista as pv
import trimesh

from .apes import (
    generate_seeder_lua,
    inspect_apes_environment,
    load_boundary_conditions,
    parse_mesh_header,
    run_wsl_tool,
)
from .config import load_cfd_flow_config
from .exact_link_flux import reconstruct_boundary
from .geometry import (
    BOUNDARY_LABELS,
    BoundingCube,
    find_seed_point,
    load_frozen_surface_partition,
    partition_surface,
)
from .io import FlowError, load_flow_inputs, read_json, sha256_file, write_json
from .port_flux_audit import (
    extract_boundary_property_indices,
    parse_bnd_header,
    parse_boundary_property_header,
    read_boundary_ids,
)
from .restart_decode import D3Q19_DIRECTIONS, read_treelm_elemlist


SOURCE_BRANCH = "codex/cfd-flow-musubi-recovery-20260828"
CURRENT_SYNCED_BASE_COMMIT = "b37d5391cdb186068411357a11943729347899d6"
FROZEN_SEEDER_RUN = "musubi_recovery_anchor003274_20260828_162530"
OLD_EXACT_AUDIT_RUN = "exact_link_flux_audit_anchor003274_20260829_105727"
OUTPUT_PREFIX = "axis_aligned_inlet_geometry_anchor003274"
REVISION = "AXIS_ALIGNED_INLET_RIGID_ROTATION_SEEDER_PREFLIGHT_V1"

MODEL_NAME = "AXIS_ALIGNED_INLET_BASELINE"
OLD_BASELINE_NAME = "OBLIQUE_INLET_DIAGNOSTIC_BASELINE"
OLD_BASELINE_STATUS = "MFR_EQ_OBLIQUE_PORT_TARGET_MISMATCH_CONFIRMED"
SUCCESS_STATUS = "CFD_FLOW_AXIS_ALIGNED_INLET_MESH_PREFLIGHT_PASS"
FAILED_STATUS = "CFD_FLOW_AXIS_ALIGNED_INLET_MESH_PREFLIGHT_FAILED"
SOURCE_NOT_PLANAR = "CFD_FLOW_AXIS_ALIGNED_INLET_SOURCE_CAP_NOT_PLANAR"
AREA_PROXY_FAILED = "CFD_FLOW_AXIS_ALIGNED_INLET_AREA_PROXY_FAILED"
NORMALIND_NOT_CLEAN = "CFD_FLOW_AXIS_ALIGNED_INLET_NORMALIND_NOT_CLEAN"
SUCCESS_NEXT = "RUN ONE NEW AXIS_ALIGNED MUSUBI BASELINE"
FAILED_NEXT = "REVIEW FIRST AXIS-ALIGNED INLET PREFLIGHT HARD-GATE FAILURE"

TARGET_INWARD_NORMAL = np.asarray((0.0, 0.0, -1.0), dtype=np.float64)
TARGET_OUTWARD_NORMAL = -TARGET_INWARD_NORMAL
TARGET_CARDINAL_INDEX_ZERO_BASED = 2
EXPECTED_DX_M = 2.0e-7
EXPECTED_Q_M3_S = 7.693508475538942e-16
EXPECTED_MASS_FLOW_KG_S = 8.124344950169123e-13
EXPECTED_DENSITY_KG_M3 = 1056.0
EXPECTED_KINEMATIC_VISCOSITY_M2_S = 3.27e-6
EXPECTED_BULK_VISCOSITY_M2_S = 2.18e-6
EXPECTED_GAUGE_PRESSURES_PA = (
    14.544978101274268,
    132.20454922317552,
    -13.700626673311461,
)
EXPECTED_PRESSURE_REFERENCE_PA = 23622.32012800001
EXPECTED_SMOOTH_INLET_AREA_M2 = 7.819752111687344e-12
EXPECTED_SMOOTH_MEAN_VELOCITY_M_S = 9.838558007536173e-5

GEOMETRY_RELATIVE_TOLERANCE = 1.0e-12
NORMAL_DOT_TOLERANCE = 1.0e-12
GRID_PHASE_TOLERANCE = 1.0e-10
AREA_PROXY_HARD_RELATIVE_ERROR = 0.05
AREA_PROXY_PREFERRED_RELATIVE_ERROR = 0.02
SOURCE_PLANARITY_TOLERANCE_IN_DX = 1.0e-4


@dataclass(frozen=True, slots=True)
class AxisAlignedLayout:
    root: Path
    geometry: Path
    geometry_solver_m: Path
    transform: Path
    seeder: Path
    seeder_mesh: Path
    qc: Path
    input: Path


def minimal_rotation_matrix(
    source: np.ndarray, target: np.ndarray
) -> tuple[np.ndarray, float]:
    """Return the shortest proper rotation mapping source onto target."""

    a = np.asarray(source, dtype=np.float64).reshape(3)
    b = np.asarray(target, dtype=np.float64).reshape(3)
    a /= np.linalg.norm(a)
    b /= np.linalg.norm(b)
    cosine = float(np.clip(np.dot(a, b), -1.0, 1.0))
    if cosine >= 1.0 - 1.0e-15:
        return np.eye(3), 0.0
    if cosine <= -1.0 + 1.0e-15:
        basis = np.eye(3)[int(np.argmin(np.abs(a)))]
        axis = np.cross(a, basis)
        axis /= np.linalg.norm(axis)
        rotation = -np.eye(3) + 2.0 * np.outer(axis, axis)
        return rotation, math.pi
    cross = np.cross(a, b)
    sine = float(np.linalg.norm(cross))
    skew = np.asarray(
        (
            (0.0, -cross[2], cross[1]),
            (cross[2], 0.0, -cross[0]),
            (-cross[1], cross[0], 0.0),
        )
    )
    rotation = np.eye(3) + skew + skew @ skew * ((1.0 - cosine) / sine**2)
    return rotation, math.acos(cosine)


def apply_rigid_transform(
    points: np.ndarray, rotation: np.ndarray, pivot: np.ndarray
) -> np.ndarray:
    values = np.asarray(points, dtype=np.float64)
    matrix = np.asarray(rotation, dtype=np.float64).reshape(3, 3)
    center = np.asarray(pivot, dtype=np.float64).reshape(3)
    return (matrix @ (values - center).T).T + center


def homogeneous_transform(
    rotation: np.ndarray, pivot: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    matrix = np.asarray(rotation, dtype=np.float64).reshape(3, 3)
    center = np.asarray(pivot, dtype=np.float64).reshape(3)
    forward = np.eye(4)
    forward[:3, :3] = matrix
    forward[:3, 3] = center - matrix @ center
    inverse = np.eye(4)
    inverse[:3, :3] = matrix.T
    inverse[:3, 3] = center - matrix.T @ center
    return forward, inverse


def inlet_plane_deviation(
    points: np.ndarray,
    faces: np.ndarray,
    face_indices: np.ndarray,
    centroid: np.ndarray,
    normal: np.ndarray,
) -> float:
    vertex_indices = np.unique(
        np.asarray(faces, dtype=np.int64)[face_indices].reshape(-1)
    )
    offsets = np.asarray(points, dtype=np.float64)[vertex_indices] - np.asarray(
        centroid, dtype=np.float64
    )
    return float(np.max(np.abs(offsets @ np.asarray(normal, dtype=np.float64))))


def compute_phase_aligned_bounding_cube(
    bounds_um: np.ndarray,
    inlet_plane_z_um: float,
    *,
    dx_m: float,
    margin_cells: int,
) -> tuple[BoundingCube, dict[str, Any]]:
    """Find the smallest 2^level cube whose z-origin makes the inlet a grid face."""

    bounds_m = np.asarray(bounds_um, dtype=np.float64).reshape(3, 2) * 1.0e-6
    lower = bounds_m[:, 0]
    upper = bounds_m[:, 1]
    plane_z_m = float(inlet_plane_z_um) * 1.0e-6
    margin_m = int(margin_cells) * float(dx_m)
    required = float(np.max(upper - lower + 2.0 * margin_m))
    first_level = max(0, int(math.ceil(math.log2(required / dx_m))))
    for level in range(first_level, first_level + 16):
        side = 2**level * float(dx_m)
        origin = (lower + upper) / 2.0 - side / 2.0
        minimum_k = math.ceil((plane_z_m - (lower[2] - margin_m)) / dx_m - 1.0e-12)
        maximum_k = math.floor(
            (plane_z_m + side - (upper[2] + margin_m)) / dx_m + 1.0e-12
        )
        if minimum_k > maximum_k:
            continue
        centered_origin_z = (lower[2] + upper[2]) / 2.0 - side / 2.0
        centered_k = int(round((plane_z_m - centered_origin_z) / dx_m))
        selected_k = min(max(centered_k, minimum_k), maximum_k)
        origin[2] = plane_z_m - selected_k * dx_m
        lower_margin = lower - origin
        upper_margin = origin + side - upper
        minimum_margin_cells = float(
            np.min(np.minimum(lower_margin, upper_margin)) / dx_m
        )
        if minimum_margin_cells + 1.0e-10 < margin_cells:
            continue
        phase = (plane_z_m - origin[2]) / dx_m
        phase_error = abs(phase - round(phase))
        cube = BoundingCube(origin, side, level, 2**level, minimum_margin_cells)
        return cube, {
            "inlet_plane_z_m": plane_z_m,
            "phase_coordinate": phase,
            "nearest_integer_face_index": int(round(phase)),
            "inlet_plane_grid_phase_error_over_dx": phase_error,
            "phase_status": "PASS" if phase_error <= GRID_PHASE_TOLERANCE else "FAIL",
        }
    raise FlowError(
        FAILED_STATUS, "No phase-aligned bounding cube satisfies the required margin"
    )


def relative_difference(before: float, after: float) -> float:
    return abs(float(after) - float(before)) / max(
        abs(float(before)), np.finfo(float).tiny
    )


def pairwise_distances(points: np.ndarray) -> np.ndarray:
    values = np.asarray(points, dtype=np.float64)
    return np.asarray(
        [
            np.linalg.norm(values[first] - values[second])
            for first in range(len(values))
            for second in range(first + 1, len(values))
        ]
    )


def cardinal_normal_gate(normal_indices: np.ndarray) -> dict[str, Any]:
    values = np.asarray(normal_indices, dtype=np.int64).reshape(-1)
    target_count = int(np.count_nonzero(values == TARGET_CARDINAL_INDEX_ZERO_BASED))
    diagonal = np.linalg.norm(D3Q19_DIRECTIONS[values].astype(np.float64), axis=1) > 1.0
    unique, counts = np.unique(values, return_counts=True)
    distribution = [
        {
            "direction_index": int(index + 1),
            "cx": int(D3Q19_DIRECTIONS[index, 0]),
            "cy": int(D3Q19_DIRECTIONS[index, 1]),
            "cz": int(D3Q19_DIRECTIONS[index, 2]),
            "element_count": int(count),
        }
        for index, count in zip(unique, counts, strict=True)
    ]
    fraction = target_count / len(values) if len(values) else 0.0
    return {
        "status": "PASS" if len(values) and target_count == len(values) else "FAIL",
        "distribution": distribution,
        "target_direction_index": 3,
        "target_direction": [0, 0, -1],
        "target_cardinal_count": target_count,
        "target_cardinal_normal_fraction": fraction,
        "diagonal_normal_ind_count": int(np.count_nonzero(diagonal)),
    }


def area_proxy_qc(
    element_count: int, smooth_area_m2: float, dx_m: float
) -> dict[str, Any]:
    proxy = int(element_count) * float(dx_m) ** 2
    ratio = proxy / float(smooth_area_m2)
    return {
        "status": "PASS"
        if abs(ratio - 1.0) <= AREA_PROXY_HARD_RELATIVE_ERROR
        else "FAIL",
        "preferred_status": "PASS"
        if abs(ratio - 1.0) <= AREA_PROXY_PREFERRED_RELATIVE_ERROR
        else "OUTSIDE_PREFERRED_2_PERCENT",
        "d3q19_inlet_element_count": int(element_count),
        "dx_m": float(dx_m),
        "smooth_inlet_area_m2": float(smooth_area_m2),
        "smooth_area_over_dx_squared": float(smooth_area_m2) / float(dx_m) ** 2,
        "mfr_eq_area_proxy_m2": proxy,
        "area_proxy_over_smooth_area": ratio,
        "absolute_ratio_error": abs(ratio - 1.0),
        "hard_maximum_ratio_error": AREA_PROXY_HARD_RELATIVE_ERROR,
    }


def _create_layout(output_root: Path) -> AxisAlignedLayout:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    root = Path(output_root) / f"{OUTPUT_PREFIX}_{stamp}"
    if root.exists():
        raise FlowError(FAILED_STATUS, f"Output already exists: {root}")
    directories = {
        name: root / name for name in ("geometry", "transform", "seeder", "qc", "input")
    }
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=False)
    geometry_solver_m = directories["geometry"] / "geometry_solver_m"
    geometry_solver_m.mkdir()
    seeder_mesh = directories["seeder"] / "mesh"
    seeder_mesh.mkdir()
    return AxisAlignedLayout(
        root,
        directories["geometry"],
        geometry_solver_m,
        directories["transform"],
        directories["seeder"],
        seeder_mesh,
        directories["qc"],
        directories["input"],
    )


def _git_value(repository: Path, *arguments: str) -> str:
    process = subprocess.run(
        ["git", "-C", str(repository), *arguments],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    return process.stdout.strip()


def _file_manifest(paths: Iterable[Path]) -> dict[str, Any]:
    return {
        str(path.resolve()): {"bytes": path.stat().st_size, "sha256": sha256_file(path)}
        for path in paths
    }


def _directory_manifest(directory: Path) -> dict[str, Any]:
    paths = sorted(path for path in directory.rglob("*") if path.is_file())
    return {
        "root": str(directory),
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


def _write_rotated_tagged_surface(
    source_path: Path, destination: Path, rotated_points_um: np.ndarray
) -> None:
    surface = pv.read(source_path).triangulate()
    surface.points = np.asarray(rotated_points_um, dtype=np.float64)
    surface.save(destination, binary=True)


def _write_rotated_full_stl(
    points_um: np.ndarray, faces: np.ndarray, destination: Path
) -> None:
    mesh = trimesh.Trimesh(
        np.asarray(points_um) * 1.0e-6, np.asarray(faces), process=False
    )
    mesh.export(destination, file_type="stl")


def _geometry_qc(
    source: Any, rotated: Any, rotation: np.ndarray, source_plane_deviation_um: float
) -> dict[str, Any]:
    source_patches = {patch.label: patch for patch in source.patches}
    rotated_patches = {patch.label: patch for patch in rotated.patches}
    labels = BOUNDARY_LABELS[1:]
    port_areas = {
        label: {
            "before_um2": source_patches[label].area_um2,
            "after_um2": rotated_patches[label].area_um2,
            "relative_difference": relative_difference(
                source_patches[label].area_um2, rotated_patches[label].area_um2
            ),
        }
        for label in labels
    }
    before_centers = np.asarray([source_patches[label].center_um for label in labels])
    after_centers = np.asarray([rotated_patches[label].center_um for label in labels])
    before_distances = pairwise_distances(before_centers)
    after_distances = pairwise_distances(after_centers)
    distance_errors = np.abs(after_distances - before_distances) / before_distances
    source_total_area = float(source.mesh_um.area)
    rotated_total_area = float(rotated.mesh_um.area)
    inlet_after = rotated_patches["inlet"]
    inlet_plane_after_um = inlet_plane_deviation(
        rotated.points_um,
        rotated.faces,
        inlet_after.face_indices,
        inlet_after.center_um,
        TARGET_OUTWARD_NORMAL,
    )
    checks = {
        "triangle_count_exact": len(source.faces) == len(rotated.faces),
        "faces_exact": bool(np.array_equal(source.faces, rotated.faces)),
        "cell_entity_ids_exact": bool(
            np.array_equal(source.entity_ids, rotated.entity_ids)
        ),
        "watertight_same": bool(
            source.mesh_um.is_watertight == rotated.mesh_um.is_watertight
        ),
        "total_area_relative_difference": relative_difference(
            source_total_area, rotated_total_area
        )
        <= GEOMETRY_RELATIVE_TOLERANCE,
        "all_port_area_relative_differences": all(
            item["relative_difference"] <= GEOMETRY_RELATIVE_TOLERANCE
            for item in port_areas.values()
        ),
        "pairwise_port_center_distances": bool(
            np.all(distance_errors <= GEOMETRY_RELATIVE_TOLERANCE)
        ),
        "rotation_orthogonal": bool(
            np.allclose(rotation.T @ rotation, np.eye(3), rtol=0.0, atol=1.0e-14)
        ),
        "rotation_determinant_one": math.isclose(
            float(np.linalg.det(rotation)), 1.0, rel_tol=0.0, abs_tol=1.0e-14
        ),
        "inlet_outward_normal_plus_z": float(
            np.dot(inlet_after.outward_normal, TARGET_OUTWARD_NORMAL)
        )
        >= 1.0 - NORMAL_DOT_TOLERANCE,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "transform_kind": "GLOBAL_RIGID_ROTATION_ONLY",
        "scale": 1.0,
        "shear": False,
        "local_deformation": False,
        "remeshing": False,
        "triangle_count_before": int(len(source.faces)),
        "triangle_count_after": int(len(rotated.faces)),
        "watertight_before": bool(source.mesh_um.is_watertight),
        "watertight_after": bool(rotated.mesh_um.is_watertight),
        "total_area_before_um2": source_total_area,
        "total_area_after_um2": rotated_total_area,
        "total_area_relative_difference": relative_difference(
            source_total_area, rotated_total_area
        ),
        "port_areas": port_areas,
        "pairwise_port_center_distances_before_um": before_distances.tolist(),
        "pairwise_port_center_distances_after_um": after_distances.tolist(),
        "maximum_pairwise_distance_relative_difference": float(np.max(distance_errors)),
        "source_inlet_plane_maximum_deviation_m": source_plane_deviation_um * 1.0e-6,
        "rotated_inlet_plane_maximum_deviation_m": inlet_plane_after_um * 1.0e-6,
        "inlet_outward_normal_after": inlet_after.outward_normal.tolist(),
        "inlet_outward_normal_plus_z_dot": float(
            np.dot(inlet_after.outward_normal, TARGET_OUTWARD_NORMAL)
        ),
        "cell_entity_ids_unchanged": bool(
            np.array_equal(source.entity_ids, rotated.entity_ids)
        ),
        "port_ids_unchanged": True,
        "wall_and_cap_labels_unchanged": True,
    }


def _bc_physics_qc(config: Any, boundary_conditions: Any) -> dict[str, Any]:
    pressure_factor = EXPECTED_DENSITY_KG_M3 * EXPECTED_DX_M**2 / (2.44140625e-8) ** 2
    pressure_reference = pressure_factor / 3.0
    bulk = (2.0 / 3.0) * boundary_conditions.kinematic_viscosity_m2_s
    checks = {
        "inlet_q": boundary_conditions.inlet_flow_m3_s == EXPECTED_Q_M3_S,
        "inlet_mass_flow": math.isclose(
            boundary_conditions.density_kg_m3 * boundary_conditions.inlet_flow_m3_s,
            EXPECTED_MASS_FLOW_KG_S,
            rel_tol=1.0e-15,
        ),
        "outlet_gauge_pressures": tuple(boundary_conditions.outlet_gauge_pressures_pa)
        == EXPECTED_GAUGE_PRESSURES_PA,
        "pressure_reference": math.isclose(
            pressure_reference,
            EXPECTED_PRESSURE_REFERENCE_PA,
            rel_tol=0.0,
            abs_tol=1.0e-9,
        ),
        "density": boundary_conditions.density_kg_m3 == EXPECTED_DENSITY_KG_M3,
        "kinematic_viscosity": boundary_conditions.kinematic_viscosity_m2_s
        == EXPECTED_KINEMATIC_VISCOSITY_M2_S,
        "bulk_viscosity": math.isclose(
            bulk, EXPECTED_BULK_VISCOSITY_M2_S, rel_tol=0.0, abs_tol=1.0e-20
        ),
        "dx": config.mesh.dx_target_m == EXPECTED_DX_M,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "inlet_target_q_m3_s": boundary_conditions.inlet_flow_m3_s,
        "inlet_target_mass_flow_kg_s": boundary_conditions.density_kg_m3
        * boundary_conditions.inlet_flow_m3_s,
        "outlet_gauge_pressures_pa": list(
            boundary_conditions.outlet_gauge_pressures_pa
        ),
        "pressure_reference_pa": pressure_reference,
        "density_kg_m3": boundary_conditions.density_kg_m3,
        "kinematic_viscosity_m2_s": boundary_conditions.kinematic_viscosity_m2_s,
        "bulk_viscosity_m2_s": bulk,
        "dx_m": config.mesh.dx_target_m,
    }


def _mesh_boundary_preflight(
    mesh_dir: Path, smooth_area_m2: float, dx_m: float
) -> dict[str, Any]:
    mesh = parse_mesh_header(mesh_dir)
    header_path = Path(mesh["header"])
    property_header = parse_boundary_property_header(
        header_path.read_text(encoding="utf-8")
    )
    boundary_header = parse_bnd_header(
        (mesh_dir / "bnd.lua").read_text(encoding="utf-8")
    )
    tree_ids, property_bits, elemlist_contract = read_treelm_elemlist(
        mesh_dir / "elemlist.lsb",
        n_elems=int(mesh["fluid_element_count"]),
    )
    property_indices = extract_boundary_property_indices(
        property_bits, property_header.bit_position
    )
    if len(property_indices) != property_header.element_count:
        raise FlowError(
            FAILED_STATUS, "Boundary property mapping count differs from header"
        )
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
    normal_qc = cardinal_normal_gate(inlet.normal_indices)
    proxy_qc = area_proxy_qc(len(inlet.cell_indices), smooth_area_m2, dx_m)
    return {
        "mesh_summary": mesh,
        "property_header": asdict(property_header),
        "boundary_header": asdict(boundary_header),
        "elemlist_contract": elemlist_contract,
        "tree_id_count": int(len(tree_ids)),
        "d3q19_inlet": {
            "element_count": int(len(inlet.cell_indices)),
            "modified_link_count": int(np.count_nonzero(inlet.incoming_masks)),
            "normal_ind": normal_qc,
            "area_proxy": proxy_qc,
        },
    }


def run_axis_aligned_inlet_preflight(project_root: Path) -> dict[str, Any]:
    """Create the rotated geometry and consume at most one Seeder call."""

    root = Path(project_root).resolve()
    branch = _git_value(root, "branch", "--show-current")
    head = _git_value(root, "rev-parse", "HEAD")
    if branch != SOURCE_BRANCH or head != CURRENT_SYNCED_BASE_COMMIT:
        raise FlowError(
            FAILED_STATUS,
            f"Expected {SOURCE_BRANCH}@{CURRENT_SYNCED_BASE_COMMIT}, found {branch}@{head}",
        )
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    inputs = load_flow_inputs(config.paths.source_surface_run)
    frozen_patch_dir = (
        config.paths.output_root / FROZEN_SEEDER_RUN / "geometry" / "geometry_solver_m"
    )
    frozen_paths = (
        inputs.tagged_surface_vtp,
        inputs.meter_surface_stl,
        *(frozen_patch_dir / f"{label}.stl" for label in BOUNDARY_LABELS),
    )
    frozen_before = _file_manifest(frozen_paths)
    layout = _create_layout(config.paths.output_root)
    manifest_path = layout.qc / "axis_aligned_inlet_preflight_manifest.json"
    summary: dict[str, Any] = {
        "status": FAILED_STATUS,
        "next": FAILED_NEXT,
        "revision": REVISION,
        "model_name": MODEL_NAME,
        "old_baseline_name": OLD_BASELINE_NAME,
        "run_root": str(layout.root),
        "branch": branch,
        "actual_head": head,
        "production_pipeline_modified": False,
        "source_frozen_surface_modified": False,
        "seeder_run_count": 0,
        "seeder_return_code": "NOT_RUN",
        "musubi_run_count": 0,
        "harvester_run_count": 0,
        "grid_convergence": "NOT_RUN",
        "proteus_runtime": "NOT_RUN",
        "microbubble": "NOT_RUN",
        "ulm": "NOT_RUN",
        "frozen_files_before": frozen_before,
        "started_at": datetime.now().isoformat(),
    }
    write_json(manifest_path, summary)
    first_failure: str | None = None
    try:
        old_audit = (
            config.paths.output_root
            / OLD_EXACT_AUDIT_RUN
            / "qc"
            / "exact_link_flux_audit_manifest.json"
        )
        old_audit_manifest = read_json(old_audit)
        if old_audit_manifest.get("status") != OLD_BASELINE_STATUS:
            raise FlowError(
                FAILED_STATUS, "Old oblique diagnostic baseline status changed"
            )
        write_json(
            layout.qc / "old_oblique_baseline_reference.json",
            {
                "status": OLD_BASELINE_STATUS,
                "model_name": OLD_BASELINE_NAME,
                "read_only": True,
                "audit_manifest": str(old_audit),
                "audit_manifest_sha256": sha256_file(old_audit),
            },
        )

        source_partition = load_frozen_surface_partition(inputs, frozen_patch_dir)
        source_inlet = source_partition.patch("inlet")
        source_plane_deviation_um = inlet_plane_deviation(
            source_partition.points_um,
            source_partition.faces,
            source_inlet.face_indices,
            source_inlet.center_um,
            source_inlet.outward_normal,
        )
        planarity_tolerance_m = SOURCE_PLANARITY_TOLERANCE_IN_DX * EXPECTED_DX_M
        if source_plane_deviation_um * 1.0e-6 > planarity_tolerance_m:
            raise FlowError(
                SOURCE_NOT_PLANAR,
                f"Source cap deviation {source_plane_deviation_um * 1.0e-6:.17g} m exceeds {planarity_tolerance_m:.17g} m",
            )

        source_inward = -source_inlet.outward_normal
        rotation, angle_rad = minimal_rotation_matrix(
            source_inward, TARGET_INWARD_NORMAL
        )
        pivot_um = source_inlet.center_um.copy()
        rotated_points_um = apply_rigid_transform(
            source_partition.points_um, rotation, pivot_um
        )
        pivot_m = pivot_um * 1.0e-6
        forward_m, inverse_m = homogeneous_transform(rotation, pivot_m)
        forward_um, inverse_um = homogeneous_transform(rotation, pivot_um)
        transform_path = layout.transform / "anatomical_to_cfd_transform.json"
        transform = {
            "status": "PASS",
            "transform_kind": "GLOBAL_MINIMUM_RIGID_ROTATION",
            "source_coordinate_system": "anatomical_fMOST",
            "target_coordinate_system": "axis_aligned_CFD",
            "rotation_target": {
                "inlet_inward": [0.0, 0.0, -1.0],
                "inlet_outward": [0.0, 0.0, 1.0],
            },
            "source_inlet_outward_normal": source_inlet.outward_normal.tolist(),
            "source_inlet_inward_normal": source_inward.tolist(),
            "rotation_angle_rad": angle_rad,
            "rotation_angle_deg": math.degrees(angle_rad),
            "rotation_matrix_3x3": rotation.tolist(),
            "pivot_um": pivot_um.tolist(),
            "pivot_m": pivot_m.tolist(),
            "forward_homogeneous_transform_4x4": forward_m.tolist(),
            "inverse_homogeneous_transform_4x4": inverse_m.tolist(),
            "homogeneous_transform_length_unit": "m",
            "forward_homogeneous_transform_4x4_um": forward_um.tolist(),
            "inverse_homogeneous_transform_4x4_um": inverse_um.tolist(),
            "inverse_usage": "Apply inverse_homogeneous_transform_4x4 to CFD/bubble/ULM coordinates in meters to recover anatomical/fMOST coordinates.",
            "scale": 1.0,
            "shear": False,
            "local_deformation": False,
        }
        write_json(transform_path, transform)

        rotated_vtp = layout.geometry / "cfd_surface_axis_aligned_inlet_um.vtp"
        rotated_stl = layout.geometry / "cfd_surface_axis_aligned_inlet_m.stl"
        _write_rotated_tagged_surface(
            inputs.tagged_surface_vtp, rotated_vtp, rotated_points_um
        )
        _write_rotated_full_stl(rotated_points_um, source_partition.faces, rotated_stl)
        rotated_inputs = replace(
            inputs,
            run_root=layout.root,
            tagged_surface_vtp=rotated_vtp,
            meter_surface_stl=rotated_stl,
        )
        rotated_partition = partition_surface(rotated_inputs, layout.geometry_solver_m)
        geometry_qc = _geometry_qc(
            source_partition, rotated_partition, rotation, source_plane_deviation_um
        )
        write_json(layout.qc / "geometry_rigid_transform_qc.json", geometry_qc)
        if geometry_qc["status"] != "PASS":
            raise FlowError(
                FAILED_STATUS, "Rigid-transform geometry invariance QC failed"
            )

        boundary_conditions = load_boundary_conditions(inputs.boundary_conditions)
        bc_qc = _bc_physics_qc(config, boundary_conditions)
        write_json(layout.qc / "bc_physics_invariance_qc.json", bc_qc)
        if bc_qc["status"] != "PASS":
            raise FlowError(FAILED_STATUS, "Frozen BC/physics invariance QC failed")

        bounds_um = np.column_stack(
            (rotated_points_um.min(axis=0), rotated_points_um.max(axis=0))
        )
        rotated_inlet = rotated_partition.patch("inlet")
        cube, phase_qc = compute_phase_aligned_bounding_cube(
            bounds_um,
            float(rotated_inlet.center_um[2]),
            dx_m=EXPECTED_DX_M,
            margin_cells=config.mesh.bounding_margin_cells,
        )
        phase_qc.update(
            {
                "status": phase_qc["phase_status"],
                "bounding_cube_origin_m": cube.origin_m.tolist(),
                "bounding_cube_side_m": cube.side_m,
                "root_level": cube.level,
                "cells_per_axis": cube.cells_per_axis,
                "minimum_margin_cells": cube.margin_cells_minimum,
                "dx_m": EXPECTED_DX_M,
                "geometry_moved_for_phase_alignment": False,
            }
        )
        write_json(layout.qc / "inlet_grid_phase_qc.json", phase_qc)
        if phase_qc["status"] != "PASS":
            raise FlowError(FAILED_STATUS, "Inlet plane grid-phase alignment failed")

        seed = find_seed_point(rotated_partition)
        seeder_lua = generate_seeder_lua(rotated_partition, seed, cube)
        (layout.seeder / "seeder.lua").write_text(seeder_lua, encoding="utf-8")
        shutil.copy2(
            inputs.boundary_conditions_json,
            layout.input / "source_boundary_conditions.json",
        )
        shutil.copy2(
            inputs.boundary_manifest_csv, layout.input / "source_boundary_manifest.csv"
        )
        write_json(
            layout.input / "source_provenance.json",
            {
                "frozen_source_run": str(inputs.run_root),
                "frozen_files": frozen_before,
                "rotated_tagged_vtp": str(rotated_vtp),
                "rotated_meter_stl": str(rotated_stl),
                "transform": str(transform_path),
                "source_geometry_read_only": True,
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
        predicted_fluid = int(
            math.ceil(
                abs(rotated_partition.mesh_um.volume) / config.mesh.dx_target_um**3
            )
        )
        predicted_ram = (
            predicted_fluid * config.resources.estimated_bytes_per_fluid_cell
        )
        ram_limit = int(
            environment.available_ram_bytes
            * config.resources.maximum_available_ram_fraction
        )
        resource_preflight = {
            "status": "PASS" if predicted_ram <= ram_limit else "FAIL",
            "lumen_volume_um3": float(abs(rotated_partition.mesh_um.volume)),
            "predicted_fluid_element_count": predicted_fluid,
            "estimated_ram_bytes": predicted_ram,
            "available_ram_bytes": environment.available_ram_bytes,
            "maximum_allowed_ram_bytes": ram_limit,
            "maximum_available_ram_fraction": config.resources.maximum_available_ram_fraction,
            "estimated_bytes_per_fluid_cell": config.resources.estimated_bytes_per_fluid_cell,
        }
        write_json(layout.qc / "resource_preflight_qc.json", resource_preflight)
        if resource_preflight["status"] != "PASS":
            raise FlowError(
                FAILED_STATUS, "Predicted mesh RAM exceeds 60% of available RAM"
            )

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
            raise FlowError(
                FAILED_STATUS, "Generated Seeder Lua failed syntax preflight"
            )

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
            raise FlowError(
                FAILED_STATUS, f"Seeder return code {seeder_run.returncode}"
            )

        mesh_preflight = _mesh_boundary_preflight(
            layout.seeder_mesh,
            rotated_inlet.area_um2 * 1.0e-12,
            EXPECTED_DX_M,
        )
        mesh_summary = mesh_preflight["mesh_summary"]
        actual_ram = (
            int(mesh_summary["fluid_element_count"])
            * config.resources.estimated_bytes_per_fluid_cell
        )
        mesh_preflight["actual_estimated_ram_bytes"] = actual_ram
        mesh_preflight["available_ram_bytes"] = environment.available_ram_bytes
        mesh_preflight["maximum_allowed_ram_bytes"] = ram_limit
        mesh_preflight["ram_status"] = "PASS" if actual_ram <= ram_limit else "FAIL"
        mesh_preflight["dx_m"] = EXPECTED_DX_M
        mesh_preflight["root_level"] = cube.level
        mesh_preflight["uniform_level_matches_root"] = (
            mesh_summary["minimum_level"] == cube.level
            and mesh_summary["maximum_level"] == cube.level
        )
        write_json(layout.qc / "seeder_mesh_qc.json", mesh_preflight)
        area_qc = mesh_preflight["d3q19_inlet"]["area_proxy"]
        normal_qc = mesh_preflight["d3q19_inlet"]["normal_ind"]
        hard_gate_failures: list[dict[str, str]] = []
        if area_qc["status"] != "PASS":
            hard_gate_failures.append(
                {
                    "status": AREA_PROXY_FAILED,
                    "reason": f"Area proxy ratio is {area_qc['area_proxy_over_smooth_area']:.17g}",
                }
            )
        if normal_qc["status"] != "PASS":
            hard_gate_failures.append(
                {
                    "status": NORMALIND_NOT_CLEAN,
                    "reason": f"normalInd distribution is {normal_qc['distribution']}",
                }
            )
        if actual_ram > ram_limit:
            hard_gate_failures.append(
                {
                    "status": FAILED_STATUS,
                    "reason": "Actual Seeder mesh exceeds 60% RAM gate",
                }
            )
        if not mesh_preflight["uniform_level_matches_root"]:
            hard_gate_failures.append(
                {
                    "status": FAILED_STATUS,
                    "reason": "Seeder mesh is not uniform at the selected root level",
                }
            )

        velocity_new = EXPECTED_MASS_FLOW_KG_S / (
            EXPECTED_DENSITY_KG_M3 * area_qc["mfr_eq_area_proxy_m2"]
        )
        velocity_qc = {
            "status": "PASS",
            "predicted_u_mfr_new_m_s": velocity_new,
            "smooth_expected_mean_velocity_m_s": EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
            "velocity_ratio": velocity_new / EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
            "prediction_only_no_cfd_run": True,
        }
        write_json(layout.qc / "mfr_eq_static_velocity_prediction_qc.json", velocity_qc)
        mesh_manifest = _directory_manifest(layout.seeder_mesh)
        mesh_manifest.update(
            {
                "status": (
                    "AXIS_ALIGNED_INLET_SEEDER_MESH_PREFLIGHT_PASS"
                    if not hard_gate_failures
                    else "AXIS_ALIGNED_INLET_SEEDER_MESH_PREFLIGHT_FAILED"
                ),
                "source_rotated_surface_sha256": sha256_file(rotated_vtp),
                "actual_dx_m": EXPECTED_DX_M,
                "root_level": cube.level,
                "fluid_element_count": int(mesh_summary["fluid_element_count"]),
                "boundary_element_count": int(mesh_summary["boundary_element_count"]),
                "qval_element_count": int(
                    mesh_summary["property_element_counts"].get("has qVal", 0)
                ),
                "d3q19_inlet_element_count": area_qc["d3q19_inlet_element_count"],
            }
        )
        write_json(layout.seeder / "mesh_manifest.json", mesh_manifest)

        frozen_after = _file_manifest(frozen_paths)
        frozen_unchanged = frozen_before == frozen_after
        write_json(
            layout.qc / "frozen_source_read_only_sha_qc.json",
            {
                "status": "PASS" if frozen_unchanged else "FAIL",
                "source_frozen_surface_modified": not frozen_unchanged,
                "before": frozen_before,
                "after": frozen_after,
            },
        )
        if not frozen_unchanged:
            hard_gate_failures.append(
                {
                    "status": FAILED_STATUS,
                    "reason": "Frozen source surface or patch changed",
                }
            )

        final_status = SUCCESS_STATUS if not hard_gate_failures else FAILED_STATUS
        final_next = SUCCESS_NEXT if not hard_gate_failures else FAILED_NEXT

        summary.update(
            {
                "status": final_status,
                "next": final_next,
                "hard_gate_failures": hard_gate_failures,
                "first_failure_status": (
                    hard_gate_failures[0]["status"] if hard_gate_failures else None
                ),
                "first_failure_reason": (
                    hard_gate_failures[0]["reason"] if hard_gate_failures else None
                ),
                "old_oblique_baseline_status": OLD_BASELINE_STATUS,
                "new_derived_geometry_path": str(layout.geometry),
                "transform_path": str(transform_path),
                "rotation_target": "INLET_INWARD_TO_NEGATIVE_Z",
                "rotation_angle_rad": angle_rad,
                "rotation_angle_deg": math.degrees(angle_rad),
                "rotation_matrix_3x3": rotation.tolist(),
                "geometry_qc": geometry_qc,
                "grid_phase_qc": phase_qc,
                "dx_m": EXPECTED_DX_M,
                "root_level": cube.level,
                "fluid_element_count": int(mesh_summary["fluid_element_count"]),
                "boundary_element_count": int(mesh_summary["boundary_element_count"]),
                "qval_element_count": int(
                    mesh_summary["property_element_counts"].get("has qVal", 0)
                ),
                "new_inlet_d3q19_globbc_count": area_qc["d3q19_inlet_element_count"],
                "smooth_inlet_area_m2": rotated_inlet.area_um2 * 1.0e-12,
                "mfr_eq_area_proxy_m2": area_qc["mfr_eq_area_proxy_m2"],
                "area_proxy_over_smooth_area": area_qc["area_proxy_over_smooth_area"],
                "normal_ind_distribution": normal_qc["distribution"],
                "target_cardinal_normal_fraction": normal_qc[
                    "target_cardinal_normal_fraction"
                ],
                "diagonal_normal_ind_count": normal_qc["diagonal_normal_ind_count"],
                "predicted_u_mfr_new_m_s": velocity_new,
                "smooth_expected_mean_velocity_m_s": EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
                "velocity_ratio": velocity_new / EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
                "ram_estimate_bytes": actual_ram,
                "ram_limit_bytes": ram_limit,
                "source_frozen_surface_modified": False,
                "frozen_files_after": frozen_after,
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
    except Exception as error:
        first_failure = str(error)
        status = error.status if isinstance(error, FlowError) else FAILED_STATUS
        if status not in {SOURCE_NOT_PLANAR, AREA_PROXY_FAILED, NORMALIND_NOT_CLEAN}:
            status = FAILED_STATUS
        frozen_after = _file_manifest(frozen_paths)
        summary.update(
            {
                "status": status,
                "next": FAILED_NEXT,
                "first_failure_reason": first_failure,
                "source_frozen_surface_modified": frozen_before != frozen_after,
                "frozen_files_after": frozen_after,
                "musubi_run_count": 0,
                "harvester_run_count": 0,
                "grid_convergence": "NOT_RUN",
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary


def finalize_existing_axis_aligned_preflight(
    project_root: Path,
    run_root: Path,
) -> dict[str, Any]:
    """Finish QC from an already completed single Seeder run without rerunning it."""

    root = Path(project_root).resolve()
    branch = _git_value(root, "branch", "--show-current")
    head = _git_value(root, "rev-parse", "HEAD")
    if branch != SOURCE_BRANCH or head != CURRENT_SYNCED_BASE_COMMIT:
        raise FlowError(
            FAILED_STATUS,
            f"Expected {SOURCE_BRANCH}@{CURRENT_SYNCED_BASE_COMMIT}, found {branch}@{head}",
        )
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    run = Path(run_root).resolve()
    output_root = config.paths.output_root.resolve()
    if run.parent != output_root or not run.name.startswith(f"{OUTPUT_PREFIX}_"):
        raise FlowError(FAILED_STATUS, f"Refusing non-axis-aligned run: {run}")
    manifest_path = run / "qc" / "axis_aligned_inlet_preflight_manifest.json"
    summary = read_json(manifest_path)
    if (
        summary.get("seeder_run_count") != 1
        or summary.get("seeder_return_code") != 0
        or summary.get("musubi_run_count") != 0
        or summary.get("harvester_run_count") != 0
    ):
        raise FlowError(
            FAILED_STATUS,
            "Existing run is not one successful Seeder with zero Musubi/Harvester calls",
        )

    inputs = load_flow_inputs(config.paths.source_surface_run)
    frozen_patch_dir = (
        config.paths.output_root / FROZEN_SEEDER_RUN / "geometry" / "geometry_solver_m"
    )
    frozen_paths = (
        inputs.tagged_surface_vtp,
        inputs.meter_surface_stl,
        *(frozen_patch_dir / f"{label}.stl" for label in BOUNDARY_LABELS),
    )
    frozen_before = summary["frozen_files_before"]
    frozen_after = _file_manifest(frozen_paths)
    frozen_unchanged = frozen_before == frozen_after

    geometry_qc = read_json(run / "qc" / "geometry_rigid_transform_qc.json")
    phase_qc = read_json(run / "qc" / "inlet_grid_phase_qc.json")
    mesh_preflight = read_json(run / "qc" / "seeder_mesh_qc.json")
    transform = read_json(run / "transform" / "anatomical_to_cfd_transform.json")
    environment = read_json(run / "input" / "apes_environment.json")
    mesh_summary = mesh_preflight["mesh_summary"]
    area_qc = mesh_preflight["d3q19_inlet"]["area_proxy"]
    normal_qc = mesh_preflight["d3q19_inlet"]["normal_ind"]
    actual_ram = int(mesh_preflight["actual_estimated_ram_bytes"])
    ram_limit = int(mesh_preflight["maximum_allowed_ram_bytes"])

    hard_gate_failures: list[dict[str, str]] = []
    if area_qc["status"] != "PASS":
        hard_gate_failures.append(
            {
                "status": AREA_PROXY_FAILED,
                "reason": f"Area proxy ratio is {area_qc['area_proxy_over_smooth_area']:.17g}",
            }
        )
    if normal_qc["status"] != "PASS":
        hard_gate_failures.append(
            {
                "status": NORMALIND_NOT_CLEAN,
                "reason": f"normalInd distribution is {normal_qc['distribution']}",
            }
        )
    if actual_ram > ram_limit:
        hard_gate_failures.append(
            {
                "status": FAILED_STATUS,
                "reason": "Actual Seeder mesh exceeds 60% RAM gate",
            }
        )
    if not mesh_preflight["uniform_level_matches_root"]:
        hard_gate_failures.append(
            {
                "status": FAILED_STATUS,
                "reason": "Seeder mesh is not uniform at the selected root level",
            }
        )
    if geometry_qc["status"] != "PASS":
        hard_gate_failures.append(
            {"status": FAILED_STATUS, "reason": "Rigid-transform geometry QC failed"}
        )
    if phase_qc["status"] != "PASS":
        hard_gate_failures.append(
            {"status": FAILED_STATUS, "reason": "Grid-phase alignment QC failed"}
        )
    if not frozen_unchanged:
        hard_gate_failures.append(
            {
                "status": FAILED_STATUS,
                "reason": "Frozen source surface or patch changed",
            }
        )

    velocity_new = EXPECTED_MASS_FLOW_KG_S / (
        EXPECTED_DENSITY_KG_M3 * area_qc["mfr_eq_area_proxy_m2"]
    )
    velocity_qc = {
        "status": "STATIC_PREDICTION_COMPLETE",
        "predicted_u_mfr_new_m_s": velocity_new,
        "smooth_expected_mean_velocity_m_s": EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
        "velocity_ratio": velocity_new / EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
        "prediction_only_no_cfd_run": True,
    }
    write_json(run / "qc" / "mfr_eq_static_velocity_prediction_qc.json", velocity_qc)
    write_json(
        run / "qc" / "frozen_source_read_only_sha_qc.json",
        {
            "status": "PASS" if frozen_unchanged else "FAIL",
            "source_frozen_surface_modified": not frozen_unchanged,
            "before": frozen_before,
            "after": frozen_after,
        },
    )
    mesh_manifest = _directory_manifest(run / "seeder" / "mesh")
    mesh_manifest.update(
        {
            "status": (
                "AXIS_ALIGNED_INLET_SEEDER_MESH_PREFLIGHT_PASS"
                if not hard_gate_failures
                else "AXIS_ALIGNED_INLET_SEEDER_MESH_PREFLIGHT_FAILED"
            ),
            "source_rotated_surface_sha256": sha256_file(
                run / "geometry" / "cfd_surface_axis_aligned_inlet_um.vtp"
            ),
            "actual_dx_m": EXPECTED_DX_M,
            "root_level": int(phase_qc["root_level"]),
            "fluid_element_count": int(mesh_summary["fluid_element_count"]),
            "boundary_element_count": int(mesh_summary["boundary_element_count"]),
            "qval_element_count": int(
                mesh_summary["property_element_counts"].get("has qVal", 0)
            ),
            "d3q19_inlet_element_count": int(area_qc["d3q19_inlet_element_count"]),
        }
    )
    write_json(run / "seeder" / "mesh_manifest.json", mesh_manifest)

    final_status = SUCCESS_STATUS if not hard_gate_failures else FAILED_STATUS
    final_next = SUCCESS_NEXT if not hard_gate_failures else FAILED_NEXT
    summary.update(
        {
            "status": final_status,
            "next": final_next,
            "hard_gate_failures": hard_gate_failures,
            "first_failure_status": (
                hard_gate_failures[0]["status"] if hard_gate_failures else None
            ),
            "first_failure_reason": (
                hard_gate_failures[0]["reason"] if hard_gate_failures else None
            ),
            "post_seeder_continuation": "PYTHON_READ_ONLY_TOPOLOGY_ANALYSIS_ONLY",
            "post_seeder_additional_seeder_calls": 0,
            "old_oblique_baseline_status": OLD_BASELINE_STATUS,
            "new_derived_geometry_path": str(run / "geometry"),
            "transform_path": str(
                run / "transform" / "anatomical_to_cfd_transform.json"
            ),
            "rotation_target": "INLET_INWARD_TO_NEGATIVE_Z",
            "rotation_angle_rad": transform["rotation_angle_rad"],
            "rotation_angle_deg": transform["rotation_angle_deg"],
            "rotation_matrix_3x3": transform["rotation_matrix_3x3"],
            "geometry_qc": geometry_qc,
            "grid_phase_qc": phase_qc,
            "dx_m": EXPECTED_DX_M,
            "root_level": int(phase_qc["root_level"]),
            "fluid_element_count": int(mesh_summary["fluid_element_count"]),
            "boundary_element_count": int(mesh_summary["boundary_element_count"]),
            "qval_element_count": int(
                mesh_summary["property_element_counts"].get("has qVal", 0)
            ),
            "new_inlet_d3q19_globbc_count": int(area_qc["d3q19_inlet_element_count"]),
            "smooth_inlet_area_m2": area_qc["smooth_inlet_area_m2"],
            "mfr_eq_area_proxy_m2": area_qc["mfr_eq_area_proxy_m2"],
            "area_proxy_over_smooth_area": area_qc["area_proxy_over_smooth_area"],
            "normal_ind_distribution": normal_qc["distribution"],
            "target_cardinal_normal_fraction": normal_qc[
                "target_cardinal_normal_fraction"
            ],
            "diagonal_normal_ind_count": normal_qc["diagonal_normal_ind_count"],
            "predicted_u_mfr_new_m_s": velocity_new,
            "smooth_expected_mean_velocity_m_s": EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
            "velocity_ratio": velocity_new / EXPECTED_SMOOTH_MEAN_VELOCITY_M_S,
            "ram_estimate_bytes": actual_ram,
            "ram_limit_bytes": ram_limit,
            "available_ram_bytes": int(environment["available_ram_bytes"]),
            "source_frozen_surface_modified": not frozen_unchanged,
            "frozen_files_after": frozen_after,
            "musubi_run_count": 0,
            "harvester_run_count": 0,
            "grid_convergence": "NOT_RUN",
            "completed_at": datetime.now().isoformat(),
        }
    )
    write_json(manifest_path, summary)
    return summary
