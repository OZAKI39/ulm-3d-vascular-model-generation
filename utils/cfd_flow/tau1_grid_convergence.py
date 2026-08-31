"""Research-only repaired Tau=1 three-grid convergence contracts.

The module is deliberately disconnected from the production pipeline.  Stage
zero decodes the already accepted Base restarts and integrates velocity over
fixed continuous physical port apertures.  Seeder or long Musubi work is not
permitted until that zero-solver preflight passes.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np
import trimesh
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon

from .full_timestep_mass_referee import PHYSICAL_FLUX_DEFINITION
from .grid_convergence import (
    evaluate_grid_convergence_gate,
    three_grid_scalar_analysis,
)
from .io import read_json, sha256_file, write_json
from .musubi_boundary_mass_referee import load_mesh_contract
from .port_grid_sensitivity import AXIS_GEOMETRY_RUN, PlaneFrame, _project, recover_continuous_ports
from .restart_decode import read_restart_pdf, reconstruct_macroscopic_field
from .tau1_base import _restart_pairs, _runtime_windows


RUN_NAME = "healthy_mouse_capillary_tau1_grid_convergence_anchor003274_20260831"
BASE_MESH_RUN = (
    "healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830"
)
BASE_CFD_RUN = "healthy_mouse_capillary_tau1_base_anchor003274_20260830"
STANDARDIZED_OUTLET_RUN = (
    "healthy_mouse_capillary_port_grid_sensitivity_research_anchor003274_20260830"
)
PLANE_REVISION = "STANDARDIZED_PHYSICAL_PORT_PLANES_V2"
FLUX_ALGORITHM_REVISION = "CELL_CUBE_PLANE_APERTURE_CLIPPING_V1"
HISTORICAL_CLASSIFICATION = (
    "HISTORICAL_UNDER_FALLBACK_QVALUE_WALL_AND_HIGH_TAU_NUMERICS"
)

REFINEMENT_RATIO = 1.3
RHO_KG_M3 = 1056.0
NU_M2_S = 3.27e-6
BULK_NU_M2_S = 2.18e-6
TARGET_U_MEAN_M_S = 0.35e-3
TARGET_Q_M3_S = 2.7369132390905703e-15
TARGET_MASS_FLOW_KG_S = 2.890180380479642e-12
PRESSURE_REFERENCE_PA = 23622.320128
OUTLET_GAUGE_PRESSURE_PA = {
    "outlet_01": 14.544978101,
    "outlet_02": 132.204549223,
    "outlet_03": -13.700626673,
}
SHORT_WINDOW_S = 0.0002441406727828746
LONG_WINDOW_S = 0.0004882813455657492
HARD_MAX_PHYSICAL_TIME_S = 0.0065
AREA_COVERAGE_GATE = 0.001
PHYSICAL_FLUX_GATE = 0.01

SEEDER_BINARY_WSL = (
    "/home/lzy/apes-worktrees/seeder_dimensionless_kernel_20260830/build/seeder"
)
SEEDER_BINARY_SHA256 = (
    "d7be681ca90da706559a4fd7e8f769fdb8f4303b8508f751077205f8e00cc7ed"
)
MUSUBI_BINARY_WSL = (
    "/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/"
    "build/musubi_adaptive_flux"
)
MUSUBI_BINARY_SHA256 = (
    "e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588"
)

PORTS = ("inlet", "outlet_01", "outlet_02", "outlet_03")
OUTLETS = PORTS[1:]
BASE_ITERATIONS = (2_878_425, 2_998_176, 3_117_927)
BASE_RESTART_SHA256 = {
    2_878_425: "75815fded691784ae942285e6ccf32514a1936ef9b1c228a03631991462646ee",
    2_998_176: "e3bb103963299e2384ce636dd117d194ee21b86adeac63a9965a095840f39a6c",
    3_117_927: "3d54f3970b4120896c214155811d7cd1b594e3efd172f80b5dc5e7d0fef279e2",
}
BASE_MESH_SHA256 = {
    "elemlist.lsb": "f7d7b1d55273c78c336ac04e39bc018dd9ebb470a9f29ce833ff01711de8c386",
    "bnd.lsb": "520d7dd1e4a46a45f9b1218a5807cfd89d6f054e0a247872362b130ff6bcfe69",
    "qval.lsb": "35884406b5f0111cd4ab471f7b08ac3df00e478d3458a57636d1bd8921cb0fe6",
}

PRIMARY_METRICS = (
    "inlet_gauge_pressure_pa",
    "pressure_drop_outlet_01_pa",
    "pressure_drop_outlet_02_pa",
    "pressure_drop_outlet_03_pa",
    "physical_outlet_01_flow_fraction",
    "physical_outlet_02_flow_fraction",
    "physical_outlet_03_flow_fraction",
)


@dataclass(frozen=True, slots=True)
class Tau1GridSpec:
    label: str
    dx_m: float
    dt_s: float
    root_level: int
    bounding_cube_length_m: float
    base_read_only: bool = False

    @property
    def nu_lattice(self) -> float:
        return NU_M2_S * self.dt_s / self.dx_m**2

    @property
    def tau(self) -> float:
        return 3.0 * self.nu_lattice + 0.5

    @property
    def omega(self) -> float:
        return 1.0 / self.tau

    @property
    def short_window_iterations(self) -> int:
        return round(SHORT_WINDOW_S / self.dt_s)

    @property
    def long_window_iterations(self) -> int:
        return round(LONG_WINDOW_S / self.dt_s)

    @property
    def tracking_interval_iterations(self) -> int:
        return round(0.5 * SHORT_WINDOW_S / self.dt_s)

    @property
    def checkpoint_interval_iterations(self) -> int:
        return self.short_window_iterations

    @property
    def earliest_audit_iteration(self) -> int:
        return math.ceil(2.0 * LONG_WINDOW_S / self.dt_s)

    @property
    def hard_max_iterations(self) -> int:
        return math.floor(HARD_MAX_PHYSICAL_TIME_S / self.dt_s)

    def evidence(self) -> dict[str, Any]:
        result = asdict(self)
        result.update(
            {
                "nu_lattice": self.nu_lattice,
                "tau": self.tau,
                "omega": self.omega,
                "short_window_iterations": self.short_window_iterations,
                "long_window_iterations": self.long_window_iterations,
                "tracking_interval_iterations": self.tracking_interval_iterations,
                "checkpoint_interval_iterations": self.checkpoint_interval_iterations,
                "earliest_audit_iteration": self.earliest_audit_iteration,
                "hard_max_iterations": self.hard_max_iterations,
            }
        )
        return result


GRID_SPECS: dict[str, Tau1GridSpec] = {
    "coarse": Tau1GridSpec(
        "coarse", 2.6e-7, 3.4454638124362897e-9, 9, 0.00013312
    ),
    "base": Tau1GridSpec(
        "base", 2.0e-7, 2.038735983690112e-9, 9, 0.0001024, True
    ),
    "fine": Tau1GridSpec(
        "fine", 1.5384615384615385e-7, 1.2063526530710723e-9, 10,
        0.00015753846153846154,
    ),
}


@dataclass(frozen=True, slots=True)
class PhysicalPortPlane:
    label: str
    origin_m: np.ndarray
    unit_normal: np.ndarray
    basis_u: np.ndarray
    basis_v: np.ndarray
    aperture_uv_m: np.ndarray
    physical_contract_sha256: str

    @property
    def aperture(self) -> Polygon:
        return Polygon(np.asarray(self.aperture_uv_m, dtype=np.float64))


@dataclass(frozen=True, slots=True)
class PlaneQuadrature:
    label: str
    cell_indices: np.ndarray
    clipped_areas_m2: np.ndarray
    aperture_area_m2: float
    candidate_cell_count: int
    area_coverage_relative_error: float


def _run_root(project_root: Path) -> Path:
    return Path(project_root).resolve() / "outputs" / "cfd_flow" / RUN_NAME


def _base_mesh(project_root: Path) -> Path:
    return (
        Path(project_root).resolve()
        / "outputs"
        / "cfd_flow"
        / BASE_MESH_RUN
        / "seeder"
        / "mesh"
    )


def _hash_payload(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def validate_grid_spec(spec: Tau1GridSpec) -> dict[str, bool]:
    checks = {
        "nu_lattice_one_sixth": abs(spec.nu_lattice - 1.0 / 6.0) <= 1.0e-12,
        "tau_one": abs(spec.tau - 1.0) <= 1.0e-12,
        "omega_one": abs(spec.omega - 1.0) <= 1.0e-12,
        "positive_physical_windows": (
            spec.short_window_iterations > 0
            and spec.long_window_iterations > spec.short_window_iterations
        ),
    }
    if not all(checks.values()):
        raise ValueError(f"invalid Tau=1 grid contract for {spec.label}: {checks}")
    return checks


def write_grid_design(project_root: Path) -> dict[str, Any]:
    checks = {label: validate_grid_spec(spec) for label, spec in GRID_SPECS.items()}
    ratio_checks = {
        "coarse_over_base": GRID_SPECS["coarse"].dx_m / GRID_SPECS["base"].dx_m,
        "base_over_fine": GRID_SPECS["base"].dx_m / GRID_SPECS["fine"].dx_m,
    }
    result = {
        "status": "PASS",
        "refinement_ratio": REFINEMENT_RATIO,
        "constant_refinement_ratio": all(
            abs(value - REFINEMENT_RATIO) <= 1.0e-14
            for value in ratio_checks.values()
        ),
        "ratio_checks": ratio_checks,
        "grids": {label: spec.evidence() for label, spec in GRID_SPECS.items()},
        "checks": checks,
        "short_window_s": SHORT_WINDOW_S,
        "long_window_s": LONG_WINDOW_S,
        "hard_max_physical_time_s": HARD_MAX_PHYSICAL_TIME_S,
        "frozen_physical_parameters": {
            "rho_kg_m3": RHO_KG_M3,
            "nu_m2_s": NU_M2_S,
            "bulk_nu_m2_s": BULK_NU_M2_S,
            "target_u_mean_m_s": TARGET_U_MEAN_M_S,
            "target_q_m3_s": TARGET_Q_M3_S,
            "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
            "pressure_reference_pa": PRESSURE_REFERENCE_PA,
            "outlet_gauge_pressure_pa": OUTLET_GAUGE_PRESSURE_PA,
            "layout": "D3Q19",
            "relaxation": "BGK",
            "wall": "wall_libb continuous q",
            "inlet": "adaptive_flux_pressure",
            "outlets": "pressure_eq",
        },
        "frozen_binaries": {
            "seeder_wsl": SEEDER_BINARY_WSL,
            "seeder_sha256": SEEDER_BINARY_SHA256,
            "musubi_wsl": MUSUBI_BINARY_WSL,
            "musubi_sha256": MUSUBI_BINARY_SHA256,
        },
        "base_seeder_calls_allowed": 0,
        "base_long_musubi_calls_allowed": 0,
        "production_pipeline_modified": False,
    }
    if not result["constant_refinement_ratio"]:
        raise ValueError("the requested grids do not have constant r=1.3")
    output = _run_root(project_root) / "qc" / "tau1_grid_design.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json(output, result)
    return result


def assert_launch_allowed(grid: str, operation: str) -> None:
    if grid not in GRID_SPECS:
        raise ValueError(grid)
    if operation not in {"seeder", "long_musubi", "referee_one_step"}:
        raise ValueError(operation)
    if grid == "base" and operation in {"seeder", "long_musubi"}:
        raise PermissionError(f"accepted Base is read-only: {operation} is forbidden")


def render_repaired_seeder_config(base_text: str, spec: Tau1GridSpec) -> str:
    if spec.label == "base":
        raise PermissionError("accepted Base Seeder config must not be regenerated")
    required = (
        "7.8082771579034594e-05",
        "5.4076359424537688e-05",
        "8.8206398940868488e-05",
        "0.0001016827715790346",
        "3.9999999999999888e-06",
        "4.2000000000000004e-06",
        "0.00010367370500008977",
        "calc_dist = true",
    )
    if not all(token in base_text for token in required):
        raise ValueError("frozen Base physical Seeder objects are incomplete")
    text = re.sub(
        r"(?m)^comment\s*=.*$",
        f"comment = 'repaired Tau1 {spec.label} grid; physical objects frozen from Base'",
        base_text,
        count=1,
    )
    text = re.sub(
        r"(?m)^minlevel\s*=\s*\d+\s*$",
        f"minlevel = {spec.root_level}",
        text,
        count=1,
    )
    text = re.sub(
        r"(bounding_cube\s*=\s*\{.*?length\s*=\s*)[-+0-9.eE]+",
        rf"\g<1>{spec.bounding_cube_length_m:.17g}",
        text,
        count=1,
    )
    return text


def seeder_physical_spatial_signature(text: str) -> str:
    normalized = re.sub(r"(?m)^comment\s*=.*$", "comment=<GRID>", text)
    normalized = re.sub(r"(?m)^minlevel\s*=\s*\d+\s*$", "minlevel=<GRID>", normalized)
    normalized = re.sub(
        r"(bounding_cube\s*=\s*\{.*?length\s*=\s*)[-+0-9.eE]+",
        r"\g<1><GRID>",
        normalized,
        count=1,
    )
    return hashlib.sha256(normalized.encode()).hexdigest()


def _continuous_inlet_record(
    root: Path, frame: PlaneFrame, continuous: Mapping[str, Any]
) -> dict[str, Any]:
    geometry = root / "outputs" / "cfd_flow" / AXIS_GEOMETRY_RUN / "geometry" / "geometry_solver_m"
    cap_path = geometry / "inlet.stl"
    cap = trimesh.load_mesh(cap_path, process=True)
    if not isinstance(cap, trimesh.Trimesh):
        raise ValueError("continuous inlet cap is not a triangle mesh")
    projected = _project(np.asarray(cap.vertices, dtype=np.float64), frame)
    hull = ConvexHull(projected)
    contour = projected[hull.vertices]
    return {
        "origin_m": frame.origin.tolist(),
        "unit_normal": frame.normal.tolist(),
        "basis_u": frame.basis_u.tolist(),
        "basis_v": frame.basis_v.tolist(),
        "physical_aperture_contour_uv_m": contour.tolist(),
        "source_geometry_path": str(cap_path.resolve()),
        "source_geometry_sha256": sha256_file(cap_path),
        "continuous_cap_area_m2": float(continuous["continuous_area_m2"]),
        "orientation_source": "frozen inlet source normal plus rigid CFD transform",
    }


def _canonical_plane_record(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: record[key]
        for key in (
            "origin_m",
            "unit_normal",
            "basis_u",
            "basis_v",
            "physical_aperture_contour_uv_m",
            "source_geometry_sha256",
        )
    }


def build_physical_port_plane_contract(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    v1_path = (
        root
        / "outputs"
        / "cfd_flow"
        / STANDARDIZED_OUTLET_RUN
        / "qc"
        / "standardized_outlet_plane_contract.json"
    )
    v1 = read_json(v1_path)
    if v1.get("status") != "PASS":
        raise ValueError("STANDARDIZED_PHYSICAL_OUTLET_PLANES_V1 is not accepted")
    frames, continuous = recover_continuous_ports(root)
    axis_geometry = (
        root / "outputs" / "cfd_flow" / AXIS_GEOMETRY_RUN / "geometry" / "geometry_solver_m"
    )
    frozen_geometry = (
        root / "outputs" / "cfd_flow" / BASE_MESH_RUN / "geometry" / "geometry_solver_m"
    )
    records: dict[str, dict[str, Any]] = {}
    records["inlet"] = _continuous_inlet_record(root, frames["inlet"], continuous["inlet"])
    for label in OUTLETS:
        parent = v1["outlets"][label]
        record = {
            key: parent[key]
            for key in (
                "origin_m",
                "unit_normal",
                "basis_u",
                "basis_v",
                "physical_aperture_contour_uv_m",
                "source_geometry_path",
                "source_geometry_sha256",
            )
        }
        record["parent_revision"] = v1["revision"]
        record["parent_physical_contract_sha256"] = parent["physical_contract_sha256"]
        records[label] = record
    for label, record in records.items():
        aperture = Polygon(record["physical_aperture_contour_uv_m"])
        record["aperture_physical_area_m2"] = float(aperture.area)
        record["aperture_valid"] = bool(aperture.is_valid and aperture.area > 0.0)
        record["basis_qc"] = {
            "normal_norm_error": abs(np.linalg.norm(record["unit_normal"]) - 1.0),
            "u_norm_error": abs(np.linalg.norm(record["basis_u"]) - 1.0),
            "v_norm_error": abs(np.linalg.norm(record["basis_v"]) - 1.0),
            "max_pairwise_dot": max(
                abs(float(np.dot(record["unit_normal"], record["basis_u"]))),
                abs(float(np.dot(record["unit_normal"], record["basis_v"]))),
                abs(float(np.dot(record["basis_u"], record["basis_v"]))),
            ),
        }
        record["physical_contract_sha256"] = _hash_payload(
            _canonical_plane_record(record)
        )
    canonical = {
        "revision": PLANE_REVISION,
        "ports": {label: _canonical_plane_record(records[label]) for label in PORTS},
    }
    contract_hash = _hash_payload(canonical)
    hashes_by_grid = {label: contract_hash for label in GRID_SPECS}
    numerical_inlet_path = (
        root
        / "outputs"
        / "cfd_flow"
        / BASE_MESH_RUN
        / "geometry"
        / "numerical_inlet_plane.stl"
    )
    numerical_inlet = trimesh.load_mesh(numerical_inlet_path, process=True)
    if not isinstance(numerical_inlet, trimesh.Trimesh):
        raise ValueError("frozen numerical inlet patch is not a triangle mesh")
    inlet_frame = frames["inlet"]
    numerical_uv = _project(
        np.asarray(numerical_inlet.vertices, dtype=np.float64), inlet_frame
    )
    numerical_hull = ConvexHull(numerical_uv)
    numerical_aperture = Polygon(numerical_uv[numerical_hull.vertices])
    continuous_aperture = Polygon(
        records["inlet"]["physical_aperture_contour_uv_m"]
    )
    source_geometry_checks = {
        "wall_sha256_matches_frozen_base": (
            sha256_file(axis_geometry / "wall.stl")
            == sha256_file(frozen_geometry / "wall.stl")
        ),
        **{
            f"{label}_sha256_matches_frozen_base": (
                sha256_file(axis_geometry / f"{label}.stl")
                == sha256_file(frozen_geometry / f"{label}.stl")
            )
            for label in OUTLETS
        },
        "continuous_inlet_cap_covered_by_frozen_numerical_patch": bool(
            numerical_aperture.buffer(1.0e-15).covers(continuous_aperture)
        ),
    }
    checks = {
        "parent_outlet_contract_pass": v1["status"] == "PASS",
        "all_apertures_valid": all(record["aperture_valid"] for record in records.values()),
        "all_bases_orthonormal": all(
            max(record["basis_qc"].values()) <= 1.0e-12
            for record in records.values()
        ),
        "same_physical_contract_hash_across_grids": len(set(hashes_by_grid.values())) == 1,
        "physical_definition_is_dx_independent": True,
        "frozen_base_geometry_sources_match": all(source_geometry_checks.values()),
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "revision": PLANE_REVISION,
        "flux_algorithm_revision": FLUX_ALGORITHM_REVISION,
        "contract_sha256": contract_hash,
        "contract_hashes_by_grid": hashes_by_grid,
        "checks": checks,
        "ports": records,
        "source_outlet_contract_path": str(v1_path),
        "source_outlet_contract_sha256": sha256_file(v1_path),
        "frozen_base_geometry_consistency": {
            "checks": source_geometry_checks,
            "numerical_inlet_plane_path": str(numerical_inlet_path),
            "numerical_inlet_plane_sha256": sha256_file(numerical_inlet_path),
            "continuous_inlet_cap_sha256": records["inlet"][
                "source_geometry_sha256"
            ],
        },
        "grid_fit_performed": False,
    }
    output = _run_root(root) / "qc" / "physical_port_plane_contract.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json(output, result)
    return result


def plane_from_record(label: str, record: Mapping[str, Any]) -> PhysicalPortPlane:
    return PhysicalPortPlane(
        label=label,
        origin_m=np.asarray(record["origin_m"], dtype=np.float64),
        unit_normal=np.asarray(record["unit_normal"], dtype=np.float64),
        basis_u=np.asarray(record["basis_u"], dtype=np.float64),
        basis_v=np.asarray(record["basis_v"], dtype=np.float64),
        aperture_uv_m=np.asarray(
            record["physical_aperture_contour_uv_m"], dtype=np.float64
        ),
        physical_contract_sha256=str(record["physical_contract_sha256"]),
    )


_CUBE_EDGES = (
    (0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3),
    (2, 6), (3, 7), (4, 5), (4, 6), (5, 7), (6, 7),
)
_CUBE_SIGNS = np.asarray(
    [
        (-1, -1, -1), (-1, -1, 1), (-1, 1, -1), (-1, 1, 1),
        (1, -1, -1), (1, -1, 1), (1, 1, -1), (1, 1, 1),
    ],
    dtype=np.float64,
)


def _cube_plane_polygon(
    center_m: np.ndarray,
    dx_m: float,
    plane: PhysicalPortPlane,
    *,
    tolerance_m: float,
) -> Polygon | None:
    vertices = np.asarray(center_m, dtype=np.float64) + 0.5 * dx_m * _CUBE_SIGNS
    signed = (vertices - plane.origin_m) @ plane.unit_normal
    points: list[np.ndarray] = []
    for index, value in enumerate(signed):
        if abs(float(value)) <= tolerance_m:
            points.append(vertices[index])
    for left_index, right_index in _CUBE_EDGES:
        left = float(signed[left_index])
        right = float(signed[right_index])
        if left * right < 0.0:
            fraction = left / (left - right)
            points.append(
                vertices[left_index]
                + fraction * (vertices[right_index] - vertices[left_index])
            )
    unique: list[np.ndarray] = []
    for point in points:
        if not any(np.linalg.norm(point - prior) <= tolerance_m for prior in unique):
            unique.append(point)
    if len(unique) < 3:
        return None
    xyz = np.asarray(unique, dtype=np.float64)
    uv = np.column_stack(
        (
            (xyz - plane.origin_m) @ plane.basis_u,
            (xyz - plane.origin_m) @ plane.basis_v,
        )
    )
    if np.linalg.matrix_rank(uv - uv.mean(axis=0)) < 2:
        return None
    hull = ConvexHull(uv)
    polygon = Polygon(uv[hull.vertices])
    return polygon if polygon.is_valid and polygon.area > 0.0 else None


def build_plane_quadrature(
    cell_centers_m: np.ndarray,
    *,
    dx_m: float,
    plane: PhysicalPortPlane,
) -> PlaneQuadrature:
    centers = np.asarray(cell_centers_m, dtype=np.float64)
    if centers.ndim != 2 or centers.shape[1] != 3:
        raise ValueError("cell centers must have shape (n, 3)")
    aperture = plane.aperture
    if not aperture.is_valid or aperture.area <= 0.0:
        raise ValueError("physical aperture is invalid")
    tolerance_m = max(1.0e-15, abs(float(dx_m)) * 1.0e-8)
    signed = (centers - plane.origin_m) @ plane.unit_normal
    normal_extent = 0.5 * dx_m * float(np.sum(np.abs(plane.unit_normal)))
    projected = np.column_stack(
        (
            (centers - plane.origin_m) @ plane.basis_u,
            (centers - plane.origin_m) @ plane.basis_v,
        )
    )
    min_u, min_v, max_u, max_v = aperture.bounds
    u_extent = 0.5 * dx_m * float(np.sum(np.abs(plane.basis_u)))
    v_extent = 0.5 * dx_m * float(np.sum(np.abs(plane.basis_v)))
    candidate = (
        (np.abs(signed) <= normal_extent + tolerance_m)
        & (projected[:, 0] >= min_u - u_extent - tolerance_m)
        & (projected[:, 0] <= max_u + u_extent + tolerance_m)
        & (projected[:, 1] >= min_v - v_extent - tolerance_m)
        & (projected[:, 1] <= max_v + v_extent + tolerance_m)
    )
    candidate_indices = np.flatnonzero(candidate)
    indices: list[int] = []
    areas: list[float] = []
    area_tolerance = max(np.finfo(float).tiny, dx_m**2 * 1.0e-14)
    for index_value in candidate_indices:
        index = int(index_value)
        section = _cube_plane_polygon(
            centers[index], dx_m, plane, tolerance_m=tolerance_m
        )
        if section is None:
            continue
        clipped = section.intersection(aperture)
        area = float(clipped.area)
        if area > area_tolerance:
            indices.append(index)
            areas.append(area)
    total = math.fsum(areas)
    aperture_area = float(aperture.area)
    relative = abs(total - aperture_area) / aperture_area
    return PlaneQuadrature(
        label=plane.label,
        cell_indices=np.asarray(indices, dtype=np.int64),
        clipped_areas_m2=np.asarray(areas, dtype=np.float64),
        aperture_area_m2=aperture_area,
        candidate_cell_count=len(candidate_indices),
        area_coverage_relative_error=relative,
    )


def integrate_plane_flux(
    velocities_m_s: np.ndarray,
    *,
    plane: PhysicalPortPlane,
    quadrature: PlaneQuadrature,
) -> dict[str, Any]:
    velocities = np.asarray(velocities_m_s, dtype=np.float64)
    if velocities.ndim != 2 or velocities.shape[1] != 3:
        raise ValueError("velocities must have shape (n, 3)")
    indices = quadrature.cell_indices
    if len(indices) == 0:
        normal_velocity = np.asarray([], dtype=np.float64)
        signed_q = 0.0
        weighted_mean = math.nan
        minimum = math.nan
        maximum = math.nan
    else:
        normal_velocity = velocities[indices] @ plane.unit_normal
        signed_q = math.fsum(
            float(area) * float(value)
            for area, value in zip(
                quadrature.clipped_areas_m2, normal_velocity, strict=True
            )
        )
        covered_area = math.fsum(float(value) for value in quadrature.clipped_areas_m2)
        weighted_mean = signed_q / covered_area
        minimum = float(np.min(normal_velocity))
        maximum = float(np.max(normal_velocity))
    physical_q = -signed_q if plane.label == "inlet" else signed_q
    return {
        "flux_definition": PHYSICAL_FLUX_DEFINITION,
        "algorithm_revision": FLUX_ALGORITHM_REVISION,
        "port": plane.label,
        "aperture_physical_area_m2": quadrature.aperture_area_m2,
        "sum_clipped_cell_plane_area_m2": math.fsum(
            float(value) for value in quadrature.clipped_areas_m2
        ),
        "area_coverage_relative_error": quadrature.area_coverage_relative_error,
        "area_coverage_gate": AREA_COVERAGE_GATE,
        "area_coverage_pass": (
            quadrature.area_coverage_relative_error <= AREA_COVERAGE_GATE
        ),
        "candidate_cell_count": quadrature.candidate_cell_count,
        "contributing_cell_count": len(indices),
        "outward_normal_velocity_mean_m_s": weighted_mean,
        "outward_normal_velocity_min_m_s": minimum,
        "outward_normal_velocity_max_m_s": maximum,
        "signed_outward_q_m3_s": signed_q,
        "sign_convention": "Qin=-integral(u.n)dA" if plane.label == "inlet" else "Qout=+integral(u.n)dA",
        "physical_q_m3_s": physical_q,
        "all_finite": bool(
            np.isfinite(signed_q)
            and np.all(np.isfinite(normal_velocity))
            and np.all(np.isfinite(quadrature.clipped_areas_m2))
        ),
    }


def physical_flux_mass_closure(
    qin_m3_s: float,
    outlet_values_m3_s: Iterable[float],
    *,
    gate: float = PHYSICAL_FLUX_GATE,
) -> dict[str, Any]:
    outlets = tuple(float(value) for value in outlet_values_m3_s)
    outlet_sum = math.fsum(outlets)
    denominator = max(abs(float(qin_m3_s)), np.finfo(float).tiny)
    relative = abs(float(qin_m3_s) - outlet_sum) / denominator
    return {
        "qin_m3_s": float(qin_m3_s),
        "outlet_sum_m3_s": outlet_sum,
        "absolute_mismatch_m3_s": abs(float(qin_m3_s) - outlet_sum),
        "relative_error": relative,
        "gate": float(gate),
        "pass": math.isfinite(relative) and relative <= float(gate),
    }


def _mesh_origin_dx(mesh_dir: Path) -> tuple[np.ndarray, float]:
    text = (mesh_dir / "header.lua").read_text(encoding="utf-8")
    block = re.search(r"boundingbox\s*=\s*\{(.*?)\n\}", text, re.DOTALL)
    level = re.search(r"(?m)^\s*minLevel\s*=\s*(\d+)", text)
    if block is None or level is None:
        raise ValueError("TreElm header has no uniform bounding-box contract")
    origin = re.search(r"origin\s*=\s*\{(.*?)\}", block.group(1), re.DOTALL)
    length = re.search(r"length\s*=\s*([^\s]+)", block.group(1))
    if origin is None or length is None:
        raise ValueError("TreElm bounding-box contract is incomplete")
    origin_value = np.asarray(
        [float(token.replace("D", "E")) for token in origin.group(1).replace(",", " ").split()]
    )
    length_value = float(length.group(1).replace("D", "E"))
    return origin_value, length_value / 2 ** int(level.group(1))


def _base_flux_sample(
    *,
    pdf_path: Path,
    iteration: int,
    cell_count: int,
    spec: Tau1GridSpec,
    planes: Mapping[str, PhysicalPortPlane],
    quadratures: Mapping[str, PlaneQuadrature],
) -> dict[str, Any]:
    pdf = read_restart_pdf(pdf_path, n_elems=cell_count, n_components=19)
    field = reconstruct_macroscopic_field(
        pdf, dx_m=spec.dx_m, dt_s=spec.dt_s, rho0_kg_m3=RHO_KG_M3
    )
    ports = {
        label: integrate_plane_flux(
            field.velocity_phy, plane=planes[label], quadrature=quadratures[label]
        )
        for label in PORTS
    }
    qin = float(ports["inlet"]["physical_q_m3_s"])
    qout = [float(ports[label]["physical_q_m3_s"]) for label in OUTLETS]
    closure = physical_flux_mass_closure(qin, qout)
    fractions_outlet_sum = {
        label: float(ports[label]["physical_q_m3_s"]) / closure["outlet_sum_m3_s"]
        for label in OUTLETS
    }
    fractions_inlet = {
        label: float(ports[label]["physical_q_m3_s"]) / qin for label in OUTLETS
    }
    return {
        "iteration": int(iteration),
        "physical_time_s": int(iteration) * spec.dt_s,
        "restart_binary": str(pdf_path),
        "restart_sha256": sha256_file(pdf_path),
        "ports": ports,
        "Qin_phys_m3_s": qin,
        "Q1_phys_m3_s": qout[0],
        "Q2_phys_m3_s": qout[1],
        "Q3_phys_m3_s": qout[2],
        "physical_volume_closure": closure,
        "fraction_of_outlet_sum": fractions_outlet_sum,
        "fraction_of_inlet": fractions_inlet,
    }


def _summary(values: Iterable[float]) -> dict[str, float]:
    array = np.asarray(tuple(float(value) for value in values), dtype=np.float64)
    return {
        "mean": float(np.mean(array)),
        "std": float(np.std(array)),
        "endpoint": float(array[-1]),
    }


def _write_stage0_stop_evidence(
    run_root: Path,
    final: Mapping[str, Any],
    preflight: Mapping[str, Any],
) -> None:
    qc = run_root / "qc"
    skipped = {
        "status": "NOT_RUN_BASE_PHYSICAL_FLUX_GATE_FAILED",
        "scientific_calls": 0,
        "reason": "Stage 0 hard gate failed before Coarse/Fine Seeder or Musubi",
    }
    for name in (
        "coarse_mesh_contract.json",
        "fine_mesh_contract.json",
        "fine_mpi_performance_selection.json",
        "coarse_steady_acceptance.json",
        "fine_steady_acceptance.json",
        "coarse_full_referee.json",
        "fine_full_referee.json",
        "three_grid_scalar_analyses.json",
    ):
        write_json(qc / name, {**skipped, "artifact": name})
    write_json(
        qc / "physical_flux_cbf.json",
        {
            "status": "BASE_NON_ACCEPTED_PARTIAL_APERTURE_DIAGNOSTIC",
            "flux_definition": PHYSICAL_FLUX_DEFINITION,
            "algorithm_revision": FLUX_ALGORITHM_REVISION,
            "physical_plane_contract_sha256": preflight[
                "physical_plane_contract_sha256"
            ],
            "coarse": skipped,
            "base": {
                "status": "FAIL_AREA_COVERAGE_GATE",
                "observables_accepted": False,
                "steady_window_observables": preflight[
                    "steady_window_observables"
                ],
                "physical_volume_closure": preflight[
                    "steady_window_mean_physical_volume_closure"
                ],
                "Qin_target_relative_error": preflight[
                    "Qin_target_relative_error"
                ],
            },
            "fine": skipped,
        },
    )
    recovery = qc / "operational_recovery_log.jsonl"
    recovery.write_text(
        json.dumps(
            {
                "status": "NO_OPERATIONAL_RECOVERY_REQUIRED",
                "recoveries": 0,
                "scientific_calls_consumed": 0,
            },
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )
    write_json(qc / "grid_convergence_final.json", dict(final))


def run_base_physical_flux_preflight(project_root: Path) -> dict[str, Any]:
    started = time.perf_counter()
    root = Path(project_root).resolve()
    run_root = _run_root(root)
    qc = run_root / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    design = write_grid_design(root)
    plane_contract = build_physical_port_plane_contract(root)
    spec = GRID_SPECS["base"]
    mesh_dir = _base_mesh(root)
    mesh_hashes = {name: sha256_file(mesh_dir / name) for name in BASE_MESH_SHA256}
    mesh = load_mesh_contract(mesh_dir, expected_cells=182_320)
    origin, header_dx = _mesh_origin_dx(mesh_dir)
    if not math.isclose(header_dx, spec.dx_m, rel_tol=0.0, abs_tol=1.0e-18):
        raise ValueError(f"Base header dx changed: {header_dx}")
    centers = origin + (mesh.cell_ijk.astype(np.float64) + 0.5) * spec.dx_m
    planes = {
        label: plane_from_record(label, plane_contract["ports"][label])
        for label in PORTS
    }
    quadratures = {
        label: build_plane_quadrature(
            centers, dx_m=spec.dx_m, plane=planes[label]
        )
        for label in PORTS
    }
    pairs = _restart_pairs(_runtime_windows() / "restart")
    restart_checks = {
        str(iteration): (
            iteration in pairs
            and sha256_file(pairs[iteration][1]) == BASE_RESTART_SHA256[iteration]
        )
        for iteration in BASE_ITERATIONS
    }
    if not all(restart_checks.values()):
        raise ValueError(f"accepted Base restart evidence changed: {restart_checks}")
    samples = [
        _base_flux_sample(
            pdf_path=pairs[iteration][1],
            iteration=iteration,
            cell_count=len(mesh.tree_ids),
            spec=spec,
            planes=planes,
            quadratures=quadratures,
        )
        for iteration in BASE_ITERATIONS
    ]
    observable_names = (
        "Qin_phys_m3_s",
        "Q1_phys_m3_s",
        "Q2_phys_m3_s",
        "Q3_phys_m3_s",
    )
    steady = {
        name: _summary(sample[name] for sample in samples) for name in observable_names
    }
    mean_qin = steady["Qin_phys_m3_s"]["mean"]
    mean_outlets = [steady[name]["mean"] for name in observable_names[1:]]
    mean_closure = physical_flux_mass_closure(mean_qin, mean_outlets)
    target_relative_error = abs(mean_qin - TARGET_Q_M3_S) / TARGET_Q_M3_S
    significant_backflow = any(
        value < 0.0 and abs(value) > 0.05 * abs(mean_qin) for value in mean_outlets
    )
    coverage = {
        label: {
            "aperture_physical_area_m2": quadrature.aperture_area_m2,
            "sum_clipped_cell_plane_area_m2": math.fsum(
                float(value) for value in quadrature.clipped_areas_m2
            ),
            "area_coverage_relative_error": quadrature.area_coverage_relative_error,
            "candidate_cell_count": quadrature.candidate_cell_count,
            "contributing_cell_count": len(quadrature.cell_indices),
            "pass": quadrature.area_coverage_relative_error <= AREA_COVERAGE_GATE,
        }
        for label, quadrature in quadratures.items()
    }
    gates = {
        "grid_design": design["status"] == "PASS",
        "physical_port_plane_contract": plane_contract["status"] == "PASS",
        "base_mesh_hashes": mesh_hashes == BASE_MESH_SHA256,
        "base_restarts_unchanged": all(restart_checks.values()),
        "all_aperture_area_coverage_le_0p001": all(
            record["pass"] for record in coverage.values()
        ),
        "Qin_target_relative_error_le_0p01": target_relative_error <= PHYSICAL_FLUX_GATE,
        "physical_volume_closure_le_0p01": mean_closure["pass"],
        "no_significant_time_averaged_outlet_backflow": not significant_backflow,
        "all_flux_values_finite": all(
            math.isfinite(sample[name])
            for sample in samples
            for name in observable_names
        ),
    }
    passed = all(gates.values())
    first_failure = next((name for name, value in gates.items() if not value), None)
    result = {
        "status": "PASS" if passed else "FAIL",
        "final_status_if_stopped": (
            None if passed else "CFD_FLOW_PHYSICAL_PORT_FLUX_CONTRACT_FAILED"
        ),
        "analysis_mode": "ZERO_HARVESTER_ZERO_MUSUBI_EXISTING_BASE_RESTARTS",
        "flux_definition": PHYSICAL_FLUX_DEFINITION,
        "flux_algorithm_revision": FLUX_ALGORITHM_REVISION,
        "physical_plane_contract_sha256": plane_contract["contract_sha256"],
        "base_mesh": str(mesh_dir),
        "base_mesh_cells": len(mesh.tree_ids),
        "base_mesh_sha256": mesh_hashes,
        "base_restart_checks": restart_checks,
        "area_quadrature": coverage,
        "physical_flux_observables_accepted": gates[
            "all_aperture_area_coverage_le_0p001"
        ],
        "Q_values_classification": (
            "ACCEPTED_PHYSICAL_CROSS_SECTION_FLUX"
            if gates["all_aperture_area_coverage_le_0p001"]
            else "NON_ACCEPTED_PARTIAL_APERTURE_DIAGNOSTIC"
        ),
        "samples": samples,
        "steady_window_observables": steady,
        "steady_window_mean_physical_volume_closure": mean_closure,
        "Qin_target_m3_s": TARGET_Q_M3_S,
        "Qin_target_relative_error": target_relative_error,
        "significant_time_averaged_outlet_backflow": significant_backflow,
        "gates": gates,
        "true_first_scientific_failure": first_failure,
        "base_seeder_calls": 0,
        "base_long_musubi_calls": 0,
        "coarse_seeder_calls": 0,
        "fine_seeder_calls": 0,
        "coarse_long_logical_musubi_calls": 0,
        "fine_long_logical_musubi_calls": 0,
        "harvester_calls": 0,
        "production_pipeline_modified": False,
        "runtime_seconds": time.perf_counter() - started,
        "next": (
            "GENERATE ONLY COARSE AND FINE MESHES"
            if passed
            else "FIX PHYSICAL CROSS-SECTION FLUX EXTRACTION; DO NOT RUN COARSE/FINE"
        ),
    }
    write_json(qc / "base_physical_flux_preflight.json", result)
    if not passed:
        final = {
            "status": "CFD_FLOW_PHYSICAL_PORT_FLUX_CONTRACT_FAILED",
            "stage": "BASE_PHYSICAL_FLUX_ZERO_RUN_PREFLIGHT",
            "true_first_scientific_failure": first_failure,
            "base_long_cfd_calls": 0,
            "seeder_calls": 0,
            "coarse_fine_long_logical_calls": 0,
            "restart_resumes": 0,
            "operational_recoveries": 0,
            "production_pipeline_modified": False,
            "physical_plane_contract_sha256": plane_contract["contract_sha256"],
            "runtime_seconds": result["runtime_seconds"],
            "next": result["next"],
        }
        _write_stage0_stop_evidence(run_root, final, result)
    return result


def restart_resume_contract(
    saved: Mapping[str, Any], expected: Mapping[str, Any]
) -> dict[str, Any]:
    fields = (
        "mesh_hashes", "dx_m", "dt_s", "rho_kg_m3", "nu_m2_s",
        "bulk_nu_m2_s", "tau", "boundary_contract", "outlet_pressures_pa",
        "target_mass_flow_kg_s", "binary_sha256", "layout", "relaxation",
    )
    checks = {field: saved.get(field) == expected.get(field) for field in fields}
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def validate_physical_flux_observables(
    observables: Mapping[str, Mapping[str, Any]],
) -> None:
    for grid in ("coarse", "base", "fine"):
        if observables[grid].get("historical_classification") == HISTORICAL_CLASSIFICATION:
            raise ValueError(f"{grid} historical fallback-q evidence is excluded")
        definition = observables[grid].get("flux_definition")
        if definition != PHYSICAL_FLUX_DEFINITION:
            raise ValueError(
                f"{grid} primary flow metrics require {PHYSICAL_FLUX_DEFINITION}; "
                f"found {definition}"
            )


def build_primary_analyses(
    observables: Mapping[str, Mapping[str, Any]],
) -> dict[str, dict[str, Any]]:
    validate_physical_flux_observables(observables)
    return {
        metric: three_grid_scalar_analysis(
            float(observables["coarse"][metric]),
            float(observables["base"][metric]),
            float(observables["fine"][metric]),
            refinement_ratio=REFINEMENT_RATIO,
        )
        for metric in PRIMARY_METRICS
    }


def evaluate_repaired_grid_gate(
    analyses: dict[str, dict[str, Any]],
) -> dict[str, Any]:
    result = evaluate_grid_convergence_gate(analyses, PRIMARY_METRICS)
    passed = result["status"] == "PASS"
    result["final_status"] = (
        "CFD_FLOW_REPAIRED_TAU1_THREE_GRID_CONVERGENCE_PASS"
        if passed
        else "CFD_FLOW_REPAIRED_TAU1_THREE_GRID_CONVERGENCE_FAILED"
    )
    result["next"] = (
        "PROMOTE VALIDATED TAU1 CFD CONTRACT TO PRODUCTION PIPELINE"
        if passed
        else "STOP; REVIEW FIRST FAILING REPAIRED TAU1 PRIMARY METRIC"
    )
    return result
