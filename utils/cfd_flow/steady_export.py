"""One-shot field export from the frozen project-steady Musubi restart."""

from __future__ import annotations

import hashlib
import math
import re
import shutil
import subprocess
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

from .apes import (
    inspect_apes_environment,
    load_boundary_conditions,
    run_wsl_tool,
    windows_to_wsl,
)
from .config import load_cfd_flow_config
from .diagnostics import validate_continuation_restart
from .geometry import SurfacePartition, load_frozen_surface_partition
from .io import FlowError, load_flow_inputs, read_json, sha256_file, write_json
from .qc import numerical_port_fluxes, reynolds_diagnostics


SOURCE_COMMIT = "ebed073222a84d2971ff50c2f490c805e7c7d59f"
SOURCE_RUN_NAME = "musubi_project_steady_confirmation_anchor003274_20260828_225334"
FROZEN_SEEDER_NAME = "musubi_recovery_anchor003274_20260828_162530"
SOLVER_WORKDIR_NAME = "musubi_solver_recovery_anchor003274_20260828_164912"
EXPORT_PREFIX = "musubi_steady_field_export_anchor003274"
EXPORT_REVISION = "FROZEN_PROJECT_STEADY_ASCII_SPATIAL_EXPORT_V1"
SUCCESS_STATUS = "CFD_FLOW_STEADY_FIELD_EXPORT_PASS_PENDING_GRID_CONVERGENCE"
SUCCESS_NEXT = "RUN MUSUBI GRID-SPACING CONVERGENCE STUDY"
HARVEST_FAILED = "CFD_FLOW_STEADY_ASCII_HARVEST_FAILED"
SCHEMA_INVALID = "CFD_FLOW_STEADY_ASCII_SCHEMA_INVALID"
PATH_INVALID = "CFD_FLOW_STEADY_EXPORT_ASCII_PATH_INVALID"
FLUX_FAILED = "CFD_FLOW_STEADY_FIELD_FLUX_QC_FAILED"
SOURCE_STATE = "FROZEN_PROJECT_STEADY_RESTART"
FROZEN_ITERATION = 198_064
EXPECTED_CELL_COUNT = 221_109
PRESSURE_REFERENCE_PA = 23_622.32012800001
PROJECT_PRESSURE_THRESHOLD_PA = 0.145905175896487
PROJECT_VELOCITY_THRESHOLD_M_S = 9.838558007536173e-8
PROJECT_PRESSURE_RATIO = 0.6126756145994807
PROJECT_VELOCITY_RATIO = 0.9987856073836835
ASCII_FOLDER_WSL = "/tmp/u3d/x/"
ASCII_LABEL = "f"
ASCII_MAX_PREDICTED_LENGTH = 60
HARVEST_TIMEOUT_S = 600
ALIGNMENT_TOLERANCE_FRACTION = 1.0e-6


@dataclass(frozen=True, slots=True)
class ExportLayout:
    root: Path
    input: Path
    harvest: Path
    flow: Path
    qc: Path
    proteus: Path
    figures: Path


@dataclass(frozen=True, slots=True)
class AsciiSpatialField:
    coordinates_m: np.ndarray
    pressure_pa: np.ndarray
    velocity_m_s: np.ndarray
    columns: tuple[str, ...]
    result_path: Path
    companion_path: Path


@dataclass(frozen=True, slots=True)
class LatticeMapping:
    cell_indices: np.ndarray
    maximum_alignment_error_m: float
    duplicate_cell_count: int
    unique_cell_count: int


_NUMBER = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"
_HEX_CORNERS = np.asarray(
    (
        (0, 0, 0),
        (1, 0, 0),
        (1, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0, 1, 1),
    ),
    dtype=np.int64,
)


def _create_layout(output_root: Path) -> ExportLayout:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    root = output_root / f"{EXPORT_PREFIX}_{stamp}"
    if root.exists():
        raise FlowError("CFD_FLOW_STEADY_FIELD_EXPORT_INVALID", f"Output exists: {root}")
    directories = {name: root / name for name in ("input", "harvest", "flow", "qc", "proteus", "figures")}
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=False)
    return ExportLayout(root=root, **directories)


def _git_value(root: Path, *arguments: str) -> str:
    process = subprocess.run(
        ["git", *arguments],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )
    if process.returncode != 0:
        raise FlowError("CFD_FLOW_STEADY_FIELD_EXPORT_INVALID", process.stderr.strip())
    return process.stdout.strip()


def _extract_lua_string(text: str, key: str) -> str:
    match = re.search(rf"(?m)^\s*{re.escape(key)}\s*=\s*['\"]([^'\"]+)['\"]", text)
    if match is None:
        raise FlowError("CFD_FLOW_STEADY_RESTART_INVALID", f"Missing {key} in Lua")
    return match.group(1)


def _extract_lua_integer(text: str, key: str) -> int:
    match = re.search(rf"(?m)^\s*{re.escape(key)}\s*=\s*(\d+)\s*,?\s*$", text)
    if match is None:
        raise FlowError("CFD_FLOW_STEADY_RESTART_INVALID", f"Missing {key} in Lua")
    return int(match.group(1))


def _wsl_to_windows(path: str) -> Path:
    match = re.fullmatch(r"/mnt/([A-Za-z])/(.+)", path)
    if match is None:
        raise FlowError("CFD_FLOW_STEADY_RESTART_INVALID", f"Unsupported WSL path: {path}")
    drive, relative = match.groups()
    return Path(f"{drive.upper()}:/{relative}").resolve()


def _lua_quote(value: str) -> str:
    return value.replace("\\", "/").replace("'", "\\'")


def generate_ascii_spatial_harvester_lua(
    *,
    solver_config_wsl: str,
    restart_header_wsl: str,
) -> str:
    """Require the exact solver config, then override only restart/tracking."""

    solver_path = Path(solver_config_wsl)
    module = solver_path.stem
    package_pattern = solver_path.parent.as_posix().rstrip("/") + "/?.lua"
    return f"""-- One-shot full field export from the frozen project-steady restart.
package.path = '{_lua_quote(package_pattern)};' .. package.path
require '{_lua_quote(module)}'

restart = {{
  read = '{_lua_quote(restart_header_wsl)}'
}}

tracking = {{
  label = '{ASCII_LABEL}',
  folder = '{ASCII_FOLDER_WSL}',
  variable = {{ 'pressure_phy', 'velocity_phy' }},
  shape = {{ kind = 'all' }},
  output = {{ format = 'asciiSpatial', use_get_point = false }}
}}
"""


def harvester_lua_contract(
    text: str,
    *,
    solver_config_wsl: str,
    restart_header_wsl: str,
) -> dict[str, Any]:
    solver_path = Path(solver_config_wsl)
    checks = {
        "requires_exact_solver_module": f"require '{solver_path.stem}'" in text,
        "solver_package_path_exact": solver_path.parent.as_posix() in text,
        "restart_exact": restart_header_wsl in text,
        "tracking_variables_exact": "variable = { 'pressure_phy', 'velocity_phy' }" in text,
        "shape_all": "shape = { kind = 'all' }" in text,
        "ascii_spatial_exact": "format = 'asciiSpatial'" in text,
        "short_folder_exact": f"folder = '{ASCII_FOLDER_WSL}'" in text,
        "short_label_exact": f"label = '{ASCII_LABEL}'" in text,
        "element_sampling": "use_get_point = false" in text,
        "no_vtk_output": "format = 'vtk'" not in text,
        "no_extra_derived_fields": not any(
            name in text for name in ("vel_mag_phy", "vorticity", "q_criterion", "wss")
        ),
    }
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def ascii_path_preflight(simulation_name: str, timestamp: str) -> dict[str, Any]:
    actual = (
        f"{ASCII_FOLDER_WSL}{simulation_name}_{ASCII_LABEL}_p00000_t{timestamp}.res"
    )
    lengths = {
        "actual_result": len(actual),
        "companion_header": len(
            f"{ASCII_FOLDER_WSL}{simulation_name}_{ASCII_LABEL}.lua"
        ),
    }
    maximum = max(lengths.values())
    return {
        "status": "PASS" if maximum <= ASCII_MAX_PREDICTED_LENGTH else "FAIL",
        "treelm_label_len": 80,
        "project_maximum": ASCII_MAX_PREDICTED_LENGTH,
        "predicted_paths": {"actual": actual},
        "lengths": lengths,
        "maximum_predicted_filename_length": maximum,
    }


def parse_uniform_mesh_lattice(header_path: Path) -> dict[str, Any]:
    text = header_path.read_text(encoding="utf-8", errors="replace")
    box = re.search(
        rf"boundingbox\s*=\s*\{{\s*origin\s*=\s*\{{(.*?)\}}\s*,?\s*length\s*=\s*({_NUMBER})",
        text,
        re.DOTALL,
    )
    if box is None:
        raise FlowError("CFD_FLOW_FROZEN_SEEDER_MESH_INVALID", "Missing boundingbox")
    level_match = re.search(r"(?m)^\s*minLevel\s*=\s*(\d+)", text)
    max_level_match = re.search(r"(?m)^\s*maxLevel\s*=\s*(\d+)", text)
    if level_match is None or max_level_match is None:
        raise FlowError("CFD_FLOW_FROZEN_SEEDER_MESH_INVALID", "Incomplete lattice header")
    origin_values = [float(value) for value in re.findall(_NUMBER, box.group(1))]
    if len(origin_values) != 3:
        raise FlowError("CFD_FLOW_FROZEN_SEEDER_MESH_INVALID", "Invalid lattice origin")
    minimum_level = int(level_match.group(1))
    maximum_level = int(max_level_match.group(1))
    side = float(box.group(2))
    dx = side / (2**minimum_level)
    n_elems = _extract_lua_integer(text, "nElems")
    return {
        "status": "PASS" if minimum_level == maximum_level else "FAIL",
        "origin_m": origin_values,
        "side_m": side,
        "minimum_level": minimum_level,
        "maximum_level": maximum_level,
        "dx_m": dx,
        "n_elems": n_elems,
    }


def _companion_variables(text: str) -> dict[str, int]:
    matches = re.findall(
        r"name\s*=\s*['\"]([^'\"]+)['\"].*?ncomponents\s*=\s*(\d+)",
        text,
        re.DOTALL,
    )
    return {name: int(components) for name, components in matches}


def parse_ascii_spatial_schema(result_path: Path, companion_path: Path) -> dict[str, Any]:
    companion = companion_path.read_text(encoding="utf-8", errors="replace")
    format_match = re.search(r"format\s*=\s*['\"]([^'\"]+)['\"]", companion)
    count_match = re.search(r"(?m)^\s*nElems\s*=\s*(\d+)", companion)
    variables = _companion_variables(companion)
    header_line = ""
    with result_path.open(encoding="utf-8", errors="replace") as stream:
        for _ in range(8):
            line = stream.readline()
            if not line:
                break
            if line.lstrip().startswith("#") and "coordX" in line:
                header_line = line
                break
    if not header_line:
        raise FlowError(SCHEMA_INVALID, "asciiSpatial result has no coordinate header")
    columns = tuple(header_line.lstrip().lstrip("#").split())
    required_variables = {"pressure_phy": 1, "velocity_phy": 3}
    expected_field_columns = {
        "pressure_phy",
        "velocity_phy_01",
        "velocity_phy_02",
        "velocity_phy_03",
    }
    checks = {
        "companion_format_ascii_spatial": bool(
            format_match and format_match.group(1).lower() == "asciispatial"
        ),
        "companion_n_elems_present": count_match is not None,
        "companion_variables_exact": variables == required_variables,
        "coordinate_columns_exact": {"coordX", "coordY", "coordZ"}.issubset(columns),
        "field_columns_exact": set(columns) - {"coordX", "coordY", "coordZ"}
        == expected_field_columns,
        "columns_unique": len(columns) == len(set(columns)),
        "seven_columns_exact": len(columns) == 7,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "columns": list(columns),
        "companion_variables": variables,
        "companion_n_elems": int(count_match.group(1)) if count_match else None,
    }
    if result["status"] != "PASS":
        raise FlowError(SCHEMA_INVALID, f"Ambiguous asciiSpatial metadata: {result}")
    return result


def read_ascii_spatial_field(result_path: Path, companion_path: Path) -> AsciiSpatialField:
    schema = parse_ascii_spatial_schema(result_path, companion_path)
    columns = tuple(schema["columns"])
    values = np.loadtxt(result_path, comments="#", dtype=np.float64, ndmin=2)
    if values.shape[1] != len(columns):
        raise FlowError(SCHEMA_INVALID, f"Data shape {values.shape} does not match header")
    positions = {name: columns.index(name) for name in columns}
    coordinates = values[:, [positions["coordX"], positions["coordY"], positions["coordZ"]]]
    pressure = values[:, positions["pressure_phy"]]
    velocity = values[
        :,
        [
            positions["velocity_phy_01"],
            positions["velocity_phy_02"],
            positions["velocity_phy_03"],
        ],
    ]
    return AsciiSpatialField(
        coordinates_m=coordinates,
        pressure_pa=pressure,
        velocity_m_s=velocity,
        columns=columns,
        result_path=result_path,
        companion_path=companion_path,
    )


def quantize_cell_centers(
    coordinates_m: np.ndarray,
    *,
    origin_m: np.ndarray,
    dx_m: float,
) -> LatticeMapping:
    coordinates = np.asarray(coordinates_m, dtype=float)
    origin = np.asarray(origin_m, dtype=float)
    scaled = (coordinates - origin[None, :]) / dx_m - 0.5
    indices = np.rint(scaled).astype(np.int64)
    reconstructed = origin[None, :] + (indices.astype(float) + 0.5) * dx_m
    errors = np.linalg.norm(coordinates - reconstructed, axis=1)
    maximum = float(np.max(errors)) if len(errors) else math.inf
    unique = np.unique(indices, axis=0)
    duplicates = int(len(indices) - len(unique))
    return LatticeMapping(
        cell_indices=indices,
        maximum_alignment_error_m=maximum,
        duplicate_cell_count=duplicates,
        unique_cell_count=int(len(unique)),
    )


def reconstruct_hexahedral_field(
    *,
    mapping: LatticeMapping,
    origin_m: np.ndarray,
    dx_m: float,
    pressure_pa: np.ndarray,
    velocity_m_s: np.ndarray,
    pressure_reference_pa: float = PRESSURE_REFERENCE_PA,
) -> pv.UnstructuredGrid:
    indices = np.asarray(mapping.cell_indices, dtype=np.int64)
    corner_indices = indices[:, None, :] + _HEX_CORNERS[None, :, :]
    unique_corners, inverse = np.unique(
        corner_indices.reshape(-1, 3), axis=0, return_inverse=True
    )
    connectivity = inverse.reshape(-1, 8)
    cells = np.column_stack((np.full(len(indices), 8, dtype=np.int64), connectivity)).reshape(-1)
    cell_types = np.full(len(indices), int(pv.CellType.HEXAHEDRON), dtype=np.uint8)
    points = np.asarray(origin_m, dtype=float)[None, :] + unique_corners * dx_m
    grid = pv.UnstructuredGrid(cells, cell_types, points)
    grid.cell_data["velocity_phy"] = np.asarray(velocity_m_s, dtype=float)
    grid.cell_data["pressure_phy"] = np.asarray(pressure_pa, dtype=float)
    grid.cell_data["pressure_gauge_pa"] = np.asarray(pressure_pa, dtype=float) - float(
        pressure_reference_pa
    )
    return grid


def _percentiles(values: np.ndarray, requested: tuple[int, ...]) -> dict[str, float]:
    array = np.asarray(values, dtype=float).reshape(-1)
    result = {"min": float(np.min(array)), "max": float(np.max(array))}
    for percentile in requested:
        result[f"p{percentile:02d}"] = float(np.percentile(array, percentile))
    return result


def build_proteus_metadata(
    *,
    flow_vtu: Path,
    inlet_area_m2: float,
    dx_m: float,
) -> dict[str, Any]:
    return {
        "status": "PROTEUS_FIELD_CONTRACT_PASS",
        "lengthUnit": 1.0,
        "velocityUnit": 1.0,
        "velocityField": "velocity_phy",
        "pressureField": "pressure_phy",
        "gaugePressureField": "pressure_gauge_pa",
        "coordinateUnit": "m",
        "velocityUnitName": "m/s",
        "pressureUnitName": "Pa",
        "flowFieldVTU": str(flow_vtu),
        "inletDiameter": float(math.sqrt(4.0 * inlet_area_m2 / math.pi)),
        "inletNormal": None,
        "inlet_normal_policy": "AUTO_DETECT_BY_BACKPROPAGATION_LATER",
        "topology": "UNIFORM_CARTESIAN_VTK_HEXAHEDRON",
        "uniformDxM": float(dx_m),
        "gridConvergenceCompleted": False,
        "microbubbleSimulationRun": False,
        "backpropagationRun": False,
    }


def _file_manifest(paths: list[Path]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
        result[str(path.resolve())] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
    return result


def _historical_inventory_digest(output_root: Path) -> dict[str, Any]:
    rows: list[str] = []
    total_bytes = 0
    for path in sorted(item for item in output_root.rglob("*") if item.is_file()):
        relative = path.relative_to(output_root).as_posix()
        if relative.split("/", 1)[0].startswith(EXPORT_PREFIX):
            continue
        stat = path.stat()
        rows.append(f"{relative}|{stat.st_size}|{stat.st_mtime_ns}")
        total_bytes += stat.st_size
    digest = hashlib.sha256("\n".join(rows).encode("utf-8")).hexdigest()
    return {"file_count": len(rows), "total_bytes": total_bytes, "stat_digest_sha256": digest}


def _tail(text: str, lines: int = 40) -> list[str]:
    return text.splitlines()[-lines:]


def _copy_ascii_outputs(distribution: str, destination: Path) -> tuple[Path, Path]:
    temporary = Path(rf"\\wsl.localhost\{distribution}\tmp\u3d\x")
    companions = sorted(temporary.glob(f"*_{ASCII_LABEL}.lua"))
    results = sorted(temporary.glob(f"*_{ASCII_LABEL}_p*_t*.res"))
    if len(companions) != 1 or len(results) != 1:
        raise FlowError(
            HARVEST_FAILED,
            f"Expected one companion and one result, found companions={len(companions)} results={len(results)}",
        )
    companion = companions[0]
    destination.mkdir(parents=True, exist_ok=True)
    copied_companion = destination / companion.name
    copied_result = destination / results[0].name
    shutil.copy2(companion, copied_companion)
    shutil.copy2(results[0], copied_result)
    return copied_result, copied_companion


def _scatter_figure(
    coordinates_m: np.ndarray,
    values: np.ndarray,
    path: Path,
    *,
    title: str,
    colorbar_label: str,
) -> None:
    stride = max(1, len(coordinates_m) // 60_000)
    points = coordinates_m[::stride] * 1.0e6
    scalars = np.asarray(values)[::stride]
    figure = plt.figure(figsize=(9, 7), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    artist = axis.scatter(
        points[:, 0], points[:, 1], points[:, 2], c=scalars, s=1.0, cmap="viridis"
    )
    figure.colorbar(artist, ax=axis, shrink=0.7, label=colorbar_label)
    axis.set_title(title)
    axis.set_xlabel("x (um)")
    axis.set_ylabel("y (um)")
    axis.set_zlabel("z (um)")
    figure.savefig(path, dpi=170)
    plt.close(figure)


def _port_diagnostics(
    fluxes: dict[str, float],
    partition: SurfacePartition,
    kinematic_viscosity_m2_s: float,
) -> dict[str, Any]:
    reynolds = reynolds_diagnostics(fluxes, partition, kinematic_viscosity_m2_s)
    result: dict[str, Any] = {}
    for patch in partition.patches:
        if patch.label == "wall":
            continue
        area = patch.area_um2 * 1.0e-12
        diameter = 2.0 * patch.equivalent_radius_um * 1.0e-6
        result[patch.label] = {
            "q_m3_s": float(fluxes[patch.label]),
            "area_m2": float(area),
            "mean_velocity_m_s": float(abs(fluxes[patch.label]) / area),
            "equivalent_diameter_m": float(diameter),
            "reynolds": float(reynolds[patch.label]),
        }
    return result


def _failure_manifest(
    path: Path,
    summary: dict[str, Any],
    error: Exception,
) -> dict[str, Any]:
    if isinstance(error, FlowError):
        status = error.status
    else:
        status = "CFD_FLOW_STEADY_FIELD_EXPORT_INVALID"
    summary.update(
        {
            "status": status,
            "failure": str(error),
            "next": "REVIEW STEADY FIELD EXPORT EVIDENCE WITHOUT RERUN",
            "completed_at": datetime.now().isoformat(),
        }
    )
    write_json(path, summary)
    return summary


def _existing_revision(output_root: Path) -> dict[str, Any] | None:
    for run_root in sorted(output_root.glob(f"{EXPORT_PREFIX}_*"), reverse=True):
        reference = run_root / "input" / "source_reference.json"
        manifest = run_root / "qc" / "steady_field_export_manifest.json"
        if not reference.is_file() or not manifest.is_file():
            continue
        if read_json(reference).get("export_revision") == EXPORT_REVISION:
            return read_json(manifest)
    return None


def run_steady_field_export(project_root: Path) -> dict[str, Any]:
    """Run no solver and at most one pinned single-process ASCII harvester."""

    root = Path(project_root).resolve()
    config = load_cfd_flow_config(root / "configs" / "cfd_flow.yaml", project_root=root)
    output_root = config.paths.output_root
    previous = _existing_revision(output_root)
    if previous is not None:
        return previous

    branch = _git_value(root, "branch", "--show-current")
    head = _git_value(root, "rev-parse", "HEAD")
    if head != SOURCE_COMMIT:
        raise FlowError(
            "CFD_FLOW_STEADY_FIELD_EXPORT_SOURCE_INVALID",
            f"Expected source commit {SOURCE_COMMIT}, found {head}",
        )
    source_run = output_root / SOURCE_RUN_NAME
    frozen_seeder = output_root / FROZEN_SEEDER_NAME
    solver_workdir = output_root / SOLVER_WORKDIR_NAME / "musubi"
    restart_header = source_run / "restart" / "roi003274_steady_lbm_lastHeader.lua"
    source_summary = read_json(source_run / "qc" / "project_steady_confirmation.json")
    steady_restart = read_json(source_run / "qc" / "steady_restart_manifest.json")
    if (
        source_summary.get("status") != "CFD_FLOW_PROJECT_STEADY_0P1PCT_CONFIRMED"
        or not source_summary["official_steady_termination"]["official_steady_termination"]
        or int(source_summary["official_steady_termination"]["confirmation_iteration"])
        != FROZEN_ITERATION
        or steady_restart.get("status") != SOURCE_STATE
    ):
        raise FlowError("CFD_FLOW_STEADY_RESTART_INVALID", "Frozen steady evidence mismatch")
    restart_qc = validate_continuation_restart(
        restart_header,
        solver_workdir=solver_workdir,
        frozen_seeder_run=frozen_seeder,
        expected_iteration=FROZEN_ITERATION,
    )
    header_text = restart_header.read_text(encoding="utf-8", errors="replace")
    solver_config_wsl = _extract_lua_string(header_text, "solver_configFile")
    solver_config = _wsl_to_windows(solver_config_wsl)
    if not solver_config.is_file():
        raise FlowError("CFD_FLOW_STEADY_RESTART_INVALID", "Header solver_configFile is missing")
    restart_header_wsl = windows_to_wsl(restart_header, config.apes.wsl_distribution)
    binary_wsl = re.search(r"binary_name\s*=\s*\{\s*['\"]([^'\"]+)", header_text, re.DOTALL)
    timestamp_token = "4.836E-03"
    if binary_wsl:
        match = re.search(r"_([0-9.]+E[-+][0-9]+)\.lsb$", binary_wsl.group(1))
        if match:
            timestamp_token = match.group(1)
    simulation_name = _extract_lua_string(solver_config.read_text(encoding="utf-8"), "simulation_name")
    path_qc = ascii_path_preflight(simulation_name, timestamp_token)
    if path_qc["status"] != "PASS":
        raise FlowError(PATH_INVALID, f"ASCII path length failed: {path_qc}")
    mesh_dir = frozen_seeder / "seeder" / "mesh"
    lattice = parse_uniform_mesh_lattice(mesh_dir / "header.lua")
    scaling_source = read_json(frozen_seeder / "qc" / "lbm_scaling_qc.json")
    dt_s = float(scaling_source["dt_s"])
    pressure_reference_pa = float(scaling_source["pressure_reference_pa"])
    if (
        lattice["status"] != "PASS"
        or lattice["n_elems"] != EXPECTED_CELL_COUNT
        or lattice["minimum_level"] != 9
        or not math.isclose(lattice["dx_m"], 2.0e-7, rel_tol=0.0, abs_tol=1.0e-21)
        or not math.isclose(float(scaling_source["dx_m"]), lattice["dx_m"], rel_tol=0.0, abs_tol=1.0e-21)
        or not math.isclose(pressure_reference_pa, PRESSURE_REFERENCE_PA, rel_tol=0.0, abs_tol=1.0e-10)
    ):
        raise FlowError(
            "CFD_FLOW_FROZEN_SEEDER_MESH_INVALID",
            f"Lattice/scaling mismatch: lattice={lattice}, scaling={scaling_source}",
        )

    critical_paths = [
        restart_header,
        Path(restart_qc["binary"]),
        solver_config,
        source_run / "qc" / "project_steady_confirmation.json",
        source_run / "qc" / "steady_restart_manifest.json",
        frozen_seeder / "seeder" / "mesh_manifest.json",
        *(path for path in sorted(mesh_dir.iterdir()) if path.is_file()),
    ]
    critical_before = _file_manifest(critical_paths)
    historical_before = _historical_inventory_digest(output_root)
    layout = _create_layout(output_root)
    manifest_path = layout.qc / "steady_field_export_manifest.json"
    summary: dict[str, Any] = {
        "status": "CFD_FLOW_STEADY_FIELD_EXPORT_PREFLIGHT",
        "export_revision": EXPORT_REVISION,
        "branch": branch,
        "source_commit": head,
        "run_root": str(layout.root),
        "source_state": SOURCE_STATE,
        "source_run": str(source_run),
        "frozen_iteration": FROZEN_ITERATION,
        "frozen_n_elems": EXPECTED_CELL_COUNT,
        "frozen_restart": str(restart_header),
        "frozen_restart_validation": restart_qc,
        "steady_evidence": {
            "status": "CFD_FLOW_PROJECT_STEADY_0P1PCT_CONFIRMED",
            "official_steady_termination": True,
            "confirmation_iteration": FROZEN_ITERATION,
            "temporal_steady_state": "FROZEN_ACCEPTED",
            "criterion_policy_name": "PROJECT_CHARACTERISTIC_SCALE_STEADY_0P1_PERCENT",
            "legacy_heuristic_threshold": {"pressure_pa": 0.001, "velocity_m_s": 1.0e-9},
            "project_steady_threshold": {
                "pressure_pa": PROJECT_PRESSURE_THRESHOLD_PA,
                "velocity_m_s": PROJECT_VELOCITY_THRESHOLD_M_S,
            },
            "project_steady_ratio": {
                "pressure": PROJECT_PRESSURE_RATIO,
                "velocity": PROJECT_VELOCITY_RATIO,
            },
        },
        "solver_config_from_restart_header": str(solver_config),
        "solver_config_sha256": sha256_file(solver_config),
        "seeder_run_count": 0,
        "musubi_run_count": 0,
        "harvester_run_count": 0,
        "grid_sweep_performed": False,
        "physics_or_bc_modified": False,
        "source_steady_restart_modified": False,
        "tetrahedral_mesh_created": False,
        "microbubble_simulation_run": False,
        "final_cfd_solution": False,
        "grid_convergence_pending": True,
        "harvester_output_format": "asciiSpatial",
        "ascii_path_preflight": path_qc,
        "lattice": lattice,
        "critical_source_manifest_before": critical_before,
        "historical_output_inventory_before": historical_before,
        "next": "RUN ONE MUS_HARVESTING ASCII-SPATIAL EXPORT",
        "started_at": datetime.now().isoformat(),
    }
    write_json(manifest_path, summary)
    write_json(
        layout.input / "source_reference.json",
        {
            "status": "PASS",
            "export_revision": EXPORT_REVISION,
            "source_state": SOURCE_STATE,
            "source_commit": head,
            "restart_header": str(restart_header),
            "restart_header_sha256": sha256_file(restart_header),
            "solver_config_from_restart_header": str(solver_config),
            "solver_config_sha256": sha256_file(solver_config),
            "mesh_header": str(mesh_dir / "header.lua"),
            "mesh_header_sha256": sha256_file(mesh_dir / "header.lua"),
        },
    )
    write_json(layout.qc / "ascii_path_preflight.json", path_qc)
    write_json(layout.qc / "frozen_restart_validation.json", restart_qc)
    write_json(layout.qc / "uniform_lattice_source.json", lattice)

    try:
        inputs = load_flow_inputs(config.paths.source_surface_run)
        partition = load_frozen_surface_partition(
            inputs,
            frozen_seeder / "geometry" / "geometry_solver_m",
        )
        bc = load_boundary_conditions(inputs.boundary_conditions)
        environment = inspect_apes_environment(config.apes)
        if environment.status != "PASS" or not environment.binaries.get("mus_harvesting"):
            raise FlowError("CFD_FLOW_ENVIRONMENT_BLOCKED", "Pinned mus_harvesting unavailable")
        write_json(layout.input / "environment.json", asdict(environment))
        harvester_lua = generate_ascii_spatial_harvester_lua(
            solver_config_wsl=solver_config_wsl,
            restart_header_wsl=restart_header_wsl,
        )
        harvester_config = layout.input / "steady_ascii_harvester.lua"
        harvester_config.write_text(harvester_lua, encoding="utf-8")
        contract = harvester_lua_contract(
            harvester_lua,
            solver_config_wsl=solver_config_wsl,
            restart_header_wsl=restart_header_wsl,
        )
        write_json(layout.input / "harvester_lua_contract.json", contract)
        if contract["status"] != "PASS":
            raise FlowError("CFD_FLOW_STEADY_HARVEST_CONFIG_INVALID", "Lua contract failed")
        luac = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=layout.root,
            command=[
                str(environment.binaries["lua_compiler"]),
                "-p",
                windows_to_wsl(harvester_config, config.apes.wsl_distribution),
            ],
            stdout_path=layout.input / "luac_stdout.log",
            stderr_path=layout.input / "luac_stderr.log",
            timeout_s=30,
        )
        if luac.returncode != 0:
            raise FlowError("CFD_FLOW_STEADY_HARVEST_CONFIG_INVALID", "Lua syntax failed")
        cleanup_before = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=layout.root,
            command=["rm", "-rf", "--", ASCII_FOLDER_WSL.rstrip("/")],
            stdout_path=layout.input / "tmp_cleanup_before_stdout.log",
            stderr_path=layout.input / "tmp_cleanup_before_stderr.log",
            timeout_s=30,
        )
        mkdir = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=layout.root,
            command=["mkdir", "-p", "--", ASCII_FOLDER_WSL],
            stdout_path=layout.input / "tmp_mkdir_stdout.log",
            stderr_path=layout.input / "tmp_mkdir_stderr.log",
            timeout_s=30,
        )
        if cleanup_before.returncode != 0 or mkdir.returncode != 0:
            raise FlowError("CFD_FLOW_STEADY_EXPORT_TEMP_INVALID", "Could not prepare /tmp/u3d/x")
        summary.update(
            {
                "status": "CFD_FLOW_STEADY_ASCII_HARVEST_RUNNING",
                "harvester_run_count": 1,
                "harvester_execution": "PINNED_OFFICIAL_SINGLE_PROCESS",
                "next": "NO_AUTOMATIC_RETRY",
            }
        )
        write_json(manifest_path, summary)
        harvest_run = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=solver_workdir,
            command=[
                str(environment.binaries["mus_harvesting"]),
                windows_to_wsl(harvester_config, config.apes.wsl_distribution),
            ],
            stdout_path=layout.harvest / "harvester_stdout.log",
            stderr_path=layout.harvest / "harvester_stderr.log",
            timeout_s=HARVEST_TIMEOUT_S,
        )
        stdout = harvest_run.stdout_path.read_text(encoding="utf-8", errors="replace")
        stderr = harvest_run.stderr_path.read_text(encoding="utf-8", errors="replace")
        temp_dir = Path(rf"\\wsl.localhost\{config.apes.wsl_distribution}\tmp\u3d\x")
        temp_results = sorted(temp_dir.glob("*")) if temp_dir.exists() else []
        restart_read = "Done reading restart." in stdout
        ascii_started = any(path.is_file() for path in temp_results)
        preserved_partial_files: list[str] = []
        if harvest_run.returncode != 0:
            for path in temp_results:
                if path.is_file():
                    target = layout.harvest / path.name
                    shutil.copy2(path, target)
                    preserved_partial_files.append(str(target))
        signal = None
        if harvest_run.returncode in (-11, 139) or "SIGSEGV" in stdout + stderr:
            signal = "SIGSEGV"
        harvester_evidence = {
            "status": "PASS" if harvest_run.returncode == 0 else "FAIL",
            "return_code": harvest_run.returncode,
            "signal": signal,
            "wall_time_s": harvest_run.wall_time_s,
            "restart_read_completed": restart_read,
            "ascii_spatial_file_started_writing": ascii_started,
            "temporary_files": [path.name for path in temp_results],
            "preserved_partial_files": preserved_partial_files,
            "last_stdout_lines": _tail(stdout),
            "last_stderr_lines": _tail(stderr),
        }
        write_json(layout.qc / "harvester_run_qc.json", harvester_evidence)
        summary["harvester"] = harvester_evidence
        write_json(manifest_path, summary)
        if harvest_run.returncode != 0:
            raise FlowError(HARVEST_FAILED, f"mus_harvesting return code {harvest_run.returncode}")
        if not restart_read:
            raise FlowError(HARVEST_FAILED, "Harvester did not confirm restart read completion")
        result_path, companion_path = _copy_ascii_outputs(
            config.apes.wsl_distribution,
            layout.harvest,
        )
        cleanup_after = run_wsl_tool(
            distribution=config.apes.wsl_distribution,
            workdir=layout.root,
            command=["rm", "-rf", "--", ASCII_FOLDER_WSL.rstrip("/")],
            stdout_path=layout.input / "tmp_cleanup_after_stdout.log",
            stderr_path=layout.input / "tmp_cleanup_after_stderr.log",
            timeout_s=30,
        )
        if cleanup_after.returncode != 0:
            raise FlowError("CFD_FLOW_STEADY_EXPORT_TEMP_INVALID", "Could not clean /tmp/u3d/x")

        schema = parse_ascii_spatial_schema(result_path, companion_path)
        write_json(layout.qc / "ascii_schema_qc.json", schema)
        field = read_ascii_spatial_field(result_path, companion_path)
        sample_count = len(field.coordinates_m)
        finite = {
            "coordinates": bool(np.all(np.isfinite(field.coordinates_m))),
            "pressure": bool(np.all(np.isfinite(field.pressure_pa))),
            "velocity": bool(np.all(np.isfinite(field.velocity_m_s))),
            "nan_count": int(
                np.isnan(field.coordinates_m).sum()
                + np.isnan(field.pressure_pa).sum()
                + np.isnan(field.velocity_m_s).sum()
            ),
            "inf_count": int(
                np.isinf(field.coordinates_m).sum()
                + np.isinf(field.pressure_pa).sum()
                + np.isinf(field.velocity_m_s).sum()
            ),
        }
        if sample_count != EXPECTED_CELL_COUNT or not all(
            finite[key] for key in ("coordinates", "pressure", "velocity")
        ):
            raise FlowError(
                "CFD_FLOW_STEADY_FIELD_IDENTITY_INVALID",
                f"Samples={sample_count}, finite={finite}",
            )
        origin = np.asarray(lattice["origin_m"], dtype=float)
        dx_m = float(lattice["dx_m"])
        root_lower = origin
        root_upper = origin + float(lattice["side_m"])
        coordinate_bounds = np.column_stack(
            (np.min(field.coordinates_m, axis=0), np.max(field.coordinates_m, axis=0))
        )
        coordinates_in_root = bool(
            np.all(field.coordinates_m >= root_lower[None, :])
            and np.all(field.coordinates_m <= root_upper[None, :])
        )
        if not coordinates_in_root or float(np.max(np.abs(field.coordinates_m))) >= 1.0e-2:
            raise FlowError("CFD_FLOW_STEADY_FIELD_COORDINATE_UNIT_INVALID", "Coordinates not meters")
        mapping = quantize_cell_centers(field.coordinates_m, origin_m=origin, dx_m=dx_m)
        alignment_limit = ALIGNMENT_TOLERANCE_FRACTION * dx_m
        lattice_qc = {
            "status": "PASS"
            if mapping.maximum_alignment_error_m <= alignment_limit
            and mapping.duplicate_cell_count == 0
            and mapping.unique_cell_count == EXPECTED_CELL_COUNT
            else "FAIL",
            "coordinate_unit": "m",
            "coordinates_in_frozen_root_bounding_box": coordinates_in_root,
            "coordinate_bounds_m": coordinate_bounds.tolist(),
            "origin_m": origin.tolist(),
            "dx_m": dx_m,
            "level": int(lattice["minimum_level"]),
            "maximum_center_alignment_error_m": mapping.maximum_alignment_error_m,
            "alignment_limit_m": alignment_limit,
            "unique_cartesian_cell_count": mapping.unique_cell_count,
            "duplicate_cell_index_count": mapping.duplicate_cell_count,
        }
        write_json(layout.qc / "lattice_alignment_qc.json", lattice_qc)
        if lattice_qc["status"] != "PASS":
            raise FlowError("CFD_FLOW_STEADY_FIELD_LATTICE_INVALID", f"Lattice QC: {lattice_qc}")

        grid = reconstruct_hexahedral_field(
            mapping=mapping,
            origin_m=origin,
            dx_m=dx_m,
            pressure_pa=field.pressure_pa,
            velocity_m_s=field.velocity_m_s,
            pressure_reference_pa=pressure_reference_pa,
        )
        flow_vtu = layout.flow / "flow_field.vtu"
        grid.save(flow_vtu, binary=True)
        reloaded = pv.read(flow_vtu).cast_to_unstructured_grid()
        cell_types, counts = np.unique(reloaded.celltypes, return_counts=True)
        cell_type_counts = {
            str(int(cell_type)): int(count)
            for cell_type, count in zip(cell_types, counts, strict=True)
        }
        speed = np.linalg.norm(np.asarray(reloaded.cell_data["velocity_phy"]), axis=1)
        gauge = np.asarray(reloaded.cell_data["pressure_gauge_pa"], dtype=float)
        xml_prefix = flow_vtu.read_bytes()[:8192]
        field_qc = {
            "status": "PASS",
            "source_ascii_spatial": str(result_path),
            "output_vtu": str(flow_vtu),
            "vtu_binary": b'format="binary"' in xml_prefix or b'format="appended"' in xml_prefix,
            "coordinate_unit": "m",
            "coordinate_bounds_m": list(reloaded.bounds),
            "dx_m": dx_m,
            "cell_count": int(reloaded.n_cells),
            "point_count": int(reloaded.n_points),
            "cell_type_counts": cell_type_counts,
            "cell_type": "VTK_HEXAHEDRON",
            "all_cells_vtk_hexahedron": bool(
                len(cell_types) == 1 and int(cell_types[0]) == int(pv.CellType.HEXAHEDRON)
            ),
            "pressure_phy_shape": list(np.asarray(reloaded.cell_data["pressure_phy"]).shape),
            "velocity_phy_shape": list(np.asarray(reloaded.cell_data["velocity_phy"]).shape),
            "pressure_gauge_pa_shape": list(gauge.shape),
            "velocity_phy_cell_data": "velocity_phy" in reloaded.cell_data,
            "pressure_phy_cell_data": "pressure_phy" in reloaded.cell_data,
            "pressure_gauge_pa_cell_data": "pressure_gauge_pa" in reloaded.cell_data,
            "finite": finite,
            "velocity_m_s": _percentiles(speed, (50, 95, 99)),
            "pressure_gauge_pa": _percentiles(gauge, (1, 50, 99)),
        }
        identity_pass = (
            field_qc["vtu_binary"]
            and reloaded.n_cells == EXPECTED_CELL_COUNT
            and field_qc["all_cells_vtk_hexahedron"]
            and field_qc["velocity_phy_shape"] == [EXPECTED_CELL_COUNT, 3]
            and field_qc["pressure_phy_shape"] == [EXPECTED_CELL_COUNT]
            and field_qc["pressure_gauge_pa_shape"] == [EXPECTED_CELL_COUNT]
        )
        field_qc["status"] = "PASS" if identity_pass else "FAIL"
        write_json(layout.qc / "field_identity_qc.json", field_qc)
        if not identity_pass:
            raise FlowError("CFD_FLOW_STEADY_FIELD_IDENTITY_INVALID", "VTU identity QC failed")
        actual_mach = float(np.max(speed) * dt_s / dx_m / math.sqrt(1.0 / 3.0))
        mach_qc = {
            "status": "PASS" if actual_mach < 0.05 else "FAIL",
            "maximum_velocity_m_s": float(np.max(speed)),
            "dx_m": dx_m,
            "dt_s": dt_s,
            "d3q19_cs": math.sqrt(1.0 / 3.0),
            "actual_lattice_mach": actual_mach,
            "maximum_allowed": 0.05,
        }
        write_json(layout.qc / "mach_qc.json", mach_qc)
        if mach_qc["status"] != "PASS":
            raise FlowError("CFD_FLOW_STEADY_FIELD_MACH_INVALID", "Actual Mach >= 0.05")

        try:
            fluxes, measured_pressures = numerical_port_fluxes(reloaded, partition, dx_m)
        except FlowError as error:
            flux_failure = {
                "status": "FAIL",
                "failure_kind": "PORT_INTEGRATION_OR_SAMPLING_FAILURE",
                "failure": str(error),
            }
            write_json(layout.qc / "port_flux_qc.json", flux_failure)
            raise FlowError(FLUX_FAILED, str(error)) from error
        q_in = float(fluxes["inlet"])
        q_out = [float(fluxes[f"outlet_{index:02d}"]) for index in range(1, 4)]
        inlet_error = abs(q_in - bc.inlet_flow_m3_s) / abs(bc.inlet_flow_m3_s)
        mass_error = abs(abs(q_in) - sum(abs(value) for value in q_out)) / abs(q_in)
        directions_pass = bool(q_in > 0.0 and all(value > 0.0 for value in q_out))
        port_details = _port_diagnostics(fluxes, partition, bc.kinematic_viscosity_m2_s)
        flux_pass = inlet_error <= 0.01 and mass_error <= 0.01 and directions_pass
        flux_qc = {
            "status": "PASS" if flux_pass else "FAIL",
            "failure_kind": None if flux_pass else "FIELD_ACTUAL_FLOW_FAILURE",
            "method": "existing numerical_port_fluxes: cap triangles translated inward 2dx with three-point barycentric quadrature",
            "q_in_m3_s": q_in,
            "q_out_m3_s": {
                f"outlet_{index:02d}": value for index, value in enumerate(q_out, start=1)
            },
            "inlet_target_m3_s": bc.inlet_flow_m3_s,
            "inlet_relative_error": inlet_error,
            "mass_conservation_error": mass_error,
            "flow_directions_pass": directions_pass,
            "maximum_allowed_inlet_relative_error": 0.01,
            "maximum_allowed_mass_conservation_error": 0.01,
            "port_diagnostics": port_details,
            "expected_1d_vs_3d_outlet_flow_diagnostic_only": [
                {
                    "label": f"outlet_{index:02d}",
                    "expected_1d_m3_s": expected,
                    "measured_3d_m3_s": q_out[index - 1],
                    "role": "DIAGNOSTIC_ONLY",
                }
                for index, expected in enumerate(
                    bc.outlet_expected_1d_flows_m3_s, start=1
                )
            ],
        }
        write_json(layout.qc / "port_flux_qc.json", flux_qc)
        pressure_qc = {
            "status": "DIAGNOSTIC",
            "method": "boundary-adjacent internal cap plane from numerical_port_fluxes",
            "outlets": [
                {
                    "label": f"outlet_{index:02d}",
                    "target_gauge_pa": target,
                    "measured_boundary_adjacent_gauge_pa": measured_pressures[
                        f"outlet_{index:02d}"
                    ],
                    "difference_pa": measured_pressures[f"outlet_{index:02d}"] - target,
                    "role": "DIAGNOSTIC_ONLY",
                }
                for index, target in enumerate(bc.outlet_gauge_pressures_pa, start=1)
            ],
        }
        write_json(layout.qc / "outlet_pressure_diagnostic.json", pressure_qc)
        if not flux_pass:
            raise FlowError(FLUX_FAILED, f"Flux gates failed: {flux_qc}")

        metadata = build_proteus_metadata(
            flow_vtu=flow_vtu,
            inlet_area_m2=partition.patch("inlet").area_um2 * 1.0e-12,
            dx_m=dx_m,
        )
        metadata_path = layout.proteus / "proteus_flow_metadata.json"
        write_json(metadata_path, metadata)
        proteus_qc = {
            "status": "PASS",
            "coordinates_meter": True,
            "velocity_phy_cell_data": field_qc["velocity_phy_cell_data"],
            "velocity_components": 3,
            "velocity_unit": "m/s",
            "pressure_phy_cell_data": field_qc["pressure_phy_cell_data"],
            "pressure_unit": "Pa",
            "cartesian_hexahedral_topology": field_qc["all_cells_vtk_hexahedron"],
            "uniform_dx_m": dx_m,
            "finite_values": finite["nan_count"] == 0 and finite["inf_count"] == 0,
            "inlet_normal": None,
            "inlet_normal_policy": "AUTO_DETECT_BY_BACKPROPAGATION_LATER",
        }
        write_json(layout.qc / "proteus_compatibility_qc.json", proteus_qc)
        velocity_figure = layout.figures / "velocity_magnitude_review.png"
        pressure_figure = layout.figures / "gauge_pressure_review.png"
        _scatter_figure(
            field.coordinates_m,
            speed,
            velocity_figure,
            title="Steady velocity magnitude",
            colorbar_label="|u| (m/s)",
        )
        _scatter_figure(
            field.coordinates_m,
            gauge,
            pressure_figure,
            title="Steady gauge pressure",
            colorbar_label="gauge pressure (Pa)",
        )

        critical_after = _file_manifest(critical_paths)
        historical_after = _historical_inventory_digest(output_root)
        source_unchanged = critical_before == critical_after
        historical_unchanged = historical_before == historical_after
        source_read_only_qc = {
            "status": "PASS" if source_unchanged and historical_unchanged else "FAIL",
            "critical_source_files_unchanged": source_unchanged,
            "historical_output_inventory_unchanged": historical_unchanged,
            "critical_source_manifest_after": critical_after,
            "historical_output_inventory_after": historical_after,
        }
        write_json(layout.qc / "source_read_only_qc.json", source_read_only_qc)
        if source_read_only_qc["status"] != "PASS":
            raise FlowError("CFD_FLOW_HISTORICAL_OUTPUT_MODIFIED", "Historical outputs changed")
        summary.update(
            {
                "status": SUCCESS_STATUS,
                "next": SUCCESS_NEXT,
                "harvester_output_format": "asciiSpatial",
                "ascii_schema": schema,
                "spatial_sample_count": sample_count,
                "unique_lattice_cell_count": mapping.unique_cell_count,
                "lattice_alignment": lattice_qc,
                "flow_vtu": str(flow_vtu),
                "field_identity": field_qc,
                "mach": mach_qc,
                "port_flux": flux_qc,
                "outlet_pressure": pressure_qc,
                "proteus_compatibility": proteus_qc,
                "proteus_metadata": str(metadata_path),
                "figures": [str(velocity_figure), str(pressure_figure)],
                "critical_source_manifest_after": critical_after,
                "historical_output_inventory_after": historical_after,
                "source_steady_restart_modified": False,
                "physics_or_bc_modified": False,
                "grid_sweep_performed": False,
                "microbubble_simulation_run": False,
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(manifest_path, summary)
        return summary
    except Exception as error:
        try:
            critical_after = _file_manifest(critical_paths)
            historical_after = _historical_inventory_digest(output_root)
            source_unchanged = critical_before == critical_after
            historical_unchanged = historical_before == historical_after
            summary.update(
                {
                    "critical_source_manifest_after": critical_after,
                    "historical_output_inventory_after": historical_after,
                    "source_steady_restart_modified": not source_unchanged,
                    "historical_outputs_modified": not historical_unchanged,
                }
            )
            write_json(
                layout.qc / "source_read_only_qc.json",
                {
                    "status": "PASS"
                    if source_unchanged and historical_unchanged
                    else "FAIL",
                    "critical_source_files_unchanged": source_unchanged,
                    "historical_output_inventory_unchanged": historical_unchanged,
                },
            )
        except Exception as audit_error:
            summary["source_read_only_audit_failure"] = str(audit_error)
        return _failure_manifest(manifest_path, summary, error)
