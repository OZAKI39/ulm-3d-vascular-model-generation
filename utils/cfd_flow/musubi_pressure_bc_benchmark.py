"""Small source-proven Musubi Hagen--Poiseuille pressure-BC benchmark.

The benchmark is isolated research code.  It creates a fixed physical pipe,
uses the accepted adaptive-flux inlet and wall_libb wall, and compares at most
the current pressure outlet plus two alternatives already present in the
pinned Musubi source.  It never launches the vascular model.
"""

from __future__ import annotations

import json
import math
import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import trimesh

from .io import sha256_file, write_json
from .musubi_boundary_mass_referee import load_mesh_contract
from .port_grid_sensitivity import RESEARCH_RUN, recover_continuous_ports


RHO_KG_M3 = 1056.0
NU_M2_S = 3.27e-6
BULK_NU_M2_S = 2.18e-6
MEAN_VELOCITY_M_S = 0.35e-3
TARGET_Q_M3_S = 2.7369132390905703e-15
TARGET_MASS_FLOW_KG_S = 2.890180380479642e-12
PRESSURE_REFERENCE_PA = 23622.32012800001
BASE_DX_M = 2.0e-7
BASE_DT_S = 2.44140625e-8
PIPE_AREA_M2 = TARGET_Q_M3_S / MEAN_VELOCITY_M_S
PIPE_RADIUS_M = math.sqrt(PIPE_AREA_M2 / math.pi)
PIPE_DIAMETER_M = 2.0 * PIPE_RADIUS_M
PIPE_LENGTH_M = 6.0 * PIPE_DIAMETER_M
CELLS_ACROSS = (12, 16, 20, 27)
ROOT_LEVEL = 8
MAXIMUM_ITERATIONS = 8_000
STEADY_WINDOW = 1_000

SEEDER_WSL = "/home/lzy/.local/bin/seeder"
MPIRUN_WSL = "/home/lzy/.local/bin/mpirun"
MUSUBI_WSL = (
    "/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/"
    "build/musubi_adaptive_flux"
)
MUSUBI_SHA256 = "e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588"
RUNTIME_ROOT_WSL = "/home/lzy/u3da/pressure_bc_benchmark_20260830_v2"

CANDIDATES = {
    "pressure_eq": "pressure_eq",
    "pressure_nonEqExpol": "pressure_noneq_expol",
    # The loader compares its normalized API token in lower case even though
    # the Fortran subroutine and tutorial heading use camel case.
    "pressure_antiBounceBack": "pressure_antibounceback",
}

CONTROLLER_PATTERN = re.compile(
    r"ADAPTIVE_FLUX_PRESSURE\s+iter=(?P<iteration>\d+)\s+"
    r"target_lattice=\s*(?P<target>[-+0-9.Ee]+)\s+"
    r"controlled_lattice=\s*(?P<controlled>[-+0-9.Ee]+)\s+"
    r"relative_error=\s*(?P<relative>[-+0-9.Ee]+)\s+"
    r"rho_boundary=\s*(?P<rho>[-+0-9.Ee]+)\s+"
    r"pressure_pa=\s*(?P<pressure>[-+0-9.Ee]+)\s+"
    r"max_lattice_velocity=\s*(?P<speed>[-+0-9.Ee]+)\s+"
    r"minimum_pdf=\s*(?P<minimum_pdf>[-+0-9.Ee]+)\s+"
    r"globBC_count=(?P<count>\d+)"
)


@dataclass(frozen=True, slots=True)
class PipeGrid:
    cells_across_diameter: int
    dx_m: float
    dt_s: float


def hagen_poiseuille_delta_p(
    *,
    length_m: float = PIPE_LENGTH_M,
    flow_m3_s: float = TARGET_Q_M3_S,
    radius_m: float = PIPE_RADIUS_M,
    dynamic_viscosity_pa_s: float = RHO_KG_M3 * NU_M2_S,
) -> float:
    if min(length_m, flow_m3_s, radius_m, dynamic_viscosity_pa_s) <= 0.0:
        raise ValueError("Hagen--Poiseuille inputs must be positive")
    return 8.0 * dynamic_viscosity_pa_s * length_m * flow_m3_s / (math.pi * radius_m**4)


def diffusive_time_step(dx_m: float) -> float:
    if dx_m <= 0.0:
        raise ValueError("dx must be positive")
    return BASE_DT_S * (float(dx_m) / BASE_DX_M) ** 2


def benchmark_grids() -> tuple[PipeGrid, ...]:
    return tuple(
        PipeGrid(count, PIPE_DIAMETER_M / count, diffusive_time_step(PIPE_DIAMETER_M / count))
        for count in CELLS_ACROSS
    )


def pipe_frame(direction: Iterable[float]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    axis = np.asarray(tuple(direction), dtype=np.float64).reshape(3)
    axis /= np.linalg.norm(axis)
    helper = np.asarray((1.0, 0.0, 0.0))
    if abs(float(np.dot(axis, helper))) > 0.85:
        helper = np.asarray((0.0, 1.0, 0.0))
    u = np.cross(axis, helper)
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    v /= np.linalg.norm(v)
    return axis, u, v


def worst_real_outlet_direction(project_root: Path) -> dict[str, Any]:
    frames, _ = recover_continuous_ports(project_root)
    candidates = {}
    for label in ("outlet_01", "outlet_02", "outlet_03"):
        direction = frames[label].normal
        # Lower maximum Cartesian component means farther from every axis.
        score = 1.0 - float(np.max(np.abs(direction)))
        candidates[label] = {"direction": direction.tolist(), "obliquity_score": score}
    label = max(candidates, key=lambda name: candidates[name]["obliquity_score"])
    return {"selected_port": label, **candidates[label], "all_outlets": candidates}


def generate_pipe_wall(
    path: Path,
    direction: Iterable[float],
    *,
    radius_m: float = PIPE_RADIUS_M,
    length_m: float = PIPE_LENGTH_M,
    sections: int = 128,
) -> dict[str, Any]:
    axis, u, v = pipe_frame(direction)
    inlet = -0.5 * length_m * axis
    outlet = 0.5 * length_m * axis
    angles = np.linspace(0.0, 2.0 * math.pi, int(sections), endpoint=False)
    ring = np.cos(angles)[:, None] * u + np.sin(angles)[:, None] * v
    vertices = np.vstack((inlet + radius_m * ring, outlet + radius_m * ring))
    faces: list[tuple[int, int, int]] = []
    for index in range(int(sections)):
        nxt = (index + 1) % int(sections)
        faces.append((index, nxt, int(sections) + nxt))
        faces.append((index, int(sections) + nxt, int(sections) + index))
    mesh = trimesh.Trimesh(vertices=vertices, faces=np.asarray(faces), process=False)
    path.parent.mkdir(parents=True, exist_ok=True)
    mesh.export(path)
    return {
        "path": str(path.resolve()),
        "sha256": sha256_file(path),
        "triangle_count": len(faces),
        "sections": int(sections),
        "axis": axis.tolist(),
        "inlet_center_m": inlet.tolist(),
        "outlet_center_m": outlet.tolist(),
        "radius_m": radius_m,
        "length_m": length_m,
    }


def _lua_vector(vector: Iterable[float]) -> str:
    return "{ " + ", ".join(f"{float(value):.17g}" for value in vector) + " }"


def generate_seeder_lua(grid: PipeGrid, direction: Iterable[float]) -> str:
    axis, u, v = pipe_frame(direction)
    cube_length = grid.dx_m * (2**ROOT_LEVEL)
    cube_origin = np.full(3, -0.5 * cube_length)
    inlet = -0.5 * PIPE_LENGTH_M * axis
    outlet = 0.5 * PIPE_LENGTH_M * axis
    halfwidth = 1.10 * PIPE_RADIUS_M
    inlet_origin = inlet - halfwidth * u - halfwidth * v
    outlet_origin = outlet - halfwidth * u - halfwidth * v
    inlet_u = 2.0 * halfwidth * v
    inlet_v = 2.0 * halfwidth * u  # v x u = -axis (outward at inlet)
    outlet_u = 2.0 * halfwidth * u
    outlet_v = 2.0 * halfwidth * v  # u x v = +axis (outward at outlet)
    return f"""-- Isolated pressure-BC pipe benchmark; not production CFD.
folder = 'mesh/'
comment = 'fixed physical Hagen-Poiseuille pipe, N={grid.cells_across_diameter}'
debug = {{ debugMode = false, debugFiles = false, debugMesh = 'debug/' }}
minlevel = {ROOT_LEVEL}
bounding_cube = {{ origin = {_lua_vector(cube_origin)}, length = {cube_length:.17g} }}
spatial_object = {{
  {{ attribute = {{ kind = 'boundary', label = 'wall', level = minlevel, calc_dist = true }}, geometry = {{ kind = 'stl', object = {{ filename = 'geometry/wall.stl' }} }} }},
  {{ attribute = {{ kind = 'boundary', label = 'inlet', level = minlevel }}, geometry = {{ kind = 'canoND', object = {{ origin = {_lua_vector(inlet_origin)}, vec = {{ {_lua_vector(inlet_u)}, {_lua_vector(inlet_v)} }} }} }} }},
  {{ attribute = {{ kind = 'boundary', label = 'outlet', level = minlevel }}, geometry = {{ kind = 'canoND', object = {{ origin = {_lua_vector(outlet_origin)}, vec = {{ {_lua_vector(outlet_u)}, {_lua_vector(outlet_v)} }} }} }} }},
  {{ attribute = {{ kind = 'seed' }}, geometry = {{ kind = 'canoND', object = {{ origin = {{ 0, 0, 0 }} }} }} }}
}}
"""


def generate_musubi_lua(grid: PipeGrid, candidate: str) -> str:
    if candidate not in CANDIDATES:
        raise ValueError(f"Unknown source-proven candidate: {candidate}")
    initial_pressure = PRESSURE_REFERENCE_PA + 0.5 * hagen_poiseuille_delta_p()
    return f"""-- Isolated fixed-geometry pressure-BC benchmark; not production CFD.
simulation_name = 'pipe_{candidate}_{grid.cells_across_diameter}'
printRuntimeInfo = true
timing_file = 'timing.res'
mesh = 'mesh/'
scaling = 'diffusive'
logging = {{ level = 5 }}
maximum_iterations = {MAXIMUM_ITERATIONS}
dx = {grid.dx_m:.17g}
dt = {grid.dt_s:.17g}
rho0_phy = {RHO_KG_M3:.17g}
nu_phy = {NU_M2_S:.17g}
bulk_viscosity_phy = {BULK_NU_M2_S:.17g}
pressure_reference_phy = {PRESSURE_REFERENCE_PA:.17g}
function outlet_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA:.17g} end
physics = {{ dt = dt, rho0 = rho0_phy }}
identify = {{ label = 'pressure_bc_pipe', kind = 'fluid', layout = 'd3q19', relaxation = 'bgk' }}
fluid = {{ kinematic_viscosity = nu_phy, bulk_viscosity = bulk_viscosity_phy }}
initial_condition = {{ pressure = {initial_pressure:.17g}, velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0 }}
boundary_condition = {{
  {{ label = 'wall', kind = 'wall_libb' }},
  {{ label = 'inlet', kind = 'adaptive_flux_pressure', mass_flowrate = {TARGET_MASS_FLOW_KG_S:.17g} }},
  {{ label = 'outlet', kind = '{CANDIDATES[candidate]}', pressure = outlet_pressure }}
}}
sim_control = {{
  time_control = {{ max = {{ iter = maximum_iterations }}, interval = {{ iter = 100 }} }},
  abort_criteria = {{ stop_file = 'stop' }}
}}
"""


def parse_controller_output(text: str) -> list[dict[str, float | int]]:
    rows: list[dict[str, float | int]] = []
    for match in CONTROLLER_PATTERN.finditer(text):
        rows.append(
            {
                "iteration": int(match.group("iteration")),
                "target_lattice": float(match.group("target")),
                "controlled_lattice": float(match.group("controlled")),
                "relative_error": float(match.group("relative")),
                "rho_boundary": float(match.group("rho")),
                "pressure_pa": float(match.group("pressure")),
                "max_lattice_velocity": float(match.group("speed")),
                "minimum_pdf": float(match.group("minimum_pdf")),
                "globBC_count": int(match.group("count")),
            }
        )
    return rows


def observed_order(dx: Sequence[float], error: Sequence[float]) -> dict[str, Any]:
    spacing = np.asarray(dx, dtype=np.float64)
    values = np.asarray(error, dtype=np.float64)
    valid = np.isfinite(spacing) & np.isfinite(values) & (spacing > 0.0) & (values > 0.0)
    if np.count_nonzero(valid) < 3:
        return {"available": False, "value": None, "reason": "fewer than three positive finite errors"}
    slope, intercept = np.polyfit(np.log(spacing[valid]), np.log(values[valid]), 1)
    return {
        "available": True,
        "value": float(slope),
        "fit_intercept": float(intercept),
        "points": int(np.count_nonzero(valid)),
        "method": "least-squares log(error) versus log(dx)",
    }


def invalid_gci_guard(values: Sequence[float], spacings: Sequence[float]) -> dict[str, Any]:
    """Reject GCI for non-monotone/oscillatory or non-positive-order data."""

    y = np.asarray(values, dtype=np.float64)
    h = np.asarray(spacings, dtype=np.float64)
    if len(y) < 3 or len(y) != len(h) or not np.all(np.isfinite(y)):
        return {"available": False, "reason": "insufficient finite values"}
    differences = np.diff(y)
    monotone = bool(np.all(differences > 0) or np.all(differences < 0))
    if not monotone:
        return {"available": False, "reason": "oscillatory/non-monotone sequence"}
    order = observed_order(h, np.abs(y - y[-1]) + np.finfo(float).eps)
    if not order["available"] or float(order["value"]) <= 0.0:
        return {"available": False, "reason": "positive observed order not established"}
    return {"available": True, "reason": None}


def compare_candidates(cases: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    grouped: dict[str, dict[str, list[Mapping[str, Any]]]] = {
        name: {} for name in CANDIDATES
    }
    for candidate in CANDIDATES:
        for orientation in ("axis_aligned", "worst_real_outlet"):
            grouped[candidate][orientation] = sorted(
                (row for row in cases if row["candidate"] == candidate and row["orientation"] == orientation),
                key=lambda row: row["dx_m"],
                reverse=True,
            )
    summaries: dict[str, Any] = {}
    for candidate, orientations in grouped.items():
        per_orientation: dict[str, Any] = {}
        fine_errors = []
        base_fine = []
        safe = True
        monotone = True
        for orientation, rows in orientations.items():
            by_resolution = {
                int(row["cells_across_diameter"]): row for row in rows
            }
            required = {12, 16, 20, 27}
            missing = sorted(required - set(by_resolution))
            if missing:
                raise ValueError(
                    f"Missing pressure benchmark resolutions for {candidate}/"
                    f"{orientation}: {missing}"
                )
            errors = [float(row["relative_delta_p_error"]) for row in rows]
            dx = [float(row["dx_m"]) for row in rows]
            base_delta = float(by_resolution[16]["delta_p_numerical_pa"])
            fine_delta = float(by_resolution[20]["delta_p_numerical_pa"])
            base_to_fine = abs(fine_delta - base_delta) / abs(base_delta)
            per_orientation[orientation] = {
                "resolution_semantics": {
                    "coarse_like": 12,
                    "base_like": 16,
                    "fine_like": 20,
                    "extra_fine_diagnostic": 27,
                },
                "analytic_relative_errors": errors,
                "observed_order": observed_order(dx, errors),
                "monotone_error_decrease": bool(np.all(np.diff(errors) <= 0.0)),
                "base_to_fine_delta_p_relative_difference": base_to_fine,
                "fine_analytic_relative_error": float(
                    by_resolution[20]["relative_delta_p_error"]
                ),
                "extra_fine_diagnostic": {
                    "cells_across_diameter": 27,
                    "delta_p_numerical_pa": float(
                        by_resolution[27]["delta_p_numerical_pa"]
                    ),
                    "analytic_relative_error": float(
                        by_resolution[27]["relative_delta_p_error"]
                    ),
                    "excluded_from_base_to_fine": True,
                },
            }
            fine_errors.append(per_orientation[orientation]["fine_analytic_relative_error"])
            base_fine.append(per_orientation[orientation]["base_to_fine_delta_p_relative_difference"])
            monotone &= per_orientation[orientation]["monotone_error_decrease"]
            safe &= all(bool(row["numerical_safety_pass"]) for row in rows)
        fine_axis = next(
            row["delta_p_numerical_pa"]
            for row in orientations["axis_aligned"]
            if int(row["cells_across_diameter"]) == 20
        )
        fine_oblique = next(
            row["delta_p_numerical_pa"]
            for row in orientations["worst_real_outlet"]
            if int(row["cells_across_diameter"]) == 20
        )
        spread = abs(float(fine_axis) - float(fine_oblique)) / hagen_poiseuille_delta_p()
        candidate_pass = (
            max(fine_errors) <= 0.02
            and max(base_fine) <= 0.02
            and spread <= 0.02
            and safe
            and monotone
        )
        summaries[candidate] = {
            "orientations": per_orientation,
            "worst_fine_analytic_error": max(fine_errors),
            "worst_base_to_fine_difference": max(base_fine),
            "orientation_spread": spread,
            "numerical_safety_pass": safe,
            "monotone_error_decrease": monotone,
            "candidate_pass": candidate_pass,
        }
    current = summaries["pressure_eq"]
    alternatives = [name for name in CANDIDATES if name != "pressure_eq"]
    passing = [name for name in alternatives if summaries[name]["candidate_pass"]]
    winner = None
    strength = "INCONCLUSIVE"
    if current["candidate_pass"]:
        better = [
            name for name in passing
            if summaries[name]["worst_fine_analytic_error"]
            <= 0.5 * current["worst_fine_analytic_error"]
        ]
        if better:
            winner = min(better, key=lambda name: summaries[name]["worst_fine_analytic_error"])
            strength = "STRONG_ALTERNATIVE_AT_LEAST_2X_LOWER_WORST_FINE_ERROR"
        else:
            winner = "pressure_eq"
            strength = "CURRENT_PASSES_NO_ALTERNATIVE_MEETS_REPLACEMENT_THRESHOLD"
    elif passing:
        winner = min(passing, key=lambda name: summaries[name]["worst_fine_analytic_error"])
        strength = "STRONG_CURRENT_FAILS_ALTERNATIVE_PASSES_ALL_2_PERCENT_GATES"
    return {
        "candidates": summaries,
        "winner": winner,
        "winner_evidence_strength": strength,
        "benchmark_conclusive": winner is not None,
    }


def _windows_to_wsl(path: Path) -> str:
    resolved = Path(path).resolve()
    drive = resolved.drive.rstrip(":").lower()
    tail = resolved.as_posix().split(":", 1)[1]
    return f"/mnt/{drive}{tail}"


def _run_wsl(command: str, *, timeout: int = 1800) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["wsl", "bash", "-lc", command],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=timeout,
    )


def _first_failure(stderr: str, stdout: str) -> str | None:
    for source in (stderr, stdout):
        for line in source.splitlines():
            stripped = line.strip()
            lowered = stripped.lower()
            if "fatal error occured in loading scalar constants" in lowered:
                # Musubi emits this known warning while filling optional strain
                # initial conditions with zero; it is not the terminating error.
                continue
            if stripped and any(token in lowered for token in ("error", "abort", "failed", "invalid", "unknown boundary")):
                return stripped
    return None


def _mesh_case(
    benchmark_root: Path,
    *,
    orientation: str,
    direction: np.ndarray,
    grid: PipeGrid,
) -> tuple[Path, dict[str, Any], int]:
    case = benchmark_root / "meshes" / orientation / f"n{grid.cells_across_diameter:02d}"
    mesh = case / "mesh"
    summary_path = case / "mesh_summary.json"
    if summary_path.is_file() and mesh.is_dir():
        prior = json.loads(summary_path.read_text(encoding="utf-8"))
        if prior.get("status") == "PASS":
            return mesh, prior, 0
    geometry = case / "geometry"
    mesh.mkdir(parents=True, exist_ok=True)
    wall = generate_pipe_wall(geometry / "wall.stl", direction)
    (case / "seeder.lua").write_text(generate_seeder_lua(grid, direction), encoding="utf-8")
    started = time.perf_counter()
    process = _run_wsl(f"cd '{_windows_to_wsl(case)}' && '{SEEDER_WSL}' seeder.lua")
    elapsed = time.perf_counter() - started
    (case / "seeder_stdout.log").write_text(process.stdout, encoding="utf-8")
    (case / "seeder_stderr.log").write_text(process.stderr, encoding="utf-8")
    mesh_complete = all(
        (mesh / name).is_file()
        for name in ("header.lua", "bnd.lua", "bnd.lsb", "elemlist.lsb", "qval.lsb")
    )
    result = {
        "status": "PASS" if process.returncode == 0 and mesh_complete else "FAIL",
        "seeder_returncode": process.returncode,
        "seeder_wall_time_s": elapsed,
        "first_failure": (
            _first_failure(process.stderr, process.stdout)
            if process.returncode != 0 or not mesh_complete
            else None
        ),
        "orientation": orientation,
        "direction": direction.tolist(),
        "cells_across_diameter": grid.cells_across_diameter,
        "dx_m": grid.dx_m,
        "dt_s": grid.dt_s,
        "wall": wall,
    }
    if process.returncode == 0 and mesh_complete:
        contract = load_mesh_contract(mesh, allow_zero_normals=True, require_runtime_order=False)
        result.update(
            {
                "fluid_cell_count": len(contract.tree_ids),
                "boundary_cell_counts": {
                    label: len(boundary.cell_indices) for label, boundary in contract.boundaries.items()
                },
                "zero_normal_counts": {
                    label: int(np.count_nonzero(boundary.normal_indices < 0))
                    for label, boundary in contract.boundaries.items()
                },
            }
        )
    write_json(summary_path, result)
    if process.returncode != 0 or not mesh_complete:
        raise RuntimeError(f"Seeder failed for {orientation} N={grid.cells_across_diameter}: {result['first_failure']}")
    return mesh, result, 1


def _run_case(
    benchmark_root: Path,
    *,
    mesh: Path,
    orientation: str,
    direction: np.ndarray,
    grid: PipeGrid,
    candidate: str,
) -> tuple[dict[str, Any], int]:
    case = benchmark_root / "runs" / candidate / orientation / f"n{grid.cells_across_diameter:02d}"
    summary_path = case / "case_summary.json"
    if summary_path.is_file():
        prior = json.loads(summary_path.read_text(encoding="utf-8"))
        if prior.get("status") == "PASS":
            return prior, 0
    case.mkdir(parents=True, exist_ok=True)
    lua = generate_musubi_lua(grid, candidate)
    (case / "musubi.lua").write_text(lua, encoding="utf-8")
    runtime = f"{RUNTIME_ROOT_WSL}/{candidate}_{orientation}_n{grid.cells_across_diameter:02d}"
    source_mesh = _windows_to_wsl(mesh)
    source_lua = _windows_to_wsl(case / "musubi.lua")
    stage = _run_wsl(
        f"mkdir -p '{RUNTIME_ROOT_WSL}' && test ! -e '{runtime}' && mkdir '{runtime}' "
        f"&& cp -a '{source_mesh}' '{runtime}/mesh' && cp '{source_lua}' '{runtime}/musubi.lua'"
    )
    if stage.returncode != 0:
        raise RuntimeError(f"WSL benchmark stage failed: {_first_failure(stage.stderr, stage.stdout) or stage.stderr.strip()}")
    started = time.perf_counter()
    process = _run_wsl(
        f"cd '{runtime}' && '{MPIRUN_WSL}' --bind-to core --map-by core -np 1 "
        f"'{MUSUBI_WSL}' musubi.lua",
        timeout=3600,
    )
    elapsed = time.perf_counter() - started
    (case / "musubi_stdout.log").write_text(process.stdout, encoding="utf-8")
    (case / "musubi_stderr.log").write_text(process.stderr, encoding="utf-8")
    rows = parse_controller_output(process.stdout)
    analytic = hagen_poiseuille_delta_p()
    if process.returncode == 0 and len(rows) >= STEADY_WINDOW:
        window = rows[-STEADY_WINDOW:]
        pressures = np.asarray([float(row["pressure_pa"]) for row in window])
        numerical = float(np.median(pressures) - PRESSURE_REFERENCE_PA)
        pressure_span = float(np.ptp(pressures) / max(abs(float(np.mean(pressures))), 1.0))
        relative_error = abs(numerical - analytic) / analytic
        mass_residual = float(max(abs(float(row["relative_error"])) for row in window))
        minimum_pdf = float(min(float(row["minimum_pdf"]) for row in window))
        max_speed = float(max(float(row["max_lattice_velocity"]) for row in window))
        all_finite = bool(np.all(np.isfinite(pressures)))
        steady = pressure_span <= 5.0e-4
    else:
        numerical = math.nan
        pressure_span = math.inf
        relative_error = math.inf
        mass_residual = math.inf
        minimum_pdf = -math.inf
        max_speed = math.inf
        all_finite = False
        steady = False
    safety = process.returncode == 0 and minimum_pdf > 0.0 and max_speed < 0.05 and all_finite
    result = {
        "status": "PASS" if process.returncode == 0 and steady and safety else "FAIL",
        "candidate": candidate,
        "musubi_kind": CANDIDATES[candidate],
        "orientation": orientation,
        "direction": direction.tolist(),
        "cells_across_diameter": grid.cells_across_diameter,
        "dx_m": grid.dx_m,
        "dt_s": grid.dt_s,
        "musubi_returncode": process.returncode,
        "musubi_wall_time_s": elapsed,
        "controller_rows": len(rows),
        "last_iteration": int(rows[-1]["iteration"]) if rows else None,
        "delta_p_numerical_pa": numerical if math.isfinite(numerical) else None,
        "delta_p_analytic_pa": analytic,
        "relative_delta_p_error": relative_error if math.isfinite(relative_error) else None,
        "target_q_m3_s": TARGET_Q_M3_S,
        "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
        "controller_mass_target_relative_error": mass_residual if math.isfinite(mass_residual) else None,
        "steady_window_pressure_relative_span": pressure_span if math.isfinite(pressure_span) else None,
        "steady_window_pass": steady,
        "minimum_pdf": minimum_pdf if math.isfinite(minimum_pdf) else None,
        "maximum_lattice_speed": max_speed if math.isfinite(max_speed) else None,
        "all_finite": all_finite,
        "numerical_safety_pass": safety,
        "first_failure": _first_failure(process.stderr, process.stdout) if process.returncode != 0 else (None if steady and safety else "steady-window or numerical-safety gate failed"),
        "runtime_root_wsl": runtime,
        "q_value_used_by_pressure_bc": False,
    }
    write_json(summary_path, result)
    if result["status"] != "PASS":
        raise RuntimeError(f"Benchmark case failed: {candidate}/{orientation}/N{grid.cells_across_diameter}: {result['first_failure']}")
    return result, 1


def run_pressure_bc_benchmark(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / RESEARCH_RUN
    benchmark_root = run_root / "pressure_bc_benchmark"
    qc = run_root / "qc"
    benchmark_root.mkdir(parents=True, exist_ok=True)
    qc.mkdir(parents=True, exist_ok=True)
    binary_unc = Path("//wsl.localhost/Ubuntu" + MUSUBI_WSL)
    if not binary_unc.is_file() or sha256_file(binary_unc) != MUSUBI_SHA256:
        raise RuntimeError("Pinned adaptive Musubi binary is missing or changed")
    worst = worst_real_outlet_direction(root)
    orientations = {
        "axis_aligned": np.asarray((0.0, 0.0, 1.0)),
        "worst_real_outlet": np.asarray(worst["direction"], dtype=np.float64),
    }
    grids = benchmark_grids()
    seeder_calls = 0
    musubi_calls = 0
    mesh_records: dict[str, Any] = {}
    cases: list[dict[str, Any]] = []
    try:
        meshes: dict[tuple[str, int], Path] = {}
        for orientation, direction in orientations.items():
            mesh_records[orientation] = {}
            for grid in grids:
                mesh, mesh_record, calls = _mesh_case(
                    benchmark_root, orientation=orientation, direction=direction, grid=grid
                )
                seeder_calls += calls
                meshes[(orientation, grid.cells_across_diameter)] = mesh
                mesh_records[orientation][str(grid.cells_across_diameter)] = mesh_record
        for candidate in CANDIDATES:
            for orientation, direction in orientations.items():
                for grid in grids:
                    case, calls = _run_case(
                        benchmark_root,
                        mesh=meshes[(orientation, grid.cells_across_diameter)],
                        orientation=orientation,
                        direction=direction,
                        grid=grid,
                        candidate=candidate,
                    )
                    musubi_calls += calls
                    cases.append(case)
        comparison = compare_candidates(cases)
        status = "PASS" if comparison["benchmark_conclusive"] else "CFD_FLOW_PRESSURE_BC_BENCHMARK_INCONCLUSIVE"
        first_failure = None
    except Exception as error:
        comparison = None
        status = "FAIL"
        first_failure = str(error)
    result = {
        "status": status,
        "first_failure": first_failure,
        "fixed_physical_geometry": {
            "diameter_m": PIPE_DIAMETER_M,
            "radius_m": PIPE_RADIUS_M,
            "length_m": PIPE_LENGTH_M,
            "area_m2": PIPE_AREA_M2,
            "analytic_delta_p_pa": hagen_poiseuille_delta_p(),
        },
        "physics": {
            "rho_kg_m3": RHO_KG_M3,
            "nu_m2_s": NU_M2_S,
            "bulk_nu_m2_s": BULK_NU_M2_S,
            "mean_velocity_m_s": MEAN_VELOCITY_M_S,
            "q_m3_s": TARGET_Q_M3_S,
            "mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
            "wall": "wall_libb with calc_dist=true",
            "inlet": "adaptive_flux_pressure",
        },
        "grids": [grid.__dict__ if hasattr(grid, "__dict__") else {"cells_across_diameter": grid.cells_across_diameter, "dx_m": grid.dx_m, "dt_s": grid.dt_s} for grid in grids],
        "orientations": {
            key: value.tolist() for key, value in orientations.items()
        },
        "worst_real_outlet_selection": worst,
        "candidates": list(CANDIDATES),
        "mesh_records": mesh_records,
        "case_results": cases,
        "comparison": comparison,
        "seeder_calls_this_invocation": seeder_calls,
        "musubi_calls_this_invocation": musubi_calls,
        "harvester_calls": 0,
        "vascular_cfd_calls": 0,
        "extra_fine_vascular_cfd_calls": 0,
        "pinned_binary_wsl": MUSUBI_WSL,
        "pinned_binary_sha256": MUSUBI_SHA256,
        "historical_failed_candidate_preflight": {
            "calls": 1,
            "candidate": "pressure_antiBounceBack",
            "incorrect_token": "pressure_antiBounceBack",
            "source_proven_correction": "pressure_antibounceback",
            "failure": "Unknown Boundary condition kind",
            "failed_configuration_evidence": str(
                benchmark_root / "failure_evidence" /
                "anti_camelcase_preflight" / "musubi.lua"
            ),
            "failed_stdout_preserved": False,
            "note": "The same Windows case log path was overwritten by the corrected-token run; the rejected configuration itself was recovered from its isolated WSL runtime.",
        },
    }
    write_json(qc / "pressure_bc_benchmark.json", result)
    return result


def refresh_pressure_benchmark_derived_summary(project_root: Path) -> dict[str, Any]:
    """Correct derived N16->N20 semantics without touching raw case results."""

    root = Path(project_root).resolve()
    path = (
        root
        / "outputs"
        / "cfd_flow"
        / RESEARCH_RUN
        / "qc"
        / "pressure_bc_benchmark.json"
    )
    result = json.loads(path.read_text(encoding="utf-8"))
    original_cases = result["case_results"]
    result["comparison"] = compare_candidates(original_cases)
    result["resolution_semantics"] = {
        "N12": "COARSE_LIKE",
        "N16": "BASE_LIKE",
        "N20": "FINE_LIKE",
        "N27": "EXTRA_FINE_DIAGNOSTIC",
        "base_to_fine_definition": "N16_TO_N20",
        "raw_case_results_modified": False,
    }
    result["derived_summary_revision"] = "PRESSURE_BENCHMARK_RESOLUTION_SEMANTICS_V2"
    write_json(path, result)
    return result


def finalize_research_decision(project_root: Path) -> dict[str, Any]:
    """Combine the three independent evidence layers without launching CFD."""

    root = Path(project_root).resolve()
    qc = root / "outputs" / "cfd_flow" / RESEARCH_RUN / "qc"
    forensic = json.loads((qc / "port_grid_sensitivity_forensics.json").read_text(encoding="utf-8"))
    planes = json.loads((qc / "standardized_outlet_plane_contract.json").read_text(encoding="utf-8"))
    benchmark = json.loads((qc / "pressure_bc_benchmark.json").read_text(encoding="utf-8"))
    production = (
        root / "cfd_flow.py",
        root / "configs" / "cfd_flow.yaml",
        root / "utils" / "cfd_flow" / "pipeline.py",
    )
    production_hashes = {str(path): sha256_file(path) for path in production}
    comparison = benchmark.get("comparison") or {}
    winner = comparison.get("winner")
    conclusive = bool(comparison.get("benchmark_conclusive"))
    planes_pass = planes.get("status") == "PASS"
    if conclusive and planes_pass:
        remediation = {
            "pressure_boundary": winner,
            "outlets": "STANDARDIZED_MATHEMATICAL_PLANES",
        }
        status = "REMEDIATION_ELIGIBLE_BUT_NOT_EXECUTED"
        root_cause = "GEOMETRY_AND_BC"
        next_step = "RUN ONE ISOLATED VASCULAR COARSE_BASE_FINE REMEDIATION VALIDATION"
    else:
        remediation = None
        status = "CFD_FLOW_GRID_SENSITIVITY_ROOT_CAUSE_UNRESOLVED"
        root_cause = "UNRESOLVED"
        next_step = (
            "ISOLATE WALL_AND_BIFURCATION_GEOMETRY_ERROR; DO NOT RUN VASCULAR EXTRA_FINE "
            "OR SWITCH PRESSURE_BC WITHOUT A CONCLUSIVE BENCHMARK"
        )
    result = {
        "status": status,
        "next": next_step,
        "production_pipeline_modified": False,
        "production_file_sha256": production_hashes,
        "extra_fine_vascular_cfd_calls": 0,
        "vascular_remediation_cfd_executed": False,
        "new_vascular_grid_results": None,
        "root_cause_before_benchmark": forensic["comparison"]["root_cause_classification_before_benchmark"],
        "root_cause_final": root_cause,
        "mathematical_outlet_plane_qc": planes["status"],
        "benchmark_status": benchmark["status"],
        "benchmark_winner": winner,
        "winner_evidence_strength": comparison.get("winner_evidence_strength", "INCONCLUSIVE"),
        "selected_remediation_config": remediation,
        "selected_pressure_bc_uses_q_value": None if winner is None else False,
        "vmtk_used_at_runtime": False,
        "flow_extension_used": False,
        "first_failure": (
            "No pressure BC candidate met the <=2% fine analytic-error, <=2% Base-to-Fine, "
            "and <=2% orientation-spread gates; no unique evidence-backed remediation exists."
            if winner is None
            else None
        ),
        "vascular_calls": {"seeder": 0, "musubi": 0, "harvester": 0},
        "source_frozen_files_modified": False,
    }
    write_json(qc / "root_cause_decision.json", result)
    return result
