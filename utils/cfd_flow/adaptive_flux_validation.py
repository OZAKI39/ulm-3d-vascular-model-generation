"""Isolated canonical and vascular smoke validation for adaptive flux Musubi."""

from __future__ import annotations

import math
import re
import shlex
import subprocess
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from .adaptive_flux_pressure_audit import (
    LBPM_COMMIT,
    MUSUBI_COMMIT,
    MUS_COMMIT,
    TREELM_COMMIT,
    physical_mass_factor,
)
from .exact_link_flux import (
    EXPECTED_DX_M,
    FROZEN_SEEDER_RUN,
    FROZEN_STEADY_RUN,
    REFERENCE_DENSITY_KG_M3,
    TARGET_MASS_FLOW_KG_S,
    TARGET_Q_M3_S,
    _file_manifest,
)
from .io import FlowError, sha256_file, write_json
from .mcclure_adaptive_flux_reference import physical_volume_flux_to_lattice


WORKTREE_WSL = "/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300"
BINARY_WSL = f"{WORKTREE_WSL}/build/musubi_adaptive_flux"
EXPECTED_BINARY_SHA256 = "e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588"
SEEDER_WSL = "/home/lzy/.local/bin/seeder"
MPIRUN_WSL = "/home/lzy/.local/bin/mpirun"
EXPECTED_DT_S = 2.441406249999999e-08
PRESSURE_REFERENCE_PA = 23622.32012800001
NU_M2_S = 3.27e-6
BULK_NU_M2_S = 2.18e-6
CANONICAL_TIMESTEPS = 20
VASCULAR_TIMESTEPS = 30
AXIS_MESH_RUN = "axis_aligned_ideal_plane_inlet_preflight_anchor003274_20260829_120444"
VALIDATION_PREFIX = "mcclure_adaptive_flux_validation_anchor003274"

CANONICAL_FAILED = "CFD_FLOW_MCCLURE_ADAPTIVE_FLUX_CANONICAL_FAILED"
SMOKE_PASS = "CFD_FLOW_MCCLURE_ADAPTIVE_FLUX_SMOKE_PASS"
NEXT_FIX_CANONICAL = "FIX ISOLATED ADAPTIVE FLUX CANONICAL TEST"
NEXT_STEADY = "RUN ONE NEW ADAPTIVE-FLUX MUSUBI STEADY BASELINE"

CONTROLLER_PATTERN = re.compile(
    r"ADAPTIVE_FLUX_PRESSURE\s+iter=\s*(?P<iteration>\d+)"
    r"\s+target_lattice=\s*(?P<target>[+\-0-9.Ee]+)"
    r"\s+controlled_lattice=\s*(?P<controlled>[+\-0-9.Ee]+)"
    r"\s+relative_error=\s*(?P<error>[+\-0-9.Ee]+)"
    r"\s+rho_boundary=\s*(?P<rho>[+\-0-9.Ee]+)"
    r"\s+pressure_pa=\s*(?P<pressure>[+\-0-9.Ee]+)"
    r"\s+max_lattice_velocity=\s*(?P<speed>[+\-0-9.Ee]+)"
    r"\s+minimum_pdf=\s*(?P<minimum_pdf>[+\-0-9.Ee]+)"
    r"\s+globBC_count=\s*(?P<count>\d+)"
)


@dataclass(frozen=True, slots=True)
class ExternalResult:
    returncode: int
    command: tuple[str, ...]
    stdout: Path
    stderr: Path


def _wsl_path(path: Path) -> str:
    resolved = Path(path).resolve()
    drive = resolved.drive.rstrip(":").lower()
    tail = resolved.as_posix().split(":", 1)[1]
    return f"/mnt/{drive}{tail}"


def _run_wsl(
    arguments: list[str], *, cwd: Path, stdout: Path, stderr: Path
) -> ExternalResult:
    shell_command = f"cd {shlex.quote(_wsl_path(cwd))} && {shlex.join(arguments)}"
    completed = subprocess.run(
        ["wsl.exe", "-e", "bash", "-lc", shell_command],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    stdout.write_text(completed.stdout, encoding="utf-8")
    stderr.write_text(completed.stderr, encoding="utf-8")
    return ExternalResult(
        returncode=completed.returncode,
        command=tuple(arguments),
        stdout=stdout,
        stderr=stderr,
    )


def _controller_rows(stdout: Path, stderr: Path) -> list[dict[str, float | int]]:
    text = stdout.read_text(encoding="utf-8", errors="replace") + "\n" + stderr.read_text(
        encoding="utf-8", errors="replace"
    )
    rows: list[dict[str, float | int]] = []
    for match in CONTROLLER_PATTERN.finditer(text):
        rows.append(
            {
                "iteration": int(match.group("iteration")),
                "target_lattice": float(match.group("target")),
                "controlled_lattice": float(match.group("controlled")),
                "relative_error": float(match.group("error")),
                "rho_boundary": float(match.group("rho")),
                "pressure_pa": float(match.group("pressure")),
                "max_lattice_velocity": float(match.group("speed")),
                "minimum_pdf": float(match.group("minimum_pdf")),
                "globBC_count": int(match.group("count")),
            }
        )
    return rows


def _summarize_controller(
    rows: list[dict[str, float | int]], *, maximum_error: float, expected_count: int | None
) -> dict[str, Any]:
    if not rows:
        raise FlowError(CANONICAL_FAILED, "No adaptive controller records were emitted")
    numeric = np.asarray(
        [
            [
                row["target_lattice"],
                row["controlled_lattice"],
                row["relative_error"],
                row["rho_boundary"],
                row["pressure_pa"],
                row["max_lattice_velocity"],
                row["minimum_pdf"],
            ]
            for row in rows
        ],
        dtype=np.float64,
    )
    counts = sorted({int(row["globBC_count"]) for row in rows})
    if not np.all(np.isfinite(numeric)):
        raise FlowError(CANONICAL_FAILED, "Adaptive controller emitted NaN/Inf")
    if expected_count is not None and counts != [expected_count]:
        raise FlowError(
            CANONICAL_FAILED, f"Adaptive inlet globBC changed: {counts}, expected {expected_count}"
        )
    result = {
        "status": "PASS",
        "records": len(rows),
        "iterations": [int(row["iteration"]) for row in rows],
        "globbc_counts": counts,
        "rho_range": [float(np.min(numeric[:, 3])), float(np.max(numeric[:, 3]))],
        "pressure_range_pa": [
            float(np.min(numeric[:, 4])),
            float(np.max(numeric[:, 4])),
        ],
        "controlled_lattice_flux_range": [
            float(np.min(numeric[:, 1])),
            float(np.max(numeric[:, 1])),
        ],
        "maximum_relative_error": float(np.max(numeric[:, 2])),
        "final_relative_error": float(numeric[-1, 2]),
        "maximum_lattice_speed": float(np.max(numeric[:, 5])),
        "minimum_pdf": float(np.min(numeric[:, 6])),
        "all_finite": True,
        "error_gate": maximum_error,
    }
    if result["maximum_relative_error"] > maximum_error:
        raise FlowError(
            CANONICAL_FAILED,
            f"controlled flux error {result['maximum_relative_error']} > {maximum_error}",
        )
    if result["maximum_lattice_speed"] >= 0.05:
        raise FlowError(CANONICAL_FAILED, "maximum lattice speed is not below 0.05")
    return result


def _canonical_seeder_lua() -> str:
    return """-- Isolated tiny Cartesian z-channel for adaptive flux validation.
folder = 'mesh/'
logging = { level = 4 }
debug = { debugMode = false, debugFiles = false }
dx = 2.0e-7
width = 8*dx
length = 16*dx
level = 5
cube = (2^level)*dx
eps = cube/(2^22)
bounding_cube = { origin = {-dx, -dx, -dx}, length = cube }
minlevel = level
spatial_object = {
  { attribute = {kind='seed'}, geometry = {kind='canoND', object={origin={width/2,width/2,length/2}}}},
  { attribute = {kind='boundary',label='inlet'}, geometry={kind='canoND',object={origin={-eps,-eps,-eps},vec={{width+2*eps,0,0},{0,width+2*eps,0}}}}},
  { attribute = {kind='boundary',label='outlet'}, geometry={kind='canoND',object={origin={-eps,-eps,length+eps},vec={{width+2*eps,0,0},{0,width+2*eps,0}}}}},
  { attribute = {kind='boundary',label='wall_x0'}, geometry={kind='canoND',object={origin={-eps,-eps,-eps},vec={{0,width+2*eps,0},{0,0,length+2*eps}}}}},
  { attribute = {kind='boundary',label='wall_x1'}, geometry={kind='canoND',object={origin={width+eps,-eps,-eps},vec={{0,width+2*eps,0},{0,0,length+2*eps}}}}},
  { attribute = {kind='boundary',label='wall_y0'}, geometry={kind='canoND',object={origin={-eps,-eps,-eps},vec={{width+2*eps,0,0},{0,0,length+2*eps}}}}},
  { attribute = {kind='boundary',label='wall_y1'}, geometry={kind='canoND',object={origin={-eps,width+eps,-eps},vec={{width+2*eps,0,0},{0,0,length+2*eps}}}}}
}
"""


def _common_musubi_lua(*, mesh: str, maximum_iterations: int, boundaries: str) -> str:
    return f"""-- Isolated adaptive flux validation; no restart and no field export.
simulation_name = 'adaptive_flux_validation'
printRuntimeInfo = true
timing_file = 'musubi_timing.res'
mesh = '{mesh}'
scaling = 'diffusive'
logging = {{ level = 5 }}
dx = {EXPECTED_DX_M:.17g}
dt = {EXPECTED_DT_S:.17g}
rho0_phy = {REFERENCE_DENSITY_KG_M3:.17g}
nu_phy = {NU_M2_S:.17g}
bulk_viscosity_phy = {BULK_NU_M2_S:.17g}
pressure_reference_phy = {PRESSURE_REFERENCE_PA:.17g}
maximum_iterations = {maximum_iterations}
sim_control = {{
  time_control = {{ max = {{iter=maximum_iterations}}, interval={{iter=1}} }},
  abort_criteria = {{ stop_file = 'stop' }}
}}
physics = {{ dt=dt, rho0=rho0_phy }}
identify = {{ label='adaptive_flux_test', kind='fluid', layout='d3q19', relaxation='bgk' }}
fluid = {{ kinematic_viscosity=nu_phy, bulk_viscosity=bulk_viscosity_phy }}
initial_condition = {{ pressure=pressure_reference_phy, velocityX=0.0, velocityY=0.0, velocityZ=0.0 }}
{boundaries}
"""


def _canonical_musubi_lua() -> str:
    boundaries = f"""boundary_condition = {{
  {{label='inlet', kind='adaptive_flux_pressure', mass_flowrate={TARGET_MASS_FLOW_KG_S:.17g}}},
  {{label='outlet', kind='pressure_eq', pressure=pressure_reference_phy}},
  {{label='wall_x0', kind='wall'}}, {{label='wall_x1', kind='wall'}},
  {{label='wall_y0', kind='wall'}}, {{label='wall_y1', kind='wall'}}
}}"""
    return _common_musubi_lua(
        mesh="./mesh/", maximum_iterations=CANONICAL_TIMESTEPS, boundaries=boundaries
    )


def _vascular_musubi_lua(mesh_wsl: str) -> str:
    boundaries = f"""function outlet_01_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA + 14.544978101274268:.17g} end
function outlet_02_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA + 132.20454922317552:.17g} end
function outlet_03_pressure(x,y,z,t) return {PRESSURE_REFERENCE_PA - 13.700626673311461:.17g} end
boundary_condition = {{
  {{label='wall', kind='wall_libb'}},
  {{label='inlet', kind='adaptive_flux_pressure', mass_flowrate={TARGET_MASS_FLOW_KG_S:.17g}}},
  {{label='outlet_01', kind='pressure_eq', pressure=outlet_01_pressure}},
  {{label='outlet_02', kind='pressure_eq', pressure=outlet_02_pressure}},
  {{label='outlet_03', kind='pressure_eq', pressure=outlet_03_pressure}}
}}"""
    return _common_musubi_lua(
        mesh=mesh_wsl, maximum_iterations=VASCULAR_TIMESTEPS, boundaries=boundaries
    )


def _range(values: list[float]) -> list[float]:
    return [float(min(values)), float(max(values))]


def run_adaptive_flux_validation(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    head = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    ).stdout.strip()
    binary_windows = Path(
        r"\\wsl.localhost\Ubuntu\home\lzy\apes-worktrees\musubi_mcclure_adaptive_flux_20260829_1300\build\musubi_adaptive_flux"
    )
    if sha256_file(binary_windows) != EXPECTED_BINARY_SHA256:
        raise FlowError(CANONICAL_FAILED, "Adaptive Musubi binary hash changed")
    output_root = root / "outputs" / "cfd_flow"
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_root = output_root / f"{VALIDATION_PREFIX}_{stamp}"
    canonical = run_root / "canonical"
    vascular = run_root / "vascular"
    qc = run_root / "qc"
    for directory in (canonical, vascular, qc):
        directory.mkdir(parents=True, exist_ok=False)

    source_paths = (
        root / "cfd_flow.py",
        root / "configs" / "cfd_flow.yaml",
        root / "utils" / "cfd_flow" / "pipeline.py",
        output_root / FROZEN_STEADY_RUN / "diagnostic_musubi.lua",
        output_root / FROZEN_SEEDER_RUN / "seeder" / "mesh" / "header.lua",
        output_root / AXIS_MESH_RUN / "seeder" / "mesh" / "header.lua",
    )
    source_before = _file_manifest(source_paths)
    summary: dict[str, Any] = {
        "status": CANONICAL_FAILED,
        "next": NEXT_FIX_CANONICAL,
        "run_root": str(run_root),
        "actual_head": head,
        "production_pipeline_modified": False,
        "pinned_musubi_source_modified": False,
        "canonical_seeder_calls": 0,
        "canonical_musubi_calls": 0,
        "vascular_seeder_calls": 0,
        "vascular_musubi_calls": 0,
        "harvester_calls": 0,
        "grid_convergence": "NOT_RUN",
    }
    write_json(qc / "adaptive_flux_validation_manifest.json", summary)

    try:
        patch = root / "patches" / "musubi" / "adaptive_flux_pressure.patch"
        provenance = {
            "status": "PASS",
            "source_base_commit": MUSUBI_COMMIT,
            "mus_submodule_base_commit": MUS_COMMIT,
            "treelm_submodule_commit": TREELM_COMMIT,
            "lbpm_reference_commit": LBPM_COMMIT,
            "worktree_path": WORKTREE_WSL,
            "patch_path": str(patch),
            "patch_sha256": sha256_file(patch),
            "compiler": "/home/lzy/.local/bin/mpif90 (GNU Fortran 13.4.0)",
            "mpi": "/home/lzy/.local/bin/mpirun",
            "build_command": "FC=/home/lzy/.local/bin/mpif90 ./bin/waf configure --notests; ./bin/waf build --notests --targets=musubi -j8",
            "binary_path": BINARY_WSL,
            "binary_sha256": EXPECTED_BINARY_SHA256,
            "pinned_executable_overwritten": False,
        }
        write_json(qc / "adaptive_flux_musubi_build_provenance.json", provenance)

        (canonical / "seeder.lua").write_text(_canonical_seeder_lua(), encoding="utf-8")
        (canonical / "musubi.lua").write_text(_canonical_musubi_lua(), encoding="utf-8")
        summary["canonical_seeder_calls"] = 1
        seed_result = _run_wsl(
            [SEEDER_WSL, "seeder.lua"],
            cwd=canonical,
            stdout=canonical / "seeder_stdout.log",
            stderr=canonical / "seeder_stderr.log",
        )
        if seed_result.returncode != 0:
            raise FlowError(CANONICAL_FAILED, "Canonical Seeder returned nonzero")
        summary["canonical_musubi_calls"] = 1
        canonical_result = _run_wsl(
            [MPIRUN_WSL, "-np", "1", BINARY_WSL, "musubi.lua"],
            cwd=canonical,
            stdout=canonical / "musubi_stdout.log",
            stderr=canonical / "musubi_stderr.log",
        )
        if canonical_result.returncode != 0:
            raise FlowError(CANONICAL_FAILED, "Canonical Musubi returned nonzero")
        canonical_rows = _controller_rows(
            canonical_result.stdout, canonical_result.stderr
        )
        canonical_qc = _summarize_controller(
            canonical_rows, maximum_error=1.0e-8, expected_count=None
        )
        canonical_qc.update(
            {
                "seeder_calls": 1,
                "musubi_calls": 1,
                "harvester_calls": 0,
                "seeder_returncode": seed_result.returncode,
                "musubi_returncode": canonical_result.returncode,
                "target_flux_relative_error_gate": 1.0e-8,
                "preferred_gate": 1.0e-10,
            }
        )
        write_json(qc / "canonical_adaptive_flux_qc.json", canonical_qc)

        axis_mesh = output_root / AXIS_MESH_RUN / "seeder" / "mesh"
        (vascular / "musubi.lua").write_text(
            _vascular_musubi_lua(_wsl_path(axis_mesh) + "/"), encoding="utf-8"
        )
        summary["vascular_musubi_calls"] = 1
        vascular_result = _run_wsl(
            [MPIRUN_WSL, "-np", "1", BINARY_WSL, "musubi.lua"],
            cwd=vascular,
            stdout=vascular / "musubi_stdout.log",
            stderr=vascular / "musubi_stderr.log",
        )
        if vascular_result.returncode != 0:
            raise FlowError(CANONICAL_FAILED, "Vascular smoke Musubi returned nonzero")
        vascular_rows = _controller_rows(vascular_result.stdout, vascular_result.stderr)
        vascular_qc = _summarize_controller(
            vascular_rows, maximum_error=1.0e-6, expected_count=287
        )
        if vascular_qc["minimum_pdf"] <= 0.0:
            raise FlowError(CANONICAL_FAILED, "Vascular smoke minimum PDF is not positive")
        factor = physical_mass_factor(
            density_kg_m3=REFERENCE_DENSITY_KG_M3,
            dx_m=EXPECTED_DX_M,
            dt_s=EXPECTED_DT_S,
        )
        masses = [float(row["controlled_lattice"]) * factor for row in vascular_rows]
        vascular_qc.update(
            {
                "seeder_calls": 0,
                "musubi_calls": 1,
                "harvester_calls": 0,
                "musubi_returncode": vascular_result.returncode,
                "fresh_initialization": True,
                "requested_timesteps": VASCULAR_TIMESTEPS,
                "target_q_physical_m3_s": TARGET_Q_M3_S,
                "target_q_lattice": physical_volume_flux_to_lattice(
                    TARGET_Q_M3_S, dx_m=EXPECTED_DX_M, dt_s=EXPECTED_DT_S
                ),
                "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
                "actual_controlled_mass_flow_range_kg_s": _range(masses),
                "final_controlled_mass_flow_kg_s": masses[-1],
                "minimum_pdf_positive": True,
                "area_proxy_used_by_controller": False,
            }
        )
        write_json(qc / "vascular_adaptive_flux_smoke_qc.json", vascular_qc)

        source_after = _file_manifest(source_paths)
        if source_before != source_after:
            raise FlowError(CANONICAL_FAILED, "Frozen source files changed during validation")
        write_json(
            qc / "source_frozen_files_unchanged_qc.json",
            {
                "status": "PASS",
                "source_frozen_files_unchanged": True,
                "before": source_before,
                "after": source_after,
            },
        )
        summary.update(
            {
                "status": SMOKE_PASS,
                "next": NEXT_STEADY,
                "adaptive_musubi_binary": BINARY_WSL,
                "adaptive_musubi_binary_sha256": EXPECTED_BINARY_SHA256,
                "canonical": canonical_qc,
                "vascular": vascular_qc,
                "vascular_mesh_path": str(axis_mesh),
                "source_frozen_files_unchanged": True,
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(qc / "adaptive_flux_validation_manifest.json", summary)
        return summary
    except Exception as error:
        summary.update(
            {
                "status": CANONICAL_FAILED,
                "next": NEXT_FIX_CANONICAL,
                "first_failure": str(error),
                "source_frozen_files_unchanged": source_before
                == _file_manifest(source_paths),
                "completed_at": datetime.now().isoformat(),
            }
        )
        write_json(qc / "adaptive_flux_validation_manifest.json", summary)
        return summary

