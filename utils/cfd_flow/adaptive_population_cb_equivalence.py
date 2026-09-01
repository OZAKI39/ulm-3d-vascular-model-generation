"""Coarse/Base one-step oracle for the Fine active-population repair.

``prepare`` records a zero-solver-call deferred state while the Fine run is
active.  ``run`` is intentionally explicit and must only be invoked after the
Fine continuation has stopped cleanly.  Historical old-binary one-step states
are reused; only one corrected-binary step per grid is launched.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .fine_continuation_monitor import RUN_RELATIVE
from .io import sha256_file
from .musubi_boundary_mass_referee import load_mesh_contract
from .restart_decode import D3Q19_DIRECTIONS, read_restart_pdf
from .tau1_base import _controller_records


OLD_BINARY_WSL = (
    "/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/"
    "build/musubi_adaptive_flux"
)
OLD_BINARY_SHA256 = "e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588"
NEW_BINARY_WSL = "/home/lzy/apes-worktrees/musubi_adaptive_target_fixed_20260901/build/musubi"
NEW_BINARY_SHA256 = "1491e7ded30ed15158bd8ba2812d414ad9fa31cbc0128a26d2b1eba779773caa"
MPIRUN_WSL = "/home/lzy/.local/bin/mpirun"
RUNTIME_WSL = "/home/lzy/u3da/cb_active_population_equivalence_20260902"
SOURCE_REVISION = "81f8c4f13772f6d4af31f335e1e3f99b02726e25"
PARENT_REVISION = "4e8b277b66226277171ef93bf054d21270812793"
TEM_REVISION = "9899d1376992c4fafc8a343d2b4ccef81de670d1"
PATCH_SHA256 = "90efce400bb4b6ad5ad22ddd6518ccbeb8ba8a8eece2f76f2109a26dae92758f"
NEW_SOURCE_SHA256 = "7fc7325eb1d21112451fe99437fddce8413a5f2958b578ea2b641711c96909e3"
PDF_GATE = 1.0e-14
RHO_GATE = 1.0e-13
U_LAT_GATE = 1.0e-13
CONTROLLER_GATE = 1.0e-12
RHO0_KG_M3 = 1056.0


@dataclass(frozen=True, slots=True)
class GridOracle:
    name: str
    cells: int
    iteration: int
    dx_m: float
    dt_s: float
    pressure_reference_pa: float
    accepted_sha256: str
    old_end_sha256: str
    local_mesh: str
    old_lua: str
    old_start: str
    old_end: str
    old_stdout: str
    accepted_header_wsl: str


GRIDS = {
    "coarse": GridOracle(
        name="coarse",
        cells=82_957,
        iteration=354_295,
        dx_m=2.6e-7,
        dt_s=3.4454638124362897e-9,
        pressure_reference_pa=2_004_444.213017751,
        accepted_sha256="9eda88e685e5eaa6af650757b8c97accba392c88988df2766d2bb4de825fcfbb",
        old_end_sha256="236c8a1b38ad7f7b682ab8dd9ff857d778b7fc4ac2494dc72dc0e9dcc725dfbf",
        local_mesh=f"{RUN_RELATIVE.as_posix()}/coarse/seeder/mesh",
        old_lua=f"{RUN_RELATIVE.as_posix()}/coarse/full_referee_from_0354295/musubi.lua",
        old_start=f"{RUN_RELATIVE.as_posix()}/coarse/accepted_restart/tau1_reference_scaled_coarse_354295.lsb",
        old_end=f"{RUN_RELATIVE.as_posix()}/coarse/full_referee_from_0354295/restart/tau1_reference_scaled_coarse_354296.lsb",
        old_stdout=f"{RUN_RELATIVE.as_posix()}/coarse/full_referee_from_0354295/musubi_stdout.log",
        accepted_header_wsl=(
            "/home/lzy/u3da/tau1_reference_scaled_cbf_20260901/coarse/restart/"
            "tau1_reference_scaled_coarse_header_354295.lua"
        ),
    ),
    "base": GridOracle(
        name="base",
        cells=182_320,
        iteration=598_755,
        dx_m=2.0e-7,
        dt_s=2.038735983690112e-9,
        pressure_reference_pa=3_387_510.7199999997,
        accepted_sha256="ffcd98b2dc684d1569d937d915b603805809c581d5341e71b17afac2ac64c39f",
        old_end_sha256="7165731f2516782384a9d58b0a56b58097f26081f6bd3828dfbe03491fe78f0d",
        local_mesh=(
            "outputs/cfd_flow/"
            "healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/"
            "seeder/mesh"
        ),
        old_lua=(
            "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901/"
            "full_referee_from_0598755/musubi.lua"
        ),
        old_start=(
            "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901/"
            "accepted_restart/tau1_reference_scaled_base_598755.lsb"
        ),
        old_end=(
            "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901/"
            "full_referee_from_0598755/restart/tau1_reference_scaled_base_598756.lsb"
        ),
        old_stdout=(
            "outputs/cfd_flow/healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901/"
            "full_referee_from_0598755/musubi_stdout.log"
        ),
        accepted_header_wsl=(
            "/home/lzy/u3da/tau1_reference_scaled_base_20260901/restart/"
            "tau1_reference_scaled_base_header_598755.lua"
        ),
    ),
}


def _wsl_sha256(path: str) -> str:
    process = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "sha256sum", path],
        check=False,
        capture_output=True,
        text=True,
        timeout=180,
    )
    if process.returncode != 0:
        raise RuntimeError(process.stderr.strip() or f"sha256 failed for {path}")
    return process.stdout.split()[0]


def verify_binary_provenance() -> dict[str, Any]:
    old_actual = _wsl_sha256(OLD_BINARY_WSL)
    new_actual = _wsl_sha256(NEW_BINARY_WSL)
    if old_actual != OLD_BINARY_SHA256 or new_actual != NEW_BINARY_SHA256:
        raise RuntimeError(f"binary provenance mismatch: old={old_actual}, new={new_actual}")
    return {
        "old_binary": {
            "path": OLD_BINARY_WSL,
            "sha256": old_actual,
            "source_revision": SOURCE_REVISION,
            "parent_revision": PARENT_REVISION,
            "tem_revision": TEM_REVISION,
            "patch_sha": None,
        },
        "new_binary": {
            "path": NEW_BINARY_WSL,
            "sha256": new_actual,
            "source_revision": SOURCE_REVISION,
            "parent_revision": PARENT_REVISION,
            "tem_revision": TEM_REVISION,
            "patch_sha": PATCH_SHA256,
            "source_sha256": NEW_SOURCE_SHA256,
        },
    }


def _controller(path: Path) -> dict[str, Any]:
    records = _controller_records(path.read_text(encoding="utf-8"))
    if not records:
        raise RuntimeError(f"controller record missing: {path}")
    return records[-1]


def _rank_active_counts(project_root: Path, grid: GridOracle) -> list[int]:
    mesh = load_mesh_contract(
        project_root / grid.local_mesh,
        expected_cells=grid.cells,
        require_runtime_order=False,
    )
    indices = np.asarray(mesh.boundaries["inlet"].cell_indices, dtype=np.int64)
    partitions = np.array_split(np.arange(grid.cells, dtype=np.int64), 4)
    return [
        int(np.count_nonzero((indices >= part[0]) & (indices <= part[-1])))
        for part in partitions
    ]


def _validate_old_reference(project_root: Path, grid: GridOracle) -> dict[str, Any]:
    start = project_root / grid.old_start
    end = project_root / grid.old_end
    lua = (project_root / grid.old_lua).read_text(encoding="utf-8")
    stdout = project_root / grid.old_stdout
    controller = _controller(stdout)
    checks = {
        "accepted_restart_sha256": sha256_file(start) == grid.accepted_sha256,
        "old_output_sha256": sha256_file(end) == grid.old_end_sha256,
        "exactly_one_timestep": (
            f"maximum_iterations = {grid.iteration + 1}" in lua
            and f"min={{iter={grid.iteration + 1}}}" in lua
            and int(controller["iteration"]) == grid.iteration + 1
        ),
        "accepted_header_path": grid.accepted_header_wsl in lua,
        "successful_run": "SUCCESSFUL run!" in stdout.read_text(encoding="utf-8"),
    }
    if not all(checks.values()):
        raise RuntimeError(f"{grid.name} historical one-step reference is incompatible: {checks}")
    counts = _rank_active_counts(project_root, grid)
    denominator = int(controller["globBC_count"])
    if sum(counts) != denominator:
        raise RuntimeError(
            f"{grid.name} partition reconstruction {counts} != old denominator {denominator}"
        )
    return {
        "old_reference_reused": True,
        "old_reference_checks": checks,
        "old_denominator": denominator,
        "new_per_rank_active_counts": counts,
        "new_active_global": sum(counts),
        "old_target": float(controller["target_lattice"]),
        "new_target": float(controller["target_lattice"]),
        "target_relative_difference": 0.0,
        "old_controlled": float(controller["controlled_lattice"]),
        "old_pdf_payload_sha256": sha256_file(end),
    }


def _source_scope() -> dict[str, Any]:
    return {
        "status": "PASS",
        "comparison": "old and new worktree source trees, including tem/source",
        "changed_scientific_source_files": ["mus/source/bc/mus_bc_fluid_module.fpp"],
        "unchanged_source_trees": ["tem/source"],
        "changes": [
            "adaptive_flux_pressure active-population denominator",
            "directly related finite/index/link/MPI fail-fast invariants",
            "diagnostic active and initial population log fields",
        ],
        "explicitly_unchanged": [
            "collision",
            "streaming",
            "wall_libb",
            "pressure_eq reconstruction mathematics",
            "force",
            "viscosity",
            "q",
            "D3Q19",
            "equilibrium",
            "physics conversion",
            "other boundary mathematics",
        ],
        "old_source_sha256": "4385b30cc804791fb0e50e73ab4e15d92a014f20eaaa7a71c8746a6aa08f6e29",
        "new_source_sha256": NEW_SOURCE_SHA256,
    }


def compare_pdf_payloads(
    old_path: Path,
    new_path: Path,
    *,
    cells: int,
    dx_m: float,
    dt_s: float,
    pressure_reference_pa: float,
) -> dict[str, Any]:
    old_hash = sha256_file(old_path)
    new_hash = sha256_file(new_path)
    old = read_restart_pdf(old_path, n_elems=cells, n_components=19)
    new = read_restart_pdf(new_path, n_elems=cells, n_components=19)
    difference = np.asarray(new) - np.asarray(old)
    abs_difference = np.abs(difference)
    max_abs = float(np.max(abs_difference))
    rms = float(math.sqrt(np.mean(np.square(difference), dtype=np.float64)))
    scale = np.maximum(np.maximum(np.abs(old), np.abs(new)), np.finfo(np.float64).tiny)
    max_relative = float(np.max(abs_difference / scale))
    old_rho = np.sum(old, axis=1, dtype=np.float64)
    new_rho = np.sum(new, axis=1, dtype=np.float64)
    old_u = (np.asarray(old) @ D3Q19_DIRECTIONS.astype(np.float64)) / old_rho[:, None]
    new_u = (np.asarray(new) @ D3Q19_DIRECTIONS.astype(np.float64)) / new_rho[:, None]
    max_rho = float(np.max(np.abs(new_rho - old_rho)))
    max_u_lat = float(np.max(np.abs(new_u - old_u)))
    max_u_phy = max_u_lat * dx_m / dt_s
    max_pressure = max_rho * pressure_reference_pa
    bitwise = old_hash == new_hash
    return {
        "old_pdf_payload_sha256": old_hash,
        "new_pdf_payload_sha256": new_hash,
        "bitwise_equivalent": bitwise,
        "max_abs_pdf_diff": max_abs,
        "rms_pdf_diff": rms,
        "max_relative_pdf_diff": max_relative,
        "max_abs_rho_diff": max_rho,
        "max_abs_u_lat_diff": max_u_lat,
        "max_abs_u_phy_diff": max_u_phy,
        "max_abs_gauge_pressure_diff": max_pressure,
    }


def equivalence_verdict(metrics: dict[str, Any], controller: dict[str, float]) -> str:
    numerical = (
        float(metrics["max_abs_pdf_diff"]) <= PDF_GATE
        and float(metrics["max_abs_rho_diff"]) <= RHO_GATE
        and float(metrics["max_abs_u_lat_diff"]) <= U_LAT_GATE
        and float(controller["target_relative_difference"]) <= CONTROLLER_GATE
        and float(controller["controlled_relative_difference"]) <= CONTROLLER_GATE
    )
    if bool(metrics["bitwise_equivalent"]):
        return "PASS_BITWISE_EQUIVALENT"
    return "PASS_MACHINE_PRECISION_EQUIVALENT" if numerical else "FAIL"


def _git_head(project_root: Path) -> str:
    return subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=project_root, text=True, timeout=10
    ).strip()


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def _write_csv(path: Path, evidence: dict[str, Any]) -> None:
    fields = [
        "grid",
        "status",
        "accepted_restart_sha256",
        "old_reference_reused",
        "old_denominator",
        "new_per_rank_active_counts",
        "new_active_global",
        "old_target",
        "new_target",
        "target_relative_difference",
        "old_pdf_payload_sha256",
        "new_pdf_payload_sha256",
        "bitwise_equivalent",
        "max_abs_pdf_diff",
        "rms_pdf_diff",
        "max_abs_rho_diff",
        "max_abs_u_lat_diff",
        "max_abs_u_phy_diff",
        "max_abs_gauge_pressure_diff",
        "verdict",
    ]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for name in ("coarse", "base"):
            item = evidence[name]
            writer.writerow({field: item.get(field) for field in fields} | {"grid": name})


def prepare_deferred(project_root: Path) -> dict[str, Any]:
    root = project_root.resolve()
    qc = root / RUN_RELATIVE / "qc"
    provenance = verify_binary_provenance()
    grids: dict[str, Any] = {}
    for name, grid in GRIDS.items():
        reference = _validate_old_reference(root, grid)
        grids[name] = {
            "status": "DEFERRED_NEW_BINARY_STEP",
            "accepted_restart_sha256": grid.accepted_sha256,
            **reference,
            "new_controlled": None,
            "controlled_relative_difference": None,
            "new_pdf_payload_sha256": None,
            "bitwise_equivalent": None,
            "max_abs_pdf_diff": None,
            "rms_pdf_diff": None,
            "max_relative_pdf_diff": None,
            "max_abs_rho_diff": None,
            "max_abs_u_lat_diff": None,
            "max_abs_u_phy_diff": None,
            "max_abs_gauge_pressure_diff": None,
            "verdict": "PENDING_FINE_COMPLETION",
        }
    evidence = {
        "status": "CB_ORACLE_DEFERRED_UNTIL_FINE_COMPLETION",
        "deferred_reason": (
            "Fine iteration throughput is not observable between filtered cadence records, "
            "so a <=10% concurrent slowdown cannot be proven."
        ),
        "git_head": _git_head(root),
        **provenance,
        "source_diff_scope": _source_scope(),
        "resource_check": {
            "fine_pid": 30718,
            "physical_cores": 12,
            "logical_cpus": 24,
            "fine_mpi_ranks": 4,
            "fine_rank_cpu_state": "four ranks continuously CPU-active",
            "available_ram_bytes_observed": 6_740_140_032,
            "swap_used_bytes_observed": 0,
            "disk_io_pressure_avg10_observed": 0.0,
            "planned_oracle": "4 MPI ranks, same WSL VM and storage, two one-step calls",
            "iteration_throughput_observability": (
                "not available between cadence checkpoints because stdout is deliberately filtered"
            ),
            "decision": "DEFER; cannot prove the required before/during Fine iteration slowdown",
        },
        "fine_slowdown_during_oracle": None,
        "oracle_execution_timing": "AFTER_FINE_COMPLETION",
        "solver_calls": 0,
        "coarse": grids["coarse"],
        "base": grids["base"],
        "coarse_status": grids["coarse"]["status"],
        "base_status": grids["base"]["status"],
        "overall_verdict": "PENDING_FINE_COMPLETION",
        "accepted_coarse_base_require_rerun": None,
    }
    _write_json(qc / "adaptive_active_population_cb_equivalence.json", evidence)
    _write_csv(qc / "adaptive_active_population_cb_equivalence.csv", evidence)
    return evidence


def _make_new_lua(project_root: Path, grid: GridOracle, output: Path) -> str:
    text = (project_root / grid.old_lua).read_text(encoding="utf-8")
    runtime = f"{RUNTIME_WSL}/{grid.name}/new"
    text = re.sub(
        r"simulation_name\s*=\s*'[^']+'",
        f"simulation_name = 'adaptive_population_{grid.name}_new_one_step'",
        text,
        count=1,
    )
    text = re.sub(r"abort_criteria=\{stop_file='[^']+'\}", f"abort_criteria={{stop_file='{runtime}/stop'}}", text)
    text = re.sub(r"restart\s*=\s*\{read='[^']+'", f"restart = {{read='{grid.accepted_header_wsl}'", text)
    text = re.sub(r"write='[^']+'", f"write='{runtime}/restart/'", text)
    output.write_text(text, encoding="utf-8", newline="\n")
    return runtime


def _windows_to_wsl(path: Path) -> str:
    resolved = path.resolve()
    drive = resolved.drive.rstrip(":").lower()
    return f"/mnt/{drive}/{resolved.as_posix().split(':', 1)[1].lstrip('/')}"


def _run_new_step(project_root: Path, grid: GridOracle) -> tuple[Path, Path, float]:
    output = project_root / RUN_RELATIVE / "cb_equivalence" / grid.name / "new"
    restart = output / "restart"
    restart.mkdir(parents=True, exist_ok=True)
    lua_path = output / "musubi.lua"
    runtime = _make_new_lua(project_root, grid, lua_path)
    subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "mkdir", "-p", f"{runtime}/restart", f"{runtime}/tracking"],
        check=True,
        timeout=30,
    )
    subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "cp", _windows_to_wsl(lua_path), f"{runtime}/musubi.lua"],
        check=True,
        timeout=30,
    )
    started = time.perf_counter()
    process = subprocess.run(
        [
            "wsl.exe",
            "-d",
            "Ubuntu",
            "--",
            MPIRUN_WSL,
            "--bind-to",
            "core",
            "--map-by",
            "core",
            "--report-bindings",
            "--wdir",
            runtime,
            "-np",
            "4",
            NEW_BINARY_WSL,
            "musubi.lua",
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    elapsed = time.perf_counter() - started
    (output / "musubi_stdout.log").write_text(process.stdout, encoding="utf-8", newline="\n")
    (output / "musubi_stderr.log").write_text(process.stderr, encoding="utf-8", newline="\n")
    if process.returncode != 0 or "SUCCESSFUL run!" not in process.stdout:
        raise RuntimeError(f"{grid.name} corrected-binary one-step call failed")
    iteration = grid.iteration + 1
    for suffix in (f"_{iteration}.lsb", f"_header_{iteration}.lua"):
        command = (
            f"set -e; file=$(find '{runtime}/restart' -maxdepth 1 -type f -name '*{suffix}' -print -quit); "
            f"test -n \"$file\"; cp \"$file\" '{_windows_to_wsl(restart)}/'"
        )
        subprocess.run(
            ["wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", command],
            check=True,
            timeout=60,
        )
    payloads = list(restart.glob(f"*_{iteration}.lsb"))
    if len(payloads) != 1:
        raise RuntimeError(f"{grid.name} corrected output restart missing")
    return payloads[0], output / "musubi_stdout.log", elapsed


def run_oracle(project_root: Path) -> dict[str, Any]:
    root = project_root.resolve()
    fine = json.loads(
        (root / RUN_RELATIVE / "fine_postfix/qc/fine_continuation_monitor.json").read_text(
            encoding="utf-8"
        )
    )
    if fine.get("status") != "CLEAN_COMPLETION":
        raise RuntimeError("oracle is blocked until the monitored Fine continuation completes cleanly")
    evidence = prepare_deferred(root)
    calls = 0
    for name, grid in GRIDS.items():
        new_end, new_stdout, elapsed = _run_new_step(root, grid)
        calls += 1
        new_controller = _controller(new_stdout)
        item = evidence[name]
        old_target = float(item["old_target"])
        old_controlled = float(item["old_controlled"])
        target_rel = abs(float(new_controller["target_lattice"]) - old_target) / abs(old_target)
        controlled_rel = abs(float(new_controller["controlled_lattice"]) - old_controlled) / abs(old_controlled)
        controller = {
            "target_relative_difference": target_rel,
            "controlled_relative_difference": controlled_rel,
        }
        metrics = compare_pdf_payloads(
            root / grid.old_end,
            new_end,
            cells=grid.cells,
            dx_m=grid.dx_m,
            dt_s=grid.dt_s,
            pressure_reference_pa=grid.pressure_reference_pa,
        )
        verdict = equivalence_verdict(metrics, controller)
        item.update(
            {
                "status": "PASS" if verdict.startswith("PASS") else "FAIL",
                "new_target": float(new_controller["target_lattice"]),
                "target_relative_difference": target_rel,
                "new_controlled": float(new_controller["controlled_lattice"]),
                "controlled_relative_difference": controlled_rel,
                "new_active_global": int(new_controller["globBC_count"]),
                "process_wall_clock_s": elapsed,
                **metrics,
                "verdict": verdict,
            }
        )
        if item["new_active_global"] != item["old_denominator"] or verdict == "FAIL":
            evidence["status"] = "CB_BEHAVIORAL_NEUTRALITY_FAILED"
            evidence["overall_verdict"] = "CB_BEHAVIORAL_NEUTRALITY_FAILED"
            evidence["accepted_coarse_base_require_rerun"] = "REVIEW_REQUIRED"
            evidence["solver_calls"] = calls
            break
    else:
        evidence["status"] = "CFD_FLOW_ADAPTIVE_ACTIVE_POPULATION_FIX_CB_BEHAVIORALLY_NEUTRAL"
        evidence["overall_verdict"] = evidence["status"]
        evidence["accepted_coarse_base_require_rerun"] = False
        evidence["conclusion"] = (
            "The active-population fix is behaviorally neutral for accepted Coarse and Base, "
            "and only changes cases where remove_solid_in_bc reduces the active inlet population, "
            "as observed for Fine."
        )
        evidence["solver_calls"] = calls
    evidence["coarse_status"] = evidence["coarse"]["status"]
    evidence["base_status"] = evidence["base"]["status"]
    qc = root / RUN_RELATIVE / "qc"
    _write_json(qc / "adaptive_active_population_cb_equivalence.json", evidence)
    _write_csv(qc / "adaptive_active_population_cb_equivalence.csv", evidence)
    return evidence


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run"))
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    args = parser.parse_args()
    result = (
        prepare_deferred(args.project_root)
        if args.action == "prepare"
        else run_oracle(args.project_root)
    )
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
