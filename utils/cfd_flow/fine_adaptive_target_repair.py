"""Finalize the source-proven Fine adaptive-target repair study.

This research-only helper consumes the already generated Fine mesh and the
separate ``fine_postfix`` logical simulation.  It never launches Seeder,
Coarse, or Base.  Expensive calls are deliberately split: the long Fine run
must already be complete before ``audit``, and Full V2 is a separate explicit
``full-v2`` action allowed only after both steady audits pass.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import time
from pathlib import Path
from typing import Any
from xml.etree import ElementTree

import numpy as np

from .full_timestep_mass_referee import FULL_IDENTITY_GATE, public_step_record, replay_full_timestep
from .io import sha256_file, write_json
from .musubi_boundary_mass_referee import load_mesh_contract
from .restart_decode import read_restart_pdf
from .tau1_base import MPIRUN_WSL, _controller_records, _restart_pairs
from .tau1_grid_convergence import (
    GRID_SPECS,
    PLANE_CONTRACT_SHA256,
    evaluate_repaired_grid_gate,
    run_root,
)
from .tau1_reference_scaled_grid import (
    GridCheckpointAuditor,
    _window_observables,
    finalize_grid_convergence,
    generate_member_lua,
    steady_window_audit,
)


POSTFIX_WSL = "/home/lzy/u3da/tau1_reference_scaled_cbf_20260901/fine_postfix"
POSTFIX_WINDOWS = Path(
    r"\\wsl.localhost\Ubuntu\home\lzy\u3da\tau1_reference_scaled_cbf_20260901\fine_postfix"
)
CANDIDATE_BINARY_WSL = (
    "/home/lzy/apes-worktrees/musubi_adaptive_target_fixed_20260901/build/musubi"
)
CANDIDATE_BINARY_WINDOWS = Path(
    r"\\wsl.localhost\Ubuntu\home\lzy\apes-worktrees\musubi_adaptive_target_fixed_20260901\build\musubi"
)
EXPECTED_CANDIDATE_SHA256 = "1491e7ded30ed15158bd8ba2812d414ad9fa31cbc0128a26d2b1eba779773caa"
ACCEPTED_ITERATION = 1_011_895
AUDIT_ITERATIONS = (404_758, 607_137, 809_516, ACCEPTED_ITERATION)
EXPECTED_MESH_HASHES = {
    "elemlist.lsb": "52fccada55dc83257b7d6f699470137c41ca3f4afeb90e27c7cedcfccd77a043",
    "bnd.lsb": "d9918517d006c3a7fc639df6f4d2927ceb2a97badb80b14fc065e78b44c053a4",
    "qval.lsb": "c9d2f28c6532d9698d49d6158e1619f1cba2872012b2fa7e7def038da200002a",
}


def _postfix_qc(root: Path) -> Path:
    return run_root(root) / "fine_postfix" / "qc"


def _controller_history() -> dict[int, dict[str, Any]]:
    records: dict[int, dict[str, Any]] = {}
    paths = [
        POSTFIX_WINDOWS / "musubi_0000000_0005000_stdout.log",
        POSTFIX_WINDOWS
        / "segments"
        / "segment_0005000_to_1011895"
        / "musubi_stdout.log",
    ]
    for path in paths:
        if path.is_file():
            for record in _controller_records(path.read_text(encoding="utf-8")):
                records[int(record["iteration"])] = record
    return records


def _mesh_hashes(root: Path) -> dict[str, str]:
    mesh = run_root(root) / "fine" / "seeder" / "mesh"
    return {name: sha256_file(mesh / name) for name in EXPECTED_MESH_HASHES}


def audit_postfix_fine(project_root: Path) -> dict[str, Any]:
    """Decode cadence checkpoints and apply both existing steady windows."""

    root = Path(project_root).resolve()
    qc = _postfix_qc(root)
    qc.mkdir(parents=True, exist_ok=True)
    if sha256_file(CANDIDATE_BINARY_WINDOWS) != EXPECTED_CANDIDATE_SHA256:
        raise RuntimeError("candidate binary hash changed")
    mesh_hashes = _mesh_hashes(root)
    if mesh_hashes != EXPECTED_MESH_HASHES:
        raise RuntimeError(f"Fine mesh hash changed: {mesh_hashes}")
    pairs = _restart_pairs(POSTFIX_WINDOWS / "restart")
    missing = [iteration for iteration in AUDIT_ITERATIONS if iteration not in pairs]
    if missing:
        raise RuntimeError(f"Fine cadence checkpoints are incomplete: {missing}")
    controllers = _controller_history()
    auditor = GridCheckpointAuditor(root, "fine")
    samples: dict[int, dict[str, Any]] = {}
    for iteration in AUDIT_ITERATIONS:
        candidates = [value for value in controllers if value <= iteration]
        if not candidates:
            raise RuntimeError(f"controller record missing at {iteration}")
        samples[iteration] = auditor.snapshot(
            iteration, pairs[iteration][1], controllers[max(candidates)]
        )
    ordered = [samples[value] for value in AUDIT_ITERATIONS]
    rho_pass = all(0.9 <= float(item["rho_lattice"]["mean"]) <= 1.1 for item in ordered)
    first = steady_window_audit(ordered[:3], GRID_SPECS["fine"], all_checkpoint_rho_pass=rho_pass)
    second = steady_window_audit(ordered[1:], GRID_SPECS["fine"], all_checkpoint_rho_pass=rho_pass)
    passed = first["status"] == second["status"] == "PASS_NON_REFEREE"
    observables = _window_observables(ordered[1:]) if passed else None
    accepted_binary = pairs[ACCEPTED_ITERATION][1]
    result = {
        "status": "PASS_NON_REFEREE_CONFIRMED" if passed else "FAIL",
        "scientific_status": (
            "CFD_FLOW_FINE_TAU1_STEADY_PASS_PENDING_FULL_V2"
            if passed
            else "CFD_FLOW_TAU1_GRID_MEMBER_STEADY_FAILED"
        ),
        "grid": "fine",
        "logical_simulations": 1,
        "process_launches_before_full_v2": 2,
        "restart_resumes_before_full_v2": 1,
        "operational_recoveries": 1,
        "accepted_iteration": ACCEPTED_ITERATION if passed else None,
        "accepted_physical_time_s": ACCEPTED_ITERATION * GRID_SPECS["fine"].dt_s if passed else None,
        "accepted_restart": (
            {
                "header": str(pairs[ACCEPTED_ITERATION][0]),
                "binary": str(accepted_binary),
                "sha256": sha256_file(accepted_binary),
            }
            if passed
            else None
        ),
        "cell_count": auditor.cells,
        "mesh_hashes": mesh_hashes,
        "physical_plane_contract_sha256": PLANE_CONTRACT_SHA256,
        "steady_audits": [first, second],
        "steady_metrics": second,
        "primary_observables": observables,
        "full_v2_residual": None,
    }
    write_json(qc / "fine_checkpoint_history.json", {"samples": ordered})
    write_json(qc / "fine_steady_acceptance.json", result)
    write_json(
        qc / "fine_physical_plane_metrics.json",
        {
            "status": "PASS" if passed else "BLOCKED_FINE_NOT_ACCEPTED",
            "physical_plane_contract_sha256": PLANE_CONTRACT_SHA256,
            "accepted_window": observables,
            "samples": ordered,
        },
    )
    if not passed:
        raise RuntimeError("CFD_FLOW_TAU1_GRID_MEMBER_STEADY_FAILED")
    return result


def run_postfix_full_v2(project_root: Path) -> dict[str, Any]:
    """Run exactly one candidate-binary step after confirmed Fine steady PASS."""

    root = Path(project_root).resolve()
    qc = _postfix_qc(root)
    acceptance_path = qc / "fine_steady_acceptance.json"
    acceptance = json.loads(acceptance_path.read_text(encoding="utf-8"))
    if acceptance["status"] != "PASS_NON_REFEREE_CONFIRMED":
        raise RuntimeError("Fine steady gates must pass before Full V2")
    pairs = _restart_pairs(POSTFIX_WINDOWS / "restart")
    header, start_binary = pairs[ACCEPTED_ITERATION]
    runtime = POSTFIX_WINDOWS / "full_v2"
    restart = runtime / "restart"
    restart.mkdir(parents=True, exist_ok=True)
    lua = generate_member_lua(
        root,
        "fine",
        maximum_iteration=ACCEPTED_ITERATION + 1,
        restart_header_wsl=f"{POSTFIX_WSL}/restart/{header.name}",
        restart_first_iteration=ACCEPTED_ITERATION + 1,
        restart_interval=1,
        restart_write_wsl=f"{POSTFIX_WSL}/full_v2/restart/",
    )
    lua = lua.replace(
        "simulation_name = 'tau1_reference_scaled_fine'",
        "simulation_name = 'tau1_reference_scaled_fine_postfix_full_v2'",
    ).replace(f"{POSTFIX_WSL.rsplit('/', 1)[0]}/fine/stop", f"{POSTFIX_WSL}/stop")
    (runtime / "musubi.lua").write_text(lua, encoding="utf-8", newline="\n")
    started = time.perf_counter()
    with (runtime / "musubi_stdout.log").open("w", encoding="utf-8", newline="\n") as stdout:
        with (runtime / "musubi_stderr.log").open("w", encoding="utf-8", newline="\n") as stderr:
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
                    "-np",
                    "4",
                    CANDIDATE_BINARY_WSL,
                    f"{POSTFIX_WSL}/full_v2/musubi.lua",
                ],
                cwd=runtime,
                stdout=stdout,
                stderr=stderr,
                check=False,
            )
    elapsed = time.perf_counter() - started
    end_pairs = _restart_pairs(restart)
    end_iteration = ACCEPTED_ITERATION + 1
    if process.returncode != 0 or end_iteration not in end_pairs:
        raise RuntimeError("Fine postfix Full V2 Musubi step failed")
    cells = 400_949
    start_pdf = np.asarray(read_restart_pdf(start_binary, n_elems=cells, n_components=19))
    end_binary = end_pairs[end_iteration][1]
    end_pdf = np.asarray(read_restart_pdf(end_binary, n_elems=cells, n_components=19))
    mesh = load_mesh_contract(
        run_root(root) / "fine" / "seeder" / "mesh",
        expected_cells=cells,
        require_runtime_order=False,
    )
    spec = GRID_SPECS["fine"]
    replay = replay_full_timestep(
        start_pdf,
        end_pdf,
        mesh,
        dx_m=spec.dx_m,
        dt_s=spec.dt_s,
        density_kg_m3=1056.0,
        target_mass_flow_kg_s=2.8901803804796421e-12,
        outlet_pressures_pa=spec.outlet_absolute_pressure_pa,
    )
    residual = float(replay["R_full_one_step_identity"])
    result = {
        "status": "PASS" if residual <= FULL_IDENTITY_GATE else "FAIL",
        "scientific_status": (
            "CFD_FLOW_FINE_TAU1_STEADY_PASS"
            if residual <= FULL_IDENTITY_GATE
            else "CFD_FLOW_TAU1_GRID_MEMBER_STEADY_FAILED"
        ),
        "iteration_start": ACCEPTED_ITERATION,
        "iteration_end": end_iteration,
        "hard_gate": FULL_IDENTITY_GATE,
        "raw_residual": residual,
        "process_wall_clock_s": elapsed,
        "candidate_binary_sha256": sha256_file(CANDIDATE_BINARY_WINDOWS),
        "restart_sha256": {
            str(ACCEPTED_ITERATION): sha256_file(start_binary),
            str(end_iteration): sha256_file(end_binary),
        },
        "referee": public_step_record(replay),
    }
    write_json(qc / "fine_full_v2_referee.json", result)
    destination = run_root(root) / "fine_postfix" / "full_v2"
    destination.mkdir(parents=True, exist_ok=True)
    for name in ("musubi.lua", "musubi_stdout.log", "musubi_stderr.log"):
        shutil.copy2(runtime / name, destination / name)
    acceptance["status"] = "PASS" if result["status"] == "PASS" else "FAIL"
    acceptance["scientific_status"] = result["scientific_status"]
    acceptance["full_v2_residual"] = residual
    acceptance["process_launches"] = 3
    acceptance["restart_resumes"] = 2
    write_json(acceptance_path, acceptance)
    if result["status"] != "PASS":
        raise RuntimeError("CFD_FLOW_TAU1_GRID_MEMBER_STEADY_FAILED")
    return result


def _verification(qc: Path) -> dict[str, Any]:
    result: dict[str, Any] = {}
    pytest_path = qc / "postfix_targeted_pytest.xml"
    if pytest_path.is_file():
        xml = ElementTree.parse(pytest_path).getroot()
        suite = xml if xml.tag == "testsuite" else xml.find("testsuite")
        if suite is not None:
            result["targeted_pytest"] = {
                "status": "PASS"
                if int(suite.attrib.get("failures", 0)) == int(suite.attrib.get("errors", 0)) == 0
                else "FAIL",
                "tests": int(suite.attrib.get("tests", 0)),
                "seconds": float(suite.attrib.get("time", 0.0)),
                "sha256": sha256_file(pytest_path),
            }
    ruff_path = qc / "postfix_targeted_ruff.json"
    if ruff_path.is_file():
        findings = json.loads(ruff_path.read_text(encoding="utf-8-sig"))
        result["targeted_ruff"] = {
            "status": "PASS" if not findings else "FAIL",
            "findings": len(findings),
            "sha256": sha256_file(ruff_path),
        }
    return result


def finalize_postfix_grid(project_root: Path) -> dict[str, Any]:
    """Reuse accepted C/B evidence and run the unchanged three-grid gate."""

    root = Path(project_root).resolve()
    qc = run_root(root) / "qc"
    postfix = json.loads((_postfix_qc(root) / "fine_steady_acceptance.json").read_text(encoding="utf-8"))
    if postfix["status"] != "PASS" or float(postfix["full_v2_residual"]) > FULL_IDENTITY_GATE:
        raise RuntimeError("Fine and Full V2 must pass before C/B/F analysis")
    shutil.copy2(_postfix_qc(root) / "fine_steady_acceptance.json", qc / "fine_steady_acceptance.json")
    result = finalize_grid_convergence(root)
    analyses = result["analyses"]
    gate = evaluate_repaired_grid_gate(analyses)
    result["final_status"] = (
        "CFD_FLOW_TAU1_GRID_CONVERGENCE_PASS"
        if result["status"] == gate["status"] == "PASS"
        else "CFD_FLOW_TAU1_GRID_CONVERGENCE_FAILED"
    )
    result["scientific_status"] = result["final_status"]
    result["verification"] = _verification(qc)
    result["seeder_calls"] = 0
    result["coarse_cfd_launches"] = 0
    result["base_cfd_launches"] = 0
    result["fine_pref_fix_instrumentation_calls"] = 1
    result["fine_postfix_target_oracle_calls"] = 1
    result["fine_logical_simulations"] = 1
    result["fine_process_launches"] = 3
    result["fine_restart_resumes"] = 2
    result["operational_recoveries"] = 1
    result["production_pipeline_modified"] = False
    result["WSS_status"] = "DEFERRED_TO_POST_GRID_PRODUCTION_VALIDATION"
    result["next"] = (
        "PROMOTE VALIDATED TAU1 + CONTINUOUS-Q CFD CONTRACT TO PRODUCTION AND RUN POST-GRID WSS VALIDATION"
        if result["final_status"] == "CFD_FLOW_TAU1_GRID_CONVERGENCE_PASS"
        else "STOP AND REVIEW THE SPECIFIC TRUE GRID-CONVERGENCE METRIC"
    )
    write_json(qc / "grid_convergence_final.json", result)
    return result


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("audit", "full-v2", "finalize"))
    parser.add_argument("project_root", nargs="?", default=Path.cwd())
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    actions = {
        "audit": audit_postfix_fine,
        "full-v2": run_postfix_full_v2,
        "finalize": finalize_postfix_grid,
    }
    result = actions[args.action](Path(args.project_root))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
