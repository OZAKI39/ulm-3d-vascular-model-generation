"""Research-only repaired Tau=1 Coarse/Base/Fine grid study runner.

Only Coarse and Fine are generated and solved.  The accepted Base is consumed
read-only.  Sparse PDFs are decoded offline with the same physical-plane flux
extractor and steady gates used by the accepted Base runner.
"""

from __future__ import annotations

import csv
import json
import math
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence
from xml.etree import ElementTree

import numpy as np

from .apes import parse_mesh_header, windows_to_wsl
from .final_base_geometry_validation import validate_final_base, write_compact_stage_files
from .full_timestep_mass_referee import (
    FULL_IDENTITY_GATE,
    public_step_record,
    replay_full_timestep,
)
from .grid_convergence import mesh_file_contract
from .io import sha256_file, write_json
from .musubi_boundary_mass_referee import (
    load_mesh_contract,
    replay_boundary_step,
    runtime_solid_cells,
)
from .physical_port_flux import (
    FLUX_ALGORITHM_REVISION,
    FLUX_DEFINITION,
    MLS_CONDITION_GATE,
    _evaluate_prepared_plane,
    _mesh_origin_dx,
    _prepare_plane_numerics,
    plane_from_v3_record,
)
from .qvalue_contract_forensics import vascular_wall_qvalue_distribution
from .restart_decode import read_restart_pdf, reconstruct_macroscopic_field
from .tau1_base import (
    BULK_NU_M2_S,
    MPI_RANKS,
    MPIRUN_WSL,
    MUSUBI_SHA256,
    MUSUBI_WSL,
    NU_M2_S,
    OUTLET_GAUGE_PRESSURE_PA,
    PROJECT_WSL,
    RHO_KG_M3,
    TARGET_MASS_FLOW_KG_S,
    TARGET_Q_M3_S,
    Tau1BaseRuntimeContract,
    _controller_records,
    _physics_and_boundaries_lua,
    _restart_pairs,
)
from .tau1_grid_convergence import (
    BASE_ACCEPTED_ITERATION,
    BASE_ACCEPTED_RESTART_SHA256,
    BASE_CFD_RUN,
    BASE_EXPECTED_CELLS,
    BASE_FULL_REFEREE_RESIDUAL,
    BASE_MESH_RUN,
    BASE_MESH_SHA256,
    CBF_RUN_NAME as RUN_NAME,
    GRID_SPECS,
    PLANE_CONTRACT_FILE_SHA256,
    PLANE_CONTRACT_SHA256,
    PLATEAU_AUDITS,
    PLATEAU_MIN_PHYSICAL_TIME_S,
    PLATEAU_RELATIVE_CHANGE,
    PRIMARY_METRICS,
    SEEDER_BINARY_SHA256,
    SEEDER_BINARY_WSL,
    build_primary_analyses,
    evaluate_repaired_grid_gate,
    render_repaired_seeder_config,
    restart_resume_contract,
    run_root,
    seeder_physical_spatial_signature,
    write_grid_design,
    Tau1GridSpec,
)


PORTS = ("inlet", "outlet_01", "outlet_02", "outlet_03")
OUTLETS = PORTS[1:]
RUNTIME_ROOT_WSL = "/home/lzy/u3da/tau1_reference_scaled_cbf_20260901"
RUNTIME_ROOT_WINDOWS = Path(
    r"\\wsl.localhost\Ubuntu\home\lzy\u3da\tau1_reference_scaled_cbf_20260901"
)
PLANE_CONTRACT_RELATIVE = Path(
    "outputs/cfd_flow/healthy_mouse_capillary_tau1_grid_convergence_"
    "anchor003274_20260831/qc/physical_port_flux_plane_contract_v3.json"
)
BASE_FINAL_RELATIVE = Path("outputs/cfd_flow") / BASE_CFD_RUN / "qc" / (
    "reference_scaled_base_final.json"
)
MASS_GATE = 0.01
VELOCITY_GATE = 0.01
PRESSURE_GATE = 0.005
INLET_GATE = 0.01
FRACTION_DRIFT_GATE = 0.01
Q_DENSITY_GATE = 0.01
CONTROLLER_GATE = 1.0e-8
MAX_LATTICE_SPEED = 0.05
BACKFLOW_FRACTION = 0.05
RHO_GATE = (0.9, 1.1)
INITIAL_SAFETY_ITERATIONS = 5_000
MAX_OPERATIONAL_RECOVERIES = 5
EXPECTED_START_HEAD = "4edbc1120ed7d4d961118004ea91759847ef2762"


@dataclass(frozen=True, slots=True)
class GridPaths:
    project_root: Path
    label: str

    @property
    def run_root(self) -> Path:
        return run_root(self.project_root)

    @property
    def member(self) -> Path:
        return self.run_root / self.label

    @property
    def mesh(self) -> Path:
        return self.member / "seeder" / "mesh"

    @property
    def qc(self) -> Path:
        return self.member / "qc"

    @property
    def segments(self) -> Path:
        return self.member / "segments"

    @property
    def runtime_wsl(self) -> str:
        return f"{RUNTIME_ROOT_WSL}/{self.label}"

    @property
    def runtime(self) -> Path:
        return RUNTIME_ROOT_WINDOWS / self.label


def _git(root: Path, *arguments: str) -> str:
    return subprocess.run(
        ["git", "-C", str(root), *arguments], check=True,
        capture_output=True, text=True, encoding="utf-8",
    ).stdout.strip()


def _base_mesh(root: Path) -> Path:
    return root / "outputs" / "cfd_flow" / BASE_MESH_RUN / "seeder" / "mesh"


def _base_seeder(root: Path) -> Path:
    return root / "outputs" / "cfd_flow" / BASE_MESH_RUN / "seeder" / "seeder.lua"


def _plane_contract(root: Path) -> Path:
    return root / PLANE_CONTRACT_RELATIVE


def _binary_windows(wsl_path: str) -> Path:
    return Path(r"\\wsl.localhost\Ubuntu" + wsl_path.replace("/", "\\"))


def _state_path(root: Path) -> Path:
    return run_root(root) / "qc" / "cbf_run_state.json"


def _update_state(root: Path, **values: Any) -> dict[str, Any]:
    path = _state_path(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    state = json.loads(path.read_text(encoding="utf-8")) if path.is_file() else {
        "stage": "STAGE_0_PREFLIGHT",
        "members": {},
        "seeder_calls": {"coarse": 0, "base": 0, "fine": 0},
        "logical_cfd": {"coarse": 0, "base": 0, "fine": 0},
        "process_launches": 0,
        "restart_resumes": 0,
        "operational_recoveries": 0,
        "next_action": "ZERO_SOLVER_PREFLIGHT",
    }
    state.update(values)
    write_json(path, state)
    return state


def _production_modified(root: Path) -> bool:
    return bool(_git(root, "diff", "--", "cfd_flow.py", "configs/cfd_flow.yaml", "utils/cfd_flow/pipeline.py"))


def preflight(project_root: Path) -> dict[str, Any]:
    """Prove all frozen inputs before the first Seeder or Musubi call."""

    root = Path(project_root).resolve()
    qc = run_root(root) / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    design = write_grid_design(root)
    base_final_path = root / BASE_FINAL_RELATIVE
    base_final = json.loads(base_final_path.read_text(encoding="utf-8"))
    accepted = Path(base_final["accepted_restart"]["binary"])
    plane = _plane_contract(root)
    plane_record = json.loads(plane.read_text(encoding="utf-8"))
    base_text = _base_seeder(root).read_text(encoding="utf-8")
    rendered = {
        label: render_repaired_seeder_config(base_text, GRID_SPECS[label])
        for label in ("coarse", "fine")
    }
    signatures = {
        "base": seeder_physical_spatial_signature(base_text),
        **{
            label: seeder_physical_spatial_signature(text)
            for label, text in rendered.items()
        },
    }
    mesh_hashes = {
        name: sha256_file(_base_mesh(root) / name) for name in BASE_MESH_SHA256
    }
    checks = {
        "branch": _git(root, "branch", "--show-current")
        == "codex/cfd-wall-force-numerics-validated-sync-20260830",
        "reference_head": _git(root, "rev-parse", "HEAD") == EXPECTED_START_HEAD,
        "base_status": base_final["status"]
        == "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_PASS",
        "base_iteration": int(base_final["accepted_restart"]["iteration"])
        == BASE_ACCEPTED_ITERATION,
        "base_restart_exists": accepted.is_file(),
        "base_restart_hash": accepted.is_file()
        and sha256_file(accepted) == BASE_ACCEPTED_RESTART_SHA256,
        "base_mesh_hashes": mesh_hashes == BASE_MESH_SHA256,
        "musubi_binary_hash": sha256_file(_binary_windows(MUSUBI_WSL))
        == MUSUBI_SHA256,
        "seeder_binary_hash": sha256_file(_binary_windows(SEEDER_BINARY_WSL))
        == SEEDER_BINARY_SHA256,
        "grid_design": design["status"] == "PASS",
        "physical_plane_contract": plane_record["contract_sha256"]
        == PLANE_CONTRACT_SHA256,
        "physical_plane_file_hash": sha256_file(plane)
        == PLANE_CONTRACT_FILE_SHA256,
        "seeder_geometry_signature_identical": len(set(signatures.values())) == 1,
        "production_pipeline_unmodified": not _production_modified(root),
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "actual_head": _git(root, "rev-parse", "HEAD"),
        "checks": checks,
        "base_final": str(base_final_path),
        "base_accepted_restart_sha256": sha256_file(accepted) if accepted.is_file() else None,
        "base_mesh_hashes": mesh_hashes,
        "seeder_geometry_signatures": signatures,
        "physical_plane_contract_sha256": PLANE_CONTRACT_SHA256,
        "physical_flux_extractor_revision": FLUX_ALGORITHM_REVISION,
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "musubi_calls": 0,
    }
    write_json(qc / "stage0_preflight.json", result)
    if result["status"] != "PASS":
        _update_state(root, stage="STOPPED", next_action="STOP_PREFLIGHT_FAILURE", preflight=result)
        raise RuntimeError(f"CBF zero-solver preflight failed: {checks}")
    _update_state(root, stage="STAGE_0B_BASE_READONLY", preflight=result, next_action="BASE_READONLY_REAUDIT")
    return result


def _trapezoid_mean(samples: Sequence[Mapping[str, Any]], getter: Any) -> float:
    values = [float(getter(item)) for item in samples]
    times = [int(item["iteration"]) for item in samples]
    integral = math.fsum(
        0.5 * (values[index] + values[index + 1])
        * (times[index + 1] - times[index])
        for index in range(len(values) - 1)
    )
    return integral / (times[-1] - times[0])


def _window_observables(samples: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    means = {
        "inlet_gauge_pressure_pa": _trapezoid_mean(
            samples, lambda item: item["inlet_gauge_pressure_pa"]
        )
    }
    for label in OUTLETS:
        means[f"DeltaP{label[-2:]}_pa"] = _trapezoid_mean(
            samples, lambda item, label=label: item["pressure_drops_pa"][label]
        )
    q = {
        label: _trapezoid_mean(
            samples,
            lambda item, label=label: item["ports"][label]["Q_velocity_m3_s"],
        )
        for label in PORTS
    }
    qout = math.fsum(q[label] for label in OUTLETS)
    means.update(
        {
            "Qin_m3_s": q["inlet"],
            "Q1_m3_s": q["outlet_01"],
            "Q2_m3_s": q["outlet_02"],
            "Q3_m3_s": q["outlet_03"],
            "Qout_sum_m3_s": qout,
            "outlet_01_flow_fraction": q["outlet_01"] / qout,
            "outlet_02_flow_fraction": q["outlet_02"] / qout,
            "outlet_03_flow_fraction": q["outlet_03"] / qout,
            "flux_definition": FLUX_DEFINITION,
            "window_sample_iterations": [int(item["iteration"]) for item in samples],
            "window_statistic": "PHYSICAL_TIME_TRAPEZOIDAL_MEAN",
        }
    )
    return means


def base_readonly_reaudit(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    source = root / "outputs" / "cfd_flow" / BASE_CFD_RUN / "qc"
    final = json.loads((source / "reference_scaled_base_final.json").read_text(encoding="utf-8"))
    history_record = json.loads(
        (source / "reference_scaled_base_checkpoint_history.json").read_text(encoding="utf-8")
    )
    history = {int(item["iteration"]): item for item in history_record["checkpoint_history"]}
    iterations = [int(value) for value in final["steady_audit"]["window_iterations"]]
    samples = [history[value] for value in iterations]
    endpoint = samples[-1]
    accepted = Path(final["accepted_restart"]["binary"])
    checks = {
        "base_final_pass": final["status"] == "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_PASS",
        "accepted_iteration": int(final["accepted_restart"]["iteration"])
        == BASE_ACCEPTED_ITERATION,
        "accepted_hash": sha256_file(accepted) == BASE_ACCEPTED_RESTART_SHA256,
        "history_window_present": len(samples) == 3,
        "endpoint_inlet_reproduced": math.isclose(
            float(endpoint["inlet_gauge_pressure_pa"]), 531.431946845226,
            rel_tol=1.0e-12, abs_tol=1.0e-9,
        ),
        "endpoint_qin_reproduced": math.isclose(
            float(endpoint["ports"]["inlet"]["Q_velocity_m3_s"]),
            2.728393297831303e-15, rel_tol=1.0e-12, abs_tol=0.0,
        ),
        "steady_qc_pass": final["steady_audit"]["status"] == "PASS_NON_REFEREE",
        "full_referee_pass": json.loads(
            (source / "reference_scaled_base_full_referee.json").read_text(encoding="utf-8")
        )["status"] == "PASS",
        "no_new_musubi": True,
        "no_new_seeder": True,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "accepted_iteration": BASE_ACCEPTED_ITERATION,
        "accepted_physical_time_s": BASE_ACCEPTED_ITERATION * GRID_SPECS["base"].dt_s,
        "accepted_restart_sha256": BASE_ACCEPTED_RESTART_SHA256,
        "cell_count": BASE_EXPECTED_CELLS,
        "mesh_hashes": BASE_MESH_SHA256,
        "steady_metrics": final["steady_audit"],
        "full_v2_residual": BASE_FULL_REFEREE_RESIDUAL,
        "endpoint": endpoint,
        "primary_observables": _window_observables(samples),
        "checks": checks,
        "seeder_calls": 0,
        "musubi_calls": 0,
        "read_only": True,
    }
    qc = run_root(root) / "qc"
    write_json(qc / "base_readonly_acceptance.json", result)
    write_json(
        qc / "base_full_referee.json",
        {
            "status": "PASS",
            "reused_read_only": True,
            "R_full_one_step_identity": BASE_FULL_REFEREE_RESIDUAL,
            "hard_gate": FULL_IDENTITY_GATE,
            "source": str(source / "reference_scaled_base_full_referee.json"),
        },
    )
    if result["status"] != "PASS":
        _update_state(root, stage="STOPPED", next_action="CFD_FLOW_CBF_BASE_REPRODUCTION_FAILED")
        raise RuntimeError("CFD_FLOW_CBF_BASE_REPRODUCTION_FAILED")
    _update_state(root, stage="STAGE_1_MESHES", base_readonly=result, next_action="GENERATE_COARSE_FINE_MESHES")
    return result


def _copy_geometry(root: Path, destination: Path) -> None:
    source = root / "outputs" / "cfd_flow" / BASE_MESH_RUN / "geometry" / "geometry_solver_m"
    destination.mkdir(parents=True, exist_ok=False)
    for label in ("wall", "outlet_01", "outlet_02", "outlet_03"):
        shutil.copy2(source / f"{label}.stl", destination / f"{label}.stl")


def _physical_plane_support(root: Path, mesh_dir: Path, expected_cells: int) -> dict[str, Any]:
    contract = load_mesh_contract(
        mesh_dir, expected_cells=expected_cells, allow_zero_normals=True,
        require_runtime_order=False,
    )
    origin, dx = _mesh_origin_dx(mesh_dir)
    centers = origin + (contract.cell_ijk.astype(np.float64) + 0.5) * dx
    solid = runtime_solid_cells(contract)
    fluid_mask = np.ones(expected_cells, dtype=bool)
    fluid_mask[np.asarray(sorted(solid), dtype=np.int64)] = False
    plane_record = json.loads(_plane_contract(root).read_text(encoding="utf-8"))
    ports: dict[str, Any] = {}
    for label in PORTS:
        central = plane_record["ports"][label]["planes"]["central"]
        plane = plane_from_v3_record(label, "central", central)
        prepared = _prepare_plane_numerics(plane, centers[fluid_mask], dx_m=dx)
        stencils = prepared["stencil_qc"]
        ports[label] = {
            "origin_m": central["origin_m"],
            "unit_normal": central["unit_normal"],
            "aperture_area_relative_error": prepared["quadrature"].area_relative_error,
            "stencils": stencils,
            "invalid_stencils": sum(int(value["invalid_stencil_count"]) for value in stencils.values()),
            "maximum_condition_number": max(float(value["max_condition_number"]) for value in stencils.values()),
        }
    checks = {
        "contract_hash": plane_record["contract_sha256"] == PLANE_CONTRACT_SHA256,
        "all_apertures_valid": max(
            float(item["aperture_area_relative_error"]) for item in ports.values()
        ) <= 1.0e-10,
        "invalid_interpolation_stencils_zero": sum(
            int(item["invalid_stencils"]) for item in ports.values()
        ) == 0,
        "condition_number_le_1e8": max(
            float(item["maximum_condition_number"]) for item in ports.values()
        ) <= MLS_CONDITION_GATE,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "contract_sha256": PLANE_CONTRACT_SHA256,
        "extractor_revision": FLUX_ALGORITHM_REVISION,
        "ports": ports,
        "checks": checks,
    }


def seed_and_audit_member(project_root: Path, label: str) -> dict[str, Any]:
    root = Path(project_root).resolve()
    if label not in {"coarse", "fine"}:
        raise ValueError(label)
    paths = GridPaths(root, label)
    manifest_path = run_root(root) / "qc" / f"{label}_mesh_contract.json"
    if manifest_path.is_file():
        existing = json.loads(manifest_path.read_text(encoding="utf-8"))
        if existing.get("status") == "PASS" and mesh_file_contract(paths.mesh)["status"] == "PASS":
            return existing
        raise RuntimeError(f"existing {label} mesh evidence is not reusable")
    spec = GRID_SPECS[label]
    state = _update_state(root)
    own_complete = (
        int(state["seeder_calls"].get(label, 0)) == 1
        and mesh_file_contract(paths.mesh)["status"] == "PASS"
    )
    text = render_repaired_seeder_config(
        _base_seeder(root).read_text(encoding="utf-8"), spec
    )
    lua = paths.member / "seeder" / "seeder.lua"
    if own_complete:
        process_returncode = 0
        stdout = (lua.parent / "seeder_stdout.log").read_text(encoding="utf-8")
        match = re.search(r"Done with Seeder in\s+([0-9.]+)\s+s", stdout)
        elapsed = float(match.group(1)) if match else 0.0
    else:
        paths.qc.mkdir(parents=True, exist_ok=False)
        paths.mesh.mkdir(parents=True, exist_ok=False)
        _copy_geometry(root, paths.member / "geometry" / "geometry_solver_m")
        lua.write_text(text, encoding="utf-8", newline="\n")
        command = (
            f"cd '{windows_to_wsl(lua.parent, 'Ubuntu')}' && "
            f"'{SEEDER_BINARY_WSL}' seeder.lua"
        )
        started = time.perf_counter()
        process = subprocess.run(
            ["wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", command],
            capture_output=True, text=True, encoding="utf-8", errors="replace",
            timeout=1800, check=False,
        )
        elapsed = time.perf_counter() - started
        process_returncode = process.returncode
        (lua.parent / "seeder_stdout.log").write_text(process.stdout, encoding="utf-8")
        (lua.parent / "seeder_stderr.log").write_text(process.stderr, encoding="utf-8")
        state["seeder_calls"][label] = 1
        write_json(_state_path(root), state)
    files = mesh_file_contract(paths.mesh)
    if process_returncode != 0 or files["status"] != "PASS":
        _update_state(root, stage="STOPPED", next_action="CFD_FLOW_TAU1_GRID_MESH_CONTRACT_FAILED")
        raise RuntimeError(f"{label} Seeder failed: returncode={process_returncode}")
    header = parse_mesh_header(paths.mesh)
    cells = int(header["fluid_element_count"])
    origin, actual_dx = _mesh_origin_dx(paths.mesh)
    geometry = validate_final_base(
        base_run=paths.member,
        continuous_surface=(
            root / "outputs" / "cfd_flow"
            / "axis_aligned_inlet_geometry_anchor003274_20260829_111451"
            / "geometry" / "cfd_surface_axis_aligned_inlet_m.stl"
        ),
    )
    geometry["solver_calls"] = {
        "seeder_semantic_calls": 1,
        "small_musubi_semantic_calls": 0,
        "vascular_musubi_semantic_calls": 0,
        "launch_failures": 0,
        "preflight_failures": 0,
    }
    write_json(paths.qc / "final_base_geometry_validation.json", geometry)
    write_compact_stage_files(paths.member, geometry)
    q_distribution = vascular_wall_qvalue_distribution(paths.mesh)
    plane_support = _physical_plane_support(root, paths.mesh, cells)
    mesh_hashes = {
        name: sha256_file(paths.mesh / name)
        for name in ("elemlist.lsb", "bnd.lsb", "qval.lsb")
    }
    checks = {
        "seeder_returncode": process_returncode == 0,
        "all_files": files["status"] == "PASS",
        "positive_cells": cells > 0,
        "actual_dx": math.isclose(actual_dx, spec.dx_m, rel_tol=0.0, abs_tol=1.0e-18),
        "uniform_level": header["minimum_level"] == header["maximum_level"] == spec.root_level,
        "all_boundary_labels": all(
            int(header["boundary_cell_counts"].get(name, 0)) > 0
            for name in ("wall", "inlet", "outlet_01", "outlet_02", "outlet_03")
        ),
        "continuous_q": q_distribution["status"] == "PASS",
        "geometry_oracle": geometry["status"] == "PASS",
        "connectivity": geometry["connectivity"]["status"] == "PASS",
        "physical_plane_support": plane_support["status"] == "PASS",
        "physical_geometry_signature": seeder_physical_spatial_signature(text)
        == seeder_physical_spatial_signature(_base_seeder(root).read_text(encoding="utf-8")),
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "grid": label,
        "dx_m": actual_dx,
        "dt_s": spec.dt_s,
        "cell_count": cells,
        "mesh_hashes": mesh_hashes,
        "mesh_origin_m": origin.tolist(),
        "boundary_cell_counts": header["boundary_cell_counts"],
        "seeder_binary_sha256": SEEDER_BINARY_SHA256,
        "seeder_wall_clock_s": elapsed,
        "seeder_calls": 1,
        "geometry_validation": geometry,
        "physical_plane_support": plane_support,
        "checks": checks,
    }
    write_json(manifest_path, result)
    if result["status"] != "PASS":
        _update_state(root, stage="STOPPED", next_action="CFD_FLOW_TAU1_GRID_MESH_CONTRACT_FAILED")
        raise RuntimeError(f"CFD_FLOW_TAU1_GRID_MESH_CONTRACT_FAILED: {label}")
    state = _update_state(root)
    state["members"].setdefault(label, {})["mesh"] = result
    state["next_action"] = "GENERATE_FINE_MESH" if label == "coarse" else "RUN_COARSE_CFD"
    write_json(_state_path(root), state)
    return result


def _runtime_contract(root: Path, label: str) -> dict[str, Any]:
    paths = GridPaths(root, label)
    spec = GRID_SPECS[label]
    return {
        "mesh_hashes": {
            name: sha256_file(paths.mesh / name)
            for name in ("elemlist.lsb", "bnd.lsb", "qval.lsb")
        },
        "dx_m": spec.dx_m,
        "dt_s": spec.dt_s,
        "rho0_kg_m3": RHO_KG_M3,
        "nu_m2_s": NU_M2_S,
        "bulk_nu_m2_s": BULK_NU_M2_S,
        "tau": spec.tau,
        "omega": spec.omega,
        "pressure_reference_pa": spec.pressure_reference_pa,
        "outlet_gauge_pressure_pa": dict(OUTLET_GAUGE_PRESSURE_PA),
        "outlet_absolute_pressure_pa": spec.outlet_absolute_pressure_pa,
        "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
        "binary_sha256": MUSUBI_SHA256,
        "physical_plane_contract_sha256": PLANE_CONTRACT_SHA256,
        "layout": "d3q19",
        "collision": "bgk",
        "boundary_contract": {
            "wall": "wall_libb_continuous_q",
            "inlet": "adaptive_flux_pressure",
            "outlets": "pressure_eq",
        },
    }


def _contract_object(spec: Tau1GridSpec) -> Tau1BaseRuntimeContract:
    return Tau1BaseRuntimeContract(
        dx_m=spec.dx_m,
        dt_s=spec.dt_s,
        pressure_reference_pa=spec.pressure_reference_pa,
    )


def generate_member_lua(
    root: Path,
    label: str,
    *,
    maximum_iteration: int,
    restart_header_wsl: str | None,
    restart_first_iteration: int,
    restart_interval: int,
    restart_write_wsl: str | None = None,
) -> str:
    paths = GridPaths(root, label)
    spec = GRID_SPECS[label]
    read = f"read='{restart_header_wsl}', " if restart_header_wsl else ""
    write = restart_write_wsl or f"{paths.runtime_wsl}/restart/"
    mesh_wsl = windows_to_wsl(paths.mesh, "Ubuntu").rstrip("/")
    sim_interval = max(1, spec.short_window_iterations // 100)
    return f"""-- Fresh/resumable repaired Tau1 {label} grid member.
simulation_name = 'tau1_reference_scaled_{label}'
printRuntimeInfo = true
timing_file = 'tracking/timing.res'
mesh = '{mesh_wsl}/'
scaling = 'diffusive'
logging = {{level=5}}
maximum_iterations = {int(maximum_iteration)}
{_physics_and_boundaries_lua(_contract_object(spec))}
sim_control = {{
  time_control={{max={{iter=maximum_iterations}}, interval={{iter={sim_interval}}}}},
  abort_criteria={{stop_file='{paths.runtime_wsl}/stop'}}
}}
restart = {{{read}write='{write}', timeformat={{use_iter=true}},
  time_control={{min={{iter={int(restart_first_iteration)}}}, max={{iter=maximum_iterations}},
    interval={{iter={int(restart_interval)}}}}}
}}
"""


def member_lua_contract(
    text: str,
    root: Path,
    label: str,
    *,
    maximum_iteration: int,
    restart_header_wsl: str | None,
) -> dict[str, Any]:
    paths = GridPaths(root, label)
    spec = GRID_SPECS[label]
    mesh_wsl = windows_to_wsl(paths.mesh, "Ubuntu").rstrip("/")
    checks = {
        "member": f"simulation_name = 'tau1_reference_scaled_{label}'" in text,
        "mesh": f"mesh = '{mesh_wsl}/'" in text,
        "maximum": f"maximum_iterations = {maximum_iteration}" in text,
        "dx": f"dx = {spec.dx_m:.17g}" in text,
        "dt": f"dt = {spec.dt_s:.17g}" in text,
        "dynamic_reference": f"pressure_reference_phy = {spec.pressure_reference_pa:.17g}" in text,
        "tau_one": math.isclose(spec.tau, 1.0, abs_tol=1.0e-14),
        "fresh_or_own_restart": (
            f"read='{restart_header_wsl}'" in text if restart_header_wsl else "read='" not in text
        ),
        "target": f"mass_flowrate={TARGET_MASS_FLOW_KG_S:.17g}" in text,
        "boundaries": "kind='wall_libb'" in text
        and "kind='adaptive_flux_pressure'" in text
        and text.count("kind='pressure_eq'") == 3,
        "iteration_filenames": "timeformat={use_iter=true}" in text,
        "no_full_field_tracking": all(
            token not in text for token in ("asciiSpatial", "format='vtk'", "shape={kind='all'}")
        ),
        "four_ranks": MPI_RANKS == 4,
    }
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def _launcher_text(short_interval: int) -> str:
    return f"""#!/usr/bin/env bash
set -euo pipefail
SEGMENT_DIR="${{1:?segment directory required}}"
cd "$SEGMENT_DIR"
mkdir -p tracking
set +e
'{MPIRUN_WSL}' --bind-to core --map-by core --report-bindings -np {MPI_RANKS} \\
  '{MUSUBI_WSL}' musubi.lua 2> musubi_stderr.log | awk -v interval={short_interval} '
/ADAPTIVE_FLUX_PRESSURE/ {{
  last=$0
  if (match($0,/iter=[0-9]+/)) {{
    value=substr($0,RSTART+5,RLENGTH-5)+0
    if (value % interval == 0) {{ print $0; fflush() }}
  }}
  next
}}
{{ print $0; fflush() }}
END {{ if (last != "") print last }}
' > musubi_stdout.log
rc=${{PIPESTATUS[0]}}
set -e
[[ "$rc" -eq 0 ]]
grep -q 'SUCCESSFUL run' musubi_stdout.log
grep -q 'ADAPTIVE_FLUX_PRESSURE' musubi_stdout.log
printf 'SEGMENT_SEMANTIC_SUCCESS=PASS\n' > semantic_status.log
"""


def _controller_history(paths: GridPaths) -> dict[int, dict[str, Any]]:
    records: dict[int, dict[str, Any]] = {}
    for root in (paths.runtime / "segments", paths.segments):
        for path in sorted(root.glob("*/musubi_stdout.log")):
            for record in _controller_records(path.read_text(encoding="utf-8")):
                records[int(record["iteration"])] = record
    return records


class GridCheckpointAuditor:
    def __init__(self, root: Path, label: str) -> None:
        self.root = root
        self.label = label
        self.paths = GridPaths(root, label)
        self.spec = GRID_SPECS[label]
        self.cells = int(parse_mesh_header(self.paths.mesh)["fluid_element_count"])
        self.mesh = load_mesh_contract(
            self.paths.mesh,
            expected_cells=self.cells,
            require_runtime_order=False,
        )
        solid = runtime_solid_cells(self.mesh)
        self.fluid_mask = np.ones(self.cells, dtype=bool)
        self.fluid_mask[np.asarray(sorted(solid), dtype=np.int64)] = False
        origin, dx = _mesh_origin_dx(self.paths.mesh)
        centers = origin + (self.mesh.cell_ijk.astype(np.float64) + 0.5) * dx
        contract = json.loads(_plane_contract(root).read_text(encoding="utf-8"))
        self.prepared = {}
        for port in PORTS:
            record = contract["ports"][port]["planes"]["central"]
            plane = plane_from_v3_record(port, "central", record)
            self.prepared[port] = _prepare_plane_numerics(
                plane, centers[self.fluid_mask], dx_m=self.spec.dx_m
            )

    def snapshot(
        self, iteration: int, binary: Path, controller: Mapping[str, Any]
    ) -> dict[str, Any]:
        spec = self.spec
        pdf = np.asarray(
            read_restart_pdf(binary, n_elems=self.cells, n_components=19),
            dtype=np.float64,
        )
        field = reconstruct_macroscopic_field(
            pdf, dx_m=spec.dx_m, dt_s=spec.dt_s, rho0_kg_m3=RHO_KG_M3
        )
        density = np.asarray(field.density_lattice, dtype=np.float64)
        velocity = np.asarray(field.velocity_phy, dtype=np.float64)
        speed = np.linalg.norm(velocity, axis=1)
        ports: dict[str, Any] = {}
        for label in PORTS:
            value = _evaluate_prepared_plane(
                self.prepared[label], velocity[self.fluid_mask], density[self.fluid_mask]
            )
            q = float(value["physical_q_m3_s"])
            q_rho = float(value["mass_normalized_volume_flux_m3_s"])
            ports[label] = {
                "Q_velocity_m3_s": q,
                "Q_rho_u_over_rho0_m3_s": q_rho,
                "area_weighted_rho_lattice": value["area_weighted_density_lattice"],
                "R_quadrature_Q": value["R_quadrature_Q"],
            }
        boundary = replay_boundary_step(
            pdf, self.mesh, dx_m=spec.dx_m, dt_s=spec.dt_s,
            density_kg_m3=RHO_KG_M3,
            target_mass_flow_kg_s=TARGET_MASS_FLOW_KG_S,
            outlet_pressures_pa=spec.outlet_absolute_pressure_pa,
        )
        inlet_rho = float(boundary["details"]["inlet"]["rho"])
        inlet_gauge = (inlet_rho - 1.0) * spec.pressure_reference_pa
        qin = float(ports["inlet"]["Q_velocity_m3_s"])
        qout = math.fsum(float(ports[label]["Q_velocity_m3_s"]) for label in OUTLETS)
        fraction_denominator = qout if abs(qout) > np.finfo(float).tiny else np.finfo(float).tiny
        q_scale = max(abs(qin), np.finfo(float).tiny)
        q_density_checks = {}
        for label in PORTS:
            q = float(ports[label]["Q_velocity_m3_s"])
            q_rho = float(ports[label]["Q_rho_u_over_rho0_m3_s"])
            denominator = max(abs(q_rho), np.finfo(float).tiny) if (
                label == "inlet" or abs(q_rho) >= 0.05 * q_scale
            ) else q_scale
            residual = abs(q - q_rho) / denominator
            q_density_checks[label] = {"residual": residual, "pass": residual <= Q_DENSITY_GATE}
        return {
            "iteration": int(iteration),
            "physical_time_s": int(iteration) * spec.dt_s,
            "restart_sha256": sha256_file(binary),
            "rho_lattice": {
                "mean": float(np.mean(density)), "median": float(np.median(density)),
                "p1": float(np.percentile(density, 1.0)),
                "p99": float(np.percentile(density, 99.0)),
            },
            "mean_speed_m_s": float(np.mean(speed)),
            "inlet_gauge_pressure_pa": inlet_gauge,
            "pressure_drops_pa": {
                label: inlet_gauge - gauge
                for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
            },
            "ports": ports,
            "Qout_sum_m3_s": qout,
            "physical_volume_closure": abs(qin - qout) / q_scale,
            "flow_fractions": {
                label: float(ports[label]["Q_velocity_m3_s"]) / fraction_denominator
                for label in OUTLETS
            },
            "Q_density_consistency": q_density_checks,
            "minimum_pdf": float(np.min(pdf)),
            "maximum_lattice_speed": float(np.max(speed) * spec.dt_s / spec.dx_m),
            "all_finite": bool(
                np.all(np.isfinite(pdf)) and np.all(np.isfinite(density))
                and np.all(np.isfinite(velocity))
            ),
            "controller": dict(controller),
        }


def steady_window_audit(
    samples: Sequence[Mapping[str, Any]], spec: Tau1GridSpec,
    *, all_checkpoint_rho_pass: bool,
) -> dict[str, Any]:
    if len(samples) != 3:
        raise ValueError("three long/short/end checkpoints required")
    iterations = [int(item["iteration"]) for item in samples]
    if any(iterations[index + 1] - iterations[index] != spec.short_window_iterations for index in (0, 1)):
        raise ValueError(f"non-cadenced physical windows: {iterations}")

    def q(item: Mapping[str, Any], label: str) -> float:
        return float(item["ports"][label]["Q_velocity_m3_s"])

    def flows(window: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
        inlet = _trapezoid_mean(window, lambda item: q(item, "inlet"))
        outlets = {
            label: _trapezoid_mean(window, lambda item, label=label: q(item, label))
            for label in OUTLETS
        }
        total = math.fsum(outlets.values())
        return {
            "mean_Qin_m3_s": inlet,
            "mean_outlets_m3_s": outlets,
            "mean_Qout_sum_m3_s": total,
            "R_mass": abs(inlet - total) / max(abs(inlet), np.finfo(float).tiny),
            "significant_backflow_by_outlet": {
                label: value < 0.0 and abs(value) > BACKFLOW_FRACTION * abs(inlet)
                for label, value in outlets.items()
            },
        }

    a, b, c = samples
    short = flows((b, c))
    long = flows((a, b, c))
    r_velocity = abs(float(c["mean_speed_m_s"]) - float(b["mean_speed_m_s"])) / max(
        abs(float(c["mean_speed_m_s"])), 1.0e-12
    )
    pressure_residuals = {
        "inlet_gauge_pressure_pa": abs(
            float(c["inlet_gauge_pressure_pa"]) - float(b["inlet_gauge_pressure_pa"])
        ) / max(abs(float(c["inlet_gauge_pressure_pa"])), 1.0),
        **{
            label: abs(
                float(c["pressure_drops_pa"][label]) - float(b["pressure_drops_pa"][label])
            ) / max(abs(float(c["pressure_drops_pa"][label])), 1.0)
            for label in OUTLETS
        },
    }
    r_pressure = max(pressure_residuals.values())
    r_inlet = abs(short["mean_Qin_m3_s"] - TARGET_Q_M3_S) / TARGET_Q_M3_S
    fraction_drift = {
        label: max(float(item["flow_fractions"][label]) for item in samples)
        - min(float(item["flow_fractions"][label]) for item in samples)
        for label in OUTLETS
    }
    controller = c["controller"]
    target_error = abs(float(controller["target_lattice"]) - spec.target_lattice) / spec.target_lattice
    gates = {
        "R_mass_short": short["R_mass"] <= MASS_GATE,
        "R_mass_long": long["R_mass"] <= MASS_GATE,
        "physical_volume_closure": float(c["physical_volume_closure"]) <= MASS_GATE,
        "R_velocity": r_velocity <= VELOCITY_GATE,
        "R_pressure": r_pressure <= PRESSURE_GATE,
        "R_inlet": r_inlet <= INLET_GATE,
        "flow_fraction_drift": max(fraction_drift.values()) <= FRACTION_DRIFT_GATE,
        "Q_density_consistency": all(
            bool(value["pass"]) for value in c["Q_density_consistency"].values()
        ),
        "rho_sanity_all_checkpoints": bool(all_checkpoint_rho_pass),
        "no_significant_averaged_backflow": not any(
            (*short["significant_backflow_by_outlet"].values(),
             *long["significant_backflow_by_outlet"].values())
        ),
        "minimum_pdf_positive": float(c["minimum_pdf"]) > 0.0,
        "maximum_lattice_speed": float(c["maximum_lattice_speed"]) < MAX_LATTICE_SPEED,
        "all_finite": bool(c["all_finite"]),
        "controller_target": target_error <= CONTROLLER_GATE,
        "controller_controlled_flux": float(controller["relative_error"]) <= CONTROLLER_GATE,
    }
    failed = [name for name, passed in gates.items() if not passed]
    return {
        "status": "PASS_NON_REFEREE" if not failed else "FAIL",
        "iteration": iterations[-1],
        "physical_time_s": iterations[-1] * spec.dt_s,
        "window_iterations": iterations,
        "nominal_short_window_iterations": spec.short_window_iterations,
        "nominal_long_window_iterations": spec.long_window_iterations,
        "actual_long_window_iterations": 2 * spec.short_window_iterations,
        "short_window": short, "long_window": long,
        "R_mass_short": short["R_mass"], "R_mass_long": long["R_mass"],
        "physical_volume_closure": float(c["physical_volume_closure"]),
        "R_velocity": r_velocity, "R_pressure": r_pressure,
        "pressure_residuals": pressure_residuals, "R_inlet": r_inlet,
        "flow_fraction_drift": fraction_drift,
        "maximum_flow_fraction_drift": max(fraction_drift.values()),
        "significant_averaged_backflow": not gates["no_significant_averaged_backflow"],
        "controller_target_expected": spec.target_lattice,
        "controller_target_observed": float(controller["target_lattice"]),
        "controller_target_error": target_error,
        "gates": gates, "failed_gates": failed,
    }


def plateau_failure(audits: Sequence[Mapping[str, Any]]) -> dict[str, Any] | None:
    if len(audits) < PLATEAU_AUDITS:
        return None
    recent = list(audits[-PLATEAU_AUDITS:])
    if float(recent[-1]["physical_time_s"]) < PLATEAU_MIN_PHYSICAL_TIME_S:
        return None
    failed_sets = [tuple(item["failed_gates"]) for item in recent]
    if not failed_sets[0] or any(value != failed_sets[0] for value in failed_sets[1:]):
        return None
    metric_names = {
        "R_mass_short": "R_mass_short", "R_mass_long": "R_mass_long",
        "physical_volume_closure": "physical_volume_closure",
        "R_velocity": "R_velocity", "R_pressure": "R_pressure",
        "R_inlet": "R_inlet", "flow_fraction_drift": "maximum_flow_fraction_drift",
        "controller_target": "controller_target_error",
    }
    changes = {}
    for gate in failed_sets[0]:
        metric = metric_names.get(gate)
        if metric is None:
            return None
        values = [float(item[metric]) for item in recent]
        changes[gate] = (max(values) - min(values)) / max(
            abs(float(np.mean(values))), np.finfo(float).tiny
        )
    if any(value >= PLATEAU_RELATIVE_CHANGE for value in changes.values()):
        return None
    return {
        "failure_mode": "SCIENTIFIC_PLATEAU_FAILURE",
        "failed_gates": list(failed_sets[0]),
        "relative_changes": changes,
        "audit_iterations": [int(item["iteration"]) for item in recent],
    }


def _archive_segment(paths: GridPaths, name: str) -> None:
    source = paths.runtime / "segments" / name
    destination = paths.segments / name
    destination.mkdir(parents=True, exist_ok=True)
    for filename in ("musubi.lua", "musubi_stdout.log", "musubi_stderr.log", "semantic_status.log"):
        if (source / filename).is_file():
            shutil.copy2(source / filename, destination / filename)


def _write_member_histories(paths: GridPaths, state: Mapping[str, Any]) -> None:
    write_json(paths.qc / "checkpoint_history.json", {
        "status": state["status"],
        "checkpoint_history": state["checkpoint_history"],
        "audits": state["audits"],
    })
    write_json(paths.qc / "physical_flux_history.json", {
        "flux_definition": FLUX_DEFINITION,
        "physical_plane_contract_sha256": PLANE_CONTRACT_SHA256,
        "samples": state["checkpoint_history"],
    })


def _launch_segment(paths: GridPaths, name: str) -> tuple[int, float]:
    launcher_wsl = f"{PROJECT_WSL}/outputs/cfd_flow/{RUN_NAME}/run_{paths.label}_segment.sh"
    started = time.perf_counter()
    process = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "/bin/bash", launcher_wsl,
         f"{paths.runtime_wsl}/segments/{name}"],
        check=False,
    )
    return process.returncode, time.perf_counter() - started


def _one_step_referee(
    root: Path, label: str, iteration: int, header: Path, binary: Path,
) -> dict[str, Any]:
    paths = GridPaths(root, label)
    name = f"full_referee_from_{iteration:07d}"
    runtime = paths.runtime / name
    restart = runtime / "restart"
    restart.mkdir(parents=True, exist_ok=True)
    lua = generate_member_lua(
        root, label, maximum_iteration=iteration + 1,
        restart_header_wsl=f"{paths.runtime_wsl}/restart/{header.name}",
        restart_first_iteration=iteration + 1, restart_interval=1,
        restart_write_wsl=f"{paths.runtime_wsl}/{name}/restart/",
    )
    (runtime / "musubi.lua").write_text(lua, encoding="utf-8", newline="\n")
    launcher_wsl = f"{PROJECT_WSL}/outputs/cfd_flow/{RUN_NAME}/run_{label}_segment.sh"
    started = time.perf_counter()
    process = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "/bin/bash", launcher_wsl,
         f"{paths.runtime_wsl}/{name}"], check=False,
    )
    elapsed = time.perf_counter() - started
    pairs = _restart_pairs(restart)
    if process.returncode != 0 or iteration + 1 not in pairs:
        raise RuntimeError(f"{label} one-step referee failed")
    cells = int(parse_mesh_header(paths.mesh)["fluid_element_count"])
    start_pdf = read_restart_pdf(binary, n_elems=cells, n_components=19)
    end_binary = pairs[iteration + 1][1]
    end_pdf = read_restart_pdf(end_binary, n_elems=cells, n_components=19)
    spec = GRID_SPECS[label]
    mesh = load_mesh_contract(
        paths.mesh,
        expected_cells=cells,
        require_runtime_order=False,
    )
    replay = replay_full_timestep(
        start_pdf, end_pdf, mesh, dx_m=spec.dx_m, dt_s=spec.dt_s,
        density_kg_m3=RHO_KG_M3,
        target_mass_flow_kg_s=TARGET_MASS_FLOW_KG_S,
        outlet_pressures_pa=spec.outlet_absolute_pressure_pa,
    )
    residual = float(replay["R_full_one_step_identity"])
    result = {
        "status": "PASS" if residual <= FULL_IDENTITY_GATE else "FAIL",
        "iteration_start": iteration, "iteration_end": iteration + 1,
        "process_wall_clock_s": elapsed, "hard_gate": FULL_IDENTITY_GATE,
        "restart_sha256": {
            str(iteration): sha256_file(binary),
            str(iteration + 1): sha256_file(end_binary),
        },
        "referee": public_step_record(replay),
    }
    write_json(run_root(root) / "qc" / f"{label}_full_referee.json", result)
    destination = paths.member / name
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(runtime, destination)
    return result


def run_member(project_root: Path, label: str) -> dict[str, Any]:
    """Run/resume one fresh grid member to confirmation or a specified stop."""

    root = Path(project_root).resolve()
    if label not in {"coarse", "fine"}:
        raise ValueError(label)
    paths = GridPaths(root, label)
    spec = GRID_SPECS[label]
    mesh_evidence = json.loads(
        (run_root(root) / "qc" / f"{label}_mesh_contract.json").read_text(encoding="utf-8")
    )
    if mesh_evidence["status"] != "PASS":
        raise RuntimeError("mesh contract must pass before Musubi")
    paths.qc.mkdir(parents=True, exist_ok=True)
    paths.segments.mkdir(parents=True, exist_ok=True)
    (paths.runtime / "restart").mkdir(parents=True, exist_ok=True)
    (paths.runtime / "segments").mkdir(parents=True, exist_ok=True)
    launcher = run_root(root) / f"run_{label}_segment.sh"
    launcher.write_text(_launcher_text(spec.short_window_iterations), encoding="utf-8", newline="\n")
    expected_contract = _runtime_contract(root, label)
    contract_path = paths.runtime / "runtime_contract.json"
    if contract_path.is_file():
        saved = json.loads(contract_path.read_text(encoding="utf-8"))
        compatible = restart_resume_contract(saved, expected_contract)
        if compatible["status"] != "PASS":
            raise RuntimeError(f"incompatible {label} restart: {compatible}")
    else:
        write_json(contract_path, expected_contract)
    state_path = paths.qc / "run_state.json"
    if state_path.is_file():
        state = json.loads(state_path.read_text(encoding="utf-8"))
    else:
        state = {
            "status": "IN_PROGRESS", "logical_simulations": 1,
            "process_launches": 0, "restart_resumes": 0,
            "operational_recoveries": 0, "solver_wall_clock_s": 0.0,
            "candidate_iteration": None, "accepted_iteration": None,
            "checkpoint_history": [], "audits": [], "safety": None,
        }
        global_state = _update_state(root)
        global_state["logical_cfd"][label] = 1
        write_json(_state_path(root), global_state)
    auditor = GridCheckpointAuditor(root, label)
    history_by_iteration = {
        int(item["iteration"]): item for item in state["checkpoint_history"]
    }

    def audit_pairs(pairs: Mapping[int, tuple[Path, Path]]) -> None:
        controllers = _controller_history(paths)
        for iteration in sorted(value for value in pairs if value not in history_by_iteration):
            candidates = [value for value in controllers if value <= iteration]
            if not candidates:
                raise RuntimeError(f"controller record absent for {label} iteration {iteration}")
            snapshot = auditor.snapshot(iteration, pairs[iteration][1], controllers[max(candidates)])
            history_by_iteration[iteration] = snapshot
            state["checkpoint_history"].append(snapshot)
            rho_ok = RHO_GATE[0] <= float(snapshot["rho_lattice"]["mean"]) <= RHO_GATE[1]
            if iteration == INITIAL_SAFETY_ITERATIONS:
                target_error = abs(
                    float(snapshot["controller"]["target_lattice"]) - spec.target_lattice
                ) / spec.target_lattice
                gates = {
                    "rho_sanity": rho_ok,
                    "controller_target": target_error <= CONTROLLER_GATE,
                    "controller_controlled_flux": float(snapshot["controller"]["relative_error"]) <= CONTROLLER_GATE,
                    "minimum_pdf_positive": float(snapshot["minimum_pdf"]) > 0.0,
                    "maximum_lattice_speed": float(snapshot["maximum_lattice_speed"]) < MAX_LATTICE_SPEED,
                    "all_finite": bool(snapshot["all_finite"]),
                }
                state["safety"] = {"status": "PASS" if all(gates.values()) else "FAIL", "iteration": iteration, "gates": gates}
                write_json(paths.qc / "initial_5000_safety.json", state["safety"])
                if not all(gates.values()):
                    state["status"] = "SAFETY_FAILURE"
                    break
            if iteration < spec.earliest_audit_iteration or iteration % spec.short_window_iterations != 0:
                continue
            needed = (
                iteration - 2 * spec.short_window_iterations,
                iteration - spec.short_window_iterations,
                iteration,
            )
            if not all(value in history_by_iteration for value in needed):
                continue
            all_rho = all(
                RHO_GATE[0] <= float(item["rho_lattice"]["mean"]) <= RHO_GATE[1]
                for item in state["checkpoint_history"]
            )
            audit = steady_window_audit(
                [history_by_iteration[value] for value in needed], spec,
                all_checkpoint_rho_pass=all_rho,
            )
            state["audits"].append(audit)
            write_json(paths.qc / "latest_steady_audit.json", audit)
            if audit["status"] == "PASS_NON_REFEREE":
                candidate = state.get("candidate_iteration")
                if candidate is not None and iteration - int(candidate) == spec.short_window_iterations:
                    state["accepted_iteration"] = iteration
                    state["status"] = "NON_REFEREE_CONFIRMED"
                    break
                state["candidate_iteration"] = iteration
            else:
                state["candidate_iteration"] = None
            plateau = plateau_failure(state["audits"])
            if plateau is not None:
                state["status"] = "SCIENTIFIC_PLATEAU_FAILURE"
                state["plateau"] = plateau
                break
        state["checkpoint_history"].sort(key=lambda item: int(item["iteration"]))
        write_json(state_path, state)
        _write_member_histories(paths, state)

    pairs = _restart_pairs(paths.runtime / "restart")
    if pairs and any(value not in history_by_iteration for value in pairs):
        audit_pairs(pairs)
    current = max(pairs) if pairs else 0
    hard_regular = spec.hard_max_iterations // spec.short_window_iterations * spec.short_window_iterations
    while state["status"] == "IN_PROGRESS" and current < hard_regular:
        if current == 0:
            maximum = INITIAL_SAFETY_ITERATIONS
            first = INITIAL_SAFETY_ITERATIONS
            interval = INITIAL_SAFETY_ITERATIONS
        elif current == INITIAL_SAFETY_ITERATIONS:
            maximum = 2 * spec.short_window_iterations
            first = spec.short_window_iterations
            interval = spec.short_window_iterations
        elif state.get("candidate_iteration") is not None:
            maximum = current + spec.short_window_iterations
            first = maximum
            interval = spec.short_window_iterations
        else:
            maximum = min(current + 2 * spec.short_window_iterations, hard_regular)
            first = current + spec.short_window_iterations
            interval = spec.short_window_iterations
        segment_name = f"segment_{current:07d}_to_{maximum:07d}"
        segment = paths.runtime / "segments" / segment_name
        segment.mkdir(parents=True, exist_ok=True)
        pairs = _restart_pairs(paths.runtime / "restart")
        restart = pairs.get(current) if current else None
        restart_header_wsl = f"{paths.runtime_wsl}/restart/{restart[0].name}" if restart else None
        lua = generate_member_lua(
            root, label, maximum_iteration=maximum,
            restart_header_wsl=restart_header_wsl,
            restart_first_iteration=first, restart_interval=interval,
        )
        contract = member_lua_contract(
            lua, root, label, maximum_iteration=maximum,
            restart_header_wsl=restart_header_wsl,
        )
        if contract["status"] != "PASS":
            raise RuntimeError(f"{label} Lua contract failed: {contract}")
        (segment / "musubi.lua").write_text(lua, encoding="utf-8", newline="\n")
        (paths.runtime / "stop").unlink(missing_ok=True)
        returncode, elapsed = _launch_segment(paths, segment_name)
        state["process_launches"] += 1
        state["restart_resumes"] += int(current > 0)
        state["solver_wall_clock_s"] += elapsed
        _archive_segment(paths, segment_name)
        pairs = _restart_pairs(paths.runtime / "restart")
        newest = max(pairs) if pairs else current
        if newest <= current:
            state["operational_recoveries"] += 1
            write_json(state_path, state)
            if state["operational_recoveries"] >= MAX_OPERATIONAL_RECOVERIES:
                state["status"] = "CFD_FLOW_OPERATIONAL_INFRASTRUCTURE_BLOCKED"
                break
            continue
        if returncode != 0:
            state["operational_recoveries"] += 1
        current = newest
        audit_pairs(pairs)
        global_state = _update_state(root)
        global_state["process_launches"] += 1
        global_state["restart_resumes"] += int(restart is not None)
        global_state["operational_recoveries"] = sum(
            int(item.get("operational_recoveries", 0))
            for item in global_state["members"].values()
        ) + int(state["operational_recoveries"])
        global_state["members"].setdefault(label, {}).update({
            "last_completed_iteration": current,
            "physical_time_s": current * spec.dt_s,
            "restart_sha256": sha256_file(pairs[current][1]),
            "failed_gates": state["audits"][-1]["failed_gates"] if state["audits"] else [],
            "operational_recoveries": state["operational_recoveries"],
        })
        global_state["next_action"] = f"CONTINUE_{label.upper()}_OR_AUDIT"
        write_json(_state_path(root), global_state)
    if state["status"] == "IN_PROGRESS":
        state["status"] = "HARD_MAX_REACHED"
    accepted_iteration = state.get("accepted_iteration")
    full: dict[str, Any] | None = None
    accepted: dict[str, Any] | None = None
    pairs = _restart_pairs(paths.runtime / "restart")
    if accepted_iteration is not None:
        header, binary = pairs[int(accepted_iteration)]
        full = _one_step_referee(root, label, int(accepted_iteration), header, binary)
        state["process_launches"] += 1
        state["solver_wall_clock_s"] += float(full["process_wall_clock_s"])
        if full["status"] != "PASS":
            state["status"] = "FULL_TIMESTEP_REFEREE_FAILURE"
        else:
            state["status"] = "PASS"
            destination = paths.member / "accepted_restart"
            destination.mkdir(exist_ok=True)
            shutil.copy2(header, destination / header.name)
            shutil.copy2(binary, destination / binary.name)
            accepted = {
                "iteration": int(accepted_iteration),
                "header": str(destination / header.name),
                "binary": str(destination / binary.name),
                "sha256": sha256_file(destination / binary.name),
            }
    write_json(state_path, state)
    _write_member_histories(paths, state)
    final_audit = state["audits"][-1] if state["audits"] else None
    history = {int(item["iteration"]): item for item in state["checkpoint_history"]}
    observables = None
    if accepted_iteration is not None:
        selected = [
            int(accepted_iteration) - 2 * spec.short_window_iterations,
            int(accepted_iteration) - spec.short_window_iterations,
            int(accepted_iteration),
        ]
        observables = _window_observables([history[value] for value in selected])
    status = (
        "PASS" if state["status"] == "PASS"
        else "CFD_FLOW_OPERATIONAL_INFRASTRUCTURE_BLOCKED"
        if state["status"] == "CFD_FLOW_OPERATIONAL_INFRASTRUCTURE_BLOCKED"
        else "CFD_FLOW_TAU1_GRID_MEMBER_STEADY_FAILED"
    )
    result = {
        "status": status, "grid": label,
        "accepted_iteration": accepted_iteration,
        "accepted_physical_time_s": int(accepted_iteration) * spec.dt_s if accepted_iteration else None,
        "accepted_restart": accepted,
        "cell_count": int(parse_mesh_header(paths.mesh)["fluid_element_count"]),
        "mesh_hashes": _runtime_contract(root, label)["mesh_hashes"],
        "runtime_contract": expected_contract,
        "steady_metrics": final_audit,
        "primary_observables": observables,
        "full_v2_residual": (
            float(full["referee"]["R_full_one_step_identity"]) if full else None
        ),
        "logical_simulations": 1, "process_launches": state["process_launches"],
        "restart_resumes": state["restart_resumes"],
        "operational_recoveries": state["operational_recoveries"],
        "solver_wall_clock_s": state["solver_wall_clock_s"],
        "plateau": state.get("plateau"),
    }
    write_json(run_root(root) / "qc" / f"{label}_steady_acceptance.json", result)
    manifest = {
        "status": "PASS" if accepted else "PRESERVED_NOT_ACCEPTED",
        "runtime_directory": paths.runtime_wsl,
        "available_complete_restarts": {
            str(value): {"header": str(pair[0]), "binary": str(pair[1]),
                         "binary_sha256": sha256_file(pair[1])}
            for value, pair in sorted(pairs.items())
        },
        "accepted_restart": accepted,
    }
    write_json(paths.qc / "restart_manifest.json", manifest)
    if status != "PASS":
        _update_state(root, stage="STOPPED", next_action=status)
        raise RuntimeError(status)
    global_state = _update_state(root)
    global_state["process_launches"] += 1
    global_state["restart_resumes"] += 1
    global_state["members"].setdefault(label, {})["steady"] = result
    global_state["next_action"] = "RUN_FINE_CFD" if label == "coarse" else "ANALYZE_CBF"
    write_json(_state_path(root), global_state)
    return result


def finalize_grid_convergence(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    qc = run_root(root) / "qc"
    coarse = json.loads((qc / "coarse_steady_acceptance.json").read_text(encoding="utf-8"))
    base = json.loads((qc / "base_readonly_acceptance.json").read_text(encoding="utf-8"))
    fine = json.loads((qc / "fine_steady_acceptance.json").read_text(encoding="utf-8"))
    records = {
        "coarse": coarse["primary_observables"],
        "base": base["primary_observables"],
        "fine": fine["primary_observables"],
    }
    for value in records.values():
        value["historical_classification"] = None
    write_json(qc / "physical_flux_cbf.json", {
        "status": "PASS", "flux_definition": FLUX_DEFINITION,
        "physical_plane_contract_sha256": PLANE_CONTRACT_SHA256,
        "extractor_revision": FLUX_ALGORITHM_REVISION,
        "grids": records,
    })
    analyses = build_primary_analyses(records)
    gate = evaluate_repaired_grid_gate(analyses)
    write_json(qc / "three_grid_scalar_analyses.json", {
        "status": gate["status"], "refinement_ratio": 1.3,
        "safety_factor": 1.25, "metrics": analyses,
    })
    csv_path = qc / "grid_convergence_primary_metrics.csv"
    fields = (
        "metric", "coarse", "base", "fine", "relative_C_B", "relative_B_F",
        "trend", "observed_order_p", "richardson_extrapolated",
        "GCI_CB_percent", "GCI_BF_percent", "pass",
    )
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for metric in PRIMARY_METRICS:
            item = analyses[metric]
            writer.writerow({
                "metric": metric, "coarse": item["coarse"], "base": item["base"],
                "fine": item["fine"], "relative_C_B": item["relative_C_B"],
                "relative_B_F": item["relative_B_F"], "trend": item["classification"],
                "observed_order_p": item["observed_order_p"],
                "richardson_extrapolated": item["richardson_extrapolated"],
                "GCI_CB_percent": item["GCI_CB_percent"],
                "GCI_BF_percent": item["GCI_BF_percent"], "pass": item["pass"],
            })
    member_pass = coarse["status"] == base["status"] == fine["status"] == "PASS"
    referee_pass = all(
        float(value) <= FULL_IDENTITY_GATE
        for value in (coarse["full_v2_residual"], base["full_v2_residual"], fine["full_v2_residual"])
    )
    mesh_pass = all(
        json.loads((qc / f"{label}_mesh_contract.json").read_text(encoding="utf-8"))["status"] == "PASS"
        for label in ("coarse", "fine")
    )
    final_pass = member_pass and referee_pass and mesh_pass and gate["status"] == "PASS"
    result = {
        **gate,
        "status": "PASS" if final_pass else "FAIL",
        "final_status": (
            "CFD_FLOW_TAU1_CBF_GRID_CONVERGENCE_PASS"
            if final_pass else "CFD_FLOW_TAU1_CBF_GRID_CONVERGENCE_FAILED"
        ),
        "actual_head_before_final_commit": _git(root, "rev-parse", "HEAD"),
        "production_pipeline_modified": _production_modified(root),
        "member_steady_pass": member_pass,
        "full_v2_pass": referee_pass,
        "mesh_contracts_pass": mesh_pass,
        "all_primary_B_to_F_le_5_percent": gate["base_fine_primary_within_5_percent"],
        "all_primary_B_to_F_le_2_percent": gate["base_fine_primary_within_preferred_2_percent"],
        "analyses": analyses,
        "WSS_status": "DEFERRED_TO_POST_GRID_PRODUCTION_VALIDATION",
        "scientific_wall_clock_s": coarse["solver_wall_clock_s"] + fine["solver_wall_clock_s"],
        "next": (
            "PROMOTE VALIDATED TAU1 CFD CONTRACT TO PRODUCTION PIPELINE AND RUN FINAL PRODUCTION REGRESSION"
            if final_pass else "STOP AND REVIEW THE SPECIFIC FAILED GRID-CONVERGENCE METRIC"
        ),
    }
    write_json(qc / "grid_convergence_final.json", result)
    _update_state(root, stage="COMPLETE", final_status=result["final_status"], next_action=result["next"])
    return result


def finalize_grid_member_failure(project_root: Path, label: str) -> dict[str, Any]:
    """Close the study without C/B/F fabrication after a member safety failure."""

    root = Path(project_root).resolve()
    if label != "fine":
        raise ValueError("the current forensic closure is specific to Fine")
    qc = run_root(root) / "qc"
    paths = GridPaths(root, label)
    spec = GRID_SPECS[label]
    member = json.loads((qc / "fine_steady_acceptance.json").read_text(encoding="utf-8"))
    safety = json.loads((paths.qc / "initial_5000_safety.json").read_text(encoding="utf-8"))
    run_state = json.loads((paths.qc / "run_state.json").read_text(encoding="utf-8"))
    sample = next(
        item for item in run_state["checkpoint_history"]
        if int(item["iteration"]) == INITIAL_SAFETY_ITERATIONS
    )
    controller = sample["controller"]
    observed = float(controller["target_lattice"])
    expected = spec.target_lattice
    count = int(controller["globBC_count"])
    implied_nonzero_samples = observed / expected * count
    source = _binary_windows(
        "/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/"
        "mus/source/bc/mus_bc_fluid_module.fpp"
    )
    forensic = {
        "status": "CONFIRMED_FINE_SAFETY_TARGET_MISMATCH",
        "classification": "FINE_ADAPTIVE_TARGET_SAMPLE_DEFICIT_AT_SAFETY_GATE",
        "iteration": INITIAL_SAFETY_ITERATIONS,
        "expected_target_lattice_formula": (
            "target_mass_flow * dt / (rho0 * dx^3)"
        ),
        "expected_target_lattice": expected,
        "observed_target_lattice": observed,
        "relative_error": abs(observed - expected) / expected,
        "observed_over_expected": observed / expected,
        "globBC_count": count,
        "implied_nonzero_mass_flow_samples": implied_nonzero_samples,
        "exact_count_ratio_hypothesis": (count - 1) / count,
        "ratio_matches_375_over_376": math.isclose(
            observed / expected, (count - 1) / count,
            rel_tol=0.0, abs_tol=1.0e-14,
        ),
        "source_accounting": {
            "path": str(source),
            "sha256": sha256_file(source),
            "sequence": [
                "massLocal = sum(massFlowRate)",
                "MPI allreduce massLocal -> massGlobal",
                "targetMassFlow = massGlobal / globBC%nElems_total",
                "targetFlux = targetMassFlow / rho0 * dtLvl / dxLvl**3",
            ],
            "inference": (
                "The printed target and count imply 375 constant-valued samples "
                "plus one zero/missing point-index evaluation before division by 376."
            ),
        },
        "other_safety_gates": {
            name: passed for name, passed in safety["gates"].items()
            if name != "controller_target"
        },
        "controlled_flux_relative_error": float(controller["relative_error"]),
        "mesh_contract_status": json.loads(
            (qc / "fine_mesh_contract.json").read_text(encoding="utf-8")
        )["status"],
        "actions": {
            "long_fine_cfd_started": False,
            "physical_parameters_changed": False,
            "controller_changed": False,
            "mesh_changed": False,
            "fine_5000_restart_preserved": True,
        },
    }
    write_json(qc / "fine_target_lattice_failure_forensic.json", forensic)
    write_json(qc / "fine_full_referee.json", {
        "status": "NOT_RUN_SAFETY_GATE_FAILED",
        "reason": "controller_target_relative_error_gt_1e-8_at_iteration_5000",
        "musubi_timesteps": 0,
    })
    coarse = json.loads((qc / "coarse_steady_acceptance.json").read_text(encoding="utf-8"))
    base = json.loads((qc / "base_readonly_acceptance.json").read_text(encoding="utf-8"))
    write_json(qc / "physical_flux_cbf.json", {
        "status": "BLOCKED_FINE_NOT_ACCEPTED",
        "flux_definition": FLUX_DEFINITION,
        "physical_plane_contract_sha256": PLANE_CONTRACT_SHA256,
        "extractor_revision": FLUX_ALGORITHM_REVISION,
        "coarse": coarse["primary_observables"],
        "base": base["primary_observables"],
        "fine_primary": None,
        "fine_5000_transient_diagnostic_only": sample,
    })
    unavailable = {
        metric: {
            "status": "UNAVAILABLE",
            "reason": "FINE_SAFETY_TARGET_GATE_FAILED_BEFORE_STEADY_ACCEPTANCE",
            "coarse": coarse["primary_observables"][metric],
            "base": base["primary_observables"][metric],
            "fine": None,
            "relative_C_B": None,
            "relative_B_F": None,
            "classification": None,
            "observed_order_p": None,
            "richardson_extrapolated": None,
            "GCI_CB_percent": None,
            "GCI_BF_percent": None,
            "pass": False,
        }
        for metric in PRIMARY_METRICS
    }
    write_json(qc / "three_grid_scalar_analyses.json", {
        "status": "UNAVAILABLE", "metrics": unavailable,
    })
    csv_path = qc / "grid_convergence_primary_metrics.csv"
    fields = (
        "metric", "coarse", "base", "fine", "relative_C_B", "relative_B_F",
        "trend", "observed_order_p", "richardson_extrapolated",
        "GCI_CB_percent", "GCI_BF_percent", "pass",
    )
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for metric, item in unavailable.items():
            writer.writerow({
                "metric": metric, "coarse": item["coarse"], "base": item["base"],
                "fine": None, "relative_C_B": None, "relative_B_F": None,
                "trend": None, "observed_order_p": None,
                "richardson_extrapolated": None, "GCI_CB_percent": None,
                "GCI_BF_percent": None, "pass": False,
            })
    state = _update_state(root)
    result = {
        "status": "FAIL",
        "final_status": "CFD_FLOW_TAU1_GRID_MEMBER_STEADY_FAILED",
        "failed_grid": label,
        "first_failed_gate": "controller_target",
        "failure_iteration": INITIAL_SAFETY_ITERATIONS,
        "failure_forensic": forensic,
        "production_pipeline_modified": _production_modified(root),
        "coarse_status": coarse["status"],
        "base_status": base["status"],
        "fine_status": member["status"],
        "three_grid_analysis": "UNAVAILABLE_FINE_NOT_ACCEPTED",
        "WSS_status": "DEFERRED_TO_POST_GRID_PRODUCTION_VALIDATION",
        "scientific_wall_clock_s": (
            float(coarse["solver_wall_clock_s"])
            + float(member["solver_wall_clock_s"])
        ),
        "seeder_calls": state["seeder_calls"],
        "logical_cfd": state["logical_cfd"],
        "process_launches": state["process_launches"],
        "restart_resumes": state["restart_resumes"],
        "operational_recoveries": state["operational_recoveries"],
        "next": "STOP AND REVIEW THE SPECIFIC FAILED GRID-CONVERGENCE METRIC",
    }
    write_json(qc / "grid_convergence_final.json", result)
    _update_state(
        root, stage="STOPPED", final_status=result["final_status"],
        next_action=result["next"],
    )
    return result


def attach_verification(project_root: Path) -> dict[str, Any]:
    """Attach final targeted-test evidence without changing the verdict."""

    root = Path(project_root).resolve()
    qc = run_root(root) / "qc"
    final_path = qc / "grid_convergence_final.json"
    result = json.loads(final_path.read_text(encoding="utf-8"))
    pytest_path = qc / "targeted_pytest.xml"
    ruff_path = qc / "targeted_ruff.json"
    verification: dict[str, Any] = {}
    if pytest_path.is_file():
        xml_root = ElementTree.parse(pytest_path).getroot()
        suite = xml_root if xml_root.tag == "testsuite" else xml_root.find("testsuite")
        if suite is not None:
            failures = int(suite.attrib.get("failures", 0))
            errors = int(suite.attrib.get("errors", 0))
            verification["targeted_pytest"] = {
                "status": "PASS" if failures == errors == 0 else "FAIL",
                "tests": int(suite.attrib.get("tests", 0)),
                "failures": failures,
                "errors": errors,
                "seconds": float(suite.attrib.get("time", 0.0)),
                "sha256": sha256_file(pytest_path),
            }
    if ruff_path.is_file():
        findings = json.loads(ruff_path.read_text(encoding="utf-8-sig"))
        verification["targeted_ruff"] = {
            "status": "PASS" if not findings else "FAIL",
            "findings": len(findings),
            "sha256": sha256_file(ruff_path),
        }
    result["verification"] = verification
    write_json(final_path, result)
    return result


def run_all(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    preflight(root)
    base_readonly_reaudit(root)
    seed_and_audit_member(root, "coarse")
    seed_and_audit_member(root, "fine")
    _update_state(root, stage="STAGE_2_COARSE_CFD", next_action="RUN_COARSE_CFD")
    run_member(root, "coarse")
    _update_state(root, stage="STAGE_3_FINE_CFD", next_action="RUN_FINE_CFD")
    run_member(root, "fine")
    _update_state(root, stage="STAGE_4_5_ANALYSIS", next_action="ANALYZE_CBF")
    return finalize_grid_convergence(root)
