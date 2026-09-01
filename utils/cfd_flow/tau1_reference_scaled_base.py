"""Fresh reference-scaled Tau1 Base runner and physical steady auditor.

This research-only module orchestrates one segmented Musubi calculation.  It
does not alter the production CFD pipeline, regenerate the mesh, or tune any
physical/numerical parameter.  Full fields are reconstructed offline from the
sparse restart PDFs and evaluated with the frozen physical-plane MLS V2 code.
"""

from __future__ import annotations

import csv
import json
import math
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .full_timestep_mass_referee import (
    FULL_IDENTITY_GATE,
    public_step_record,
    replay_full_timestep,
)
from .io import sha256_file, write_json
from .musubi_boundary_mass_referee import (
    load_mesh_contract,
    replay_boundary_step,
    runtime_solid_cells,
)
from .restart_decode import read_restart_pdf, reconstruct_macroscopic_field
from .tau1_base import (
    BULK_NU_M2_S,
    DX_M,
    EXPECTED_CELLS,
    MESH_HASHES,
    MPI_RANKS,
    MPIRUN_WSL,
    MUSUBI_SHA256,
    MUSUBI_WSL,
    NU_M2_S,
    OUTLET_GAUGE_PRESSURE_PA,
    PRESSURE_REFERENCE_PA,
    PROJECT_WSL,
    RHO_KG_M3,
    TAU1_DT_S,
    TARGET_MASS_FLOW_KG_S,
    TARGET_Q_M3_S,
    Tau1BaseRuntimeContract,
    _controller_records,
    _mesh_path,
    _physics_and_boundaries_lua,
    _restart_pairs,
)
from .tau1_reference_pressure import (
    PHYSICAL_FLUX_RUN,
    PLANE_CONTRACT_SHA256,
    _flux_snapshot,
    _prepare_flux_evaluator,
)


RUN_NAME = (
    "healthy_mouse_capillary_tau1_reference_scaled_base_"
    "anchor003274_20260901"
)
RUNTIME_WSL = "/home/lzy/u3da/tau1_reference_scaled_base_20260901"
RUNTIME_WINDOWS = Path(
    r"\\wsl.localhost\Ubuntu\home\lzy\u3da"
    r"\tau1_reference_scaled_base_20260901"
)
SHORT_WINDOW_ITERATIONS = 119_751
LONG_WINDOW_ITERATIONS = 239_502
EARLIEST_AUDIT_ITERATION = 479_004
SEGMENT_ITERATIONS = 239_502
HARD_MAX_ITERATIONS = 6_131_250
LAST_REGULAR_CHECKPOINT = (
    HARD_MAX_ITERATIONS // SHORT_WINDOW_ITERATIONS * SHORT_WINDOW_ITERATIONS
)
PLATEAU_MIN_PHYSICAL_TIME_S = 0.009
PLATEAU_AUDITS = 4
MASS_GATE = 0.01
VELOCITY_GATE = 0.01
PRESSURE_GATE = 0.005
INLET_GATE = 0.01
FRACTION_DRIFT_GATE = 0.01
Q_DENSITY_GATE = 0.01
RHO_GATE = (0.9, 1.1)
CONTROLLER_GATE = 1.0e-8
MAX_LATTICE_SPEED = 0.05
BACKFLOW_FRACTION = 0.05
PLATEAU_RELATIVE_CHANGE = 0.01
PLANE_CONTRACT_RELATIVE = (
    f"outputs/cfd_flow/{PHYSICAL_FLUX_RUN}/qc/"
    "physical_port_flux_plane_contract_v3.json"
)
PORTS = ("inlet", "outlet_01", "outlet_02", "outlet_03")
OUTLETS = PORTS[1:]


@dataclass(frozen=True, slots=True)
class RunPaths:
    project_root: Path

    @property
    def run_root(self) -> Path:
        return self.project_root / "outputs" / "cfd_flow" / RUN_NAME

    @property
    def qc(self) -> Path:
        return self.run_root / "qc"

    @property
    def segments(self) -> Path:
        return self.run_root / "segments"

    @property
    def runtime(self) -> Path:
        return RUNTIME_WINDOWS


def outlet_absolute_pressures() -> dict[str, float]:
    """Derive, rather than independently pin, all absolute outlet pressures."""

    return {
        label: PRESSURE_REFERENCE_PA + float(gauge)
        for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
    }


def runtime_resume_manifest(project_root: Path) -> dict[str, Any]:
    contract = Tau1BaseRuntimeContract()
    plane_path = Path(project_root).resolve() / PLANE_CONTRACT_RELATIVE
    return {
        "mesh_hashes": dict(MESH_HASHES),
        "dx_m": DX_M,
        "dt_s": TAU1_DT_S,
        "rho0_kg_m3": RHO_KG_M3,
        "nu_m2_s": NU_M2_S,
        "bulk_nu_m2_s": BULK_NU_M2_S,
        "nu_lattice": contract.nu_lattice,
        "tau": contract.tau,
        "omega": contract.omega,
        "pressure_reference_pa": PRESSURE_REFERENCE_PA,
        "outlet_gauge_pressure_pa": dict(OUTLET_GAUGE_PRESSURE_PA),
        "outlet_absolute_pressure_pa": outlet_absolute_pressures(),
        "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
        "target_volume_flow_m3_s": TARGET_Q_M3_S,
        "layout": "d3q19",
        "collision": "bgk",
        "boundary_contract": {
            "wall": "wall_libb_continuous_q",
            "inlet": "adaptive_flux_pressure",
            "outlets": "pressure_eq",
        },
        "binary_sha256": MUSUBI_SHA256,
        "physical_plane_contract_sha256": PLANE_CONTRACT_SHA256,
        "physical_plane_file_sha256": sha256_file(plane_path),
    }


def restart_compatibility(
    saved: Mapping[str, Any], expected: Mapping[str, Any]
) -> dict[str, Any]:
    """Require every state-defining field, including dynamic P_ref."""

    checks = {key: saved.get(key) == value for key, value in expected.items()}
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def _runtime_contract(paths: RunPaths) -> dict[str, Any]:
    root = paths.project_root
    contract = Tau1BaseRuntimeContract()
    mesh = _mesh_path(root)
    mesh_hashes = {name: sha256_file(mesh / name) for name in MESH_HASHES}
    binary = Path(
        r"\\wsl.localhost\Ubuntu\home\lzy\apes-worktrees"
        r"\musubi_mcclure_adaptive_flux_20260829_1300\build"
        r"\musubi_adaptive_flux"
    )
    plane = root / PLANE_CONTRACT_RELATIVE
    smoke_referee = (
        root
        / "outputs"
        / "cfd_flow"
        / "healthy_mouse_capillary_tau1_reference_scaled_smoke_anchor003274_20260901"
        / "qc"
        / "tau1_reference_scaled_smoke.json"
    )
    smoke = json.loads(smoke_referee.read_text(encoding="utf-8"))
    manifest = runtime_resume_manifest(root)
    fresh_lua = generate_segment_lua(
        maximum_iteration=SEGMENT_ITERATIONS,
        restart_header_wsl=None,
        restart_first_iteration=SHORT_WINDOW_ITERATIONS,
    )
    fresh_lua_contract = segment_lua_contract(
        fresh_lua,
        maximum_iteration=SEGMENT_ITERATIONS,
        restart_header_wsl=None,
        restart_first_iteration=SHORT_WINDOW_ITERATIONS,
    )
    checks = {
        "mesh_files_and_hashes": mesh_hashes == MESH_HASHES,
        "mesh_cells": "nElems = 182320" in (mesh / "header.lua").read_text(
            encoding="utf-8"
        ),
        "binary_hash": binary.is_file() and sha256_file(binary) == MUSUBI_SHA256,
        "dx": DX_M == 2.0e-7,
        "dt": TAU1_DT_S == 2.038735983690112e-9,
        "rho0": RHO_KG_M3 == 1056.0,
        "viscosities": NU_M2_S == 3.27e-6 and BULK_NU_M2_S == 2.18e-6,
        "nu_lattice": math.isclose(contract.nu_lattice, 1.0 / 6.0),
        "tau": math.isclose(contract.tau, 1.0),
        "omega": math.isclose(contract.omega, 1.0),
        "dynamic_pressure_formula": math.isclose(
            PRESSURE_REFERENCE_PA,
            RHO_KG_M3 * (1.0 / 3.0) * DX_M**2 / TAU1_DT_S**2,
            rel_tol=2.0e-15,
        ),
        "initial_rho_lattice": math.isclose(
            PRESSURE_REFERENCE_PA
            / (RHO_KG_M3 * (1.0 / 3.0) * DX_M**2 / TAU1_DT_S**2),
            1.0,
        ),
        "outlet_gauges_frozen": manifest["outlet_gauge_pressure_pa"]
        == dict(OUTLET_GAUGE_PRESSURE_PA),
        "outlet_absolute_derived": all(
            math.isclose(
                manifest["outlet_absolute_pressure_pa"][label]
                - PRESSURE_REFERENCE_PA,
                gauge,
                abs_tol=5.0e-10,
            )
            for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
        ),
        "target_mass_flow": TARGET_MASS_FLOW_KG_S == 2.890180380479642e-12,
        "fresh_restart_read_none": "read='" not in fresh_lua,
        "old_restart_not_in_read_path": all(
            token not in fresh_lua for token in ("2878425", "2998176", "3117927")
        ),
        "plane_contract": (
            plane.is_file()
            and json.loads(plane.read_text(encoding="utf-8"))["contract_sha256"]
            == PLANE_CONTRACT_SHA256
        ),
        "full_timestep_v2_prevalidated": (
            smoke.get("gates", {}).get("full_timestep_referee_le_1e8") is True
        ),
        "production_pipeline_unchanged": True,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "contract": manifest,
        "checks": checks,
        "fresh_initial_condition": {
            "velocity_m_s": [0.0, 0.0, 0.0],
            "pressure_pa": PRESSURE_REFERENCE_PA,
            "rho_lattice": 1.0,
            "restart_read": None,
        },
        "fresh_segment_lua_contract": fresh_lua_contract,
        "checkpoint_contract": {
            "short_window_iterations": SHORT_WINDOW_ITERATIONS,
            "short_window_s": SHORT_WINDOW_ITERATIONS * TAU1_DT_S,
            "long_window_iterations": LONG_WINDOW_ITERATIONS,
            "long_window_s": LONG_WINDOW_ITERATIONS * TAU1_DT_S,
            "earliest_audit_iteration": EARLIEST_AUDIT_ITERATION,
            "hard_max_iterations": HARD_MAX_ITERATIONS,
            "last_regular_checkpoint": LAST_REGULAR_CHECKPOINT,
        },
        "production_pipeline_modified": False,
        "seeder_calls": 0,
    }
    if result["status"] != "PASS":
        raise RuntimeError(f"reference-scaled Base preflight failed: {checks}")
    return result


def generate_segment_lua(
    *,
    maximum_iteration: int,
    restart_header_wsl: str | None,
    restart_first_iteration: int,
    restart_interval: int = SHORT_WINDOW_ITERATIONS,
    restart_write_wsl: str = f"{RUNTIME_WSL}/restart/",
    stop_file_wsl: str = f"{RUNTIME_WSL}/stop",
) -> str:
    read = f"read='{restart_header_wsl}', " if restart_header_wsl else ""
    first = int(restart_first_iteration)
    if first > maximum_iteration or restart_interval <= 0:
        raise ValueError("invalid restart schedule")
    contract = Tau1BaseRuntimeContract()
    return f"""-- One fresh/resumable reference-scaled Tau1 Base simulation.
simulation_name = 'tau1_reference_scaled_base'
printRuntimeInfo = true
timing_file = 'tracking/timing.res'
mesh = '{PROJECT_WSL}/outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/seeder/mesh/'
scaling = 'diffusive'
logging = {{level=5}}
maximum_iterations = {maximum_iteration}
{_physics_and_boundaries_lua(contract)}
sim_control = {{
  time_control={{max={{iter=maximum_iterations}}, interval={{iter=1198}}}},
  abort_criteria={{stop_file='{stop_file_wsl}'}}
}}
restart = {{{read}write='{restart_write_wsl}',
  timeformat={{use_iter=true}},
  time_control={{min={{iter={first}}}, max={{iter=maximum_iterations}},
    interval={{iter={restart_interval}}}}}
}}
"""


def segment_lua_contract(
    text: str,
    *,
    maximum_iteration: int,
    restart_header_wsl: str | None,
    restart_first_iteration: int,
    restart_interval: int = SHORT_WINDOW_ITERATIONS,
) -> dict[str, Any]:
    contract = Tau1BaseRuntimeContract()
    checks = {
        "maximum": f"maximum_iterations = {maximum_iteration}" in text,
        "dynamic_reference": (
            f"pressure_reference_phy = {PRESSURE_REFERENCE_PA:.17g}" in text
        ),
        "initial_zero_velocity": "velocityX=0.0" in text,
        "tau_one_inputs": (
            f"dt = {TAU1_DT_S:.17g}" in text
            and f"nu_phy = {NU_M2_S:.17g}" in text
            and math.isclose(contract.tau, 1.0)
        ),
        "fresh_or_own_restart": (
            f"read='{restart_header_wsl}'" in text
            if restart_header_wsl
            else "read='" not in text
        ),
        "restart_first": f"min={{iter={restart_first_iteration}}}" in text,
        "checkpoint_interval": f"interval={{iter={restart_interval}}}" in text,
        "iteration_filenames": "timeformat={use_iter=true}" in text,
        "no_full_field_tracking": all(
            token not in text for token in ("asciiSpatial", "format='vtk'", "shape={kind='all'}")
        ),
        "boundaries": (
            "kind='wall_libb'" in text
            and "kind='adaptive_flux_pressure'" in text
            and text.count("kind='pressure_eq'") == 3
        ),
        "four_rank_external_launcher": MPI_RANKS == 4,
        "target": f"mass_flowrate={TARGET_MASS_FLOW_KG_S:.17g}" in text,
        "no_old_restart": all(
            token not in text for token in ("2878425", "2998176", "3117927")
        ),
    }
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def _launcher_script() -> str:
    return f"""#!/usr/bin/env bash
set -euo pipefail
SEGMENT_DIR="${{1:?segment directory required}}"
cd "$SEGMENT_DIR"
MUSUBI='{MUSUBI_WSL}'
MPIRUN='{MPIRUN_WSL}'
[[ -s musubi.lua && -x "$MUSUBI" && -x "$MPIRUN" ]]
[[ "$(sha256sum "$MUSUBI" | awk '{{print $1}}')" == '{MUSUBI_SHA256}' ]]
mkdir -p tracking '{RUNTIME_WSL}/restart'
set +e
"$MPIRUN" --bind-to core --map-by core --report-bindings -np {MPI_RANKS} \
  "$MUSUBI" musubi.lua 2> musubi_stderr.log | awk -v interval={SHORT_WINDOW_ITERATIONS} '
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


def prepare(project_root: Path) -> dict[str, Any]:
    paths = RunPaths(Path(project_root).resolve())
    paths.qc.mkdir(parents=True, exist_ok=True)
    paths.segments.mkdir(exist_ok=True)
    paths.runtime.mkdir(parents=True, exist_ok=True)
    (paths.runtime / "restart").mkdir(exist_ok=True)
    (paths.runtime / "segments").mkdir(exist_ok=True)
    contract = _runtime_contract(paths)
    write_json(paths.qc / "reference_scaled_base_runtime_contract.json", contract)
    expected = contract["contract"]
    runtime_manifest = paths.runtime / "runtime_contract.json"
    if runtime_manifest.is_file():
        saved = json.loads(runtime_manifest.read_text(encoding="utf-8"))
        compatible = restart_compatibility(saved, expected)
        if compatible["status"] != "PASS":
            raise RuntimeError(f"incompatible Base runtime: {compatible}")
    else:
        write_json(runtime_manifest, expected)
    launcher = paths.run_root / "run_segment.sh"
    launcher.write_text(_launcher_script(), encoding="utf-8", newline="\n")
    (paths.run_root / "fresh_initial_segment.lua").write_text(
        generate_segment_lua(
            maximum_iteration=SEGMENT_ITERATIONS,
            restart_header_wsl=None,
            restart_first_iteration=SHORT_WINDOW_ITERATIONS,
        ),
        encoding="utf-8",
        newline="\n",
    )
    return contract


class CheckpointAuditor:
    """Cache the frozen plane numerics and audit sparse restart endpoints."""

    def __init__(self, project_root: Path) -> None:
        self.root = Path(project_root).resolve()
        self.mesh = load_mesh_contract(_mesh_path(self.root), expected_cells=EXPECTED_CELLS)
        solid = runtime_solid_cells(self.mesh)
        self.fluid_mask = np.ones(EXPECTED_CELLS, dtype=bool)
        self.fluid_mask[np.asarray(sorted(solid), dtype=np.int64)] = False
        self.prepared, self.plane_evidence = _prepare_flux_evaluator(
            self.root, self.mesh, self.fluid_mask
        )

    def snapshot(
        self,
        iteration: int,
        binary: Path,
        controller: Mapping[str, Any],
    ) -> dict[str, Any]:
        pdf = np.asarray(
            read_restart_pdf(binary, n_elems=EXPECTED_CELLS, n_components=19),
            dtype=np.float64,
        )
        field = reconstruct_macroscopic_field(
            pdf, dx_m=DX_M, dt_s=TAU1_DT_S, rho0_kg_m3=RHO_KG_M3
        )
        density = np.asarray(field.density_lattice, dtype=np.float64)
        velocity = np.asarray(field.velocity_phy, dtype=np.float64)
        speed = np.linalg.norm(velocity, axis=1)
        ports = _flux_snapshot(self.prepared, velocity, density, self.fluid_mask)
        boundary = replay_boundary_step(
            pdf,
            self.mesh,
            dx_m=DX_M,
            dt_s=TAU1_DT_S,
            density_kg_m3=RHO_KG_M3,
            target_mass_flow_kg_s=TARGET_MASS_FLOW_KG_S,
            outlet_pressures_pa=outlet_absolute_pressures(),
        )
        inlet_rho = float(boundary["details"]["inlet"]["rho"])
        inlet_gauge = (inlet_rho - 1.0) * PRESSURE_REFERENCE_PA
        drops = {
            label: inlet_gauge - gauge
            for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
        }
        qin = float(ports["inlet"]["Q_velocity_m3_s"])
        qout = math.fsum(
            float(ports[label]["Q_velocity_m3_s"]) for label in OUTLETS
        )
        q_scale = max(abs(qin), np.finfo(float).tiny)
        fraction_denominator = (
            qout if abs(qout) > np.finfo(float).tiny else np.finfo(float).tiny
        )
        q_density_checks: dict[str, Any] = {}
        for label in PORTS:
            q = float(ports[label]["Q_velocity_m3_s"])
            q_rho = float(ports[label]["Q_rho_u_over_rho0_m3_s"])
            if label == "inlet" or abs(q_rho) >= 0.05 * q_scale:
                denominator = max(abs(q_rho), np.finfo(float).tiny)
                denominator_kind = "abs_Q_rho_u_over_rho0"
            else:
                denominator = q_scale
                denominator_kind = "abs_Qin_near_zero_outlet"
            residual = abs(q - q_rho) / denominator
            q_density_checks[label] = {
                "residual": residual,
                "denominator_kind": denominator_kind,
                "pass": residual <= Q_DENSITY_GATE,
            }
        return {
            "iteration": int(iteration),
            "physical_time_s": int(iteration) * TAU1_DT_S,
            "restart_sha256": sha256_file(binary),
            "rho_lattice": {
                "mean": float(np.mean(density)),
                "median": float(np.median(density)),
                "p1": float(np.percentile(density, 1.0)),
                "p99": float(np.percentile(density, 99.0)),
            },
            "mean_speed_m_s": float(np.mean(speed)),
            "inlet_gauge_pressure_pa": inlet_gauge,
            "pressure_drops_pa": drops,
            "ports": ports,
            "Qout_sum_m3_s": qout,
            "physical_volume_closure": abs(qin - qout) / q_scale,
            "flow_fractions": {
                label: float(ports[label]["Q_velocity_m3_s"])
                / fraction_denominator
                for label in OUTLETS
            },
            "Q_density_consistency": q_density_checks,
            "minimum_pdf": float(np.min(pdf)),
            "maximum_lattice_speed": float(np.max(speed) * TAU1_DT_S / DX_M),
            "all_finite": bool(
                np.all(np.isfinite(pdf))
                and np.all(np.isfinite(density))
                and np.all(np.isfinite(velocity))
            ),
            "controller": dict(controller),
        }


def _trapezoid_mean(samples: Sequence[Mapping[str, Any]], getter: Any) -> float:
    if len(samples) < 2:
        raise ValueError("a physical-time mean requires at least two samples")
    values = [float(getter(item)) for item in samples]
    times = [int(item["iteration"]) for item in samples]
    integral = math.fsum(
        0.5 * (values[index] + values[index + 1])
        * (times[index + 1] - times[index])
        for index in range(len(values) - 1)
    )
    return integral / (times[-1] - times[0])


def steady_window_audit(
    samples: Sequence[Mapping[str, Any]], *, all_checkpoint_rho_pass: bool
) -> dict[str, Any]:
    """Apply the frozen short/long physical-window acceptance definitions."""

    if len(samples) != 3:
        raise ValueError("steady audit requires long-start/short-start/end samples")
    a, b, c = samples
    iterations = [int(item["iteration"]) for item in samples]
    if (
        iterations[1] - iterations[0] != SHORT_WINDOW_ITERATIONS
        or iterations[2] - iterations[1] != SHORT_WINDOW_ITERATIONS
    ):
        raise ValueError(f"non-physical audit windows: {iterations}")

    def q(sample: Mapping[str, Any], label: str) -> float:
        return float(sample["ports"][label]["Q_velocity_m3_s"])

    def flow_metrics(window: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
        inlet = _trapezoid_mean(window, lambda item: q(item, "inlet"))
        outlets = {
            label: _trapezoid_mean(window, lambda item, label=label: q(item, label))
            for label in OUTLETS
        }
        outlet_sum = math.fsum(outlets.values())
        closure = abs(inlet - outlet_sum) / max(abs(inlet), np.finfo(float).tiny)
        backflow = {
            label: value < 0.0 and abs(value) > BACKFLOW_FRACTION * abs(inlet)
            for label, value in outlets.items()
        }
        return {
            "mean_Qin_m3_s": inlet,
            "mean_outlets_m3_s": outlets,
            "mean_Qout_sum_m3_s": outlet_sum,
            "R_mass": closure,
            "significant_backflow_by_outlet": backflow,
        }

    short = flow_metrics((b, c))
    long = flow_metrics((a, b, c))
    r_velocity = abs(float(c["mean_speed_m_s"]) - float(b["mean_speed_m_s"])) / max(
        abs(float(c["mean_speed_m_s"])), 1.0e-12
    )
    pressure_names = ("inlet_gauge_pressure_pa", *OUTLETS)

    def pressure(sample: Mapping[str, Any], name: str) -> float:
        if name == "inlet_gauge_pressure_pa":
            return float(sample[name])
        return float(sample["pressure_drops_pa"][name])

    pressure_residuals = {
        name: abs(pressure(c, name) - pressure(b, name))
        / max(abs(pressure(c, name)), 1.0)
        for name in pressure_names
    }
    r_pressure = max(pressure_residuals.values())
    r_inlet = abs(short["mean_Qin_m3_s"] - TARGET_Q_M3_S) / TARGET_Q_M3_S
    fraction_drift = {
        label: max(float(item["flow_fractions"][label]) for item in samples)
        - min(float(item["flow_fractions"][label]) for item in samples)
        for label in OUTLETS
    }
    max_fraction_drift = max(fraction_drift.values())
    controller = c["controller"]
    controller_target_error = abs(
        float(controller["target_lattice"])
        - Tau1BaseRuntimeContract().target_lattice_flux
    ) / Tau1BaseRuntimeContract().target_lattice_flux
    q_density_pass = all(
        bool(value["pass"]) for value in c["Q_density_consistency"].values()
    )
    gates = {
        "R_mass_short": short["R_mass"] <= MASS_GATE,
        "R_mass_long": long["R_mass"] <= MASS_GATE,
        "physical_volume_closure": float(c["physical_volume_closure"]) <= MASS_GATE,
        "R_velocity": r_velocity <= VELOCITY_GATE,
        "R_pressure": r_pressure <= PRESSURE_GATE,
        "R_inlet": r_inlet <= INLET_GATE,
        "flow_fraction_drift": max_fraction_drift <= FRACTION_DRIFT_GATE,
        "Q_density_consistency": q_density_pass,
        "rho_sanity_all_checkpoints": bool(all_checkpoint_rho_pass),
        "no_significant_averaged_backflow": not any(
            (*short["significant_backflow_by_outlet"].values(),
             *long["significant_backflow_by_outlet"].values())
        ),
        "minimum_pdf_positive": float(c["minimum_pdf"]) > 0.0,
        "maximum_lattice_speed": float(c["maximum_lattice_speed"]) < MAX_LATTICE_SPEED,
        "all_finite": bool(c["all_finite"]),
        "controller_target": controller_target_error <= CONTROLLER_GATE,
        "controller_controlled_flux": float(controller["relative_error"])
        <= CONTROLLER_GATE,
    }
    failed = [name for name, passed in gates.items() if not passed]
    return {
        "status": "PASS_NON_REFEREE" if not failed else "FAIL",
        "iteration": iterations[-1],
        "physical_time_s": iterations[-1] * TAU1_DT_S,
        "window_iterations": iterations,
        "short_window": short,
        "long_window": long,
        "R_mass_short": short["R_mass"],
        "R_mass_long": long["R_mass"],
        "physical_volume_closure": float(c["physical_volume_closure"]),
        "R_velocity": r_velocity,
        "R_pressure": r_pressure,
        "pressure_residuals": pressure_residuals,
        "R_inlet": r_inlet,
        "flow_fraction_drift": fraction_drift,
        "maximum_flow_fraction_drift": max_fraction_drift,
        "significant_averaged_backflow": any(
            (*short["significant_backflow_by_outlet"].values(),
             *long["significant_backflow_by_outlet"].values())
        ),
        "controller_target_expected": Tau1BaseRuntimeContract().target_lattice_flux,
        "controller_target_observed": float(controller["target_lattice"]),
        "controller_target_error": controller_target_error,
        "gates": gates,
        "failed_gates": failed,
    }


def plateau_failure(audits: Sequence[Mapping[str, Any]]) -> dict[str, Any] | None:
    """Detect the requested four-audit, one-gate physical plateau."""

    if len(audits) < PLATEAU_AUDITS:
        return None
    recent = list(audits[-PLATEAU_AUDITS:])
    if float(recent[-1]["physical_time_s"]) < PLATEAU_MIN_PHYSICAL_TIME_S:
        return None
    failed = [item.get("failed_gates", []) for item in recent]
    if not all(len(names) == 1 and names == failed[0] for names in failed):
        return None
    gate = failed[0][0]
    metric_names = {
        "R_mass_short": "R_mass_short",
        "R_mass_long": "R_mass_long",
        "physical_volume_closure": "physical_volume_closure",
        "R_velocity": "R_velocity",
        "R_pressure": "R_pressure",
        "R_inlet": "R_inlet",
        "flow_fraction_drift": "maximum_flow_fraction_drift",
        "controller_target": "controller_target_error",
    }
    metric_name = metric_names.get(gate)
    if metric_name is None:
        return None
    values = [float(item[metric_name]) for item in recent]
    relative_change = (max(values) - min(values)) / max(
        abs(float(np.mean(values))), np.finfo(float).tiny
    )
    if relative_change >= PLATEAU_RELATIVE_CHANGE:
        return None
    return {
        "failure_mode": "SCIENTIFIC_PLATEAU_FAILURE",
        "failed_gate": gate,
        "metric": metric_name,
        "values": values,
        "relative_change": relative_change,
        "audit_iterations": [int(item["iteration"]) for item in recent],
    }


def acceptance_transition(
    *, candidate_iteration: int | None, current_audit_pass: bool, iteration: int
) -> tuple[str, int | None]:
    """First PASS becomes a candidate; one further short checkpoint accepts."""

    if not current_audit_pass:
        return "CONTINUE", None
    if candidate_iteration is None:
        return "CANDIDATE", int(iteration)
    if int(iteration) - int(candidate_iteration) == SHORT_WINDOW_ITERATIONS:
        return "CONFIRMED", int(candidate_iteration)
    return "CANDIDATE", int(iteration)


def _all_controller_records(paths: RunPaths) -> dict[int, dict[str, Any]]:
    records: dict[int, dict[str, Any]] = {}
    roots = (paths.runtime / "segments", paths.segments)
    for root in roots:
        for path in sorted(root.glob("*/musubi_stdout.log")):
            for record in _controller_records(path.read_text(encoding="utf-8")):
                records[int(record["iteration"])] = record
    return records


def _write_histories(paths: RunPaths, state: Mapping[str, Any]) -> None:
    history = list(state.get("checkpoint_history", []))
    audits = list(state.get("audits", []))
    write_json(paths.qc / "reference_scaled_base_checkpoint_history.json", {
        "status": "IN_PROGRESS" if state.get("status") == "IN_PROGRESS" else state.get("status"),
        "checkpoint_history": history,
        "audits": audits,
    })
    fieldnames = (
        "iteration", "physical_time_s", "rho_mean", "rho_median", "rho_p1",
        "rho_p99", "mean_speed_m_s", "Qin_m3_s", "Q1_m3_s", "Q2_m3_s",
        "Q3_m3_s", "Qout_sum_m3_s", "physical_volume_closure",
        "inlet_gauge_pressure_pa", "minimum_pdf", "maximum_lattice_speed",
    )
    with (paths.qc / "reference_scaled_base_checkpoint_history.csv").open(
        "w", encoding="utf-8", newline=""
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for item in history:
            writer.writerow({
                "iteration": item["iteration"],
                "physical_time_s": item["physical_time_s"],
                "rho_mean": item["rho_lattice"]["mean"],
                "rho_median": item["rho_lattice"]["median"],
                "rho_p1": item["rho_lattice"]["p1"],
                "rho_p99": item["rho_lattice"]["p99"],
                "mean_speed_m_s": item["mean_speed_m_s"],
                "Qin_m3_s": item["ports"]["inlet"]["Q_velocity_m3_s"],
                "Q1_m3_s": item["ports"]["outlet_01"]["Q_velocity_m3_s"],
                "Q2_m3_s": item["ports"]["outlet_02"]["Q_velocity_m3_s"],
                "Q3_m3_s": item["ports"]["outlet_03"]["Q_velocity_m3_s"],
                "Qout_sum_m3_s": item["Qout_sum_m3_s"],
                "physical_volume_closure": item["physical_volume_closure"],
                "inlet_gauge_pressure_pa": item["inlet_gauge_pressure_pa"],
                "minimum_pdf": item["minimum_pdf"],
                "maximum_lattice_speed": item["maximum_lattice_speed"],
            })
    write_json(paths.qc / "reference_scaled_base_physical_flux_history.json", {
        "flux_definition": "PHYSICAL_INTERIOR_CROSS_SECTION_VELOCITY_FLUX",
        "density_weighted_role": "SCALING_DIAGNOSTIC_ONLY",
        "plane_contract_sha256": PLANE_CONTRACT_SHA256,
        "samples": [
            {
                "iteration": item["iteration"],
                "physical_time_s": item["physical_time_s"],
                "ports": item["ports"],
                "Qout_sum_m3_s": item["Qout_sum_m3_s"],
                "physical_volume_closure": item["physical_volume_closure"],
                "flow_fractions": item["flow_fractions"],
                "Q_density_consistency": item["Q_density_consistency"],
            }
            for item in history
        ],
    })


def _archive_segment(paths: RunPaths, name: str) -> None:
    source = paths.runtime / "segments" / name
    destination = paths.segments / name
    destination.mkdir(parents=True, exist_ok=True)
    for filename in (
        "musubi.lua", "musubi_stdout.log", "musubi_stderr.log", "semantic_status.log"
    ):
        if (source / filename).is_file():
            shutil.copy2(source / filename, destination / filename)


def _run_process(script_wsl: str, segment_wsl: str) -> tuple[int, float]:
    started = time.perf_counter()
    process = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "/bin/bash", script_wsl, segment_wsl],
        check=False,
    )
    return process.returncode, time.perf_counter() - started


def _one_step_referee(
    paths: RunPaths, iteration: int, header: Path, binary: Path
) -> dict[str, Any]:
    name = f"full_referee_from_{iteration:07d}"
    runtime = paths.runtime / name
    restart = runtime / "restart"
    restart.mkdir(parents=True, exist_ok=True)
    header_wsl = f"{RUNTIME_WSL}/restart/{header.name}"
    lua = generate_segment_lua(
        maximum_iteration=iteration + 1,
        restart_header_wsl=header_wsl,
        restart_first_iteration=iteration + 1,
        restart_interval=1,
        restart_write_wsl=f"{RUNTIME_WSL}/{name}/restart/",
        stop_file_wsl=f"{RUNTIME_WSL}/{name}/stop",
    )
    # The diagnostic is deliberately a single continuation timestep in a
    # separate directory; it never advances or mutates the primary Base.
    (runtime / "musubi.lua").write_text(lua, encoding="utf-8", newline="\n")
    launcher_wsl = f"{PROJECT_WSL}/outputs/cfd_flow/{RUN_NAME}/run_segment.sh"
    returncode, elapsed = _run_process(launcher_wsl, f"{RUNTIME_WSL}/{name}")
    if returncode != 0:
        raise RuntimeError("one-step full-timestep referee process failed")
    pairs = _restart_pairs(restart)
    if iteration + 1 not in pairs:
        raise RuntimeError("one-step full-timestep referee restart is absent")
    start_pdf = read_restart_pdf(binary, n_elems=EXPECTED_CELLS, n_components=19)
    end_binary = pairs[iteration + 1][1]
    end_pdf = read_restart_pdf(end_binary, n_elems=EXPECTED_CELLS, n_components=19)
    mesh = load_mesh_contract(_mesh_path(paths.project_root), expected_cells=EXPECTED_CELLS)
    replay = replay_full_timestep(
        start_pdf,
        end_pdf,
        mesh,
        dx_m=DX_M,
        dt_s=TAU1_DT_S,
        density_kg_m3=RHO_KG_M3,
        target_mass_flow_kg_s=TARGET_MASS_FLOW_KG_S,
        outlet_pressures_pa=outlet_absolute_pressures(),
    )
    result = {
        "status": "PASS"
        if replay["R_full_one_step_identity"] <= FULL_IDENTITY_GATE
        else "FAIL",
        "iteration_start": iteration,
        "iteration_end": iteration + 1,
        "process_wall_clock_s": elapsed,
        "hard_gate": FULL_IDENTITY_GATE,
        "restart_sha256": {
            str(iteration): sha256_file(binary),
            str(iteration + 1): sha256_file(end_binary),
        },
        "referee": public_step_record(replay),
    }
    write_json(paths.qc / "reference_scaled_base_full_referee.json", result)
    destination = paths.run_root / name
    if destination.exists():
        shutil.rmtree(destination)
    shutil.copytree(runtime, destination)
    return result


def run(project_root: Path) -> dict[str, Any]:
    """Run/resume the sole fresh Base until confirmation, plateau, or hard max."""

    paths = RunPaths(Path(project_root).resolve())
    prepare(paths.project_root)
    state_path = paths.qc / "reference_scaled_base_run_state.json"
    state: dict[str, Any]
    if state_path.is_file():
        state = json.loads(state_path.read_text(encoding="utf-8"))
    else:
        state = {
            "status": "IN_PROGRESS",
            "logical_simulations": 1,
            "process_launches": 0,
            "restart_resumes": 0,
            "operational_recoveries": 0,
            "total_solver_wall_clock_s": 0.0,
            "total_primary_iterations": 0,
            "candidate_iteration": None,
            "checkpoint_history": [],
            "audits": [],
        }
    auditor = CheckpointAuditor(paths.project_root)
    pairs = _restart_pairs(paths.runtime / "restart")
    current = max(pairs) if pairs else 0
    state["total_primary_iterations"] = max(
        int(state.get("total_primary_iterations", 0)), current
    )
    history_by_iteration = {
        int(item["iteration"]): item for item in state["checkpoint_history"]
    }
    launcher_wsl = f"{PROJECT_WSL}/outputs/cfd_flow/{RUN_NAME}/run_segment.sh"

    def audit_new_checkpoints(
        available_pairs: Mapping[int, tuple[Path, Path]],
    ) -> None:
        controllers = _all_controller_records(paths)
        new_iterations = sorted(
            iteration
            for iteration in available_pairs
            if iteration not in history_by_iteration
        )
        for iteration in new_iterations:
            controller = controllers.get(iteration)
            if controller is None:
                candidates = [key for key in controllers if key <= iteration]
                if not candidates:
                    raise RuntimeError(f"controller record absent at {iteration}")
                controller = controllers[max(candidates)]
            snapshot = auditor.snapshot(
                iteration, available_pairs[iteration][1], controller
            )
            history_by_iteration[iteration] = snapshot
            state["checkpoint_history"].append(snapshot)
            rho_mean = float(snapshot["rho_lattice"]["mean"])
            if not (RHO_GATE[0] <= rho_mean <= RHO_GATE[1]):
                state["status"] = "SCIENTIFIC_RHO_SANITY_FAILURE"
                break
            if iteration < EARLIEST_AUDIT_ITERATION:
                continue
            needed = (
                iteration - LONG_WINDOW_ITERATIONS,
                iteration - SHORT_WINDOW_ITERATIONS,
                iteration,
            )
            if not all(item in history_by_iteration for item in needed):
                continue
            all_rho_pass = all(
                RHO_GATE[0] <= float(item["rho_lattice"]["mean"]) <= RHO_GATE[1]
                for item in state["checkpoint_history"]
            )
            audit = steady_window_audit(
                [history_by_iteration[item] for item in needed],
                all_checkpoint_rho_pass=all_rho_pass,
            )
            state["audits"].append(audit)
            write_json(paths.qc / "reference_scaled_base_steady_audit.json", audit)
            transition, candidate = acceptance_transition(
                candidate_iteration=state.get("candidate_iteration"),
                current_audit_pass=audit["status"] == "PASS_NON_REFEREE",
                iteration=iteration,
            )
            state["candidate_iteration"] = candidate
            if transition == "CONFIRMED":
                full = _one_step_referee(
                    paths,
                    iteration,
                    available_pairs[iteration][0],
                    available_pairs[iteration][1],
                )
                state["process_launches"] += 1
                state["total_solver_wall_clock_s"] += full["process_wall_clock_s"]
                state["full_referee_iterations"] = 1
                if full["status"] == "PASS":
                    state["status"] = "PASS"
                    state["accepted_iteration"] = iteration
                else:
                    state["status"] = "FULL_TIMESTEP_REFEREE_FAILURE"
                break
            plateau = plateau_failure(state["audits"])
            if plateau is not None:
                state["status"] = "SCIENTIFIC_PLATEAU_FAILURE"
                state["plateau"] = plateau
                break
        state["checkpoint_history"].sort(key=lambda item: int(item["iteration"]))
        _write_histories(paths, state)
        write_json(state_path, state)

    # An interrupted Python monitor may leave solver-complete restart pairs.
    # Audit those first and never relaunch their already-computed interval.
    if pairs and any(iteration not in history_by_iteration for iteration in pairs):
        audit_new_checkpoints(pairs)

    while state["status"] == "IN_PROGRESS" and current < LAST_REGULAR_CHECKPOINT:
        maximum = min(current + SEGMENT_ITERATIONS, LAST_REGULAR_CHECKPOINT)
        if maximum - current not in {SHORT_WINDOW_ITERATIONS, SEGMENT_ITERATIONS}:
            maximum = current + SHORT_WINDOW_ITERATIONS
        segment_name = f"segment_{current:07d}_to_{maximum:07d}"
        segment_runtime = paths.runtime / "segments" / segment_name
        segment_runtime.mkdir(parents=True, exist_ok=True)
        pairs = _restart_pairs(paths.runtime / "restart")
        latest = pairs.get(current) if current else None
        restart_header_wsl = (
            f"{RUNTIME_WSL}/restart/{latest[0].name}" if latest else None
        )
        lua = generate_segment_lua(
            maximum_iteration=maximum,
            restart_header_wsl=restart_header_wsl,
            restart_first_iteration=current + SHORT_WINDOW_ITERATIONS,
        )
        lua_check = segment_lua_contract(
            lua,
            maximum_iteration=maximum,
            restart_header_wsl=restart_header_wsl,
            restart_first_iteration=current + SHORT_WINDOW_ITERATIONS,
        )
        if lua_check["status"] != "PASS":
            raise RuntimeError(f"segment Lua contract failed: {lua_check}")
        (segment_runtime / "musubi.lua").write_text(
            lua, encoding="utf-8", newline="\n"
        )
        (paths.runtime / "stop").unlink(missing_ok=True)
        returncode, elapsed = _run_process(
            launcher_wsl, f"{RUNTIME_WSL}/segments/{segment_name}"
        )
        state["process_launches"] += 1
        state["restart_resumes"] += int(current > 0)
        state["total_solver_wall_clock_s"] += elapsed
        _archive_segment(paths, segment_name)
        pairs = _restart_pairs(paths.runtime / "restart")
        newest = max(pairs) if pairs else current
        if newest <= current:
            state["operational_recoveries"] += 1
            write_json(state_path, state)
            if state["operational_recoveries"] >= 5:
                state["status"] = "CFD_FLOW_OPERATIONAL_INFRASTRUCTURE_BLOCKED"
                break
            continue
        if returncode != 0:
            # A wall-clock/launcher failure after an intact restart is a
            # recoverable process failure, not authority to recompute it.
            state["operational_recoveries"] += 1
        current = newest
        state["total_primary_iterations"] = current
        audit_new_checkpoints(pairs)
        keep = sorted(pairs)[-3:]
        for iteration, (header, binary) in pairs.items():
            if iteration not in keep:
                header.unlink(missing_ok=True)
                binary.unlink(missing_ok=True)

    if state["status"] == "IN_PROGRESS":
        state["status"] = "HARD_MAX_REACHED"
    write_json(state_path, state)
    _write_histories(paths, state)
    return finalize(paths.project_root)


def finalize(project_root: Path) -> dict[str, Any]:
    paths = RunPaths(Path(project_root).resolve())
    state = json.loads(
        (paths.qc / "reference_scaled_base_run_state.json").read_text(encoding="utf-8")
    )
    accepted_iteration = state.get("accepted_iteration")
    pairs = _restart_pairs(paths.runtime / "restart")
    accepted: dict[str, Any] | None = None
    if accepted_iteration is not None:
        header, binary = pairs[int(accepted_iteration)]
        destination = paths.run_root / "accepted_restart"
        destination.mkdir(exist_ok=True)
        shutil.copy2(header, destination / header.name)
        shutil.copy2(binary, destination / binary.name)
        accepted = {
            "iteration": int(accepted_iteration),
            "header": str(destination / header.name),
            "binary": str(destination / binary.name),
            "sha256": sha256_file(destination / binary.name),
        }
    manifest = {
        "status": "PASS" if accepted else "PRESERVED_NOT_ACCEPTED",
        "runtime_directory": RUNTIME_WSL,
        "available_complete_restarts": {
            str(iteration): {
                "header": str(header),
                "header_sha256": sha256_file(header),
                "binary": str(binary),
                "binary_sha256": sha256_file(binary),
            }
            for iteration, (header, binary) in sorted(pairs.items())
        },
        "accepted_restart": accepted,
    }
    write_json(paths.qc / "reference_scaled_base_restart_manifest.json", manifest)
    audits = state.get("audits", [])
    final_audit = audits[-1] if audits else None
    status = (
        "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_PASS"
        if state["status"] == "PASS"
        else "CFD_FLOW_OPERATIONAL_INFRASTRUCTURE_BLOCKED"
        if state["status"] == "CFD_FLOW_OPERATIONAL_INFRASTRUCTURE_BLOCKED"
        else "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_FAILED"
    )
    final = {
        "status": status,
        "failure_mode": (
            "SCIENTIFIC_PLATEAU_FAILURE"
            if state["status"] == "SCIENTIFIC_PLATEAU_FAILURE"
            else None
        ),
        "next": (
            "RUN REPAIRED TAU1 COARSE/BASE/FINE GRID CONVERGENCE"
            if status == "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_PASS"
            else "STOP"
        ),
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "fresh_base_logical_simulations": 1,
        "musubi_process_launches": state["process_launches"],
        "restart_resumes": state["restart_resumes"],
        "operational_recoveries": state["operational_recoveries"],
        "total_primary_iterations": state["total_primary_iterations"],
        "total_referee_iterations": state.get("full_referee_iterations", 0),
        "total_solver_wall_clock_s": state["total_solver_wall_clock_s"],
        "accepted_restart": accepted,
        "steady_audit": final_audit,
        "plateau": state.get("plateau"),
        "runtime_contract": runtime_resume_manifest(paths.project_root),
        "first_failed_scientific_gate": (
            final_audit.get("failed_gates", [None])[0]
            if final_audit and final_audit.get("failed_gates")
            else None
        ),
    }
    write_json(paths.qc / "reference_scaled_base_final.json", final)
    return final
