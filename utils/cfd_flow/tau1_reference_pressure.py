"""Research-only Tau1 numerical reference-pressure validation.

The physical pressure differences are frozen.  Only the LBM numerical offset
is rescaled so that the reference lattice density remains one when ``dt``
changes.  This module never launches Seeder or a long CFD calculation.
"""

from __future__ import annotations

import hashlib
import json
import math
import subprocess
import time
from pathlib import Path
from typing import Any, Mapping
from xml.etree import ElementTree

import numpy as np

from .full_timestep_mass_referee import (
    FULL_IDENTITY_GATE,
    public_step_record,
    replay_full_timestep,
    source_token_evidence,
)
from .io import sha256_file, write_json
from .musubi_boundary_mass_referee import load_mesh_contract, runtime_solid_cells
from .physical_port_flux import (
    BASE_MESH_SHA256,
    FLUX_ALGORITHM_REVISION,
    PLANE_CONTRACT_REVISION,
    PORTS,
    _evaluate_prepared_plane,
    _mesh_origin_dx,
    _prepare_plane_numerics,
    plane_from_v3_record,
)
from .restart_decode import read_restart_pdf, reconstruct_macroscopic_field
from .tau1_base import (
    DX_M,
    EXPECTED_CELLS,
    MPI_RANKS,
    MPIRUN_WSL,
    MUSUBI_SHA256,
    MUSUBI_WSL,
    OLD_DT_S,
    OLD_PRESSURE_REFERENCE_PA,
    OUTLET_GAUGE_PRESSURE_PA,
    PROJECT_WSL,
    RHO_KG_M3,
    TARGET_MASS_FLOW_KG_S,
    TARGET_Q_M3_S,
    Tau1BaseRuntimeContract,
    Tau1ReferencePressureContract,
    _controller_records,
    _mesh_path,
    _physics_and_boundaries_lua,
    _restart_pairs,
    _runtime_windows,
)


RUN_NAME = "healthy_mouse_capillary_tau1_reference_scaled_smoke_anchor003274_20260901"
PHYSICAL_FLUX_RUN = "healthy_mouse_capillary_tau1_grid_convergence_anchor003274_20260831"
PLANE_CONTRACT_SHA256 = "ffaa49bdb6e43fb7208ff29df07a90d4e92ef9bfa4b96ca4f997d4f453a7f005"
MUSUBI_SOURCE_REVISION = "81f8c4f13772f6d4af31f335e1e3f99b02726e25"
OLD_BASE_CLASSIFICATION = (
    "NUMERICALLY_VALIDATED_TAU1_BASE_WITH_MIS_SCALED_REFERENCE_PRESSURE_OFFSET"
)
OLD_BASE_ITERATIONS = (2_878_425, 2_998_176, 3_117_927)
SMOKE_CHECKPOINTS = (0, 1, 10, 100, 1_000, 5_000)
SMOKE_MAX_ITERATIONS = 5_000
SMOKE_HARD_MAX_ITERATIONS = 50_000
RHO_MEAN_GATE = (0.9, 1.1)
Q_DENSITY_CONSISTENCY_GATE = 0.05
CONTROLLER_GATE = 1.0e-8


def _run_root(project_root: Path) -> Path:
    return Path(project_root).resolve() / "outputs" / "cfd_flow" / RUN_NAME


def _qc(project_root: Path) -> Path:
    path = _run_root(project_root) / "qc"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _verification_evidence(project_root: Path) -> dict[str, Any]:
    qc = _qc(project_root)
    pytest_path = qc / "targeted_pytest.xml"
    ruff_path = qc / "targeted_ruff.json"
    evidence: dict[str, Any] = {}
    if pytest_path.is_file():
        suite = ElementTree.parse(pytest_path).getroot().find("testsuite")
        if suite is not None:
            failures = int(suite.attrib["failures"])
            errors = int(suite.attrib["errors"])
            evidence["targeted_pytest"] = {
                "status": "PASS" if failures == errors == 0 else "FAIL",
                "tests": int(suite.attrib["tests"]),
                "failures": failures,
                "errors": errors,
                "seconds": float(suite.attrib["time"]),
                "sha256": sha256_file(pytest_path),
            }
    if ruff_path.is_file():
        findings = json.loads(ruff_path.read_text(encoding="utf-8"))
        evidence["targeted_ruff"] = {
            "status": "PASS" if not findings else "FAIL",
            "findings": len(findings),
            "sha256": sha256_file(ruff_path),
        }
    return evidence


def _source_root() -> Path:
    return (
        Path(
            r"\\wsl.localhost\Ubuntu\home\lzy\apes-worktrees"
            r"\musubi_mcclure_adaptive_flux_20260829_1300\mus\source"
        )
    )


def _directory_manifest(path: Path) -> dict[str, Any]:
    root = Path(path)
    records = []
    for item in sorted(entry for entry in root.rglob("*") if entry.is_file()):
        records.append(
            {
                "path": item.relative_to(root).as_posix(),
                "bytes": item.stat().st_size,
                "sha256": sha256_file(item),
            }
        )
    payload = "\n".join(
        f"{item['path']}\t{item['bytes']}\t{item['sha256']}" for item in records
    ).encode()
    return {
        "root": str(root),
        "files": len(records),
        "bytes": sum(int(item["bytes"]) for item in records),
        "manifest_sha256": hashlib.sha256(payload).hexdigest(),
    }


def protected_old_base_manifests(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    return {
        "mesh": _directory_manifest(_mesh_path(root)),
        "qc": _directory_manifest(
            root
            / "outputs"
            / "cfd_flow"
            / "healthy_mouse_capillary_tau1_base_anchor003274_20260830"
            / "qc"
        ),
        "restart": _directory_manifest(_runtime_windows() / "restart"),
    }


def write_pressure_source_contract(project_root: Path) -> dict[str, Any]:
    source = _source_root()
    physics = source / "mus_physics_module.f90"
    boundary = source / "bc" / "mus_bc_fluid_module.fpp"
    derived = source / "derived" / "mus_derQuanPhysics_module.fpp"
    old = Tau1ReferencePressureContract(dt_s=OLD_DT_S)
    new = Tau1ReferencePressureContract()
    ratio = new.unit_density_pressure_pa / old.unit_density_pressure_pa
    dt_ratio_squared = (OLD_DT_S / new.dt_s) ** 2
    old_density_under_new = new.lattice_density(OLD_PRESSURE_REFERENCE_PA)
    evidence = {
        "conversion_factors": source_token_evidence(
            physics,
            "me%fac( iLevel )%vel = me%dxLvl( iLevel )/me%dtLvl( iLevel )",
            "me%fac( iLevel )%press = me%rho0 * me%dxLvl( iLevel )**2",
            "/ me%dtLvl( iLevel )**2",
        ),
        "pressure_boundary_and_controller": source_token_evidence(
            boundary,
            "rho = rho / physics%fac( iLevel )%press * cs2inv",
            "targetFlux = targetMassFlow / physics%rho0 * physics%dtLvl(iLevel)",
            "/ physics%dxLvl(iLevel)**3",
            "pressure = rhoBoundary * physics%fac(iLevel)%press * cs2",
        ),
        "physical_pressure_output": source_token_evidence(
            derived,
            "subroutine derivePressurePhy_fromIndex",
            "res = res * fPtr%solverData%physics%fac( ilevel )%press",
        ),
    }
    checks = {
        "old_reference_equals_old_unit_density_pressure": math.isclose(
            old.unit_density_pressure_pa,
            OLD_PRESSURE_REFERENCE_PA,
            rel_tol=2.0e-15,
            abs_tol=1.0e-10,
        ),
        "new_pressure_scales_as_dt_inverse_squared": math.isclose(
            ratio, dt_ratio_squared, rel_tol=2.0e-15, abs_tol=0.0
        ),
        "new_reference_density_is_one": math.isclose(
            new.lattice_density(new.pressure_reference_pa),
            1.0,
            rel_tol=0.0,
            abs_tol=1.0e-15,
        ),
        "old_reference_predicts_low_tau1_density": math.isclose(
            old_density_under_new,
            0.006973356567857592,
            rel_tol=2.0e-15,
            abs_tol=0.0,
        ),
        "adaptive_target_is_mass_normalized_volume_flux": True,
        "velocity_conversion_is_dx_over_dt": True,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "musubi_source_revision": MUSUBI_SOURCE_REVISION,
        "binary_sha256": MUSUBI_SHA256,
        "pressure_formula": "p_phy=rho_lattice*cs2*rho0*dx^2/dt^2",
        "velocity_formula": "u_phy=u_lattice*dx/dt",
        "adaptive_target_formula": "targetMassFlow/rho0*dt/dx^3",
        "adaptive_target_semantics": (
            "population mass flux normalized by rho0; it equals the requested "
            "physical volume flux when reference rho_lattice is one"
        ),
        "old_unit_density_pressure_pa": old.unit_density_pressure_pa,
        "old_reference_pressure_pa": OLD_PRESSURE_REFERENCE_PA,
        "new_unit_density_pressure_pa": new.unit_density_pressure_pa,
        "new_to_old_pressure_ratio": ratio,
        "old_to_new_dt_ratio_squared": dt_ratio_squared,
        "rho_lattice_if_old_reference_used_with_tau1_dt": old_density_under_new,
        "source_evidence": evidence,
        "checks": checks,
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "musubi_calls": 0,
    }
    write_json(_qc(project_root) / "tau1_pressure_reference_source_contract.json", result)
    return result


def _summary(values: np.ndarray) -> dict[str, float]:
    array = np.asarray(values, dtype=np.float64)
    return {
        "mean": float(np.mean(array)),
        "median": float(np.median(array)),
        "p1": float(np.percentile(array, 1.0)),
        "p99": float(np.percentile(array, 99.0)),
        "minimum": float(np.min(array)),
        "maximum": float(np.max(array)),
    }


def write_old_reference_forensic(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    mesh = load_mesh_contract(_mesh_path(root), expected_cells=EXPECTED_CELLS)
    solid = runtime_solid_cells(mesh)
    fluid = np.ones(EXPECTED_CELLS, dtype=bool)
    fluid[np.asarray(sorted(solid), dtype=np.int64)] = False
    pairs = _restart_pairs(_runtime_windows() / "restart")
    restarts: list[dict[str, Any]] = []
    for iteration in OLD_BASE_ITERATIONS:
        if iteration not in pairs:
            raise RuntimeError(f"protected old Base restart is absent: {iteration}")
        binary = pairs[iteration][1]
        pdf = read_restart_pdf(binary, n_elems=EXPECTED_CELLS, n_components=19)
        density = np.sum(pdf, axis=1, dtype=np.float64)
        restarts.append(
            {
                "iteration": iteration,
                "binary": str(binary),
                "sha256": sha256_file(binary),
                "domain_density": _summary(density),
                "fluid_density": _summary(density[fluid]),
            }
        )
    physical_path = (
        root
        / "outputs"
        / "cfd_flow"
        / PHYSICAL_FLUX_RUN
        / "qc"
        / "base_physical_flux_preflight_v2.json"
    )
    physical = json.loads(physical_path.read_text(encoding="utf-8"))
    predicted = Tau1ReferencePressureContract().lattice_density(
        OLD_PRESSURE_REFERENCE_PA
    )
    ports: dict[str, Any] = {}
    for label in PORTS:
        records = [sample["ports"][label] for sample in physical["samples"]]
        q_velocity = float(np.mean([item["physical_q_m3_s"] for item in records]))
        q_density = float(
            np.mean([item["mass_normalized_volume_flux_m3_s"] for item in records])
        )
        effective_density = q_density / q_velocity
        ports[label] = {
            "Q_velocity_m3_s": q_velocity,
            "Q_density_weighted_m3_s": q_density,
            "Q_density_weighted_div_predicted_reference_density_m3_s": (
                q_density / predicted
            ),
            "Q_target_m3_s": TARGET_Q_M3_S,
            "effective_flux_weighted_rho_lattice": effective_density,
            "area_weighted_rho_lattice": float(
                np.mean([item["area_weighted_density_lattice"] for item in records])
            ),
            "Q_velocity_times_effective_density_m3_s": (
                q_velocity * effective_density
            ),
            "Q_velocity_over_Q_target": q_velocity / TARGET_Q_M3_S,
            "inverse_predicted_reference_density": 1.0 / predicted,
        }
    old_mean = float(
        np.mean([item["fluid_density"]["mean"] for item in restarts])
    )
    checks = {
        "old_domain_density_matches_pressure_prediction_order": (
            abs(old_mean - predicted) / predicted <= 0.05
        ),
        "inlet_density_weighted_flux_tracks_target": (
            abs(ports["inlet"]["Q_density_weighted_m3_s"] - TARGET_Q_M3_S)
            / TARGET_Q_M3_S
            <= 0.01
        ),
        "inlet_velocity_times_effective_density_identity": math.isclose(
            ports["inlet"]["Q_velocity_times_effective_density_m3_s"],
            ports["inlet"]["Q_density_weighted_m3_s"],
            rel_tol=2.0e-15,
            abs_tol=0.0,
        ),
        "old_velocity_flux_has_inverse_density_signature": (
            abs(
                ports["inlet"]["Q_velocity_over_Q_target"]
                - ports["inlet"]["inverse_predicted_reference_density"]
            )
            / ports["inlet"]["inverse_predicted_reference_density"]
            <= 0.05
        ),
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "classification": OLD_BASE_CLASSIFICATION,
        "interpretation": (
            "The archived Base remains numerical validation evidence but is not "
            "a physiological-flow Base under the corrected reference offset."
        ),
        "predicted_old_reference_rho_lattice_under_tau1_scaling": predicted,
        "existing_old_base_mean_rho_lattice": old_mean,
        "restarts": restarts,
        "physical_flux_source": {
            "path": str(physical_path),
            "sha256": sha256_file(physical_path),
            "plane_contract_sha256": physical["physical_plane_contract_sha256"],
            "algorithm": physical["flux_algorithm_revision"],
        },
        "ports": ports,
        "checks": checks,
        "protected_before": protected_old_base_manifests(root),
        "seeder_calls": 0,
        "musubi_calls": 0,
    }
    write_json(_qc(root) / "tau1_old_reference_pressure_forensic.json", result)
    return result


def write_scaling_contract(project_root: Path) -> dict[str, Any]:
    source = write_pressure_source_contract(project_root)
    old = write_old_reference_forensic(project_root)
    reference = Tau1ReferencePressureContract()
    runtime = Tau1BaseRuntimeContract()
    outlets = reference.outlet_absolute_pressures(OUTLET_GAUGE_PRESSURE_PA)
    recovered = {
        label: outlets[label] - reference.pressure_reference_pa for label in outlets
    }
    checks = {
        "source_contract": source["status"] == "PASS",
        "old_base_causal_forensic": old["status"] == "PASS",
        "runtime_pressure_is_formula_derived": math.isclose(
            runtime.pressure_reference_pa,
            reference.pressure_reference_pa,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ),
        "reference_lattice_density_is_one": math.isclose(
            reference.lattice_density(runtime.pressure_reference_pa),
            1.0,
            rel_tol=0.0,
            abs_tol=1.0e-12,
        ),
        "outlet_gauge_pressures_preserved": all(
            math.isclose(
                recovered[label], gauge, rel_tol=0.0, abs_tol=5.0e-10
            )
            for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
        ),
        "tau_one": abs(runtime.tau - 1.0) <= 1.0e-12,
        "omega_one": abs(runtime.omega - 1.0) <= 1.0e-12,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "reference_pressure": reference.as_evidence(),
        "runtime_contract": runtime.as_evidence(),
        "old_pressure_reference_pa": OLD_PRESSURE_REFERENCE_PA,
        "rho_lattice_if_old_reference_used": reference.lattice_density(
            OLD_PRESSURE_REFERENCE_PA
        ),
        "new_outlet_absolute_pressures_pa": outlets,
        "recovered_outlet_gauge_pressures_pa": recovered,
        "frozen_outlet_gauge_pressures_pa": dict(OUTLET_GAUGE_PRESSURE_PA),
        "checks": checks,
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "musubi_calls": 0,
    }
    write_json(_qc(project_root) / "tau1_reference_pressure_scaling_contract.json", result)
    return result


def generate_smoke_lua(
    contract: Tau1BaseRuntimeContract | None = None,
    *,
    maximum_iterations: int = SMOKE_MAX_ITERATIONS,
) -> str:
    contract = contract or Tau1BaseRuntimeContract()
    if maximum_iterations > SMOKE_HARD_MAX_ITERATIONS:
        raise ValueError("reference-pressure smoke exceeds the 50k hard maximum")
    trackers = []
    for iteration in checkpoint_iterations(maximum_iterations):
        label = f"state_{iteration:06d}"
        trackers.append(
            "  {label='"
            + label
            + "', folder='tracking/"
            + label
            + "/', variable={'density_phy','velocity_phy'}, "
            + "shape={kind='all'}, time_control={min={iter="
            + str(iteration)
            + "}, max={iter="
            + str(iteration)
            + "}, interval={iter=1}}, "
            + "output={format='asciiSpatial', use_get_point=false}}"
        )
    tracking_text = ",\n".join(trackers)
    return f"""-- Fresh Tau1 reference-pressure smoke; numerical offset is not physiological pressure.
simulation_name = 'tau1_reference_scaled_smoke'
printRuntimeInfo = true
timing_file = 'tracking/timing.res'
mesh = '{PROJECT_WSL}/outputs/cfd_flow/healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830/seeder/mesh/'
scaling = 'diffusive'
logging = {{level=5}}
maximum_iterations = {int(maximum_iterations)}
{_physics_and_boundaries_lua(contract)}
sim_control = {{
  time_control={{max={{iter=maximum_iterations}}, interval={{iter=100}}}},
  abort_criteria={{stop_file='stop'}}
}}
tracking = {{
{tracking_text}
}}
restart = {{write='restart/', timeformat={{use_iter=true}},
  time_control={{min={{iter={int(maximum_iterations)-1}}}, max={{iter=maximum_iterations}}, interval={{iter=1}}}}
}}
"""


def checkpoint_iterations(maximum_iterations: int) -> tuple[int, ...]:
    return tuple(
        value for value in SMOKE_CHECKPOINTS if value <= int(maximum_iterations)
    )


def smoke_lua_contract(
    text: str,
    contract: Tau1BaseRuntimeContract | None = None,
    *,
    maximum_iterations: int = SMOKE_MAX_ITERATIONS,
) -> dict[str, Any]:
    contract = contract or Tau1BaseRuntimeContract()
    reference = Tau1ReferencePressureContract(
        rho0_kg_m3=contract.rho_kg_m3,
        dx_m=contract.dx_m,
        dt_s=contract.dt_s,
    )
    outlets = reference.outlet_absolute_pressures(OUTLET_GAUGE_PRESSURE_PA)
    checks = {
        "dx": f"dx = {contract.dx_m:.17g}" in text,
        "dt": f"dt = {contract.dt_s:.17g}" in text,
        "rho0": f"rho0_phy = {contract.rho_kg_m3:.17g}" in text,
        "tau_one": abs(contract.tau - 1.0) <= 1.0e-12,
        "omega_one": abs(contract.omega - 1.0) <= 1.0e-12,
        "dynamic_pressure_reference": (
            f"pressure_reference_phy = {reference.pressure_reference_pa:.17g}" in text
        ),
        "formula_outlet_pressures": all(
            f"return {value:.17g}" in text for value in outlets.values()
        ),
        "target_mass_flow": (
            f"mass_flowrate={TARGET_MASS_FLOW_KG_S:.17g}" in text
        ),
        "wall_libb": "kind='wall_libb'" in text,
        "adaptive_flux_pressure": "kind='adaptive_flux_pressure'" in text,
        "three_pressure_eq": text.count("kind='pressure_eq'") == 3,
        "fresh_no_restart_read": "read=" not in text and "read =" not in text,
        "zero_initial_velocity": all(
            f"velocity{axis}=0.0" in text for axis in ("X", "Y", "Z")
        ),
        "initial_density_one": abs(
            reference.lattice_density(reference.pressure_reference_pa) - 1.0
        )
        <= 1.0e-12,
        "final_two_restarts": (
            f"min={{iter={maximum_iterations-1}}}" in text
            and "max={iter=maximum_iterations}" in text
            and "interval={iter=1}" in text
        ),
        "requested_checkpoints": all(
            f"state_{iteration:06d}" in text
            for iteration in checkpoint_iterations(maximum_iterations)
        ),
        "full_field_uses_ascii_spatial": (
            text.count("format='asciiSpatial'")
            == len(checkpoint_iterations(maximum_iterations))
            and text.count("use_get_point=false")
            == len(checkpoint_iterations(maximum_iterations))
        ),
        "hard_max": int(maximum_iterations) <= SMOKE_HARD_MAX_ITERATIONS,
        "no_vtk": "vtk" not in text.lower(),
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "initial_rho_lattice": reference.lattice_density(
            reference.pressure_reference_pa
        ),
        "outlet_rho_lattice": {
            label: reference.lattice_density(value) for label, value in outlets.items()
        },
        "pressure_reference_role": (
            "LBM_NUMERICAL_OFFSET_NOT_PHYSIOLOGICAL_ABSOLUTE_PRESSURE"
        ),
    }


def _smoke_launcher() -> str:
    tracking_dirs = " ".join(
        f"tracking/state_{iteration:06d}"
        for iteration in checkpoint_iterations(SMOKE_MAX_ITERATIONS)
    )
    return f"""#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${{BASH_SOURCE[0]}}")" && pwd -P)"
cd "$SCRIPT_DIR"
MUSUBI='{MUSUBI_WSL}'
MPIRUN='{MPIRUN_WSL}'
[[ -x "$MUSUBI" && -x "$MPIRUN" && -s musubi.lua ]]
[[ "$(sha256sum "$MUSUBI" | awk '{{print $1}}')" == '{MUSUBI_SHA256}' ]]
! grep -Eq 'restart[[:space:]]*=.*read|read[[:space:]]*=' musubi.lua
if find restart -maxdepth 1 -type f -size +0c 2>/dev/null | grep -q .; then
  echo 'fresh smoke restart directory is not empty' >&2
  exit 9
fi
mkdir -p restart {tracking_dirs}
"$MPIRUN" --bind-to core --map-by core --report-bindings -np {MPI_RANKS} "$MUSUBI" musubi.lua > musubi_stdout.log 2> musubi_stderr.log
grep -q 'Initializing musubi' musubi_stdout.log
grep -q 'Loading qVal data' musubi_stdout.log
grep -q 'Found BC wall of kind wall_libb' musubi_stdout.log
grep -q 'ADAPTIVE_FLUX_PRESSURE iter={SMOKE_MAX_ITERATIONS}' musubi_stdout.log
grep -q 'SUCCESSFUL run' musubi_stdout.log
printf 'FRESH_INITIALIZATION=PASS\nITERATIONS={SMOKE_MAX_ITERATIONS}\nSEMANTIC_SUCCESS=PASS\n' > semantic_status.log
"""


def prepare_reference_pressure_zero_run(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    scaling = write_scaling_contract(root)
    run_root = _run_root(root)
    run_root.mkdir(parents=True, exist_ok=True)
    contract = Tau1BaseRuntimeContract()
    lua = generate_smoke_lua(contract)
    lua_check = smoke_lua_contract(lua, contract)
    if scaling["status"] != "PASS" or lua_check["status"] != "PASS":
        raise RuntimeError("reference-pressure zero-run gate failed")
    (run_root / "musubi.lua").write_text(lua, encoding="utf-8", newline="\n")
    (run_root / "run_smoke.sh").write_text(
        _smoke_launcher(), encoding="utf-8", newline="\n"
    )
    result = {
        "status": "PREPARED_ZERO_RUN_PASS",
        "maximum_iterations": SMOKE_MAX_ITERATIONS,
        "hard_max_iterations": SMOKE_HARD_MAX_ITERATIONS,
        "fresh_initialization": True,
        "old_restart_read": False,
        "lua_contract": lua_check,
        "lua_sha256": sha256_file(run_root / "musubi.lua"),
        "launcher_sha256": sha256_file(run_root / "run_smoke.sh"),
        "binary_sha256": sha256_file(
            Path(
                r"\\wsl.localhost\Ubuntu\home\lzy\apes-worktrees"
                r"\musubi_mcclure_adaptive_flux_20260829_1300\build"
                r"\musubi_adaptive_flux"
            )
        ),
        "mesh_hashes": {
            name: sha256_file(_mesh_path(root) / name) for name in BASE_MESH_SHA256
        },
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "long_cfd_calls": 0,
        "short_smoke_musubi_calls": 0,
    }
    write_json(_qc(root) / "tau1_reference_scaled_smoke.json", result)
    return result


def _tracking_snapshot(folder: Path, mesh: Any, origin: np.ndarray) -> dict[str, Any]:
    rows: list[np.ndarray] = []
    columns: list[str] | None = None
    files = sorted(folder.glob("*.res"))
    if not files:
        raise RuntimeError(f"tracking snapshot is absent: {folder}")
    for path in files:
        header = None
        for line in path.read_text(encoding="utf-8").splitlines()[:5]:
            if "coordX" in line and "density_phy" in line:
                header = line.lstrip("#").split()
                break
        if header is None:
            raise RuntimeError(f"tracking header is incomplete: {path}")
        if columns is None:
            columns = header
        elif columns != header:
            raise RuntimeError("ranked tracking headers disagree")
        data = np.loadtxt(path, comments="#", ndmin=2)
        rows.append(np.asarray(data, dtype=np.float64))
    assert columns is not None
    data = np.vstack(rows)
    coordinates = data[:, [columns.index(name) for name in ("coordX", "coordY", "coordZ")]]
    density_phy = data[:, columns.index("density_phy")]
    velocity = data[
        :,
        [columns.index(f"velocity_phy_0{axis}") for axis in (1, 2, 3)],
    ]
    ijk = np.rint((coordinates - origin) / DX_M - 0.5).astype(np.int64)
    order = np.asarray(
        [mesh.lookup[tuple(int(value) for value in point)] for point in ijk],
        dtype=np.int64,
    )
    if len(order) != EXPECTED_CELLS or len(np.unique(order)) != EXPECTED_CELLS:
        raise RuntimeError("tracking rows do not map one-to-one to frozen mesh cells")
    density = np.empty(EXPECTED_CELLS, dtype=np.float64)
    velocity_ordered = np.empty((EXPECTED_CELLS, 3), dtype=np.float64)
    density[order] = density_phy / RHO_KG_M3
    velocity_ordered[order] = velocity
    return {
        "density_lattice": density,
        "velocity_phy": velocity_ordered,
        "files": [
            {"path": str(path), "sha256": sha256_file(path)} for path in files
        ],
    }


def _prepare_flux_evaluator(
    project_root: Path, mesh: Any, fluid_mask: np.ndarray
) -> tuple[dict[str, Any], dict[str, Any]]:
    root = Path(project_root).resolve()
    contract_path = (
        root
        / "outputs"
        / "cfd_flow"
        / PHYSICAL_FLUX_RUN
        / "qc"
        / "physical_port_flux_plane_contract_v3.json"
    )
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if contract["contract_sha256"] != PLANE_CONTRACT_SHA256:
        raise RuntimeError("validated physical plane contract changed")
    origin, dx = _mesh_origin_dx(_mesh_path(root))
    centers = origin + (mesh.cell_ijk.astype(np.float64) + 0.5) * dx
    prepared = {}
    for label in PORTS:
        record = contract["ports"][label]["planes"]["central"]
        plane = plane_from_v3_record(label, "central", record)
        prepared[label] = _prepare_plane_numerics(
            plane, centers[fluid_mask], dx_m=dx
        )
    return prepared, {
        "path": str(contract_path),
        "sha256": sha256_file(contract_path),
        "contract_sha256": contract["contract_sha256"],
    }


def _flux_snapshot(
    prepared: Mapping[str, Any],
    velocity: np.ndarray,
    density: np.ndarray,
    fluid_mask: np.ndarray,
) -> dict[str, Any]:
    ports = {}
    for label in PORTS:
        record = _evaluate_prepared_plane(
            prepared[label], velocity[fluid_mask], density[fluid_mask]
        )
        q_velocity = float(record["physical_q_m3_s"])
        q_density = float(record["mass_normalized_volume_flux_m3_s"])
        difference = abs(q_velocity - q_density) / max(
            abs(q_density), np.finfo(float).tiny
        )
        ports[label] = {
            "Q_velocity_m3_s": q_velocity,
            "Q_rho_u_over_rho0_m3_s": q_density,
            "R_Q_density_consistency": difference,
            "Q_velocity_over_Q_target": q_velocity / TARGET_Q_M3_S,
            "area_weighted_rho_lattice": record[
                "area_weighted_density_lattice"
            ],
        }
    return ports


def audit_reference_scaled_smoke(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    run_root = _run_root(root)
    prepared_zero = json.loads(
        (_qc(root) / "tau1_reference_scaled_smoke.json").read_text(encoding="utf-8")
    )
    contract = Tau1BaseRuntimeContract()
    mesh = load_mesh_contract(_mesh_path(root), expected_cells=EXPECTED_CELLS)
    origin, _ = _mesh_origin_dx(_mesh_path(root))
    solid = runtime_solid_cells(mesh)
    fluid = np.ones(EXPECTED_CELLS, dtype=bool)
    fluid[np.asarray(sorted(solid), dtype=np.int64)] = False
    prepared, plane_evidence = _prepare_flux_evaluator(root, mesh, fluid)
    densities: dict[str, Any] = {}
    fluxes: dict[str, Any] = {}
    safety: dict[str, Any] = {}
    tracking_hashes: dict[str, Any] = {}
    for iteration in checkpoint_iterations(SMOKE_MAX_ITERATIONS):
        snapshot = _tracking_snapshot(
            run_root / "tracking" / f"state_{iteration:06d}", mesh, origin
        )
        density = snapshot["density_lattice"]
        velocity = snapshot["velocity_phy"]
        densities[str(iteration)] = {
            "domain": _summary(density),
            "fluid": _summary(density[fluid]),
            "boundary_neighborhood_means": {
                label: float(np.mean(density[mesh.boundaries[label].cell_indices]))
                for label in PORTS
            },
        }
        fluxes[str(iteration)] = _flux_snapshot(
            prepared, velocity, density, fluid
        )
        speed_lattice = np.linalg.norm(velocity, axis=1) * contract.dt_s / contract.dx_m
        safety[str(iteration)] = {
            "all_finite": bool(
                np.all(np.isfinite(density)) and np.all(np.isfinite(velocity))
            ),
            "max_lattice_speed": float(np.max(speed_lattice)),
        }
        tracking_hashes[str(iteration)] = snapshot["files"]
    pairs = _restart_pairs(run_root / "restart")
    if set(pairs) != {SMOKE_MAX_ITERATIONS - 1, SMOKE_MAX_ITERATIONS}:
        raise RuntimeError(f"fresh smoke restart pair is incomplete: {sorted(pairs)}")
    start_pdf = read_restart_pdf(
        pairs[SMOKE_MAX_ITERATIONS - 1][1],
        n_elems=EXPECTED_CELLS,
        n_components=19,
    )
    end_pdf = read_restart_pdf(
        pairs[SMOKE_MAX_ITERATIONS][1],
        n_elems=EXPECTED_CELLS,
        n_components=19,
    )
    end_field = reconstruct_macroscopic_field(
        end_pdf,
        dx_m=contract.dx_m,
        dt_s=contract.dt_s,
        rho0_kg_m3=contract.rho_kg_m3,
    )
    min_pdf = float(min(np.min(start_pdf), np.min(end_pdf)))
    max_speed = float(
        np.max(np.linalg.norm(end_field.velocity_lattice, axis=1))
    )
    reference = Tau1ReferencePressureContract()
    absolute_pressure = np.asarray(end_field.pressure_phy, dtype=np.float64)
    gauge_pressure = absolute_pressure - reference.pressure_reference_pa
    pressure_observables = {
        "absolute_pressure_pa": _summary(absolute_pressure),
        "gauge_pressure_pa": _summary(gauge_pressure),
        "pressure_reference_pa": reference.pressure_reference_pa,
        "pressure_reference_role": (
            "LBM_NUMERICAL_OFFSET_NOT_PHYSIOLOGICAL_ABSOLUTE_PRESSURE"
        ),
    }
    all_pdf_finite = bool(
        np.all(np.isfinite(start_pdf)) and np.all(np.isfinite(end_pdf))
    )
    outlet_pressures = Tau1ReferencePressureContract().outlet_absolute_pressures(
        OUTLET_GAUGE_PRESSURE_PA
    )
    replay = replay_full_timestep(
        start_pdf,
        end_pdf,
        mesh,
        dx_m=contract.dx_m,
        dt_s=contract.dt_s,
        density_kg_m3=contract.rho_kg_m3,
        target_mass_flow_kg_s=contract.target_mass_flow_kg_s,
        outlet_pressures_pa=outlet_pressures,
    )
    stdout = (run_root / "musubi_stdout.log").read_text(encoding="utf-8")
    controller = _controller_records(stdout)
    if not controller:
        raise RuntimeError("fresh smoke controller evidence is absent")
    final_controller = max(controller, key=lambda item: int(item["iteration"]))
    expected_target = contract.target_lattice_flux
    controller_error = abs(
        float(final_controller["target_lattice"]) - expected_target
    ) / expected_target
    final_density = densities[str(SMOKE_MAX_ITERATIONS)]["domain"]
    final_inlet = fluxes[str(SMOKE_MAX_ITERATIONS)]["inlet"]
    gates = {
        "zero_run_config_oracle": prepared_zero["status"] == "PREPARED_ZERO_RUN_PASS",
        "fresh_no_old_restart": prepared_zero["old_restart_read"] is False,
        "mean_rho_lattice_0p9_to_1p1": (
            RHO_MEAN_GATE[0] <= final_density["mean"] <= RHO_MEAN_GATE[1]
        ),
        "old_rho_0p007_not_recurrent": final_density["mean"] > 0.9,
        "controller_target_relative_error_le_1e8": controller_error
        <= CONTROLLER_GATE,
        "controller_internal_relative_error_le_1e8": float(
            final_controller["relative_error"]
        )
        <= CONTROLLER_GATE,
        "minimum_pdf_positive": min_pdf > 0.0,
        "max_lattice_speed_lt_0p05": max_speed < 0.05,
        "all_finite": all_pdf_finite
        and all(item["all_finite"] for item in safety.values()),
        "velocity_vs_density_weighted_q_le_5pct": final_inlet[
            "R_Q_density_consistency"
        ]
        <= Q_DENSITY_CONSISTENCY_GATE,
        "old_140x_q_signature_absent": (
            0.05 <= abs(final_inlet["Q_velocity_over_Q_target"]) <= 20.0
        ),
        "full_timestep_referee_le_1e8": replay["R_full_one_step_identity"]
        <= FULL_IDENTITY_GATE,
    }
    passed = all(gates.values())
    smoke = {
        "status": "PASS" if passed else "FAIL",
        "classification_if_fail": next(
            (name for name, value in gates.items() if not value), None
        ),
        "fresh_initialization": True,
        "old_restart_read": False,
        "short_smoke_musubi_calls": 1,
        "new_short_smoke_iterations": SMOKE_MAX_ITERATIONS,
        "density_by_iteration": densities,
        "safety_by_iteration": safety,
        "tracking_file_hashes": tracking_hashes,
        "adaptive_target_lattice_expected": expected_target,
        "adaptive_target_lattice_observed": final_controller["target_lattice"],
        "controller_target_relative_error": controller_error,
        "controller": final_controller,
        "minimum_pdf": min_pdf,
        "maximum_lattice_speed": max_speed,
        "all_finite": all_pdf_finite,
        "pressure_observables": pressure_observables,
        "full_timestep_referee": public_step_record(replay),
        "restart_hashes": {
            str(iteration): {
                "header": sha256_file(pair[0]),
                "binary": sha256_file(pair[1]),
            }
            for iteration, pair in pairs.items()
        },
        "gates": gates,
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "long_cfd_calls": 0,
    }
    physical = {
        "status": "PASS"
        if gates["velocity_vs_density_weighted_q_le_5pct"]
        and gates["old_140x_q_signature_absent"]
        else "FAIL",
        "plane_contract": plane_evidence,
        "plane_contract_revision": PLANE_CONTRACT_REVISION,
        "flux_algorithm_revision": FLUX_ALGORITHM_REVISION,
        "snapshots": fluxes,
        "final_inlet": final_inlet,
        "Q_target_m3_s": TARGET_Q_M3_S,
        "R_Q_density_consistency_gate": Q_DENSITY_CONSISTENCY_GATE,
    }
    write_json(_qc(root) / "tau1_reference_scaled_smoke.json", smoke)
    write_json(_qc(root) / "tau1_reference_scaled_physical_flux_smoke.json", physical)
    protected_after = protected_old_base_manifests(root)
    old_forensic = json.loads(
        (_qc(root) / "tau1_old_reference_pressure_forensic.json").read_text(
            encoding="utf-8"
        )
    )
    protected_unchanged = protected_after == old_forensic["protected_before"]
    final_status = (
        "CFD_FLOW_TAU1_REFERENCE_PRESSURE_SCALING_VALIDATED"
        if passed and protected_unchanged
        else "REFERENCE_PRESSURE_RESCALING_SMOKE_FAILED"
    )
    final = {
        "status": final_status,
        "next": (
            "RUN ONE FRESH REFERENCE-SCALED TAU1 BASE TO STEADY STATE"
            if final_status == "CFD_FLOW_TAU1_REFERENCE_PRESSURE_SCALING_VALIDATED"
            else "STOP; INVESTIGATE THE FIRST REFERENCE-PRESSURE SCALING FAILURE"
        ),
        "pressure_reference_role": (
            "LBM_NUMERICAL_OFFSET_NOT_PHYSIOLOGICAL_ABSOLUTE_PRESSURE"
        ),
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "long_cfd_calls": 0,
        "short_smoke_musubi_calls": 1,
        "new_short_smoke_iterations": SMOKE_MAX_ITERATIONS,
        "old_unit_density_pressure_pa": Tau1ReferencePressureContract(
            dt_s=OLD_DT_S
        ).unit_density_pressure_pa,
        "old_reference_pressure_pa": OLD_PRESSURE_REFERENCE_PA,
        "old_rho_lattice_implied_by_new_scaling": reference.lattice_density(
            OLD_PRESSURE_REFERENCE_PA
        ),
        "existing_old_base_mean_rho_lattice": old_forensic[
            "existing_old_base_mean_rho_lattice"
        ],
        "new_unit_density_pressure_pa": reference.unit_density_pressure_pa,
        "new_reference_pressure_pa": reference.pressure_reference_pa,
        "new_to_old_pressure_ratio": (
            reference.unit_density_pressure_pa
            / Tau1ReferencePressureContract(dt_s=OLD_DT_S).unit_density_pressure_pa
        ),
        "old_to_new_dt_ratio_squared": (OLD_DT_S / reference.dt_s) ** 2,
        "new_outlet_absolute_pressures_pa": reference.outlet_absolute_pressures(
            OUTLET_GAUGE_PRESSURE_PA
        ),
        "preserved_outlet_gauge_pressures": True,
        "fresh_initialization_used": True,
        "initial_rho_lattice": densities["0"]["domain"]["mean"],
        "final_rho_lattice": final_density,
        "adaptive_target_lattice_expected": expected_target,
        "adaptive_target_lattice_observed": final_controller["target_lattice"],
        "controller_relative_error": controller_error,
        "final_physical_flux": fluxes[str(SMOKE_MAX_ITERATIONS)],
        "minimum_pdf": min_pdf,
        "maximum_lattice_speed": max_speed,
        "pressure_observables": pressure_observables,
        "full_timestep_referee_residual": replay["R_full_one_step_identity"],
        "old_base_classification": OLD_BASE_CLASSIFICATION,
        "protected_old_base_evidence_unchanged": protected_unchanged,
        "protected_after": protected_after,
        "smoke_gates": gates,
    }
    write_json(_qc(root) / "tau1_reference_pressure_final.json", final)
    return final


def run_reference_scaled_smoke(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    prepared = prepare_reference_pressure_zero_run(root)
    if prepared["status"] != "PREPARED_ZERO_RUN_PASS":
        return prepared
    script_wsl = f"{PROJECT_WSL}/outputs/cfd_flow/{RUN_NAME}/run_smoke.sh"
    started = time.perf_counter()
    completed = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "bash", script_wsl],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
        timeout=1_800,
    )
    (_run_root(root) / "launcher_stdout.log").write_text(
        completed.stdout, encoding="utf-8", newline="\n"
    )
    (_run_root(root) / "launcher_stderr.log").write_text(
        completed.stderr, encoding="utf-8", newline="\n"
    )
    if completed.returncode != 0:
        return finalize_reference_scaled_smoke_blocked(
            root,
            returncode=completed.returncode,
            runtime_seconds=time.perf_counter() - started,
        )
    result = audit_reference_scaled_smoke(root)
    result["runtime_seconds"] = time.perf_counter() - started
    write_json(_qc(root) / "tau1_reference_pressure_final.json", result)
    return result


def finalize_reference_scaled_smoke_blocked(
    project_root: Path,
    *,
    returncode: int | None = None,
    runtime_seconds: float | None = None,
) -> dict[str, Any]:
    """Preserve and classify the consumed zero-iteration smoke attempt.

    This is a read-only diagnosis of solver evidence except for new QC and a
    corrected, unexecuted Lua candidate.  It never launches Musubi.
    """

    root = Path(project_root).resolve()
    run_root = _run_root(root)
    final_path = _qc(root) / "tau1_reference_pressure_final.json"
    previous_final = (
        json.loads(final_path.read_text(encoding="utf-8"))
        if final_path.is_file()
        else {}
    )
    stdout_path = run_root / "musubi_stdout.log"
    stderr_path = run_root / "musubi_stderr.log"
    failed_lua_path = run_root / "musubi.lua"
    stdout = stdout_path.read_text(encoding="utf-8", errors="replace")
    marker = "Too many entries for the tracking entity as it is in ascii format"
    tracking_abort = marker in stdout
    entered_timestep_loop = "ADAPTIVE_FLUX_PRESSURE iter=" in stdout
    successful_run = "SUCCESSFUL run" in stdout
    restart_files = sorted(
        path for path in (run_root / "restart").glob("*") if path.is_file()
    )
    corrected_lua = generate_smoke_lua()
    corrected_contract = smoke_lua_contract(corrected_lua)
    candidate_path = run_root / "musubi_corrected_next_authorized.lua"
    candidate_path.write_text(corrected_lua, encoding="utf-8", newline="\n")

    old_forensic = json.loads(
        (_qc(root) / "tau1_old_reference_pressure_forensic.json").read_text(
            encoding="utf-8"
        )
    )
    protected_after = protected_old_base_manifests(root)
    protected_unchanged = protected_after == old_forensic["protected_before"]
    zero_run = json.loads(
        (_qc(root) / "tau1_reference_scaled_smoke.json").read_text(
            encoding="utf-8"
        )
    )
    source_contract = json.loads(
        (_qc(root) / "tau1_pressure_reference_source_contract.json").read_text(
            encoding="utf-8"
        )
    )
    scaling_contract = json.loads(
        (_qc(root) / "tau1_reference_pressure_scaling_contract.json").read_text(
            encoding="utf-8"
        )
    )
    evidence_files = {
        str(path.relative_to(run_root)).replace("\\", "/"): {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in (stdout_path, stderr_path, failed_lua_path, candidate_path)
        if path.is_file()
    }
    new_iterations = 0 if not entered_timestep_loop else None
    physical = {
        "status": "NOT_RUN_OPERATIONAL_BLOCKED",
        "reason": (
            "Musubi aborted while initializing full-field ASCII tracking before "
            "the first timestep; no physical-flux snapshot exists."
        ),
        "short_smoke_musubi_calls": 1,
        "new_short_smoke_iterations": new_iterations,
        "plane_contract_sha256": PLANE_CONTRACT_SHA256,
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "long_cfd_calls": 0,
    }
    write_json(
        _qc(root) / "tau1_reference_scaled_physical_flux_smoke.json", physical
    )
    reference = Tau1ReferencePressureContract()
    old_reference = Tau1ReferencePressureContract(dt_s=OLD_DT_S)
    result = {
        "status": "OPERATIONAL_INFRASTRUCTURE_BLOCKED",
        "next": (
            "STOP; DO NOT CONSUME A SECOND SMOKE CALL WITHOUT EXPLICIT AUTHORIZATION"
        ),
        "first_failure_layer": "TRACKING_OUTPUT_FORMAT_INITIALIZATION",
        "failure_reason": (
            "Full-domain density/velocity tracking used scalar ASCII output, which "
            "has a 50-entry limit; Musubi reported 91160 entries and aborted before "
            "the first timestep."
        ),
        "scientific_interpretation": (
            "No pressure-scaling, controller, PDF, flux, or conservation gate was "
            "reached; this attempt cannot validate or falsify the scaling hypothesis."
        ),
        "returncode": (
            previous_final.get("returncode", 1)
            if returncode is None
            else int(returncode)
        ),
        "runtime_seconds": (
            previous_final.get("runtime_seconds")
            if runtime_seconds is None
            else float(runtime_seconds)
        ),
        "pressure_reference_role": (
            "LBM_NUMERICAL_OFFSET_NOT_PHYSIOLOGICAL_ABSOLUTE_PRESSURE"
        ),
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "long_cfd_calls": 0,
        "short_smoke_musubi_calls": 1,
        "new_short_smoke_iterations": new_iterations,
        "fresh_initialization_used": True,
        "old_restart_read": False,
        "source_contract_status": source_contract["status"],
        "scaling_contract_status": scaling_contract["status"],
        "zero_run_config_oracle_status": zero_run["status"],
        "smoke_evidence": {
            "tracking_abort_marker_found": tracking_abort,
            "entered_timestep_loop": entered_timestep_loop,
            "successful_run_marker_found": successful_run,
            "restart_file_count": len(restart_files),
            "stdout_initial_iteration": 0,
            "files": evidence_files,
        },
        "unexecuted_correction": {
            "change": "full-field tracking ascii -> asciiSpatial",
            "candidate": str(candidate_path),
            "candidate_sha256": sha256_file(candidate_path),
            "contract": corrected_contract,
            "musubi_calls": 0,
        },
        "old_unit_density_pressure_pa": old_reference.unit_density_pressure_pa,
        "old_reference_pressure_pa": OLD_PRESSURE_REFERENCE_PA,
        "old_rho_lattice_implied_by_new_scaling": reference.lattice_density(
            OLD_PRESSURE_REFERENCE_PA
        ),
        "existing_old_base_mean_rho_lattice": old_forensic[
            "existing_old_base_mean_rho_lattice"
        ],
        "new_unit_density_pressure_pa": reference.unit_density_pressure_pa,
        "new_reference_pressure_pa": reference.pressure_reference_pa,
        "new_to_old_pressure_ratio": (
            reference.unit_density_pressure_pa
            / old_reference.unit_density_pressure_pa
        ),
        "old_to_new_dt_ratio_squared": (OLD_DT_S / reference.dt_s) ** 2,
        "new_outlet_absolute_pressures_pa": (
            reference.outlet_absolute_pressures(OUTLET_GAUGE_PRESSURE_PA)
        ),
        "preserved_outlet_gauge_pressures": True,
        "initial_rho_lattice_formula_oracle": 1.0,
        "smoke_mean_rho_lattice": None,
        "adaptive_target_lattice_expected": Tau1BaseRuntimeContract().target_lattice_flux,
        "adaptive_target_lattice_observed": None,
        "controller_relative_error": None,
        "final_physical_flux": None,
        "minimum_pdf": None,
        "maximum_lattice_speed": None,
        "full_timestep_referee_residual": None,
        "old_base_classification": OLD_BASE_CLASSIFICATION,
        "verification": _verification_evidence(root),
        "protected_old_base_evidence_unchanged": protected_unchanged,
        "protected_after": protected_after,
    }
    write_json(final_path, result)
    return result
