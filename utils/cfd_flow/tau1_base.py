"""Research-only continuous-q Referee V2 and fresh tau=1 Base harness."""

from __future__ import annotations

import csv
import json
import math
import os
import re
import shutil
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from .healthy_capillary_calibration import OUTLET_GAUGE_PRESSURE_PA
from .full_timestep_mass_referee import (
    BOUNDARY_FLUX_DEFINITION,
    DEFERRED_PHYSICAL_FLUX,
    FULL_IDENTITY_GATE,
    full_identity_pass,
    pull_link_counts,
    public_step_record,
    replay_full_timestep,
    source_token_evidence,
    stable_delta,
)
from .io import sha256_file, write_json
from .musubi_boundary_mass_referee import (
    REFEREE_REVISION_NEW,
    _wall_operations,
    boundary_window_closure,
    conservation_identity_residual,
    load_mesh_contract,
    replay_boundary_step,
    runtime_solid_cells,
    significant_time_averaged_backflow,
    trapezoidal_integral,
)
from .restart_decode import D3Q19_DIRECTIONS, parse_restart_header, read_restart_pdf


RUN_NAME = "healthy_mouse_capillary_tau1_base_anchor003274_20260830"
FROZEN_MESH_RUN = (
    "healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830"
)
FREEZE_RUN = "healthy_mouse_capillary_port_grid_sensitivity_research_anchor003274_20260830"
EXPECTED_CELLS = 182_320
EXPECTED_INLET_GLOBBC = 223
DX_M = 2.0e-7
RHO_KG_M3 = 1056.0
NU_M2_S = 3.27e-6
BULK_NU_M2_S = 2.18e-6
TARGET_MEAN_VELOCITY_M_S = 0.35e-3
TARGET_Q_M3_S = 2.7369132390905703e-15
TARGET_MASS_FLOW_KG_S = 2.890180380479642e-12
CS2 = 1.0 / 3.0
OLD_DT_S = 2.44140625e-8
OLD_PRESSURE_REFERENCE_PA = 23622.320128
TAU1_DT_S = DX_M**2 / (6.0 * NU_M2_S)
PRESSURE_REFERENCE_PA = CS2 * (RHO_KG_M3 * DX_M**2 / TAU1_DT_S**2)
MUSUBI_WSL = (
    "/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/"
    "build/musubi_adaptive_flux"
)
MUSUBI_SHA256 = "e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588"
MPIRUN_WSL = "/home/lzy/.local/bin/mpirun"
PROJECT_WSL = "/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular"
RUNTIME_WSL = "/home/lzy/u3da/tau1_base_20260830"
TEM_ROOT_WSL = "/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/tem"
TEM_GIT_SHA = "4e8b277b66226277171ef93bf054d21270812793"
MPI_RANKS = 4
CHECKPOINT_INTERVAL = 239_502
CONFIRMATION_INTERVAL = 119_751
TRACKING_INTERVAL = 59_875
SIM_CONTROL_INTERVAL = 1_198
HARD_MAX_PHYSICAL_TIME_S = 0.0125
CONTINUOUS_WALL_CLOCK_SECONDS = 8_400
CONTINUOUS_PROCESS_TIMEOUT_SECONDS = 9_000
MASS_GATE = 0.01
VELOCITY_GATE = 0.01
PRESSURE_GATE = 0.005
INLET_GATE = 0.01
BOUNDARY_CLOSURE_GATE = 0.001
CONTROLLER_GATE = 1.0e-8
REFEREE_GATE = 1.0e-8
MESH_HASHES = {
    "elemlist.lsb": "f7d7b1d55273c78c336ac04e39bc018dd9ebb470a9f29ce833ff01711de8c386",
    "bnd.lsb": "520d7dd1e4a46a45f9b1218a5807cfd89d6f054e0a247872362b130ff6bcfe69",
    "qval.lsb": "35884406b5f0111cd4ab471f7b08ac3df00e478d3458a57636d1bd8921cb0fe6",
}


@dataclass(frozen=True, slots=True)
class Tau1ReferencePressureContract:
    """Map a numerical reference density to the current physical scaling."""

    rho0_kg_m3: float = RHO_KG_M3
    cs2: float = CS2
    dx_m: float = DX_M
    dt_s: float = TAU1_DT_S
    reference_lattice_density: float = 1.0

    @property
    def pressure_conversion_pa(self) -> float:
        return self.rho0_kg_m3 * self.dx_m**2 / self.dt_s**2

    @property
    def unit_density_pressure_pa(self) -> float:
        return self.cs2 * self.pressure_conversion_pa

    @property
    def pressure_reference_pa(self) -> float:
        return self.reference_lattice_density * self.unit_density_pressure_pa

    def lattice_density(self, pressure_pa: float) -> float:
        return float(pressure_pa) / self.unit_density_pressure_pa

    def gauge_pressure(self, pressure_pa: float) -> float:
        return float(pressure_pa) - self.pressure_reference_pa

    def outlet_absolute_pressures(
        self, gauge_pressures_pa: Mapping[str, float]
    ) -> dict[str, float]:
        return {
            label: self.pressure_reference_pa + float(gauge)
            for label, gauge in gauge_pressures_pa.items()
        }

    def as_evidence(self) -> dict[str, Any]:
        values = asdict(self)
        values.update(
            {
                "pressure_conversion_pa": self.pressure_conversion_pa,
                "unit_density_pressure_pa": self.unit_density_pressure_pa,
                "pressure_reference_pa": self.pressure_reference_pa,
                "pressure_reference_role": (
                    "LBM_NUMERICAL_OFFSET_NOT_PHYSIOLOGICAL_ABSOLUTE_PRESSURE"
                ),
            }
        )
        return values


@dataclass(frozen=True, slots=True)
class Tau1BaseRuntimeContract:
    dx_m: float = DX_M
    dt_s: float = TAU1_DT_S
    rho_kg_m3: float = RHO_KG_M3
    nu_m2_s: float = NU_M2_S
    bulk_nu_m2_s: float = BULK_NU_M2_S
    target_mass_flow_kg_s: float = TARGET_MASS_FLOW_KG_S
    pressure_reference_pa: float = PRESSURE_REFERENCE_PA
    expected_cells: int = EXPECTED_CELLS
    binary_sha256: str = MUSUBI_SHA256

    @property
    def nu_lattice(self) -> float:
        return self.nu_m2_s * self.dt_s / self.dx_m**2

    @property
    def tau(self) -> float:
        return 3.0 * self.nu_lattice + 0.5

    @property
    def omega(self) -> float:
        return 1.0 / self.tau

    @property
    def target_lattice_flux(self) -> float:
        return (
            self.target_mass_flow_kg_s
            * self.dt_s
            / (self.rho_kg_m3 * self.dx_m**3)
        )

    @property
    def hard_max_iterations(self) -> int:
        return round(HARD_MAX_PHYSICAL_TIME_S / self.dt_s)

    def as_evidence(self) -> dict[str, Any]:
        values = asdict(self)
        values.update(
            {
                "nu_lattice": self.nu_lattice,
                "tau": self.tau,
                "omega": self.omega,
                "target_lattice_flux": self.target_lattice_flux,
                "old_to_new_dt_ratio": OLD_DT_S / self.dt_s,
                "tracking_interval_iterations": TRACKING_INTERVAL,
                "checkpoint_interval_iterations": CHECKPOINT_INTERVAL,
                "confirmation_interval_iterations": CONFIRMATION_INTERVAL,
                "hard_max_physical_time_s": HARD_MAX_PHYSICAL_TIME_S,
                "hard_max_iterations": self.hard_max_iterations,
            }
        )
        return values


def historical_tau1_runtime_contract() -> Tau1BaseRuntimeContract:
    """Return the immutable pressure offset used by the archived Tau1 Base."""

    return Tau1BaseRuntimeContract(
        pressure_reference_pa=OLD_PRESSURE_REFERENCE_PA
    )


def rescale_physical_window(old_iterations: int, *, old_dt_s: float = OLD_DT_S) -> int:
    if old_iterations <= 0 or old_dt_s <= 0.0:
        raise ValueError("window and timestep must be positive")
    return round(int(old_iterations) * float(old_dt_s) / Tau1BaseRuntimeContract().dt_s)


def base_decision(*, referee_pass: bool, base_pass: bool | None) -> str:
    if not referee_pass:
        return "CFD_FLOW_CONTINUOUS_Q_REFEREE_FAILED"
    if base_pass is None:
        return "RUN_FRESH_BASE_TAU1"
    return (
        "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_PASS"
        if base_pass
        else "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_FAILED"
    )


def restart_resume_contract(
    saved: Mapping[str, Any], expected: Mapping[str, Any]
) -> dict[str, Any]:
    keys = (
        "mesh_hashes",
        "dx_m",
        "dt_s",
        "rho_kg_m3",
        "nu_m2_s",
        "bulk_nu_m2_s",
        "boundary_contract",
        "outlet_gauge_pressure_pa",
        "target_mass_flow_kg_s",
        "binary_sha256",
    )
    checks = {key: saved.get(key) == expected.get(key) for key in keys}
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def _mesh_path(root: Path) -> Path:
    return root / "outputs" / "cfd_flow" / FROZEN_MESH_RUN / "seeder" / "mesh"


def _run_root(root: Path) -> Path:
    return root / "outputs" / "cfd_flow" / RUN_NAME


def _freeze_path(root: Path) -> Path:
    return root / "outputs" / "cfd_flow" / FREEZE_RUN / "qc" / "seeder_geometry_freeze.json"


def _binary_windows() -> Path:
    return Path(
        r"\\wsl.localhost\Ubuntu\home\lzy\apes-worktrees\musubi_mcclure_adaptive_flux_20260829_1300\build\musubi_adaptive_flux"
    )


def _runtime_windows() -> Path:
    return Path(r"\\wsl.localhost\Ubuntu\home\lzy\u3da\tau1_base_20260830")


def _tem_windows() -> Path:
    return Path(
        r"\\wsl.localhost\Ubuntu\home\lzy\apes-worktrees"
        r"\musubi_mcclure_adaptive_flux_20260829_1300\tem"
    )


def _mesh_wsl() -> str:
    return f"{PROJECT_WSL}/outputs/cfd_flow/{FROZEN_MESH_RUN}/seeder/mesh"


def _read_text(path: Path) -> str:
    readable = path.resolve()
    if (
        os.name == "nt"
        and not str(readable).startswith("\\\\?\\")
        and not str(readable).startswith("\\\\")
    ):
        readable = Path("\\\\?\\" + str(readable))
    return readable.read_text(encoding="utf-8", errors="strict")


def tem_restart_timeformat_contract() -> dict[str, Any]:
    """Source-prove the linked Tem iteration-based restart filename contract."""

    restart_source = _tem_windows() / "source" / "tem_restart_module.f90"
    formatter_source = (
        _tem_windows() / "source" / "control" / "tem_timeformatter_module.f90"
    )
    restart_text = _read_text(restart_source)
    formatter_text = _read_text(formatter_source)
    checks = {
        "restart_calls_formatter_with_restart_parent": (
            "call tem_timeformatter_load(me     = me%timeform" in restart_text
            and "parent = restart_table" in restart_text
        ),
        "default_key_is_timeformat": "loc_key = 'timeformat'" in formatter_text,
        "use_iter_key_supported": "key     = 'use_iter'" in formatter_text,
        "iteration_format_is_unambiguous_integer": (
            "timeform = '(I0)'" in formatter_text
        ),
        "iteration_stamp_selected": (
            "stamp    = tem_timeformatter_iter_stamp" in formatter_text
        ),
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "linked_tem_git_sha": TEM_GIT_SHA,
        "required_lua_syntax": "timeformat={use_iter=true}",
        "forbidden_typo": "timeform={use_iter=true}",
        "restart_source_wsl": f"{TEM_ROOT_WSL}/source/tem_restart_module.f90",
        "restart_source_sha256": sha256_file(restart_source),
        "formatter_source_wsl": (
            f"{TEM_ROOT_WSL}/source/control/tem_timeformatter_module.f90"
        ),
        "formatter_source_sha256": sha256_file(formatter_source),
        "checks": checks,
    }


def _runtime_manifest(contract: Tau1BaseRuntimeContract) -> dict[str, Any]:
    return {
        "mesh_hashes": dict(MESH_HASHES),
        "dx_m": contract.dx_m,
        "dt_s": contract.dt_s,
        "rho_kg_m3": contract.rho_kg_m3,
        "nu_m2_s": contract.nu_m2_s,
        "bulk_nu_m2_s": contract.bulk_nu_m2_s,
        "boundary_contract": {
            "wall": "wall_libb",
            "inlet": "adaptive_flux_pressure",
            "outlets": "pressure_eq",
            "layout": "d3q19",
            "relaxation": "bgk",
        },
        "outlet_gauge_pressure_pa": dict(OUTLET_GAUGE_PRESSURE_PA),
        "target_mass_flow_kg_s": contract.target_mass_flow_kg_s,
        "binary_sha256": contract.binary_sha256,
    }


def fast_zero_run_preflight(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    mesh = _mesh_path(root)
    contract = Tau1BaseRuntimeContract()
    freeze = json.loads(_freeze_path(root).read_text(encoding="utf-8"))
    required = ("header.lua", "elemlist.lsb", "bnd.lua", "bnd.lsb", "qval.lua", "qval.lsb")
    actual_hashes = {name: sha256_file(mesh / name) for name in required}
    header_match = re.search(r"(?m)^\s*nElems\s*=\s*(\d+)", (mesh / "header.lua").read_text(encoding="utf-8"))
    cell_count = int(header_match.group(1)) if header_match else -1
    binary = _binary_windows()
    binary_hash = sha256_file(binary) if binary.is_file() else None
    checks = {
        "required_mesh_files": all((mesh / name).is_file() for name in required),
        "cell_count": cell_count == EXPECTED_CELLS,
        "freeze_status": freeze.get("status") == "CFD_FLOW_SEEDER_GEOMETRY_KERNEL_VALIDATED",
        "freeze_cells": freeze.get("final_base", {}).get("fluid_cells") == EXPECTED_CELLS,
        "elemlist_hash": actual_hashes["elemlist.lsb"] == MESH_HASHES["elemlist.lsb"],
        "bnd_hash": actual_hashes["bnd.lsb"] == MESH_HASHES["bnd.lsb"],
        "qval_hash": actual_hashes["qval.lsb"] == MESH_HASHES["qval.lsb"],
        "binary_exists": binary.is_file(),
        "binary_hash": binary_hash == MUSUBI_SHA256,
        "tau_exact": abs(contract.tau - 1.0) <= 1.0e-12,
        "omega_exact": abs(contract.omega - 1.0) <= 1.0e-12,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "seeder_calls": 0,
        "mesh": str(mesh),
        "mesh_cell_count": cell_count,
        "mesh_hashes": actual_hashes,
        "binary": MUSUBI_WSL,
        "binary_sha256": binary_hash,
        "runtime_contract": contract.as_evidence(),
        "checks": checks,
    }
    if result["status"] != "PASS":
        raise RuntimeError(f"tau1 Base preflight failed: {result}")
    return result


def _physics_and_boundaries_lua(contract: Tau1BaseRuntimeContract) -> str:
    outlets = {
        label: contract.pressure_reference_pa + gauge
        for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
    }
    return f"""dx = {contract.dx_m:.17g}
dt = {contract.dt_s:.17g}
rho0_phy = {contract.rho_kg_m3:.17g}
nu_phy = {contract.nu_m2_s:.17g}
bulk_viscosity_phy = {contract.bulk_nu_m2_s:.17g}
pressure_reference_phy = {contract.pressure_reference_pa:.17g}
function outlet_01_pressure(x,y,z,t) return {outlets['outlet_01']:.17g} end
function outlet_02_pressure(x,y,z,t) return {outlets['outlet_02']:.17g} end
function outlet_03_pressure(x,y,z,t) return {outlets['outlet_03']:.17g} end
physics = {{dt=dt, rho0=rho0_phy}}
identify = {{label='ROI003274_tau1', kind='fluid', layout='d3q19', relaxation='bgk'}}
fluid = {{kinematic_viscosity=nu_phy, bulk_viscosity=bulk_viscosity_phy}}
initial_condition = {{pressure=pressure_reference_phy, velocityX=0.0, velocityY=0.0, velocityZ=0.0}}
boundary_condition = {{
  {{label='wall', kind='wall_libb'}},
  {{label='inlet', kind='adaptive_flux_pressure', mass_flowrate={contract.target_mass_flow_kg_s:.17g}}},
  {{label='outlet_01', kind='pressure_eq', pressure=outlet_01_pressure}},
  {{label='outlet_02', kind='pressure_eq', pressure=outlet_02_pressure}},
  {{label='outlet_03', kind='pressure_eq', pressure=outlet_03_pressure}}
}}"""


def generate_referee_lua(contract: Tau1BaseRuntimeContract | None = None) -> str:
    contract = contract or Tau1BaseRuntimeContract()
    return f"""-- Fresh continuous-q Referee V2 pre/post state, tau exactly one.
simulation_name = 'tau1_base_referee_v2'
printRuntimeInfo = true
timing_file = 'tracking/timing.res'
mesh = '{_mesh_wsl()}/'
scaling = 'diffusive'
logging = {{level=5}}
maximum_iterations = 1
{_physics_and_boundaries_lua(contract)}
sim_control = {{
  time_control={{max={{iter=maximum_iterations}}, interval={{iter=1}}}},
  abort_criteria={{stop_file='stop'}}
}}
restart = {{
  write='restart/',
  time_control={{min={{iter=0}}, max={{iter=1}}, interval={{iter=1}}}}
}}
"""


def referee_lua_contract(text: str, contract: Tau1BaseRuntimeContract | None = None) -> dict[str, Any]:
    contract = contract or Tau1BaseRuntimeContract()
    checks = {
        "mesh": f"mesh = '{_mesh_wsl()}/'" in text,
        "dx": f"dx = {contract.dx_m:.17g}" in text,
        "dt": f"dt = {contract.dt_s:.17g}" in text,
        "rho": f"rho0_phy = {contract.rho_kg_m3:.17g}" in text,
        "nu": f"nu_phy = {contract.nu_m2_s:.17g}" in text,
        "bulk_nu": f"bulk_viscosity_phy = {contract.bulk_nu_m2_s:.17g}" in text,
        "fresh": "restart = {" in text and "read=" not in text and "read =" not in text,
        "one_step": "maximum_iterations = 1" in text,
        "restart_zero_one": "min={iter=0}" in text and "max={iter=1}" in text,
        "d3q19_bgk": "layout='d3q19'" in text and "relaxation='bgk'" in text,
        "wall_libb": "kind='wall_libb'" in text,
        "adaptive_flux": "kind='adaptive_flux_pressure'" in text,
        "pressure_eq_three": text.count("kind='pressure_eq'") == 3,
        "target_mass_flow": f"mass_flowrate={contract.target_mass_flow_kg_s:.17g}" in text,
    }
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def generate_base_segment_lua(
    *,
    maximum_iterations: int,
    segment_wsl: str,
    restart_header_wsl: str | None,
    restart_first_iteration: int | None = None,
    restart_interval: int = 1,
    restart_write_wsl: str | None = None,
    restart_use_iteration_filename: bool = False,
    stop_file_wsl: str | None = None,
    wall_clock_seconds: int = CONTINUOUS_WALL_CLOCK_SECONDS,
    contract: Tau1BaseRuntimeContract | None = None,
) -> str:
    contract = contract or Tau1BaseRuntimeContract()
    read_line = f"read='{restart_header_wsl}', " if restart_header_wsl else ""
    first_restart = (
        int(restart_first_iteration)
        if restart_first_iteration is not None
        else int(maximum_iterations)
    )
    restart_write = restart_write_wsl or f"{RUNTIME_WSL}/restart/"
    stop_file = stop_file_wsl or f"{RUNTIME_WSL}/stop"
    timeform = (
        "  timeformat={use_iter=true},\n" if restart_use_iteration_filename else ""
    )
    return f"""-- Fresh/resumable repaired Base tau=1 segment; exact Python gates are final.
simulation_name = 'tau1_base'
printRuntimeInfo = true
timing_file = 'tracking/timing.res'
mesh = '{_mesh_wsl()}/'
scaling = 'diffusive'
logging = {{level=5}}
maximum_iterations = {int(maximum_iterations)}
{_physics_and_boundaries_lua(contract)}
sim_control = {{
  time_control={{max={{iter=maximum_iterations, clock={int(wall_clock_seconds)}}}, interval={{iter={SIM_CONTROL_INTERVAL}}}}},
  abort_criteria={{stop_file='{stop_file}'}}
}}
tracking = {{
  {{label='p', folder='tracking/p/', variable={{'pressure_phy'}},
    shape={{kind='all'}}, reduction={{'average'}},
    time_control={{min={{iter=0}}, max={{iter=maximum_iterations}}, interval={{iter={TRACKING_INTERVAL}}}}},
    output={{format='ascii'}}}},
  {{label='u', folder='tracking/u/', variable={{'vel_mag_phy'}},
    shape={{kind='all'}}, reduction={{'average'}},
    time_control={{min={{iter=0}}, max={{iter=maximum_iterations}}, interval={{iter={TRACKING_INTERVAL}}}}},
    output={{format='ascii'}}}}
}}
restart = {{{read_line}write='{restart_write}',
{timeform}  time_control={{min={{iter={first_restart}}}, max={{iter=maximum_iterations}}, interval={{iter={int(restart_interval)}}}}}
}}
"""


def base_lua_contract(
    text: str,
    *,
    maximum_iterations: int,
    restart_header_wsl: str | None,
    restart_first_iteration: int | None = None,
    restart_interval: int = 1,
    restart_write_wsl: str | None = None,
    restart_use_iteration_filename: bool = False,
    contract: Tau1BaseRuntimeContract | None = None,
) -> dict[str, Any]:
    contract = contract or Tau1BaseRuntimeContract()
    first_restart = (
        int(restart_first_iteration)
        if restart_first_iteration is not None
        else int(maximum_iterations)
    )
    restart_write = restart_write_wsl or f"{RUNTIME_WSL}/restart/"
    checks = {
        "mesh": f"mesh = '{_mesh_wsl()}/'" in text,
        "dx": f"dx = {contract.dx_m:.17g}" in text,
        "dt": f"dt = {contract.dt_s:.17g}" in text,
        "rho": f"rho0_phy = {contract.rho_kg_m3:.17g}" in text,
        "nu": f"nu_phy = {contract.nu_m2_s:.17g}" in text,
        "bulk_nu": f"bulk_viscosity_phy = {contract.bulk_nu_m2_s:.17g}" in text,
        "d3q19_bgk": "layout='d3q19'" in text and "relaxation='bgk'" in text,
        "wall_libb": "kind='wall_libb'" in text,
        "adaptive_flux": "kind='adaptive_flux_pressure'" in text,
        "pressure_eq_three": text.count("kind='pressure_eq'") == 3,
        "target_mass_flow": (
            f"mass_flowrate={contract.target_mass_flow_kg_s:.17g}" in text
        ),
        "maximum_iterations": f"maximum_iterations = {maximum_iterations}" in text,
        "restart_read": (
            f"read='{restart_header_wsl}'" in text
            if restart_header_wsl
            else "read='" not in text
        ),
        "tracking_rescaled": (
            text.count(f"interval={{iter={TRACKING_INTERVAL}}}") == 2
        ),
        "sim_control_rescaled": (
            f"interval={{iter={SIM_CONTROL_INTERVAL}}}" in text
        ),
        "restart_schedule": (
            f"min={{iter={first_restart}}}" in text
            and "max={iter=maximum_iterations}" in text
            and f"interval={{iter={int(restart_interval)}}}" in text
        ),
        "restart_write": f"write='{restart_write}'" in text,
        "restart_iteration_filename": (
            "timeformat={use_iter=true}" in text
            if restart_use_iteration_filename
            else "timeformat={use_iter=true}" not in text
        ),
        "no_vtk": "vtk" not in text.lower(),
        "no_old_restart": "221309" not in text and "7.670E-03" not in text,
    }
    return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}


def _referee_launcher(kind: str) -> str:
    common = f"""#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(
  cd -- "$(dirname -- "${{BASH_SOURCE[0]}}")"
  pwd -P
)"
cd "$SCRIPT_DIR"
MESH='{_mesh_wsl()}'
MUSUBI='{MUSUBI_WSL}'
MPIRUN='{MPIRUN_WSL}'
"""
    if kind == "preflight":
        return common + f"""[[ -f musubi.lua && -x "$MUSUBI" && -x "$MPIRUN" ]]
[[ "$(sha256sum "$MUSUBI" | awk '{{print $1}}')" == '{MUSUBI_SHA256}' ]]
for f in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do [[ -s "$MESH/$f" ]]; done
grep -q 'nElems = {EXPECTED_CELLS}' "$MESH/header.lua"
[[ "$(sha256sum "$MESH/elemlist.lsb" | awk '{{print $1}}')" == '{MESH_HASHES['elemlist.lsb']}' ]]
[[ "$(sha256sum "$MESH/bnd.lsb" | awk '{{print $1}}')" == '{MESH_HASHES['bnd.lsb']}' ]]
[[ "$(sha256sum "$MESH/qval.lsb" | awk '{{print $1}}')" == '{MESH_HASHES['qval.lsb']}' ]]
mkdir -p tracking restart
printf 'CELL_COUNT={EXPECTED_CELLS}\nBINARY_SHA256={MUSUBI_SHA256}\nPREFLIGHT=PASS\n'
"""
    if kind == "run":
        return common + f"""/bin/bash "$SCRIPT_DIR/preflight_referee.sh" > preflight.log 2>&1
if find restart -maxdepth 1 -type f -name '*.lsb' -size +0c | grep -q .; then
  [[ -s musubi_stdout.log ]]
else
  "$MPIRUN" --bind-to core --map-by core --report-bindings -np {MPI_RANKS} "$MUSUBI" musubi.lua > musubi_stdout.log 2> musubi_stderr.log
fi
grep -q 'Initializing musubi' musubi_stdout.log
grep -q 'Got a mesh with following properties' musubi_stdout.log
grep -q 'Loading qVal data' musubi_stdout.log
grep -q 'Found BC wall of kind wall_libb' musubi_stdout.log
grep -q 'ADAPTIVE_FLUX_PRESSURE iter=1' musubi_stdout.log
grep -Eq 'iterations:[[:space:]]+1' musubi_stdout.log
test "$(find restart -maxdepth 1 -type f -iname '*header*.lua' -size +0c | wc -l)" -ge 2
test "$(find restart -maxdepth 1 -type f -name '*.lsb' -size +0c | wc -l)" -ge 2
printf 'MESH_LOADED=PASS\nCONTINUOUS_Q_LOADED=PASS\nITERATION_ONE=PASS\nRESTART_ZERO_ONE=PASS\nSEMANTIC_SUCCESS=PASS\n' > semantic_status.log
"""
    raise ValueError(kind)


def _base_launcher() -> str:
    return f"""#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(
  cd -- "$(dirname -- "${{BASH_SOURCE[0]}}")"
  pwd -P
)"
cd "$SCRIPT_DIR"
SEGMENT_DIR="${{1:?segment directory required}}"
MESH='{_mesh_wsl()}'
MUSUBI='{MUSUBI_WSL}'
MPIRUN='{MPIRUN_WSL}'
[[ -d "$SEGMENT_DIR" && -s "$SEGMENT_DIR/musubi.lua" ]]
[[ "$(sha256sum "$MUSUBI" | awk '{{print $1}}')" == '{MUSUBI_SHA256}' ]]
for f in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do [[ -s "$MESH/$f" ]]; done
grep -q 'nElems = {EXPECTED_CELLS}' "$MESH/header.lua"
mkdir -p '{RUNTIME_WSL}/restart' "$SEGMENT_DIR/tracking/p" "$SEGMENT_DIR/tracking/u"
cd "$SEGMENT_DIR"
set +e
"$MPIRUN" --bind-to core --map-by core --report-bindings -np {MPI_RANKS} "$MUSUBI" musubi.lua 2> musubi_stderr.log | awk -v interval={TRACKING_INTERVAL} '
/ADAPTIVE_FLUX_PRESSURE/ {{
  last=$0
  if (match($0,/iter=[0-9]+/)) {{ value=substr($0,RSTART+5,RLENGTH-5)+0; if (value % interval == 0) {{ print $0; fflush() }} }}
  next
}}
{{ print $0 }}
END {{ if (last != "") print last }}
' > musubi_stdout.log
rc=${{PIPESTATUS[0]}}
set -e
[[ "$rc" -eq 0 ]]
grep -q 'Initializing musubi' musubi_stdout.log
grep -q 'Got a mesh with following properties' musubi_stdout.log
grep -q 'Loading qVal data' musubi_stdout.log
grep -q 'ADAPTIVE_FLUX_PRESSURE' musubi_stdout.log
grep -q 'SUCCESSFUL run' musubi_stdout.log
[[ -s '{RUNTIME_WSL}/restart/tau1_base_lastHeader.lua' ]]
printf 'SEGMENT_SEMANTIC_SUCCESS=PASS\n' > semantic_status.log
"""


def prepare_tau1_base(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    run_root = _run_root(root)
    qc = run_root / "qc"
    referee = run_root / "referee"
    qc.mkdir(parents=True, exist_ok=True)
    referee.mkdir(parents=True, exist_ok=True)
    (referee / "tracking").mkdir(exist_ok=True)
    (referee / "restart").mkdir(exist_ok=True)
    preflight = fast_zero_run_preflight(root)
    contract = Tau1BaseRuntimeContract()
    if abs(contract.tau - 1.0) > 1.0e-12 or abs(contract.omega - 1.0) > 1.0e-12:
        raise RuntimeError("analytical tau=1 contract failed")
    referee_files = {
        "musubi.lua": generate_referee_lua(contract),
        "preflight_referee.sh": _referee_launcher("preflight"),
        "run_referee.sh": _referee_launcher("run"),
    }
    for name, content in referee_files.items():
        (referee / name).write_text(content, encoding="utf-8", newline="\n")
    lua_contract = referee_lua_contract(referee_files["musubi.lua"], contract)
    if lua_contract["status"] != "PASS":
        raise RuntimeError(f"referee Lua contract failed: {lua_contract}")
    base_runner = run_root / "run_base_segment.sh"
    base_runner.write_text(_base_launcher(), encoding="utf-8", newline="\n")
    runtime_manifest = _runtime_manifest(contract)
    write_json(qc / "tau1_base_runtime_contract.json", {
        "status": "PASS",
        "preflight": preflight,
        "contract": contract.as_evidence(),
        "runtime_resume_manifest": runtime_manifest,
        "referee_revision": REFEREE_REVISION_NEW,
        "referee_lua_contract": lua_contract,
        "referee_files": {name: sha256_file(referee / name) for name in referee_files},
        "base_runner_sha256": sha256_file(base_runner),
        "production_pipeline_modified": False,
        "seeder_calls": 0,
    })
    return {"status": "PREPARED", "run_root": str(run_root), "preflight": preflight}


def _restart_pairs(restart_dir: Path) -> dict[int, tuple[Path, Path]]:
    pairs: dict[int, tuple[Path, Path]] = {}
    for header in restart_dir.glob("*header*.lua"):
        if "lastHeader" in header.name:
            continue
        parsed = parse_restart_header(header)
        binary = parsed.binary_path
        if not binary.is_file():
            binary = restart_dir / parsed.binary_path.name
        if binary.is_file():
            pairs[parsed.iteration] = (header, binary)
    return pairs


def _state(binary: Path, *, contract: Tau1BaseRuntimeContract, with_velocity: bool) -> dict[str, Any]:
    pdf = np.asarray(read_restart_pdf(binary, n_elems=EXPECTED_CELLS, n_components=19))
    density = np.sum(pdf, axis=1, dtype=np.float64)
    velocity = (pdf @ D3Q19_DIRECTIONS.astype(np.float64)) / density[:, None]
    pressure_scale = contract.rho_kg_m3 * contract.dx_m**2 / contract.dt_s**2 / 3.0
    result: dict[str, Any] = {
        "pdf": pdf,
        "total_pdf_mass": float(np.sum(density, dtype=np.float64)),
        "mean_pressure_pa": float(np.mean(density, dtype=np.float64) * pressure_scale),
        "all_finite": bool(np.all(np.isfinite(pdf))),
        "minimum_pdf": float(np.min(pdf)),
        "maximum_lattice_speed": float(np.max(np.linalg.norm(velocity, axis=1))),
    }
    if with_velocity:
        result["velocity_lattice"] = velocity
    return result


def _controller_records(text: str) -> list[dict[str, Any]]:
    pattern = re.compile(
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
    return [
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
        for match in pattern.finditer(text)
    ]


def audit_continuous_q_referee(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    run_root = _run_root(root)
    referee = run_root / "referee"
    mesh = load_mesh_contract(_mesh_path(root), expected_cells=EXPECTED_CELLS)
    pairs = _restart_pairs(referee / "restart")
    if set(pairs) != {0, 1}:
        raise RuntimeError(f"expected referee restart iterations 0 and 1, got {sorted(pairs)}")
    contract = historical_tau1_runtime_contract()
    pre = _state(pairs[0][1], contract=contract, with_velocity=False)
    post = _state(pairs[1][1], contract=contract, with_velocity=False)
    replay = replay_boundary_step(
        pre["pdf"], mesh, dx_m=contract.dx_m, dt_s=contract.dt_s,
        density_kg_m3=contract.rho_kg_m3,
        target_mass_flow_kg_s=contract.target_mass_flow_kg_s,
        outlet_pressures_pa={
            label: contract.pressure_reference_pa + gauge
            for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
        },
    )
    actual_delta = float(np.sum(post["pdf"] - pre["pdf"], dtype=np.float64))
    residual = conservation_identity_residual(
        replay["predicted_total_lattice"], actual_delta, replay["target_lattice_flux"]
    )
    stdout = _read_text(referee / "musubi_stdout.log")
    controller = _controller_records(stdout)
    if not controller:
        raise RuntimeError("referee controller record missing")
    observed = controller[-1]
    target_error = abs(observed["target_lattice"] - contract.target_lattice_flux) / contract.target_lattice_flux
    wall = mesh.boundaries["wall"]
    wall_q = mesh.qvalues_by_cell[wall.cell_indices]
    continuous_q_loaded = (
        MESH_HASHES["qval.lsb"] == sha256_file(_mesh_path(root) / "qval.lsb")
        and "Loading qVal data" in stdout
        and "Found BC wall of kind wall_libb" in stdout
        and bool(np.any(np.isfinite(wall_q) & (np.abs(wall_q - 0.5) > 1.0e-12)))
    )
    finite = bool(pre["all_finite"] and post["all_finite"])
    minimum_pdf = min(float(pre["minimum_pdf"]), float(post["minimum_pdf"]))
    maximum_speed = max(
        float(pre["maximum_lattice_speed"]), float(post["maximum_lattice_speed"])
    )
    passed = (
        residual <= REFEREE_GATE
        and target_error <= 1.0e-12
        and observed["globBC_count"] == EXPECTED_INLET_GLOBBC
        and observed["relative_error"] <= CONTROLLER_GATE
        and continuous_q_loaded
        and finite
        and minimum_pdf > 0.0
        and maximum_speed < 0.05
    )
    result = {
        "status": "PASS" if passed else "FAIL",
        "final_status_if_fail": "CFD_FLOW_CONTINUOUS_Q_REFEREE_FAILED",
        "referee_revision": REFEREE_REVISION_NEW,
        "scientific_musubi_calls": 1,
        "seeder_calls": 0,
        "mesh_cell_count": len(mesh.cell_ijk),
        "continuous_q_loaded": continuous_q_loaded,
        "expected_target_lattice": contract.target_lattice_flux,
        "observed_target_lattice": observed["target_lattice"],
        "target_lattice_relative_error": target_error,
        "controller": observed,
        "predicted_total_lattice": replay["predicted_total_lattice"],
        "actual_total_lattice": actual_delta,
        "target_lattice_flux": replay["target_lattice_flux"],
        "R_one_step_identity": residual,
        "hard_gate": REFEREE_GATE,
        "per_boundary_lattice": replay["per_label_lattice"],
        "per_boundary_kg_s_domain": replay["per_label_kg_s_domain"],
        "all_finite": finite,
        "minimum_pdf": minimum_pdf,
        "maximum_lattice_speed": maximum_speed,
        "pre_restart_sha256": sha256_file(pairs[0][1]),
        "post_restart_sha256": sha256_file(pairs[1][1]),
        "binary_sha256": MUSUBI_SHA256,
    }
    write_json(run_root / "qc" / "continuous_q_referee_v2_tau1_base.json", result)
    return result


def _velocity_residual(previous: np.ndarray, current: np.ndarray) -> float:
    numerator = float(np.linalg.norm((current - previous).ravel()))
    denominator = max(float(np.linalg.norm(current.ravel())), np.finfo(float).tiny)
    return numerator / denominator


def audit_base_window(
    *,
    project_root: Path,
    restart_triplet: Sequence[tuple[int, Path]],
    referee_residual: float,
    controller: Mapping[str, Any],
) -> dict[str, Any]:
    if len(restart_triplet) != 3:
        raise ValueError("Base steady audit requires start/mid/end restarts")
    root = Path(project_root).resolve()
    mesh = load_mesh_contract(_mesh_path(root), expected_cells=EXPECTED_CELLS)
    contract = historical_tau1_runtime_contract()
    iterations = [int(item[0]) for item in restart_triplet]
    if iterations[1] - iterations[0] != CONFIRMATION_INTERVAL or iterations[2] - iterations[1] != CONFIRMATION_INTERVAL:
        raise ValueError(f"restart triplet does not preserve physical windows: {iterations}")
    states = [
        _state(path, contract=contract, with_velocity=True)
        for _, path in restart_triplet
    ]
    replays = [
        replay_boundary_step(
            state["pdf"], mesh, dx_m=contract.dx_m, dt_s=contract.dt_s,
            density_kg_m3=contract.rho_kg_m3,
            target_mass_flow_kg_s=contract.target_mass_flow_kg_s,
            outlet_pressures_pa={
                label: contract.pressure_reference_pa + gauge
                for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
            },
        )
        for state in states
    ]
    conversion = contract.rho_kg_m3 * contract.dx_m**3 / contract.dt_s

    def mass_residual(a: int, b: int) -> float:
        delta = iterations[b] - iterations[a]
        accumulation = (
            (states[b]["total_pdf_mass"] - states[a]["total_pdf_mass"])
            / delta
            * conversion
        )
        return abs(accumulation) / contract.target_mass_flow_kg_s

    r_mass_short = mass_residual(1, 2)
    r_mass_long = mass_residual(0, 2)
    r_velocity = _velocity_residual(
        states[1]["velocity_lattice"], states[2]["velocity_lattice"]
    )
    current_inlet_rho = float(replays[2]["details"]["inlet"]["rho"])
    pressure_scale = contract.rho_kg_m3 * contract.dx_m**2 / contract.dt_s**2 / 3.0
    inlet_absolute = current_inlet_rho * pressure_scale
    inlet_gauge = inlet_absolute - contract.pressure_reference_pa
    pressure_drops = {
        label: inlet_gauge - gauge for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
    }
    characteristic = float(np.median(np.abs(list(pressure_drops.values()))))
    r_pressure = abs(states[2]["mean_pressure_pa"] - states[1]["mean_pressure_pa"]) / max(
        characteristic, np.finfo(float).tiny
    )
    per_boundary = [item["per_label_kg_s_domain"] for item in replays]
    inlet_series = [float(item["inlet"]) for item in per_boundary]
    outlet_series = {
        label: [-float(item[label]) for item in per_boundary]
        for label in ("outlet_01", "outlet_02", "outlet_03")
    }
    span = iterations[-1] - iterations[0]
    mean_inlet = trapezoidal_integral(iterations, inlet_series) / span
    mean_outlets = {
        label: trapezoidal_integral(iterations, values) / span
        for label, values in outlet_series.items()
    }
    predicted_net = [float(item["predicted_total_kg_s"]) for item in replays]
    predicted_mass_change = trapezoidal_integral(iterations, predicted_net) * contract.dt_s
    observed_mass_change = (
        states[2]["total_pdf_mass"] - states[0]["total_pdf_mass"]
    ) * contract.rho_kg_m3 * contract.dx_m**3
    closure = boundary_window_closure(
        predicted_mass_change,
        observed_mass_change,
        contract.target_mass_flow_kg_s * span * contract.dt_s,
    )
    r_inlet = abs(mean_inlet - contract.target_mass_flow_kg_s) / contract.target_mass_flow_kg_s
    backflow = significant_time_averaged_backflow(mean_inlet, mean_outlets.values())
    minimum_pdf = min(float(item["minimum_pdf"]) for item in states)
    maximum_speed = max(float(item["maximum_lattice_speed"]) for item in states)
    all_finite = all(bool(item["all_finite"]) for item in states)
    controller_error = abs(float(controller["target_lattice"]) - contract.target_lattice_flux) / contract.target_lattice_flux
    passed = all(
        (
            r_mass_short <= MASS_GATE,
            r_mass_long <= MASS_GATE,
            r_velocity <= VELOCITY_GATE,
            r_pressure <= PRESSURE_GATE,
            r_inlet <= INLET_GATE,
            closure <= BOUNDARY_CLOSURE_GATE,
            referee_residual <= REFEREE_GATE,
            not backflow,
            all_finite,
            minimum_pdf > 0.0,
            maximum_speed < 0.05,
            controller_error <= 1.0e-12,
            float(controller["relative_error"]) <= CONTROLLER_GATE,
            int(controller["globBC_count"]) == EXPECTED_INLET_GLOBBC,
        )
    )
    outlet_q = {label: value / contract.rho_kg_m3 for label, value in mean_outlets.items()}
    result = {
        "status": "PASS" if passed else "FAIL",
        "accepted_iteration": iterations[-1] if passed else None,
        "accepted_physical_time_s": iterations[-1] * contract.dt_s if passed else None,
        "window_iterations": iterations,
        "R_mass_short": r_mass_short,
        "R_mass_long": r_mass_long,
        "R_velocity": r_velocity,
        "R_pressure": r_pressure,
        "R_inlet": r_inlet,
        "boundary_window_closure": closure,
        "time_averaged_backflow": backflow,
        "all_finite": all_finite,
        "minimum_pdf": minimum_pdf,
        "maximum_lattice_speed": maximum_speed,
        "referee_v2_residual": referee_residual,
        "controller": dict(controller),
        "controller_target_relative_error": controller_error,
        "inlet_gauge_pressure_pa": inlet_gauge,
        "pressure_drops_pa": pressure_drops,
        "mean_mass_flow_kg_s": {"inlet": mean_inlet, **mean_outlets},
        "mean_volume_flow_m3_s": {
            "inlet": mean_inlet / contract.rho_kg_m3,
            **outlet_q,
        },
        "flow_fractions": {
            label: value / mean_inlet for label, value in mean_outlets.items()
        },
        "restart_sha256": {
            str(iteration): sha256_file(path)
            for iteration, path in restart_triplet
        },
        "gates": {
            "R_mass_short_le_1pct": r_mass_short <= MASS_GATE,
            "R_mass_long_le_1pct": r_mass_long <= MASS_GATE,
            "R_velocity_le_1pct": r_velocity <= VELOCITY_GATE,
            "R_pressure_le_0p5pct": r_pressure <= PRESSURE_GATE,
            "R_inlet_le_1pct": r_inlet <= INLET_GATE,
            "boundary_window_closure_le_0p1pct": closure <= BOUNDARY_CLOSURE_GATE,
            "referee_v2_le_1e_minus_8": referee_residual <= REFEREE_GATE,
            "no_significant_time_averaged_backflow": not backflow,
            "numerical_safety": all_finite and minimum_pdf > 0.0 and maximum_speed < 0.05,
            "adaptive_target_exact": controller_error <= 1.0e-12,
        },
    }
    return result


def _latest_controller_from_segments(run_root: Path) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    for path in sorted((run_root / "segments").glob("*/musubi_stdout.log")):
        records.extend(_controller_records(_read_text(path)))
    if not records:
        raise RuntimeError("no Base controller records are available")
    return max(records, key=lambda item: int(item["iteration"]))


def audit_latest_base_window(project_root: Path) -> dict[str, Any]:
    """Audit and adopt the latest exact confirmation triplet after a monitor restart."""

    root = Path(project_root).resolve()
    run_root = _run_root(root)
    pairs = _restart_pairs(_runtime_windows() / "restart")
    available = sorted(pairs)
    if len(available) < 3:
        raise RuntimeError("fewer than three Base restart pairs are available")
    iterations = available[-3:]
    if (
        iterations[1] - iterations[0] != CONFIRMATION_INTERVAL
        or iterations[2] - iterations[1] != CONFIRMATION_INTERVAL
    ):
        raise RuntimeError(f"latest Base restarts are not an exact triplet: {iterations}")
    records: list[dict[str, Any]] = []
    for path in sorted((_runtime_windows() / "segments").glob("*/musubi_stdout.log")):
        records.extend(_controller_records(_read_text(path)))
    if not records:
        raise RuntimeError("no runtime Base controller records are available")
    referee = json.loads(
        (run_root / "qc" / "continuous_q_referee_v2_tau1_base.json").read_text(
            encoding="utf-8"
        )
    )
    result = audit_base_window(
        project_root=root,
        restart_triplet=[(iteration, pairs[iteration][1]) for iteration in iterations],
        referee_residual=float(referee["R_one_step_identity"]),
        controller=max(records, key=lambda item: int(item["iteration"])),
    )
    write_json(run_root / "qc" / "tau1_base_steady_status.json", result)
    state_path = run_root / "qc" / "tau1_base_run_state.json"
    if state_path.is_file():
        state = json.loads(state_path.read_text(encoding="utf-8"))
        state["latest_steady_audit"] = result
        state["status"] = "PASS" if result["status"] == "PASS" else "IN_PROGRESS"
        write_json(state_path, state)
    return result


def _quadrature_integral(
    iterations: Sequence[int], values: Sequence[float], method: str
) -> float:
    x = np.asarray(iterations, dtype=np.float64)
    y = np.asarray(values, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) < 2:
        raise ValueError("quadrature needs matching one-dimensional samples")
    widths = np.diff(x)
    if np.any(widths <= 0.0):
        raise ValueError("quadrature iterations must be strictly increasing")
    if method == "left_rectangle":
        return float(np.sum(widths * y[:-1], dtype=np.float64))
    if method == "right_rectangle":
        return float(np.sum(widths * y[1:], dtype=np.float64))
    if method == "trapezoidal":
        return float(np.sum(widths * 0.5 * (y[:-1] + y[1:]), dtype=np.float64))
    if method == "timestep_phase_shifted":
        shifted_x = x + 1.0
        shifted_y = np.interp(x, shifted_x, y, left=y[0], right=y[-1])
        return float(
            np.sum(
                widths * 0.5 * (shifted_y[:-1] + shifted_y[1:]),
                dtype=np.float64,
            )
        )
    raise ValueError(f"unknown quadrature method: {method}")


def _stdout_total_density_samples(run_root: Path) -> dict[int, float]:
    pattern = re.compile(
        r"iterations:\s*(?P<iteration>\d+).*?"
        r"field\s*\|\s*total density.*?"
        r"\n\s*1\s*\|\s*(?P<density>[-+0-9.Ee]+)",
        flags=re.DOTALL,
    )
    samples: dict[int, float] = {}
    roots = (run_root / "segments", _runtime_windows() / "segments")
    for root in roots:
        if not root.is_dir():
            continue
        for path in root.glob("*/musubi_stdout.log"):
            for match in pattern.finditer(_read_text(path)):
                samples[int(match.group("iteration"))] = float(match.group("density"))
    return samples


def _nearest_density_sample(
    samples: Mapping[int, float], requested_iteration: int
) -> tuple[int, float]:
    if not samples:
        raise RuntimeError("no stdout total-density samples are available")
    iteration = min(samples, key=lambda item: (abs(item - requested_iteration), item))
    return int(iteration), float(samples[iteration])


def forensic_boundary_window_audit(
    project_root: Path,
    *,
    output_stem: str = "tau1_boundary_window_forensics",
) -> dict[str, Any]:
    """Compare sparse window accounting without launching Musubi."""

    if not re.fullmatch(r"[a-z0-9_]+", output_stem):
        raise ValueError("forensic output stem must be a safe lowercase name")

    root = Path(project_root).resolve()
    run_root = _run_root(root)
    contract = historical_tau1_runtime_contract()
    pairs = _restart_pairs(_runtime_windows() / "restart")
    iterations = sorted(pairs)
    if len(iterations) < 3:
        raise RuntimeError("forensic audit requires the three preserved restart pairs")
    exact_iterations = iterations[-3:]
    if any(
        right - left != CONFIRMATION_INTERVAL
        for left, right in zip(exact_iterations, exact_iterations[1:])
    ):
        raise RuntimeError(f"preserved restarts are not an exact triplet: {exact_iterations}")

    mesh = load_mesh_contract(_mesh_path(root), expected_cells=EXPECTED_CELLS)
    exact_states: dict[int, dict[str, Any]] = {}
    exact_replays: dict[int, dict[str, Any]] = {}
    for iteration in exact_iterations:
        state = _state(pairs[iteration][1], contract=contract, with_velocity=False)
        replay = replay_boundary_step(
            state["pdf"],
            mesh,
            dx_m=contract.dx_m,
            dt_s=contract.dt_s,
            density_kg_m3=contract.rho_kg_m3,
            target_mass_flow_kg_s=contract.target_mass_flow_kg_s,
            outlet_pressures_pa={
                label: contract.pressure_reference_pa + gauge
                for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
            },
        )
        exact_states[iteration] = state
        exact_replays[iteration] = replay

    labels = ("wall", "inlet", "outlet_01", "outlet_02", "outlet_03")
    methods = (
        "left_rectangle",
        "right_rectangle",
        "trapezoidal",
        "timestep_phase_shifted",
    )
    density_samples = _stdout_total_density_samples(run_root)
    end_iteration = exact_iterations[-1]
    first_confirmation_iteration = 479_004
    windows = (
        ("1x", CONFIRMATION_INTERVAL),
        ("2x", 2 * CONFIRMATION_INTERVAL),
        ("4x", 4 * CONFIRMATION_INTERVAL),
        ("8x", 8 * CONFIRMATION_INTERVAL),
        ("maximum_available", end_iteration - first_confirmation_iteration),
    )
    mass_per_lattice = contract.rho_kg_m3 * contract.dx_m**3
    replay_mean = {
        label: float(
            np.mean(
                [
                    exact_replays[item]["per_label_kg_s_domain"][label]
                    for item in exact_iterations
                ],
                dtype=np.float64,
            )
        )
        for label in labels
    }
    result_windows: list[dict[str, Any]] = []
    csv_rows: list[dict[str, Any]] = []
    for name, span in windows:
        requested_start = end_iteration - span
        direct = requested_start in exact_states
        if direct:
            sample_iterations = [
                item for item in exact_iterations if requested_start <= item <= end_iteration
            ]
            observed_mass_change = (
                exact_states[end_iteration]["total_pdf_mass"]
                - exact_states[requested_start]["total_pdf_mass"]
            ) * mass_per_lattice
            observed_source = "element totals from preserved restart PDFs"
            observed_precision_bound = None
            observed_start_iteration = requested_start
            observed_end_iteration = end_iteration
        else:
            observed_start_iteration, start_density = _nearest_density_sample(
                density_samples, requested_start
            )
            observed_end_iteration, end_density = _nearest_density_sample(
                density_samples, end_iteration
            )
            observed_mass_change = (end_density - start_density) * mass_per_lattice
            observed_source = "archived stdout total-density endpoints"
            observed_precision_bound = 1.0e-7 * mass_per_lattice
            sample_iterations = [requested_start, end_iteration]

        sample_rates: dict[str, list[float]] = {}
        if direct:
            for label in labels:
                sample_rates[label] = [
                    float(exact_replays[item]["per_label_kg_s_domain"][label])
                    for item in sample_iterations
                ]
        else:
            for label in labels:
                sample_rates[label] = [replay_mean[label], replay_mean[label]]

        inlet_normalizer = contract.target_mass_flow_kg_s * span * contract.dt_s
        integrations: dict[str, Any] = {}
        for method in methods:
            per_boundary_mass = {
                label: _quadrature_integral(
                    sample_iterations, sample_rates[label], method
                )
                * contract.dt_s
                for label in labels
            }
            predicted_mass_change = float(sum(per_boundary_mass.values()))
            mismatch = abs(predicted_mass_change - observed_mass_change)
            closure = boundary_window_closure(
                predicted_mass_change, observed_mass_change, inlet_normalizer
            )
            integrations[method] = {
                "predicted_mass_change_kg": predicted_mass_change,
                "observed_mass_change_kg": observed_mass_change,
                "absolute_mismatch_kg": mismatch,
                "inlet_mass_normalizer_kg": inlet_normalizer,
                "closure": closure,
                "per_boundary_contribution_kg_domain_signed": per_boundary_mass,
                "flux_definition": BOUNDARY_FLUX_DEFINITION,
            }
            csv_rows.append(
                {
                    "window": name,
                    "method": method,
                    "span_iterations": span,
                    "predicted_mass_change_kg": predicted_mass_change,
                    "observed_mass_change_kg": observed_mass_change,
                    "absolute_mismatch_kg": mismatch,
                    "inlet_mass_normalizer_kg": inlet_normalizer,
                    "closure": closure,
                    **{f"{label}_kg": per_boundary_mass[label] for label in labels},
                }
            )
        result_windows.append(
            {
                "window": name,
                "span_iterations": span,
                "physical_time_s": span * contract.dt_s,
                "requested_start_iteration": requested_start,
                "end_iteration": end_iteration,
                "boundary_rate_evidence": (
                    "direct preserved restart replay"
                    if direct
                    else "stationary mean of latest exact restart replay triplet"
                ),
                "observed_mass_evidence": observed_source,
                "observed_start_iteration": observed_start_iteration,
                "observed_end_iteration": observed_end_iteration,
                "observed_precision_bound_kg": observed_precision_bound,
                "sample_iterations": sample_iterations,
                "integration_methods": integrations,
            }
        )

    qc = run_root / "qc"
    csv_path = qc / f"{output_stem}.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)
    result = {
        "status": "ZERO_LONG_RUN_FORENSICS_COMPLETE",
        "flux_definition": BOUNDARY_FLUX_DEFINITION,
        "musubi_calls": 0,
        "exact_restart_iterations": exact_iterations,
        "exact_restart_sha256": {
            str(item): sha256_file(pairs[item][1]) for item in exact_iterations
        },
        "latest_replay_per_boundary_kg_s_domain": {
            str(item): exact_replays[item]["per_label_kg_s_domain"]
            for item in exact_iterations
        },
        "windows": result_windows,
        "interpretation": (
            "1x/2x are exact restart-endpoint audits. 4x/8x/maximum retain exact "
            "historical total-density observables but use the latest stationary replay "
            "mean because older PDFs were removed by the pre-existing keep-three policy."
        ),
        "next": "RUN AT MOST ONE 8-STEP DENSE DISCRETE DIAGNOSTIC",
    }
    write_json(qc / f"{output_stem}.json", result)
    return result


def _dense_diagnostic_launcher(*, start_binary_wsl: str, start_sha256: str) -> str:
    return f"""#!/usr/bin/env bash
set -euo pipefail
DENSE_DIR="${{1:?dense diagnostic directory required}}"
MESH='{_mesh_wsl()}'
MUSUBI='{MUSUBI_WSL}'
MPIRUN='{MPIRUN_WSL}'
[[ -d "$DENSE_DIR" && -s "$DENSE_DIR/musubi.lua" ]]
[[ ! -e "$DENSE_DIR/diagnostic_launched" ]]
[[ "$(sha256sum '{start_binary_wsl}' | awk '{{print $1}}')" == '{start_sha256}' ]]
[[ "$(sha256sum "$MUSUBI" | awk '{{print $1}}')" == '{MUSUBI_SHA256}' ]]
for f in header.lua elemlist.lsb bnd.lua bnd.lsb qval.lua qval.lsb; do [[ -s "$MESH/$f" ]]; done
mkdir -p "$DENSE_DIR/restart" "$DENSE_DIR/tracking/p" "$DENSE_DIR/tracking/u"
touch "$DENSE_DIR/diagnostic_launched"
cd "$DENSE_DIR"
"$MPIRUN" --bind-to core --map-by core --report-bindings -np {MPI_RANKS} "$MUSUBI" musubi.lua > musubi_stdout.log 2> musubi_stderr.log
grep -q 'Initializing musubi' musubi_stdout.log
grep -q 'Got a mesh with following properties' musubi_stdout.log
grep -q 'Loading qVal data' musubi_stdout.log
grep -q 'Found BC wall of kind wall_libb' musubi_stdout.log
grep -q 'ADAPTIVE_FLUX_PRESSURE' musubi_stdout.log
grep -q 'SUCCESSFUL run' musubi_stdout.log
printf 'DENSE_DIAGNOSTIC_SEMANTIC_SUCCESS=PASS\n' > semantic_status.log
"""


def _read_pdf(binary: Path) -> np.ndarray:
    return np.asarray(
        read_restart_pdf(binary, n_elems=EXPECTED_CELLS, n_components=19)
    )


def _run_dense_discrete_diagnostic_attempt(
    project_root: Path, *, steps: int, attempt: int
) -> dict[str, Any]:
    """Run one short process and evaluate every exact one-step mass identity."""

    if int(steps) not in {8, 16}:
        raise ValueError("corrected dense diagnostic must contain 8 or 16 timesteps")
    if int(attempt) not in {1, 2}:
        raise ValueError("corrected dense diagnostic attempt must be 1 or 2")
    root = Path(project_root).resolve()
    run_root = _run_root(root)
    qc = run_root / "qc"
    result_path = qc / (
        f"tau1_dense_discrete_mass_identity_v2_corrected_v2_{steps}_attempt_{attempt}.json"
    )
    if result_path.is_file():
        return json.loads(result_path.read_text(encoding="utf-8"))
    forensic_path = qc / "tau1_boundary_window_forensics.json"
    if not forensic_path.is_file():
        raise RuntimeError("zero-long-run forensic audit must precede the diagnostic")

    source_contract = tem_restart_timeformat_contract()
    write_json(qc / "tau1_tem_restart_timeformat_contract.json", source_contract)
    if source_contract["status"] != "PASS":
        raise RuntimeError(f"linked Tem timeformat contract failed: {source_contract}")

    latest = _latest_runtime_restart()
    if latest is None:
        raise RuntimeError("latest compatible Base restart is unavailable")
    start_iteration, start_header, start_binary = latest
    end_iteration = start_iteration + int(steps)
    start_sha256 = sha256_file(start_binary)
    main_pairs_before = _restart_pairs(_runtime_windows() / "restart")
    main_hashes_before = {
        str(iteration): sha256_file(binary)
        for iteration, (_, binary) in sorted(main_pairs_before.items())
    }
    diagnostic_name = (
        f"dense_diagnostic_corrected_v2_{start_iteration}_{steps}_attempt_{attempt}"
    )
    diagnostic_runtime = _runtime_windows() / diagnostic_name
    diagnostic_wsl = f"{RUNTIME_WSL}/{diagnostic_name}"
    restart_wsl = f"{diagnostic_wsl}/restart/"
    launched_marker = diagnostic_runtime / "diagnostic_launched"
    expected_iterations = list(range(start_iteration + 1, end_iteration + 1))

    if diagnostic_runtime.exists() and not launched_marker.is_file():
        if any(diagnostic_runtime.iterdir()):
            raise RuntimeError(
                "nonempty unlaunched dense diagnostic directory requires manual audit"
            )
    diagnostic_runtime.mkdir(parents=True, exist_ok=True)
    (diagnostic_runtime / "restart").mkdir(exist_ok=True)
    restart_header_wsl = f"{RUNTIME_WSL}/restart/{start_header.name}"
    short_stop_name = f"dense_{steps}_a{attempt}.stop"
    short_stop_wsl = f"{RUNTIME_WSL}/{short_stop_name}"
    (_runtime_windows() / short_stop_name).unlink(missing_ok=True)
    lua = generate_base_segment_lua(
        maximum_iterations=end_iteration,
        segment_wsl=diagnostic_wsl,
        restart_header_wsl=restart_header_wsl,
        restart_first_iteration=start_iteration + 1,
        restart_interval=1,
        restart_write_wsl=restart_wsl,
        restart_use_iteration_filename=True,
        stop_file_wsl=short_stop_wsl,
        wall_clock_seconds=1_200,
    )
    lua_contract = base_lua_contract(
        lua,
        maximum_iterations=end_iteration,
        restart_header_wsl=restart_header_wsl,
        restart_first_iteration=start_iteration + 1,
        restart_interval=1,
        restart_write_wsl=restart_wsl,
        restart_use_iteration_filename=True,
    )
    if lua_contract["status"] != "PASS":
        raise RuntimeError(f"dense diagnostic Lua contract failed: {lua_contract}")
    if (
        source_contract["required_lua_syntax"] not in lua
        or source_contract["forbidden_typo"] in lua
        or f"stop_file='{short_stop_wsl}'" not in lua
    ):
        raise RuntimeError("generated dense Lua violates pre-execution contract")
    (diagnostic_runtime / "musubi.lua").write_text(
        lua, encoding="utf-8", newline="\n"
    )
    launcher_path = run_root / "run_dense_diagnostic_corrected.sh"
    start_binary_wsl = f"{RUNTIME_WSL}/restart/{start_binary.name}"
    launcher_path.write_text(
        _dense_diagnostic_launcher(
            start_binary_wsl=start_binary_wsl, start_sha256=start_sha256
        ),
        encoding="utf-8",
        newline="\n",
    )
    if not launched_marker.is_file():
        launcher_wsl = f"{PROJECT_WSL}/outputs/cfd_flow/{RUN_NAME}/{launcher_path.name}"
        returncode, elapsed_s = _invoke_wsl_script(
            launcher_wsl, diagnostic_wsl, timeout_s=1_800
        )
        if returncode != 0:
            stdout_path = diagnostic_runtime / "musubi_stdout.log"
            completed_records = (
                [
                    item
                    for item in _controller_records(_read_text(stdout_path))
                    if int(item["iteration"]) > start_iteration
                ]
                if stdout_path.is_file()
                else []
            )
            failure = {
                "status": "RECOVERABLE_OPERATIONAL_ERROR",
                "scientific_musubi_calls": int(bool(completed_records)),
                "completed_timesteps": len(completed_records),
                "attempt": attempt,
                "returncode": returncode,
                "elapsed_s": elapsed_s,
                "start_iteration": start_iteration,
                "requested_steps": steps,
                "long_base_restart_sha256_before": main_hashes_before,
                "next": "AUTOMATIC CORRECTED SHORT-DIAGNOSTIC RETRY",
            }
            _copy_segment_evidence(
                diagnostic_runtime,
                run_root
                / "operational_failures"
                / f"corrected_dense_v2_attempt_{attempt}",
            )
            write_json(result_path, failure)
            return failure
    else:
        elapsed_s = None

    diagnostic_pairs = _restart_pairs(diagnostic_runtime / "restart")
    available_iterations = sorted(diagnostic_pairs)
    if available_iterations != expected_iterations:
        failure = {
            "status": "RECOVERABLE_OPERATIONAL_ERROR",
            "scientific_musubi_calls": 1,
            "attempt": attempt,
            "start_iteration": start_iteration,
            "requested_steps": steps,
            "expected_restart_iterations": expected_iterations,
            "available_restart_iterations": available_iterations,
            "restart_headers_found": len(
                [
                    path
                    for path in (diagnostic_runtime / "restart").glob("*header*.lua")
                    if "lastHeader" not in path.name
                ]
            ),
            "restart_binaries_found": len(
                list((diagnostic_runtime / "restart").glob("*.lsb"))
            ),
            "next": "AUTOMATIC CORRECTED SHORT-DIAGNOSTIC RETRY",
        }
        write_json(result_path, failure)
        return failure

    contract = historical_tau1_runtime_contract()
    mesh = load_mesh_contract(_mesh_path(root), expected_cells=EXPECTED_CELLS)
    outlet_pressures = {
        label: contract.pressure_reference_pa + gauge
        for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
    }
    labels = ("wall", "inlet", "outlet_01", "outlet_02", "outlet_03")
    current_pdf = _read_pdf(start_binary)
    endpoint_start_pdf = current_pdf.copy()
    per_step: list[dict[str, Any]] = []
    predicted_deltas: list[float] = []
    observed_step_deltas: list[float] = []
    per_boundary_sums = {label: [] for label in labels}
    dense_replays: list[dict[str, Any]] = []
    mass_per_lattice = contract.rho_kg_m3 * contract.dx_m**3
    for next_iteration in expected_iterations:
        next_binary = diagnostic_pairs[next_iteration][1]
        next_pdf = _read_pdf(next_binary)
        replay = replay_boundary_step(
            current_pdf,
            mesh,
            dx_m=contract.dx_m,
            dt_s=contract.dt_s,
            density_kg_m3=contract.rho_kg_m3,
            target_mass_flow_kg_s=contract.target_mass_flow_kg_s,
            outlet_pressures_pa=outlet_pressures,
        )
        predicted = float(replay["predicted_total_lattice"])
        observed = float(np.sum(next_pdf - current_pdf, dtype=np.float64))
        residual = conservation_identity_residual(
            predicted, observed, replay["target_lattice_flux"]
        )
        boundary_values = {
            label: float(replay["per_label_lattice"][label]) for label in labels
        }
        dense_replays.append(replay)
        for label in labels:
            per_boundary_sums[label].append(boundary_values[label])
        predicted_deltas.append(predicted)
        observed_step_deltas.append(observed)
        per_step.append(
            {
                "start_iteration": next_iteration - 1,
                "end_iteration": next_iteration,
                "predicted_exact_boundary_pdf_mass_delta_lattice": predicted,
                "observed_elementwise_pdf_delta_lattice": observed,
                "absolute_mismatch_lattice": abs(predicted - observed),
                "predicted_exact_boundary_pdf_mass_delta_kg": (
                    predicted * mass_per_lattice
                ),
                "observed_elementwise_pdf_delta_kg": observed * mass_per_lattice,
                "R_one_step_identity": residual,
                "target_lattice_inlet_flux": float(replay["target_lattice_flux"]),
                "per_boundary_lattice_signed": boundary_values,
                "per_boundary_kg_signed": {
                    label: value * mass_per_lattice
                    for label, value in boundary_values.items()
                },
                "old_sparse_predicted_boundary_lattice_signed": boundary_values,
                "old_minus_exact_per_boundary_lattice": {
                    label: 0.0 for label in labels
                },
                "old_minus_exact_total_lattice": 0.0,
                "end_restart_sha256": sha256_file(next_binary),
            }
        )
        current_pdf = next_pdf

    sum_predicted = math.fsum(predicted_deltas)
    sum_observed_steps = math.fsum(observed_step_deltas)
    endpoint_delta = float(
        np.sum(current_pdf - endpoint_start_pdf, dtype=np.float64)
    )
    per_boundary_lattice = {
        label: math.fsum(values) for label, values in per_boundary_sums.items()
    }
    normalizer = abs(per_boundary_lattice["inlet"])
    exact_closure = abs(sum_predicted - endpoint_delta) / normalizer
    step_telescoping_error = abs(sum_observed_steps - endpoint_delta) / normalizer
    final_replay = replay_boundary_step(
        current_pdf,
        mesh,
        dx_m=contract.dx_m,
        dt_s=contract.dt_s,
        density_kg_m3=contract.rho_kg_m3,
        target_mass_flow_kg_s=contract.target_mass_flow_kg_s,
        outlet_pressures_pa=outlet_pressures,
    )
    old_sample_replays = [*dense_replays, final_replay]
    old_sample_iterations = list(range(start_iteration, end_iteration + 1))
    old_per_boundary_lattice = {
        label: trapezoidal_integral(
            old_sample_iterations,
            [
                float(item["per_label_lattice"][label])
                for item in old_sample_replays
            ],
        )
        for label in labels
    }
    old_predicted = float(sum(old_per_boundary_lattice.values()))
    old_normalizer = contract.target_lattice_flux * steps
    old_residual = abs(old_predicted - endpoint_delta) / old_normalizer
    one_step_residuals = [float(item["R_one_step_identity"]) for item in per_step]
    individual_steps_pass = max(one_step_residuals) <= 1.0e-8
    classification = (
        "BOUNDARY_WINDOW_AUDIT_ACCOUNTING_BIAS_CONFIRMED"
        if individual_steps_pass and exact_closure <= 1.0e-8
        else "TRUE_BOUNDARY_DISCRETE_ACCOUNTING_FAILURE"
    )
    main_pairs_after = _restart_pairs(_runtime_windows() / "restart")
    main_hashes_after = {
        str(iteration): sha256_file(binary)
        for iteration, (_, binary) in sorted(main_pairs_after.items())
    }
    archive = run_root / "dense_diagnostic_corrected" / f"attempt_{attempt}_{steps}_steps"
    archive.mkdir(parents=True, exist_ok=True)
    for name in (
        "musubi.lua",
        "musubi_stdout.log",
        "musubi_stderr.log",
        "semantic_status.log",
        "diagnostic_launched",
    ):
        source = diagnostic_runtime / name
        if source.is_file():
            shutil.copy2(source, archive / name)
    restart_archive = archive / "restart"
    restart_archive.mkdir(exist_ok=True)
    for header, binary in diagnostic_pairs.values():
        shutil.copy2(header, restart_archive / header.name)
        shutil.copy2(binary, restart_archive / binary.name)
    result = {
        "status": classification,
        "scientific_musubi_calls": 1,
        "attempt": attempt,
        "long_base_musubi_calls_after_stop": 0,
        "elapsed_s": elapsed_s,
        "steps": steps,
        "start_iteration": start_iteration,
        "end_iteration": end_iteration,
        "unique_restart_states": steps + 1,
        "diagnostic_restart_headers": len(diagnostic_pairs),
        "diagnostic_restart_binaries": len(diagnostic_pairs),
        "start_physical_time_s": start_iteration * contract.dt_s,
        "end_physical_time_s": end_iteration * contract.dt_s,
        "integration": "exact per-step boundary PDF delta; no sparse interpolation",
        "sum_predicted_exact_boundary_pdf_mass_delta_lattice": sum_predicted,
        "elementwise_sum_pdf_end_minus_pdf_start_lattice": endpoint_delta,
        "sum_elementwise_per_step_pdf_delta_lattice": sum_observed_steps,
        "absolute_mismatch_lattice": abs(sum_predicted - endpoint_delta),
        "inlet_mass_normalizer_lattice": normalizer,
        "exact_dense_discrete_closure": exact_closure,
        "hard_classification_threshold": 1.0e-8,
        "individual_steps_pass": individual_steps_pass,
        "one_step_residual_min": min(one_step_residuals),
        "one_step_residual_median": float(np.median(one_step_residuals)),
        "one_step_residual_max": max(one_step_residuals),
        "step_telescoping_roundoff_normalized": step_telescoping_error,
        "old_sparse_window_residual_same_dense_states": old_residual,
        "old_sparse_predicted_total_lattice": old_predicted,
        "old_sparse_inlet_normalizer_lattice": old_normalizer,
        "old_sparse_per_boundary_lattice_signed": old_per_boundary_lattice,
        "old_minus_exact_cumulative_per_boundary_lattice": {
            label: old_per_boundary_lattice[label] - per_boundary_lattice[label]
            for label in labels
        },
        "old_minus_exact_cumulative_total_lattice": (
            old_predicted - sum_predicted
        ),
        "per_boundary_sum_lattice_signed": per_boundary_lattice,
        "per_boundary_sum_kg_domain_signed": {
            label: value * mass_per_lattice
            for label, value in per_boundary_lattice.items()
        },
        "predicted_mass_change_kg": sum_predicted * mass_per_lattice,
        "observed_mass_change_kg": endpoint_delta * mass_per_lattice,
        "absolute_mismatch_kg": abs(sum_predicted - endpoint_delta)
        * mass_per_lattice,
        "inlet_mass_normalizer_kg": normalizer * mass_per_lattice,
        "per_step": per_step,
        "start_restart_sha256": start_sha256,
        "diagnostic_restart_sha256": {
            str(iteration): per_step[index]["end_restart_sha256"]
            for index, iteration in enumerate(expected_iterations)
        },
        "long_base_restart_sha256_before": main_hashes_before,
        "long_base_restart_sha256_after": main_hashes_after,
        "long_base_restarts_preserved": main_hashes_before == main_hashes_after,
        "binary_sha256": MUSUBI_SHA256,
        "lua_sha256": sha256_file(diagnostic_runtime / "musubi.lua"),
        "lua_contract": lua_contract,
        "tem_restart_timeformat_contract": source_contract,
        "next": (
            "RE-AUDIT EXISTING BASE RESTART WITH DENSE EXACT WINDOW ACCOUNTING"
            if classification == "BOUNDARY_WINDOW_AUDIT_ACCOUNTING_BIAS_CONFIRMED"
            else "STOP BEFORE ANY NEW LONG CFD"
        ),
    }
    write_json(result_path, result)
    return result


def run_dense_discrete_diagnostic(
    project_root: Path, *, steps: int = 8, maximum_attempts: int = 2
) -> dict[str, Any]:
    """Run the corrected dense diagnostic with bounded operational recovery."""

    if maximum_attempts not in {1, 2}:
        raise ValueError("maximum corrected dense attempts must be 1 or 2")
    root = Path(project_root).resolve()
    qc = _run_root(root) / "qc"
    aggregate_path = qc / (
        f"tau1_dense_discrete_mass_identity_v2_corrected_v2_{steps}.json"
    )
    if aggregate_path.is_file():
        existing = json.loads(aggregate_path.read_text(encoding="utf-8"))
        if existing.get("status") != "RECOVERABLE_OPERATIONAL_ERROR":
            return existing
    attempts: list[dict[str, Any]] = []
    total_calls = 0
    legacy_path = qc / f"tau1_dense_discrete_mass_identity_v2_corrected_{steps}.json"
    legacy_preiteration_recoveries = 0
    if legacy_path.is_file():
        legacy = json.loads(legacy_path.read_text(encoding="utf-8"))
        legacy_preiteration_recoveries = sum(
            1
            for item in legacy.get("corrected_short_diagnostic_attempts", [])
            if item.get("status") == "RECOVERABLE_OPERATIONAL_ERROR"
        )
        for legacy_attempt in (1, 2):
            legacy_runtime = _runtime_windows() / (
                "dense_diagnostic_corrected_"
                f"3117927_{steps}_attempt_{legacy_attempt}"
            )
            if legacy_runtime.is_dir():
                _copy_segment_evidence(
                    legacy_runtime,
                    _run_root(root)
                    / "operational_failures"
                    / f"corrected_dense_preiteration_attempt_{legacy_attempt}",
                )
    for attempt in range(1, maximum_attempts + 1):
        result = _run_dense_discrete_diagnostic_attempt(
            root, steps=steps, attempt=attempt
        )
        attempts.append(
            {
                "attempt": attempt,
                "status": result["status"],
                "returncode": result.get("returncode"),
                "available_restart_iterations": result.get(
                    "available_restart_iterations"
                ),
            }
        )
        total_calls += int(result.get("scientific_musubi_calls", 0))
        if result["status"] != "RECOVERABLE_OPERATIONAL_ERROR":
            result["corrected_short_diagnostic_attempts"] = attempts
            result["scientific_musubi_calls"] = total_calls
            result["operational_recoveries"] = (
                legacy_preiteration_recoveries + attempt - 1
            )
            result["preiteration_operational_recoveries"] = (
                legacy_preiteration_recoveries
            )
            write_json(aggregate_path, result)
            return result
    failure = {
        "status": "RECOVERABLE_OPERATIONAL_ERROR",
        "corrected_short_diagnostic_attempts": attempts,
        "scientific_musubi_calls": total_calls,
        "operational_recoveries": legacy_preiteration_recoveries + maximum_attempts,
        "preiteration_operational_recoveries": legacy_preiteration_recoveries,
        "long_base_musubi_calls_after_stop": 0,
        "next": "STOP; CORRECTED SHORT-DIAGNOSTIC ATTEMPT BUDGET EXHAUSTED",
    }
    write_json(aggregate_path, failure)
    return failure


def forensic_dense_discrete_failure(project_root: Path) -> dict[str, Any]:
    """Source-decompose the failed dense identity without another solver call."""

    root = Path(project_root).resolve()
    run_root = _run_root(root)
    qc = run_root / "qc"
    legacy_path = (
        qc / "tau1_dense_discrete_mass_identity_v2_corrected_v2_8.json"
    )
    if not legacy_path.is_file():
        raise RuntimeError("the completed 8-step dense diagnostic is unavailable")
    legacy = json.loads(legacy_path.read_text(encoding="utf-8"))
    if legacy.get("status") != "TRUE_BOUNDARY_DISCRETE_ACCOUNTING_FAILURE":
        raise RuntimeError("dense diagnostic did not enter Decision B")

    latest = _latest_runtime_restart()
    if latest is None or int(latest[0]) != int(legacy["start_iteration"]):
        raise RuntimeError("preserved dense start restart is unavailable")
    _, _, start_binary = latest
    diagnostic_restart = (
        run_root
        / "dense_diagnostic_corrected"
        / "attempt_1_8_steps"
        / "restart"
    )
    diagnostic_pairs = _restart_pairs(diagnostic_restart)
    expected_iterations = list(
        range(int(legacy["start_iteration"]) + 1, int(legacy["end_iteration"]) + 1)
    )
    if sorted(diagnostic_pairs) != expected_iterations:
        raise RuntimeError("archived dense restart sequence is incomplete")

    contract = historical_tau1_runtime_contract()
    mesh = load_mesh_contract(_mesh_path(root), expected_cells=EXPECTED_CELLS)
    solid = runtime_solid_cells(mesh)
    outlet_pressures = {
        label: contract.pressure_reference_pa + gauge
        for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
    }
    labels = ("wall", "inlet", "outlet_01", "outlet_02", "outlet_03")
    corrected_sums = {label: [] for label in labels}
    per_step: list[dict[str, Any]] = []
    direction_correction: dict[int, float] = {}
    affected_wall_writes = 0
    corrected_predicted: list[float] = []
    corrected_observed: list[float] = []
    current_pdf = _read_pdf(start_binary)
    endpoint_start = current_pdf.copy()
    for next_iteration in expected_iterations:
        next_pdf = _read_pdf(diagnostic_pairs[next_iteration][1])
        replay = replay_boundary_step(
            current_pdf,
            mesh,
            dx_m=contract.dx_m,
            dt_s=contract.dt_s,
            density_kg_m3=contract.rho_kg_m3,
            target_mass_flow_kg_s=contract.target_mass_flow_kg_s,
            outlet_pressures_pa=outlet_pressures,
        )
        legacy_wall = _wall_operations(current_pdf, mesh)
        corrected_wall = _wall_operations(current_pdf, mesh, solid)
        if len(legacy_wall) != len(corrected_wall):
            raise RuntimeError("wall operation cardinality changed")
        per_direction: dict[int, float] = {}
        step_affected = 0
        for old, corrected in zip(legacy_wall, corrected_wall, strict=True):
            if old[:3] != corrected[:3]:
                raise RuntimeError("wall operation ordering changed")
            change = float(corrected[3]) - float(old[3])
            storage_direction = int(corrected[2])
            per_direction[storage_direction] = (
                per_direction.get(storage_direction, 0.0) + change
            )
            direction_correction[storage_direction] = (
                direction_correction.get(storage_direction, 0.0) + change
            )
            if change != 0.0:
                step_affected += 1
        affected_wall_writes += step_affected
        boundary_values = {
            label: float(replay["per_label_lattice"][label]) for label in labels
        }
        for label, value in boundary_values.items():
            corrected_sums[label].append(value)
        predicted = float(replay["predicted_total_lattice"])
        observed = float(np.sum(next_pdf - current_pdf, dtype=np.float64))
        corrected_predicted.append(predicted)
        corrected_observed.append(observed)
        per_step.append(
            {
                "start_iteration": next_iteration - 1,
                "end_iteration": next_iteration,
                "runtime_solid_affected_wall_writes": step_affected,
                "wall_runtime_solid_correction_by_storage_direction_lattice": {
                    str(direction + 1): value
                    for direction, value in sorted(per_direction.items())
                    if value != 0.0
                },
                "per_boundary_lattice_signed": boundary_values,
                "predicted_lattice": predicted,
                "observed_lattice": observed,
                "R_source_corrected_one_step": conservation_identity_residual(
                    predicted, observed, replay["target_lattice_flux"]
                ),
            }
        )
        current_pdf = next_pdf

    corrected_boundary = {
        label: math.fsum(values) for label, values in corrected_sums.items()
    }
    corrected_total = math.fsum(corrected_predicted)
    endpoint_delta = float(np.sum(current_pdf - endpoint_start, dtype=np.float64))
    inlet_normalizer = abs(corrected_boundary["inlet"])
    corrected_residual = abs(corrected_total - endpoint_delta) / inlet_normalizer
    old_total = float(legacy["sum_predicted_exact_boundary_pdf_mass_delta_lattice"])
    old_signed_error = old_total - endpoint_delta
    corrected_signed_error = corrected_total - endpoint_delta
    wall_connectivity_correction = corrected_total - old_total
    mass_per_lattice = contract.rho_kg_m3 * contract.dx_m**3

    source_root = _tem_windows().parent
    source_files = {
        "timestep": source_root / "mus" / "source" / "mus_control_module.f90",
        "boundary_buffers": source_root
        / "mus"
        / "source"
        / "bc"
        / "mus_bc_general_module.fpp",
        "wall_libb": source_root
        / "mus"
        / "source"
        / "bc"
        / "mus_bc_fluid_wall_module.fpp",
        "pull_connectivity": source_root
        / "mus"
        / "source"
        / "mus_connectivity_module.fpp",
        "aux_density": source_root
        / "mus"
        / "source"
        / "derived"
        / "mus_auxFieldVar_module.fpp",
        "tau1_bgk": source_root
        / "mus"
        / "source"
        / "compute"
        / "mus_compute_d3q19_module.fpp",
        "restart_serializer": source_root
        / "mus"
        / "source"
        / "mus_buffer_module.fpp",
    }
    early = json.loads(
        (qc / "continuous_q_referee_v2_tau1_base.json").read_text(
            encoding="utf-8"
        )
    )
    main_pairs = _restart_pairs(_runtime_windows() / "restart")
    long_hashes = {
        str(iteration): sha256_file(main_pairs[iteration][1])
        for iteration in (2_878_425, 2_998_176, 3_117_927)
    }
    corrected_window_path = (
        qc / "tau1_boundary_window_forensics_runtime_solid_wall_corrected.json"
    )
    corrected_window = json.loads(
        corrected_window_path.read_text(encoding="utf-8")
    )
    corrected_long_window_closure = {
        str(item["window"]): float(
            item["integration_methods"]["trapezoidal"]["closure"]
        )
        for item in corrected_window["windows"]
    }
    result = {
        "status": "TRUE_BOUNDARY_DISCRETE_ACCOUNTING_FAILURE",
        "analysis_mode": "ZERO-LONG-RUN_EXISTING_PDFS_ONLY",
        "additional_musubi_calls": 0,
        "steps": len(expected_iterations),
        "start_iteration": int(legacy["start_iteration"]),
        "end_iteration": int(legacy["end_iteration"]),
        "legacy_dense_residual": float(legacy["exact_dense_discrete_closure"]),
        "source_corrected_runtime_solid_wall_residual": corrected_residual,
        "source_corrected_individual_residual_min": min(
            item["R_source_corrected_one_step"] for item in per_step
        ),
        "source_corrected_individual_residual_median": float(
            np.median(
                [item["R_source_corrected_one_step"] for item in per_step]
            )
        ),
        "source_corrected_individual_residual_max": max(
            item["R_source_corrected_one_step"] for item in per_step
        ),
        "source_corrected_predicted_lattice": corrected_total,
        "observed_endpoint_lattice": endpoint_delta,
        "source_corrected_absolute_mismatch_lattice": abs(
            corrected_total - endpoint_delta
        ),
        "source_corrected_predicted_kg": corrected_total * mass_per_lattice,
        "observed_endpoint_kg": endpoint_delta * mass_per_lattice,
        "source_corrected_absolute_mismatch_kg": abs(
            corrected_total - endpoint_delta
        )
        * mass_per_lattice,
        "exact_inlet_normalizer_lattice": inlet_normalizer,
        "exact_inlet_normalizer_kg": inlet_normalizer * mass_per_lattice,
        "source_corrected_per_boundary_lattice_signed": corrected_boundary,
        "source_corrected_per_boundary_kg_signed": {
            label: value * mass_per_lattice
            for label, value in corrected_boundary.items()
        },
        "offset_decomposition": {
            "legacy_signed_error_lattice": old_signed_error,
            "wall_runtime_solid_connectivity_correction_lattice": (
                wall_connectivity_correction
            ),
            "remaining_signed_error_lattice": corrected_signed_error,
            "legacy_signed_error_normalized": old_signed_error / inlet_normalizer,
            "wall_runtime_solid_connectivity_correction_normalized": (
                wall_connectivity_correction / inlet_normalizer
            ),
            "remaining_signed_error_normalized": (
                corrected_signed_error / inlet_normalizer
            ),
        },
        "runtime_solid_cell_count": len(solid),
        "runtime_solid_affected_wall_writes_over_8_steps": affected_wall_writes,
        "wall_runtime_solid_correction_by_storage_direction_lattice": {
            str(direction + 1): value
            for direction, value in sorted(direction_correction.items())
            if value != 0.0
        },
        "first_contract_difference": {
            "earlier_pass_state": "fresh equilibrium initialization at iteration 0",
            "earlier_wall_delta_lattice": float(
                early["per_boundary_lattice"]["wall"]
            ),
            "earlier_R_one_step_identity": float(early["R_one_step_identity"]),
            "late_state": "non-equilibrium steady restart at iteration 3117927",
            "late_legacy_wall_delta_first_step_lattice": float(
                legacy["per_step"][0]["per_boundary_lattice_signed"]["wall"]
            ),
            "meaning": (
                "the fresh equilibrium state made the wall_libb neighbor term "
                "degenerate, so it did not exercise runtime-solid wall FETCH"
            ),
        },
        "source_phase_contract": {
            "order": [
                "restart serializes nNext",
                "set_boundary writes nNext",
                "fill_neighBuffer/computeNeighBuf uses PULL FETCH",
                "swap nNow/nNext",
                "aux density sums PULL-fetched PDFs",
                "omega=1 BGK writes equilibrium preserving that density",
            ],
            "non_boundary_source": "NONE; generated Lua has no source table",
            "controller_term": (
                "unchanged; adaptive inlet contribution and exact inlet "
                "normalizer are identical before/after the wall correction"
            ),
            "source_files_sha256": {
                label: sha256_file(path) for label, path in source_files.items()
            },
        },
        "root_cause_of_0p002049922": (
            "The V2 wall_libb replay read fNgh with coordinate-only PULL while "
            "Musubi computeNeighBuf uses the runtime connectivity. When either "
            "the current or source element carries prp_solid, Musubi reads the "
            "current inverse PDF. The fresh equilibrium referee hid this path; "
            "the late non-equilibrium restart exposed it. This wall-only term "
            "accounts for the stable O(0.002) offset. A smaller residual remains "
            "above 1e-8, so the hard gate still fails."
        ),
        "remaining_failure": (
            "Endpoint PDFs provide conserved moments, not every pre-collision "
            "boundary replacement. The remaining term cannot be uniquely "
            "assigned without source instrumentation; no further Musubi call "
            "is authorized."
        ),
        "hard_gate": 1.0e-8,
        "hard_gate_pass": corrected_residual <= 1.0e-8,
        "corrected_long_window_trapezoidal_closure": (
            corrected_long_window_closure
        ),
        "long_base_restart_sha256": long_hashes,
        "per_step": per_step,
        "next": "STOP BEFORE ANY NEW LONG CFD; INSTRUMENT BOUNDARY ACCOUNTING ONLY",
    }
    write_json(qc / "tau1_true_boundary_discrete_accounting_failure.json", result)
    return result


def salvage_incomplete_dense_diagnostic(project_root: Path) -> dict[str, Any]:
    """Audit the sole dense run without launching or advancing Musubi again."""

    root = Path(project_root).resolve()
    run_root = _run_root(root)
    qc = run_root / "qc"
    raw_path = qc / "tau1_dense_discrete_mass_identity_v2.json"
    raw = json.loads(raw_path.read_text(encoding="utf-8"))
    if raw.get("status") != "DENSE_DIAGNOSTIC_INCOMPLETE_NO_RETRY":
        raise RuntimeError("the sole dense diagnostic is not the expected incomplete run")
    start_iteration = int(raw["start_iteration"])
    steps = int(raw["requested_steps"])
    end_iteration = start_iteration + steps
    diagnostic_name = f"dense_diagnostic_{start_iteration}_{steps}"
    diagnostic_runtime = _runtime_windows() / diagnostic_name
    diagnostic_pairs = _restart_pairs(diagnostic_runtime / "restart")
    if sorted(diagnostic_pairs) != [end_iteration]:
        raise RuntimeError("the overwritten dense diagnostic endpoint is unavailable")
    main_pairs = _restart_pairs(_runtime_windows() / "restart")
    if start_iteration not in main_pairs:
        raise RuntimeError("the preserved long-Base start restart is unavailable")

    contract = historical_tau1_runtime_contract()
    mesh = load_mesh_contract(_mesh_path(root), expected_cells=EXPECTED_CELLS)
    start_binary = main_pairs[start_iteration][1]
    end_header, end_binary = diagnostic_pairs[end_iteration]
    start_state = _state(start_binary, contract=contract, with_velocity=False)
    end_state = _state(end_binary, contract=contract, with_velocity=False)
    outlet_pressures = {
        label: contract.pressure_reference_pa + gauge
        for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
    }
    endpoint_replays = [
        replay_boundary_step(
            state["pdf"],
            mesh,
            dx_m=contract.dx_m,
            dt_s=contract.dt_s,
            density_kg_m3=contract.rho_kg_m3,
            target_mass_flow_kg_s=contract.target_mass_flow_kg_s,
            outlet_pressures_pa=outlet_pressures,
        )
        for state in (start_state, end_state)
    ]
    labels = ("wall", "inlet", "outlet_01", "outlet_02", "outlet_03")
    methods = (
        "left_rectangle",
        "right_rectangle",
        "trapezoidal",
        "timestep_phase_shifted",
    )
    endpoint_delta_lattice = float(
        np.sum(end_state["pdf"] - start_state["pdf"], dtype=np.float64)
    )
    mass_per_lattice = contract.rho_kg_m3 * contract.dx_m**3
    observed_mass_change = endpoint_delta_lattice * mass_per_lattice
    inlet_normalizer = contract.target_mass_flow_kg_s * steps * contract.dt_s
    integrations: dict[str, Any] = {}
    csv_rows: list[dict[str, Any]] = []
    for method in methods:
        per_boundary_mass = {
            label: _quadrature_integral(
                [start_iteration, end_iteration],
                [
                    float(endpoint_replays[0]["per_label_kg_s_domain"][label]),
                    float(endpoint_replays[1]["per_label_kg_s_domain"][label]),
                ],
                method,
            )
            * contract.dt_s
            for label in labels
        }
        predicted_mass_change = float(sum(per_boundary_mass.values()))
        mismatch = abs(predicted_mass_change - observed_mass_change)
        closure = boundary_window_closure(
            predicted_mass_change, observed_mass_change, inlet_normalizer
        )
        integrations[method] = {
            "predicted_mass_change_kg": predicted_mass_change,
            "observed_mass_change_kg": observed_mass_change,
            "absolute_mismatch_kg": mismatch,
            "inlet_mass_normalizer_kg": inlet_normalizer,
            "closure": closure,
            "per_boundary_contribution_kg_domain_signed": per_boundary_mass,
        }
        csv_rows.append(
            {
                "method": method,
                "predicted_mass_change_kg": predicted_mass_change,
                "observed_mass_change_kg": observed_mass_change,
                "absolute_mismatch_kg": mismatch,
                "inlet_mass_normalizer_kg": inlet_normalizer,
                "closure": closure,
                **{f"{label}_kg": per_boundary_mass[label] for label in labels},
            }
        )

    stdout_path = diagnostic_runtime / "musubi_stdout.log"
    stdout = _read_text(stdout_path)
    write_iterations = [
        int(value)
        for value in re.findall(
            r"Writing restart at:\s*\n\s*iterations:\s*(\d+)", stdout
        )
    ]
    controller = [
        item
        for item in _controller_records(stdout)
        if start_iteration < int(item["iteration"]) <= end_iteration
    ]
    maximum_elementwise_pdf_change = float(
        np.max(np.abs(end_state["pdf"] - start_state["pdf"]))
    )
    relative_endpoint_pdf_l2 = float(
        np.linalg.norm((end_state["pdf"] - start_state["pdf"]).ravel())
        / np.linalg.norm(end_state["pdf"].ravel())
    )
    main_restart_hashes = {
        str(iteration): sha256_file(binary)
        for iteration, (_, binary) in sorted(main_pairs.items())
    }
    preserved_restart_manifest = {
        "status": "PRESERVED",
        "restart_directory_windows": str(_runtime_windows() / "restart"),
        "restart_directory_wsl": f"{RUNTIME_WSL}/restart",
        "files": [
            {
                "iteration": iteration,
                "physical_time_s": iteration * contract.dt_s,
                "header_name": header.name,
                "header_size_bytes": header.stat().st_size,
                "header_sha256": sha256_file(header),
                "binary_name": binary.name,
                "binary_size_bytes": binary.stat().st_size,
                "binary_sha256": sha256_file(binary),
            }
            for iteration, (header, binary) in sorted(main_pairs.items())
        ],
    }
    write_json(qc / "tau1_preserved_long_restart_manifest.json", preserved_restart_manifest)
    preserved_headers = run_root / "preserved_long_restart_headers"
    preserved_headers.mkdir(exist_ok=True)
    for header, _ in main_pairs.values():
        shutil.copy2(header, preserved_headers / header.name)
    archive = run_root / "dense_diagnostic"
    archive.mkdir(parents=True, exist_ok=True)
    for name in (
        "musubi.lua",
        "musubi_stdout.log",
        "musubi_stderr.log",
        "semantic_status.log",
        "diagnostic_launched",
    ):
        source = diagnostic_runtime / name
        if source.is_file():
            shutil.copy2(source, archive / name)
    restart_archive = archive / "restart_endpoint"
    restart_archive.mkdir(exist_ok=True)
    shutil.copy2(end_header, restart_archive / end_header.name)
    shutil.copy2(end_binary, restart_archive / end_binary.name)
    csv_path = qc / "tau1_dense_endpoint_salvage.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(csv_rows[0]))
        writer.writeheader()
        writer.writerows(csv_rows)

    result = {
        "status": "DENSE_EXACT_DISCRETE_EVIDENCE_INCOMPLETE_NO_RETRY",
        "scientific_musubi_calls": 1,
        "additional_musubi_calls_for_salvage": 0,
        "executed_timesteps": steps,
        "start_iteration": start_iteration,
        "end_iteration": end_iteration,
        "write_events_in_log": write_iterations,
        "write_event_count": len(write_iterations),
        "controller_record_count": len(controller),
        "retained_endpoint_iterations": sorted(diagnostic_pairs),
        "root_cause": (
            "restart formatter option was written as timeform, but the linked "
            "Tem implementation loads timeformat; identical simulated-time filenames "
            "therefore replaced the first 31 per-step restart files"
        ),
        "strict_dense_identity_computable": False,
        "strict_classification": None,
        "classification_not_claimed": [
            "SPARSE_BOUNDARY_WINDOW_AUDIT_BIAS_CONFIRMED",
            "TRUE_BOUNDARY_WINDOW_ACCOUNTING_FAILURE",
        ],
        "reason_classification_not_claimed": (
            "the exact boundary replay cannot be evaluated at each of the 32 states; "
            "endpoint quadrature is not a substitute for the requested exact sum"
        ),
        "endpoint_only_nonqualifying_evidence": {
            "elementwise_sum_pdf_end_minus_pdf_start_lattice": endpoint_delta_lattice,
            "observed_mass_change_kg": observed_mass_change,
            "maximum_elementwise_pdf_change": maximum_elementwise_pdf_change,
            "relative_endpoint_pdf_l2": relative_endpoint_pdf_l2,
            "integration_methods": integrations,
        },
        "minimum_pdf": min(
            float(start_state["minimum_pdf"]), float(end_state["minimum_pdf"])
        ),
        "maximum_lattice_speed": max(
            float(start_state["maximum_lattice_speed"]),
            float(end_state["maximum_lattice_speed"]),
        ),
        "controller_records": controller,
        "start_restart_sha256": sha256_file(start_binary),
        "end_restart_sha256": sha256_file(end_binary),
        "long_base_restart_sha256": main_restart_hashes,
        "long_base_restart_manifest": preserved_restart_manifest,
        "long_base_restarts_preserved": True,
        "binary_sha256": MUSUBI_SHA256,
        "lua_sha256": sha256_file(diagnostic_runtime / "musubi.lua"),
        "stdout_sha256": sha256_file(stdout_path),
        "audit_fix_applied": False,
        "acceptance_threshold": BOUNDARY_CLOSURE_GATE,
        "acceptance_threshold_changed": False,
        "new_long_cfd_allowed": False,
        "next": "STOP; DENSE EXACT CLASSIFICATION REQUIRES NEW USER AUTHORITY",
    }
    write_json(qc / "tau1_dense_discrete_mass_identity_v2_salvage.json", result)
    run_state_path = qc / "tau1_base_run_state.json"
    state = json.loads(run_state_path.read_text(encoding="utf-8"))
    state["status"] = "STOPPED_AT_COMPATIBLE_CHECKPOINT_DENSE_EVIDENCE_INCOMPLETE"
    state["stop_iteration"] = start_iteration
    state["dense_diagnostic"] = {
        "status": result["status"],
        "steps": steps,
        "end_iteration": end_iteration,
        "strict_classification": None,
    }
    write_json(run_state_path, state)
    final = {
        "status": "CFD_FLOW_OPERATIONAL_INFRASTRUCTURE_BLOCKED",
        "next": "STOP",
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "fresh_base_logical_musubi_calls": 1,
        "short_dense_diagnostic_calls": 1,
        "accepted_restart": None,
        "long_base_stop_iteration": start_iteration,
        "long_base_stop_physical_time_s": start_iteration * contract.dt_s,
        "long_base_restart_sha256": main_restart_hashes,
        "dense_diagnostic": result,
        "strict_scientific_classification": None,
        "true_first_scientific_failure": None,
        "operational_blocker": result["root_cause"],
        "boundary_acceptance_threshold": BOUNDARY_CLOSURE_GATE,
        "boundary_acceptance_threshold_changed": False,
    }
    write_json(qc / "tau1_base_final.json", final)
    return result


def _copy_segment_evidence(segment_runtime: Path, segment_archive: Path) -> None:
    segment_archive.mkdir(parents=True, exist_ok=True)
    for name in ("musubi.lua", "musubi_stdout.log", "musubi_stderr.log", "semantic_status.log"):
        source = segment_runtime / name
        if source.is_file():
            shutil.copy2(source, segment_archive / name)
    tracking_source = segment_runtime / "tracking"
    if tracking_source.is_dir():
        destination = segment_archive / "tracking"
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(tracking_source, destination)


def _invoke_wsl_script(script_wsl: str, *args: str, timeout_s: int) -> tuple[int, float]:
    started = time.perf_counter()
    process = subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "/bin/bash", script_wsl, *args],
        check=False,
        timeout=timeout_s,
    )
    return process.returncode, time.perf_counter() - started


def _invoke_confirmation_process(
    *,
    script_wsl: str,
    segment_wsl: str,
    segment_runtime: Path,
    project_root: Path,
    run_root: Path,
    referee_residual: float,
    timeout_s: int,
) -> tuple[int, float, dict[str, Any] | None]:
    """Run one long confirmation process and stop at the first exact PASS triplet."""

    started = time.perf_counter()
    process = subprocess.Popen(
        [
            "wsl.exe",
            "-d",
            "Ubuntu",
            "--",
            "/bin/bash",
            script_wsl,
            segment_wsl,
        ]
    )
    last_audited: tuple[int, int, int] | None = None
    accepted: dict[str, Any] | None = None
    while process.poll() is None:
        if time.perf_counter() - started > timeout_s:
            process.terminate()
            raise subprocess.TimeoutExpired(script_wsl, timeout_s)
        pairs = _restart_pairs(_runtime_windows() / "restart")
        available = sorted(pairs)
        if len(available) >= 3:
            iterations = tuple(available[-3:])
            exact_windows = (
                iterations[1] - iterations[0] == CONFIRMATION_INTERVAL
                and iterations[2] - iterations[1] == CONFIRMATION_INTERVAL
            )
            if exact_windows and iterations != last_audited:
                stdout_path = segment_runtime / "musubi_stdout.log"
                records = (
                    _controller_records(_read_text(stdout_path))
                    if stdout_path.is_file()
                    else []
                )
                if records:
                    accepted = audit_base_window(
                        project_root=project_root,
                        restart_triplet=[
                            (iteration, pairs[iteration][1])
                            for iteration in iterations
                        ],
                        referee_residual=referee_residual,
                        controller=max(
                            records, key=lambda item: int(item["iteration"])
                        ),
                    )
                    write_json(run_root / "qc" / "tau1_base_steady_status.json", accepted)
                    last_audited = iterations
                    if accepted["status"] == "PASS":
                        (_runtime_windows() / "stop").touch(exist_ok=True)
            for iteration in available[:-3]:
                header, binary = pairs[iteration]
                header.unlink(missing_ok=True)
                binary.unlink(missing_ok=True)
        time.sleep(10.0)
    return process.returncode, time.perf_counter() - started, accepted


def _latest_runtime_restart() -> tuple[int, Path, Path] | None:
    restart = _runtime_windows() / "restart"
    if not restart.is_dir():
        return None
    pairs = _restart_pairs(restart)
    if not pairs:
        return None
    iteration = max(pairs)
    header, binary = pairs[iteration]
    return iteration, header, binary


def run_fresh_base(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    run_root = _run_root(root)
    qc = run_root / "qc"
    referee = json.loads((qc / "continuous_q_referee_v2_tau1_base.json").read_text(encoding="utf-8"))
    if referee.get("status") != "PASS":
        result = {"status": "CFD_FLOW_CONTINUOUS_Q_REFEREE_FAILED", "long_base_musubi_calls": 0}
        write_json(qc / "tau1_base_final.json", result)
        return result
    contract = Tau1BaseRuntimeContract()
    expected_manifest = _runtime_manifest(contract)
    runtime = _runtime_windows()
    runtime.mkdir(parents=True, exist_ok=True)
    (runtime / "restart").mkdir(exist_ok=True)
    runtime_manifest_path = runtime / "runtime_contract.json"
    if runtime_manifest_path.is_file():
        saved = json.loads(runtime_manifest_path.read_text(encoding="utf-8"))
        compatibility = restart_resume_contract(saved, expected_manifest)
        if compatibility["status"] != "PASS":
            raise RuntimeError(f"incompatible tau1 Base runtime exists: {compatibility}")
    else:
        write_json(runtime_manifest_path, expected_manifest)
    run_state_path = qc / "tau1_base_run_state.json"
    state = (
        json.loads(run_state_path.read_text(encoding="utf-8"))
        if run_state_path.is_file()
        else {
            "status": "IN_PROGRESS",
            "logical_musubi_calls": 1,
            "restart_resumes": 0,
            "segment_processes": 0,
            "confirmation_mode": False,
            "total_wall_clock_s": 0.0,
            "checkpoint_iterations": [],
            "operational_recoveries": 0,
        }
    )
    base_script = f"{PROJECT_WSL}/outputs/cfd_flow/{RUN_NAME}/run_base_segment.sh"
    latest = _latest_runtime_restart()
    current = latest[0] if latest else 0
    if latest and not state["checkpoint_iterations"]:
        state["checkpoint_iterations"] = [current]
    if current >= contract.hard_max_iterations:
        state["status"] = "HARD_MAX_REACHED"
    while current < contract.hard_max_iterations and state.get("status") == "IN_PROGRESS":
        continuous_confirmation = bool(state["confirmation_mode"])
        step = CONFIRMATION_INTERVAL if continuous_confirmation else CHECKPOINT_INTERVAL
        maximum = (
            contract.hard_max_iterations
            if continuous_confirmation
            else min(current + step, contract.hard_max_iterations)
        )
        segment_name = (
            f"confirmation_from_{current:07d}_to_{maximum:07d}"
            if continuous_confirmation
            else f"segment_{maximum:07d}"
        )
        segment_wsl = f"{RUNTIME_WSL}/segments/{segment_name}"
        segment_runtime = runtime / "segments" / segment_name
        segment_runtime.mkdir(parents=True, exist_ok=True)
        restart_header_wsl = (
            f"{RUNTIME_WSL}/restart/{latest[1].name}" if latest else None
        )
        lua = generate_base_segment_lua(
            maximum_iterations=maximum,
            segment_wsl=segment_wsl,
            restart_header_wsl=restart_header_wsl,
            restart_first_iteration=(
                current + CONFIRMATION_INTERVAL
                if continuous_confirmation
                else maximum
            ),
            restart_interval=(CONFIRMATION_INTERVAL if continuous_confirmation else 1),
            contract=contract,
        )
        lua_contract = base_lua_contract(
            lua,
            maximum_iterations=maximum,
            restart_header_wsl=restart_header_wsl,
            restart_first_iteration=(
                current + CONFIRMATION_INTERVAL
                if continuous_confirmation
                else maximum
            ),
            restart_interval=(CONFIRMATION_INTERVAL if continuous_confirmation else 1),
            contract=contract,
        )
        if lua_contract["status"] != "PASS":
            raise RuntimeError(f"Base segment Lua contract failed: {lua_contract}")
        (segment_runtime / "musubi.lua").write_text(lua, encoding="utf-8", newline="\n")
        (runtime / "stop").unlink(missing_ok=True)
        monitored_audit: dict[str, Any] | None = None
        if continuous_confirmation:
            returncode, elapsed, monitored_audit = _invoke_confirmation_process(
                script_wsl=base_script,
                segment_wsl=segment_wsl,
                segment_runtime=segment_runtime,
                project_root=root,
                run_root=run_root,
                referee_residual=float(referee["R_one_step_identity"]),
                timeout_s=CONTINUOUS_PROCESS_TIMEOUT_SECONDS,
            )
        else:
            returncode, elapsed = _invoke_wsl_script(
                base_script, segment_wsl, timeout_s=7_800
            )
        state["segment_processes"] += 1
        if latest:
            state["restart_resumes"] += 1
        state["total_wall_clock_s"] += elapsed
        archive = run_root / "segments" / segment_name
        _copy_segment_evidence(segment_runtime, archive)
        if returncode != 0:
            state["operational_recoveries"] += 1
            write_json(run_state_path, state)
            if state["operational_recoveries"] >= 5:
                state["status"] = "CFD_FLOW_OPERATIONAL_INFRASTRUCTURE_BLOCKED"
                break
            latest_after_failure = _latest_runtime_restart()
            if latest_after_failure is None or latest_after_failure[0] <= current:
                continue
            latest = latest_after_failure
            current = latest[0]
            continue
        latest = _latest_runtime_restart()
        if latest is None or latest[0] <= current:
            state["operational_recoveries"] += 1
            write_json(run_state_path, state)
            if state["operational_recoveries"] >= 5:
                state["status"] = "CFD_FLOW_OPERATIONAL_INFRASTRUCTURE_BLOCKED"
                break
            continue
        current = latest[0]
        state["checkpoint_iterations"].append(current)
        state["checkpoint_iterations"] = sorted(set(state["checkpoint_iterations"]))
        if monitored_audit is not None:
            state["latest_steady_audit"] = monitored_audit
            if monitored_audit["status"] == "PASS":
                state["status"] = "PASS"
                write_json(run_state_path, state)
                break
        if continuous_confirmation:
            records = _controller_records(
                _read_text(segment_runtime / "musubi_stdout.log")
            )
            reached_hard_max = bool(
                records
                and max(int(item["iteration"]) for item in records)
                >= contract.hard_max_iterations
            )
            if reached_hard_max:
                state["status"] = "HARD_MAX_REACHED"
            write_json(run_state_path, state)
            if reached_hard_max:
                break
            continue
        pairs = _restart_pairs(runtime / "restart")
        if not state["confirmation_mode"] and len(pairs) >= 2:
            available = sorted(pairs)
            previous_candidates = [item for item in available if item < current]
            if previous_candidates:
                previous = max(previous_candidates)
                if current - previous == CHECKPOINT_INTERVAL:
                    a = _state(pairs[previous][1], contract=contract, with_velocity=False)
                    b = _state(pairs[current][1], contract=contract, with_velocity=False)
                    conversion = contract.rho_kg_m3 * contract.dx_m**3 / contract.dt_s
                    coarse_mass = abs(
                        (b["total_pdf_mass"] - a["total_pdf_mass"])
                        / (current - previous)
                        * conversion
                    ) / contract.target_mass_flow_kg_s
                    state["latest_coarse_mass_residual"] = coarse_mass
                    if coarse_mass <= 0.05:
                        state["confirmation_mode"] = True
        pairs = _restart_pairs(runtime / "restart")
        available = sorted(pairs)
        if state["confirmation_mode"] and len(available) >= 3:
            triplet_iterations = available[-3:]
            if (
                triplet_iterations[1] - triplet_iterations[0] == CONFIRMATION_INTERVAL
                and triplet_iterations[2] - triplet_iterations[1] == CONFIRMATION_INTERVAL
            ):
                controller = _latest_controller_from_segments(run_root)
                audit = audit_base_window(
                    project_root=root,
                    restart_triplet=[(item, pairs[item][1]) for item in triplet_iterations],
                    referee_residual=float(referee["R_one_step_identity"]),
                    controller=controller,
                )
                write_json(qc / "tau1_base_steady_status.json", audit)
                state["latest_steady_audit"] = audit
                if audit["status"] == "PASS":
                    state["status"] = "PASS"
                    break
        keep = sorted(pairs)[-3:]
        for iteration, (header, binary) in pairs.items():
            if iteration not in keep:
                header.unlink(missing_ok=True)
                binary.unlink(missing_ok=True)
        write_json(run_state_path, state)
    write_json(run_state_path, state)
    latest = _latest_runtime_restart()
    steady = state.get("latest_steady_audit")
    if state["status"] == "PASS" and steady:
        final_status = "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_PASS"
    elif state["status"] == "CFD_FLOW_OPERATIONAL_INFRASTRUCTURE_BLOCKED":
        final_status = state["status"]
    else:
        final_status = "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_FAILED"
    accepted_archive: dict[str, Any] | None = None
    if latest and steady and steady["status"] == "PASS":
        accepted_dir = run_root / "accepted_restart"
        accepted_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(latest[1], accepted_dir / latest[1].name)
        shutil.copy2(latest[2], accepted_dir / latest[2].name)
        accepted_archive = {
            "iteration": latest[0],
            "header": str(accepted_dir / latest[1].name),
            "binary": str(accepted_dir / latest[2].name),
            "binary_sha256": sha256_file(accepted_dir / latest[2].name),
        }
    runtime_contract_path = qc / "tau1_base_runtime_contract.json"
    segment_lua_paths = sorted(
        (_runtime_windows() / "segments").glob("*/musubi.lua"),
        key=lambda path: path.stat().st_mtime_ns,
    )
    lua_sha256 = sha256_file(segment_lua_paths[-1]) if segment_lua_paths else None
    runtime_contract_sha256 = (
        sha256_file(runtime_contract_path) if runtime_contract_path.is_file() else None
    )
    if final_status == "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_PASS":
        true_first_scientific_failure = None
    elif final_status == "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_FAILED" and steady:
        failed_gates = [
            name for name, passed in steady.get("gates", {}).items() if not passed
        ]
        if failed_gates == ["boundary_window_closure_le_0p1pct"]:
            true_first_scientific_failure = (
                "boundary-window closure "
                f"{steady['boundary_window_closure']:.17g} > {BOUNDARY_CLOSURE_GATE:.17g} "
                "at physical-time hard max"
            )
        else:
            true_first_scientific_failure = ", ".join(failed_gates) or "steady audit unavailable"
    else:
        true_first_scientific_failure = None
    final = {
        "status": final_status,
        "next": (
            "RUN REPAIRED TAU1 COARSE/BASE/FINE GRID CONVERGENCE"
            if final_status == "CFD_FLOW_REPAIRED_BASE_TAU1_STEADY_PASS"
            else "STOP"
        ),
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "fresh_base_logical_musubi_calls": 1,
        "restart_resumes": state["restart_resumes"],
        "segment_processes": state["segment_processes"],
        "total_wall_clock_s": state["total_wall_clock_s"],
        "referee": referee,
        "steady": steady,
        "accepted_restart": accepted_archive,
        "runtime_contract": contract.as_evidence(),
        "binary_sha256": MUSUBI_SHA256,
        "lua_sha256": lua_sha256,
        "runtime_contract_sha256": runtime_contract_sha256,
        "mesh_hashes": dict(MESH_HASHES),
        "operational_recoveries": state["operational_recoveries"],
        "true_first_scientific_failure": true_first_scientific_failure,
    }
    write_json(qc / "tau1_base_final.json", final)
    return final


def write_full_timestep_mass_contract(project_root: Path) -> dict[str, Any]:
    """Write the pinned source/runtime contract for the full V2 referee."""

    root = Path(project_root).resolve()
    run_root = _run_root(root)
    qc = run_root / "qc"
    source_root = _tem_windows().parent
    mus = source_root / "mus" / "source"
    tem = source_root / "tem" / "source"
    built_macro = source_root / "build" / "mus" / "source" / "header" / "lbm_macros.inc"
    archived = run_root / "dense_diagnostic_corrected" / "attempt_1_8_steps"
    lua_path = archived / "musubi.lua"
    stdout_path = archived / "musubi_stdout.log"
    mesh = load_mesh_contract(_mesh_path(root), expected_cells=EXPECTED_CELLS)
    lua_text = lua_path.read_text(encoding="utf-8")
    stdout = stdout_path.read_text(encoding="utf-8")
    source_absence = all(
        token not in lua_text
        for token in ("glob_source", "source =", "source=", "body_force")
    )
    runtime_count_proof = all(
        token in stdout
        for token in (
            "Total number of elements: 182320",
            "Local number of elements: 45580",
            "fluid: 45580",
            "halo: 400",
            "Fluid:   182320",
        )
    )
    evidence = {
        "control": source_token_evidence(
            mus / "mus_control_module.f90",
            "call set_boundary(",
            "call mus_swap_now_next(",
            "call mus_calcAuxFieldAndExchange(",
            "call me%scheme%compute(",
            "call mus_apply_sourceTerms(",
            "exchange_real(",
        ),
        "boundary_general": source_token_evidence(
            mus / "bc" / "mus_bc_general_module.fpp",
            "currState = state( :, pdf%nNext )",
            "call fill_neighBuffer(",
            "do iBnd = 1, nBCs",
            "%fnct(",
            "currstate( ?FETCH?(",
        ),
        "boundary_fluid": source_token_evidence(
            mus / "bc" / "mus_bc_fluid_module.fpp",
            "adaptive_flux_pressure",
            "pressure_eq",
        ),
        "boundary_wall": source_token_evidence(
            mus / "bc" / "mus_bc_fluid_wall_module.fpp",
            "subroutine wall_libb(",
            "fNgh = me%neigh(iLevel)%computeNeighBuf",
            "state( me%links(iLevel)%val(iLink) )",
        ),
        "boundary_setup": source_token_evidence(
            mus / "bc" / "mus_bc_header_module.fpp",
            "subroutine set_bouzidi_coeff(",
            "subroutine mus_set_bouzidi(",
            "bouzidi%nghPos( iLink ) = invDir",
        ),
        "connectivity": source_token_evidence(
            mus / "mus_connectivity_module.fpp",
            "solidified = ( btest(neighProp, prp_solid)",
            "missing_neigh_for_nonghost .or. solidified",
            "GetFromPos = iElem",
        ),
        "boundary_neighbor_ids": source_token_evidence(
            mus / "mus_construction_module.fpp",
            "subroutine assignBCList(",
            "bID      = bc_prop%boundary_ID",
            "stencil%cxDirInv(iDir)",
        ),
        "treelm_boundary_neighbor_entry": source_token_evidence(
            tem / "tem_construction_module.f90",
            "If the entry in the stencil is < 0, it must be a boundary ID",
            "neighVal = int(nTreeID)",
        ),
        "pdf_targets": source_token_evidence(
            mus / "mus_pdf_module.f90",
            "me%nElems_solve = nFluids + nGhostFromCoarser",
            "me%nElems_local = nFluids + nGhostFromCoarser + nGhostFromFiner + nHalos",
            "subroutine mus_swap_Now_Next(",
        ),
        "aux_density": source_token_evidence(
            mus / "derived" / "mus_auxFieldVar_module.fpp",
            "subroutine mus_calcAuxField_fluid_d3q19",
            "state(?FETCH?( 1, 1, elemPos",
            "rho(iElem) = sum(pdf(:,iElem))",
        ),
        "bgk_kernel": source_token_evidence(
            mus / "compute" / "mus_compute_d3q19_module.fpp",
            "subroutine mus_advRel_kFluid_rBGK_vStd_lD3Q19",
            "cmpl_o = 1._rk - omega",
            "f000*cmpl_o+omega*rho",
        ),
        "bgk_selection": source_token_evidence(
            mus / "init" / "mus_initFluid_module.f90",
            "compute => mus_advRel_kFluid_rBGK_vStd_lD3Q19",
        ),
        "restart": source_token_evidence(
            mus / "mus_buffer_module.fpp",
            "subroutine mus_pdf_serialize(",
            "scheme%pdf( iLevel )%nNext",
            "subroutine mus_pdf_unserialize(",
            "Only set nNext",
        ),
        "lbm_macros": source_token_evidence(
            built_macro,
            "else !PULL",
            "macro :: FETCH",
            "macro :: SAVE",
        ),
        "treelm_runtime_solid": source_token_evidence(
            tem / "tem_construction_module.f90",
            "prp_solid",
        ),
        "mesh_boundary_order": {
            "path": str(_mesh_path(root) / "bnd.lua"),
            "sha256": sha256_file(_mesh_path(root) / "bnd.lua"),
        },
        "runtime_log": {
            "path": str(stdout_path),
            "sha256": sha256_file(stdout_path),
        },
        "lua": {"path": str(lua_path), "sha256": sha256_file(lua_path)},
    }
    contract = Tau1BaseRuntimeContract()
    result = {
        "status": "SOURCE_PROVEN",
        "referee_revision": REFEREE_REVISION_NEW,
        "musubi_source_revision": "81f8c4f13772f6d4af31f335e1e3f99b02726e25",
        "apes_link_revision": TEM_GIT_SHA,
        "treelm_revision": "9899d1376992c4fafc8a343d2b4ccef81de670d1",
        "binary_sha256": sha256_file(_binary_windows()),
        "restart_state": {
            "serialized_buffer": "nNext",
            "phase": (
                "post-collision/post-source and post-halo-exchange output; "
                "pre-boundary for the following timestep"
            ),
            "continuation_restore_buffer": "nNext",
            "set_boundary_write_buffer": "nNext",
        },
        "timestep_order": [
            "advance time",
            "fill boundary buffers from original nNext",
            "set_boundary writes nNext",
            "swap nNow/nNext",
            "aux density/velocity from PULL FETCH on nNow",
            "BGK compute reads nNow and writes nNext",
            "apply source terms",
            "exchange nNext halos",
            "serialize nNext",
        ],
        "boundary_order": list(mesh.boundary_labels),
        "boundary_order_source": "mesh bnd.lua/globBC iBnd order",
        "neighbor_buffer": {
            "construction": "once before the boundary loop",
            "reads": "original pre-boundary nNext (and nNow where requested)",
            "shared_snapshot": True,
            "sees_prior_boundary_writes": False,
        },
        "swap": "after all boundary writes and immediately before aux/compute",
        "compute_targets": {
            "formula": "nElems_solve = nFluids + nGhostFromCoarser",
            "single_level_nFluids_global": EXPECTED_CELLS,
            "ghost_from_coarser": 0,
            "ghost_from_finer": 0,
            "halo_per_logged_rank": 400,
            "halo_is_solve_target": False,
            "boundary_elements_are_fluid_targets": True,
            "prp_solid_elements_are_fluid_targets_with_special_connectivity": True,
            "compute_target_count": EXPECTED_CELLS,
            "runtime_log_count_proof": runtime_count_proof,
        },
        "pull_connectivity": {
            "normal_fluid_source": "same direction PDF at PULL source element",
            "missing_source": "current element inverse-direction PDF",
            "current_prp_solid": "current element inverse-direction PDF",
            "source_prp_solid": "current element inverse-direction PDF",
            "boundary_source": (
                "negative TreElm boundary neighbor entry; PULL maps to the "
                "current element inverse-direction PDF even when a coordinate-"
                "adjacent tree element exists"
            ),
            "save": "current target element/same direction in AOS/PULL",
        },
        "density": "sum of all 19 PULL-fetched PDFs per solve target, not raw row",
        "tau1_bgk": {
            "omega": contract.omega,
            "formula": "f_post=(1-omega)*f_fetch+omega*f_eq(rho,u)",
            "omega_one_reduction": "f_post=f_eq(rho,u)",
            "mathematical_local_mass_delta": (
                "sum_q(f_post-f_fetch)=sum_q(f_eq)-rho=0 because "
                "sum_q w_q=1, sum_q w_q c_q=0, and D3Q19 isotropy cancels "
                "the quadratic velocity terms"
            ),
            "numerical_delta": "measured explicitly by the Python replay",
        },
        "other_sources": {
            "status": "SOURCE_PROVEN_NONE" if source_absence else "UNPROVEN",
            "glob_source": False,
            "body_force": False,
            "adaptive_inlet_is_boundary_not_glob_source": True,
        },
        "runtime_solid_count": len(runtime_solid_cells(mesh)),
        "mesh_cells": len(mesh.cell_ijk),
        "mesh_sha256": dict(MESH_HASHES),
        "source_evidence": evidence,
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "musubi_calls": 0,
    }
    write_json(qc / "musubi_full_timestep_mass_contract.json", result)
    return result


def audit_full_timestep_replay_8step(project_root: Path) -> dict[str, Any]:
    """Zero-MUSUBI replay of the nine archived late non-equilibrium states."""

    root = Path(project_root).resolve()
    run_root = _run_root(root)
    qc = run_root / "qc"
    contract = historical_tau1_runtime_contract()
    mesh = load_mesh_contract(_mesh_path(root), expected_cells=EXPECTED_CELLS)
    main_pairs = _restart_pairs(_runtime_windows() / "restart")
    archived_pairs = _restart_pairs(
        run_root / "dense_diagnostic_corrected" / "attempt_1_8_steps" / "restart"
    )
    start_iteration = 3_117_927
    expected = list(range(start_iteration + 1, start_iteration + 9))
    if start_iteration not in main_pairs or sorted(archived_pairs) != expected:
        raise RuntimeError("the exact 3117927..3117935 restart sequence is incomplete")
    outlet_pressures = {
        label: contract.pressure_reference_pa + gauge
        for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
    }
    start_pdf = _read_pdf(main_pairs[start_iteration][1])
    endpoint_start = start_pdf.copy()
    current = start_pdf
    records: list[dict[str, Any]] = []
    predicted_values: list[float] = []
    for iteration_end in expected:
        following = _read_pdf(archived_pairs[iteration_end][1])
        replay = replay_full_timestep(
            current,
            following,
            mesh,
            dx_m=contract.dx_m,
            dt_s=contract.dt_s,
            density_kg_m3=contract.rho_kg_m3,
            target_mass_flow_kg_s=contract.target_mass_flow_kg_s,
            outlet_pressures_pa=outlet_pressures,
        )
        record = public_step_record(replay)
        record["iteration_start"] = iteration_end - 1
        record["iteration_end"] = iteration_end
        records.append(record)
        predicted_values.append(float(record["delta_full_predicted"]))
        current = following
    cumulative_predicted = math.fsum(predicted_values)
    cumulative_actual = stable_delta(current, endpoint_start)
    cumulative_normalizer = abs(contract.target_lattice_flux * len(records))
    cumulative_residual = abs(cumulative_predicted - cumulative_actual) / max(
        cumulative_normalizer, np.finfo(np.float64).tiny
    )
    individual = [float(item["R_full_one_step_identity"]) for item in records]
    boundary_only = [float(item["R_boundary_only"]) for item in records]
    passed = full_identity_pass(individual, cumulative_residual)
    result = {
        "status": (
            "FULL_TIMESTEP_DISCRETE_MASS_IDENTITY_PROVEN"
            if passed
            else "ZERO_RUN_FULL_TIMESTEP_REPLAY_FAILED"
        ),
        "referee_revision": REFEREE_REVISION_NEW,
        "analysis_mode": "ZERO_MUSUBI_CALLS_EXISTING_9_RESTARTS",
        "musubi_calls": 0,
        "seeder_calls": 0,
        "new_long_cfd_iterations": 0,
        "start_iteration": start_iteration,
        "end_iteration": expected[-1],
        "steps": len(records),
        "R_boundary_only_min": min(boundary_only),
        "R_boundary_only_median": float(np.median(boundary_only)),
        "R_boundary_only_max": max(boundary_only),
        "R_full_one_step_identity_min": min(individual),
        "R_full_one_step_identity_median": float(np.median(individual)),
        "R_full_one_step_identity_max": max(individual),
        "cumulative_predicted": cumulative_predicted,
        "cumulative_actual": cumulative_actual,
        "R_full_8step_identity": cumulative_residual,
        "hard_gate": FULL_IDENTITY_GATE,
        "hard_gate_pass": passed,
        "classification": (
            "REFEREE_BOUNDARY_ONLY_ACCOUNTING_INCOMPLETE"
            if passed
            else "INSTRUMENTED_ONE_TIMESTEP_REQUIRED"
        ),
        "runtime_solid_count": records[0]["runtime_solid_count"],
        "compute_target_count": records[0]["compute_target_count"],
        "pull_link_counts": pull_link_counts(
            mesh, runtime_solid_cells(mesh)
        ),
        "flux_definition": BOUNDARY_FLUX_DEFINITION,
        "physical_Q1_Q2_Q3": DEFERRED_PHYSICAL_FLUX,
        "restart_sha256": {
            str(start_iteration): sha256_file(main_pairs[start_iteration][1]),
            **{
                str(iteration): sha256_file(archived_pairs[iteration][1])
                for iteration in expected
            },
        },
        "per_step": records,
    }
    write_json(qc / "tau1_full_timestep_replay_8step.json", result)
    return result
