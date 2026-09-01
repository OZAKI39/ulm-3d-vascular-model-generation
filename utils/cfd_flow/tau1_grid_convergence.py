"""Research-only repaired Tau=1 three-grid convergence contracts.

Physical port-flux extraction lives in :mod:`physical_port_flux`; this module
contains only grid-design and convergence-gate contracts.
"""

from __future__ import annotations

import hashlib
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

from .full_timestep_mass_referee import PHYSICAL_FLUX_DEFINITION
from .grid_convergence import (
    evaluate_grid_convergence_gate,
    three_grid_scalar_analysis,
)
from .io import write_json
from .tau1_base import PRESSURE_REFERENCE_PA


RUN_NAME = "healthy_mouse_capillary_tau1_grid_convergence_anchor003274_20260831"
BASE_MESH_RUN = (
    "healthy_mouse_capillary_dimensionless_qvalue_base_preflight_anchor003274_20260830"
)
BASE_CFD_RUN = "healthy_mouse_capillary_tau1_base_anchor003274_20260830"
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
OUTLET_GAUGE_PRESSURE_PA = {
    "outlet_01": 14.544978101,
    "outlet_02": 132.204549223,
    "outlet_03": -13.700626673,
}
SHORT_WINDOW_S = 0.0002441406727828746
LONG_WINDOW_S = 0.0004882813455657492
HARD_MAX_PHYSICAL_TIME_S = 0.0065

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
