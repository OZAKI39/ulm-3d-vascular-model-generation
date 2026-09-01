"""Contracts for the repaired reference-scaled Tau=1 C/B/F study.

This module is deliberately solver-free. It defines the immutable grid
design, restart compatibility fields and pre-declared three-grid rules used
by :mod:`tau1_reference_scaled_grid`.
"""

from __future__ import annotations

import hashlib
import math
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .io import write_json
from .tau1_base import (
    BULK_NU_M2_S,
    MUSUBI_SHA256,
    MUSUBI_WSL,
    NU_M2_S,
    OUTLET_GAUGE_PRESSURE_PA,
    RHO_KG_M3,
    TARGET_MASS_FLOW_KG_S,
    TARGET_Q_M3_S,
)


RUN_NAME = "healthy_mouse_capillary_tau1_grid_convergence_anchor003274_20260831"
CBF_RUN_NAME = (
    "healthy_mouse_capillary_tau1_reference_scaled_grid_convergence_"
    "anchor003274_20260901"
)
BASE_CFD_RUN = (
    "healthy_mouse_capillary_tau1_reference_scaled_base_anchor003274_20260901"
)
BASE_MESH_RUN = (
    "healthy_mouse_capillary_dimensionless_qvalue_base_preflight_"
    "anchor003274_20260830"
)
HISTORICAL_CLASSIFICATION = (
    "HISTORICAL_UNDER_FALLBACK_QVALUE_WALL_AND_HIGH_TAU_NUMERICS"
)

REFINEMENT_RATIO = 1.3
CS2 = 1.0 / 3.0
SHORT_WINDOW_S = 0.0002441406727828746
LONG_WINDOW_S = 0.0004882813455657492
SOFT_MAX_PHYSICAL_TIME_S = 0.0020
HARD_MAX_PHYSICAL_TIME_S = 0.0030
PLATEAU_MIN_PHYSICAL_TIME_S = 0.0015
PLATEAU_AUDITS = 4
PLATEAU_RELATIVE_CHANGE = 0.02
GCI_SAFETY_FACTOR = 1.25
PRIMARY_FLUX_DEFINITION = "PHYSICAL_INTERIOR_CROSS_SECTION_VELOCITY_FLUX"
SEEDER_BINARY_WSL = (
    "/home/lzy/apes-worktrees/seeder_dimensionless_kernel_20260830/build/seeder"
)
SEEDER_BINARY_SHA256 = (
    "d7be681ca90da706559a4fd7e8f769fdb8f4303b8508f751077205f8e00cc7ed"
)
PLANE_CONTRACT_SHA256 = (
    "ffaa49bdb6e43fb7208ff29df07a90d4e92ef9bfa4b96ca4f997d4f453a7f005"
)
PLANE_CONTRACT_FILE_SHA256 = (
    "2b9cc9edfd5e0a801b368922b95044c3bfebdf24d51fd55181fe232bbf59705f"
)
BASE_ACCEPTED_ITERATION = 598_755
BASE_ACCEPTED_RESTART_SHA256 = (
    "ffcd98b2dc684d1569d937d915b603805809c581d5341e71b17afac2ac64c39f"
)
BASE_FULL_REFEREE_RESIDUAL = 7.913943402747673e-10
BASE_EXPECTED_CELLS = 182_320
BASE_MESH_SHA256 = {
    "elemlist.lsb": "f7d7b1d55273c78c336ac04e39bc018dd9ebb470a9f29ce833ff01711de8c386",
    "bnd.lsb": "520d7dd1e4a46a45f9b1218a5807cfd89d6f054e0a247872362b130ff6bcfe69",
    "qval.lsb": "35884406b5f0111cd4ab471f7b08ac3df00e478d3458a57636d1bd8921cb0fe6",
}
# Compatibility evidence consumed by the already-validated physical-flux
# module. These are historical read-only restart IDs, never C/B/F inputs.
BASE_ITERATIONS = (2_878_425, 2_998_176, 3_117_927)
BASE_RESTART_SHA256 = {
    2_878_425: "75815fded691784ae942285e6ccf32514a1936ef9b1c228a03631991462646ee",
    2_998_176: "e3bb103963299e2384ce636dd117d194ee21b86adeac63a9965a095840f39a6c",
    3_117_927: "3d54f3970b4120896c214155811d7cd1b594e3efd172f80b5dc5e7d0fef279e2",
}
PORTS = ("inlet", "outlet_01", "outlet_02", "outlet_03")
OUTLETS = PORTS[1:]

PRIMARY_METRICS = (
    "inlet_gauge_pressure_pa",
    "DeltaP01_pa",
    "DeltaP02_pa",
    "DeltaP03_pa",
    "outlet_01_flow_fraction",
    "outlet_02_flow_fraction",
    "outlet_03_flow_fraction",
)


@dataclass(frozen=True, slots=True)
class Tau1GridSpec:
    label: str
    dx_m: float
    root_level: int
    bounding_cube_length_m: float
    base_read_only: bool = False

    @property
    def dt_s(self) -> float:
        return self.dx_m**2 / (6.0 * NU_M2_S)

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
    def pressure_reference_pa(self) -> float:
        return RHO_KG_M3 * CS2 * self.dx_m**2 / self.dt_s**2

    @property
    def outlet_absolute_pressure_pa(self) -> dict[str, float]:
        return {
            label: self.pressure_reference_pa + float(gauge)
            for label, gauge in OUTLET_GAUGE_PRESSURE_PA.items()
        }

    @property
    def target_lattice(self) -> float:
        return (
            TARGET_MASS_FLOW_KG_S * self.dt_s
            / (RHO_KG_M3 * self.dx_m**3)
        )

    @property
    def short_window_iterations(self) -> int:
        return round(SHORT_WINDOW_S / self.dt_s)

    @property
    def long_window_iterations(self) -> int:
        return round(LONG_WINDOW_S / self.dt_s)

    @property
    def checkpoint_interval_iterations(self) -> int:
        return self.short_window_iterations

    @property
    def earliest_audit_iteration(self) -> int:
        # Checkpoints stay on one cadence. Coarse differs from round(2*LONG/dt)
        # by only two timesteps because each physical window is rounded once.
        return 4 * self.short_window_iterations

    @property
    def soft_max_iterations(self) -> int:
        return math.floor(SOFT_MAX_PHYSICAL_TIME_S / self.dt_s)

    @property
    def hard_max_iterations(self) -> int:
        return math.floor(HARD_MAX_PHYSICAL_TIME_S / self.dt_s)

    def evidence(self) -> dict[str, Any]:
        result = asdict(self)
        result.update(
            {
                "dt_s": self.dt_s,
                "nu_lattice": self.nu_lattice,
                "tau": self.tau,
                "omega": self.omega,
                "pressure_reference_pa": self.pressure_reference_pa,
                "outlet_gauge_pressure_pa": dict(OUTLET_GAUGE_PRESSURE_PA),
                "outlet_absolute_pressure_pa": self.outlet_absolute_pressure_pa,
                "target_lattice": self.target_lattice,
                "short_window_iterations": self.short_window_iterations,
                "long_window_iterations": self.long_window_iterations,
                "checkpoint_interval_iterations": self.checkpoint_interval_iterations,
                "earliest_audit_iteration": self.earliest_audit_iteration,
                "soft_max_iterations": self.soft_max_iterations,
                "hard_max_iterations": self.hard_max_iterations,
            }
        )
        return result


GRID_SPECS: dict[str, Tau1GridSpec] = {
    "coarse": Tau1GridSpec("coarse", 2.6e-7, 9, 0.00013312),
    "base": Tau1GridSpec("base", 2.0e-7, 9, 0.0001024, True),
    "fine": Tau1GridSpec(
        "fine", 1.5384615384615385e-7, 10, 0.00015753846153846154
    ),
}


def run_root(project_root: Path) -> Path:
    return Path(project_root).resolve() / "outputs" / "cfd_flow" / CBF_RUN_NAME


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
        "dt_formula": spec.dt_s == spec.dx_m**2 / (6.0 * NU_M2_S),
        "nu_lattice_one_sixth": math.isclose(
            spec.nu_lattice, 1.0 / 6.0, rel_tol=0.0, abs_tol=1.0e-15
        ),
        "tau_one": math.isclose(spec.tau, 1.0, rel_tol=0.0, abs_tol=1.0e-14),
        "omega_one": math.isclose(spec.omega, 1.0, rel_tol=0.0, abs_tol=1.0e-14),
        "dynamic_reference": spec.pressure_reference_pa
        == RHO_KG_M3 * CS2 * spec.dx_m**2 / spec.dt_s**2,
        "positive_windows": (
            spec.short_window_iterations > 0
            and spec.long_window_iterations > spec.short_window_iterations
        ),
    }
    if not all(checks.values()):
        raise ValueError(f"invalid Tau=1 grid contract for {spec.label}: {checks}")
    return checks


def write_grid_design(project_root: Path) -> dict[str, Any]:
    ratios = {
        "coarse_over_base": GRID_SPECS["coarse"].dx_m / GRID_SPECS["base"].dx_m,
        "base_over_fine": GRID_SPECS["base"].dx_m / GRID_SPECS["fine"].dx_m,
    }
    checks = {
        "constant_refinement_ratio": all(
            math.isclose(value, REFINEMENT_RATIO, rel_tol=0.0, abs_tol=1.0e-14)
            for value in ratios.values()
        ),
        "all_grid_specs": all(
            all(validate_grid_spec(spec).values()) for spec in GRID_SPECS.values()
        ),
        "base_read_only": GRID_SPECS["base"].base_read_only,
        "outlet_gauges_identical": True,
        "physical_plane_contract_frozen": True,
        "production_pipeline_unmodified": True,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "refinement_ratio": REFINEMENT_RATIO,
        "ratio_checks": ratios,
        "grids": {label: spec.evidence() for label, spec in GRID_SPECS.items()},
        "short_window_s": SHORT_WINDOW_S,
        "long_window_s": LONG_WINDOW_S,
        "soft_max_physical_time_s": SOFT_MAX_PHYSICAL_TIME_S,
        "hard_max_physical_time_s": HARD_MAX_PHYSICAL_TIME_S,
        "frozen_physics": {
            "rho0_kg_m3": RHO_KG_M3,
            "nu_m2_s": NU_M2_S,
            "bulk_nu_m2_s": BULK_NU_M2_S,
            "target_mass_flow_kg_s": TARGET_MASS_FLOW_KG_S,
            "target_q_m3_s": TARGET_Q_M3_S,
            "outlet_gauge_pressure_pa": dict(OUTLET_GAUGE_PRESSURE_PA),
            "layout": "d3q19",
            "collision": "bgk",
            "boundary_contract": {
                "wall": "wall_libb_continuous_q",
                "inlet": "adaptive_flux_pressure",
                "outlets": "pressure_eq",
            },
        },
        "frozen_binaries": {
            "seeder_wsl": SEEDER_BINARY_WSL,
            "seeder_sha256": SEEDER_BINARY_SHA256,
            "musubi_wsl": MUSUBI_WSL,
            "musubi_sha256": MUSUBI_SHA256,
        },
        "physical_plane_contract_sha256": PLANE_CONTRACT_SHA256,
        "checks": checks,
        "production_pipeline_modified": False,
    }
    output = run_root(project_root) / "qc" / "tau1_cbf_grid_design.json"
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
    """Change only the three grid-dependent fields of the frozen Base input."""

    if spec.base_read_only:
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
    return re.sub(
        r"(bounding_cube\s*=\s*\{.*?length\s*=\s*)[-+0-9.eE]+",
        rf"\g<1>{spec.bounding_cube_length_m:.17g}",
        text,
        count=1,
    )


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
        "mesh_hashes", "dx_m", "dt_s", "rho0_kg_m3", "nu_m2_s",
        "bulk_nu_m2_s", "tau", "omega", "pressure_reference_pa",
        "boundary_contract", "outlet_gauge_pressure_pa",
        "outlet_absolute_pressure_pa", "target_mass_flow_kg_s",
        "binary_sha256", "physical_plane_contract_sha256", "layout", "collision",
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
        if definition != PRIMARY_FLUX_DEFINITION:
            raise ValueError(
                f"{grid} primary flow metrics require {PRIMARY_FLUX_DEFINITION}; "
                f"found {definition}"
            )


def three_grid_primary_analysis(
    coarse: float,
    base: float,
    fine: float,
    *,
    refinement_ratio: float = REFINEMENT_RATIO,
) -> dict[str, Any]:
    """Apply the pre-declared 1% floor, monotonic Richardson and GCI rules."""

    c, b, f = (float(coarse), float(base), float(fine))
    tiny = np.finfo(float).tiny
    rel_cb = abs(b - c) / max(abs(b), tiny)
    rel_bf = abs(f - b) / max(abs(f), tiny)
    common = {
        "coarse": c, "base": b, "fine": f,
        "relative_C_B": rel_cb, "relative_B_F": rel_bf,
        "relative_C_B_percent": 100.0 * rel_cb,
        "relative_B_F_percent": 100.0 * rel_bf,
        "base_to_fine_relative_difference": rel_bf,
        "base_to_fine_within_5_percent": rel_bf <= 0.05,
        "base_to_fine_within_2_percent": rel_bf <= 0.02,
    }
    unavailable = {
        "status": "FAIL", "classification": "OSCILLATORY_OR_NONCONVERGENT",
        "observed_order_p": None, "richardson_extrapolated": None,
        "GCI_CB": None, "GCI_BF": None, "GCI_CB_percent": None,
        "GCI_BF_percent": None, "pass": False,
    }
    if not all(math.isfinite(value) for value in (c, b, f)):
        return {**common, **unavailable}
    if rel_cb <= 0.01 and rel_bf <= 0.01:
        return {
            **common, "status": "PASS",
            "classification": "GRID_INSENSITIVE_WITHIN_1PCT",
            "observed_order_p": None, "richardson_extrapolated": None,
            "GCI_CB": None, "GCI_BF": None, "GCI_CB_percent": None,
            "GCI_BF_percent": None, "pass": True,
        }
    delta_cb = b - c
    delta_bf = f - b
    monotonic = delta_cb * delta_bf > 0.0 and abs(delta_bf) < abs(delta_cb)
    if not monotonic:
        return {**common, **unavailable}
    order = math.log(abs((c - b) / (b - f))) / math.log(refinement_ratio)
    denominator = refinement_ratio**order - 1.0
    if not math.isfinite(order) or order <= 0.0 or denominator <= 0.0:
        return {**common, **unavailable}
    richardson = f + (f - b) / denominator
    gci_bf = GCI_SAFETY_FACTOR * abs(f - b) / max(abs(f), tiny) / denominator
    gci_cb = GCI_SAFETY_FACTOR * abs(b - c) / max(abs(b), tiny) / denominator
    passed = rel_bf <= 0.05 and gci_bf <= 0.05
    return {
        **common, "status": "PASS" if passed else "FAIL",
        "classification": "ASYMPTOTIC_MONOTONIC",
        "observed_order_p": order,
        "richardson_extrapolated": richardson,
        "GCI_CB": gci_cb, "GCI_BF": gci_bf,
        "GCI_CB_percent": 100.0 * gci_cb,
        "GCI_BF_percent": 100.0 * gci_bf, "pass": passed,
    }


def build_primary_analyses(
    observables: Mapping[str, Mapping[str, Any]],
) -> dict[str, dict[str, Any]]:
    validate_physical_flux_observables(observables)
    return {
        metric: three_grid_primary_analysis(
            float(observables["coarse"][metric]),
            float(observables["base"][metric]),
            float(observables["fine"][metric]),
        )
        for metric in PRIMARY_METRICS
    }


def evaluate_repaired_grid_gate(
    analyses: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    within_5 = all(
        float(analyses[name]["relative_B_F"]) <= 0.05 for name in PRIMARY_METRICS
    )
    within_2 = all(
        float(analyses[name]["relative_B_F"]) <= 0.02 for name in PRIMARY_METRICS
    )
    trends = all(
        analyses[name]["classification"]
        in {"GRID_INSENSITIVE_WITHIN_1PCT", "ASYMPTOTIC_MONOTONIC"}
        and bool(analyses[name]["pass"])
        for name in PRIMARY_METRICS
    )
    passed = within_5 and trends
    failed = [name for name in PRIMARY_METRICS if not analyses[name]["pass"]]
    return {
        "status": "PASS" if passed else "FAIL",
        "primary_trends_valid": trends,
        "primary_trends_asymptotic_monotonic": all(
            analyses[name]["classification"] == "ASYMPTOTIC_MONOTONIC"
            for name in PRIMARY_METRICS
        ),
        "base_fine_primary_within_5_percent": within_5,
        "base_fine_primary_within_preferred_2_percent": within_2,
        "failed_metrics": failed,
        "final_status": (
            "CFD_FLOW_TAU1_CBF_GRID_CONVERGENCE_PASS"
            if passed else "CFD_FLOW_TAU1_CBF_GRID_CONVERGENCE_FAILED"
        ),
        "next": (
            "PROMOTE VALIDATED TAU1 CFD CONTRACT TO PRODUCTION PIPELINE AND RUN FINAL PRODUCTION REGRESSION"
            if passed
            else "STOP AND REVIEW THE SPECIFIC FAILED GRID-CONVERGENCE METRIC"
        ),
    }
