"""Solver-free optional three-grid analysis for the validated Tau=1 contract.

Fine steady evidence was not completed, so this utility cannot report a
formal current C/B/F result. It is excluded from the default integration
context and may only analyze future complete evidence.
"""

from __future__ import annotations

import hashlib
import math
import re
from dataclasses import dataclass
from typing import Any, Mapping

import numpy as np

from .validated_contract import (
    OUTLET_GAUGE_PRESSURES_PA,
    TARGET_MASS_FLOW_KG_S,
    pressure_reference_pa,
    relaxation_from_physical,
    restart_compatibility as validated_restart_compatibility,
    target_lattice_flux,
    tau1_time_step_s,
)


REFINEMENT_RATIO = 1.3
GCI_SAFETY_FACTOR = 1.25
PRIMARY_FLUX_DEFINITION = "PHYSICAL_INTERIOR_CROSS_SECTION_VELOCITY_FLUX"
HISTORICAL_CLASSIFICATION = "HISTORICAL_UNDER_FALLBACK_QVALUE_WALL_AND_HIGH_TAU_NUMERICS"
PRIMARY_METRICS = (
    "inlet_gauge_pressure_pa", "DeltaP01_pa", "DeltaP02_pa", "DeltaP03_pa",
    "outlet_01_flow_fraction", "outlet_02_flow_fraction", "outlet_03_flow_fraction",
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
        return tau1_time_step_s(self.dx_m)

    @property
    def nu_lattice(self) -> float:
        return relaxation_from_physical(self.dx_m, self.dt_s)[0]

    @property
    def tau(self) -> float:
        return relaxation_from_physical(self.dx_m, self.dt_s)[1]

    @property
    def omega(self) -> float:
        return relaxation_from_physical(self.dx_m, self.dt_s)[2]

    @property
    def pressure_reference_pa(self) -> float:
        return pressure_reference_pa(self.dx_m, self.dt_s)

    @property
    def outlet_absolute_pressure_pa(self) -> dict[str, float]:
        return {label: self.pressure_reference_pa + gauge for label, gauge in OUTLET_GAUGE_PRESSURES_PA.items()}

    @property
    def target_lattice(self) -> float:
        return target_lattice_flux(TARGET_MASS_FLOW_KG_S, self.dx_m, self.dt_s)


GRID_SPECS: dict[str, Tau1GridSpec] = {
    "coarse": Tau1GridSpec("coarse", 2.6e-7, 9, 0.00013312),
    "base": Tau1GridSpec("base", 2.0e-7, 9, 0.0001024, True),
    "fine": Tau1GridSpec("fine", 1.5384615384615385e-7, 10, 0.00015753846153846154),
}


def validate_grid_spec(spec: Tau1GridSpec) -> dict[str, bool]:
    checks = {
        "dt_formula": spec.dt_s == tau1_time_step_s(spec.dx_m),
        "nu_lattice_one_sixth": math.isclose(spec.nu_lattice, 1.0 / 6.0, rel_tol=0.0, abs_tol=1.0e-15),
        "tau_one": math.isclose(spec.tau, 1.0, rel_tol=0.0, abs_tol=1.0e-14),
        "omega_one": math.isclose(spec.omega, 1.0, rel_tol=0.0, abs_tol=1.0e-14),
        "dynamic_reference": spec.pressure_reference_pa == pressure_reference_pa(spec.dx_m, spec.dt_s),
    }
    if not all(checks.values()):
        raise ValueError(f"invalid Tau=1 grid contract for {spec.label}: {checks}")
    return checks


def assert_launch_allowed(grid: str, operation: str) -> None:
    if grid not in GRID_SPECS:
        raise ValueError(grid)
    if operation not in {"seeder", "long_musubi", "referee_one_step"}:
        raise ValueError(operation)
    if grid == "base" and operation in {"seeder", "long_musubi"}:
        raise PermissionError(f"accepted Base is read-only: {operation} is forbidden")


def render_repaired_seeder_config(base_text: str, spec: Tau1GridSpec) -> str:
    """Change only grid metadata while preserving frozen physical objects."""

    if spec.base_read_only:
        raise PermissionError("accepted Base Seeder config must not be regenerated")
    required = (
        "7.8082771579034594e-05", "5.4076359424537688e-05",
        "8.8206398940868488e-05", "0.0001016827715790346", "calc_dist = true",
    )
    if not all(token in base_text for token in required):
        raise ValueError("frozen Base physical Seeder objects are incomplete")
    text = re.sub(
        r"(?m)^comment\s*=.*$",
        f"comment = 'repaired Tau1 {spec.label} grid; physical objects frozen from Base'",
        base_text,
        count=1,
    )
    text = re.sub(r"(?m)^minlevel\s*=\s*\d+\s*$", f"minlevel = {spec.root_level}", text, count=1)
    return re.sub(
        r"(bounding_cube\s*=\s*\{.*?length\s*=\s*)[-+0-9.eE]+",
        rf"\g<1>{spec.bounding_cube_length_m:.17g}", text, count=1,
    )


def seeder_physical_spatial_signature(text: str) -> str:
    normalized = re.sub(r"(?m)^comment\s*=.*$", "comment=<GRID>", text)
    normalized = re.sub(r"(?m)^minlevel\s*=\s*\d+\s*$", "minlevel=<GRID>", normalized)
    normalized = re.sub(
        r"(bounding_cube\s*=\s*\{.*?length\s*=\s*)[-+0-9.eE]+",
        r"\g<1><GRID>", normalized, count=1,
    )
    return hashlib.sha256(normalized.encode()).hexdigest()


def restart_resume_contract(saved: Mapping[str, Any], expected: Mapping[str, Any]) -> dict[str, Any]:
    legacy_fields = (
        "mesh_hashes", "dx_m", "dt_s", "rho_kg_m3", "nu_m2_s", "bulk_nu_m2_s",
        "tau", "boundary_contract", "outlet_pressures_pa", "target_mass_flow_kg_s",
        "binary_sha256", "layout", "relaxation",
    )
    if any(field in saved or field in expected for field in ("rho_kg_m3", "outlet_pressures_pa", "relaxation")):
        checks = {field: saved.get(field) == expected.get(field) for field in legacy_fields}
        return {"status": "PASS" if all(checks.values()) else "FAIL", "checks": checks}
    return validated_restart_compatibility(saved, expected)


def validate_physical_flux_observables(observables: Mapping[str, Mapping[str, Any]]) -> None:
    for grid in ("coarse", "base", "fine"):
        if observables[grid].get("historical_classification") == HISTORICAL_CLASSIFICATION:
            raise ValueError(f"{grid} historical fallback-q evidence is excluded")
        definition = observables[grid].get("flux_definition")
        if definition != PRIMARY_FLUX_DEFINITION:
            raise ValueError(f"{grid} primary flow metrics require {PRIMARY_FLUX_DEFINITION}; found {definition}")


def three_grid_primary_analysis(
    coarse: float, base: float, fine: float, *, refinement_ratio: float = REFINEMENT_RATIO
) -> dict[str, Any]:
    """Apply the pre-declared 1% floor and monotonic Richardson/GCI rules."""

    c, b, f = float(coarse), float(base), float(fine)
    tiny = np.finfo(float).tiny
    rel_cb = abs(b - c) / max(abs(b), tiny)
    rel_bf = abs(f - b) / max(abs(f), tiny)
    common = {
        "coarse": c, "base": b, "fine": f, "relative_C_B": rel_cb, "relative_B_F": rel_bf,
        "relative_C_B_percent": 100.0 * rel_cb, "relative_B_F_percent": 100.0 * rel_bf,
        "base_to_fine_relative_difference": rel_bf,
        "base_to_fine_within_5_percent": rel_bf <= 0.05,
        "base_to_fine_within_2_percent": rel_bf <= 0.02,
    }
    unavailable = {
        "status": "FAIL", "classification": "OSCILLATORY_OR_NONCONVERGENT",
        "observed_order_p": None, "richardson_extrapolated": None,
        "GCI_CB": None, "GCI_BF": None, "GCI_CB_percent": None, "GCI_BF_percent": None,
        "pass": False,
    }
    if not all(math.isfinite(value) for value in (c, b, f)):
        return {**common, **unavailable}
    if rel_cb <= 0.01 and rel_bf <= 0.01:
        return {
            **common, "status": "PASS", "classification": "GRID_INSENSITIVE_WITHIN_1PCT",
            "observed_order_p": None, "richardson_extrapolated": None,
            "GCI_CB": None, "GCI_BF": None, "GCI_CB_percent": None, "GCI_BF_percent": None,
            "pass": True,
        }
    delta_cb, delta_bf = b - c, f - b
    if not (delta_cb * delta_bf > 0.0 and abs(delta_bf) < abs(delta_cb)):
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
        **common, "status": "PASS" if passed else "FAIL", "classification": "ASYMPTOTIC_MONOTONIC",
        "observed_order_p": order, "richardson_extrapolated": richardson,
        "GCI_CB": gci_cb, "GCI_BF": gci_bf,
        "GCI_CB_percent": 100.0 * gci_cb, "GCI_BF_percent": 100.0 * gci_bf,
        "pass": passed,
    }


def build_primary_analyses(observables: Mapping[str, Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    validate_physical_flux_observables(observables)
    return {
        metric: three_grid_primary_analysis(
            float(observables["coarse"][metric]),
            float(observables["base"][metric]),
            float(observables["fine"][metric]),
        )
        for metric in PRIMARY_METRICS
    }


def evaluate_repaired_grid_gate(analyses: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    within_5 = all(float(analyses[name]["relative_B_F"]) <= 0.05 for name in PRIMARY_METRICS)
    within_2 = all(float(analyses[name]["relative_B_F"]) <= 0.02 for name in PRIMARY_METRICS)
    trends = all(
        analyses[name]["classification"] in {"GRID_INSENSITIVE_WITHIN_1PCT", "ASYMPTOTIC_MONOTONIC"}
        and bool(analyses[name]["pass"])
        for name in PRIMARY_METRICS
    )
    passed = within_5 and trends
    return {
        "status": "PASS" if passed else "FAIL", "primary_trends_valid": trends,
        "primary_trends_asymptotic_monotonic": all(
            analyses[name]["classification"] == "ASYMPTOTIC_MONOTONIC" for name in PRIMARY_METRICS
        ),
        "base_fine_primary_within_5_percent": within_5,
        "base_fine_primary_within_preferred_2_percent": within_2,
        "failed_metrics": [name for name in PRIMARY_METRICS if not analyses[name]["pass"]],
    }
