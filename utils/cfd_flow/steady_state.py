"""Pure, generic acceptance audit for validated physical-time CFD samples."""

from __future__ import annotations

import math
from collections.abc import Callable, Mapping, Sequence
from typing import Any

import numpy as np

from .validated_contract import (
    CONTROLLER_GATE,
    FLOW_FRACTION_DRIFT_GATE,
    INLET_GATE,
    MASS_GATE,
    MAXIMUM_LATTICE_SPEED,
    OUTLETS,
    PRESSURE_GATE,
    Q_DENSITY_GATE,
    SIGNIFICANT_BACKFLOW_FRACTION,
    TARGET_VOLUME_FLOW_M3_S,
    VELOCITY_GATE,
    ValidatedTau1Contract,
)


def trapezoidal_mean(
    samples: Sequence[Mapping[str, Any]],
    getter: Callable[[Mapping[str, Any]], float],
) -> float:
    if len(samples) < 2:
        raise ValueError("a physical-time mean requires at least two samples")
    values = [float(getter(item)) for item in samples]
    iterations = [int(item["iteration"]) for item in samples]
    duration = iterations[-1] - iterations[0]
    if duration <= 0 or any(right <= left for left, right in zip(iterations, iterations[1:])):
        raise ValueError("sample iterations must be strictly increasing")
    integral = math.fsum(
        0.5 * (values[index] + values[index + 1])
        * (iterations[index + 1] - iterations[index])
        for index in range(len(values) - 1)
    )
    return integral / duration


def audit_steady_window(
    samples: Sequence[Mapping[str, Any]],
    *,
    all_checkpoint_rho_pass: bool,
    expected_short_window_iterations: int | None = None,
    contract: ValidatedTau1Contract | None = None,
) -> dict[str, Any]:
    """Audit long-start/short-start/end samples without any solver I/O."""

    if len(samples) != 3:
        raise ValueError("steady audit requires long-start/short-start/end samples")
    contract = contract or ValidatedTau1Contract()
    a, b, c = samples
    iterations = [int(item["iteration"]) for item in samples]
    if iterations[1] - iterations[0] != iterations[2] - iterations[1]:
        raise ValueError(f"unequal physical audit windows: {iterations}")
    if expected_short_window_iterations is not None and iterations[1] - iterations[0] != expected_short_window_iterations:
        raise ValueError(f"unexpected physical audit window: {iterations}")

    def q(sample: Mapping[str, Any], label: str) -> float:
        return float(sample["ports"][label]["Q_velocity_m3_s"])

    def flow_metrics(window: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
        inlet = trapezoidal_mean(window, lambda item: q(item, "inlet"))
        outlets = {
            label: trapezoidal_mean(window, lambda item, label=label: q(item, label))
            for label in OUTLETS
        }
        outlet_sum = math.fsum(outlets.values())
        scale = max(abs(inlet), np.finfo(float).tiny)
        backflow = {
            label: value < 0.0 and abs(value) > SIGNIFICANT_BACKFLOW_FRACTION * abs(inlet)
            for label, value in outlets.items()
        }
        return {
            "mean_Qin_m3_s": inlet,
            "mean_outlets_m3_s": outlets,
            "mean_Qout_sum_m3_s": outlet_sum,
            "R_mass": abs(inlet - outlet_sum) / scale,
            "significant_backflow_by_outlet": backflow,
        }

    short = flow_metrics((b, c))
    long = flow_metrics((a, b, c))
    r_velocity = abs(float(c["mean_speed_m_s"]) - float(b["mean_speed_m_s"])) / max(
        abs(float(c["mean_speed_m_s"])), 1.0e-12
    )
    pressure_names = ("inlet_gauge_pressure_pa", *OUTLETS)

    def pressure(sample: Mapping[str, Any], name: str) -> float:
        return float(sample[name] if name == "inlet_gauge_pressure_pa" else sample["pressure_drops_pa"][name])

    pressure_residuals = {
        name: abs(pressure(c, name) - pressure(b, name)) / max(abs(pressure(c, name)), 1.0)
        for name in pressure_names
    }
    r_pressure = max(pressure_residuals.values())
    r_inlet = abs(short["mean_Qin_m3_s"] - TARGET_VOLUME_FLOW_M3_S) / TARGET_VOLUME_FLOW_M3_S
    fraction_drift = {
        label: max(float(item["flow_fractions"][label]) for item in samples)
        - min(float(item["flow_fractions"][label]) for item in samples)
        for label in OUTLETS
    }
    controller = c["controller"]
    target_error = abs(float(controller["target_lattice"]) - contract.target_lattice_flux) / contract.target_lattice_flux
    q_density_pass = all(
        float(value["residual"]) <= Q_DENSITY_GATE and bool(value.get("pass", True))
        for value in c["Q_density_consistency"].values()
    )
    significant_backflow = any(
        (*short["significant_backflow_by_outlet"].values(), *long["significant_backflow_by_outlet"].values())
    )
    gates = {
        "R_mass_short": short["R_mass"] <= MASS_GATE,
        "R_mass_long": long["R_mass"] <= MASS_GATE,
        "physical_volume_closure": float(c["physical_volume_closure"]) <= MASS_GATE,
        "R_velocity": r_velocity <= VELOCITY_GATE,
        "R_pressure": r_pressure <= PRESSURE_GATE,
        "R_inlet": r_inlet <= INLET_GATE,
        "flow_fraction_drift": max(fraction_drift.values()) <= FLOW_FRACTION_DRIFT_GATE,
        "Q_density_consistency": q_density_pass,
        "rho_sanity_all_checkpoints": bool(all_checkpoint_rho_pass),
        "no_significant_averaged_backflow": not significant_backflow,
        "minimum_pdf_positive": float(c["minimum_pdf"]) > 0.0,
        "maximum_lattice_speed": float(c["maximum_lattice_speed"]) < MAXIMUM_LATTICE_SPEED,
        "all_finite": bool(c["all_finite"]),
        "controller_target": target_error <= CONTROLLER_GATE,
        "controller_controlled_flux": float(controller["relative_error"]) <= CONTROLLER_GATE,
    }
    failed = [name for name, passed in gates.items() if not passed]
    return {
        "status": "PASS_NON_REFEREE" if not failed else "FAIL",
        "iteration": iterations[-1],
        "physical_time_s": iterations[-1] * contract.dt_s,
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
        "maximum_flow_fraction_drift": max(fraction_drift.values()),
        "significant_averaged_backflow": significant_backflow,
        "controller_target_expected": contract.target_lattice_flux,
        "controller_target_observed": float(controller["target_lattice"]),
        "controller_target_error": target_error,
        "gates": gates,
        "failed_gates": failed,
    }


def acceptance_transition(
    *, candidate_iteration: int | None, current_audit_pass: bool, iteration: int
) -> tuple[str, int | None]:
    if not current_audit_pass:
        return "CONTINUE", None
    if candidate_iteration is None:
        return "CANDIDATE", int(iteration)
    return "CONFIRMED", int(candidate_iteration)
