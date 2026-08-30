"""Contracts for the isolated Musubi periodic wall-geometry benchmark."""

from __future__ import annotations

from typing import Any, Iterable, Sequence

import numpy as np


def poiseuille_pressure_gradient_pa_m(
    *, mean_velocity_m_s: float, radius_m: float, rho_kg_m3: float, nu_m2_s: float
) -> float:
    """Pressure-gradient magnitude giving the requested circular-pipe mean."""

    if min(mean_velocity_m_s, radius_m, rho_kg_m3, nu_m2_s) <= 0.0:
        raise ValueError("Poiseuille inputs must be positive")
    return 8.0 * rho_kg_m3 * nu_m2_s * mean_velocity_m_s / radius_m**2


def source_force_vector(pressure_gradient_pa_m: float, axis: Iterable[float]) -> np.ndarray:
    """Musubi PIP_Force contract: glob_source.force is the physical gradient."""

    direction = np.asarray(tuple(axis), dtype=np.float64).reshape(3)
    direction /= np.linalg.norm(direction)
    return float(pressure_gradient_pa_m) * direction


def poiseuille_velocity_profile(
    radial_distance_m: np.ndarray | Sequence[float],
    *,
    radius_m: float,
    mean_velocity_m_s: float,
) -> np.ndarray:
    radial = np.asarray(radial_distance_m, dtype=np.float64)
    return 2.0 * float(mean_velocity_m_s) * (1.0 - (radial / float(radius_m)) ** 2)


def relative_l2(numerical: np.ndarray, analytic: np.ndarray) -> float:
    observed = np.asarray(numerical, dtype=np.float64)
    reference = np.asarray(analytic, dtype=np.float64)
    return float(np.linalg.norm(observed - reference) / np.linalg.norm(reference))


def effective_radius_ratio(flow_numerical: float, flow_exact: float) -> float:
    """Infer R_eff/R from Q proportional to R**4 at fixed gradient/viscosity."""

    if flow_numerical <= 0.0 or flow_exact <= 0.0:
        raise ValueError("Flow rates must be positive")
    return float(flow_numerical / flow_exact) ** 0.25


def early_stop_decision(
    samples: Sequence[dict[str, float]],
    *,
    window_iterations: int = 200,
    mean_drift_limit: float = 1.0e-4,
    profile_drift_limit: float = 1.0e-4,
) -> dict[str, Any]:
    if len(samples) < 2:
        return {"stop": False, "reason": "insufficient samples"}
    last_iteration = int(samples[-1]["iteration"])
    window = [row for row in samples if int(row["iteration"]) >= last_iteration - window_iterations]
    if len(window) < 2 or int(window[-1]["iteration"]) - int(window[0]["iteration"]) < window_iterations:
        return {"stop": False, "reason": "window shorter than required iterations"}

    def drift(key: str) -> float:
        values = np.asarray([float(row[key]) for row in window])
        return float(np.ptp(values) / max(abs(float(np.mean(values))), np.finfo(float).tiny))

    mean_drift = drift("mean_axial_velocity")
    profile_drift = drift("profile_l2_error")
    safety = (
        all(bool(row.get("all_finite", False)) for row in window)
        and min(float(row["minimum_pdf"]) for row in window) > 0.0
        and max(float(row["maximum_lattice_speed"]) for row in window) < 0.05
    )
    return {
        "stop": bool(mean_drift <= mean_drift_limit and profile_drift <= profile_drift_limit and safety),
        "mean_relative_drift": mean_drift,
        "profile_relative_drift": profile_drift,
        "numerical_safety": safety,
    }


def wall_stage_decision(axis_n27_error: float, *, oblique_n27_error: float | None = None) -> str:
    if axis_n27_error <= 0.02 and oblique_n27_error is not None and oblique_n27_error <= 0.02:
        return "WALL_GEOMETRY_ISOLATED_PASS"
    if axis_n27_error > 0.05:
        return "CFD_FLOW_WALL_GEOMETRY_ERROR_IDENTIFIED"
    return "CFD_FLOW_WALL_GEOMETRY_ERROR_MODERATE"


def observed_order_three(dx: Sequence[float], error: Sequence[float]) -> float | None:
    spacing = np.asarray(dx, dtype=np.float64)
    values = np.asarray(error, dtype=np.float64)
    if len(spacing) != 3 or len(values) != 3 or np.any(spacing <= 0.0) or np.any(values <= 0.0):
        return None
    return float(np.polyfit(np.log(spacing), np.log(values), 1)[0])
