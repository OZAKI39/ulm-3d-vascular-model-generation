"""Pure analysis contracts for the gated interior-pressure pipe stage."""

from __future__ import annotations

from typing import Any, Iterable

import numpy as np


def interior_probe_positions(
    center_m: Iterable[float], axis: Iterable[float], diameter_m: float
) -> dict[str, Any]:
    center = np.asarray(tuple(center_m), dtype=np.float64).reshape(3)
    direction = np.asarray(tuple(axis), dtype=np.float64).reshape(3)
    direction /= np.linalg.norm(direction)
    diameter = float(diameter_m)
    if diameter <= 0.0:
        raise ValueError("diameter must be positive")
    a = center - diameter * direction
    b = center + diameter * direction
    return {
        "A_m": a,
        "B_m": b,
        "separation_m": 2.0 * diameter,
        "coordinates_along_axis_in_diameters": (-1.0, 1.0),
    }


def analytic_interior_delta_p(pressure_gradient_pa_m: float, diameter_m: float) -> float:
    return float(pressure_gradient_pa_m) * 2.0 * float(diameter_m)


def pressure_stage_decision(
    *, n16_delta_p: float, n20_delta_p: float, analytic_delta_p: float
) -> dict[str, Any]:
    error_n20 = abs(float(n20_delta_p) - float(analytic_delta_p)) / abs(float(analytic_delta_p))
    base_to_fine = abs(float(n20_delta_p) - float(n16_delta_p)) / abs(float(n16_delta_p))
    passed = error_n20 <= 0.02 and base_to_fine <= 0.02
    return {
        "status": "PRESSURE_EQ_BULK_EFFECT_PASS" if passed else "PRESSURE_EQ_NEEDS_N27",
        "n20_interior_delta_p_relative_error": error_n20,
        "n16_to_n20_delta_p_relative_difference": base_to_fine,
    }


def distinguish_endpoint_offset(
    *, boundary_delta_p_error: float, interior_delta_p_error: float
) -> str:
    if boundary_delta_p_error > 0.05 and interior_delta_p_error <= 0.02:
        return "BOUNDARY_ENDPOINT_PRESSURE_OFFSET"
    return "PRESSURE_BC_BULK_EFFECT_NOT_CLEARED"
