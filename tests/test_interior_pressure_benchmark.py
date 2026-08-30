from __future__ import annotations

import math

import numpy as np

from utils.cfd_flow.pipe_interior_pressure_benchmark import (
    analytic_interior_delta_p,
    distinguish_endpoint_offset,
    interior_probe_positions,
    pressure_stage_decision,
)


def test_interior_probe_positions_are_minus_and_plus_one_diameter() -> None:
    probes = interior_probe_positions((1.0, 2.0, 3.0), (0.0, 0.0, 2.0), 0.5)
    assert np.allclose(probes["A_m"], (1.0, 2.0, 2.5))
    assert np.allclose(probes["B_m"], (1.0, 2.0, 3.5))
    assert math.isclose(probes["separation_m"], 1.0)


def test_interior_pressure_semantics_and_decision() -> None:
    assert math.isclose(analytic_interior_delta_p(10.0, 0.5), 10.0)
    result = pressure_stage_decision(n16_delta_p=9.9, n20_delta_p=10.0, analytic_delta_p=10.0)
    assert result["status"] == "PRESSURE_EQ_BULK_EFFECT_PASS"
    assert distinguish_endpoint_offset(
        boundary_delta_p_error=0.35, interior_delta_p_error=0.01
    ) == "BOUNDARY_ENDPOINT_PRESSURE_OFFSET"


def test_n16_to_n20_is_base_to_fine() -> None:
    result = pressure_stage_decision(n16_delta_p=100.0, n20_delta_p=110.0, analytic_delta_p=110.0)
    assert math.isclose(result["n16_to_n20_delta_p_relative_difference"], 0.1)
    assert result["status"] == "PRESSURE_EQ_NEEDS_N27"
