from __future__ import annotations

import math

import numpy as np

from utils.cfd_flow.musubi_pressure_bc_benchmark import (
    CANDIDATES,
    PIPE_DIAMETER_M,
    TARGET_Q_M3_S,
    benchmark_grids,
    compare_candidates,
    diffusive_time_step,
    hagen_poiseuille_delta_p,
    invalid_gci_guard,
    observed_order,
    parse_controller_output,
    pipe_frame,
)


def test_hagen_poiseuille_formula_and_fixed_geometry() -> None:
    radius = 2.0
    assert math.isclose(
        hagen_poiseuille_delta_p(length_m=3.0, flow_m3_s=5.0, radius_m=radius, dynamic_viscosity_pa_s=7.0),
        8.0 * 7.0 * 3.0 * 5.0 / (math.pi * radius**4),
    )
    assert PIPE_DIAMETER_M > 0.0 and TARGET_Q_M3_S > 0.0


def test_diffusive_scaling_and_registered_grids() -> None:
    grids = benchmark_grids()
    assert [grid.cells_across_diameter for grid in grids] == [12, 16, 20, 27]
    ratios = [grid.dt_s / grid.dx_m**2 for grid in grids]
    assert np.allclose(ratios, ratios[0], rtol=1e-14)
    assert diffusive_time_step(grids[1].dx_m) == grids[1].dt_s


def test_rotated_pipe_frame() -> None:
    axis, u, v = pipe_frame((0.2, 0.5, 0.8))
    matrix = np.vstack((axis, u, v))
    assert np.allclose(matrix @ matrix.T, np.eye(3), atol=1e-14)
    assert np.allclose(np.cross(u, v), axis, atol=1e-14)


def test_bc_source_contract_is_bounded_and_uses_actual_loader_tokens() -> None:
    assert list(CANDIDATES) == [
        "pressure_eq",
        "pressure_nonEqExpol",
        "pressure_antiBounceBack",
    ]
    assert CANDIDATES["pressure_nonEqExpol"] == "pressure_noneq_expol"
    assert CANDIDATES["pressure_antiBounceBack"] == "pressure_antibounceback"


def test_controller_parser() -> None:
    line = (
        "ADAPTIVE_FLUX_PRESSURE iter=12 target_lattice= 1.0E-3 controlled_lattice= 1.0E-3 "
        "relative_error= 2.0E-13 rho_boundary= 1.1 pressure_pa= 2.3E4 "
        "max_lattice_velocity= 1.2E-4 minimum_pdf= 4.0E-2 globBC_count=123"
    )
    row = parse_controller_output(line)[0]
    assert row["iteration"] == 12
    assert row["globBC_count"] == 123
    assert row["relative_error"] == 2e-13


def test_observed_order_and_invalid_gci() -> None:
    h = np.asarray((0.3, 0.2, 0.1))
    order = observed_order(h, h**2)
    assert order["available"] and math.isclose(order["value"], 2.0, rel_tol=1e-12)
    assert invalid_gci_guard((1.0, 0.9, 0.95), h)["available"] is False


def test_candidate_comparison_requires_all_gates() -> None:
    cases = []
    for candidate in CANDIDATES:
        for orientation in ("axis_aligned", "worst_real_outlet"):
            for index, grid in enumerate(benchmark_grids()):
                error = 0.08 / (index + 1)
                cases.append(
                    {
                        "candidate": candidate,
                        "orientation": orientation,
                        "dx_m": grid.dx_m,
                        "relative_delta_p_error": error,
                        "delta_p_numerical_pa": 100.0 * (1.0 - error),
                        "numerical_safety_pass": True,
                    }
                )
    report = compare_candidates(cases)
    assert report["benchmark_conclusive"] is False
    assert report["winner"] is None

