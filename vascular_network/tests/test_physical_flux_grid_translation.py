from __future__ import annotations

from utils.cfd_flow.physical_port_flux import run_analytic_flux_oracles


def test_base_grid_translation_sensitivity() -> None:
    record = run_analytic_flux_oracles()["grid_translation"]
    assert set(record["cases"]) == {"0.0", "0.25", "0.5"}
    assert record["relative_spread"] <= 0.005
