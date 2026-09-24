from __future__ import annotations

from utils.cfd_flow.physical_port_flux import run_analytic_flux_oracles


def test_oblique_poiseuille_cbf_oracle() -> None:
    records = run_analytic_flux_oracles()["poiseuille_oblique"]
    assert records["coarse"]["relative_error"] <= 0.01
    assert records["base"]["relative_error"] <= 0.005
    assert records["fine"]["relative_error"] <= 0.005
    assert records["fine"]["relative_error"] <= records["coarse"]["relative_error"]
