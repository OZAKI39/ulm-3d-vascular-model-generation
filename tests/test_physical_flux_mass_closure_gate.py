from __future__ import annotations

from utils.cfd_flow.tau1_grid_convergence import physical_flux_mass_closure


def test_physical_flux_mass_closure_uses_unrenormalized_flows() -> None:
    assert physical_flux_mass_closure(10.0, (2.0, 3.0, 5.0))["pass"]
    failed = physical_flux_mass_closure(10.0, (2.0, 3.0, 4.8))
    assert not failed["pass"]
    assert failed["relative_error"] > 0.01
