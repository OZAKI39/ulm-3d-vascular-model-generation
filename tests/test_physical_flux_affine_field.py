from __future__ import annotations

from utils.cfd_flow.physical_port_flux import run_analytic_flux_oracles


def test_affine_oblique_field_oracle() -> None:
    record = run_analytic_flux_oracles()["affine_field"]
    assert record["relative_error"] <= 1.0e-8
    assert record["stencil_qc"]["high"]["invalid_stencil_count"] == 0
