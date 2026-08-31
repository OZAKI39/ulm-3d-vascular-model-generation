from __future__ import annotations

from utils.cfd_flow.tau1_grid_convergence import restart_resume_contract


def test_resume_requires_every_frozen_numerical_field() -> None:
    expected = {
        "mesh_hashes": {"a": "b"}, "dx_m": 1.0, "dt_s": 2.0,
        "rho_kg_m3": 3.0, "nu_m2_s": 4.0, "bulk_nu_m2_s": 5.0,
        "tau": 1.0, "boundary_contract": {"wall": "wall_libb"},
        "outlet_pressures_pa": {"o": 1.0}, "target_mass_flow_kg_s": 6.0,
        "binary_sha256": "sha", "layout": "d3q19", "relaxation": "bgk",
    }
    assert restart_resume_contract(expected, expected)["status"] == "PASS"
    changed = dict(expected)
    changed["dt_s"] = 3.0
    assert restart_resume_contract(changed, expected)["status"] == "FAIL"
