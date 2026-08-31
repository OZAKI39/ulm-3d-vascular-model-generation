from __future__ import annotations

from utils.cfd_flow.tau1_base import restart_resume_contract


def test_restart_resumes_only_when_entire_contract_matches() -> None:
    expected = {
        "mesh_hashes": {"qval.lsb": "continuous-q"},
        "dx_m": 2.0e-7,
        "dt_s": 2.038735983690112e-9,
        "rho_kg_m3": 1056.0,
        "nu_m2_s": 3.27e-6,
        "bulk_nu_m2_s": 2.18e-6,
        "boundary_contract": {"wall": "wall_libb"},
        "outlet_gauge_pressure_pa": {"outlet_01": 14.544978101},
        "target_mass_flow_kg_s": 2.890180380479642e-12,
        "binary_sha256": "e801",
    }
    assert restart_resume_contract(expected, expected)["status"] == "PASS"
    stale = dict(expected, dt_s=2.44140625e-8)
    result = restart_resume_contract(stale, expected)
    assert result["status"] == "FAIL"
    assert result["checks"]["dt_s"] is False
