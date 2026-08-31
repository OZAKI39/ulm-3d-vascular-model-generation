from __future__ import annotations

import json
from pathlib import Path

from utils.cfd_flow.physical_port_flux import RUN_NAME


ROOT = Path(__file__).resolve().parents[1]


def test_base_gate_reports_first_real_failure_without_solver_calls() -> None:
    path = ROOT / "outputs" / "cfd_flow" / RUN_NAME / "qc" / "base_physical_flux_preflight_v2.json"
    record = json.loads(path.read_text(encoding="utf-8"))
    assert record["gates"]["analytic_flux_oracles"]
    assert record["gates"]["physical_aperture_area_error_le_1e-10"]
    assert record["gates"]["quadrature_convergence_le_1e-3"]
    assert record["gates"]["all_interpolation_stencils_valid"]
    assert record["true_first_scientific_failure"] == "BASE_QIN_TARGET_FAILED"
    assert record["seeder_calls"] == record["normal_musubi_calls"] == 0
    assert record["instrumented_musubi_calls"] == record["long_cfd_iterations"] == 0
