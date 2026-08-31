from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN_NAME = "healthy_mouse_capillary_tau1_grid_convergence_anchor003274_20260831"


def test_old_cell_cube_flux_has_evidence_but_no_live_code() -> None:
    path = ROOT / "outputs" / "cfd_flow" / RUN_NAME / "qc" / "cell_cube_plane_aperture_clipping_v1_forensic.json"
    record = json.loads(path.read_text(encoding="utf-8"))
    assert record["status"] == "HISTORICAL_FAILED_FLUX_EXTRACTOR"
    assert record["classification"] == "RETIRED_FOR_PHYSICAL_PORT_FLUX"
    assert record["root_cause"] == "LBM_NODE_IS_NOT_PHYSICAL_CONTROL_VOLUME"
    assert not record["historical_code_preserved"]
    assert record["retirement_evidence_preserved"]
    assert all(port["coverage_fraction"] < 0.20 for port in record["ports"].values())
