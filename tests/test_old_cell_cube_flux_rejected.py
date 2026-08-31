from __future__ import annotations

import json
from pathlib import Path

from utils.cfd_flow.physical_port_flux import (
    LEGACY_CLASSIFICATION,
    LEGACY_RETIREMENT,
    LEGACY_ROOT_CAUSE,
    RUN_NAME,
)


ROOT = Path(__file__).resolve().parents[1]


def test_old_cell_cube_flux_is_preserved_but_rejected() -> None:
    path = ROOT / "outputs" / "cfd_flow" / RUN_NAME / "qc" / "cell_cube_plane_aperture_clipping_v1_forensic.json"
    record = json.loads(path.read_text(encoding="utf-8"))
    assert record["status"] == LEGACY_CLASSIFICATION
    assert record["classification"] == LEGACY_RETIREMENT
    assert record["root_cause"] == LEGACY_ROOT_CAUSE
    assert record["historical_code_preserved"]
    assert all(port["coverage_fraction"] < 0.20 for port in record["ports"].values())
