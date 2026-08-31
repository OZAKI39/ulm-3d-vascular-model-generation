from __future__ import annotations

from pathlib import Path

from utils.cfd_flow.tau1_grid_convergence import build_physical_port_plane_contract


ROOT = Path(__file__).resolve().parents[1]


def test_physical_plane_hash_is_identical_for_all_three_grids() -> None:
    result = build_physical_port_plane_contract(ROOT)
    hashes = result["contract_hashes_by_grid"]
    assert result["status"] == "PASS"
    assert hashes["coarse"] == hashes["base"] == hashes["fine"]
