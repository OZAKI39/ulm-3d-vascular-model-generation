from __future__ import annotations

from pathlib import Path

from utils.cfd_flow.standardized_outlet_planes import (
    build_standardized_outlet_plane_contract,
)


ROOT = Path(__file__).resolve().parents[1]


def test_standardized_outlet_planes_are_safe_and_grid_independent() -> None:
    report = build_standardized_outlet_plane_contract(ROOT)
    assert report["status"] == "PASS"
    assert all(report["checks"].values())
    hashes = report["contract_hashes_by_grid"]
    assert hashes["coarse"] == hashes["base"] == hashes["fine"]
    for outlet in report["outlets"].values():
        assert outlet["plane_covers_only_corresponding_opening"]
        assert outlet["coverage_area_over_aperture_area"] < 1.5
        assert max(outlet["basis_qc"].values()) < 1e-12


def test_plane_contract_does_not_enable_extensions_or_vmtk_runtime() -> None:
    report = build_standardized_outlet_plane_contract(ROOT)
    assert report["flow_extension_used"] is False
    assert report["vmtk_used_at_runtime"] is False

