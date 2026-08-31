from __future__ import annotations

from pathlib import Path

from utils.cfd_flow.physical_port_flux import (
    MINIMUM_CLEARANCE_M,
    build_interior_plane_contract,
)


ROOT = Path(__file__).resolve().parents[1]


def test_v3_planes_are_geometry_only_single_lumen_sections() -> None:
    contract = build_interior_plane_contract(ROOT)
    assert contract["status"] == "PASS"
    assert len(set(contract["contract_hashes_by_grid"].values())) == 1
    for port in contract["ports"].values():
        central = port["planes"]["central"]
        assert central["centerline_arclength_from_port_m"] >= MINIMUM_CLEARANCE_M
        assert central["distance_to_nearest_bifurcation_m"] >= MINIMUM_CLEARANCE_M
        for plane in port["planes"].values():
            assert plane["slice_qc"]["nontrivial_local_secondary_component_count"] == 0
            assert plane["slice_qc"]["polygon_valid"]
