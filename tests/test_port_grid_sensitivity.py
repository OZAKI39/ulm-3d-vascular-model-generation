from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from utils.cfd_flow.port_grid_sensitivity import (
    FLOW_SPLITS,
    RESEARCH_RUN,
    compare_grid_ports,
    orthonormal_plane_basis,
    recover_continuous_ports,
)


ROOT = Path(__file__).resolve().parents[1]


def test_plane_basis_is_orthonormal() -> None:
    normal = np.asarray((0.21, -0.62, 0.71))
    normal /= np.linalg.norm(normal)
    u, v = orthonormal_plane_basis(normal)
    assert np.isclose(np.linalg.norm(u), 1.0)
    assert np.isclose(np.linalg.norm(v), 1.0)
    assert abs(float(np.dot(normal, u))) < 1e-14
    assert abs(float(np.dot(normal, v))) < 1e-14
    assert abs(float(np.dot(u, v))) < 1e-14


def test_continuous_port_recovery_uses_real_geometry() -> None:
    frames, ports = recover_continuous_ports(ROOT)
    assert set(frames) == {"inlet", "outlet_01", "outlet_02", "outlet_03"}
    assert 7.7e-12 < ports["inlet"]["continuous_area_m2"] < 7.9e-12
    assert all(ports[label]["perimeter_m"] > 0.0 for label in ports)
    assert all(ports[label]["branch_distance"]["distance_to_nearest_bifurcation_um"] > 0.0 for label in ports)


def test_existing_discrete_port_report_contains_pressure_valid_metrics() -> None:
    path = ROOT / "outputs" / "cfd_flow" / RESEARCH_RUN / "qc" / "port_grid_sensitivity_forensics.json"
    report = json.loads(path.read_text(encoding="utf-8"))
    for grid in ("coarse", "base", "fine"):
        for port in ("inlet", "outlet_01", "outlet_02", "outlet_03"):
            row = report["grids"][grid]["ports"][port]
            assert 0 < row["pressure_valid_cell_count"] <= row["boundary_cell_count"]
            assert row["seeder_boundary_cell_count_all_lattice_sides"] >= row["d3q19_globBC_count"]
            assert 0.0 < row["pressure_valid_fraction"] <= 1.0
            assert row["pressure_valid_count_dx2_proxy_m2"] > 0.0
            assert row["true_plane_signed_cell_center_distance_over_dx"]["maximum_absolute"] >= 0.0
    assert report["seeder_calls"] == report["musubi_calls"] == report["harvester_calls"] == 0
    assert report["extra_fine_vascular_cfd_calls"] == 0


def test_flow_split_values_are_frozen() -> None:
    assert FLOW_SPLITS["coarse"]["outlet_02"] == 0.3455414
    assert FLOW_SPLITS["base"]["outlet_03"] == 0.4400351
    assert FLOW_SPLITS["fine"]["outlet_02"] == 0.5525517
