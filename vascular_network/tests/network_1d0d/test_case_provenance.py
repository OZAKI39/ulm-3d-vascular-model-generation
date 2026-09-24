"""Integration checks against the saved original A and current production mesh."""
from pathlib import Path
import csv
import json
import numpy as np
import pytest
from network_1d0d.audit import sha256
from network_1d0d.hydraulic_resistance import linear_radius_resistance

DATA = Path(__file__).resolve().parents[2]/'reports/a_network_1d0d_boundary_v1/data'


@pytest.fixture
def summary():
    if not (DATA/'final_summary.json').exists(): pytest.skip('Run the case provenance audit first')
    return json.loads((DATA/'final_summary.json').read_text())


def test_protected_source_and_production_hashes(summary):
    hashes = json.loads((DATA/'protected_input_hashes.json').read_text())
    assert all(sha256(path) == digest for path, digest in hashes.items())


def test_case_exact_cuts_preserve_every_raw_edge_resistance(summary):
    graph = np.load(DATA/'a_graph_with_exact_cuts_si.npz')
    e, x, r = graph['edges'], graph['xyz_m'], graph['radius_m']
    resistance = linear_radius_resistance(np.linalg.norm(x[e[:, 1]]-x[e[:, 0]], axis=1), r[e[:, 0]], r[e[:, 1]])
    with open(DATA/'a_edge_hydraulic_resistance.csv') as f:
        raw = np.array([float(row['R_Pa_s_m3']) for row in csv.DictReader(f)])
    grouped = np.bincount(graph['original_edge'], weights=resistance, minlength=len(raw))
    np.testing.assert_allclose(grouped, raw, rtol=1e-12, atol=0)


def test_blocked_case_never_exports_numeric_bc(summary):
    assert summary['status'] == 'A_NETWORK_SOURCE_AMBIGUOUS'
    bc = json.loads((DATA/'roi_fixed_pressure_bc_from_A.json').read_text())
    assert bc['status'] == 'BLOCKED_DO_NOT_USE'
    assert bc['raw_pressure_Pa'] is None and bc['source_A_solution_SHA'] is None
    assert summary['3D_validation_split'] is None
    assert summary['full_A_total_inflow'] is None


def test_all_saved_roi_cut_points_match_current_cap_planes(summary):
    for p in summary['roi_ports']:
        assert p['position_error_um'] < 1e-8 and p['radius_error_um'] < 1e-8
        assert p['fem_cap']['plane_offset_error_m'] < 2e-11
        assert p['fem_cap']['real_cut_to_cap_axial_um'] > 0
    mapping = json.loads((DATA/'roi_ports_in_a.json').read_text())
    assert mapping['extra_saved_roi_ports'] == 0 and not mapping['missing_incident_edges']


def test_uncertain_terminal_not_automatically_grounded(summary):
    with open(DATA/'a_terminal_audit.csv') as f: terminals = list(csv.DictReader(f))
    t = next(t for t in terminals if int(t['original_id']) == summary['roi_ports'][-1]['original_node_id'])
    assert t['terminal_class'] == 'UNKNOWN_TERMINAL'
    assert t['deep_interior_candidate'] == 'True'
    assert not summary['terminal_sensitivity']['MODEL_A']['valid']
