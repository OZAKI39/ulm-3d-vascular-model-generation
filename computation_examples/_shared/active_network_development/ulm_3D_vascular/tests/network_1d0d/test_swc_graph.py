"""Permanent analytic validation for steady Newtonian networks."""
import numpy as np
import pytest
from network_1d0d.network_solver import solve_network
from network_1d0d.swc_graph import load_swc, audit_graph


def test_cycles_preserved(tmp_path):
    path = tmp_path/'cycle.swc'
    path.write_text('10 0 0 0 0 1 30\n20 0 1 0 0 1 10\n30 0 1 1 0 1 20\n')
    graph = load_swc(path, [1, 1, 2]); a = audit_graph(graph)
    assert a['cycle_rank'] == 1 and a['edges'] == 3
    s = solve_network(3, graph.edges, np.ones(3), {0: 1., 2: 0.})
    assert s.pressure[1] == pytest.approx(.5)


def test_swc_coordinate_transform_and_invalid_ids(tmp_path):
    path = tmp_path/'small.swc'; path.write_text('10 0 1 2 3 4 -1\n20 0 2 3 4 5 10\n')
    g = load_swc(path, [1, 1, 2]); np.testing.assert_array_equal(g.xyz_um[0], [1, 2, 6])
    assert g.radius_um[0] == 4
    path.write_text('10 0 0 0 0 1 -1\n10 0 1 0 0 1 10\n')
    with pytest.raises(ValueError, match='Duplicate'): load_swc(path, [1, 1, 1])
