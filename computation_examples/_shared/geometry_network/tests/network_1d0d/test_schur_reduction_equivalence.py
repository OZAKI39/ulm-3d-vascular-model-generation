"""Permanent analytic validation for steady Newtonian networks."""
import numpy as np
from network_1d0d.network_solver import solve_network
from network_1d0d.network_reduction import reduce_boundary


def test_schur_reduction_equivalence():
    # Two ports, an internal junction, grounded terminal and nonzero-pressure source.
    edges = [[0, 2], [1, 2], [2, 3], [2, 4], [0, 1]]
    R = np.array([2., 3., 7., 5., 11.])*1e17
    reduced = reduce_boundary(5, edges, R, [0, 1], {3: 0., 4: 13.})
    for p0, p1 in [(1., 2.), (9., -3.), (20., 15.)]:
        full = solve_network(5, edges, R, {0: p0, 1: p1, 3: 0., 4: 13.})
        np.testing.assert_allclose(full.node_outflow[:2], reduced['Y'] @ [p0, p1]+reduced['forcing'], rtol=1e-13, atol=1e-30)
    assert reduced['passive'] and reduced['symmetry_error'] < 1e-14
    np.testing.assert_allclose(reduced['Z'] @ reduced['Y'], np.eye(2), atol=1e-14)
    assert np.max(abs(reduced['forcing'])) > 0


def test_singular_passive_operator_has_no_resistance_inverse():
    result = reduce_boundary(2, [[0, 1]], np.array([1e17]), [0, 1], {})
    assert result['passive'] and result['Z'] is None
