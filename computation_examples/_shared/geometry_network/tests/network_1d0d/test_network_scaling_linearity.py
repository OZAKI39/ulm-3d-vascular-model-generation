"""Permanent analytic validation for steady Newtonian networks."""
import numpy as np
import pytest
from network_1d0d.network_solver import solve_network, operating_point_scale


def test_network_scaling_linearity():
    edges, R = [[0, 1], [1, 2]], np.array([2., 3.])*1e17
    a = solve_network(3, edges, R, {0: 1., 2: 0.})
    scale = operating_point_scale(a.flow[1], 1.551359160440232e-14)
    b = solve_network(3, edges, R, {0: scale, 2: 0.})
    np.testing.assert_allclose(b.pressure, scale*a.pressure, rtol=1e-13)
    np.testing.assert_allclose(b.flow, scale*a.flow, rtol=1e-13, atol=0)


@pytest.mark.parametrize('q', [0., -1., float('nan')])
def test_scaling_does_not_hide_reversed_inlet(q):
    with pytest.raises(ValueError): operating_point_scale(q, 1.)
