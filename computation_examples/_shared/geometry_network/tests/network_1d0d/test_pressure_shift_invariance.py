"""Permanent analytic validation for steady Newtonian networks."""
import numpy as np
from network_1d0d.network_solver import solve_network
from network_1d0d.boundary_conditions import shift_pressure_gauge


def test_pressure_shift_invariance():
    edges, R = [[0, 1], [1, 2]], np.array([2., 3.])*1e17
    a = solve_network(3, edges, R, {0: 10., 2: 0.})
    b = solve_network(3, edges, R, {0: 51., 2: 41.})
    np.testing.assert_allclose(a.flow, b.flow, rtol=2e-14, atol=0)
    p = np.array([2., 10., -5.]); shifted = shift_pressure_gauge(p)
    np.testing.assert_allclose(p[:, None]-p, shifted[:, None]-shifted)
