"""Permanent analytic validation for steady Newtonian networks."""
import numpy as np
from network_1d0d.network_solver import solve_network


def test_series_resistance():
    s = solve_network(3, [[0, 1], [1, 2]], np.array([2e17, 7e17]), {0: 90., 2: 0.})
    np.testing.assert_allclose(s.flow, 1e-16, rtol=1e-13, atol=0)
