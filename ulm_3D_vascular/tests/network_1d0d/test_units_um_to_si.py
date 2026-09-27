"""Permanent analytic validation for steady Newtonian networks."""
import numpy as np
import pytest
from network_1d0d.hydraulic_resistance import linear_radius_resistance, um_to_m


def test_units_um_to_si():
    r = linear_radius_resistance(um_to_m(100), um_to_m(2), um_to_m(2))
    assert r == pytest.approx(8*.00345312*100e-6/(np.pi*(2e-6)**4), rel=1e-14)
    np.testing.assert_allclose(um_to_m([1, 2, 3]), [1e-6, 2e-6, 3e-6], atol=0)
