"""Permanent analytic validation for steady Newtonian networks."""
import numpy as np
import pytest
from scipy.integrate import quad
from network_1d0d.hydraulic_resistance import linear_radius_resistance


@pytest.mark.parametrize('ratio', [1., 1+1e-12, .3, 2., 10.])
def test_linear_radius_segment(ratio):
    r, L, mu = 1.4e-6, 37e-6, .00345312
    dimensionless, _ = quad(lambda t: (1+(ratio-1)*t)**-4, 0, 1, epsabs=1e-13, epsrel=1e-13)
    expected = 8*mu*L/(np.pi*r**4)*dimensionless
    assert linear_radius_resistance(L, r, ratio*r, mu) == pytest.approx(expected, rel=2e-13)


def test_global_radius_resistance_scaling():
    for scale in [.90, .95, 1., 1.05, 1.10]:
        a = linear_radius_resistance(50e-6, 1.1e-6, 2.4e-6)
        b = linear_radius_resistance(50e-6, scale*1.1e-6, scale*2.4e-6)
        assert b/a == pytest.approx(scale**-4, rel=1e-14)
