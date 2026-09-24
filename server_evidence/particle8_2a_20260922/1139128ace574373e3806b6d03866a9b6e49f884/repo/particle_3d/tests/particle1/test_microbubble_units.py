import numpy as np
import pytest
from particle_3d.sonovue_adapter import diameter_um_to_radius_m
from particle_3d.microbubble import stokes_drag_force_n


def test_diameter_to_radius_si_and_stokes_dimensions():
    assert diameter_um_to_radius_m(2.) == np.float64(1e-6)
    force=stokes_drag_force_n([1.,-2.,3.],[0.,0.,0.],.003,1e-6)
    np.testing.assert_allclose(force, -6*np.pi*.003*1e-6*np.array([1.,-2.,3.]),rtol=8*np.finfo(float).eps,atol=0)


@pytest.mark.parametrize('diameter',[0.,-1.,np.nan,np.inf])
def test_invalid_diameter_rejected(diameter):
    with pytest.raises(ValueError): diameter_um_to_radius_m(diameter)
