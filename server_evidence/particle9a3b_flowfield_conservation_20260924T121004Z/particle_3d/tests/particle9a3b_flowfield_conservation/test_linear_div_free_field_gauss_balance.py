import numpy as np
from conftest import balance

def test_linear_div_free_field_gauss_balance(cube):
    x=cube['unit'];u=.001*np.column_stack([.4+x[:,0],-x[:,1],np.zeros(8)])
    r=balance(cube,u);scale=.001*cube['scale']**2
    assert np.isclose(r['Q_section_numpy_m3_s'],.77*scale,rtol=2e-13,atol=0)
    assert abs(r['volume_integral_divergence_m3_s'])<scale*2e-13
    assert abs(r['closure_residual_m3_s'])<scale*2e-13
    assert np.isclose(r['surface_fluxes']['WALL']['Q_m3_s'],-.37*scale,rtol=2e-13,atol=0)
