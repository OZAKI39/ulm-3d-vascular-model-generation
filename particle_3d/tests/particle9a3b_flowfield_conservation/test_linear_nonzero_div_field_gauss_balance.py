import numpy as np
from conftest import balance

def test_linear_nonzero_div_field_gauss_balance(cube):
    u=.001*cube['unit']*[1,2,3];r=balance(cube,u);scale=.001*cube['scale']**2
    assert np.isclose(r['volume_integral_divergence_m3_s'],6*.37*scale,rtol=2e-13,atol=0)
    assert np.isclose(r['Q_section_numpy_m3_s'],.37*scale,rtol=2e-13,atol=0)
    assert abs(r['closure_residual_m3_s'])<scale*2e-13
