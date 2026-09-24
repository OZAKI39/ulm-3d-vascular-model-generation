import numpy as np
from conftest import balance

def test_constant_field_internal_flux(cube):
    u=np.tile([.002,0,0],(8,1));r=balance(cube,u);q=.002*cube['scale']**2
    assert np.isclose(r['Q_section_numpy_m3_s'],q,rtol=2e-13,atol=0)
    assert abs(r['closure_residual_m3_s'])<q*2e-13
    assert r['volume_integral_divergence_m3_s']==0
