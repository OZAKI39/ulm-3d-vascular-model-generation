import numpy as np
from network_1d0d.boundary_conditions import shift_pressure_gauge

def test_signed_raw_cap_pressure_common_shift(bc):
    p=np.array([r['pressure_cap_raw_Pa'] for r in bc['ports']])
    assert p[2]<0
    shifted=shift_pressure_gauge(p)
    assert shifted.min()==0
    np.testing.assert_allclose(shifted[:,None]-shifted,p[:,None]-p,rtol=1e-15,atol=1e-12)
