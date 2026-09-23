import numpy as np
from audit_math import candidate

def test_square_root_shear_linear_slip():
    f=candidate(1e-6,1e-5,100,.003,1000)
    np.testing.assert_allclose(candidate(1e-6,1e-5,np.array([.25,1,4])*100,.003,1000)/f,[.5,1,2])
    assert np.isclose(candidate(1e-6,2e-5,100,.003,1000),2*f,atol=0)
