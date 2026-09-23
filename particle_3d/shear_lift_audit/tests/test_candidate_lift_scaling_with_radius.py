import numpy as np
from audit_math import candidate

def test_radius_squared():
    f=candidate(1e-6,1e-5,100,.003,1000)
    np.testing.assert_allclose(candidate(np.array([.5,1,2])*1e-6,1e-5,100,.003,1000)/f,[.25,1,4])
