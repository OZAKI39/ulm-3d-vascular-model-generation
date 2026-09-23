import numpy as np
from audit_math import dimensionless,local_tensors

def test_radius_reynolds_and_shear_tensor():
    g=np.zeros((2,3,3));g[0,0,1]=20;g[1,0,0]=3;g[1,1,1]=-3
    e,w,curl,gamma=local_tensors(g)
    np.testing.assert_allclose(gamma,[20,6])
    np.testing.assert_allclose(curl[0],[0,0,-20])
    np.testing.assert_allclose(e+w,g)
    z=dimensionless(1e-6,1e-4,100.,1e-6,1e-6)
    assert np.isclose(z['Re_p'],1e-4)
    assert np.isclose(z['Re_G'],1e-4)
    assert np.isclose(z['wall_margin'],.02)
