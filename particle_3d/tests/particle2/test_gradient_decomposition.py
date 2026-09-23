import numpy as np
from particle_3d.rbc_orientation import decompose_gradient


def test_symmetric_and_skew_parts_reconstruct_gradient():
    g=np.array([[2.,3.,-4.],[5.,-6.,7.],[8.,-9.,10.]])
    e,w=decompose_gradient(g)
    np.testing.assert_array_equal(e,e.T);np.testing.assert_array_equal(w,-w.T)
    np.testing.assert_array_equal(e+w,g)
    p=np.array([.2,.3,.4])
    curl=np.array([g[2,1]-g[1,2],g[0,2]-g[2,0],g[1,0]-g[0,1]])
    np.testing.assert_allclose(w@p,np.cross(.5*curl,p),rtol=0,atol=16*np.finfo(float).eps)
