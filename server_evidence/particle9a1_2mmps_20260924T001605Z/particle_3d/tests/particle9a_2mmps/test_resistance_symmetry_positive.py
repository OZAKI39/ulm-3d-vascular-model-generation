import numpy as np
from particle_3d.planar_wall_hydrodynamics import planar_wall_affine_block

def test_resistance_symmetry_positive():
    a=1e-6;mu=.00345312;diag=np.r_[np.full(3,6*np.pi*mu*a),np.full(3,8*np.pi*mu*a**3)]
    D=np.diag(1/np.sqrt(diag))
    for xi in [.001,.003,.01,.1,.5,1.]:
        r,_,_=planar_wall_affine_block(a,mu,xi*a,[.3,.4,np.sqrt(.75)],np.zeros((3,3)),np.zeros(6))
        R=D@(np.diag(diag)+r)@D
        np.testing.assert_allclose(R,R.T,rtol=0,atol=2e-14)
        assert np.linalg.eigvalsh(R).min()>0
