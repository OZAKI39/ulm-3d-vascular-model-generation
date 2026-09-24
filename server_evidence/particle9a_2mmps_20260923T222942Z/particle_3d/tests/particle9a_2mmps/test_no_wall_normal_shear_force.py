import numpy as np
from particle_3d.planar_wall_hydrodynamics import planar_wall_affine_block

def test_no_wall_normal_shear_force():
    for n in [[0,0,1],[.6,0,.8],[.3,.4,np.sqrt(.75)]]:
        _,_,d=planar_wall_affine_block(1e-6,.003,1e-8,n,np.arange(9).reshape(3,3)*100,np.zeros(6))
        f=np.array(d['shear_force_N'])
        assert abs(f@np.array(n))<=32*np.finfo(float).eps*np.linalg.norm(f)
