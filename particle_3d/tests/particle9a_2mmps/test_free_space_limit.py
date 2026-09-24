import numpy as np
from particle_3d.planar_wall_hydrodynamics import planar_wall_affine_block,local_mobility_4d

def test_free_space_limit():
    for xi in [1.,2.]:
        r,b,d=planar_wall_affine_block(1e-6,.00345312,xi*1e-6,[0,0,1],np.eye(3),np.arange(6)*.001)
        assert np.array_equal(r,np.zeros((6,6))) and np.array_equal(b,np.zeros(6))
        m,free,_,_=local_mobility_4d(1e-6,.00345312,xi)
        np.testing.assert_array_equal(m,free)
