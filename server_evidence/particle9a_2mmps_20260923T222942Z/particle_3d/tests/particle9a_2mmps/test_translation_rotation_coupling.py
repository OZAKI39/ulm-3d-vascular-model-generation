import numpy as np
from particle_3d.planar_wall_hydrodynamics import local_mobility_4d,stable_wall_frame

def test_translation_rotation_coupling():
    n=np.array([.6,0,.8]);t=np.array([-.8,0,.6]);T=stable_wall_frame(n)[:,:2]
    m,_,_,_=local_mobility_4d(1e-6,.00345312,.01)
    q=m@np.r_[T.T@t*1e-12,[0,0]]
    v=T@q[:2];omega=T@q[2:]/1e-6
    assert v@t>0 and omega[1]>0
    assert abs(omega@n)<1e-12
