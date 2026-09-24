import numpy as np
from particle_3d.planar_wall_hydrodynamics import local_mobility_4d,stable_wall_frame

def test_rotation_translation_coupling():
    n=np.array([.6,0,.8]);t=np.array([-.8,0,.6]);T=stable_wall_frame(n)[:,:2]
    m,_,_,_=local_mobility_4d(1e-6,.00345312,.01)
    q=m@np.r_[[0,0],T.T@np.array([0,1e-18,0])/1e-6]
    assert (T@q[:2])@t>0 and (T@q[2:])[1]>0
