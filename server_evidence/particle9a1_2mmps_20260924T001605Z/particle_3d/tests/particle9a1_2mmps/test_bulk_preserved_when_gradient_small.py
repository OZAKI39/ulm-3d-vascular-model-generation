import numpy as np
from particle_3d.planar_wall_hydrodynamics import planar_wall_affine_block

def test_bulk_survives_small_gradient():
 for gamma in [0.,1e-3]:
  g=np.zeros((3,3));g[0,2]=gamma
  _,_,d=planar_wall_affine_block(1e-6,.003,1e-8,[0,0,1],g,[.002,0,0,0,0,0])
  assert np.linalg.norm(d['target_tangential_velocity_xyz'])>.002/100
