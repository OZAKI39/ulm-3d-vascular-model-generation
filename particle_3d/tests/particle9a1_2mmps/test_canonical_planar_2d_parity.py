import numpy as np
from particle_3d.planar_wall_hydrodynamics import planar_wall_affine_block

def test_canonical_matches_original(old_planar):
 n=np.array([.6,0,.8]);t=np.array([-.8,0,.6]);g=500*np.outer(t,n)
 for xi in [.001,.003,.01,.1,.5,1.,2.]:
  u=np.r_[500e-6*(1+xi)*t,250*np.cross(n,t)]
  old=old_planar.planar_wall_affine_block(1e-6,.00345312,xi*1e-6,n,g,u)
  new=planar_wall_affine_block(1e-6,.00345312,xi*1e-6,n,g,u)
  np.testing.assert_array_equal(old[0],new[0])
  np.testing.assert_allclose(old[1],new[1],rtol=2e-14,atol=1e-28)
  for key in ['target_tangential_velocity_xyz','target_tangential_omega_xyz']:
   np.testing.assert_allclose(old[2][key],new[2][key],rtol=2e-14,atol=1e-18)
