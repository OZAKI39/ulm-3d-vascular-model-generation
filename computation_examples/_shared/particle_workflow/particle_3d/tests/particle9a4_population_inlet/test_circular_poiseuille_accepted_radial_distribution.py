import numpy as np

def test_source_and_conditional_radial_cdfs(poiseuille):
 rad=np.array([np.linalg.norm(r['position_m'][:2])/2e-6 for r in poiseuille.rows])
 for x in [.2,.4,.6,.8]:assert abs(np.mean(rad<=x)-(2*x*x-x**4))<.015
 for d in [.6e-6,2.4e-6]:
  rows=[r for r in poiseuille.rows if r['diameter_m']==d and r['particle_id'] is not None]
  radius=np.array([np.linalg.norm(r['position_m'][:2])/2e-6 for r in rows]);xc=(2e-6-d/2-2e-9)/2e-6
  for t in [.25,.5,.75]:
   x=xc*t;expected=(2*x*x-x**4)/(2*xc*xc-xc**4)
   assert abs(np.mean(radius<=x)-expected)<.025
