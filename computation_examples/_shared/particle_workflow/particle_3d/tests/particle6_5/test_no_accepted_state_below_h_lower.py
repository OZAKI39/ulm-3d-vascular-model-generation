import numpy as np,pytest
from particle_3d.nearfield_handoff import require_admissible_initial,InitialBelowContinuumHandoff,pair_handoff_certificate,wall_handoff_certificate
from particle_3d.particle_shapes import Sphere
from particle_3d.particle65_cases import synthetic

def test_accepted_and_continuous(wall_run,pair_run,real_saved):
 for d in [wall_run,pair_run,real_saved]:
  assert all(r['g_nf_m']>=-r['geometry_roundoff_m'] for r in d['interactions'])
  for s in d['states']:
   for p in s['projection']['continuous_certificates']:
    for k in ['minimum_g_nf_m','minimum_g_nf_bound_m']:
     if k in p:assert p[k]>=-p['roundoff_m']

def test_initial_below_rejected_without_movement(policy):
 shapes,_,wall,_=synthetic('wall');s=shapes[17].moved([0,0,1.001e-6]);before=s.center_m.copy()
 with pytest.raises(InitialBelowContinuumHandoff):require_admissible_initial({17:s},wall,policy,.00345312)
 assert np.array_equal(before,s.center_m)

def test_pair_positive_endpoints_can_cross():
 a=Sphere([-3e-6,0,0],1e-6);b=Sphere([0,0,0],1e-6);end=a.moved([3e-6,0,0])
 ok,p=pair_handoff_certificate(a,b,end,b,2e-9);assert not ok and p['fraction_of_minimum']==.5

def test_wall_positive_endpoints_can_cross():
 _,_,wall,_=synthetic('wall');a=Sphere([0,0,2e-6],1e-6);b=a.moved([0,0,-2e-6])
 ok,p=wall_handoff_certificate(a,b,wall,2e-9);assert not ok
