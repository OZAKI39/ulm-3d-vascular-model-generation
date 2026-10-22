import numpy as np,pytest
from mouse_vessel_flow.particles.nearfield_handoff import require_admissible_initial,InitialBelowContinuumHandoff,pair_handoff_certificate,wall_handoff_certificate
from mouse_vessel_flow.particles.particle_shapes import Sphere
from mouse_vessel_flow.particles.particle65_cases import synthetic

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
