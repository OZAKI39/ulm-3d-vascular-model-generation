import numpy as np
from particle_3d.particle_shapes import Sphere
from particle_3d.hydrodynamic_resistance import PhysicalNearField

def test_transition_limits(policy):
 for x in [.01,.05]:
  values=[policy.evaluate(PhysicalNearField(1,None,(x+d)*1e-6,[0,0,1],0),{1:Sphere([0,0,0],1e-6)},.00345312)[0] for d in [-1e-9,0,1e-9]]
  assert max(values)-min(values)<2e-12
