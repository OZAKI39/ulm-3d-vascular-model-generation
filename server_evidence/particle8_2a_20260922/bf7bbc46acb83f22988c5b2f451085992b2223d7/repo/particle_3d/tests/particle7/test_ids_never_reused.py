from particle_3d.lammps_state import BridgeParticle
from particle_3d.particle_shapes import Sphere
import pytest
def test_tombstone(engine):
 p=BridgeParticle.from_shape(7,Sphere([0,0,0],1e-6)); engine.bridge.insert(p); engine.bridge.remove(7)
 with pytest.raises(ValueError,match='NEVER_REUSED'): engine.bridge.insert(p)
