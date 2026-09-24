from particle_3d.lammps_state import BridgeParticle
from particle_3d.particle_shapes import Sphere

def test_removed_from_pairs(engine):
 for tag in [2,9]: engine.bridge.insert(BridgeParticle.from_shape(tag,Sphere([tag*1e-6,0,0],1e-6)))
 assert engine.bridge.raw_pairs==[(2,9)]
 engine.bridge.remove(2); assert engine.bridge.raw_pairs==[] and [p.particle_id for p in engine.bridge.read()]==[9]
