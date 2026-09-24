from particle_3d.lammps_state import BridgeParticle
from particle_3d.particle_shapes import Sphere

def test_not_local_index(engine):
 for tag in [100,7,305]: engine.bridge.insert(BridgeParticle.from_shape(tag,Sphere([tag*1e-8,0,0],1e-7)))
 engine.bridge.remove(7); assert [p.particle_id for p in engine.bridge.read()]==[100,305]
 engine.bridge.remove(100); engine.bridge.remove(305); assert engine.bridge.read()==[]
 assert 'delete_atoms group p7_delete compress no' in engine.bridge.commands
