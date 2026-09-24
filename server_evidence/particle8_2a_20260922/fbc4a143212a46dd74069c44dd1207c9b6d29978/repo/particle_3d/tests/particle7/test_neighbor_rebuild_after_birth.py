from particle_3d.lammps_state import BridgeParticle
from particle_3d.particle_shapes import Sphere
from particle_3d.lammps_neighbors import standalone_candidates

def test_pairs(engine):
 for tag,x in [(91,0),(5,5e-6),(200,100e-6)]:
  engine.bridge.insert(BridgeParticle.from_shape(tag,Sphere([x,0,0],1e-6)))
  shapes={p.particle_id:p.shape() for p in engine.bridge.read()}
  assert engine.bridge.raw_pairs==standalone_candidates(shapes,engine.bridge.policy)
 assert engine.bridge.raw_pairs==[(5,91)]
