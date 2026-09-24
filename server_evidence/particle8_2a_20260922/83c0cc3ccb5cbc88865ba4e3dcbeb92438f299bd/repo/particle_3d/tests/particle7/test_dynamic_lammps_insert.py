from particle_3d.lammps_state import BridgeParticle
from particle_3d.particle_shapes import Ellipsoid,Capsule
from particle_3d.injection_admission import event_shape
import numpy as np
def test_properties(engine,source):
 e=source.next_rbc(); shape=event_shape(e,[0,0,0]); p=BridgeParticle.from_shape(801,shape,q=e['q'],velocity=[.1,.2,.3],omega=[1,2,3]); engine.bridge.insert(p)
 a=engine.bridge.read()[0]
 for k in p.to_dict(): assert np.array_equal(getattr(p,k),getattr(a,k))
 cap=BridgeParticle.from_shape(17,Capsule([20e-6,0,0],[1,0,0],1e-6,2e-6)); engine.bridge.insert(cap)
 assert [x.particle_id for x in engine.bridge.read()]==[17,801]
