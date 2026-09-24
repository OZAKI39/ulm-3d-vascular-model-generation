from copy import deepcopy
from particle_3d.particle_shapes import Sphere
from particle_3d.lammps_state import BridgeParticle

def test_block_then_unblock(engine):
 e=engine.scheduler.pop(); frozen=deepcopy(e); engine.events.append(deepcopy(e)); engine.pending['RBC'].append(e); engine.time_s=e['scheduled_time_s']
 blocker=BridgeParticle.from_shape(100000,Sphere([0,0,0],1e-3)); engine.active[100000]=blocker
 # Blocker is a validation exclusion object, never inserted into LAMMPS.
 engine.admission.guard=3; engine.retry(); assert len(engine.pending['RBC'])==1
 for k,v in frozen.items():
  if k!='attempt_count': assert e[k]==v
 del engine.active[100000]; engine.retry(); assert not engine.pending['RBC']
 assert e['particle_id'] in engine.active and engine.births[e['particle_id']]['q']==frozen['q']

def test_cannot_admit_before_scheduled_time(engine):
 e=engine.scheduler.pop(); engine.pending['RBC'].append(e); engine.retry()
 assert list(engine.pending['RBC'])==[e] and not engine.active
 engine.time_s=e['scheduled_time_s']; engine.retry()
 assert not engine.pending['RBC'] and e['particle_id'] in engine.active
