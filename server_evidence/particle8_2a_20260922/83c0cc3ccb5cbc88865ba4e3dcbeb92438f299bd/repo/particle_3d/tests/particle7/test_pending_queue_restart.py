from copy import deepcopy
from particle_3d.particle7_cases import channel
from particle_3d.injection_admission import FiniteSizeAdmission,PlugPathAdmission
from particle_3d.particle7_checkpoint import write_checkpoint,read_checkpoint

def test_pending_exact(engine,tmp_path):
 e=engine.scheduler.pop(); e['attempt_count']=17; e['last_candidate']={'position_m':[0,0,0],'status':'DEFORMATION_INFEASIBLE'}
 engine.events.append(deepcopy(e)); engine.pending['RBC'].append(e); engine.time_s=e['scheduled_time_s']
 write_checkpoint(tmp_path/'r',engine,{'test':'pending'})
 s,w,c,v=channel(); b=read_checkpoint(tmp_path/'r',lambda source:PlugPathAdmission(s,source,width=200e-6,velocity=v),c,{'test':'pending'})
 try:
  assert list(b.pending['RBC'])==[e]
  assert b.scheduler.next_rbc==engine.scheduler.next_rbc
  b.retry(); assert b.births[e['particle_id']]['geometry']==e['geometry']
 finally: b.bridge.close()

def test_restart_rejects_changed_inlet_model(engine,tmp_path):
 import pytest
 write_checkpoint(tmp_path/'r',engine,{'test':'model'})
 s,w,c,v=channel()
 with pytest.raises(ValueError,match='inlet/admission'):
  read_checkpoint(tmp_path/'r',lambda source:PlugPathAdmission(s,source,width=100e-6,velocity=v),c,{'test':'model'})
