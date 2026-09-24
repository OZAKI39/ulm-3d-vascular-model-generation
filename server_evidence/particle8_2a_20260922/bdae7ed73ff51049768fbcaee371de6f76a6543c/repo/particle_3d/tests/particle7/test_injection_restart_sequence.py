from particle_3d.particle7_cases import synthetic_engine,channel
from particle_3d.injection_admission import FiniteSizeAdmission,PlugPathAdmission
from particle_3d.particle7_checkpoint import write_checkpoint,read_checkpoint
from particle_3d.particle3_cases import json_safe
import json

def test_binary_exact(tmp_path):
 a=synthetic_engine(); b=synthetic_engine()
 try:
  a.step_to(.002); b.step_to(.002); write_checkpoint(tmp_path/'restart',b,{'test':'P7'}); b.bridge.close()
  s,w,c,v=channel()
  b=read_checkpoint(tmp_path/'restart',lambda source:PlugPathAdmission(s,source,width=200e-6,velocity=v),c,{'test':'P7'})
  a.step_to(.004); b.step_to(.004)
  for key in ['events','births','exits','pending','active','scheduler']:
   assert json.dumps(a.state()[key],sort_keys=True,default=json_safe)==json.dumps(b.state()[key],sort_keys=True,default=json_safe),key
 finally: a.bridge.close(); b.bridge.close()
