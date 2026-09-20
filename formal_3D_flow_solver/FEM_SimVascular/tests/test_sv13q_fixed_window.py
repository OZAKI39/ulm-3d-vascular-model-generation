from sv13q_support import *
def test_exact_native_checkpoint_and_ten_steps():
 for d in cases():
  if d['mode']!='window':continue
  assert d['initial_checkpoint_sha256']==policy()['checkpoint_sha256']
  if d['status']=='PASS':assert d['steps']==10 and d['stop_step']==70 and d['start_step']==60
  assert case(d['candidate'],'SMOKE')['status']=='PASS'
def test_final_fields_and_reload_healthy():
 for d in cases():
  if d['status']=='PASS':
   assert d['linear_failures']==d['nonlinear_failures']==0
   assert d['reload']['velocity_finite'] and d['reload']['pressure_finite']
   assert d['measurement']['wall_noslip_pass'] and d['normal_exit']
