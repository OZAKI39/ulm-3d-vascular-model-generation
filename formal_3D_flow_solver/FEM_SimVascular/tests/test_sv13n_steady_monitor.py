from sv13n_support import *
from sv_validation.sv13 import stop_gate
def policy():return json.loads((ROOT/'configs/sv1_3/policy.json').read_text())
def state():return dict(epsilon_Q=0.,epsilon_mass=0.,velocity_finite=True,pressure_finite=True,wall_noslip_pass=True)
def test_single_interval_is_not_steady():
 assert not stop_gate([dict(E_u=0.,E_Q=0.)],state(),0,0,policy())
def test_full_consecutive_requirement_and_failures():
 rows=[dict(E_u=0.,E_Q=0.) for _ in range(policy()['steady_last_intervals'])]
 assert stop_gate(rows,state(),0,0,policy())
 assert not stop_gate(rows,state(),1,0,policy())
 assert not stop_gate(rows,state(),0,1,policy())
 rows[-2]['E_u']=policy()['velocity_change_limit']*1.01
 assert not stop_gate(rows,state(),0,0,policy())
@pytest.mark.parametrize('key',['velocity_change_limit','flow_change_limit','mass_limit','steady_last_intervals','save_interval_steps'])
def test_changed_production_criterion_rejected(key):
 p=policy();p[key]*=2
 with pytest.raises(GateError):production_monitor_policy_gate(p,policy())
def test_nan_field_cannot_stop():
 rows=[dict(E_u=0.,E_Q=0.) for _ in range(5)];bad=state();bad['velocity_finite']=False
 assert not stop_gate(rows,bad,0,0,policy())
def test_actual_first_full_steady_point():
 d=accepted('gpu_steady_history');passed=[]
 for i,s in enumerate(d['states']):
  if stop_gate(d['intervals'][:i],s,0,0,policy()):passed.append(s['step'])
 assert passed and min(passed)==d['first_full_steady_step']
