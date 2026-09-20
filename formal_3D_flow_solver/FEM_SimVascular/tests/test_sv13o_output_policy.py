from sv13o_support import *
def test_actual_VTU_cadence_and_forced_final():
 for p in R.glob('*_acceptance.json'):
  d=json.loads(p.read_text())
  if d['status']=='PASS':output_policy_gate(d['VTU_steps'],d['stop_step'],d['start_step'])
 policy=json.loads((C/'policy.json').read_text());assert policy['VTU_cadence']==policy['restart_cadence']==10
def test_per_step_output_and_missing_final_are_rejected():
 with pytest.raises(GateError):output_policy_gate(list(range(1,12)),11)
 with pytest.raises(GateError):output_policy_gate([10],11)
 assert output_policy_gate([10,11],11)

