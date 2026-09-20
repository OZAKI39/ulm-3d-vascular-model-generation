from sv13p_support import *
def test_accepted_cases_have_positive_reasons_and_finite_fields():
 good=[d for d in cases() if d['status']=='PASS'];assert good
 for d in good:
  assert d['linear_failures']==d['nonlinear_failures']==0 and d['normal_exit'] and not d['errors']
  assert d['Mat']=='seqaijcusparse' and d['Vec']=='seqcuda'
  assert all(v>0 for v in d['KSP_reason_values'].values())
  assert d['reload']['velocity_finite'] and d['reload']['pressure_finite']
  for s in d['runtime_semantics']:semantics_gate(s,d['PETSC_OPTIONS'])
def test_relaxed_ksp_rejected():
 d=accepted('P1_SMOKE_acceptance');s=dict(d['runtime_semantics'][0],rtol=1e-6)
 with pytest.raises(GateError):semantics_gate(s,d['PETSC_OPTIONS'])
def test_failed_smokes_remain_rejected():
 for name,iterations in [('P2',200),('P3',400),('P4',200)]:
  d=read(name+'_SMOKE_acceptance');e=read('remote/'+name+'_SMOKE_execution')
  assert d['status']=='FAIL' and d['errors'] and e['monitor_stop']['reason']=='DIVERGED_BREAKDOWN'
  assert any(r['diverged'] and r['iterations']==iterations for r in e['history']['petsc_reasons'])
  assert not (R/(name+'_WINDOW_acceptance.json')).exists()
