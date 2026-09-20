from sv13o_support import *
def test_all_accepted_solvers_have_positive_reasons_and_finite_fields():
 rows=[json.loads(p.read_text()) for p in R.glob('*_acceptance.json')]
 passed=[d for d in rows if d['status']=='PASS'];assert len(passed)>=4
 for d in passed:
  assert not d['errors'] and d['linear_failures']==d['nonlinear_failures']==0 and d['normal_exit']
  assert d['Mat']=='seqaijcusparse' and d['Vec']=='seqcuda'
  assert d['KSP_reason_values'] and all(v>0 for v in d['KSP_reason_values'].values())
  assert d['reload']['velocity_finite'] and d['reload']['pressure_finite']
  if d['measurement']:assert all(math.isfinite(d['measurement'][k]) for k in ('epsilon_Q','epsilon_mass','Q_in_m3_s','Q_out_total_m3_s'))
  for s in d['runtime_semantics']:semantics_gate(s,d['PETSC_OPTIONS'])
def test_relaxed_tolerance_is_rejected():
 d=accepted('PERF_A_acceptance');s=dict(d['runtime_semantics'][0],rtol=1e-6)
 with pytest.raises(GateError):semantics_gate(s,d['PETSC_OPTIONS'])

