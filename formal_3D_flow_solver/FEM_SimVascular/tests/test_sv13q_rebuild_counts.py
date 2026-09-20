from sv13q_support import *
def test_counters_agree_with_actual_petsc_factorizations():
 for d in cases():
  if d['status']!='PASS':continue
  r=d['reuse'];events=d['profile']['events']
  assert r['ILU_rebuild_count']==events['MatLUFactorNum']['count']
  assert r['KSP_attempts']==events['KSPSolve']['count']
  assert r['ILU_reuse_count']+r['ILU_rebuild_count']==r['KSP_attempts']
  assert d['total_attempt_iterations']==sum(a['outcome']['iterations'] for a in r['trace'])
  assert d['profile']['nested_not_additive']
