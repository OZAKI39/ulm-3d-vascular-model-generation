from sv13p_support import *
def test_lifecycle_audited_before_reuse():
 d=accepted('ksp_lifecycle_audit');assert d['KSP_PC_persist_across_linear_solves']
 assert d['target_policy']=='REUSE_WITHIN_TIMESTEP'
 for name in ('P1_SMOKE','P1_WINDOW'):
  d=accepted(name+'_acceptance');trace=d['reuse']['trace']
  assert len({r['ksp'] for r in trace})==len({r['pc'] for r in trace})==1
