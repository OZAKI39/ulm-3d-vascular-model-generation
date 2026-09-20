from sv13q_support import *
def test_actual_gpu_helper_recovery_probe():
 d=read('remote/adaptive_recovery_probe');assert d['status']=='PASS' and d['run']['exit_code']==0
 assert d['synthetic_fault_injection'] and d['not_CFD'] and d['not_performance_evidence']
 assert d['production_header_sha256']==digest(ROOT/'external/sv13q/svMultiPhysics-reuse/Code/Source/solver/sv13q_reuse.h')
 log=(ROOT/'logs/sv1_3q/remote'/Path(d['run']['log']).name).read_text()
 assert 'old_reason=-5' in log and 'original_rhs_restored=1' in log and 'recovery=STALE_ILU_RECOVERED' in log
 assert 'recovery=FRESH_RETRY_FAILED' in log and 'no_third_attempt=1' in log
 assert float(re.search(r'solution_error_inf=([^ ]+)',log)[1])<1e-9
 assert 'synthetic_failure' not in (ROOT/'external/sv13q/svMultiPhysics-reuse/Code/Source/solver/sv13q_reuse.h').read_text()
def test_real_recoveries_are_never_hidden():
 for d in cases():
  h=read('remote/'+d['name']+'_execution')['history']
  assert d['observed_KSP_reasons']==h['all_petsc_attempt_reasons']
  if d['reuse']:
   assert d['reuse']['recovery_count']==raw(d).count('recovery=STALE_ILU_RECOVERED')
