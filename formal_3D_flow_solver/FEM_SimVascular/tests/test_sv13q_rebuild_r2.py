from sv13q_support import *
def test_r2_observed_policy_and_health():
 check_candidate('R2')
 for d in cases():
  if d['candidate']=='R2' and d['status']=='PASS':
   validate_reuse(raw(d),d['profile'],read('remote/'+d['name']+'_execution')['history'],'R2')
