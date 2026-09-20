from sv13q_support import *
def test_r5_observed_policy_and_health():
 check_candidate('R5')
 for d in cases():
  if d['candidate']=='R5' and d['status']=='PASS':
   validate_reuse(raw(d),d['profile'],read('remote/'+d['name']+'_execution')['history'],'R5')
