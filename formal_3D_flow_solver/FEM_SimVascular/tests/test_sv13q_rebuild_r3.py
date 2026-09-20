from sv13q_support import *
def test_r3_observed_policy_and_health():
 check_candidate('R3')
 for d in cases():
  if d['candidate']=='R3' and d['status']=='PASS':
   validate_reuse(raw(d),d['profile'],read('remote/'+d['name']+'_execution')['history'],'R3')
