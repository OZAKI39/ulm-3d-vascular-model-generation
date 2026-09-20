from sv13q_support import *
def test_adaptive_rule_frozen_and_actual_trace_valid():
 p=policy();assert p['adaptive_threshold']==1.5 and p['max_reuse_age']==5
 check_candidate('RA')
 for d in cases():
  if d['candidate']=='RA' and d['status']=='PASS':validate_reuse(raw(d),d['profile'],read('remote/'+d['name']+'_execution')['history'],'RA')
def test_adaptive_timeline_retains_reference_age_and_reason():
 d=case('RA','SMOKE');rows=attempts(raw(d));assert rows
 for r in rows:assert all(k in r for k in ('age','ref','rebuild_reason','reuse','step','outcome'))
