from sv13q_support import *
def test_wall_based_selection_and_full_gate():
 w=read('winner');eligible=w['eligible'];assert w['full_required']==bool(eligible and eligible[0]['improvement']>=.10)
 if eligible:
  assert w['candidate']==min(eligible,key=lambda x:x['observed_total_wall_s'])['candidate']
  assert w['observed_total_wall_s']==w['early_wall_s']+w['late_wall_s']
  assert abs(w['improvement']-(1-w['observed_total_wall_s']/w['baseline_total_wall_s']))<1e-12
 for e in eligible:
  assert case(e['candidate'],'EARLY')['status']=='PASS'
  assert e['early_wall_s']<=1.05*read('baseline_r1')['early']['wall_time_s']
