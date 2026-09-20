from sv13q_support import *
def test_only_late_top_two_get_early_runs():
 rank=read('late_ranking');assert len(rank['top_two'])<=2
 actual=[d['candidate'] for d in cases() if d['mode']=='early'];assert set(actual)==set(rank['top_two'])
 assert rank['top_two']==[r['candidate'] for r in rank['ranking'][:2]]
 for d in cases():
  if d['mode']=='early':
   assert d['start_step']==10 and d['initial_checkpoint_sha256']==policy()['early_checkpoint_sha256']
   if d['status']=='PASS':assert d['stop_step']==20 and d['steps']==10
 assert read('early_checkpoint')['status']=='PASS'
