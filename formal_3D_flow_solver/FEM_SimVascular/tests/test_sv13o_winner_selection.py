from sv13o_support import *
def test_health_and_wall_time_control_selection():
 assert selection_action(100,dict(status='PASS',wall_time_s=89))['action']=='ACCEPT'
 assert selection_action(100,dict(status='PASS',wall_time_s=90))['action']=='ACCEPT'
 assert selection_action(100,dict(status='PASS',wall_time_s=97))['action']=='KEEP_INCUMBENT'
 assert selection_action(100,dict(status='PASS',wall_time_s=93))['action']=='CONFIRM_ONCE'
 assert selection_action(100,dict(status='PASS',wall_time_s=95))['action']=='CONFIRM_ONCE'
 assert selection_action(100,dict(status='FAIL',wall_time_s=1))['action']=='REJECT_UNHEALTHY'
def test_only_ambiguous_timing_may_repeat():
 first=dict(status='PASS',wall_time_s=93)
 assert confirm_selection(100,first,dict(status='PASS',wall_time_s=92))['accepted']
 assert not confirm_selection(100,first,dict(status='FAIL',wall_time_s=10))['accepted']
 with pytest.raises(GateError):confirm_selection(100,dict(status='PASS',wall_time_s=50),first)
def test_actual_winner_and_repeat_budget():
 w=accepted('winner');assert w['candidate'] in 'ABCD'
 assert w['fixed_window_reduction']==pytest.approx(1-w['selection_wall_time_s']/w['baseline_A_wall_s'])
 for p in (C/'runplans').glob('*_CONFIRM.json'):
  candidate=json.loads(p.read_text())['candidate'];d=next(x for x in w['decisions'] if x['candidate']==candidate)
  assert d['action']=='CONFIRM_ONCE' and .05<=d['reduction']<.10
 assert not any('median' in p.name.lower() for p in (C/'runplans').iterdir())
 profile=accepted('profile_runs');assert profile['runs']<=2 and profile['excluded_from_winner_selection']
