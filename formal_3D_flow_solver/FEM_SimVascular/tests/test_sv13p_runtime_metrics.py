from sv13p_support import *
def test_actual_window_metrics_are_consistent():
 for d in cases():
  if d['status']!='PASS':continue
  st=d['statistics'];assert st['mean_iterations']==pytest.approx(st['total_iterations']/st['KSP_solves'])
  assert d['wall_time_s']>0 and d['sampled_device_memory_max_MiB']>0
  if d['mode']=='window':assert d['profile']['events']['MatMult']['time_s'] is not None
  assert d['profile']['nested_not_additive']
def test_na_event_counts_retained_and_stages_aggregated():
 p=profile_events('--- Event Stage 0: A\nMatMult   3 1.0 n/a n/a 0\n--- Event Stage 1: B\nMatMult   4 1.0 2.0 1.0 0\n')
 assert p['events']['MatMult']['count']==7 and p['events']['MatMult']['time_s'] is None
def test_p1_cost_reduction_is_observed_not_requested():
 d=accepted('P1_WINDOW_acceptance');e=d['profile']['events']
 assert e['MatLUFactorNum']['count']==10 and e['KSPSolve']['count']==20
 assert e['PCSetUp']['count']==10 and e['PCSetUpOnBlocks']['count']==20
 assert e['MatLUFactorNum']['time_s']<e['PCSetUpOnBlocks']['time_s']
 assert d['statistics']['total_iterations']>read('reference_freeze')['window_baseline']['statistics']['total_iterations']
 assert d['wall_time_s']<.9*read('reference_freeze')['window_baseline']['wall_time_s']
