from sv13q_support import *
def test_r1_is_frozen_p_winner_without_rerun():
 b=read('baseline_r1');assert b['no_baseline_rerun'];assert b['late']['wall_time_s']==213.3558020569617
 assert b['late']['statistics']['total_iterations']==6720
 assert b['late']['reuse']['PC_rebuild_count']==10
 assert not any(p.stem.startswith('R1_') for p in (C/'runplans').glob('*.json'))
def test_r1_early_observed_interval_is_from_existing_run():
 e=read('baseline_r1')['early'];assert e['wall_time_s']==e['end_row']['reported_elapsed_s']-e['start_row']['reported_elapsed_s']
 assert e['start_row']['step']==10 and e['end_row']['step']==20
 assert not e['standalone_restart_baseline'] and e['limitation']
