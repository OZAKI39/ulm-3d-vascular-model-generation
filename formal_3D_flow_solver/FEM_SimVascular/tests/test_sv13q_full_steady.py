from sv13q_support import *
def test_at_most_one_full_run_and_complete_acceptance():
 w=read('winner');plans=[p for p in (C/'runplans').glob('REAL_*.json')]
 assert len(plans)==int(w['full_required'])
 if not w['full_required']:assert read('full_steady_decision')['status']=='NOT_REQUIRED';return
 d=read('winner_steady_candidate');assert d['status']=='PASS' and d['safe_stop'] and d['normal_exit']
 assert d['scientific_equivalence']=='DEFERRED' and not d['production_changed']
 assert d['reload']['velocity_finite'] and d['reload']['pressure_finite']
 assert max(d['measurement']['epsilon_Q'],d['measurement']['epsilon_mass'])<=policy()['production_policy']['mass_limit']
 assert d['reuse']['ILU_rebuild_count']==d['profile']['events']['MatLUFactorNum']['count']
