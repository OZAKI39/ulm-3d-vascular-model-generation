from particle_3d.particle65_cases import run_case

def test_crossing_trial_is_split():
 d=run_case('wall',.032,.032)
 assert d['attempts'][0]['minimum_g_nf_m']<0 and not d['attempts'][0]['accepted']
 assert d['summary']['rejected_trials']>0 and max(x['depth'] for x in d['ledger'])>0
 assert d['summary']['handoff_time_s']<.032 and d['summary']['final_time_s']==.032
 assert d['summary']['tangential_displacement_m']>6.39e-8
