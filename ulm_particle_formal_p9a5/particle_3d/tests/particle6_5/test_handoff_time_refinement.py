from particle_3d.particle65_cases import run_case

def test_crossing_trial_is_split():
 d=run_case('wall',.032,.032)
 # P9-A.1 resolves the first crossing as a certified real-time event. The
 # old rejected midpoint staircase is intentionally no longer required.
 assert d['attempts'][0]['accepted']
 assert d['summary']['all_accepted_above_lower']
 assert d['summary']['rejected_trials']==0 and max(x['depth'] for x in d['ledger'])==0
 assert d['ledger'][0]['handoff_event']['kind']=='WALL'
 assert d['summary']['handoff_time_s']<.032 and d['summary']['final_time_s']==.032
 assert d['summary']['tangential_displacement_m']>6.39e-8
