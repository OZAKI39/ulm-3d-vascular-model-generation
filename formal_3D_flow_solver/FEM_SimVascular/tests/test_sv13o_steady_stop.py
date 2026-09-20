from sv13o_support import *
def test_original_steady_rule_first_eligible_point():
 h=accepted('steady_history');c=accepted('optimized_steady_candidate');policy=json.loads((ROOT/'configs/sv1_3/policy.json').read_text())
 eligible=[s['step'] for s in h['states'] if stop_gate([i for i in h['intervals'] if i['step']<=s['step']],s,0,0,policy)]
 assert eligible and eligible[0]==c['first_full_steady_step']
 assert c['required_consecutive_intervals']==policy['steady_last_intervals'] and h['formal_check_cadence']==policy['save_interval_steps']
 r=read('steady_stop_request');assert r['kind']=='STEADY' and r['detail']['first_full_steady_step']==eligible[0]
 assert c['safe_stop'] and c['normal_exit'] and c['first_full_steady_step']<=c['stop_step']<=c['first_full_steady_step']+1
def test_no_extra_full_CFD():
 plans=[json.loads(p.read_text()) for p in (C/'runplans').glob('*.json')]
 assert [p['name'] for p in plans if p['mode']=='full']==['REAL_VASCULAR_GPU_PERF']
 assert all(p['MPI_ranks']==p['GPUs']==p['OMP_NUM_THREADS']==1 for p in plans)
 assert all('CPU' not in p['name'] for p in plans)

