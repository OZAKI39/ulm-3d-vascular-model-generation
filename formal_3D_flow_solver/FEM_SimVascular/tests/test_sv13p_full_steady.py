from sv13p_support import *
def test_full_run_is_unique_and_only_for_large_gain():
 w=read('winner');paths=list((ROOT/'configs/sv1_3p/runplans').glob('REAL_VASCULAR_GPU_PC_WINNER.json'))
 assert len(paths)==int(w['full_required'])
 if not w['full_required']:assert read('full_steady_decision')['status']=='NOT_REQUIRED';return
 d=accepted('winner_steady_candidate');h=accepted('steady_history')
 assert d['safe_stop'] and d['normal_exit'] and d['required_consecutive_intervals']==5
 assert d['first_full_steady_step']==h['first_full_steady_step']
 assert d['scientific_equivalence']=='DEFERRED' and not d['production_changed']
def test_actual_stop_latency_is_preserved_and_reported():
 d=read('stop_latency');c=accepted('winner_steady_candidate')
 assert d['first_qualifying_saved_step']==c['first_full_steady_step']==70
 assert d['actual_stop_step']==c['stop_step']==72
 assert d['steps_after_first_qualifying_state']==2 and d['extra_full_runs']==0
 assert d['native_stop_request_unix_s']==read('steady_stop_request')['request_unix_s']
 assert d['final_process_wall_s']==c['wall_time_s']
 assert c['VTU_steps']==[10,20,30,40,50,60,70,72]
