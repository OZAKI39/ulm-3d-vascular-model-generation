from sv13q_support import *
def test_frozen_monitor_and_recorded_native_stop():
 if not read('winner')['full_required']:assert read('full_steady_decision')['status']=='NOT_REQUIRED';return
 s=read('stop_latency');assert s['status']=='PASS'
 assert s['solver_ack_timestamp']>=s['stop_request_timestamp']
 assert s['actual_stop_step']>=s['first_qualifying_saved_step']
 h=read('steady_history');assert h['status']=='PASS' and h['stop_eligible']
 assert h['policy_sha256']==policy()['production_policy_sha256']
 assert s['additional_full_runs']==0
 assert read('remote/native_solver_stopped')['active_solver_pids']==[]
