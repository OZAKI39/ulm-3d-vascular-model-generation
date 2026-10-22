"""Reassess only a completed CFD run after documented monitor classification fix."""
from pathlib import Path
import json,sys,hashlib,time,argparse
V=Path(__file__).resolve().parents[1];sys.path.insert(0,str(V/'scripts/production'))
from log_acceptance import classify
p=argparse.ArgumentParser();p.add_argument('--case',type=Path,required=True);a=p.parse_args();c=a.case.resolve();assert c.is_relative_to(V)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def put(p,obj):p.write_text(json.dumps(obj,indent=2)+'\n')
e=c/'reports/execution.json';old=json.loads(e.read_text());policy=json.loads((c/'policy.json').read_text());assert old['status']=='FAIL' and old['exit_code']==0 and old['stop']['consecutive']>=policy['steady_consecutive_intervals']
assert all(sha(c/k)==h for k,h in json.loads((c/'input_hashes.json').read_text()).items())
manifest=json.loads((c/'frozen_flow/manifest.json').read_text());assert manifest['input_configuration_sha256']==sha(c/'run/solver.xml')
for k,v in manifest['files'].items():assert sha(c/'frozen_flow'/k)==v['sha256']
quality=json.loads((c/'reports/flow_quality.json').read_text());assert quality['accepted_final_and_log_checks'] and quality['last_nonlinear_relative']<1e-10
assert old['solver_log_sha256']==sha(c/'run/solver.log')
last=old['states'][-1]
assert max(last[k] for k in ['velocity_relative_change','pressure_relative_change','wss_area_L2_relative_change'])<policy['steady_change_limit']
assert last['epsilon_mass']<=policy['mass_limit'] and last['epsilon_Q']<=1e-6 and last['wall_max_speed_m_s']<=1e-12
accept=classify((c/'run/solver.log').read_text(),policy['dt_s']);assert accept['accepted'] and len(accept['verified_recovered_attempts'])==1
original=c/'reports/execution_original_monitor.json';assert not original.exists();original.write_bytes(e.read_bytes())
put(c/'reports/linear_attempt_acceptance.json',accept)
receipt=dict(reassessed_unix=time.time(),original_status=old['status'],original_execution_sha256=sha(original),raw_log_sha256=sha(c/'run/solver.log'),classifier_sha256=sha(V/'scripts/log_acceptance.py'),reassessment_script_sha256=sha(Path(__file__)),CFD_rerun=False,inputs_outputs_logs_unchanged=True,reason='Original monitor rejected any DIVERGED marker, including verified discarded stale-preconditioner attempt followed by a successful fresh solve. Independent field, accepted correction and restored-RHS checks pass. Failed attempt remains recorded.',independent_flow_quality_sha256=sha(c/'reports/flow_quality.json'))
put(c/'reports/monitor_reassessment.json',receipt);old.update(status='PASS',linear_attempt_acceptance=accept,monitor_reassessment=receipt);put(e,old);print(json.dumps(receipt,indent=2))
