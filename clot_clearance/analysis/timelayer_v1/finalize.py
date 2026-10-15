"""Finalize evidence identities only. No native simulation is launched."""
from common import *
import time,subprocess
begin=time.perf_counter()
before=json.loads((O/'HISTORICAL_PROTECTION_BEFORE.json').read_text())
mismatch=[]
for p,rec in before.items():
 q=Path(p)
 if not q.is_file() or q.stat().st_size!=rec['bytes'] or sha(q)!=rec['sha256']: mismatch.append(p)
dump(O/'HISTORICAL_PROTECTION_AFTER.json',{'status':'PASS' if not mismatch else 'FAIL','checked_files':len(before),'checked_bytes':sum(v['bytes'] for v in before.values()),'mismatches':mismatch,'before_manifest_sha256':sha(O/'HISTORICAL_PROTECTION_BEFORE.json'),'scope':'All pre-registered historical bugfix/V3/macro result files and listed protected native libraries/source/env; excluded build/source trees except explicitly listed identity files','check_wall_s':time.perf_counter()-begin})
assert not mismatch
ident=json.loads((O/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text());assert all(sha(p)==v for p,v in ident['files_sha256'].items())
ident['status']='PASS_ON_DISK_AND_ACTUAL_PROCESS_MAPS';ident.pop('new_native_compute_count',None)
ident['formal_completed_dynamic_cases']=2;ident['native_instance_attempts']=3;ident['setup_only_failed_attempts']=1;ident['setup_failure_accepted_dynamic_steps']=0
ident['actual_runtime']={}
for name in ['dt_0p5us','dt_0p25us']:
 d=json.loads((O/name/'LOADED_LIBRARY.json').read_text());libs={p:v for p,v in d['actual_maps'].items() if 'liblammps.so' in p};assert len(libs)==1 and next(iter(libs.values()))==EXPECTED
 ident['actual_runtime'][name]=d
ident['end_of_round_recheck']='All listed source/library/env SHA unchanged';dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',ident)
sp=json.loads((O/'SMALL_SPECIMEN_IDENTITY.json').read_text());sp['status']='PASS_HISTORICAL_STATIC_IDENTITY_WITH_AUTHORIZED_DYNAMIC_ADAPTATION';sp['dynamic_adaptation_authorization']='ACCEPTED: same nodes/bonds/material/15 um amplitude; freeze continuous time history before two dt runs';sp['authorization_record']='DYNAMIC_ADAPTATION_AUTHORIZATION.json';sp['design_freeze']='TIMESTEP_DESIGN_FREEZE.json';dump(O/'SMALL_SPECIMEN_IDENTITY.json',sp)
code="""import json,xmlrpc.client,supervisor.xmlrpc
r=xmlrpc.client.ServerProxy('http://127.0.0.1',transport=supervisor.xmlrpc.SupervisorTransport(None,None,'unix:///var/run/supervisor.sock'))
print(json.dumps([r.supervisor.getProcessInfo(n) for n in NAMES]))
""".replace('NAMES',repr(['clot_timelayer_'+C['run'],'clot_timelayer_'+C['run']+'_formal']))
p=subprocess.run(['/usr/bin/python3','-c',code],capture_output=True,text=True,check=True);managed=json.loads(p.stdout);assert managed[1]['statename']=='EXITED' and managed[1]['exitstatus']==0
dump(O/'MANAGED_COMPLETION.json',managed)
freeze=json.loads((O/'TIMESTEP_DESIGN_FREEZE.json').read_text());fcheck={n:sha(R/'analysis/timelayer_v1'/n)==v for n,v in freeze['diagnostic_source_sha256'].items()};assert all(fcheck.values())
assert sha(O/'TIMESTEP_DESIGN_FREEZE.json')==json.loads((O/'CAMPAIGN_LAUNCHED.json').read_text())['freeze_sha256']
dump(O/'FINAL_EXECUTION_FREEZE_CHECK.json',{'frozen_execution_code':fcheck,'design_freeze_sha256':sha(O/'TIMESTEP_DESIGN_FREEZE.json'),'freeze_precedes_formal_launch':freeze['implementation_revision']['epoch']<json.loads((O/'CAMPAIGN_LAUNCHED.json').read_text())['epoch'],'no_additional_native_runs':True})
print(json.dumps({'protected_files':len(before),'protected_bytes':sum(v['bytes'] for v in before.values()),'identity':ident['status'],'formal_process_state':managed[1]['statename'],'exit_status':managed[1]['exitstatus'],'wall_s':time.perf_counter()-begin}))
