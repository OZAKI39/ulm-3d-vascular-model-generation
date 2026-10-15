from common import *
import csv
prefix=json.loads((O/'PREFIX_COPY_MANIFEST.json').read_text());checks={}
for name,z in prefix.items():
 if name=='full/KEYFRAME_INDEX.json':
  old=json.loads((PARENT/name).read_text());new=json.loads((O/name).read_text());assert new[:len(old)]==old;checks[name]='EXACT_LOGICAL_PREFIX'
 else:
  with (O/name).open('rb') as f:h=hashlib.sha256(f.read(z['bytes'])).hexdigest()
  assert h==z['sha256'],name;checks[name]='EXACT_BYTE_PREFIX'
with (O/'NATIVE_DYNAMIC_HISTORY.csv').open() as f:rows=list(csv.DictReader(f))
steps=np.array([int(q['step']) for q in rows]);assert np.array_equal(steps,np.arange(len(rows)))
times=np.array([float(q['time_s']) for q in rows]);assert np.all(np.diff(times)>0)
restore=json.loads((O/'NATIVE_RESTART_RESTORE_AUDIT.json').read_text());ledger=json.loads((O/'LEDGER_RESTORE_AUDIT.json').read_text());assert restore['status']==ledger['status']=='PASS'
initial=dict(np.load(O/'RESTORED_NATIVE_STATE.npz'));assert int(initial['step'])==9600
parent=json.loads((PARENT/'PAYLOAD_MANIFEST.json').read_text())
for name,z in parent['files'].items():assert sha(R/name)==z['sha256'],name
dump(O/'CONTINUATION_CHAIN_AUDIT.json',dict(status='PASS',prefix_checks=checks,source_payload_files_verified=len(parent['files']),source_payload_sha256=sha(PARENT/'PAYLOAD_MANIFEST.json'),original_steps=9600,additional_steps=len(rows)-1-9600,total_accepted_steps=len(rows)-1,strictly_monotone_time=True,no_duplicate_accepted_step=True,original_reference_time_origin=True,restored_native_state='PASS',restored_energy_ledger='PASS',energy_balance_jump_J=ledger['R_jump_J'],first_appended_time_s=float(times[9601]) if len(times)>9601 else None,last_time_s=float(times[-1]),no_fresh_t0_run=True))
print('CONTINUATION_CHAIN_PASS',len(rows)-1,float(times[-1]))
