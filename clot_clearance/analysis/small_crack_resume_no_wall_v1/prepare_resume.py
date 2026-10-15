from common import *
import shutil,copy,difflib,ast
assert not O.exists(), 'Never overwrite an archived run'
oldmanifest=json.loads((PARENT/'PAYLOAD_MANIFEST.json').read_text())
for name,z in oldmanifest['files'].items():assert sha(R/name)==z['sha256'],name
for p,h in json.loads((PARENT/'PROTECTED_BASELINES.json').read_text()).items():assert sha(p)==h,p
cp=json.loads((PARENT/'full/LATEST_NATIVE_CHECKPOINT.json').read_text())
assert cp['step']==9600 and abs(cp['time_s']-.0024)<1e-15
for k in ['state','restart']:assert sha(cp[k])==cp[k+'_sha256']
s=np.load(cp['state']);assert s['step']==9600 and np.array_equal(s['x'],np.load(PARENT/'full/last_accepted.npz')['x'])
trial=json.loads((PARENT/'full/TRIAL_ACCEPT_REJECT_AUDIT.json').read_text());assert trial['rejected']==0 and trial['base_intervals']==9600
O.mkdir(parents=True);F.mkdir(parents=True);(O/'full').mkdir();(O/'full/rejected').mkdir();(O/'full/checkpoints').mkdir()
shutil.copytree(PARENT/'provenance',O/'provenance')
for name in ['GEOMETRY_AND_INITIAL_CUT_AUDIT.json','MINIMAL_IMPLEMENTATION_CHECK.json','MIXED_DOF_LEDGER_CHECK.json','PROTECTED_BASELINES.json','PREFIX_50US_COMPARISON.json']:
 shutil.copy2(PARENT/name,O/name)
shutil.copy2(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md',O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md')
prefix={}
for name in ['native_frames.bin','native_frames.csv','damage_updates.bin.gz','contact_pairs.bin.gz','events.jsonl','trials.jsonl','KEYFRAME_INDEX.json','STREAM_LAYOUT.json']:
 shutil.copy2(PARENT/'full'/name,O/'full'/name);prefix['full/'+name]=dict(bytes=(PARENT/'full'/name).stat().st_size,sha256=sha(PARENT/'full'/name))
shutil.copytree(PARENT/'full/keyframes',O/'full/keyframes')
shutil.copy2(PARENT/'NATIVE_DYNAMIC_HISTORY.csv',O/'NATIVE_DYNAMIC_HISTORY.csv');prefix['NATIVE_DYNAMIC_HISTORY.csv']=dict(bytes=(PARENT/'NATIVE_DYNAMIC_HISTORY.csv').stat().st_size,sha256=sha(PARENT/'NATIVE_DYNAMIC_HISTORY.csv'))
freeze=copy.deepcopy(json.loads((PARENT/'SMALL_CRACK_DESIGN_FREEZE.json').read_text()))
freeze['run_id']=RUN;freeze['budget']['main_wall_s']=None;freeze['budget']['RSS_stop_GiB']=None
freeze['budget']['memory_policy']='No application RSS or wall-clock ceiling; explicit user authorization. Original output limit and all numerical gates retained.'
freeze['budget']['estimate']='Only resume 2.4–2.5ms: 400 nominal dt=0.25us steps, no new geometry or parameter case. No fixed wall time ceiling.'
freeze['continuation']=dict(parent=str(PARENT),checkpoint=cp,start_s=cp['time_s'],planned_end_s=END,prior_steps=9600,original_case_continuation=True,original_t0_replayed=False,user_instruction='不设置 1800 秒墙钟预算，继续进行运行')
runtime=['common.py','native.py','checks.py','spatial_model.py','mechanics.py','journal.py','integrator.py','run_protocol.py','launch.py']
freeze['driver_sha256']={n:sha(A/n) for n in runtime}
for name in runtime:ast.parse((A/name).read_text())
dump(O/'SMALL_CRACK_DESIGN_FREEZE.json',freeze);(O/'SMALL_CRACK_DESIGN_FREEZE.sha256').write_text(sha(O/'SMALL_CRACK_DESIGN_FREEZE.json')+'  SMALL_CRACK_DESIGN_FREEZE.json\n')
dump(O/'PREFIX_COPY_MANIFEST.json',prefix)
dump(O/'CONTINUATION_AUTHORIZATION.json',dict(user_instruction=freeze['continuation']['user_instruction'],wall_limit_s=None,RSS_limit_GiB=None,no_replacement_time_or_memory_cap=True,output_limit_GiB=2,physical_window_s=[cp['time_s'],END],source_checkpoint=cp,source_payload_sha256=sha(PARENT/'PAYLOAD_MANIFEST.json'),source_payload_files_verified=len(oldmanifest['files']),all_physics_and_numerical_gates_unchanged=True))
identity=json.loads((PARENT/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text());identity.update(run_id=RUN,parent_archive=str(PARENT),no_CPP_build_this_run=True,reused_library=str(LIB),library_sha256=sha(LIB),continuation_from_native_restart=True);dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',identity)
(O/'driver_diff').mkdir()
for n in runtime:
 diff=''.join(difflib.unified_diff((PARENT_A/n).read_text().splitlines(True),(A/n).read_text().splitlines(True),fromfile='parent/'+n,tofile='resume/'+n))
 if diff:(O/'driver_diff'/(n+'.diff')).write_text(diff)
print('FROZEN_NATIVE_RESUME',RUN,sha(O/'SMALL_CRACK_DESIGN_FREEZE.json'))
