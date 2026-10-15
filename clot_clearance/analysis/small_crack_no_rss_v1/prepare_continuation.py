from common import *
import shutil,difflib,ast
from datetime import datetime,timezone
assert not O.exists(),'Never overwrite an archived run'
prior=json.loads((PRIOR_SMALL/'PAYLOAD_MANIFEST.json').read_text())
for name,item in prior['files'].items():assert sha(R/name)==item['sha256'],name
for folder in [O,F,O/'provenance',O/'full',O/'source_diffs']:folder.mkdir(parents=True,exist_ok=True)
shutil.copy2(R/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.md',O/'ACTIVE_DEVELOPMENT_ENTRYPOINTS.before.md')
for p in (PRIOR_SMALL/'provenance').iterdir():shutil.copy2(p,O/'provenance'/p.name)
for name in ['MINIMAL_IMPLEMENTATION_CHECK.json','MIXED_DOF_LEDGER_CHECK.json','GEOMETRY_AND_INITIAL_CUT_AUDIT.json','PROTECTED_BASELINES.json']:shutil.copy2(PRIOR_SMALL/name,O/name)
old=json.loads((PRIOR_SMALL/'SMALL_CRACK_DESIGN_FREEZE.json').read_text());new=json.loads(json.dumps(old));new.update(run_id=RUN,frozen_utc=datetime.now(timezone.utc).isoformat(),source_run=str(PRIOR_SMALL),source_freeze_sha256=sha(PRIOR_SMALL/'SMALL_CRACK_DESIGN_FREEZE.json'),user_instruction='不设置 4 GiB 预算，继续进行运行',physical_protocol_changed=False,CPP_recompiled=False,library_reused=str(LIB))
new['budget']['RSS_stop_GiB']=None
new['budget']['memory_policy']='No application-level RSS ceiling; measure actual peak RSS. User explicitly withdrew 4 GiB limit. Original wall/output and all numerical gates unchanged.'
new['driver_sha256']={name:sha(A/name) for name in old['driver_sha256']}
for name,h in new['reference_sha256'].items():assert sha(O/'provenance'/name)==h,name
assert sha(LIB)==old['library_sha256']
for name in old['driver_sha256']:
 src=R/'analysis/small_crack_v1'/name;dst=A/name
 if sha(src)!=sha(dst):(O/'source_diffs'/(name+'.diff')).write_text(''.join(difflib.unified_diff(src.read_text().splitlines(True),dst.read_text().splitlines(True),fromfile='prior/'+name,tofile='no_rss_limit/'+name)))
for name in ['native.py','checks.py','spatial_model.py','mechanics.py','journal.py','integrator.py','freeze.py','launch.py']:assert sha(A/name)==sha(R/'analysis/small_crack_v1'/name),name
for p in A.glob('*.py'):ast.parse(p.read_text())
dump(O/'SMALL_CRACK_DESIGN_FREEZE.json',new);(O/'SMALL_CRACK_DESIGN_FREEZE.sha256').write_text(sha(O/'SMALL_CRACK_DESIGN_FREEZE.json')+'  SMALL_CRACK_DESIGN_FREEZE.json\n')
identity=json.loads((PRIOR_SMALL/'SOURCE_AND_LIBRARY_IDENTITY.json').read_text());identity.update(continuation_run=RUN,CPP_changed_in_this_continuation=False,CPP_recompiled_in_this_continuation=False,reused_library=True,previous_run=str(PRIOR_SMALL));dump(O/'SOURCE_AND_LIBRARY_IDENTITY.json',identity)
dump(O/'CONTINUATION_AUTHORIZATION.json',dict(user_instruction=new['user_instruction'],memory_stop_removed=True,replacement_memory_cap=None,all_other_budget_and_numerical_gates_unchanged=True,geometry_material_loading_contact_CPP_unchanged=True,start='Fresh original precut t=0; last50us FP64 is not a complete native restart. No trajectory prefix splicing.',minimal_checks_reused_from_previous_run=True,minimal_check_sha256=sha(PRIOR_SMALL/'MINIMAL_IMPLEMENTATION_CHECK.json'),previous_157_file_manifest_verified=True,source_manifest_sha256=sha(PRIOR_SMALL/'PAYLOAD_MANIFEST.json'),additional_cases_authorized=1))
(O/'REPRODUCE.txt').write_text(f'Source /workspace/y_clot/env.sh. Driver {A}/launch.py full. Reused library {LIB}, SHA {sha(LIB)}. Same frozen geometry, material, loading, event integrator, numerical gates; only RSS stop disabled by explicit user instruction. No new C++ build or P1 physical rerun. Native runtime requires LD_LIBRARY_PATH={LIB.parent}:existing and PYTHONPATH={NATIVE_PREV}/source/python:existing, one thread. Existing archive guards prohibit overwriting. Full protocol is 2.5ms; do not assume completion without full/DYNAMIC_AUDIT.json. Prior t0-50us is a separate archived run, never spliced into this one.\n')
print('CONTINUATION_READY',RUN,sha(O/'SMALL_CRACK_DESIGN_FREEZE.json'))
