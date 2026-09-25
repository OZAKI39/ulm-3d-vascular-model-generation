from pathlib import Path
import json,hashlib,shutil,difflib,subprocess
R=Path('/workspace/microbubble_lammps/results/lammps2026_bubble_interaction_20260916_091244');W=Path('/workspace/lammps_migration/new_20260916_091244')
assert json.loads((R/'MIGRATION_EXECUTION_STATE.json').read_text())['state']=='FAIL'
failed=R/'failed_attempts';failed.mkdir();src=R/'migration/coupling';dst=failed/'coupling_before_explicit_finalize';src.rename(dst);src.mkdir()
for line in (dst/'FROZEN_INPUT_SHA256SUMS').read_text().splitlines():
 h,n=line.split('  ',1);f=dst/n;assert hashlib.sha256(f.read_bytes()).hexdigest()==h;g=src/n;g.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,g)
for f in ['MIGRATION_EXECUTION_STATE.json','MIGRATION_REGRESSION_STARTED.json']:
 shutil.copy2(R/f,failed/f)
for stage,name in [('coupling','src/coupled_lammps_main.cpp'),('passive','src/passive_main.cpp')]:
 root=R/'migration'/stage;f=root/name;before=f.read_text();after=before.replace('#include "lammps.h"','#include "lammps.h"\n#include "library.h"').replace('delete lmp;MPI_Finalize();','delete lmp;lammps_kokkos_finalize();lammps_python_finalize();lammps_plugin_finalize();MPI_Finalize();').replace('delete l;}catch','delete l;lammps_kokkos_finalize();lammps_python_finalize();lammps_plugin_finalize();}catch')
 assert before!=after;f.write_text(after)
 patch={'stage':stage,'scope':'EXPLICIT_LIFECYCLE_API_COMPATIBILITY_ONLY','math_changed':False,'reference':'frozen target src/main.cpp finalize() sequence after deleting LAMMPS and before MPI_Finalize','before_sha256':hashlib.sha256(before.encode()).hexdigest(),'after_sha256':hashlib.sha256(after.encode()).hexdigest(),'diff':''.join(difflib.unified_diff(before.splitlines(True),after.splitlines(True),fromfile='before/'+name,tofile='after/'+name))}
 (R/'provenance'/(stage+'_LIFECYCLE_API_PATCH.json')).write_text(json.dumps(patch,indent=2)+'\n')
 # Only an API registration/lifecycle entrypoint changes. No scientific input or gate changes.
 if stage=='passive':shutil.copy2(root/'FROZEN_INPUT_SHA256SUMS',failed/'PASSIVE_PRE_API_PATCH_SHA256SUMS')
 lines=[hashlib.sha256(q.read_bytes()).hexdigest()+'  '+str(q.relative_to(root)) for q in sorted(root.rglob('*')) if q.is_file() and q.name!='FROZEN_INPUT_SHA256SUMS']
 (root/'FROZEN_INPUT_SHA256SUMS').write_text('\n'.join(lines)+'\n')
shutil.copy2(W/'upstream/lammps/src/main.cpp',R/'provenance/TARGET_UPSTREAM_MAIN.cpp')
(R/'provenance/EXPLICIT_API_REQUALIFICATION.json').write_text(json.dumps({'reason':'Kokkos finalization missing in legacy external embedding entrypoints; source confirms explicit API cleanup required','first_attempt_preserved':str(dst),'first_attempt_exit_code':134,'first_attempt_completed_steps':1000,'automatic_same_binary_retry':False,'requalification':'Run coupling suite once with corrected adapter entrypoint, then first passive suite','authorization':'User sections 15-18 authorize API compatibility patches and regression validation; no physics changes','source_correction_gate':'FAIL_UNCHANGED','phase_B':'PROHIBITED','phase_C':'PROHIBITED'},indent=2)+'\n')
script='''from pathlib import Path
import subprocess,json,os,traceback
R=Path(__R__);W=Path(__W__);env=dict(os.environ,HWLOC_COMPONENTS='-gl',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0',PYTHONDONTWRITEBYTECODE='1')
py='/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019/venv/bin/python'
def state(s,**k):(R/'MIGRATION_EXECUTION_STATE_V2.json').write_text(json.dumps(dict(state=s,**k),indent=2)+'\\n')
try:
 for stage in ['coupling','passive']:
  state('BUILD_API_LIFECYCLE_PATCH',stage=stage)
  with (R/'provenance'/(stage+'_build_v2.log')).open('w') as f:p=subprocess.run(['cmake','--build',str(W/stage/'build'),'--parallel','4'],stdout=f,stderr=subprocess.STDOUT,env=env)
  assert p.returncode==0
  state('RUN_REGRESSION',stage=stage)
  with (R/'provenance'/(stage+'_regression_v2.log')).open('w') as f:p=subprocess.run([py,str(R/'scripts/regression_driver.py'),stage],stdout=f,stderr=subprocess.STDOUT,env=env,timeout=1800)
  assert p.returncode==0,stage
 state('REGRESSIONS_PASS_SOURCE_GATE_FAIL',promotion_allowed=False)
except Exception as e:state('FAIL',error=repr(e),traceback=traceback.format_exc());raise
'''.replace('__R__',repr(str(R))).replace('__W__',repr(str(W)))
(R/'scripts/migration_api_requalification.py').write_text(script)
config=W/'supervisord.conf'
with config.open('a') as f:f.write('\n[program:api_requalification]\ncommand=/usr/bin/python3 '+str(R/'scripts/migration_api_requalification.py')+'\nautostart=true\nautorestart=false\nstartsecs=0\nstdout_logfile='+str(W/'requalification.stdout')+'\nstderr_logfile='+str(W/'requalification.stderr')+'\n')
for op in ['reread','update','status']:print(subprocess.check_output(['supervisorctl','-c',str(config),op],text=True))
