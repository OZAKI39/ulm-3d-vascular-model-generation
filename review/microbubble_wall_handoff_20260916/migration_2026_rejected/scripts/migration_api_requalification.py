from pathlib import Path
import subprocess,json,os,traceback
R=Path('/workspace/microbubble_lammps/results/lammps2026_bubble_interaction_20260916_091244');W=Path('/workspace/lammps_migration/new_20260916_091244');env=dict(os.environ,HWLOC_COMPONENTS='-gl',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0',PYTHONDONTWRITEBYTECODE='1')
py='/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019/venv/bin/python'
def state(s,**k):(R/'MIGRATION_EXECUTION_STATE_V2.json').write_text(json.dumps(dict(state=s,**k),indent=2)+'\n')
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
