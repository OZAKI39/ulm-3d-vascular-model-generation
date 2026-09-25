"""Independent migration candidate qualification. Never promotes or removes an installation."""
from pathlib import Path
import os,sys,json,subprocess,time,traceback,hashlib
R=Path(__file__).resolve().parents[1];W=Path('/workspace/lammps_migration/new_20260916_091244')
PYTHON='/workspace/microbubble_lammps/work/palabos_lammps_coupling_v0_20260915_224019/venv/bin/python'
ENV=dict(os.environ,HWLOC_COMPONENTS='-gl',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0',PYTHONDONTWRITEBYTECODE='1')
def save(p,d):p.write_text(json.dumps(d,indent=2)+'\n')
def state(s,**k):save(R/'MIGRATION_EXECUTION_STATE.json',dict(state=s,**k));print(s,k,flush=True)
def command(args,log,timeout=900):
 with log.open('w') as f:
  f.write('COMMAND_JSON='+json.dumps(args)+'\n');f.flush();p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,env=ENV,timeout=timeout)
 assert p.returncode==0,str(log)
def main():
 with (R/'MIGRATION_REGRESSION_STARTED.json').open('x') as f:json.dump({'time':time.time(),'no_automatic_solver_retries':True},f)
 state('WAIT_BUILD')
 while True:
  s=json.loads((R/'provenance/BUILD_STATUS.json').read_text())
  if s['state']=='FAIL':raise RuntimeError('candidate build failed')
  if s['state']=='PASS':break
  time.sleep(10)
 for stage in ['coupling','passive']:
  state('BUILD_CUSTOM',stage=stage);b=W/stage/'build'
  command(['cmake','-S',str(R/'migration'/stage),'-B',str(b),'-G','Ninja','-DCMAKE_BUILD_TYPE=Release','-DENGINE_ROOT='+str(W)],R/'provenance'/(stage+'_configure.log'))
  command(['cmake','--build',str(b),'--parallel','4'],R/'provenance'/(stage+'_build.log'))
 # Engine's old contract remains an immutable regression definition; this is actual new build provenance.
 p=R/'migration/engine/build_provenance';p.mkdir();build=json.loads((R/'provenance/LAMMPS_2026_BUILD_PROVENANCE.json').read_text())
 save(p/'LAMMPS_BUILD_PROVENANCE.json',{'status':'PASS','builds':build,'actual_target':json.loads((R/'provenance/TARGET_RELEASE_FROZEN.json').read_text()),'role':'2026 migration candidate; old contract source denotes historical regression provenance'})
 state('ENGINE_REGRESSION')
 command([PYTHON,str(R/'migration/engine/scripts/run_stage.py'),'--root',str(R/'migration/engine'),'--work',str(W)],R/'provenance/engine_regression.log',1800)
 for stage in ['coupling','passive']:
  state('CUSTOM_REGRESSION',stage=stage)
  command([PYTHON,str(R/'scripts/regression_driver.py'),stage],R/'provenance'/(stage+'_regression.log'),1800)
 state('REGRESSIONS_PASS_SOURCE_GATE_SEPARATE',promotion_performed=False)
if __name__=='__main__':
 try:main()
 except Exception as e:state('FAIL',error=repr(e),traceback=traceback.format_exc());raise
