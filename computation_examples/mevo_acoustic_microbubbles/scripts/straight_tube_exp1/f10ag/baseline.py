"""Reproduce frozen AF before admitting any new scientific dynamics."""
import time,shutil,subprocess,os,importlib.metadata as md
import numpy as np,h5py
from .common import *
def main():
 start=time.perf_counter();BASE.mkdir(parents=True,exist_ok=True)
 protected=read(AF/'PREREGISTRATION.json')['protected_sha256']
 for folder in ['scripts/straight_tube_exp1/f10af','scripts/straight_tube_exp1/f10aa','src/mevo_mb/fast_sim','src/mevo_mb/acoustics']:
  for p in (ROOT/folder).rglob('*.py'):protected[str(p.relative_to(ROOT))]=sha(p)
 for p in [AF/'INITIAL_STATES.h5',AF/'ACTIONS.json',AF/'SCIENTIFIC_RESULTS.json',AF/'FINAL_GATE.json']:
  protected[str(p.relative_to(ROOT))]=sha(p)
 for p in (AF/'fields').rglob('*.h5'):protected[str(p.relative_to(ROOT))]=sha(p)
 write(BASE/'PREREGISTRATION.json',dict(protocol=read(CONFIG),protected_sha256=protected,baseline_commit='4eb54462b91cfc9d81dc08d8edd6803175724523',identity=identity()))
 for name in ['INITIAL_STATES.h5','ACTIONS.json','SOURCE_ORDER.json']:shutil.copy2(AF/name,BASE/name)
 (BASE/'fields').symlink_to(AF/'fields',target_is_directory=True)
 seal=read(AF/'FINAL_EXECUTION_SEAL.json');matched={p:sha(ROOT/p)==v for p,v in seal['source_sha256'].items()};assert all(matched.values())
 env=dict(cpu=subprocess.check_output(['lscpu'],text=True),load=os.getloadavg(),gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,memory.total,memory.used,utilization.gpu,driver_version','--format=csv'],text=True),versions={n:md.version(n) for n in ['numpy','scipy','warp-lang','h5py']},history_source_match=matched,identity=identity())
 from threadpoolctl import threadpool_info
 env['blas']=threadpool_info();write(BASE/'ENVIRONMENT.json',env)
 run,slow,audit,paired=bind();names=['long50_X_PAIR','pair_r0_ROT_PLUS','pair_r0_ROT_MINUS','pair_contact','pair_compact'];plan=[]
 for name in names:
  old=read(AF/'hf'/name/'SUMMARY.json');t=dict(old['task'],id='baseline_'+name);plan.append((old,t))
 write(BASE/'BASELINE_PLAN.json',dict(tasks=[t for _,t in plan],identity=identity()))
 checks=[]
 for old,t in plan:
  r=run.perform(t);check=dict(id=t['id'],status=r['status'],old_status=old['status'],completed_s=r.get('completed_s'),old_completed_s=old.get('completed_s'),task_wall_s=r['task_wall_s'])
  with h5py.File(AF/old['trajectory_path']) as h,h5py.File(BASE/r['trajectory_path']) as z:
   assert h['positions_m'].shape==z['positions_m'].shape
   maximum=0.
   for j in range(0,len(h['time_s']),500):maximum=max(maximum,float(abs(h['positions_m'][j:j+500]-z['positions_m'][j:j+500]).max()))
   check['max_abs_position_difference_m']=maximum;check['PASS']=maximum==0 and r['status']==old['status'] and r['completed_s']==old['completed_s']
  check['audit']=audit.one(r);assert check['PASS'];checks.append(check)
  if name in ['long50_X_PAIR','pair_contact','pair_compact']:
   slow.integrate(t);check['paired']=paired.compare(t)
  write(BASE/'BASELINE_RESULTS.json',dict(rows=checks,PASS=all(z['PASS'] for z in checks),wall_s=time.perf_counter()-start,identity=identity()))
 speed=checks[0]['task_wall_s']/.05
 write(BASE/'BUDGET_FREEZE.json',dict(total_scientific_wall_s=6000,per_run_wall_s=600,workers=3,baseline_s_per_physical_s=speed,baseline_wall_s=time.perf_counter()-start,proposed_ladder_s=[.05,.075,.1,.15,.25],max_LF_s=.25,HF_admission='Predicted cost below 540 s; stop original gates or 600 s; refinement separately bounded',identity=identity()))
 print('BASELINE_ALL_PASS',time.perf_counter()-start,flush=True)
if __name__=='__main__':main()
