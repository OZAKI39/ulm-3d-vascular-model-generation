"""Finite scientific stages under one Supervisor label, never restart old batches."""
import time,sys,subprocess
from .common import *
def main():
 assert read(BASE/'COMPARISON_V2.json')['fair_complete']
 phases=[('tests','f10ac_qd.check',[],120),('winner','f10ac_qd.extend',['winner'],650),('frequency_dynamics','f10ad_multifrequency.dynamics',[],1200),('frequency_search','f10ac_qd.extend',['frequencies'],1250),('independent_elites','f10ac_qd.validate',[],900)]
 start=time.perf_counter();records=[]
 consumed=read(BASE/'COMPARISON.json')['wall_s']+read(BASE/'COMPARISON_V2.json')['wall_s']+read(BASE/'PILOT.json')['wall_s']+read(AD/'BUILD_SUMMARY.json')['wall_s']
 write(BASE/'WORKFLOW_START.json',dict(phases=phases,prior_recorded_s=consumed,global_scientific_budget_s=read(CONFIG)['total_scientific_wall_budget_s'],identity=identity()))
 for name,module,args,budget in phases:
  remaining=read(CONFIG)['total_scientific_wall_budget_s']-consumed-(time.perf_counter()-start)-600
  if remaining<30:raise RuntimeError('Registered global budget reached before '+name)
  print('STAGE_BEGIN',name,flush=True);begin=time.perf_counter();p=subprocess.run([sys.executable,'-m','scripts.straight_tube_exp1.'+module,*args],timeout=min(budget,remaining))
  records.append(dict(name=name,returncode=p.returncode,wall_s=time.perf_counter()-begin,command=[sys.executable,'-m','scripts.straight_tube_exp1.'+module,*args]));write(BASE/'WORKFLOW_PROGRESS.json',records)
  if p.returncode:raise RuntimeError('Scientific stage failed; completed outputs preserved: '+name)
 write(BASE/'WORKFLOW_COMPLETE.json',dict(stages=records,wall_s=time.perf_counter()-start,identity=identity()))
if __name__=='__main__':main()
