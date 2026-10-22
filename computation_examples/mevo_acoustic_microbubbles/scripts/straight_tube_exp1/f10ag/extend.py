"""Optional 250 ms HF capture confirmation, cost-admitted after actual 150 ms."""
import time
from .common import *
from .experiment import batch
from .paired import compare

def main():
 begin=time.perf_counter();verify();assert (BASE/'EXPERIMENT_DONE.json').exists();run,slow,_,_=bind();budget=read(BASE/'BUDGET_FREEZE.json');used=read(BASE/'EXPERIMENT_DONE.json')['wall_s']+budget['baseline_wall_s'];tasks=[];decisions=[]
 for state in ['near','random_0']:
  r=read(BASE/'hf'/('long150_'+state)/'SUMMARY.json');prediction=r['task_wall_s']*.25/.15*1.1;allow=r['status']=='PASS' and prediction<540 and used+prediction<4500
  decisions.append(dict(state=state,predicted_task_wall_s=prediction,admitted=allow))
  if allow:tasks.append(run.task('long250_'+state,state='AG32_'+state,duration=.25));used+=prediction
 write(BASE/'HF250_PLAN.json',dict(tasks=tasks,decisions=decisions,identity=identity()));rows=batch(tasks,kinds=('HF','LF'),workers=2) if tasks else [];pairs=[compare(t) for t in tasks]
 # A serial, same-state timing probe avoids treating shared-resource contention
 # as an intrinsic step cost. Its new cost rule is frozen before the probe.
 probe=run.task('serial_DT_probe',state='AG32_random_0',duration=.0025,variant='DT05')
 write(BASE/'SERIAL_DT_PROBE_PLAN.json',dict(task=probe,admission='Measured serial integration cost ×60×1.15 <540 s, otherwise keep 50 ms reference',identity=identity()))
 measured=run.perform(probe);pred=measured.get('integration_s',600)/.0025*.15*1.15
 prior=read(BASE/'REFINEMENT_RESULTS.json')['rows'];refs=[]
 if not any(r['status']=='PASS' and r['task'].get('variant')=='DT05' for r in prior):
  if measured['status']=='PASS' and pred<540:
   refs=[run.task('long150_random_0_DT05',state='AG32_random_0',duration=.15,variant='DT05')]
  else:
   refs=[run.task('dtcheck50',state='AG32_random_0',duration=.05),run.task('dtcheck50_DT05',state='AG32_random_0',duration=.05,variant='DT05'),run.task('dtcheck50_FINE',state='AG32_random_0',duration=.05,variant='FINE')]
 write(BASE/'DT_FALLBACK_PLAN.json',dict(tasks=refs,serial_probe_prediction_s=pred,identity=identity()));rr=[run.perform(t) for t in refs]
 # Spatial refinement of the continued historical N64 state has its own scope.
 extra=run.task('long75_AF_X_FINE',state='N64_wide_r0',duration=.075,variant='FINE')
 write(BASE/'AF75_FINE_PLAN.json',dict(task=extra,identity=identity()));rr.append(run.perform(extra))
 write(BASE/'HF250_RESULTS.json',dict(rows=rows,paired=pairs,extra_refinements=rr,wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
