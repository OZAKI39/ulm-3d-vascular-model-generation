"""Bounded winner continuation and discrete-frequency QD; no GP between frequencies."""
import time,signal,numpy as np
from scipy.stats import qmc
from ribs.emitters import BayesianOptimizationEmitter
from ribs.schedulers import BayesianOptimizationScheduler
from .common import *
from .run import task,batch
from .search import archive,add,state,quality,eligible,measured,action,make_cma,alarm

def campaign(method,label,count,frequency=2500000.,prior=None):
 start=time.perf_counter();ar=archive(method=='CMA_MAE');passive=archive();connected=archive();records=[];xs_all=[];overhead=0.;skips=[]
 sob=qmc.Sobol(4,scramble=True,seed=301818)
 if prior:
  seed_x=np.array(prior['solutions']);seed_r=prior['rows'];sob.fast_forward(len(seed_r));add(ar,seed_x,seed_r);add(passive,seed_x,seed_r)
 else:
  seed_x=sob.random_base2(3);seed_r=batch([task(label+'_init_'+str(i),amps=action(x)[0],phase=action(x)[1],frequency=frequency) for i,x in enumerate(seed_x)],4);records+=seed_r;xs_all+=seed_x.tolist();add(ar,seed_x,seed_r);add(passive,seed_x,seed_r);add(connected,seed_x,seed_r,True)
 sched=None
 if method=='BOP':
  mask=np.array([measured(r) for r in seed_r]);good=seed_x[mask]
  if not len(good):return dict(rows=records,solutions=xs_all,stop='NO_VALID_SEED',wall_s=time.perf_counter()-start)
  em=BayesianOptimizationEmitter(ar,lower_bounds=np.zeros(4),upper_bounds=np.ones(4),initial_solutions=good,batch_size=1,search_nrestarts=1,upscale_schedule=None,seed=301819);sched=BayesianOptimizationScheduler(ar,[em]);sched.ask();sched.tell([quality(r) for r,m in zip(seed_r,mask) if m],[r['descriptor'] for r,m in zip(seed_r,mask) if m])
 elif method=='CMA_MAE':sched=make_cma(ar,ar.best_elite['solution'] if len(ar) else None,301819)
 stop=None
 while len(records)<count:
  if time.perf_counter()-start>600:stop='TIME_BUDGET_STOP';break
  begin=time.perf_counter()
  try:
   signal.signal(signal.SIGALRM,alarm);signal.alarm(45);xs=sob.random(4) if method=='SOBOL' else sched.ask()
  except Exception as e:stop=repr(e);break
  finally:signal.alarm(0)
  overhead+=time.perf_counter()-begin
  rs=batch([task(label+'_'+str(len(records)+i),amps=action(x)[0],phase=action(x)[1],frequency=frequency) for i,x in enumerate(xs)],4);records+=rs;xs_all+=np.asarray(xs).tolist();add(passive,xs,rs);add(connected,xs,rs,True)
  begin=time.perf_counter();mask=np.array([measured(r) for r in rs])
  if method=='BOP':
   if mask.any():sched._cur_solutions=np.array(xs)[mask];sched._num_emitted=[int(mask.sum())];sched.tell([quality(r) for r,m in zip(rs,mask) if m],[r['descriptor'] for r,m in zip(rs,mask) if m])
   else:sched._last_called='tell';skips.append(len(records))
  elif method=='CMA_MAE':
   if mask.all():sched.tell([quality(r) for r in rs],[r['descriptor'] for r in rs])
   else:add(ar,xs,rs);sched=make_cma(ar,ar.best_elite['solution'] if len(ar) else None,301819+len(records));skips.append(len(records))
  overhead+=time.perf_counter()-begin
  write(BASE/'extensions'/label/'PROGRESS.json',dict(method=method,frequency=frequency,rows=records,solutions=xs_all,archive=state(passive),connected_archive=state(connected),wall_s=time.perf_counter()-start,optimizer_wall_s=overhead))
 result=dict(method=method,frequency=frequency,rows=records,solutions=xs_all,archive=state(passive),connected_archive=state(connected),wall_s=time.perf_counter()-start,optimizer_wall_s=overhead,invalid_update_skips=skips,stop=stop,continuation='Sobol continuation; BOP refits historical observations; CMA official reseed from known best, no reconstructed private covariance' if prior else 'independent fixed-frequency optimizer',identity=identity())
 write(BASE/'extensions'/label/'RESULT.json',result);return result

def main():
 import sys
 mode=sys.argv[1] if len(sys.argv)>1 else 'winner';comparison=read(BASE/'COMPARISON_V2.json');winner=comparison['winner']
 if winner is None:write(BASE/'extensions/NOT_RUN.json',dict(reason='No complete fair optimizer comparison'));return
 if mode=='winner':campaign(winner,'winner',32,prior=read(BASE/'search_v2'/winner/'RESULT.json'))
 elif mode=='frequencies':
  admit=read(AD/'SHORT_ADMISSION.json');plans=[]
  for f in [2400000.,2600000.]:
   if any(r['task'].get('frequency')==f and r.get('numerical_admitted') for r in admit):plans.append(f)
  write(BASE/'FREQUENCY_SEARCH_PREREGISTRATION.json',dict(frequencies=plans,evaluations_each=16,method=winner,independent_models=True,unchanged_archive=True,identity=identity()))
  for f in plans:campaign(winner,'frequency_'+str(int(f)),16,frequency=f)
if __name__=='__main__':main()
