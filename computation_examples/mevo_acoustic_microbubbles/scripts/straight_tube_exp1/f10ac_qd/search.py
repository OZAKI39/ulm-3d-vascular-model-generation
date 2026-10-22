"""Thin ask/evaluate/tell integration of SciPy Sobol and official pyribs optimizers."""
import time,signal,warnings,json
import numpy as np
from scipy.stats import qmc
from ribs.archives import GridArchive
from ribs.emitters import BayesianOptimizationEmitter,EvolutionStrategyEmitter
from ribs.schedulers import BayesianOptimizationScheduler,Scheduler
from .common import *
from .run import task,batch
SEED=301818

class BoundedGridArchive(GridArchive):
 """Only archive admission is filtered; official emitters receive true observations.

 pyribs clips out-of-bounds measures into boundary cells by default. This adapter
 explicitly excludes those entries and returns no-add metadata, without altering
 the official GP, acquisition, CMA update, or threshold update.
 """
 def add(self,solution,objective,measures,**fields):
  x=np.asarray(solution);y=np.asarray(objective);m=np.asarray(measures);inside=np.all((m>=self.lower_bounds)&(m<self.upper_bounds),axis=1)
  result=dict(status=np.zeros(len(x),np.int32),value=np.zeros(len(x),float))
  if inside.any():
   info=super().add(x[inside],y[inside],m[inside],**{k:np.asarray(v)[inside] for k,v in fields.items()})
   for key,values in info.items():
    if key not in result:result[key]=np.zeros(len(x),dtype=values.dtype)
    result[key][inside]=values
  return result

def archive(cma=False):
 cfg=read(CONFIG)
 return BoundedGridArchive(solution_dim=4,dims=cfg['grid'],ranges=cfg['ranges'],seed=SEED,**(dict(learning_rate=.01,threshold_min=0) if cma else {}))
def measured(r):return r['scientific_status']=='VALID_COMPLETE_TRAJECTORY' and r['descriptor'] is not None
def action(x):
 a,phase,y,yp=np.asarray(x,float)
 if np.any(np.asarray(x)<0) or np.any(np.asarray(x)>1):raise ValueError('Action outside registered bounds')
 return [1100*y,1100*y,200+900*a,200+900*a],[0,2*np.pi*yp,0,np.deg2rad(150+90*phase)]
def eligible(r):
 if r['scientific_status']!='VALID_COMPLETE_TRAJECTORY':return False
 d=np.array(r['descriptor']);limits=np.array(read(CONFIG)['ranges'])
 return bool(np.all(d>=limits[:,0]) and np.all(d<limits[:,1]))
def quality(r):return 100/(1+np.linalg.norm(r['displacement_m'])/1e-6)
def add(ar,xs,rows,connected=False):
 mask=np.array([eligible(r) and (not connected or r['connected_at_end']) for r in rows])
 if mask.any():ar.add(np.asarray(xs)[mask],np.array([quality(r) for r in rows])[mask],np.array([r['descriptor'] if r['descriptor'] is not None else [np.nan,np.nan] for r in rows])[mask])
def state(ar):return dict(coverage=float(ar.stats.coverage),qd_score=float(ar.stats.qd_score),cells=len(ar),elites=ar.data())
def alarm(signum,frame):raise TimeoutError('OPTIMIZER_ASK_45_SECOND_LIMIT')
def make_cma(ar,x0=None,seed=SEED):
 e=EvolutionStrategyEmitter(ar,x0=np.full(4,.5) if x0 is None else x0,sigma0=.22,ranker='imp',selection_rule='mu',restart_rule='basic',lower_bounds=np.zeros(4),upper_bounds=np.ones(4),batch_size=4,seed=seed)
 return Scheduler(ar,[e])
def main():
 verify();assert all(r['PASS'] for r in read(BASE/'REGRESSION.json'))
 cfg=read(CONFIG);start=time.perf_counter();workers=read(BASE/'PILOT.json')['workers_selected'];workers=min(workers,4)
 common=qmc.Sobol(4,scramble=True,seed=SEED).random_base2(3)
 write(BASE/'SEARCH_V2_PREREGISTRATION.json',dict(cfg=cfg,initial_sobol=common,workers=workers,ask_budget_s=45,invalid_policy='No fake descriptor or reward. BOP learns all real completed descriptors, including out-of-archive bounds; archive admission alone excludes out-of-bounds. CMA physically failed batch has no adaptation and is restarted via official emitter; all valid trajectories still enter passive archive.',seed=SEED,pre_search=True,identity=identity()))
 t=time.perf_counter();initial=batch([task('common_v2_'+str(i),amps=action(x)[0],phase=action(x)[1]) for i,x in enumerate(common)],workers);common_wall=time.perf_counter()-t
 write(BASE/'COMMON_V2.json',dict(solutions=common,rows=initial,wall_s=common_wall))
 summaries={}
 for method in ['SOBOL','BOP','CMA_MAE']:
  begun=time.perf_counter();ar=archive(method=='CMA_MAE');passive=archive();connected=archive();records=[];optimizer_wall=0.;updates_skipped=[];suggestions=[]
  add(ar,common,initial);add(passive,common,initial);add(connected,common,initial,True)
  records.extend(initial);suggestions.extend(common.tolist());history=[dict(evaluations=8,elapsed_s=common_wall,**{k:v for k,v in state(passive).items() if k!='elites'})]
  if method=='BOP':
   emitter=BayesianOptimizationEmitter(ar,lower_bounds=np.zeros(4),upper_bounds=np.ones(4),initial_solutions=common,batch_size=1,search_nrestarts=1,upscale_schedule=None,seed=SEED)
   sched=BayesianOptimizationScheduler(ar,[emitter]);asked=sched.ask();mask=np.array([measured(r) for r in initial]);sched._cur_solutions=np.asarray(asked)[mask];sched._num_emitted=[int(mask.sum())]
   if not mask.any():raise RuntimeError('No valid initial samples for BOP')
   tm=time.perf_counter();sched.tell([quality(r) for r,m in zip(initial,mask) if m],[r['descriptor'] for r,m in zip(initial,mask) if m]);optimizer_wall+=time.perf_counter()-tm
  elif method=='CMA_MAE':sched=make_cma(ar)
  sobol=qmc.Sobol(4,scramble=True,seed=SEED);sobol.fast_forward(8)
  stop=None
  while len(records)<cfg['comparison_evaluations_per_method']:
   if time.perf_counter()-start>cfg['optimizer_comparison_s']:stop='TIME_BUDGET_STOP';break
   tm=time.perf_counter()
   try:
    signal.signal(signal.SIGALRM,alarm);signal.alarm(45)
    xs=sobol.random(4) if method=='SOBOL' else sched.ask()
   except Exception as e:stop=repr(e);break
   finally:signal.alarm(0)
   optimizer_wall+=time.perf_counter()-tm
   tasks=[task(method.lower()+'_v2_'+str(len(records)+i),amps=action(x)[0],phase=action(x)[1]) for i,x in enumerate(xs)]
   rows=batch(tasks,workers);records.extend(rows);suggestions.extend(np.asarray(xs).tolist());add(passive,xs,rows);add(connected,xs,rows,True)
   tm=time.perf_counter();mask=np.array([measured(r) for r in rows])
   if method=='BOP':
    if mask.any():
     sched._cur_solutions=np.asarray(xs)[mask];sched._num_emitted=[int(mask.sum())];sched.tell([quality(r) for r,m in zip(rows,mask) if m],[r['descriptor'] for r,m in zip(rows,mask) if m])
    else:
     # No invented observation: official GP remains unchanged; move ask/tell handshake only.
     sched._last_called='tell';updates_skipped.append(dict(at=len(records),reason='no valid response; unchanged official GP'))
   elif method=='CMA_MAE':
    if mask.all():sched.tell([quality(r) for r in rows],[r['descriptor'] for r in rows])
    else:
     add(ar,xs,rows);updates_skipped.append(dict(at=len(records),reason='failed generation; official restart without fabricated descriptors'))
     center=ar.best_elite['solution'] if len(ar) else np.full(4,.5);sched=make_cma(ar,center,SEED+len(records))
   else:add(ar,xs,rows)
   optimizer_wall+=time.perf_counter()-tm
   history.append(dict(evaluations=len(records),elapsed_s=common_wall+time.perf_counter()-begun,**{k:v for k,v in state(passive).items() if k!='elites'}))
   write(BASE/'search_v2'/method/'PROGRESS.json',dict(rows=records,solutions=suggestions,history=history,optimizer_wall_s=optimizer_wall,updates_skipped=updates_skipped,archive=state(passive),connected_archive=state(connected)))
   print(method,'evaluations',len(records),'cells',len(passive),'optimizer_s',round(optimizer_wall,2),flush=True)
  result=dict(method=method,rows=records,solutions=suggestions,history=history,optimizer_wall_s=optimizer_wall,updates_skipped=updates_skipped,archive=state(passive),connected_archive=state(connected),common_cached=8,true_new_calls=len(records)-8,total_charged_wall_s=common_wall+time.perf_counter()-begun,stop=stop,complete=len(records)==cfg['comparison_evaluations_per_method'],identity=identity())
  write(BASE/'search_v2'/method/'RESULT.json',result);summaries[method]={k:v for k,v in result.items() if k not in ['rows','solutions']}
 candidates=[m for m,r in summaries.items() if r['complete']]
 winner=max(candidates,key=lambda m:(summaries[m]['archive']['coverage']/summaries[m]['total_charged_wall_s'],summaries[m]['archive']['qd_score']/summaries[m]['total_charged_wall_s'])) if len(candidates)==3 else None
 write(BASE/'COMPARISON_V2.json',dict(methods=summaries,winner=winner,fair_complete=len(candidates)==3,wall_s=time.perf_counter()-start,common_wall_s=common_wall,scope='single optimizer seed; fixed S2 initial, 2ms, 2.5MHz; no universal optimizer ranking'))
if __name__=='__main__':main()
