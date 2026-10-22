"""Original F10AA trajectory integration, recording and contact gates, unchanged."""
import time,traceback,concurrent.futures,multiprocessing
import numpy as np,h5py
from scripts.straight_tube_exp1.f10aa.trajectory import Campaign as Original
from scripts.straight_tube_exp1.f10ab.run import AuditEngine
from scripts.straight_tube_exp1.f10ab.metrics import measure,geometry,series
from scripts.straight_tube_exp1.f10ad_multifrequency.interface import context,make_schedule
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from .common import *
class Campaign(Original):
 def __init__(self,out,task):
  verify();self.started=time.perf_counter();self.out=out;out.mkdir(parents=True,exist_ok=False)
  self.cfg=dict(max_trajectories=1,per_run_budget_s=120,remote_budget_s=task.get('remaining_s',120))
  self.backend='B1';self.state=task['state'];self.p=policy();self.c=context(task.get('frequency',2500000.),task.get('shell',1.))
  self.ident=identity();self.ident['execution_commit']=self.ident['source_commit'];self.records=[]
  f=task.get('frequency',2500000.);fine=OLD/'cache/fine.h5' if f==2500000. else AD/str(round(f))/'fine.h5'
  self.bases={'coarse':self.c['basis'],'fine':CachedWaterBasis(fine)}
  self.input_library=BASE/task.get('initial_library','INITIAL_STATES.h5')
  if not self.input_library.resolve().is_relative_to(BASE.resolve()):raise ValueError('Initial library must be inside this campaign')
  with h5py.File(self.input_library) as h:x=h[self.state][:]
  if task.get('offset_m'):x=x+np.asarray(task['offset_m'])
  self.states={len(x):x};self.schedules={'CONTROL':make_schedule(task['amplitudes'],task['phases'],task['duration'],f)}
 def engine(self,n,cache='coarse'):
  e=AuditEngine(self.states[n],self.c['physics'],self.c['response'],self.bases[cache],self.c['cavity'],policy=self.p)
  e.gate_records=[];self.last_engine=e;return e

def perform(task):
 start=time.perf_counter();out=BASE/'runs'/task['id']
 try:
  c=Campaign(out,task);r=c.run(next(iter(c.states)),'CONTROL',task['duration'],task.get('variant','DT1'))
  folder=out/'runs'/r['name'];path=folder/'trajectory.h5'
  # Retain original INPUT; add explicit authoritative child protocol and physical inputs.
  write(folder/'F10AC_INPUT.json',dict(task=task,identity=identity(),protocol_sha256=sha(BASE/'PREREGISTRATION.json'),initial_states_sha256=sha(c.input_library),field=c.bases[r['cache']].metadata,response=c.c['response'].metadata))
  write(folder/'SOLVER_GATES.json',c.last_engine.gate_records)
  with h5py.File(path) as h:
   x=h['positions_m'][:];t=h['time_s'][:];m=series(t,x,10);initial=measure(x[0],x[0]);final=measure(x[0],x[-1])
  for z in [initial,final]:z.pop('aligned_positions_m',None)
  np.savez_compressed(out/'METRICS.npz',**{k:v for k,v in m.items() if k!='aligned_positions_m'})
  descriptor=[float(np.log(final['Rg_3d_m']/initial['Rg_3d_m'])),float(np.log(final['aspect_3d']/initial['aspect_3d']))]
  status='VALID_COMPLETE_TRAJECTORY' if r['status']=='PASS' else ('TIME_BUDGET_STOP' if any('BUDGET' in s for s in r['reasons']) else ('OUT_OF_MODEL_DOMAIN' if any('APECSS' in s or 'cache' in s or 'Outside' in s for s in r['reasons']) else 'NUMERICAL_FAILURE'))
  r.update(id=task['id'],task=task,trajectory_path=str(path.relative_to(BASE)),initial=initial,final=final,descriptor=descriptor if status=='VALID_COMPLETE_TRAJECTORY' else None,scientific_status=status,connected_at_end=final['proximity_components']==1,connected_sampled_path=bool(np.all(m['proximity_components']==1)))
 except Exception as e:r=dict(id=task['id'],task=task,status='FAIL',scientific_status='OUT_OF_MODEL_DOMAIN' if isinstance(e,ValueError) else 'NUMERICAL_FAILURE',reason=repr(e),traceback=traceback.format_exc(),descriptor=None)
 r['task_wall_s']=time.perf_counter()-start;write(out/'SUMMARY.json',r);return clean(r)
def batch(tasks,workers=2):
 with concurrent.futures.ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as ex:return list(ex.map(perform,tasks))
def task(id,state='S2_r0',amps=(0,0,1036.2016581041792,1036.2016581041792),phase=(0,0,0,np.pi),duration=.002,**kw):
 return dict(id=id,state=state,amplitudes=list(amps),phases=list(phase),duration=duration,**kw)
