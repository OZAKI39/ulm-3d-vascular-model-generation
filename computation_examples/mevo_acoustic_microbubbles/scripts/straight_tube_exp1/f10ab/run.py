"""Reuse F10AA Campaign.run verbatim; only supply frozen shapes and schedules."""
import argparse,time,traceback,concurrent.futures,multiprocessing,importlib.metadata,os
import numpy as np,h5py
from scipy.linalg import cho_factor,cho_solve
from mevo_mb.fast_sim.accelerated_water import WaterAccelerated
from mevo_mb.fast_sim.accelerated_contact import resistance
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from mevo_mb.control.scheduled_pzt import Schedule
from scripts.straight_tube_exp1.f10aa.trajectory import Campaign as OriginalCampaign
from .metrics import deformation,measure
from .common import *

class AuditEngine(WaterAccelerated):
 def step(self,dt):
  value=super().step(dt)
  self.gate_records.append(value['gates'].copy())
  return value

class Campaign(OriginalCampaign):
 def __init__(self,out,state_id,action_record,remaining_s):
  verify();self.started=time.perf_counter();self.out=out;out.mkdir(parents=True,exist_ok=False)
  self.cfg=dict(max_trajectories=1,per_run_budget_s=240,remote_budget_s=remaining_s)
  self.backend='B1';self.state=state_id;self.p=policy();self.c=context();self.ident=identity();self.ident['execution_commit']=self.ident['source_commit'];self.records=[]
  self.bases={'coarse':self.c['basis'],'fine':CachedWaterBasis(OLD/'cache/fine.h5')}
  with h5py.File(BASE/'INITIAL_STATES.h5') as h:x=h[state_id][:]
  self.states={len(x):x};self.schedules={action_record['name']:Schedule.from_record(action_record['schedule'])}
 def engine(self,n,cache='coarse'):
  e=AuditEngine(self.states[n],self.c['physics'],self.c['response'],self.bases[cache],self.c['cavity'],policy=self.p)
  e.gate_records=[];self.last_engine=e;return e

def perform(task):
 begin=time.perf_counter();out=BASE/task['stage']/task['run_id']
 action=dict(name=task['action'],schedule=task.get('schedule',read(BASE/'ACTIONS.json').get(task['action'])))
 try:
  c=Campaign(out,task['state_id'],action,task.get('remaining_s',240))
  n=next(iter(c.states));r=c.run(n,task['action'],task['duration_s'],task.get('variant','DT1'))
  folder=out/'runs'/r['name'];write(folder/'SOLVER_GATES.json',c.last_engine.gate_records)
  r.update(run_id=task['run_id'],state_id=task['state_id'],stage=task['stage'],trajectory_path=str((folder/'trajectory.h5').relative_to(BASE)),task=task)
  with h5py.File(folder/'trajectory.h5') as h:
   xx=h['positions_m'];r['initial_shape']=measure(xx[0],xx[0]);r['final_shape']=measure(xx[0],xx[-1])
  # Do not embed positions in tiny index; retained completely in original HDF5.
  for key in ['initial_shape','final_shape']:r[key].pop('aligned_positions_m',None)
  write(out/'SUMMARY.json',r)
 except Exception as e:
  r=dict(run_id=task['run_id'],state_id=task['state_id'],stage=task['stage'],status='EXCEPTION',reason=repr(e),traceback=traceback.format_exc(),task=task)
  write(out/'EXCEPTION.json',r)
 r['total_task_wall_s']=time.perf_counter()-begin
 return clean(r)

def static():
 verify();out=BASE/'static';out.mkdir(exist_ok=False);begin=time.perf_counter();c=context();p=policy();actions=read(BASE/'ACTIONS.json');records=[]
 with h5py.File(BASE/'INITIAL_STATES.h5') as h:
  for state_id in [k for k in h if k.endswith('r0')]:
   x=h[state_id][:];op=resistance(x,R);chol=cho_factor(op['A']);rows={}
   for name,rec in actions.items():
    try:
     sim=WaterAccelerated(x,c['physics'],c['response'],c['basis'],c['cavity'],policy=p);sim.control(Schedule.from_record(rec).segments[0].event())
     f=sim.fields(x);v=cho_solve(chol,f['total_force_N'].ravel()/sim.zeta).reshape(x.shape)
     primary_v=cho_solve(chol,f['primary_force_N'].ravel()/sim.zeta).reshape(x.shape)
     secondary_v=cho_solve(chol,f['secondary_force_N'].ravel()/sim.zeta).reshape(x.shape)
     step=sim.step(1e-6);u=c['basis'].vector(rec['segments'][0]['amplitude_Pa_equivalent'],rec['segments'][0]['phase_rad'],rec['segments'][0]['enabled'])
     pp,gg=c['basis'].query_basis(x);grad=np.einsum('i,inj->nj',u,gg)
     rates=deformation(x,x+1e-6*v);r=dict(state_id=state_id,action=name,status='PASS',scope='FIXED_STATE_FORCE_DIAGNOSTIC',
       source_l1_Pa_equivalent=float(abs(u).sum()),source_l2_squared_Pa2=float(np.vdot(u,u).real),u=u,
       pressure_min_Pa=float(abs(f['pressure_Pa']).min()),pressure_max_Pa=float(abs(f['pressure_Pa']).max()),gradient_Pa_per_m=grad,
       center_force_N=f['total_force_N'].mean(0),center_primary_force_N=f['primary_force_N'].mean(0),center_secondary_force_N=f['secondary_force_N'].mean(0),
       force_relative_rms_N=float(np.sqrt(np.mean(np.sum((f['total_force_N']-f['total_force_N'].mean(0))**2,axis=1)))),
       velocity_free_m_s=v,primary_free_velocity_m_s=primary_v,secondary_free_velocity_m_s=secondary_v,
       velocity_constrained_m_s=step['result']['velocity_m_s'],contact_force_N=step['result']['contact_force_N'],
       free_deviatoric_rate_s1=rates['deviatoric_stretch_3d']/1e-6,free_signed_xy_rate_s1=rates['directional_xy_stretch_3d']/1e-6,
       single_step_gates=step['gates'],resistance_min_eigenvalue=float(np.linalg.eigvalsh(op['A'])[0]),**f)
     # Store F1/F2, relative force, resistance, exact state and gradients separately.
     with h5py.File(out/f'{state_id}__{name}.h5','w') as dest:
      for k,val in {**f,'positions_m':x,'velocity_free_m_s':v,'resistance_matrix':op['A'],'gradient_Pa_per_m':grad,
        'relative_force_N':f['total_force_N']-f['total_force_N'].mean(0),'velocity_constrained_m_s':step['result']['velocity_m_s']}.items():dest[k]=val
      dest.attrs['scope']='FIXED_STATE_FORCE_DIAGNOSTIC'
    except Exception as e:r=dict(state_id=state_id,action=name,status='FAIL',reason=repr(e),traceback=traceback.format_exc(),scope='FIXED_STATE_FORCE_DIAGNOSTIC')
    rows[name]=r;records.append(r);print('STATIC',state_id,name,r['status'],flush=True)
   write(out/(state_id+'.json'),rows)
 write(out/'SUMMARY.json',dict(records=records,wall_s=time.perf_counter()-begin,identity=identity()))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('stage');args=ap.parse_args()
 if args.stage=='media':
  from .figures import main as figures
  from .videos import main as videos
  figures();videos();return
 if args.stage in ['target_evaluate','extension_evaluate']:
  from .target import evaluate,extension_evaluate
  (evaluate if args.stage=='target_evaluate' else extension_evaluate)();return
 if args.stage=='static':static();return
 if args.stage=='sensitivity_analysis':
  from .sensitivity import analyze
  analyze();return
 if args.stage=='confirm_evaluate':
  from .evaluate import confirmation
  confirmation();return
 if args.stage=='verify_sensitivity':
  from .sensitivity import verify_checks
  verify_checks();return
 if args.stage=='evaluate':
  from .evaluate import numerical_admission
  numerical_admission();return
 if args.stage.startswith('analyze_'):
  from .analyze import analyze_stage
  analyze_stage(args.stage[len('analyze_'):]);return
 import fcntl,datetime
 lockfile=(BASE/'evidence/physical_stage.lock');lockfile.parent.mkdir(exist_ok=True)
 lock=lockfile.open('a');fcntl.flock(lock,fcntl.LOCK_EX)
 tasks=read(BASE/'plans'/f'{args.stage}.json');cfg=read(CONFIG);stage_dir=BASE/args.stage;stage_dir.mkdir(exist_ok=False)
 write(stage_dir/'START.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),identity=identity(),serialized_physical_stages=True))
 start=time.perf_counter();results=[]
 consumed=sum(read(p).get('wall_s',0) for p in BASE.glob('*/COMPLETE.json'))
 if (BASE/'static/SUMMARY.json').exists():consumed+=read(BASE/'static/SUMMARY.json')['wall_s']
 deadline=start+max(0,min(tasks['stage_budget_s'],cfg['total_server_budget_s']-consumed))
 with concurrent.futures.ProcessPoolExecutor(max_workers=tasks.get('workers',2),mp_context=multiprocessing.get_context('spawn')) as ex:
  pending={};queue=list(tasks['tasks'])
  while queue or pending:
   while queue and len(pending)<tasks.get('workers',2):
    item=queue.pop(0);remaining=deadline-time.perf_counter()
    if remaining<30:
     results.append(dict(run_id=item['run_id'],status='NOT_RUN',reason='STAGE_WALL_BUDGET'));continue
    item['remaining_s']=min(remaining,240);pending[ex.submit(perform,item)]=item
   if not pending:break
   done,_=concurrent.futures.wait(pending,timeout=1,return_when=concurrent.futures.FIRST_COMPLETED)
   for future in done:
    item=pending.pop(future)
    try:r=future.result()
    except Exception as e:r=dict(run_id=item['run_id'],status='EXCEPTION',reason=repr(e),task=item)
    results.append(r);write(stage_dir/'INDEX.json',results)
    print('FINISHED',item['run_id'],r['status'],'elapsed',round(time.perf_counter()-start,2),flush=True)
 write(stage_dir/'COMPLETE.json',dict(wall_s=time.perf_counter()-start,identity=identity(),results=len(results),plan_sha256=sha(BASE/'plans'/f'{args.stage}.json')))
if __name__=='__main__':main()
