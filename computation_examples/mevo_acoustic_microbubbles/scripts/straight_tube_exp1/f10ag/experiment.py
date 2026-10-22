"""Bounded full original dynamics; task arrays preregistered before each stage."""
import time,json,concurrent.futures,multiprocessing,shutil
import numpy as np,h5py
from scipy.spatial.distance import pdist
from scripts.straight_tube_exp1.f10af.states import dispersed
from scripts.straight_tube_exp1.f10af.field import response_context
from .common import *
from .paired import compare

def worker(job):
 run,slow,_,_=bind();t,kind=job
 return run.perform(t) if kind=='HF' else slow.integrate(t)
def batch(tasks,kinds=('HF','LF'),workers=3):
 jobs=[(t,k) for t in tasks for k in kinds]
 with concurrent.futures.ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as ex:return list(ex.map(worker,jobs))
def states():
 c=response_context();ref=next(v for v in read(BASE/'STATIC_RESULTS.json')['rows'] if v['frequency_Hz']==2500000 and v['phase_deg']==0);roots=np.array([r['position_m'][0] for r in ref['roots']]);out={};rng=np.random.default_rng(463020)
 for name,d in [('near',15e-6),('mid',45e-6),('far',100e-6)]:
  q=[]
  while len(q)<32:
   i=len(q);v=np.array([roots[i%len(roots)]+(-1 if i%2 else 1)*d,rng.uniform(-350e-6,350e-6),rng.uniform(-60e-6,60e-6)])
   if not q or np.linalg.norm(np.array(q)-v,axis=1).min()>8*R:q.append(v)
  out['AG32_'+name]=np.array(q)
 for j in range(3):out[f'AG32_random_{j}']=dispersed(32,463101+j)
 ang=np.arange(32)*2*np.pi/32;out['AG32_ring']=np.c_[150e-6*np.cos(ang),150e-6*np.sin(ang),np.zeros(32)]
 out['AG32_local']=np.stack(np.meshgrid(np.arange(4),np.arange(4),np.arange(2),indexing='ij'),-1).reshape(-1,3)*(2.4*R);out['AG32_local']-=out['AG32_local'].mean(0)
 rows=[]
 shutil.copy2(BASE/'INITIAL_STATES.h5',BASE/'evidence/BASELINE_INITIAL_STATES.h5')
 with h5py.File(BASE/'INITIAL_STATES.h5','a') as h:
  for name,x in out.items():
   assert name not in h;gap=float(pdist(x).min()-2*R);assert gap>=0 and all(c['cavity'].query(z,R)['valid'] for z in x);assert np.all(abs(x)<[600e-6,500e-6,80e-6]);h[name]=x
   rows.append(dict(id=name,N=len(x),minimum_gap_m=gap,xyz_extent_m=np.ptp(x,axis=0),source='PRESCRIBED_RING_INITIAL_CONDITION' if name.endswith('ring') else ('PRESCRIBED_NEAR_BAND_INITIAL_CONDITION' if name.endswith(('near','mid','far')) else 'Frozen random cloud' if 'random' in name else 'PRESCRIBED_LOCAL_COMPACT_INITIAL_CONDITION')))
 write(BASE/'NEW_STATES.json',dict(rows=rows,sha256=sha(BASE/'INITIAL_STATES.h5'),identity=identity()))
def main():
 start=time.perf_counter();verify();assert read(BASE/'BASELINE_RESULTS.json')['PASS'];run,slow,_,paired=bind()
 if not (BASE/'NEW_STATES.json').exists():states()
 tasks=[run.task('pilot_'+s,state='AG32_'+s,duration=.02) for s in ['near','mid','far','random_0','random_1','random_2']]
 tasks += [run.task('pilot_ring_'+a,action=a,state='AG32_ring',duration=.02) for a in ['ROT_PLUS','ROT_MINUS']]
 tasks += [run.task('pilot_local_'+a,action=a,state='AG32_local',duration=.002) for a in ['ROT_PLUS','ROT_MINUS']]
 tasks += [run.task('pilot_highgradient',action='X_180',state='AG32_random_0',duration=.02)]
 write(BASE/'PILOT_PLAN.json',dict(tasks=tasks,identity=identity()));rows=batch(tasks);checks=[compare(t) for t in tasks];write(BASE/'PILOT_RESULTS.json',dict(rows=rows,paired=checks,wall_s=time.perf_counter()-start,identity=identity()))
 # Only numerically successful complete pairs progress. No pattern-based early rejection.
 admitted={r['id'] for r in checks if r['PASS']};long=[run.task('long75_AF_X',state='N64_wide_r0',duration=.075)]
 long += [run.task('long150_'+s,state='AG32_'+s,duration=.15) for s in ['near','mid','far','random_0','random_1','random_2'] if 'pilot_'+s in admitted]
 long += [run.task('ring50_'+a,action=a,state='AG32_ring',duration=.05) for a in ['ROT_PLUS','ROT_MINUS'] if 'pilot_ring_'+a in admitted]
 long += [run.task('random50_'+a+f'_{j}',action=a,state=f'AG32_random_{j}',duration=.05) for a in ['ROT_PLUS','ROT_MINUS'] for j in range(3)]
 # Prioritize capture; unique actions and random seeds have their own real trajectories.
 write(BASE/'LONG_PLAN.json',dict(tasks=long,stop='Original gate or per-run 600 s',identity=identity()));begin=time.perf_counter();rr=batch(long);cc=[compare(t) for t in long];write(BASE/'LONG_RESULTS.json',dict(rows=rr,paired=cc,wall_s=time.perf_counter()-begin,identity=identity()))
 # Longest affordable refinement and high-gradient reference, exact original dt/cache.
 refs=[t for t in long if t['id'] in ['long150_random_0','long150_mid'] and next(v['PASS'] for v in cc if v['id']==t['id'])]
 ref_tasks=[dict(t,id=t['id']+'_'+v,variant=v) for t in refs for v in ['DT05','FINE']]
 # Refinements predicted from actual HF integration wall, with margin.
 ref_tasks=[t for t in ref_tasks if read(BASE/'hf'/t['id'].rsplit('_',1)[0]/'SUMMARY.json')['wall_s']*(2.1 if t['variant']=='DT05' else 1.2)<540]
 write(BASE/'REFINEMENT_PLAN.json',dict(tasks=ref_tasks,identity=identity()));begin=time.perf_counter();rr=batch(ref_tasks,kinds=('HF',));write(BASE/'REFINEMENT_RESULTS.json',dict(rows=rr,wall_s=time.perf_counter()-begin,identity=identity()))
 write(BASE/'EXPERIMENT_DONE.json',dict(wall_s=time.perf_counter()-start,identity=identity()))
 print('EXPERIMENT_DONE',time.perf_counter()-start,flush=True)
if __name__=='__main__':main()
