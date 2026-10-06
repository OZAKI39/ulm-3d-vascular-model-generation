"""Read-only Stage4L evidence audit and fixed-state profiling, no new trajectory."""
from pathlib import Path
import json, hashlib, time, collections, os, platform, importlib, subprocess
import numpy as np
import h5py
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
ROOT=Path('/workspace/mevo_acoustic_microbubbles_stage5a')
OUT=ROOT/'results/stage5/motion_solver/audit_v2'
OUT.mkdir(parents=True,exist_ok=False)
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False,default=lambda a:a.item())+'\n')
def sha(p):return hashlib.file_digest(open(p,'rb'),'sha256').hexdigest()
def dist(a):
 a=np.asarray(a);return dict(count=len(a),min=float(np.min(a)),p05=float(np.quantile(a,.05)),median=float(np.median(a)),p95=float(np.quantile(a,.95)),max=float(np.max(a))) if len(a) else None
paths=json.load(open(ROOT/'baseline_paths.json'));folders=sorted(set(str(Path(v['original_server_path']).parent) for v in paths));rows={};identity=[]
for folder in folders:
 p=Path(folder);s=json.load(open(p/'SUMMARY.json'))
 if s['status']!='COMPLETE':continue
 key=f"{s['case']}_{s['requested_duration_s']*1000:g}ms"
 inputs=[]
 for record in paths:
  q=Path(record['original_server_path'])
  if q.parent!=p:continue
  digest=sha(q);expected=record.get('sha256')
  if expected and digest!=expected:raise RuntimeError(f'HASH_MISMATCH {q}')
  inputs.append(dict(path=str(q),sha256=digest,expected_sha256=expected,verified=expected is None or digest==expected,prior_digest_available=expected is not None))
 steps=[json.loads(line) for line in open(p/'steps.jsonl')];solves=[json.loads(line) for line in open(p/'TRANSLATION_SOLVES.jsonl')] if (p/'TRANSLATION_SOLVES.jsonl').exists() else []
 n=len(steps);ident=json.load(open(p/'IDENTITY.json'));nominal=ident['config']['particles']['dt_s'];t0=ident['initial_time_s'];before=np.r_[t0,[v['time_s'] for v in steps[:-1]]];rem=t0+(np.floor((before-t0+1e-14)/nominal)+1)*nominal-before
 rem=np.minimum(rem,t0+s['requested_duration_s']-before)
 ref=collections.Counter(v.get('contact_refinement','none') for v in solves)
 row=dict(run=str(p),inputs=inputs,summary=s,accepted_steps=n,evaluate_mobility_calls=sum('solver' in v for v in solves),contact_calls=len(solves),projection_calls=sum(v['projection_calls'] for v in steps),dt_s=dist([v['dt_s'] for v in steps]),limiter_trigger_fraction=float(np.mean(np.array([v['dt_s'] for v in steps])<rem*(1-1e-7))),constraint_count=dist([v['contact_constraints'] for v in solves]),matrix_size=dist([v['matrix_size'] for v in solves if 'matrix_size' in v]),recorded_contact_s=sum(v.get('contact_solve_time_s',0) for v in solves),recorded_assembly_and_solve_s=sum(v.get('assembly_and_solve_time_s',0) for v in solves),recorded_linear_solve_s=sum(v.get('solve_time_s',0) for v in solves),qr_refinement_counts=dict(ref),lbfgs_iterations_historical=None,unrecorded_historical=['wall-query time','response time','pair-force time','L-BFGS-B iterations','component distribution'])
 row['historical_translation_log_available']=bool(solves)
 if not solves:
  for k in ['evaluate_mobility_calls','contact_calls','recorded_contact_s','recorded_assembly_and_solve_s','recorded_linear_solve_s']:row[k]=None
 rows[key]=row
 identity.append(dict(key=key,run=str(p),identity_sha256=sha(p/'IDENTITY.json')))
write(OUT/'HISTORICAL_AUDIT.json',rows)
write(OUT/'BASELINE_SELECTION.json',identity)
# Fresh resource/dependency observation. No package installation.
env=dict(python=platform.python_version(),cpu_quota=Path('/sys/fs/cgroup/cpu.max').read_text().strip(),memory_limit_bytes=int(Path('/sys/fs/cgroup/memory.max').read_text()),cpuset=Path('/sys/fs/cgroup/cpuset.cpus.effective').read_text().strip(),logical_cpus=os.cpu_count(),threads={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']})
for name in ['numpy','scipy','cupy','osqp','proxsuite']:
 try: m=importlib.import_module(name);env[name]=m.__version__
 except ImportError:env[name]='not installed'
env['gpu']=subprocess.check_output(['nvidia-smi','--query-gpu=name,memory.total,memory.used,utilization.gpu,driver_version','--format=csv'],text=True)
write(OUT/'ENVIRONMENT.json',env)
from scripts.stage4l.dense_probe import load_case
from mevo_mb.coupling.integrators import GlobalSynchronousIntegrator
from mevo_mb.bubbles.ordered_contact_sweep import OrderedContactSweep
import mevo_mb.bubbles.tracker as tracking
import mevo_mb.bubbles.translation_backends as tb
from mevo_mb.bubbles.translation_backends import _pairs,_gap_vectors
profile=[]
for case in ['B2','B3']:
 p=Path(rows[case+'_5ms']['run']);config,field,flow,gold,inputs=load_case(Path('/workspace/mevo_acoustic_microbubbles'),case)
 config['particles'].update(translation_backend='lubrication_coupled',lubrication=dict(activation=.1,floor=.001,residual_tolerance=1e-11))
 tracker=GlobalSynchronousIntegrator(config,field,flow);sweep=OrderedContactSweep(OUT/'ordered_cache');tracking.project_contacts=sweep.project
 stat=collections.defaultdict(lambda:dict(calls=0,seconds=0.));opt=[]
 def wrap(obj,name,key):
  orig=getattr(obj,name)
  def call(*a,**kw):
   start=time.perf_counter()
   try:return orig(*a,**kw)
   finally:stat[key]['calls']+=1;stat[key]['seconds']+=time.perf_counter()-start
  setattr(obj,name,call)
 for obj,name,key in [(tracker,'evaluate','evaluate'),(tracker.translation,'solve','mobility'),(tracker.translation,'contact','contact'),(tracker,'_project','projection'),(tracker.cylinder,'wall_distance','wall_distance'),(tracker.cylinder,'wall_gap_and_normal','wall_gap'),(tracker.cylinder,'project','wall_project')]:wrap(obj,name,key)
 original_min=tb.minimize
 def mini(*a,**kw):
  fit=original_min(*a,**kw);opt.append(dict(nit=int(fit.nit),nfev=int(fit.nfev)));return fit
 tb.minimize=mini
 with h5py.File(p/'trajectories.h5','r') as h:
  times=h['time_s'][:]
  for index in np.unique(np.linspace(0,len(times)-2,9,dtype=int)):
   mask=h['active'][index];x=h['positions_m'][index][mask];r=h['radii_m'][:][mask];t=float(times[index]);dt=float(times[index+1]-t)
   stat.clear();opt.clear();tracker.timings={k:0 for k in tracker.timings};tracker.translation_diagnostics.clear();start=time.perf_counter()
   state=tracker.evaluate(x,r,t);initial=tracker.last_translation
   # Exact union of lubrication and the baseline finite-horizon broad phase.
   vel=state['total_velocity_m_s'];reach=2*np.linalg.norm(vel-vel.mean(axis=0),axis=1).max()*dt
   pairs=_pairs(x,max((2+.1)*r[0],2*r[0]+.001*r[0]+reach));g=coo_matrix((np.ones(len(pairs)),(pairs[:,0],pairs[:,1])),shape=(len(x),len(x)));count,labels=connected_components(g,directed=False);sizes=np.bincount(labels)
   new,_,_,_,_=tracker.step(x,r,t,dt,state);tracker.evaluate(new,r,t+dt)
   profile.append(dict(case=case,source_run=str(p),frame_index=int(index),time_s=t,dt_s=dt,wall_s=time.perf_counter()-start,calls=dict(stat),tracker_timings=tracker.timings.copy(),lbfgs=list(opt),component_sizes=sizes.tolist(),largest_component_fraction=float(max(sizes)/len(x)),constraint_counts=[v.get('contact_constraints') for v in tracker.translation_diagnostics],refinements=[v.get('contact_refinement') for v in tracker.translation_diagnostics]))
   write(OUT/'FIXED_STATE_PROFILES.json',profile)
 tb.minimize=original_min
 print(case,'profile complete',flush=True)
write(OUT/'AUDIT_COMPLETE.json',dict(status='COMPLETE',scientific_time_advanced=False,scope='Archived files read and hashed; independent single-step replays only, not concatenated trajectories',baseline_identity=identity))
print('AUDIT_COMPLETE',flush=True)
