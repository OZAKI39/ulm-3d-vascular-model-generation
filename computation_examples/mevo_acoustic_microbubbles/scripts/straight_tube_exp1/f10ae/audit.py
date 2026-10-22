"""Full accepted-step audit and self-contained compact datasets, preserving failures."""
import time,concurrent.futures,multiprocessing,json,h5py,numpy as np,importlib.metadata as md
from scripts.straight_tube_exp1.f10aa.analyze import independent
from .common import *

def one(r):
 path=BASE/r['trajectory_path'];folder=path.parent;dest=BASE/'dataset'/(r['id']+'.h5');dest.parent.mkdir(exist_ok=True)
 check=independent(folder) if r['status']=='PASS' else dict(PASS=False,checked_steps=0,reason='Incomplete run: original attempted-step gates and rejected proposal retained, not admitted as full trajectory')
 with h5py.File(path) as h,h5py.File(dest,'w') as d:
  t=h['time_s'][:];ix=np.unique(np.r_[np.arange(0,len(t),10),len(t)-1]).astype(int);xx=h['positions_m'][:]
  d.attrs['failure_proposal_semantics']='Rejected integration proposal only when attempted_steps > accepted_steps; for a field-query exception the legacy proposal may repeat the previous accepted state.'
  d['time_s']=t[ix];d['positions_m']=xx[ix];d['particle_id']=h['particle_id'][:];d['raw_frame_index']=ix
  d['all_attempt_metrics']=h['step_metrics'][:];d['all_attempt_metrics'].attrs['columns']=h['step_metrics'].attrs['columns'];d['all_attempt_failed_gates']=h['step_failed_gates'][:]
  if 'failure_proposal_m' in h:d['failure_proposal_m']=h['failure_proposal_m'][:]
  if 'time_s' in h['force_samples']:
   f=h['force_samples'];jf=np.unique(np.r_[np.arange(0,len(f['time_s']),10),len(f['time_s'])-1]).astype(int)
   for k in f:d.create_dataset('force_samples/'+k,data=f[k][:][jf],compression='gzip')
   d['all_attempt_force_time_s']=f['time_s'][:];d['all_attempt_gate_pass']=f['gate_pass'][:]
  for k,v in h.attrs.items():d.attrs[k]=v
  d.attrs.update(raw_relative_path=r['trajectory_path'],raw_sha256=sha(path),derived_sampling='Exact stored frames every 10 accepted steps plus end; no interpolation',task_json=json.dumps(r['task']),scientific_status=r['scientific_status'])
  # Every accepted step legality is also inspected independently of saved success flags.
  ii,jj=np.triu_indices(xx.shape[1],1);mingap=np.inf
  for block in np.array_split(xx,max(1,len(xx)//128)):mingap=min(mingap,float(np.linalg.norm(block[:,ii]-block[:,jj],axis=-1).min()-2*R))
  check.update(minimum_actual_gap_m=mingap,actual_ids_unique=bool(np.array_equal(h['particle_id'][:],np.arange(xx.shape[1]))),no_time_padding=bool(abs(t[-1]-r['completed_s'])<1e-14))

 return dict(id=r['id'],scientific_status=r['scientific_status'],trajectory_path=r['trajectory_path'],raw_sha256=sha(path),raw_bytes=path.stat().st_size,dataset_path=str(dest.relative_to(BASE)),dataset_sha256=sha(dest),audit=check)

def main():
 begin=time.perf_counter();verify();rows=[read(p) for p in sorted((BASE/'runs').glob('*/SUMMARY.json'))];withdata=[r for r in rows if 'trajectory_path' in r]
 old=read(BASE/'AUDIT.json')['rows'] if (BASE/'AUDIT.json').exists() else [];known={r['id']:r for r in old};todo=[r for r in withdata if r['id'] not in known]
 with concurrent.futures.ProcessPoolExecutor(max_workers=4,mp_context=multiprocessing.get_context('spawn')) as pool:
  for r in pool.map(one,todo):known[r['id']]=r;write(BASE/'AUDIT_PROGRESS.json',list(known.values()))
 versions={n:dict(version=md.version(n),license=md.metadata(n).get('License-Expression',md.metadata(n).get('License'))) for n in ['ax-platform','botorch','torch','nevergrad','scipy','numpy','warp-lang','osqp']}
 write(BASE/'evidence/DEPENDENCIES.json',versions)
 write(BASE/'AUDIT.json',dict(rows=list(known.values()),attempts=len(rows),complete=sum(r['status']=='PASS' for r in rows),failures=[r for r in rows if r['status']!='PASS'],all_completed_independent_pass=all(r['audit']['PASS'] for r in known.values() if r['scientific_status']=='VALID_COMPLETE_TRAJECTORY'),checked_steps=sum(r['audit']['checked_steps'] for r in known.values()),wall_s=time.perf_counter()-begin,identity=identity()))
 write(BASE/'ALL_RUNS.json',rows)
if __name__=='__main__':main()
