"""Freeze representative elites before unseen seeds and numerical refinements."""
import time,h5py,numpy as np
from scripts.straight_tube_exp1.f10z_scaling.analyze import compare
from .common import *
from .run import task,batch
from .search import eligible

def main():
 begin=time.perf_counter();rows=[]
 for p in (BASE/'search_v2').glob('*/RESULT.json'):rows+=read(p)['rows']
 for p in (BASE/'extensions').glob('*/RESULT.json'):rows+=read(p)['rows']
 uniq={r['id']:r for r in rows if eligible(r) and r['final']['similarity_rmsd_3d_m']>1e-9};candidates=list(uniq.values());chosen=[]
 # Pick descriptor extremes in separate directions, without inspecting held-out outcomes.
 for direction in [np.array([0,1]),np.array([0,-1]),np.array([-1,0])]:
  pool=[r for r in candidates if r['id'] not in [c['id'] for c in chosen]]
  if pool:chosen.append(max(pool,key=lambda r:np.dot(r['descriptor'],direction)))
 plans=[]
 for k,r in enumerate(chosen):
  for state in ['S2_r0','S2_unseen_0','S2_unseen_1','S2_unseen_2']:
   for variant in ['DT1','DT05','FINE']:
    t=r['task'].copy();t.update(id=f'elite{k}_{state}_{variant}',state=state,variant=variant);plans.append(t)
 write(BASE/'ELITE_VALIDATION_PREREGISTRATION.json',dict(selected=chosen,plans=plans,criterion=read(CONFIG)['numerical_floor'],no_heldout_selection=True,identity=identity()))
 outputs=batch(plans,4);write(BASE/'ELITE_VALIDATION_RUNS.json',outputs);checks=[];groups=[]
 for k,original in enumerate(chosen):
  rep=[]
  for state in ['S2_r0','S2_unseen_0','S2_unseen_1','S2_unseen_2']:
   a=next(r for r in outputs if r['id']==f'elite{k}_{state}_DT1');cc=[];noise=np.zeros(2);res_noise=0
   for v in ['DT05','FINE']:
    b=next(r for r in outputs if r['id']==f'elite{k}_{state}_{v}')
    c=compare((BASE/a['trajectory_path']).parent,(BASE/b['trajectory_path']).parent) if 'trajectory_path' in a and 'trajectory_path' in b else dict(admitted=False,reason='INITIALIZATION_FAILURE')
    c.update(reference=a['id'],comparison=b['id']);checks.append(c);cc.append(c)
    if a['descriptor'] is not None and b['descriptor'] is not None:
     noise=np.maximum(noise,abs(np.array(a['descriptor'])-b['descriptor']));res_noise=max(res_noise,abs(a['final']['similarity_rmsd_3d_m']-b['final']['similarity_rmsd_3d_m']))
   numeric=all(c['admitted'] for c in cc);d=np.array(a['descriptor']) if a['descriptor'] is not None else np.zeros(2)
   rep.append(dict(id=a['id'],state=state,descriptor=a['descriptor'],noise_bound=noise,residual_noise_m=res_noise,numerical_pass=numeric,resolved=numeric and a['final']['similarity_rmsd_3d_m']>max(1e-9,5*res_noise),descriptor_resolved=(abs(d)>np.maximum(1e-4,5*noise)).tolist(),connected=a.get('connected_at_end',False),components=a.get('final',{}).get('proximity_components'),similarity_nm=a.get('final',{}).get('similarity_rmsd_3d_m',0)*1e9))
  signs=[np.sign(r['descriptor']).tolist() if r['descriptor'] is not None else None for r in rep]
  repeated=len(rep)==4 and all(r['resolved'] for r in rep) and all(s==signs[0] for s in signs)
  groups.append(dict(elite=k,source=original['id'],frequency=original['task'].get('frequency',2500000.),replicates=rep,signed_descriptor_repeated=repeated,independent_seeds=3,connected_all=all(r['connected'] for r in rep),scope='Free-distribution deformation; connected group evaluated separately; geometric breakup is not material fracture'))
 write(BASE/'VALIDATION.json',dict(groups=groups,comparisons=checks,wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
