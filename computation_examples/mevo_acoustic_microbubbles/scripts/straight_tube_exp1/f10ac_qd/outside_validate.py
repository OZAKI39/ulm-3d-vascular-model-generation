"""One retained out-of-archive candidate; no archive-bound or comparison changes."""
import time,h5py,numpy as np
from scripts.straight_tube_exp1.f10ab.prepare import shape
from scripts.straight_tube_exp1.f10z_scaling.analyze import compare
from .common import *
from .run import batch

def main():
 start=time.perf_counter();candidates=[]
 for p in (BASE/'search_v2').glob('*/RESULT.json'):
  candidates += [r for r in read(p)['rows'] if r['status']=='PASS' and r['descriptor'] is not None and r['descriptor'][1]<-.01 and r['descriptor'][0]<-.08]
 if not candidates:write(BASE/'OUTSIDE_VALIDATION.json',dict(status='NOT_RUN',reason='No completed negative-aspect result outside the frozen archive'));return
 chosen=min(candidates,key=lambda r:r['descriptor'][1]);seeds=[401801,401802,401803];library=BASE/'OUTSIDE_INITIAL_STATES.h5'
 with h5py.File(BASE/'INITIAL_STATES.h5') as src,h5py.File(library,'x') as h:
  h['OD_search']=src['S2_r0'][:]
  for i,seed in enumerate(seeds):h['OD_unseen_'+str(i)]=shape('S2',seed)
 plans=[]
 for state in ['OD_search','OD_unseen_0','OD_unseen_1','OD_unseen_2']:
  for variant in ['DT1','DT05','FINE']:
   t=chosen['task'].copy();t.update(id='outside_'+state+'_'+variant,state=state,variant=variant,initial_library=library.name);plans.append(t)
 write(BASE/'OUTSIDE_PREREGISTRATION.json',dict(chosen=chosen,seeds=seeds,plans=plans,input_sha256=sha(library),selection='Most negative completed aspect response outside the original Rg archive bound; one candidate, no new optimization, no archive resizing',criterion=read(CONFIG)['numerical_floor'],budget_s=180,pre_outside_validation=True,identity=identity()))
 rs=batch(plans,4);write(BASE/'OUTSIDE_RUNS.json',rs);checks=[];rep=[]
 for state in ['OD_search','OD_unseen_0','OD_unseen_1','OD_unseen_2']:
  a=next(r for r in rs if r['id']=='outside_'+state+'_DT1');cc=[];noise=np.zeros(2);resnoise=0
  for v in ['DT05','FINE']:
   b=next(r for r in rs if r['id']=='outside_'+state+'_'+v);z=compare((BASE/a['trajectory_path']).parent,(BASE/b['trajectory_path']).parent) if 'trajectory_path' in a and 'trajectory_path' in b else dict(admitted=False,reason='Initialization failed');z.update(reference=a['id'],comparison=b['id']);cc.append(z);checks.append(z)
   if a['descriptor'] is not None and b['descriptor'] is not None:
    noise=np.maximum(noise,abs(np.array(a['descriptor'])-b['descriptor']));resnoise=max(resnoise,abs(a['final']['similarity_rmsd_3d_m']-b['final']['similarity_rmsd_3d_m']))
  numeric=all(c['admitted'] for c in cc);d=np.array(a['descriptor']) if a['descriptor'] is not None else np.zeros(2)
  rep.append(dict(id=a['id'],state=state,descriptor=a['descriptor'],noise_bound=noise,residual_noise_m=resnoise,numerical_pass=numeric,resolved=numeric and a['final']['similarity_rmsd_3d_m']>max(1e-9,5*resnoise),descriptor_resolved=(abs(d)>np.maximum(1e-4,5*noise)).tolist(),connected=a.get('connected_at_end',False),components=a.get('final',{}).get('proximity_components'),similarity_nm=a.get('final',{}).get('similarity_rmsd_3d_m',0)*1e9))
 repeated=all(r['resolved'] and all(r['descriptor_resolved']) and all(np.array(r['descriptor'])<0) for r in rep)
 group=dict(elite='OD',source=chosen['id'],frequency=chosen['task'].get('frequency',2500000.),replicates=rep,signed_descriptor_repeated=repeated,independent_seeds=3,connected_all=all(r['connected'] for r in rep),scope='Retained outside-archive distribution response; never counted as QD coverage or an archived elite')
 write(BASE/'OUTSIDE_VALIDATION.json',dict(status='COMPLETE',group=group,comparisons=checks,wall_s=time.perf_counter()-start,identity=identity()))
if __name__=='__main__':main()
