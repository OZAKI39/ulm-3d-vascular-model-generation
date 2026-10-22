"""Preregister positions and bounded actions without inspecting future outcomes."""
import numpy as np,h5py,json,subprocess
from scipy.spatial.distance import pdist
from mevo_mb.control.scheduled_pzt import Segment,Schedule,one_hot,axis_pair
from mevo_mb.coupling.pdms_cavity import apparatus
from .common import *
from .state import series

def cloud(seed,half,min_distance,n=32,prior=None):
 rng=np.random.default_rng(seed);q=[] if prior is None else list(np.array(prior,copy=True));start=len(q)
 for _ in range(100000):
  a=rng.uniform(-half,half,3)
  if not q or min(np.linalg.norm(a-b) for b in q)>=min_distance:q.append(a)
  if len(q)-start==n:return np.array(q)
 raise ValueError('Registered rejection sampler exhausted')
def generate():
 original=np.load(OLD/'initial_states.npz')['aggregate_m'];states={'original':original};meta={'original':dict(seed=None,split='historical_baseline',kind='prescribed connected demo, not acoustic formation')}
 for name,seed in [('low_train',18001),('medium_train',18002),('medium_validation',28002),('medium_test',38002)]:
  x=cloud(seed,40e-6 if name.startswith('low') else 18e-6,12e-6 if name.startswith('low') else 5e-6);states[name]=x-x.mean(0);meta[name]=dict(seed=seed,split=name.split('_')[-1],kind='rejection sampled dispersed')
 dense=cloud(18003,7e-6,3.85e-6,n=16);dense=cloud(18004,30e-6,6e-6,n=16,prior=dense);states['local_dense']=dense-dense.mean(0);meta['local_dense']=dict(seed=[18003,18004],split='train',kind='16+16 prescribed density contrast; no contacts')
 groups=[]
 for i,(cx,cy) in enumerate([(-17,-14),(-17,14),(17,-14),(17,14)]):
  q=np.array([[a,b,c] for a in [-1,1] for b in [-1,1] for c in [-1,1]])*2.2e-6+np.random.default_rng(18010+i).uniform(-.08e-6,.08e-6,(8,3));groups.append(q+np.array([cx,cy,0])*1e-6)
 states['multi_cluster']=np.concatenate(groups);states['multi_cluster']-=states['multi_cluster'].mean(0);meta['multi_cluster']=dict(seed=list(range(18010,18014)),split='train',kind='four local groups with 4.4 um spacing, no near-contact graph edges')
 q=original*1.055+np.random.default_rng(18005).uniform(-.012,.012,original.shape)*R;states['near_contact']=q-q.mean(0);meta['near_contact']=dict(seed=18005,split='train',kind='prescribed nearly touching layout, not formed by acoustics')
 return states,meta

def action(name,duration=.002):
 amp=1036.2016581041792;phase=210.
 if name=='OFF':return Schedule((Segment(0,0,duration,'ONE_HOT_PZT',(0.,)*4,(0.,)*4,(False,)*4),))
 if name in ['XM','XP','YM','YP']:return Schedule((one_hot({'XM':'X-','XP':'X+','YM':'Y-','YP':'Y+'}[name],400 if name[0]=='X' else 1000,duration_s=duration),))
 if name.startswith('AMP'):amp*=1+(float(name[3:])/100)
 if name.startswith('PHASE'):phase+=float(name[5:])
 if name=='PAIR_OFF':return Schedule((axis_pair(x_amplitude=amp,x_delta_deg=phase,duration_s=duration/2),Segment(1,duration/2,duration/2,'ONE_HOT_PZT',(0.,)*4,(0.,)*4,(False,)*4)))
 return Schedule((axis_pair(x_amplitude=amp,x_delta_deg=phase,duration_s=duration),))

def main():
 BASE.mkdir(parents=True,exist_ok=True)
 if (BASE/'PREREGISTRATION.json').exists():raise RuntimeError('Preregistration immutable; do not overwrite')
 states,meta=generate();cavity,_=apparatus(read(OLD/'CONFIG_SNAPSHOT.json')['geometry']);stats={}
 with h5py.File(BASE/'INITIAL_STATES.h5','w') as h:
  for name,x in states.items():
   d=pdist(x);g=np.array([cavity.closest_wall(a)[0] for a in x])-R;m=series(np.array([0.]),x[None])
   assert d.min()>=2*R and g.min()>0 and np.abs(x).max()<54e-6
   h.create_dataset(name,data=x);stats[name]=dict(**meta[name],min_distance_m=float(d.min()),min_surface_gap_m=float(d.min()-2*R),wall_clearance_m=float(g.min()),nn_m=float(m['nn_m'][0]),rg_m=float(m['rg_m'][0]),covariance_m2=m['tensor_m2'][0],principal_axis=m['principal_axis'][0],axis_identifiable=bool(m['axis_identifiable'][0]),components=int(m['components'][0]),near_edges=int(m['edge_count'][0]),contact_count=int(m['contact_count'][0]),legal=True)
 write(BASE/'INITIAL_STATE_AUDIT.json',stats)
 schedules={name:action(name,.020).record() for name in ['OFF','XM','XP','PAIR','YM','YP','AMP-5','AMP+5','PHASE-3','PHASE+3','PHASE-6','PHASE+6']};schedules['PAIR_OFF']=action('PAIR_OFF',.002).record()
 write(BASE/'RUN_SCHEDULES.json',schedules)
 protocol=dict(labels=LABELS,training_ready=False,frequency_Hz=2.5e6,source_limit_Pa_equivalent=1500,flow='STRICT_ZERO',primary_states=['original','low_train','medium_train','local_dense','multi_cluster','near_contact'],heldout_states=['medium_validation','medium_test'],screen_s=.0005,formal_s=.002,dt_s=1e-6,half_dt_s=.5e-6,extension_s=[.003,.005,.010,.020],extension_state='medium_train',extension_action='PAIR',extension_requires='Previous level complete with unchanged numerical gates, dt/cache full-window comparison, no wall correction, <5% Rg cost, predicted next triple <180s and remaining campaign budget. Fresh integration from t=0 at every level; no stitching.',remote_budget_s=1500,per_run_budget_s=90,max_trajectories=110,initial_graph='d <= 2.1R; contact d <= 2R+1nm; proximity not fusion',signal_noise_multiplier=5,formation_rule='new edges persistent >=0.2ms (exploratory) and separately original 2ms persistence, lower components, relative Rg reduction>=1%, endpoints and edge events agree under dt/cache; not physical fusion',world_model_gate='at least 3 independent training initial states plus held-out initial states and held-out actions; converged, distinct centroid AND shape responses >=5x measured numerical differences; local normalized response rank>=3 for joint position-shape model; no random time-sample split',derivative_parameters={'pair_amplitude_percent':[-5,5],'X_plus_phase_deg':[-3,3,-6,6]},Y_scope='HOMOGENEOUS_WATER_PROXY_SPECIFIC',radius_m=R)
 write(CONFIG,protocol)
 protected={p:sha(ROOT/p) for p in PROTECTED};inputs={str((OLD/p).relative_to(ROOT)):sha(OLD/p) for p in INPUTS}
 for p in ['INITIAL_STATES.h5','RUN_SCHEDULES.json','INITIAL_STATE_AUDIT.json']:inputs[str((BASE/p).relative_to(ROOT))]=sha(BASE/p)
 inputs[str(CONFIG.relative_to(ROOT))]=sha(CONFIG)
 write(BASE/'PREREGISTRATION.json',dict(protocol=protocol,pre_results=True,input_sha256=inputs,protected_source_sha256=protected))
 # Verify old completed files directly rather than copying report numbers.
 rows=[]
 for p in sorted((ROOT/'results/straight_tube_exp1/f10s').glob('*/runs/*/RESULT.json')):
  v=read(p)
  if not (p.parent/'trajectory.h5').exists():continue
  with h5py.File(p.parent/'trajectory.h5') as h:t=h['time_s'][:];x=h['positions_m'][:];ids=h['particle_id'][:]
  assert sha(p.parent/'trajectory.h5')==v['trajectory_sha256'];assert np.array_equal(ids,np.arange(32));assert np.array_equal(x[0],states['original']);assert np.all(np.diff(t)>0)
  rows.append(dict(name=p.parent.name,complete=abs(t[-1]-v['requested_s'])<1e-13,status=v['status'],completed_s=float(t[-1]),COM_displacement_m=x[-1].mean(0)-x[0].mean(0),trajectory_sha256=sha(p.parent/'trajectory.h5')))
 assert len(rows)==11 and all(r['complete'] and r['status']=='PASS' for r in rows)
 write(BASE/'evidence/F10S_BASELINE_AUDIT.json',dict(trajectories=rows,protected_source_sha256=protected,latest_base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),verified=True))
 print('Preregistered',len(states),'states and',len(schedules),'schedules. F10S 11 trajectories verified.')
if __name__=='__main__':main()
