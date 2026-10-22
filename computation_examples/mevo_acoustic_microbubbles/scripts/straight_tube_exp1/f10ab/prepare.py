"""Freeze finite shape/action design before any F10AB physical integrations."""
import itertools,h5py,numpy as np
from scipy.spatial.distance import pdist
from scipy.spatial.transform import Rotation
from mevo_mb.control.scheduled_pzt import Segment,Schedule,one_hot,axis_pair
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from .metrics import geometry
from .common import *

def grow(seed,n,axes):
 rng=np.random.default_rng(seed);axes=np.array(axes)*1e-6;pts=[np.zeros(3)]
 for trial in range(300000):
  direction=rng.normal(size=3);direction/=np.linalg.norm(direction)
  p=pts[rng.integers(len(pts))]+direction*rng.uniform(3.62e-6,3.70e-6)
  if np.sum((p/axes)**2)>1 or np.min(np.linalg.norm(np.asarray(pts)-p,axis=1))<3.60e-6:continue
  pts.append(p)
  if len(pts)==n:
   x=np.array(pts);return x-x.mean(0)
 raise RuntimeError('Initial geometry construction exhausted attempts')

def shape(kind,seed,n=64):
 if kind=='S1':x=grow(seed,n,[11,11,11])
 elif kind in ['S2','S3','S4']:
  x=grow(seed,n,[24,8,7])
  if kind=='S3':x=x@Rotation.from_euler('z',90,degrees=True).as_matrix().T
  if kind=='S4':x=x@Rotation.from_euler('zy',[35,12],degrees=True).as_matrix().T
 elif kind=='S5':
  a=grow(seed,n//2,[9,8,7]);b=grow(seed+10000,n//2,[9,8,7]);x=np.r_[a+[-14e-6,0,0],b+[14e-6,0,0]]
 else:
  from scripts.straight_tube_exp1.f10z.prepare import cloud
  x=cloud(seed,18e-6*(n/32)**(1/3),5e-6,n=n)
 return x-x.mean(0)

def schedules():
 A=1036.2016581041792;t=.01
 s={'A0_OFF':Segment(0,0,t,'ONE_HOT_PZT',(0,)*4,(0,)*4,(False,)*4)}
 s.update(A1_XM=one_hot('X-',A,duration_s=t),A2_XP=one_hot('X+',A,duration_s=t))
 for name,deg in [('A3_X0',0),('A4_X90',90),('A5_X180',180),('A6_X270',270),('A7_X210',210)]:s[name]=axis_pair(x_amplitude=A,x_delta_deg=deg,duration_s=t)
 s.update(A8_YM=one_hot('Y-',A,duration_s=t),A9_YP=one_hot('Y+',A,duration_s=t),A10_Y90=axis_pair(y_amplitude=A,y_delta_deg=90,x_enabled=False,y_enabled=True,duration_s=t))
 s['A11_FOUR']=Segment(0,0,t,'FOUR_CHANNEL_INDEPENDENT',(A,)*4,(0,np.pi/2,0,7*np.pi/6),(True,)*4)
 # Fixed L2 source budget: sum |u|^2 equals a single A-amplitude source.
 for name,deg in [('B3_X0_L2',0),('B4_X90_L2',90),('B7_X210_L2',210)]:s[name]=axis_pair(x_amplitude=A/np.sqrt(2),x_delta_deg=deg,duration_s=t)
 return {k:Schedule((v,)).record() for k,v in s.items()}

def main():
 assert not (BASE/'PREREGISTRATION.json').exists()
 cfg=dict(radius_m=R,frequency_Hz=2.5e6,flow='STRICT_ZERO',channel_order=['Y-','Y+','X-','X+'],source_limit_Pa_equivalent=1500,
  stages=['static','pilot','discovery','confirmation','sensitivity','target_if_admitted'],total_server_budget_s=3600,per_run_budget_s=240,max_workers=2,
  dt_s=1e-6,half_dt_s=.5e-6,cache_spacing_m=[3e-6,1.5e-6],short_s=.0005,formal_s=.002,max_extension_s=.005,
  metric_stride_steps=10,physical_output='every accepted step',discovery_replicates=[0],confirmation_replicates=[1,2],
  shape_response_rule='Similarity residual > max(5*dt/cache residual uncertainty,1 nm); deviatoric/shear response > max(5*uncertainty,1e-4); repeat signed response on two held-out realizations; require full windows and unchanged numerical gates.',
  action_difference_rule='Pairwise anisotropic/shear output difference > max(5*sum numerical uncertainties,1e-4); centroid similarity <=20 nm is diagnostic, not target success.',
  orientation_eigen_gap_relative=.03,affine_max_condition=100,affine_max_relative_residual=.1,
  discovery_selection='Static diagnostic selects at most four non-OFF actions: historical A7, max free deviatoric rate, opposite/min rate, and L2-normalized control when legal. Fixed action library retained even if not integrated.',
  extension_rule='Only a selected numerically resolved repeated case; predicted triple cost <180 s, remaining budget >300 s, no unresolved contact failure; max 5ms.',
  sensitivity_rule='Central differences in pair amplitude (5%,2.5%) and relative X+ phase (6deg,3deg); at two windows or initial shapes; base A7. Two explored parameters cannot prove four-channel global controllability.',
  target_rule='Only after full-window dt/cache and held-out signed-response checks. Freeze target and <=4 candidate schedules before target runs; max eight legal contiguous segments.',
  matching_priority='Exact N/radius/center; then requested geometry and legal connectivity; Rg/density/minimum distances measured, not forced to match.',labels=LABELS,training_ready=False)
 write(CONFIG,cfg);c=context();bases=[c['basis'],CachedWaterBasis(OLD/'cache/fine.h5')];meta={};states={'S0_r0':c['x0']}
 for k in range(1,7):
  for rep in range(3):states[f'S{k}_r{rep}']=shape(f'S{k}',281000+k*100+rep)
 with h5py.File(BASE/'INITIAL_STATES.h5','w') as h:
  for name,x in states.items():
   h[name]=x;m=geometry(x);wall=min(c['cavity'].closest_wall(a)[0]-R for a in x)
   legal=pdist(x).min()>=2*R-1e-12 and wall>100e-6
   for b in bases:
    q=(x-b.origin)/b.spacing;legal=legal and bool(np.all(q>=1) and np.all(q<np.array(b.shape)-2))
   assert legal,name
   k,rep=name.split('_r');meta[name]=dict(shape=k,replicate=int(rep),seed=None if k=='S0' else 281000+int(k[1:])*100+int(rep),N=len(x),role='historical' if k=='S0' else ('discovery' if rep=='0' else 'held_out_confirmation'),construction='prescribed initial geometry, not acoustically formed',legal=legal,wall_clearance_m=wall,**m)
 write(BASE/'INITIAL_STATE_AUDIT.json',meta);write(BASE/'ACTIONS.json',schedules())
 protected={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'src/mevo_mb').rglob('*.py')}
 for p in [ROOT/'configs/straight_tube_exp1/f10q/protocol.json',ROOT/'configs/straight_tube_exp1/f10q/algorithm_revision_1.json']:protected[str(p.relative_to(ROOT))]=sha(p)
 inputs={str(p.relative_to(ROOT)):sha(p) for p in [CONFIG,BASE/'INITIAL_STATES.h5',BASE/'INITIAL_STATE_AUDIT.json',BASE/'ACTIONS.json',OLD/'cache/coarse.h5',OLD/'cache/fine.h5',OLD/'inputs/response.h5',OLD/'CONFIG_SNAPSHOT.json',OLD/'initial_states.npz']}
 history=[]
 for stage in ['f10s','f10z','f10z_scaling','f10aa']:
  for p in (ROOT/'results/straight_tube_exp1'/stage).rglob('trajectory.h5'):
   with h5py.File(p) as h:
    history.append(dict(path=str(p.relative_to(ROOT)),sha256=sha(p),N=h['positions_m'].shape[1],samples=len(h['time_s']),completed_s=float(h['time_s'][-1]),has_fixed_id='particle_id' in h))
 write(BASE/'BASELINE_AUDIT.json',dict(baseline_commit='6028598d0012912cc19aa766fa786c4afb0a91a5',historical_trajectories=history,source_sha256=protected))
 write(BASE/'PREREGISTRATION.json',dict(protocol=cfg,model_sha256=protected,input_sha256=inputs,pre_physics_results=True,identity=identity()))
 print({k:(v['proximity_components'],round(v['Rg_3d_m']*1e6,3),round(v['aspect_3d'],3)) for k,v in meta.items()})
if __name__=='__main__':main()
