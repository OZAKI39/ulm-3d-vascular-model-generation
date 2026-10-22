"""Freeze historical states and deterministic mechanical fixtures BEFORE optimization."""
import itertools, h5py, numpy as np
from scipy.spatial.distance import pdist
from .common import *

def main():
 BASE.mkdir(parents=True,exist_ok=True); assert not (BASE/'PREREGISTRATION.json').exists()
 cfg=dict(N=[32,64,128],dt_s=1e-6,repeat_fixed=9,warmup=2,BLAS_threads=[1,2,4],
  full_windows_s={'32':.002,'64':.002,'128':.001},near_contact_window_s=.0002,
  concurrent_workers=[1,2,4],parallel_jobs=4,parallel_window_s=.0002,
  total_compute_budget_s=2400,per_trajectory_budget_s=240,max_extension_s=.002,
  extension_rule='Only if numerical and dt/cache gates pass, N128 projected triple <120s; 2ms maximum; no change of action or dt',
  fixed_fixture_rule='cubic and FCC lattices of touching-scale spheres with positive 1e-5 relative separation; manufactured compression u=-20*x with F=zeta*A*u, mechanics-only, never called acoustic trajectory',
  labels=LABELS,training_ready=False,baseline_commit='df76b1e9d311024398bc2bca40b6cf0f6382ba64')
 write(CONFIG,cfg);meta={};states={};audit=[]
 # Hash every available historical trajectory and its RESULT declaration.
 for stage in ['f10s','f10z','f10z_scaling']:
  for p in (ROOT/'results/straight_tube_exp1'/stage).glob('*/runs/*/RESULT.json'):
   raw=p.parent/'trajectory.h5'
   if raw.exists():
    d=read(p);assert sha(raw)==d['trajectory_sha256'],raw
    audit.append(dict(path=str(raw.relative_to(ROOT)),sha256=sha(raw),status=d['status']))
 with h5py.File(SCALING/'INITIAL_STATES.h5') as h:
  for n in cfg['N']:
   x=h[str(n)][:];states[f'N{n}_dispersed']=x;meta[f'N{n}_dispersed']=dict(kind='historical',force='acoustic')
   # Historical endpoint has active contact for 32/64; 128 is explicitly no-hard-contact.
   p=next((ROOT/'results/straight_tube_exp1/f10z').glob('*/runs/medium_train__PAIR__2000us__DT1/trajectory.h5')) if n==32 else next(SCALING.glob(f'*/runs/N{n}_PAIR_*us_DT1/trajectory.h5'))
   if n!=32:
    p=max(SCALING.glob(f'*/runs/N{n}_PAIR_*us_DT1/trajectory.h5'),key=lambda p:read(p.parent/'RESULT.json')['completed_s'])
   with h5py.File(p) as q:
    states[f'N{n}_historical_end']=q['positions_m'][-1]
   meta[f'N{n}_historical_end']=dict(kind='historical accepted endpoint',force='acoustic',source_sha256=sha(p))
   for kind in ['pair','cubic','fcc']:
    if kind=='pair':
     a=x.copy();a[:2]=[[-45e-6-R*(1+1e-5),0,0],[-45e-6+R*(1+1e-5),0,0]]
    elif kind=='cubic':
     a=np.array(list(itertools.product(range(8),range(4),range(4)))[:n],float)*2*R*(1+1e-5);a-=a.mean(0)
    else:
     motif=np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0]])
     a=np.concatenate([motif+[i,j,k] for i,j,k in itertools.product(range(4),range(4),range(4))])[:n]*2*np.sqrt(2)*R*(1+1e-5);a-=a.mean(0)
    assert pdist(a).min()>=2*R
    key=f'N{n}_{kind}';states[key]=a;meta[key]=dict(kind='manufactured legal mechanics fixture',force='compression',seed=None,rule=cfg['fixed_fixture_rule'])
 original=np.load(OLD/'initial_states.npz')['aggregate_m'];states['N32_original']=original;meta['N32_original']=dict(kind='F10S prescribed near-contact',force='acoustic')
 with h5py.File(BASE/'STATES.h5','w') as h:
  for name,x in states.items():
   h[name]=x;meta[name].update(N=len(x),min_gap_m=float(pdist(x).min()-2*R),max_abs_m=float(abs(x).max()))
 write(BASE/'STATE_RULES.json',meta)
 write(BASE/'BASELINE_AUDIT.json',dict(trajectories=audit,verified=True,identity=identity()))
 frozen={str(p.relative_to(ROOT)):sha(p) for p in (ROOT/'src/mevo_mb').rglob('*.py')}
 for p in [ROOT/'configs/straight_tube_exp1/f10q/protocol.json',ROOT/'configs/straight_tube_exp1/f10q/algorithm_revision_1.json']:frozen[str(p.relative_to(ROOT))]=sha(p)
 inputs={}
 for p in [CONFIG,BASE/'STATES.h5',BASE/'STATE_RULES.json',SCALING/'RUN_SCHEDULES.json',OLD/'cache/coarse.h5',OLD/'cache/fine.h5',OLD/'inputs/response.h5',OLD/'CONFIG_SNAPSHOT.json',OLD/'initial_states.npz']:
  inputs[str(p.relative_to(ROOT))]=sha(p)
 write(BASE/'PREREGISTRATION.json',dict(protocol=cfg,input_sha256=inputs,frozen_source_sha256=frozen))
 print('Frozen',len(states),'states and verified',len(audit),'historical trajectories')
if __name__=='__main__':main()
