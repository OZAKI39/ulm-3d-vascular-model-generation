"""Correlate actual immutable accepted-time tracks; no particle integration."""
from pathlib import Path
import argparse,gzip,json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from audit_field import js,sha
ap=argparse.ArgumentParser();ap.add_argument('--population-root',type=Path,required=True);a=ap.parse_args();root=a.population_root;out=root/'cfd_wall_audit'
f=pd.read_csv(out/'all_cells.csv');N=len(f);eta=f.eta_div.to_numpy();layer=f.layer.to_numpy();vol=f.volume_m3.to_numpy()
trackroot=root/'analysis/realization_0000_contactfix/tracks';raw=root/'runs/realization_0000_contactfix'
terminals={r['particle_id']:r for r in map(json.loads,gzip.open(raw/'terminals.jsonl.gz','rt'))}
states=np.zeros(N,dtype=np.int64);res=np.zeros(N);unique=np.zeros(N,dtype=np.int64);rows=[];acc=[];invalid=0;duplicates=0;position_mismatch=0
z=np.load(out/'recomputed_arrays.npz');c=z['connectivity'];xyz=z['xyz'];X=xyz[c];mat=(X[:,1:]-X[:,:1]).transpose(0,2,1);inv=np.linalg.inv(mat)
speedup=pd.read_csv(out/'inputs/visible_speedup_candidates.csv').set_index('particle_id');hashes={}
for k,path in enumerate(sorted(trackroot.glob('mb_*.npz'))):
 pid=int(path.stem.split('_')[-1]);hashes[path.name]=sha(path)
 with np.load(path) as d:t=d['physical_time_s'];tid=d['tetra_id'].astype(int);x=d['x_m'];gap=d['wall_gap_m']
 valid=(tid>=0)&(tid<N);invalid+=int((~valid).sum());assert valid.all()
 assert np.all(np.diff(t)>0);ci=tid;dtime=np.diff(t)
 # Check real positions against the recorded canonical tetra; interface roundoff has tolerance 1e-8 barycentric.
 b123=np.einsum('nij,nj->ni',inv[ci],x-X[ci,0]);b0=1-b123.sum(1);outside=(b123.min(1)<-1e-8)|(b0<-1e-8)|(b123.max(1)>1+1e-8)
 position_mismatch+=int(outside.sum())
 states+=np.bincount(ci,minlength=N);res+=np.bincount(ci[:-1],weights=dtime,minlength=N);unique[np.unique(ci)]+=1
 terminal=terminals[pid];entry=dict(particle_id=pid,state_count=len(t),residence_s=float(t[-1]-t[0]),terminal_status=terminal['status'],outlet=terminal.get('outlet'),last_tetra_id=int(ci[-1]),near_wall_state_count=int((layer[ci]<=2).sum()),W0_state_count=int((layer[ci]==0).sum()),eta_gt_0p1_state_count=int((eta[ci]>.1).sum()),eta_gt_0p5_state_count=int((eta[ci]>.5).sum()),eta_gt_0p1_residence_s=float(dtime[eta[ci[:-1]]>.1].sum()),eta_gt_0p5_residence_s=float(dtime[eta[ci[:-1]]>.5].sum()),eta_max=float(eta[ci].max()))
 rows.append(entry)
 if pid in speedup.index:
  r=speedup.loc[pid];window=(t>=r.physical_time_s-r.window_s)&(t<=r.physical_time_s+r.window_s);idx=int(np.argmin(abs(t-r.physical_time_s)))
  before=(t>=r.physical_time_s-r.window_s)&(t<r.physical_time_s)
  acc.append(dict(particle_id=pid,physical_time_s=r.physical_time_s,video_time_s=r.video_time_s,onset_tetra=int(ci[idx]),onset_eta=float(eta[ci[idx]]),window_tetra_count=len(np.unique(ci[window])),window_any_eta_gt_0p1=bool(np.any(eta[ci[window]]>.1)),window_any_eta_gt_0p5=bool(np.any(eta[ci[window]]>.5)),before_eta_median=float(np.median(eta[ci[before]])) if before.any() else None,window_W0_state_fraction=float(np.mean(layer[ci[window]]==0))))
 if k%200==0:print(k,'tracks processed',flush=True)
pd.DataFrame(dict(tetra_id=np.arange(N),layer=layer,eta_div=eta,state_count=states,residence_s_start_cell_assignment=res,unique_particles=unique,visited=unique>0)).to_csv(out/'particle_overlap.csv',index=False)
particles=pd.DataFrame(rows);particles.to_csv(out/'particle_exposure.csv',index=False);acc=pd.DataFrame(acc);acc.to_csv(out/'acceleration_overlap.csv',index=False)
groups=[]
for l in range(3):
 mask=layer==l;visited=mask&(unique>0)
 for threshold in [.1,.5]:
  bad=eta>threshold
  groups.append(dict(layer=f'W{l}',eta_strictly_greater_than=threshold,all_cell_count=int(mask.sum()),all_flag_fraction=float(np.mean(bad[mask])),visited_cell_count=int(visited.sum()),visited_flag_fraction=float(np.mean(bad[visited])),particle_state_count=int(states[mask].sum()),flag_state_count=int(states[mask&bad].sum()),flag_state_fraction=float(states[mask&bad].sum()/states[mask].sum()),residence_s=res[mask].sum(),flag_residence_fraction=float(res[mask&bad].sum()/res[mask].sum())))
pd.DataFrame(groups).to_csv(out/'particle_layer_exposure.csv',index=False)
assoc=[]
for (status,outlet),df in particles.fillna({'outlet':'NONE'}).groupby(['terminal_status','outlet']):
 assoc.append(dict(status=status,outlet=outlet,particle_count=len(df),median_residence_s=df.residence_s.median(),median_fraction_time_eta_gt_0p5=float(np.median(df.eta_gt_0p5_residence_s/df.residence_s)),ever_eta_gt_0p5_count=int((df.eta_gt_0p5_state_count>0).sum())))
pd.DataFrame(assoc).to_csv(out/'outcome_associations.csv',index=False)
summary=dict(track_count=len(rows),state_count=int(states.sum()),visited_cell_count=int((unique>0).sum()),near_wall_visited_cells=int(((unique>0)&(layer<=2)).sum()),state_eta_gt_0p1=int(states[eta>.1].sum()),state_eta_gt_0p5=int(states[eta>.5].sum()),state_eta_gt_0p1_fraction=float(states[eta>.1].sum()/states.sum()),state_eta_gt_0p5_fraction=float(states[eta>.5].sum()/states.sum()),acceleration_cases=len(acc),acceleration_window_eta_gt_0p5=int(acc.window_any_eta_gt_0p5.sum()),acceleration_onset_eta_gt_0p5=int((acc.onset_eta>.5).sum()),failure_rows=particles[particles.terminal_status=='SOLVER_FAILURE'].replace({np.nan:None}).to_dict('records'),residence_vs_exposure_fraction_spearman=float(spearmanr(particles.residence_s,particles.eta_gt_0p5_residence_s/particles.residence_s).statistic),position_outside_recorded_tetra_barycentric_1e_minus8=position_mismatch,invalid_tetra_ids=invalid,residence_definition='Sum each accepted interval duration assigned to its START tetra; a path may traverse tetra within an interval. This is an endpoint-cell occupancy approximation, not exact geometric residence integration.',association_not_causation=True,old_experiment_outcome_dependent_stop=True)
js(out/'PARTICLE_OVERLAP_SUMMARY.json',summary);js(out/'validation/trajectory_hashes.json',dict(hashes=hashes,manifest_sha256=sha(raw/'MANIFEST.json')))
print(json.dumps(summary,default=lambda x:x.item()),flush=True)
