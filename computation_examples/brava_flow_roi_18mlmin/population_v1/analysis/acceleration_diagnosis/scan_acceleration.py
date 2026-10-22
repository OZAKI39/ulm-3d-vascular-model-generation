"""Read-only scan of physical-time animation and accepted trajectory endpoints."""
from pathlib import Path
import numpy as np,json,csv
root=Path(__file__).resolve().parents[2];out=Path(__file__).resolve().parent
base=root/'analysis/realization_0000_contactfix';d=np.load(base/'visualization_data.npz')
t=d['frame_time_s'];x=d['frame_x_m'];ids=d['particle_id'];dt=np.diff(t);speed=np.linalg.norm(np.diff(x,axis=0),axis=2)/dt[:,None]
rows=[];w=3
for j,pid in enumerate(ids):
 z=speed[:,j];best=None
 for k in range(w,len(z)-w):
  before=z[k-w:k];after=z[k:k+w]
  if not(np.isfinite(before).all() and np.isfinite(after).all()):continue
  a=float(before.mean());b=float(after.mean());ratio=b/max(a,1e-15)
  if a<.002 and b>.01 and ratio>=5:
   row=dict(particle_id=int(pid),frame_index=k,physical_time_s=float(t[k]),video_time_s=k/24,
    before_m_s=a,after_m_s=b,ratio=ratio,delta_m_s=b-a,window_s=float(t[k]-t[k-w]))
   if best is None or row['ratio']>best['ratio']:best=row
 if best:rows.append(best)
rows.sort(key=lambda r:-r['ratio'])
with (out/'visible_speedup_candidates.csv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
print('ANIMATION_TIMING',dict(frames=len(t),step_min=dt.min(),step_max=dt.max(),physical_per_video_s=float(dt.mean()*24)),flush=True)
print('SLOW_FAST_TRACKS',len(rows),'TOP',json.dumps(rows[:12]),flush=True)
jumps=[];max_single_step_increment=0.
for n,p in enumerate(sorted((base/'tracks').glob('*.npz'))):
 with np.load(p) as z:
  tt=z['physical_time_s'];v=z['velocity_m_s'];s=np.linalg.norm(v,axis=1);delta=np.diff(s)
  max_single_step_increment=max(max_single_step_increment,float(delta[1:].max(initial=0.)))
  keep=(delta>.005)&(s[1:]>2*np.maximum(s[:-1],1e-14));keep[:1]=False
  ii=np.flatnonzero(keep)
  if len(ii):
   i=int(ii[np.argmax(delta[ii])])+1
   jumps.append(dict(particle_id=int(p.stem[3:]),index=i,evaluated_time_s=float(tt[i-1]),
    before_m_s=float(s[i-1]),after_m_s=float(s[i]),ratio=float(s[i]/max(s[i-1],1e-15)),
    delta_m_s=float(s[i]-s[i-1]),time_between_evaluations_s=float(tt[i-1]-tt[i-2]),
    displacement_between_evaluations_m=float(np.linalg.norm(z['x_m'][i-1]-z['x_m'][i-2])),
    gap_before_m=float(z['wall_gap_m'][i-2]),gap_after_m=float(z['wall_gap_m'][i-1]),
    tetra_before=int(z['tetra_id'][i-2]),tetra_after=int(z['tetra_id'][i-1])))
 if n%300==0:print('SCANNED',n,flush=True)
jumps.sort(key=lambda r:-r['delta_m_s'])
with (out/'adjacent_step_jumps.csv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(jumps[0]) if jumps else ['particle_id','index']);wr.writeheader();wr.writerows(jumps)
result=dict(animation_timing=dict(frames=len(t),fps=24,minimum_frame_physical_step_s=float(dt.min()),maximum_frame_physical_step_s=float(dt.max()),physical_per_video_second=float(dt.mean()*24)),
 selection='3 adjacent video intervals before/after: before<2 mm/s, after>10 mm/s, ratio>=5; one largest ratio per ID',
 visible_slow_fast_particle_count=len(rows),visible_top=rows[:20],adjacent_step_jump_particle_count=len(jumps),adjacent_jump_top=jumps[:20],maximum_single_step_speed_increase_m_s=max_single_step_increment)
(out/'SCAN_SUMMARY.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result),flush=True)
