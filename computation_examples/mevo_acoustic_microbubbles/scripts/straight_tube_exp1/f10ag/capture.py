"""Capture metrics with bracketed entries, hold intervals and right censoring."""
import numpy as np

def capture(times,x,roots,threshold=20e-6,hold=.01,full_3d=False):
 times=np.asarray(times);x=np.asarray(x);roots=np.asarray(roots);dist=abs(x[:,:,0,None]-roots) if not full_3d else np.sqrt((x[:,:,0,None]-roots)**2+np.sum(x[:,:,1:]**2,axis=2)[:,:,None]);assignment=dist.argmin(2);distance=dist.min(2);rows=[]
 for i in range(x.shape[1]):
  inside=distance[:,i]<=threshold;ids=np.flatnonzero(inside);first=int(ids[0]) if len(ids) else None;held=None;reexit=False;stretches=[];start=None
  # Maximal contiguous intervals with a fixed band assignment, vectorized.
  cuts=np.r_[True,(~inside[:-1]) | (assignment[1:,i]!=assignment[:-1,i])]
  ends=np.r_[(~inside[1:]) | (assignment[1:,i]!=assignment[:-1,i]),True]
  stretches=list(zip(np.flatnonzero(inside & cuts),np.flatnonzero(inside & ends)))
  for s,e in stretches:
   if times[e]-times[s]>=hold-1e-12 and held is None:held=s
  if first is not None:reexit=bool((~inside[first:]).any())
  rows.append(dict(particle_id=i,initial_distance_m=float(distance[0,i]),final_distance_m=float(distance[-1,i]),initial_band=int(assignment[0,i]),final_band=int(assignment[-1,i]),first_entry_s=None if first is None else float(times[first]),entry_bracket_s=None if first is None else [float(times[max(0,first-1)]),float(times[first])],first_held_entry_s=None if held is None else float(times[held]),right_censored=held is None,reexit=reexit,completed_s=float(times[-1]),threshold_m=threshold,hold_s=hold,final_hold_s=float(times[-1]-times[stretches[-1][0]]) if stretches and stretches[-1][1]==len(times)-1 else 0.))
 return rows
