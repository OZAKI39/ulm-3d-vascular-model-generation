"""Audit real ladder prefixes and stop redundant unpaired work after cost evidence."""
import time
import numpy as np,h5py
from .common import *
from .experiment import batch

def main():
 begin=time.perf_counter();verify();assert (BASE/'HF250_RESULTS.json').exists();run,slow,_,_=bind();rows=[];checkpoints=[];tasks=[];decisions=[]
 for state in ['random_0','mid','near']:
  r=read(BASE/'lf'/('long150_'+state)/'SUMMARY.json');rows.append(r)
  with h5py.File(BASE/r['trajectory_path']) as h:
   times=h['time_s'][:]
   for T in [.05,.075,.1,.15]:
    ix=min(max(np.searchsorted(times,T+1e-13,side='right')-1,0),len(times)-1);x=h['positions_m'][ix]
    checkpoints.append(dict(state=state,target_checkpoint_s=T,actual_accepted_time_s=float(times[ix]),raw_node_index=int(ix),source_trajectory=r['trajectory_path'],source_sha256=r['trajectory_sha256'],minimum_wall_gap_m=float((.0008-np.linalg.norm(x[:,1:],axis=1)-R).min()),scope='Actual accepted prefix of completed 150 ms LF; not a separately rerun trajectory, no padding or interpolation'))
  existing=BASE/'lf'/('long250_'+state)/'SUMMARY.json'
  if existing.exists():rows.append(read(existing))
  else:
   fallback_fraction=r['HF_fallback_steps']*1e-6/r['completed_s'];predicted=r['wall_s']*.25/.15
   # No further unpaired extension: actual paired 250 ms cases already test the
   # longest registered window and spend >90% of physical time in original HF.
   decisions.append(dict(state=state,status='NOT_RUN',candidate_s=.25,actual_150ms_fallback_fraction=fallback_fraction,predicted_cost_s=predicted,reason='Existing 250 ms HF/LF pairs suffice; near-contact fallback removes the long-time cost benefit. Preserve 150 ms mid state as last paired evidence.'))
 write(BASE/'LF_LADDER_PLAN.json',dict(checkpoints=checkpoints,new_tasks=tasks,additional_candidates=decisions,reason='Reuse already completed matching fixed-action trajectories; compute only genuinely new windows',identity=identity()));rows+=batch(tasks,kinds=('LF',),workers=1)
 write(BASE/'LF_LADDER_RESULTS.json',dict(rows=rows,checkpoints=checkpoints,wall_s=time.perf_counter()-begin,identity=identity()))
 print('LF_LADDER_DONE',time.perf_counter()-begin,flush=True)
if __name__=='__main__':main()
