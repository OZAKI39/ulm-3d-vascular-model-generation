"""Paired observational comparison: original 0.25 ms vs requested 1.0 ms, same births."""
from pathlib import Path
import csv,json,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[4];sys.path.insert(0,str(ROOT/'particle_3d/src'))
from particle_3d.formal_cohort_p9a5 import REL,DT,canonical,digest
R=ROOT/REL;old=ROOT/'particle_3d/reports/particle9a4_population_inlet/outputs/smoke30/trajectories'
rows=[]
for i in range(1,25):
 a=np.load(old/f'mb_{i:06d}.npz')['samples'];folder=R/'outputs/benchmark/workers_01/tracks'/f'mb_{i:06d}';b=np.load(folder/'trajectory.npz')['samples'];m=json.loads((folder/'metrics.json').read_text())
 old_m=json.loads((old/f'mb_{i:06d}.json').read_text())
 end=min(a[-1,0],b[-1,0]);t=np.linspace(0,end,501)
 xyz_a=np.column_stack([np.interp(t,a[:,0],a[:,j]) for j in [1,2,3]]);xyz_b=np.column_stack([np.interp(t,b[:,0],b[:,j]) for j in [1,2,3]])
 distance=np.linalg.norm(xyz_a-xyz_b,axis=1)
 old_o='O'+str(int(old_m['exit_outlet'].split('_')[-1])) if old_m.get('completed') else None
 rows.append(dict(particle_id=i,source_event_id=m['source_event_id'],old_dt_s=.00025,new_dt_s=DT,old_outlet=old_o,new_outlet=m['outlet'],old_end_reason=old_m['end_reason'],new_status=m['status'],old_age_s=float(a[-1,0]),new_age_s=float(b[-1,0]),age_change_s=float(b[-1,0]-a[-1,0]),shared_age_s=float(end),max_interpolated_position_difference_m=float(distance.max()),rms_interpolated_position_difference_m=float(np.sqrt(np.mean(distance**2))),old_samples_sha256=digest(old/f'mb_{i:06d}.npz'),new_samples_sha256=digest(folder/'trajectory.npz')))
summary=dict(old_dt_s=.00025,new_dt_s=DT,count=24,identical_frozen_births=True,outlet_changes=sum(r['old_outlet']!=r['new_outlet'] for r in rows),interpretation='Paired change audit only, not a time-convergence proof; positions linearly interpolated at 501 shared physical ages; no birth selected by outcome',rows=rows)
(R/'data/dt_comparison_24.json').write_bytes(canonical(summary))
with (R/'data/dt_comparison_24.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print(json.dumps({k:v for k,v in summary.items() if k!='rows'},indent=2))
