"""Read-only audit of the completed P9-A.4 cohort; no trajectory advancement."""
from pathlib import Path
import importlib.util,json,sys,time
from collections import Counter
import numpy as np
ROOT=Path(__file__).resolve().parents[4]
R=ROOT/'particle_3d/reports/particle9a4_population_inlet';D=R/'data'
sys.path.insert(0,str(ROOT/'particle_3d/src'))
sys.path.insert(0,str(ROOT/'particle_3d/reports/network_derived_flow_mb_validation_v1/scripts'))
from particle_3d.population_inlet_p9a4 import load_new_environment
from particle_3d.particle82a_geometry import InletGeometryAudit,lower_gap
from particle_3d.validation_boundary import ValidationBoundaryClassifier
import analyze as original

def read(p):return json.loads(p.read_text())
def write(p,value):
 with p.open('x') as f:json.dump(value,f,indent=2,default=original.safe,allow_nan=False);f.write('\n')

def main():
 start=time.perf_counter();env=load_new_environment(ROOT);original.R=R
 geo=InletGeometryAudit(env);inlet=ValidationBoundaryClassifier({'INLET':env.boundaries['INLET']})
 events=read(D/'inlet100k/smoke30_births.json')['events'];rows=[];audits=[]
 for e in events:
  row,audit=original.one('smoke30',e,env)
  a=np.load(R/'outputs/smoke30/trajectories'/f"mb_{e['particle_id']:06d}.npz")['samples']
  gaps=np.array([env.wall.nearest_center_triangle(p)[1]-e['radius_m'] for p in a[:,1:4]])
  row.update(independent_wall_penetration_count=int((gaps < -env.wall.roundoff_m).sum()),
   independent_handoff_violation_count=int((gaps-lower_gap(e['radius_m']) < -env.wall.roundoff_m).sum()),
   independent_minimum_wall_gap_m=float(gaps.min()),saved_sample_count=len(a))
  escapes=[]
  for i,(x,y) in enumerate(zip(a[:-1,1:4],a[1:,1:4])):
   # A birth on the open cap moving inward is not an escape. Test all
   # outward-moving segments against the actual (slightly nonplanar) cap.
   if (y-x)@geo.normal >= 0:continue
   hit=inlet.first_event(x,y)
   if hit is not None:escapes.append(dict(segment=i,position_m=hit.position_m,triangle_id=hit.triangle_id))
  row['inlet_escape_count']=len(escapes);audit['all_saved_segments_outward_inlet_crossings']=escapes
  rows.append(row);audits.append(audit)
  print(e['particle_id'],row['status'],'stationary_supported',row['stationary_three_contact_supported'],flush=True)
 counts=Counter(r['status'] for r in rows)
 keys=['penetration_count','handoff_violation_count','continuous_certificate_violations','nan_inf_count',
  'accepted_diagnostic_nan_inf_count','unclassified_solver_corruption','independent_wall_penetration_count',
  'independent_handoff_violation_count','inlet_escape_count']
 safety={k:sum(r[k] for r in rows) for k in keys}
 point=read(R/'outputs/smoke30/point_completed.json')['rows']
 point_counts={o:sum(p['outlet']=='OUTLET_'+o[-1].zfill(2) for p in point) for o in ['O1','O2','O3']}
 point_counts['noexit']=sum(p['outlet'] is None for p in point)
 unsupported=[r['bubble_id'] for r in rows if r['status']=='STATIONARY' and not r['stationary_three_contact_supported']]
 passed=not any(safety.values()) and not unsupported and all(r['status'] in ['COMPLETED','STATIONARY'] for r in rows)
 summary=dict(N=30,status_counts=dict(counts),completed=counts['COMPLETED'],stationary=counts['STATIONARY'],solver_fail=counts['SOLVER_FAIL'],
  stationary_ids=[r['bubble_id'] for r in rows if r['status']=='STATIONARY'],unsupported_stationary_ids=unsupported,
  safety=safety,outlet_counts={o:sum(r['outlet']==o for r in rows) for o in ['O1','O2','O3']},point_basin_counts=point_counts,
  point_noexit_ids=[p['bubble_id'] for p in point if p['outlet'] is None],
  continuous_certificate_count=sum(a['continuous_certificate_count'] for a in audits),
  saved_sample_count=sum(r['saved_sample_count'] for r in rows),
  minimum_gap_m=min(r['independent_minimum_wall_gap_m'] for r in rows),
  minimum_g_nf_m=min(r['minimum_g_nf_m'] for r in rows),
  smoke_runtime_s=read(R/'outputs/smoke30/completed.json')['runtime_s'],
  analysis_runtime_s=time.perf_counter()-start,passed=bool(passed),
  scope='Fixed first 30 accepted births; stationary means supported three-contact constrained equilibrium, not physiological trapping proof; point basins descriptive only')
 write(D/'smoke30_metrics.json',rows);write(D/'smoke30_contact_audit.json',audits);write(D/'smoke30_summary.json',summary)
 print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
