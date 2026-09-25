"""Cross-link raw solve, constraint and position correction logs by stage/ID."""
from pathlib import Path
import json,csv,numpy as np
from flow_geometry import Geometry
S=Path(__file__).resolve().parents[1];g=Geometry(S/'geometry/GEOMETRY_ARRAYS.npz');results=[]
def rows(p):return list(csv.DictReader(p.open()))
for case,D in [('LOW',S/'runs/LOW_release_mpi1'),('MEDIUM',S/'runs/MEDIUM_release_mpi1'),('LONG_PREFIX',S/'validation/abort_prefix_mpi1')]:
 stage={(int(r['step']),int(r['stage']),int(r['particle_id'])):r for r in rows(D/'INTEGRATION_STAGES.csv')};radii={int(r['particle_id']):float(r['radius_m']) for r in stage.values()};events=0
 for r in rows(D/'wall_constraint_timeseries.csv'):
  key=(int(r['step']),int(r['stage']),int(r['particle_id']));q=stage[key];a=float(q['radius_m'])
  assert all(float(r['v'+k+'_raw'])==float(q['v'+k+'_m_s']) for k in 'xyz')
  assert all(float(r['omega_raw_'+k])==float(q['omega_'+k]) for k in 'xyz')
  assert abs(float(r['gap_over_radius'])*a-float(r['gap']))<1e-20
  x=np.array([float(q[k+'_m']) for k in 'xyz']);assert all(x[i]==float(r[k+'_eval']) for i,k in enumerate('xyz'))
  n=np.array([float(r['n'+k]) for k in 'xyz']);raw=np.array([float(r['v'+k+'_raw']) for k in 'xyz']);used=np.array([float(r['v'+k+'_used']) for k in 'xyz'])
  assert abs(raw@n-float(r['vn_raw']))<1e-18 and abs(used@n-float(r['vn_used']))<1e-18
 for r in rows(D/'GEOMETRIC_PROJECTION_EVENTS.csv'):
  a=radii[int(r['particle_id'])];before=np.array([float(r['before_'+k]) for k in 'xyz']);after=np.array([float(r['after_'+k]) for k in 'xyz'])
  gb=g.distance(before)[0]-a;ga=g.distance(after)[0]-a
  assert abs(gb-float(r['gap_before']))<2e-14 and abs(ga-float(r['gap_after']))<2e-14
  assert gb<1e-10+2e-14
  # Correction targets the inherited safety surface, not an arbitrary clearance.
  if int(r['accepted']):assert abs(ga-1e-10)<2e-14
  assert 0<float(r['correction_distance'])<=2.5e-11
  events+=1
 results.append(dict(case=case,status='PASS',raw_stage_rows_linked=len(stage),correction_attempts_linked=events))
(S/'validation/LOG_LINKAGE_VALIDATION.json').write_text(json.dumps(dict(status='PASS',cases=results),indent=2)+'\n');print(results)
