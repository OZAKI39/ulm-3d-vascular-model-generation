"""Explicit CPU pair-physics/contact fallback stage cost, no GPU claim."""
from pathlib import Path
import argparse,sys,time,json
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True);args=ap.parse_args()
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(args.source_root/'particle_3d/src'))
from particle_3d.particle_shapes import Sphere
from particle_3d.population_v1.components import assemble_components,solve_components
from particle_3d.population_v1.pairs import candidates,sphere_gap
from particle_3d.nearfield_handoff import handoff_constraints
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.particle6_stepper import bind_query_dependency
from particle_3d.population_v1.provenance import write_json
policy=NearFieldRegularizationV1();rows=[]
for n in [2,40]:
 shapes={i:Sphere([i*(2e-6+2e-9),0,0],1e-6) for i in range(n)}
 free={i:np.array([.001*(n/2-i)/n,2e-4,0,0,0,0]) for i in shapes}
 times=[]
 for repeat in range(6):
  a=time.perf_counter();system=assemble_components(shapes,free,.00345312,None,policy=policy,gradient_at=lambda x:np.zeros((3,3)))
  b=time.perf_counter();constraints,_=bind_query_dependency(handoff_constraints,dict(all_pairs=candidates,pair_gap=sphere_gap))(shapes,None,policy)
  c=time.perf_counter();result=solve_components(system,particles=shapes,constraints=constraints);d=time.perf_counter()
  if repeat:times.append([b-a,c-b,d-c])
  assert all(c['max_constraint_violation_m_s']<1e-10 for c in result.record['components'])
 med=np.median(times,axis=0);rows.append(dict(particles=n,contacts=len(constraints),backend='CPU_EXPLICIT_COMPONENT_FALLBACK',
  pair_resistance_assembly_s=float(med[0]),contact_geometry_s=float(med[1]),coupled_contact_solve_s=float(med[2]),
  total_stage_median_sum_s=float(sum(med)),replicates=5,wall='none; isolate pair mechanism',GPU_pair_compute_claimed=False))
write_json(root/'validation/phase9_pair_component_cost.json',dict(rows=rows,scope='Separate CPU timing; production remains hybrid',GPU_concurrent_production=True,
 limitation='CPU host contention may affect timing; not a full production speedup benchmark'))
print(json.dumps(rows),flush=True)
