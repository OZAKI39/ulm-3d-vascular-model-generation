"""Read-only interim failure evidence; never changes production state."""
from pathlib import Path
import sys,json,zlib,argparse
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True);ap.add_argument('--run-name',required=True);args=ap.parse_args()
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(args.source_root/'particle_3d/src'))
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.provenance import write_json
from particle_3d.wall_gap import wall_gap
from particle_3d.particle_shapes import Sphere
run=root/'runs'/args.run_name
def complete_rows(name):
 text=zlib.decompressobj(31).decompress((run/(name+'.jsonl.gz')).read_bytes()).decode();rows=[]
 for line in text.splitlines():
  try:rows.append(json.loads(line))
  except json.JSONDecodeError:pass
 return rows
term=[r for r in complete_rows('terminals') if r['status']=='SOLVER_FAILURE'];source={r['particle_id']:r for r in complete_rows('source') if r['particle_id'] is not None}
env=load_environment(json.loads((root/'data/environment.json').read_text()));rows=[]
for r in term:
 e=r.get('support_evidence',[])
 if not e:continue
 x=np.array(e[-1]['position_m']);radius=source[r['particle_id']]['radius_m'];sample=env.field.sample(x);tet,_=env.field.locate(x);g=wall_gap(Sphere(x,radius),env.wall)
 rows.append(dict(particle_id=r['particle_id'],trigger=r['trigger'],physical_time_s=r['physical_time_s'],position_m=x.tolist(),tetra_id=tet,
  FEM_velocity_m_s=sample.velocity_m_s.tolist(),tetra_vertex_speeds_m_s=np.linalg.norm(env.field.velocity_nodes_m_s[env.field.tetra[tet]],axis=1).tolist(),
  nearest_wall_triangle_id=g.wall_triangle_id,wall_gap_m=g.gap_m,contact_rank=e[-1]['solver'].get('contact_redundancy',{}).get('rank_after',0),
  maximum_GPU_CPU_velocity_difference_m_s=max(float(np.max(abs(np.array(v['velocity'])-v['independently_recomputed_velocity_m_s']))) for v in e),support_verified=r['support_verified']))
write_json(root/'validation/live_failure_probe.json',dict(scope='INTERIM_SUBSET_NOT_FINAL_RESULT',flow_sha256=env.flow_sha256,run=args.run_name,rows=rows))
print(json.dumps(rows,default=lambda x:x.item()),flush=True)
