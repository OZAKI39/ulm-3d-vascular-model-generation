"""Extra parity probe at actual failed-state history, distinct from unit fixtures."""
from pathlib import Path
import argparse,json,zlib,sys
import numpy as np
import torch
ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True);ap.add_argument('--run-name',required=True);args=ap.parse_args()
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(args.source_root/'particle_3d/src'))
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.gpu_field import GPUField
from particle_3d.population_v1.gpu_wall import GPUWall
from particle_3d.population_v1.gpu_resistance import assemble_isolated,solve_isolated
from particle_3d.population_v1.provenance import write_json
from particle_3d.particle_shapes import Sphere
from particle_3d.wall_gap import wall_gap
run=root/'runs'/args.run_name
def read(name):
 text=zlib.decompressobj(31).decompress((run/(name+'.jsonl.gz')).read_bytes()).decode();rows=[]
 for l in text.splitlines():
  try:rows.append(json.loads(l))
  except json.JSONDecodeError:pass
 return rows
source={r['particle_id']:r for r in read('source') if r['particle_id'] is not None};terms=[r for r in read('terminals') if r['status']=='SOLVER_FAILURE' and r.get('support_evidence')]
env=load_environment(json.loads((root/'data/environment.json').read_text()));gf=GPUField(env.field,cache_path=root/'data/tetra_acceleration.npz');gw=GPUWall(env.wall);result=[]
for term in terms:
 pid=term['particle_id'];e=term['support_evidence'];x=np.array([r['position_m'] for r in e]);radius=source[pid]['radius_m'];tx=torch.tensor(x,device='cuda');tr=torch.full((len(x),),radius,device='cuda',dtype=torch.float64)
 f=gf.sample(tx);w=gw.query(tx,tr);free=torch.cat([f[1][:,:3],.5*f[1][:,13:16]],dim=1)
 matrix,rhs,diagonal,assembly=assemble_isolated(tr,env.mu,w[1],f[1][:,4:13].reshape(-1,3,3),free)
 velocity,hydro,audit=solve_isolated(matrix,rhs,diagonal,free,w[3],w[4])
 cpu=[env.field.sample(p) for p in x];cpuwall=[wall_gap(Sphere(p,radius),env.wall) for p in x]
 ref=np.array([r['independently_recomputed_velocity_m_s'] for r in e]);saved=np.array([r['velocity'] for r in e]);new=velocity[:,:3].cpu().numpy()
 record=dict(particle_id=pid,n=len(e),
  field_velocity_max_error_m_s=float(np.max(abs(np.array([s.velocity_m_s for s in cpu])-f[1][:,:3].cpu().numpy()))),
  gradient_max_error_s_inv=float(np.max(abs(np.array([s.velocity_gradient_s_inv for s in cpu])-f[1][:,4:13].cpu().numpy().reshape(-1,3,3)))),
  wall_gap_max_error_m=float(np.max(abs(np.array([s.gap_m for s in cpuwall])-w[1][:,0].cpu().numpy()))),
  wall_nearest_triangle_mismatches=int(np.count_nonzero(np.array([s.wall_triangle_id for s in cpuwall])!=w[0].cpu().numpy())),
  fresh_GPU_vs_CPU_velocity_max_error_m_s=float(np.max(abs(new-ref))),fresh_GPU_vs_recorded_GPU_max_error_m_s=float(np.max(abs(new-saved))),
  recorded_GPU_vs_CPU_max_error_m_s=float(np.max(abs(saved-ref))),
  original_CPU_contact_counts=[s['solver']['contact_count'] for s in e],fresh_GPU_contact_counts=w[4].cpu().tolist(),
  wall_roundoff_m=env.wall.roundoff_m,positions_m=x.tolist(),fresh_GPU_velocity_m_s=new.tolist(),CPU_velocity_m_s=ref.tolist(),recorded_GPU_velocity_m_s=saved.tolist())
 result.append(record);print({k:v for k,v in record.items() if k not in ['positions_m','fresh_GPU_velocity_m_s','CPU_velocity_m_s','recorded_GPU_velocity_m_s']},flush=True)
write_json(root/'validation/actual_failed_state_parity.json',dict(scope='ACTUAL_FAILED_STATE_INTERIM_PROBE_NOT_GENERIC_PARITY_PASS',flow_sha256=env.flow_sha256,results=result))
