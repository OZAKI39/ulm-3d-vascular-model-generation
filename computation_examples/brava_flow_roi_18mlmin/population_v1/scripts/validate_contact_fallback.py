"""Actual failed-state regression through the engine's effective CPU fallback."""
from pathlib import Path
import argparse,json,sys,gzip
import numpy as np
import torch
ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True);args=ap.parse_args()
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(args.source_root/'particle_3d/src'))
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.gpu_engine import GPUPopulationEngine
from particle_3d.population_v1.gpu_resistance import assemble_isolated,solve_isolated
from particle_3d.population_v1.provenance import write_json
env=load_environment(json.loads((root/'data/environment.json').read_text()));engine=GPUPopulationEngine(env,pairs=True,cache_path=root/'data/tetra_acceleration.npz')
probe=json.loads((root/'validation/actual_failed_state_parity.json').read_text());r=next(r for r in probe['results'] if r['particle_id']==1131)
with gzip.open(root/'runs/realization_0000_final/source.jsonl.gz','rt') as f:source=next(json.loads(l) for l in f if json.loads(l)['particle_id']==1131)
x=np.array(r['positions_m']);radius=source['radius_m'];n=len(x);p=engine.pool;p.n=n
p.particle_id[:n]=torch.arange(1,n+1,device='cuda');p.x[:n]=torch.tensor(x,device='cuda');p.radius[:n]=radius
sampled=engine.field.sample(p.x[:n]);wall=engine.wall.query(p.x[:n],p.radius[:n]);free=torch.cat([sampled[1][:,:3],.5*sampled[1][:,13:16]],dim=1)
matrix,rhs,diagonal,assembly=assemble_isolated(p.radius[:n],env.mu,wall[1],sampled[1][:,4:13].reshape(-1,3,3),free)
velocity,_,audit=solve_isolated(matrix,rhs,diagonal,free,wall[3],wall[4]);candidate=velocity[:,:3].cpu().numpy().copy();ref=np.array(r['CPU_velocity_m_s'])
assert bool(audit['fallback'].any()),audit['fallback']
for i in torch.nonzero(audit['fallback']).flatten().cpu().tolist():
 ix=torch.tensor([i],device='cuda');v,_=engine._cpu_components(ix,free);velocity[ix]=v
error=float(np.max(abs(velocity[:,:3].cpu().numpy()-ref)))
assert error<1e-14,error
write_json(root/'validation/actual_contact_fallback_regression.json',dict(PASS=True,particle_id=1131,actual_history_states=n,
 before_effective_GPU_vs_original_CPU_max_error_m_s=float(np.max(abs(candidate-ref))),after_effective_engine_vs_original_CPU_max_error_m_s=error,
 marked_for_original_CPU_fallback=int(audit['fallback'].sum()),policy='Unresolved dual-Gram singular values use preserved original positive-cone reduction; GPU candidate discarded',
 scientific_tolerance_not_relaxed=True,flow_sha256=env.flow_sha256))
print('ACTUAL_CONTACT_FALLBACK_REGRESSION',error,flush=True)
