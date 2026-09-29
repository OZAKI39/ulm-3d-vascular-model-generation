"""Sparse independent realization batching: actual physics, unchanged clocks."""
from pathlib import Path
import argparse,json,sys,time,resource,subprocess
import numpy as np
import torch
ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True);args=ap.parse_args()
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(args.source_root/'particle_3d/src'))
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.gpu_engine import GPUPopulationEngine
from particle_3d.population_v1.gpu_field import GPUField
from particle_3d.population_v1.gpu_wall import GPUWall
from particle_3d.population_v1.gpu_batch import prefetch
from particle_3d.population_v1.provenance import write_json
from particle_3d.lammps_state import BridgeParticle
from particle_3d.particle_shapes import Sphere
class Null:
 def snapshot(self,*a):pass
 def terminal(self,*a):pass
 def _row(self,*a):pass
 def audit(self,*a):pass
env=load_environment(json.loads((root/'data/environment.json').read_text()));gf=GPUField(env.field,cache_path=root/'data/tetra_acceleration.npz');gw=GPUWall(env.wall)
data=np.load(root/'data/benchmark_input.npz');results=[];coordinates={}
for mode in ['sequential','batched']:
 engines=[]
 for rid in range(16):
  e=GPUPopulationEngine(env,pairs=True,realization_id=rid,gpu_field=gf,gpu_wall=gw,recorder=Null())
  x=data['x'][rid];v=data['field_values'][rid];p=BridgeParticle.from_shape(1,Sphere(x,1e-6),velocity=v[:3],omega=.5*v[13:16])
  e.insert_accepted(p,dict(particle_id=1,source_event_id=1,realization_id=rid,birth_time_s=0.,clearance_m=float(data['gap'][rid])))
  engines.append(e)
 torch.cuda.synchronize();start=time.perf_counter();cpu=time.process_time();prefix=0.
 for step in range(1,21):
  if mode=='batched':
   a=time.perf_counter();prefetch(engines);torch.cuda.synchronize();prefix+=time.perf_counter()-a
  for e in engines:e.advance(step*.0005)
 torch.cuda.synchronize();elapsed=time.perf_counter()-start;cpuused=time.process_time()-cpu
 coordinates[mode]=np.stack([e.pool.x[0].cpu().numpy() for e in engines])
 results.append(dict(mode=mode,realizations=16,particles_per_realization=1,nominal_steps_per_realization=20,
  total_particle_updates=sum(e.backend_counts['accepted_particle_updates'] for e in engines),wall_s=elapsed,
  particle_updates_s=sum(e.backend_counts['accepted_particle_updates'] for e in engines)/elapsed,
  CPU_utilization_percent_one_core_100=100*cpuused/elapsed,GPU_peak_allocated_bytes=torch.cuda.max_memory_allocated(),
  batched_prefix_wall_s=prefix,nvidia_smi_snapshot=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used','--format=csv,noheader,nounits'],text=True).strip(),
  IO='disabled in this compute benchmark; separate compressed-IO benchmark supplied',
  scope='16 independent clocks; physical FEM/wall/contact stages; no shared realization collisions'))
 print(results[-1],flush=True)
error=float(np.max(abs(coordinates['sequential']-coordinates['batched'])))
assert error<1e-13,error
write_json(root/'validation/phase9_sparse_batch_benchmark.json',dict(PASS=True,rows=results,
 maximum_position_difference_m=error,batched_vs_sequential_speedup=results[0]['wall_s']/results[1]['wall_s'],
 limitation='One timing per mode, small sparse workloads; launch and per-realization certification dominate.'))
