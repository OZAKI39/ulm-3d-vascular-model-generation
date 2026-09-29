"""Separated stage timings and actual certified particle updates on real BraVa.

These synthetic test populations are not the physical infusion realization.
They deliberately select noninteracting interior states so CPU/GPU execute an
identical full 0.5 ms interval without unrelated shared-clock handoff events.
"""
from pathlib import Path
import argparse,json,sys,time,resource,subprocess,csv
import numpy as np
ap=argparse.ArgumentParser();ap.add_argument('--source-root',type=Path,required=True)
ap.add_argument('--mode',choices=['gpu','cpu'],required=True);args=ap.parse_args()
sys.path.insert(0,str(args.source_root/'particle_3d/src'))
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.provenance import write_json
from particle_3d.population_v1.engine import PopulationEngine
from particle_3d.lammps_state import BridgeParticle
from particle_3d.particle_shapes import Sphere
root=Path(__file__).resolve().parents[1];env=load_environment(json.loads((root/'data/environment.json').read_text()))
counts=[1,10,100,1000,10000];dt=.0005

class NullOutput:
    def snapshot(self,r):pass
    def terminal(self,r):pass
    def audit(self,r):pass
    def _row(self,*a):pass

def usage():
    r=resource.getrusage(resource.RUSAGE_SELF);return r.ru_utime+r.ru_stime

if args.mode=='gpu':
    import torch
    from particle_3d.population_v1.gpu_field import GPUField
    from particle_3d.population_v1.gpu_wall import GPUWall
    from particle_3d.population_v1.gpu_neighbors import neighbors
    from particle_3d.population_v1.gpu_resistance import assemble_isolated,solve_isolated
    from particle_3d.population_v1.gpu_engine import GPUPopulationEngine
    from particle_3d.population_v1.gpu_state import GPUActivePool
    from particle_3d.population_v1.output import PopulationOutput
    gf=GPUField(env.field,cache_path=root/'data/tetra_acceleration.npz');gw=GPUWall(env.wall)
    rng=np.random.default_rng(2026092991);accepted=[]
    while sum(len(r) for r in accepted)<12000:
        ids=rng.integers(0,len(env.field.tetra),50000);w=rng.dirichlet(np.ones(4),len(ids))
        x=np.einsum('ni,nij->nj',w,env.field.points_m[env.field.tetra[ids]])
        xx=torch.tensor(x,device='cuda');r=torch.full((len(x),),1e-6,device='cuda',dtype=torch.float64)
        f=gf.sample(xx);g=gw.query(xx,r)
        clear=(g[1][:,0]>1.1*dt*torch.linalg.vector_norm(f[1][:,:3],dim=1)+2*r)&(f[0]>=0)
        accepted.append(x[clear.cpu().numpy()]);print('BENCHMARK_INTERIOR_STATES',sum(len(r) for r in accepted),flush=True)
    x=np.concatenate(accepted)[:12000];xx=torch.tensor(x,device='cuda');r=torch.full((len(x),),1e-6,device='cuda',dtype=torch.float64)
    f=gf.sample(xx);v=f[1][:,:3];g=gw.query(xx,r);real=torch.zeros(len(x),device='cuda',dtype=torch.int64)
    pairs,geom,_=neighbors(xx,r,real,delta=dt*v)
    keep=torch.ones(len(x),device='cuda',dtype=torch.bool)
    if len(pairs):
        i,j=pairs[:,0],pairs[:,1];rel=xx[i]-xx[j];dv=dt*(v[i]-v[j]);den=(dv*dv).sum(1)
        t=torch.where(den>0,-(rel*dv).sum(1)/den,torch.zeros_like(den)).clamp(0.,1.)
        h=torch.linalg.vector_norm(rel+t[:,None]*dv,dim=1)-2e-6
        bad=pairs[h<.05e-6];keep[bad.flatten()]=False
    # Keep each segment away from all open-cap AABBs too.
    end=xx+dt*v;lo=torch.minimum(xx,end);hi=torch.maximum(xx,end)
    for role in sorted(set(env.classifier.roles)):
        tri=env.classifier.xyz[env.classifier.roles==role]
        bl=torch.tensor(tri.min((0,1)),device='cuda');bh=torch.tensor(tri.max((0,1)),device='cuda')
        keep&=~((lo<=bh)&(hi>=bl)).all(1)
    ix=torch.nonzero(keep).flatten()[:10000]
    if len(ix)<10000:raise ValueError('Insufficient collision-free benchmark interior states')
    x=xx[ix].cpu().numpy();values=f[1][ix].cpu().numpy();gap=g[1][ix,0].cpu().numpy()
    np.savez_compressed(root/'data/benchmark_input.npz',x=x,field_values=values,gap=gap,radius=np.full(len(x),1e-6))
    def timed(fn):
        times=[]
        for _ in range(4):
            a=torch.cuda.Event(enable_timing=True);b=torch.cuda.Event(enable_timing=True)
            a.record();value=fn();b.record();b.synchronize();times.append(a.elapsed_time(b)/1000)
        return float(np.median(times[1:])),value
else:
    data=np.load(root/'data/benchmark_input.npz');x=data['x'];values=data['field_values'];gap=data['gap']

rows=[]
for n in counts:
    print('BENCHMARK_START',args.mode,n,flush=True)
    if args.mode=='gpu':
        tx=torch.tensor(x[:n],device='cuda');radius=torch.full((n,),1e-6,device='cuda',dtype=torch.float64)
        real=torch.zeros(n,device='cuda',dtype=torch.int64)
        ts,field=timed(lambda:gf.sample(tx));tw,wall=timed(lambda:gw.query(tx,radius))
        tn,nb=timed(lambda:neighbors(tx,radius,real))
        free=torch.cat([field[1][:,:3],.5*field[1][:,13:16]],dim=1)
        ta,assembled=timed(lambda:assemble_isolated(radius,env.mu,wall[1],field[1][:,4:13].reshape(-1,3,3),free))
        matrix,rhs,diagonal,_=assembled
        tc,solved=timed(lambda:solve_isolated(matrix,rhs,diagonal,free,wall[3],wall[4]))
        tp,proof=timed(lambda:gw.path(tx,solved[0][:,:3],radius,dt))
        engine=GPUPopulationEngine(env,pairs=True,recorder=NullOutput(),gpu_field=gf,gpu_wall=gw)
        p=GPUActivePool(n);p.n=n;p.used_ids={(0,i+1) for i in range(n)}
        p.particle_id[:n]=torch.arange(1,n+1,device='cuda');p.source_event_id[:n]=p.particle_id[:n]
        p.x[:n]=tx;p.radius[:n]=radius;p.diameter[:n]=2*radius;p.q[:n,0]=1.
        p.tetra_id[:n]=field[0]
        p.v[:n]=field[1][:,:3];p.omega[:n]=.5*field[1][:,13:16];p.wall_gap[:n]=wall[1][:,0]
        engine.pool=p;engine.safety._resize(n);engine.safety.n=n
        engine.safety.arrays['count'][:n]=1;engine.safety.arrays['history'][:n]=tx[:,None,:]
        engine.safety.arrays['previous_nominal'][:n]=tx
        engine.counts['accepted_birth_count']=n;engine.maximum_active=n
        engine.births={i+1:dict(source_event_id=i+1,birth_time_s=0.) for i in range(n)}
        torch.cuda.synchronize();start=time.perf_counter();cpu_start=usage();engine.advance(dt);torch.cuda.synchronize();elapsed=time.perf_counter()-start
        updates=engine.backend_counts['accepted_particle_updates']
        gpu_memory=torch.cuda.max_memory_allocated();query=subprocess.check_output(['nvidia-smi','--query-gpu=utilization.gpu,memory.used','--format=csv,noheader,nounits'],text=True).strip()
        row=dict(N_active=n,backend='CUDA_HYBRID',wall_s=elapsed,cpu_s=usage()-cpu_start,
            accepted_particle_updates=updates,steps_s=engine.backend_counts['accepted_intervals']/elapsed,particle_updates_s=updates/elapsed,
            field_sampling_s=ts,wall_query_s=tw,neighbor_build_and_exact_gap_s=tn,
            resistance_assembly_s=ta,contact_solve_s=tc,continuous_wall_certificate_s=tp,
            pair_candidate_count=len(nb[0]),gpu_peak_allocated_bytes=gpu_memory,nvidia_smi_snapshot=query,
            backend_counts=dict(engine.backend_counts))
        out=root/'validation_runs'/f'benchmark_io_N{n}'
        writer=PopulationOutput(out,dict(role='BENCHMARK_IO_ONLY',flow_sha256=env.flow_sha256))
        start=time.perf_counter();engine.recorder=writer;engine.flush_output();writer.close(dict(benchmark=True))
        row['batched_state_transfer_and_compressed_IO_s']=time.perf_counter()-start
        row['CPU_utilization_percent_one_core_100']=100*row['cpu_s']/elapsed
        del engine
    else:
        engine=PopulationEngine(env,pairs=True,recorder=NullOutput());engine.record=lambda:None
        for i in range(n):
            p=BridgeParticle.from_shape(i+1,Sphere(x[i],1e-6),velocity=values[i,:3],omega=.5*values[i,13:16])
            engine.insert_accepted(p,dict(particle_id=i+1,source_event_id=i+1,realization_id=0,birth_time_s=0.,clearance_m=gap[i]))
        start=time.perf_counter();cpu_start=usage();engine.advance(dt);elapsed=time.perf_counter()-start
        row=dict(N_active=n,backend='CPU_COMPONENT_REFERENCE',wall_s=elapsed,cpu_s=usage()-cpu_start,
            accepted_particle_updates=n,steps_s=1/elapsed,particle_updates_s=n/elapsed,
            peak_rss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            CPU_utilization_percent_one_core_100=100*(usage()-cpu_start)/elapsed)
        if engine.pool.n!=n:raise AssertionError('Benchmark event unexpectedly removed a particle')
    rows.append(row);print(json.dumps(row),flush=True)
    write_json(root/'validation'/f'phase9_benchmark_{args.mode}.json',dict(phase=9,flow_sha256=env.flow_sha256,
        role='SYNTHETIC_NONINTERACTING_BENCHMARK_STATES_IN_ACTUAL_FLOW_NOT_PHYSICAL_INFUSION',nominal_dt_s=dt,rows=rows))
