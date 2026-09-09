"""Native Mirheo SDPD worker. API adapted from fluid_comparison/gpu_worker.py.

Mirheo API examples Copyright ETH Zurich, MIT; this file is not HemoCell code.
One coordinator per rank, one initialization, no restart or thermal rescaling.
"""
import argparse
import csv
import ctypes
import json
import math
import os
from pathlib import Path
import sys
import time
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--spec',required=True)
    args=parser.parse_args()  # --help must not import CUDA/MPI.
    s=json.loads(Path(args.spec).read_text())
    if os.environ.get('OMPI_COMM_WORLD_SIZE')!='2' or os.environ.get('SOLVER_BENCHMARK_AUTHORIZED')!=s['plan_sha256']:
        raise RuntimeError('AUTHORIZED_LAUNCH_REQUIRED')
    sys.path.insert(0,s['project_root'])
    from py_scripts.fluid_comparison.gpu_worker import thermal_bins
    from native_mpi import NativeMPI as MPI
    import mirheo as mir
    comm=MPI();compute=comm.rank==0
    start=time.perf_counter();d=Path.cwd();p=s['physics'];c=p['candidate'];m=p['mapping'];dt=c['dt_star']
    runtime=ctypes.CDLL('libcudart.so.12') if compute else None
    if compute:
        runtime.cudaMemcpy.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_size_t,ctypes.c_int]
        runtime.cudaMemcpy.restype=ctypes.c_int
    def sync():
        if compute and runtime.cudaDeviceSynchronize()!=0:raise RuntimeError('CUDA_SYNC_FAILED')
    def density_channel(pv):
        a=pv.local.per_particle['saved_density'].__cuda_array_interface__
        out=np.empty(a['shape'],dtype=np.dtype(a['typestr']))
        if a.get('strides') not in (None,out.strides):raise RuntimeError('DENSITY_STRIDE_MISMATCH')
        if runtime.cudaMemcpy(out.ctypes.data,int(a['data'][0]),out.nbytes,2)!=0:raise RuntimeError('CUDA_COPY_FAILED')
        return out.reshape(-1).astype(float)
    u=mir.Mirheo((1,1,1),tuple(p['domain_star']),debug_level=1,log_filename='mirheo',no_splash=True)
    pv=mir.ParticleVectors.ParticleVector('pv',mass=m['m_star'])
    u.registerParticleVector(pv,mir.InitialConditions.Uniform(number_density=m['n_star']))
    if compute:
        N=len(pv.getCoordinates())
        v=np.random.default_rng(20260908).normal(0,math.sqrt(m['kBT_star']/m['m_star']),(N,3))
        v-=v.mean(axis=0);pv.setVelocities(v.tolist())
    den=mir.Interactions.Pairwise('density',c['rc_star'],kind='Density',density_kernel='WendlandC2')
    interaction=mir.Interactions.Pairwise('sdpd',c['rc_star'],kind='SDPD',viscosity=c['viscosity_mu_star'],kBT=m['kBT_star'],
        EOS=c['EOS'],sound_speed=c['sound_speed_star'],rho_0=c['rho_0_star'],density_kernel='WendlandC2',stress=False)
    for inter in (den,interaction):u.registerInteraction(inter);u.setInteraction(inter,pv,pv)
    rest=mir.Integrators.VelocityVerlet('rest')
    flow=mir.Integrators.VelocityVerlet_withPeriodicForce('flow',force=p['particle_force_star'],direction='x')
    u.registerIntegrator(rest);u.registerIntegrator(flow);u.setIntegrator(rest,pv)
    u.registerPlugins(mir.Plugins.createParticleChecker('checker',check_every=s['sample_steps']))
    u.registerPlugins(mir.Plugins.createParticleChannelSaver('save_density',pv,'densities','saved_density'))
    streams=[]
    def writer(name,header):
        f=(d/name).open('x',buffering=1);streams.append(f);w=csv.writer(f);w.writerow(header);return w
    if compute:
        profiles=writer('profiles.csv',['step','time_si','bin','y_si','count','ux_si','uy_si','uz_si','rho_si'])
        moments=writer('moments.csv',['step','time_si','N','temperature_K','temperature_coarse_K','COM_temperature_K','kernel_rho_mean_si','kernel_rho_cv','max_speed_si'])
        timings=writer('timings.csv',['step','time_si','chunk_steps','compute_max_rank_s','sampling_max_rank_s','compute_cumulative_s','sampling_cumulative_s','elapsed_worker_s'])
    velocity_scale=m['si_per_star']['velocity'];mass_scale=m['si_per_star']['mass_density'];T=p['temperature_K']/m['kBT_star']
    def sample(step):
        if not compute:return
        pos=np.asarray(pv.getCoordinates(),float);vel=np.asarray(pv.getVelocities(),float)
        if not np.isfinite(pos).all() or not np.isfinite(vel).all():raise RuntimeError('NONFINITE_STATE')
        thermal,counts,av,v2=thermal_bins(pos,vel,p['domain_star'][1],p['bins'],m['m_star'])
        fine=thermal_bins(pos,vel,p['domain_star'][1],p['bins']*2,m['m_star'])[0]
        for b in range(p['bins']):
            profiles.writerow([step,step*s['dt_si'],b,(b+.5)*p['bin_width_si'],int(counts[b]),*(av[b]*velocity_scale),counts[b]*m['m_star']*m['M0']/(np.prod(p['box_si'])/p['bins'])])
        com=m['m_star']*np.sum((vel-vel.mean(axis=0))**2)/(3*(len(vel)-1))
        kernel=density_channel(pv)*m['m_star']*mass_scale if step else None
        moments.writerow([step,step*s['dt_si'],len(vel),fine*T,thermal*T,com*T,
                          float(kernel.mean()) if kernel is not None else '',float(kernel.std()/kernel.mean()) if kernel is not None else '',float(np.linalg.norm(vel,axis=1).max()*velocity_scale)])
    sync();setup=comm.allreduce(time.perf_counter()-start,op=MPI.MAX)
    sample(0);done=0;compute_s=sample_s=0.
    while done<s['steps']:
        stop=comm.bcast((d/'STOP_REQUESTED').exists() if compute else None,root=0)
        if stop:break
        if done==s['warmup_steps']:
            u.deregisterIntegrator(rest)
            u.setIntegrator(flow,pv)
        n=min(s['sample_steps'],s['steps']-done)
        if done<s['warmup_steps']:n=min(n,s['warmup_steps']-done)
        comm.Barrier();before=time.perf_counter();u.run(n,dt=dt);sync()
        comp=comm.allreduce(time.perf_counter()-before,op=MPI.MAX);compute_s+=comp;done+=n
        before=time.perf_counter();sample(done);samp=comm.allreduce(time.perf_counter()-before,op=MPI.MAX);sample_s+=samp
        if compute:timings.writerow([done,done*s['dt_si'],n,comp,samp,compute_s,sample_s,time.perf_counter()-start])
    if compute:
        for f in streams:f.close()
        (d/'completion.json').write_text(json.dumps(dict(actual_steps=done,actual_time_si=done*s['dt_si'],dt_si=s['dt_si'],rank_count=2,compute_ranks=1,
            postprocess_ranks=1,fluid_particles=N,wall_particles=0,cell_count=0,setup_s=setup,geometry_s=None,compute_s=compute_s,sampling_s=sample_s,
            worker_elapsed_s=time.perf_counter()-start,completed=done==s['steps'],coordinator_count_per_rank=1,velocity_initializations=1,
            stress_output=False,pressure_output=False,density_phase='pre-integration, one dt before sampled end velocity',
            synchronization='cudaDeviceSynchronize then MPI_MAX across both ranks'),allow_nan=False,indent=2))
    del u


if __name__=='__main__':main()
