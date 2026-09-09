"""Small native liquid experiment; adapted from fluid_physics/gpu_worker.py.

Mirheo Density/SDPD API wiring follows tests/sdpd/rest.py and double_poiseuille.py,
Copyright ETH Zurich, MIT (MIRHEO_LICENSE.txt saved in every attempt).
No CPU integrator, alternate thermostat, or simultaneous DPD/SDPD force pair.
"""
import argparse
import csv
import ctypes
import json
import math
import os
from pathlib import Path
import time
import numpy as np
try:
    from .observations import cuda_layout,phase_audit,velocity_statistics
except ImportError:  # Frozen standalone copies in historical diagnostic attempts.
    from observations import cuda_layout,phase_audit,velocity_statistics


def continuous_chunks(u,d,compute,steps,every,dt,sync):
    """Keep one coordinator and its native state; both MPI ranks read one sequence."""
    if steps<=0 or every<=0 or dt<=0:raise ValueError('INVALID_CONTINUOUS_PLAN')
    done=0;index=0
    while True:
        control=d/'control'/f'chunk_{index:07d}.json'
        if compute:
            n=0 if (d/'STOP_REQUESTED').exists() else min(every,steps-done)
            tmp=control.with_suffix('.tmp');tmp.write_text(json.dumps({'steps':n}));tmp.rename(control)
        else:
            while not control.exists():time.sleep(.001)
        n=json.loads(control.read_text())['steps']
        if not n:break
        start=time.monotonic();u.run(n,dt=dt);sync();done+=n
        yield index,done,time.monotonic()-start
        index+=1


def thermal_bins(pos,vel,Ly,bins,mass):
    ids=np.minimum(bins-1,np.floor(np.mod(pos[:,1],Ly)/Ly*bins).astype(int))
    counts=np.bincount(ids,minlength=bins)
    sums=np.column_stack([np.bincount(ids,weights=vel[:,j],minlength=bins) for j in range(3)])
    mean=np.divide(sums,counts[:,None],out=np.zeros_like(sums),where=counts[:,None]>0)
    v2=np.bincount(ids,weights=np.sum(vel**2,axis=1),minlength=bins)
    kBT=mass*np.sum((vel-mean[ids])**2)/(3*(len(vel)-np.count_nonzero(counts)))
    return float(kBT),counts,mean,v2


def main(observer_factory=None):
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--spec',required=True);args=parser.parse_args()
    if os.environ.get('OMPI_COMM_WORLD_SIZE')!='2':raise RuntimeError('REQUIRES_TWO_RANK_AUTHORIZED_LAUNCHER')
    setup_started=time.monotonic();spec=json.loads(Path(args.spec).read_text());d=Path.cwd()
    rank=int(os.environ['OMPI_COMM_WORLD_RANK']);compute=rank==0
    import mirheo as mir
    runtime=ctypes.CDLL('libcudart.so.12') if compute else None
    if compute:
        runtime.cudaMemcpy.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_size_t,ctypes.c_int]
        runtime.cudaMemcpy.restype=ctypes.c_int
    def sync():
        if compute and runtime.cudaDeviceSynchronize()!=0:raise RuntimeError('CUDA_SYNC_FAILED')
    def channel(pv,name):
        a=pv.local.per_particle[name].__cuda_array_interface__;shape,dtype,strides,size=cuda_layout(a)
        if size==0:return np.empty(shape,dtype=dtype)
        raw=np.empty(size,dtype=np.uint8)
        if runtime.cudaMemcpy(raw.ctypes.data,int(a['data'][0]),size,2)!=0:raise RuntimeError('CUDA_CHANNEL_COPY_FAILED '+name)
        return np.ndarray(shape,dtype=dtype,buffer=raw,strides=tuple(strides)).copy()
    cand=spec['candidate'];task=spec['task'];domain=tuple(spec['domain_star']);dt=spec['dt_star'];mass=spec['m_star']
    u=mir.Mirheo((1,1,1),domain,debug_level=1,log_filename='mirheo',no_splash=True)
    pv=mir.ParticleVectors.ParticleVector('pv',mass=mass)
    u.registerParticleVector(pv,mir.InitialConditions.Uniform(number_density=task['n_star']))
    if compute:
        N=len(pv.getCoordinates());v=np.random.default_rng(spec['initial_velocity_seed']).normal(0,math.sqrt(spec['kBT_star']/mass),(N,3))
        v-=v.mean(axis=0);pv.setVelocities(v.tolist())
        sync();initial_v=np.asarray(pv.getVelocities(),float);initial_x=np.asarray(pv.getCoordinates(),float)
        np.save(d/'initial_positions_global.npy',initial_x,allow_pickle=False)
        np.save(d/'initial_velocities.npy',initial_v,allow_pickle=False)
        (d/'initial_state.json').write_text(json.dumps({'time_star':0.,'phase':'after Uniform and Gaussian COM subtraction, before first force',
            'Uniform_seed_rule':'rank + std::hash(pv name), pinned library','velocity_seed':spec['initial_velocity_seed'],
            **velocity_statistics(initial_v,mass)},indent=2,allow_nan=False)+'\n')
    observer=observer_factory(spec,d,pv,channel) if compute and observer_factory else None
    density_enabled=cand['method']=='SDPD' or task.get('passive_density',False)
    if density_enabled:
        den=mir.Interactions.Pairwise('density',spec['rc_star'],kind='Density',density_kernel='WendlandC2')
        u.registerInteraction(den);u.setInteraction(den,pv,pv)
    if cand['method']=='SDPD':
        if not cand.get('density_interaction'):raise RuntimeError('SDPD_REQUIRES_DENSITY')
        interaction=mir.Interactions.Pairwise('sdpd',spec['rc_star'],kind='SDPD',viscosity=cand['viscosity_mu_star'],
            kBT=spec['kBT_star'],EOS='Linear',sound_speed=cand['sound_speed_star'],rho_0=cand['rho_0_star'],
            density_kernel='WendlandC2',stress=True,stress_period=dt/2)
    else:
        interaction=mir.Interactions.Pairwise('dpd',spec['rc_star'],kind='DPD',a=cand['a'],gamma=cand['gamma'],
            power=cand['power'],kBT=spec['kBT_star'],stress=True,stress_period=dt/2)
    u.registerInteraction(interaction);u.setInteraction(interaction,pv,pv)
    vv=mir.Integrators.VelocityVerlet_withPeriodicForce('vv',force=task['force_star'],direction='x') if task['force_star'] else mir.Integrators.VelocityVerlet('vv')
    u.registerIntegrator(vv);u.setIntegrator(vv,pv)
    every=spec['snapshot_every'];bins=int(round(domain[1]/spec['y_bin_star']))
    u.registerPlugins(mir.Plugins.createStats('stats',every=every,filename='native_stats.csv'))
    u.registerPlugins(mir.Plugins.createVirialPressurePlugin('pressure',pv,lambda r:1.,(1.,1.,1.),every,'pressure'))
    u.registerPlugins(mir.Plugins.createParticleChecker('checker',check_every=every))
    if density_enabled:
        saved_channels=[('densities','saved_density'),('positions','saved_positions'),('velocities','saved_velocities'),('__forces','saved_forces')]
        if observer_factory:saved_channels.append(('stresses','saved_stresses'))
        for source,saved in saved_channels:
            u.registerPlugins(mir.Plugins.createParticleChannelSaver('save_'+source,pv,source,saved))
    streams=[]
    def writer(name,header):
        f=(d/name).open('x',buffering=1);streams.append(f);w=csv.writer(f);w.writerow(header);return w
    if compute:
        moments=writer('moments.csv',['step','time_star','N','mass_star','kBT_raw_star','kBT_COM_star','kBT_thermal_star','kBT_bin_coarse_star','mean_vx','mean_vy','mean_vz','max_speed_star','compute_wall_s','snapshot_wall_s'])
        profiles=writer('profile_samples.csv',['step','time_star','bin','y_star','count','ux','uy','uz','sum_v2'])
        density=writer('density_samples.csv',['step','density_time_star','N','global_n_star','global_rho_star','kernel_number_mean','kernel_number_std','kernel_number_min','kernel_number_max','kernel_mass_mean','EOS_model_mean_star']) if density_enabled else None
        timings=writer('chunk_timing.csv',['step','time_star','compute_sync_s','sampling_transfer_output_s'])
    sync();setup_s=time.monotonic()-setup_started;done=0;index=0;compute_s=0.;sample_s=0.;first_compute=None;snapshot_count=0
    marks=sorted(set(spec.get('snapshot_marks_steps',[every]+[max(every,int(spec['steps']*f)//every*every) for f in [.25,.5,.75,1.]])));observations=[]
    for index,done,chunk_s in continuous_chunks(u,d,compute,spec['steps'],every,dt,sync):
        compute_s+=chunk_s
        if first_compute is None:first_compute=chunk_s
        if compute:
            start=time.monotonic();pos=np.asarray(pv.getCoordinates(),float);vel=np.asarray(pv.getVelocities(),float)
            if not np.isfinite(pos).all() or not np.isfinite(vel).all():raise RuntimeError('NONFINITE_PARTICLE_STATE')
            N=len(vel);mean=vel.mean(axis=0);raw=mass*np.sum(vel**2)/(3*N);com=mass*np.sum((vel-mean)**2)/(3*(N-1))
            coarse,counts,av,v2=thermal_bins(pos,vel,domain[1],bins,mass)
            fine=thermal_bins(pos,vel,domain[1],bins*2,mass)[0]
            measured_temperature=fine if task['kind']=='flow' else com
            profiles.writerows([done,done*dt,j,(j+.5)*spec['y_bin_star'],counts[j],*av[j],v2[j]] for j in range(bins))
            if density_enabled:
                den_host=channel(pv,'saved_density').reshape(-1).astype(float)
                if len(den_host)!=N or not np.isfinite(den_host).all() or np.min(den_host)<=0:raise RuntimeError('INVALID_NATIVE_DENSITY_FIELD')
                eos=float(np.mean(cand['sound_speed_star']**2*(mass*den_host-cand['rho_0_star']))) if cand['method']=='SDPD' else ''
                density.writerow([done,(done-1)*dt,N,N/np.prod(domain),N*mass/np.prod(domain),np.mean(den_host),np.std(den_host),np.min(den_host),np.max(den_host),mass*np.mean(den_host),eos])
                if marks and done>=marks[0]:
                    while marks and done>=marks[0]:marks.pop(0)
                    x4=channel(pv,'saved_positions');x=x4[:,:3];v4=channel(pv,'saved_velocities');f=channel(pv,'saved_forces')[:,:3]
                    current_v4=channel(pv,'velocities');current_x4=channel(pv,'positions')
                    ids_equal=bool(np.array_equal(x4[:,3].view(np.uint32),current_x4[:,3].view(np.uint32)) and np.array_equal(v4[:,3].view(np.uint32),current_v4[:,3].view(np.uint32)))
                    coordinate_error=pos-current_x4[:,:3]-.5*np.array(domain);coordinate_error-=np.array(domain)*np.rint(coordinate_error/np.array(domain))
                    checks=phase_audit(x,v4[:,:3],f,current_x4[:,:3],current_v4[:,:3],dt,mass,domain)
                    checks.update(step=done,pre_time_star=(done-1)*dt,post_time_star=done*dt,packed_ID_words_preserved=ids_equal,
                        getter_velocity_max_abs=float(np.max(abs(vel-current_v4[:,:3]))),getter_global_coordinate_PBC_max_abs=float(np.max(abs(coordinate_error))),
                        velocity=velocity_statistics(vel,mass),force_RMS_star=float(np.sqrt(np.mean(f.astype(float)**2))))
                    if not ids_equal:raise RuntimeError('PARTICLE_ID_ORDER_MISMATCH')
                    observations.append(checks)
                    for label,values in [('pre_velocities',v4[:,:3]),('pre_forces',f),('post_velocities',current_v4[:,:3]),('post_positions',current_x4[:,:3])]:
                        np.save(d/f'snapshot_{snapshot_count:02d}_{label}.npy',values,allow_pickle=False)
                    if snapshot_count==0:
                        layouts={}
                        for name in ['positions','velocities','saved_density','saved_positions','saved_velocities','saved_forces']:
                            interface=pv.local.per_particle[name].__cuda_array_interface__;shape,dtype,strides,size=cuda_layout(interface)
                            layouts[name]={'shape':shape,'dtype':dtype.str,'strides':strides,'byte_span':size,'owned_only':True}
                        (d/'channel_layouts.json').write_text(json.dumps(layouts,indent=2)+'\n')
                    np.save(d/f'snapshot_{snapshot_count:02d}_positions.npy',x,allow_pickle=False)
                    np.save(d/f'snapshot_{snapshot_count:02d}_kernel_number.npy',den_host,allow_pickle=False)
                    (d/f'snapshot_{snapshot_count:02d}.json').write_text(json.dumps({'step':done,'time_star':(done-1)*dt,'phase':'beforeIntegration, matching native density','coordinate_frame':'native local frame; periodic distances are translation invariant'}))
                    if observer:observer.snapshot(snapshot_count,done)
                    snapshot_count+=1
            if observer:observer.sample(done,vel,den_host,compute_s,sample_s,start,chunk_s)
            spent=time.monotonic()-start;sample_s+=spent
            moments.writerow([done,done*dt,N,mass,raw,com,measured_temperature,coarse,*mean,np.max(np.linalg.norm(vel,axis=1)),compute_s,sample_s])
            timings.writerow([done,done*dt,chunk_s,spent])
    if compute:
        for f in streams:f.close()
        (d/'phase_observations.json').write_text(json.dumps(observations,indent=2,allow_nan=False)+'\n')
        result={'steps':done,'actual_N':N,'time_star':done*dt,'setup_wall_s':setup_s,'compute_wall_s':compute_s,
                'first_chunk_compute_s':first_compute,'steady_compute_wall_s':compute_s-(first_compute or 0),
                'steady_steps':max(0,done-every),'snapshot_wall_s':sample_s,'density_snapshots':snapshot_count,
                'density_enabled':density_enabled,'synchronization':'u.run followed by cudaDeviceSynchronize before timing end',
                'status':'COMPLETED_PLANNED_STEPS' if done==spec['steps'] else 'COOPERATIVE_PARTIAL'}
        (d/'worker_completion.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        if observer:observer.finish(result)


if __name__=='__main__':main()
