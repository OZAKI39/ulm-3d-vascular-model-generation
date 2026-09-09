"""Authorized native SDPD restart diagnostic/segment worker; imports are CPU safe.

Uses the existing continuous_chunks and EquilibrationObserver reductions. Only
fresh runs initialize velocities. Diagnostic restores intentionally expose the
native missing-RNG limitation and are never accepted into a formal chain.
"""
import argparse
import csv
import ctypes
import math
import os
from pathlib import Path
import time
import numpy as np
from py_scripts.fluid_physics.common import read_json,write_json,fingerprint
from .gpu_worker import continuous_chunks
from .equilibration_worker import EquilibrationObserver,RAW_COLUMNS,mechanical_pressure
from .observations import cuda_layout
from .extended_restart import restore_coordinator,verify_checkpoint,require_readable_checkpoint_channels


class SegmentObserver(EquilibrationObserver):
    def __init__(self,spec,directory,pv,channel,absolute_start):
        self.spec,self.d,self.pv,self.channel=spec,Path(directory),pv,channel
        self.volume=float(np.prod(spec['domain_star']));self.initial_ids=np.asarray(pv.get_indices(),np.int64);self.n=len(self.initial_ids)
        if self.n!=spec['expected_N'] or len(np.unique(self.initial_ids))!=self.n:raise ValueError('PARTICLE_ID_OR_COUNT_ERROR')
        self.stream=(self.d/'raw_statistics.csv').open('x',buffering=1);self.writer=csv.DictWriter(self.stream,fieldnames=RAW_COLUMNS);self.writer.writeheader()
        self.snapshots=[]
        if absolute_start==0:
            row=self.velocity_row(0,np.asarray(pv.getVelocities(),float));self.writer.writerow(row)
        else:
            write_json(self.d/'restore_boundary.json',{'step':absolute_start,'time_star':absolute_start*spec['dt_star'],
                       'is_new_sample':False,'velocity_initialization_count':0,'COM_subtraction_count':0})

    def observation(self,step):
        vel=np.asarray(self.pv.getVelocities(),float);den=self.channel(self.pv,'saved_density').reshape(-1).astype(float)
        row=self.velocity_row(step,vel)
        row['pressure_star']=mechanical_pressure(self.channel(self.pv,'saved_stresses'),vel,self.spec['m_star'],self.volume,self.spec['dt_star'],step)[0]
        row['kernel_number_mean']=float(den.mean())
        return row

    def save_snapshot(self,index,step):
        prefix=self.d/f'snapshot_{index:02d}'
        for label,values in [('positions',self.channel(self.pv,'saved_positions')[:,:3]),
            ('pre_velocities',self.channel(self.pv,'saved_velocities')[:,:3]),
            ('pre_forces',self.channel(self.pv,'saved_forces')[:,:3]),
            ('post_positions',self.channel(self.pv,'positions')[:,:3]),
            ('post_velocities',self.channel(self.pv,'velocities')[:,:3]),
            ('kernel_number',self.channel(self.pv,'saved_density').reshape(-1))]:
            np.save(str(prefix)+'_'+label+'.npy',values,allow_pickle=False)
        self.snapshot(index,step)
        write_json(self.d/f'snapshot_{index:02d}.json',{'step':step,'pre_time_star':(step-1)*self.spec['dt_star'],
            'post_time_star':step*self.spec['dt_star'],'phase':'saved pre-integration positions/density/stress, post-integration x/v; not a restart checkpoint'})


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--spec',required=True);args=parser.parse_args()
    if os.environ.get('OMPI_COMM_WORLD_SIZE')!='2':raise RuntimeError('REQUIRES_TWO_RANK_AUTHORIZED_LAUNCHER')
    spec=read_json(Path(args.spec));d=Path.cwd();plan=spec['extended_plan'];task=spec['extended_task'];mode=task['mode']
    rank=int(os.environ['OMPI_COMM_WORLD_RANK']);compute=rank==0;restore=bool(spec.get('restore_directory'))
    diagnostic=mode=='diagnostic_restore';dt=spec['dt_star'];cand=spec['candidate'];start_clock=time.monotonic()
    if dt!=1e-6 or spec['task']['force_star']!=0 or cand['method']!='SDPD':raise ValueError('FIXED_UNFORCED_SDPD_REQUIRED')
    if restore:
        checkpoint=verify_checkpoint(spec['restore_directory'],require_rng=not diagnostic,expected_manifest=spec['restore_manifest_sha256'])
        require_readable_checkpoint_channels(spec['restore_directory'])  # CPU gate before importing Mirheo/CUDA
        if (checkpoint['frozen_parameters_sha256']!=fingerprint(plan['parameters'])
                or checkpoint['source_contract_sha256']!=plan['restart_contract_sha256']
                or checkpoint['dt_star']!=dt or checkpoint['mass_star']!=spec['m_star']):
            raise ValueError('CHECKPOINT_FROZEN_PARAMETERS_OR_SOURCE_CONTRACT_MISMATCH')
    import mirheo as mir
    runtime=ctypes.CDLL('libcudart.so.12') if compute else None
    if compute:
        runtime.cudaMemcpy.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_size_t,ctypes.c_int];runtime.cudaMemcpy.restype=ctypes.c_int
    def sync():
        if compute and runtime.cudaDeviceSynchronize()!=0:raise RuntimeError('CUDA_SYNC_FAILED')
    def channel(pv,name):
        interface=pv.local.per_particle[name].__cuda_array_interface__;shape,dtype,strides,size=cuda_layout(interface)
        raw=np.empty(size,dtype=np.uint8)
        if size and runtime.cudaMemcpy(raw.ctypes.data,int(interface['data'][0]),size,2)!=0:raise RuntimeError('CUDA_COPY_FAILED')
        return np.ndarray(shape,dtype=dtype,buffer=raw,strides=tuple(strides)).copy()
    checkpoint_every=plan['restart_test']['N1'] if task['kind']=='restart_diagnostic' else plan['segment_steps']
    u=mir.Mirheo((1,1,1),tuple(spec['domain_star']),debug_level=1,log_filename='mirheo',no_splash=True,
                 checkpoint_every=checkpoint_every,checkpoint_folder=str(d/'native_checkpoint'),checkpoint_mode=plan['checkpoint_mode'])
    pv=mir.ParticleVectors.ParticleVector('pv',mass=spec['m_star'])
    initial=mir.InitialConditions.FromArray([],[]) if restore else mir.InitialConditions.Uniform(number_density=spec['task']['n_star'])
    u.registerParticleVector(pv,initial)
    if compute and not restore:
        count=len(pv.getCoordinates());vel=np.random.default_rng(spec['initial_velocity_seed']).normal(0,math.sqrt(spec['kBT_star']/spec['m_star']),(count,3))
        vel-=vel.mean(axis=0);pv.setVelocities(vel.tolist())
    den=mir.Interactions.Pairwise('density',spec['rc_star'],kind='Density',density_kernel='WendlandC2')
    u.registerInteraction(den);u.setInteraction(den,pv,pv)
    interaction=mir.Interactions.Pairwise('sdpd',spec['rc_star'],kind='SDPD',viscosity=cand['viscosity_mu_star'],
        kBT=spec['kBT_star'],EOS='Linear',sound_speed=cand['sound_speed_star'],rho_0=cand['rho_0_star'],
        density_kernel='WendlandC2',stress=True,stress_period=dt/2)
    u.registerInteraction(interaction);u.setInteraction(interaction,pv,pv)
    vv=mir.Integrators.VelocityVerlet('vv');u.registerIntegrator(vv);u.setIntegrator(vv,pv)
    every=spec['snapshot_every']
    u.registerPlugins(mir.Plugins.createStats('stats',every=every,filename='native_stats.csv'))
    u.registerPlugins(mir.Plugins.createVirialPressurePlugin('pressure',pv,lambda r:1.,(1.,1.,1.),every,'pressure'))
    u.registerPlugins(mir.Plugins.createParticleChecker('checker',check_every=every))
    for source,saved in [('densities','saved_density'),('positions','saved_positions'),('velocities','saved_velocities'),('__forces','saved_forces'),('stresses','saved_stresses')]:
        u.registerPlugins(mir.Plugins.createParticleChannelSaver('save_'+source,pv,source,saved))
    start=0
    if restore:
        m=restore_coordinator(u,spec['restore_directory'],diagnostic_only=diagnostic,
                              expected_manifest=spec['restore_manifest_sha256'],compute=compute);start=m['absolute_step']
    sync()
    observer=SegmentObserver(spec,d,pv,channel,start) if compute else None
    def save_state(label):
        state=u.getState()
        np.savez(d/(label+'.npz'),ids=np.asarray(pv.get_indices(),np.int64),positions=np.asarray(pv.getCoordinates(),float),
                 velocities=np.asarray(pv.getVelocities(),float),step=state.current_step,time_star=state.current_time,mass_star=spec['m_star'])
    if compute:
        save_state('restored_before_advance' if restore else 'initial_state')
        write_json(d/'initialization_record.json',{'absolute_start_step':start,'time_star':start*dt,
             'velocity_initialization_count':0 if restore else 1,'COM_subtraction_count':0 if restore else 1,
             'restart_mode':'DIAGNOSTIC_INCOMPLETE_RNG' if diagnostic else 'FULL_CHECKPOINT' if restore else 'COLD_START',
             'coordinator_restart_called':restore,'extra_thermostat':False,'rank_layout':[1,1,1],'world_size':2,
             'native_stats_file_role':'Secondary: native Stats restore copies prior CSV before appending; never interpreted as newly sampled rows.',
             'raw_statistics_file_role':'Exclusive new segment measurements, absolute steps; restored boundary is a separate JSON, not a sample.'})
    N1,N2=plan['restart_test']['N1'],plan['restart_test']['N2'];cuts=[N1,N1+1,N1+N2] if task['kind']=='restart_diagnostic' else [start+1,start+200000]
    marks=[round(t/dt) for t in plan['snapshot_times_star']];done=0;compute_s=0.;sample_s=0.;snapshot_count=0
    committed_end=task.get('committed_end_step',start+spec['steps'])
    for index,done,chunk_s in continuous_chunks(u,d,compute,spec['steps'],every,dt,sync,start_step=start,stop_steps=cuts):
        absolute=start+done;compute_s+=chunk_s
        if compute:
            state=u.getState()
            if state.current_step!=absolute or not math.isclose(state.current_time,absolute*dt,abs_tol=2e-10,rel_tol=0):raise RuntimeError('NATIVE_ABSOLUTE_TIME_DIVERGED')
            ids=np.asarray(pv.get_indices(),np.int64)
            if not np.array_equal(np.sort(ids),np.sort(observer.initial_ids)):raise RuntimeError('PARTICLE_ID_CHANGED')
            tick=time.monotonic()
            if absolute%every==0 and absolute<=committed_end:
                observer.sample(absolute,np.asarray(pv.getVelocities(),float),channel(pv,'saved_density').reshape(-1),compute_s,sample_s,tick,chunk_s)
            if absolute in cuts or done==spec['steps']:
                save_state(f'state_{absolute}');write_json(d/f'observation_{absolute}.json',observer.observation(absolute))
            if absolute in marks:
                observer.save_snapshot(snapshot_count,absolute);snapshot_count+=1
            sample_s+=time.monotonic()-tick
    if compute:
        observer.stream.close()
        write_json(d/'worker_completion.json',{'status':'COMPLETED_PLANNED_STEPS' if done==spec['steps'] else 'COOPERATIVE_PARTIAL',
            'steps':start+done,'segment_advanced_steps':done,'segment_start_step':start,'committed_end_step':committed_end if done==spec['steps'] else None,
            'time_star':u.getState().current_time,'compute_wall_s':compute_s,'sampling_wall_s':sample_s,
            'worker_elapsed_s':time.monotonic()-start_clock,'snapshots':observer.snapshots,'diagnostic_only':task['kind']=='restart_diagnostic'})


if __name__=='__main__':main()
