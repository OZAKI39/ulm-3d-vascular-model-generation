"""Native DPD + one deformable membrane. Adapted from Mirheo MIT examples.

Sources: membranes_solvents.py, walls/analytic/couette.py, membrane/wlc.py.
Copyright ETH Zurich; see MIRHEO_LICENSE. Continuous diagnostic: one native run per phase with built-in HDF5 observation.
All checkpoint-free relaxation and in-memory wall-switch work is timed.
"""
import argparse
import csv
import ctypes as C
import gc
import json
import math
import os
from pathlib import Path
import time
import numpy as np
from common import write, sha, read
from rbc_checks import scan_frames, preparation, probe_checks

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',required=True);a=p.parse_args()
    s=json.loads(Path(a.spec).read_text());c=s['config'];d=Path.cwd();role=s['role'];dt=s['dt']
    if not (s.get('continuous_observation') and s.get('cloud_full_protocol')=='CLOUD_FULL_RBC_V1' and not s.get('stop_after_relaxation') and role=='main' and s['steps']==400000 and s['prep_steps']==60000 and dt==.0005):raise RuntimeError('CLOUD_FULL_RBC_V1_REQUIRED')
    if os.environ.get('OMPI_COMM_WORLD_SIZE')!='2' or os.environ.get('CLOUD_JOB_AUTHORIZED')!=s['cloud_plan_sha256']:raise RuntimeError('CLOUD_AUTHORIZED_TWO_RANK_LAUNCH_REQUIRED')
    from py_scripts.solver_benchmark.native_mpi import NativeMPI
    from cloud_geometry import read_off, write_off, order_vertices, geometry
    mpi=NativeMPI();mpi.check(mpi.lib.MPI_Init(None,None));compute=mpi.rank==0
    if s.get("native_library"):
        import hashlib, sys
        native=Path(s["native_library"])
        if hashlib.sha256(native.read_bytes()).hexdigest()!=s["native_library_sha256"]:raise RuntimeError("NATIVE_LIBRARY_HASH_MISMATCH")
        sys.path.insert(0,str(native.parent))
    import mirheo as mir
    import libmirheo
    if str(Path(libmirheo.__file__).resolve())!=str(native.resolve()) or sha(libmirheo.__file__)!=s['native_library_sha256']:raise RuntimeError('LOADED_LIBRARY_MISMATCH')
    write(d/f'rank_{mpi.rank}.json',dict(rank=mpi.rank,pid=os.getpid(),mpi_size=int(os.environ['OMPI_COMM_WORLD_SIZE']),library=libmirheo.__file__,sha256=sha(libmirheo.__file__)))
    if compute:write(d/'loaded_library.json',dict(path=libmirheo.__file__,sha256=sha(libmirheo.__file__),rank=mpi.rank))
    runtime=C.CDLL('libcudart.so.12') if compute else None
    if compute:runtime.cudaMemcpy.argtypes=[C.c_void_p,C.c_void_p,C.c_size_t,C.c_int]
    def sync():
        if compute and runtime.cudaDeviceSynchronize()!=0:raise RuntimeError('CUDA_SYNC_FAILED')
    def channel(pv,name):
        interface=pv.local.per_particle[name].__cuda_array_interface__;out=np.empty(interface['shape'],dtype=interface['typestr'])
        if runtime.cudaMemcpy(out.ctypes.data,interface['data'][0],out.nbytes,2):raise RuntimeError('CUDA_COPY_FAILED')
        return out
    before=time.perf_counter();dp=c['dpd'];geo=c['geometry'];L=geo['periodic_length'];H=geo['gap'];pad=geo['wall_layer']
    vertices,faces=read_off(s['mesh']);shape=geometry(vertices,faces);has_cell=role not in ('material','empty')
    material=role=='material';probe=role=='membrane_probe'
    domain=tuple(c['material_check']['domain']) if material else (L,L,H+2*pad)
    handles=[]
    def membrane(u,prior=None):
        initial=s['mesh']
        if prior is not None and compute:
            initial=str(d/'same_process_relaxed_geometry.off')
            write_off(initial,prior['rbc']['positions']-[L/2,L/2,H/2+pad],faces)
        # Native MembraneIC initializes oldPositions consistently with the deformed
        # geometry; the second OFF retains the original stress-free reference.
        mesh=mir.ParticleVectors.MembraneMesh(initial,s['mesh']);rbc=mir.ParticleVectors.MembraneVector('rbc',mass=c['mirheo_membrane']['mass'],mesh=mesh)
        u.registerParticleVector(rbc,mir.InitialConditions.Membrane([[L/2,L/2,H/2+pad,1.,0.,0.,0.]]))
        mp=c['mirheo_membrane'];params={k:mp[k] for k in ('x0','mpow','ka','ka_tot','kv_tot','theta','gammaC','kBT')}
        params.update(ks=s.get('ks',mp['ks_candidate']),kb=s.get('kb',mp['kb_candidate']),tot_area=shape['area'],tot_volume=shape['volume'])
        force=mir.Interactions.MembraneForces('membrane',mp['model'],mp['bending'],stress_free=mp['stress_free'],**params)
        u.registerInteraction(force);u.setInteraction(force,rbc,rbc);handles.extend([mesh,force]);return rbc
    def dpd(u,name,first,second,mem=False):
        interaction=mir.Interactions.Pairwise(name,dp['rc'],kind='DPD',a=0. if mem else dp['a'],gamma=dp['membrane_gamma'] if mem else dp['gamma'],kBT=dp['kBT'],power=dp['power'])
        u.registerInteraction(interaction);u.setInteraction(interaction,first,second);handles.append(interaction);return interaction
    def coordinator(label):
        os.environ['RBC_REPAIR_PHASE']=label
        return mir.Mirheo((1,1,1),domain,debug_level=s.get('debug_level',1),log_filename=label,no_splash=True)
    def capture(pv):
        if not compute:return None
        ids=np.asarray(pv.get_indices(),np.int64);order=np.argsort(ids)
        return dict(ids=ids[order],positions=np.asarray(pv.getCoordinates(),float)[order],velocities=np.asarray(pv.getVelocities(),float)[order])
    def make_pv(u,name,state=None):
        pv=mir.ParticleVectors.ParticleVector(name,mass=dp['mass'])
        ic=mir.InitialConditions.FromArray(state['positions'].tolist(),state['velocities'].tolist()) if state is not None and compute else mir.InitialConditions.Uniform(dp['number_density'])
        u.registerParticleVector(pv,ic);return pv
    streams=[]
    def writer(name,cols):
        if not compute:return None
        f=(d/name).open('x',buffering=1);streams.append(f);out=csv.writer(f);out.writerow(cols);return out
    profile=writer('profiles.csv',['step','phase','time_star','strain','bin','z','count','ux','uy','uz','density'])
    localflow=writer('local_flow.csv',['step','phase','time_star','strain','ix','iy','iz','x','y','z','count','ux','uy','uz','density']) if not material else None
    moments=writer('moments.csv',['step','phase','time_star','strain','N','inner_N','temperature','max_speed','wall_crossings'])
    membrane_out=writer('vertices.csv',['step','phase','time_star','strain','vertex','x','y','z','fx','fy','fz'])
    timings=writer('timings.csv',['step','phase','chunk_steps','compute_s','output_s'])
    timings_total=dict(setup_s=0.,relaxation_s=0.,coupled_s=0.,output_s=0.,wall_preparation_s=0.,state_handoff_s=0.)
    # A coordinator cannot change wall velocity after initialization in this version.
    # Positions/velocities persist in RAM. IDs, classification, oldPositions, forces and RNG are regenerated; this is NOT a full-state checkpoint.
    # No external prepared states are accepted; both stages are part of E2E.
    state=None;total_steps=0;initialN=round(L*L*H*dp['number_density']);initial_innerN=None;full=True
    def stage(label,count,sample,wall_speed,prior):
        nonlocal total_steps,initialN,initial_innerN,full
        t=time.perf_counter();u=coordinator(label);walls=[];wall_pvs=[];probe_ids={}
        fluid=make_pv(u,'outer',prior['fluid'] if prior is not None else None)
        if compute and prior is None:
            rng=np.random.default_rng(s['seed']);v=rng.normal(0,math.sqrt(dp['kBT']/dp['mass']),(len(fluid.getCoordinates()),3));v-=v.mean(axis=0);fluid.setVelocities(v.tolist())
        interaction=mir.Interactions.Pairwise('solvent',dp['rc'],kind='DPD',a=dp['a'],gamma=dp['gamma'],kBT=dp['kBT'],power=dp['power']);u.registerInteraction(interaction)
        vv=mir.Integrators.VelocityVerlet('vv');u.registerIntegrator(vv)
        if not material:
            for index,sign,z in ((0,-1,pad),(1,1,H+pad)):
                wall=mir.Walls.MovingPlane('wall'+str(index),normal=(0,0,sign),pointThrough=(0,0,z),velocity=(sign*wall_speed,0,0));u.registerWall(wall,1000);walls.append(wall)
                t0=time.perf_counter()
                if prior is None:
                    wp=u.makeFrozenWallParticles(pvName='wallpv'+str(index),walls=[wall],interactions=[interaction],integrator=mir.Integrators.VelocityVerlet('wallrelax'+str(index)),number_density=dp['number_density'],mass=dp['mass'],dt=dt,nsteps=dp['wall_relax_steps'])
                else:
                    wp=make_pv(u,'wallpv'+str(index),prior['walls'][index] if compute else None);wall.attachFrozenParticles(wp)
                if compute:wp.setVelocities(np.tile([sign*wall_speed,0.,0.],(len(wp.getCoordinates()),1)).tolist())
                sync();timings_total['wall_preparation_s']+=mpi.allreduce(time.perf_counter()-t0)
                move=mir.Integrators.Translate('move'+str(index),velocity=(sign*wall_speed,0,0));u.registerIntegrator(move);u.setIntegrator(move,wp);handles.extend([move,wall,wp]);wall_pvs.append(wp)
        rbc=None;inner=None
        if has_cell:
            rbc=membrane(u,prior)
            if prior is not None and compute:rbc.setVelocities(prior['rbc']['velocities'].tolist())
            checker=mir.BelongingCheckers.Mesh('classifier');u.registerObjectBelongingChecker(checker,rbc)
            inner=u.applyObjectBelongingChecker(checker,fluid,correct_every=0,inside='inner');handles.append(checker)
            u.setInteraction(interaction,inner,inner);u.setInteraction(interaction,inner,fluid);u.setIntegrator(vv,inner);u.setIntegrator(vv,rbc)
            dpd(u,'outer_membrane',fluid,rbc,True);dpd(u,'inner_membrane',inner,rbc,True)
            from py_scripts.single_rbc_repair.coupling import bind_bouncers
            handles.extend(bind_bouncers(mir,u,rbc,fluid,inner,s.get('bouncer_policy','shared')))
            u.registerPlugins(mir.Plugins.createForceSaver('membrane_force',rbc))
        u.setInteraction(interaction,fluid,fluid);u.setIntegrator(vv,fluid)
        for w,wp in zip(walls,wall_pvs):
            for pv in ([fluid,inner] if has_cell else [fluid]):u.setWall(w,pv);u.setInteraction(interaction,wp,pv)
            # Detect RBC-wall crossings in saved meshes; do not mask a failure by bouncing membrane nodes.
        if material and label=='shear':
            driven=mir.Integrators.VelocityVerlet_withPeriodicForce('driven',force=c['material_check']['force_per_particle'],direction='x');u.registerIntegrator(driven);u.deregisterIntegrator(vv);u.setIntegrator(driven,fluid);handles.append(driven)
        u.registerPlugins(mir.Plugins.createParticleChecker('finite_check',check_every=1000))
        sync();timings_total['setup_s']+=mpi.allreduce(time.perf_counter()-t)
        def sample_state(done):
            if not compute:return
            states=[capture(fluid)]+([capture(inner)] if has_cell else [])
            pos=np.concatenate([x['positions'] for x in states]);vel=np.concatenate([x['velocities'] for x in states]);ids=np.concatenate([x['ids'] for x in states])
            if not np.isfinite(pos).all() or not np.isfinite(vel).all() or len(np.unique(ids))!=len(ids):raise RuntimeError('BAD_FLUID_STATE_OR_DUPLICATE_ID')
            if initialN is not None and len(ids)!=initialN:raise RuntimeError('FLUID_MASS_CHANGED')
            time_star=done*dt;strain=done*dt*c['protocol']['shear_rate'] if label=='shear' else 0.
            axis=1 if material else 2;offset=0 if material else pad;length=domain[axis] if material else H;bins=c['material_check']['bins'] if material else round(H)
            bin_id=np.clip(((pos[:,axis]-offset)/length*bins).astype(int),0,bins-1);n=np.bincount(bin_id,minlength=bins);means=np.array([np.bincount(bin_id,weights=vel[:,k],minlength=bins)/np.maximum(n,1) for k in range(3)]).T
            thermal=float(np.sum((vel-means[bin_id])**2)*dp['mass']/(3*max(1,len(vel)-np.count_nonzero(n))))
            vol=np.prod(domain) if material else L*L*H
            for b in range(bins):profile.writerow([total_steps+done,label,time_star,strain,b,(b+.5)*length/bins,int(n[b]),*means[b],n[b]*dp['mass']/(vol/bins)])
            if not material:
                xyz=pos-[0,0,pad];grid=np.clip((xyz/np.array([L,L,H])*8).astype(int),0,7);cell=(grid[:,0]*8+grid[:,1])*8+grid[:,2]
                cn=np.bincount(cell,minlength=512);cu=np.array([np.bincount(cell,weights=vel[:,k],minlength=512)/np.maximum(cn,1) for k in range(3)]).T
                for ix in range(8):
                    for iy in range(8):
                        for iz in range(8):
                            j=(ix*8+iy)*8+iz
                            localflow.writerow([total_steps+done,label,time_star,strain,ix,iy,iz,(ix+.5)*L/8,(iy+.5)*L/8,(iz+.5)*H/8,int(cn[j]),*cu[j],cn[j]*dp['mass']/(L*L*H/512)])
            crossings=0 if material else int(np.count_nonzero((pos[:,2]<pad-1e-5)|(pos[:,2]>H+pad+1e-5)))
            moments.writerow([total_steps+done,label,time_star,strain,len(pos),len(states[1]['ids']) if has_cell else 0,thermal,float(np.linalg.norm(vel,axis=1).max()),crossings])
            if has_cell:
                rm=capture(rbc);v=rm['positions']-[0,0,pad]
                force=order_vertices(rbc.get_indices(),channel(rbc,'forces'),len(vertices)) if done else np.zeros_like(v)
                for i,(x,f) in enumerate(zip(v,force)):membrane_out.writerow([total_steps+done,label,time_star,strain,i,*x,*f])
                # Deterministic identity-aware probes; analysis reports limited leakage evidence separately.
                chosen=[]
                for membership,st in enumerate(states):
                    if membership not in probe_ids:
                        probe_ids[membership]=st['ids'][np.arange(0,len(st['ids']),max(1,len(st['ids'])//96))[:96]].copy()
                    take=np.searchsorted(st['ids'],probe_ids[membership])
                    if np.any(take>=len(st['ids'])) or not np.array_equal(st['ids'][take],probe_ids[membership]):raise RuntimeError('PROBE_IDS_CHANGED_WITHIN_PHASE')
                    chosen.extend([[membership,int(st['ids'][i]),*(st['positions'][i]-[0,0,pad])] for i in take])
                np.savez_compressed(d/f'fluid_probes_{label}_{done:08d}.npz',probes=np.array(chosen),strain=strain,phase=label,global_step=total_steps+done,membrane_vertex_ids=rm['ids'])
                if s.get('diagnostic_state',False):
                    np.savez_compressed(d/f'membrane_state_{label}_{done:08d}.npz',ids=rm['ids'],positions=v,velocities=rm['velocities'],phase=label,global_step=total_steps+done,strain=strain)
        # Native dump callbacks observe positions before forces at the recorded
        # step, without returning through Simulation::run initialization.
        dump_root=d/f'native_{label}'
        if compute:dump_root.mkdir(exist_ok=False)
        mpi.Barrier()
        # Membrane at fixed cadence; full fluid only at phase step zero for cold-state and handoff audits.
        for pv,name in ((rbc,'rbc'),(fluid,'outer'),(inner,'inner')):
            channels=['forces'] if name=='rbc' else []
            u.registerPlugins(mir.Plugins.createDumpParticles('observe_'+name,pv,dump_every=sample if name=='rbc' else count,channel_names=channels,path=str(dump_root/name)))
        for name,pvs in [('fluid',[fluid,inner]),('outer',[fluid]),('inner',[inner]),('membrane',[rbc])]:
            u.registerPlugins(mir.Plugins.createStats('stats_'+name,every=1000,pvs=pvs,filename=str(d/f'native_stats_{label}_{name}.csv')))
        u.registerPlugins(mir.Plugins.createDumpAverage('local_response',[fluid,inner],sample_every=sample,dump_every=sample,
            bin_size=(L/8,L/8,(H+2*pad)/8),channels=['velocities'],path=str(dump_root/'local_flow')))
        if compute:
            write(d/'phase_timings.json',dict(**timings_total,last_phase=label,native_returned=False))
            write(d/f'phase_{label}.json',dict(phase=label,native_returned=False,planned_steps=count,successful_returned_steps=None,dt=dt,started_at=time.time()))
            write(d/'progress.json',dict(phase=label,native_returned=False,planned_steps=count,successful_returned_steps=None,dt=dt,prep_returned_steps=total_steps))
            if prior is not None:
                actual=capture(rbc);ordered_old=order_vertices(rbc.get_indices(),channel(rbc,'old_positions')[:,:3],len(vertices))
                expected_local=actual['positions']-np.array(domain)/2
                error=float(np.max(np.abs(ordered_old-expected_local)))
                write(d/'handoff_membrane_initialization.json',dict(oldPositions_reinitialized_to_current=True,maximum_coordinate_error=error,
                    reference_mesh_sha256=sha(s['mesh']),positions_max_error=float(np.max(np.abs(actual['positions']-prior['rbc']['positions']))),
                    velocities_max_error=float(np.max(np.abs(actual['velocities']-prior['rbc']['velocities']))),state_continuity='partial: coordinates/velocities only'))
                if error>1e-5:raise RuntimeError('MEMBRANE_OLD_POSITIONS_INITIALIZATION_MISMATCH')
        from py_scripts.single_rbc_repair.continuous_protocol import evolve_stage
        mpi.Barrier();t=time.perf_counter()
        done=evolve_stage(u,count,dt,sample);sync()
        if compute and int(u.getState().current_step)!=count:raise RuntimeError('NATIVE_STEP_COUNT_MISMATCH')
        if compute:
            actual=dict(phase=label,native_returned=True,successful_returned_steps=int(u.getState().current_step),dt=dt,native_time_star=float(u.getState().current_time),
                        strain=float(u.getState().current_time)*c['protocol']['shear_rate'] if label=='shear' else 0.)
            write(d/f'phase_{label}.json',actual);write(d/'progress.json',actual)
        cost=mpi.allreduce(time.perf_counter()-t)
        timings_total['relaxation_s' if label=='relaxation' else 'coupled_s']+=cost
        t=time.perf_counter();sample_state(done)
        outcost=mpi.allreduce(time.perf_counter()-t);timings_total['output_s']+=outcost
        if compute:timings.writerow([total_steps+done,label,count,cost,outcost])
        transfer_start=time.perf_counter();state=None
        if compute:
            fs=[capture(fluid)]+([capture(inner)] if has_cell else []);ids=np.concatenate([x['ids'] for x in fs]);order=np.argsort(ids)
            state=dict(fluid=dict(ids=ids[order],positions=np.concatenate([x['positions'] for x in fs])[order],velocities=np.concatenate([x['velocities'] for x in fs])[order]),walls=[capture(x) for x in wall_pvs],rbc=capture(rbc) if has_cell else None)
            if initial_innerN is None:initial_innerN=len(fs[1]['ids'])
            if label=='relaxation':
                state['fluid']['membership']=np.concatenate([np.full(len(st['ids']),member) for member,st in enumerate(fs)])[order]
                np.savez_compressed(d/'handoff_saved_arrays.npz',**{'fluid_'+k:v for k,v in state['fluid'].items()},
                    **{'rbc_'+k:v for k,v in state['rbc'].items()},**{f'wall{j}_'+k:v for j,wall in enumerate(state['walls']) for k,v in wall.items()},
                    membrane_old_positions_local=order_vertices(rbc.get_indices(),channel(rbc,'old_positions')[:,:3],len(vertices)))
        total_steps+=done;del u;handles.clear();gc.collect();mpi.Barrier();timings_total['state_handoff_s']+=mpi.allreduce(time.perf_counter()-transfer_start)
        if compute:write(d/'phase_timings.json',dict(**timings_total,last_phase=label,native_returned=True))
        return state,done
    if material:prep=round(c['material_check']['relaxation_time']/dt);steps=round(c['material_check']['drive_time']/dt);sample=round(c['material_check']['sample_time']/dt)
    else:prep=s['prep_steps'];steps=s['steps'];sample=s['sample_steps']
    state,done=stage('relaxation',prep,s['prep_sample_steps'],0.,None)
    # Save a checksum of the in-memory handoff to prove provenance; no reload path exists.
    if compute:
        import hashlib
        (d/'handoff.json').write_text(json.dumps(dict(source='same process freshly computed relaxation',checkpoints=0,fluid_particles=initialN,inner_particles=initial_innerN,positions_sha256=hashlib.sha256(state['fluid']['positions'].tobytes()).hexdigest(),old_ids_sha256=hashlib.sha256(state['fluid']['ids'].tobytes()).hexdigest(),ids_remapping_planned_by_sorted_old_id=True,execution_state='prepared arrays only; second coordinator subject to preparation gate')))
    relaxation_done=done
    gate_pass=None
    if compute:
        t=time.perf_counter();cache=scan_frames(d,s,final=True);gate=preparation(list(cache.values()),vertices,s)
        members=probe_checks(d,s,list(cache.values()));gate['membership']=members
        with (d/'moments.csv').open() as f:observed=list(csv.DictReader(f))
        last=observed[-1];gate['endpoint_temperature']=float(last['temperature'])
        gate['fluid_endpoint_pass']=int(last['N'])==initialN and int(last['wall_crossings'])==0 and abs(float(last['temperature'])/c['dpd']['kBT']-1)<=c['screen']['temperature_relative_error']
        gate['passed']=gate['passed'] and gate['fluid_endpoint_pass'] and members['tested_point_frames']>0 and members['confirmed_mismatches']==0 and all(q['status']=='CHECKED' for q in members['checks'])
        gate['status']='PASS' if gate['passed'] else 'PREPARATION_QUALITY_FAILED_OR_INCOMPLETE'
        gate['cpu_analysis_s']=time.perf_counter()-t
        write(d/'preparation_gate.json',gate);gate_pass=gate['passed']
    gate_pass=mpi.bcast(gate_pass)
    if not gate_pass:
        for f in streams:f.close()
        if compute:write(d/'stop_reason.json',dict(reason='PREPARATION_QUALITY_FAILED_OR_INCOMPLETE',prep_successful_returned_steps=relaxation_done,shear_successful_returned_steps=0,strain=0.,no_retry=True))
        mpi.Barrier();mpi.check(mpi.lib.MPI_Finalize());return

    if s.get('stop_after_relaxation',False):
        done=0;steps=0;full=full and relaxation_done==prep
    elif full:state,done=stage('shear',steps,sample,c['protocol']['shear_rate']*H/2,state)
    if compute:
        for f in streams:f.close()
        (d/'completion.json').write_text(json.dumps(dict(completed=full and done==steps,role=role,cell_count=int(has_cell),vertices=len(vertices) if has_cell else 0,actual_steps=total_steps,prep_steps=prep,shear_steps=done,strain_end=done*dt*c['protocol']['shear_rate'] if not material else None,dt_star=dt,fluid_particles=initialN,rank_count=2,precision_bits=32,**timings_total,worker_elapsed_s=time.perf_counter()-before,timing_note='setup includes wall_preparation; native output is included in relaxation/coupled; handoff and endpoint output separately timed; do not sum contained categories',coordinators_sequential=1 if s.get('stop_after_relaxation',False) else 2,wall_hidden_steps=2*dp['wall_relax_steps'] if not material else 0,checkpoints_loaded=0,membrane_updates=total_steps if has_cell else 0,interior_viscosity_ratio=1.,bouncer='bounce_back',membrane_thermostat_gammaC=0.,classifier_corrections=0,initial_classifier_calls=2,run_api_mode='one_continuous_call_per_phase',native_dump_sampling='beforeForces at the indexed step; saved forces are from the preceding force-saver update',native_dump_output_included_in_relaxation_s=True),indent=2))
    mpi.Barrier();mpi.check(mpi.lib.MPI_Finalize())

if __name__=='__main__':main()
