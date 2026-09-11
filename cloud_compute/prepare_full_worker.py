"""Derive a bounded full worker from the already deployed continuous smoke adapter."""
from pathlib import Path
import difflib,hashlib,json
b=Path(__file__).resolve().parent
source=b/'cloud_repair_worker.py';old=source.read_text();s=old

def change(a,z):
    global s
    assert s.count(a)==1,(a,s.count(a))
    s=s.replace(a,z)

change("import numpy as np", "import numpy as np\nfrom common import write, sha, read\nfrom rbc_checks import scan_frames, preparation, probe_checks")
change("if not (s.get('continuous_observation') and s.get('stop_after_relaxation') and role=='main' and s['steps']==0):raise RuntimeError('CONTINUOUS_ZERO_SHEAR_DIAGNOSTIC_ONLY')", "if not (s.get('continuous_observation') and s.get('cloud_full_protocol')=='CLOUD_FULL_RBC_V1' and not s.get('stop_after_relaxation') and role=='main' and s['steps']==400000 and s['prep_steps']==60000 and dt==.0005):raise RuntimeError('CLOUD_FULL_RBC_V1_REQUIRED')")
change('import mirheo as mir', '''import mirheo as mir
    import libmirheo
    if str(Path(libmirheo.__file__).resolve())!=str(native.resolve()) or sha(libmirheo.__file__)!=s['native_library_sha256']:raise RuntimeError('LOADED_LIBRARY_MISMATCH')
    write(d/f'rank_{mpi.rank}.json',dict(rank=mpi.rank,pid=os.getpid(),mpi_size=int(os.environ['OMPI_COMM_WORLD_SIZE']),library=libmirheo.__file__,sha256=sha(libmirheo.__file__)))
    if compute:write(d/'loaded_library.json',dict(path=libmirheo.__file__,sha256=sha(libmirheo.__file__),rank=mpi.rank))''')
# Remove unreachable material-probe loops, so the sole numerical entry is evolve_stage.
a=s.index("    response=writer('native_response.csv'");z=s.index("    profile=writer('profiles.csv'");s=s[:a]+s[z:]
change("state=None;total_steps=0;initialN=None;initial_innerN=None;full=True", "state=None;total_steps=0;initialN=round(L*L*H*dp['number_density']);initial_innerN=None;full=True")
change("# Preserve all state in RAM inside this cold process, then rebuild only the native handles.", "# Positions/velocities persist in RAM. IDs, classification, oldPositions, forces and RNG are regenerated; this is NOT a full-state checkpoint.")
change("time_star=done*dt if label=='shear' else done*dt-c['protocol']['relaxation_time'];strain=max(0,time_star)*c['protocol']['shear_rate'] if not material else 0.", "time_star=done*dt;strain=done*dt*c['protocol']['shear_rate'] if label=='shear' else 0.")
change("for pv,name in ((rbc,'rbc'),):", "for pv,name in ((rbc,'rbc'),(fluid,'outer'),(inner,'inner')):")
change("dump_every=sample,channel_names=channels", "dump_every=sample if name=='rbc' else count,channel_names=channels")
change("# Cloud smoke observes native membrane frames only; no full-fluid frame stream.", "# Membrane at fixed cadence; full fluid only at phase step zero for cold-state and handoff audits.")
change("        from py_scripts.single_rbc_repair.continuous_protocol import evolve_stage", '''        for name,pvs in [('fluid',[fluid,inner]),('outer',[fluid]),('inner',[inner]),('membrane',[rbc])]:
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
        from py_scripts.single_rbc_repair.continuous_protocol import evolve_stage''')
change("        cost=mpi.allreduce(time.perf_counter()-t)", """        if compute:
            actual=dict(phase=label,native_returned=True,successful_returned_steps=int(u.getState().current_step),dt=dt,native_time_star=float(u.getState().current_time),
                        strain=float(u.getState().current_time)*c['protocol']['shear_rate'] if label=='shear' else 0.)
            write(d/f'phase_{label}.json',actual);write(d/'progress.json',actual)
        cost=mpi.allreduce(time.perf_counter()-t)""")
change("            if initialN is None:initialN=len(ids);initial_innerN=len(fs[1]['ids']) if has_cell else 0", """            if initial_innerN is None:initial_innerN=len(fs[1]['ids'])
            if label=='relaxation':
                state['fluid']['membership']=np.concatenate([np.full(len(st['ids']),member) for member,st in enumerate(fs)])[order]
                np.savez_compressed(d/'handoff_saved_arrays.npz',**{'fluid_'+k:v for k,v in state['fluid'].items()},
                    **{'rbc_'+k:v for k,v in state['rbc'].items()},**{f'wall{j}_'+k:v for j,wall in enumerate(state['walls']) for k,v in wall.items()},
                    membrane_old_positions_local=order_vertices(rbc.get_indices(),channel(rbc,'old_positions')[:,:3],len(vertices)))""")
change("state,done=stage('relaxation',prep,min(sample,max(1,prep//10)),0.,None)", "state,done=stage('relaxation',prep,s['prep_sample_steps'],0.,None)")
change("total_steps+=done;del u;handles.clear();gc.collect();mpi.Barrier();timings_total['state_handoff_s']+=mpi.allreduce(time.perf_counter()-transfer_start);return state,done", "total_steps+=done;del u;handles.clear();gc.collect();mpi.Barrier();timings_total['state_handoff_s']+=mpi.allreduce(time.perf_counter()-transfer_start)\n        if compute:write(d/'phase_timings.json',dict(**timings_total,last_phase=label,native_returned=True))\n        return state,done")
change("ids_remapped_by_sorted_old_id=True", "ids_remapping_planned_by_sorted_old_id=True,execution_state='prepared arrays only; second coordinator subject to preparation gate'")
change("    relaxation_done=done", """    relaxation_done=done
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
""")
change("initial_classifier_calls=1", "initial_classifier_calls=2")
# Explicit timing identities avoid summing the contained wall/native-output costs twice.
change("worker_elapsed_s=time.perf_counter()-before,coordinators_sequential", "worker_elapsed_s=time.perf_counter()-before,timing_note='setup includes wall_preparation; native output is included in relaxation/coupled; handoff and endpoint output separately timed; do not sum contained categories',coordinators_sequential")
(b/'rbc_full_worker.py').write_text(s)
patch=''.join(difflib.unified_diff(old.splitlines(True),s.splitlines(True),fromfile='cloud_repair_worker.py',tofile='rbc_full_worker.py'))
(b/'rbc_full_worker.patch').write_text(patch)
(b/'rbc_full_worker_derivation.json').write_text(json.dumps(dict(source=str(source),source_sha256=hashlib.sha256(old.encode()).hexdigest(),worker_sha256=hashlib.sha256(s.encode()).hexdigest(),patch_sha256=hashlib.sha256(patch.encode()).hexdigest(),physics_parameters_changed=False,native_sources_changed=False,scope='one continuous call per phase; full endpoint conditional on unchanged preparation quality; CPU observations and explicit partial handoff audit'),indent=2)+'\n')
