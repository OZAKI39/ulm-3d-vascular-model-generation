"""P8.2 recorder: identical P8.1 stepping, parameterized computation budgets only.

Baseline factor 1 retains the original guard ordering/arithmetic. Factors 4/16
increase call, no-accept, stationary and accepted-progress-window budgets. No
force, wall gap, radius, admission candidate or acceptance rule is changed.
"""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import json,os,time
import numpy as np
from .particle81_simulation import (DT,HORIZON,MAX_PROVIDER_CALLS,STATIONARY_STEPS,
    COLUMNS,environment,SavedTrajectoryStepper,IntegrationSafetyStop,dump)
from .particle8_replay import read,digest,canonical_hash
from .particle3_cases import json_safe
from .injection_admission import FiniteSizeAdmission
from .lammps_neighbors import ValidationNeighborPolicy
from .physical_time_refinement import PhysicalTimeRefinementError
from .nearfield_regularization import STATES

def integrate_one(event,dt=DT,output=None,force=False,guard_factor=1,checkpoint=None):
    if output is None: raise ValueError("Explicit P8.2 output required")
    if guard_factor not in [1,4,16]: raise ValueError("Guard factor must be predeclared: 1, 4 or 16")
    max_calls=MAX_PROVIDER_CALLS*guard_factor
    stationary_steps=STATIONARY_STEPS*guard_factor
    progress_window=64*guard_factor
    output=Path(output);folder=output/'trajectories';folder.mkdir(parents=True,exist_ok=True)
    pid=event['particle_id'];stem=f'mb_{pid:06d}';meta_path=folder/(stem+'.json');array_path=folder/(stem+'.npz')
    config=dict(dt_s=dt,horizon_s=HORIZON,maximum_provider_calls=max_calls,
        roundoff_progress_window_accepted_steps=progress_window,roundoff_progress_budget_multiplier=16,
        minimum_64_substep_progress_nominal_dt_fraction=.01,
        progress_check_provider_interval=128,maximum_rejected_provider_calls_without_acceptance=512*guard_factor,
        stationary_unchanged_nominal_steps=stationary_steps,nearfield_contract='NEAR_FIELD_REGULARIZATION_V1')
    before=None;original=None
    if checkpoint is not None:
        checkpoint=Path(checkpoint);original=read(checkpoint);before=np.load(checkpoint.with_suffix('.npz'))['samples']
        if original['end_reason']!='INTEGRATION_SAFETY_STOP' or not len(before):raise ValueError('Expected saved original safety-stop checkpoint')
        if original['birth_metadata_sha256']!=canonical_hash(event) or original['integration_config']['dt_s']!=dt:
            raise ValueError('Checkpoint initial state or nominal dt changed')
        if original['samples_sha256']!=digest(checkpoint.with_suffix('.npz')):raise ValueError('Historical checkpoint changed')
        config.update(checkpoint_metadata_sha256=digest(checkpoint),checkpoint_samples_sha256=original['samples_sha256'],
            continuation='EXACT_SAVED_ACCEPTED_ENDPOINT_AND_ORIGINAL_PENDING_BINARY_INTERVALS')
    if meta_path.exists() and not force:
        existing=read(meta_path)
        if existing['birth_metadata_sha256']!=canonical_hash(event) or existing['integration_config']!=config or digest(array_path)!=existing['samples_sha256']:
            raise ValueError('Resume data mismatch for stable ID '+str(pid))
        return existing
    env=environment();metadata=deepcopy(event);start=time.perf_counter()
    rng=np.random.default_rng(event['position_seed']);rng_initial=deepcopy(rng.bit_generator.state)
    source=SimpleNamespace(rng={'MB_POSITION':rng,'ADMISSION_RETRY':rng})
    admission=FiniteSizeAdmission(env.sampler,source,wall=env.wall,field=env.field)
    if checkpoint is None:p=admission.attempt(metadata,{})
    else:
        from .particle82_checkpoint import particle_at_saved_endpoint
        p=particle_at_saved_endpoint(original,before)
    result=dict(particle_id=pid,birth_metadata=event,birth_metadata_sha256=canonical_hash(event),integration_config=config,
        birth_time_s=event['birth_time_s'],radius_m=event['radius_m'],diameter_um=event['diameter_um'],initial_q=event['q'],
        position_rng_initial=rng_initial,position_rng_final=deepcopy(rng.bit_generator.state),admission_candidates=admission.candidates,
        original_geometry_preserved=True,artificial_inlet_offset_m=0.,independently_integrated=True,
        dynamics='UNCHANGED_PARTICLE65STEPPER_NO_LAMMPS_NO_MB_MB_OR_RBC_COUPLING',
        timestep_role='P81_DATASET_VALIDATION_ONLY_NOT_PRODUCTION_FROZEN',columns=COLUMNS,
        velocity_role='HELD_ON_PRECEDING_ACCEPTED_INTERVAL; FIRST_ROW_FREE_INITIAL_DIAGNOSTIC',
        terminal_event=None,completed=False,exit_outlet=None,exit_time_s=None,residence_time_s=None,path_length_m=0.,failure_detail=None)
    if checkpoint is not None:
        result=deepcopy(original)
        result.update(integration_config=config,failure_detail=None,
            continuation_checkpoint=dict(metadata_sha256=digest(checkpoint),samples_sha256=original['samples_sha256'],
                saved_prefix_count=len(before),restart_elapsed_time_s=float(before[-1,0]),
                admission_not_redrawn=True,historical_states_not_reintegrated=True,
                retry_policy='Unaccepted trials may be retried; cumulative provider count restored; rejected-call streak restarts. Nominal grid and pending depths restored.'))
    calls=0;stepper=None;last_checked_count=0;no_acceptance_checks=0
    last_position=None;last_sample=None
    if p is None:
        result['end_reason']='INLET_ADMISSION_GUARD_EXHAUSTED_UNRESOLVED'
        samples=np.empty((0,len(COLUMNS)))
    else:
        if checkpoint is None:
            xyz=env.sampler.triangles[metadata['last_candidate']['triangle_id']]
            uv=np.linalg.lstsq((xyz[1:]-xyz[0]).T,p.position-xyz[0],rcond=None)[0]
            result.update(inlet_face=metadata['last_candidate']['triangle_id'],inlet_barycentric=[1.-float(uv.sum()),*uv.tolist()],
                inlet_position_m=p.position.tolist(),admit_time_s=event['birth_time_s'])
        def provider(i,shape,t):
            nonlocal calls,last_position,last_sample,last_checked_count,no_acceptance_checks
            calls+=1
            if calls>max_calls:raise IntegrationSafetyStop('DETERMINISTIC_PROVIDER_CALL_GUARD')
            if stepper is not None and calls%128==0:
                count=len(stepper.samples)
                no_acceptance_checks=no_acceptance_checks+1 if count==last_checked_count else 0
                last_checked_count=count
                if no_acceptance_checks>=4*guard_factor:raise IntegrationSafetyStop('REJECTED_TRIAL_PROGRESS_GUARD_512_EVALUATIONS_WITHOUT_ACCEPTED_STATE')
                if count>=progress_window:
                    recent=np.array(stepper.samples[-progress_window:]);positions=recent[:,1:4]
                    if np.linalg.norm(np.ptp(positions,axis=0))<=16*env.wall.roundoff_m:
                        raise IntegrationSafetyStop('ROUNDOFF_SCALE_STAGNATION_64_ACCEPTED_SUBSTEPS_NOT_PHYSIOLOGICAL_TRAPPING_PROOF')
                    if recent[-1,0]-recent[0,0]<dt*.01:
                        raise IntegrationSafetyStop('EXCESSIVE_SUBDIVISION_PROGRESS_GUARD_64_STEPS_BELOW_1_PERCENT_NOMINAL_DT')
            key=tuple(shape.center_m)
            if key!=last_position:last_position=key;last_sample=env.field.sample(shape.center_m)
            s=last_sample
            if not s.inside_lumen:raise IntegrationSafetyStop('CENTER_OUTSIDE_FROZEN_LUMEN')
            return s.velocity_m_s,.5*s.vorticity_s_inv
        query=ValidationNeighborPolicy(20e-6,.5e-6,'P81_INDEPENDENT_SINGLE_MB_NO_PAIR_QUERY')
        try:
            stepper=SavedTrajectoryStepper([p],query,provider,env.mu,wall=env.wall,boundary_classifier=env.classifier,
                physical_time_s=0. if checkpoint is None else float(before[-1,0]))
            unchanged=0;previous_position=p.position.copy()
            first_k=1
            if checkpoint is not None:
                from .particle82_checkpoint import restore_stepper_recording,pending_intervals,preceding_nominal_state,advance_pending_interval
                restore_stepper_recording(stepper,original,before)
                calls=original['provider_calls'];last_checked_count=len(before)
                pending,first_k=pending_intervals(before,dt)
                previous_position,unchanged=preceding_nominal_state(before,dt)
                result['continuation_checkpoint']['pending_intervals']=[list(row) for row in pending]
                result['continuation_checkpoint']['first_future_nominal_index']=first_k
                for begin,end,depth in pending:
                    if stepper.time_s!=begin:raise ValueError('Pending checkpoint interval starts at wrong saved time')
                    advance_pending_interval(stepper,end,depth)
                    if stepper.boundary_event.startswith('OUTLET_'):break
                if pending:
                    position=stepper.read()[0].position
                    unchanged=unchanged+1 if np.array_equal(position,previous_position) else 0
                    previous_position=position.copy()
            for k in range(first_k,int(np.ceil(HORIZON/dt))+1):
                if stepper.boundary_event.startswith('OUTLET_'):break
                stepper.step_to(min(k*dt,HORIZON))
                if stepper.boundary_event.startswith('OUTLET_'):break
                position=stepper.read()[0].position
                unchanged=unchanged+1 if np.array_equal(position,previous_position) else 0
                previous_position=position.copy()
                if unchanged>=stationary_steps:
                    raise IntegrationSafetyStop('EXACT_STATIONARY_COORDINATES_32_NOMINAL_STEPS_NOT_PHYSIOLOGICAL_TRAPPING_PROOF')
            result['end_reason']=stepper.boundary_event if stepper.boundary_event!='ACTIVE' else 'PHYSICAL_RESIDENCE_HORIZON_REACHED'
        except (IntegrationSafetyStop,PhysicalTimeRefinementError,ValueError) as error:
            result['end_reason']='INTEGRATION_SAFETY_STOP' if isinstance(error,IntegrationSafetyStop) else 'NUMERICAL_REFINEMENT_FAILURE'
            result['failure_detail']=json.loads(json.dumps(getattr(error,'record',str(error)),default=json_safe))
        samples=np.array(stepper.samples,dtype=float) if stepper is not None else np.empty((0,len(COLUMNS)))
        if checkpoint is not None and (len(samples)<len(before) or samples[:len(before)].tobytes()!=before.tobytes()):
            raise ValueError('True checkpoint continuation failed exact historical prefix preservation')
        if stepper is not None:
            result.update(accepted_steps=stepper.accepted_count,rejected_trials=stepper.rejected_count,
                maximum_refinement_depth=stepper.maximum_depth,minimum_original_wall_gap_m=stepper.minimum_gap,
                maximum_position_identity_error_m=stepper.maximum_position_identity_error,nearfield_states=dict(stepper.states))
            result['query_cache']=dict(stepper.cache_info()._asdict(),role='EXACT_INPUT_MEMOIZATION_NO_NUMERICAL_CHANGE',adapter_revision='V2_EXACT_BVH_AND_WALL_GAP')
        if len(samples):
            result.update(last_elapsed_time_s=float(samples[-1,0]),last_physical_time_s=float(event['birth_time_s']+samples[-1,0]),
                path_length_m=float(np.linalg.norm(np.diff(samples[:,1:4],axis=0),axis=1).sum()))
        if result['end_reason'].startswith('OUTLET_'):
            # Recheck original final Euler segment against the official outlet triangles.
            hit=env.classifier.first_event(samples[-2,1:4],samples[-1,1:4])
            if hit is None or hit.role!=result['end_reason']:raise ValueError('Terminal boundary could not be independently reclassified')
            result.update(completed=True,exit_outlet=hit.role,exit_time_s=event['birth_time_s']+float(samples[-1,0]),
                residence_time_s=float(samples[-1,0]),terminal_event=dict(role=hit.role,role_triangle_id=hit.role_triangle_id,
                    position_m=samples[-1,1:4].tolist(),segment_start_m=samples[-2,1:4].tolist(),
                    segment_fraction=hit.segment_fraction,semantics='FIRST_CENTER_CROSSING_OF_OFFICIAL_OPEN_OUTLET; NO_FINITE_BODY_CLEARANCE_CLAIM'))
    result.update(sample_count=len(samples),provider_calls=calls,wall_seconds=time.perf_counter()-start,
        computational_guard_factor=guard_factor,
        nominal_no_motion_count=locals().get('unchanged',0),
        rejected_provider_check_count=no_acceptance_checks,
        progress_window_accepted_states=progress_window)
    if len(samples):
        last=env.field.sample(samples[-1,1:4])
        recent=samples[-min(progress_window,len(samples)):]
        result['stop_diagnostic']=dict(position_m=samples[-1,1:4].tolist(),tetra_id=last.tetra_id,
            elapsed_time_s=float(samples[-1,0]),physical_time_s=float(event['birth_time_s']+samples[-1,0]),
            nearest_wall_gap_m=float(samples[-1,14]),g_nf_m=float(samples[-1,15]),h_lower_m=float(samples[-1,16]),
            speed_m_s=float(np.linalg.norm(samples[-1,4:7])),nearfield_state=STATES[int(samples[-1,17])],
            maximum_subdivision_depth=result.get('maximum_refinement_depth'),
            last_accepted_displacement_m=float(np.linalg.norm(samples[-1,1:4]-samples[-2,1:4])) if len(samples)>1 else 0.,
            recent_position_span_m=float(np.linalg.norm(np.ptp(recent[:,1:4],axis=0))),
            nominal_no_motion_count=locals().get('unchanged',0),provider_calls=calls,
            roundoff_motion_count=sum(np.linalg.norm(np.diff(recent[:,1:4],axis=0),axis=1)<=16*env.wall.roundoff_m))
    temp=array_path.with_suffix('.tmp.npz');np.savez_compressed(temp,samples=samples);os.replace(temp,array_path)
    result['samples_sha256']=digest(array_path);result['samples_path']=str(array_path.relative_to(output))
    temp=meta_path.with_suffix('.tmp.json');dump(temp,result);os.replace(temp,meta_path)
    return result
