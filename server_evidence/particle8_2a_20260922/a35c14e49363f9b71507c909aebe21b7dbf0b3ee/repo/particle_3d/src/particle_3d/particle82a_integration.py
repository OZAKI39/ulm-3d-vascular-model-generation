"""P8.2A explicit accepted birth state, unchanged P8.1 integration body.

Only initialization/metadata differ. Physics uses the original SavedTrajectoryStepper,
P6.5 force/geometry/refinement and identical nominal timestep/guard ordering.
The source-generation diff and a common-state parity test audit this separation.
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

def integrate_admitted(event,dt=DT,output=None,force=False):
    if output is None: raise ValueError("Explicit P8.2A output required")
    output=Path(output);folder=output/'trajectories';folder.mkdir(parents=True,exist_ok=True)
    pid=event['particle_id'];stem=f'mb_{pid:06d}';meta_path=folder/(stem+'.json');array_path=folder/(stem+'.npz')
    config=dict(dt_s=dt,horizon_s=HORIZON,maximum_provider_calls=MAX_PROVIDER_CALLS,
        roundoff_progress_window_accepted_steps=64,roundoff_progress_budget_multiplier=16,
        minimum_64_substep_progress_nominal_dt_fraction=.01,
        progress_check_provider_interval=128,maximum_rejected_provider_calls_without_acceptance=512,
        stationary_unchanged_nominal_steps=STATIONARY_STEPS,nearfield_contract='NEAR_FIELD_REGULARIZATION_V1')
    if meta_path.exists() and not force:
        existing=read(meta_path)
        if existing['birth_metadata_sha256']!=canonical_hash(event) or existing['integration_config']!=config or digest(array_path)!=existing['samples_sha256']:
            raise ValueError('Resume data mismatch for stable ID '+str(pid))
        return existing
    env=environment();metadata=deepcopy(event);start=time.perf_counter()
    checker=FiniteSizeAdmission(None,None,wall=env.wall,field=env.field)
    p,status,detail=checker.check(metadata,np.asarray(event['birth_center_m']),{})
    if p is None: raise ValueError('Saved admission failed original checker at trajectory start: '+status)
    admission=SimpleNamespace(candidates=[])
    rng_initial=None; rng=None
    result=dict(particle_id=pid,birth_metadata=event,birth_metadata_sha256=canonical_hash(event),integration_config=config,
        birth_time_s=event['birth_time_s'],radius_m=event['radius_m'],diameter_um=event['diameter_um'],initial_q=event['q'],
        position_rng_initial=None,position_rng_final=None,admission_candidates=[],
        admission_strategy=event['admission_strategy'],anchor_m=event['anchor_m'],
        entry_transition_time_s=event.get('entry_transition_time_s',0.),
        entry_transition_role='POINT_STREAMLINE_REPRESENTATION_TIME_NOT_FINITE_SIZE_ENTRY_DYNAMICS',
        original_geometry_preserved=True,admission_input_preserved=True,independently_integrated=True,
        dynamics='UNCHANGED_PARTICLE65STEPPER_NO_LAMMPS_NO_MB_MB_OR_RBC_COUPLING',
        timestep_role='P81_DATASET_VALIDATION_ONLY_NOT_PRODUCTION_FROZEN',columns=COLUMNS,
        velocity_role='HELD_ON_PRECEDING_ACCEPTED_INTERVAL; FIRST_ROW_FREE_INITIAL_DIAGNOSTIC',
        terminal_event=None,completed=False,exit_outlet=None,exit_time_s=None,residence_time_s=None,path_length_m=0.,failure_detail=None)
    calls=0;stepper=None;last_checked_count=0;no_acceptance_checks=0
    last_position=None;last_sample=None
    if p is None:
        result['end_reason']='INLET_ADMISSION_GUARD_EXHAUSTED_UNRESOLVED'
        samples=np.empty((0,len(COLUMNS)))
    else:
        result.update(inlet_face=event['anchor_triangle'],
            inlet_position_m=event['anchor_m'],birth_center_m=event['birth_center_m'],
            admit_time_s=event['birth_time_s'])
        def provider(i,shape,t):
            nonlocal calls,last_position,last_sample,last_checked_count,no_acceptance_checks
            calls+=1
            if calls>MAX_PROVIDER_CALLS:raise IntegrationSafetyStop('DETERMINISTIC_PROVIDER_CALL_GUARD')
            if stepper is not None and calls%128==0:
                count=len(stepper.samples)
                no_acceptance_checks=no_acceptance_checks+1 if count==last_checked_count else 0
                last_checked_count=count
                if no_acceptance_checks>=4:raise IntegrationSafetyStop('REJECTED_TRIAL_PROGRESS_GUARD_512_EVALUATIONS_WITHOUT_ACCEPTED_STATE')
                if count>=64:
                    recent=np.array(stepper.samples[-64:]);positions=recent[:,1:4]
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
            stepper=SavedTrajectoryStepper([p],query,provider,env.mu,wall=env.wall,boundary_classifier=env.classifier)
            unchanged=0;previous_position=p.position.copy()
            for k in range(1,int(np.ceil(HORIZON/dt))+1):
                stepper.step_to(min(k*dt,HORIZON))
                if stepper.boundary_event.startswith('OUTLET_'):break
                position=stepper.read()[0].position
                unchanged=unchanged+1 if np.array_equal(position,previous_position) else 0
                previous_position=position.copy()
                if unchanged>=STATIONARY_STEPS:
                    raise IntegrationSafetyStop('EXACT_STATIONARY_COORDINATES_32_NOMINAL_STEPS_NOT_PHYSIOLOGICAL_TRAPPING_PROOF')
            result['end_reason']=stepper.boundary_event if stepper.boundary_event!='ACTIVE' else 'PHYSICAL_RESIDENCE_HORIZON_REACHED'
        except (IntegrationSafetyStop,PhysicalTimeRefinementError,ValueError) as error:
            result['end_reason']='INTEGRATION_SAFETY_STOP' if isinstance(error,IntegrationSafetyStop) else 'NUMERICAL_REFINEMENT_FAILURE'
            result['failure_detail']=json.loads(json.dumps(getattr(error,'record',str(error)),default=json_safe))
        samples=np.array(stepper.samples,dtype=float) if stepper is not None else np.empty((0,len(COLUMNS)))
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
    result.update(sample_count=len(samples),provider_calls=calls,wall_seconds=time.perf_counter()-start)
    temp=array_path.with_suffix('.tmp.npz');np.savez_compressed(temp,samples=samples);os.replace(temp,array_path)
    result['samples_sha256']=digest(array_path);result['samples_path']=str(array_path.relative_to(output))
    temp=meta_path.with_suffix('.tmp.json');dump(temp,result);os.replace(temp,meta_path)
    return result

