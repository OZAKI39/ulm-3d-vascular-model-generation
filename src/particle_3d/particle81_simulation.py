"""Independent full-vessel MB integrations using the unchanged P6.5 stepper.

No streamline substitution, position projection, altered radii or new forces.
Physical integration uses elapsed time; absolute acquisition clocks are saved
separately to avoid loss of substep precision after hours of acquisition.
"""
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from collections import Counter,namedtuple
import json,os,subprocess,time,multiprocessing
import numpy as np
from .particle8_replay import REPO,read,write,digest,canonical_hash
from .particle3_cases import json_safe,write_json
from .particle7_cases import real_inlet,SONOVUE
from .injection_population import PopulationSource,FluxClock,LinearProfile,ConstantMBConcentrationV0,H_D,C_MB
from .injection_admission import FiniteSizeAdmission
from .validation_boundary import ValidationBoundaryClassifier
from .particle65_motion import Particle65Stepper
from .lammps_neighbors import ValidationNeighborPolicy
from .nearfield_regularization import NearFieldRegularizationV1,STATES
from .nearfield_handoff import exact_interactions
from .hydrodynamic_resistance import viscosity_from_frozen
from .physical_time_refinement import PhysicalTimeRefinementError
from .particle81_cache import cached_step_to

OUTPUT=REPO/'particle_3d/outputs/particle8_1'
SEED=2026092181
DT=.00025 # P7 single-MB validation step; not a production timestep decision.
HORIZON=1.5 # Per-MB elapsed physical-time safety limit, not acquisition duration.
MAX_PROVIDER_CALLS=16000 # Deterministic computation guard; incomplete paths retained.
STATIONARY_STEPS=32 # Exact unchanged centers, not a physical speed cutoff.
_ENV=None


def dump(path,value):write_json(path,value)


def environment():
    global _ENV
    if _ENV is None:
        provenance,mesh,field,boundaries,audit,samplers,wall=real_inlet()
        classifier=ValidationBoundaryClassifier({k:v for k,v in boundaries.items() if k.startswith('OUTLET_')})
        mu,viscosity=viscosity_from_frozen(REPO/'formal_3D_flow_solver/FEM_SimVascular')
        _ENV=SimpleNamespace(provenance=provenance,mesh=mesh,field=field,boundaries=boundaries,audit=audit,
            sampler=samplers['INLET'],wall=wall,classifier=classifier,mu=mu,viscosity=viscosity)
    return _ENV


def lock_upstream():
    path=OUTPUT/'data/upstream_lock.json'
    if path.exists():
        data=read(path)
        for name,sha in data['sha256'].items():
            if digest(REPO/name)!=sha:raise ValueError('Upstream changed: '+name)
        return len(data['sha256'])
    names=subprocess.check_output(['git','ls-files','particle_3d'],cwd=REPO,text=True).splitlines()
    names=[p for p in names if 'particle81' not in p and 'PARTICLE8_1' not in p and '/particle8_1/' not in p]
    dump(path,dict(baseline_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
        sha256={name:digest(REPO/name) for name in names}))
    return len(names)


def prepare_population(count):
    env=environment();lock_upstream();source=PopulationSource(SONOVUE,SEED)
    clock=FluxClock(LinearProfile([0],[env.sampler.Q_m3_s]),ConstantMBConcentrationV0())
    events=[]
    for pid in range(1,count+1):
        event=source.next_mb()
        # P7 Haar SO(3) sampler; no RBC are created in this independent MB dataset.
        event['q']=source.orientation()[0].tolist()
        event.update(particle_id=pid,attempt_count=0,birth_time_s=clock.time_at(pid),scheduled_time_s=clock.time_at(pid),
            orientation_role='ISOTROPIC_RANDOM_ORIENTATION_V0_MODEL_ASSUMPTION_NOT_MEASURED',
            position_seed=[SEED,pid,81])
        events.append(event)
    result=dict(schema='PARTICLE81_BIRTH_LEDGER_V1',events=events,seed=SEED,source_rng_final=source.state(),
        Q_in_m3_s=env.sampler.Q_m3_s,C_MB_m3=C_MB,H_D_feed=H_D,
        H_D_role='RETAINED_FEED_PARAMETER_NO_RBC_SIMULATION_OR_FORCED_TUBE_HCT',
        scheduler='P7_FLUX_CLOCK_DETERMINISTIC_CUMULATIVE_INTEGER_THRESHOLDS',
        acquisition_birth_window_s=events[-1]['birth_time_s'],nominal_expected_births=count,
        independent_MB_position_streams='PCG64 SeedSequence([master_seed, stable_id, 81]); same P7 exact flux sampler',
        geometry_sampling='UNCHANGED_SONOVUE_INVERSE_CDF_ONCE_PER_ID_NO_DIAMETER_REJECTION_RESAMPLING')
    old=OUTPUT/'data/birth_ledger.json'
    if old.exists():
        before=read(old)['events']
        if before!=events[:len(before)]:raise ValueError('Extending the acquisition must preserve every existing birth')
    dump(old,result)
    dump(OUTPUT/'data/frozen_provenance.json',dict(provenance=env.provenance,viscosity=env.viscosity,flux_audit=env.audit))
    dump(OUTPUT/'data/inlet_mesh.json',dict(triangles_m=env.sampler.triangles,normal_velocity_m_s=env.sampler.q,
        weights_m3_s=env.sampler.weights,expected_mean_position_m=env.sampler.expectation()))
    arrays={}
    for role,s in env.boundaries.items():
        arrays[role+'_points_m']=np.asarray(s.points,dtype=float);arrays[role+'_faces']=s.faces.reshape(-1,4)[:,1:]
    np.savez_compressed(OUTPUT/'data/full_frozen_geometry.npz',**arrays)
    return events


class IntegrationSafetyStop(RuntimeError):pass


class SavedTrajectoryStepper(Particle65Stepper):
    """Record every upstream accepted substep; discard only verbose duplicate audits."""
    def __init__(self,*args,**kwargs):
        memoize=kwargs.pop('memoize',True)
        super().__init__(*args,**kwargs)
        self.memoize=memoize
        if memoize:self.advance_cached,self.cache_info=cached_step_to(self)
        else:
            self.advance_cached=None
            info=namedtuple('CacheInfo','hits misses maxsize currsize')
            self.cache_info=lambda:info(0,0,0,0)
        self.samples=[];self.accepted_count=0;self.rejected_count=0;self.maximum_depth=0
        self.minimum_gap=np.inf;self.maximum_position_identity_error=0.;self.states=Counter()
        p=self.read()[0]
        rows=exact_interactions({p.particle_id:p.shape()},self.wall,self.regularization,self.mu)[1]
        self.record_sample(p,rows[0],0.,0.)

    def record_sample(self,p,interaction,dt,evaluation_time):
        state=interaction['interaction_state'];self.states[state]+=1
        gap=interaction['h_geom_m'];self.minimum_gap=min(self.minimum_gap,gap)
        self.samples.append([self.time_s,*p.position,*p.velocity,*p.omega,*p.q,gap,
            interaction['g_nf_m'],interaction['h_lower_m'],STATES.index(state),dt,evaluation_time])

    def _accepted(self,world):
        previous=self.read()[0];dt=world.time_s-self.time_s
        super()._accepted(world)
        p=self.read()[0];interaction=world.projection['accepted_interactions'][0]
        residual=np.max(np.abs(p.position-(previous.position+dt*p.velocity)))
        self.maximum_position_identity_error=max(self.maximum_position_identity_error,float(residual))
        self.record_sample(p,interaction,dt,world.projection['evaluated_at_time_s'])
        self.accepted_count+=1;self.accepted_worlds.clear()

    def step_to(self,end_time):
        try:return self.advance_cached(end_time) if self.memoize else super().step_to(end_time)
        finally:
            self.rejected_count+=sum(not a['accepted'] for a in self.attempts)
            self.maximum_depth=max([self.maximum_depth]+[a['depth'] for a in self.ledger])
            self.attempts.clear();self.ledger.clear()


COLUMNS=['elapsed_time_s','x_m','y_m','z_m','vx_m_s','vy_m_s','vz_m_s','omega_x_s_inv','omega_y_s_inv','omega_z_s_inv',
         'q_w','q_x','q_y','q_z','h_geom_m','g_nf_m','h_lower_m','nearfield_state_code','accepted_dt_s','velocity_evaluated_elapsed_s']


def integrate_one(event,dt=DT,output=OUTPUT,force=False):
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
    rng=np.random.default_rng(event['position_seed']);rng_initial=deepcopy(rng.bit_generator.state)
    source=SimpleNamespace(rng={'MB_POSITION':rng,'ADMISSION_RETRY':rng})
    admission=FiniteSizeAdmission(env.sampler,source,wall=env.wall,field=env.field)
    p=admission.attempt(metadata,{})
    result=dict(particle_id=pid,birth_metadata=event,birth_metadata_sha256=canonical_hash(event),integration_config=config,
        birth_time_s=event['birth_time_s'],radius_m=event['radius_m'],diameter_um=event['diameter_um'],initial_q=event['q'],
        position_rng_initial=rng_initial,position_rng_final=deepcopy(rng.bit_generator.state),admission_candidates=admission.candidates,
        original_geometry_preserved=True,artificial_inlet_offset_m=0.,independently_integrated=True,
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
        xyz=env.sampler.triangles[metadata['last_candidate']['triangle_id']]
        uv=np.linalg.lstsq((xyz[1:]-xyz[0]).T,p.position-xyz[0],rcond=None)[0]
        result.update(inlet_face=metadata['last_candidate']['triangle_id'],inlet_barycentric=[1.-float(uv.sum()),*uv.tolist()],
            inlet_position_m=p.position.tolist(),admit_time_s=event['birth_time_s'])
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


def run(count,workers=8,first_id=1):
    if not 1<=first_id<=count:raise ValueError("Invalid scheduled ID range")
    events=prepare_population(count)[first_id-1:];completed=Counter();results=[]
    # Each worker owns its particle and RNG; immutable Frozen inputs inherited.
    # Process parallelism changes neither stable IDs nor birth/acceptance ordering.
    context=multiprocessing.get_context('fork')
    with context.Pool(workers) as pool:
        for result in pool.imap_unordered(integrate_one,events,chunksize=1):
            results.append(result);completed[result['end_reason']]+=1
            if len(results)%10==0 or len(results)==len(events):
                print('PROGRESS',len(results),'/',len(events),dict(completed),flush=True)
    results.sort(key=lambda r:r['particle_id'])
    summary_name='run_summary.json' if first_id==1 else f'run_partition_{first_id}_{count}.json'
    dump(OUTPUT/'data'/summary_name,dict(scheduled_total=count,first_id=first_id,scheduled=len(events),end_reasons=dict(completed),completed=sum(r['completed'] for r in results),
        workers=workers,algorithm='UNCHANGED_P65_STEPPER_INDEPENDENT_TRAJECTORIES',
        acquisition_birth_window_s=events[-1]['birth_time_s'],maximum_provider_calls=MAX_PROVIDER_CALLS,
        maximum_particle_residence_horizon_s=HORIZON,nominal_dt_s=DT))
    return results
