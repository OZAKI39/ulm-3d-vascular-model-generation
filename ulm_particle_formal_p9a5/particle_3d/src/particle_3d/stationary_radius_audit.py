"""Bounded virtual-radius diagnostics. No production caller or admission change."""
from pathlib import Path
import time,hashlib
import numpy as np
from .routing_stationary_audit import read,sha,downstream_feasibility
from .particle81_simulation import environment,dump,DT,STATIONARY_STEPS,MAX_PROVIDER_CALLS,IntegrationSafetyStop
from .nearfield_regularization import NearFieldRegularizationV1
from .nearfield_handoff import handoff_constraints
from .particle_shapes import Sphere
from .particle65_motion import assemble_v1
from .particle9a_motion import augment_planar_system,Particle9AStepper
from .wall_gap import wall_gap
from .resistance_solver import solve_resistance
from .lammps_state import BridgeParticle
from .lammps_neighbors import ValidationNeighborPolicy
from .physical_time_refinement import PhysicalTimeRefinementError
from .convex_triangle import triangle_closest_many

DIAGNOSTIC_ROLE='DIAGNOSTIC_ONLY_VIRTUAL_RADIUS_NOT_REAL_BUBBLE_SHRINKAGE'

def state_audit(pid,radius,x):
    env=environment();field=env.field.sample(x);shapes={pid:Sphere(x,radius)};policy=NearFieldRegularizationV1()
    free={pid:np.r_[field.velocity_m_s,.5*field.vorticity_s_inv]}
    base=assemble_v1(shapes,free,env.mu,env.wall,policy=policy)
    system=augment_planar_system(base,shapes,env.mu,env.wall,lambda _:field.velocity_gradient_s_inv,wall_gap)
    cs,records=handoff_constraints(shapes,env.wall,policy);sol=solve_resistance(system,constraints=cs);unconstrained=solve_resistance(system)
    normals=np.array([c.normal for c in cs]);sv=np.linalg.svd(normals,compute_uv=False)
    tol=np.sqrt(max(normals.shape)*np.finfo(float).eps)*sv[0] if len(sv) else 0.
    red=sol.record['contact_redundancy']
    return dict(particle_id=pid,radius_m=radius,diameter_um=radius*2e6,final_position_m=np.asarray(x).tolist(),
        wall_gap_m=wall_gap(shapes[pid],env.wall).gap_m,contacts=records,contact_triangle_ids=[c.canonical_id[3] for c in cs],
        contact_feature_types=[c.canonical_id[4] for c in cs],contact_normals=normals.tolist(),contact_count=len(cs),
        normal_jacobian_singular_values=sv.tolist(),normal_jacobian_rank=int(np.sum(sv>tol)),normal_rank_tolerance=tol,
        contact_condition=sol.record['contact_kkt']['condition'],contact_rank=red['rank_after'],contact_redundancy=red,
        retained_rows_independent=red['rank_after']==red['retained_constraint_count'],multipliers=sol.record['multipliers'],
        free_fem_velocity_m_s=field.velocity_m_s.tolist(),unconstrained_particle_velocity=unconstrained.velocity.tolist(),constrained_velocity=sol.velocity.tolist(),
        solver=sol.record,downstream=downstream_feasibility(normals,field.velocity_m_s))

def checkpoint(samples,radius,roundoff,wall=None,pid=None):
    displacement=np.linalg.norm(np.diff(samples[:,1:4],axis=0),axis=1)
    # A directly observed outgoing movement, larger than the original stagnation budget.
    threshold=max(16*roundoff,np.sqrt(np.finfo(float).eps)*radius)
    moving=np.flatnonzero(displacement>threshold)
    if wall is not None:
        for candidate in moving[::-1]:
            cs,_=handoff_constraints({pid:Sphere(samples[candidate,1:4],radius)},wall,NearFieldRegularizationV1())
            normals=np.array([c.normal for c in cs])
            if len(cs)<3 or np.linalg.matrix_rank(normals)<3:
                moving=moving[moving<=candidate];break
        else:raise ValueError('No movable checkpoint before rank-three contacts')
    if not len(moving):raise ValueError('No demonstrably moving checkpoint')
    k=int(moving[-1]);row=samples[k].copy()
    return k,row,dict(rule='LAST_SAVED_STATE_BEFORE_RANK_THREE_CONTACTS_WITH_OBSERVED_NEXT_DISPLACEMENT_GT_MAX_16_ROUNDOFF_SQRT_EPS_RADIUS',index=k,
        observed_next_displacement_m=float(displacement[k]),movement_threshold_m=threshold,
        elapsed_time_s=float(row[0]),position_m=row[1:4].tolist(),quaternion=row[10:14].tolist(),
        prefix_sha256=hashlib.sha256(samples[:k+1].tobytes()).hexdigest(),original_radius_m=radius)

def virtual_trial(pid,radius,ratio,row,final_position,folder):
    env=environment();start=time.time();a=radius*ratio;shape=Sphere(row[1:4],a)
    p=BridgeParticle.from_shape(pid,shape,q=row[10:14],velocity=row[4:7],omega=row[7:10])
    calls=0;stepper=None;last_checked=0;no_accept=0
    def provider(i,s,t):
        nonlocal calls,last_checked,no_accept
        calls+=1
        if calls>MAX_PROVIDER_CALLS:raise IntegrationSafetyStop('DETERMINISTIC_PROVIDER_CALL_GUARD')
        if stepper is not None and calls%128==0:
            n=len(stepper.samples);no_accept=no_accept+1 if n==last_checked else 0;last_checked=n
            if no_accept>=4:raise IntegrationSafetyStop('REJECTED_TRIAL_PROGRESS_GUARD_512_EVALUATIONS_WITHOUT_ACCEPTED_STATE')
            if n>=64:
                recent=np.array(stepper.samples[-64:])
                if np.linalg.norm(np.ptp(recent[:,1:4],axis=0))<=16*env.wall.roundoff_m:raise IntegrationSafetyStop('ROUNDOFF_SCALE_STAGNATION_64_ACCEPTED_SUBSTEPS')
                if recent[-1,0]-recent[0,0]<DT*.01:raise IntegrationSafetyStop('EXCESSIVE_SUBDIVISION_PROGRESS_GUARD')
        f=env.field.sample(s.center_m)
        if not f.inside_lumen:raise IntegrationSafetyStop('CENTER_OUTSIDE_FROZEN_LUMEN')
        return f.velocity_m_s,.5*f.vorticity_s_inv
    field=env.field.sample(final_position);uhat=field.velocity_m_s/np.linalg.norm(field.velocity_m_s)
    # Local region: three original radii around final jam; use accumulated local-FEM
    # downstream displacement through bends, not a fixed global-direction projection.
    duration=min(.1,max(64*DT,20*radius/np.linalg.norm(field.velocity_m_s)))
    original_contacts=handoff_constraints({pid:Sphere(final_position,radius)},env.wall,NearFieldRegularizationV1())[0]
    original_triangles=env.wall.triangles[[c.canonical_id[3] for c in original_contacts]]
    criteria=dict(local_time_budget_s=duration,hotspot_radius_m=3*radius,downstream_distance_m=3*radius,
        passage_criterion_revision='LOCAL_CURVED_FLOW_PASSAGE_V2',
        downstream_measure='SUM_OF_NOMINAL_DISPLACEMENT_DOT_LOCAL_FEM_UNIT_DIRECTION',
        require_all_original_contact_triangles_cleared=True,
        consecutive_downstream_steps=5,local_only=True,ratios_order=[1.,.99,.95,.90],original_nominal_grid_preserved=True,
        original_provider_and_progress_guards=True)
    status='LOCAL_TIME_BUDGET_UNRESOLVED';reason=None;unchanged=0;downstream_steps=0;cumulative_downstream=0.
    try:
        stepper=Particle9AStepper([p],ValidationNeighborPolicy(20e-6,.5e-6,'DIAGNOSTIC_ONLY_SINGLE_MB'),provider,env.mu,
            wall=env.wall,boundary_classifier=env.classifier,physical_time_s=float(row[0]),
            gradient_provider=lambda x:env.field.sample(x).velocity_gradient_s_inv)
        previous=p.position.copy();end=float(row[0])+duration
        first=int(np.floor(float(row[0])/DT))+1
        for k in range(first,int(np.ceil(end/DT))+1):
            target=min(k*DT,end)
            if target<=stepper.time_s:continue
            stepper.step_to(target);x=stepper.read()[0].position;delta=x-previous
            unchanged=unchanged+1 if np.array_equal(x,previous) else 0
            s=env.field.sample(previous);downstream_steps=downstream_steps+1 if float(delta@s.velocity_m_s)>0 and np.linalg.norm(delta)>16*env.wall.roundoff_m else 0
            cumulative_downstream+=float(delta@s.velocity_m_s/np.linalg.norm(s.velocity_m_s)) if np.linalg.norm(s.velocity_m_s) else 0.
            previous=x.copy()
            offset=x-final_position;cleared=False
            if cumulative_downstream>3*radius:
                near,_=triangle_closest_many(x,original_triangles)
                cleared=bool(np.all(np.linalg.norm(x-near,axis=1)-a>2e-9+16*env.wall.roundoff_m))
            if cleared and np.linalg.norm(offset)>3*radius and downstream_steps>=5:
                status='PASSED_LOCAL_HOTSPOT';break
            if stepper.boundary_event.startswith('OUTLET_'):
                status='PASSED_LOCAL_HOTSPOT' if cleared else 'OUTLET_BEFORE_LOCAL_CRITERION_UNRESOLVED';break
            if unchanged>=STATIONARY_STEPS:status='STATIONARY';reason='EXACT_STATIONARY_COORDINATES_32_NOMINAL_STEPS';break
    except (IntegrationSafetyStop,PhysicalTimeRefinementError,ValueError) as err:
        reason=str(err);status='STATIONARY' if 'STAGNATION' in reason else 'NUMERICAL_OR_GEOMETRIC_UNRESOLVED'
    samples=np.array(stepper.samples if stepper is not None else [row]);folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    file=folder/f'mb_{pid:06d}_ratio_{round(ratio*100):03d}.npz';np.savez_compressed(file,samples=samples)
    x=samples[-1,1:4];result=dict(role=DIAGNOSTIC_ROLE,particle_id=pid,radius_ratio=ratio,original_radius_m=radius,virtual_radius_m=a,
        status=status,reason=reason,criteria=criteria,checkpoint_time_s=float(row[0]),last_time_s=float(samples[-1,0]),
        final_position_m=x.tolist(),distance_from_original_jam_m=float(np.linalg.norm(x-final_position)),downstream_distance_m=float((x-final_position)@uhat),
        cumulative_local_downstream_m=cumulative_downstream,provider_calls=calls,wall_seconds=time.time()-start,trajectory_sha256=sha(file),trajectory_path=str(file),production_result=False)
    if status=='STATIONARY':
        try:result['final_contact_audit']=state_audit(pid,a,x)
        except Exception as err:result['final_contact_audit_error']=str(err)
    dump(file.with_suffix('.json'),result);return result

def radius_job(args):
    meta_path,folder=args;meta_path=Path(meta_path);m=read(meta_path);pid=m['particle_id'];radius=m['radius_m'];samples=np.load(meta_path.with_suffix('.npz'))['samples'];final=samples[-1,1:4]
    env=environment();audit=state_audit(pid,radius,final);k,row,cp=checkpoint(samples,radius,env.wall.roundoff_m,env.wall,pid)
    cp['contacts_at_checkpoint']=handoff_constraints({pid:Sphere(row[1:4],radius)},env.wall,NearFieldRegularizationV1())[1]
    cp['absolute_birth_time_s']=m['birth_time_s'];cp['absolute_checkpoint_time_s']=m['birth_time_s']+row[0]
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True);np.savez_compressed(folder/f'prefix_{pid:06d}.npz',samples=samples[:k+1])
    trials=[]
    for ratio in [1.,.99,.95,.90]:
        trial=virtual_trial(pid,radius,ratio,row,final,folder);trials.append(trial)
        if ratio==1.:
            trial['reference_jam_reproduced']=trial['status']=='STATIONARY' and trial['distance_from_original_jam_m']<=16*env.wall.roundoff_m
        elif trial['status']=='PASSED_LOCAL_HOTSPOT':break
    reference_ok=trials[0]['reference_jam_reproduced'];passed=[t for t in trials if t['radius_ratio']<1 and t['status']=='PASSED_LOCAL_HOTSPOT']
    unresolved=any('UNRESOLVED' in t['status'] for t in trials)
    critical=None;description='UNRESOLVED';category='NUMERICAL_OR_GEOMETRIC_UNRESOLVED'
    if passed and reference_ok and not unresolved:
        passing=passed[0]['radius_ratio'];previous=trials[trials.index(passed[0])-1]['radius_ratio'];critical=dict(highest_tested_passing=passing,lowest_tested_failing=previous,interval=[passing,previous],role='LOCAL_PASSAGE_BRACKET_NOT_EXACT_CRITICAL_RADIUS; MONOTONICITY_NOT_PROVEN')
        description='HIGHLY_RADIUS_SENSITIVE' if passing>=.95 else 'INTERMEDIATE'
        if audit['retained_rows_independent'] and audit['downstream']['feasible_downstream_direction'] is False:
            category='RIGID_MODEL_HIGH_SENSITIVITY' if passing>=.95 else 'RIGID_GEOMETRIC_JAM'
    elif reference_ok and not unresolved and all(t['status']=='STATIONARY' and t.get('final_contact_audit',{}).get('downstream',{}).get('feasible_downstream_direction') is False for t in trials):
        description='ROBUST_RIGID_SIZE_EXCLUSION';critical=dict(highest_tested_passing=None,lowest_tested_failing=.9,interval=None,role='NO_PASSAGE_OBSERVED_AT_RATIOS_1_099_095_090')
        if audit['retained_rows_independent'] and audit['downstream']['feasible_downstream_direction'] is False:category='RIGID_GEOMETRIC_JAM'
    result=dict(audit,role=DIAGNOSTIC_ROLE,checkpoint=cp,trials=trials,reference_jam_reproduced=reference_ok,
        critical_radius_ratio=critical,radius_sensitivity=description,classification=category)
    dump(folder/f'stationary_{pid:06d}.json',result);return result
