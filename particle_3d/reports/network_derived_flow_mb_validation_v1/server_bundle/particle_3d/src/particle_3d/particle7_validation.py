"""Reproducible P7 evidence generation. All populations are validation-only."""
from pathlib import Path
from copy import deepcopy
from dataclasses import replace
import json,math
import numpy as np
from .particle7_cases import *
from .particle3_cases import write_json,write_rows,json_safe
from .injection_population import C_MB,H_D,InjectionScheduler,PopulationSource,LinearProfile
from .injection_admission import FiniteSizeAdmission,PlugPathAdmission,event_shape
from .rbc_orientation import short_axis
from .particle_shapes import Sphere
from .lammps_state import BridgeParticle
from .particle7_checkpoint import write_checkpoint,read_checkpoint
from .particle65_motion import Particle65Stepper
from .hydrodynamic_resistance import viscosity_from_frozen


def plain(x): return json.loads(json.dumps(x,default=json_safe))


def scheduler_validation(folder,Q):
    source=PopulationSource(SONOVUE); s=InjectionScheduler(source,LinearProfile([0],[Q])); start=s.state()
    horizon=2.25/(C_MB*Q); events=[]; rows=[]
    for t in np.linspace(0,horizon,301):
        events.extend(s.through(float(t)))
        rows.append(dict(time_s=t,mb_expected=s.mb_clock.cumulative(t),mb_scheduled=s.mb_count,
            mb_error=s.mb_count-s.mb_clock.cumulative(t),rbc_target_m3=s.rbc_clock.cumulative(t),rbc_scheduled_m3=s.rbc_volume,
            rbc_residual_m3=s.rbc_clock.cumulative(t)-s.rbc_volume,next_rbc_volume_m3=s.next_rbc['volume_m3']))
    write_rows(folder/'03_04_scheduler.csv',rows); write_json(folder/'03_04_scheduled_events.json',events)
    partitions=[]
    for count in [1,2,4]:
        p=InjectionScheduler.restore(start); these=[]
        for k in range(1,count+1): these.extend(p.through(horizon*k/count))
        partitions.append(dict(dt_s=horizon/count,event_count=len(these),exact_equal=these==events))
    assert all(v['exact_equal'] for v in partitions)
    write_json(folder/'08_timestep_parity.json',dict(title='NOT PRODUCTION TIMESTEP SELECTION',partitions=partitions,events=events))
    return rows


def sampler_validation(folder,sampler):
    examples=[]
    for name,values in [('constant',[1,1,1]),('linear',[0,1,0]),('sign_changing',[-1,1,1])]:
        s=InletFluxSampler([[[0,0,0],[1,0,0],[0,1,0]]],[values]); p,ids=s.sample(np.random.default_rng(2026092107),100000)
        np.savetxt(folder/f'02_synthetic_{name}.csv',p,delimiter=',',header='x,y,z',comments='')
        examples.append(dict(case=name,Q=s.Q_m3_s,expected_mean=s.expectation(),sample_mean=p.mean(0),n=len(p),
            sampling_method='EXACT_LINEAR_DENSITY_DIRICHLET_MIXTURE',proposal_acceptance=1.))
    actual=np.loadtxt(folder/'02_actual_inlet_samples.csv',delimiter=',',skiprows=1)
    counts=np.bincount(actual[:,3].astype(int),minlength=len(sampler.weights)); expected=sampler.weights/sampler.Q_m3_s
    from scipy.stats import chi2
    sel=expected*len(actual)>=5; statistic=float(np.sum((counts[sel]-expected[sel]*len(actual))**2/(expected[sel]*len(actual))))
    write_json(folder/'02_sampling_statistics.json',dict(synthetic=examples,actual=dict(n=len(actual),expected_mean_m=sampler.expectation(),sample_mean_m=actual[:,:3].mean(0),
        chi_square=statistic,chi_square_dof=int(sel.sum()-1),chi_square_p_value=float(chi2.sf(statistic,sel.sum()-1)),acceptance=1.)))
    write_rows(folder/'02_actual_triangle_counts.csv',[dict(triangle_id=i,expected_fraction=expected[i],sampled_fraction=counts[i]/len(actual),sample_count=counts[i]) for i in range(len(expected))])
    source=PopulationSource(SONOVUE); q=source.orientation(100000); axes=np.array([short_axis(x) for x in q])
    np.savetxt(folder/'05_orientation.csv',np.column_stack([q,axes]),delimiter=',',header='qw,qx,qy,qz,px,py,pz',comments='')
    write_json(folder/'05_orientation_statistics.json',dict(n=len(q),mean_short_axis=axes.mean(0),second_moments=(axes**2).mean(0),
        model='ISOTROPIC_RANDOM_SO3_V0',role='PROVISIONAL_NOT_MEASURED_C57BL6_ORIENTATION',algorithm='NORMALIZED_4D_GAUSSIAN_HAAR_S3_QUOTIENT'))


def admission_validation(folder):
    s,w,c,v=channel(); src=PopulationSource(SONOVUE); a=FiniteSizeAdmission(s,src,wall=w,velocity=v)
    mb=src.next_mb(); mb.update(particle_id=1,attempt_count=0); rbc=src.next_rbc(); rbc.update(particle_id=2,attempt_count=0)
    records=[]
    for e in [mb,rbc]:
        for label,point,active in [('ACCEPT',[0,20e-6,0],{}),('WALL',[99.999e-6,0,0],{}),('PAIR',[0,0,0],{99:Sphere([0,0,0],3e-6)})]:
            p,status,detail=a.check(e,point,active)
            records.append(dict(species=e['species'],case=label,status=status,position_m=point,particle=e,accepted=p.to_dict() if p else None,**detail))
    write_json(folder/'06_admission.json',dict(width_m=200e-6,blocker_center_m=[0,0,0],blocker_radius_m=3e-6,cases=records))
    engine=synthetic_engine(bridge=False,guard=1); e=engine.scheduler.pop(); original=deepcopy(e)
    engine.events.append(deepcopy(e)); engine.pending['RBC'].append(e); engine.time_s=e['scheduled_time_s']
    peek=np.random.default_rng(); peek.bit_generator.state=deepcopy(engine.scheduler.source.rng['RBC_POSITION'].bit_generator.state)
    blocker_center=engine.admission.sampler.sample(peek)[0][0]
    blocker_radius=min(2e-6,float(np.min(100e-6-abs(blocker_center[:2])))/2)
    engine.active[99999]=BridgeParticle.from_shape(99999,Sphere(blocker_center,blocker_radius))
    engine.retry(); blocked=deepcopy(list(engine.pending['RBC'])); engine.active.clear(); engine.time_s=.01; engine.retry()
    admitted=engine.births[e['particle_id']]
    for k in ['particle_id','species','scheduled_time_s','geometry','q','provenance']: assert admitted[k]==original[k]
    write_json(folder/'07_pending_queue.json',dict(original=original,blocked=blocked,admitted=admitted,queue_sizes=[1,1,0],times_s=[e['scheduled_time_s'],.005,.01],
        candidate_records=engine.admission.candidates,geometry_orientation_id_preserved=True,blocker_center_m=blocker_center,blocker_radius_m=blocker_radius,blocker_fits_side_walls=True))


def lifecycle_restart_validation(folder):
    a=synthetic_engine(); b=synthetic_engine(); comparisons=[]; checkpoint=folder/'checkpoint'
    if checkpoint.exists(): raise FileExistsError('Evidence checkpoint already exists; select fresh run directory')
    try:
        for k in range(1,21):
            t=.25*k/20; a.step_to(t); b.step_to(t)
            if k==10:
                write_checkpoint(checkpoint,b,{'role':'P7_SYNTHETIC_COMMON_PLUG_VALIDATION','P65':'57dfe7bd15b29aca9a941ab4fd74dae081d9f161'})
                b.bridge.close(); assert b.bridge.closed
                s,w,c,v=channel()
                b=read_checkpoint(checkpoint,lambda source:PlugPathAdmission(s,source,width=200e-6,velocity=v),c,
                    {'role':'P7_SYNTHETIC_COMMON_PLUG_VALIDATION','P65':'57dfe7bd15b29aca9a941ab4fd74dae081d9f161'})
            errors={key:plain(a.state()[key])==plain(b.state()[key]) for key in ['events','births','exits','pending','active','scheduler']}
            assert all(errors.values()),errors
            comparisons.append(dict(time_s=t,event_count=len(a.events),**errors))
        write_json(folder/'09_lifecycle.json',dict(role='P4_MIXED_COMMON_PLUG_KINEMATIC_VALIDATION_ONLY',events=a.events,births=a.births,exits=a.exits,
            active=[p.to_dict() for p in a.active.values()],accounting=a.accounting(),commands=a.bridge.commands,
            force_audits=a.bridge.force_audits,neighbor_comparisons=a.neighbor_comparisons,neighbor_mismatches=a.neighbor_mismatches))
        write_rows(folder/'09_lifecycle_balance.csv',a.timeline)
        write_json(folder/'10_restart_parity.json',dict(status='PASS',destroyed_actual_lammps_before_restart=True,comparisons=comparisons,
            final_event_count=len(a.events),checkpoint_event_count=comparisons[9]['event_count'],all_sequence_fields_exact=True,
            active_state_exact=True,next_presampled_rbc_exact=True,full_rng_states_exact=True,uses_seed_replay=False))
        write_json(folder/'10_continuous_final.json',a.state()); write_json(folder/'10_restarted_final.json',b.state())
    finally: a.bridge.close(); b.bridge.close()


class RealAdmission(FiniteSizeAdmission):
    def attempt(self,event,active):
        self.guard=512 if event['species']=='MB' else 1
        return super().attempt(event,active)


def real_smoke(folder):
    provenance,mesh,field,boundaries,audit,samplers,wall=real_inlet()
    Q=samplers['INLET'].Q_m3_s; src=PopulationSource(SONOVUE); scheduler=InjectionScheduler(src,LinearProfile([0],[Q]))
    classifier=ValidationBoundaryClassifier({k:v for k,v in boundaries.items() if k.startswith('OUTLET_')})
    admission=RealAdmission(samplers['INLET'],src,wall=wall,field=field,guard=1)
    query=ValidationNeighborPolicy(20e-6,.5e-6,'P7_REAL_SHORT_ADMISSION_SMOKE_ONLY')
    bridge=DynamicParticleBridge([],query); paths=[]; motion_records=[]
    mu,viscosity=viscosity_from_frozen(REPO/'formal_3D_flow_solver/FEM_SimVascular')
    def provider(i,shape,t):
        sample=field.sample(shape.center_m)
        if not sample.inside_lumen: raise ValueError('Real smoke center outside frozen lumen')
        return np.asarray(sample.velocity_m_s),.5*np.asarray(sample.vorticity_s_inv)
    def mover(active,start,end):
        # This real fixture measures entry congestion with admitted shapes held
        # as exact obstacles. Real mixed-RBC transport is not an accepted model.
        # Physical transport and outlet removal are independently exercised by
        # the long plug fixture and the actual-LAMMPS lifecycle/restart fixture.
        return dict(active)
    engine=PopulationLifecycle(scheduler,admission,classifier,bridge=bridge,mover=mover)
    # First actual MB event plus 1 ms of remaining physical time; fixed dose and Q.
    horizon=scheduler.mb_clock.time_at(1)+.001
    try:
        row=engine.step_to(horizon)
        tets=mesh.cells.reshape(-1,5)[:,1:]; xyz=np.asarray(mesh.points,dtype=float)[tets]
        volume=float(np.sum(np.abs(np.linalg.det(np.stack([xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0],xyz[:,3]-xyz[:,0]],axis=2))))/6)
        result=dict(role='REAL_INLET_INJECTION_SMOKE',horizon_s=horizon,accounting=row,rejection_counts=dict(admission.counts),
            rbc_admission_status='REAL_RBC_ADMISSION_LIMITED_BY_UPSTREAM_MODEL' if row['pending_rbc_count'] else 'NO_BACKLOG_IN_THIS_SHORT_WINDOW',
            real_rbc_passage='NOT_ESTABLISHED',candidate_guard_per_attempt=dict(MB=512,RBC=1),guard_role='IMPLEMENTATION_SAFETY_GUARD',
            retry_schedule='EACH_MOVEMENT_OUTLET_DELETION_AND_NEW_BIRTH; PER_SPECIES_FIFO',
            observed_local_tube_hct=row['active_rbc_volume_m3']/volume,control_volume_m3=volume,control_volume_role='ENTIRE_FROZEN_LUMEN_CENTER_ASSIGNED_DIAGNOSTIC',
            neighbor_mismatches=engine.neighbor_mismatches,motion_role='ADMISSION_ONLY_STATIC_ACTIVE_OBSTACLES_NO_PASSAGE_CLAIM',
            model_for_admitted_MB='P65_SPHERE_ADMISSION_HANDOFF_CHECKS',
            model_for_RBC_admission='UNCHANGED_P3_FREE_OBLATE_THEN_CAPILLARY_SURROGATE',events=engine.events,births=engine.births,
            pending={k:list(v) for k,v in engine.pending.items()},active=[p.to_dict() for p in engine.active.values()],exits=engine.exits,
            candidates=admission.candidates,mb_trajectory=paths,mb_motion_attempts=motion_records)
        write_json(folder/'13_real_smoke.json',result)
        # Save actual official inlet and nearby wall, not a reconstructed cylinder.
        center=samplers['INLET'].triangles.reshape(-1,3).mean(0)
        ids=wall.candidates(center,8e-6)
        write_json(folder/'13_real_geometry.json',dict(inlet_triangles_m=samplers['INLET'].triangles,nearby_wall_triangles_m=wall.triangles[ids]))
        return result
    finally: bridge.close()


def isolated_real_mb_transport(folder):
    """Independent single-MB entry/remaining-time check with unchanged P6.5 dynamics."""
    _,mesh,field,boundaries,audit,samplers,wall=real_inlet()
    source=PopulationSource(SONOVUE); event=source.next_mb(); event.update(particle_id=1,attempt_count=0,
        scheduled_time_s=1/(C_MB*samplers['INLET'].Q_m3_s))
    admission=FiniteSizeAdmission(samplers['INLET'],source,wall=wall,field=field)
    p=admission.attempt(event,{})
    if p is None: raise AssertionError('Isolated first original MB failed the 512-draw admission guard')
    classifier=ValidationBoundaryClassifier({k:v for k,v in boundaries.items() if k.startswith('OUTLET_')})
    query=ValidationNeighborPolicy(20e-6,.5e-6,'P7_ISOLATED_MB_REMAINING_TIME_VALIDATION_ONLY')
    mu,viscosity=viscosity_from_frozen(REPO/'formal_3D_flow_solver/FEM_SimVascular')
    def provider(i,shape,t):
        sample=field.sample(shape.center_m)
        if not sample.inside_lumen: raise ValueError('P0 membership failed')
        return np.asarray(sample.velocity_m_s),.5*np.asarray(sample.vorticity_s_inv)
    with DynamicParticleBridge([p],query) as bridge:
        stepper=Particle65Stepper([p],query,provider,mu,wall=wall,bridge=bridge,boundary_classifier=classifier,physical_time_s=event['scheduled_time_s'])
        states=[dict(time_s=stepper.time_s,particle=p.to_dict())]
        for k in range(1,5):
            stepper.step_to(event['scheduled_time_s']+k*.00025)
            states.append(dict(time_s=stepper.time_s,particle=stepper.read()[0].to_dict()))
        distance=float(np.linalg.norm(stepper.read()[0].position-p.position))
        assert distance>0
        result=dict(status='PASS',role='INDEPENDENT_SINGLE_MB_FIXTURE_NOT_COMBINED_REAL_SMOKE_POPULATION',geometry=event,
            birth_position_m=p.position,artificial_inlet_offset_m=0.,elapsed_after_birth_s=.001,displacement_m=distance,
            model='UNCHANGED_P65_SPHERE_NORMAL_NEARFIELD',viscosity_provenance=viscosity,states=states,candidates=admission.candidates,attempts=stepper.attempts)
        write_json(Path(folder)/'13_isolated_mb_transport.json',result)
        return result


def inlet_audit_validation(folder):
    """Rebuild official-face audit and the complete 100,000-point saved sample."""
    folder=Path(folder); folder.mkdir(parents=True,exist_ok=True)
    provenance,mesh,field,boundaries,audit,samplers,wall=real_inlet()
    write_json(folder/'01_flux_audit.json',audit); write_json(folder/'frozen_provenance.json',provenance)
    sampler=samplers['INLET']; points,ids=sampler.sample(np.random.default_rng(2026092107),100000)
    np.savetxt(folder/'02_actual_inlet_samples.csv',np.column_stack([points,ids]),delimiter=',',header='x_m,y_m,z_m,triangle_id',comments='')
    write_json(folder/'01_inlet_mesh.json',dict(triangles_m=sampler.triangles,q_m_s=sampler.q,weights_m3_s=sampler.weights,expectation_m=sampler.expectation()))
    return audit,samplers
