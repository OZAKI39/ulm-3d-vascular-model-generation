"""Executable bridge experiments. Comparisons always evolve independent states."""
from pathlib import Path
from dataclasses import replace
import json,platform,sys,importlib.metadata
import numpy as np
from .particle6_cases import *
from .lammps_bridge import LammpsParticleBridge,PROPERTY_COMMAND
from .lammps_neighbors import exact_physics_pairs
from .particle6_checkpoint import write_checkpoint,read_checkpoint,frozen_provenance,MODEL_FLAGS
from .particle6_stepper import Particle6Stepper
from .audit import read_frozen,sha256
from .field import FrozenFEMField
from .hydrodynamic_resistance import viscosity_from_frozen
from .validation_boundary import ValidationBoundaryClassifier
from .particle5_motion import simulate_resistance


def environment(bridge):
    import lammps
    lmp=bridge._lmp
    return dict(lammps_version=bridge.version,lammps_distribution=importlib.metadata.version('lammps'),
        python_version=sys.version,python_executable=sys.executable,platform=platform.platform(),
        lammps_python_module=lammps.__file__,lammps_library_path=str(lmp.lib._name),
        lammps_mpi_available=lmp.has_mpi_support,mpi_world_size=lmp.extract_setting('world_size'),
        mpich_distribution=importlib.metadata.version('mpich'),mpich_library=str(Path(sys.prefix)/'lib/libmpi.so.12'),
        installed_packages=lmp.installed_packages,capabilities=dict(python_module=True,pair_zero=lmp.has_style('pair','zero'),
        property_atom=lmp.has_style('fix','property/atom'),python_neighbor_access=lmp.find_pair_neighlist('zero')>=0,
        binary_restart='VERIFIED_BY_07_CHECKPOINT_AND_PERMANENT_TESTS'),
        installation_scope='PROJECT_VENV_ONLY',fem_venv_modified=False,system_python_modified=False,
        initial_lammps_commands={'lmp':None,'lmp_serial':None,'lmp_mpi':None},
        initial_python_import='ModuleNotFoundError',parallel_lammps_validated=False,
        dependency_reuse='READ_ONLY .pth to pre-existing FEM numerical packages; no installation into FEM venv')


def neighbor_case(mu):
    shapes,groups,policy=six_pair_shapes();particles=[BridgeParticle.from_shape(i,s) for i,s in shapes.items()]
    with LammpsParticleBridge(particles[::-1],policy) as bridge:
        a=standalone_candidates(shapes,policy);b=bridge.rebuild()
        ea,details=exact_physics_pairs(shapes,a,mu);eb,_=exact_physics_pairs(shapes,b,mu)
        return dict(policy=policy.to_dict(),groups=groups,particles=[p.to_dict() for p in particles],
            all_possible_pairs=all_pairs(shapes),standalone_candidates=a,lammps_candidates=b,
            standalone_exact_pairs=ea,lammps_exact_pairs=eb,exact_gap_records=details,
            false_negatives=len(set(a)-set(b)),mismatch_count=len(set(a)^set(b)),
            extra_candidates=len(set(b)-set(eb)),filtered_mismatch=len(set(ea)^set(eb)))


def scenario(mixed=False):
    if mixed:
        p,policy=mixed_particles();return p,policy,mixed_provider,None,'P4_MIXED_KINEMATIC_VALIDATION_ONLY'
    p,policy,wall=sphere_particles();return p,policy,sphere_provider,wall,'P5_SPHERE_RESISTANCE'


def snapshot(stepper):
    return dict(step=stepper.step_index,time_s=stepper.time_s,particles=[p.to_dict() for p in stepper.read()],
        pairs=standalone_candidates({p.particle_id:p.shape() for p in stepper.read()},stepper.policy) if stepper.bridge is None else stepper.bridge.raw_pairs)


def compare(a,b,steps):
    e=state_errors(a.read(),b.read(),steps=steps)
    pairs=standalone_candidates({p.particle_id:p.shape() for p in a.read()},a.policy)
    e.update(step=steps,time_s=a.time_s,neighbor_mismatch_count=len(set(pairs)^set(b.bridge.raw_pairs)))
    if not e['exact_equal'] or e['neighbor_mismatch_count']:raise AssertionError(e)
    return e


def run_parity(mu,mixed=False,steps=100):
    p,policy,provider,wall,model=scenario(mixed)
    with LammpsParticleBridge(p[::-1],policy) as bridge:
        a=Particle6Stepper(p,policy,provider,mu,wall=wall,model=model)
        b=Particle6Stepper(bridge.read(),policy,provider,mu,wall=wall,model=model,bridge=bridge)
        errors=[];states_a=[snapshot(a)];states_b=[snapshot(b)]
        for step in range(1,steps+1):
            a.step_to(step*.001);b.step_to(step*.001);errors.append(compare(a,b,step))
            states_a.append(snapshot(a));states_b.append(snapshot(b))
        return dict(model=model,steps=steps,dt_s=.001,policy=policy.to_dict(),independently_evolved=True,
            reference='ORIGINAL_UNMODIFIED_P4_OR_P5_TRIAL; ORIGINAL_P3_REFINEMENT; ORIGINAL_P2_Q_INTEGRATOR',
            errors=errors,standalone=states_a,bridge=states_b,force_audits=bridge.force_audits,commands=bridge.commands,
            accepted_substeps_standalone=len(a.accepted_worlds),accepted_substeps_bridge=len(b.accepted_worlds))


def resistance_case(mu):
    p,policy,wall=sphere_particles()
    with LammpsParticleBridge(p[::-1],policy) as bridge:
        stepper=Particle6Stepper(bridge.read(),policy,sphere_provider,mu,wall=wall,bridge=bridge)
        a=resistance_snapshot(p,sphere_provider,mu,wall)
        b=resistance_snapshot(bridge.read(),sphere_provider,mu,wall,stepper=stepper)
        errors={k:float(np.max(np.abs(a[k]-b[k]),initial=0.)) for k in ['R','b','U','J']}
        assert all(e==0 for e in errors.values()) and a['active_constraints']==b['active_constraints']
        return dict(standalone=a,bridge=b,max_errors=errors,sparsity_equal=np.array_equal(a['R']!=0,b['R']!=0),
            active_constraints_equal=a['active_constraints']==b['active_constraints'])


def restart_case(mu,folder,repo,provenance,mixed=True,steps=100,checkpoint_step=40):
    p,policy,provider,wall,model=scenario(mixed)
    continuous=LammpsParticleBridge(p,policy);split=LammpsParticleBridge(p[::-1],policy)
    a=Particle6Stepper(continuous.read(),policy,provider,mu,wall=wall,model=model,bridge=continuous)
    b=Particle6Stepper(split.read(),policy,provider,mu,wall=wall,model=model,bridge=split)
    errors=[];sa=[snapshot(a)];sb=[snapshot(b)];audits=[];commands=[]
    try:
        for step in range(1,steps+1):
            a.step_to(step*.001);b.step_to(step*.001)
            errors.append(compare(a,b,step));sa.append(snapshot(a));sb.append(snapshot(b))
            if step==checkpoint_step:
                before=b.read();before_global=b.global_metadata();old_pairs=b.bridge.raw_pairs
                write_checkpoint(folder,b,provenance,repo)
                audits.extend(split.force_audits);commands.extend(split.commands)
                split.close();assert split.closed
                # No particle constructor/fixture/sampling on the restart path.
                b=read_checkpoint(folder,provenance,provider,mu,wall=wall)
                metadata_error=state_errors(before,b.read(),steps=step)
                assert metadata_error['exact_equal'] and before_global==b.global_metadata() and old_pairs==b.bridge.raw_pairs
                metadata=dict(before=[p.to_dict() for p in before],after=[p.to_dict() for p in b.read()],
                    errors=metadata_error,global_before=before_global,global_after=b.global_metadata(),
                    destroyed_before_new_instance=True,retained_quaternion_role='RETAINED_QUATERNION_NOT_CAPSULE_ORIENTATION')
        audits.extend(a.bridge.force_audits+b.bridge.force_audits);commands.extend(a.bridge.commands+b.bridge.commands)
        assert a.time_s==b.time_s and a.step_index==b.step_index
        return dict(status='PASS',model=model,steps=steps,checkpoint_step=checkpoint_step,
            checkpoint_path=str(Path(folder).relative_to(repo)) if Path(folder).is_relative_to(repo) else str(folder),
            final_errors=errors[-1],errors=errors,continuous=sa,restarted=sb,metadata=metadata,
            physical_time_equal=a.time_s==b.time_s,step_index_equal=a.step_index==b.step_index,
            force_audits=audits,commands=commands)
    finally:
        continuous.close();split.close()
        if b.bridge is not split:b.bridge.close()


def rebuild_stress():
    policy=ValidationNeighborPolicy(5e-6,.5e-6,'CUTOFF_CROSSING_STORAGE_UPDATE_STRESS')
    a=BridgeParticle.from_shape(107,Sphere([0,0,0],1e-6));rows=[]
    with LammpsParticleBridge([a,BridgeParticle.from_shape(701,Sphere([7e-6,0,0],1e-6))],policy) as bridge:
        for step,distance in enumerate([7e-6,6e-6,5.4e-6,4e-6,6e-6,3e-6,8e-6,5.4e-6]):
            prior=bridge.raw_pairs.copy();ps=[a,BridgeParticle.from_shape(701,Sphere([distance,0,0],1e-6))]
            bridge.write(ps);current=bridge.rebuild();expected=standalone_candidates({p.particle_id:p.shape() for p in ps},policy)
            assert current==expected
            rows.append(dict(step=step,distance_m=distance,query_radius_m=policy.query_radius_m,
                previous_list=prior,standalone=expected,lammps=current,standalone_candidate=bool(expected),
                lammps_candidate=bool(current),previous_stale=prior!=expected,mismatch_count=len(set(current)^set(expected))))
        return dict(policy=policy.to_dict(),rows=rows,force_audits=bridge.force_audits,commands=bridge.commands,
            interpretation='EXPLICIT_STORAGE_UPDATES_FOR_QUERY_TEST; NOT_PHYSICAL_INTEGRATION')


def real_case(repo,mu):
    repo=Path(repo);fem=repo/'formal_3D_flow_solver/FEM_SimVascular'
    _,mesh,flow,boundaries=read_frozen(fem);field=FrozenFEMField.from_grids(mesh,flow)
    wall=WallGeometry.from_frozen(fem);classifier=ValidationBoundaryClassifier(boundaries)
    source=repo/'particle_3d/reports/particle5/data/10_initialization.json';p5init=json.loads(source.read_text())
    init=p5init['original_p4_initialization'];shapes={p['particle_id']:Sphere(p['center_m'],p['radius_m']) for p in init['initial_state']['particles']}
    def provider(i,shape,time):
        sample=field.sample(shape.center_m)
        if not sample.inside_lumen:raise ValueError('PARTICLE_CENTER_OUTSIDE_FROZEN_LUMEN')
        return np.asarray(sample.velocity_m_s),.5*np.asarray(sample.vorticity_s_inv)
    dt=init['validation_dt_s'];horizon=init['horizon_s']
    # Run the complete original P5 entrypoint again, independent of the bridge adapter.
    end,original,ledger=simulate_resistance(shapes,provider,mu,dt,horizon,wall=wall,boundary_classifier=classifier)
    policy=ValidationNeighborPolicy(20e-6,.5e-6,'UNCHANGED_P5_REAL_TWO_MB_REPLAY_QUERY_ONLY')
    particles=[BridgeParticle.from_shape(i,s,velocity=provider(i,s,0)[0],omega=provider(i,s,0)[1]) for i,s in shapes.items()]
    with LammpsParticleBridge(particles[::-1],policy) as bridge:
        a=Particle6Stepper(particles,policy,provider,mu,wall=wall,boundary_classifier=classifier)
        b=Particle6Stepper(bridge.read(),policy,provider,mu,wall=wall,boundary_classifier=classifier,bridge=bridge)
        errors=[];sa=[snapshot(a)];sb=[snapshot(b)]
        # Match original P5's requested interval arithmetic exactly.
        t=0.;step=0
        while t<horizon:
            step+=1;target=min(step*dt,horizon);a.step_to(target);b.step_to(target);t=target
            errors.append(compare(a,b,a.step_index));sa.append(snapshot(a));sb.append(snapshot(b))
        bridged=[original[0]]+b.accepted_worlds;adapted=[original[0]]+a.accepted_worlds
        assert len(bridged)==len(original)==len(adapted)
        gap_pair=gap_wall=residual_error=0.;activation_equal=True;outlet_equal=True;original_state_error=0.
        for ref,got,ad in zip(original,bridged,adapted):
            assert ref['time_s']==got['time_s']==ad['time_s']
            for x,y in zip(ref['particles'],got['particles']):
                for k in ['center_m','velocity_m_s','omega_s_inv']:
                    original_state_error=max(original_state_error,float(np.max(np.abs(np.asarray(x[k])-y[k]))))
                gap_wall=max(gap_wall,abs(x['wall_gap_m']-y['wall_gap_m']))
            gap_pair=max(gap_pair,max((abs(x['gap_m']-y['gap_m']) for x,y in zip(ref['pair_gaps'],got['pair_gaps'])),default=0.))
            if ref['projection']:
                activation_equal &= ref['projection']['potential_blocks']==got['projection']['potential_blocks']
                residual_error=max(residual_error,abs(ref['projection']['relative_backward_residual']-got['projection']['relative_backward_residual']))
            outlet_equal &= ref['boundary_event']==got['boundary_event']
        assert original_state_error==gap_wall==gap_pair==residual_error==0 and activation_equal and outlet_equal
        return dict(status='PASS',initialization=p5init,source_path=str(source.relative_to(repo)),source_sha256=sha256(source),
            policy=policy.to_dict(),dt_s=dt,horizon_s=horizon,requested_steps=a.step_index,accepted_substeps=len(b.accepted_worlds),
            errors=errors,standalone=sa,bridge=sb,original_p5_states=original,bridge_world_states=bridged,
            original_p5_state_max_error=original_state_error,pair_gap_max_error_m=gap_pair,wall_gap_max_error_m=gap_wall,
            residual_max_error=residual_error,nearfield_activation_equal=activation_equal,outlet_events_equal=outlet_equal,
            force_audits=bridge.force_audits,commands=bridge.commands,real_rbc_lammps_dynamics=MODEL_FLAGS['real_rbc_lammps_dynamics'])
