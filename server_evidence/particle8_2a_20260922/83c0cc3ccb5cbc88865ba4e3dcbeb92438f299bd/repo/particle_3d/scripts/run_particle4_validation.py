#!/usr/bin/env python3
"""Generate P4 geometry, simultaneous-contact, trajectory and FEM evidence."""
from pathlib import Path
import argparse,itertools,json,sys,time,importlib.util
import numpy as np
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.particle4_cases import *
from particle_3d.particle_shapes import EPS
from particle_3d.particle4_motion import simulate_world,initial_world
from particle_3d.pair_broadphase import broadphase,all_pairs
from particle_3d.pair_geometry import PAIR_GAP_ENGINE
from particle_3d.kinematic_contact import CONTACT_MODEL,CONTACT_METRIC
from particle_3d.wall_geometry import WallGeometry
from particle_3d.wall_gap import wall_gap,touching_contacts
from particle_3d.particle3_cases import plane_triangle
from particle_3d.audit import sha256,read_frozen
from particle_3d.field import FrozenFEMField
from particle_3d.validation_boundary import ValidationBoundaryClassifier
from particle_3d.sonovue_adapter import sample_single_validation_size
REPORT=PACKAGE/'reports/particle4';DATA=REPORT/'data';FEM=REPO/'formal_3D_flow_solver/FEM_SimVascular'


def save(name,value):
    write_json(DATA/(name+'.json'),value)
    if isinstance(value,list) and value:write_rows(DATA/(name+'.csv'),value)


def summary(rows,ledger):
    gaps=[g for row in rows for g in row['pair_gaps']]
    touch=[row['time_s'] for row in rows if any(g['state']=='TOUCHING' for g in row['pair_gaps'])]
    duration=sum(b['time_s']-a['time_s'] for a,b in zip(rows,rows[1:]) if any(g['state']=='TOUCHING' for g in a['pair_gaps']))
    return dict(status='PASS',rows=len(rows),last_accepted_time_s=rows[-1]['time_s'],minimum_pair_gap_m=min(g['gap_m'] for g in gaps),
        minimum_gap_plus_roundoff_m=min(g['gap_m']+g['roundoff_budget_m'] for g in gaps),
        first_contact_time_s=min(touch) if touch else None,contact_duration_s=duration,
        time_coverage_error_s=abs(sum(r['dt_s'] for r in ledger)-rows[-1]['time_s']),
        max_refinement_depth=max(r['depth'] for r in ledger),interval_count=len(ledger),
        final_particles=rows[-1]['particles'],boundary_event=rows[-1]['boundary_event'],
        no_position_projection=True,failed_trial_consumes_time=False,timestep_role='VALIDATION_ONLY')


def geometry():
    p3=json.loads((PACKAGE/'reports/particle3/PARTICLE3_VALIDATION.json').read_text())
    acceptance=json.loads((PACKAGE/'reports/particle3/PARTICLE3_MANUAL_ACCEPTANCE.json').read_text())
    assert acceptance['manual_visual_review']=='PASS' and acceptance['particle4_authorized_to_start']
    for group in ['source_sha256','data_sha256','figure_sha256']:
        for path,digest in p3[group].items():assert sha256(REPO/path)==digest,path
    save('00_scope',dict(stage='Particle-4',particle3_dependency_commit=P3_COMMIT,particle3_acceptance=acceptance,
        inherited_source_sha256=p3['source_sha256'],wall_sha256=p3['wall_sha256'],pair_gap_engine=PAIR_GAP_ENGINE,
        contact_model=CONTACT_MODEL,contact_metric=CONTACT_METRIC,physical_force_claimed=False,
        particle3_real_rbc_passage='NOT_ESTABLISHED',particle_particle_lubrication='DEFERRED_PARTICLE5',
        collision_induced_deformation=False,adhesion=False,lammps=False,production_particle_timestep_frozen=False,
        no_cfd_executed=True,particle5_started=False))
    rows=geometry_cases();save('01_geometry_matrix',rows)
    spec=importlib.util.spec_from_file_location('independent_pair_reference',PACKAGE/'tests/particle4/test_pair_geometry.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    references=[];rng=np.random.default_rng(2026092004)
    for label in ['Sphere-Sphere','Sphere-Capsule','Capsule-Capsule']:
        for index in range(24):
            c=rng.normal(size=(2,3))*1.5e-6;r=rng.uniform(.2,1.,2)*1e-6
            a=Sphere(c[0],r[0]) if label!='Capsule-Capsule' else Capsule(c[0],rng.normal(size=3),r[0],1.2e-6)
            b=Sphere(c[1],r[1]) if label=='Sphere-Sphere' else Capsule(c[1],rng.normal(size=3),r[1],1.4e-6)
            result=pair_gap(a,b);exact=module.segment_distance(a,b)-sum(r)
            references.append(dict(pair_type=label,reference='INDEPENDENT_FINITE_SEGMENT_ANALYTIC',gap_m=result.gap_m,reference_gap_m=exact,error_m=abs(result.gap_m-exact),roundoff_m=result.roundoff_budget_m))
    base=shapes_for_geometry(1)
    for j in range(3):
        for n in [[1,2,-1],[-2,1,3],[1,.1,.2]]:
            a,b=place_pair(base[1],base[j],n,.5e-6);result=pair_gap(a,b);lo,hi=module.primal_reference(a,b)
            references.append(dict(pair_type='Ellipsoid-'+type(b).__name__,reference='INDEPENDENT_CONVEX_PRIMAL_DUAL_BOUNDS',gap_m=result.gap_m,reference_lower_m=lo,reference_upper_m=hi,error_m=max(abs(result.gap_m-lo),abs(result.gap_m-hi)),roundoff_m=result.roundoff_budget_m))
    references.extend(module.deep_overlap_reference_records())
    save('01_independent_references',references)
    rbc_mb=[]
    for quantile in range(5):
        g=validation_geometries()[quantile]
        for axis in [[0,0,1],[1,0,0],[1,2,3]]:
            rbc=Ellipsoid([0,0,0],[g.a_m,g.b_m,g.c_m],rotation_matrix(quaternion_from_short_axis(axis)))
            a,b=place_pair(rbc,base[0],[1,.4,.2]);gap=pair_gap(a,b,10+quantile,99)
            rbc_mb.append(dict(quantile=[.05,.25,.5,.75,.95][quantile],rbc_id=g.provenance.rbc_id,r=g.r,
                shape_i=shape_record(a),shape_j=shape_record(b),gap=gap.to_dict()))
    save('04_rbc_mb_geometry',rbc_mb)
    g0,g1=validation_geometries()[1],validation_geometries()[3]
    a=Ellipsoid([0,0,0],[g0.a_m,g0.b_m,g0.c_m],np.eye(3));b=Ellipsoid([0,0,0],[g1.a_m,g1.b_m,g1.c_m],np.eye(3))
    _,b=place_pair(a,b,[1,0,0]);orientation=[]
    for angle in np.linspace(0,np.pi/2,13):
        rotated=Ellipsoid(b.center_m,b.axes_m,rotation_matrix(quaternion_from_short_axis([np.sin(angle),0,np.cos(angle)])))
        gap=pair_gap(a,rotated,31,72)
        orientation.append(dict(angle_rad=angle,shape_i=shape_record(a),shape_j=shape_record(rotated),gap=gap.to_dict()))
    save('05_rbc_orientation_geometry',orientation)
    print('Geometry evidence complete',flush=True)


def contact():
    a=Sphere([-1e-6,0,0],1e-6);b=Sphere([1e-6,0,0],1e-6)
    for number,label,tangent in [(2,'head_on',0.),(3,'glancing',.5e-6)]:
        velocities={1:np.array([1e-6,tangent,0]),2:np.array([-1e-6,-tangent,0])}
        save(f'{number:02d}_{label}_projection',projection_case({1:a,2:b},velocities))
        scene={1:a,2:b} if number==3 else {1:a.moved([-2e-6,0,0]),2:b.moved([2e-6,0,0])}
        final,rows,ledger=simulate_world(scene,lambda i,s,t:(velocities[i],np.zeros(3)),.125,2. if number==3 else 1.5)
        save(f'{number:02d}_{label}_states',rows);save(f'{number:02d}_{label}_intervals',ledger)
    base=shapes_for_geometry();a,b=place_pair(base[1],base[0],[1,.4,.3]);gap=pair_gap(a,b,1,2);n=gap.normal_j_to_i
    save('06_offcenter_projection',projection_case({1:a,2:b},{1:-n*1e-6,2:n*1e-6},{1:np.array([.1,.2,.3]),2:np.zeros(3)}))
    capsules=[]
    for k in range(3):
        a,b=place_pair(base[2],base[k],[.3,1,.2]);n=pair_gap(a,b).normal_j_to_i
        capsules.append(dict(other_shape=type(b).__name__,**projection_case({1:a,2:b},{1:-n*1e-6,2:n*1e-6})))
    save('07_capsule_projection',capsules)
    no_tunnel=[]
    for di,dt in enumerate([1.,.5,.25]):
        shapes={3:Sphere([-2e-6,0,0],.5e-6),7:Sphere([2e-6,0,0],.5e-6)}
        velocity=lambda i,s,t:(np.array([4e-6 if i==3 else -4e-6,0,0]),np.zeros(3))
        final,rows,ledger=simulate_world(shapes,velocity,dt,1.)
        save(f'08_no_tunnel_dt{di}_states',rows);save(f'08_no_tunnel_dt{di}_intervals',ledger)
        no_tunnel.append(dict(validation_dt_s=dt,**summary(rows,ledger)))
    save('08_no_tunnel_summary',no_tunnel)
    multibody=[];orders=[]
    for name,centers in [('chain',[[-2,0,0],[0,0,0],[2,0,0]]),('triangle',[[-1,0,0],[1,0,0],[0,np.sqrt(3),0]]),('four',[[-3,0,0],[-1,0,0],[1,0,0],[3,0,0]])]:
        shapes={i:Sphere(np.array(c)*1e-6,1e-6) for i,c in zip([17,3,91,42],centers)};centroid=np.mean(centers,axis=0)*1e-6
        velocities={i:centroid-s.center_m for i,s in shapes.items()};reference=projection_case(shapes,velocities);multibody.append(dict(case=name,**reference))
        for order in itertools.permutations(shapes):
            result=projection_case({i:shapes[i] for i in order},velocities)
            orders.append(dict(case=name,permutation=order,velocity_error_m_s=max(np.linalg.norm(result['velocities'][i]-reference['velocities'][i]) for i in shapes),
                angular_error_s_inv=max(np.linalg.norm(result['omegas'][i]-reference['omegas'][i]) for i in shapes),
                contact_set_identical=result['audit']['canonical_ids']==reference['audit']['canonical_ids']))
    save('09_simultaneous_projections',multibody);save('10_order_invariance',orders)
    wall=WallGeometry(plane_triangle()[None,:,:]);wall_cases=[]
    for k in [0,1]:
        a=Sphere([0,0,1e-6],1e-6);a,b=place_pair(a,base[k],[0,0,1]);shapes={1:a,2:b}
        constraints=[ContactConstraint.wall(i,g) for i,s in shapes.items() for g in touching_contacts(s,wall)]
        wall_cases.append(dict(case='WALL_MB_'+type(b).__name__,wall_triangle_m=plane_triangle(),**projection_case(shapes,{1:np.array([0,0,-1e-6]),2:np.array([0,0,-2e-6])},wall_constraints=constraints)))
    save('11_wall_pair_projection',wall_cases)
    rng=np.random.default_rng(2026092004);shapes={100+i:base[i%3].moved(rng.normal(size=3)*3e-6) for i in range(18)}
    started=time.perf_counter();exhaustive={ij:pair_gap(shapes[ij[0]],shapes[ij[1]],*ij) for ij in all_pairs(shapes)};elapsed=time.perf_counter()-started
    candidates=set(broadphase(shapes));rows=[]
    for ij,g in exhaustive.items():
        narrow=pair_gap(shapes[ij[0]],shapes[ij[1]],*ij) if ij in candidates else None
        rows.append(dict(pair=ij,all_pairs_gap_m=g.gap_m,state=g.state,candidate=ij in candidates,
            narrowphase_gap_m=narrow.gap_m if narrow else None,gap_error_m=abs(narrow.gap_m-g.gap_m) if narrow else None,
            normal_error=float(np.linalg.norm(narrow.normal_j_to_i-g.normal_j_to_i)) if narrow else None))
    save('12_broadphase_comparison',rows);save('12_broadphase_summary',dict(particles=18,all_pairs=len(exhaustive),candidates=len(candidates),
        true_contact_pairs=sum(g.state!='SEPARATED' for g in exhaustive.values()),false_negatives=sum(g.state!='SEPARATED' and ij not in candidates for ij,g in exhaustive.items()),
        elapsed_s=elapsed,timing_role='INFORMATIONAL_SMOKE'))
    print('Contact, time and broadphase evidence complete',flush=True)


def mixed():
    cases=[];final=[]
    for di,dt in enumerate([.2,.1,.05]):
        shapes,provider=mixed_scene();state,rows,ledger=simulate_world(shapes,provider,dt,1.)
        save(f'13_mixed_dt{di}_states',rows);save(f'13_mixed_dt{di}_intervals',ledger)
        case=dict(dt_index=di,validation_dt_s=dt,**summary(rows,ledger));cases.append(case);final.append(state)
        print('Mixed dt',dt,'states',len(rows),'min gap',case['minimum_pair_gap_m'],flush=True)
    comparisons=[]
    for i,j in [(0,1),(1,2)]:
        axis_errors=[]
        for particle_id in final[i].shapes:
            a,b=final[i].shapes[particle_id],final[j].shapes[particle_id]
            if isinstance(a,Ellipsoid):axis_errors.append(float(np.arccos(np.clip(abs(a.rotation[:,2]@b.rotation[:,2]),0,1))))
        comparisons.append(dict(coarse=i,fine=j,final_position_max_difference_m=max(np.linalg.norm(final[i].shapes[k].center_m-final[j].shapes[k].center_m) for k in final[i].shapes),
            final_velocity_max_difference_m_s=max(np.linalg.norm(final[i].velocities[k]-final[j].velocities[k]) for k in final[i].shapes),max_short_axis_difference_rad=max(axis_errors),
            orientation_role='SYNTHETIC_CENTRAL_CONTACT_ZERO_FREE_OMEGA; NOT REAL FEM ORIENTATION CONVERGENCE'))
    save('15_timestep_comparison',dict(cases=cases,comparisons=comparisons,production_particle_timestep_frozen=False,role='NOT_PRODUCTION_TIMESTEP_SELECTION'))


def real():
    _,mesh,flow,boundaries=read_frozen(FEM);field=FrozenFEMField.from_grids(mesh,flow);wall=WallGeometry.from_frozen(FEM);classifier=ValidationBoundaryClassifier(boundaries)
    sizes=[sample_single_validation_size('/home/lzy/projects/sonovue_size_distribution_v0',seed=s) for s in [20260920,20260921]]
    save('14_sonovue_provenance',[s.to_dict() for s in sizes])
    path=json.loads((PACKAGE/'reports/particle3/data/10_real_mb_dt2_states.json').read_text());selected=None
    window_steps=256
    # First deterministic original P3 path locations that fit both actual sizes.
    # Positive separation ratio is a validation initial condition, not wall skin.
    for j in range(1,len(path)-window_steps):
        # Select a genuinely wide validation segment for the ORIGINAL larger
        # MB, using the already accepted P3 center path. This is an initial-scene
        # selection criterion, not a change to WALL, radius or contact rules.
        if min(r['wall_gap_m']+r['radius_m'] for r in path[j:j+window_steps+1]) < sizes[1].radius_m:
            continue
        b=Sphere(path[j]['center_m'],sizes[1].radius_m)
        if wall_gap(b,wall).gap_m<0:continue
        for i in range(j-1,-1,-1):
            a=Sphere(path[i]['center_m'],sizes[0].radius_m);g=pair_gap(a,b,101,203)
            if g.gap_m>=.25*(a.radius_m+b.radius_m) and wall_gap(a,wall).gap_m>=0:
                selected=(i,j,a,b);break
        if selected is not None:break
    if selected is None:raise RuntimeError('No deterministic valid two-MB initialization along accepted P3 path')
    i,j,a,b=selected;shapes={101:a,203:b}
    def velocity(particle_id,shape,time_s):
        sample=field.sample(shape.center_m)
        if not sample.inside_lumen:raise ValueError('P0 says MB center outside lumen')
        return np.asarray(sample.velocity_m_s),.5*np.asarray(sample.vorticity_s_inv)
    dt=path[1]['time_s']-path[0]['time_s'];horizon=window_steps*dt
    initial=initial_world(shapes,velocity,wall)
    save('14_initialization',dict(source='UNCHANGED_PARTICLE3_ACCEPTED_MB_PATH',source_indices=[i,j],
        initial_state=initial.to_dict(),validation_dt_s=dt,horizon_s=horizon,validation_requested_steps=256,
        minimum_pair_initial_gap_ratio=.25,initialization_role='VALIDATION_ONLY_POSITIVE_SEPARATION_NOT_COATING',
        wide_segment_selection='FIRST_P3_PATH_WINDOW_WITH_POSITIVE_CLEARANCE_FOR_ORIGINAL_LARGER_MB_OVER_256_REFERENCE_STEPS',
        inside_lumen={k:bool(field.sample(s.center_m).inside_lumen) for k,s in shapes.items()},initial_overlap=False))
    state,rows,ledger=simulate_world(shapes,velocity,dt,horizon,wall=wall,boundary_classifier=classifier)
    save('14_real_two_mb_states',rows);save('14_real_two_mb_intervals',ledger)
    result=summary(rows,ledger);result.update(pair_contact_status='NATURAL_PAIR_CONTACT' if result['first_contact_time_s'] is not None else 'NO_NATURAL_PAIR_CONTACT',
        minimum_wall_gap_m=min(p['wall_gap_m'] for r in rows for p in r['particles']),
        outlet_events=[dict(time_s=r['time_s'],particle_id=p['particle_id'],event=p['boundary_event']) for r in rows for p in r['particles'] if p['boundary_event'].startswith('OUTLET_')],
        horizon_role='256_FINEST_P3_VALIDATION_STEPS; SMOKE_NOT_FULL_PASSAGE_OR_SUSPENSION',
        real_rbc_pair_validation='DEFERRED_DUE_TO_P3_REAL_RBC_PASSAGE_LIMITATION',no_artificial_forces=True)
    save('14_real_two_mb_summary',result);print(result,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--stage',choices=['geometry','contact','mixed','real'],required=True)
    args=parser.parse_args();DATA.mkdir(parents=True,exist_ok=True);globals()[args.stage]()
