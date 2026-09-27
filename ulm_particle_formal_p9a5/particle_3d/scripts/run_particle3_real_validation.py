#!/usr/bin/env python3
"""Actual WALL patch diagnostics followed by unchanged-inlet full-WALL FEM replays."""
from pathlib import Path
import argparse,json,sys,time
from dataclasses import replace
import numpy as np
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.particle3_cases import *
from particle_3d.particle_shapes import unit
from particle_3d.particle1_cases import AffineValidationField
from particle_3d.particle2_cases import initial_orientation
from particle_3d.sonovue_adapter import sample_single_validation_size
from particle_3d.rbc_capillary_surrogate import area_feasible_interval,capsule_length
from particle_3d.wall_gap import wall_gap
from particle_3d.particle3_motion import initial_wall_state,contact_trial,SurrogateCaseStopped
from particle_3d.physical_time_refinement import refine_interval,PhysicalTimeRefinementError
from particle_3d.field import FrozenFEMField
from particle_3d.audit import read_frozen
from particle_3d.validation_boundary import ValidationBoundaryClassifier
REPORT=PACKAGE/'reports/particle3';DATA=REPORT/'data'
FEM=REPO/'formal_3D_flow_solver/FEM_SimVascular'


def save(name,rows):
    write_json(DATA/(name+'.json'),rows)
    if isinstance(rows,list) and rows:write_rows(DATA/(name+'.csv'),rows)


def controlled_stage():
    wall=WallGeometry.from_frozen(FEM);field=FrozenFEMField.from_frozen(FEM)
    _,geometries=population_and_geometries()
    sizes=[sample_single_validation_size('/home/lzy/projects/sonovue_size_distribution_v0',seed=s) for s in [20260920,20260921,20260922]]
    save('08_sonovue_samples',[s.to_dict() for s in sizes])
    ids=representative_wall_ids(wall);static=[];paths=[];intervals=[];summaries=[]
    for index,tid in enumerate([ids[i] for i in [1,4,8,12,16]]):
        tri=wall.triangles[tid];point=tri.mean(axis=0);n=wall.normal_in[tid]
        tangent=unit(tri[1]-tri[0]);inradius=2*wall.areas_m2[tid]/sum(np.linalg.norm(tri[(i+1)%3]-tri[i]) for i in range(3))
        patch=WallGeometry(tri[None,:,:],provenance=dict(role='ACTUAL_WALL_FINITE_TRIANGLE_PATCH_DIAGNOSTIC',original_triangle_id=tid))
        g=geometries[index];q=quaternion_from_short_axis(n);r=sizes[index%3].radius_m
        rmin,rmax,_=area_feasible_interval(g);rc=(rmin+rmax)/2;length=capsule_length(g.volume_m3,rc)
        shapes=[Sphere(point+n*(1.25*r),r),Ellipsoid.from_rbc(point+n*(1.25*g.c_m),g,q),
                Capsule(point+n*1.25*(rc+length/2),n,rc,length)]
        for shape in shapes:
            local=wall_gap(shape,patch);global_gap=wall_gap(shape,wall,inside_lumen=bool(field.sample(shape.center_m).inside_lumen))
            static.append(dict(particle_type='MB' if isinstance(shape,Sphere) else 'RBC',shape_mode=shape.mode,
                triangle_id=tid,triangle_m=tri,center_m=shape.center_m,wall_point_m=local.wall_point_m,particle_point_m=local.particle_point_m,
                normal=local.normal_inward,gap_m=local.gap_m,full_wall_gap_m=global_gap.gap_m,center_inside_lumen=global_gap.inside_lumen,
                geometry=dict(radius_m=r) if isinstance(shape,Sphere) else dict(axes_m=[g.a_m,g.b_m,g.c_m],q=q,R_cap_m=rc,L_cap_m=length,axis=n,V0_m3=g.volume_m3),
                scope='PATCH_ONLY_CONTROLLED_GEOMETRY; full_wall_gap is separately reported, not a global-fit claim',
                g0_ratio=.25,g0_role='VALIDATION_ONLY_NOT_COATING',analytic_patch_gap_m=.25*(r if isinstance(shape,Sphere) else g.c_m if isinstance(shape,Ellipsoid) else rc+length/2)))
        for ptype,size in [('MB',r),('RBC',g.c_m)]:
            for direction in ['APPROACH','OBLIQUE','TANGENT','SEPARATING']:
                g0=.25*size if direction in ['APPROACH','OBLIQUE'] else 0.
                velocity=(-2*g0*n if direction in ['APPROACH','OBLIQUE'] else .25*size*n if direction=='SEPARATING' else np.zeros(3))
                if direction in ['OBLIQUE','TANGENT']:velocity+=.1*inradius*tangent
                synthetic=AffineValidationField(velocity,np.zeros((3,3)))
                state=initial_wall_state(point+n*(size+g0),q,synthetic,patch,radius=r if ptype=='MB' else None,
                    geometry=g if ptype=='RBC' else None,allow_deformation=False)
                rows=[state.to_dict()];ledger=[];trial=contact_trial(synthetic,patch,allow_deformation=False,on_accept=lambda s:rows.append(s.to_dict()))
                for endpoint in [.25,.5,.75,1.]:state,_=refine_interval(state,endpoint,trial,ledger=ledger)
                case=f'patch{tid}_{ptype}_{direction}'
                paths.extend(dict(case=case,particle_type=ptype,direction=direction,original_wall_triangle_id=tid,**row) for row in rows)
                intervals.extend(dict(case=case,**row) for row in ledger)
                corrections=np.array([row['corrected_velocity_m_s'] for row in rows])-velocity
                tangent_error=float(np.max(np.linalg.norm(corrections-(corrections@n)[:,None]*n,axis=1)))
                summaries.append(dict(case=case,particle_type=ptype,direction=direction,triangle_id=tid,scope='ACTUAL_WALL_PATCH_ONLY',
                    min_gap_m=min(row['wall_gap_m'] for row in rows),roundoff_m=state.gap.roundoff_m,
                    max_normal_error_m_s=max(row['normal_constraint_error_m_s'] for row in rows),max_tangential_error_m_s=tangent_error,
                    time_coverage_error_s=abs(sum(row['dt_s'] for row in ledger)-1),accepted_intervals=len(ledger),max_depth=max(row['depth'] for row in ledger),
                    free_velocity_m_s=velocity,normal_inward=n,full_wall_initial_gap_m=wall_gap(shapes[0 if ptype=='MB' else 1],wall).gap_m))
        print(f'Actual WALL patch {tid} controlled cases complete',flush=True)
    save('08_real_wall_static',static);save('09_controlled_paths',paths);save('09_controlled_intervals',intervals);save('09_controlled_summary',summaries)


def replay(particle,geometry_index,dt_index,field,wall,classifier,initial,size,geometries):
    dt=initial['validation_timesteps_s'][-1]/2**dt_index
    name=f'{particle}_g{geometry_index}_dt{dt_index}' if particle=='rbc' else f'mb_dt{dt_index}'
    prefix='11_real_'+name if particle=='rbc' else '10_real_'+name
    geometry=geometries[geometry_index] if particle=='rbc' else None
    rows=[];ledger=[];failure=None;started=time.perf_counter();status='VALIDATION_HORIZON_REACHED'
    state=None
    try:
        state=initial_wall_state(initial['initial_position_m'],initial_orientation(),field,wall,radius=size['radius_m'],geometry=geometry)
        rows.append(state.to_dict())
        trial=contact_trial(field,wall,boundary_classifier=classifier,on_accept=lambda s:rows.append(s.to_dict()))
        for step in range(1,int(np.ceil(initial['horizon_s']/dt))+1):
            state,_=refine_interval(state,min(step*dt,initial['horizon_s']),trial,ledger=ledger)
            if step%50==0:print(f'{name} step={step} t={state.time_s:.8g} mode={state.shape_mode} gap={state.gap.gap_m:.5g}',flush=True)
            if state.boundary_event.startswith('OUTLET_'):status=state.boundary_event;break
    except SurrogateCaseStopped as error:status=error.record['status'];failure=error.record
    except PhysicalTimeRefinementError as error:status='PhysicalTimeRefinementError';failure=error.record
    # If a right child fails, all certified left children remain in rows/ledger.
    last_time=rows[-1]['time_s'] if rows else 0.
    events=[];previous=None
    for row in rows:
        if row['shape_mode']!=previous:
            events.append(dict(time_s=row['time_s'],from_mode=previous,to_mode=row['shape_mode']));previous=row['shape_mode']
    contact_rows=[r for r in rows if r['contact_state']=='TOUCHING']
    contact_duration=sum(b['time_s']-a['time_s'] for a,b in zip(rows,rows[1:]) if a['contact_state']=='TOUCHING')
    contact_episodes=sum(r['contact_state']=='TOUCHING' and (i==0 or rows[i-1]['contact_state']!='TOUCHING') for i,r in enumerate(rows))
    summary=dict(case=name,particle_type=particle.upper(),geometry_index=geometry_index if geometry else None,
        rbc_id=geometry.provenance.rbc_id if geometry else None,r=geometry.r if geometry else None,
        validation_dt_s=dt,dt_index=dt_index,timestep_role='VALIDATION_ONLY',status=status,rows=len(rows),
        initial_center_m=initial['initial_position_m'],initial_q=initial_orientation(),radius_m=size['radius_m'] if geometry is None else None,
        last_accepted_time_s=last_time,minimum_gap_m=min((r['wall_gap_m'] for r in rows),default=None),
        accepted_states_finite=all(np.isfinite(r['center_m']+r['q']+[r['wall_gap_m']]).all() for r in rows),
        accepted_gap_constraint_pass=all(r['wall_gap_m']>=-r['roundoff_m'] for r in rows),
        min_gap_minus_roundoff_m=min((r['wall_gap_m']+r['roundoff_m'] for r in rows),default=None),
        first_contact_time_s=contact_rows[0]['time_s'] if contact_rows else None,contact_count=contact_episodes,
        contact_state_count=len(contact_rows),contact_duration_s=contact_duration,shape_mode_events=events,
        accepted_intervals=len(ledger),max_refinement_depth=max((r['depth'] for r in ledger),default=0),
        time_coverage_error_s=abs(sum(r['dt_s'] for r in ledger)-last_time),
        failed_interval_consumed_time=False,failure=failure,elapsed_wall_time_s=time.perf_counter()-started,
        upstream_orientation_status='NO_UPSTREAM_FREE_OBLATE_SEGMENT' if geometry and not any(r['shape_mode']=='FREE_OBLATE' for r in rows) else 'AVAILABLE' if geometry else 'NOT_APPLICABLE',
        production_particle_timestep_frozen=False)
    save(prefix+'_states',rows);save(prefix+'_intervals',ledger);save(prefix+'_summary',summary)
    print(json.dumps(summary,default=json_safe,allow_nan=False),flush=True)
    return summary


def fem_stage(particle,gi,di):
    _,mesh,flow,boundaries=read_frozen(FEM)
    field=FrozenFEMField.from_grids(mesh,flow);wall=WallGeometry.from_frozen(FEM);classifier=ValidationBoundaryClassifier(boundaries)
    initial=json.loads((PACKAGE/'reports/particle1/data/05_real_initialization.json').read_text())
    size=json.loads((PACKAGE/'reports/particle1/data/01_single_mb_size_provenance.json').read_text())
    _,geometries=population_and_geometries()
    summaries=[]
    for geometry in ([gi] if gi is not None else range(5) if particle=='rbc' else [0]):
        for dt in ([di] if di is not None else range(3)):
            summaries.append(replay(particle,geometry,dt,field,wall,classifier,initial,size,geometries))
    if any(s['status']=='PhysicalTimeRefinementError' for s in summaries):raise SystemExit(1)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--stage',choices=['controlled','mb','rbc'],required=True)
    p.add_argument('--geometry',type=int,choices=range(5));p.add_argument('--dt-index',type=int,choices=range(3))
    args=p.parse_args();DATA.mkdir(parents=True,exist_ok=True)
    if args.stage=='controlled':controlled_stage()
    else:fem_stage(args.stage,args.geometry,args.dt_index)


if __name__=='__main__':main()
