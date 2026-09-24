#!/usr/bin/env python3
"""Particle-3 staged validation; frozen inputs read-only, no CFD or next stage."""
from pathlib import Path
import argparse,json,sys,time
import numpy as np
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.particle3_cases import *
from particle_3d.audit import sha256,read_frozen
from particle_3d.field import FrozenFEMField
from particle_3d.particle2_audit import check_dependencies
from particle_3d.wall_gap import wall_gap
from particle_3d.wall_contact import contact_velocity
from particle_3d.particle1_cases import AffineValidationField
from particle_3d.particle3_motion import initial_wall_state,contact_trial,SurrogateCaseStopped
from particle_3d.physical_time_refinement import refine_interval,PhysicalTimeRefinementError,MAX_REFINEMENT_DEPTH
from particle_3d.rbc_capillary_surrogate import select_rbc_shape,area_feasible_interval
REPORT=PACKAGE/'reports/particle3';DATA=REPORT/'data'
FEM=REPO/'formal_3D_flow_solver/FEM_SimVascular'


def save(name,records):
    write_json(DATA/(name+'.json'),records)
    if isinstance(records,list) and records:write_rows(DATA/(name+'.csv'),records)


def geometry_stage():
    check_dependencies(REPO)
    p2=json.loads((PACKAGE/'reports/particle2/PARTICLE2_VALIDATION.json').read_text())
    if p2['manual_visual_review']!='PASS' or p2['git_commit']!=P2_COMMIT:raise ValueError('Accepted P2 required')
    for group in ['source_sha256','data_sha256','figure_sha256']:
        for p,h in p2[group].items():
            if sha256(REPO/p)!=h:raise ValueError('P2 dependency changed: '+p)
    wall=WallGeometry.from_frozen(FEM);field=FrozenFEMField.from_frozen(FEM)
    scope=dict(stage='Particle-3',particle2_commit=P2_COMMIT,wall=wall.provenance,
        wall_gap_engine='SUPPORT_FEATURE_STATIONARY_SIGNED_DISTANCE_V0',contact_model='HARD_FRICTIONLESS_KINEMATIC_V0',
        rbc_deformation_model='REDUCED_ORDER_CAPSULE_V0',roundoff_rule='512*float64_eps*actual_coordinate_or_shape_scale; no physical coating',
        maximum_refinement_depth=MAX_REFINEMENT_DEPTH,wall_lubrication='DEFERRED_PARTICLE5',adhesion='OFF',
        transit_time_penalty='DEFERRED_PARTICLE5',rbc_mb_lateral_displacement='DEFERRED_PARTICLE4_5',
        production_particle_timestep_frozen=False,particle4_started=False,no_cfd_executed=True)
    save('00_particle3_scope_and_wall_contract',scope)
    ids=representative_wall_ids(wall);save('01_wall_normals',normal_audit(wall,field,ids))
    sphere,ellipsoid,capsule=geometry_validation_records()
    save('02_sphere_gap',sphere);save('03_ellipsoid_gap',ellipsoid);save('03_capsule_gap',capsule)
    print('Geometry evidence saved: 18 actual WALL normals and synthetic signed-gap cases',flush=True)


def contact_stage():
    wall=WallGeometry(plane_triangle()[None,:,:]);records=[]
    for label,v in [('INCOMING',[2e-4,1e-4,-3e-4]),('TANGENT',[2e-4,1e-4,0]),('SEPARATING',[2e-4,1e-4,3e-4]),('NO_CONTACT',[2e-4,1e-4,-3e-4])]:
        sphere=Sphere([0,0,2e-6 if label=='NO_CONTACT' else 1e-6],1e-6);gap=wall_gap(sphere,wall)
        free=np.array(v);omega=np.array([3.,-4.,5.]);corrected,w=contact_velocity(free,omega,sphere.center_m,gap)
        n=gap.normal_inward;normal=float((corrected+np.cross(w,gap.particle_point_m-sphere.center_m))@n)
        tangent=float(np.linalg.norm((corrected-free)-((corrected-free)@n)*n))
        records.append(dict(case=label,particle_type='MB',center_m=sphere.center_m,radius_m=sphere.radius_m,free_velocity=free,
            corrected_velocity=corrected,omega=omega,normal=n,gap_m=gap.gap_m,contact_point_m=gap.particle_point_m,
            normal_constraint_error_m_s=max(0.,-normal) if gap.state=='TOUCHING' else 0.,tangential_velocity_error_m_s=tangent))
    rot=rotation_matrix(quaternion_from_short_axis([1,0,1]));axes=np.array([2,2,.5])*1e-6;height=np.linalg.norm(axes*(rot.T@np.array([0,0,1])))
    e=Ellipsoid([0,0,height],axes,rot);gap=wall_gap(e,wall);free=np.array([1e-4,2e-4,-1e-4]);omega=np.array([0,20,0.])
    corrected,w=contact_velocity(free,omega,e.center_m,gap);n=gap.normal_inward
    records.append(dict(case='ROTATING_ELLIPSOID',particle_type='RBC',center_m=e.center_m,axes_m=axes,rotation=rot,free_velocity=free,
        corrected_velocity=corrected,omega=omega,normal=n,gap_m=gap.gap_m,contact_point_m=gap.particle_point_m,
        normal_constraint_error_m_s=abs(float((corrected+np.cross(w,gap.particle_point_m-e.center_m))@n)),
        tangential_velocity_error_m_s=float(np.linalg.norm((corrected-free)-((corrected-free)@n)*n))))
    save('04_contact_velocity',records)
    trajectories=[];intervals=[];summaries=[]
    field=AffineValidationField(np.array([1e-6,0,-4e-6]),np.zeros((3,3)))
    for pieces in [1,2,4]:
        state=initial_wall_state([0,0,2e-6],[1,0,0,0],field,wall,radius=.5e-6)
        rows=[state.to_dict()];ledger=[];trial=contact_trial(field,wall,on_accept=lambda s:rows.append(s.to_dict()))
        for i in range(pieces):state,_=refine_interval(state,(i+1)/pieces,trial,ledger=ledger)
        trajectories.extend(dict(requested_dt_s=1/pieces,**r) for r in rows)
        intervals.extend(dict(requested_dt_s=1/pieces,**r) for r in ledger)
        summaries.append(dict(requested_dt_s=1/pieces,end_time_s=state.time_s,accepted_subintervals=len(ledger),
                              max_depth=max(r['depth'] for r in ledger),time_coverage_error=abs(sum(r['dt_s'] for r in ledger)-1),
                              final_center_m=state.shape.center_m,min_gap_m=min(r['wall_gap_m'] for r in rows)))
    save('05_time_trajectories',trajectories);save('05_time_intervals',intervals);save('05_time_summary',summaries)
    print('Hard contact and complete physical-time subdivision evidence saved',flush=True)


def deformation_stage():
    population,representatives=population_and_geometries()
    indices=stratified_indices(population.samples,P3_SELECTION_SEED,n=128)
    chosen=[RBCGeometry.from_population(population,int(i)) for i in indices]
    # Additional legal accepted samples near the shape-rejection edge.
    extra=[RBCGeometry.from_population(population,int(np.argmax(population.samples['r']))),
           RBCGeometry.from_population(population,int(np.argmin(population.samples['D_um'])))]
    records=[];examples=[];q=np.array([1.,0.,0.,0.]);velocity=np.array([0,0,1e-4])
    for tube_index,radius in enumerate(P3_TUBE_RADII_M):
        wall=cylinder_wall(float(radius))
        for group,geometries in [('distribution128',chosen),('representatives_and_extremes',representatives+extra)]:
            for g in geometries:
                d=select_rbc_shape(np.zeros(3),g,q,velocity,wall)
                row=dict(group=group,tube_index=tube_index,tube_radius_m=radius,tube_apothem_m=wall.provenance['tube_apothem_m'],
                    tube_role='VALIDATION_GEOMETRY_PARAMETERS_NOT_MOUSE_CAPILLARY_DISTRIBUTION',rbc_id=g.provenance.rbc_id,
                    D_um=g.provenance.D_um,V_fL=g.provenance.V_fL,a_m=g.a_m,c_m=g.c_m,r=g.r,
                    status=d.status,original_oblate_conflict_gap_m=d.original_oblate_gap_m,area_budget_m2=d.area_budget_m2,
                    radius_interval_m=d.radius_interval_m,search_evaluations=d.search_evaluations,radius_optimality_bound_m=d.radius_optimality_bound_m,
                    volume_m3=g.volume_m3,reason=d.reason)
                if d.shape is not None:row['gap_m']=d.gap.gap_m
                if isinstance(d.shape,Capsule):
                    row.update(R_cap_m=d.shape.radius_m,L_cap_m=d.shape.cylindrical_length_m,
                        capsule_axis_world=d.shape.axis_world,area_ratio=d.shape.area_m2/d.area_budget_m2,
                        volume_relative_error=abs(d.shape.volume_m3/g.volume_m3-1),
                        area_violation_m2=max(0.,d.shape.area_m2-d.area_budget_m2),
                        analytic_max_radius_m=min(area_feasible_interval(g)[1],wall.provenance['tube_apothem_m']))
                (records if group=='distribution128' else examples).append(row)
        print(f'Tube {tube_index+1}/{len(P3_TUBE_RADII_M)} complete',flush=True)
    save('06_deformation_examples',examples);save('07_feasibility_map',records)
    save('07_selection',dict(population_seed=2026092002,selection_seed=P3_SELECTION_SEED,selection_N=128,indices=indices,
         population_array_sha256=population.metadata['sample_structured_array_sha256'],radius_role='VALIDATION_GEOMETRY_PARAMETERS',tube_radii_m=P3_TUBE_RADII_M))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--stage',choices=['geometry','contact','deformation'],required=True)
    args=parser.parse_args();DATA.mkdir(parents=True,exist_ok=True)
    {'geometry':geometry_stage,'contact':contact_stage,'deformation':deformation_stage}[args.stage]()


if __name__=='__main__':main()
