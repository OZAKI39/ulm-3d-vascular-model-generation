#!/usr/bin/env python3
"""Generate P5 scientific evidence; original FEM is read-only, no CFD execution."""
from pathlib import Path
import argparse,json,sys
import numpy as np
PACKAGE=Path(__file__).resolve().parents[1];REPO=PACKAGE.parent
sys.path.insert(0,str(PACKAGE/'src'))
from particle_3d.audit import sha256,read_frozen
from particle_3d.particle3_cases import write_json,write_rows
from particle_3d.particle_shapes import Sphere,EPS
from particle_3d.hydrodynamic_resistance import viscosity_from_frozen,PhysicalNearField
from particle_3d.resistance_assembly import assemble_resistance_system
from particle_3d.resistance_solver import solve_resistance
from particle_3d.particle5_cases import benchmarks,algebraic_pair_block,sparse_scene,contact_cases,synthetic_motion
from particle_3d.particle5_motion import simulate_resistance
from particle_3d.field import FrozenFEMField
from particle_3d.wall_geometry import WallGeometry
from particle_3d.wall_gap import wall_gap
from particle_3d.validation_boundary import ValidationBoundaryClassifier
FEM=REPO/'formal_3D_flow_solver/FEM_SimVascular';REPORT=PACKAGE/'reports/particle5';DATA=REPORT/'data'


def save(name,value):
    write_json(DATA/(name+'.json'),value)
    write_rows(DATA/(name+'.csv'),value if isinstance(value,list) else [value])


def synthetic(mu,provenance):
    scope=dict(stage='Particle-5',resistance_formulation='AFFINE_BLOCK_RESISTANCE_BALANCE_V0',
        generalized_velocity='[V1,Omega1,...,VN,OmegaN]',matrix_dimension='6N',rhs='R_self @ U_free; stationary WALL; pair relative velocity; F_nonhydro=0',
        dynamic_viscosity_pa_s=mu,viscosity_source=provenance,
        sphere_wall_lubrication_model='LEADING_NORMAL_ASYMPTOTIC_V0',sphere_pair_lubrication_model='LEADING_NORMAL_ASYMPTOTIC_V0',
        nonspherical_lubrication_model='NOT_FROZEN',nonspherical_rejection='NONSPHERICAL_LUBRICATION_NOT_FROZEN',
        production_lubrication_cutoff_frozen=False,production_pair_lubrication_cutoff=None,production_wall_lubrication_cutoff=None,
        validation_nearfield_ratio_max=.01,eligibility_role='VALIDATION_ONLY',synthetic_ratios_role='VALIDATION_GAPS_ONLY',
        production_particle_timestep_frozen=False,real_rbc_hydrodynamics='DEFERRED',particle3_real_rbc_passage='NOT_ESTABLISHED',
        manual_visual_review='PENDING_USER_REVIEW',particle6_started=False,no_cfd_executed=True,
        forbidden_added_physics=['mass','inertia','Brownian','gravity','buoyancy','lift','added mass','Basset history','adhesion','acoustics','friction','restitution','springs','RBC membrane','collision deformation','LAMMPS','hematocrit','continuous injection'])
    save('00_scope',scope)
    stokes,walls,pairs=benchmarks(mu)
    save('01_stokes',stokes);save('02_wall',walls);save('03_pair',pairs)
    r=algebraic_pair_block();permutation=np.r_[np.arange(6,12),np.arange(6)]
    rng=np.random.default_rng(2026092504);x=rng.normal(size=(512,12));diss=np.einsum('bi,ij,bj->b',x,r,x)
    from particle_3d.particle4_cases import shapes_for_geometry,shape_record
    original=shapes_for_geometry()
    save('04_algebraic',dict(role='ALGEBRAIC ASSEMBLY TEST ONLY',warning='NOT PHYSICAL RBC LUBRICATION MODEL',
        coefficient_test=1.,coefficient_units='INTERNALLY_NORMALIZED_TEST_UNITS',physical_solver_allowed=False,
        source_rbc=shape_record(original[1]),source_mb=shape_record(original[0]),matrix=r,
        symmetry_error=float(np.max(np.abs(r-r.T))),swap_error=float(np.max(np.abs(r-r[np.ix_(permutation,permutation)]))),
        minimum_eigenvalue=float(np.linalg.eigvalsh(r).min()),eigenvalue_roundoff_bound=64*EPS*np.linalg.norm(r,2),
        minimum_dissipation=float(diss.min()),sample_count=len(x)))
    shapes,s,candidates=sparse_scene(mu);_,all_s,all_pairs=sparse_scene(mu,candidates=False)
    solved=solve_resistance(s);reference=solve_resistance(all_s)
    scaled=s.matrix.toarray()/np.sqrt(np.outer(s.self_diagonal,s.self_diagonal))
    save('05_sparse',dict(particle_count=len(shapes),ids=s.ids,centers_m=[shapes[i].center_m for i in s.ids],
        radii_m=[shapes[i].radius_m for i in s.ids],matrix=s.matrix.toarray(),self_scaled_matrix=scaled,
        candidate_pairs=candidates,all_pairs_count=len(all_pairs),blocks=s.blocks,
        candidate_query_role='VALIDATION_QUERY_BOUND_NOT_PHYSICAL_MOTION',
        matrix_difference=float(np.max(np.abs((s.matrix-all_s.matrix).toarray()))),
        velocity_difference=float(np.max(np.abs(solved.velocity-reference.velocity))),solver=solved.record,
        sparse_all_pairs_equivalence=bool(np.array_equal(s.matrix.toarray(),all_s.matrix.toarray()))))
    rng=np.random.default_rng(2026092506);rows=[]
    for k in range(1024):
        u=rng.normal(size=len(s.free))*np.tile([1e-5]*3+[10.]*3,len(s.ids))
        d=s.dissipation(u)
        y=np.sqrt(s.self_diagonal)*u
        rows.append(dict(sample=k,**d,quadratic_form=float(u@(s.matrix@u)),
            roundoff_lower_bound=-64*EPS*len(u)*np.linalg.norm(scaled,2)*float(y@y),
            roundoff_role='SELF_SCALED_QUADRATIC_FORM_OPERATION_BOUND'))
    save('06_dissipation',rows)
    save('07_contact',contact_cases(mu))
    summaries=[]
    for kind in ['wall','pair']:
        for divisor in [1,2,4]:
            name=f'11_{kind}_dt{divisor}';summary,states,ledger=synthetic_motion(mu,kind,divisor)
            save(name+'_states',states);save(name+'_intervals',ledger);summaries.append(summary)
    p4_summary,p4_states,p4_ledger=synthetic_motion(mu,'wall',p4=True)
    save('08_p4_wall_states',p4_states);save('08_p4_wall_intervals',p4_ledger)
    save('08_delay',dict(p4=p4_summary,p5=summaries[0],p5_source='11_wall_dt1_states.json',constant_multiplier_used=False))
    # Independent implicit exact ODE: h+r*log(h) decreases at free relative speed.
    for summary in summaries:
        r=1e-6 if summary['kind']=='wall' else .5e-6
        h0=.1*r;speed=1e-6 if summary['kind']=='wall' else 2e-6
        h=summary['minimum_gap_m']
        summary['implicit_ode_defect_m']=float(h-h0+r*np.log(h/h0)+speed*summary['horizon_s'])
        summary['comparison_role']='VALIDATION_ONLY; NOT PRODUCTION TIMESTEP SELECTION'
    save('11_timestep',summaries)
    print('synthetic evidence complete',flush=True)


def real(mu):
    _,mesh,flow,boundaries=read_frozen(FEM)
    field=FrozenFEMField.from_grids(mesh,flow);wall=WallGeometry.from_frozen(FEM)
    source=PACKAGE/'reports/particle4/data/14_real_two_mb_states.json'
    initial_source=PACKAGE/'reports/particle4/data/14_initialization.json'
    states=json.loads(source.read_text());init=json.loads(initial_source.read_text())
    candidates=[]
    for row_index,state in enumerate(states):
        for p in state['particles']:
            if 0<p['wall_gap_m']/p['radius_m']<=.01:
                shape=Sphere(p['center_m'],p['radius_m']);gap=wall_gap(shape,wall);sample=field.sample(shape.center_m)
                vn=float(np.asarray(sample.velocity_m_s)@gap.normal_inward)
                candidates.append((gap.gap_m,row_index,p,shape,gap,sample,vn,state['time_s']))
    minimum=min(candidates,key=lambda r:r[0]);approach=min((c for c in candidates if c[6]<0),key=lambda r:r[0])
    replays=[]
    for label,entry in [('MINIMUM_ORIGINAL_GAP_SEPARATING',minimum),('MINIMUM_ORIGINAL_APPROACHING_GAP',approach)]:
        _,row_index,p,shape,gap,sample,vn,time=entry
        free=np.r_[sample.velocity_m_s,.5*np.asarray(sample.vorticity_s_inv)]
        spec=PhysicalNearField(p['particle_id'],None,gap.gap_m,gap.normal_inward,gap.roundoff_m)
        system=assemble_resistance_system({p['particle_id']:shape},{p['particle_id']:free},[spec],mu=mu)
        solved=solve_resistance(system);n=gap.normal_inward;v=solved.velocity[:3]
        tangent_delta=(np.eye(3)-np.outer(n,n))@(v-free[:3])
        distance=np.linalg.norm(wall.triangle_centers-shape.center_m,axis=1);local=np.argsort(distance)[:48]
        replays.append(dict(case=label,source_path=str(source.relative_to(REPO)),source_sha256=sha256(source),source_row_index=row_index,
            source_time_s=time,original_particle=p,center_m=shape.center_m,radius_m=shape.radius_m,
            original_center_unchanged=True,original_radius_unchanged=True,gap_m=gap.gap_m,roundoff_m=gap.roundoff_m,
            gap_ratio=gap.gap_m/shape.radius_m,normal=n,wall_triangle_id=gap.wall_triangle_id,
            wall_point_m=gap.wall_point_m,wall_triangle_m=wall.triangles[gap.wall_triangle_id],wall_patch_m=wall.triangles[local],
            inside_lumen=bool(sample.inside_lumen),free=free,lubricated=solved.velocity,
            free_normal_velocity_m_s=vn,lubricated_normal_velocity_m_s=float(v@n),attenuation_ratio=float((v@n)/vn),
            analytic_attenuation=gap.gap_m/(gap.gap_m+shape.radius_m),tangential_difference_m_s=float(np.linalg.norm(tangent_delta)),
            angular_difference_s_inv=float(np.linalg.norm(solved.velocity[3:]-free[3:])),eligibility=system.blocks[0],solver=solved.record,
            interpretation='ONE_WAY_FROZEN_FLOW_SPHERE_NORMAL_LEADING_CORRECTION_ONLY'))
    save('09_real_near_wall',replays)
    print('real frozen-state wall replays complete',flush=True)
    shapes={p['particle_id']:Sphere(p['center_m'],p['radius_m']) for p in init['initial_state']['particles']}
    def provider(i,shape,time):
        sample=field.sample(shape.center_m)
        if not sample.inside_lumen:raise ValueError('PARTICLE_CENTER_OUTSIDE_FROZEN_LUMEN')
        return np.asarray(sample.velocity_m_s),.5*np.asarray(sample.vorticity_s_inv)
    save('10_initialization',dict(source=str(initial_source.relative_to(REPO)),source_sha256=sha256(initial_source),
        original_p4_initialization=init,unchanged_dt_s=init['validation_dt_s'],unchanged_horizon_s=init['horizon_s'],
        identical_initial_centers=True,identical_original_sonovue_radii=True,wall_provenance=wall.provenance))
    end,rows,ledger=simulate_resistance(shapes,provider,mu,init['validation_dt_s'],init['horizon_s'],wall=wall,
                                     boundary_classifier=ValidationBoundaryClassifier(boundaries))
    save('10_real_two_mb_states',rows);save('10_real_two_mb_intervals',ledger)
    pair_records=[r for state in rows[1:] for r in state['projection']['pair_eligibility']]
    save('10_pair_eligibility',pair_records)
    reff=np.prod([s.radius_m for s in shapes.values()])/sum(s.radius_m for s in shapes.values())
    active=sum(r['active'] for r in pair_records)
    summary=dict(status='PASS',rows=len(rows),final_time_s=end.time_s,horizon_s=init['horizon_s'],
        requested_dt_s=init['validation_dt_s'],time_coverage_error_s=abs(sum(l['dt_s'] for l in ledger)-end.time_s),
        minimum_pair_gap_m=min(g['gap_m'] for r in rows for g in r['pair_gaps']),relevant_pair_radius_m=reff,
        minimum_pair_gap_ratio=min(g['gap_m']/reff for r in rows for g in r['pair_gaps']),
        minimum_wall_gap_m=min(p['wall_gap_m'] for r in rows for p in r['particles']),pair_active_block_count=active,
        pair_lubrication_status='PAIR_LUBRICATION_NOT_ACTIVE_IN_THIS_REAL_SMOKE' if not active else 'VALIDATION_NEAR_FIELD_ACTIVE',
        wall_active_block_count=sum(b['active'] for r in rows[1:] for b in r['projection']['potential_blocks'] if b['kind']=='WALL'),
        maximum_condition=max(r['projection']['condition_estimate'] for r in rows[1:]),
        maximum_residual=max(r['projection']['relative_backward_residual'] for r in rows[1:]),
        outlet_events=[r['boundary_event'] for r in rows if r['boundary_event'].startswith('OUTLET_')],
        fixed_original_mb_radii=True,initial_state_moved=False,real_rbc_hydrodynamics='DEFERRED',
        deferral_reasons=['REAL_RBC_PASSAGE_NOT_ESTABLISHED','NONSPHERICAL_LUBRICATION_NOT_FROZEN'])
    save('10_real_summary',summary)
    print(json.dumps(summary,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--part',choices=['synthetic','real','all'],default='all');args=parser.parse_args()
    mu,provenance=viscosity_from_frozen(FEM)
    if args.part in ['all','synthetic']:synthetic(mu,provenance)
    if args.part in ['all','real']:real(mu)


if __name__=='__main__':main()
