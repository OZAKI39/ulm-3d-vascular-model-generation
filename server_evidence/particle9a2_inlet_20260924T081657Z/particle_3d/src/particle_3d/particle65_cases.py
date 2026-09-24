"""Fixed, saved-state validation fixtures. All dt/neighbor settings VALIDATION_ONLY."""
from pathlib import Path
from functools import lru_cache
import json
import numpy as np
from .particle_shapes import Sphere
from .particle3_cases import plane_triangle
from .wall_geometry import WallGeometry
from .particle6_validation import ValidationBoundaryClassifier
from .lammps_state import BridgeParticle
from .lammps_neighbors import ValidationNeighborPolicy
from .nearfield_regularization import NearFieldRegularizationV1
from .nearfield_handoff import exact_interactions
from .particle65_motion import Particle65Stepper
from .particle6_stepper import Particle6Stepper
from .audit import read_frozen,sha256
from .field import FrozenFEMField

MU=.00345312
REPO=Path(__file__).resolve().parents[3]


def synthetic(case):
    if case=='wall':
        shapes={17:Sphere([0,0,1.01e-6],1e-6)};wall=WallGeometry([plane_triangle()])
        def provider(i,s,t):return np.array([2e-6,0,-1e-4]),np.array([.1,.2,.3])
    elif case=='pair':
        shapes={17:Sphere([-1.005e-6,0,0],1e-6),203:Sphere([1.005e-6,0,0],1e-6)};wall=None
        def provider(i,s,t):return np.array([5e-5 if i==17 else -5e-5,2e-6,0]),np.array([.1,.2,.3])
    else:raise ValueError(case)
    return shapes,provider,wall,None


@lru_cache(maxsize=1)
def real_fixture():
    fem=REPO/'formal_3D_flow_solver/FEM_SimVascular'
    _,mesh,flow,boundaries=read_frozen(fem);field=FrozenFEMField.from_grids(mesh,flow)
    wall=WallGeometry.from_frozen(fem);classifier=ValidationBoundaryClassifier(boundaries)
    source=REPO/'particle_3d/reports/particle5/data/09_real_near_wall.json'
    row=json.loads(source.read_text())[1]
    assert row['case']=='MINIMUM_ORIGINAL_APPROACHING_GAP'
    p=row['original_particle'];shapes={p['particle_id']:Sphere(p['center_m'],p['radius_m'])}
    def provider(i,s,t):
        sample=field.sample(s.center_m)
        if not sample.inside_lumen:raise ValueError('PARTICLE_CENTER_OUTSIDE_FROZEN_LUMEN')
        return np.asarray(sample.velocity_m_s),.5*np.asarray(sample.vorticity_s_inv)
    provenance=dict(source_path=str(source.relative_to(REPO)),source_sha256=sha256(source),source_index=1,
        selected_saved_state=row['case'],original_particle=p,original_center_unchanged=True,original_radius_unchanged=True,
        original_source_path=row['source_path'],original_source_sha256=row['source_sha256'],original_source_time_s=row['source_time_s'],
        time_origin='REPLAY_ELAPSED_TIME_FROM_EXACT_SAVED_STATE',free_velocity_source='UNCHANGED_FROZEN_FEM_VELOCITY_AND_HALF_VORTICITY')
    return shapes,provider,wall,classifier,provenance


def records_for(shapes,provider):
    return [BridgeParticle.from_shape(i,s,velocity=provider(i,s,0)[0],omega=provider(i,s,0)[1]) for i,s in sorted(shapes.items())]


def interaction_rows(shapes,provider,wall,policy,time,*,projection=None):
    _,rows,_=exact_interactions(shapes,wall,policy,MU);projection=projection or {};ids=sorted(shapes);idx={i:6*k for k,i in enumerate(ids)}
    free=np.concatenate([np.r_[provider(i,shapes[i],time)[0],provider(i,shapes[i],time)[1]] for i in ids])
    hydro=np.asarray(projection.get('hydro_generalized_velocity',free));solved=np.asarray(projection.get('constrained_generalized_velocity',hydro))
    # At accepted endpoints these velocities are held over the preceding interval;
    # their force/constraint evaluation time is recorded separately, not mislabeled.
    evaluated_free=projection.get('free_generalized_velocities')
    if evaluated_free:free=np.concatenate([evaluated_free.get(i,evaluated_free.get(str(i))) for i in ids])
    for r in rows:
        n=np.asarray(r['normal']);i=r['particle_i_id'];j=r['particle_j_id'];row=np.zeros(len(free));row[idx[i]:idx[i]+3]=n
        if j is not None:row[idx[j]:idx[j]+3]=-n
        r.update(time_s=time,free_vn_m_s=float(row@free),lubricated_vn_m_s=float(row@hydro),constrained_vn_m_s=float(row@solved),
            velocity_evaluated_at_time_s=projection.get('evaluated_at_time_s',time),
            constraint_count=len(projection.get('handoff_constraints',[])))
    return rows


def run_case(case,dt,horizon,*,floor=2e-9,old=False,bridge=None):
    if case=='real':shapes,provider,wall,classifier,provenance=real_fixture()
    else:shapes,provider,wall,classifier=synthetic(case);provenance={'fixture':case}
    policy=NearFieldRegularizationV1(floor,'USER_APPROVED_V1' if floor==2e-9 else 'SENSITIVITY_ONLY')
    query=ValidationNeighborPolicy(20e-6,.5e-6,'PARTICLE65_SMALL_SCENE_QUERY_VALIDATION_ONLY')
    records=records_for(shapes,provider);cls=Particle6Stepper if old else Particle65Stepper
    kw={} if old else {'policy':policy}
    stepper=cls(records,query,provider,MU,wall=wall,boundary_classifier=classifier,bridge=bridge,**kw)
    initial_rows=interaction_rows(shapes,provider,wall,policy,0.)
    count=0
    while stepper.time_s<horizon and stepper.boundary_event=='ACTIVE':
        count+=1;stepper.step_to(min(count*dt,horizon))
    rows=initial_rows.copy();position_rows=[dict(time_s=0.,particle_id=p.particle_id,center_m=p.position.tolist(),velocity_m_s=p.velocity.tolist(),radius_m=p.radius) for p in records]
    for world in stepper.accepted_worlds:
        current={p['particle_id']:Sphere(p['center_m'],p['radius_m']) for p in world['particles']}
        projection=world['projection']
        if old:
            projection=dict(projection,constrained_generalized_velocity=np.concatenate([np.r_[p['velocity_m_s'],p['omega_s_inv']] for p in world['particles']]))
            projection['hydro_generalized_velocity']=projection['constrained_generalized_velocity']
        rows.extend(interaction_rows(current,provider,wall,policy,world['time_s'],projection=projection))
        position_rows.extend(dict(time_s=world['time_s'],particle_id=p['particle_id'],center_m=p['center_m'],velocity_m_s=p['velocity_m_s'],radius_m=p['radius_m']) for p in world['particles'])
    contacts=[r for r in rows if r['interaction_state']=='CONTINUUM_HANDOFF_CONTACT'];entry=[];was=False
    for r in rows:
        now=r['interaction_state']=='CONTINUUM_HANDOFF_CONTACT'
        if now and not was:entry.append(r['time_s'])
        was=now
    final=stepper.read();delta=final[0].position-records[0].position;n=np.asarray(initial_rows[0]['normal'])
    covered=sum(r['dt_s'] for r in stepper.ledger);target=stepper.time_s if stepper.boundary_event!='ACTIVE' else horizon
    summary=dict(case=case,model='HISTORICAL_P5_DEFAULT' if old else 'NEAR_FIELD_REGULARIZATION_V1',dt_s=dt,dt_role='VALIDATION_ONLY',
        horizon_s=horizon,h_molecular_floor_m=floor,parameter_role=policy.role,minimum_h_geom_m=min(r['h_geom_m'] for r in rows),
        minimum_g_nf_m=min(r['g_nf_m'] for r in rows),h_lower_m=initial_rows[0]['h_lower_m'],
        handoff_time_s=entry[0] if entry else None,handoff_event_count=len(entry),
        handoff_time_resolution='GEOMETRY_ROUNDOFF_BAND; NO_ENDPOINT_CLIP',
        tangential_displacement_m=float(np.linalg.norm(delta-n*(n@delta))),tangential_reference='INITIAL_EXACT_WALL_OR_PAIR_NORMAL',
        final_position_m=[p.position.tolist() for p in final],final_velocity_m_s=[p.velocity.tolist() for p in final],
        final_time_s=stepper.time_s,covered_time_s=covered,time_coverage_error_s=abs(covered-target),
        outlet_event=stepper.boundary_event,requested_steps=count,accepted_intervals=len(stepper.ledger),
        rejected_trials=sum(not r['accepted'] for r in getattr(stepper,'attempts',[])),
        all_accepted_above_lower=all(r['continuum_state_admissible'] for r in rows),
        finite=bool(all(np.isfinite(p.position).all() and np.isfinite(p.velocity).all() for p in final)))
    return dict(summary=summary,initialization=provenance,interactions=rows,positions=position_rows,
        states=stepper.accepted_worlds,ledger=stepper.ledger,attempts=getattr(stepper,'attempts',[]))
