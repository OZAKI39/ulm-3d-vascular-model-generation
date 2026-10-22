"""Explicit V1 motion policy; unchanged P5 resistance metric and P3 time subdivision.

No position projection, no geometry offsets and no hidden time consumption.
A rejected trial is bisected by the ORIGINAL P3 refine_interval, including the
remaining time after first continuum handoff. Tangential motion remains free.
"""
from dataclasses import replace
import numpy as np
from scipy import sparse
from .nearfield_regularization import NearFieldRegularizationV1,CONTRACT_NAME
from .nearfield_handoff import (exact_interactions,require_admissible_initial,handoff_constraints,
    pair_handoff_certificate,wall_handoff_certificate,first_handoff_event,partitioned_wall_handoff_certificate)
from .resistance_assembly import ResistanceSystem
from .resistance_solver import solve_resistance,contact_jacobian,ResistanceSystemIllConditioned
from .hydrodynamic_resistance import sphere_self_diagonal
from .pair_broadphase import all_pairs
from .pair_geometry import pair_gap
from .wall_gap import wall_gap
from .kinematic_contact import MultiContactInfeasible
from .particle4_motion import ParticleWorldState,initial_world,pair_interval_certificate
from .physical_time_refinement import refine_interval,TrialNeedsSubdivision,swept_clearance_certificate
from .particle6_stepper import Particle6Stepper


def assemble_v1(shapes,free,mu,wall=None,*,policy=None,raw_candidates=None):
    policy=policy or NearFieldRegularizationV1();ids=tuple(sorted(shapes));index={i:6*k for k,i in enumerate(ids)}
    specs,records,_=exact_interactions(shapes,wall,policy,mu)
    required={(r['particle_i_id'],r['particle_j_id']) for r in records if r['kind']=='PAIR' and r['active']}
    if raw_candidates is not None and not required.issubset(raw_candidates):raise ValueError('V1_QUERY_MISSES_PHYSICS_ELIGIBLE_PAIRS')
    diag=np.concatenate([sphere_self_diagonal(shapes[i],mu) for i in ids]);u=np.concatenate([free[i] for i in ids])
    if u.shape!=diag.shape or not np.isfinite(u).all():raise ValueError('Finite 6N free velocity required')
    matrix=sparse.diags(diag,format='csr');blocks=[];rows=[]
    for spec in specs:
        z,n,record=policy.evaluate(spec,shapes,mu)
        # A LAMMPS candidate is only a broadphase candidate, never a physics decision.
        if spec.particle_j_id is not None and raw_candidates is not None and (spec.particle_i_id,spec.particle_j_id) not in raw_candidates:continue
        row=np.zeros(len(diag));row[index[spec.particle_i_id]:index[spec.particle_i_id]+3]=n
        if spec.particle_j_id is not None:row[index[spec.particle_j_id]:index[spec.particle_j_id]+3]=-n
        if z:
            r=sparse.csr_matrix(row.reshape(1,-1));matrix+=z*(r.T@r)
        blocks.append(record);rows.append(row)
    return ResistanceSystem(ids,matrix,diag*u,diag,u,blocks,rows)


def v1_trial(provider,mu,wall,policy,*,boundary_classifier=None,on_accept=None,attempts=None,query=None,matrix_audit=None):
    attempts=[] if attempts is None else attempts
    def trial(old,dt):
        requested_dt=dt;free={i:np.concatenate(provider(i,s,old.time_s)) for i,s in old.shapes.items()}
        raw=None if query is None else query(old.shapes)
        system=assemble_v1(old.shapes,free,mu,wall,policy=policy,raw_candidates=raw)
        contacts,contact_rows=handoff_constraints(old.shapes,wall,policy)
        try:solved=solve_resistance(system,particles=old.shapes,constraints=contacts)
        except (ResistanceSystemIllConditioned,MultiContactInfeasible) as error:raise TrialNeedsSubdivision(str(error)) from error
        v={i:solved.velocity[6*k:6*k+3] for k,i in enumerate(system.ids)};w={i:solved.velocity[6*k+3:6*k+6] for k,i in enumerate(system.ids)}
        end={i:s.moved(s.center_m+dt*v[i]) for i,s in old.shapes.items()}
        events={};event='ACTIVE'
        if boundary_classifier is not None:
            hits=[(hit.segment_fraction,i,hit) for i,s in old.shapes.items() if (hit:=boundary_classifier.first_event(s.center_m,end[i].center_m)) is not None]
            if hits:
                fraction,i,hit=min(hits,key=lambda x:(x[0],x[1]))
                if not hit.role.startswith('OUTLET_'):raise TrialNeedsSubdivision('Center boundary event '+hit.role)
                dt*=fraction;events[i]=hit.role;event=hit.role;end={i:s.moved(s.center_m+dt*v[i]) for i,s in old.shapes.items()}
        handoff_event=first_handoff_event(old.shapes,v,dt,wall,policy)
        if handoff_event is not None:
            event_time=old.time_s+handoff_event['time_from_start_s']
            if event_time<=old.time_s:raise TrialNeedsSubdivision('HANDOFF_EVENT_TIME_NOT_REPRESENTABLE')
            dt=event_time-old.time_s
            end={i:s.moved(s.center_m+dt*v[i]) for i,s in old.shapes.items()}
            events={};event='ACTIVE'
        _,end_records,_=exact_interactions(end,wall,policy,mu)
        attempt=dict(t0_s=old.time_s,requested_dt_s=requested_dt,trial_t1_s=old.time_s+dt,
            minimum_g_nf_m=min((r['g_nf_m'] for r in end_records),default=None),accepted=False)
        proofs=[]
        def reject(reason):
            attempts.append(dict(attempt,reason=reason));raise TrialNeedsSubdivision(reason)
        for i,j in all_pairs(end):
            lower=policy.lower_handoff_gap(policy.reference_length(old.shapes[i],old.shapes[j]))['h_lower_m']
            safe,record=pair_handoff_certificate(old.shapes[i],old.shapes[j],end[i],end[j],lower)
            if not safe:reject(record['proof'])
            geometric,proof=pair_interval_certificate(old.shapes[i],old.shapes[j],end[i],end[j],w[i],w[j],dt)
            if not geometric:reject(proof)
            proofs.append(record)
        wall_gaps={}
        if wall is not None:
            for i,s in end.items():
                lower=policy.lower_handoff_gap(policy.reference_length(s))['h_lower_m']
                safe,record=wall_handoff_certificate(old.shapes[i],s,wall,lower)
                if not safe and record['proof']=='HANDOFF_SWEEP_NOT_CERTIFIED':
                    safe,record=partitioned_wall_handoff_certificate(old.shapes[i],s,wall,lower)
                if not safe:reject(record['proof'])
                geometric,proof,triangle=swept_clearance_certificate(old.shapes[i],s,wall)
                if not geometric:reject(proof)
                proofs.append(record);wall_gaps[i]=wall_gap(s,wall)
        if any(not r['continuum_state_admissible'] for r in end_records):reject('HANDOFF_END_STATE_BELOW_LOWER')
        gaps=[pair_gap(end[i],end[j],i,j) for i,j in all_pairs(end)]
        j,ordered=contact_jacobian(system,old.shapes,contacts)
        audit=dict(solved.record,model=CONTRACT_NAME,potential_blocks=system.blocks,accepted_interactions=end_records,
            handoff_constraints=contact_rows,continuous_certificates=proofs,
            free_generalized_velocities={i:free[i] for i in system.ids},hydro_generalized_velocity=solved.unconstrained,
            constrained_generalized_velocity=solved.velocity,evaluated_at_time_s=old.time_s,position_projection=False,
            handoff_event=handoff_event)
        if matrix_audit is not None:matrix_audit.append(dict(time_s=old.time_s,R=system.matrix.toarray(),b=system.rhs,U=solved.velocity,J=j,
            eligible_pairs=[(r['particle_i_id'],r['particle_j_id']) for r in system.blocks if r['kind']=='PAIR' and r['active']],
            constraint_ids=[c.canonical_id for c in ordered]))
        new=ParticleWorldState(old.time_s+dt,end,v,w,gaps,wall_gaps,audit,events,event)
        attempts.append(dict(attempt,accepted=True,reason='CERTIFIED_G_NF_AND_ORIGINAL_GEOMETRY'))
        if on_accept is not None:on_accept(new)
        return new,'V1_CONTINUOUS_G_NF_AND_ORIGINAL_GEOMETRY; NO_POSITION_PROJECTION'
    return trial


class Particle65Stepper(Particle6Stepper):
    def __init__(self,particles,neighbor_policy,provider,mu,*,policy=None,wall=None,bridge=None,boundary_classifier=None,physical_time_s=0.,step_index=0):
        super().__init__(particles,neighbor_policy,provider,mu,wall=wall,bridge=bridge,boundary_classifier=boundary_classifier,
            physical_time_s=physical_time_s,step_index=step_index)
        self.model=CONTRACT_NAME;self.regularization=policy or NearFieldRegularizationV1();self.attempts=[];self.matrix_audit=None
        require_admissible_initial({p.particle_id:p.shape() for p in self.read()},wall,self.regularization,mu)

    def step_to(self,end_time):
        records=self.read();shapes={p.particle_id:p.shape() for p in records};old=initial_world(shapes,self.provider,self.wall)
        old=replace(old,time_s=self.time_s,velocities={p.particle_id:p.velocity for p in records},omegas={p.particle_id:p.omega for p in records},boundary_event=self.boundary_event,boundary_events=self.boundary_events)
        trial=v1_trial(self.provider,self.mu,self.wall,self.regularization,boundary_classifier=self.boundary_classifier,
            on_accept=self._accepted,attempts=self.attempts,query=None if self.bridge is None else self._raw,matrix_audit=self.matrix_audit)
        world,_=refine_interval(old,end_time,trial,ledger=self.ledger);self.step_index+=1
        return world

    def global_metadata(self):
        return dict(super().global_metadata(),near_field_contract=CONTRACT_NAME,h_molecular_floor_m=self.regularization.h_molecular_floor_m,parameter_role=self.regularization.role)
