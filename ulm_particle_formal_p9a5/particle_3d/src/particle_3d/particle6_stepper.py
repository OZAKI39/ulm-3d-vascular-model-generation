"""Adapter around the unchanged P4/P5 trial bytecode and P3 time integrator.

Private function globals supply the neighbor-query dependency. No historical
module, function code object, formula, or global binding is mutated. This avoids
process-wide monkeypatching while reusing the exact previously tested integrator.
"""
from dataclasses import replace
from types import FunctionType
import numpy as np
from . import particle4_motion as p4,particle5_motion as p5
from .physical_time_refinement import refine_interval
from .pair_broadphase import broadphase as original_broadphase
from .wall_gap import wall_gap
from .hydrodynamic_resistance import PhysicalNearField
from .resistance_assembly import assemble_resistance_system
from .lammps_state import BridgeParticle
from .lammps_neighbors import standalone_candidates
from .rbc_orientation import advance_orientation


def bind_query_dependency(function,overrides):
    namespace=dict(function.__globals__);namespace.update(overrides)
    result=FunctionType(function.__code__,namespace,function.__name__,function.__defaults__,function.__closure__)
    result.__kwdefaults__=dict(function.__kwdefaults__ or {})
    return result


class Particle6Stepper:
    def __init__(self,particles,policy,velocity_provider,mu,*,bridge=None,model='P5_SPHERE_RESISTANCE',
                 wall=None,boundary_classifier=None,physical_time_s=0.,step_index=0,boundary_event='ACTIVE',boundary_events=None):
        if model not in ['P5_SPHERE_RESISTANCE','P4_MIXED_KINEMATIC_VALIDATION_ONLY']:raise ValueError('Unapproved upstream model')
        self.particles=list(particles);self.policy=policy;self.provider=velocity_provider;self.mu=mu
        self.bridge=bridge;self.model=model;self.wall=wall;self.boundary_classifier=boundary_classifier
        self.time_s=float(physical_time_s);self.step_index=int(step_index);self.boundary_event=boundary_event
        self.boundary_events={} if boundary_events is None else dict(boundary_events)
        self.ledger=[];self.accepted_worlds=[];self.neighbor_queries=[];self.resistance_systems=[]
        self.last_projection={}
        if bridge is not None:self.particles=bridge.read() # LAMMPS is authoritative storage

    def read(self):
        return self.bridge.read() if self.bridge is not None else list(self.particles)

    def _raw(self,shapes):
        if self.bridge is None:pairs=standalone_candidates(shapes,self.policy)
        else:
            current={p.particle_id:p for p in self.bridge.read()}
            if set(current)!=set(shapes) or any(not np.array_equal(current[i].position,s.center_m) for i,s in shapes.items()):
                raise ValueError('BRIDGE_STORAGE_NOT_SYNCHRONIZED_WITH_ACCEPTED_STATE')
            pairs=self.bridge.rebuild()
        self.neighbor_queries.append(dict(time_s=self.time_s,pairs=pairs))
        return pairs

    @staticmethod
    def _require_coverage(required,raw):
        if not set(required).issubset(raw):raise ValueError('VALIDATION_NEIGHBOR_QUERY_DOES_NOT_COVER_UPSTREAM_BROADPHASE')
        return sorted(set(raw).intersection(required))

    def _broadphase(self,shapes,**kwargs):
        return self._require_coverage(original_broadphase(shapes,**kwargs),self._raw(shapes))

    def _assemble(self,shapes,free,mu,wall=None,*,use_broadphase=True,role='VALIDATION_ONLY'):
        if role!='VALIDATION_ONLY':raise ValueError('P6 dynamics preserves P5 real-validation eligibility')
        pairs=self._require_coverage(p5.validation_candidates(shapes),self._raw(shapes))
        walls=[]
        if wall is not None:
            for i,s in shapes.items():
                g=wall_gap(s,wall);walls.append(PhysicalNearField(i,None,g.gap_m,g.normal_inward,g.roundoff_m,role))
        system=assemble_resistance_system(shapes,free,walls,[p5.pair_spec(shapes,i,j,role=role) for i,j in pairs],mu=mu)
        self.resistance_systems.append(system)
        return system

    def _accepted(self,world):
        previous={p.particle_id:p for p in self.read()};dt=world.time_s-self.time_s
        if dt<=0:raise ValueError('Accepted state must advance physical time')
        records=[]
        for i,s in sorted(world.shapes.items()):
            p=previous[i]
            q=p.q if p.mode_code==3 else advance_orientation(p.q,world.omegas[i],dt)
            records.append(BridgeParticle.from_shape(i,s,q=q,velocity=world.velocities[i],omega=world.omegas[i]))
        if self.bridge is not None:
            self.bridge.write(records);self.bridge.rebuild()
        self.particles=records;self.time_s=world.time_s;self.boundary_event=world.boundary_event
        self.boundary_events=world.boundary_events;self.last_projection=world.projection
        self.accepted_worlds.append(world.to_dict())

    def step_to(self,end_time):
        records=self.read();shapes={p.particle_id:p.shape() for p in records}
        old=p4.initial_world(shapes,self.provider,self.wall)
        old=replace(old,time_s=self.time_s,velocities={p.particle_id:p.velocity for p in records},
            omegas={p.particle_id:p.omega for p in records},boundary_event=self.boundary_event,boundary_events=self.boundary_events)
        if self.model=='P5_SPHERE_RESISTANCE':
            if self.bridge is None:
                # Reference calls original P5 directly, without query substitution.
                factory=p5.resistance_trial
            else:factory=bind_query_dependency(p5.resistance_trial,dict(assemble_scene=self._assemble))
            trial=factory(self.provider,self.mu,wall=self.wall,boundary_classifier=self.boundary_classifier,on_accept=self._accepted)
        else:
            factory=p4.world_trial if self.bridge is None else bind_query_dependency(p4.world_trial,dict(broadphase=self._broadphase))
            trial=factory(self.provider,wall=self.wall,boundary_classifier=self.boundary_classifier,on_accept=self._accepted)
        world,_=refine_interval(old,end_time,trial,ledger=self.ledger)
        self.step_index+=1
        return world

    def global_metadata(self):
        return dict(physical_time_s=self.time_s,particle_step_index=self.step_index,model=self.model,
            boundary_event=self.boundary_event,boundary_events=self.boundary_events,
            rng_states=None,rng_role='DETERMINISTIC_CONTINUATION; NO_RESAMPLING_ON_RESTART')
