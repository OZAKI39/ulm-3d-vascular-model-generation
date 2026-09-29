"""CPU global-time reference using the original P6.5/P9-A.1 trial functions."""
import numpy as np
from .. import particle65_motion as p65
from ..particle6_stepper import bind_query_dependency
from ..particle9a_motion import augment_planar_system
from ..nearfield_regularization import NearFieldRegularizationV1
from ..wall_gap import wall_gap
from ..nearfield_handoff import exact_interactions, handoff_constraints, first_handoff_event
from ..particle4_motion import ParticleWorldState
from ..physical_time_refinement import refine_interval


class ReferenceDynamics:
    def __init__(self, environment, *, pairs=False):
        self.env = environment
        self.pairs = bool(pairs)
        self.policy = NearFieldRegularizationV1()
        self.last_ledger = []

    def provider(self, pid, sphere, time):
        if hasattr(self,'monitor'):self.monitor.provider(pid)
        sample = self.env.field.sample(sphere.center_m)
        if not sample.inside_lumen:
            raise ValueError('CENTER_OUTSIDE_AUTHORITATIVE_FEM')
        return sample.velocity_m_s, .5*sample.vorticity_s_inv

    def _trial_factory(self):
        if self.pairs:
            return self._pair_trial_factory()
        pair_list = {} if self.pairs else {'all_pairs': lambda shapes: []}
        interactions = bind_query_dependency(exact_interactions, pair_list)
        contacts = bind_query_dependency(handoff_constraints, pair_list)
        handoff = bind_query_dependency(first_handoff_event, pair_list)
        from .components import assemble_components,solve_components,sparse_jacobian
        def assembly(shapes, free, mu, wall=None, **options):
            return assemble_components(shapes,free,mu,wall,pairs_enabled=False,
                gradient_at=lambda x:self.env.field.sample(x).velocity_gradient_s_inv,**options)
        return bind_query_dependency(p65.v1_trial, dict(
            assemble_v1=assembly, exact_interactions=interactions,
            handoff_constraints=contacts, first_handoff_event=handoff,
            solve_resistance=solve_components,contact_jacobian=sparse_jacobian,**pair_list))

    def _pair_trial_factory(self):
        from .pairs import candidates, sphere_gap
        from .components import assemble_components, solve_components, sparse_jacobian
        from ..particle4_motion import pair_interval_certificate
        from ..nearfield_handoff import pair_handoff_certificate
        swept=[]
        interactions=bind_query_dependency(exact_interactions,dict(all_pairs=candidates,pair_gap=sphere_gap))
        contacts=bind_query_dependency(handoff_constraints,dict(all_pairs=candidates,pair_gap=sphere_gap))
        def handoff(shapes,velocities,dt,wall,policy):
            swept[:]=candidates(shapes,velocities=velocities,dt=dt)
            function=bind_query_dependency(first_handoff_event,dict(all_pairs=lambda s:list(swept)))
            return function(shapes,velocities,dt,wall,policy)
        def assembly(shapes,free,mu,wall=None,**options):
            return assemble_components(shapes,free,mu,wall,
                gradient_at=lambda x:self.env.field.sample(x).velocity_gradient_s_inv,**options)
        # Exact minimum-distance certificate for spheres is stronger than the
        # generic convex support-plane certificate and needs no GJK witness.
        def geometry(a,b,ea,eb,wa,wb,dt):
            safe,record=pair_handoff_certificate(a,b,ea,eb,0.)
            return safe,record['proof']
        return bind_query_dependency(p65.v1_trial,dict(assemble_v1=assembly,
            exact_interactions=interactions,handoff_constraints=contacts,first_handoff_event=handoff,
            all_pairs=lambda s:list(swept),pair_gap=sphere_gap,pair_interval_certificate=geometry,
            solve_resistance=solve_components,contact_jacobian=sparse_jacobian))

    def advance(self, pool, start, end, on_accept):
        ids = pool.particle_id[:pool.n]
        old = ParticleWorldState(start, pool.shapes(),
            {int(pid): pool.v[k].copy() for k, pid in enumerate(ids)},
            {int(pid): pool.omega[k].copy() for k, pid in enumerate(ids)}, [], {}, {}, {})
        trial = self._trial_factory()(self.provider, self.env.mu, self.env.wall, self.policy,
                                      boundary_classifier=self.env.classifier, on_accept=on_accept)
        self.last_ledger = []
        return refine_interval(old, end, trial, ledger=self.last_ledger)[0]
