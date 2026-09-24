"""DIAGNOSTIC_ONLY observers. Never imported by the production path.

Private dependency binding observes the existing trial, solve and recursion;
it returns the very same solver objects and never edits their arrays. Logging
and the optional checkpoint model switch both default OFF. SI throughout.
"""
from copy import deepcopy
from functools import lru_cache
from types import MethodType
import numpy as np
from .particle6_stepper import bind_query_dependency
from .planar_wall_hydrodynamics import (stable_wall_frame, local_mobility_4d,
    planar_wall_affine_block, SHEAR_FORCE_COEFFICIENT, SHEAR_TORQUE_COEFFICIENT)
from .wall_gap import wall_gap
from .particle_shapes import Sphere
from .particle82a_geometry import perimeter_edges, segment_distance


def closure_metrics(velocity, gradient, normal, radius, gap):
    """Both compared vectors are m/s. Epsilon is numerical, never a gate."""
    frame = stable_wall_frame(normal); n = frame[:, 2]
    p = np.eye(3)-np.outer(n, n)
    ut = p@np.asarray(velocity); shear = p@np.asarray(gradient)@n
    H = radius+max(gap, 0.); linear = H*shear
    b, l = np.linalg.norm(ut), np.linalg.norm(linear)
    # Relative floating-point floor referenced to the supplied velocity scales.
    reference = max(np.finfo(float).tiny, np.finfo(float).eps*
                    max(np.linalg.norm(velocity), b, l))
    angle = (float(np.degrees(np.arccos(np.clip(ut@linear/(b*l), -1., 1.))))
             if b > reference and l > reference else None)
    return dict(tangential_bulk_velocity_xyz=ut.tolist(), local_shear_vector_xyz=shear.tolist(),
        H_m=float(H), Hs_velocity_scale_xyz=linear.tolist(), bulk_t_speed=float(b),
        linear_shear_speed=float(l), closure_speed_ratio=float(l/max(b, reference)),
        closure_relative_mismatch=float(np.linalg.norm(ut-linear)/max(b, l, reference)),
        closure_angle_deg=angle, small_reference_m_s=float(reference),
        closure_angle_unavailable_reason=None if angle is not None else 'ONE_VECTOR_AT_NUMERICAL_RESOLUTION')


def decomposition(u, omega, gradient, normal, a, h, mu):
    result = closure_metrics(u, gradient, normal, a, h)
    t = stable_wall_frame(normal)[:, :2]
    dr, db, d = planar_wall_affine_block(a, mu, h, normal, gradient, np.r_[u, omega])
    near = local_mobility_4d(a, mu, h/a)[2]
    s = np.array(result['local_shear_vector_xyz'])
    force = SHEAR_FORCE_COEFFICIENT*6*np.pi*mu*a*result['H_m']*s
    torque = SHEAR_TORQUE_COEFFICIENT*4*np.pi*mu*a**3*np.cross(s, stable_wall_frame(normal)[:, 2])
    bulk = np.r_[t.T@u, a*t.T@omega]
    wall = near@np.r_[t.T@force, t.T@torque/a]
    w = d['wall_weight']; target = (1-w)*bulk+w*wall
    result.update(wall_weight=w, gap_ratio=h/a, q_bulk=bulk.tolist(), q_wall=wall.tolist(),
        q_target=target.tolist(), q_bulk_norm=float(np.linalg.norm(bulk)),
        q_wall_norm=float(np.linalg.norm(wall)), q_target_norm=float(np.linalg.norm(target)),
        wall_target_velocity_xyz=(t@wall[:2]).tolist(), wall_target_omega_xyz=(t@wall[2:]/a).tolist(),
        blended_target_velocity_xyz=(t@target[:2]).tolist(), blended_target_omega_xyz=(t@target[2:]/a).tolist(),
        target_tangential_speed=float(np.linalg.norm(target[:2])),
        p9a_candidate_delta_R_norm=float(np.linalg.norm(dr)), p9a_candidate_delta_b_norm=float(np.linalg.norm(db)))
    return result


def inlet_distance(point, origin, inward_normal):
    return float((np.asarray(point)-origin)@inward_normal)


def triangle_change(previous, current):
    if previous is None:
        return dict(wall_triangle_changed=False, normal_change_angle_from_previous_deg=0.)
    a = np.asarray(previous['wall_normal_xyz']); b = np.asarray(current['wall_normal_xyz'])
    # atan2 gives exact zero for repeated identical normals, unlike acos(dot).
    angle = float(np.degrees(np.arctan2(np.linalg.norm(np.cross(a, b)), a@b)))
    return dict(wall_triangle_changed=previous['nearest_wall_triangle_id'] != current['nearest_wall_triangle_id'],
                normal_change_angle_from_previous_deg=angle)


class Observer:
    def __init__(self, env=None, event=None, model='P9A'):
        self.env = env; self.event = deepcopy(event or {}); self.model = model
        self.trials = []; self.states = []; self.depth = 0; self.current = None
        self.nominal_dt = None; self.previous = None; self.checkpoint = None
        if env is None:
            return
        surf = env.boundaries['INLET']; points = np.asarray(surf.points, float)
        faces = surf.faces.reshape(-1, 4)[:, 1:]
        self.origin = points.mean(axis=0)
        self.inward = np.linalg.svd(points-self.origin, full_matrices=False)[2][2]
        gids = np.asarray(surf.point_data['GlobalNodeID'], int)-1
        if self.inward@env.field.velocity_nodes_m_s[gids].mean(axis=0) < 0:
            self.inward *= -1
        self.perimeter, edges = perimeter_edges(points, faces)
        self.rim_ids = set(gids[np.unique(edges)].tolist())
        self.plane = dict(origin_m=self.origin.tolist(), inward_normal=self.inward.tolist(),
                          maximum_nonplanarity_m=float(np.max(abs((points-self.origin)@self.inward))))
        self.geometry = lru_cache(maxsize=512)(self._geometry)

    def _geometry(self, center, radius):
        env = self.env; shape = Sphere(center, radius); g = wall_gap(shape, env.wall)
        x = np.asarray(center); field = env.field.sample(x); n = np.asarray(g.normal_inward)
        lower = max(2e-9, .001*radius)
        result = dict(position_xyz=list(center), radius_m=radius, tetra_id=int(field.tetra_id),
            inside_lumen=bool(field.inside_lumen), wall_gap_m=float(g.gap_m), gap_ratio=float(g.gap_m/radius),
            nearest_wall_triangle_id=int(g.wall_triangle_id), nearest_wall_feature=g.wall_feature,
            wall_normal_xyz=n.tolist(), wall_point_m=np.asarray(g.wall_point_m).tolist(),
            handoff_gap_m=lower, g_nf_m=float(g.gap_m-lower), geometry_roundoff_m=float(g.roundoff_m),
            distance_to_inlet_plane_m=inlet_distance(x, self.origin, self.inward),
            wall_inlet_normal_angle_deg=float(np.degrees(np.arccos(np.clip(n@self.inward, -1., 1.)))),
            wall_point_distance_to_inlet_rim_m=segment_distance(np.asarray(g.wall_point_m), self.perimeter))
        gids = getattr(env.wall, 'global_node_ids', None)
        result['wall_triangle_touches_inlet_rim'] = (bool(self.rim_ids.intersection(gids[g.wall_triangle_id]))
                                                     if gids is not None else None)
        if field.inside_lumen:
            u = field.velocity_m_s; omega = .5*field.vorticity_s_inv
            result.update(FEM_velocity_xyz=u.tolist(), FEM_gradient_3x3=field.velocity_gradient_s_inv.tolist(),
                          FEM_inward_velocity_m_s=float(u@self.inward), bulk_omega_xyz=omega.tolist())
            result.update(decomposition(u, omega, field.velocity_gradient_s_inv, n, radius, g.gap_m, env.mu))
        else:
            result['FEM_unavailable_reason'] = 'CENTER_OUTSIDE_FROZEN_LUMEN; NO_EXTRAPOLATION'
        return result

    def snapshot(self, shape, time):
        result = dict(particle_id=self.event.get('particle_id'), model=self.model, time_s=float(time),
                      physical_time_s=float(time+self.event.get('birth_time_s', 0.)))
        if self.env is not None:
            result.update(self.geometry(tuple(shape.center_m), shape.radius_m))
        return result

    def record_initial(self, stepper):
        if self.env is None:
            return
        p = stepper.read()[0]; row = self.snapshot(p.shape(), stepper.time_s)
        row.update(velocity_xyz=p.velocity.tolist(), omega_xyz=p.omega.tolist(), quaternion=p.q.tolist(),
            state_role='INITIAL_FREE_DIAGNOSTIC_NOT_SOLVED_VELOCITY', trial_accepted=None,
            solver_unavailable_reason='NO_TRIAL_AT_INITIAL_ROW', accepted_dt_s=0., boundary_role='ACTIVE')
        self.states.append(row)


def instrument(stepper, observer, *, enabled=False):
    """Bind a private logging copy; disabled returns the untouched instance."""
    if not enabled:
        return stepper
    method = stepper.advance_cached.__func__
    factory = method.__globals__['v1_trial']
    original_solve = factory.__globals__['solve_resistance']
    original_certificate = factory.__globals__['wall_handoff_certificate']
    original_refine = method.__globals__['refine_interval']

    def solve(system, **kwargs):
        solved = original_solve(system, **kwargs)
        row = observer.current
        row.update(solver_record=deepcopy(solved.record), velocity_xyz=solved.velocity[:3].tolist(),
            omega_xyz=solved.velocity[3:6].tolist(), unconstrained_velocity_xyz=solved.unconstrained[:3].tolist(),
            resistance_condition_estimate=solved.record['condition_estimate'],
            handoff_active=bool(solved.record['contact_count']),
            constraint_multipliers=solved.record.get('multipliers', []),
            delta_R_norm=float(np.linalg.norm(system.planar_matrix.toarray())) if hasattr(system, 'planar_matrix') else 0.,
            delta_b_norm=float(np.linalg.norm(system.planar_rhs)) if hasattr(system, 'planar_rhs') else 0.)
        if observer.env is not None:
            n = np.asarray(row['wall_normal_xyz']); v = solved.velocity[:3]
            row.update(normal_velocity_before_constraint=float(solved.unconstrained[:3]@n),
                normal_velocity_after_constraint=float(v@n), tangential_speed=float(np.linalg.norm(v-(v@n)*n)),
                omega_magnitude=float(np.linalg.norm(solved.velocity[3:6])),
                final_velocity_speed=float(np.linalg.norm(v)),
                velocity_dot_inlet_inward_normal=float(v@observer.inward))
        return solved

    def certificate(*args, **kwargs):
        result = original_certificate(*args, **kwargs)
        observer.current['handoff_certificate'] = deepcopy(result[1])
        observer.current['handoff_certificate_safe'] = bool(result[0])
        return result

    bound_factory = bind_query_dependency(factory, {'solve_resistance': solve,
                                                    'wall_handoff_certificate': certificate})

    def trial_factory(*args, **kwargs):
        trial = bound_factory(*args, **kwargs)
        def observed(old, dt):
            shape = next(iter(old.shapes.values()))
            row = observer.snapshot(shape, old.time_s)
            row.update(trial_index=len(observer.trials), requested_dt_s=observer.nominal_dt,
                trial_dt_s=float(dt), accepted_dt_s=0., subdivision_depth=observer.depth,
                trial_accepted=False, trial_rejection_reason=None, boundary_role=old.boundary_event,
                state_role='TRIAL_START_ACTUAL_VELOCITY_EVALUATION_POINT')
            if observer.env is not None:
                row.update(triangle_change(observer.previous, row))
                observer.previous = row
            observer.current = row
            try:
                result = trial(old, dt)
                new = result[0]
                row.update(trial_accepted=True, accepted_dt_s=float(new.time_s-old.time_s),
                    effective_trial_dt_s=float(new.time_s-old.time_s), boundary_role=new.boundary_event)
                if observer.env is not None:
                    end = observer.snapshot(next(iter(new.shapes.values())), new.time_s)
                    p = stepper.read()[0]
                    end.update(velocity_xyz=p.velocity.tolist(), omega_xyz=p.omega.tolist(), quaternion=p.q.tolist(),
                        trial_index=row['trial_index'], trial_accepted=True, trial_rejection_reason=None,
                        requested_dt_s=observer.nominal_dt, trial_dt_s=float(dt),
                        accepted_dt_s=row['accepted_dt_s'], subdivision_depth=observer.depth,
                        boundary_role=new.boundary_event, velocity_evaluated_time_s=float(old.time_s),
                        state_role='ACCEPTED_ENDPOINT_VELOCITY_HELD_FROM_PRECEDING_TRIAL_START',
                        velocity_dot_inlet_inward_normal=float(p.velocity@observer.inward),
                        handoff_active=row.get('handoff_active'),
                        delta_R_norm=row.get('delta_R_norm'), delta_b_norm=row.get('delta_b_norm'),
                        resistance_condition_estimate=row.get('resistance_condition_estimate'))
                    observer.states.append(end)
                return result
            except Exception as error:
                row['trial_rejection_reason'] = str(error)
                row['exception_type'] = type(error).__name__
                if 'solver_record' not in row:
                    row['solver_unavailable_reason'] = 'TRIAL_TERMINATED_BEFORE_SUCCESSFUL_SOLVE'
                raise
            finally:
                observer.trials.append(row)
                observer.current = None
        return observed

    def recursion(*args, **kwargs):
        old_depth = observer.depth; observer.depth = kwargs.get('depth', 0)
        try:
            return bound_refine(*args, **kwargs)
        finally:
            observer.depth = old_depth
    bound_refine = bind_query_dependency(original_refine, {'refine_interval': recursion})
    stepper.advance_cached = MethodType(bind_query_dependency(method,
        {'v1_trial': trial_factory, 'refine_interval': recursion}), stepper)
    return stepper


def wrap_nominal_steps(stepper, observer, *, enabled=False, checkpoint_p65_time=None):
    """DIAGNOSTIC_ONLY exact-prefix replay; no new checkpoint integrator/guards.

    The optional switch is at an existing nominal boundary, preserving the
    admission, accepted history, provider counters, and full 1.5 s horizon.
    """
    if not enabled:
        return stepper
    original = stepper.step_to
    observer.record_initial(stepper)
    def step(self, end_time):
        observer.nominal_dt = float(end_time-self.time_s)
        if checkpoint_p65_time is not None and observer.checkpoint is None and self.time_s == checkpoint_p65_time:
            from .particle81_cache import cached_step_to
            p = self.read()[0]
            observer.checkpoint = dict(time_s=float(self.time_s), position=p.position.tolist(),
                velocity=p.velocity.tolist(), omega=p.omega.tolist(), q=p.q.tolist(),
                accepted_prefix_samples=len(self.samples), switch='DIAGNOSTIC_ONLY_P9A_TO_P65')
            self.advance_cached, self.cache_info = cached_step_to(self)
            observer.model = 'P65_AFTER_P9A_CHECKPOINT'
            instrument(self, observer, enabled=True)
        return original(end_time)
    stepper.step_to = MethodType(step, stepper)
    return instrument(stepper, observer, enabled=True)
