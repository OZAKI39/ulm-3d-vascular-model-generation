"""Paired A/B/C diagnostics from one immutable event, with retained retries."""
from copy import deepcopy
from types import SimpleNamespace
import hashlib
import numpy as np
from .particle82a_geometry import (InletGeometryAudit, classify_failure,
    first_admissible, lower_gap)
from .particle82_point_native import NativePointTracer
from .sonovue_adapter import read_sonovue
from .particle7_cases import SONOVUE

MASTER_SEED = 2026092282
GUARD = 512
_AUDIT = None


def context():
    global _AUDIT
    if _AUDIT is None:
        from .particle81_simulation import environment
        env = environment()
        contract, distribution, receipt = read_sonovue(SONOVUE)
        _AUDIT = SimpleNamespace(env=env, geometry=InletGeometryAudit(env),
            native=NativePointTracer(env), distribution=distribution, contract=contract)
    return _AUDIT


def common_event(event_id, master_seed=MASTER_SEED, ctx=None):
    c = ctx or context()
    size_rng = np.random.default_rng([master_seed, event_id, 0])
    position_rng = np.random.default_rng([master_seed, event_id, 1])
    size_state = deepcopy(size_rng.bit_generator.state)
    position_state = deepcopy(position_rng.bit_generator.state)
    diameter = float(c.distribution.inverse_cdf(size_rng.random()))
    position, triangle = c.env.sampler.sample(position_rng)
    event = dict(event_id=int(event_id), particle_id=int(event_id), master_seed=master_seed,
        anchor_m=position[0].tolist(), anchor_triangle=int(triangle[0]),
        first_diameter_um=diameter, first_radius_m=diameter*.5e-6,
        size_rng_before_first=size_state, position_rng_before_anchor=position_state,
        position_rng_after_anchor=deepcopy(position_rng.bit_generator.state),
        B_retry_seed=[master_seed, event_id, 2])
    return event


def sphere_check(point, radius, ctx, *, known=None, owner=None):
    """Exact sphere distance inequality; no change to original P7 checker.

    This avoids constructing a bridge particle at millions of rejected trials.
    Production audits independently re-run original FiniteSizeAdmission.check.
    """
    env = ctx.env
    cell = env.field.locate(point)[0] if owner is None else owner
    if cell < 0:
        return 'NO_OWNER_TETRA', dict(owner_tetra=-1)
    wall, perimeter, _ = known if known is not None else ctx.geometry.distances(point)
    cause = classify_failure(radius=radius, wall_distance=wall,
        perimeter_distance=perimeter, tolerance=env.wall.roundoff_m)
    # A/B use precisely the original WALL + nearfield rule. The perimeter
    # label only refines an already-rejected wall collision; it cannot add a
    # new rejection to the unchanged current rule.
    accepted = wall-radius-lower_gap(radius) >= -env.wall.roundoff_m
    if accepted:
        cause = 'ACCEPTED'
    elif cause == 'ACCEPTED':
        cause = 'OTHER'
    return cause, dict(owner_tetra=int(cell), wall_distance_m=wall,
        perimeter_distance_m=perimeter, wall_gap_m=wall-radius,
        handoff_margin_m=float(wall-radius-lower_gap(radius)))


def method_a(event, ctx=None, guard=GUARD):
    c = ctx or context(); e = deepcopy(event)
    rng = np.random.default_rng(); rng.bit_generator.state = deepcopy(e['position_rng_after_anchor'])
    point = np.asarray(e['anchor_m']); triangle = e['anchor_triangle']; radius = e['first_radius_m']
    attempts = []; accepted = False
    for draw in range(guard):
        if draw:
            points, ids = c.env.sampler.sample(rng); point = points[0]; triangle = int(ids[0])
        cause, detail = sphere_check(point, radius, c)
        attempts.append(dict(draw=draw, position_m=point.tolist(), diameter_um=e['first_diameter_um'],
                             triangle=triangle, cause=cause, **detail))
        if cause == 'ACCEPTED':
            accepted = True; break
    return dict(method='A', event_id=e['event_id'], accepted=accepted,
        status='ACCEPTED' if accepted else 'ADMISSION_GUARD_EXHAUSTED',
        anchor_m=e['anchor_m'], birth_center_m=point.tolist() if accepted else None,
        diameter_um=e['first_diameter_um'], radius_m=radius,
        attempts=attempts, attempt_count=len(attempts),
        position_rng_final=deepcopy(rng.bit_generator.state),
        entry_transition_time_s=0., s_birth_m=0.)


def method_b(event, ctx=None, guard=GUARD):
    c = ctx or context(); e = deepcopy(event); point = np.asarray(e['anchor_m'])
    rng = np.random.default_rng(e['B_retry_seed'])
    diameters = np.r_[e['first_diameter_um'], c.distribution.inverse_cdf(rng.random(guard-1))]
    known = c.geometry.distances(point); owner = c.env.field.locate(point)[0]
    attempts = []; accepted = False
    for draw, diameter in enumerate(diameters):
        radius = float(diameter*.5e-6)
        cause, detail = sphere_check(point, radius, c, known=known, owner=owner)
        attempts.append(dict(draw=draw, position_m=point.tolist(), diameter_um=float(diameter),
                             triangle=e['anchor_triangle'], cause=cause, **detail))
        if cause == 'ACCEPTED':
            accepted = True; break
    return dict(method='B', event_id=e['event_id'], accepted=accepted,
        status='ACCEPTED' if accepted else 'DIAGNOSTIC_SIZE_RETRY_GUARD_EXHAUSTED',
        anchor_m=e['anchor_m'], birth_center_m=point.tolist() if accepted else None,
        diameter_um=float(diameters[len(attempts)-1]) if accepted else None,
        radius_m=float(diameters[len(attempts)-1]*.5e-6) if accepted else None,
        attempts=attempts, attempt_count=len(attempts),
        size_distribution_role='POSITION_CONDITIONED_ACCEPTED_DISTRIBUTION_NOT_ORIGINAL_SONOVUE',
        B_retry_rng_generated_count=guard-1, entry_transition_time_s=0., s_birth_m=0.)


def method_c(event, point_path, ctx=None, *, tolerance=1e-9, horizon_factor=1.):
    c = ctx or context(); e = deepcopy(event); radius = e['first_radius_m']
    aperture_radius = c.geometry.aperture_radius(e['anchor_m'])
    result = first_admissible(point_path, lambda p: c.geometry.full_margin(p, radius),
        horizon=c.geometry.search_horizon*horizon_factor, tolerance=tolerance)
    accepted = result['status'] == 'ACCEPTED'
    return dict(method='C', event_id=e['event_id'], accepted=accepted,
        anchor_m=e['anchor_m'], birth_center_m=result.get('center_m'),
        diameter_um=e['first_diameter_um'], radius_m=radius,
        attempt_count=1, attempts=[], path_definition='FROZEN_P1_VELOCITY_STREAMLINE_POLYLINE',
        normal_fallback=False, anchor_aperture_passable=bool(radius <= aperture_radius),
        anchor_aperture_max_radius_m=aperture_radius,
        transition_role='ENTRY_REPRESENTATION_TRANSITION_NOT_A_FINITE_SIZE_DYNAMICAL_TRAJECTORY',
        physical_crossing_certified=False, **result)


def audit_event(event, ctx=None, *, point_step=.2e-6, search_tolerance=1e-9, horizon_factor=1.):
    c = ctx or context(); e = deepcopy(event)
    trace = c.native.trace(e['anchor_m'], step_m=point_step)
    path = trace.pop('path'); e['point_tracer_basin'] = trace['outlet'] or 'UNRESOLVED_POINT_PATH'
    e['point_tracer_end_reason'] = trace['end_reason']
    wall, perimeter, cap = c.geometry.distances(e['anchor_m'], full=True)
    cause, detail = sphere_check(e['anchor_m'], e['first_radius_m'], c, known=(wall, perimeter, cap))
    e.update(wall_distance_m=wall, perimeter_distance_m=perimeter, cap_distance_m=cap,
        Dmax_um=2*min(wall, perimeter)*1e6,
        geometry_pass_probability=float(c.distribution.cdf(2*min(wall, perimeter)*1e6)),
        aperture_passable=bool(e['first_radius_m'] <= min(wall, perimeter)),
        center_on_plane_accepted=cause == 'ACCEPTED', current_first_cause=cause,
        full_domain_at_anchor=bool(cause == 'ACCEPTED' and cap >= e['first_radius_m']))
    A = method_a(e, c); B = method_b(e, c)
    C = method_c(e, path, c, tolerance=search_tolerance, horizon_factor=horizon_factor)
    for result in (A, B, C):
        result['anchor_basin'] = e['point_tracer_basin']
        if not result['accepted']:
            result['birth_point_basin'] = None
        elif result['method'] == 'B' or np.array_equal(result['birth_center_m'], e['anchor_m']):
            result['birth_point_basin'] = e['point_tracer_basin']
        else:
            birth_trace = c.native.trace(result['birth_center_m'], step_m=point_step)
            result['birth_point_basin'] = birth_trace['outlet'] or 'UNRESOLVED_POINT_PATH'
        result['common_event_sha256'] = hashlib.sha256(__import__('json').dumps(e, sort_keys=True).encode()).hexdigest()
    # Retain only the inlet-neighborhood part for C replay; full MB trajectories
    # are separately integrated with the original P6.5 stepper.
    arc = np.r_[0., np.cumsum(np.linalg.norm(np.diff(path[:, 1:], axis=0), axis=1))]
    n = min(len(path), int(np.searchsorted(arc, c.geometry.search_horizon*horizon_factor))+2)
    return e, (A, B, C), path[:n]
