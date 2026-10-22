"""Simultaneous frictionless minimum kinematic correction, never a force law."""
from dataclasses import dataclass, asdict
import numpy as np
from .particle_shapes import Ellipsoid, EPS, unit

CONTACT_MODEL = 'HARD_FRICTIONLESS_PARTICLE_CONTACT_V0'
CONTACT_METRIC = 'GEOMETRY_SCALED_KINEMATIC_MINIMUM_CORRECTION'


class MultiContactInfeasible(RuntimeError):
    def __init__(self, reason, record=None):
        self.record = dict(status='MULTI_CONTACT_INFEASIBLE', reason=reason, diagnostics=record)
        super().__init__(str(self.record))


@dataclass(frozen=True)
class ContactConstraint:
    canonical_id: tuple
    particle_i_id: int
    particle_j_id: int | None
    normal: np.ndarray
    point_i_m: np.ndarray
    point_j_m: np.ndarray | None = None
    offset_m_s: float = 0.

    @classmethod
    def pair(cls, gap):
        if gap.state != 'TOUCHING':
            raise ValueError('Only touching pairs can enter velocity projection')
        return cls(('PAIR', *gap.canonical_pair_id), gap.particle_i_id, gap.particle_j_id,
                   gap.normal_j_to_i, gap.point_i_m, gap.point_j_m)

    @classmethod
    def wall(cls, particle_id, gap):
        if isinstance(gap, tuple):
            triangle_id, gap = gap
        else:
            triangle_id = gap.wall_triangle_id
        if abs(gap.gap_m) > gap.roundoff_m:
            raise ValueError('Only touching WALL features can enter projection')
        return cls(('WALL', particle_id, triangle_id, gap.wall_feature), particle_id,
                   None, gap.normal_inward, gap.particle_point_m)


@dataclass
class ContactProjection:
    velocities: dict
    omegas: dict
    translation_corrections: dict
    angular_corrections: dict
    record: dict


def _dual_active_set(jacobian, speeds, velocity_scale):
    count = len(speeds)
    if not count:
        return np.empty(0), np.empty(0), dict(iterations=0, condition=1., velocity_budget_m_s=0., multiplier_budget=0., complementarity_budget=0.)
    gram = jacobian @ jacobian.T
    singular = np.linalg.svd(jacobian, compute_uv=False)
    retained = singular[singular > max(jacobian.shape)*EPS*singular[0]]
    condition = max(1., float(singular[0]/retained[-1])) if len(retained) else 1.
    # Explicit operation-count / conditioning roundoff budget, no physical floor.
    budget = 512*EPS*max(1, count, jacobian.shape[1])*condition*velocity_scale
    multiplier_budget = budget/max(float(np.max(np.diag(gram))), np.finfo(float).tiny)
    multipliers = np.zeros(count); active = []; seen = set()
    max_iterations = 16*(count+1)**2  # implementation guard, not a model parameter
    for iteration in range(max_iterations):
        speeds_new = speeds + gram@multipliers
        inactive = [k for k in range(count) if k not in active and speeds_new[k] < -budget]
        if not inactive:
            break
        entering = min(inactive, key=lambda k: (speeds_new[k], k)); active.append(entering); active.sort()
        for _ in range(max_iterations):
            candidate = np.zeros(count)
            matrix = gram[np.ix_(active, active)]
            candidate[active] = np.linalg.lstsq(matrix, -speeds[active], rcond=max(matrix.shape)*EPS)[0]
            residual = matrix@candidate[active]+speeds[active]
            if np.max(np.abs(residual), initial=0.) > budget:
                raise MultiContactInfeasible('Rank-deficient active constraints are incompatible', dict(active=active, residual=residual.tolist()))
            negative = [k for k in active if candidate[k] < -multiplier_budget]
            if not negative:
                multipliers = np.maximum(candidate, 0.); break
            alpha = min(multipliers[k]/(multipliers[k]-candidate[k]) for k in negative)
            multipliers += alpha*(candidate-multipliers)
            leaving = [k for k in active if multipliers[k] <= multiplier_budget and candidate[k] < 0]
            if not leaving:
                raise MultiContactInfeasible('Active-set step did not remove a negative multiplier')
            for k in leaving:
                active.remove(k); multipliers[k] = 0.
        signature = (tuple(active), tuple(np.flatnonzero(multipliers > multiplier_budget)))
        if signature in seen and np.min(speeds+gram@multipliers) < -budget:
            raise MultiContactInfeasible('Active-set cycling guard')
        seen.add(signature)
    else:
        raise MultiContactInfeasible('Active-set iteration guard')
    speeds_new = speeds+gram@multipliers
    complementarity_budget = budget*max(velocity_scale, np.max(np.abs(multipliers), initial=0.))
    if np.min(speeds_new) < -budget or np.max(np.abs(speeds_new*multipliers)) > complementarity_budget:
        raise MultiContactInfeasible('Final KKT audit failed')
    return multipliers, speeds_new, dict(iterations=iteration+1, condition=condition,
        velocity_budget_m_s=budget, multiplier_budget=multiplier_budget, complementarity_budget=complementarity_budget)


def project_contacts(shapes, free_velocities, free_omegas, constraints):
    """All WALL and pair constraints together; stable IDs never depend on position."""
    ids = sorted(shapes); constraints = sorted(constraints, key=lambda c: c.canonical_id)
    if len(set(c.canonical_id for c in constraints)) != len(constraints):
        raise ValueError('Duplicate canonical contact ID')
    slices = {}; values = []; cursor = 0
    for i in ids:
        v = np.asarray(free_velocities[i], dtype=float); omega = np.asarray(free_omegas[i], dtype=float)
        if v.shape != (3,) or omega.shape != (3,) or not np.isfinite(np.r_[v, omega]).all():
            raise ValueError('Finite vector velocities required')
        active = isinstance(shapes[i], Ellipsoid); length = shapes[i].bounding_radius_m
        slices[i] = (cursor, active, length); values.extend(v)
        if active:
            values.extend(length*omega)
        cursor += 6 if active else 3
    free = np.asarray(values); jacobian = np.zeros((len(constraints), cursor)); offsets = []
    for row, c in enumerate(constraints):
        n = unit(c.normal); offsets.append(c.offset_m_s)
        for particle_id, sign, point in [(c.particle_i_id, 1., c.point_i_m), (c.particle_j_id, -1., c.point_j_m)]:
            if particle_id is None:
                continue
            start, active, length = slices[particle_id]
            jacobian[row, start:start+3] = sign*n
            if active:
                jacobian[row, start+3:start+6] = sign*np.cross(point-shapes[particle_id].center_m, n)/length
    speeds = jacobian@free + offsets
    scale = max(np.linalg.norm(free), np.max(np.abs(speeds), initial=0.), np.finfo(float).tiny)
    multipliers, corrected_speeds, audit = _dual_active_set(jacobian, speeds, scale)
    correction = jacobian.T@multipliers; corrected = free+correction
    velocities = {}; omegas = {}; dv = {}; dw = {}
    for i in ids:
        start, active, length = slices[i]
        velocities[i] = corrected[start:start+3]; dv[i] = correction[start:start+3]
        dw[i] = correction[start+3:start+6]/length if active else np.zeros(3)
        omegas[i] = np.asarray(free_omegas[i])+dw[i]
    pair_contributions = {i:np.zeros(3) for i in ids}; contributions = []
    for k, c in enumerate(constraints):
        value = multipliers[k]*unit(c.normal)
        if c.particle_j_id is not None:
            pair_contributions[c.particle_i_id] += value
            pair_contributions[c.particle_j_id] -= value
        contributions.append(dict(canonical_id=c.canonical_id, multiplier=multipliers[k],
                                  particle_i_translation=value, particle_j_translation=-value if c.particle_j_id is not None else None))
    # Audit the actual corrected point velocities, independently of Gram-matrix
    # multiplication used inside the dual solve.
    actual_speeds=[]
    for c in constraints:
        speed=float(c.offset_m_s)
        for i,sign,point in [(c.particle_i_id,1.,c.point_i_m),(c.particle_j_id,-1.,c.point_j_m)]:
            if i is None:continue
            v=velocities[i].copy()
            if isinstance(shapes[i],Ellipsoid):v+=np.cross(omegas[i],point-shapes[i].center_m)
            speed+=sign*float(unit(c.normal)@v)
        actual_speeds.append(speed)
    actual_speeds=np.asarray(actual_speeds)
    if np.min(actual_speeds,initial=0.) < -audit['velocity_budget_m_s']:
        raise MultiContactInfeasible('Actual contact-point velocity audit failed')
    pair_balance=sum(pair_contributions.values(),np.zeros(3))
    record = dict(status='PASS', solver='DETERMINISTIC_DUAL_ACTIVE_SET_V0', metric=CONTACT_METRIC,
        multiplier_role='KINEMATIC_CONTACT_MULTIPLIER', physical_force_claimed=False,
        canonical_ids=[c.canonical_id for c in constraints], particle_ids=ids,
        free_normal_speeds_m_s=speeds, multipliers=multipliers, corrected_normal_speeds_m_s=actual_speeds,
        dual_matrix_normal_speeds_m_s=corrected_speeds,
        max_primal_violation_m_s=float(max(0., -np.min(actual_speeds, initial=0.))),
        max_dual_violation=float(max(0., -np.min(multipliers, initial=0.))),
        max_complementarity_residual=float(np.max(np.abs(multipliers*actual_speeds), initial=0.)),
        stationarity_residual=float(np.linalg.norm(corrected-free-jacobian.T@multipliers)),
        pair_translation_balance_error_m_s=float(np.linalg.norm(pair_balance)),
        total_translation_correction_m_s=sum(dv.values(), np.zeros(3)), constraint_contributions=contributions, **audit)
    return ContactProjection(velocities, omegas, dv, dw, record)
