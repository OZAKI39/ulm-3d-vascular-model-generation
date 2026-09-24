"""P4 Minkowski support geometry using the immutable P3 convex shapes.

GJK supplies the separating direction. Smooth support-stationary refinement
removes polytope approximation error; penetration enumerates deterministic
stationary basins and every capsule segment sign / kink subspace. No effective
RBC radius or point-cloud geometry is used. Failure to close the witness
identity is an explicit geometry error, never an invented signed distance.
"""
from dataclasses import dataclass, asdict
from itertools import combinations, product
import numpy as np
from scipy.optimize import minimize, lsq_linear
from scipy.linalg import null_space
from .particle_shapes import Sphere, Ellipsoid, Capsule, EPS, roundoff_length, unit

PAIR_GAP_ENGINE = 'GJK_MINKOWSKI_SUPPORT_STATIONARY_V0'


class PairGeometryError(RuntimeError):
    pass


@dataclass(frozen=True)
class ParticlePairGapResult:
    particle_i_id: int
    particle_j_id: int
    shape_i: str
    shape_j: str
    gap_m: float
    state: str
    point_i_m: np.ndarray
    point_j_m: np.ndarray
    normal_j_to_i: np.ndarray
    distance_m: float
    geometry_iterations: int
    roundoff_budget_m: float
    canonical_pair_id: tuple

    def to_dict(self):
        return asdict(self)


def _simplex(points):
    """Closest convex combination, exhaustive faces of a <=4 vertex simplex."""
    best = None
    for size in range(1, min(4, len(points)) + 1):
        for ids in combinations(range(len(points)), size):
            p = np.array([points[k] for k in ids])
            gram = p @ p.T
            system = np.block([[gram, np.ones((size, 1))], [np.ones((1, size)), np.zeros((1, 1))]])
            weights = np.linalg.lstsq(system, np.r_[np.zeros(size), 1.], rcond=None)[0][:size]
            if min(weights) < -64 * EPS:
                continue
            weights = np.maximum(weights, 0); weights /= weights.sum()
            closest = weights @ p; score = closest @ closest
            if best is None or score < best[0]:
                keep = [k for k, w in zip(ids, weights) if w > 64 * EPS]
                best = score, closest, [points[k] for k in keep]
    if best is None:
        raise PairGeometryError('GJK simplex barycentric failure')
    return best[1:]


def _gjk(a, b, scale, pad):
    delta = (a.center_m - b.center_m) / scale
    def support(n):
        return delta + ((a.support(n)-a.center_m) - (b.support(-n)-b.center_m))/scale
    direction = delta if np.linalg.norm(delta) else np.array([1., 0, 0])
    points = [support(-direction)]; closest = points[0]
    for iteration in range(1, 129):
        length = np.linalg.norm(closest)
        if length <= pad/scale:
            return unit(direction), iteration
        direction = closest/length; point = support(-direction)
        if length - point@direction <= pad/scale:
            return direction, iteration
        if any(np.linalg.norm(point-p) <= 32*EPS for p in points):
            return direction, iteration  # stationary refinement still mandatory
        closest, points = _simplex(points + [point])
    return unit(direction), 128


def _components(shape, scale):
    if isinstance(shape, Ellipsoid):
        return shape.quadratic/scale**2, 0., None
    if isinstance(shape, (Sphere, Capsule)):
        segment = (shape.axis_world, shape.cylindrical_length_m/(2*scale)) if isinstance(shape, Capsule) and shape.cylindrical_length_m else None
        return None, shape.radius_m/scale, segment
    raise ValueError('Only SPHERE_MB, FREE_OBLATE, CAPILLARY_DEFORMED states are simulatable')


def _stationary(d, matrices, radius, basis, seeds):
    """Stationary support directions on a sphere or a capsule kink subspace."""
    dimension = basis.shape[1]
    db = basis.T @ d; qs = [basis.T @ q @ basis for q in matrices]
    if dimension == 1:
        return [basis[:, 0], -basis[:, 0]], 0
    if not qs:
        if np.linalg.norm(db) > 64*EPS:
            return [basis @ (db/np.linalg.norm(db)), -basis @ (db/np.linalg.norm(db))], 0
        return [basis @ np.eye(dimension)[k]*sign for k in range(dimension) for sign in [-1, 1]], 0
    def evaluate(x):
        value = radius - db@x; gradient = -db.copy(); hessian = np.zeros((dimension, dimension))
        for q in qs:
            z = q@x; h = np.sqrt(x@z)
            value += h; gradient += z/h
            hessian += q/h - np.outer(z, z)/h**3
        return value, gradient, hessian
    candidates = []; evaluations = 0
    for seed in seeds:
        x = basis.T @ seed
        if np.linalg.norm(x) <= 64*EPS:
            continue
        x = x/np.linalg.norm(x)
        result = minimize(lambda v: evaluate(v)[0], x, jac=lambda v: evaluate(v)[1],
                          method='SLSQP', constraints={'type': 'eq', 'fun': lambda v: v@v-1,
                          'jac': lambda v: 2*v}, options={'ftol': 8*EPS, 'maxiter': 100})
        evaluations += result.nit; x = result.x/np.linalg.norm(result.x)
        # Newton on the tangent plane: function-value stopping alone loses
        # sqrt(eps) accuracy in contact normals and witness coordinates.
        for _ in range(16):
            _, gradient, hessian = evaluate(x); tangent = null_space(x[None, :])
            residual = tangent.T @ gradient
            if np.linalg.norm(residual) <= 16*EPS*max(1., np.linalg.norm(db)):
                break
            reduced = tangent.T @ (hessian - (x@gradient)*np.eye(dimension)) @ tangent
            step = np.linalg.lstsq(reduced, -residual, rcond=64*EPS)[0]
            if np.linalg.norm(step) > .5:
                break
            x = x + tangent@step; x /= np.linalg.norm(x); evaluations += 1
        candidates.append(basis@x)
    return candidates, evaluations


def _witness(a, b, n, gap, scale, segments):
    points = []
    kink = []
    for index, (shape, direction) in enumerate([(a, -n), (b, n)]):
        point = shape.support(direction)
        if isinstance(shape, Capsule) and shape.cylindrical_length_m:
            dot = direction@shape.axis_world
            if abs(dot) <= 512*EPS:
                point = shape.center_m + shape.radius_m*direction
                kink.append((index, shape.axis_world*shape.cylindrical_length_m/2))
        points.append(point.copy())
    if kink:
        matrix = np.column_stack([(1 if i == 0 else -1)*v/scale for i, v in kink])
        target = (gap*n - (points[0]-points[1]))/scale
        solution = lsq_linear(matrix, target, bounds=(-1., 1.), tol=32*EPS, lsmr_tol=32*EPS, max_iter=1000).x
        for weight, (i, vector) in zip(solution, kink):
            points[i] += weight*vector
    return points


def _ordered_gap(a, b, i, j):
    scale = max(a.bounding_radius_m, b.bounding_radius_m, np.linalg.norm(a.center_m-b.center_m))
    pad = roundoff_length(a.center_m, b.center_m, a.bounding_radius_m, b.bounding_radius_m)
    seed, iterations = _gjk(a, b, scale, pad)
    components = [_components(s, scale) for s in [a, b]]
    matrices = [q for q, _, _ in components if q is not None]
    radius = sum(r for _, r, _ in components)
    segments = [s for _, _, s in components if s is not None]
    d = (a.center_m-b.center_m)/scale
    seeds = [seed]
    full_seeds = [unit(np.array(v, dtype=float)) for v in product([-1, 0, 1], repeat=3) if any(v)]
    full_seeds += [seed, -seed]
    best = None
    for exhaustive in [False, True]:
        for kink_mask in product([False, True], repeat=len(segments)):
            axes = [s[0] for s, kink in zip(segments, kink_mask) if kink]
            basis = null_space(np.array(axes)) if axes else np.eye(3)
            if not basis.shape[1]:
                continue
            free_indices = [k for k, kink in enumerate(kink_mask) if not kink]
            for signs in product([-1, 1], repeat=len(free_indices)):
                adjusted = d.copy()
                for k, sign in zip(free_indices, signs):
                    adjusted -= sign*segments[k][1]*segments[k][0]
                normals, count = _stationary(adjusted, matrices, radius, basis, full_seeds if exhaustive else seeds)
                iterations += count
                for n in normals:
                    if any(sign*(n@segments[k][0]) < -512*EPS for k, sign in zip(free_indices, signs)):
                        continue
                    value = d@n - radius - sum(np.sqrt(n@q@n) for q in matrices) - sum(h*abs(n@axis) for axis, h in segments)
                    if best is None or value > best[0]:
                        best = value, n
        if best is not None:
            gap = best[0]*scale; n = best[1]
            pi, pj = _witness(a, b, n, gap, scale, segments)
            error = np.linalg.norm(pi-pj-gap*n)
            if error <= pad and (gap >= -pad or exhaustive or not matrices):
                break
    if best is None or error > pad:
        raise PairGeometryError(f'Support witness did not close: residual={error if best else None}, budget={pad}')
    state = 'SEPARATED' if gap > pad else 'PENETRATING' if gap < -pad else 'TOUCHING'
    return ParticlePairGapResult(i, j, a.mode, b.mode, float(gap), state, pi, pj, n,
                                float(abs(gap)), iterations, pad, (i, j))


def pair_gap(shape_i, shape_j, particle_i_id=0, particle_j_id=1):
    """Stable IDs canonicalize arithmetic; reversal swaps witnesses and normal."""
    i, j = int(particle_i_id), int(particle_j_id)
    if i == j:
        raise ValueError('Distinct stable particle IDs required')
    if i < j:
        return _ordered_gap(shape_i, shape_j, i, j)
    r = _ordered_gap(shape_j, shape_i, j, i)
    return ParticlePairGapResult(i, j, shape_i.mode, shape_j.mode, r.gap_m, r.state,
        r.point_j_m, r.point_i_m, -r.normal_j_to_i, r.distance_m, r.geometry_iterations,
        r.roundoff_budget_m, r.canonical_pair_id)
