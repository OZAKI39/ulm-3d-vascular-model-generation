"""CPU-only, opt-in F10G references. No production solver imports or defaults changed.

A is the *existing* dimensionless, pair-normal resistance. zeta*A has units
N s/m. Its inverse / zeta is used only as a <=96-DOF diagnostic mobility.
"""

from time import perf_counter
import numpy as np
from scipy import linalg
from scipy.spatial.distance import pdist, squareform


def normals(x):
    x = np.asarray(x, dtype=float)
    i, j = np.triu_indices(len(x), 1)
    dx = x[j] - x[i]
    d = np.linalg.norm(dx, axis=1)
    if np.any(d <= 0):
        raise ValueError("Coincident centers are not a valid contact normal")
    return i, j, d, dx / d[:, None]


def resistance(x, radius, activation=0.1, floor=0.001):
    i, j, d, n = normals(x)
    e = (d - 2 * radius) / radius
    mask = e < activation
    i, j, d, n, e = i[mask], j[mask], d[mask], n[mask], e[mask]

    def phi(v):
        return 0.5 / v - 0.45 * np.log(v) - (3 / 56) * v * np.log(v)

    c = 0.5 * (phi(np.maximum(e, floor)) - phi(activation))
    a = np.eye(3 * len(x))
    coupling = np.zeros((len(x), len(x), 3, 3))
    parent = np.arange(len(x))

    def root(k):
        while parent[k] != k:
            k = parent[k]
        return k

    for u, v, normal, coef in zip(i, j, n, c):
        b = coef * np.outer(normal, normal)
        s, t = slice(3 * u, 3 * u + 3), slice(3 * v, 3 * v + 3)
        a[s, s] += b
        a[t, t] += b
        a[s, t] -= b
        a[t, s] -= b
        coupling[u, v] = coupling[v, u] = b
        ru, rv = root(u), root(v)
        parent[max(ru, rv)] = min(ru, rv)
    groups = [np.array([j for j in range(len(x)) if root(j) == k]) for k in range(len(x)) if root(k) == k]
    blocks = np.array([a[3 * k : 3 * k + 3, 3 * k : 3 * k + 3] for k in range(len(x))])
    return dict(
        A=a,
        blocks=blocks,
        coupling=coupling,
        groups=groups,
        pairs=np.c_[i, j],
        normals=n,
        gaps_m=d - 2 * radius,
        coefficients=c,
        floored=e < floor,
    )


def jacobi(op, free, count):
    u = np.array(free, copy=True)
    for _ in range(count):
        rhs = free + np.einsum("ijkl,jl->ik", op["coupling"], u)
        u = 0.5 * (u + np.linalg.solve(op["blocks"], rhs[..., None])[..., 0])
    return u


def balance(op, free, u):
    v = u.copy()
    for ids in op["groups"]:
        v[ids] += np.mean(free[ids] - v[ids], axis=0)
    return v


def overlaps(x, radius):
    return np.maximum(0, 2 * radius - pdist(x))


def pbd(prediction, radius, rounds=12, relaxation=0.5):
    x = prediction.copy()
    records = []
    for k in range(rounds):
        i, j, d, n = normals(x)
        impulse = 0.5 * relaxation * np.maximum(0, 2 * radius - d)[:, None] * n
        dx = np.zeros_like(x)
        np.add.at(dx, i, -impulse)
        np.add.at(dx, j, impulse)
        x += dx
        records.append(
            dict(
                round=k + 1,
                positions_m=x.copy(),
                correction_m=dx,
                residual_m=float(overlaps(x, radius).max(initial=0)),
                max_correction_m=float(np.linalg.norm(dx, axis=1).max()),
            )
        )
    return x, records


def measure(fn, repeats=25):
    fn()  # explicit one untimed warmup; no trajectory is advanced
    samples = []
    for _ in range(repeats):
        start = perf_counter()
        fn()
        samples.append(perf_counter() - start)
    return dict(
        repeats=repeats,
        median_s=float(np.median(samples)),
        min_s=min(samples),
        max_s=max(samples),
        samples_s=samples,
    )


def shape(x, radius):
    y = x - np.mean(x, axis=0)
    tensor = y.T @ y / len(x)
    eig = np.linalg.eigvalsh(tensor)
    adjacency = squareform(pdist(x) <= 2.1 * radius)
    np.fill_diagonal(adjacency, True)
    members = []
    remaining = set(range(len(x)))
    while remaining:
        todo = [min(remaining)]
        visited = set()
        while todo:
            i = todo.pop()
            if i in visited:
                continue
            visited.add(i)
            todo.extend(int(j) for j in np.flatnonzero(adjacency[i]) if j not in visited)
        remaining -= visited
        members.append(sorted(visited))
    return dict(
        centroid_m=x.mean(axis=0),
        Rg_m=float(np.sqrt(np.trace(tensor))),
        shape_tensor_m2=tensor,
        aspect_ratio=float(np.sqrt(eig[-1] / eig[0])) if eig[0] > 0 else None,
        members=members,
        edges=np.argwhere(np.triu(adjacency, 1)),
    )


def compare(x, ref):
    a = x - x.mean(axis=0)
    b = ref - ref.mean(axis=0)
    u, _, v = np.linalg.svd(a.T @ b)
    rot = u @ np.diag([1, 1, np.linalg.det(u @ v)]) @ v
    return dict(
        centroid_delta_m=x.mean(axis=0) - ref.mean(axis=0),
        centroid_norm_m=float(np.linalg.norm(x.mean(axis=0) - ref.mean(axis=0))),
        max_fixed_ID_m=float(np.linalg.norm(x - ref, axis=1).max()),
        Kabsch_RMSD_m=float(np.sqrt(np.mean(np.sum((a @ rot - b) ** 2, axis=1)))),
    )


def contact_geometry(x, radius):
    """All pairs: no missed neighbor that another contact could activate.

    Rows of B are distance gradients (-n_i,+n_j); positive B dx separates.
    Normals are frozen at the initial state, not claimed nonlinear BE normals.
    """
    i, j, d, n = normals(x)
    B = np.zeros((len(i), 3 * len(x)))
    for k, (a, b) in enumerate(zip(i, j)):
        B[k, 3 * a : 3 * a + 3] = -n[k]
        B[k, 3 * b : 3 * b + 3] = n[k]
    return B, np.c_[i, j], d - 2 * radius


def soft_contact(x, free, radius, zeta, dt, kc, A=None, max_iterations=24, atol_m=1e-13):
    """Bounded CPU semismooth Newton prototype for a frozen-normal spring.

    Minimizes .5 q^T A q + .5 alpha ||(b-Bq)_+||^2 in q=dx_contact/R,
    b=predicted signed overlap/R, alpha=dt*kc/zeta. A=I for B1; actual
    resistance for B2. Neither an inertial XPBD update nor production QP.
    """
    start = perf_counter()
    if dt <= 0 or kc <= 0 or radius <= 0 or zeta <= 0:
        raise ValueError("Positive physical parameters required")
    x = np.asarray(x, dtype=float)
    free = np.asarray(free, dtype=float)
    A = np.eye(x.size) if A is None else np.array(A)
    if (
        not np.all(np.isfinite(x))
        or not np.all(np.isfinite(free))
        or not np.allclose(A, A.T, rtol=0, atol=1e-12)
    ):
        raise ValueError("Finite inputs and symmetric resistance required")
    factor = linalg.cho_factor(A)
    u = linalg.cho_solve(factor, free.ravel())
    B, pairs, gaps = contact_geometry(x, radius)
    b = (-gaps - dt * (B @ u)) / radius
    alpha = dt * kc / zeta
    q = np.zeros(x.size)
    init_s = perf_counter() - start

    def objective(q):
        overlap = np.maximum(0, b - B @ q)
        return 0.5 * (q @ A @ q + alpha * (overlap @ overlap))

    records = []
    converged = False
    begin = perf_counter()
    for iteration in range(max_iterations + 1):
        overlap = np.maximum(0, b - B @ q)
        gradient = A @ q - alpha * B.T @ overlap
        residual = linalg.cho_solve(factor, gradient) * radius
        residual_m = float(np.linalg.norm(residual.reshape(-1, 3), axis=1).max())
        if residual_m <= atol_m:
            converged = True
            break
        if iteration == max_iterations:
            break
        active = overlap > 0
        H = A + alpha * B[active].T @ B[active]
        direction = linalg.solve(H, -gradient, assume_a="pos")
        old = objective(q)
        rate = 1.0
        for _ in range(20):
            if objective(q + rate * direction) <= old + 1e-4 * rate * np.dot(gradient, direction):
                break
            rate *= 0.5
        else:
            break
        dx = rate * direction * radius
        q += rate * direction
        records.append(
            dict(
                iteration=iteration + 1,
                residual_before_m=residual_m,
                line_search_rate=rate,
                max_correction_m=float(np.linalg.norm(dx.reshape(-1, 3), axis=1).max()),
                correction_m=dx.reshape(-1, 3),
                objective=objective(q),
            )
        )
    contact_dx = q.reshape(-1, 3) * radius
    positions = x + dt * u.reshape(-1, 3) + contact_dx
    delta = np.maximum(0, b - B @ q) * radius
    reaction = kc * delta
    force = (B.T @ reaction).reshape(-1, 3)
    velocity = u.reshape(-1, 3) + linalg.cho_solve(factor, force.ravel() / zeta).reshape(-1, 3)
    consistent = positions - x - dt * velocity
    return dict(
        positions_m=positions,
        prediction_m=x + dt * u.reshape(-1, 3),
        free_m_s=free,
        mobility_free_m_s=u.reshape(-1, 3),
        corrected_velocity_m_s=velocity,
        contact_displacement_m=contact_dx,
        contact_force_N=force,
        pair_force_N=reaction,
        pairs=pairs,
        frozen_normal_overlap_m=delta,
        geometric_overlap_m=overlaps(positions, radius),
        relative_normal_velocity_m_s=B @ velocity.ravel(),
        mobility_equation_residual_m=residual_m,
        force_position_consistency_m=float(np.max(np.abs(consistent))),
        net_contact_force_N=np.sum(force, axis=0),
        centroid_contact_shift_m=contact_dx.mean(axis=0),
        converged=converged,
        iterations=records,
        max_single_iteration_m=max((r["max_correction_m"] for r in records), default=0.0),
        max_total_contact_correction_m=float(np.linalg.norm(contact_dx, axis=1).max()),
        initialization_s=init_s,
        contact_s=perf_counter() - begin,
        total_s=perf_counter() - start,
        kc_N_m=kc,
        dt_s=dt,
        model="frozen initial normals and mobility; unilateral overdamped BE spring",
    )


def enumerated_reference(x, free, radius, zeta, dt, kc, A=None):
    """Independent <=4-bubble dual active-set enumeration (<=64 subsets).

    Solve f/kc + dt B M B^T f = b on each guessed active set, then
    check the unilateral complementarity conditions on ALL pairs.
    This independently solves force multipliers rather than primal Newton.
    """
    from itertools import product

    if len(x) > 4:
        raise ValueError("Exhaustive reference limited to four bubbles")
    A = np.eye(x.size) if A is None else A
    B, pairs, gaps = contact_geometry(x, radius)
    D = np.linalg.inv(A)
    pred = dt * (D @ free.ravel())
    b = (-gaps - B @ pred) / radius
    S = dt * kc / zeta * (B @ D @ B.T)
    for flags in product([False, True], repeat=len(pairs)):
        active = np.array(flags)
        f = np.zeros(len(pairs))  # f scaled by kc*radius
        if np.any(active):
            f[active] = np.linalg.solve(np.eye(active.sum()) + S[np.ix_(active, active)], b[active])
        slack = b - S @ f
        if np.all(f[active] >= -1e-12) and np.all(slack[~active] <= 1e-12):
            force = B.T @ (f * kc * radius)
            positions = x + (pred + dt / zeta * D @ force).reshape(-1, 3)
            return dict(
                positions_m=positions,
                pair_force_N=f * kc * radius,
                unilateral_max_error=float(np.max(np.abs(f - np.maximum(slack, 0)))),
            )
    raise RuntimeError("Independent complementarity reference found no valid set")
