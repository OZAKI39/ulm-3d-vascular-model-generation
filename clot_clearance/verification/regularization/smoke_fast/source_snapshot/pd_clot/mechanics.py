"""Energy-consistent NOSB correspondence and provisional non-affine penalty.

Unordered ij bond: xi=Xj-Xi, eta=xj-xi. Positive pair force is on i;
the exact negative is added to j. Forces are N, energy density is Pa.
"""
import numpy as np
from numba import njit


@njit(cache=True)
def shape_data(n, pairs, xi, weight, integrity, volume):
    K = np.zeros((n, 3, 3))
    degree = np.zeros(n, dtype=np.int64)
    total = np.zeros(n); intact = np.zeros(n)
    for k in range(len(pairs)):
        i, j = pairs[k]
        q = weight[k] * integrity[k]
        for a in range(3):
            for b in range(3):
                z = q * xi[k, a] * xi[k, b]
                K[i, a, b] += z * volume[j]
                K[j, a, b] += z * volume[i]
        total[i] += volume[j]; total[j] += volume[i]
        intact[i] += integrity[k] * volume[j]
        intact[j] += integrity[k] * volume[i]
        if integrity[k] > 0:
            degree[i] += 1; degree[j] += 1
    return K, degree, intact / total


def prepare_shape(cloud, integrity, safety):
    K, degree, remaining = shape_data(len(cloud.X), cloud.pairs, cloud.xi,
                                      cloud.weight, integrity, cloud.volume)
    eig = np.linalg.eigvalsh(K)
    active = degree > 0
    bad = active & ((degree < safety['minimum_neighbors']) |
                    (eig[:, 0] <= eig[:, 2] / safety['maximum_shape_condition']))
    if bad.any():
        ids = np.flatnonzero(bad)
        raise FloatingPointError(f'Insufficient 3-D shape support at particles {ids[:20].tolist()}; '
                                 'no pseudoinverse, erosion or hidden bond deletion applied')
    inv = np.zeros_like(K)
    inv[active] = np.linalg.inv(K[active])
    trace = np.trace(K, axis1=1, axis2=2)
    return inv, trace, remaining, active


@njit(cache=True)
def evaluate(x, X, pairs, xi, weight, integrity, volume, invK, traceK,
             remaining, active, shear, bulk, stabilization, J_min):
    n = len(x); A = np.zeros((n, 3, 3))
    for k in range(len(pairs)):
        i, j = pairs[k]; q = weight[k] * integrity[k]
        for a in range(3):
            for b in range(3):
                t = q * (x[j, a] - x[i, a]) * xi[k, b]
                A[i, a, b] += t * volume[j]
                A[j, a, b] += t * volume[i]
    F = np.zeros_like(A); J = np.full(n, np.nan)
    PK = np.zeros_like(A); energy = np.zeros(n); hg = np.zeros(n)
    coefficient = np.zeros(n)
    lame = bulk - 2 * shear / 3
    for i in range(n):
        if not active[i]:
            continue  # isolated point: no continuum F or material stress exists
        F[i] = A[i] @ invK[i]
        det = np.linalg.det(F[i]); J[i] = det
        if not np.isfinite(det) or det <= J_min:
            raise ValueError('Invalid deformation Jacobian: stop, inspect FAILED.json')
        logJ = np.log(det)
        P = shear * (F[i] - np.linalg.inv(F[i]).T) + lame * logJ * np.linalg.inv(F[i]).T
        PK[i] = remaining[i] * P @ invK[i]
        energy[i] = remaining[i] * (.5 * shear * (np.sum(F[i]**2) - 3)
                                     - shear * logJ + .5 * lame * logJ**2)
        coefficient[i] = stabilization * shear * remaining[i] / traceK[i]
    force = np.zeros((n, 3))
    for k in range(len(pairs)):
        i, j = pairs[k]; q = weight[k] * integrity[k]
        if q == 0:
            continue
        eta = x[j] - x[i]
        ri = eta - F[i] @ xi[k]; rj = eta - F[j] @ xi[k]
        fij = q * volume[i] * volume[j] * (
            (PK[i] + PK[j]) @ xi[k] + coefficient[i] * ri + coefficient[j] * rj)
        force[i] += fij; force[j] -= fij
        hg[i] += .5 * coefficient[i] * q * np.dot(ri, ri) * volume[j]
        hg[j] += .5 * coefficient[j] * q * np.dot(rj, rj) * volume[i]
    return force, F, J, energy, hg


def evaluate_cloud(cloud, x, integrity, material, safety, shape=None):
    if shape is None:
        shape = prepare_shape(cloud, integrity, safety)
    return evaluate(x, cloud.X, cloud.pairs, cloud.xi, cloud.weight, integrity,
                    cloud.volume, *shape, material['shear_modulus_Pa'],
                    material['bulk_modulus_Pa'], material['stabilization_alpha'],
                    safety['minimum_J'])


def timestep_limit(cloud, c):
    m = c['material']; rho = c['clot']['density_kg_m3']
    speed = np.sqrt((m['bulk_modulus_Pa'] + (4/3 + 2*m['stabilization_alpha']) *
                     m['shear_modulus_Pa']) / rho)
    return c['safety']['dt_safety_factor'] * cloud.spacing / speed
