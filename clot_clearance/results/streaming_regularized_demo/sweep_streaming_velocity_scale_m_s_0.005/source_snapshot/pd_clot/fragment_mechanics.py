"""Explicit supported-subspace extension; original full-rank kernel is unchanged.

Fracture may leave line/surface families. Their constitutive energy is defined
on the reference subspace spanned by SURVIVING bonds. This is a disclosed
prototype approximation, not a pseudoinverse pretending to be a 3-D F.
"""
import numpy as np
from numba import njit
from .mechanics import shape_data, evaluate_cloud


def prepare_supported_shape(cloud, integrity, safety):
    K, degree, remaining = shape_data(len(cloud.X), cloud.pairs, cloud.xi,
                                      cloud.weight, integrity, cloud.volume)
    eig, basis = np.linalg.eigh(K)
    tolerance = safety['support_rank_relative_tolerance']
    supported = eig > tolerance * np.maximum(eig[:, -1:], np.finfo(float).tiny)
    ranks = supported.sum(axis=1).astype(np.int64)
    if safety['rank_deficient_policy'] == 'error' and np.any((degree > 0) & (ranks < 3)):
        raise FloatingPointError('Reduced-dimensional active family; explicit extension disabled')
    if safety['rank_deficient_policy'] not in ['error', 'intrinsic_correspondence']:
        raise ValueError('Unknown rank-deficient policy')
    for i in np.flatnonzero(ranks):
        e = eig[i, supported[i]]
        if e[-1] / e[0] > safety['maximum_shape_condition']:
            raise FloatingPointError(f'Ill-conditioned supported family at particle {i}')
    reciprocal = np.zeros_like(eig)
    np.divide(1, eig, out=reciprocal, where=supported)
    inverse = np.einsum('nij,nj,nkj->nik', basis, reciprocal, basis)
    trace = np.trace(K, axis1=1, axis2=2)
    return inverse, trace, remaining, ranks, basis


@njit(cache=True)
def evaluate_intrinsic(x, pairs, xi, weight, integrity, volume, inverse, trace,
                       remaining, ranks, basis, shear, bulk, alpha, minimum_J):
    n = len(x); A = np.zeros((n, 3, 3))
    for k in range(len(pairs)):
        i, j = pairs[k]; q = weight[k] * integrity[k]
        for a in range(3):
            for b in range(3):
                t = q * (x[j, a] - x[i, a]) * xi[k, b]
                A[i, a, b] += t * volume[j]
                A[j, a, b] += t * volume[i]
    F = np.zeros_like(A); J = np.full(n, np.nan)
    PK = np.zeros_like(A); energy = np.zeros(n); hg = np.zeros(n); coefficient = np.zeros(n)
    lame = bulk - 2 * shear / 3
    for i in range(n):
        rank = ranks[i]
        if rank == 0: continue
        F[i] = A[i] @ inverse[i]
        if rank == 3:
            jac = np.linalg.det(F[i]); J[i] = jac
            if not np.isfinite(jac) or jac <= minimum_J:
                raise ValueError('Invalid full-rank deformation Jacobian')
            logJ = np.log(jac); invT = np.linalg.inv(F[i]).T
            P = shear * (F[i] - invT) + lame * logJ * invT
            energy[i] = remaining[i] * (.5*shear*(np.sum(F[i]**2)-3) - shear*logJ + .5*lame*logJ**2)
        else:
            U = np.ascontiguousarray(basis[i, :, 3-rank:])
            B = F[i] @ U
            C = B.T @ B
            jac = np.sqrt(np.linalg.det(C)); J[i] = jac
            if not np.isfinite(jac) or jac <= minimum_J:
                raise ValueError('Collapsed intrinsic line/area stretch')
            logJ = np.log(jac)
            Psub = shear*B + (lame*logJ-shear)*(B @ np.linalg.inv(C))
            P = Psub @ U.T
            energy[i] = remaining[i] * (.5*shear*(np.sum(B**2)-rank) - shear*logJ + .5*lame*logJ**2)
        PK[i] = remaining[i] * P @ inverse[i]
        coefficient[i] = alpha * shear * remaining[i] / trace[i]
    force = np.zeros((n,3))
    for k in range(len(pairs)):
        i, j = pairs[k]; q = weight[k]*integrity[k]
        if q == 0: continue
        eta = x[j]-x[i]
        ri = eta-F[i]@xi[k]; rj = eta-F[j]@xi[k]
        fij = q*volume[i]*volume[j]*((PK[i]+PK[j])@xi[k] + coefficient[i]*ri + coefficient[j]*rj)
        force[i] += fij; force[j] -= fij
        hg[i] += .5*coefficient[i]*q*np.dot(ri,ri)*volume[j]
        hg[j] += .5*coefficient[j]*q*np.dot(rj,rj)*volume[i]
    return force,F,J,energy,hg


def evaluate_supported(cloud, x, integrity, material, safety, shape):
    inverse, trace, remaining, ranks, basis = shape
    if np.all(ranks == 3):
        # Exact original implementation, not a replacement continuum kernel.
        return evaluate_cloud(cloud,x,integrity,material,safety,
                              (inverse,trace,remaining,ranks>0))
    return evaluate_intrinsic(x,cloud.pairs,cloud.xi,cloud.weight,integrity,cloud.volume,
        inverse,trace,remaining,ranks,basis,material['shear_modulus_Pa'],material['bulk_modulus_Pa'],
        material['stabilization_alpha'],safety['minimum_J'])
