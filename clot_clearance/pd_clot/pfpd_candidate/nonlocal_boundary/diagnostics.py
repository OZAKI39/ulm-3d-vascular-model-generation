"""Independent continuum targets, affine work, support and field diagnostics."""
import numpy as np
from ..svk_spectral_split import svk_split_energy, green_lagrange


def virtual_variations(length):
    result = []
    for i in range(3):
        for j in range(3):
            A = np.zeros((3, 3)); A[i, j] = 1.
            result.append((f'A_{"xyz"[i]}{"xyz"[j]}', A, np.zeros(3)))
    for i in range(3):
        c = np.eye(3)[i]*length
        result.append((f'translation_{"xyz"[i]}', np.zeros((3, 3)), c))
    return result


def load_checks(geometry, X, load, surrogate, gate=1e-6, P_override=None):
    P = surrogate.response()['P'] if P_override is None else np.asarray(P_override); L = max(geometry.lengths)
    scale = np.linalg.norm(P)*L**2; mscale = scale*L
    faces = surrogate.exact_face_integrals(geometry, P)
    target = sum(r['force'] for r in faces); moment_target = sum(r['moment'] for r in faces)
    resultant = np.sum(load, axis=0); moment = np.cross(X-geometry.center, load).sum(axis=0)
    work = []
    for name, A, c in virtual_variations(L):
        variation = (X-geometry.center)@A.T+c
        actual = float(np.sum(load*variation))
        continuum = float(geometry.volume*np.sum(P*A))
        stable = mscale*max(np.linalg.norm(A), np.linalg.norm(c)/L)
        nonzero = abs(continuum) > 100*np.finfo(float).eps*stable
        denom = abs(continuum) if nonzero else stable
        work.append(dict(variation=name, discrete_J=actual, target_J=continuum,
                         error_J=abs(actual-continuum), normalized_error=abs(actual-continuum)/denom,
                         nonzero_target=bool(nonzero), normalization='abs(target)' if nonzero else 'stable stress-volume variation scale',
                         PASS=bool(abs(actual-continuum)/denom <= gate)))
    # First moments are work-equivalent paired-face loads, not a spatial face allocation.
    first_moment = load.T@(X-geometry.center)
    return dict(target_resultant_N=target, discrete_resultant_N=resultant,
                resultant_absolute_error_N=float(np.linalg.norm(resultant-target)),
                resultant_error=float(np.linalg.norm(resultant-target)/scale),
                target_moment_Nm=moment_target, discrete_moment_Nm=moment,
                moment_absolute_error_Nm=float(np.linalg.norm(moment-moment_target)),
                moment_error=float(np.linalg.norm(moment-moment_target)/mscale),
                work=work, virtual_work_error=max(r['normalized_error'] for r in work),
                force_scale_N=scale, moment_scale_Nm=mscale,
                first_moment_Pa=first_moment/geometry.volume,
                first_moment_relative_error=float(np.linalg.norm(first_moment/geometry.volume-P)/np.linalg.norm(P)),
                resultant_PASS=bool(np.linalg.norm(resultant-target)/scale <= gate),
                moment_PASS=bool(np.linalg.norm(moment-moment_target)/mscale <= gate),
                virtual_work_PASS=all(r['PASS'] for r in work))


def support_profile(geometry, family, load, h):
    distance = geometry.distance(family.X); norm = np.linalg.norm(load, axis=1)
    layer = np.rint(distance/h).astype(int); rows = []
    cumulative = 0.; total = norm.sum()
    for k in np.unique(layer):
        mask = layer == k; cumulative += norm[mask].sum()
        rows.append(dict(depth_h=int(k), depth_m=float(k*h), depth_over_delta=float(k*h/family.horizon),
                         nodes=int(mask.sum()), load_L1_N=float(norm[mask].sum()),
                         load_mean_N=float(norm[mask].mean()), load_max_N=float(norm[mask].max()),
                         cumulative_L1_fraction=float(cumulative/total)))
    quantiles = {f'd{q}': next(dict(depth_m=r['depth_m'], depth_h=r['depth_h'], depth_over_delta=r['depth_over_delta'])
                              for r in rows if r['cumulative_L1_fraction'] >= q/100) for q in [50, 90, 95, 99]}
    return rows, quantiles


def fields(geometry, family, state, surrogate):
    err = np.linalg.norm(state.Fbar-surrogate.F, axis=(1, 2))
    scale = np.linalg.norm(surrogate.F-np.eye(3)); core = geometry.distance(family.X) >= family.horizon
    plus, _ = svk_split_energy(green_lagrange(state.Ftilde), surrogate.material.lame, surrogate.material.mu)
    # Bond-quadrature energies averaged to nodes by existing energy weights.
    i, j = family.pairs.T; weights = family.wbar*family.volume[i]*family.volume[j]
    nodal_plus = np.zeros(len(family.X)); denom = np.zeros(len(family.X))
    for ids in [i, j]:
        np.add.at(nodal_plus, ids, weights*plus); np.add.at(denom, ids, weights)
    nodal_plus /= denom
    return dict(F_error_max=float(err.max()/scale),
                F_error_core_RMS=float(np.sqrt(np.mean(err[core]**2))/scale) if core.any() else None,
                F_error_boundary_RMS=float(np.sqrt(np.mean(err[~core]**2))/scale),
                Psi_plus_max_Pa=float(plus.max()), Psi_plus_P95_Pa=float(np.percentile(plus, 95)),
                Psi_total_J=float(state.energy), core_nodes=int(core.sum()),
                core_definition='Geometric d>=delta; not a fully closed multi-hop bulk stencil'), nodal_plus
