"""Independent, read-only convex quadratic stationary-state certificate.

FEM velocity direction is not the force driving the affine wall model. A
geometrically feasible FEM-downstream direction need not lower its quadratic
objective. Never change the production matrix, solver, or recorded trajectory.
"""
from itertools import combinations
import numpy as np


def independent_contact_equilibrium(matrix, rhs, jacobian, recorded_velocity):
    R, b, J, U = map(lambda x: np.asarray(x, float), (matrix, rhs, jacobian, recorded_velocity))
    scale = np.sqrt(np.diag(R))
    A = R/scale[:, None]/scale[None, :]
    beta = b/scale; force_scale = np.linalg.norm(beta)
    if force_scale == 0:
        raise ValueError('Zero-load audit requires a separate certificate')
    beta /= force_scale
    H = J/scale[None, :]
    row_scale = np.linalg.norm(H, axis=1)
    if np.any(row_scale == 0):
        raise ValueError('Zero constraint row')
    H /= row_scale[:, None]
    spectrum = np.linalg.eigvalsh(A)
    if spectrum.min() <= 0:
        raise ValueError('Original scaled matrix is not positive definite')
    candidates = []
    for n in range(min(len(H), len(b))+1):
        for active in combinations(range(len(H)), n):
            h = H[list(active)]
            if n:
                singular = np.linalg.svd(h, compute_uv=False)
                if singular[-1] <= np.sqrt(max(h.shape)*np.finfo(float).eps)*singular[0]:
                    continue  # Independent active sets only; all rows checked below.
            K = np.block([[A, -h.T], [h, np.zeros((n,n))]])
            try:
                sol = np.linalg.solve(K, np.r_[beta, np.zeros(n)])
            except np.linalg.LinAlgError:
                continue
            z, lam = sol[:len(b)], np.zeros(len(H))
            lam[list(active)] = sol[len(b):]
            tolerance = 512*np.finfo(float).eps*max(K.shape)*np.linalg.cond(K)*max(1., np.linalg.norm(sol,np.inf))
            primal = H@z; residual = A@z-beta-H.T@lam
            if primal.min() < -tolerance or lam.min() < -tolerance or np.linalg.norm(residual,np.inf)>tolerance:
                continue
            candidates.append(dict(active=active,z=z,lam=lam,tolerance=tolerance,
                objective=float(.5*z@A@z-beta@z),minimum_primal=float(primal.min()),
                minimum_multiplier=float(lam.min()),residual=float(np.linalg.norm(residual,np.inf)),
                complementarity=float(np.max(np.abs(lam*primal)))))
    if not candidates:
        raise ValueError('No independently certified KKT active set')
    best = min(candidates, key=lambda c:c['objective'])
    z0 = U*scale/force_scale
    error = float(np.linalg.norm(z0-best['z'],np.inf))
    velocity = best['z']*force_scale/scale
    budget = best['tolerance']*force_scale/scale
    return dict(method='INDEPENDENT_EXHAUSTIVE_ACTIVE_SET_ORIGINAL_6D_SPD_QUADRATIC',
        original_matrix_unchanged=True,independent_velocity=velocity.tolist(),
        velocity_roundoff_budget=budget.tolist(),minimum_scaled_eigenvalue=float(spectrum.min()),
        active_constraints=list(best['active']),feasible_active_sets=len(candidates),
        KKT_residual_scaled=best['residual'],minimum_primal_scaled=best['minimum_primal'],
        minimum_multiplier_scaled=best['minimum_multiplier'],complementarity_scaled=best['complementarity'],
        recorded_velocity_scaled_error=error,scaled_roundoff_budget=best['tolerance'],
        recorded_solution_verified=bool(error<=best['tolerance']),
        independent_translation_zero_within_roundoff=bool(np.all(np.abs(velocity[:3])<=budget[:3])),
        conclusion='CURRENT_MODEL_MECHANICAL_EQUILIBRIUM; NOT_GEOMETRIC_NO_PATH_OR_PHYSIOLOGICAL_TRAPPING')
