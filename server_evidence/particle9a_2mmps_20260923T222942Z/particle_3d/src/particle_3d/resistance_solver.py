"""Self-scaled sparse solve and simultaneous R-metric hard contact; no inverse."""
from dataclasses import dataclass
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scipy.linalg import solve_triangular
from .kinematic_contact import _dual_active_set
from .particle_shapes import EPS


class ResistanceSystemIllConditioned(RuntimeError):
    def __init__(self, record):
        self.record = dict(status='RESISTANCE_SYSTEM_ILL_CONDITIONED', **record)
        super().__init__(str(self.record))


@dataclass
class ResistanceSolution:
    velocity: np.ndarray
    unconstrained: np.ndarray
    record: dict


def contact_jacobian(system, particles, constraints):
    """P4 contact-point Jacobian lifted to full sphere 6N coordinates.

    A sphere's normal through its center gives exactly zero rotational row,
    consistent with P4's inactive sphere contact angular correction.
    """
    ids = {i: 6*k for k, i in enumerate(system.ids)}
    contacts = sorted(constraints, key=lambda c: c.canonical_id)
    rows = np.zeros((len(contacts), len(system.free)))
    for k, c in enumerate(contacts):
        if c.offset_m_s != 0:
            raise ValueError('P5 V0 requires stationary nonpenetration constraints')
        i = ids[c.particle_i_id]; rows[k, i:i+3] = c.normal
        if c.particle_j_id is not None:
            j = ids[c.particle_j_id]; rows[k, j:j+3] = -c.normal
    return rows, contacts


def solve_resistance(system, *, particles=None, constraints=(), scaled=True):
    r, b = system.matrix, system.rhs
    d = 1/np.sqrt(system.self_diagonal)
    a = (sparse.diags(d)@r@sparse.diags(d)).tocsc()
    dense = a.toarray()
    symmetry = float(np.max(np.abs(dense-dense.T)))
    if symmetry > 64*EPS*max(1., np.max(np.abs(dense))):
        raise ValueError('NONSYMMETRIC_RESISTANCE_SYSTEM')
    eig = np.linalg.eigvalsh(dense)
    condition = float(eig[-1]/eig[0]) if eig[0] > 0 else float('inf')
    numerical_bound = EPS*len(b)*condition
    if not np.isfinite(condition) or numerical_bound >= 1/512:
        raise ResistanceSystemIllConditioned(dict(condition_estimate=condition, floating_error_bound=numerical_bound,
            minimum_scaled_eigenvalue=float(eig[0]), action='REQUEST_PHYSICAL_TIME_SUBDIVISION; NO_GAP_FLOOR'))
    lu = splu(a, permc_spec='NATURAL')
    u = d*lu.solve(d*b) if scaled else splu(r.tocsc(), permc_spec='NATURAL').solve(b)
    rhs_scaled = d*b
    residual = float(np.linalg.norm(a@(u/d)-rhs_scaled)/(np.linalg.norm(a.toarray(), ord=2)*np.linalg.norm(u/d)+np.linalg.norm(rhs_scaled))) if np.any(b) else 0.
    if not np.isfinite(u).all() or residual > 512*EPS*len(b):
        raise ResistanceSystemIllConditioned(dict(condition_estimate=condition, relative_backward_residual=residual))
    hydro = u.copy()
    audit = dict(condition_estimate=condition, condition_role='SELF_DIAGONALLY_SCALED_6N_SYSTEM',
        minimum_scaled_eigenvalue=float(eig[0]), eigenvalue_roundoff_bound=float(64*EPS*len(b)*np.linalg.norm(dense, 2)),
        matrix_symmetry_error=float(np.max(np.abs((r-r.T).toarray()))),
        scaled_matrix_symmetry_error=symmetry, relative_backward_residual=residual, scaled=scaled,
        multiplier_role='KINEMATIC_CONSTRAINT_MULTIPLIER', physical_contact_force_claimed=False,
        contact_count=0, max_constraint_violation_m_s=0., resistance_correction_objective=0.)
    if constraints:
        j, contacts = contact_jacobian(system, particles, constraints)
        # B B^T = J R^-1 J^T using triangular solves, never an explicit inverse.
        chol = np.linalg.cholesky(dense)
        whitened = solve_triangular(chol, (j*d).T, lower=True).T
        scale = max(float(np.max(np.abs(hydro.reshape(-1, 6)[:, :3]))),
                    float(np.max(np.abs(system.free.reshape(-1, 6)[:, :3]))), np.finfo(float).tiny)
        lam, _, kkt = _dual_active_set(whitened, j@hydro, scale)
        u = hydro+d*lu.solve(d*(j.T@lam))
        speed = j@u
        if np.min(speed) < -kkt['velocity_budget_m_s']:
            raise ResistanceSystemIllConditioned(dict(reason='CONTACT_PRIMAL_RESIDUAL', speeds=speed.tolist()))
        delta = u-hydro
        audit.update(contact_count=len(contacts), contact_ids=[c.canonical_id for c in contacts],
            multipliers=lam.tolist(), normal_speeds_m_s=speed.tolist(), contact_kkt=kkt,
            max_constraint_violation_m_s=max(0., -float(np.min(speed))),
            resistance_correction_objective=float(.5*delta@(r@delta)),
            max_complementarity_residual=float(np.max(np.abs(lam*speed))))
    audit['dissipation'] = system.dissipation(u)
    return ResistanceSolution(u, hydro, audit)
