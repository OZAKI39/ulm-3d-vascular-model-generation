"""Sparse elimination; only the small retained boundary matrix is dense."""
import numpy as np
from scipy.sparse.linalg import splu
from .network_solver import conductance_matrix


def reduce_boundary(n, edges, resistance, ports, fixed_pressure):
    """Return q_ports = Y p_ports + forcing with q directed INTO the exterior.

    Nonzero source pressures give affine forcing. Ports already clamped to a
    reservoir are constraints and must be handled separately, not inverted.
    """
    ports = np.asarray(ports, dtype=int); fixed = np.asarray(sorted(fixed_pressure), dtype=int)
    if len(ports) > 128:
        raise ValueError('This API is for a small boundary, not dense full-network elimination')
    if len(set(ports)) != len(ports) or set(ports) & set(fixed) or any(i < 0 or i >= n for i in np.r_[ports, fixed]):
        raise ValueError('Distinct valid ports and disjoint fixed boundary required')
    G = conductance_matrix(n, edges, resistance)
    internal = np.setdiff1d(np.arange(n), np.r_[ports, fixed])
    Y = G[ports][:, ports].toarray()
    pf = np.asarray([fixed_pressure[int(i)] for i in fixed])
    forcing = G[ports][:, fixed] @ pf
    if len(internal):
        A = G[internal][:, internal].tocsc(); scale = abs(A).max()
        try:
            lu = splu(A/scale)
        except RuntimeError as exc:
            raise ValueError('Unreferenced external component or singular elimination') from exc
        coupling = G[internal][:, ports].tocsc()
        # Solve columns individually: never construct a dense internal NxN.
        for j in range(len(ports)):
            response = lu.solve(coupling[:, [j]].toarray().ravel()/scale)
            Y[:, j] -= G[ports][:, internal] @ response
        forcing = forcing-G[ports][:, internal] @ lu.solve((G[internal][:, fixed] @ pf)/scale)
    norm = max(float(np.linalg.norm(Y)), np.finfo(float).tiny)
    symmetry = float(np.linalg.norm(Y-Y.T)/norm)
    eig = np.linalg.eigvalsh((Y+Y.T)/2)
    passive = bool(eig.min() >= -1e-10*norm)
    cond = float(np.linalg.cond(Y))
    Z = np.linalg.inv(Y) if eig.min() > 1e-12*norm and cond < 1e12 else None
    return dict(Y=Y, Z=Z, forcing=forcing, symmetry_error=symmetry,
                min_eigenvalue=float(eig.min()), passive=passive,
                condition_number=cond if np.isfinite(cond) else None,
                resistance_form_valid=Z is not None)
