"""Sparse conductance assembly and grounded Dirichlet solve, including cycles."""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import splu, onenormest, LinearOperator


def conductance_matrix(n, edges, resistance):
    edges, resistance = np.asarray(edges, dtype=int), np.asarray(resistance, dtype=float)
    if (edges.shape != (len(resistance), 2) or np.any(edges < 0) or np.any(edges >= n)
            or np.any(edges[:, 0] == edges[:, 1]) or not np.all(np.isfinite(resistance)) or np.any(resistance <= 0)):
        raise ValueError('Invalid edges or positive finite resistance required')
    u, v = edges.T; g = 1/resistance
    return coo_matrix((np.r_[g, g, -g, -g], (np.r_[u, v, u, v], np.r_[u, v, v, u])), shape=(n, n)).tocsr()


@dataclass
class Solution:
    pressure: np.ndarray
    flow: np.ndarray
    node_outflow: np.ndarray
    audit: dict


def solve_network(n, edges, resistance, fixed_pressure):
    G = conductance_matrix(n, edges, resistance)
    b = np.array(sorted(fixed_pressure), dtype=int)
    if len(b) == 0 or np.any(b < 0) or np.any(b >= n):
        raise ValueError('Valid pressure reference required')
    count, comp = connected_components(G, directed=False)
    if len(set(comp[b])) != count:
        raise ValueError('Unreferenced component: pressure gauge is undetermined')
    p = np.zeros(n); p[b] = [fixed_pressure[int(i)] for i in b]
    if not np.all(np.isfinite(p)):
        raise ValueError('Nonfinite pressure')
    free = np.setdiff1d(np.arange(n), b)
    cond = 1.
    if len(free):
        A = G[free][:, free].tocsc()
        scale = abs(A).max()
        A = A/scale
        rhs = -(G[free][:, b] @ p[b])/scale
        lu = splu(A)
        p[free] = lu.solve(rhs)
        inverse = LinearOperator(A.shape, matvec=lu.solve, rmatvec=lambda x: lu.solve(x, trans='T'), matmat=lu.solve)
        cond = float(onenormest(A)*onenormest(inverse))
    e = np.asarray(edges); q = (p[e[:, 0]]-p[e[:, 1]])/resistance
    out = np.bincount(e[:, 0], weights=q, minlength=n)-np.bincount(e[:, 1], weights=q, minlength=n)
    scale_q = max(float(np.max(abs(q))), np.finfo(float).tiny)
    residual = float(np.max(abs(out[free]))/scale_q) if len(free) else 0.
    boundary_mass = float(abs(out[b].sum())/scale_q)
    if not np.all(np.isfinite(p)) or max(residual, boundary_mass) > 1e-9:
        raise ValueError('Network mass balance failed')
    return Solution(p, q, out, dict(matrix_size=n, nnz=G.nnz, free_nodes=len(free),
                    solver='scipy.sparse.linalg.splu; scaled sparse system',
                    free_matrix_condition_1_estimate=cond, max_internal_relative_residual=residual,
                    boundary_mass_relative_residual=boundary_mass))


def operating_point_scale(unit_roi_inlet_q, target_q):
    if not np.isfinite(unit_roi_inlet_q) or unit_roi_inlet_q <= 0 or not np.isfinite(target_q) or target_q <= 0:
        raise ValueError('ROI inlet flow must be finite positive in the declared orientation; never take abs')
    return target_q/unit_roi_inlet_q
