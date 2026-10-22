"""Deterministic frictionless no-overlap projection for the demo, not coalescence."""
import numpy as np
from scipy.spatial.distance import pdist, squareform

from mevo_mb.coupling.geometry import Cylinder


def project_contacts(positions_m, radii_m, cylinder: Cylinder, tolerance_m: float, max_iterations: int):
    x, wall_hit = cylinder.project(positions_m, radii_m)
    events = set()
    min_distance = radii_m[:, None] + radii_m[None, :]
    for iteration in range(max_iterations):
        distances = squareform(pdist(x))
        overlap = min_distance - distances
        pairs = np.argwhere(np.triu(overlap > tolerance_m, 1))
        if not len(pairs):
            return x, wall_hit, sorted(events), iteration
        for i, j in pairs:
            delta = x[j]-x[i]
            d = np.linalg.norm(delta)
            if d == 0:
                raise RuntimeError("Coincident centers are invalid; reduce timestep")
            amount = max(0., radii_m[i]+radii_m[j]-d)
            # mobility-weighted constraint correction: mu_i proportional to 1/R_i
            fraction = radii_m[j] / (radii_m[i]+radii_m[j])
            x[i] -= fraction*amount*delta/d
            x[j] += (1-fraction)*amount*delta/d
            events.add((int(i), int(j)))
        x, hit = cylinder.project(x, radii_m)
        wall_hit |= hit
    raise RuntimeError("Contact projection did not converge; no unconverged step accepted")


def constrain_pair_velocities(x, radii, velocity, tolerance_m: float, max_iterations: int):
    """Mobility-metric projection onto all touching-pair constraints together.

    With A v >= 0 and R=diag(radius), v=v_free+R^-1 A.T lambda.
    The dual nonnegative QP is exactly NNLS with C=R^-1/2 A.T and
    d=-R^1/2 v_free. The positive common 6*pi*mu factor cancels.
    """
    from scipy.optimize import nnls
    distances = squareform(pdist(x))
    pairs = np.argwhere(np.triu(distances <= radii[:, None]+radii[None, :]+tolerance_m, 1))
    if not len(pairs):
        return velocity.copy()
    matrix = np.zeros((len(pairs), 3*len(x)))
    for row, (i, j) in enumerate(pairs):
        normal = (x[j]-x[i])/distances[i, j]
        matrix[row, 3*i:3*i+3] = -normal
        matrix[row, 3*j:3*j+3] = normal
    sqrt_r = np.repeat(np.sqrt(radii/np.median(radii)), 3)
    scale = max(float(np.linalg.norm(velocity, axis=1).max()), np.finfo(float).tiny)
    free = velocity.ravel()/scale
    multiplier, _ = nnls((matrix/sqrt_r).T, -sqrt_r*free, maxiter=max_iterations*max(1, len(pairs)))
    constrained = (free + matrix.T@multiplier/(sqrt_r*sqrt_r))*scale
    gap_rates = matrix@constrained
    if gap_rates.min() < -max(1e-12, 1e-10*scale):
        raise RuntimeError("Simultaneous contact projection violates normal feasibility")
    return constrained.reshape(velocity.shape)
