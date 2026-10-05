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
    """Project onto non-closing normal relative velocities at touching pairs."""
    distances = squareform(pdist(x))
    pairs = np.argwhere(np.triu(distances <= radii[:, None]+radii[None, :]+tolerance_m, 1))
    v = velocity.copy()
    for _ in range(max_iterations):
        largest = 0.
        for i, j in pairs:
            normal = (x[j]-x[i])/distances[i, j]
            closing = float(np.dot(v[j]-v[i], normal))
            if closing < 0:
                fraction = radii[j]/(radii[i]+radii[j])
                v[i] += fraction*closing*normal
                v[j] -= (1-fraction)*closing*normal
                largest = max(largest, -closing)
        if largest <= 1e-12:
            return v
    # Position-level no-overlap still enforces feasibility; do not silently claim velocity convergence.
    raise RuntimeError("Contact velocity projection failed to converge")
