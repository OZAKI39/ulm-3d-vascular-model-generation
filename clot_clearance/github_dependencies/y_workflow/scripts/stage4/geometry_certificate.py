"""Optional exact no-op certificate for a Stage 4 wall projection.

The experiment runner enables this only through an explicit option. A known interior point and
the 1-Lipschitz Euclidean distance bound can prove that a new bubble is still
strictly clear of the same closed surface. No approximate distances are returned
or inserted into the exact query cache. All uncertain cases use the old method.
"""
import numpy as np


def install(geometry, *, per_particle=False):
    original = geometry.project
    counts = dict(calls=0, certified_inactive=0, partial_batches=0,
                  certified_particles=0, minimum_certified_margin_m=None)
    anchor = None
    anchor_mesh = geometry._mesh_distance

    def project(positions, radii):
        nonlocal anchor, anchor_mesh
        counts["calls"] += 1
        x = np.asarray(positions)
        if geometry._mesh_distance is not anchor_mesh:
            anchor, anchor_mesh = None, geometry._mesh_distance
        previous = geometry._distance_cache
        if previous is not None and previous[0].shape == x.shape:
            anchor = previous
        elif per_particle:
            previous = anchor
        if len(x) and geometry._mesh_distance is not None and previous is not None and previous[0].shape == x.shape:
            margin = -previous[1]-np.linalg.norm(x-previous[0], axis=1)-np.asarray(radii)
            unit_gradient = np.abs(np.linalg.norm(previous[2], axis=1)-1.) < 1e-10
            clear = (previous[1] < 0) & unit_gradient & (margin > 1e-10)
            if clear.all() or (per_particle and clear.any()):
                counts["certified_particles"] += int(clear.sum())
                smallest = float(margin[clear].min())
                prior = counts["minimum_certified_margin_m"]
                counts["minimum_certified_margin_m"] = smallest if prior is None else min(prior, smallest)
                result, hits = x.copy(), np.zeros(len(x), bool)
                if clear.all():
                    counts["certified_inactive"] += 1
                else:
                    counts["partial_batches"] += 1
                    # Independent pointwise wall corrections retain the original
                    # exact method for every row not proven to be a no-op.
                    r = np.broadcast_to(radii, (len(x),))
                    result[~clear], hits[~clear] = original(x[~clear], r[~clear])
                return result, hits
        return original(positions, radii)

    geometry.project = project
    return counts
