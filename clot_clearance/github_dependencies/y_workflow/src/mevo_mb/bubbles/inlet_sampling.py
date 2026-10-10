"""Paired finite-bubble cross-section initialization with explicit area/flux measures."""
import numpy as np
from mevo_mb.bubbles.population import BubblePopulation


def initialize_inlet(geometry, flow, config):
    count, seed = int(config["count"]), int(config["seed"])
    rng = np.random.default_rng(seed)
    radius = float(config["radius_m"])
    limit = (geometry.parent_radius-radius)*config["cross_section_radius_fraction"]
    axial = config["initialization_x_m"]
    policy = config["sampling_policy"]
    minimum = max(2*radius, config["minimum_separation_m"])
    if policy == "uniform_area":
        i = np.arange(count)
        r = limit*np.sqrt((i+.5)/count)
        theta = i*np.pi*(3-np.sqrt(5)) + rng.uniform(0, 2*np.pi)
        positions = np.column_stack([np.full(count, axial), r*np.cos(theta), r*np.sin(theta)])
    elif policy == "flux_weighted":
        if hasattr(flow, "maximum_velocity"):
            bound = flow.maximum_velocity*1.01
        elif hasattr(flow, "fields"):
            source = flow.fields[0]
            bound = float(np.max(np.asarray(source.grid.point_data[source.array])[:, 0]))*source.scale*1.01
        else:
            raise ValueError("Flux sampling needs a certified upper velocity bound")
        if bound <= 0:
            raise ValueError("No positive inlet flux")
        accepted = []
        for _ in range(10000):
            n = max(128, 4*(count-len(accepted)))
            r = limit*np.sqrt(rng.random(n))
            theta = rng.uniform(0, 2*np.pi, n)
            proposals = np.column_stack([np.full(n, axial), r*np.cos(theta), r*np.sin(theta)])
            speed = np.maximum(flow.velocity(proposals, 0)[:, 0], 0)
            if np.any(speed > bound):
                raise ValueError("Certified flux rejection-sampling bound exceeded")
            selected = proposals[rng.random(n)*bound < speed]
            for p in selected:
                if not accepted or np.min(np.linalg.norm(np.asarray(accepted)-p, axis=1)) >= minimum:
                    accepted.append(p)
                if len(accepted) == count:
                    break
            if len(accepted) == count:
                break
        if len(accepted) != count:
            raise RuntimeError("Insufficient non-overlapping flux samples")
        positions = np.asarray(accepted)
    else:
        raise ValueError("Sampling must be uniform_area or flux_weighted")
    if not geometry.contains(positions, radius).all():
        raise ValueError("Initial positions outside finite-radius lumen")
    return BubblePopulation(positions, np.full(count, radius), np.arange(count))
