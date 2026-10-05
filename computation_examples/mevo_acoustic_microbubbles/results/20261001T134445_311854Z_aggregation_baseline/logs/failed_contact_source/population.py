"""Demo size distributions and deterministic non-overlapping initial positions."""
from dataclasses import dataclass

import numpy as np

from mevo_mb.coupling.geometry import Cylinder


@dataclass
class BubbleSizeDistribution:
    config: dict

    def sample(self, count: int, rng: np.random.Generator) -> np.ndarray:
        c = self.config
        kind = c["kind"]
        if kind == "monodisperse":
            out = np.full(count, c["radius_m"])
        elif kind == "samples":
            values = np.asarray(c["samples_m"], float)
            values = values[(values >= c["minimum_radius_m"]) & (values <= c["maximum_radius_m"])]
            if len(values) == 0:
                raise ValueError("No provided radii inside the truncation range")
            out = rng.choice(values, size=count)
        elif kind == "lognormal":
            from scipy.special import ndtr, ndtri
            sigma = c["lognormal_sigma"]
            if sigma <= 0:
                raise ValueError("Positive lognormal sigma required")
            lo, hi = np.log(np.array([c["minimum_radius_m"], c["maximum_radius_m"]])/c["radius_m"])/sigma
            quantile = rng.uniform(ndtr(lo), ndtr(hi), count)
            out = c["radius_m"] * np.exp(sigma*ndtri(quantile))
        else:
            raise ValueError(f"Unknown size distribution {kind}")
        if np.any(out < c["minimum_radius_m"]) or np.any(out > c["maximum_radius_m"]):
            raise ValueError("Radii outside declared truncation range")
        return out


@dataclass
class BubblePopulation:
    positions_m: np.ndarray
    radii_m: np.ndarray
    ids: np.ndarray

    @classmethod
    def initialize(cls, config: dict):
        c = config["particles"]
        rng = np.random.default_rng(config["seed"])
        radii = BubbleSizeDistribution(c["size_distribution"]).sample(c["count"], rng)
        cylinder = Cylinder(**config["vessel"])
        init = c["initialization"]
        positions = []
        for i, radius in enumerate(radii):
            for _ in range(100000):
                if init["kind"] == "cloud":
                    direction = rng.normal(size=3)
                    x = np.asarray(init["cloud_center_m"]) + direction/np.linalg.norm(direction)*init["cloud_radius_m"]*rng.random()**(1/3)
                elif init["kind"] == "cylinder":
                    angle = rng.uniform(0, 2*np.pi)
                    r = np.sqrt(rng.random())*init["radial_fraction"]*(cylinder.radius_m-radius)
                    x = np.array([rng.uniform(-init["x_half_width_m"], init["x_half_width_m"]), r*np.cos(angle), r*np.sin(angle)])
                else:
                    raise ValueError("Unknown spatial initialization")
                if cylinder.contains(x[None], radius)[0] and (not positions or np.all(np.linalg.norm(np.asarray(positions)-x, axis=1) >= np.maximum(radius+radii[:i], c["minimum_separation_m"]))):
                    positions.append(x)
                    break
            else:
                raise RuntimeError("Unable to initialize cloud without overlaps; reduce number or increase volume")
        return cls(np.asarray(positions), radii, np.arange(c["count"]))
