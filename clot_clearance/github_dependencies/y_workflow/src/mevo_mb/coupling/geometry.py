"""Cartesian acoustic coordinates and finite-radius cylindrical lumen boundaries."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def grid_axes(shape: list[int], spacing_m: list[float]) -> tuple[np.ndarray, ...]:
    """Match kWaveGrid: even grids have a sample exactly at zero."""
    return tuple((np.arange(n) - n // 2) * d for n, d in zip(shape, spacing_m, strict=True))


@dataclass(frozen=True)
class Cylinder:
    radius_m: float
    length_m: float

    def wall_distance(self, positions_m):
        return np.minimum(self.radius_m-np.linalg.norm(positions_m[:, 1:], axis=1),
                          self.length_m/2-np.abs(positions_m[:, 0]))

    def contains(self, positions_m: np.ndarray, radii_m: np.ndarray | float = 0.0) -> np.ndarray:
        x = np.asarray(positions_m)
        return (np.linalg.norm(x[..., 1:], axis=-1) <= self.radius_m - radii_m) & (np.abs(x[..., 0]) <= self.length_m / 2 - radii_m)

    def mask(self, axes: tuple[np.ndarray, ...]) -> np.ndarray:
        x, y, z = np.meshgrid(*axes, indexing="ij", sparse=True)
        return (np.abs(x) <= self.length_m / 2) & (y*y + z*z <= self.radius_m**2)

    def project(self, positions_m: np.ndarray, radii_m: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Demo hard projection, including end caps; no wall hydrodynamic claim."""
        x = positions_m.copy()
        radial = np.linalg.norm(x[:, 1:], axis=1)
        limit = self.radius_m - radii_m
        hit = radial > limit
        x[hit, 1:] *= (limit[hit] / radial[hit])[:, None]
        cap = np.abs(x[:, 0]) > self.length_m / 2 - radii_m
        x[:, 0] = np.clip(x[:, 0], -self.length_m / 2 + radii_m, self.length_m / 2 - radii_m)
        return x, hit | cap

    def constrain_velocity(self, positions_m: np.ndarray, radii_m: np.ndarray, velocity_m_s: np.ndarray) -> np.ndarray:
        v = velocity_m_s.copy()
        radial = np.linalg.norm(positions_m[:, 1:], axis=1)
        normal = positions_m[:, 1:] / np.maximum(radial[:, None], np.finfo(float).tiny)
        outward = np.sum(v[:, 1:] * normal, axis=1)
        active = (radial >= self.radius_m - radii_m - 1e-12) & (outward > 0)
        v[active, 1:] -= outward[active, None] * normal[active]
        cap = np.abs(positions_m[:, 0]) >= self.length_m / 2 - radii_m - 1e-12
        v[cap & (v[:, 0] * positions_m[:, 0] > 0), 0] = 0
        return v


def geometry_from_config(config):
    if config.get("geometry", {}).get("kind") == "bifurcation":
        from mevo_mb.coupling.bifurcation import BifurcationGeometry
        geometry = BifurcationGeometry({k:v for k,v in config["geometry"].items()
                                       if k not in ("kind", "mesh_surface")})
        if config["geometry"].get("mesh_surface"):
            import pyvista as pv
            geometry.bind_mesh_surface(pv.read(config["geometry"]["mesh_surface"]))
        return geometry
    return Cylinder(**config["vessel"])
