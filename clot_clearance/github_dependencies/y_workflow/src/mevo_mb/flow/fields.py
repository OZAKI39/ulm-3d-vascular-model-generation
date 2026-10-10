"""Velocity providers with the same SI interface for analytic and VTK fields."""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path

import numpy as np


class FlowField(ABC):
    @abstractmethod
    def velocity(self, positions_m: np.ndarray, time_s: float) -> np.ndarray:
        """Return N by 3 velocity in m/s at physical time_s."""


class ZeroFlowField(FlowField):
    def velocity(self, positions_m, time_s):
        return np.zeros_like(positions_m, dtype=float)


@dataclass
class PoiseuilleCylinderFlowField(FlowField):
    radius_m: float
    maximum_velocity_m_s: float

    def velocity(self, positions_m, time_s):
        r2 = np.sum(np.asarray(positions_m)[:, 1:] ** 2, axis=1)
        v = np.zeros_like(positions_m, dtype=float)
        v[:, 0] = self.maximum_velocity_m_s * np.clip(1 - r2 / self.radius_m**2, 0, 1)
        return v


class SimVascularFlowField(FlowField):
    """Static VTK point velocity interpolation; no nearest-neighbor extrapolation.

    Explicit scale parameters convert file coordinates and velocity into SI.
    Outside-domain probes raise, preserving invalid classifications.
    """

    def __init__(
        self,
        path: str | Path,
        velocity_array: str = "Velocity",
        coordinate_scale_to_m: float = 1.0,
        velocity_scale_to_m_s: float = 1.0,
        spatial_method: str = "vtk",
    ):
        import pyvista as pv

        self.grid = pv.read(path)
        self.grid.points *= coordinate_scale_to_m
        if velocity_array not in self.grid.point_data:
            raise ValueError("A point-data velocity array is required")
        self.array = velocity_array
        self.scale = velocity_scale_to_m_s
        if spatial_method not in ("vtk", "tetra_p1"):
            raise ValueError("Spatial method must be vtk or tetra_p1")
        self.sampler = None
        if spatial_method == "tetra_p1":
            from mevo_mb.flow.tetra import TetraP1Sampler
            self.sampler = TetraP1Sampler(self.grid,self.grid.point_data[self.array])

    def velocity(self, positions_m, time_s):
        import pyvista as pv

        if self.sampler is not None:
            return self.sampler.sample(positions_m)*self.scale
        sampled = pv.PolyData(np.asarray(positions_m)).sample(self.grid)
        if not np.all(sampled["vtkValidPointMask"]):
            raise ValueError("Flow query outside VTK mesh; no extrapolation performed")
        return np.asarray(sampled[self.array], float) * self.scale


class VTKFlowField(FlowField):
    """Generic static or time-indexed VTU, linear interpolation, no extrapolation.

    Existing SimVascularFlowField remains unchanged for baseline compatibility.
    Snapshots may use different meshes; each spatial probe must be inside both.
    """

    def __init__(self, paths, times_s=None, **kwargs):
        if isinstance(paths, (str, Path)):
            paths = [paths]
        self.fields = [SimVascularFlowField(p, **kwargs) for p in paths]
        self.times = np.array([0.0] if times_s is None else times_s, dtype=float)
        if (
            len(self.fields) != len(self.times)
            or not len(self.times)
            or not np.isfinite(self.times).all()
            or np.any(np.diff(self.times) <= 0)
        ):
            raise ValueError("One strictly increasing finite time per VTU required")

    def velocity(self, positions_m, time_s):
        if len(self.fields) == 1:
            return self.fields[0].velocity(positions_m, time_s)
        if not np.isfinite(time_s) or not self.times[0] <= time_s <= self.times[-1]:
            raise ValueError("Flow time outside available snapshots; no extrapolation")
        j = min(int(np.searchsorted(self.times, time_s, side="right")) - 1, len(self.times) - 2)
        w = (time_s - self.times[j]) / (self.times[j + 1] - self.times[j])
        return (1 - w) * self.fields[j].velocity(positions_m, time_s) + w * self.fields[j + 1].velocity(
            positions_m, time_s
        )
