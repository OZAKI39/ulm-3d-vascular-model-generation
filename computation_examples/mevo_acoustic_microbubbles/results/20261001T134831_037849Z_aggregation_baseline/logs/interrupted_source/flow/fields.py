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
        r2 = np.sum(np.asarray(positions_m)[:, 1:]**2, axis=1)
        v = np.zeros_like(positions_m, dtype=float)
        v[:, 0] = self.maximum_velocity_m_s * np.clip(1-r2/self.radius_m**2, 0, 1)
        return v


class SimVascularFlowField(FlowField):
    """Static VTK point velocity interpolation; no nearest-neighbor extrapolation.

    Explicit scale parameters convert file coordinates and velocity into SI.
    Outside-domain probes raise, preserving invalid classifications.
    """
    def __init__(self, path: str | Path, velocity_array: str = "Velocity", coordinate_scale_to_m: float = 1., velocity_scale_to_m_s: float = 1.):
        import pyvista as pv
        self.grid = pv.read(path)
        self.grid.points *= coordinate_scale_to_m
        if velocity_array not in self.grid.point_data:
            raise ValueError("A point-data velocity array is required")
        self.array = velocity_array
        self.scale = velocity_scale_to_m_s

    def velocity(self, positions_m, time_s):
        import pyvista as pv
        sampled = pv.PolyData(np.asarray(positions_m)).sample(self.grid)
        if not np.all(sampled["vtkValidPointMask"]):
            raise ValueError("Flow query outside VTK mesh; no extrapolation performed")
        return np.asarray(sampled[self.array], float)*self.scale
