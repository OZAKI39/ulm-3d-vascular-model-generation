"""Small typed records used by Ultraliser input and quality control."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


class GeometryValidationError(RuntimeError):
    """Raised when a saved ROI cannot be serialized or validated safely."""


@dataclass(slots=True)
class BranchGeometry:
    branch_id: int
    local_node_ids: tuple[int, ...]
    source_global_nodes: tuple[int, ...]
    source_global_edges: tuple[int, ...]
    raw_points_um: np.ndarray
    raw_radius_um: np.ndarray
    points_um: np.ndarray = field(default_factory=lambda: np.empty((0, 3), dtype=float))
    radius_um: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))
    arc_length_um: np.ndarray = field(default_factory=lambda: np.empty(0, dtype=float))

    @property
    def length_um(self) -> float:
        points = self.points_um if len(self.points_um) else self.raw_points_um
        return float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())


@dataclass(frozen=True, slots=True)
class RadiusFidelitySample:
    branch_id: int
    sample_index: int
    arc_length_um: float
    center_um: tuple[float, float, float]
    tangent: tuple[float, float, float]
    source_radius_um: float
    reconstructed_radius_um: float
    relative_error: float
    section_xy_um: tuple[tuple[float, float], ...] = ()

    def report(self) -> dict[str, Any]:
        return {
            "branch_id": self.branch_id,
            "sample_index": self.sample_index,
            "arc_length_um": self.arc_length_um,
            "center_x_um": self.center_um[0],
            "center_y_um": self.center_um[1],
            "center_z_um": self.center_um[2],
            "tangent_x": self.tangent[0],
            "tangent_y": self.tangent[1],
            "tangent_z": self.tangent[2],
            "source_radius_um": self.source_radius_um,
            "reconstructed_radius_um": self.reconstructed_radius_um,
            "relative_error": self.relative_error,
            "absolute_relative_error": abs(self.relative_error),
        }
