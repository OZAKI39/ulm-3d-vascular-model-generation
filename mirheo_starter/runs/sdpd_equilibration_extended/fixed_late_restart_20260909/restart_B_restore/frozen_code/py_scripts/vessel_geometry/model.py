"""Data contracts for geometry identity (not numerical boundary conditions)."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np


class GeometryError(ValueError):
    def __init__(self, code: str, message: str, *, status: str = "FAIL", **details):
        self.code, self.status, self.details = code, status, details
        super().__init__(f"{code}：{message}")


@dataclass
class BoundaryPatch:
    display_name: str
    entity_id: int
    port_id: str | None
    boundary_index: int | None
    kind: str
    original_role: str
    boundary_origin: str
    face_ids: list[int]
    area_m2: float
    center_m: list[float]
    outward_normal: list[float] | None
    normal_status: str
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def triangle_count(self) -> int:
        return len(self.face_ids)


@dataclass
class VesselGeometry:
    points_source: np.ndarray
    points_m: np.ndarray
    triangles: np.ndarray
    original_point_ids: np.ndarray
    original_face_ids: np.ndarray
    entity_ids: np.ndarray
    patches: list[BoundaryPatch]
    source_length_unit: str
    to_meter: float
    unit_source: dict[str, Any]
    source_cell_data: dict[str, np.ndarray] = field(default_factory=dict)
    checks: dict[str, Any] = field(default_factory=dict)

    def summary(self) -> dict[str, Any]:
        used = np.unique(self.triangles)
        return {
            "point_count": len(self.points_m), "triangle_count": len(self.triangles),
            "entity_ids": [int(x) for x in np.unique(self.entity_ids)],
            "bounds_m": [self.points_m.min(0).tolist(), self.points_m.max(0).tolist()],
            "extent_m": np.ptp(self.points_m, axis=0).tolist(),
            "referenced_bounds_m": [self.points_m[used].min(0).tolist(), self.points_m[used].max(0).tolist()],
            "referenced_extent_m": np.ptp(self.points_m[used], axis=0).tolist(),
            "unused_point_count": len(self.points_m) - len(used),
            "source_length_unit": self.source_length_unit, "output_length_unit": "m",
            "to_meter": self.to_meter,
        }


@dataclass
class ImportResult:
    run_id: str
    package_path: Path
    geometry: VesselGeometry
    manifest: dict[str, Any]
