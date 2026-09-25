"""Immutable surface partitioning and Cartesian-baseline geometry helpers."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh

from .io import FlowError, FlowInputs, sha256_file


BOUNDARY_LABELS = ("wall", "inlet", "outlet_01", "outlet_02", "outlet_03")


@dataclass(frozen=True, slots=True)
class BoundaryPatch:
    label: str
    entity_id: int
    face_indices: np.ndarray
    triangle_count: int
    area_um2: float
    center_um: np.ndarray
    outward_normal: np.ndarray
    equivalent_radius_um: float
    path_m: Path
    sha256: str
    port_id: str | None


@dataclass(frozen=True, slots=True)
class SurfacePartition:
    points_um: np.ndarray
    faces: np.ndarray
    entity_ids: np.ndarray
    mesh_um: trimesh.Trimesh
    patches: tuple[BoundaryPatch, ...]
    qc: dict[str, Any]

    def patch(self, label: str) -> BoundaryPatch:
        return next(item for item in self.patches if item.label == label)


@dataclass(frozen=True, slots=True)
class BoundingCube:
    origin_m: np.ndarray
    side_m: float
    level: int
    cells_per_axis: int
    margin_cells_minimum: float


@dataclass(frozen=True, slots=True)
class SeedPoint:
    coordinates_m: np.ndarray
    coordinates_um: np.ndarray
    candidate_offset_radius: float
    seed_inside_lumen: bool
    method: str


def _triangles(points: np.ndarray, faces: np.ndarray) -> np.ndarray:
    return np.asarray(points, dtype=float)[np.asarray(faces, dtype=np.int64)]


def _triangle_area_centroid_normal(
    points: np.ndarray, faces: np.ndarray, *, mesh_signed_volume: float
) -> tuple[float, np.ndarray, np.ndarray]:
    triangles = _triangles(points, faces)
    cross = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    double_area = np.linalg.norm(cross, axis=1)
    area = 0.5 * double_area
    total = float(area.sum())
    if total <= 0.0:
        raise FlowError("CFD_FLOW_BOUNDARY_PARTITION_FAILED", "Zero-area boundary patch")
    centers = triangles.mean(axis=1)
    centroid = np.average(centers, axis=0, weights=area)
    normal = cross.sum(axis=0)
    norm = float(np.linalg.norm(normal))
    if norm <= 0.0:
        raise FlowError("CFD_FLOW_BOUNDARY_PARTITION_FAILED", "Boundary patch has no stable normal")
    normal = normal / norm
    if mesh_signed_volume < 0.0:
        normal = -normal
    return total, centroid, normal


def _manifest_mapping(inputs: FlowInputs, entity_values: set[int]) -> tuple[dict[str, int], dict[str, str]]:
    inlet_rows = [row for row in inputs.boundary_manifest if row.get("role") == "ASSUMED_INLET"]
    outlet_rows = [row for row in inputs.boundary_manifest if row.get("role") == "ASSUMED_OUTLET"]
    if len(inlet_rows) != 1 or len(outlet_rows) != 3:
        raise FlowError("CFD_FLOW_BOUNDARY_PARTITION_FAILED", "Manifest must contain 1 inlet and 3 outlets")
    mapping = {"inlet": int(inlet_rows[0]["vmtk_cap_entity_id"])}
    ports = {"inlet": inlet_rows[0]["port_id"]}
    for index, row in enumerate(outlet_rows, start=1):
        label = f"outlet_{index:02d}"
        mapping[label] = int(row["vmtk_cap_entity_id"])
        ports[label] = row["port_id"]
    cap_ids = set(mapping.values())
    wall_ids = entity_values - cap_ids
    if len(wall_ids) != 1 or len(cap_ids) != 4 or cap_ids | wall_ids != entity_values:
        raise FlowError(
            "CFD_FLOW_BOUNDARY_PARTITION_FAILED",
            f"Could not infer one wall plus four unique caps from entities {sorted(entity_values)}",
        )
    mapping["wall"] = next(iter(wall_ids))
    ports["wall"] = ""
    return mapping, ports


def partition_surface(inputs: FlowInputs, output_dir: Path) -> SurfacePartition:
    """Split CellEntityIds into five exact, unshifted meter STL patches."""

    surface = pv.read(inputs.tagged_surface_vtp).triangulate()
    if "CellEntityIds" not in surface.cell_data:
        raise FlowError("CFD_FLOW_BOUNDARY_PARTITION_FAILED", "CellEntityIds is missing")
    faces = np.asarray(surface.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    points = np.asarray(surface.points, dtype=float)
    entity_ids = np.asarray(surface.cell_data["CellEntityIds"], dtype=np.int64).reshape(-1)
    if len(faces) == 0 or len(entity_ids) != len(faces):
        raise FlowError("CFD_FLOW_BOUNDARY_PARTITION_FAILED", "Invalid tagged triangle surface")
    mesh = trimesh.Trimesh(points.copy(), faces.copy(), process=False)
    if not mesh.is_watertight:
        raise FlowError("CFD_FLOW_BOUNDARY_PARTITION_FAILED", "Frozen surface is not watertight")
    mapping, ports = _manifest_mapping(inputs, set(int(value) for value in np.unique(entity_ids)))
    output_dir.mkdir(parents=True, exist_ok=True)
    assignments = np.zeros(len(faces), dtype=np.int8)
    patches: list[BoundaryPatch] = []
    for label in BOUNDARY_LABELS:
        entity = mapping[label]
        face_indices = np.flatnonzero(entity_ids == entity)
        assignments[face_indices] += 1
        area, center, normal = _triangle_area_centroid_normal(
            points, faces[face_indices], mesh_signed_volume=float(mesh.volume)
        )
        patch_mesh = trimesh.Trimesh(
            vertices=(points * 1.0e-6).copy(), faces=faces[face_indices].copy(), process=False
        )
        patch_mesh.remove_unreferenced_vertices()
        patch_path = output_dir / f"{label}.stl"
        patch_mesh.export(patch_path, file_type="stl")
        patches.append(
            BoundaryPatch(
                label=label,
                entity_id=entity,
                face_indices=face_indices,
                triangle_count=len(face_indices),
                area_um2=area,
                center_um=center,
                outward_normal=normal,
                equivalent_radius_um=math.sqrt(area / math.pi),
                path_m=patch_path,
                sha256=sha256_file(patch_path),
                port_id=ports[label] or None,
            )
        )
    missing = int(np.count_nonzero(assignments == 0))
    duplicate = int(np.count_nonzero(assignments > 1))
    expected_counts = {
        label: int(next(row["triangle_count"] for row in inputs.boundary_manifest if row["port_id"] == ports[label]))
        for label in BOUNDARY_LABELS[1:]
    }
    actual_counts = {patch.label: patch.triangle_count for patch in patches}
    count_match = all(actual_counts[label] == expected_counts[label] for label in expected_counts)
    qc = {
        "status": "PASS" if missing == 0 and duplicate == 0 and count_match else "FAIL",
        "source_triangle_count": int(len(faces)),
        "patch_triangle_count_sum": int(sum(item.triangle_count for item in patches)),
        "missing_triangles": missing,
        "duplicate_triangles": duplicate,
        "overlap_triangles": duplicate,
        "unexpected_triangles": 0,
        "manifest_cap_counts_match": count_match,
        "coordinate_transform": "meter = micrometer * 1e-6",
        "translation_applied": False,
        "surface_geometry_modified": False,
        "patches": {
            item.label: {
                "entity_id": item.entity_id,
                "triangle_count": item.triangle_count,
                "area_um2": item.area_um2,
                "path_m": str(item.path_m),
                "sha256": item.sha256,
            }
            for item in patches
        },
    }
    if qc["status"] != "PASS":
        raise FlowError("CFD_FLOW_BOUNDARY_PARTITION_FAILED", f"Partition QC failed: {qc}")
    return SurfacePartition(points, faces, entity_ids, mesh, tuple(patches), qc)


def load_frozen_surface_partition(inputs: FlowInputs, patch_dir: Path) -> SurfacePartition:
    """Load existing solver patches and rebuild only in-memory port metadata."""

    surface = pv.read(inputs.tagged_surface_vtp).triangulate()
    if "CellEntityIds" not in surface.cell_data:
        raise FlowError("CFD_FLOW_FROZEN_SEEDER_MESH_INVALID", "CellEntityIds is missing")
    faces = np.asarray(surface.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    points = np.asarray(surface.points, dtype=float)
    entity_ids = np.asarray(surface.cell_data["CellEntityIds"], dtype=np.int64).reshape(-1)
    mesh = trimesh.Trimesh(points.copy(), faces.copy(), process=False)
    mapping, ports = _manifest_mapping(inputs, set(int(value) for value in np.unique(entity_ids)))
    patches: list[BoundaryPatch] = []
    assignments = np.zeros(len(faces), dtype=np.int8)
    for label in BOUNDARY_LABELS:
        path = Path(patch_dir) / f"{label}.stl"
        if not path.is_file():
            raise FlowError("CFD_FLOW_FROZEN_SEEDER_MESH_INVALID", f"Missing frozen patch: {path}")
        entity = mapping[label]
        face_indices = np.flatnonzero(entity_ids == entity)
        assignments[face_indices] += 1
        area, center, normal = _triangle_area_centroid_normal(
            points,
            faces[face_indices],
            mesh_signed_volume=float(mesh.volume),
        )
        patches.append(
            BoundaryPatch(
                label=label,
                entity_id=entity,
                face_indices=face_indices,
                triangle_count=len(face_indices),
                area_um2=area,
                center_um=center,
                outward_normal=normal,
                equivalent_radius_um=math.sqrt(area / math.pi),
                path_m=path,
                sha256=sha256_file(path),
                port_id=ports[label] or None,
            )
        )
    if np.any(assignments != 1):
        raise FlowError("CFD_FLOW_FROZEN_SEEDER_MESH_INVALID", "Frozen partition is not an exact union")
    qc = {
        "status": "PASS",
        "mode": "READ_ONLY_FROZEN_PATCH_REFERENCE",
        "source_triangle_count": int(len(faces)),
        "patch_triangle_count_sum": int(sum(item.triangle_count for item in patches)),
        "surface_geometry_modified": False,
    }
    return SurfacePartition(points, faces, entity_ids, mesh, tuple(patches), qc)


def robust_inside(mesh_um: trimesh.Trimesh, points_um: np.ndarray) -> np.ndarray:
    """Return deterministic closed-mesh containment, with signed-distance fallback."""

    points = np.atleast_2d(np.asarray(points_um, dtype=float))
    try:
        result = np.asarray(mesh_um.contains(points), dtype=bool)
    except (ModuleNotFoundError, ValueError):
        result = np.asarray(trimesh.proximity.signed_distance(mesh_um, points) > 0.0, dtype=bool)
    return result


def find_seed_point(partition: SurfacePartition) -> SeedPoint:
    """Try only 0.5R, 1R, and 2R inward from the inlet area centroid."""

    inlet = partition.patch("inlet")
    factors = (0.5, 1.0, 2.0)
    candidates = np.asarray(
        [inlet.center_um - factor * inlet.equivalent_radius_um * inlet.outward_normal for factor in factors]
    )
    inside = robust_inside(partition.mesh_um, candidates)
    for factor, point, accepted in zip(factors, candidates, inside, strict=True):
        if bool(accepted):
            return SeedPoint(point * 1.0e-6, point, factor, True, "closed_mesh_contains")
    raise FlowError(
        "CFD_FLOW_SEED_POINT_FAILED",
        "The deterministic 0.5R, 1R, and 2R inward candidates are outside the closed lumen",
    )


def compute_bounding_cube(bounds_um: tuple[float, ...], dx_m: float, margin_cells: int) -> BoundingCube:
    """Return the minimum centered 2^level cube with the requested margin."""

    bounds = np.asarray(bounds_um, dtype=float).reshape(3, 2) * 1.0e-6
    lower, upper = bounds[:, 0], bounds[:, 1]
    span = upper - lower
    required_side = float(np.max(span) + 2.0 * margin_cells * dx_m)
    level = int(math.ceil(math.log2(required_side / dx_m)))
    cells = 2**level
    side = cells * dx_m
    center = (lower + upper) / 2.0
    origin = center - side / 2.0
    actual_margin = float(np.min(np.minimum(lower - origin, origin + side - upper)) / dx_m)
    if actual_margin + 1.0e-12 < margin_cells:
        raise FlowError("CFD_FLOW_MESH_RESOURCE_LIMIT", "Bounding cube margin construction failed")
    return BoundingCube(origin, side, level, cells, actual_margin)


def cells_across_diameter(partition: SurfacePartition, dx_um: float) -> dict[str, Any]:
    values = {
        item.label: 2.0 * item.equivalent_radius_um / dx_um
        for item in partition.patches
        if item.label != "wall"
    }
    array = np.asarray(list(values.values()), dtype=float)
    return {
        "basis": "four final cap equivalent diameters",
        "minimum": float(np.min(array)),
        "p05": float(np.percentile(array, 5)),
        "median": float(np.median(array)),
        "per_port": values,
    }


def parabolic_velocity(
    point: np.ndarray,
    center: np.ndarray,
    inward_normal: np.ndarray,
    equivalent_radius: float,
    maximum_velocity: float,
) -> np.ndarray:
    """Evaluate the requested clamped circular parabola for any 3-D normal."""

    normal = np.asarray(inward_normal, dtype=float)
    normal /= np.linalg.norm(normal)
    offset = np.asarray(point, dtype=float) - np.asarray(center, dtype=float)
    radial = offset - np.dot(offset, normal) * normal
    factor = max(0.0, 1.0 - float(np.dot(radial, radial)) / equivalent_radius**2)
    return normal * maximum_velocity * factor
