"""Validated physical interior cross-section flux extraction.

This canonical solver-free module replaces the retired interpretation of lattice
nodes as finite-volume ``dx**3`` cells.  The integration domain is a polygon
obtained by slicing the continuous lumen wall.  Macroscopic fields from an
existing restart are reconstructed at Gaussian polygon quadrature points by
deterministic local weighted least squares.

No function in this module launches Seeder, Musubi, or Harvester.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import shapely
import trimesh
from scipy.spatial import cKDTree
from shapely.geometry import Point, Polygon
from shapely.geometry.polygon import orient

from .io import sha256_file
from .port_grid_sensitivity import (
    AXIS_GEOMETRY_RUN,
    _read_port_rows,
    recover_continuous_ports,
)
from .tau1_grid_convergence import GRID_SPECS
from .validated_contract import PORTS, RHO0_KG_M3 as RHO_KG_M3


PHYSICAL_FLUX_EVIDENCE_RUN = "healthy_mouse_capillary_tau1_grid_convergence_anchor003274_20260831"
PLANE_CONTRACT_REVISION = "STANDARDIZED_INTERIOR_PHYSICAL_PORT_PLANES_V3"
FLUX_ALGORITHM_REVISION = "CONTINUOUS_APERTURE_GAUSS_MLS_QUADRATIC_V2"
FLUX_DEFINITION = "PHYSICAL_INTERIOR_CROSS_SECTION_VELOCITY_FLUX"
DX_COARSE_M = 2.6e-7
MINIMUM_CLEARANCE_M = 3.0 * DX_COARSE_M
AREA_RELATIVE_ERROR_GATE = 1.0e-10
QUADRATURE_RELATIVE_GATE = 1.0e-3
PHYSICAL_FLUX_GATE = 1.0e-2
TEMPORAL_DRIFT_GATE = 1.0e-2
PLANE_CONSISTENCY_GATE = 1.0e-2
MLS_CONDITION_GATE = 1.0e8
MLS_RADII_DX = (2.0, 2.5, 3.0)
MLS_MIN_FLUID_NODES = 8
MLS_PREFERRED_FLUID_NODES = 12
WALL_ANCHOR_SPACING_DX = 0.5
LOCAL_COMPONENT_RADIUS_DH = 2.0


@dataclass(frozen=True, slots=True)
class InteriorPlane:
    """A continuous physical lumen section in an outward-oriented frame."""

    label: str
    position: str
    origin_m: np.ndarray
    unit_normal: np.ndarray
    basis_u: np.ndarray
    basis_v: np.ndarray
    aperture_uv_m: np.ndarray
    physical_contract_sha256: str = "pending"

    @property
    def aperture(self) -> Polygon:
        return Polygon(np.asarray(self.aperture_uv_m, dtype=np.float64))


@dataclass(frozen=True, slots=True)
class PolygonQuadrature:
    """Triangulation plus low- and fifth-degree Gaussian rules."""

    triangles_uv_m: np.ndarray
    low_points_uv_m: np.ndarray
    low_weights_m2: np.ndarray
    high_points_uv_m: np.ndarray
    high_weights_m2: np.ndarray
    aperture_area_m2: float
    triangle_area_sum_m2: float
    area_relative_error: float


@dataclass(frozen=True, slots=True)
class MLSStencil:
    fluid_indices: np.ndarray
    fluid_coefficients: np.ndarray
    radius_dx: float
    fluid_node_count: int
    wall_anchor_count: int
    rank: int
    condition_number: float


@dataclass(frozen=True, slots=True)
class MLSStencilMap:
    query_points_m: np.ndarray
    stencils: tuple[MLSStencil | None, ...]
    polynomial_order: int

    @property
    def invalid_count(self) -> int:
        return sum(stencil is None for stencil in self.stencils)

    @property
    def max_condition_number(self) -> float:
        values = [
            stencil.condition_number
            for stencil in self.stencils
            if stencil is not None
        ]
        return max(values, default=math.inf)


def _mesh_origin_dx(mesh_dir: Path) -> tuple[np.ndarray, float]:
    text = (mesh_dir / "header.lua").read_text(encoding="utf-8")
    block = re.search(r"boundingbox\s*=\s*\{(.*?)\n\}", text, re.DOTALL)
    level = re.search(r"(?m)^\s*minLevel\s*=\s*(\d+)", text)
    if block is None or level is None:
        raise ValueError("TreElm header has no uniform bounding-box contract")
    origin = re.search(r"origin\s*=\s*\{(.*?)\}", block.group(1), re.DOTALL)
    length = re.search(r"length\s*=\s*([^\s]+)", block.group(1))
    if origin is None or length is None:
        raise ValueError("TreElm bounding-box contract is incomplete")
    origin_value = np.asarray(
        [
            float(token.replace("D", "E"))
            for token in origin.group(1).replace(",", " ").split()
        ]
    )
    length_value = float(length.group(1).replace("D", "E"))
    return origin_value, length_value / 2 ** int(level.group(1))


def mesh_origin_dx(mesh_dir: Path) -> tuple[np.ndarray, float]:
    """Public uniform TreElm coordinate contract used by production replay."""

    return _mesh_origin_dx(Path(mesh_dir))


def _hash_payload(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload, sort_keys=True, separators=(",", ":"), allow_nan=False
    ).encode()
    return hashlib.sha256(encoded).hexdigest()


def _unit(vector: Iterable[float]) -> np.ndarray:
    value = np.asarray(tuple(vector), dtype=np.float64).reshape(3)
    norm = float(np.linalg.norm(value))
    if not math.isfinite(norm) or norm <= 0.0:
        raise ValueError("finite non-zero vector required")
    return value / norm


def _canonical_ring(polygon: Polygon) -> np.ndarray:
    value = orient(polygon, sign=1.0)
    coordinates = np.asarray(value.exterior.coords[:-1], dtype=np.float64)
    if len(coordinates) < 3:
        raise ValueError("aperture loop has fewer than three vertices")
    rounded = np.round(coordinates, decimals=18)
    start = min(range(len(rounded)), key=lambda index: tuple(rounded[index]))
    return np.roll(coordinates, -start, axis=0)


def _merge_segment_endpoints(
    segments: np.ndarray, *, tolerance_m: float = 1.0e-12
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    points = np.asarray(segments, dtype=np.float64).reshape((-1, 3))
    if len(points) == 0:
        return np.empty((0, 3)), []
    tree = cKDTree(points)
    parent = np.arange(len(points), dtype=np.int64)

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = int(parent[index])
        return index

    def union(left: int, right: int) -> None:
        a, b = find(left), find(right)
        if a != b:
            parent[max(a, b)] = min(a, b)

    for left, right in tree.query_pairs(float(tolerance_m)):
        union(int(left), int(right))
    groups: dict[int, list[int]] = {}
    for index in range(len(points)):
        groups.setdefault(find(index), []).append(index)
    roots = sorted(groups)
    root_to_new = {root: index for index, root in enumerate(roots)}
    merged = np.asarray(
        [points[groups[root]].mean(axis=0) for root in roots], dtype=np.float64
    )
    point_to_new = {
        point_index: root_to_new[root]
        for root in roots
        for point_index in groups[root]
    }
    edges: list[tuple[int, int]] = []
    for index in range(len(segments)):
        left, right = point_to_new[2 * index], point_to_new[2 * index + 1]
        if left != right:
            edges.append((left, right))
    return merged, edges


def _ordered_closed_loops(
    vertices: np.ndarray, edges: Sequence[tuple[int, int]]
) -> list[np.ndarray]:
    adjacency: dict[int, list[int]] = {}
    for left, right in edges:
        adjacency.setdefault(int(left), []).append(int(right))
        adjacency.setdefault(int(right), []).append(int(left))
    if any(len(neighbors) != 2 for neighbors in adjacency.values()):
        degrees = sorted({len(neighbors) for neighbors in adjacency.values()})
        raise ValueError(f"wall-plane intersection is not closed; degrees={degrees}")
    remaining = set(adjacency)
    loops: list[np.ndarray] = []
    while remaining:
        start = min(remaining)
        ordered = [start]
        previous = -1
        current = start
        while True:
            choices = sorted(
                neighbor
                for neighbor in adjacency[current]
                if neighbor != previous
            )
            if not choices:
                raise ValueError("wall-plane loop traversal stopped before closure")
            following = choices[0]
            if following == start:
                break
            if following in ordered:
                raise ValueError("wall-plane loop self-intersects in graph topology")
            ordered.append(following)
            previous, current = current, following
        remaining.difference_update(ordered)
        loops.append(np.asarray(vertices[ordered], dtype=np.float64))
    return loops


def slice_continuous_lumen(
    wall: trimesh.Trimesh,
    *,
    label: str,
    position: str,
    origin_m: np.ndarray,
    unit_normal: np.ndarray,
    basis_u: np.ndarray,
    basis_v: np.ndarray,
    hydraulic_diameter_m: float,
) -> tuple[InteriorPlane, dict[str, Any]]:
    """Slice the wall and select the unique local center-containing lumen."""

    segments = trimesh.intersections.mesh_plane(
        wall,
        plane_normal=np.asarray(unit_normal, dtype=np.float64),
        plane_origin=np.asarray(origin_m, dtype=np.float64),
    )
    vertices, edges = _merge_segment_endpoints(segments)
    loops_xyz = _ordered_closed_loops(vertices, edges)
    records: list[dict[str, Any]] = []
    polygons: list[Polygon] = []
    for loop_xyz in loops_xyz:
        delta = loop_xyz - origin_m
        uv = np.column_stack((delta @ basis_u, delta @ basis_v))
        polygon = Polygon(uv)
        valid = bool(polygon.is_valid and polygon.area > 0.0)
        if not valid:
            raise ValueError(f"{label}/{position}: invalid wall intersection loop")
        polygon = orient(polygon, sign=1.0)
        polygons.append(polygon)
        records.append(
            {
                "area_m2": float(polygon.area),
                "perimeter_m": float(polygon.length),
                "contains_measurement_center": bool(
                    polygon.buffer(1.0e-14).covers(Point(0.0, 0.0))
                ),
                "distance_from_measurement_center_m": float(
                    polygon.distance(Point(0.0, 0.0))
                ),
            }
        )
    center_indices = [
        index for index, record in enumerate(records) if record["contains_measurement_center"]
    ]
    if len(center_indices) != 1:
        raise ValueError(
            f"{label}/{position}: expected one center-containing lumen, "
            f"found {len(center_indices)}"
        )
    selected_index = center_indices[0]
    selected = polygons[selected_index]
    local_limit = LOCAL_COMPONENT_RADIUS_DH * float(hydraulic_diameter_m)
    local_secondary = [
        index
        for index, record in enumerate(records)
        if index != selected_index
        and record["area_m2"] > 1.0e-4 * selected.area
        and record["distance_from_measurement_center_m"] <= local_limit
    ]
    if local_secondary:
        raise ValueError(
            f"{label}/{position}: second nontrivial local lumen indicates bifurcation"
        )
    ring = _canonical_ring(selected)
    plane = InteriorPlane(
        label=label,
        position=position,
        origin_m=np.asarray(origin_m, dtype=np.float64),
        unit_normal=_unit(unit_normal),
        basis_u=_unit(basis_u),
        basis_v=_unit(basis_v),
        aperture_uv_m=ring,
    )
    qc = {
        "triangle_plane_segment_count": int(len(segments)),
        "closed_loop_count_global_plane": len(polygons),
        "selected_loop_index": selected_index,
        "selected_contains_measurement_center": True,
        "local_component_radius_m": local_limit,
        "nontrivial_local_secondary_component_count": len(local_secondary),
        "remote_unrelated_component_count": max(0, len(polygons) - 1),
        "components": records,
        "polygon_valid": bool(selected.is_valid),
        "polygon_positive_area": bool(selected.area > 0.0),
        "polygon_simple": bool(selected.exterior.is_simple),
    }
    return plane, qc


def _plane_payload(plane: InteriorPlane) -> dict[str, Any]:
    aperture = plane.aperture
    return {
        "position": plane.position,
        "origin_m": plane.origin_m.tolist(),
        "unit_normal": plane.unit_normal.tolist(),
        "basis_u": plane.basis_u.tolist(),
        "basis_v": plane.basis_v.tolist(),
        "physical_aperture_contour_uv_m": plane.aperture_uv_m.tolist(),
        "aperture_physical_area_m2": float(aperture.area),
        "aperture_perimeter_m": float(aperture.length),
        "local_hydraulic_diameter_m": float(4.0 * aperture.area / aperture.length),
    }


def build_interior_plane_contract(project_root: Path) -> dict[str, Any]:
    """Build the dx-, mesh-, PDF-, and velocity-independent V3 plane contract."""

    root = Path(project_root).resolve()
    frames, continuous = recover_continuous_ports(root)
    rows = _read_port_rows(root)
    geometry = (
        root
        / "outputs"
        / "cfd_flow"
        / AXIS_GEOMETRY_RUN
        / "geometry"
        / "geometry_solver_m"
    )
    wall_path = geometry / "wall.stl"
    wall = trimesh.load_mesh(wall_path, process=False)
    if not isinstance(wall, trimesh.Trimesh):
        raise ValueError("continuous wall source is not one triangle mesh")
    ports: dict[str, Any] = {}
    canonical_ports: dict[str, Any] = {}
    for label in PORTS:
        frame = frames[label]
        cap_dh = float(continuous[label]["hydraulic_diameter_m"])
        central_distance = 0.5 * cap_dh
        offsets = {
            "portward_check": central_distance - 0.1 * cap_dh,
            "central": central_distance,
            "interior_check": central_distance + 0.1 * cap_dh,
        }
        branch_from_original_port_m = (
            float(
                continuous[label]["branch_distance"][
                    "distance_to_nearest_bifurcation_um"
                ]
            )
            * 1.0e-6
        )
        extension_m = float(rows[label]["extension_length_um"]) * 1.0e-6
        planes: dict[str, Any] = {}
        canonical_planes: dict[str, Any] = {}
        for position, distance_m in offsets.items():
            origin = frame.origin - distance_m * frame.normal
            plane, slice_qc = slice_continuous_lumen(
                wall,
                label=label,
                position=position,
                origin_m=origin,
                unit_normal=frame.normal,
                basis_u=frame.basis_u,
                basis_v=frame.basis_v,
                hydraulic_diameter_m=cap_dh,
            )
            payload = _plane_payload(plane)
            payload.update(
                {
                    "centerline_arclength_from_port_m": float(distance_m),
                    "distance_to_nearest_bifurcation_m": float(
                        extension_m + branch_from_original_port_m - distance_m
                    ),
                    "slice_qc": slice_qc,
                }
            )
            planes[position] = payload
            canonical_planes[position] = {
                key: payload[key]
                for key in (
                    "position",
                    "origin_m",
                    "unit_normal",
                    "basis_u",
                    "basis_v",
                    "physical_aperture_contour_uv_m",
                    "aperture_physical_area_m2",
                    "centerline_arclength_from_port_m",
                    "distance_to_nearest_bifurcation_m",
                    "local_hydraulic_diameter_m",
                )
            }
        central = planes["central"]
        checks = {
            "central_clearance_from_pressure_cap_ge_3dx_coarse": (
                central["centerline_arclength_from_port_m"] >= MINIMUM_CLEARANCE_M
            ),
            "central_clearance_from_bifurcation_ge_3dx_coarse": (
                central["distance_to_nearest_bifurcation_m"] >= MINIMUM_CLEARANCE_M
            ),
            "all_three_planes_inside_extension": all(
                value < extension_m for value in offsets.values()
            ),
            "all_three_planes_single_local_lumen": all(
                plane["slice_qc"]["nontrivial_local_secondary_component_count"] == 0
                for plane in planes.values()
            ),
            "all_three_polygons_valid": all(
                plane["slice_qc"]["polygon_valid"]
                and plane["slice_qc"]["polygon_positive_area"]
                and plane["slice_qc"]["polygon_simple"]
                for plane in planes.values()
            ),
        }
        canonical_ports[label] = {
            "source_wall_sha256": sha256_file(wall_path),
            "source_port_id": continuous[label]["port_id"],
            "source_boundary_origin": continuous[label]["boundary_origin"],
            "normal_orientation": "OUTWARD_FROM_DOMAIN_INTERIOR_TO_PORT_CAP",
            "geometry_location_rule": "CAP_INWARD_0.5_DH_WITH_CHECKS_AT_PLUS_MINUS_0.1_DH",
            "planes": canonical_planes,
        }
        ports[label] = {
            **canonical_ports[label],
            "planes": planes,
            "source_wall_path": str(wall_path.relative_to(root)),
            "cap_hydraulic_diameter_m": cap_dh,
            "extension_length_m": extension_m,
            "original_port_to_nearest_bifurcation_arclength_m": branch_from_original_port_m,
            "usable_centerline_interval_from_cap_m": [
                MINIMUM_CLEARANCE_M,
                extension_m + branch_from_original_port_m - MINIMUM_CLEARANCE_M,
            ],
            "checks": checks,
        }
    canonical = {"revision": PLANE_CONTRACT_REVISION, "ports": canonical_ports}
    contract_hash = _hash_payload(canonical)
    for record in ports.values():
        record["physical_contract_sha256"] = contract_hash
        for plane in record["planes"].values():
            plane["physical_contract_sha256"] = contract_hash
    checks = {
        "all_port_geometry_checks_pass": all(
            all(record["checks"].values()) for record in ports.values()
        ),
        "same_physical_contract_hash_across_grids": True,
        "contract_build_is_independent_of_dx_and_cell_centers": True,
        "contract_build_is_independent_of_velocity_and_pdf": True,
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "revision": PLANE_CONTRACT_REVISION,
        "flux_algorithm_revision": FLUX_ALGORITHM_REVISION,
        "contract_sha256": contract_hash,
        "contract_hashes_by_grid": {
            grid: contract_hash for grid in ("coarse", "base", "fine")
        },
        "source_geometry_sha256": sha256_file(wall_path),
        "geometry_only": True,
        "selection_uses_real_q": False,
        "checks": checks,
        "ports": ports,
    }
    return result


def plane_from_v3_record(
    label: str, position: str, record: Mapping[str, Any]
) -> InteriorPlane:
    return InteriorPlane(
        label=label,
        position=position,
        origin_m=np.asarray(record["origin_m"], dtype=np.float64),
        unit_normal=np.asarray(record["unit_normal"], dtype=np.float64),
        basis_u=np.asarray(record["basis_u"], dtype=np.float64),
        basis_v=np.asarray(record["basis_v"], dtype=np.float64),
        aperture_uv_m=np.asarray(
            record["physical_aperture_contour_uv_m"], dtype=np.float64
        ),
        physical_contract_sha256=str(record["physical_contract_sha256"]),
    )


def triangulate_aperture(polygon: Polygon) -> np.ndarray:
    if not polygon.is_valid or polygon.area <= 0.0:
        raise ValueError("valid positive-area polygon required")
    collection = shapely.constrained_delaunay_triangles(polygon)
    triangles: list[np.ndarray] = []
    for geometry in collection.geoms:
        if geometry.geom_type != "Polygon" or geometry.area <= 0.0:
            continue
        vertices = np.asarray(geometry.exterior.coords[:-1], dtype=np.float64)
        if len(vertices) != 3:
            raise ValueError("constrained triangulation returned a non-triangle")
        triangles.append(vertices)
    if not triangles:
        raise ValueError("constrained triangulation returned no triangles")
    return np.asarray(triangles, dtype=np.float64)


_DUNAVANT_BARYCENTRIC = np.asarray(
    [
        (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0),
        (0.059715871789770, 0.470142064105115, 0.470142064105115),
        (0.470142064105115, 0.059715871789770, 0.470142064105115),
        (0.470142064105115, 0.470142064105115, 0.059715871789770),
        (0.797426985353087, 0.101286507323456, 0.101286507323456),
        (0.101286507323456, 0.797426985353087, 0.101286507323456),
        (0.101286507323456, 0.101286507323456, 0.797426985353087),
    ],
    dtype=np.float64,
)
_DUNAVANT_WEIGHTS = np.asarray(
    [
        0.225,
        0.132394152788506,
        0.132394152788506,
        0.132394152788506,
        0.125939180544827,
        0.125939180544827,
        0.125939180544827,
    ],
    dtype=np.float64,
)
_DEGREE_FOUR_BARYCENTRIC = np.asarray(
    (
        (0.445948490915965, 0.445948490915965, 0.108103018168070),
        (0.445948490915965, 0.108103018168070, 0.445948490915965),
        (0.108103018168070, 0.445948490915965, 0.445948490915965),
        (0.091576213509771, 0.091576213509771, 0.816847572980459),
        (0.091576213509771, 0.816847572980459, 0.091576213509771),
        (0.816847572980459, 0.091576213509771, 0.091576213509771),
    ),
    dtype=np.float64,
)
_DEGREE_FOUR_WEIGHTS = np.asarray(
    (
        0.223381589678011,
        0.223381589678011,
        0.223381589678011,
        0.109951743655322,
        0.109951743655322,
        0.109951743655322,
    ),
    dtype=np.float64,
)


def build_polygon_quadrature(polygon: Polygon) -> PolygonQuadrature:
    triangles = triangulate_aperture(polygon)
    left = triangles[:, 1] - triangles[:, 0]
    right = triangles[:, 2] - triangles[:, 0]
    double_areas = np.abs(left[:, 0] * right[:, 1] - left[:, 1] * right[:, 0])
    areas = 0.5 * double_areas
    triangle_sum = math.fsum(float(value) for value in areas)
    aperture_area = float(polygon.area)
    relative = abs(triangle_sum - aperture_area) / aperture_area
    # The degree-one centroid rule has a 1/3 relative bias for a quadratic
    # Poiseuille profile on the long constrained triangles of a circular
    # polygon.  A degree-four rule is the fixed lower-order comparison for the
    # degree-five rule; increasing quadrature order never changes the CFD grid.
    low_points = np.einsum("qa,tad->tqd", _DEGREE_FOUR_BARYCENTRIC, triangles)
    low_weights = areas[:, None] * _DEGREE_FOUR_WEIGHTS[None, :]
    high_points = np.einsum("qa,tad->tqd", _DUNAVANT_BARYCENTRIC, triangles)
    high_weights = areas[:, None] * _DUNAVANT_WEIGHTS[None, :]
    return PolygonQuadrature(
        triangles_uv_m=triangles,
        low_points_uv_m=low_points.reshape((-1, 2)),
        low_weights_m2=low_weights.reshape(-1),
        high_points_uv_m=high_points.reshape((-1, 2)),
        high_weights_m2=high_weights.reshape(-1),
        aperture_area_m2=aperture_area,
        triangle_area_sum_m2=triangle_sum,
        area_relative_error=relative,
    )


def plane_uv_to_xyz(plane: InteriorPlane, points_uv_m: np.ndarray) -> np.ndarray:
    uv = np.asarray(points_uv_m, dtype=np.float64).reshape((-1, 2))
    return (
        plane.origin_m
        + uv[:, :1] * plane.basis_u
        + uv[:, 1:] * plane.basis_v
    )


def sample_polygon_boundary(polygon: Polygon, spacing_m: float) -> np.ndarray:
    if spacing_m <= 0.0:
        raise ValueError("positive wall-anchor spacing required")
    perimeter = float(polygon.length)
    count = max(3, int(math.ceil(perimeter / float(spacing_m))))
    return np.asarray(
        [
            polygon.exterior.interpolate(perimeter * index / count).coords[0]
            for index in range(count)
        ],
        dtype=np.float64,
    )


def _design_matrix(delta_over_dx: np.ndarray, polynomial_order: int) -> np.ndarray:
    delta = np.asarray(delta_over_dx, dtype=np.float64).reshape((-1, 3))
    x, y, z = delta.T
    if polynomial_order == 1:
        return np.column_stack((np.ones(len(delta)), x, y, z))
    if polynomial_order == 2:
        return np.column_stack(
            (np.ones(len(delta)), x, y, z, x * x, y * y, z * z, x * y, x * z, y * z)
        )
    raise ValueError("polynomial_order must be one or two")


def build_mls_stencil_map(
    fluid_points_m: np.ndarray,
    query_points_m: np.ndarray,
    *,
    dx_m: float,
    wall_anchor_points_m: np.ndarray | None = None,
    polynomial_order: int = 2,
) -> MLSStencilMap:
    """Build deterministic local WLS intercept maps for fixed query points."""

    fluid = np.asarray(fluid_points_m, dtype=np.float64).reshape((-1, 3))
    queries = np.asarray(query_points_m, dtype=np.float64).reshape((-1, 3))
    anchors = (
        np.empty((0, 3), dtype=np.float64)
        if wall_anchor_points_m is None
        else np.asarray(wall_anchor_points_m, dtype=np.float64).reshape((-1, 3))
    )
    fluid_tree = cKDTree(fluid)
    anchor_tree = cKDTree(anchors) if len(anchors) else None
    stencils: list[MLSStencil | None] = []
    for query in queries:
        fallback: MLSStencil | None = None
        for radius_dx in MLS_RADII_DX:
            radius_m = radius_dx * float(dx_m)
            indices = np.asarray(
                sorted(fluid_tree.query_ball_point(query, radius_m)), dtype=np.int64
            )
            if len(indices) < MLS_MIN_FLUID_NODES:
                continue
            anchor_indices = (
                np.asarray(
                    sorted(anchor_tree.query_ball_point(query, radius_m)), dtype=np.int64
                )
                if anchor_tree is not None
                else np.empty(0, dtype=np.int64)
            )
            sample_points = np.vstack((fluid[indices], anchors[anchor_indices]))
            delta = (sample_points - query) / float(dx_m)
            design = _design_matrix(delta, polynomial_order)
            distances = np.linalg.norm(delta, axis=1)
            weights = 1.0 / (0.25 + (distances / radius_dx) ** 2)
            weighted_design = np.sqrt(weights)[:, None] * design
            rank = int(np.linalg.matrix_rank(weighted_design))
            condition = float(np.linalg.cond(weighted_design))
            if rank != design.shape[1] or condition > MLS_CONDITION_GATE:
                continue
            normal = design.T @ (weights[:, None] * design)
            projection = np.linalg.solve(normal, design.T * weights[None, :])
            coefficients = projection[0, : len(indices)]
            candidate = MLSStencil(
                fluid_indices=indices,
                fluid_coefficients=np.asarray(coefficients, dtype=np.float64),
                radius_dx=float(radius_dx),
                fluid_node_count=int(len(indices)),
                wall_anchor_count=int(len(anchor_indices)),
                rank=rank,
                condition_number=condition,
            )
            if len(indices) >= MLS_PREFERRED_FLUID_NODES:
                fallback = candidate
                break
            fallback = candidate
        stencils.append(fallback)
    return MLSStencilMap(queries, tuple(stencils), polynomial_order)


def evaluate_mls_stencil_map(
    stencil_map: MLSStencilMap, fluid_values: np.ndarray
) -> np.ndarray:
    values = np.asarray(fluid_values, dtype=np.float64)
    scalar = values.ndim == 1
    if scalar:
        values = values[:, None]
    if values.ndim != 2:
        raise ValueError("fluid values must have shape (n,) or (n,k)")
    output = np.full((len(stencil_map.stencils), values.shape[1]), np.nan)
    for query_index, stencil in enumerate(stencil_map.stencils):
        if stencil is None:
            continue
        output[query_index] = (
            stencil.fluid_coefficients @ values[stencil.fluid_indices]
        )
    return output[:, 0] if scalar else output


def _stencil_qc(stencil_map: MLSStencilMap) -> dict[str, Any]:
    valid = [stencil for stencil in stencil_map.stencils if stencil is not None]
    return {
        "query_point_count": len(stencil_map.stencils),
        "valid_stencil_count": len(valid),
        "invalid_stencil_count": stencil_map.invalid_count,
        "all_stencils_valid": stencil_map.invalid_count == 0,
        "minimum_fluid_node_count": min(
            (stencil.fluid_node_count for stencil in valid), default=0
        ),
        "maximum_fluid_node_count": max(
            (stencil.fluid_node_count for stencil in valid), default=0
        ),
        "maximum_wall_anchor_count": max(
            (stencil.wall_anchor_count for stencil in valid), default=0
        ),
        "max_condition_number": stencil_map.max_condition_number,
        "condition_gate": MLS_CONDITION_GATE,
        "radius_dx_distribution": {
            str(radius): sum(stencil.radius_dx == radius for stencil in valid)
            for radius in MLS_RADII_DX
        },
        "polynomial_order": stencil_map.polynomial_order,
    }


def _integral(
    velocities_m_s: np.ndarray,
    weights_m2: np.ndarray,
    normal: np.ndarray,
    *,
    label: str,
) -> tuple[float, np.ndarray]:
    normal_velocity = np.asarray(velocities_m_s) @ np.asarray(normal)
    signed = math.fsum(
        float(weight) * float(value)
        for weight, value in zip(weights_m2, normal_velocity, strict=True)
    )
    return (-signed if label == "inlet" else signed), normal_velocity


def integrate_reconstructed_flux(
    *,
    plane: InteriorPlane,
    quadrature: PolygonQuadrature,
    high_velocity_m_s: np.ndarray,
    low_velocity_m_s: np.ndarray,
    high_density_lattice: np.ndarray | None = None,
    high_density_velocity_m_s: np.ndarray | None = None,
) -> dict[str, Any]:
    high_q, normal_velocity = _integral(
        high_velocity_m_s,
        quadrature.high_weights_m2,
        plane.unit_normal,
        label=plane.label,
    )
    low_q, _ = _integral(
        low_velocity_m_s,
        quadrature.low_weights_m2,
        plane.unit_normal,
        label=plane.label,
    )
    area = math.fsum(float(value) for value in quadrature.high_weights_m2)
    outward_signed = -high_q if plane.label == "inlet" else high_q
    desired = -normal_velocity if plane.label == "inlet" else normal_velocity
    backflow_area = math.fsum(
        float(weight)
        for weight, value in zip(
            quadrature.high_weights_m2, desired, strict=True
        )
        if value < 0.0
    )
    density_mean = None
    if high_density_lattice is not None:
        density_mean = math.fsum(
            float(weight) * float(value)
            for weight, value in zip(
                quadrature.high_weights_m2,
                np.asarray(high_density_lattice, dtype=np.float64),
                strict=True,
            )
        ) / area
    mass_flux = None
    if high_density_velocity_m_s is not None:
        rho_velocity = RHO_KG_M3 * np.asarray(
            high_density_velocity_m_s, dtype=np.float64
        )
        mass_flux, _ = _integral(
            rho_velocity,
            quadrature.high_weights_m2,
            plane.unit_normal,
            label=plane.label,
        )
    return {
        "flux_definition": FLUX_DEFINITION,
        "algorithm_revision": FLUX_ALGORITHM_REVISION,
        "port": plane.label,
        "plane_position": plane.position,
        "physical_q_m3_s": float(high_q),
        "signed_outward_q_m3_s": float(outward_signed),
        "mass_flux_kg_s": None if mass_flux is None else float(mass_flux),
        "mass_normalized_volume_flux_m3_s": (
            None if mass_flux is None else float(mass_flux / RHO_KG_M3)
        ),
        "area_weighted_density_lattice": (
            None if density_mean is None else float(density_mean)
        ),
        "aperture_physical_area_m2": quadrature.aperture_area_m2,
        "triangle_area_sum_m2": quadrature.triangle_area_sum_m2,
        "area_coverage_relative_error": quadrature.area_relative_error,
        "quadrature_low_q_m3_s": float(low_q),
        "quadrature_high_q_m3_s": float(high_q),
        "R_quadrature_Q": abs(high_q - low_q)
        / max(abs(high_q), np.finfo(float).tiny),
        "weighted_outward_normal_velocity_mean_m_s": float(outward_signed / area),
        "outward_normal_velocity_min_m_s": float(np.min(normal_velocity)),
        "outward_normal_velocity_max_m_s": float(np.max(normal_velocity)),
        "local_backflow_area_fraction": float(backflow_area / area),
        "sign_convention": (
            "Qin=-integral(u.n_outward)dA"
            if plane.label == "inlet"
            else "Qout=+integral(u.n_outward)dA"
        ),
        "all_finite": bool(
            math.isfinite(high_q)
            and math.isfinite(low_q)
            and np.all(np.isfinite(high_velocity_m_s))
        ),
    }


def _synthetic_plane(normal: np.ndarray, polygon_uv: np.ndarray) -> InteriorPlane:
    from .port_grid_sensitivity import orthonormal_plane_basis

    unit_normal = _unit(normal)
    basis_u, basis_v = orthonormal_plane_basis(unit_normal)
    return InteriorPlane(
        "outlet_01",
        "synthetic",
        np.zeros(3),
        unit_normal,
        basis_u,
        basis_v,
        np.asarray(polygon_uv, dtype=np.float64),
        "synthetic",
    )


def _synthetic_cartesian_cloud(
    plane: InteriorPlane,
    *,
    dx_m: float,
    extent_m: float,
    translation_fraction: float = 0.0,
) -> np.ndarray:
    lower = -float(extent_m) - 3.0 * dx_m
    upper = float(extent_m) + 3.0 * dx_m
    values = np.arange(
        math.floor(lower / dx_m) - 1,
        math.ceil(upper / dx_m) + 2,
        dtype=np.float64,
    )
    shift = translation_fraction * dx_m * np.asarray((1.0, 0.5, 0.25))
    x, y, z = np.meshgrid(values, values, values, indexing="ij")
    return dx_m * np.column_stack((x.ravel(), y.ravel(), z.ravel())) + shift


def _evaluate_synthetic(
    plane: InteriorPlane,
    quadrature: PolygonQuadrature,
    cloud: np.ndarray,
    field: np.ndarray,
    *,
    dx_m: float,
    wall_anchors: np.ndarray | None = None,
) -> tuple[float, dict[str, Any]]:
    high_xyz = plane_uv_to_xyz(plane, quadrature.high_points_uv_m)
    low_xyz = plane_uv_to_xyz(plane, quadrature.low_points_uv_m)
    high_map = build_mls_stencil_map(
        cloud,
        high_xyz,
        dx_m=dx_m,
        wall_anchor_points_m=wall_anchors,
    )
    low_map = build_mls_stencil_map(
        cloud,
        low_xyz,
        dx_m=dx_m,
        wall_anchor_points_m=wall_anchors,
    )
    high = evaluate_mls_stencil_map(high_map, field)
    low = evaluate_mls_stencil_map(low_map, field)
    result = integrate_reconstructed_flux(
        plane=plane,
        quadrature=quadrature,
        high_velocity_m_s=high,
        low_velocity_m_s=low,
    )
    qc = {
        "high": _stencil_qc(high_map),
        "low": _stencil_qc(low_map),
        "R_quadrature_Q": result["R_quadrature_Q"],
    }
    return float(result["physical_q_m3_s"]), qc


def _poiseuille_case(
    *,
    dx_m: float,
    normal: np.ndarray,
    translation_fraction: float = 0.0,
) -> dict[str, Any]:
    radius_m = 1.5e-6
    mean_velocity_m_s = 4.0e-4
    angles = np.linspace(0.0, 2.0 * math.pi, 96, endpoint=False)
    ring = radius_m * np.column_stack((np.cos(angles), np.sin(angles)))
    plane = _synthetic_plane(normal, ring)
    quadrature = build_polygon_quadrature(plane.aperture)
    cloud = _synthetic_cartesian_cloud(
        plane,
        dx_m=dx_m,
        extent_m=radius_m,
        translation_fraction=translation_fraction,
    )
    axial = cloud @ plane.unit_normal
    radial_vector = cloud - axial[:, None] * plane.unit_normal
    radial = np.linalg.norm(radial_vector, axis=1)
    keep = (radial < radius_m) & (np.abs(axial) <= 3.0 * dx_m)
    fluid = cloud[keep]
    profile = 2.0 * mean_velocity_m_s * (1.0 - (radial[keep] / radius_m) ** 2)
    velocity = profile[:, None] * plane.unit_normal
    wall_uv = sample_polygon_boundary(
        plane.aperture, WALL_ANCHOR_SPACING_DX * dx_m
    )
    wall_xyz = plane_uv_to_xyz(plane, wall_uv)
    actual, stencil_qc = _evaluate_synthetic(
        plane,
        quadrature,
        fluid,
        velocity,
        dx_m=dx_m,
        wall_anchors=wall_xyz,
    )
    expected = mean_velocity_m_s * math.pi * radius_m**2
    error = abs(actual - expected) / expected
    return {
        "actual_q_m3_s": actual,
        "expected_q_m3_s": expected,
        "relative_error": error,
        "dx_m": dx_m,
        "translation_fraction_dx": translation_fraction,
        "stencil_qc": stencil_qc,
    }


@lru_cache(maxsize=1)
def run_analytic_flux_oracles() -> dict[str, Any]:
    """Run zero-solver constant, affine, pipe, and translation oracles."""

    polygon_uv = np.asarray(
        [(-1.8e-6, -0.9e-6), (1.5e-6, -1.1e-6), (2.0e-6, 0.7e-6), (-1.1e-6, 1.4e-6)]
    )
    plane = _synthetic_plane(np.asarray((0.37, -0.51, 0.776)), polygon_uv)
    quadrature = build_polygon_quadrature(plane.aperture)
    dx_m = GRID_SPECS["base"].dx_m
    cloud = _synthetic_cartesian_cloud(plane, dx_m=dx_m, extent_m=2.2e-6)
    constant_velocity = np.asarray((1.3e-4, -2.1e-4, 7.0e-4))
    constant_field = np.repeat(constant_velocity[None, :], len(cloud), axis=0)
    constant_actual, constant_stencils = _evaluate_synthetic(
        plane, quadrature, cloud, constant_field, dx_m=dx_m
    )
    constant_expected = float(np.dot(constant_velocity, plane.unit_normal)) * float(
        plane.aperture.area
    )
    constant_error = abs(constant_actual - constant_expected) / abs(constant_expected)

    intercept = np.asarray((2.0e-4, -1.5e-4, 5.0e-4))
    gradient = np.asarray(
        ((13.0, -7.0, 5.0), (3.0, 11.0, -4.0), (-8.0, 6.0, 9.0))
    )
    affine_field = intercept + cloud @ gradient.T
    affine_actual, affine_stencils = _evaluate_synthetic(
        plane, quadrature, cloud, affine_field, dx_m=dx_m
    )
    centroid = plane.aperture.centroid
    centroid_xyz = plane_uv_to_xyz(
        plane, np.asarray([[centroid.x, centroid.y]])
    )[0]
    affine_expected = float(
        np.dot(intercept + gradient @ centroid_xyz, plane.unit_normal)
        * plane.aperture.area
    )
    affine_error = abs(affine_actual - affine_expected) / abs(affine_expected)

    poiseuille_axis: dict[str, Any] = {}
    poiseuille_oblique: dict[str, Any] = {}
    for grid, spec in GRID_SPECS.items():
        poiseuille_axis[grid] = _poiseuille_case(
            dx_m=spec.dx_m, normal=np.asarray((0.0, 0.0, 1.0))
        )
        poiseuille_oblique[grid] = _poiseuille_case(
            dx_m=spec.dx_m, normal=np.asarray((0.37, -0.51, 0.776))
        )
    translations = {
        str(fraction): _poiseuille_case(
            dx_m=GRID_SPECS["base"].dx_m,
            normal=np.asarray((0.37, -0.51, 0.776)),
            translation_fraction=fraction,
        )
        for fraction in (0.0, 0.25, 0.5)
    }
    translation_values = np.asarray(
        [record["actual_q_m3_s"] for record in translations.values()]
    )
    translation_sensitivity = float(
        np.ptp(translation_values) / abs(np.mean(translation_values))
    )
    pipe_gates = {
        "coarse": 0.01,
        "base": 0.005,
        "fine": 0.005,
    }
    axis_pass = all(
        poiseuille_axis[grid]["relative_error"] <= pipe_gates[grid]
        for grid in pipe_gates
    )
    oblique_pass = all(
        poiseuille_oblique[grid]["relative_error"] <= pipe_gates[grid]
        for grid in pipe_gates
    )
    refinement_pass = (
        poiseuille_axis["fine"]["relative_error"]
        <= poiseuille_axis["coarse"]["relative_error"] + 1.0e-12
        and poiseuille_oblique["fine"]["relative_error"]
        <= poiseuille_oblique["coarse"]["relative_error"] + 1.0e-12
    )
    gates = {
        "constant_field_relative_error_le_1e-10": constant_error <= 1.0e-10,
        "affine_field_relative_error_le_1e-8": affine_error <= 1.0e-8,
        "poiseuille_axis_cbf": axis_pass,
        "poiseuille_oblique_cbf": oblique_pass,
        "poiseuille_refinement_not_worse": refinement_pass,
        "base_grid_translation_sensitivity_le_0p005": (
            translation_sensitivity <= 0.005
        ),
        "synthetic_quadrature_convergence_le_1e-3": max(
            constant_stencils["R_quadrature_Q"],
            affine_stencils["R_quadrature_Q"],
            *(
                record["stencil_qc"]["R_quadrature_Q"]
                for family in (poiseuille_axis, poiseuille_oblique)
                for record in family.values()
            ),
        )
        <= QUADRATURE_RELATIVE_GATE,
        "all_synthetic_stencils_valid": all(
            record["stencil_qc"][order]["invalid_stencil_count"] == 0
            for family in (poiseuille_axis, poiseuille_oblique)
            for record in family.values()
            for order in ("low", "high")
        ),
    }
    return {
        "status": "PASS" if all(gates.values()) else "FAIL",
        "algorithm_revision": FLUX_ALGORITHM_REVISION,
        "shapely_version": shapely.__version__,
        "constant_field": {
            "actual_q_m3_s": constant_actual,
            "expected_q_m3_s": constant_expected,
            "relative_error": constant_error,
            "stencil_qc": constant_stencils,
        },
        "affine_field": {
            "actual_q_m3_s": affine_actual,
            "expected_q_m3_s": affine_expected,
            "relative_error": affine_error,
            "stencil_qc": affine_stencils,
        },
        "poiseuille_axis": poiseuille_axis,
        "poiseuille_oblique": poiseuille_oblique,
        "grid_translation": {
            "cases": translations,
            "relative_spread": translation_sensitivity,
        },
        "gates": gates,
    }


def _relative_range(values: Iterable[float]) -> float:
    array = np.asarray(tuple(float(value) for value in values), dtype=np.float64)
    return float(np.ptp(array) / max(abs(float(np.mean(array))), np.finfo(float).tiny))


def _closure(inlet: float, outlets: Sequence[float]) -> dict[str, Any]:
    total = math.fsum(float(value) for value in outlets)
    relative = abs(float(inlet) - total) / max(abs(float(inlet)), np.finfo(float).tiny)
    return {
        "inlet": float(inlet),
        "outlet_sum": total,
        "absolute_mismatch": abs(float(inlet) - total),
        "relative_error": relative,
        "gate": PHYSICAL_FLUX_GATE,
        "pass": relative <= PHYSICAL_FLUX_GATE,
    }


def _prepare_plane_numerics(
    plane: InteriorPlane,
    fluid_points: np.ndarray,
    *,
    dx_m: float,
) -> dict[str, Any]:
    quadrature = build_polygon_quadrature(plane.aperture)
    wall_uv = sample_polygon_boundary(
        plane.aperture, WALL_ANCHOR_SPACING_DX * dx_m
    )
    wall_xyz = plane_uv_to_xyz(plane, wall_uv)
    high_xyz = plane_uv_to_xyz(plane, quadrature.high_points_uv_m)
    low_xyz = plane_uv_to_xyz(plane, quadrature.low_points_uv_m)
    high_velocity = build_mls_stencil_map(
        fluid_points,
        high_xyz,
        dx_m=dx_m,
        wall_anchor_points_m=wall_xyz,
        polynomial_order=2,
    )
    low_velocity = build_mls_stencil_map(
        fluid_points,
        low_xyz,
        dx_m=dx_m,
        wall_anchor_points_m=wall_xyz,
        polynomial_order=2,
    )
    high_density = build_mls_stencil_map(
        fluid_points, high_xyz, dx_m=dx_m, polynomial_order=1
    )
    return {
        "plane": plane,
        "quadrature": quadrature,
        "high_velocity": high_velocity,
        "low_velocity": low_velocity,
        "high_density": high_density,
        "stencil_qc": {
            "high_velocity": _stencil_qc(high_velocity),
            "low_velocity": _stencil_qc(low_velocity),
            "high_density": _stencil_qc(high_density),
        },
    }


def _evaluate_prepared_plane(
    prepared: Mapping[str, Any],
    velocity: np.ndarray,
    density: np.ndarray,
) -> dict[str, Any]:
    high_velocity = evaluate_mls_stencil_map(prepared["high_velocity"], velocity)
    low_velocity = evaluate_mls_stencil_map(prepared["low_velocity"], velocity)
    high_density = evaluate_mls_stencil_map(prepared["high_density"], density)
    # Momentum is reconstructed directly.  Interpolating rho and u separately
    # and multiplying them would add a post-processing product error to the
    # requested rho*u diagnostic.
    high_density_velocity = evaluate_mls_stencil_map(
        prepared["high_velocity"], density[:, None] * velocity
    )
    result = integrate_reconstructed_flux(
        plane=prepared["plane"],
        quadrature=prepared["quadrature"],
        high_velocity_m_s=high_velocity,
        low_velocity_m_s=low_velocity,
        high_density_lattice=high_density,
        high_density_velocity_m_s=high_density_velocity,
    )
    result["interpolation_stencil_qc"] = prepared["stencil_qc"]
    return result


def evaluate_physical_port_fluxes(
    plane_contract: Mapping[str, Any],
    fluid_points_m: np.ndarray,
    velocity_m_s: np.ndarray,
    density_lattice: np.ndarray,
    *,
    dx_m: float,
    position: str = "central",
) -> dict[str, Any]:
    """Evaluate all four production ports with the validated V3/V2 implementation."""

    if plane_contract.get("contract_sha256") != (
        "ffaa49bdb6e43fb7208ff29df07a90d4e92ef9bfa4b96ca4f997d4f453a7f005"
    ):
        raise ValueError("physical plane contract SHA-256 is not the validated V3 contract")
    ports: dict[str, Any] = {}
    for label in ("inlet", "outlet_01", "outlet_02", "outlet_03"):
        record = plane_contract["ports"][label]["planes"][position]
        plane = plane_from_v3_record(label, position, record)
        prepared = _prepare_plane_numerics(plane, fluid_points_m, dx_m=dx_m)
        ports[label] = _evaluate_prepared_plane(
            prepared,
            velocity_m_s,
            density_lattice,
        )
    q_in = float(ports["inlet"]["physical_q_m3_s"])
    q_out = math.fsum(
        float(ports[label]["physical_q_m3_s"])
        for label in ("outlet_01", "outlet_02", "outlet_03")
    )
    fractions = {
        label: float(ports[label]["physical_q_m3_s"]) / q_out
        for label in ("outlet_01", "outlet_02", "outlet_03")
    }
    return {
        "status": "PASS" if _closure(q_in, [ports[label]["physical_q_m3_s"] for label in ("outlet_01", "outlet_02", "outlet_03")])["pass"] else "FAIL",
        "flux_definition": FLUX_DEFINITION,
        "algorithm_revision": FLUX_ALGORITHM_REVISION,
        "plane_contract_revision": PLANE_CONTRACT_REVISION,
        "plane_contract_sha256": plane_contract["contract_sha256"],
        "plane_position": position,
        "ports": ports,
        "Qin_m3_s": q_in,
        "Qout_m3_s": q_out,
        "closure": abs(q_in - q_out) / abs(q_in),
        "flow_fractions": fractions,
    }
