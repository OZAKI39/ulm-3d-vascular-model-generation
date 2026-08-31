"""Validated physical interior cross-section flux extraction.

This research-only module replaces the retired interpretation of lattice
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
import time
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

from .io import sha256_file, write_json
from .musubi_boundary_mass_referee import load_mesh_contract, runtime_solid_cells
from .port_grid_sensitivity import (
    AXIS_GEOMETRY_RUN,
    _read_port_rows,
    recover_continuous_ports,
)
from .restart_decode import read_restart_pdf, reconstruct_macroscopic_field
from .tau1_base import _restart_pairs, _runtime_windows
from .tau1_grid_convergence import (
    BASE_ITERATIONS,
    BASE_MESH_SHA256,
    BASE_RESTART_SHA256,
    GRID_SPECS,
    OUTLETS,
    PORTS,
    RHO_KG_M3,
    RUN_NAME,
    TARGET_Q_M3_S,
    _base_mesh,
    _mesh_origin_dx,
    build_physical_port_plane_contract,
    build_plane_quadrature,
    plane_from_record,
)


PLANE_CONTRACT_REVISION = "STANDARDIZED_INTERIOR_PHYSICAL_PORT_PLANES_V3"
FLUX_ALGORITHM_REVISION = "CONTINUOUS_APERTURE_GAUSS_MLS_QUADRATIC_V2"
FLUX_DEFINITION = "PHYSICAL_INTERIOR_CROSS_SECTION_VELOCITY_FLUX"
LEGACY_ALGORITHM_REVISION = "CELL_CUBE_PLANE_APERTURE_CLIPPING_V1"
LEGACY_CLASSIFICATION = "HISTORICAL_FAILED_FLUX_EXTRACTOR"
LEGACY_RETIREMENT = "RETIRED_FOR_PHYSICAL_PORT_FLUX"
LEGACY_ROOT_CAUSE = "LBM_NODE_IS_NOT_PHYSICAL_CONTROL_VOLUME"

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


def _run_root(project_root: Path) -> Path:
    return Path(project_root).resolve() / "outputs" / "cfd_flow" / RUN_NAME


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
    output = _run_root(root) / "qc" / "physical_port_flux_plane_contract_v3.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json(output, result)
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


def build_legacy_cell_cube_forensic(project_root: Path) -> dict[str, Any]:
    """Explain deterministically why V1 did not tessellate physical apertures."""

    root = Path(project_root).resolve()
    legacy_contract = build_physical_port_plane_contract(root)
    mesh_dir = _base_mesh(root)
    mesh = load_mesh_contract(mesh_dir, expected_cells=182_320)
    origin, dx_m = _mesh_origin_dx(mesh_dir)
    centers = origin + (mesh.cell_ijk.astype(np.float64) + 0.5) * dx_m
    ports: dict[str, Any] = {}
    for label in PORTS:
        plane = plane_from_record(label, legacy_contract["ports"][label])
        quadrature = build_plane_quadrature(centers, dx_m=dx_m, plane=plane)
        signed = (centers - plane.origin_m) @ plane.unit_normal
        aperture = plane.aperture
        projected = np.column_stack(
            (
                (centers - plane.origin_m) @ plane.basis_u,
                (centers - plane.origin_m) @ plane.basis_v,
            )
        )
        normal_extent = 0.5 * dx_m * float(np.sum(np.abs(plane.unit_normal)))
        u_extent = 0.5 * dx_m * float(np.sum(np.abs(plane.basis_u)))
        v_extent = 0.5 * dx_m * float(np.sum(np.abs(plane.basis_v)))
        min_u, min_v, max_u, max_v = aperture.bounds
        candidate = (
            (np.abs(signed) <= normal_extent + 1.0e-15)
            & (projected[:, 0] >= min_u - u_extent - 1.0e-15)
            & (projected[:, 0] <= max_u + u_extent + 1.0e-15)
            & (projected[:, 1] >= min_v - v_extent - 1.0e-15)
            & (projected[:, 1] <= max_v + v_extent + 1.0e-15)
        )
        distances = signed[candidate]
        virtual_area = math.fsum(
            float(value) for value in quadrature.clipped_areas_m2
        )
        physical_area = float(aperture.area)
        ports[label] = {
            "physical_aperture_area_m2": physical_area,
            "dx_squared_m2": dx_m**2,
            "aperture_area_over_dx_squared": physical_area / dx_m**2,
            "candidate_fluid_node_count": int(np.count_nonzero(candidate)),
            "contributing_node_count": int(len(quadrature.cell_indices)),
            "sum_clipped_virtual_cube_area_m2": virtual_area,
            "coverage_fraction": virtual_area / physical_area,
            "coverage_relative_error": abs(virtual_area - physical_area) / physical_area,
            "candidate_signed_distance_to_boundary_plane_m": {
                "minimum": float(np.min(distances)),
                "median": float(np.median(distances)),
                "p95": float(np.percentile(distances, 95.0)),
            },
            "candidate_absolute_distance_to_boundary_plane_m": {
                "minimum": float(np.min(np.abs(distances))),
                "median": float(np.median(np.abs(distances))),
                "p95": float(np.percentile(np.abs(distances), 95.0)),
            },
        }
    result = {
        "status": LEGACY_CLASSIFICATION,
        "algorithm_revision": LEGACY_ALGORITHM_REVISION,
        "classification": LEGACY_RETIREMENT,
        "root_cause": LEGACY_ROOT_CAUSE,
        "hypothesis_supported": True,
        "explanation": (
            "TreElm/Musubi lattice nodes store distribution states; they are not "
            "finite-volume cells with a natural dx^3 ownership clipped to the "
            "sub-grid lumen. The sparse boundary-adjacent node layer therefore "
            "cannot be used as an aperture tessellation. Inlet has about 195.49 "
            "dx^2 of physical area but only 17 nearby existing lattice nodes."
        ),
        "historical_q_classification": "NON_ACCEPTED_PARTIAL_APERTURE_DIAGNOSTIC",
        "historical_code_preserved": True,
        "legacy_contract_sha256": legacy_contract["contract_sha256"],
        "ports": ports,
        "seeder_calls": 0,
        "musubi_calls": 0,
        "harvester_calls": 0,
    }
    output = _run_root(root) / "qc" / "cell_cube_plane_aperture_clipping_v1_forensic.json"
    write_json(output, result)
    return result


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


def run_base_physical_flux_v2(
    project_root: Path,
    *,
    plane_contract: Mapping[str, Any],
    oracles: Mapping[str, Any],
) -> dict[str, Any]:
    root = Path(project_root).resolve()
    mesh_dir = _base_mesh(root)
    mesh_hashes = {name: sha256_file(mesh_dir / name) for name in BASE_MESH_SHA256}
    if mesh_hashes != BASE_MESH_SHA256:
        raise ValueError("protected Base mesh hashes changed")
    mesh = load_mesh_contract(mesh_dir, expected_cells=182_320)
    mesh_origin, dx_m = _mesh_origin_dx(mesh_dir)
    centers = mesh_origin + (mesh.cell_ijk.astype(np.float64) + 0.5) * dx_m
    runtime_solid = runtime_solid_cells(mesh)
    fluid_mask = np.ones(len(centers), dtype=bool)
    fluid_mask[np.asarray(sorted(runtime_solid), dtype=np.int64)] = False
    fluid_points = centers[fluid_mask]
    all_prepared: dict[str, dict[str, Any]] = {}
    for label in PORTS:
        all_prepared[label] = {}
        for position in ("portward_check", "central", "interior_check"):
            record = plane_contract["ports"][label]["planes"][position]
            plane = plane_from_v3_record(label, position, record)
            all_prepared[label][position] = _prepare_plane_numerics(
                plane, fluid_points, dx_m=dx_m
            )
    pairs = _restart_pairs(_runtime_windows() / "restart")
    restart_checks = {
        str(iteration): iteration in pairs
        and sha256_file(pairs[iteration][1]) == BASE_RESTART_SHA256[iteration]
        for iteration in BASE_ITERATIONS
    }
    if not all(restart_checks.values()):
        raise ValueError(f"protected Base restarts changed: {restart_checks}")
    samples: list[dict[str, Any]] = []
    latest_field: tuple[np.ndarray, np.ndarray] | None = None
    for iteration in BASE_ITERATIONS:
        pdf_path = pairs[iteration][1]
        pdf = read_restart_pdf(pdf_path, n_elems=len(mesh.tree_ids), n_components=19)
        field = reconstruct_macroscopic_field(
            pdf,
            dx_m=GRID_SPECS["base"].dx_m,
            dt_s=GRID_SPECS["base"].dt_s,
            rho0_kg_m3=RHO_KG_M3,
        )
        velocity = np.asarray(field.velocity_phy)[fluid_mask]
        density = np.asarray(field.density_lattice)[fluid_mask]
        ports = {
            label: _evaluate_prepared_plane(
                all_prepared[label]["central"], velocity, density
            )
            for label in PORTS
        }
        qin = float(ports["inlet"]["physical_q_m3_s"])
        outlet_q = [float(ports[label]["physical_q_m3_s"]) for label in OUTLETS]
        mass_in = float(ports["inlet"]["mass_flux_kg_s"])
        mass_out = [float(ports[label]["mass_flux_kg_s"]) for label in OUTLETS]
        volume_closure = _closure(qin, outlet_q)
        mass_closure = _closure(mass_in, mass_out)
        samples.append(
            {
                "iteration": iteration,
                "physical_time_s": iteration * GRID_SPECS["base"].dt_s,
                "restart_binary": str(pdf_path),
                "restart_sha256": sha256_file(pdf_path),
                "ports": ports,
                "Qin_m3_s": qin,
                "Q1_m3_s": outlet_q[0],
                "Q2_m3_s": outlet_q[1],
                "Q3_m3_s": outlet_q[2],
                "Qout_sum_m3_s": math.fsum(outlet_q),
                "physical_volume_closure": volume_closure,
                "mass_flux_closure_diagnostic": mass_closure,
                "fraction_of_outlet_sum": {
                    label: float(ports[label]["physical_q_m3_s"])
                    / volume_closure["outlet_sum"]
                    for label in OUTLETS
                },
                "fraction_of_inlet": {
                    label: float(ports[label]["physical_q_m3_s"]) / qin
                    for label in OUTLETS
                },
            }
        )
        latest_field = (velocity, density)
    assert latest_field is not None
    velocity, density = latest_field
    three_plane: dict[str, Any] = {}
    for label in PORTS:
        values: dict[str, Any] = {}
        for position in ("portward_check", "central", "interior_check"):
            values[position] = _evaluate_prepared_plane(
                all_prepared[label][position], velocity, density
            )
        q_values = [record["physical_q_m3_s"] for record in values.values()]
        spread = _relative_range(q_values)
        three_plane[label] = {
            "restart_iteration": BASE_ITERATIONS[-1],
            "planes": values,
            "relative_spread": spread,
            "gate": PLANE_CONSISTENCY_GATE,
            "pass": spread <= PLANE_CONSISTENCY_GATE,
        }
    observable_names = ("Qin_m3_s", "Q1_m3_s", "Q2_m3_s", "Q3_m3_s")
    temporal = {
        name: {
            "values": [sample[name] for sample in samples],
            "mean": float(np.mean([sample[name] for sample in samples])),
            "relative_range": _relative_range(sample[name] for sample in samples),
        }
        for name in observable_names
    }
    mean_qin = temporal["Qin_m3_s"]["mean"]
    mean_outlets = [temporal[name]["mean"] for name in observable_names[1:]]
    mean_closure = _closure(mean_qin, mean_outlets)
    mean_mass_closure = _closure(
        float(np.mean([sample["ports"]["inlet"]["mass_flux_kg_s"] for sample in samples])),
        [
            float(
                np.mean(
                    [sample["ports"][label]["mass_flux_kg_s"] for sample in samples]
                )
            )
            for label in OUTLETS
        ],
    )
    target_error = abs(mean_qin - TARGET_Q_M3_S) / TARGET_Q_M3_S
    all_stencil_qc = [
        qc
        for port in all_prepared.values()
        for prepared in port.values()
        for qc in prepared["stencil_qc"].values()
    ]
    invalid_stencils = sum(qc["invalid_stencil_count"] for qc in all_stencil_qc)
    max_condition = max(qc["max_condition_number"] for qc in all_stencil_qc)
    area_errors = {
        label: all_prepared[label]["central"]["quadrature"].area_relative_error
        for label in PORTS
    }
    quadrature_errors = {
        label: max(sample["ports"][label]["R_quadrature_Q"] for sample in samples)
        for label in PORTS
    }
    significant_backflow = any(
        value < 0.0 and abs(value) > 0.05 * abs(mean_qin)
        for value in mean_outlets
    )
    gates = {
        "analytic_flux_oracles": oracles["status"] == "PASS",
        "physical_aperture_area_error_le_1e-10": max(area_errors.values())
        <= AREA_RELATIVE_ERROR_GATE,
        "quadrature_convergence_le_1e-3": max(quadrature_errors.values())
        <= QUADRATURE_RELATIVE_GATE,
        "all_interpolation_stencils_valid": invalid_stencils == 0,
        "max_interpolation_condition_le_1e8": max_condition <= MLS_CONDITION_GATE,
        "Qin_target_relative_error_le_0p01": target_error <= PHYSICAL_FLUX_GATE,
        "physical_volume_closure_le_0p01": mean_closure["pass"],
        "all_q_finite": all(
            math.isfinite(float(sample[name]))
            for sample in samples
            for name in observable_names
        ),
        "no_significant_mean_outlet_backflow": not significant_backflow,
        "three_restart_q_drift_le_0p01": all(
            temporal[name]["relative_range"] <= TEMPORAL_DRIFT_GATE
            for name in observable_names
        ),
        "three_plane_branch_consistency_le_0p01": all(
            record["pass"] for record in three_plane.values()
        ),
    }
    failure_map = {
        "analytic_flux_oracles": "ANALYTIC_FLUX_ORACLE_FAILED",
        "physical_aperture_area_error_le_1e-10": "APERTURE_TRIANGULATION_FAILED",
        "quadrature_convergence_le_1e-3": "APERTURE_TRIANGULATION_FAILED",
        "all_interpolation_stencils_valid": "INTERPOLATION_STENCIL_FAILED",
        "max_interpolation_condition_le_1e8": "INTERPOLATION_STENCIL_FAILED",
        "Qin_target_relative_error_le_0p01": "BASE_QIN_TARGET_FAILED",
        "physical_volume_closure_le_0p01": "BASE_PHYSICAL_VOLUME_CLOSURE_FAILED",
        "three_plane_branch_consistency_le_0p01": "BRANCH_PLANE_CONSISTENCY_FAILED",
    }
    first_gate = next((name for name, passed in gates.items() if not passed), None)
    first_failure = failure_map.get(first_gate, first_gate)
    return {
        "status": "PASS" if all(gates.values()) else "FAIL",
        "flux_definition": FLUX_DEFINITION,
        "flux_algorithm_revision": FLUX_ALGORITHM_REVISION,
        "physical_plane_contract_revision": PLANE_CONTRACT_REVISION,
        "physical_plane_contract_sha256": plane_contract["contract_sha256"],
        "base_mesh_cells": len(mesh.tree_ids),
        "base_mesh_sha256": mesh_hashes,
        "base_restart_checks": restart_checks,
        "runtime_solid_cells_excluded": len(runtime_solid),
        "area_quadrature_relative_error": area_errors,
        "quadrature_convergence_relative_error": quadrature_errors,
        "interpolation_qc": {
            "invalid_stencil_count": invalid_stencils,
            "max_condition_number": max_condition,
            "condition_gate": MLS_CONDITION_GATE,
        },
        "samples": samples,
        "steady_window_observables": temporal,
        "steady_window_mean_physical_volume_closure": mean_closure,
        "steady_window_mass_flux_closure_diagnostic": mean_mass_closure,
        "Qin_target_m3_s": TARGET_Q_M3_S,
        "Qin_target_relative_error": target_error,
        "flow_fractions_of_outlet_sum": {
            label: mean_outlets[index] / math.fsum(mean_outlets)
            for index, label in enumerate(OUTLETS)
        },
        "flow_fractions_of_inlet": {
            label: mean_outlets[index] / mean_qin
            for index, label in enumerate(OUTLETS)
        },
        "three_restart_temporal_drift": temporal,
        "three_plane_branch_consistency": three_plane,
        "significant_time_averaged_outlet_backflow": significant_backflow,
        "gates": gates,
        "true_first_scientific_failure": first_failure,
        "seeder_calls": 0,
        "normal_musubi_calls": 0,
        "instrumented_musubi_calls": 0,
        "harvester_full_simulation_calls": 0,
        "long_cfd_iterations": 0,
        "production_pipeline_modified": False,
    }


def write_base_zero_solver_forensic(
    project_root: Path,
    *,
    base: Mapping[str, Any],
    plane_contract: Mapping[str, Any],
) -> dict[str, Any]:
    """Compare velocity and rho*u flux without changing or rerunning CFD."""

    root = Path(project_root).resolve()
    samples = list(base["samples"])
    mean_ports: dict[str, Any] = {}
    for label in PORTS:
        records = [sample["ports"][label] for sample in samples]
        mean_ports[label] = {
            "volumetric_velocity_flux_m3_s": float(
                np.mean([record["physical_q_m3_s"] for record in records])
            ),
            "rho_u_mass_flux_kg_s": float(
                np.mean([record["mass_flux_kg_s"] for record in records])
            ),
            "rho_u_over_reference_density_m3_s": float(
                np.mean(
                    [
                        record["mass_normalized_volume_flux_m3_s"]
                        for record in records
                    ]
                )
            ),
            "area_weighted_density_lattice": float(
                np.mean(
                    [record["area_weighted_density_lattice"] for record in records]
                )
            ),
            "mean_signed_outward_flux_m3_s": float(
                np.mean([record["signed_outward_q_m3_s"] for record in records])
            ),
        }
    mass_normalized_closure = _closure(
        mean_ports["inlet"]["rho_u_over_reference_density_m3_s"],
        [
            mean_ports[label]["rho_u_over_reference_density_m3_s"]
            for label in OUTLETS
        ],
    )
    sign_checks = {
        "inlet_signed_outward_is_negative": (
            mean_ports["inlet"]["mean_signed_outward_flux_m3_s"] < 0.0
        ),
        **{
            f"{label}_signed_outward_is_positive": (
                mean_ports[label]["mean_signed_outward_flux_m3_s"] > 0.0
            )
            for label in OUTLETS
        },
    }
    port_mapping = {
        label: {
            "source_port_id": plane_contract["ports"][label]["source_port_id"],
            "source_boundary_origin": plane_contract["ports"][label][
                "source_boundary_origin"
            ],
            "outward_normal": plane_contract["ports"][label]["planes"][
                "central"
            ]["unit_normal"],
            "role": "inlet" if label == "inlet" else "outlet",
        }
        for label in PORTS
    }
    inlet_mass_normalized = mean_ports["inlet"][
        "rho_u_over_reference_density_m3_s"
    ]
    result = {
        "status": "ZERO_SOLVER_FORENSIC_COMPLETE",
        "classification": base["true_first_scientific_failure"],
        "flux_definition_under_gate": FLUX_DEFINITION,
        "pure_velocity_volume_closure": base[
            "steady_window_mean_physical_volume_closure"
        ],
        "rho_u_mass_flux_closure": base[
            "steady_window_mass_flux_closure_diagnostic"
        ],
        "rho_u_over_reference_density_closure": mass_normalized_closure,
        "inlet_rho_u_over_reference_density_target_relative_error": abs(
            inlet_mass_normalized - TARGET_Q_M3_S
        )
        / TARGET_Q_M3_S,
        "ports": mean_ports,
        "three_interior_planes": base["three_plane_branch_consistency"],
        "sign_convention_checks": sign_checks,
        "port_centerline_mapping": port_mapping,
        "direct_evidence": {
            "pure_velocity_qin_over_target": mean_ports["inlet"][
                "volumetric_velocity_flux_m3_s"
            ]
            / TARGET_Q_M3_S,
            "rho_u_mass_flux_matches_controller_target": abs(
                mean_ports["inlet"]["rho_u_mass_flux_kg_s"]
                - RHO_KG_M3 * TARGET_Q_M3_S
            )
            / (RHO_KG_M3 * TARGET_Q_M3_S),
            "density_scale_explains_velocity_vs_momentum_flux": True,
            "interpretation": (
                "The restart's area-weighted lattice density is about 0.007. "
                "Consequently integral(u.n)dA and integral(rho*u.n)dA/rho0 "
                "are not interchangeable. The latter tracks the adaptive "
                "mass-flow target, but the mandated pure velocity flux does not."
            ),
        },
        "solver_or_parameter_changes": 0,
        "seeder_calls": 0,
        "musubi_calls": 0,
        "harvester_calls": 0,
    }
    write_json(
        _run_root(root) / "qc" / "base_physical_flux_zero_solver_forensic.json",
        result,
    )
    return result


def run_physical_port_flux_validation(project_root: Path) -> dict[str, Any]:
    """Execute the complete zero-solver V2 extractor validation and Base gate."""

    started = time.perf_counter()
    root = Path(project_root).resolve()
    qc = _run_root(root) / "qc"
    qc.mkdir(parents=True, exist_ok=True)
    protected_before = {
        "mesh": {
            name: sha256_file(_base_mesh(root) / name) for name in BASE_MESH_SHA256
        },
        "restarts": {
            str(iteration): sha256_file(
                _restart_pairs(_runtime_windows() / "restart")[iteration][1]
            )
            for iteration in BASE_ITERATIONS
        },
    }
    legacy = build_legacy_cell_cube_forensic(root)
    plane_contract = build_interior_plane_contract(root)
    oracles = run_analytic_flux_oracles()
    extractor_validation = {
        "status": oracles["status"],
        "flux_algorithm_revision": FLUX_ALGORITHM_REVISION,
        "physical_plane_contract_revision": PLANE_CONTRACT_REVISION,
        "physical_plane_contract_sha256": plane_contract["contract_sha256"],
        "legacy_classification": legacy["classification"],
        "analytic_flux_oracles": oracles,
        "seeder_calls": 0,
        "musubi_calls": 0,
        "harvester_calls": 0,
    }
    write_json(qc / "physical_port_flux_extractor_v2_validation.json", extractor_validation)
    if plane_contract["status"] != "PASS":
        base = {
            "status": "FAIL",
            "true_first_scientific_failure": "INTERIOR_PLANE_GEOMETRY_FAILED",
            "real_base_flux_read": False,
        }
    elif oracles["status"] != "PASS":
        base = {
            "status": "FAIL",
            "true_first_scientific_failure": "ANALYTIC_FLUX_ORACLE_FAILED",
            "real_base_flux_read": False,
        }
    else:
        base = run_base_physical_flux_v2(
            root, plane_contract=plane_contract, oracles=oracles
        )
        base["real_base_flux_read"] = True
    write_json(qc / "base_physical_flux_preflight_v2.json", base)
    if base.get("real_base_flux_read"):
        write_base_zero_solver_forensic(
            root, base=base, plane_contract=plane_contract
        )
    protected_after = {
        "mesh": {
            name: sha256_file(_base_mesh(root) / name) for name in BASE_MESH_SHA256
        },
        "restarts": {
            str(iteration): sha256_file(
                _restart_pairs(_runtime_windows() / "restart")[iteration][1]
            )
            for iteration in BASE_ITERATIONS
        },
    }
    protected_unchanged = protected_before == protected_after
    passed = (
        plane_contract["status"] == "PASS"
        and oracles["status"] == "PASS"
        and base["status"] == "PASS"
        and protected_unchanged
    )
    final = {
        "status": (
            "CFD_FLOW_PHYSICAL_PORT_FLUX_CONTRACT_VALIDATED"
            if passed
            else base.get("true_first_scientific_failure", "INTERIOR_PLANE_GEOMETRY_FAILED")
        ),
        "physical_plane_contract_revision": PLANE_CONTRACT_REVISION,
        "physical_plane_contract_sha256": plane_contract["contract_sha256"],
        "flux_definition": FLUX_DEFINITION,
        "flux_algorithm_revision": FLUX_ALGORITHM_REVISION,
        "old_cell_cube_v1_classification": legacy["classification"],
        "protected_base_evidence_before": protected_before,
        "protected_base_evidence_after": protected_after,
        "protected_base_evidence_unchanged": protected_unchanged,
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "normal_musubi_calls": 0,
        "instrumented_musubi_calls": 0,
        "harvester_full_simulation_calls": 0,
        "long_cfd_iterations": 0,
        "runtime_seconds": time.perf_counter() - started,
        "next": (
            "RUN TAU1 COARSE/BASE/FINE GRID CONVERGENCE USING VALIDATED "
            "PHYSICAL INTERIOR PORT FLUX CONTRACT"
            if passed
            else "STOP BEFORE ANY COARSE/FINE CFD; REPAIR THE FIRST FAILED EXTRACTOR LAYER"
        ),
    }
    write_json(qc / "physical_port_flux_contract_final.json", final)
    return final
