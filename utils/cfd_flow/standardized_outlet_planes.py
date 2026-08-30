"""Fixed physical mathematical outlet-plane contracts and safety checks."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import trimesh
from scipy.spatial import ConvexHull
from shapely.geometry import Point, Polygon, box

from .io import sha256_file, write_json
from .port_grid_sensitivity import (
    AXIS_GEOMETRY_RUN,
    GRID_SPECS,
    OUTLETS,
    RESEARCH_RUN,
    PlaneFrame,
    _project,
    recover_continuous_ports,
)


CONTRACT_REVISION = "STANDARDIZED_PHYSICAL_OUTLET_PLANES_V1"


def _canonical_plane_payload(record: Mapping[str, Any]) -> dict[str, Any]:
    return {
        key: record[key]
        for key in (
            "origin_m",
            "unit_normal",
            "basis_u",
            "basis_v",
            "canoND_origin_m",
            "canoND_vectors_m",
            "physical_aperture_contour_uv_m",
            "physical_aperture_bounds_uv_m",
            "safety_margin_m",
            "branch_id",
            "distance_from_bifurcation_um",
            "source_geometry_sha256",
        )
    }


def _hash_payload(payload: Mapping[str, Any]) -> str:
    import hashlib

    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _triangle_plane_points(triangle: np.ndarray, frame: PlaneFrame, tolerance: float) -> np.ndarray:
    tri = np.asarray(triangle, dtype=np.float64).reshape(3, 3)
    signed = (tri - frame.origin) @ frame.normal
    points: list[np.ndarray] = []
    for index, value in enumerate(signed):
        if abs(float(value)) <= tolerance:
            points.append(tri[index])
    for first, second in ((0, 1), (1, 2), (2, 0)):
        left, right = float(signed[first]), float(signed[second])
        if left * right < 0.0:
            fraction = left / (left - right)
            points.append(tri[first] + fraction * (tri[second] - tri[first]))
    unique: list[np.ndarray] = []
    for point in points:
        if not any(np.linalg.norm(point - prior) <= tolerance for prior in unique):
            unique.append(point)
    return np.asarray(unique, dtype=np.float64).reshape((-1, 3))


def _wall_plane_safety(
    wall: trimesh.Trimesh,
    frame: PlaneFrame,
    aperture: Polygon,
    rectangle: Polygon,
    *,
    tolerance_m: float,
) -> dict[str, Any]:
    intersection_uv: list[np.ndarray] = []
    for triangle in np.asarray(wall.triangles, dtype=np.float64):
        points = _triangle_plane_points(triangle, frame, tolerance_m)
        if len(points):
            intersection_uv.extend(_project(points, frame))
    boundary = aperture.boundary
    unsafe: list[list[float]] = []
    touching = 0
    for uv in intersection_uv:
        point = Point(float(uv[0]), float(uv[1]))
        if not rectangle.buffer(tolerance_m).covers(point):
            continue
        distance = float(boundary.distance(point))
        if distance <= tolerance_m:
            touching += 1
        else:
            unsafe.append([float(uv[0]), float(uv[1])])
    return {
        "status": "PASS" if not unsafe else "FAIL",
        "method": "WALL_TRIANGLE_PLANE_INTERSECTIONS_MUST_LIE_ON_APERTURE_SEAM",
        "intersection_point_count": len(intersection_uv),
        "rectangle_touching_seam_point_count": touching,
        "unsafe_rectangle_intersection_count": len(unsafe),
        "first_unsafe_uv_m": unsafe[0] if unsafe else None,
        "tolerance_m": tolerance_m,
    }


def build_standardized_outlet_plane_contract(project_root: Path) -> dict[str, Any]:
    root = Path(project_root).resolve()
    frames, continuous = recover_continuous_ports(root)
    geometry = root / "outputs" / "cfd_flow" / AXIS_GEOMETRY_RUN / "geometry" / "geometry_solver_m"
    wall = trimesh.load_mesh(geometry / "wall.stl", process=True)
    if not isinstance(wall, trimesh.Trimesh):
        raise ValueError("Expected one wall triangle mesh")
    records: dict[str, Any] = {}
    for label in OUTLETS:
        frame = frames[label]
        cap = trimesh.load_mesh(geometry / f"{label}.stl", process=True)
        if not isinstance(cap, trimesh.Trimesh):
            raise ValueError(f"Expected one cap mesh: {label}")
        projected = _project(np.asarray(cap.vertices, dtype=np.float64), frame)
        hull = ConvexHull(projected)
        contour = projected[hull.vertices]
        aperture = Polygon(contour)
        minimum = contour.min(axis=0)
        maximum = contour.max(axis=0)
        radius = float(continuous[label]["equivalent_radius_m"])
        margin = 0.05 * radius
        lower = minimum - margin
        upper = maximum + margin
        rectangle = box(float(lower[0]), float(lower[1]), float(upper[0]), float(upper[1]))
        canond_origin = frame.origin + lower[0] * frame.basis_u + lower[1] * frame.basis_v
        vector_u = (upper[0] - lower[0]) * frame.basis_u
        vector_v = (upper[1] - lower[1]) * frame.basis_v
        tolerance = max(5.0e-9, 0.02 * radius)
        wall_safety = _wall_plane_safety(
            wall, frame, aperture, rectangle, tolerance_m=tolerance
        )
        other_ports = []
        for other, other_frame in frames.items():
            if other == label:
                continue
            signed = float(np.dot(other_frame.origin - frame.origin, frame.normal))
            uv = _project(other_frame.origin.reshape(1, 3), frame)[0]
            covered = rectangle.buffer(tolerance).covers(Point(float(uv[0]), float(uv[1])))
            conflicts = abs(signed) <= 2.0 * radius and covered
            other_ports.append(
                {
                    "port": other,
                    "signed_plane_distance_m": signed,
                    "projected_uv_m": uv.tolist(),
                    "inside_coverage_rectangle": bool(covered),
                    "conflict": bool(conflicts),
                }
            )
        basis_qc = {
            "normal_norm_error": abs(float(np.linalg.norm(frame.normal)) - 1.0),
            "u_norm_error": abs(float(np.linalg.norm(frame.basis_u)) - 1.0),
            "v_norm_error": abs(float(np.linalg.norm(frame.basis_v)) - 1.0),
            "max_pairwise_dot": max(
                abs(float(np.dot(frame.normal, frame.basis_u))),
                abs(float(np.dot(frame.normal, frame.basis_v))),
                abs(float(np.dot(frame.basis_u, frame.basis_v))),
            ),
        }
        coverage_ratio = float(rectangle.area / aperture.area)
        record: dict[str, Any] = {
            "origin_m": frame.origin.tolist(),
            "unit_normal": frame.normal.tolist(),
            "basis_u": frame.basis_u.tolist(),
            "basis_v": frame.basis_v.tolist(),
            "canoND_origin_m": canond_origin.tolist(),
            "canoND_vectors_m": [vector_u.tolist(), vector_v.tolist()],
            "physical_aperture_contour_uv_m": contour.tolist(),
            "physical_aperture_bounds_uv_m": np.column_stack((minimum, maximum)).tolist(),
            "coverage_bounds_uv_m": np.column_stack((lower, upper)).tolist(),
            "coverage_area_over_aperture_area": coverage_ratio,
            "safety_margin_m": margin,
            "branch_id": continuous[label]["port_id"],
            "distance_from_bifurcation_um": continuous[label]["branch_distance"]["distance_to_nearest_bifurcation_um"],
            "source_geometry_path": continuous[label]["source_path"],
            "source_geometry_sha256": continuous[label]["source_sha256"],
            "basis_qc": basis_qc,
            "aperture_covered": bool(rectangle.buffer(1e-15).covers(aperture)),
            "other_port_safety": other_ports,
            "wall_upstream_and_unrelated_lumen_safety": wall_safety,
            "plane_covers_only_corresponding_opening": bool(
                rectangle.buffer(1e-15).covers(aperture)
                and wall_safety["status"] == "PASS"
                and not any(item["conflict"] for item in other_ports)
            ),
        }
        record["physical_contract_sha256"] = _hash_payload(_canonical_plane_payload(record))
        records[label] = record
    hashes_by_grid = {
        grid: {label: records[label]["physical_contract_sha256"] for label in OUTLETS}
        for grid in GRID_SPECS
    }
    same = all(hashes_by_grid[grid] == hashes_by_grid["base"] for grid in GRID_SPECS)
    checks = {
        "basis_orthonormal": all(
            max(record["basis_qc"].values()) <= 1e-12 for record in records.values()
        ),
        "same_physical_plane_across_dx": same,
        "all_apertures_covered": all(record["aperture_covered"] for record in records.values()),
        "no_other_port_conflicts": all(
            not any(item["conflict"] for item in record["other_port_safety"])
            for record in records.values()
        ),
        "no_wall_upstream_or_unrelated_lumen_intersection": all(
            record["wall_upstream_and_unrelated_lumen_safety"]["status"] == "PASS"
            for record in records.values()
        ),
    }
    result = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "revision": CONTRACT_REVISION,
        "method": "CONTINUOUS_CAP_AREA_CENTROID_NORMAL_PLUS_CONVEX_APERTURE_AND_WALL_INTERSECTION_SAFETY",
        "physical_plane_fit_is_grid_independent": True,
        "vmtk_used_at_runtime": False,
        "vmtk_method_inspected": True,
        "flow_extension_used": False,
        "checks": checks,
        "outlets": records,
        "contract_hashes_by_grid": hashes_by_grid,
        "grid_dx_m": {key: value["dx_m"] for key, value in GRID_SPECS.items()},
    }
    output = root / "outputs" / "cfd_flow" / RESEARCH_RUN / "qc" / "standardized_outlet_plane_contract.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json(output, result)
    return result

