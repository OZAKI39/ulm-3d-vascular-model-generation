"""Topology, boundary, collision, and unchanged-core quality controls."""

from __future__ import annotations

from typing import Any, Iterable

import numpy as np
import trimesh

from utils.cfd_lumen.ultraliser_qc import _triangle_intersections

from .config import LocalCutConfig, SurfaceQCConfig
from .io import BoundaryInput
from .types import BoundarySurfaceResult, TaggedSurface


def surface_topology_qc(
    surface: TaggedSurface, config: SurfaceQCConfig
) -> tuple[dict[str, Any], list[tuple[int, int]]]:
    """Evaluate the closed combined surface without silently repairing it."""

    mesh = surface.mesh()
    sorted_edges = np.sort(np.asarray(mesh.edges, dtype=np.int64), axis=1)
    _, edge_counts = np.unique(sorted_edges, axis=0, return_counts=True)
    boundary_edges = int(np.count_nonzero(edge_counts == 1))
    nonmanifold_edges = int(np.count_nonzero(edge_counts > 2))
    diagonal = float(np.linalg.norm(np.ptp(surface.vertices, axis=0)))
    area_tolerance = max(np.finfo(float).eps * diagonal**2 * 100.0, 1.0e-18)
    areas = np.asarray(mesh.area_faces, dtype=float)
    repeated = np.asarray(
        [len(set(map(int, face))) < 3 for face in surface.faces], dtype=bool
    )
    degenerate = int(np.count_nonzero((areas <= area_tolerance) | repeated))
    intersections, candidates = _triangle_intersections(
        mesh, np.arange(len(mesh.faces), dtype=np.int64)
    )
    component_count = len(mesh.split(only_watertight=False))
    checks = {
        "single_component": not config.require_single_component
        or component_count == 1,
        "watertight": not config.require_watertight or bool(mesh.is_watertight),
        "zero_boundary_edges": not config.require_zero_boundary_edges
        or boundary_edges == 0,
        "zero_nonmanifold_edges": not config.require_zero_nonmanifold_edges
        or nonmanifold_edges == 0,
        "zero_self_intersections": not config.require_zero_self_intersections
        or not intersections,
        "zero_degenerate_triangles": not config.require_zero_degenerate_triangles
        or degenerate == 0,
    }
    report = {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "vertex_count": int(len(surface.vertices)),
        "triangle_count": int(len(surface.faces)),
        "component_count": int(component_count),
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "boundary_edge_count": boundary_edges,
        "nonmanifold_edge_count": nonmanifold_edges,
        "self_intersection_count": len(intersections),
        "self_intersection_candidate_pairs_checked": candidates,
        "degenerate_triangle_count": degenerate,
        "degenerate_area_tolerance_um2": area_tolerance,
        "surface_area_um2": float(mesh.area),
        "volume_um3": float(abs(mesh.volume)),
        "bounds_um": np.asarray(mesh.bounds, dtype=float).tolist(),
    }
    return report, intersections


def boundary_geometry_qc(
    boundaries: Iterable[BoundaryInput],
    cut_reports: Iterable[dict[str, Any]],
    results: Iterable[BoundarySurfaceResult],
    config: SurfaceQCConfig,
) -> dict[str, Any]:
    boundary_by_index = {item.index: item for item in boundaries}
    cuts = list(cut_reports)
    result_list = list(results)
    records: list[dict[str, Any]] = []
    for cut, result in zip(cuts, result_list, strict=True):
        boundary = boundary_by_index[result.boundary_index]
        checks = {
            "exactly_one_simple_closed_loop": cut["cut_loop_count"] == 1,
            "single_local_surface_component": cut[
                "candidate_surface_component_count"
            ]
            == 1,
            "loop_center_matches_boundary": cut["loop_center_distance_um"]
            <= boundary.source_radius_um,
            "loop_plane_residual": cut["plane_residual_um"]
            <= config.maximum_cap_planarity_error_um,
            "positive_finite_cap_area": bool(
                np.isfinite(result.actual_cap_area_um2)
                and result.actual_cap_area_um2 > 0
            ),
            "equivalent_radius": cut["equivalent_radius_relative_error"]
            <= config.maximum_equivalent_radius_relative_error,
            "extension_length": result.extension_length_error_um
            <= config.maximum_extension_length_error_um,
            "extension_axis": result.extension_axis_dot
            >= config.minimum_normal_dot,
            "distal_cap_planarity": result.cap_planarity_error_um
            <= config.maximum_cap_planarity_error_um,
            "distal_cap_normal": result.minimum_cap_normal_dot
            >= config.minimum_normal_dot,
        }
        records.append(
            {
                "port_id": boundary.port_id,
                "boundary_origin": boundary.boundary_origin,
                "role": boundary.role,
                "status": "PASS" if all(checks.values()) else "FAIL",
                "checks": checks,
                **cut,
                **result.report(),
            }
        )
    count_checks = {
        "expected_boundary_count": len(records) == config.expected_boundary_count,
        "expected_tagged_cap_count": len(
            {item.boundary_index for item in result_list}
        )
        == config.expected_boundary_count,
        "one_inlet": sum(item.role == "ASSUMED_INLET" for item in result_list) == 1,
        "expected_outlet_count": sum(
            item.role == "ASSUMED_OUTLET" for item in result_list
        )
        == config.expected_boundary_count - 1,
    }
    passed = all(count_checks.values()) and all(
        item["status"] == "PASS" for item in records
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "checks": count_checks,
        "boundary_count": len(records),
        "boundaries": records,
    }


def core_surface_preservation_qc(
    original: trimesh.Trimesh,
    final: trimesh.Trimesh,
    boundaries: Iterable[BoundaryInput],
    local: LocalCutConfig,
    qc: SurfaceQCConfig,
) -> dict[str, Any]:
    """Densely compare original vertices and centroids outside surgery cylinders."""

    samples = np.vstack(
        (
            np.asarray(original.vertices, dtype=float),
            np.asarray(original.triangles_center, dtype=float),
        )
    )
    keep = np.ones(len(samples), dtype=bool)
    for boundary in boundaries:
        relative = samples - boundary.center_um
        axial = relative @ boundary.outward_normal
        radial = np.linalg.norm(
            relative - np.outer(axial, boundary.outward_normal), axis=1
        )
        in_zone = (
            (radial <= local.local_radial_radius_factor * boundary.source_radius_um)
            & (
                axial
                >= -local.local_axial_back_radius_factor * boundary.source_radius_um
            )
            & (
                axial
                <= local.local_axial_forward_radius_factor
                * boundary.source_radius_um
            )
        )
        keep &= ~in_zone
    core = samples[keep]
    _, distances, _ = trimesh.proximity.closest_point(final, core)
    maximum = float(np.max(distances)) if len(distances) else float("inf")
    p95 = float(np.percentile(distances, 95)) if len(distances) else float("inf")
    checks = {
        "dense_core_samples_available": len(core) > 0,
        "maximum_distance": maximum <= qc.maximum_core_surface_distance_um,
        "p95_distance": p95 <= qc.maximum_core_surface_p95_distance_um,
    }
    return {
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "original_sample_count": int(len(samples)),
        "excluded_local_sample_count": int(np.count_nonzero(~keep)),
        "core_sample_count": int(len(core)),
        "max_core_surface_distance_um": maximum,
        "P95_core_surface_distance_um": p95,
        "maximum_allowed_um": qc.maximum_core_surface_distance_um,
        "P95_maximum_allowed_um": qc.maximum_core_surface_p95_distance_um,
    }


def extension_collision_qc(
    surface: TaggedSurface, intersections: Iterable[tuple[int, int]]
) -> dict[str, Any]:
    pairs = [
        (int(first), int(second))
        for first, second in intersections
        if surface.face_kind[first] > 0 or surface.face_kind[second] > 0
    ]
    return {
        "status": "PASS" if not pairs else "FAIL",
        "extension_collision_count": len(pairs),
        "intersection_face_pairs": pairs,
    }
