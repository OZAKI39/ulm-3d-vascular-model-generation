"""Measured local mesh targets and extension-only triangle quality controls."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pyvista as pv
import trimesh

from .config import MeshQualityConfig
from .io import BoundaryInput, SurfacePrepareError
from .types import BoundarySurfaceResult, TaggedSurface


ASPECT_RATIO_DEFINITION = (
    "(sqrt(3)/2) * longest_edge / altitude_to_longest_edge; equilateral=1"
)


@dataclass(frozen=True, slots=True)
class TriangleMetrics:
    edge_lengths: np.ndarray
    areas: np.ndarray
    minimum_angles_deg: np.ndarray
    maximum_angles_deg: np.ndarray
    aspect_ratios: np.ndarray


def triangle_metrics(vertices: np.ndarray, faces: np.ndarray) -> TriangleMetrics:
    triangles = np.asarray(vertices, dtype=float)[np.asarray(faces, dtype=np.int64)]
    edge_vectors = np.stack(
        (
            triangles[:, 1] - triangles[:, 0],
            triangles[:, 2] - triangles[:, 1],
            triangles[:, 0] - triangles[:, 2],
        ),
        axis=1,
    )
    lengths = np.linalg.norm(edge_vectors, axis=2)
    area_twice = np.linalg.norm(
        np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0]),
        axis=1,
    )
    if np.any(area_twice <= np.finfo(float).eps):
        raise SurfacePrepareError("EXTENSION_MESH_HAS_DEGENERATE_TRIANGLE")
    cosines = []
    for center, first, second in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        first_vector = triangles[:, first] - triangles[:, center]
        second_vector = triangles[:, second] - triangles[:, center]
        cosines.append(
            np.sum(first_vector * second_vector, axis=1)
            / (
                np.linalg.norm(first_vector, axis=1)
                * np.linalg.norm(second_vector, axis=1)
            )
        )
    angles = np.degrees(np.arccos(np.clip(np.stack(cosines, axis=1), -1.0, 1.0)))
    return TriangleMetrics(
        edge_lengths=lengths,
        areas=0.5 * area_twice,
        minimum_angles_deg=np.min(angles, axis=1),
        maximum_angles_deg=np.max(angles, axis=1),
        aspect_ratios=(np.sqrt(3.0) / 2.0)
        * np.max(lengths, axis=1) ** 2
        / area_twice,
    )


def measure_local_original_mesh(
    original: trimesh.Trimesh,
    boundary: BoundaryInput,
    *,
    sampling_radius_factor: float,
) -> dict[str, Any]:
    """Measure inward local wall triangles and exclude the outward rounded end."""

    vertices = np.asarray(original.vertices, dtype=float)
    faces = np.asarray(original.faces, dtype=np.int64)
    centers = vertices[faces].mean(axis=1)
    relative = centers - boundary.center_um
    axial = relative @ boundary.outward_normal
    mask = (
        np.linalg.norm(relative, axis=1)
        <= sampling_radius_factor * boundary.source_radius_um
    ) & (axial <= 0.0)
    if np.count_nonzero(mask) < 10:
        raise SurfacePrepareError(f"LOCAL_ORIGINAL_MESH_SAMPLE_INSUFFICIENT:{boundary.port_id}")
    metrics = triangle_metrics(vertices, faces[mask])
    edges = metrics.edge_lengths.ravel()
    return {
        "port_id": boundary.port_id,
        "boundary_origin": boundary.boundary_origin,
        "role": boundary.role,
        "sampling_radius_factor": sampling_radius_factor,
        "rounded_outward_cap_excluded": True,
        "triangle_count": int(np.count_nonzero(mask)),
        "edge_length_min_um": float(np.min(edges)),
        "edge_length_p25_um": float(np.percentile(edges, 25)),
        "edge_length_median_um": float(np.median(edges)),
        "edge_length_p75_um": float(np.percentile(edges, 75)),
        "edge_length_p95_um": float(np.percentile(edges, 95)),
        "triangle_area_median_um2": float(np.median(metrics.areas)),
        "triangle_area_p95_um2": float(np.percentile(metrics.areas, 95)),
        "aspect_ratio_definition": ASPECT_RATIO_DEFINITION,
        "aspect_ratio_median": float(np.median(metrics.aspect_ratios)),
        "aspect_ratio_p95": float(np.percentile(metrics.aspect_ratios, 95)),
    }


def _neighbor_area_ratios(faces: np.ndarray, areas: np.ndarray) -> np.ndarray:
    edge_faces: dict[tuple[int, int], list[int]] = {}
    for local_id, face in enumerate(faces):
        for first, second in zip(face, np.roll(face, -1)):
            edge_faces.setdefault(tuple(sorted((int(first), int(second)))), []).append(local_id)
    ratios = np.ones(len(faces), dtype=float)
    for linked in edge_faces.values():
        if len(linked) != 2:
            continue
        first, second = linked
        ratio = max(areas[first], areas[second]) / min(areas[first], areas[second])
        ratios[first] = max(ratios[first], ratio)
        ratios[second] = max(ratios[second], ratio)
    return ratios


def summarize_extension_mesh(
    vertices: np.ndarray,
    faces: np.ndarray,
    *,
    target_edge_length_um: float,
    local_original_median_edge_length_um: float,
    quality: MeshQualityConfig | None,
    quality_face_mask: np.ndarray | None = None,
) -> dict[str, Any]:
    metrics = triangle_metrics(vertices, faces)
    edges = metrics.edge_lengths.ravel()
    neighbor_ratios = _neighbor_area_ratios(faces, metrics.areas)
    maximum_edge_ratio = np.max(metrics.edge_lengths, axis=1) / target_edge_length_um
    if quality is None:
        bad = np.zeros(len(faces), dtype=bool)
    else:
        bad = (
            (
                (metrics.minimum_angles_deg < quality.minimum_triangle_angle_deg)
                & (metrics.aspect_ratios > quality.maximum_aspect_ratio)
            )
            | (maximum_edge_ratio > quality.maximum_edge_length_to_local_target_ratio)
        )
    aggregate_mask = (
        np.ones(len(faces), dtype=bool)
        if quality_face_mask is None
        else np.asarray(quality_face_mask, dtype=bool)
    )
    aggregate_bad = bad[aggregate_mask]
    return {
        "extension_triangle_count": int(len(faces)),
        "edge_length_min_um": float(np.min(edges)),
        "edge_length_median_um": float(np.median(edges)),
        "edge_length_p95_um": float(np.percentile(edges, 95)),
        "edge_length_max_um": float(np.max(edges)),
        "minimum_triangle_angle_deg": float(np.min(metrics.minimum_angles_deg)),
        "triangle_angle_p05_deg": float(np.percentile(metrics.minimum_angles_deg, 5)),
        "maximum_triangle_angle_deg": float(np.max(metrics.maximum_angles_deg)),
        "aspect_ratio_definition": ASPECT_RATIO_DEFINITION,
        "aspect_ratio_median": float(np.median(metrics.aspect_ratios)),
        "aspect_ratio_p95": float(np.percentile(metrics.aspect_ratios, 95)),
        "aspect_ratio_max": float(np.max(metrics.aspect_ratios)),
        "neighbor_area_ratio_p95": float(np.percentile(neighbor_ratios, 95)),
        "neighbor_area_ratio_max": float(np.max(neighbor_ratios)),
        "maximum_edge_length_to_local_target_ratio": float(np.max(maximum_edge_ratio)),
        "mesh_size_transition_ratio": float(
            np.median(edges) / local_original_median_edge_length_um
        ),
        "all_extension_bad_triangle_count": int(np.count_nonzero(bad)),
        "all_extension_bad_triangle_fraction": float(np.mean(bad)),
        "bad_triangle_count": int(np.count_nonzero(aggregate_bad)),
        "bad_triangle_fraction": float(np.mean(aggregate_bad)),
    }


def extension_mesh_quality_qc(
    surface: TaggedSurface,
    boundaries: Iterable[BoundaryInput],
    results: Iterable[BoundarySurfaceResult],
    quality: MeshQualityConfig,
) -> dict[str, Any]:
    result_by_index = {item.boundary_index: item for item in results}
    records: list[dict[str, Any]] = []
    for boundary in boundaries:
        result = result_by_index[boundary.index]
        mask = (surface.face_kind == 1) & (surface.extension_index == boundary.index)
        faces = surface.faces[mask]
        summary = summarize_extension_mesh(
            surface.vertices,
            faces,
            target_edge_length_um=result.target_edge_length_um,
            local_original_median_edge_length_um=result.local_original_median_edge_length_um,
            quality=quality,
            quality_face_mask=surface.extension_band[mask] > 0,
        )
        interface_mask = mask & (surface.extension_band == 0)
        interface_faces = surface.faces[interface_mask]
        interface = triangle_metrics(surface.vertices, interface_faces)
        interface_edge_ratio = float(
            np.median(interface.edge_lengths) / result.target_edge_length_um
        )
        finite = bool(
            all(
                np.isfinite(value)
                for value in summary.values()
                if isinstance(value, float)
            )
        )
        checks = {
            "finite_metrics": finite,
            "intermediate_rings_present": result.ring_count > 2,
            "bad_triangle_fraction": summary["bad_triangle_fraction"]
            <= quality.maximum_bad_triangle_fraction,
            "interface_edge_length_ratio": interface_edge_ratio
            <= quality.maximum_interface_edge_length_ratio,
            "neighbor_area_ratio_p95": summary["neighbor_area_ratio_p95"]
            <= quality.maximum_neighbor_area_ratio,
        }
        records.append(
            {
                "port_id": boundary.port_id,
                "boundary_origin": boundary.boundary_origin,
                "role": boundary.role,
                "local_original_median_edge_length_um": result.local_original_median_edge_length_um,
                "target_edge_length_um": result.target_edge_length_um,
                "ring_count": result.ring_count,
                **summary,
                "interface_triangle_count": int(len(interface_faces)),
                "interface_edge_length_ratio": interface_edge_ratio,
                "interface_triangle_aspect_ratio_p95": float(
                    np.percentile(interface.aspect_ratios, 95)
                ),
                "interface_minimum_angle_deg": float(
                    np.min(interface.minimum_angles_deg)
                ),
                "quality_aggregate_excludes_separately_checked_interface_band": True,
                "checks": checks,
                "status": "PASS" if all(checks.values()) else "FAIL",
            }
        )
    return {
        "status": "PASS"
        if all(record["status"] == "PASS" for record in records)
        else "FAIL",
        "aspect_ratio_definition": ASPECT_RATIO_DEFINITION,
        "thresholds": {
            "minimum_triangle_angle_deg": quality.minimum_triangle_angle_deg,
            "maximum_aspect_ratio": quality.maximum_aspect_ratio,
            "maximum_edge_length_to_local_target_ratio": quality.maximum_edge_length_to_local_target_ratio,
            "maximum_neighbor_area_ratio": quality.maximum_neighbor_area_ratio,
            "maximum_interface_edge_length_ratio": quality.maximum_interface_edge_length_ratio,
            "maximum_bad_triangle_fraction": quality.maximum_bad_triangle_fraction,
        },
        "boundaries": records,
    }


def _read_previous_vtp(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    data = pv.read(path).triangulate()
    faces = np.asarray(data.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    boundary_type = np.asarray(data.cell_data["boundary_type_code"], dtype=np.uint8)
    return np.asarray(data.points, dtype=float), faces, boundary_type


def compare_previous_extensions(
    previous_vtp: Path,
    refined: TaggedSurface,
    boundaries: Iterable[BoundaryInput],
    results: Iterable[BoundarySurfaceResult],
) -> list[dict[str, Any]]:
    old_vertices, old_faces, old_type = _read_previous_vtp(previous_vtp)
    old_centers = old_vertices[old_faces].mean(axis=1)
    result_by_index = {item.boundary_index: item for item in results}
    rows: list[dict[str, Any]] = []
    for boundary in boundaries:
        result = result_by_index[boundary.index]
        relative = old_centers - boundary.center_um
        axial = relative @ boundary.outward_normal
        radial = np.linalg.norm(
            relative - np.outer(axial, boundary.outward_normal), axis=1
        )
        old_mask = (
            (old_type == 0)
            & (axial > 1.0e-8)
            & (axial < boundary.extension_length_um + 1.0e-8)
            & (radial < 2.0 * boundary.source_radius_um)
        )
        new_mask = (refined.face_kind == 1) & (
            refined.extension_index == boundary.index
        )
        old_summary = summarize_extension_mesh(
            old_vertices,
            old_faces[old_mask],
            target_edge_length_um=result.target_edge_length_um,
            local_original_median_edge_length_um=result.local_original_median_edge_length_um,
            quality=None,
        )
        new_summary = summarize_extension_mesh(
            refined.vertices,
            refined.faces[new_mask],
            target_edge_length_um=result.target_edge_length_um,
            local_original_median_edge_length_um=result.local_original_median_edge_length_um,
            quality=None,
        )
        rows.append(
            {
                "port_id": boundary.port_id,
                "previous_triangle_count": old_summary["extension_triangle_count"],
                "refined_triangle_count": new_summary["extension_triangle_count"],
                "previous_median_edge_um": old_summary["edge_length_median_um"],
                "refined_median_edge_um": new_summary["edge_length_median_um"],
                "previous_aspect_ratio_p95": old_summary["aspect_ratio_p95"],
                "refined_aspect_ratio_p95": new_summary["aspect_ratio_p95"],
                "previous_aspect_ratio_max": old_summary["aspect_ratio_max"],
                "refined_aspect_ratio_max": new_summary["aspect_ratio_max"],
                "previous_minimum_angle_deg": old_summary["minimum_triangle_angle_deg"],
                "refined_minimum_angle_deg": new_summary["minimum_triangle_angle_deg"],
                "triangle_count_increased": new_summary["extension_triangle_count"]
                > old_summary["extension_triangle_count"],
                "p95_aspect_ratio_improved": new_summary["aspect_ratio_p95"]
                < old_summary["aspect_ratio_p95"],
                "minimum_angle_improved": new_summary["minimum_triangle_angle_deg"]
                > old_summary["minimum_triangle_angle_deg"],
            }
        )
    return rows
