"""Measured local targets and triangle metrics used by the formal VMTK path."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import trimesh

from .config import MeshQualityConfig
from .io import BoundaryInput, SurfacePrepareError


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
        np.cross(
            triangles[:, 1] - triangles[:, 0],
            triangles[:, 2] - triangles[:, 0],
        ),
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
    """Measure inward local wall triangles and exclude the rounded vessel end."""

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
        raise SurfacePrepareError(
            f"LOCAL_ORIGINAL_MESH_SAMPLE_INSUFFICIENT:{boundary.port_id}"
        )
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
            edge_faces.setdefault(
                tuple(sorted((int(first), int(second)))), []
            ).append(local_id)
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
