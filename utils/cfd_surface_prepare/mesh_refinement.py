"""Structured, extension-only ring refinement with locked end geometry."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .config import ExtensionMeshConfig
from .io import SurfacePrepareError
from .local_cut import orthogonal_basis


@dataclass(frozen=True, slots=True)
class RefinedRings:
    points: np.ndarray
    stations_um: np.ndarray
    target_areas_um2: np.ndarray
    regularized_loop: np.ndarray
    target_edge_length_um: float
    transition_length_um: float

    @property
    def ring_count(self) -> int:
        return int(len(self.points))


def polygon_area_centroid(points_2d: np.ndarray) -> tuple[float, np.ndarray]:
    following = np.roll(points_2d, -1, axis=0)
    cross = points_2d[:, 0] * following[:, 1] - following[:, 0] * points_2d[:, 1]
    signed_area = 0.5 * float(np.sum(cross))
    if not np.isfinite(signed_area) or abs(signed_area) <= 1.0e-15:
        raise SurfacePrepareError("EXTENSION_LOOP_AREA_INVALID")
    centroid = np.asarray(
        (
            np.sum((points_2d[:, 0] + following[:, 0]) * cross),
            np.sum((points_2d[:, 1] + following[:, 1]) * cross),
        ),
        dtype=float,
    ) / (6.0 * signed_area)
    return abs(signed_area), centroid


def _restore_area_centroid(
    points_2d: np.ndarray, target_area: float, target_centroid: np.ndarray
) -> np.ndarray:
    area, centroid = polygon_area_centroid(points_2d)
    scale = np.sqrt(target_area / area)
    return (points_2d - centroid) * scale + target_centroid


def regularize_loop(
    points: np.ndarray,
    normal: np.ndarray,
    *,
    iterations: int,
    relaxation: float,
) -> np.ndarray:
    """Cyclically smooth a planar loop while preserving area, centroid, and plane."""

    first, second = orthogonal_basis(normal)
    plane_origin = np.mean(points, axis=0)
    points_2d = np.column_stack(
        ((points - plane_origin) @ first, (points - plane_origin) @ second)
    )
    target_area, target_centroid = polygon_area_centroid(points_2d)
    current = points_2d.copy()
    for _ in range(iterations):
        neighbor_average = 0.5 * (
            np.roll(current, 1, axis=0) + np.roll(current, -1, axis=0)
        )
        current += relaxation * (neighbor_average - current)
        current = _restore_area_centroid(current, target_area, target_centroid)
    closed = np.vstack((current, current[0]))
    segment_lengths = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
    equal_stations = np.linspace(0.0, cumulative[-1], len(current), endpoint=False)
    equalized: list[np.ndarray] = []
    for station in equal_stations:
        segment = min(
            int(np.searchsorted(cumulative, station, side="right") - 1),
            len(current) - 1,
        )
        fraction = (station - cumulative[segment]) / segment_lengths[segment]
        equalized.append(
            closed[segment] + fraction * (closed[segment + 1] - closed[segment])
        )
    current = _restore_area_centroid(
        np.asarray(equalized), target_area, target_centroid
    )
    result = plane_origin + np.outer(current[:, 0], first) + np.outer(
        current[:, 1], second
    )
    residual = np.max(np.abs((result - plane_origin) @ normal))
    if residual > 1.0e-10:
        raise SurfacePrepareError("EXTENSION_LOOP_REGULARIZATION_LEFT_PLANE")
    return result


def calculate_ring_count(
    extension_length_um: float,
    target_edge_length_um: float,
    *,
    minimum: int,
    maximum: int,
) -> int:
    raw = int(np.ceil(extension_length_um / target_edge_length_um))
    return int(np.clip(raw, minimum, maximum))


def _smoothstep(value: float) -> float:
    clipped = float(np.clip(value, 0.0, 1.0))
    return clipped * clipped * (3.0 - 2.0 * clipped)


def _ring_target(
    actual: np.ndarray,
    regularized: np.ndarray,
    normal: np.ndarray,
    station_um: float,
    transition_length_um: float,
) -> np.ndarray:
    alpha = _smoothstep(station_um / transition_length_um)
    shape = actual + alpha * (regularized - actual)
    return shape + station_um * normal


def _section_center_3d(
    points: np.ndarray, first: np.ndarray, second: np.ndarray
) -> np.ndarray:
    origin = np.mean(points, axis=0)
    points_2d = np.column_stack(
        ((points - origin) @ first, (points - origin) @ second)
    )
    _, center_2d = polygon_area_centroid(points_2d)
    return origin + center_2d[0] * first + center_2d[1] * second


def _constrained_smoothing_pass(
    rings: np.ndarray,
    stations: np.ndarray,
    target_areas: np.ndarray,
    target_centroids: np.ndarray,
    normal: np.ndarray,
    factor: float,
    *,
    transition_only: bool,
    transition_length_um: float,
) -> np.ndarray:
    result = rings.copy()
    first, second = orthogonal_basis(normal)
    for ring_id in range(1, len(rings) - 1):
        if transition_only and stations[ring_id] > transition_length_um:
            continue
        current = rings[ring_id]
        average = 0.25 * (
            np.roll(current, 1, axis=0)
            + np.roll(current, -1, axis=0)
            + rings[ring_id - 1]
            + rings[ring_id + 1]
        )
        proposed = current + factor * (average - current)
        relative = proposed - target_centroids[ring_id]
        proposed_2d = np.column_stack((relative @ first, relative @ second))
        proposed_2d = _restore_area_centroid(
            proposed_2d, target_areas[ring_id], np.zeros(2, dtype=float)
        )
        result[ring_id] = (
            target_centroids[ring_id]
            + np.outer(proposed_2d[:, 0], first)
            + np.outer(proposed_2d[:, 1], second)
        )
    return result


def build_refined_rings(
    proximal: np.ndarray,
    normal: np.ndarray,
    extension_length_um: float,
    local_original_median_edge_length_um: float,
    source_radius_um: float,
    config: ExtensionMeshConfig,
) -> RefinedRings:
    """Build multiple axial rings; Ring 0 and the distal ring remain exact and locked."""

    target = (
        local_original_median_edge_length_um
        * config.refinement.target_axial_spacing_factor
    )
    if not np.isfinite(target) or target <= 0:
        raise SurfacePrepareError("LOCAL_ORIGINAL_MESH_TARGET_INVALID")
    ring_count = calculate_ring_count(
        extension_length_um,
        target,
        minimum=config.refinement.minimum_ring_count,
        maximum=config.refinement.maximum_ring_count,
    )
    if ring_count <= 2:
        raise SurfacePrepareError("EXTENSION_INTERMEDIATE_RINGS_REQUIRED")
    regularized = regularize_loop(
        proximal,
        normal,
        iterations=config.transition.loop_regularization.iterations,
        relaxation=config.transition.loop_regularization.relaxation,
    )
    transition_length = min(
        extension_length_um,
        config.transition.transition_length_diameters * 2.0 * source_radius_um,
    )
    stations = np.linspace(0.0, extension_length_um, ring_count)
    first, second = orthogonal_basis(normal)
    proximal_center = _section_center_3d(proximal, first, second)
    target_rings: list[np.ndarray] = []
    for station in stations:
        target_ring = _ring_target(
            proximal,
            regularized,
            normal,
            float(station),
            transition_length,
        )
        expected_center = proximal_center + float(station) * normal
        target_ring += expected_center - _section_center_3d(
            target_ring, first, second
        )
        target_rings.append(target_ring)
    targets = np.asarray(target_rings, dtype=float)
    target_areas = np.empty(ring_count, dtype=float)
    target_centroids = np.empty((ring_count, 3), dtype=float)
    for ring_id, ring in enumerate(targets):
        origin = np.mean(ring, axis=0)
        points_2d = np.column_stack(
            ((ring - origin) @ first, (ring - origin) @ second)
        )
        area, centroid_2d = polygon_area_centroid(points_2d)
        target_areas[ring_id] = area
        target_centroids[ring_id] = (
            origin + centroid_2d[0] * first + centroid_2d[1] * second
        )
    rings = targets.copy()
    for _ in range(config.smoothing.iterations):
        rings = _constrained_smoothing_pass(
            rings,
            stations,
            target_areas,
            target_centroids,
            normal,
            config.smoothing.lambda_factor,
            transition_only=config.smoothing.transition_only,
            transition_length_um=transition_length,
        )
        rings = _constrained_smoothing_pass(
            rings,
            stations,
            target_areas,
            target_centroids,
            normal,
            config.smoothing.mu_factor,
            transition_only=config.smoothing.transition_only,
            transition_length_um=transition_length,
        )
    rings[0] = proximal
    rings[-1] = regularized + extension_length_um * normal
    return RefinedRings(
        points=rings,
        stations_um=stations,
        target_areas_um2=target_areas,
        regularized_loop=regularized,
        target_edge_length_um=target,
        transition_length_um=transition_length,
    )


def _triangle_quality(points: np.ndarray) -> tuple[float, float]:
    edges = np.asarray(
        (
            np.linalg.norm(points[1] - points[0]),
            np.linalg.norm(points[2] - points[1]),
            np.linalg.norm(points[0] - points[2]),
        )
    )
    area_twice = np.linalg.norm(
        np.cross(points[1] - points[0], points[2] - points[0])
    )
    if area_twice <= np.finfo(float).eps:
        return float("inf"), 0.0
    aspect = float(np.max(edges) ** 2 / area_twice)
    cosines = np.asarray(
        [
            np.dot(points[(index + 1) % 3] - points[index], points[(index + 2) % 3] - points[index])
            / (
                np.linalg.norm(points[(index + 1) % 3] - points[index])
                * np.linalg.norm(points[(index + 2) % 3] - points[index])
            )
            for index in range(3)
        ]
    )
    minimum_angle = float(np.min(np.degrees(np.arccos(np.clip(cosines, -1.0, 1.0)))))
    return aspect, minimum_angle


def choose_quad_diagonal(
    vertices: np.ndarray, first_ring: np.ndarray, second_ring: np.ndarray, index: int
) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    """Choose the diagonal with lower worst aspect ratio, then higher minimum angle."""

    following = (index + 1) % len(first_ring)
    a = int(first_ring[index])
    a_next = int(first_ring[following])
    b = int(second_ring[index])
    b_next = int(second_ring[following])
    options = (
        ((a, a_next, b_next), (a, b_next, b)),
        ((a, a_next, b), (a_next, b_next, b)),
    )
    scores: list[tuple[float, float]] = []
    for option in options:
        qualities = [_triangle_quality(vertices[np.asarray(face)]) for face in option]
        scores.append(
            (
                max(item[0] for item in qualities),
                -min(item[1] for item in qualities),
            )
        )
    return options[int(scores[1] < scores[0])]
