"""Spatially correlated synthetic molecular targets inside a compact influence region."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy import ndimage

from .molecular_target_auto_selection import AutomaticInfluenceAnchorResult
from .molecular_target_candidates import MolecularTargetCandidateCatalog


RANDOM_FIELD_ALGORITHM = "squared_exponential_random_fourier_v1"


@dataclass(frozen=True)
class SpatiallyHeterogeneousTargetResult:
    """One reproducible, wall-area-weighted heterogeneous target realization."""

    anchor_candidate_id: str
    influence_center_xz_um: np.ndarray
    influence_radius_um: float
    requested_influence_wall_area_fraction: float
    achieved_influence_wall_area_fraction: float
    influence_accessible_wall_area_um2: float
    inaccessible_wall_area_inside_influence_um2: float
    requested_positive_wall_fraction_within_influence: float
    achieved_positive_wall_fraction_within_influence: float
    target_positive_wall_area_um2: float
    target_positive_network_wall_area_fraction: float
    correlation_length_um: float
    correlation_length_grid_cells: float
    random_seed: int
    random_field_modes: int
    random_field_algorithm: str
    random_wavevectors_um_inv: np.ndarray
    random_phases_rad: np.ndarray
    random_field_threshold: float
    influence_region_mask: np.ndarray
    influence_wall_mask: np.ndarray
    target_positive_wall_mask: np.ndarray
    influence_wall_flat_indices: np.ndarray
    influence_wall_random_field: np.ndarray
    patch_count: int


def build_spatially_heterogeneous_target(
    catalog: MolecularTargetCandidateCatalog,
    anchor: AutomaticInfluenceAnchorResult,
    *,
    influence_wall_area_fraction: float,
    positive_wall_fraction_within_influence: float,
    correlation_length_um: float,
    random_seed: int,
    random_field_modes: int,
) -> SpatiallyHeterogeneousTargetResult:
    """Generate a compact tissue-space region with correlated positive wall patches."""

    _validate_inputs(
        catalog,
        influence_wall_area_fraction,
        positive_wall_fraction_within_influence,
        correlation_length_um,
        random_seed,
        random_field_modes,
    )
    weights = np.asarray(catalog.wall_area_weight_um2, dtype=np.float64)
    accessible_wall = (
        np.asarray(catalog.accessible_wall_mask, dtype=bool)
        & np.asarray(catalog.solid_wall_mask, dtype=bool)
        & (weights > 0.0)
    )
    accessible_indices = np.argwhere(accessible_wall)
    accessible_flat = np.ravel_multi_index(
        (accessible_indices[:, 0], accessible_indices[:, 1]),
        catalog.shape,
    )
    accessible_coordinates = np.column_stack(
        (
            catalog.x_coordinates_um[accessible_indices[:, 0]],
            catalog.z_coordinates_um[accessible_indices[:, 1]],
        )
    )
    center = np.asarray([anchor.anchor_x_um, anchor.anchor_z_um], dtype=np.float64)
    distances = np.linalg.norm(accessible_coordinates - center[None, :], axis=1)
    accessible_weights = weights[accessible_wall]
    target_influence_area = (
        float(influence_wall_area_fraction)
        * float(catalog.network_endothelial_wall_area_um2)
    )
    available_area = float(np.sum(accessible_weights))
    if target_influence_area > available_area:
        raise ValueError(
            "The requested influence-region wall area exceeds the wall area expected to receive "
            "at least one bubble during the configured observation time."
        )
    radius = _area_matched_radius(
        distances,
        accessible_weights,
        accessible_flat,
        target_influence_area,
    )

    x_grid, z_grid = np.meshgrid(
        catalog.x_coordinates_um,
        catalog.z_coordinates_um,
        indexing="ij",
    )
    distance_squared = (x_grid - center[0]) ** 2 + (z_grid - center[1]) ** 2
    radius_tolerance = max(radius * radius, 1.0) * 32.0 * np.finfo(float).eps
    influence_region = distance_squared <= radius * radius + radius_tolerance
    influence_wall = influence_region & accessible_wall
    influence_area = float(np.sum(weights[influence_wall]))
    if influence_area <= 0.0:
        raise ValueError("The automatic influence region contains no accessible wall area.")

    influence_indices = np.argwhere(influence_wall)
    influence_flat = np.ravel_multi_index(
        (influence_indices[:, 0], influence_indices[:, 1]),
        catalog.shape,
    )
    influence_coordinates = np.column_stack(
        (
            catalog.x_coordinates_um[influence_indices[:, 0]],
            catalog.z_coordinates_um[influence_indices[:, 1]],
        )
    )
    wavevectors, phases = generate_random_fourier_coefficients(
        correlation_length_um,
        random_seed,
        random_field_modes,
    )
    field_values = evaluate_random_fourier_field(
        influence_coordinates,
        center,
        wavevectors,
        phases,
    )
    selected_local, threshold = _weighted_positive_sites(
        field_values,
        weights[influence_wall],
        influence_flat,
        float(positive_wall_fraction_within_influence),
    )
    target_wall = np.zeros(catalog.shape, dtype=bool)
    target_indices = influence_indices[selected_local]
    target_wall[target_indices[:, 0], target_indices[:, 1]] = True
    target_area = float(np.sum(weights[target_wall]))
    achieved_positive_fraction = target_area / influence_area

    closed_wall = (
        np.asarray(catalog.solid_wall_mask, dtype=bool)
        & ~np.asarray(catalog.open_boundary_mask, dtype=bool)
        & (weights > 0.0)
    )
    inaccessible_inside = influence_region & closed_wall & ~accessible_wall
    inaccessible_area = float(np.sum(weights[inaccessible_inside]))
    spacing_um = _grid_spacing_um(catalog)
    _, patch_count = ndimage.label(
        target_wall,
        structure=np.ones((3, 3), dtype=np.uint8),
    )
    return SpatiallyHeterogeneousTargetResult(
        anchor_candidate_id=anchor.anchor_candidate_id,
        influence_center_xz_um=np.ascontiguousarray(center, dtype=np.float64),
        influence_radius_um=float(radius),
        requested_influence_wall_area_fraction=float(influence_wall_area_fraction),
        achieved_influence_wall_area_fraction=(
            influence_area / float(catalog.network_endothelial_wall_area_um2)
        ),
        influence_accessible_wall_area_um2=influence_area,
        inaccessible_wall_area_inside_influence_um2=inaccessible_area,
        requested_positive_wall_fraction_within_influence=float(
            positive_wall_fraction_within_influence
        ),
        achieved_positive_wall_fraction_within_influence=float(
            achieved_positive_fraction
        ),
        target_positive_wall_area_um2=target_area,
        target_positive_network_wall_area_fraction=(
            target_area / float(catalog.network_endothelial_wall_area_um2)
        ),
        correlation_length_um=float(correlation_length_um),
        correlation_length_grid_cells=float(correlation_length_um) / spacing_um,
        random_seed=int(random_seed),
        random_field_modes=int(random_field_modes),
        random_field_algorithm=RANDOM_FIELD_ALGORITHM,
        random_wavevectors_um_inv=np.ascontiguousarray(wavevectors, dtype=np.float64),
        random_phases_rad=np.ascontiguousarray(phases, dtype=np.float64),
        random_field_threshold=float(threshold),
        influence_region_mask=np.ascontiguousarray(influence_region, dtype=bool),
        influence_wall_mask=np.ascontiguousarray(influence_wall, dtype=bool),
        target_positive_wall_mask=np.ascontiguousarray(target_wall, dtype=bool),
        influence_wall_flat_indices=np.ascontiguousarray(influence_flat, dtype=np.int64),
        influence_wall_random_field=np.ascontiguousarray(field_values, dtype=np.float64),
        patch_count=int(patch_count),
    )


def generate_random_fourier_coefficients(
    correlation_length_um: float,
    random_seed: int,
    random_field_modes: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Draw the reusable spectral coefficients of a squared-exponential field."""

    rng = np.random.Generator(np.random.PCG64(int(random_seed)))
    wavevectors = rng.normal(
        loc=0.0,
        scale=1.0 / float(correlation_length_um),
        size=(int(random_field_modes), 2),
    )
    phases = rng.uniform(0.0, 2.0 * np.pi, size=int(random_field_modes))
    return wavevectors, phases


def evaluate_random_fourier_field(
    coordinates_xz_um: np.ndarray,
    coordinate_origin_xz_um: np.ndarray,
    wavevectors_um_inv: np.ndarray,
    phases_rad: np.ndarray,
    *,
    chunk_size: int = 4096,
) -> np.ndarray:
    """Evaluate one continuous realization at physical coordinates in bounded memory."""

    coordinates = np.asarray(coordinates_xz_um, dtype=np.float64)
    origin = np.asarray(coordinate_origin_xz_um, dtype=np.float64)
    wavevectors = np.asarray(wavevectors_um_inv, dtype=np.float64)
    phases = np.asarray(phases_rad, dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 2:
        raise ValueError("Random-field coordinates must have shape (N, 2).")
    if origin.shape != (2,):
        raise ValueError("Random-field coordinate origin must have shape (2,).")
    if wavevectors.ndim != 2 or wavevectors.shape[1] != 2:
        raise ValueError("Random-field wavevectors must have shape (M, 2).")
    if phases.shape != (wavevectors.shape[0],):
        raise ValueError("Random-field phases must contain one value per wavevector.")
    values = np.empty(coordinates.shape[0], dtype=np.float64)
    scale = math.sqrt(2.0 / wavevectors.shape[0])
    relative = coordinates - origin[None, :]
    for start in range(0, coordinates.shape[0], int(chunk_size)):
        stop = min(start + int(chunk_size), coordinates.shape[0])
        angles = relative[start:stop] @ wavevectors.T
        angles += phases[None, :]
        values[start:stop] = scale * np.sum(np.cos(angles), axis=1)
    return values


def _area_matched_radius(
    distances_um: np.ndarray,
    wall_area_weights_um2: np.ndarray,
    flat_indices: np.ndarray,
    requested_area_um2: float,
) -> float:
    order = np.lexsort((flat_indices, distances_um))
    sorted_distances = np.asarray(distances_um, dtype=np.float64)[order]
    sorted_weights = np.asarray(wall_area_weights_um2, dtype=np.float64)[order]
    unique_distances, first_indices = np.unique(sorted_distances, return_index=True)
    group_areas = np.add.reduceat(sorted_weights, first_indices)
    cumulative = np.cumsum(group_areas)
    selected = int(np.argmin(np.abs(cumulative - float(requested_area_um2))))
    return float(unique_distances[selected])


def _weighted_positive_sites(
    field_values: np.ndarray,
    wall_area_weights_um2: np.ndarray,
    flat_indices: np.ndarray,
    positive_fraction: float,
) -> tuple[np.ndarray, float]:
    values = np.asarray(field_values, dtype=np.float64)
    weights = np.asarray(wall_area_weights_um2, dtype=np.float64)
    if positive_fraction == 1.0:
        return np.ones(values.size, dtype=bool), float("-inf")
    order = np.lexsort((np.asarray(flat_indices, dtype=np.int64), -values))
    cumulative = np.cumsum(weights[order])
    requested_area = float(positive_fraction) * float(np.sum(weights))
    count = int(np.argmin(np.abs(cumulative - requested_area))) + 1
    selected = np.zeros(values.size, dtype=bool)
    selected[order[:count]] = True
    threshold = float(np.min(values[selected]))
    return selected, threshold


def _grid_spacing_um(catalog: MolecularTargetCandidateCatalog) -> float:
    dx = np.diff(np.asarray(catalog.x_coordinates_um, dtype=np.float64))
    dz = np.diff(np.asarray(catalog.z_coordinates_um, dtype=np.float64))
    spacing = float(np.mean(np.concatenate((dx, dz))))
    if not math.isfinite(spacing) or spacing <= 0.0:
        raise ValueError("Candidate physical axes must have positive grid spacing.")
    return spacing


def _validate_inputs(
    catalog: MolecularTargetCandidateCatalog,
    influence_fraction: float,
    positive_fraction: float,
    correlation_length_um: float,
    random_seed: int,
    random_field_modes: int,
) -> None:
    if not catalog.automatic_metrics_available:
        raise ValueError(
            "Revised-v10 automatic target generation requires a v3 candidate catalog."
        )
    if not math.isfinite(influence_fraction) or not 0.0 < influence_fraction < 1.0:
        raise ValueError("Influence-region wall-area fraction must be in (0, 1).")
    if not math.isfinite(positive_fraction) or not 0.0 < positive_fraction <= 1.0:
        raise ValueError("Target-positive wall fraction must be in (0, 1].")
    if not math.isfinite(correlation_length_um) or correlation_length_um <= 0.0:
        raise ValueError("Target correlation length must be finite and positive.")
    spacing_um = _grid_spacing_um(catalog)
    if correlation_length_um < spacing_um:
        raise ValueError(
            "Target correlation length must span at least one physical grid cell."
        )
    if isinstance(random_seed, bool) or int(random_seed) != random_seed or random_seed < 0:
        raise ValueError("Random seed must be a non-negative integer.")
    if (
        isinstance(random_field_modes, bool)
        or int(random_field_modes) != random_field_modes
        or random_field_modes < 64
    ):
        raise ValueError("Random-field mode count must be an integer of at least 64.")
    if not math.isfinite(float(catalog.mapped_endothelial_wall_area_um2)) or (
        catalog.mapped_endothelial_wall_area_um2 <= 0.0
    ):
        raise ValueError(
            "The candidate grid contains no physically weighted endothelial wall samples."
        )
