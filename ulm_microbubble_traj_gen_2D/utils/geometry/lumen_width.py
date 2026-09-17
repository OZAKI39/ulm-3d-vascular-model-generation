"""Local lumen diameter and under-resolution metrics from graph and mask data."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage


@dataclass(frozen=True)
class LocalMaskDiameter:
    """Local mask diameter assigned to every lumen cell from the medial axis."""

    diameter_px: np.ndarray
    diameter_um: np.ndarray
    skeleton_mask: np.ndarray
    skeleton_diameter_px: np.ndarray


@dataclass(frozen=True)
class NarrowLumenResult:
    """Narrow-cell classification using graph and mask effective diameters."""

    narrow_mask: np.ndarray
    graph_diameter_px: np.ndarray
    mask_diameter_px: np.ndarray
    effective_diameter_px: np.ndarray
    mask_diameter: LocalMaskDiameter
    min_diameter_px: float


def compute_local_mask_diameter(lumen_mask: np.ndarray, spacing_um: float) -> LocalMaskDiameter:
    """Estimate local mask diameter by propagating medial-axis diameters to lumen cells."""

    lumen = np.asarray(lumen_mask, dtype=bool)
    spacing = max(float(spacing_um), np.finfo(float).eps)
    diameter_px = np.full(lumen.shape, np.nan, dtype=np.float32)
    diameter_um = np.full(lumen.shape, np.nan, dtype=np.float32)
    skeleton_diameter_px = np.full(lumen.shape, np.nan, dtype=np.float32)
    if not np.any(lumen):
        return LocalMaskDiameter(
            diameter_px=diameter_px,
            diameter_um=diameter_um,
            skeleton_mask=np.zeros(lumen.shape, dtype=bool),
            skeleton_diameter_px=skeleton_diameter_px,
        )

    skeleton, distance_cells = _medial_axis_with_distance(lumen)
    if not np.any(skeleton):
        return LocalMaskDiameter(
            diameter_px=diameter_px,
            diameter_um=diameter_um,
            skeleton_mask=skeleton,
            skeleton_diameter_px=skeleton_diameter_px,
        )

    skeleton_diameter_px[skeleton] = (2.0 * np.asarray(distance_cells[skeleton], dtype=float)).astype(np.float32)
    component_labels, component_count = ndimage.label(lumen, structure=np.ones((3, 3), dtype=bool))
    component_slices = ndimage.find_objects(component_labels)
    for label in range(1, component_count + 1):
        bounds = component_slices[label - 1]
        if bounds is None:
            continue
        local_component = component_labels[bounds] == label
        local_skeleton = skeleton[bounds] & local_component
        if not np.any(local_skeleton):
            continue
        _, nearest = ndimage.distance_transform_edt(~local_skeleton, return_indices=True)
        local_distance = np.asarray(distance_cells[bounds][nearest[0], nearest[1]], dtype=float)
        local_diameter = diameter_px[bounds]
        local_diameter[local_component] = (2.0 * local_distance[local_component]).astype(np.float32)
    diameter_um[lumen] = (diameter_px[lumen] * spacing).astype(np.float32)
    return LocalMaskDiameter(
        diameter_px=diameter_px,
        diameter_um=diameter_um,
        skeleton_mask=skeleton,
        skeleton_diameter_px=skeleton_diameter_px,
    )


def classify_narrow_lumen_cells(
    lumen_mask: np.ndarray,
    spacing_um: float,
    radius_um: np.ndarray,
    *,
    min_diameter_px: float = 8.0,
    junction_mask: np.ndarray | None = None,
    open_boundary_mask: np.ndarray | None = None,
) -> NarrowLumenResult:
    """Classify lumen cells whose effective local vessel diameter is under-resolved."""

    lumen = np.asarray(lumen_mask, dtype=bool)
    spacing = max(float(spacing_um), np.finfo(float).eps)
    radius = np.asarray(radius_um, dtype=float)
    if radius.shape != lumen.shape:
        raise ValueError("radius_um must have the same shape as lumen_mask.")

    mask_diameter = compute_local_mask_diameter(lumen, spacing)
    graph_diameter_px = np.full(lumen.shape, np.nan, dtype=np.float32)
    valid_radius = lumen & np.isfinite(radius) & (radius > 0.0)
    graph_diameter_px[valid_radius] = (2.0 * radius[valid_radius] / spacing).astype(np.float32)

    effective = _finite_min(graph_diameter_px, mask_diameter.diameter_px)
    if junction_mask is not None:
        junction = np.asarray(junction_mask, dtype=bool)
        if junction.shape != lumen.shape:
            raise ValueError("junction_mask must have the same shape as lumen_mask.")
        use_mask_only = lumen & junction & np.isfinite(mask_diameter.diameter_px)
        effective[use_mask_only] = mask_diameter.diameter_px[use_mask_only]
    if open_boundary_mask is not None:
        open_boundary = np.asarray(open_boundary_mask, dtype=bool)
        if open_boundary.shape != lumen.shape:
            raise ValueError("open_boundary_mask must have the same shape as lumen_mask.")
        open_boundary = _open_boundary_diagnostic_region(open_boundary, lumen, min_diameter_px)
        use_graph_only = lumen & open_boundary & np.isfinite(graph_diameter_px)
        effective[use_graph_only] = graph_diameter_px[use_graph_only]
    effective[~lumen] = np.nan
    narrow = lumen & np.isfinite(effective) & (effective < float(min_diameter_px))
    return NarrowLumenResult(
        narrow_mask=narrow,
        graph_diameter_px=graph_diameter_px,
        mask_diameter_px=mask_diameter.diameter_px,
        effective_diameter_px=effective.astype(np.float32),
        mask_diameter=mask_diameter,
        min_diameter_px=float(min_diameter_px),
    )


def _finite_min(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    result = np.full(first.shape, np.nan, dtype=np.float32)
    first_finite = np.isfinite(first)
    second_finite = np.isfinite(second)
    both = first_finite & second_finite
    result[both] = np.minimum(first[both], second[both])
    result[first_finite & ~second_finite] = first[first_finite & ~second_finite]
    result[second_finite & ~first_finite] = second[second_finite & ~first_finite]
    return result


def _open_boundary_diagnostic_region(open_boundary_mask: np.ndarray, lumen_mask: np.ndarray, min_diameter_px: float) -> np.ndarray:
    open_boundary = np.asarray(open_boundary_mask, dtype=bool) & np.asarray(lumen_mask, dtype=bool)
    if not np.any(open_boundary):
        return open_boundary
    iterations = max(1, int(np.ceil(0.5 * float(min_diameter_px))))
    structure = ndimage.generate_binary_structure(2, 1)
    return ndimage.binary_dilation(open_boundary, structure=structure, iterations=iterations) & np.asarray(lumen_mask, dtype=bool)


def _medial_axis_with_distance(lumen_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    try:
        from skimage.morphology import medial_axis
    except ImportError as exc:  # pragma: no cover - dependency is available in the project environment
        raise ImportError("scikit-image is required to compute local lumen diameter from a mask.") from exc

    try:
        skeleton, distance_cells = medial_axis(lumen_mask, return_distance=True, rng=0)
    except TypeError:  # pragma: no cover - older scikit-image fallback
        skeleton, distance_cells = medial_axis(lumen_mask, return_distance=True)
    return np.asarray(skeleton, dtype=bool), np.asarray(distance_cells, dtype=float)
