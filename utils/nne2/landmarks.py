"""Match NNE2 reference-image measurement landmarks to a 3-D branch graph."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage
from scipy.spatial import cKDTree
from skimage.feature import match_template
from skimage.transform import resize

from ..graph.model import HierarchicalGraphResult
from .catalog import NNE2Record
from .config import NNE2Config
from .model import AnchorMatch


def _unit_contrast(image: np.ndarray) -> np.ndarray:
    image = np.asarray(image, dtype=np.float32)
    low, high = np.percentile(image, (5.0, 99.5))
    if high <= low:
        return np.zeros_like(image, dtype=np.float32)
    output = np.clip((image - low) / (high - low), 0.0, 1.0)
    background = ndimage.gaussian_filter(output, sigma=7.0, mode="nearest")
    feature = output - background
    std = float(np.std(feature))
    return np.asarray(feature / max(std, 1.0e-6), dtype=np.float32)


def _reference_green_and_anchor(path: Path) -> tuple[np.ndarray, tuple[float, float], str]:
    with Image.open(path) as image:
        rgb = np.asarray(image.convert("RGB"), dtype=np.float32)
    green = rgb[..., 1]
    red = rgb[..., 0]
    height, width = green.shape
    yy, xx = np.indices(green.shape)
    central = (
        (xx >= 0.2 * width)
        & (xx <= 0.8 * width)
        & (yy >= 0.2 * height)
        & (yy <= 0.8 * height)
    )
    red_marker = central & (red > 1.20 * green) & (red > 45)
    green_bright = green >= np.percentile(green[central], 92.0)
    overlap = red_marker & green_bright
    if np.count_nonzero(overlap) >= 3:
        weights = red[overlap] + green[overlap]
        anchor_x = float(np.average(xx[overlap], weights=weights))
        anchor_y = float(np.average(yy[overlap], weights=weights))
        method = "red_green_measurement_overlap"
    else:
        anchor_x = (width - 1) / 2.0
        anchor_y = (height - 1) / 2.0
        method = "reference_image_center"
    return green, (anchor_x, anchor_y), method


def _register_reference(
    reference_file: Path,
    normalized_stack_zyx: np.ndarray,
    record: NNE2Record,
    spacing_xyz_um: tuple[float, float, float],
    config: NNE2Config,
) -> tuple[tuple[float, float, float], int, float, str, str]:
    assert record.ref_um_per_px is not None
    assert record.stack_index is not None
    green, anchor_xy, anchor_method = _reference_green_and_anchor(reference_file)
    height, width = green.shape
    half_px = max(16, int(round(config.registration_patch_um / record.ref_um_per_px / 2.0)))
    half_px = min(half_px, int(0.42 * min(height, width)))
    center_x, center_y = anchor_xy
    left = max(0, int(round(center_x)) - half_px)
    right = min(width, int(round(center_x)) + half_px + 1)
    top = max(0, int(round(center_y)) - half_px)
    bottom = min(height, int(round(center_y)) + half_px + 1)
    patch = green[top:bottom, left:right]
    sx, sy, sz = spacing_xyz_um
    patch_width_um = patch.shape[1] * record.ref_um_per_px
    patch_height_um = patch.shape[0] * record.ref_um_per_px
    target_shape = (
        max(12, int(round(patch_height_um / sy))),
        max(12, int(round(patch_width_um / sx))),
    )
    target_shape = (
        min(target_shape[0], normalized_stack_zyx.shape[1] - 2),
        min(target_shape[1], normalized_stack_zyx.shape[2] - 2),
    )
    template = resize(
        patch,
        target_shape,
        order=1,
        preserve_range=True,
        anti_aliasing=True,
    ).astype(np.float32)
    template = _unit_contrast(template)
    expected_z = record.stack_index - 1
    z_min = max(0, expected_z - config.registration_z_search)
    z_max = min(normalized_stack_zyx.shape[0] - 1, expected_z + config.registration_z_search)
    best_score = -np.inf
    best = (0, 0, expected_z)
    for z_index in range(z_min, z_max + 1):
        image = _unit_contrast(normalized_stack_zyx[z_index])
        response = match_template(image, template, pad_input=False)
        flat_index = int(np.nanargmax(response))
        y0, x0 = np.unravel_index(flat_index, response.shape)
        score = float(response[y0, x0])
        if score > best_score:
            best_score = score
            best = (x0, y0, z_index)
    x0, y0, z_index = best
    anchor_rel_x = (center_x - left) / max(1, patch.shape[1] - 1)
    anchor_rel_y = (center_y - top) / max(1, patch.shape[0] - 1)
    x_px = x0 + anchor_rel_x * (template.shape[1] - 1)
    y_px = y0 + anchor_rel_y * (template.shape[0] - 1)
    xyz_um = (float(x_px * sx), float(y_px * sy), float(z_index * sz))
    status = "registered" if best_score >= config.min_registration_score else "low_confidence"
    return xyz_um, z_index + 1, best_score, status, anchor_method


def match_measurement_anchors(
    records: list[NNE2Record],
    normalized_stack_zyx: np.ndarray,
    graph: HierarchicalGraphResult,
    spacing_xyz_um: tuple[float, float, float],
    config: NNE2Config,
    logger: logging.Logger | None = None,
) -> list[AnchorMatch]:
    logger = logger or logging.getLogger("ulm_3d_vascular")
    points: list[np.ndarray] = []
    branch_ids: list[np.ndarray] = []
    for branch in graph.branches:
        points.append(np.asarray(branch.points_raw_lps_um, dtype=float))
        branch_ids.append(np.full(len(branch.points_raw_lps_um), branch.branch_id, dtype=int))
    if not points:
        raise ValueError("Cannot match landmarks to an empty branch graph")
    all_points = np.concatenate(points, axis=0)
    all_branch_ids = np.concatenate(branch_ids, axis=0)
    tree = cKDTree(all_points)
    output: list[AnchorMatch] = []
    for record in records:
        assert record.reference_file is not None
        assert record.subject_id is not None
        assert record.tree_id is not None
        assert record.branching_order is not None
        assert record.depth_um is not None
        assert record.stack_index is not None
        seed, matched_index, registration_score, registration_status, anchor_method = (
            _register_reference(
                record.reference_file,
                normalized_stack_zyx,
                record,
                spacing_xyz_um,
                config,
            )
        )
        query_count = min(256, len(all_points))
        distances, point_indices = tree.query(np.asarray(seed), k=query_count)
        distances = np.atleast_1d(distances)
        point_indices = np.atleast_1d(point_indices)
        nearest_by_branch: dict[int, float] = {}
        for distance_value, point_index in zip(distances, point_indices):
            branch_id = int(all_branch_ids[int(point_index)])
            nearest_by_branch[branch_id] = min(
                float(distance_value), nearest_by_branch.get(branch_id, np.inf)
            )
        ranked_candidates = sorted(
            nearest_by_branch.items(), key=lambda item: (item[1], item[0])
        )[:5]
        matched_branch = ranked_candidates[0][0]
        distance = ranked_candidates[0][1]
        candidate_branch_ids = tuple(item[0] for item in ranked_candidates)
        candidate_distances = tuple(item[1] for item in ranked_candidates)
        close_alternatives = [
            item for item in ranked_candidates
            if item[1] <= min(config.max_anchor_distance_um, float(distance) + 5.0)
        ]
        ambiguity_status = (
            "multiple_nearby_branch_candidates"
            if len(close_alternatives) > 1
            else "single_nearest_branch_candidate"
        )
        match_status = "matched"
        if float(distance) > config.max_anchor_distance_um:
            matched_branch = None
            match_status = "too_far_from_segmented_centerline"
        elif registration_status != "registered":
            match_status = "matched_low_registration_confidence"
        output.append(
            AnchorMatch(
                record_id=record.record_id,
                subject_id=record.subject_id,
                tree_id=record.tree_id,
                branching_order=record.branching_order,
                depth_um=record.depth_um,
                expected_stack_index=record.stack_index,
                matched_stack_index=matched_index,
                seed_xyz_um=seed,
                registration_score=registration_score,
                registration_status=registration_status,
                anchor_pixel_method=anchor_method,
                matched_branch_id=matched_branch,
                branch_distance_um=float(distance),
                match_status=match_status,
                candidate_branch_ids=candidate_branch_ids,
                candidate_distances_um=candidate_distances,
                ambiguity_status=ambiguity_status,
            )
        )
        logger.info(
            "Record %d BO=%d registration=%.3f branch=%s distance=%.2f um",
            record.record_id,
            record.branching_order,
            registration_score,
            matched_branch,
            float(distance),
        )
    return output
