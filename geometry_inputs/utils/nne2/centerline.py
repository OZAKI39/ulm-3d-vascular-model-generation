"""NNE2 centerline extraction and conversion to the common branch graph."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy import ndimage

from ..config import HierarchicalGraphConfig
from ..graph.extraction import build_hierarchical_graph
from ..graph.model import HierarchicalGraphResult
from ..voxel.skeleton import SkeletonResult, extract_coarse_skeleton
from .config import NNE2Config


@dataclass(slots=True)
class NNE2CenterlineResult:
    skeleton_zyx: np.ndarray
    radius_zyx_um: np.ndarray
    skeleton_report: SkeletonResult
    graph: HierarchicalGraphResult

    def report(self) -> dict[str, Any]:
        return {
            "skeleton": self.skeleton_report.report(),
            "graph": self.graph.report(),
            "array_axis_order": "ZYX in NNE2 volumes; XYZ in graph voxel indices",
            "radius_volume": {
                "unit": "micrometer",
                "minimum_inside_mask_um": float(
                    np.min(self.radius_zyx_um[self.radius_zyx_um > 0])
                ) if np.any(self.radius_zyx_um > 0) else 0.0,
                "maximum_um": float(np.max(self.radius_zyx_um)),
            },
        }


def graph_config_from_nne2(config: NNE2Config) -> HierarchicalGraphConfig:
    """Translate NNE2 command-line settings to the common Step 3 settings."""
    return HierarchicalGraphConfig(
        neighbor_connectivity=config.graph_connectivity,
        smoothing_enabled=config.graph_smoothing_enabled,
        smoothing_window_points=config.graph_smoothing_window_points,
        resample_step_um=config.graph_resample_step_um,
        junction_direction_distance_um=config.junction_direction_distance_um,
        short_branch_warning_um=config.short_branch_warning_um,
        large_junction_warning_voxels=config.large_junction_warning_voxels,
        high_degree_warning=config.high_degree_warning,
        smoothing_deviation_warning_um=config.smoothing_deviation_warning_um,
        save_graphml=config.save_graphml,
        save_vtp=config.save_vtp,
        save_npz=config.save_npz,
    )


def extract_centerline(
    mask_zyx: np.ndarray,
    spacing_xyz_um: tuple[float, float, float],
    *,
    cached_skeleton_zyx: np.ndarray | None = None,
    logger: logging.Logger | None = None,
) -> tuple[np.ndarray, np.ndarray, SkeletonResult]:
    """Perform Step 2 only: radius estimate and coarse centerline extraction."""
    logger = logger or logging.getLogger("ulm_3d_vascular")
    mask_zyx = np.asarray(mask_zyx, dtype=bool)
    mask_xyz = np.transpose(mask_zyx, (2, 1, 0))
    if cached_skeleton_zyx is None:
        skeleton_result = extract_coarse_skeleton(mask_xyz, spacing_xyz_um, logger)
    else:
        cached_xyz = np.transpose(np.asarray(cached_skeleton_zyx, dtype=bool), (2, 1, 0))
        if cached_xyz.shape != mask_xyz.shape or np.any(cached_xyz & ~mask_xyz):
            raise ValueError("Cached NNE2 skeleton is incompatible with its vessel mask")
        component_count = int(
            ndimage.label(cached_xyz, structure=ndimage.generate_binary_structure(3, 3))[1]
        )
        skeleton_result = SkeletonResult(
            skeleton=cached_xyz,
            skeleton_voxel_count=int(np.count_nonzero(cached_xyz)),
            connected_component_count=component_count,
            voxels_outside_mask=0,
            approximate_length_um=(
                int(np.count_nonzero(cached_xyz)) * float(np.mean(spacing_xyz_um))
            ),
        )
    sx, sy, sz = spacing_xyz_um
    radius_zyx_um = ndimage.distance_transform_edt(
        mask_zyx, sampling=(sz, sy, sx)
    ).astype(np.float32)
    return (
        np.transpose(skeleton_result.skeleton, (2, 1, 0)),
        radius_zyx_um,
        skeleton_result,
    )


def build_graph_from_centerline(
    skeleton_zyx: np.ndarray,
    mask_zyx: np.ndarray,
    spacing_xyz_um: tuple[float, float, float],
    graph_config: HierarchicalGraphConfig,
    *,
    logger: logging.Logger | None = None,
) -> HierarchicalGraphResult:
    """Perform Step 3 graph extraction from completed Step 2 arrays."""
    logger = logger or logging.getLogger("ulm_3d_vascular")
    graph_config.validate()
    skeleton_xyz = np.transpose(np.asarray(skeleton_zyx, dtype=bool), (2, 1, 0))
    mask_xyz = np.transpose(np.asarray(mask_zyx, dtype=bool), (2, 1, 0))
    graph = build_hierarchical_graph(
        skeleton_xyz,
        mask_xyz,
        origin_lps_um=(0.0, 0.0, 0.0),
        spacing_um=spacing_xyz_um,
        config=graph_config,
        logger=logger,
    )
    graph.coordinate_system = "NNE2_stack_XYZ_not_anatomically_registered"
    graph.radius_source = "NNE2 segmented voxel distance transform; coarse estimate"
    return graph


def build_centerline_graph(
    mask_zyx: np.ndarray,
    spacing_xyz_um: tuple[float, float, float],
    *,
    cached_skeleton_zyx: np.ndarray | None = None,
    graph_connectivity: int = 26,
    graph_config: HierarchicalGraphConfig | None = None,
    logger: logging.Logger | None = None,
) -> NNE2CenterlineResult:
    logger = logger or logging.getLogger("ulm_3d_vascular")
    skeleton_zyx, radius_zyx_um, skeleton_result = extract_centerline(
        mask_zyx,
        spacing_xyz_um,
        cached_skeleton_zyx=cached_skeleton_zyx,
        logger=logger,
    )
    effective_config = graph_config or HierarchicalGraphConfig(
        neighbor_connectivity=graph_connectivity, smoothing_enabled=False
    )
    graph = build_graph_from_centerline(
        skeleton_zyx,
        mask_zyx,
        spacing_xyz_um,
        effective_config,
        logger=logger,
    )
    return NNE2CenterlineResult(
        skeleton_zyx=skeleton_zyx,
        radius_zyx_um=radius_zyx_um,
        skeleton_report=skeleton_result,
        graph=graph,
    )
