"""SWC-centric reference preservation and analysis-network selection.

The released, human-corrected SWC is the authoritative topology.  Segmentation
masks are deliberately limited to registration QC and display: they never add,
remove, or redirect an SWC edge in this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
from scipy import ndimage

from .swc_io import SWCData, swc_from_arrays


@dataclass(slots=True)
class SWCAnalysisResult:
    """Complete reference SWC plus the single component used downstream."""

    reference_swc: SWCData
    analysis_swc: SWCData
    component_records: list[dict[str, Any]]
    reference_only_node_ids: np.ndarray
    summary: dict[str, Any]


def _component_length_um(
    swc: SWCData,
    indices: np.ndarray,
    index_by_id: dict[int, int],
) -> float:
    node_ids = {int(swc.node_ids[index]) for index in indices}
    length = 0.0
    for index in indices:
        parent_id = int(swc.parent_ids[index])
        if parent_id == -1 or parent_id not in node_ids:
            continue
        length += float(
            np.linalg.norm(swc.points_um[index] - swc.points_um[index_by_id[parent_id]])
        )
    return length


def select_analysis_swc(
    reference_swc: SWCData,
    *,
    spacing_xyz_um: tuple[float, float, float],
    volume_shape_zyx: tuple[int, int, int] | None,
    analysis_component_id: int | None = None,
) -> SWCAnalysisResult:
    """Preserve the reference forest and select one unchanged component for analysis.

    When ``analysis_component_id`` is omitted, the component with the greatest
    physical centerline length is selected.  Selection never means that the other
    components are erroneous; they remain part of ``reference_swc`` and are marked
    ``REFERENCE_ONLY`` in the audit table.
    """

    graph = reference_swc.directed_graph()
    components = [set(nodes) for nodes in nx.weakly_connected_components(graph)]
    if not components:
        raise ValueError("Reference SWC contains no connected component")

    index_by_id = {
        int(node_id): index for index, node_id in enumerate(reference_swc.node_ids)
    }
    records: list[dict[str, Any]] = []
    component_indices: list[np.ndarray] = []
    for component_id, nodes in enumerate(components):
        indices = np.asarray(
            [index for index, node_id in enumerate(reference_swc.node_ids) if int(node_id) in nodes],
            dtype=np.int64,
        )
        component_indices.append(indices)
        points = reference_swc.points_um[indices]
        extent = np.ptp(points, axis=0) if len(points) else np.zeros(3, dtype=float)
        length_um = _component_length_um(reference_swc, indices, index_by_id)
        roots = [
            int(reference_swc.node_ids[index])
            for index in indices
            if int(reference_swc.parent_ids[index]) == -1
        ]
        records.append(
            {
                "swc_component_id": component_id,
                "node_count": int(len(indices)),
                "edge_count": int(np.count_nonzero(reference_swc.parent_ids[indices] != -1)),
                "total_length_um": length_um,
                "bbox_x_um": float(extent[0]),
                "bbox_y_um": float(extent[1]),
                "bbox_z_um": float(extent[2]),
                "bbox_diagonal_um": float(np.linalg.norm(extent)),
                "root_count": len(roots),
                "root_ids": "|".join(str(value) for value in roots),
                "source_index_min": int(indices[0]),
                "source_index_max": int(indices[-1]),
            }
        )

    ranking = sorted(
        range(len(records)),
        key=lambda component_id: (
            -float(records[component_id]["total_length_um"]),
            -int(records[component_id]["node_count"]),
            -float(records[component_id]["bbox_diagonal_um"]),
            component_id,
        ),
    )
    rank_by_component = {component_id: rank + 1 for rank, component_id in enumerate(ranking)}
    if analysis_component_id is None:
        selected_component_id = int(ranking[0])
        selection_rule = "maximum_total_centerline_length"
    else:
        selected_component_id = int(analysis_component_id)
        if selected_component_id < 0 or selected_component_id >= len(components):
            raise ValueError(
                f"analysis_component_id={selected_component_id} is outside the available "
                f"component range 0..{len(components) - 1}"
            )
        selection_rule = "explicit_component_id"

    for component_id, record in enumerate(records):
        selected = component_id == selected_component_id
        record.update(
            {
                "length_rank": rank_by_component[component_id],
                "selected_for_analysis": selected,
                "decision": "ANALYSIS_NETWORK" if selected else "REFERENCE_ONLY",
                "decision_reason": (
                    "explicitly selected component"
                    if selected and analysis_component_id is not None
                    else (
                        "greatest total centerline length; selected for the single-network workflow"
                        if selected
                        else "preserved in reference_swc but outside the current single-network analysis"
                    )
                ),
            }
        )

    selected_indices = component_indices[selected_component_id]
    analysis_swc = swc_from_arrays(
        path=Path(reference_swc.path),
        node_ids=reference_swc.node_ids[selected_indices],
        type_codes=reference_swc.type_codes[selected_indices],
        points_voxel_xyz=reference_swc.points_voxel_xyz[selected_indices],
        radius_raw_um=reference_swc.radius_raw_um[selected_indices],
        parent_ids=reference_swc.parent_ids[selected_indices],
        spacing_xyz_um=spacing_xyz_um,
        volume_shape_zyx=volume_shape_zyx,
    )
    selected_node_ids = set(analysis_swc.node_ids.astype(int).tolist())
    reference_only_node_ids = np.asarray(
        [
            int(node_id)
            for node_id in reference_swc.node_ids
            if int(node_id) not in selected_node_ids
        ],
        dtype=np.int64,
    )
    selected_record = records[selected_component_id]
    excluded_records = [
        record for component_id, record in enumerate(records) if component_id != selected_component_id
    ]
    summary = {
        "method": "swc_centric_reference_and_analysis_selection",
        "reference_swc_role": "complete_human_corrected_topology_reference",
        "analysis_swc_role": "selected_single_component_for_graph_roi_and_simulation",
        "component_selection_rule": selection_rule,
        "selected_component_id": selected_component_id,
        "reference_component_count": reference_swc.component_count,
        "analysis_component_count": analysis_swc.component_count,
        "reference_node_count": reference_swc.node_count,
        "reference_edge_count": reference_swc.edge_count,
        "reference_total_length_um": float(
            sum(float(record["total_length_um"]) for record in records)
        ),
        "analysis_node_count": analysis_swc.node_count,
        "analysis_edge_count": analysis_swc.edge_count,
        "analysis_total_length_um": float(selected_record["total_length_um"]),
        "reference_only_component_count": len(excluded_records),
        "reference_only_node_count": int(len(reference_only_node_ids)),
        "reference_only_edge_count": int(
            sum(int(record["edge_count"]) for record in excluded_records)
        ),
        "reference_only_total_length_um": float(
            sum(float(record["total_length_um"]) for record in excluded_records)
        ),
        "reference_topology_modified": False,
        "analysis_topology_modified": False,
        "new_node_count": 0,
        "new_edge_count": 0,
        "parent_relation_change_count": 0,
        "reference_only_components_are_errors": False,
        "selected_root_ids": analysis_swc.root_ids,
        "strict_qc_status": "PASS",
    }
    return SWCAnalysisResult(
        reference_swc=reference_swc,
        analysis_swc=analysis_swc,
        component_records=records,
        reference_only_node_ids=reference_only_node_ids,
        summary=summary,
    )


def _dense_segment_voxels(start_xyz: np.ndarray, end_xyz: np.ndarray) -> np.ndarray:
    maximum_delta = float(np.max(np.abs(end_xyz - start_xyz)))
    count = max(2, int(np.ceil(maximum_delta * 2.0)) + 1)
    weights = np.linspace(0.0, 1.0, count)[:, None]
    return start_xyz[None, :] + weights * (end_xyz - start_xyz)[None, :]


def _mask_support(mask_zyx: np.ndarray, swc: SWCData) -> tuple[float, float]:
    foreground = np.asarray(mask_zyx, dtype=bool)
    shape_xyz = np.asarray(foreground.shape[::-1], dtype=np.int64)
    node_indices = np.rint(swc.points_voxel_xyz).astype(np.int64)
    node_valid = np.all((node_indices >= 0) & (node_indices < shape_xyz), axis=1)
    node_supported = np.zeros(swc.node_count, dtype=bool)
    usable = node_indices[node_valid]
    node_supported[node_valid] = foreground[usable[:, 2], usable[:, 1], usable[:, 0]]

    index_by_id = {int(node_id): index for index, node_id in enumerate(swc.node_ids)}
    supported_samples = 0
    total_samples = 0
    for index, parent_raw in enumerate(swc.parent_ids):
        parent_id = int(parent_raw)
        if parent_id == -1:
            continue
        samples = _dense_segment_voxels(
            swc.points_voxel_xyz[index_by_id[parent_id]],
            swc.points_voxel_xyz[index],
        )
        sample_indices = np.rint(samples).astype(np.int64)
        valid = np.all((sample_indices >= 0) & (sample_indices < shape_xyz), axis=1)
        values = np.zeros(len(sample_indices), dtype=bool)
        usable = sample_indices[valid]
        values[valid] = foreground[usable[:, 2], usable[:, 1], usable[:, 0]]
        supported_samples += int(np.count_nonzero(values))
        total_samples += len(values)
    return (
        float(np.mean(node_supported)) if swc.node_count else 0.0,
        float(supported_samples / total_samples) if total_samples else 1.0,
    )


def evaluate_optional_mask_qc(
    mask_zyx: np.ndarray | None,
    reference_swc: SWCData,
    analysis_swc: SWCData,
) -> dict[str, Any]:
    """Measure registration support without using the mask to alter topology."""

    if mask_zyx is None:
        return {
            "role": "optional_registration_qc_and_visualization_only",
            "available": False,
            "skip_reason": "Mask TIFF was not provided; SWC analysis remains enabled.",
            "used_for_component_selection": False,
            "used_for_node_or_edge_removal": False,
            "used_for_topology_repair": False,
            "mask_modified": False,
            "foreground_voxel_count": None,
            "component_count_6": None,
            "component_count_18": None,
            "component_count_26": None,
            "reference_swc_node_support_fraction": None,
            "reference_swc_dense_edge_support_fraction": None,
            "analysis_swc_node_support_fraction": None,
            "analysis_swc_dense_edge_support_fraction": None,
        }

    foreground = np.asarray(mask_zyx) > 0
    component_counts = {}
    for name, rank in (("6", 1), ("18", 2), ("26", 3)):
        _, count = ndimage.label(
            foreground,
            structure=ndimage.generate_binary_structure(3, rank),
        )
        component_counts[name] = int(count)
    reference_node_support, reference_edge_support = _mask_support(
        foreground, reference_swc
    )
    analysis_node_support, analysis_edge_support = _mask_support(
        foreground, analysis_swc
    )
    return {
        "role": "optional_registration_qc_and_visualization_only",
        "available": True,
        "skip_reason": None,
        "used_for_component_selection": False,
        "used_for_node_or_edge_removal": False,
        "used_for_topology_repair": False,
        "mask_modified": False,
        "foreground_voxel_count": int(np.count_nonzero(foreground)),
        "component_count_6": component_counts["6"],
        "component_count_18": component_counts["18"],
        "component_count_26": component_counts["26"],
        "reference_swc_node_support_fraction": reference_node_support,
        "reference_swc_dense_edge_support_fraction": reference_edge_support,
        "analysis_swc_node_support_fraction": analysis_node_support,
        "analysis_swc_dense_edge_support_fraction": analysis_edge_support,
    }
