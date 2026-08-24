"""Step 1 equivalent: validate and clean Schmid graph records without inventing geometry."""

from __future__ import annotations

import logging

import networkx as nx
import numpy as np

from .config import SchmidPKLConfig
from .model import CleanEdge, CleanupDecision, SchmidCleanupResult, SchmidInputData, optional_float


def _polyline_length(points: np.ndarray) -> float:
    if len(points) < 2:
        return 0.0
    return float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())


def _align_geometry(
    points: np.ndarray,
    diameters: np.ndarray,
    coordinate_u: np.ndarray,
    coordinate_v: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, bool, float]:
    forward_error = float(
        np.linalg.norm(points[0] - coordinate_u) + np.linalg.norm(points[-1] - coordinate_v)
    )
    reverse_error = float(
        np.linalg.norm(points[0] - coordinate_v) + np.linalg.norm(points[-1] - coordinate_u)
    )
    if reverse_error + np.finfo(float).eps < forward_error:
        return points[::-1].copy(), diameters[::-1].copy(), True, reverse_error
    return points.copy(), diameters.copy(), False, forward_error


def _invalid_edge_reason(source: SchmidInputData, edge_id: int, valid_nodes: set[int]) -> str | None:
    u, v = (int(value) for value in source.edge_tuples[edge_id])
    if u == v:
        return "self_loop"
    if u not in valid_nodes or v not in valid_nodes:
        return "references_invalid_node"
    scalar_values = (
        source.mean_diameter_um[edge_id],
        source.flow_um3_per_ms[edge_id],
        source.source_length_um[edge_id],
    )
    if not np.all(np.isfinite(scalar_values)):
        return "non_finite_scalar"
    if source.mean_diameter_um[edge_id] <= 0:
        return "non_positive_mean_diameter"
    if source.source_length_um[edge_id] <= 0:
        return "non_positive_source_length"
    if source.flow_um3_per_ms[edge_id] < 0:
        return "negative_flow_magnitude"
    points = source.point_sequences_um[edge_id]
    diameters = source.diameter_profiles_um[edge_id]
    if len(points) < 2:
        return "fewer_than_two_geometry_points"
    if len(diameters) != len(points):
        return "diameter_profile_length_mismatch"
    if not np.all(np.isfinite(points)):
        return "non_finite_geometry_point"
    if not np.all(np.isfinite(diameters)):
        return "non_finite_diameter_profile"
    if np.any(diameters <= 0):
        return "non_positive_profile_diameter"
    return None


def clean_schmid_input(
    source: SchmidInputData,
    config: SchmidPKLConfig,
    logger: logging.Logger | None = None,
) -> SchmidCleanupResult:
    decisions: list[CleanupDecision] = []
    valid_node_mask = np.all(np.isfinite(source.coordinates_um), axis=1) & np.isfinite(
        source.pressure_mmhg
    )
    valid_nodes = set(np.flatnonzero(valid_node_mask).tolist())
    removed_node_ids = np.flatnonzero(~valid_node_mask).astype(int).tolist()
    for node_id in removed_node_ids:
        decisions.append(
            CleanupDecision("node", node_id, "removed", "non_finite_coordinate_or_pressure")
        )

    candidate_edges: list[CleanEdge] = []
    removed_edge_ids: list[int] = []
    for edge_id in range(source.edge_count):
        reason = _invalid_edge_reason(source, edge_id, valid_nodes)
        if reason is not None:
            removed_edge_ids.append(edge_id)
            decisions.append(CleanupDecision("edge", edge_id, "removed", reason))
            continue
        u, v = (int(value) for value in source.edge_tuples[edge_id])
        points, diameters, reversed_geometry, endpoint_error = _align_geometry(
            source.point_sequences_um[edge_id],
            source.diameter_profiles_um[edge_id],
            source.coordinates_um[u],
            source.coordinates_um[v],
        )
        geometric_length = _polyline_length(points)
        source_length = float(source.source_length_um[edge_id])
        length_relative_error = abs(geometric_length - source_length) / source_length
        if geometric_length <= config.endpoint_tolerance_um:
            geometry_status = "degenerate_source_geometry"
            decisions.append(
                CleanupDecision(
                    "edge",
                    edge_id,
                    "retained_with_warning",
                    "centerline points have near-zero length; topology and flow were preserved",
                )
            )
        elif endpoint_error > config.endpoint_tolerance_um:
            geometry_status = "endpoint_mismatch"
            decisions.append(
                CleanupDecision(
                    "edge",
                    edge_id,
                    "retained_with_warning",
                    f"centerline endpoint error is {endpoint_error:.6g} um",
                )
            )
        elif length_relative_error > config.length_relative_warning:
            geometry_status = "length_mismatch"
            decisions.append(
                CleanupDecision(
                    "edge",
                    edge_id,
                    "retained_with_warning",
                    f"source/geometric length relative difference is {length_relative_error:.3%}",
                )
            )
        else:
            geometry_status = "valid"
            action = "reversed_geometry" if reversed_geometry else "kept"
            reason_text = (
                "point and diameter sequences were reversed to match tuple endpoints"
                if reversed_geometry
                else "passed validation"
            )
            decisions.append(CleanupDecision("edge", edge_id, action, reason_text))

        candidate_edges.append(
            CleanEdge(
                edge_id=edge_id,
                node_u=u,
                node_v=v,
                points_u_to_v_um=points,
                diameter_u_to_v_um=diameters,
                mean_diameter_um=float(source.mean_diameter_um[edge_id]),
                source_length_um=source_length,
                geometric_length_um=geometric_length,
                flow_um3_per_ms=float(source.flow_um3_per_ms[edge_id]),
                vessel_type_code=int(source.vessel_type_code[edge_id]),
                hematocrit=optional_float(float(source.hematocrit[edge_id])),
                hematocrit_boundary=optional_float(float(source.hematocrit_boundary[edge_id])),
                red_blood_cell_count=optional_float(float(source.red_blood_cell_count[edge_id])),
                geometry_reversed=reversed_geometry,
                geometry_status=geometry_status,
                endpoint_error_um=endpoint_error,
                length_relative_error=length_relative_error,
            )
        )

    graph = nx.Graph()
    graph.add_nodes_from(valid_nodes)
    graph.add_edges_from((edge.node_u, edge.node_v) for edge in candidate_edges)
    components = list(nx.connected_components(graph))
    original_component_count = len(components)
    if config.keep_largest_component and components:
        retained_nodes = max(components, key=lambda item: (len(item), -min(item)))
    else:
        retained_nodes = set().union(*components) if components else set()
    retained_component_count = 1 if config.keep_largest_component and retained_nodes else len(components)

    edges: list[CleanEdge] = []
    for edge in candidate_edges:
        if edge.node_u in retained_nodes and edge.node_v in retained_nodes:
            edges.append(edge)
        else:
            removed_edge_ids.append(edge.edge_id)
            decisions.append(
                CleanupDecision("edge", edge.edge_id, "removed", "outside_retained_main_component")
            )
    component_removed_nodes = sorted(valid_nodes - retained_nodes)
    removed_node_ids.extend(component_removed_nodes)
    for node_id in component_removed_nodes:
        decisions.append(
            CleanupDecision("node", node_id, "removed", "outside_retained_main_component")
        )

    result = SchmidCleanupResult(
        source=source,
        valid_node_ids=np.asarray(sorted(retained_nodes), dtype=np.int64),
        edges=edges,
        decisions=decisions,
        original_component_count=original_component_count,
        retained_component_count=retained_component_count,
        removed_node_ids=sorted(set(removed_node_ids)),
        removed_edge_ids=sorted(set(removed_edge_ids)),
    )
    if logger is not None:
        report = result.report()
        logger.info(
            "Step 1 cleanup: vertices %d -> %d; edges %d -> %d; components %d -> %d",
            report["source_vertex_count"],
            report["retained_vertex_count"],
            report["source_edge_count"],
            report["retained_edge_count"],
            report["original_component_count"],
            report["retained_component_count"],
        )
        logger.warning(
            "Source edges with degenerate centerline geometry retained for topology: %d",
            report["geometry_status_counts"].get("degenerate_source_geometry", 0),
        )
    return result
