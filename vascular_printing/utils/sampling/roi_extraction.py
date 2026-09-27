"""Spatially indexed extraction of connected, globally traceable vascular ROIs."""

from __future__ import annotations

import re
from collections import Counter, defaultdict
from dataclasses import dataclass
from typing import Any, Iterable

import networkx as nx
import numpy as np
from scipy.spatial import cKDTree

from .roi_boundary import boundary_face, clip_segment_to_box, point_in_box
from .sampling_config import SamplingConfig
from .sampling_types import CutPort, GlobalEdge, GlobalVascularModel, ROIRecord


@dataclass(frozen=True, slots=True)
class ExtractionBatch:
    anchor_ids: tuple[int, ...]
    candidates: tuple[ROIRecord, ...]
    rejected_reasons: dict[str, int]


@dataclass(frozen=True, slots=True)
class _TemporaryNode:
    key: tuple[Any, ...]
    position_um: np.ndarray
    radius_um: float
    global_node_id: int | None
    boundary_face: str | None


@dataclass(frozen=True, slots=True)
class _TemporaryEdge:
    global_edge_id: int
    upstream_key: tuple[Any, ...]
    downstream_key: tuple[Any, ...]
    points_um: np.ndarray
    radius_um: np.ndarray
    length_um: float


class EdgeSpatialIndex:
    """KD-tree over global-edge midpoints followed by exact AABB filtering."""

    def __init__(self, model: GlobalVascularModel) -> None:
        self.model = model
        self.edge_min = np.asarray(
            [np.minimum(edge.upstream_position_um, edge.downstream_position_um) for edge in model.edges],
            dtype=float,
        )
        self.edge_max = np.asarray(
            [np.maximum(edge.upstream_position_um, edge.downstream_position_um) for edge in model.edges],
            dtype=float,
        )
        self.midpoints = (self.edge_min + self.edge_max) * 0.5
        half_lengths = np.linalg.norm((self.edge_max - self.edge_min) * 0.5, axis=1)
        self.maximum_half_length = float(np.max(half_lengths)) if len(half_lengths) else 0.0
        self.tree = cKDTree(self.midpoints) if len(self.midpoints) else None

    def query(self, bbox_min_um: np.ndarray, bbox_max_um: np.ndarray) -> np.ndarray:
        if self.tree is None:
            return np.empty(0, dtype=np.int64)
        lower = np.asarray(bbox_min_um, dtype=float)
        upper = np.asarray(bbox_max_um, dtype=float)
        center = (lower + upper) * 0.5
        radius = float(np.linalg.norm((upper - lower) * 0.5)) + self.maximum_half_length
        nearby = np.asarray(self.tree.query_ball_point(center, radius), dtype=np.int64)
        if not len(nearby):
            return nearby
        overlaps = np.all(self.edge_max[nearby] >= lower, axis=1) & np.all(
            self.edge_min[nearby] <= upper, axis=1
        )
        return np.sort(nearby[overlaps])


def global_model_from_swc(
    swc: Any,
    *,
    source_model_id: str,
    source_mouse_id: str,
    model_bounds_xyz_um: tuple[float, float, float, float, float, float] | None = None,
) -> GlobalVascularModel:
    """Adapt the existing SWC arrays without changing their node identities."""

    node_ids = np.asarray(swc.node_ids, dtype=np.int64)
    positions = np.asarray(swc.points_um, dtype=float)
    radii = np.asarray(swc.radius_raw_um, dtype=float)
    parents = np.asarray(swc.parent_ids, dtype=np.int64)
    index_by_id = {int(node_id): index for index, node_id in enumerate(node_ids)}
    edges: list[GlobalEdge] = []
    incident: dict[int, list[int]] = defaultdict(list)
    for child_index, (child_id, parent_id) in enumerate(zip(node_ids.tolist(), parents.tolist())):
        if int(parent_id) == -1:
            continue
        parent_index = index_by_id.get(int(parent_id))
        if parent_index is None:
            continue
        edge_id = len(edges)
        edges.append(
            GlobalEdge(
                edge_id=edge_id,
                upstream_node_id=int(parent_id),
                downstream_node_id=int(child_id),
                upstream_position_um=positions[parent_index].copy(),
                downstream_position_um=positions[child_index].copy(),
                upstream_radius_um=float(radii[parent_index]),
                downstream_radius_um=float(radii[child_index]),
            )
        )
        incident[int(parent_id)].append(edge_id)
        incident[int(child_id)].append(edge_id)
    if model_bounds_xyz_um is None:
        minimum = np.min(positions, axis=0)
        maximum = np.max(positions, axis=0)
        model_bounds_xyz_um = (
            float(minimum[0]), float(maximum[0]),
            float(minimum[1]), float(maximum[1]),
            float(minimum[2]), float(maximum[2]),
        )
    return GlobalVascularModel(
        source_model_id=source_model_id,
        source_mouse_id=source_mouse_id,
        node_ids=node_ids.copy(),
        node_positions_um=positions.copy(),
        node_radius_um=radii.copy(),
        parent_ids=parents.copy(),
        edges=tuple(edges),
        model_bounds_xyz_um=model_bounds_xyz_um,
        node_index_by_id=index_by_id,
        incident_edge_ids_by_node={key: tuple(sorted(value)) for key, value in incident.items()},
        global_degree_by_node={int(node_id): len(incident.get(int(node_id), ())) for node_id in node_ids},
    )


def _interior_anchor_candidates(model: GlobalVascularModel, roi_size_um: np.ndarray) -> np.ndarray:
    bounds_min = np.asarray(model.model_bounds_xyz_um[0::2], dtype=float)
    bounds_max = np.asarray(model.model_bounds_xyz_um[1::2], dtype=float)
    effective_size = np.minimum(roi_size_um, bounds_max - bounds_min)
    half = effective_size * 0.5
    positions = model.node_positions_um
    interior = np.all((positions >= bounds_min + half) & (positions <= bounds_max - half), axis=1)
    critical = np.asarray(
        [model.global_degree_by_node[int(node_id)] != 2 for node_id in model.node_ids],
        dtype=bool,
    )
    indices = np.flatnonzero(interior & critical)
    if not len(indices):
        indices = np.flatnonzero(interior)
    if not len(indices):
        indices = np.flatnonzero(critical)
    return indices if len(indices) else np.arange(model.node_count, dtype=np.int64)


def _greedy_minimum_distance(
    positions: np.ndarray,
    ordered_indices: Iterable[int],
    *,
    minimum_distance_um: float,
    maximum: int,
) -> list[int]:
    selected: list[int] = []
    for index in ordered_indices:
        candidate = int(index)
        if not selected or float(
            np.min(np.linalg.norm(positions[selected] - positions[candidate], axis=1))
        ) >= minimum_distance_um:
            selected.append(candidate)
        if len(selected) >= maximum:
            break
    return selected


def generate_anchor_ids(model: GlobalVascularModel, config: SamplingConfig) -> tuple[int, ...]:
    """Generate spatially non-redundant global-node anchors deterministically."""

    candidate_indices = _interior_anchor_candidates(model, np.asarray(config.roi_size_um, dtype=float))
    order = np.argsort(model.node_ids[candidate_indices], kind="stable")
    candidate_indices = candidate_indices[order]
    positions = model.node_positions_um[candidate_indices]
    rng = np.random.default_rng(config.seed)
    if config.anchor_mode in {"random", "poisson_disk"}:
        local_order = rng.permutation(len(candidate_indices)).tolist()
        selected_local = _greedy_minimum_distance(
            positions,
            local_order,
            minimum_distance_um=config.min_anchor_distance_um,
            maximum=config.max_candidate_anchors,
        )
    else:
        centroid = np.mean(positions, axis=0)
        first = int(np.argmin(np.linalg.norm(positions - centroid, axis=1)))
        selected_local = [first]
        minimum_distance = np.linalg.norm(positions - positions[first], axis=1)
        minimum_distance[first] = -np.inf
        while len(selected_local) < min(config.max_candidate_anchors, len(candidate_indices)):
            next_index = int(np.argmax(minimum_distance))
            if float(minimum_distance[next_index]) < config.min_anchor_distance_um:
                break
            selected_local.append(next_index)
            minimum_distance = np.minimum(
                minimum_distance,
                np.linalg.norm(positions - positions[next_index], axis=1),
            )
            minimum_distance[selected_local] = -np.inf
    return tuple(int(model.node_ids[candidate_indices[index]]) for index in selected_local)


def _bbox_for_anchor(
    model: GlobalVascularModel,
    anchor_position_um: np.ndarray,
    requested_size_um: tuple[float, float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    model_min = np.asarray(model.model_bounds_xyz_um[0::2], dtype=float)
    model_max = np.asarray(model.model_bounds_xyz_um[1::2], dtype=float)
    effective_size = np.minimum(np.asarray(requested_size_um, dtype=float), model_max - model_min)
    half = effective_size * 0.5
    center = np.maximum(np.asarray(anchor_position_um, dtype=float), model_min + half)
    center = np.minimum(center, model_max - half)
    return center - half, center + half, center


def _node_sort_key(key: tuple[Any, ...]) -> tuple[Any, ...]:
    return (0, int(key[1])) if key[0] == "g" else (1, int(key[1]), str(key[2]))


def _safe_identifier(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._") or "model"


def extract_connected_roi(
    model: GlobalVascularModel,
    spatial_index: EdgeSpatialIndex,
    *,
    anchor_id: int,
    roi_size_um: tuple[float, float, float],
) -> tuple[ROIRecord | None, str | None]:
    """Clip global edges, then retain only the connected component containing the anchor."""

    anchor_index = model.node_index_by_id[int(anchor_id)]
    anchor_position = model.node_positions_um[anchor_index]
    bbox_min, bbox_max, bbox_center = _bbox_for_anchor(model, anchor_position, roi_size_um)
    node_payload: dict[tuple[Any, ...], _TemporaryNode] = {}
    temporary_edges: list[_TemporaryEdge] = []
    for edge_index in spatial_index.query(bbox_min, bbox_max).tolist():
        edge = model.edges[int(edge_index)]
        clipped = clip_segment_to_box(
            edge.upstream_position_um,
            edge.downstream_position_um,
            bbox_min,
            bbox_max,
        )
        if clipped is None:
            continue
        start_is_global = clipped.start_t <= 1.0e-10 and point_in_box(
            edge.upstream_position_um, bbox_min, bbox_max
        )
        end_is_global = clipped.end_t >= 1.0 - 1.0e-10 and point_in_box(
            edge.downstream_position_um, bbox_min, bbox_max
        )
        upstream_key: tuple[Any, ...] = (
            ("g", edge.upstream_node_id)
            if start_is_global
            else ("c", edge.edge_id, "entry")
        )
        downstream_key: tuple[Any, ...] = (
            ("g", edge.downstream_node_id)
            if end_is_global
            else ("c", edge.edge_id, "exit")
        )
        start_radius = edge.upstream_radius_um + clipped.start_t * (
            edge.downstream_radius_um - edge.upstream_radius_um
        )
        end_radius = edge.upstream_radius_um + clipped.end_t * (
            edge.downstream_radius_um - edge.upstream_radius_um
        )
        node_payload.setdefault(
            upstream_key,
            _TemporaryNode(
                upstream_key,
                clipped.start_um,
                float(start_radius),
                edge.upstream_node_id if start_is_global else None,
                clipped.start_face,
            ),
        )
        node_payload.setdefault(
            downstream_key,
            _TemporaryNode(
                downstream_key,
                clipped.end_um,
                float(end_radius),
                edge.downstream_node_id if end_is_global else None,
                clipped.end_face,
            ),
        )
        length = float(np.linalg.norm(clipped.end_um - clipped.start_um))
        temporary_edges.append(
            _TemporaryEdge(
                edge.edge_id,
                upstream_key,
                downstream_key,
                np.asarray((clipped.start_um, clipped.end_um), dtype=float),
                np.asarray((start_radius, end_radius), dtype=float),
                length,
            )
        )
    if not temporary_edges:
        return None, "empty"

    raw_graph = nx.Graph()
    raw_graph.add_nodes_from(node_payload)
    for index, edge in enumerate(temporary_edges):
        raw_graph.add_edge(edge.upstream_key, edge.downstream_key, temporary_edge_index=index)
    components = list(nx.connected_components(raw_graph))
    anchor_key = ("g", int(anchor_id))
    retained_nodes = next((nodes for nodes in components if anchor_key in nodes), None)
    if retained_nodes is None:
        return None, "anchor_component_missing"
    retained_edges = [
        edge
        for edge in temporary_edges
        if edge.upstream_key in retained_nodes and edge.downstream_key in retained_nodes
    ]
    if not retained_edges:
        return None, "empty"
    if any(
        not np.all(np.isfinite(edge.radius_um)) or np.any(edge.radius_um <= 0)
        for edge in retained_edges
    ):
        return None, "invalid_radius"

    ordered_keys = sorted(retained_nodes, key=_node_sort_key)
    local_id_by_key = {key: index for index, key in enumerate(ordered_keys)}
    local_positions = np.asarray([node_payload[key].position_um for key in ordered_keys], dtype=float)
    local_radii = np.asarray([node_payload[key].radius_um for key in ordered_keys], dtype=float)
    local_global_ids = np.asarray(
        [node_payload[key].global_node_id if node_payload[key].global_node_id is not None else -1 for key in ordered_keys],
        dtype=np.int64,
    )
    local_edges = np.asarray(
        [
            (local_id_by_key[edge.upstream_key], local_id_by_key[edge.downstream_key])
            for edge in retained_edges
        ],
        dtype=np.int64,
    ).reshape((-1, 2))
    local_graph = nx.Graph()
    local_graph.add_nodes_from(range(len(ordered_keys)))
    local_graph.add_edges_from((int(a), int(b)) for a, b in local_edges)
    if not nx.is_connected(local_graph):
        return None, "disconnected"

    retained_global_edge_ids = {edge.global_edge_id for edge in retained_edges}
    true_terminal_local_ids: list[int] = []
    true_terminal_global_ids: list[int] = []
    cut_ports: list[CutPort] = []
    for local_id, key in enumerate(ordered_keys):
        if local_graph.degree(local_id) != 1:
            continue
        payload = node_payload[key]
        if payload.global_node_id is not None and model.global_degree_by_node[payload.global_node_id] == 1:
            true_terminal_local_ids.append(local_id)
            true_terminal_global_ids.append(payload.global_node_id)
            continue
        if payload.global_node_id is None:
            incident_retained = [
                edge.global_edge_id
                for edge in retained_edges
                if key in {edge.upstream_key, edge.downstream_key}
            ]
            cut_edge_ids = incident_retained[:1]
        else:
            cut_edge_ids = [
                edge_id
                for edge_id in model.incident_edge_ids_by_node.get(payload.global_node_id, ())
                if edge_id not in retained_global_edge_ids
            ]
            if not cut_edge_ids:
                cut_edge_ids = [
                    edge.global_edge_id
                    for edge in retained_edges
                    if key in {edge.upstream_key, edge.downstream_key}
                ][:1]
        for cut_edge_id in cut_edge_ids:
            cut_ports.append(
                CutPort(
                    cut_port_id="pending",
                    local_node_id=local_id,
                    global_edge_id=int(cut_edge_id),
                    intersection_position_um=tuple(float(value) for value in payload.position_um),
                    radius_at_cut_um=float(payload.radius_um),
                    boundary_face=payload.boundary_face
                    or boundary_face(payload.position_um, bbox_min, bbox_max),
                )
            )

    roi_id = f"{_safe_identifier(model.source_model_id)}__anchor_{int(anchor_id):06d}"
    numbered_ports = tuple(
        CutPort(
            cut_port_id=f"{roi_id}__cut_{index:03d}",
            local_node_id=port.local_node_id,
            global_edge_id=port.global_edge_id,
            intersection_position_um=port.intersection_position_um,
            radius_at_cut_um=port.radius_at_cut_um,
            boundary_face=port.boundary_face,
        )
        for index, port in enumerate(cut_ports)
    )
    return ROIRecord(
        roi_id=roi_id,
        source_model_id=model.source_model_id,
        source_mouse_id=model.source_mouse_id,
        anchor_id=int(anchor_id),
        anchor_position_um=tuple(float(value) for value in anchor_position),
        bbox_min_um=tuple(float(value) for value in bbox_min),
        bbox_max_um=tuple(float(value) for value in bbox_max),
        bbox_center_um=tuple(float(value) for value in bbox_center),
        bbox_size_um=tuple(float(value) for value in bbox_max - bbox_min),
        global_node_ids=tuple(sorted(int(value) for value in local_global_ids if value >= 0)),
        global_edge_ids=tuple(sorted(edge.global_edge_id for edge in retained_edges)),
        local_node_ids=np.arange(len(ordered_keys), dtype=np.int64),
        local_node_global_ids=local_global_ids,
        local_node_positions_um=local_positions,
        local_node_radius_um=local_radii,
        local_edges=local_edges,
        local_edge_ids=np.arange(len(retained_edges), dtype=np.int64),
        local_edge_global_ids=np.asarray(
            [edge.global_edge_id for edge in retained_edges], dtype=np.int64
        ),
        local_edge_points_um=np.asarray([edge.points_um for edge in retained_edges], dtype=float),
        local_edge_radius_um=np.asarray([edge.radius_um for edge in retained_edges], dtype=float),
        true_terminal_local_ids=tuple(true_terminal_local_ids),
        true_terminal_global_ids=tuple(true_terminal_global_ids),
        cut_ports=numbered_ports,
        raw_component_count=len(components),
        raw_total_vessel_length_um=float(sum(edge.length_um for edge in temporary_edges)),
        retained_component_length_um=float(sum(edge.length_um for edge in retained_edges)),
    ), None


def extract_candidate_rois(
    model: GlobalVascularModel,
    config: SamplingConfig,
    *,
    anchor_ids: tuple[int, ...] | None = None,
) -> ExtractionBatch:
    """Generate anchors and extract their connected source-derived subgraphs."""

    anchors = anchor_ids if anchor_ids is not None else generate_anchor_ids(model, config)
    index = EdgeSpatialIndex(model)
    candidates: list[ROIRecord] = []
    rejected: Counter[str] = Counter()
    for anchor_id in anchors:
        roi, reason = extract_connected_roi(
            model,
            index,
            anchor_id=anchor_id,
            roi_size_um=config.roi_size_um,
        )
        if roi is None:
            rejected[reason or "invalid_graph"] += 1
            continue
        if config.max_cut_ports is not None and roi.cut_port_count > config.max_cut_ports:
            rejected["too_many_cut_ports"] += 1
            continue
        candidates.append(roi)
    return ExtractionBatch(anchors, tuple(candidates), dict(sorted(rejected.items())))
