"""v4 separation of immutable representative CORE_ROI from an expanded CFD_DOMAIN."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np

from utils.rodent_vasculature.swc_io import load_normalized_swc
from utils.sampling.roi_extraction import global_model_from_swc
from utils.sampling.sampling_types import CutPort, GlobalVascularModel, ROIRecord
from utils.sampling.structural_features import _branch_paths

from .config import CFDLumenConfig
from .geometry_preprocess import resample_branches, validate_and_extract_branches
from .local_implicit_junction import define_junction_collars
from .types import ContextDomainResult, GeometryValidationError


@dataclass(frozen=True, slots=True)
class _RootSegment:
    owner_index: int
    core_local_node_id: int
    outward_global_node_id: int
    global_edge_id: int
    length_um: float


@dataclass(frozen=True, slots=True)
class _Frontier:
    owner_index: int
    global_node_id: int
    incoming_global_edge_id: int


def _core_signature(roi: ROIRecord) -> str:
    digest = hashlib.sha256()
    for array in (
        roi.local_node_ids,
        roi.local_node_global_ids,
        roi.local_node_positions_um,
        roi.local_node_radius_um,
        roi.local_edges,
        roi.local_edge_ids,
        roi.local_edge_global_ids,
        roi.local_edge_points_um,
        roi.local_edge_radius_um,
    ):
        contiguous = np.ascontiguousarray(array)
        digest.update(str(contiguous.dtype).encode())
        digest.update(str(contiguous.shape).encode())
        digest.update(contiguous.tobytes())
    digest.update(roi.roi_id.encode())
    digest.update(json.dumps([port.report() for port in roi.cut_ports], sort_keys=True).encode())
    return digest.hexdigest()


def _load_verified_model(
    sampling_run: Path,
    core_roi: ROIRecord,
    config: CFDLumenConfig,
) -> tuple[GlobalVascularModel, Path]:
    project_root = Path(sampling_run).resolve().parents[2]
    configured = config.context_domain.source_rodent_run
    if configured:
        roots = [Path(configured).expanduser().resolve()]
    else:
        base = project_root / "outputs" / "rodent_vasculature"
        roots = sorted(
            (path for path in base.iterdir() if path.is_dir()),
            key=lambda path: (path.stat().st_mtime_ns, path.name),
            reverse=True,
        ) if base.is_dir() else []
    manifests = [
        root / "samples" / core_roi.source_model_id / "preprocess_manifest.json"
        for root in roots
    ]
    edge_manifest = Path(sampling_run) / "manifests" / "global_edges.csv"
    expected_edges: dict[int, tuple[int, int]] = {}
    if edge_manifest.is_file():
        with edge_manifest.open(encoding="utf-8-sig", newline="") as stream:
            for row in csv.DictReader(stream):
                if row["source_model_id"] == core_roi.source_model_id:
                    expected_edges[int(row["global_edge_id"])] = (
                        int(row["upstream_global_node_id"]),
                        int(row["downstream_global_node_id"]),
                    )
    failures: list[str] = []
    for manifest_path in manifests:
        if not manifest_path.is_file():
            continue
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            swc = load_normalized_swc(
                Path(manifest["normalized_swc_path"]),
                source_path=Path(manifest["record"]["swc_path"]),
                spacing_xyz_um=tuple(map(float, manifest["spacing_xyz_um"])),
                volume_shape_zyx=tuple(map(int, manifest["image_metadata"]["shape_zyx"])),
            )
            model = global_model_from_swc(
                swc,
                source_model_id=core_roi.source_model_id,
                source_mouse_id=core_roi.source_mouse_id,
            )
            if config.context_domain.verify_source_geometry:
                actual_edges = {
                    edge.edge_id: (edge.upstream_node_id, edge.downstream_node_id)
                    for edge in model.edges
                }
                if expected_edges and actual_edges != expected_edges:
                    raise ValueError("global edge identities differ from the Sampling manifest")
                for local_id, global_id in enumerate(core_roi.local_node_global_ids):
                    if int(global_id) < 0:
                        continue
                    source_index = model.node_index_by_id.get(int(global_id))
                    if source_index is None:
                        raise ValueError(f"CORE global node {int(global_id)} is absent")
                    if not np.allclose(
                        model.node_positions_um[source_index],
                        core_roi.local_node_positions_um[local_id],
                        rtol=0.0,
                        atol=1.0e-10,
                    ) or not np.isclose(
                        model.node_radius_um[source_index],
                        core_roi.local_node_radius_um[local_id],
                        rtol=0.0,
                        atol=1.0e-10,
                    ):
                        raise ValueError(f"CORE node {int(global_id)} geometry differs from source")
            return model, manifest_path.resolve()
        except Exception as exc:
            failures.append(f"{manifest_path}: {type(exc).__name__}: {exc}")
    detail = failures[0] if failures else "no matching preprocess_manifest.json was found"
    raise GeometryValidationError(
        f"CFD_CONTEXT_SOURCE_NOT_FOUND: cannot verify global SWC for "
        f"{core_roi.source_model_id}: {detail}"
    )


def _global_graph(model: GlobalVascularModel) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(map(int, model.node_ids))
    for edge in model.edges:
        graph.add_edge(
            edge.upstream_node_id,
            edge.downstream_node_id,
            global_edge_id=edge.edge_id,
            length_um=float(
                np.linalg.norm(edge.downstream_position_um - edge.upstream_position_um)
            ),
        )
    return graph


def _global_branch_map(graph: nx.Graph) -> dict[int, int]:
    output: dict[int, int] = {}
    for branch_id, path in enumerate(_branch_paths(graph)):
        for first, second in zip(path[:-1], path[1:]):
            output[int(graph.edges[first, second]["global_edge_id"])] = branch_id
    return output


def _core_port_distances(roi: ROIRecord) -> dict[int, float | None]:
    graph = nx.Graph()
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    for first, second in np.asarray(roi.local_edges, dtype=np.int64):
        graph.add_edge(
            int(first),
            int(second),
            weight=float(np.linalg.norm(positions[first] - positions[second])),
        )
    junctions = [int(node) for node in graph if graph.degree(node) >= 3]
    output: dict[int, float | None] = {}
    for index, port in enumerate(roi.cut_ports):
        values = [
            nx.shortest_path_length(graph, int(port.local_node_id), node, weight="weight")
            for node in junctions
            if nx.has_path(graph, int(port.local_node_id), node)
        ]
        output[index] = float(min(values)) if values else None
    return output


class _DomainBuilder:
    def __init__(
        self,
        core_roi: ROIRecord,
        model: GlobalVascularModel,
        config: CFDLumenConfig,
    ) -> None:
        self.core = core_roi
        self.model = model
        self.config = config
        self.graph = _global_graph(model)
        self.root_segments: dict[int, _RootSegment] = {}
        self.added_full_edges: dict[int, int] = {}
        self.frontiers: list[_Frontier] = []
        self.extension_reasons: dict[int, str] = {}

    def _position(self, global_node_id: int) -> np.ndarray:
        return np.asarray(
            self.model.node_positions_um[self.model.node_index_by_id[global_node_id]],
            dtype=float,
        )

    def start_extension(self, owner_index: int, reason: str) -> None:
        if owner_index in self.root_segments:
            self.expand_owner(owner_index)
            return
        port = self.core.cut_ports[owner_index]
        local_id = int(port.local_node_id)
        incident = np.asarray(self.core.local_edges)[
            np.any(np.asarray(self.core.local_edges) == local_id, axis=1)
        ]
        if len(incident) != 1:
            raise GeometryValidationError(
                f"CORE port {port.cut_port_id} is not a degree-one node"
            )
        neighbor = int(
            incident[0, 1] if int(incident[0, 0]) == local_id else incident[0, 0]
        )
        cut = np.asarray(port.intersection_position_um, dtype=float)
        outward = cut - np.asarray(self.core.local_node_positions_um[neighbor], dtype=float)
        edge = self.model.edges[int(port.global_edge_id)]
        candidates = (edge.upstream_node_id, edge.downstream_node_id)
        outward_node = max(
            candidates,
            key=lambda node_id: float(np.dot(self._position(node_id) - cut, outward)),
        )
        length = float(np.linalg.norm(self._position(outward_node) - cut))
        if length <= self.config.geometry.minimum_edge_length_um:
            raise GeometryValidationError(
                f"CFD context extension from {port.cut_port_id} has a zero-length source segment"
            )
        self.root_segments[owner_index] = _RootSegment(
            owner_index,
            local_id,
            outward_node,
            int(port.global_edge_id),
            length,
        )
        self.extension_reasons[owner_index] = reason
        self.frontiers.append(_Frontier(owner_index, outward_node, int(port.global_edge_id)))
        self._expand_bifurcation_frontiers(owner_index)

    def _expand_frontier(self, frontier: _Frontier) -> list[_Frontier]:
        candidates: list[tuple[int, int]] = []
        for neighbor in sorted(self.graph.neighbors(frontier.global_node_id)):
            edge_id = int(self.graph.edges[frontier.global_node_id, neighbor]["global_edge_id"])
            if edge_id == frontier.incoming_global_edge_id or edge_id in self.added_full_edges:
                continue
            candidates.append((int(neighbor), edge_id))
        if not candidates:
            return [frontier]
        following: list[_Frontier] = []
        for neighbor, edge_id in candidates:
            self.added_full_edges[edge_id] = frontier.owner_index
            following.append(_Frontier(frontier.owner_index, neighbor, edge_id))
        return following

    def _expand_bifurcation_frontiers(self, owner_index: int) -> None:
        for _ in range(64):
            changed = False
            following: list[_Frontier] = []
            for frontier in self.frontiers:
                if (
                    frontier.owner_index == owner_index
                    and self.graph.degree(frontier.global_node_id) >= 3
                ):
                    following.extend(self._expand_frontier(frontier))
                    changed = True
                else:
                    following.append(frontier)
            self.frontiers = following
            if not changed:
                return
        raise GeometryValidationError("CFD context encountered an unresolved bifurcation chain")

    def expand_owner(self, owner_index: int) -> None:
        selected = [item for item in self.frontiers if item.owner_index == owner_index]
        if not selected:
            raise GeometryValidationError(
                f"CFD context for CORE port {owner_index} reached no expandable frontier"
            )
        retained = [item for item in self.frontiers if item.owner_index != owner_index]
        following: list[_Frontier] = []
        for frontier in selected:
            following.extend(self._expand_frontier(frontier))
        if following == selected:
            raise GeometryValidationError(
                f"CFD context for CORE port {owner_index} reached a true global terminal"
            )
        self.frontiers = [*retained, *following]
        self._expand_bifurcation_frontiers(owner_index)
        if len(self.added_full_edges) + len(self.root_segments) > (
            self.config.context_domain.max_added_global_edges
        ):
            raise GeometryValidationError("CFD_CONTEXT_MAX_EXPANSION_EXCEEDED")

    def build_roi(
        self,
    ) -> tuple[ROIRecord, dict[int, int], dict[int, int]]:
        positions = [np.asarray(point, dtype=float).copy() for point in self.core.local_node_positions_um]
        radii = list(map(float, self.core.local_node_radius_um))
        global_ids = list(map(int, self.core.local_node_global_ids))
        global_to_local = {
            int(global_id): local_id
            for local_id, global_id in enumerate(global_ids)
            if int(global_id) >= 0
        }

        def ensure_node(global_node_id: int) -> int:
            if global_node_id in global_to_local:
                return global_to_local[global_node_id]
            local_id = len(positions)
            source_index = self.model.node_index_by_id[global_node_id]
            positions.append(np.asarray(self.model.node_positions_um[source_index], dtype=float).copy())
            radii.append(float(self.model.node_radius_um[source_index]))
            global_ids.append(global_node_id)
            global_to_local[global_node_id] = local_id
            return local_id

        edges = [tuple(map(int, edge)) for edge in np.asarray(self.core.local_edges)]
        edge_global_ids = list(map(int, self.core.local_edge_global_ids))
        edge_points = [np.asarray(points, dtype=float).copy() for points in self.core.local_edge_points_um]
        edge_radii = [np.asarray(values, dtype=float).copy() for values in self.core.local_edge_radius_um]
        for root in self.root_segments.values():
            outer_local = ensure_node(root.outward_global_node_id)
            edges.append((root.core_local_node_id, outer_local))
            edge_global_ids.append(root.global_edge_id)
            edge_points.append(np.asarray((positions[root.core_local_node_id], positions[outer_local])))
            edge_radii.append(np.asarray((radii[root.core_local_node_id], radii[outer_local])))
        for edge_id in sorted(self.added_full_edges):
            source = self.model.edges[edge_id]
            first = ensure_node(source.upstream_node_id)
            second = ensure_node(source.downstream_node_id)
            if tuple(sorted((first, second))) in {tuple(sorted(edge)) for edge in edges}:
                continue
            edges.append((first, second))
            edge_global_ids.append(edge_id)
            edge_points.append(np.asarray((positions[first], positions[second])))
            edge_radii.append(np.asarray((radii[first], radii[second])))

        frontier_by_local: dict[int, int] = {}
        final_ports: list[CutPort] = []
        for index, port in enumerate(self.core.cut_ports):
            if index in self.root_segments:
                continue
            final_ports.append(
                CutPort(
                    cut_port_id=port.cut_port_id,
                    local_node_id=port.local_node_id,
                    global_edge_id=port.global_edge_id,
                    intersection_position_um=port.intersection_position_um,
                    radius_at_cut_um=port.radius_at_cut_um,
                    boundary_face=port.boundary_face,
                    boundary_role="CORE_ROI_BOUNDARY|CFD_BOUNDARY_PORT",
                )
            )
        true_terminal_local = list(map(int, self.core.true_terminal_local_ids))
        true_terminal_global = list(map(int, self.core.true_terminal_global_ids))
        for ordinal, frontier in enumerate(
            sorted(self.frontiers, key=lambda item: (item.owner_index, item.global_node_id))
        ):
            local_id = ensure_node(frontier.global_node_id)
            if self.graph.degree(frontier.global_node_id) == 1:
                true_terminal_local.append(local_id)
                true_terminal_global.append(frontier.global_node_id)
                continue
            outward_edges = [
                int(self.graph.edges[frontier.global_node_id, neighbor]["global_edge_id"])
                for neighbor in self.graph.neighbors(frontier.global_node_id)
                if int(self.graph.edges[frontier.global_node_id, neighbor]["global_edge_id"])
                != frontier.incoming_global_edge_id
            ]
            global_edge_id = outward_edges[0] if outward_edges else frontier.incoming_global_edge_id
            final_ports.append(
                CutPort(
                    cut_port_id=(
                        f"{self.core.roi_id}__cfd_cut_{frontier.owner_index:03d}_{ordinal:03d}"
                    ),
                    local_node_id=local_id,
                    global_edge_id=global_edge_id,
                    intersection_position_um=tuple(map(float, positions[local_id])),
                    radius_at_cut_um=float(radii[local_id]),
                    boundary_face="context",
                    boundary_role="CFD_BOUNDARY_PORT",
                )
            )
            frontier_by_local[local_id] = frontier.owner_index

        position_array = np.asarray(positions, dtype=float)
        radius_array = np.asarray(radii, dtype=float)
        edge_array = np.asarray(edges, dtype=np.int64).reshape((-1, 2))
        minimum = np.minimum(np.asarray(self.core.bbox_min_um), position_array.min(axis=0))
        maximum = np.maximum(np.asarray(self.core.bbox_max_um), position_array.max(axis=0))
        lengths = np.linalg.norm(position_array[edge_array[:, 1]] - position_array[edge_array[:, 0]], axis=1)
        roi = ROIRecord(
            roi_id=f"{self.core.roi_id}__cfd_domain",
            source_model_id=self.core.source_model_id,
            source_mouse_id=self.core.source_mouse_id,
            anchor_id=self.core.anchor_id,
            anchor_position_um=self.core.anchor_position_um,
            bbox_min_um=tuple(map(float, minimum)),
            bbox_max_um=tuple(map(float, maximum)),
            bbox_center_um=tuple(map(float, 0.5 * (minimum + maximum))),
            bbox_size_um=tuple(map(float, maximum - minimum)),
            global_node_ids=tuple(sorted(set(global_id for global_id in global_ids if global_id >= 0))),
            global_edge_ids=tuple(sorted(set(edge_global_ids))),
            local_node_ids=np.arange(len(positions), dtype=np.int64),
            local_node_global_ids=np.asarray(global_ids, dtype=np.int64),
            local_node_positions_um=position_array,
            local_node_radius_um=radius_array,
            local_edges=edge_array,
            local_edge_ids=np.arange(len(edges), dtype=np.int64),
            local_edge_global_ids=np.asarray(edge_global_ids, dtype=np.int64),
            local_edge_points_um=np.asarray(edge_points, dtype=float),
            local_edge_radius_um=np.asarray(edge_radii, dtype=float),
            true_terminal_local_ids=tuple(sorted(set(true_terminal_local))),
            true_terminal_global_ids=tuple(sorted(set(true_terminal_global))),
            cut_ports=tuple(final_ports),
            raw_component_count=self.core.raw_component_count,
            raw_total_vessel_length_um=self.core.raw_total_vessel_length_um,
            retained_component_length_um=float(lengths.sum()),
        )
        final_port_by_local = {
            int(port.local_node_id): index for index, port in enumerate(roi.cut_ports)
        }
        return roi, frontier_by_local, final_port_by_local

    def owner_added_edges(self, owner_index: int) -> set[int]:
        output = {
            root.global_edge_id
            for index, root in self.root_segments.items()
            if index == owner_index
        }
        output.update(
            edge_id for edge_id, owner in self.added_full_edges.items() if owner == owner_index
        )
        return output

    def owner_added_length(self, owner_index: int) -> float:
        total = (
            self.root_segments[owner_index].length_um
            if owner_index in self.root_segments
            else 0.0
        )
        total += sum(
            float(np.linalg.norm(
                self.model.edges[edge_id].downstream_position_um
                - self.model.edges[edge_id].upstream_position_um
            ))
            for edge_id, owner in self.added_full_edges.items()
            if owner == owner_index
        )
        return total


def build_cfd_context_domain(
    core_roi: ROIRecord,
    sampling_run: Path,
    config: CFDLumenConfig,
) -> ContextDomainResult:
    """Follow real global SWC paths until every final CFD port passes v3 collar checks."""

    signature_before = _core_signature(core_roi)
    if not config.context_domain.enabled:
        return ContextDomainResult(
            core_roi=core_roi,
            cfd_roi=core_roi,
            source_manifest=None,
            port_mappings=[],
            added_global_edge_ids=(),
            added_global_branch_ids=(),
            added_centerline_length_um=0.0,
            report={"status": "DISABLED", "core_roi_signature": signature_before},
        )
    model, manifest_path = _load_verified_model(sampling_run, core_roi, config)
    builder = _DomainBuilder(core_roi, model, config)
    conflict_history: list[dict[str, Any]] = []
    final_roi: ROIRecord | None = None
    for iteration in range(config.context_domain.max_added_global_edges + 1):
        domain, frontier_by_local, final_port_by_local = builder.build_roi()
        branches, _ = validate_and_extract_branches(domain, config)
        resample_branches(branches, config)
        try:
            define_junction_collars(domain, branches, config)
            final_roi = domain
            break
        except GeometryValidationError as exc:
            code = getattr(exc, "code", None)
            branch_id = getattr(exc, "branch_id", None)
            if code not in {
                "JUNCTION_PORT_REGION_CONFLICT",
                "JUNCTION_COLLAR_TOPOLOGY_CONFLICT",
                "JUNCTION_COLLAR_SEPARATION_FAILED",
            }:
                raise
            candidate_local_ids: list[int] = []
            if branch_id is not None:
                branch = next((item for item in branches if item.branch_id == branch_id), None)
                if branch is not None:
                    candidate_local_ids = [
                        node_id
                        for node_id in (branch.local_node_ids[0], branch.local_node_ids[-1])
                        if node_id in final_port_by_local
                    ]
            if not candidate_local_ids:
                junction_id = getattr(exc, "junction_node_id", None)
                for branch in branches:
                    if junction_id in {branch.local_node_ids[0], branch.local_node_ids[-1]}:
                        candidate_local_ids.extend(
                            node_id
                            for node_id in (branch.local_node_ids[0], branch.local_node_ids[-1])
                            if node_id in final_port_by_local
                        )
            if not candidate_local_ids:
                raise GeometryValidationError(
                    f"{code} cannot be resolved by extending a CFD boundary port: {exc}"
                ) from exc
            owners: set[int] = set()
            for local_id in candidate_local_ids:
                port = domain.cut_ports[final_port_by_local[local_id]]
                if local_id in frontier_by_local:
                    owners.add(frontier_by_local[local_id])
                else:
                    owner = next(
                        index
                        for index, original in enumerate(core_roi.cut_ports)
                        if original.cut_port_id == port.cut_port_id
                    )
                    owners.add(owner)
            conflict_history.append(
                {
                    "iteration": iteration,
                    "code": code,
                    "message": str(exc),
                    "owners_extended": sorted(owners),
                }
            )
            for owner in sorted(owners):
                if owner in builder.root_segments:
                    builder.expand_owner(owner)
                else:
                    builder.start_extension(owner, str(exc))
    if final_roi is None:
        raise GeometryValidationError("CFD_CONTEXT_MAX_EXPANSION_EXCEEDED")

    branch_map = _global_branch_map(builder.graph)
    core_distances = _core_port_distances(core_roi)
    mappings: list[dict[str, Any]] = []
    for owner, original in enumerate(core_roi.cut_ports):
        final_for_owner = [
            port
            for port in final_roi.cut_ports
            if (
                port.cut_port_id == original.cut_port_id
                or port.cut_port_id.startswith(f"{core_roi.roi_id}__cfd_cut_{owner:03d}_")
            )
        ]
        added_edges = builder.owner_added_edges(owner)
        mappings.append(
            {
                "original_cut_port": original.report(),
                "new_cfd_ports": [port.report() for port in final_for_owner],
                "original_distance_to_nearest_junction_um": core_distances[owner],
                "added_centerline_length_um": builder.owner_added_length(owner),
                "added_global_edge_ids": sorted(added_edges),
                "added_branch_ids": sorted(
                    {branch_map[edge_id] for edge_id in added_edges if edge_id in branch_map}
                ),
                "reason_for_extension": builder.extension_reasons.get(owner, "original port passed"),
            }
        )
    added_edges = set(builder.added_full_edges)
    added_edges.update(root.global_edge_id for root in builder.root_segments.values())
    signature_after = _core_signature(core_roi)
    if signature_after != signature_before:
        raise GeometryValidationError("CORE_ROI_MUTATION_DETECTED")
    added_length = sum(builder.owner_added_length(owner) for owner in builder.root_segments)
    report = {
        "status": "PASS",
        "sampling_feature_cluster_representative_domain": "CORE_ROI",
        "mb_ulm_statistics_default_domain": "CORE_ROI",
        "cfd_geometry_and_flow_domain": "CFD_DOMAIN",
        "core_roi_id": core_roi.roi_id,
        "cfd_domain_id": final_roi.roi_id,
        "core_roi_signature_before": signature_before,
        "core_roi_signature_after": signature_after,
        "core_roi_unchanged": signature_before == signature_after,
        "original_cut_ports_retained_in_core_record": True,
        "source_manifest": str(manifest_path),
        "source_global_node_count": model.node_count,
        "source_global_edge_count": model.edge_count,
        "core_node_count": core_roi.node_count,
        "cfd_node_count": final_roi.node_count,
        "core_edge_count": core_roi.edge_count,
        "cfd_edge_count": final_roi.edge_count,
        "core_cut_port_count": len(core_roi.cut_ports),
        "final_cfd_port_count": len(final_roi.cut_ports),
        "added_centerline_length_um": added_length,
        "added_global_edge_ids": sorted(added_edges),
        "added_global_branch_ids": sorted(
            {branch_map[edge_id] for edge_id in added_edges if edge_id in branch_map}
        ),
        "conflict_history": conflict_history,
        "final_JUNCTION_PORT_REGION_CONFLICT_count": 0,
        "port_mappings": mappings,
    }
    return ContextDomainResult(
        core_roi=core_roi,
        cfd_roi=final_roi,
        source_manifest=manifest_path,
        port_mappings=mappings,
        added_global_edge_ids=tuple(sorted(added_edges)),
        added_global_branch_ids=tuple(report["added_global_branch_ids"]),
        added_centerline_length_um=float(added_length),
        report=report,
    )
