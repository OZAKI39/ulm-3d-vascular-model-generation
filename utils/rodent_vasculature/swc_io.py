"""Strict SWC parsing with parent-to-current-node flow semantics."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np


class EmptySWCError(ValueError):
    """Raised when an SWC has no usable records."""


@dataclass(slots=True)
class SWCData:
    path: Path
    node_ids: np.ndarray
    type_codes: np.ndarray
    points_voxel_xyz: np.ndarray
    points_um: np.ndarray
    radius_raw_um: np.ndarray
    parent_ids: np.ndarray
    root_ids: list[int]
    component_count: int
    validation: dict[str, Any]

    @property
    def node_count(self) -> int:
        return len(self.node_ids)

    @property
    def edge_count(self) -> int:
        return int(np.count_nonzero(self.parent_ids != -1))

    @property
    def structurally_valid(self) -> bool:
        return not any(
            self.validation[key]
            for key in (
                "malformed_lines",
                "duplicate_node_ids",
                "missing_parent_ids",
                "self_parent_node_ids",
                "directed_cycle_node_ids",
                "noninteger_identifier_lines",
                "nonfinite_coordinate_row_indices",
            )
        ) and bool(self.validation["root_component_count_matches"])

    def directed_graph(self) -> nx.DiGraph:
        graph = nx.DiGraph()
        for index, node_id in enumerate(self.node_ids.tolist()):
            graph.add_node(
                int(node_id),
                source_index=index,
                type_code=int(self.type_codes[index]),
            )
        for node_id, parent_id in zip(self.node_ids.tolist(), self.parent_ids.tolist()):
            if parent_id != -1 and parent_id in graph:
                graph.add_edge(int(parent_id), int(node_id), relationship="parent_to_child")
        return graph

    def report(self) -> dict[str, Any]:
        return {
            "path": str(self.path),
            "node_count": self.node_count,
            "edge_count": self.edge_count,
            "root_count": len(self.root_ids),
            "root_ids": self.root_ids,
            "component_count": self.component_count,
            "flow_direction_rule": "parent_id node -> current node",
            "flow_direction_is_measured": False,
            "structurally_valid": self.structurally_valid,
            "radius_valid_fraction": float(
                np.mean(np.isfinite(self.radius_raw_um) & (self.radius_raw_um > 0))
            ),
            "validation": self.validation,
        }


def _integer(value: str, line_number: int, field: str, errors: list[int]) -> int:
    parsed = float(value)
    rounded = int(round(parsed))
    if not np.isfinite(parsed) or not np.isclose(parsed, rounded, rtol=0.0, atol=1.0e-9):
        errors.append(line_number)
        raise ValueError(f"SWC {field} is not integer-valued at line {line_number}: {value}")
    return rounded


def swc_from_arrays(
    *,
    path: Path,
    node_ids: np.ndarray,
    type_codes: np.ndarray,
    points_voxel_xyz: np.ndarray,
    radius_raw_um: np.ndarray,
    parent_ids: np.ndarray,
    spacing_xyz_um: tuple[float, float, float],
    volume_shape_zyx: tuple[int, int, int] | None,
    parser_validation: dict[str, Any] | None = None,
) -> SWCData:
    """Build and validate canonical SWC data from in-memory arrays.

    This is shared by the text parser, analysis-component selection, and the
    normalized NPZ loader so downstream stages consume an explicitly saved
    SWC data object instead of silently changing data sources.
    """

    node_ids = np.asarray(node_ids, dtype=np.int64).reshape(-1)
    type_codes = np.asarray(type_codes, dtype=np.int32).reshape(-1)
    points_voxel = np.asarray(points_voxel_xyz, dtype=np.float64).reshape((-1, 3))
    radii = np.asarray(radius_raw_um, dtype=np.float64).reshape(-1)
    parent_ids = np.asarray(parent_ids, dtype=np.int64).reshape(-1)
    lengths = {
        len(node_ids),
        len(type_codes),
        len(points_voxel),
        len(radii),
        len(parent_ids),
    }
    if len(lengths) != 1 or not node_ids.size:
        raise ValueError("SWC arrays must be non-empty and have matching row counts")

    unique_ids, counts = np.unique(node_ids, return_counts=True)
    duplicates = unique_ids[counts > 1].astype(int).tolist()
    id_set = set(node_ids.tolist())
    missing_parents = sorted(
        {int(parent) for parent in parent_ids.tolist() if parent != -1 and parent not in id_set}
    )
    self_parents = node_ids[node_ids == parent_ids].astype(int).tolist()

    graph = nx.DiGraph()
    graph.add_nodes_from(int(value) for value in node_ids.tolist())
    graph.add_edges_from(
        (int(parent), int(node))
        for node, parent in zip(node_ids.tolist(), parent_ids.tolist())
        if parent != -1 and parent in id_set and node != parent
    )
    try:
        cycle_edges = nx.find_cycle(graph, orientation="original")
    except nx.NetworkXNoCycle:
        cycle_nodes: list[int] = []
    else:
        cycle_nodes = sorted(
            {int(edge[0]) for edge in cycle_edges} | {int(edge[1]) for edge in cycle_edges}
        )
    roots = sorted(int(value) for value in node_ids[parent_ids == -1].tolist())
    component_count = nx.number_weakly_connected_components(graph) if graph else 0
    nonfinite_coordinate_rows = np.flatnonzero(
        ~np.all(np.isfinite(points_voxel), axis=1)
    ).tolist()
    outside_rows: list[int] = []
    if volume_shape_zyx is not None:
        maximum_xyz = np.asarray(volume_shape_zyx[::-1], dtype=float)
        outside_rows = np.flatnonzero(
            np.any((points_voxel < -1.0e-6) | (points_voxel > maximum_xyz + 1.0e-6), axis=1)
        ).tolist()
    nonpositive_radius_nodes = node_ids[
        ~np.isfinite(radii) | (radii <= 0)
    ].astype(int).tolist()
    parser_validation = parser_validation or {}
    validation = {
        "malformed_lines": list(parser_validation.get("malformed_lines", [])),
        "noninteger_identifier_lines": list(
            parser_validation.get("noninteger_identifier_lines", [])
        ),
        "duplicate_node_ids": duplicates,
        "missing_parent_ids": missing_parents,
        "self_parent_node_ids": self_parents,
        "directed_cycle_node_ids": cycle_nodes,
        "nonfinite_coordinate_row_indices": nonfinite_coordinate_rows,
        "outside_volume_row_indices": outside_rows,
        "nonpositive_radius_node_ids": nonpositive_radius_nodes,
        "root_component_count_matches": len(roots) == component_count,
    }
    return SWCData(
        path=path,
        node_ids=node_ids,
        type_codes=type_codes,
        points_voxel_xyz=points_voxel,
        points_um=points_voxel * np.asarray(spacing_xyz_um, dtype=float),
        radius_raw_um=radii,
        parent_ids=parent_ids,
        root_ids=roots,
        component_count=component_count,
        validation=validation,
    )


def load_swc(
    path: Path,
    *,
    spacing_xyz_um: tuple[float, float, float],
    volume_shape_zyx: tuple[int, int, int] | None,
) -> SWCData:
    rows: list[tuple[int, int, float, float, float, float, int]] = []
    noninteger_lines: list[int] = []
    malformed_lines: list[int] = []
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8", errors="replace").splitlines(), 1
    ):
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 7:
            malformed_lines.append(line_number)
            continue
        try:
            node_id = _integer(parts[0], line_number, "node id", noninteger_lines)
            type_code = _integer(parts[1], line_number, "type code", noninteger_lines)
            parent_id = _integer(parts[6], line_number, "parent id", noninteger_lines)
            x, y, z, radius = (float(parts[index]) for index in (2, 3, 4, 5))
        except ValueError:
            if line_number not in malformed_lines:
                malformed_lines.append(line_number)
            continue
        rows.append((node_id, type_code, x, y, z, radius, parent_id))
    if not rows:
        raise EmptySWCError(f"SWC contains no usable records: {path}")

    array = np.asarray(rows, dtype=float)
    return swc_from_arrays(
        path=path,
        node_ids=array[:, 0],
        type_codes=array[:, 1],
        points_voxel_xyz=array[:, 2:5],
        radius_raw_um=array[:, 5],
        parent_ids=array[:, 6],
        spacing_xyz_um=spacing_xyz_um,
        volume_shape_zyx=volume_shape_zyx,
        parser_validation={
            "malformed_lines": malformed_lines,
            "noninteger_identifier_lines": noninteger_lines,
        },
    )


def save_normalized_swc(data: SWCData, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        node_ids=data.node_ids,
        type_codes=data.type_codes,
        points_voxel_xyz=data.points_voxel_xyz,
        points_um=data.points_um,
        radius_raw_um=data.radius_raw_um,
        radius_valid=np.isfinite(data.radius_raw_um) & (data.radius_raw_um > 0),
        parent_ids=data.parent_ids,
        root_ids=np.asarray(data.root_ids, dtype=np.int64),
        component_count=np.asarray(data.component_count, dtype=np.int64),
        validation_json=np.asarray(json.dumps(data.validation, ensure_ascii=False)),
    )
    return path


def save_swc_text(data: SWCData, path: Path) -> Path:
    """Write a derived SWC representation while preserving node identifiers."""

    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# validated SWC derivative; columns: id type x_voxel y_voxel z_voxel radius_um parent_id",
        "# structural direction: parent_id node -> current node; not measured blood flow",
    ]
    for node_id, type_code, point, radius, parent_id in zip(
        data.node_ids,
        data.type_codes,
        data.points_voxel_xyz,
        data.radius_raw_um,
        data.parent_ids,
    ):
        lines.append(
            f"{int(node_id)} {int(type_code)} "
            f"{float(point[0]):.9g} {float(point[1]):.9g} {float(point[2]):.9g} "
            f"{float(radius):.9g} {int(parent_id)}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def load_normalized_swc(
    path: Path,
    *,
    source_path: Path,
    spacing_xyz_um: tuple[float, float, float],
    volume_shape_zyx: tuple[int, int, int] | None,
) -> SWCData:
    with np.load(path, allow_pickle=False) as arrays:
        node_ids = arrays["node_ids"].astype(np.int64)
        type_codes = arrays["type_codes"].astype(np.int32)
        points_voxel = arrays["points_voxel_xyz"].astype(np.float64)
        radii = arrays["radius_raw_um"].astype(np.float64)
        parent_ids = arrays["parent_ids"].astype(np.int64)
        validation_json = str(arrays["validation_json"].item()) if "validation_json" in arrays else "{}"
    try:
        stored_validation = json.loads(validation_json)
    except (TypeError, json.JSONDecodeError):
        stored_validation = {}
    return swc_from_arrays(
        path=source_path,
        node_ids=node_ids,
        type_codes=type_codes,
        points_voxel_xyz=points_voxel,
        radius_raw_um=radii,
        parent_ids=parent_ids,
        spacing_xyz_um=spacing_xyz_um,
        volume_shape_zyx=volume_shape_zyx,
        parser_validation=stored_validation,
    )
