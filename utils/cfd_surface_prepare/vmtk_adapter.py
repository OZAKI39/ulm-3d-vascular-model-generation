"""Build the small VMTK centerline adapter from saved SWC-derived 1-D data."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pyvista as pv

from .io import BoundaryInput, SurfacePrepareError


@dataclass(frozen=True, slots=True)
class CenterlineAdapterResult:
    path: Path
    records: tuple[dict[str, Any], ...]


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise SurfacePrepareError(f"Missing SWC-derived 1-D artifact: {path}")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def _point(row: dict[str, str]) -> np.ndarray:
    return np.asarray([float(row[f"{axis}_um"]) for axis in "xyz"], dtype=float)


def _unique_child(
    node_id: int,
    children: dict[int, list[int]],
    visited: set[int],
) -> int | None:
    available = [child for child in children.get(node_id, []) if child not in visited]
    return available[0] if len(available) == 1 else None


def _expand_actual_path(
    parent_id: int,
    child_id: int,
    *,
    nodes: dict[int, dict[str, str]],
    children: dict[int, list[int]],
    target_length_um: float,
) -> list[int]:
    """Extend along unique real graph neighbours; never interpolate a centreline."""

    path = [parent_id, child_id]
    visited = set(path)

    def length() -> float:
        points = np.asarray([_point(nodes[node_id]) for node_id in path])
        return float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())

    while length() < target_length_um:
        changed = False
        parent_text = nodes[path[0]].get("parent_id", "").strip()
        if parent_text and int(parent_text) >= 0 and int(parent_text) not in visited:
            ancestor = int(parent_text)
            if ancestor in nodes:
                path.insert(0, ancestor)
                visited.add(ancestor)
                changed = True
        if length() >= target_length_um:
            break
        descendant = _unique_child(path[-1], children, visited)
        if descendant is not None:
            path.append(descendant)
            visited.add(descendant)
            changed = True
        if not changed:
            break
    return path


def build_centerline_adapter(
    preprocess_run: Path,
    boundaries: Iterable[BoundaryInput],
    output_path: Path,
) -> CenterlineAdapterResult:
    """Write four local polylines composed only of saved global-1D/SWC nodes."""

    root = Path(preprocess_run)
    node_rows = _read_csv(root / "global_1d" / "nodes.csv")
    edge_rows = _read_csv(root / "global_1d" / "edges.csv")
    nodes = {int(row["node_id"]): row for row in node_rows}
    edges = {int(row["global_edge_id"]): row for row in edge_rows}
    children: dict[int, list[int]] = {}
    for row in node_rows:
        node_id = int(row["node_id"])
        parent_text = row.get("parent_id", "").strip()
        if parent_text and int(parent_text) >= 0:
            children.setdefault(int(parent_text), []).append(node_id)

    all_points: list[np.ndarray] = []
    line_cells: list[np.ndarray] = []
    point_boundary_index: list[int] = []
    point_radius_um: list[float] = []
    records: list[dict[str, Any]] = []
    boundary_list = list(boundaries)
    for boundary in boundary_list:
        if boundary.global_edge_id not in edges:
            raise SurfacePrepareError(
                f"Missing global edge {boundary.global_edge_id} for {boundary.port_id}"
            )
        edge = edges[boundary.global_edge_id]
        parent_id = int(edge["parent_node_id"])
        child_id = int(edge["child_node_id"])
        if parent_id not in nodes or child_id not in nodes:
            raise SurfacePrepareError(
                f"Global edge endpoints unavailable for {boundary.port_id}"
            )
        edge_vector = _point(nodes[child_id]) - _point(nodes[parent_id])
        edge_length = float(np.linalg.norm(edge_vector))
        if edge_length <= 0.0:
            raise SurfacePrepareError(f"Zero-length source edge for {boundary.port_id}")
        edge_direction = edge_vector / edge_length
        tangent_alignment = float(abs(np.dot(edge_direction, boundary.simulation_tangent)))
        if tangent_alignment < 0.999:
            raise SurfacePrepareError(
                f"SWC_CENTERLINE_DIRECTION_MISMATCH:{boundary.port_id}"
            )
        # Approximately two diameters are requested, but graph expansion stops at
        # a branch or terminal rather than inventing a continuation.
        target_length = 4.0 * boundary.source_radius_um
        node_path = _expand_actual_path(
            parent_id,
            child_id,
            nodes=nodes,
            children=children,
            target_length_um=target_length,
        )
        points = np.asarray([_point(nodes[node_id]) for node_id in node_path])
        arc_length = float(np.linalg.norm(np.diff(points, axis=0), axis=1).sum())
        start = len(all_points)
        all_points.extend(points)
        line_cells.append(np.arange(start, start + len(points), dtype=np.int64))
        point_boundary_index.extend([boundary.index] * len(points))
        point_radius_um.extend(float(nodes[node_id]["radius_um"]) for node_id in node_path)
        segment_start = _point(nodes[parent_id])
        projection = float(np.dot(boundary.center_um - segment_start, edge_direction))
        projected = segment_start + np.clip(projection, 0.0, edge_length) * edge_direction
        projection_error = float(np.linalg.norm(projected - boundary.center_um))
        records.append(
            {
                "boundary_index": boundary.index,
                "port_id": boundary.port_id,
                "global_edge_id": boundary.global_edge_id,
                "source_parent_node_id": parent_id,
                "source_child_node_id": child_id,
                "source_node_ids": node_path,
                "source_point_count": len(node_path),
                "source_arc_length_um": arc_length,
                "target_local_centerline_length_um": target_length,
                "source_arc_length_diameters": arc_length
                / (2.0 * boundary.source_radius_um),
                "boundary_to_source_edge_projection_error_um": projection_error,
                "source_edge_direction_abs_dot_simulation_tangent": tangent_alignment,
                "invented_centerline_points": 0,
                "source": "saved global_1d nodes/edges derived from source SWC",
            }
        )

    cells: list[int] = []
    for line in line_cells:
        cells.extend([len(line), *map(int, line)])
    polydata = pv.PolyData()
    polydata.points = np.asarray(all_points, dtype=float)
    polydata.lines = np.asarray(cells, dtype=np.int64)
    polydata.point_data["boundary_index"] = np.asarray(
        point_boundary_index, dtype=np.int32
    )
    polydata.point_data["MaximumInscribedSphereRadius"] = np.asarray(
        point_radius_um, dtype=float
    )
    polydata.cell_data["boundary_index"] = np.asarray(
        [item.index for item in boundary_list], dtype=np.int32
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    polydata.save(output_path, binary=True)
    return CenterlineAdapterResult(output_path.resolve(), tuple(records))
