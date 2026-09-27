"""Machine-readable exports for directed vascular graphs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import vtk

from ..io import write_csv, write_json
from .model import DirectedVascularGraph


def _graphml_value(value: Any) -> Any:
    if isinstance(value, (str, int, float, bool)):
        return value
    if value is None:
        return ""
    if isinstance(value, np.generic):
        return value.item()
    return ";".join(map(str, value)) if isinstance(value, (list, tuple, set)) else str(value)


def _graphml_copy(graph: nx.Graph) -> nx.Graph:
    output = graph.__class__()
    for node, attributes in graph.nodes(data=True):
        output.add_node(node, **{key: _graphml_value(value) for key, value in attributes.items()})
    if graph.is_multigraph():
        for upstream, downstream, key, attributes in graph.edges(keys=True, data=True):
            output.add_edge(
                upstream,
                downstream,
                key=key,
                **{name: _graphml_value(value) for name, value in attributes.items()},
            )
    else:
        for upstream, downstream, attributes in graph.edges(data=True):
            output.add_edge(
                upstream,
                downstream,
                **{name: _graphml_value(value) for name, value in attributes.items()},
            )
    return output


def _write_vtp(result: DirectedVascularGraph, path: Path) -> Path:
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    branch_id_array = vtk.vtkIntArray()
    branch_id_array.SetName("branch_id")
    radius_array = vtk.vtkDoubleArray()
    radius_array.SetName("radius_um")
    direction_array = vtk.vtkDoubleArray()
    direction_array.SetName("parent_to_current_direction_xyz")
    direction_array.SetNumberOfComponents(3)
    depth_array = vtk.vtkIntArray()
    depth_array.SetName("branch_depth")
    upstream_array = vtk.vtkLongLongArray()
    upstream_array.SetName("upstream_swc_node_id")
    downstream_array = vtk.vtkLongLongArray()
    downstream_array.SetName("downstream_swc_node_id")

    for branch in result.branches:
        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(len(branch.derived_points_um))
        for local_index, (point, radius, direction) in enumerate(
            zip(branch.derived_points_um, branch.derived_radius_um, branch.direction_vectors_xyz)
        ):
            point_id = points.InsertNextPoint(*(float(value) for value in point))
            polyline.GetPointIds().SetId(local_index, point_id)
            branch_id_array.InsertNextValue(branch.branch_id)
            radius_array.InsertNextValue(float(radius))
            direction_array.InsertNextTuple3(*(float(value) for value in direction))
            depth_array.InsertNextValue(branch.depth)
        lines.InsertNextCell(polyline)
        upstream_array.InsertNextValue(branch.upstream_node_id)
        downstream_array.InsertNextValue(branch.downstream_node_id)

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    for array in (branch_id_array, radius_array, direction_array, depth_array):
        polydata.GetPointData().AddArray(array)
    for array in (upstream_array, downstream_array):
        polydata.GetCellData().AddArray(array)
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(polydata)
    writer.SetDataModeToBinary()
    if writer.Write() != 1:
        raise OSError(f"Failed to write VTP: {path}")
    return path


def export_directed_graph(
    result: DirectedVascularGraph,
    output_dir: Path,
    *,
    save_graphml: bool,
    save_vtp: bool,
    save_npz: bool,
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    branch_table = output_dir / "directed_branches.csv"
    write_csv([branch.summary() for branch in result.branches], branch_table)
    paths.append(branch_table)
    node_table = output_dir / "directed_nodes.csv"
    write_csv(
        [dict(node_id=node_id, **attributes) for node_id, attributes in result.junction_graph.nodes(data=True)],
        node_table,
    )
    paths.append(node_table)
    edge_table = output_dir / "source_parent_to_current_edges.csv"
    write_csv(
        [
            {"parent_node_id": upstream, "current_node_id": downstream, "direction": "parent_to_current"}
            for upstream, downstream in result.source_graph.edges
        ],
        edge_table,
    )
    paths.append(edge_table)
    summary_path = output_dir / "directed_graph_summary.json"
    write_json(result.report(), summary_path)
    paths.append(summary_path)

    if save_graphml:
        junction_path = output_dir / "junction_graph_parent_to_current.graphml"
        branch_path = output_dir / "branch_hierarchy_parent_to_current.graphml"
        nx.write_graphml(_graphml_copy(result.junction_graph), junction_path)
        branch_attributes = nx.DiGraph()
        for branch in result.branches:
            branch_attributes.add_node(branch.branch_id, **branch.summary())
        for upstream, downstream, attributes in result.branch_graph.edges(data=True):
            branch_attributes.add_edge(upstream, downstream, **attributes)
        nx.write_graphml(_graphml_copy(branch_attributes), branch_path)
        paths.extend((junction_path, branch_path))
    if save_npz:
        offsets = [0]
        all_points: list[np.ndarray] = []
        all_directions: list[np.ndarray] = []
        all_radii: list[np.ndarray] = []
        for branch in result.branches:
            all_points.append(branch.derived_points_um)
            all_directions.append(branch.direction_vectors_xyz)
            all_radii.append(branch.derived_radius_um)
            offsets.append(offsets[-1] + len(branch.derived_points_um))
        npz_path = output_dir / "directed_branch_geometry.npz"
        np.savez_compressed(
            npz_path,
            points_um=np.concatenate(all_points) if all_points else np.empty((0, 3)),
            direction_parent_to_current_xyz=(
                np.concatenate(all_directions) if all_directions else np.empty((0, 3))
            ),
            radius_um=np.concatenate(all_radii) if all_radii else np.empty(0),
            branch_offsets=np.asarray(offsets, dtype=np.int64),
            branch_ids=np.asarray([branch.branch_id for branch in result.branches]),
            upstream_node_ids=np.asarray([branch.upstream_node_id for branch in result.branches]),
            downstream_node_ids=np.asarray([branch.downstream_node_id for branch in result.branches]),
        )
        paths.append(npz_path)
    if save_vtp:
        paths.append(_write_vtp(result, output_dir / "directed_branches_parent_to_current.vtp"))
    return paths
