"""Portable exports for hierarchical vascular graph artifacts."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import networkx as nx
import numpy as np
import vtk

from ..io import write_csv, write_json
from .model import HierarchicalGraphResult


def _graph_payload(result: HierarchicalGraphResult) -> dict[str, Any]:
    relations = []
    for branch_a, branch_b, key, data in result.branch_as_node_graph.edges(
        keys=True, data=True
    ):
        relations.append(
            {
                "relation_id": str(key),
                "branch_a": int(branch_a),
                "branch_b": int(branch_b),
                **dict(data),
            }
        )
    return {
        "schema_version": "1.0",
        "representation": {
            "name": "Hierarchical Vascular Representation",
            "quality_level": "whole-network coarse navigation graph",
            "coordinate_system": result.coordinate_system,
            "physical_unit": result.physical_unit,
            "origin_lps_um": result.origin_lps_um,
            "spacing_um": result.spacing_um,
            "radius_source": result.radius_source,
            "flow_direction_known": False,
            "parent_daughter_fields_available": False,
            "approved_for_cfd": False,
            "approved_as_final_geometry_training_truth": False,
        },
        "summary": result.report(),
        "scale_1_topology": {
            "nodes": [item.summary() for item in result.nodes],
            "branches": [
                {
                    "branch_id": item.branch_id,
                    "node_u": item.node_u,
                    "node_v": item.node_v,
                    "cycle_ids": item.cycle_ids,
                }
                for item in result.branches
            ],
            "cycles": [item.to_dict() for item in result.cycles],
            "branch_as_node_relations": relations,
        },
        "scale_2_branch_morphology": {
            "branches": [item.summary() for item in result.branches],
            "junction_local_geometry": [item.to_dict() for item in result.junctions],
        },
        "scale_3_dense_geometry": {
            "description": (
                "Raw voxel-center geometry is retained losslessly. Smoothed values are "
                "derived navigation measurements and never replace raw points."
            ),
            "branches": [item.dense_geometry() for item in result.branches],
            "node_regions": [
                {
                    "node_id": item.node_id,
                    "voxel_indices_xyz": item.voxel_indices_xyz,
                }
                for item in result.nodes
            ],
        },
    }


def _ragged_arrays(result: HierarchicalGraphResult) -> dict[str, np.ndarray]:
    raw_offsets = [0]
    smooth_offsets = [0]
    node_offsets = [0]
    for branch in result.branches:
        raw_offsets.append(raw_offsets[-1] + len(branch.points_raw_lps_um))
        smooth_offsets.append(smooth_offsets[-1] + len(branch.points_smoothed_lps_um))
    for node in result.nodes:
        node_offsets.append(node_offsets[-1] + len(node.voxel_indices_xyz))

    def concatenate(attribute: str, columns: int | None = None) -> np.ndarray:
        arrays = [np.asarray(getattr(item, attribute)) for item in result.branches]
        if arrays:
            return np.concatenate(arrays, axis=0)
        return np.empty((0, columns), dtype=float) if columns else np.empty(0, dtype=float)

    return {
        "branch_ids": np.asarray([item.branch_id for item in result.branches], dtype=np.int32),
        "branch_node_u": np.asarray([item.node_u for item in result.branches], dtype=np.int32),
        "branch_node_v": np.asarray([item.node_v for item in result.branches], dtype=np.int32),
        "raw_offsets": np.asarray(raw_offsets, dtype=np.int64),
        "raw_voxel_indices_xyz": concatenate("voxel_indices_xyz", 3).astype(np.int32),
        "raw_points_lps_um": concatenate("points_raw_lps_um", 3).astype(np.float64),
        "raw_arc_length_um": concatenate("arc_length_raw_um").astype(np.float64),
        "raw_coarse_radius_um": concatenate("coarse_radius_raw_um").astype(np.float64),
        "smoothed_offsets": np.asarray(smooth_offsets, dtype=np.int64),
        "smoothed_points_lps_um": concatenate("points_smoothed_lps_um", 3).astype(np.float64),
        "smoothed_arc_length_um": concatenate("arc_length_smoothed_um").astype(np.float64),
        "smoothed_coarse_radius_um": concatenate("coarse_radius_smoothed_um").astype(np.float64),
        "smoothed_local_direction_lps": concatenate("local_direction_smoothed", 3).astype(np.float64),
        "smoothed_curvature_per_um": concatenate("curvature_smoothed_per_um").astype(np.float64),
        "node_ids": np.asarray([item.node_id for item in result.nodes], dtype=np.int32),
        "node_offsets": np.asarray(node_offsets, dtype=np.int64),
        "node_region_voxel_indices_xyz": (
            np.concatenate([item.voxel_indices_xyz for item in result.nodes], axis=0).astype(np.int32)
            if result.nodes
            else np.empty((0, 3), dtype=np.int32)
        ),
    }


def _write_branch_vtp(result: HierarchicalGraphResult, path: Path) -> None:
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    point_radius = vtk.vtkDoubleArray()
    point_radius.SetName("coarse_radius_um")
    point_branch = vtk.vtkIntArray()
    point_branch.SetName("branch_id")
    point_order = vtk.vtkIntArray()
    point_order.SetName("point_order")
    cell_branch = vtk.vtkIntArray()
    cell_branch.SetName("branch_id")
    cell_node_u = vtk.vtkIntArray()
    cell_node_u.SetName("node_u")
    cell_node_v = vtk.vtkIntArray()
    cell_node_v.SetName("node_v")
    cell_length = vtk.vtkDoubleArray()
    cell_length.SetName("length_raw_um")

    point_id = 0
    for branch in result.branches:
        polyline = vtk.vtkPolyLine()
        polyline.GetPointIds().SetNumberOfIds(len(branch.points_raw_lps_um))
        for order, (coordinate, radius) in enumerate(
            zip(branch.points_raw_lps_um, branch.coarse_radius_raw_um)
        ):
            points.InsertNextPoint(*(float(value) for value in coordinate))
            polyline.GetPointIds().SetId(order, point_id)
            point_radius.InsertNextValue(float(radius))
            point_branch.InsertNextValue(branch.branch_id)
            point_order.InsertNextValue(order)
            point_id += 1
        lines.InsertNextCell(polyline)
        cell_branch.InsertNextValue(branch.branch_id)
        cell_node_u.InsertNextValue(branch.node_u)
        cell_node_v.InsertNextValue(branch.node_v)
        cell_length.InsertNextValue(float(branch.arc_length_raw_um[-1]))

    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    polydata.GetPointData().AddArray(point_radius)
    polydata.GetPointData().AddArray(point_branch)
    polydata.GetPointData().AddArray(point_order)
    polydata.GetCellData().AddArray(cell_branch)
    polydata.GetCellData().AddArray(cell_node_u)
    polydata.GetCellData().AddArray(cell_node_v)
    polydata.GetCellData().AddArray(cell_length)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(polydata)
    writer.SetDataModeToBinary()
    if writer.Write() != 1:
        raise OSError(f"Failed to write branch VTP: {path}")


def _write_node_vtp(result: HierarchicalGraphResult, path: Path) -> None:
    type_codes = {
        "isolated": 0,
        "terminal": 1,
        "connector": 2,
        "junction": 3,
        "complex_junction": 4,
        "cycle_anchor": 5,
    }
    points = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()
    node_ids = vtk.vtkIntArray()
    node_ids.SetName("node_id")
    degrees = vtk.vtkIntArray()
    degrees.SetName("graph_degree")
    node_types = vtk.vtkIntArray()
    node_types.SetName("node_type_code")
    voxel_counts = vtk.vtkIntArray()
    voxel_counts.SetName("node_region_voxel_count")
    for point_id, node in enumerate(result.nodes):
        points.InsertNextPoint(*node.representative_lps_um)
        vertex = vtk.vtkVertex()
        vertex.GetPointIds().SetId(0, point_id)
        vertices.InsertNextCell(vertex)
        node_ids.InsertNextValue(node.node_id)
        degrees.InsertNextValue(node.graph_degree)
        node_types.InsertNextValue(type_codes[node.node_type])
        voxel_counts.InsertNextValue(len(node.voxel_indices_xyz))
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetVerts(vertices)
    polydata.GetPointData().AddArray(node_ids)
    polydata.GetPointData().AddArray(degrees)
    polydata.GetPointData().AddArray(node_types)
    polydata.GetPointData().AddArray(voxel_counts)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(polydata)
    writer.SetDataModeToBinary()
    if writer.Write() != 1:
        raise OSError(f"Failed to write node VTP: {path}")


def export_hierarchical_graph(
    result: HierarchicalGraphResult,
    graphs_dir: Path,
    tables_dir: Path,
    *,
    save_graphml: bool,
    save_vtp: bool,
    save_npz: bool,
) -> list[Path]:
    graphs_dir.mkdir(parents=True, exist_ok=True)
    tables_dir.mkdir(parents=True, exist_ok=True)
    output: list[Path] = []

    json_path = graphs_dir / "hierarchical_vascular_graph.json"
    write_json(_graph_payload(result), json_path)
    output.append(json_path)

    if save_npz:
        npz_path = graphs_dir / "branch_geometry.npz"
        np.savez_compressed(npz_path, **_ragged_arrays(result))
        output.append(npz_path)

    if save_graphml:
        junction_graph_path = graphs_dir / "junction_branch_graph.graphml"
        branch_graph_path = graphs_dir / "branch_as_node_graph.graphml"
        nx.write_graphml(result.junction_graph, junction_graph_path)
        nx.write_graphml(result.branch_as_node_graph, branch_graph_path)
        output.extend((junction_graph_path, branch_graph_path))

    if save_vtp:
        branch_vtp_path = graphs_dir / "branch_centerlines.vtp"
        node_vtp_path = graphs_dir / "vascular_nodes.vtp"
        _write_branch_vtp(result, branch_vtp_path)
        _write_node_vtp(result, node_vtp_path)
        output.extend((branch_vtp_path, node_vtp_path))

    nodes_path = tables_dir / "nodes.csv"
    branches_path = tables_dir / "branches.csv"
    angles_path = tables_dir / "junction_angles.csv"
    cycles_path = tables_dir / "cycles.csv"
    write_csv(
        [
            {
                "node_id": node.node_id,
                "node_type": node.node_type,
                "graph_degree": node.graph_degree,
                "voxel_count": len(node.voxel_indices_xyz),
                "x_lps_um": node.representative_lps_um[0],
                "y_lps_um": node.representative_lps_um[1],
                "z_lps_um": node.representative_lps_um[2],
                "incident_branch_ids": json.dumps(node.incident_branch_ids),
                "cycle_ids": json.dumps(node.cycle_ids),
            }
            for node in result.nodes
        ],
        nodes_path,
    )
    write_csv(
        [
            {
                **{
                    key: json.dumps(value) if isinstance(value, list) else value
                    for key, value in branch.summary().items()
                },
                "geometry_level": "coarse_navigation",
            }
            for branch in result.branches
        ],
        branches_path,
    )
    write_csv(
        [
            {"node_id": junction.node_id, **angle}
            for junction in result.junctions
            for angle in junction.pairwise_angles_deg
        ],
        angles_path,
    )
    write_csv(
        [
            {
                "cycle_id": cycle.cycle_id,
                "cycle_type": cycle.cycle_type,
                "node_ids": json.dumps(cycle.node_ids),
                "branch_ids": json.dumps(cycle.branch_ids),
            }
            for cycle in result.cycles
        ],
        cycles_path,
    )
    output.extend((nodes_path, branches_path, angles_path, cycles_path))
    return output


def verify_graph_exports(
    paths: list[Path], expected_nodes: int, expected_branches: int
) -> list[str]:
    errors: list[str] = []
    by_name = {path.name: path for path in paths}
    try:
        payload = json.loads(
            by_name["hierarchical_vascular_graph.json"].read_text(encoding="utf-8")
        )
        if payload["summary"]["node_count"] != expected_nodes:
            errors.append("JSON node count differs from memory")
        if payload["summary"]["branch_count"] != expected_branches:
            errors.append("JSON branch count differs from memory")
    except Exception as exc:
        errors.append(f"JSON reopen failed: {exc}")

    if "branch_geometry.npz" in by_name:
        try:
            with np.load(by_name["branch_geometry.npz"]) as arrays:
                if len(arrays["branch_ids"]) != expected_branches:
                    errors.append("NPZ branch count differs from memory")
        except Exception as exc:
            errors.append(f"NPZ reopen failed: {exc}")

    for name, expected_node_count, expected_edge_count in (
        ("junction_branch_graph.graphml", expected_nodes, expected_branches),
        ("branch_as_node_graph.graphml", expected_branches, None),
    ):
        if name not in by_name:
            continue
        try:
            graph = nx.read_graphml(by_name[name])
            if graph.number_of_nodes() != expected_node_count:
                errors.append(f"{name} node count differs from memory")
            if expected_edge_count is not None and graph.number_of_edges() != expected_edge_count:
                errors.append(f"{name} edge count differs from memory")
        except Exception as exc:
            errors.append(f"{name} reopen failed: {exc}")

    for name, expected_points_or_cells, mode in (
        ("branch_centerlines.vtp", expected_branches, "cells"),
        ("vascular_nodes.vtp", expected_nodes, "points"),
    ):
        if name not in by_name:
            continue
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(str(by_name[name]))
        try:
            reader.Update()
            data = reader.GetOutput()
            value = data.GetNumberOfCells() if mode == "cells" else data.GetNumberOfPoints()
            if value != expected_points_or_cells:
                errors.append(f"{name} {mode} count differs from memory")
        except Exception as exc:
            errors.append(f"{name} reopen failed: {exc}")
    return errors
