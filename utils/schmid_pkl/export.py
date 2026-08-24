"""Portable tables, graph formats, VTK geometry, and preview volume exports."""

from __future__ import annotations

import json
from html import escape
from pathlib import Path
from typing import Any

import nibabel as nib
import networkx as nx
import numpy as np
import vtk
from skimage.draw import line_nd

from ..io import write_csv, write_json
from ..reporting.acceptance import AcceptanceResult
from .config import SchmidPKLConfig
from .layout import SchmidOutputLayout
from .model import DirectedGraphResult, VESSEL_TYPE_NAMES


def _offsets(arrays: list[np.ndarray]) -> np.ndarray:
    values = [0]
    for array in arrays:
        values.append(values[-1] + len(array))
    return np.asarray(values, dtype=np.int64)


def _concatenate(arrays: list[np.ndarray], shape: tuple[int, ...], dtype: Any) -> np.ndarray:
    return np.concatenate(arrays, axis=0) if arrays else np.empty(shape, dtype=dtype)


def export_normalized_source(result: DirectedGraphResult, output_dir: Path) -> Path:
    cleanup = result.cleanup
    source = cleanup.source
    edges = cleanup.edges
    points = [item.points_u_to_v_um for item in edges]
    diameters = [item.diameter_u_to_v_um for item in edges]
    path = output_dir / "cleaned_schmid_source.npz"
    np.savez_compressed(
        path,
        valid_node_ids=cleanup.valid_node_ids,
        coordinates_source_xyz_um=source.coordinates_um[cleanup.valid_node_ids],
        pressure_mmhg=source.pressure_mmhg[cleanup.valid_node_ids],
        pressure_boundary_mmhg=source.pressure_boundary_mmhg[cleanup.valid_node_ids],
        edge_ids=np.asarray([item.edge_id for item in edges], dtype=np.int64),
        edge_node_u=np.asarray([item.node_u for item in edges], dtype=np.int64),
        edge_node_v=np.asarray([item.node_v for item in edges], dtype=np.int64),
        mean_diameter_um=np.asarray([item.mean_diameter_um for item in edges]),
        source_length_um=np.asarray([item.source_length_um for item in edges]),
        flow_um3_per_ms=np.asarray([item.flow_um3_per_ms for item in edges]),
        vessel_type_code=np.asarray([item.vessel_type_code for item in edges], dtype=np.int16),
        point_offsets=_offsets(points),
        points_source_xyz_um=_concatenate(points, (0, 3), np.float64),
        diameter_offsets=_offsets(diameters),
        diameter_profiles_um=_concatenate(diameters, (0,), np.float64),
    )
    return path


def export_branch_geometry(result: DirectedGraphResult, output_dir: Path) -> Path:
    branches = result.branches
    raw_points = [item.points_raw_um for item in branches]
    raw_radius = [item.radius_raw_um for item in branches]
    smooth_points = [item.points_smoothed_um for item in branches]
    smooth_radius = [item.radius_smoothed_um for item in branches]
    smooth_arc = [item.arc_length_smoothed_um for item in branches]
    smooth_direction = [item.local_direction_smoothed for item in branches]
    smooth_curvature = [item.curvature_smoothed_per_um for item in branches]
    raw_edges = [np.asarray(item.raw_edge_ids, dtype=np.int64) for item in branches]
    path = output_dir / "branch_geometry.npz"
    np.savez_compressed(
        path,
        branch_ids=np.asarray([item.branch_id for item in branches], dtype=np.int64),
        branch_node_a=np.asarray([item.node_a for item in branches], dtype=np.int64),
        branch_node_b=np.asarray([item.node_b for item in branches], dtype=np.int64),
        branch_upstream_node=np.asarray(
            [item.upstream_node if item.upstream_node is not None else -1 for item in branches],
            dtype=np.int64,
        ),
        branch_downstream_node=np.asarray(
            [item.downstream_node if item.downstream_node is not None else -1 for item in branches],
            dtype=np.int64,
        ),
        branch_direction_known=np.asarray(
            [item.direction_status == "known" for item in branches], dtype=np.uint8
        ),
        raw_edge_offsets=_offsets(raw_edges),
        raw_edge_ids=_concatenate(raw_edges, (0,), np.int64),
        raw_offsets=_offsets(raw_points),
        raw_points_source_xyz_um=_concatenate(raw_points, (0, 3), np.float64),
        raw_radius_um=_concatenate(raw_radius, (0,), np.float64),
        smoothed_offsets=_offsets(smooth_points),
        smoothed_points_source_xyz_um=_concatenate(smooth_points, (0, 3), np.float64),
        smoothed_radius_um=_concatenate(smooth_radius, (0,), np.float64),
        smoothed_arc_length_um=_concatenate(smooth_arc, (0,), np.float64),
        smoothed_local_direction=_concatenate(smooth_direction, (0, 3), np.float64),
        smoothed_curvature_per_um=_concatenate(smooth_curvature, (0,), np.float64),
        node_ids=np.asarray([item.node_id for item in result.nodes], dtype=np.int64),
        node_coordinates_source_xyz_um=np.asarray([item.coordinates_um for item in result.nodes]),
        node_pressure_mmhg=np.asarray([item.pressure_mmhg for item in result.nodes]),
    )
    return path


def _hierarchical_payload(result: DirectedGraphResult) -> dict[str, Any]:
    topology_branches = []
    morphology_branches = []
    for item in result.branches:
        summary = item.summary()
        topology_branches.append(
            {
                key: value
                for key, value in summary.items()
                if key
                in {
                    "branch_id",
                    "raw_edge_ids",
                    "node_a",
                    "node_b",
                    "upstream_node",
                    "downstream_node",
                    "direction_status",
                    "parent_branch_ids",
                    "child_branch_ids",
                    "flow_um3_per_ms",
                    "pressure_drop_mmhg",
                    "vessel_type_codes",
                }
            }
        )
        morphology_branches.append(
            {
                key: value
                for key, value in summary.items()
                if key
                not in {
                    "raw_edge_ids",
                    "parent_branch_ids",
                    "child_branch_ids",
                    "cycle_ids",
                }
            }
        )
    return {
        "schema_version": "1.0",
        "representation": {
            "name": "Directed Schmid Hierarchical Vascular Graph",
            "coordinate_system": "source_XYZ_not_anatomically_labeled",
            "physical_unit": "micrometer",
            "directed": True,
            "direction_method": "higher endpoint pressure to lower endpoint pressure",
            "flow_values_are_magnitudes": True,
            "multiple_parents_and_children_allowed": True,
            "original_surface_available": False,
            "approved_for_cfd": False,
        },
        "summary": result.report(),
        "scale_1_directed_topology": {
            "nodes": [
                {key: value for key, value in item.summary().items() if key != "cycle_ids"}
                for item in result.nodes
            ],
            "branches": topology_branches,
            "cycle_basis": {
                "count": len(result.cycles),
                "full_records": "../reports/cycle_report.json",
                "table": "../tables/cycles.csv",
            },
        },
        "scale_2_branch_morphology": {
            "branches": morphology_branches
        },
        "scale_3_dense_geometry": {
            "description": (
                "Dense arrays are stored externally to keep this JSON practical for NW1. "
                "Coordinates are source XYZ micrometers, not asserted LPS."
            ),
            "branch_count": len(result.branches),
            "npz": "branch_geometry.npz",
            "vtp": "branch_centerlines.vtp",
        },
    }


def _write_branch_vtp(result: DirectedGraphResult, path: Path) -> None:
    points = vtk.vtkPoints()
    lines = vtk.vtkCellArray()
    radius = vtk.vtkDoubleArray()
    radius.SetName("radius_um")
    branch_ids = vtk.vtkIntArray()
    branch_ids.SetName("branch_id")
    direction_known = vtk.vtkIntArray()
    direction_known.SetName("direction_known")
    flows = vtk.vtkDoubleArray()
    flows.SetName("flow_um3_per_ms")
    pressure_drops = vtk.vtkDoubleArray()
    pressure_drops.SetName("pressure_drop_mmhg")
    vessel_types = vtk.vtkIntArray()
    vessel_types.SetName("primary_vessel_type_code")
    geometry_codes = vtk.vtkIntArray()
    geometry_codes.SetName("geometry_issue")
    for branch in result.branches:
        sequence = branch.points_smoothed_um
        radii = branch.radius_smoothed_um
        line = vtk.vtkPolyLine()
        line.GetPointIds().SetNumberOfIds(len(sequence))
        for local_index, coordinate in enumerate(sequence):
            point_id = points.InsertNextPoint(*(float(value) for value in coordinate))
            line.GetPointIds().SetId(local_index, point_id)
            radius.InsertNextValue(float(radii[local_index]))
        lines.InsertNextCell(line)
        branch_ids.InsertNextValue(branch.branch_id)
        direction_known.InsertNextValue(int(branch.direction_status == "known"))
        flows.InsertNextValue(branch.flow_um3_per_ms)
        pressure_drops.InsertNextValue(branch.pressure_drop_mmhg or 0.0)
        vessel_types.InsertNextValue(branch.vessel_type_codes[0])
        geometry_codes.InsertNextValue(int(branch.geometry_status != "valid"))
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    polydata.GetPointData().AddArray(radius)
    for array in (branch_ids, direction_known, flows, pressure_drops, vessel_types, geometry_codes):
        polydata.GetCellData().AddArray(array)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(polydata)
    writer.SetDataModeToBinary()
    if writer.Write() != 1:
        raise OSError(f"Failed to write VTP: {path}")


def _write_node_vtp(result: DirectedGraphResult, path: Path) -> None:
    points = vtk.vtkPoints()
    vertices = vtk.vtkCellArray()
    node_ids = vtk.vtkIntArray()
    node_ids.SetName("node_id")
    pressure = vtk.vtkDoubleArray()
    pressure.SetName("pressure_mmhg")
    indegree = vtk.vtkIntArray()
    indegree.SetName("directed_indegree")
    outdegree = vtk.vtkIntArray()
    outdegree.SetName("directed_outdegree")
    role_codes = {name: index for index, name in enumerate(sorted({n.node_role for n in result.nodes}))}
    roles = vtk.vtkIntArray()
    roles.SetName("node_role_code")
    for point_id, node in enumerate(result.nodes):
        points.InsertNextPoint(*node.coordinates_um)
        vertex = vtk.vtkVertex()
        vertex.GetPointIds().SetId(0, point_id)
        vertices.InsertNextCell(vertex)
        node_ids.InsertNextValue(node.node_id)
        pressure.InsertNextValue(node.pressure_mmhg)
        indegree.InsertNextValue(len(node.incoming_branch_ids))
        outdegree.InsertNextValue(len(node.outgoing_branch_ids))
        roles.InsertNextValue(role_codes[node.node_role])
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetVerts(vertices)
    for array in (node_ids, pressure, indegree, outdegree, roles):
        polydata.GetPointData().AddArray(array)
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(path))
    writer.SetInputData(polydata)
    writer.SetDataModeToBinary()
    if writer.Write() != 1:
        raise OSError(f"Failed to write VTP: {path}")
    write_json({"node_role_codes": role_codes}, path.with_name("vascular_node_role_codes.json"))


def export_graph_artifacts(result: DirectedGraphResult, layout: SchmidOutputLayout) -> list[Path]:
    outputs: list[Path] = []
    normalized = export_normalized_source(result, layout.normalized)
    outputs.append(normalized)
    graph_json = layout.graphs / "directed_hierarchical_vascular_graph.json"
    write_json(_hierarchical_payload(result), graph_json)
    outputs.append(graph_json)
    geometry = export_branch_geometry(result, layout.graphs)
    outputs.append(geometry)
    graphml_items = (
        (result.directed_junction_graph, layout.graphs / "directed_junction_branch_graph.graphml"),
        (result.directed_branch_graph, layout.graphs / "directed_branch_as_node_graph.graphml"),
        (result.all_connectivity_graph, layout.graphs / "all_connectivity_graph.graphml"),
    )
    for graph, path in graphml_items:
        nx.write_graphml(graph, path)
        outputs.append(path)
    branch_vtp = layout.graphs / "branch_centerlines.vtp"
    node_vtp = layout.graphs / "vascular_nodes.vtp"
    _write_branch_vtp(result, branch_vtp)
    _write_node_vtp(result, node_vtp)
    outputs.extend((branch_vtp, node_vtp, node_vtp.with_name("vascular_node_role_codes.json")))

    nodes_path = layout.tables / "nodes.csv"
    branches_path = layout.tables / "branches.csv"
    raw_edges_path = layout.tables / "raw_edges.csv"
    relations_path = layout.tables / "parent_child_relations.csv"
    cycles_path = layout.tables / "cycles.csv"
    unresolved_path = layout.tables / "unresolved_directions.csv"
    cleanup_path = layout.tables / "cleanup_decisions.csv"
    flow_path = layout.tables / "flow_conservation.csv"
    write_csv(
        [
            {
                **{
                    key: json.dumps(value) if isinstance(value, (list, tuple)) else value
                    for key, value in node.summary().items()
                }
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
                }
            }
            for branch in result.branches
        ],
        branches_path,
    )
    write_csv([edge.topology_summary() for edge in result.cleanup.edges], raw_edges_path)
    write_csv(
        [
            {
                "parent_branch_id": parent,
                "child_branch_id": child,
                "shared_node_id": data["shared_node_id"],
                "relationship": "parent_to_child",
            }
            for parent, child, data in result.directed_branch_graph.edges(data=True)
        ],
        relations_path,
    )
    write_csv(
        [
            {
                "cycle_id": item.cycle_id,
                "cycle_type": item.cycle_type,
                "node_ids": json.dumps(item.node_ids),
                "branch_ids": json.dumps(item.branch_ids),
            }
            for item in result.cycles
        ],
        cycles_path,
    )
    write_csv(
        [
            {
                "level": "raw_edge",
                "item_id": edge.edge_id,
                "direction_status": edge.direction_status,
                "node_u": edge.node_u,
                "node_v": edge.node_v,
                "pressure_u_mmhg": result.cleanup.source.pressure_mmhg[edge.node_u],
                "pressure_v_mmhg": result.cleanup.source.pressure_mmhg[edge.node_v],
                "flow_um3_per_ms": edge.flow_um3_per_ms,
            }
            for edge in result.cleanup.edges
            if edge.direction_status != "known"
        ]
        + [
            {
                "level": "branch",
                "item_id": branch.branch_id,
                "direction_status": branch.direction_status,
                "node_u": branch.node_a,
                "node_v": branch.node_b,
                "pressure_u_mmhg": result.cleanup.source.pressure_mmhg[branch.node_a],
                "pressure_v_mmhg": result.cleanup.source.pressure_mmhg[branch.node_b],
                "flow_um3_per_ms": branch.flow_um3_per_ms,
            }
            for branch in result.branches
            if branch.direction_status != "known"
        ],
        unresolved_path,
    )
    write_csv([item.to_dict() for item in result.cleanup.decisions], cleanup_path)
    write_csv(result.flow_conservation, flow_path)
    outputs.extend(
        (
            nodes_path,
            branches_path,
            raw_edges_path,
            relations_path,
            cycles_path,
            unresolved_path,
            cleanup_path,
            flow_path,
        )
    )
    return outputs


def write_preview_volume(
    result: DirectedGraphResult, config: SchmidPKLConfig, layout: SchmidOutputLayout
) -> tuple[Path | None, dict[str, Any]]:
    all_points = np.concatenate([item.points_raw_um for item in result.branches], axis=0)
    spacing = float(config.preview_voxel_size_um)
    origin = all_points.min(axis=0) - spacing * 2
    maximum = all_points.max(axis=0) + spacing * 2
    shape = tuple((np.ceil((maximum - origin) / spacing).astype(int) + 1).tolist())
    voxel_count = int(np.prod(shape, dtype=np.int64))
    report = {
        "requested": True,
        "voxel_size_um": spacing,
        "origin_source_xyz_um": origin,
        "shape_xyz": shape,
        "voxel_count": voxel_count,
        "max_voxel_count": config.max_preview_voxel_count,
        "used_for_graph_construction": False,
        "represents_original_lumen_surface": False,
    }
    if voxel_count > config.max_preview_voxel_count:
        report["written"] = False
        report["reason"] = "estimated grid exceeds max_preview_voxel_count"
        return None, report
    volume = np.zeros(shape, dtype=np.uint8)
    for branch in result.branches:
        indices = np.rint((branch.points_raw_um - origin) / spacing).astype(int)
        volume[tuple(indices.T)] = 1
        differences = np.abs(np.diff(indices, axis=0))
        jump_ids = np.flatnonzero(np.max(differences, axis=1) > 1)
        for jump_id in jump_ids:
            coordinates = line_nd(indices[jump_id], indices[jump_id + 1], endpoint=True)
            volume[coordinates] = 1
    affine = np.eye(4, dtype=float)
    affine[0, 0] = affine[1, 1] = affine[2, 2] = spacing
    affine[:3, 3] = origin
    image = nib.Nifti1Image(volume, affine)
    image.header.set_xyzt_units("micron")
    image.header["descrip"] = "Derived centerline preview; source XYZ axes are not anatomical LPS/RAS"[:79]
    path = layout.volumes / "derived_centerline_preview.nii.gz"
    nib.save(image, str(path))
    report["written"] = True
    report["foreground_voxel_count"] = int(np.count_nonzero(volume))
    return path, report


def write_acceptance_html(
    layout: SchmidOutputLayout,
    result: DirectedGraphResult,
    acceptance: AcceptanceResult,
    images: list[Path],
    files: list[Path],
) -> None:
    checks = "".join(
        f"<tr><td>{escape(item.status)}</td><td>{escape(item.name)}</td>"
        f"<td>{escape(item.message)}</td></tr>"
        for item in acceptance.checks
    )
    image_html = "".join(
        f'<section><h3>{escape(path.name)}</h3><img src="{escape(path.relative_to(layout.run_root).as_posix())}" '
        'style="max-width:100%;border:1px solid #ccc"></section>'
        for path in images
        if path.is_file()
    )
    file_html = "".join(
        f'<li><a href="{escape(path.relative_to(layout.run_root).as_posix())}">{escape(path.relative_to(layout.run_root).as_posix())}</a></li>'
        for path in files
        if path.is_file() and path != layout.html_report
    )
    summary = result.report()
    html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>Directed Schmid graph acceptance</title>
<style>body{{font-family:Arial,sans-serif;max-width:1400px;margin:2rem auto;padding:0 1rem;color:#222}}
table{{border-collapse:collapse;width:100%}}td,th{{border:1px solid #ccc;padding:.45rem;text-align:left}}
.status{{font-size:1.5rem;font-weight:bold}}code{{background:#eee;padding:.1rem .25rem}}</style></head><body>
<h1>Directed Schmid vascular graph acceptance</h1>
<p class="status">Overall status: {escape(acceptance.overall_status)}</p>
<p>Nodes: {summary['node_count']:,}; branches: {summary['branch_count']:,}; known directions: {summary['known_direction_branch_count']:,}; unresolved: {summary['unresolved_direction_branch_count']:,}.</p>
<p>Direction means higher simulated pressure to lower simulated pressure. This graph allows multiple parents and children and is not forced into a tree.</p>
<h2>Checks</h2><table><thead><tr><th>Status</th><th>Check</th><th>Meaning</th></tr></thead><tbody>{checks}</tbody></table>
<h2>Visual review</h2>{image_html}<h2>Files</h2><ul>{file_html}</ul></body></html>"""
    layout.html_report.write_text(html, encoding="utf-8")
