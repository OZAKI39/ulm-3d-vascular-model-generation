"""Thin configured-environment runner for official VMTK filters.

All project-specific mapping and QC remain in the main environment.  This file
only reads VTP, invokes VMTK, and writes VTP/STL exchange files.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import sys
from pathlib import Path

import vtk
from vmtk import vtkvmtk, vmtkscripts


def _read_vtp(path: str) -> vtk.vtkPolyData:
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(path)
    reader.Update()
    output = vtk.vtkPolyData()
    output.DeepCopy(reader.GetOutput())
    if output.GetNumberOfCells() == 0:
        raise RuntimeError(f"VTP has no cells: {path}")
    return output


def _write_vtp(surface: vtk.vtkPolyData, path: str) -> None:
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(path)
    writer.SetInputData(surface)
    writer.SetDataModeToBinary()
    if writer.Write() != 1:
        raise RuntimeError(f"Failed to write VTP: {path}")


def _write_stl(surface: vtk.vtkPolyData, path: str) -> None:
    triangle = vtk.vtkTriangleFilter()
    triangle.SetInputData(surface)
    triangle.PassLinesOff()
    triangle.PassVertsOff()
    triangle.Update()
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(path)
    writer.SetInputData(triangle.GetOutput())
    writer.SetFileTypeToBinary()
    if writer.Write() != 1:
        raise RuntimeError(f"Failed to write STL: {path}")


def _run(request: dict[str, object]) -> dict[str, object]:
    surface = _read_vtp(str(request["input_surface_vtp"]))
    centerlines = _read_vtp(str(request["centerlines_vtp"]))
    parameters = dict(request["parameters"])  # type: ignore[arg-type]

    flow = vtkvmtk.vtkvmtkPolyDataFlowExtensionsFilter()
    flow.SetInputData(surface)
    flow.SetCenterlines(centerlines)
    flow.SetSigma(float(parameters["sigma"]))
    flow.SetAdaptiveExtensionLength(int(bool(parameters["adaptive_extension_length"])))
    flow.SetExtensionRatio(float(parameters["extension_ratio"]))
    flow.SetAdaptiveExtensionRadius(int(bool(parameters["adaptive_extension_radius"])))
    flow.SetAdaptiveNumberOfBoundaryPoints(
        int(bool(parameters["adaptive_boundary_points"]))
    )
    flow.SetTransitionRatio(float(parameters["transition_ratio"]))
    preserve_shape_api = hasattr(flow, "SetPreserveCrossSectionShape")
    if preserve_shape_api:
        flow.SetPreserveCrossSectionShape(
            int(bool(parameters["preserve_cross_section_shape"]))
        )
    elif bool(parameters["preserve_cross_section_shape"]):
        raise RuntimeError(
            "Installed official VMTK release cannot preserve the source cross-section"
        )
    flow.SetExtensionModeToUseCenterlineDirection()
    flow.SetInterpolationModeToThinPlateSpline()
    flow.Update()
    raw = vtk.vtkPolyData()
    raw.DeepCopy(flow.GetOutput())
    if raw.GetNumberOfCells() <= surface.GetNumberOfCells():
        raise RuntimeError("VMTK flow-extension filter did not add surface cells")
    _write_vtp(raw, str(request["raw_vtp"]))
    _write_stl(raw, str(request["raw_stl"]))

    remesher = vmtkscripts.vmtkSurfaceRemeshing()
    remesher.Surface = raw
    remesher.ElementSizeMode = "edgelength"
    remesher.TargetEdgeLength = float(request["target_edge_length_um"])
    remesher.PreserveBoundaryEdges = 1
    remesher.Execute()
    remeshed_open = vtk.vtkPolyData()
    remeshed_open.DeepCopy(remesher.Surface)
    _write_vtp(remeshed_open, str(request["remeshed_open_vtp"]))

    capper = vmtkscripts.vmtkSurfaceCapper()
    capper.Surface = remeshed_open
    capper.Method = "simple"
    capper.Interactive = 0
    capper.TriangleOutput = 1
    capper.CellEntityIdsArrayName = "CellEntityIds"
    capper.CellEntityIdOffset = 1
    capper.Execute()
    capped = vtk.vtkPolyData()
    capped.DeepCopy(capper.Surface)
    _write_vtp(capped, str(request["capped_vtp"]))

    entity_ids = capped.GetCellData().GetArray("CellEntityIds")
    unique_entities: list[int] = []
    wall_entity_id: int | None = None
    cap_entity_ids: list[int] = []
    if entity_ids is not None:
        values = [int(entity_ids.GetTuple1(index)) for index in range(capped.GetNumberOfCells())]
        unique_entities = sorted(set(values))
        wall_entity_id = max(unique_entities, key=values.count)
        cap_entity_ids = [value for value in unique_entities if value != wall_entity_id]
    return {
        "status": "PASS",
        "python": sys.version,
        "python_executable": sys.executable,
        "vmtk_package_version": importlib.metadata.version("vmtk"),
        "vtk_version": vtk.vtkVersion.GetVTKVersion(),
        "official_flow_filter": "vtkvmtkPolyDataFlowExtensionsFilter",
        "official_remesher": "vmtkscripts.vmtkSurfaceRemeshing",
        "official_capper": "vmtkscripts.vmtkSurfaceCapper",
        "input_points": surface.GetNumberOfPoints(),
        "input_cells": surface.GetNumberOfCells(),
        "raw_points": raw.GetNumberOfPoints(),
        "raw_cells": raw.GetNumberOfCells(),
        "remeshed_open_points": remeshed_open.GetNumberOfPoints(),
        "remeshed_open_cells": remeshed_open.GetNumberOfCells(),
        "capped_points": capped.GetNumberOfPoints(),
        "capped_cells": capped.GetNumberOfCells(),
        "cell_entity_ids": unique_entities,
        "wall_entity_id": wall_entity_id,
        "cap_entity_ids": cap_entity_ids,
        "parameters": parameters,
        "remesh": {
            "element_size_mode": "edgelength",
            "target_edge_length_um": float(request["target_edge_length_um"]),
            "preserve_boundary_edges": True,
            "parameter_sweep_count": 0,
            "target_size_count": 1,
        },
        "cap": {
            "method": "simple",
            "cell_entity_ids_array_name": "CellEntityIds",
            "cell_entity_id_offset": 1,
        },
        "custom_tps_implementation": False,
        "preserve_cross_section_shape_api_available": preserve_shape_api,
        "preserve_cross_section_shape_effective": False,
        "preserve_cross_section_shape_compatibility": (
            "explicit official filter setter"
            if preserve_shape_api
            else "official v1.5.0 behavior always transitions the source profile to the target circle"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--request", type=Path, required=True)
    parser.add_argument("--result", type=Path, required=True)
    args = parser.parse_args()
    request = json.loads(args.request.read_text(encoding="utf-8"))
    result = _run(request)
    args.result.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
