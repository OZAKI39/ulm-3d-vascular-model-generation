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


def _runtime() -> dict[str, object]:
    return {
        "status": "PASS",
        "python": sys.version,
        "python_executable": sys.executable,
        "vmtk_package_version": importlib.metadata.version("vmtk"),
        "vtk_version": vtk.vtkVersion.GetVTKVersion(),
        "custom_tps_implementation": False,
    }


def _run_extension(request: dict[str, object]) -> dict[str, object]:
    surface = _read_vtp(str(request["input_surface_vtp"]))
    parameters = dict(request["parameters"])  # type: ignore[arg-type]
    extension_mode = str(parameters["extension_mode"])

    flow = vtkvmtk.vtkvmtkPolyDataFlowExtensionsFilter()
    flow.SetInputData(surface)
    if extension_mode == "centerlinedirection":
        if "centerlines_vtp" not in request:
            raise RuntimeError("INVALID_VMTK_EXTENSION_MODE:centerlines_required")
        centerlines = _read_vtp(str(request["centerlines_vtp"]))
        flow.SetCenterlines(centerlines)
        flow.SetExtensionModeToUseCenterlineDirection()
        centerlines_used = True
        direction_api = "SetExtensionModeToUseCenterlineDirection"
    elif extension_mode == "boundarynormal":
        flow.SetExtensionModeToUseNormalToBoundary()
        centerlines_used = False
        direction_api = "SetExtensionModeToUseNormalToBoundary"
    else:
        raise RuntimeError("INVALID_VMTK_EXTENSION_MODE")
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
    flow.SetInterpolationModeToThinPlateSpline()
    flow.Update()
    raw = vtk.vtkPolyData()
    raw.DeepCopy(flow.GetOutput())
    if raw.GetNumberOfCells() <= surface.GetNumberOfCells():
        raise RuntimeError("VMTK flow-extension filter did not add surface cells")
    _write_vtp(raw, str(request["raw_vtp"]))
    _write_stl(raw, str(request["raw_stl"]))
    return {
        **_runtime(),
        "operation": "extension",
        "official_flow_filter": "vtkvmtkPolyDataFlowExtensionsFilter",
        "extension_mode_effective": extension_mode,
        "centerlines_used_for_extension_direction": centerlines_used,
        "official_direction_api": direction_api,
        "official_interpolation_api": "SetInterpolationModeToThinPlateSpline",
        "input_points": surface.GetNumberOfPoints(),
        "input_cells": surface.GetNumberOfCells(),
        "raw_points": raw.GetNumberOfPoints(),
        "raw_cells": raw.GetNumberOfCells(),
        "parameters": parameters,
        "preserve_cross_section_shape_api_available": preserve_shape_api,
        "preserve_cross_section_shape_effective": False,
        "preserve_cross_section_shape_compatibility": (
            "explicit official filter setter"
            if preserve_shape_api
            else "official v1.5.0 behavior always transitions the source profile to the target circle"
        ),
    }


def _run_remesh_cap(request: dict[str, object]) -> dict[str, object]:
    """LEGACY_GLOBAL_REMESH_REFERENCE_ONLY."""

    raw = _read_vtp(str(request["raw_vtp"]))

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
        **_runtime(),
        "operation": "remesh_cap",
        "legacy_capability": "LEGACY_GLOBAL_REMESH_REFERENCE_ONLY",
        "surface_remesher_called": True,
        "global_surface_remeshing_performed": True,
        "official_remesher": "vmtkscripts.vmtkSurfaceRemeshing",
        "official_capper": "vmtkscripts.vmtkSurfaceCapper",
        "raw_points": raw.GetNumberOfPoints(),
        "raw_cells": raw.GetNumberOfCells(),
        "remeshed_open_points": remeshed_open.GetNumberOfPoints(),
        "remeshed_open_cells": remeshed_open.GetNumberOfCells(),
        "capped_points": capped.GetNumberOfPoints(),
        "capped_cells": capped.GetNumberOfCells(),
        "cell_entity_ids": unique_entities,
        "wall_entity_id": wall_entity_id,
        "cap_entity_ids": cap_entity_ids,
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
    }


def _run_cap_only(request: dict[str, object]) -> dict[str, object]:
    """Cap the official RAW flow-extension output without surface remeshing."""

    raw = _read_vtp(str(request["raw_vtp"]))
    capper = vmtkscripts.vmtkSurfaceCapper()
    capper.Surface = raw
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
        values = [
            int(entity_ids.GetTuple1(index))
            for index in range(capped.GetNumberOfCells())
        ]
        unique_entities = sorted(set(values))
        wall_entity_id = max(unique_entities, key=values.count)
        cap_entity_ids = [value for value in unique_entities if value != wall_entity_id]
    return {
        **_runtime(),
        "operation": "cap_only",
        "official_capper": "vmtkscripts.vmtkSurfaceCapper",
        "surface_remesher_called": False,
        "global_surface_remeshing_performed": False,
        "raw_points": raw.GetNumberOfPoints(),
        "raw_cells": raw.GetNumberOfCells(),
        "capped_points": capped.GetNumberOfPoints(),
        "capped_cells": capped.GetNumberOfCells(),
        "cell_entity_ids": unique_entities,
        "wall_entity_id": wall_entity_id,
        "cap_entity_ids": cap_entity_ids,
        "cap": {
            "method": "simple",
            "interactive": False,
            "triangle_output": True,
            "cell_entity_ids_array_name": "CellEntityIds",
            "cell_entity_id_offset": 1,
        },
    }


def _run_entity_remesh(request: dict[str, object]) -> dict[str, object]:
    """Remesh only explicitly active entities with the official VMTK filter."""

    raw = _read_vtp(str(request["raw_vtp"]))
    entity_array_name = str(request["entity_array_name"])
    entity_array = raw.GetCellData().GetArray(entity_array_name)
    if entity_array is None:
        raise RuntimeError("VMTK_ENTITY_ASSIGNMENT_FAILED:array_missing")
    input_entity_ids = sorted(
        {
            int(entity_array.GetTuple1(index))
            for index in range(raw.GetNumberOfCells())
        }
    )
    expected_entity_ids = sorted(
        int(value) for value in request["expected_entity_ids"]  # type: ignore[union-attr]
    )
    if input_entity_ids != expected_entity_ids:
        raise RuntimeError("VMTK_ENTITY_ASSIGNMENT_FAILED:unexpected_ids")
    excluded_entity_ids = sorted(
        int(value) for value in request["excluded_entity_ids"]  # type: ignore[union-attr]
    )
    active_entity_ids = sorted(
        int(value) for value in request["active_entity_ids"]  # type: ignore[union-attr]
    )
    if (
        not excluded_entity_ids
        or not active_entity_ids
        or set(excluded_entity_ids) & set(active_entity_ids)
        or sorted(set(excluded_entity_ids) | set(active_entity_ids))
        != expected_entity_ids
    ):
        raise RuntimeError("VMTK_ENTITY_ASSIGNMENT_FAILED:unsafe_exclusion")

    remesher = vmtkscripts.vmtkSurfaceRemeshing()
    remesher.Surface = raw
    remesher.CellEntityIdsArrayName = entity_array_name
    remesher.ExcludeEntityIds = excluded_entity_ids
    remesher.ElementSizeMode = str(request["element_size_mode"])
    remesher.TargetEdgeLength = float(request["target_edge_length_um"])
    remesher.PreserveBoundaryEdges = int(bool(request["preserve_boundary_edges"]))
    remesher.Execute()
    remeshed = vtk.vtkPolyData()
    remeshed.DeepCopy(remesher.Surface)
    output_entities = remeshed.GetCellData().GetArray(entity_array_name)
    if output_entities is None:
        raise RuntimeError("VMTK_ENTITY_ASSIGNMENT_FAILED:output_array_missing")
    output_entity_ids = sorted(
        {
            int(output_entities.GetTuple1(index))
            for index in range(remeshed.GetNumberOfCells())
        }
    )
    if output_entity_ids != input_entity_ids:
        raise RuntimeError("VMTK_ENTITY_ASSIGNMENT_FAILED:output_ids_changed")
    _write_vtp(remeshed, str(request["remeshed_open_vtp"]))
    _write_stl(remeshed, str(request["remeshed_open_stl"]))
    return {
        **_runtime(),
        "operation": "entity_remesh",
        "official_remesher": "vmtkscripts.vmtkSurfaceRemeshing",
        "cell_entity_ids_array": entity_array_name,
        "excluded_entity_ids": excluded_entity_ids,
        "active_entity_ids": active_entity_ids,
        "input_entity_ids": input_entity_ids,
        "output_entity_ids": output_entity_ids,
        "target_edge_length_um": float(request["target_edge_length_um"]),
        "element_size_mode": str(request["element_size_mode"]),
        "preserve_boundary_edges": bool(request["preserve_boundary_edges"]),
        "raw_points": raw.GetNumberOfPoints(),
        "raw_cells": raw.GetNumberOfCells(),
        "remeshed_open_points": remeshed.GetNumberOfPoints(),
        "remeshed_open_cells": remeshed.GetNumberOfCells(),
        "surface_remesher_called": True,
        "global_surface_remeshing_performed": False,
        "entity_aware_extension_remeshing_performed": True,
        "official_default_parameters_retained": {
            "number_of_iterations": remesher.NumberOfIterations,
            "connectivity_optimization_iterations": (
                remesher.NumberOfConnectivityOptimizationIterations
            ),
            "relaxation": remesher.Relaxation,
            "aspect_ratio_threshold": remesher.AspectRatioThreshold,
        },
    }


def _run(request: dict[str, object]) -> dict[str, object]:
    operation = str(request.get("operation", ""))
    if operation == "extension":
        return _run_extension(request)
    if operation == "remesh_cap":
        return _run_remesh_cap(request)
    if operation == "cap_only":
        return _run_cap_only(request)
    if operation == "entity_remesh":
        return _run_entity_remesh(request)
    raise RuntimeError("INVALID_VMTK_OPERATION")


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
