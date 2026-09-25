from __future__ import annotations

import vtk

from utils.config import MeshCleanupConfig
from utils.mesh.cleanup import cleanup_mesh


def _sphere(radius: float, center: tuple[float, float, float], resolution: int) -> vtk.vtkPolyData:
    source = vtk.vtkSphereSource()
    source.SetRadius(radius)
    source.SetCenter(center)
    source.SetThetaResolution(resolution)
    source.SetPhiResolution(resolution)
    source.Update()
    output = vtk.vtkPolyData()
    output.DeepCopy(source.GetOutput())
    return output


def test_cleanup_removes_only_component_below_all_thresholds() -> None:
    large = _sphere(10.0, (0.0, 0.0, 0.0), 24)
    tiny = _sphere(0.5, (25.0, 0.0, 0.0), 8)
    append = vtk.vtkAppendPolyData()
    append.AddInputData(large)
    append.AddInputData(tiny)
    append.Update()

    result = cleanup_mesh(
        append.GetOutput(),
        MeshCleanupConfig(
            min_component_faces=200,
            min_component_area_um2=50.0,
            min_component_diagonal_um=5.0,
        ),
    )

    assert sum(item.decision == "keep" for item in result.components) == 1
    assert sum(item.decision == "remove" for item in result.components) == 1
    assert sum(item.component_type == "small_fragment" for item in result.components) == 1
    assert result.removed_mesh is not None
    assert result.cleaned_quality.connected_component_count == 1
    assert result.cleaned_quality.boundary_edge_count == 0
    assert result.cleaned_quality.non_manifold_edge_count == 0


def test_component_is_kept_when_any_size_measure_is_large() -> None:
    large = _sphere(8.0, (0.0, 0.0, 0.0), 20)
    thin_but_long = _sphere(1.0, (20.0, 0.0, 0.0), 12)
    append = vtk.vtkAppendPolyData()
    append.AddInputData(large)
    append.AddInputData(thin_but_long)
    append.Update()

    result = cleanup_mesh(
        append.GetOutput(),
        MeshCleanupConfig(
            component_policy="conservative",
            min_component_faces=10_000,
            min_component_area_um2=10_000.0,
            min_component_diagonal_um=1.5,
        ),
    )
    assert len(result.components) == 2
    assert all(item.decision == "keep" for item in result.components)


def test_main_network_only_removes_a_large_disconnected_island() -> None:
    main = _sphere(10.0, (0.0, 0.0, 0.0), 28)
    branched_size_island = _sphere(4.0, (30.0, 0.0, 0.0), 20)
    append = vtk.vtkAppendPolyData()
    append.AddInputData(main)
    append.AddInputData(branched_size_island)
    append.Update()

    result = cleanup_mesh(append.GetOutput(), MeshCleanupConfig())

    assert result.cleaned_quality.connected_component_count == 1
    assert sum(item.component_type == "main_network" for item in result.components) == 1
    assert sum(item.component_type == "island_network" for item in result.components) == 1
    assert result.removed_island_networks_mesh is not None
    assert result.summary()["island_network_count"] == 1


def test_main_component_can_be_selected_explicitly() -> None:
    large = _sphere(9.0, (0.0, 0.0, 0.0), 24)
    smaller = _sphere(4.0, (25.0, 0.0, 0.0), 18)
    append = vtk.vtkAppendPolyData()
    append.AddInputData(large)
    append.AddInputData(smaller)
    append.Update()
    automatic = cleanup_mesh(append.GetOutput(), MeshCleanupConfig())
    smaller_id = next(
        item.component_id
        for item in automatic.components
        if item.component_type == "island_network"
    )

    explicit = cleanup_mesh(
        append.GetOutput(), MeshCleanupConfig(main_component_id=smaller_id)
    )
    assert explicit.main_network_selection.selected_component_id == smaller_id
    assert explicit.main_network_selection.selection_method == "explicit component ID"
    assert explicit.cleaned_quality.connected_component_count == 1
