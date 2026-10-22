"""Conservative connected-component cleanup for vascular STL surfaces."""

from __future__ import annotations

import logging
from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np
import vtk
from pymeshfix import clean_from_arrays
from vtk.util.numpy_support import vtk_to_numpy

from ..config import MeshCleanupConfig
from .components import MainNetworkSelection, is_small_fragment, select_main_network
from .quality import (
    MeshQuality,
    build_polydata,
    compact_polydata,
    duplicate_triangle_count,
    edge_topology,
    face_areas,
    measure_mesh_quality,
    polydata_arrays,
    remove_degenerate_triangles,
)


@dataclass(slots=True)
class ComponentRecord:
    component_id: int
    area_rank: int
    decision: str
    component_type: str
    decision_reason: str
    triangle_count: int
    surface_area_um2: float
    enclosed_volume_um3: float
    bbox_x_um: float
    bbox_y_um: float
    bbox_z_um: float
    bbox_diagonal_um: float
    boundary_edge_count: int
    non_manifold_edge_count: int
    duplicate_triangle_count: int
    repair_status: str = "not_needed"

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class MeshCleanupResult:
    cleaned_mesh: vtk.vtkPolyData
    removed_mesh: vtk.vtkPolyData | None
    removed_small_fragments_mesh: vtk.vtkPolyData | None
    removed_island_networks_mesh: vtk.vtkPolyData | None
    input_quality: MeshQuality
    cleaned_quality: MeshQuality
    components: list[ComponentRecord]
    main_network_selection: MainNetworkSelection
    degenerate_triangles_removed: int
    repair_failures: list[int] = field(default_factory=list)

    def summary(self) -> dict[str, Any]:
        removed = [item for item in self.components if item.decision == "remove"]
        kept = [item for item in self.components if item.decision == "keep"]
        small_fragments = [
            item for item in self.components if item.component_type == "small_fragment"
        ]
        island_networks = [
            item for item in self.components if item.component_type == "island_network"
        ]
        main_record = next(
            item for item in self.components if item.component_type == "main_network"
        )
        input_area = self.input_quality.surface_area_um2
        removed_area = sum(item.surface_area_um2 for item in removed)
        return {
            "input_quality": self.input_quality.to_dict(),
            "cleaned_quality": self.cleaned_quality.to_dict(),
            "component_count": len(self.components),
            "kept_component_count": len(kept),
            "removed_component_count": len(removed),
            "small_fragment_count": len(small_fragments),
            "island_network_count": len(island_networks),
            "main_network_component_id": main_record.component_id,
            "main_network_surface_area_um2": main_record.surface_area_um2,
            "main_network_surface_area_fraction": (
                main_record.surface_area_um2 / input_area if input_area else 0.0
            ),
            "small_fragment_surface_area_um2": sum(
                item.surface_area_um2 for item in small_fragments
            ),
            "island_network_surface_area_um2": sum(
                item.surface_area_um2 for item in island_networks
            ),
            "removed_surface_area_um2": removed_area,
            "removed_surface_area_fraction": removed_area / input_area if input_area else 0.0,
            "main_network_selection": self.main_network_selection.to_dict(),
            "degenerate_triangles_removed": self.degenerate_triangles_removed,
            "repair_failures": self.repair_failures,
        }


@dataclass(slots=True)
class _ComponentAnalysis:
    points: np.ndarray
    faces: np.ndarray
    labels: np.ndarray
    records: list[ComponentRecord]
    selection: MainNetworkSelection


def _triangulate_and_merge(mesh: vtk.vtkPolyData) -> vtk.vtkPolyData:
    triangles = vtk.vtkTriangleFilter()
    triangles.SetInputData(mesh)
    triangles.PassLinesOff()
    triangles.PassVertsOff()

    cleaner = vtk.vtkCleanPolyData()
    cleaner.SetInputConnection(triangles.GetOutputPort())
    cleaner.PointMergingOn()
    cleaner.ToleranceIsAbsoluteOn()
    cleaner.SetAbsoluteTolerance(0.0)
    cleaner.ConvertLinesToPointsOff()
    cleaner.ConvertPolysToLinesOff()
    cleaner.ConvertStripsToPolysOff()
    cleaner.Update()

    output = vtk.vtkPolyData()
    output.DeepCopy(cleaner.GetOutput())
    return output


def _analyze_components(mesh: vtk.vtkPolyData, config: MeshCleanupConfig) -> _ComponentAnalysis:
    connectivity = vtk.vtkConnectivityFilter()
    connectivity.SetInputData(mesh)
    connectivity.SetExtractionModeToAllRegions()
    connectivity.ColorRegionsOn()
    connectivity.Update()
    connected = connectivity.GetOutput()
    points, faces = polydata_arrays(connected)
    region_array = connected.GetCellData().GetArray("RegionId")
    if region_array is None:
        raise RuntimeError("VTK did not produce connected-component labels")
    labels = vtk_to_numpy(region_array).astype(np.int64, copy=False)
    component_count = int(connectivity.GetNumberOfExtractedRegions())

    areas = face_areas(points, faces)
    a = points[faces[:, 0]]
    b = points[faces[:, 1]]
    c = points[faces[:, 2]]
    signed_volumes = np.einsum("ij,ij->i", a, np.cross(b, c)) / 6.0
    counts = np.bincount(labels, minlength=component_count)
    area_sums = np.bincount(labels, weights=areas, minlength=component_count)
    volume_sums = np.abs(np.bincount(labels, weights=signed_volumes, minlength=component_count))

    minimums = np.full((component_count, 3), np.inf, dtype=np.float64)
    maximums = np.full((component_count, 3), -np.inf, dtype=np.float64)
    for axis in range(3):
        cell_min = np.minimum(np.minimum(a[:, axis], b[:, axis]), c[:, axis])
        cell_max = np.maximum(np.maximum(a[:, axis], b[:, axis]), c[:, axis])
        np.minimum.at(minimums[:, axis], labels, cell_min)
        np.maximum.at(maximums[:, axis], labels, cell_max)
    extents = maximums - minimums
    diagonals = np.linalg.norm(extents, axis=1)
    ranks = np.empty(component_count, dtype=np.int64)
    ranks[np.argsort(area_sums)[::-1]] = np.arange(1, component_count + 1)
    selection = select_main_network(
        policy=config.component_policy,
        requested_component_id=config.main_component_id,
        triangle_counts=counts,
        surface_areas=area_sums,
        diagonals=diagonals,
    )
    main_id = selection.selected_component_id

    order = np.argsort(labels, kind="stable")
    offsets = np.concatenate(([0], np.cumsum(np.bincount(labels[order], minlength=component_count))))
    records: list[ComponentRecord] = []
    for component_id in range(component_count):
        component_faces = faces[order[offsets[component_id] : offsets[component_id + 1]]]
        boundary_edges, non_manifold_edges = edge_topology(component_faces)
        duplicates = duplicate_triangle_count(component_faces)
        small_fragment = is_small_fragment(
            triangle_count=int(counts[component_id]),
            surface_area_um2=float(area_sums[component_id]),
            diagonal_um=float(diagonals[component_id]),
            min_faces=config.min_component_faces,
            min_area_um2=config.min_component_area_um2,
            min_diagonal_um=config.min_component_diagonal_um,
        )
        if component_id == main_id:
            decision = "keep"
            component_type = "main_network"
            reason = f"selected main network by {selection.selection_method}"
        elif small_fragment:
            decision = "remove"
            component_type = "small_fragment"
            reason = "small fragment below face, area, and diagonal thresholds"
        elif config.component_policy == "main_network_only":
            decision = "remove"
            component_type = "island_network"
            reason = "disconnected from the selected main network"
        else:
            decision = "keep"
            component_type = "retained_secondary_network"
            reason = "exceeds at least one conservative size threshold"
        records.append(
            ComponentRecord(
                component_id=component_id,
                area_rank=int(ranks[component_id]),
                decision=decision,
                component_type=component_type,
                decision_reason=reason,
                triangle_count=int(counts[component_id]),
                surface_area_um2=float(area_sums[component_id]),
                enclosed_volume_um3=float(volume_sums[component_id]),
                bbox_x_um=float(extents[component_id, 0]),
                bbox_y_um=float(extents[component_id, 1]),
                bbox_z_um=float(extents[component_id, 2]),
                bbox_diagonal_um=float(diagonals[component_id]),
                boundary_edge_count=boundary_edges,
                non_manifold_edge_count=non_manifold_edges,
                duplicate_triangle_count=duplicates,
            )
        )
    return _ComponentAnalysis(
        points=points,
        faces=faces,
        labels=labels,
        records=records,
        selection=selection,
    )


def _repair_component(
    points: np.ndarray,
    faces: np.ndarray,
    config: MeshCleanupConfig,
) -> tuple[vtk.vtkPolyData | None, str]:
    component = compact_polydata(points, faces)
    local_points, local_faces = polydata_arrays(component)
    original_face_count = len(local_faces)
    try:
        repaired_points, repaired_faces = clean_from_arrays(
            np.ascontiguousarray(local_points, dtype=np.float64),
            np.ascontiguousarray(local_faces, dtype=np.int32),
            verbose=False,
            joincomp=False,
            remove_smallest_components=False,
        )
    except Exception as exc:  # pragma: no cover - backend-specific failure
        return None, f"repair backend error: {exc}"
    repaired_faces = np.asarray(repaired_faces, dtype=np.int64).reshape((-1, 3))
    if len(repaired_faces) == 0:
        return None, "repair produced an empty component"
    change = abs(len(repaired_faces) - original_face_count) / max(original_face_count, 1)
    boundary, non_manifold = edge_topology(repaired_faces)
    duplicates = duplicate_triangle_count(repaired_faces)
    if change > config.max_repair_face_change_fraction:
        return None, f"repair changed {change:.1%} of faces"
    if boundary or non_manifold or duplicates:
        return None, (
            f"repair left boundary={boundary}, non_manifold={non_manifold}, "
            f"duplicates={duplicates}"
        )
    return build_polydata(repaired_points, repaired_faces), "repaired"


def _append_meshes(meshes: list[vtk.vtkPolyData]) -> vtk.vtkPolyData:
    append = vtk.vtkAppendPolyData()
    for mesh in meshes:
        if mesh.GetNumberOfCells():
            append.AddInputData(mesh)
    append.Update()
    output = vtk.vtkPolyData()
    output.DeepCopy(append.GetOutput())
    return output


def _orient_normals(mesh: vtk.vtkPolyData) -> vtk.vtkPolyData:
    normals = vtk.vtkPolyDataNormals()
    normals.SetInputData(mesh)
    normals.ConsistencyOn()
    normals.AutoOrientNormalsOn()
    normals.SplittingOff()
    normals.ComputePointNormalsOff()
    normals.ComputeCellNormalsOn()
    normals.NonManifoldTraversalOn()
    normals.Update()
    output = vtk.vtkPolyData()
    output.DeepCopy(normals.GetOutput())
    return output


def _smooth(mesh: vtk.vtkPolyData, config: MeshCleanupConfig) -> vtk.vtkPolyData:
    if config.smoothing_iterations == 0:
        return mesh
    smoother = vtk.vtkWindowedSincPolyDataFilter()
    smoother.SetInputData(mesh)
    smoother.SetNumberOfIterations(config.smoothing_iterations)
    smoother.SetPassBand(config.smoothing_pass_band)
    smoother.NormalizeCoordinatesOn()
    smoother.BoundarySmoothingOff()
    smoother.FeatureEdgeSmoothingOff()
    smoother.NonManifoldSmoothingOff()
    smoother.Update()
    output = vtk.vtkPolyData()
    output.DeepCopy(smoother.GetOutput())
    return output


def cleanup_mesh(
    input_mesh: vtk.vtkPolyData,
    config: MeshCleanupConfig,
    logger: logging.Logger | None = None,
) -> MeshCleanupResult:
    logger = logger or logging.getLogger("ulm_3d_vascular")
    logger.info("Measuring input mesh quality")
    input_quality = measure_mesh_quality(input_mesh, config.degenerate_area_epsilon_um2)

    prepared = _triangulate_and_merge(input_mesh)
    prepared, degenerate_removed = remove_degenerate_triangles(
        prepared, config.degenerate_area_epsilon_um2
    )
    logger.info("Removed %d degenerate triangles", degenerate_removed)
    analysis = _analyze_components(prepared, config)
    kept_ids = {item.component_id for item in analysis.records if item.decision == "keep"}
    removed_ids = {item.component_id for item in analysis.records if item.decision == "remove"}
    logger.info(
        "Connected components: %d total, %d kept, %d removed",
        len(analysis.records),
        len(kept_ids),
        len(removed_ids),
    )
    logger.info(
        "Main network: component %d selected by %s; rank agreement=%s",
        analysis.selection.selected_component_id,
        analysis.selection.selection_method,
        analysis.selection.ranking_agrees,
    )

    issue_ids = {
        item.component_id
        for item in analysis.records
        if item.decision == "keep"
        and (item.boundary_edge_count or item.non_manifold_edge_count or item.duplicate_triangle_count)
    }
    labels = analysis.labels
    normal_keep_mask = np.isin(labels, list(kept_ids - issue_ids))
    meshes_to_keep: list[vtk.vtkPolyData] = []
    if np.any(normal_keep_mask):
        meshes_to_keep.append(compact_polydata(analysis.points, analysis.faces[normal_keep_mask]))

    repair_failures: list[int] = []
    records_by_id = {item.component_id: item for item in analysis.records}
    for component_id in sorted(issue_ids):
        record = records_by_id[component_id]
        faces = analysis.faces[labels == component_id]
        if config.repair_non_manifold:
            logger.info("Repairing topology of component %d", component_id)
            repaired, status = _repair_component(analysis.points, faces, config)
        else:
            repaired, status = None, "repair disabled"
        if repaired is None:
            record.repair_status = f"failed: {status}"
            repair_failures.append(component_id)
            meshes_to_keep.append(compact_polydata(analysis.points, faces))
            logger.warning("Component %d repair failed: %s", component_id, status)
        else:
            record.repair_status = status
            meshes_to_keep.append(repaired)

    cleaned = _append_meshes(meshes_to_keep)
    cleaned = _orient_normals(cleaned)
    cleaned = _smooth(cleaned, config)
    cleaned, final_degenerate_removed = remove_degenerate_triangles(
        cleaned, config.degenerate_area_epsilon_um2
    )
    degenerate_removed += final_degenerate_removed

    removed_mesh: vtk.vtkPolyData | None = None
    removed_small_fragments_mesh: vtk.vtkPolyData | None = None
    removed_island_networks_mesh: vtk.vtkPolyData | None = None
    if removed_ids:
        remove_mask = np.isin(labels, list(removed_ids))
        if np.any(remove_mask):
            removed_mesh = compact_polydata(analysis.points, analysis.faces[remove_mask])
    small_fragment_ids = {
        item.component_id
        for item in analysis.records
        if item.component_type == "small_fragment"
    }
    if small_fragment_ids:
        small_mask = np.isin(labels, list(small_fragment_ids))
        if np.any(small_mask):
            removed_small_fragments_mesh = compact_polydata(
                analysis.points, analysis.faces[small_mask]
            )
    island_network_ids = {
        item.component_id
        for item in analysis.records
        if item.component_type == "island_network"
    }
    if island_network_ids:
        island_mask = np.isin(labels, list(island_network_ids))
        if np.any(island_mask):
            removed_island_networks_mesh = compact_polydata(
                analysis.points, analysis.faces[island_mask]
            )

    logger.info("Measuring cleaned mesh quality")
    cleaned_quality = measure_mesh_quality(cleaned, config.degenerate_area_epsilon_um2)
    return MeshCleanupResult(
        cleaned_mesh=cleaned,
        removed_mesh=removed_mesh,
        removed_small_fragments_mesh=removed_small_fragments_mesh,
        removed_island_networks_mesh=removed_island_networks_mesh,
        input_quality=input_quality,
        cleaned_quality=cleaned_quality,
        components=analysis.records,
        main_network_selection=analysis.selection,
        degenerate_triangles_removed=degenerate_removed,
        repair_failures=repair_failures,
    )
