"""Stable run layout and unit-explicit CFD geometry exports."""

from __future__ import annotations

import csv
import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pyvista as pv
import trimesh
import yaml

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .types import (
    BranchGeometry,
    CFDRunLayout,
    HybridBuildDetails,
    PatchResult,
    PortGeometry,
)


def write_json(path: Path, payload: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def write_csv(
    path: Path,
    rows: Iterable[dict[str, Any]],
    *,
    fieldnames: Iterable[str] | None = None,
) -> Path:
    payload = list(rows)
    names = list(fieldnames or (payload[0].keys() if payload else ()))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names, extrasaction="ignore")
        if names:
            writer.writeheader()
            writer.writerows(payload)
    return path


def create_run_layout(output_root: Path, run_id: str | None = None) -> CFDRunLayout:
    base = Path(output_root).resolve() / "model_generate"
    base.mkdir(parents=True, exist_ok=True)
    if run_id is not None:
        if not re.fullmatch(r"[A-Za-z0-9._-]+", run_id):
            raise ValueError("--run-id may contain only letters, digits, dot, underscore, and hyphen")
        run_root = base / run_id
        if run_root.exists():
            raise FileExistsError(f"Refusing to overwrite existing CFD lumen run: {run_root}")
    else:
        stem = datetime.now().strftime("run_%Y%m%d_%H%M%S")
        run_root = base / stem
        counter = 1
        while run_root.exists():
            run_root = base / f"{stem}_{counter:02d}"
            counter += 1
    folders = {name: run_root / name for name in ("logs", "config", "manifests", "figures", "report", "rois")}
    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=False)
    return CFDRunLayout(run_root, folders["logs"], folders["config"], folders["manifests"], folders["figures"], folders["report"], folders["rois"])


def create_roi_layout(run_layout: CFDRunLayout, roi_id: str) -> dict[str, Path]:
    root = run_layout.rois / roi_id
    folders = {
        name: root / name
        for name in ("source", "centerline", "geometry", "boundary", "qc", "figures", "mesh")
    }
    folders["ports"] = folders["geometry"] / "ports"
    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=False)
    folders["root"] = root
    return folders


def write_resolved_config(path: Path, config: CFDLumenConfig) -> Path:
    path.write_text(
        yaml.safe_dump(config.report(), sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )
    return path


def write_source_metadata(
    roi: ROIRecord,
    branches: list[BranchGeometry],
    directories: dict[str, Path],
    *,
    sampling_run: Path,
) -> list[Path]:
    metadata = {
        "roi_id": roi.roi_id,
        "source_model_id": roi.source_model_id,
        "source_mouse_id": roi.source_mouse_id,
        "sampling_run": str(Path(sampling_run).resolve()),
        "geometry_source": "saved_representative_SWC_ROI",
        "segmentation_mask_used": False,
        "coordinate_unit": "um",
        "radius_unit": "um",
        "anchor_id": roi.anchor_id,
        "anchor_position_um": list(roi.anchor_position_um),
        "bbox_min_um": list(roi.bbox_min_um),
        "bbox_max_um": list(roi.bbox_max_um),
        "node_count": roi.node_count,
        "edge_count": roi.edge_count,
        "branch_count": len(branches),
        "cut_port_count": len(roi.cut_ports),
        "global_node_ids": list(map(int, roi.global_node_ids)),
        "global_edge_ids": list(map(int, roi.global_edge_ids)),
    }
    metadata_path = write_json(directories["source"] / "roi_metadata.json", metadata)
    mapping_path = write_csv(
        directories["source"] / "branch_mapping.csv",
        [
            {
                "branch_id": branch.branch_id,
                "local_node_ids": ";".join(map(str, branch.local_node_ids)),
                "source_global_nodes": ";".join(map(str, branch.source_global_nodes)),
                "source_global_edges": ";".join(map(str, branch.source_global_edges)),
                "raw_point_count": len(branch.raw_points_um),
                "resampled_point_count": len(branch.points_um),
                "length_um": branch.length_um,
            }
            for branch in branches
        ],
    )
    return [metadata_path, mapping_path]


def _centerline_polydata(branches: list[BranchGeometry], *, raw: bool) -> pv.PolyData:
    point_blocks = [branch.raw_points_um if raw else branch.points_um for branch in branches]
    radius_blocks = [branch.raw_radius_um if raw else branch.radius_um for branch in branches]
    points = np.concatenate(point_blocks, axis=0)
    radii = np.concatenate(radius_blocks, axis=0)
    point_branch = np.concatenate(
        [np.full(len(block), branch.branch_id, dtype=np.int32) for branch, block in zip(branches, point_blocks)]
    )
    lines: list[int] = []
    cell_branch: list[int] = []
    offset = 0
    for branch, block in zip(branches, point_blocks):
        lines.extend((len(block), *range(offset, offset + len(block))))
        cell_branch.append(branch.branch_id)
        offset += len(block)
    polydata = pv.PolyData(points, lines=np.asarray(lines, dtype=np.int64))
    polydata.point_data["radius_um"] = radii
    polydata.point_data["branch_id"] = point_branch
    polydata.cell_data["branch_id"] = np.asarray(cell_branch, dtype=np.int32)
    return polydata


def write_centerlines(branches: list[BranchGeometry], directories: dict[str, Path]) -> list[Path]:
    raw_path = directories["centerline"] / "raw_centerline.vtp"
    resampled_path = directories["centerline"] / "resampled_centerline.vtp"
    _centerline_polydata(branches, raw=True).save(raw_path)
    _centerline_polydata(branches, raw=False).save(resampled_path)
    return [raw_path, resampled_path]


def write_constructed_centerlines(
    details: HybridBuildDetails,
    directories: dict[str, Path],
) -> list[Path]:
    """Persist v5 CFD-derived extension points separately from source geometry."""

    branches = details.constructed_branches
    if not branches:
        return []
    polydata = _centerline_polydata(branches, raw=False)
    point_type = np.concatenate(
        [np.asarray(branch.construction_point_type, dtype=np.uint8) for branch in branches]
    )
    source_edge = np.concatenate(
        [
            np.asarray(branch.construction_source_global_edge_id, dtype=np.int64)
            for branch in branches
        ]
    )
    distance = np.concatenate(
        [
            np.asarray(
                branch.construction_distance_from_core_boundary_um, dtype=float
            )
            for branch in branches
        ]
    )
    original_core = np.concatenate(
        [
            np.asarray(
                branch.construction_original_core_cut_position_um, dtype=float
            )
            for branch in branches
        ],
        axis=0,
    )
    source_cut = np.concatenate(
        [np.asarray(branch.construction_source_cut_port_id) for branch in branches]
    )
    polydata.point_data["point_type"] = point_type
    polydata.point_data["source_global_edge_id"] = source_edge
    polydata.point_data["distance_from_core_boundary_um"] = distance
    polydata.point_data["original_core_cut_position_um"] = original_core
    # VTK XML supports string arrays, but a CSV mirror below remains the
    # authoritative human-readable representation across viewers.
    polydata.point_data["source_cut_port_id"] = source_cut
    vtp_path = directories["centerline"] / "constructed_centerline.vtp"
    polydata.save(vtp_path)

    csv_path = write_csv(
        directories["centerline"] / "constructed_centerline_points.csv",
        details.port_extension_rows,
    )
    metadata_path = write_json(
        directories["centerline"] / "constructed_centerline_metadata.json",
        {
            "point_type_codes": {
                "SOURCE_CENTERLINE": 0,
                "CFD_EXTENSION_ORIGIN": 1,
                "CFD_EXTENSION": 2,
            },
            "source_geometry_modified": False,
            "core_roi_boundary_type": "CORE_ROI_BOUNDARY",
            "final_boundary_type": "CFD_BOUNDARY_PORT",
            "construction": "centerline extension before one vtkTubeFilter",
            "separate_extension_cylinder": False,
            "derived_extension_point_count": len(details.port_extension_rows),
        },
    )
    return [vtp_path, csv_path, metadata_path]


def _mesh_polydata(
    mesh: trimesh.Trimesh,
    patch: PatchResult,
    scale: float,
    *,
    face_region: np.ndarray | None = None,
) -> pv.PolyData:
    faces = np.column_stack((np.full(len(mesh.faces), 3, dtype=np.int64), mesh.faces)).ravel()
    polydata = pv.PolyData(np.asarray(mesh.vertices, dtype=float) * scale, faces)
    polydata.cell_data["patch_id"] = patch.patch_id
    polydata.cell_data["patch_type"] = patch.patch_type
    polydata.cell_data["port_id"] = patch.port_id
    if face_region is not None and len(face_region) == len(mesh.faces):
        polydata.cell_data["surface_region"] = np.asarray(face_region, dtype=np.uint8)
    return polydata


def _scaled_submesh(mesh: trimesh.Trimesh, face_mask: np.ndarray, scale: float) -> trimesh.Trimesh:
    indices = np.flatnonzero(face_mask)
    if len(indices) == 0:
        return trimesh.Trimesh()
    result = mesh.submesh([indices], append=True, repair=False)
    result.apply_scale(scale)
    return result


def write_geometry_exports(
    mesh: trimesh.Trimesh,
    patch: PatchResult,
    ports: list[PortGeometry],
    directories: dict[str, Path],
    *,
    face_region: np.ndarray | None = None,
) -> list[Path]:
    geometry = directories["geometry"]
    paths: list[Path] = []
    surface_um = geometry / "lumen_surface_um.vtp"
    surface_m_vtp = geometry / "lumen_surface_m.vtp"
    surface_m_stl = geometry / "lumen_surface_m.stl"
    cfd_polydata = _mesh_polydata(mesh, patch, 1.0, face_region=face_region)
    cfd_polydata.save(surface_um)
    _mesh_polydata(mesh, patch, 1.0e-6, face_region=face_region).save(surface_m_vtp)
    scaled = mesh.copy()
    scaled.apply_scale(1.0e-6)
    scaled.export(surface_m_stl)
    paths.extend((surface_um, surface_m_vtp, surface_m_stl))
    cfd_path = geometry / "lumen_surface_cfd.vtp"
    cfd_polydata.save(cfd_path)
    visualization_path = geometry / "lumen_surface_visualization.vtp"
    visual = cfd_polydata.compute_normals(
        cell_normals=True,
        point_normals=True,
        split_vertices=False,
        consistent_normals=True,
        auto_orient_normals=True,
        inplace=False,
    )
    visual.save(visualization_path)
    paths.extend((cfd_path, visualization_path))
    wall_path = geometry / "wall_m.stl"
    wall = _scaled_submesh(mesh, patch.patch_type == 0, 1.0e-6)
    wall.export(wall_path)
    paths.append(wall_path)
    for port in ports:
        port_path = directories["ports"] / f"port_{port.port_id:03d}_m.stl"
        port_mesh = _scaled_submesh(mesh, patch.port_id == port.port_id, 1.0e-6)
        port_mesh.export(port_path)
        paths.append(port_path)
    return paths


def write_units(path: Path) -> Path:
    return write_json(
        path,
        {
            "source_coordinate_unit": "um",
            "source_radius_unit": "um",
            "visualization_and_qc_unit": "um",
            "cfd_export_unit": "m",
            "scale_um_to_m": 1.0e-6,
            "lumen_surface_cfd.vtp": (
                "authoritative micrometre-coordinate geometry and face connectivity; "
                "no visualization-only point normals"
            ),
            "lumen_surface_visualization.vtp": (
                "same vertex coordinates with separately computed smooth point normals"
            ),
            "stl_unit_note": "STL is unitless; all *_m.stl vertex coordinates are metres.",
        },
    )


def verify_volume_mesh(
    surface_stl_m: Path,
    output_path: Path,
    *,
    minimum_radius_um: float,
    config: CFDLumenConfig,
) -> dict[str, Any]:
    if not config.volume_mesh.enabled:
        return {"enabled": False, "status": "NOT_RUN"}
    try:
        import gmsh
    except ImportError as exc:
        return {
            "enabled": True,
            "status": "FAIL",
            "failure_reason": f"gmsh is not installed: {exc}",
        }
    initialized_here = False
    try:
        if not gmsh.isInitialized():
            gmsh.initialize()
            initialized_here = True
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.clear()
        gmsh.model.add("cfd_lumen")
        gmsh.merge(str(Path(surface_stl_m).resolve()))
        gmsh.model.mesh.classifySurfaces(np.deg2rad(40.0), True, True, np.pi)
        gmsh.model.mesh.createGeometry()
        gmsh.model.geo.synchronize()
        surfaces = [tag for dimension, tag in gmsh.model.getEntities(2) if dimension == 2]
        loop = gmsh.model.geo.addSurfaceLoop(surfaces)
        gmsh.model.geo.addVolume([loop])
        gmsh.model.geo.synchronize()
        target = max(minimum_radius_um * 1.0e-6 * config.volume_mesh.characteristic_length_factor, 1.0e-9)
        gmsh.option.setNumber("Mesh.MeshSizeMin", target * 0.5)
        gmsh.option.setNumber("Mesh.MeshSizeMax", target * 2.0)
        gmsh.model.mesh.generate(3)
        node_tags, _, _ = gmsh.model.mesh.getNodes()
        element_types, element_tags, _ = gmsh.model.mesh.getElements(3)
        tags = np.concatenate([np.asarray(values, dtype=np.int64) for values in element_tags]) if element_tags else np.empty(0, dtype=np.int64)
        qualities = np.asarray(gmsh.model.mesh.getElementQualities(tags.tolist()), dtype=float) if len(tags) else np.empty(0)
        volumes = (
            np.asarray(gmsh.model.mesh.getElementQualities(tags.tolist(), "volume"), dtype=float)
            if len(tags)
            else np.empty(0)
        )
        gmsh.write(str(output_path.resolve()))
        return {
            "enabled": True,
            "status": "PASS" if len(tags) and len(node_tags) else "FAIL",
            "number_of_nodes": int(len(node_tags)),
            "number_of_cells": int(len(tags)),
            "element_types": list(map(int, element_types)),
            "minimum_cell_quality": float(qualities.min()) if len(qualities) else None,
            "negative_volume_cell_count": int(np.count_nonzero(volumes <= 0.0)),
            "boundary_patch_recovery_status": "NOT_VERIFIED_BY_OPTIONAL_TETRAHEDRALIZATION",
            "mesh_path": str(output_path),
        }
    except Exception as exc:
        return {"enabled": True, "status": "FAIL", "failure_reason": f"{type(exc).__name__}: {exc}"}
    finally:
        if initialized_here:
            gmsh.finalize()
