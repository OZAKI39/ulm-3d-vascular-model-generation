"""Geometry-only validation and figures for raw Ultraliser surfaces."""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import TYPE_CHECKING, Any

import numpy as np
import pyvista as pv
import trimesh

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .export import write_csv, write_json
from .geometry_preprocess import validate_and_extract_branches
from .mesh_defects import _triangle_intersections
from .surface_qc import evaluate_radius_fidelity

if TYPE_CHECKING:
    from .ultraliser_backend import UltraliserLayout


def _load_mesh(path: Path) -> trimesh.Trimesh:
    loaded = trimesh.load_mesh(path, process=True, validate=False)
    if isinstance(loaded, trimesh.Scene):
        loaded = loaded.to_mesh()
    if not isinstance(loaded, trimesh.Trimesh) or len(loaded.faces) == 0:
        raise ValueError(f"No triangular surface was loaded from {path}")
    return loaded


def _polydata(mesh: trimesh.Trimesh) -> pv.PolyData:
    faces = np.column_stack(
        (np.full(len(mesh.faces), 3, dtype=np.int64), np.asarray(mesh.faces, dtype=np.int64))
    ).ravel()
    return pv.PolyData(np.asarray(mesh.vertices, dtype=float), faces)


def _write_geometry_outputs(
    raw_surface: Path,
    final_directory: Path,
) -> tuple[trimesh.Trimesh, dict[str, Any]]:
    surface_um_stl = final_directory / "ultraliser_surface_um.stl"
    surface_um_vtp = final_directory / "ultraliser_surface_um.vtp"
    surface_m_stl = final_directory / "ultraliser_surface_m.stl"
    shutil.copy2(raw_surface, surface_um_stl)
    mesh = _load_mesh(surface_um_stl)
    _polydata(mesh).save(surface_um_vtp, binary=True)
    mesh_m = trimesh.Trimesh(
        vertices=np.asarray(mesh.vertices, dtype=float) * 1.0e-6,
        faces=np.asarray(mesh.faces, dtype=np.int64),
        process=False,
    )
    mesh_m.export(surface_m_stl)
    units = {
        "canonical_source_morphology_coordinate_unit": "um",
        "canonical_source_morphology_radius_unit": "um",
        "raw_ultraliser_geometry_unit": "um",
        "ultraliser_surface_um.vtp": "um; topology-only conversion, no smoothing/remeshing/vertex motion",
        "ultraliser_surface_um.stl": "um; byte-for-byte copy of official watertight STL",
        "ultraliser_surface_m.stl": "m; copied coordinates scaled by 1e-6 for future CFD",
        "scale_um_to_m": 1.0e-6,
    }
    write_json(final_directory / "units.json", units)
    return mesh, units


def _surface_qc(mesh: trimesh.Trimesh) -> dict[str, Any]:
    sorted_edges = np.sort(np.asarray(mesh.edges, dtype=np.int64), axis=1)
    _, edge_counts = np.unique(sorted_edges, axis=0, return_counts=True)
    boundary_edges = int(np.count_nonzero(edge_counts == 1))
    nonmanifold_edges = int(np.count_nonzero(edge_counts > 2))
    areas = np.asarray(mesh.area_faces, dtype=float)
    diagonal = float(np.linalg.norm(np.ptp(np.asarray(mesh.vertices), axis=0)))
    area_tolerance = max(np.finfo(float).eps * diagonal**2 * 100.0, 1.0e-18)
    repeated_index = np.asarray(
        [len(set(map(int, face))) < 3 for face in np.asarray(mesh.faces)], dtype=bool
    )
    degenerate = int(np.count_nonzero((areas <= area_tolerance) | repeated_index))
    all_faces = np.arange(len(mesh.faces), dtype=np.int64)
    intersection_pairs, candidate_pairs = _triangle_intersections(mesh, all_faces)
    components = mesh.split(only_watertight=False)
    report = {
        "status": "PASS"
        if (
            len(components) == 1
            and mesh.is_watertight
            and boundary_edges == 0
            and nonmanifold_edges == 0
            and len(intersection_pairs) == 0
            and degenerate == 0
        )
        else "FAIL",
        "vertex_count": int(len(mesh.vertices)),
        "triangle_count": int(len(mesh.faces)),
        "component_count": int(len(components)),
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "boundary_edge_count": boundary_edges,
        "nonmanifold_edge_count": nonmanifold_edges,
        "self_intersection_count": int(len(intersection_pairs)),
        "self_intersection_candidate_pairs_checked": int(candidate_pairs),
        "self_intersection_scope": "entire mesh via exact R-tree AABB candidates and vtkTriangle",
        "degenerate_triangle_count": degenerate,
        "degenerate_area_tolerance_um2": area_tolerance,
        "surface_area_um2": float(mesh.area),
        "signed_volume_um3": float(mesh.volume),
        "volume_um3": float(abs(mesh.volume)),
        "bounds_um": np.asarray(mesh.bounds, dtype=float).tolist(),
        "raw_ultraliser_geometry_checked_before_port_clipping": True,
    }
    return report


def _source_branches(roi: ROIRecord, config: CFDLumenConfig):
    branches, _ = validate_and_extract_branches(roi, config)
    for branch in branches:
        branch.points_um = np.asarray(branch.raw_points_um, dtype=float).copy()
        branch.radius_um = np.asarray(branch.raw_radius_um, dtype=float).copy()
        branch.arc_length_um = np.concatenate(
            (
                [0.0],
                np.cumsum(np.linalg.norm(np.diff(branch.points_um, axis=0), axis=1)),
            )
        )
    return branches


def _junctions(roi: ROIRecord) -> list[tuple[int, np.ndarray, float]]:
    edges = np.asarray(roi.local_edges, dtype=np.int64)
    degree = np.bincount(edges.ravel(), minlength=roi.node_count)
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float)
    return [
        (int(node), positions[node], float(radii[node]))
        for node in np.flatnonzero(degree >= 3)
    ]


def _find_previous_surface(layout: "UltraliserLayout") -> Path | None:
    outputs_root = layout.run_root.parent.parent
    candidates = list(outputs_root.glob("model_generate/**/geometry/lumen_surface_m.stl"))
    if not candidates:
        return None
    return max(candidates, key=lambda path: path.stat().st_mtime)


def _load_previous_um(path: Path) -> tuple[trimesh.Trimesh, float]:
    mesh = _load_mesh(path)
    diagonal = float(np.linalg.norm(np.ptp(np.asarray(mesh.vertices), axis=0)))
    scale = 1.0e6 if diagonal < 1.0e-2 else 1.0
    if scale != 1.0:
        mesh = trimesh.Trimesh(
            vertices=np.asarray(mesh.vertices, dtype=float) * scale,
            faces=np.asarray(mesh.faces, dtype=np.int64),
            process=True,
        )
    return mesh, scale


def _worst_previous_junction(
    previous: trimesh.Trimesh | None,
    junctions: list[tuple[int, np.ndarray, float]],
) -> tuple[int, np.ndarray, float, dict[str, Any]]:
    if not junctions:
        center = np.asarray(previous.centroid if previous is not None else (0.0, 0.0, 0.0))
        return -1, center, 5.0, {"selection_method": "surface centroid; no junction in ROI"}
    if previous is None:
        selected = max(junctions, key=lambda row: row[2])
        return (*selected, {"selection_method": "largest source junction radius; old STL unavailable"})
    centers = np.asarray(previous.triangles_center, dtype=float)
    adjacency = np.asarray(previous.face_adjacency, dtype=np.int64)
    angles = np.degrees(np.asarray(previous.face_adjacency_angles, dtype=float))
    rows: list[dict[str, Any]] = []
    for node_id, center, radius in junctions:
        local = np.linalg.norm(centers - center, axis=1) <= max(6.0 * radius, 3.0)
        local_adjacency = local[adjacency[:, 0]] & local[adjacency[:, 1]] if len(adjacency) else []
        values = angles[local_adjacency] if len(adjacency) else np.empty(0)
        rows.append(
            {
                "junction_node_id": node_id,
                "normal_jump_p99_deg": float(np.percentile(values, 99)) if len(values) else -1.0,
                "local_face_count": int(np.count_nonzero(local)),
            }
        )
    worst = max(rows, key=lambda row: row["normal_jump_p99_deg"])
    selected = next(row for row in junctions if row[0] == worst["junction_node_id"])
    return (
        *selected,
        {
            "selection_method": "maximum old-surface local face-normal jump p99",
            "old_junction_metrics": rows,
            "selected_junction_node_id": int(worst["junction_node_id"]),
        },
    )


def _camera(
    roi: ROIRecord,
    node_id: int,
    center: np.ndarray,
    extent: float,
) -> tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]:
    edges = np.asarray(roi.local_edges, dtype=np.int64)
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    neighbors = sorted(
        int(second if first == node_id else first)
        for first, second in edges
        if int(first) == node_id or int(second) == node_id
    )
    directions = []
    for neighbor in neighbors:
        vector = positions[neighbor] - center
        norm = float(np.linalg.norm(vector))
        if norm > 0:
            directions.append(vector / norm)
    view = np.asarray((1.0, 1.0, 1.0), dtype=float)
    best_norm = 0.0
    for first in range(len(directions)):
        for second in range(first + 1, len(directions)):
            candidate = np.cross(directions[first], directions[second])
            norm = float(np.linalg.norm(candidate))
            if norm > best_norm:
                view, best_norm = candidate, norm
    view /= np.linalg.norm(view)
    dominant = int(np.argmax(np.abs(view)))
    if view[dominant] < 0:
        view *= -1.0
    up = directions[0] if directions else np.asarray((0.0, 0.0, 1.0))
    up = up - np.dot(up, view) * view
    if np.linalg.norm(up) < 1.0e-8:
        up = np.cross(view, np.asarray((1.0, 0.0, 0.0)))
    up /= np.linalg.norm(up)
    position = center + view * (4.0 * extent)
    return tuple(position), tuple(center), tuple(up)


def _local_mesh(mesh: trimesh.Trimesh, center: np.ndarray, extent: float) -> trimesh.Trimesh:
    mask = np.linalg.norm(np.asarray(mesh.triangles_center) - center, axis=1) <= 1.8 * extent
    face_ids = np.flatnonzero(mask)
    if len(face_ids) == 0:
        return mesh
    return mesh.submesh([face_ids], append=True, repair=False)


def _prepare_plotter(plotter: pv.Plotter) -> None:
    plotter.set_background("white")
    plotter.enable_anti_aliasing("ssaa")


def _full_surface_figure(mesh: trimesh.Trimesh, path: Path) -> Path:
    plotter = pv.Plotter(off_screen=True, window_size=(1400, 1000))
    _prepare_plotter(plotter)
    plotter.add_mesh(_polydata(mesh), color="#C84A42", smooth_shading=True, specular=0.15)
    plotter.add_text("Raw Ultraliser watertight surface (um)", font_size=13, color="black")
    plotter.view_isometric()
    plotter.reset_camera()
    plotter.show(screenshot=str(path), auto_close=True)
    return path


def _wireframe_figure(mesh: trimesh.Trimesh, path: Path) -> Path:
    plotter = pv.Plotter(off_screen=True, window_size=(1400, 1000))
    _prepare_plotter(plotter)
    plotter.add_mesh(
        _polydata(mesh),
        color="#D97757",
        show_edges=True,
        edge_color="#263238",
        line_width=0.35,
        smooth_shading=False,
    )
    plotter.add_text("Raw Ultraliser surface wireframe", font_size=13, color="black")
    plotter.view_isometric()
    plotter.reset_camera()
    plotter.show(screenshot=str(path), auto_close=True)
    return path


def _junction_figure(
    mesh: trimesh.Trimesh,
    roi: ROIRecord,
    node_id: int,
    center: np.ndarray,
    radius: float,
    path: Path,
) -> tuple[Path, tuple[Any, Any, Any], float]:
    extent = max(10.0 * radius, 6.0)
    local = _local_mesh(mesh, center, extent)
    camera = _camera(roi, node_id, center, extent)
    plotter = pv.Plotter(off_screen=True, window_size=(1400, 1000))
    _prepare_plotter(plotter)
    plotter.add_mesh(
        _polydata(local),
        color="#C84A42",
        smooth_shading=True,
        show_edges=True,
        edge_color="#6E2D2A",
        line_width=0.4,
    )
    plotter.add_text(f"Worst previous junction close-up: node {node_id}", font_size=13, color="black")
    plotter.camera_position = camera
    plotter.show(screenshot=str(path), auto_close=True)
    return path, camera, extent


def _comparison_figure(
    previous: trimesh.Trimesh,
    current: trimesh.Trimesh,
    center: np.ndarray,
    extent: float,
    camera: tuple[Any, Any, Any],
    path: Path,
) -> Path:
    previous_local = _local_mesh(previous, center, extent)
    current_local = _local_mesh(current, center, extent)
    plotter = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(1800, 900))
    for index, (surface, title, color) in enumerate(
        (
            (previous_local, "Previous final STL", "#7586A5"),
            (current_local, "Raw Ultraliser", "#C84A42"),
        )
    ):
        plotter.subplot(0, index)
        _prepare_plotter(plotter)
        plotter.add_mesh(
            _polydata(surface),
            color=color,
            smooth_shading=True,
            show_edges=True,
            edge_color="#303030",
            line_width=0.35,
        )
        plotter.add_text(title, font_size=13, color="black")
        plotter.camera_position = camera
    plotter.link_views()
    plotter.show(screenshot=str(path), auto_close=True)
    return path


def _validation_figures(
    mesh: trimesh.Trimesh,
    roi: ROIRecord,
    layout: "UltraliserLayout",
    previous_surface: Path | None,
) -> dict[str, Any]:
    full = _full_surface_figure(mesh, layout.figures / "ultraliser_full_surface.png")
    wireframe = _wireframe_figure(mesh, layout.figures / "ultraliser_wireframe.png")
    previous_path = previous_surface or _find_previous_surface(layout)
    previous = None
    previous_scale = None
    if previous_path is not None and Path(previous_path).is_file():
        previous, previous_scale = _load_previous_um(Path(previous_path))
    node_id, center, radius, selection = _worst_previous_junction(previous, _junctions(roi))
    closeup, camera, extent = _junction_figure(
        mesh,
        roi,
        node_id,
        center,
        radius,
        layout.figures / "ultraliser_junction_closeup.png",
    )
    comparison = None
    if previous is not None:
        comparison = _comparison_figure(
            previous,
            mesh,
            center,
            extent,
            camera,
            layout.figures / "previous_vs_ultraliser.png",
        )
    else:
        placeholder = pv.Plotter(off_screen=True, window_size=(1400, 900))
        _prepare_plotter(placeholder)
        placeholder.add_text("Previous final STL was not found", position="upper_left", color="black")
        placeholder.show(
            screenshot=str(layout.figures / "previous_vs_ultraliser.png"), auto_close=True
        )
        comparison = layout.figures / "previous_vs_ultraliser.png"
    report = {
        "full_surface": str(full),
        "junction_closeup": str(closeup),
        "wireframe": str(wireframe),
        "previous_vs_ultraliser": str(comparison),
        "previous_surface": str(previous_path) if previous_path is not None else None,
        "previous_surface_scale_to_um_for_visualization": previous_scale,
        "same_camera_used_for_comparison": True,
        "camera_position": [list(map(float, row)) for row in camera],
        "closeup_extent_um": extent,
        **selection,
    }
    write_json(layout.qc / "visualization_provenance.json", report)
    return report


def finalize_ultraliser_outputs(
    roi: ROIRecord,
    config: CFDLumenConfig,
    *,
    layout: "UltraliserLayout",
    raw_surface: Path,
    previous_surface: Path | None,
) -> dict[str, Any]:
    mesh, _ = _write_geometry_outputs(raw_surface, layout.final)
    surface_report = _surface_qc(mesh)
    write_json(layout.qc / "surface_qc.json", surface_report)
    branches = _source_branches(roi, config)
    fidelity_samples, fidelity_report = evaluate_radius_fidelity(mesh, branches, roi, config)
    write_csv(layout.qc / "radius_fidelity.csv", [sample.report() for sample in fidelity_samples])
    write_json(layout.qc / "radius_fidelity.json", fidelity_report)
    figures = _validation_figures(mesh, roi, layout, previous_surface)
    return {
        "surface_qc": surface_report,
        "radius_fidelity": fidelity_report,
        "figures": figures,
    }


def write_ultraliser_report(path: Path, summary: dict[str, Any]) -> Path:
    surface = summary.get("surface_qc", {})
    radius = summary.get("radius_fidelity", {})
    command_path = Path(summary.get("run_root", ".")) / "ultraliser_command.txt"
    command = command_path.read_text(encoding="utf-8").strip() if command_path.is_file() else "N/A"
    visual = summary.get("visual_assessment", "PENDING_MANUAL_REVIEW")
    lines = [
        "# Ultraliser ROI 003274 report",
        "",
        f"1. **编译/执行**：成功；commit `{summary.get('ultraliser_git_commit')}`。",
        f"2. **正式命令**：`{command}`",
        (
            "3. **Packing algorithm**："
            f"`{summary.get('packing_algorithm')}`；{summary.get('packing_algorithm_basis')}。"
        ),
        (
            "4. **voxels_per_micron**："
            f"`{summary.get('voxels_per_micron')}`；由最小直径目标和 150,000,000 voxel 上限自动确定。"
        ),
        f"5. **正式运行耗时**：`{summary.get('final_runtime_seconds')}` s。",
        (
            "6. **Surface topology**："
            f"watertight={surface.get('watertight')}，components={surface.get('component_count')}，"
            f"self_intersections={surface.get('self_intersection_count')}，"
            f"nonmanifold_edges={surface.get('nonmanifold_edge_count')}。"
        ),
        (
            "7. **Radius P95 absolute relative error**："
            f"`{radius.get('p95_absolute_relative_error')}`。"
        ),
        f"8. **Junction 锯齿人工核验**：{visual}。",
        f"9. **相对旧模型**：{summary.get('previous_comparison', visual)}。",
        f"10. **建议**：`{summary.get('decision', 'PENDING_MANUAL_REVIEW')}`。",
        "",
        "## Compatibility note",
        "",
        (
            "本地 Ultraliser vascular reader 仅接受 H5/VMV，和论文声称的 vascular SWC 直读不一致。"
            "`roi_core.swc` 是保持 µm 坐标、半径及 parent-child 的 canonical 输入；正式 executable "
            "通过 `roi_core.h5` 官方 vascular schema 读取同一几何，未修改 Ultraliser C++ 算法。"
        ),
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def finalize_visual_assessment(
    run_root: Path,
    *,
    accepted: bool,
    junction_artifact_assessment: str,
    previous_comparison: str,
) -> dict[str, Any]:
    run_root = Path(run_root).resolve()
    summary_path = run_root / "qc" / "run_summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    topology_pass = summary.get("surface_qc", {}).get("status") == "PASS"
    accepted = bool(accepted and topology_pass)
    summary["visual_assessment"] = junction_artifact_assessment
    summary["previous_comparison"] = previous_comparison
    summary["ULTRALISER_RAW_SURFACE_ACCEPTED"] = accepted
    summary["decision"] = "ADOPT_ULTRALISER" if accepted else "ULTRALISER_GEOMETRY_NOT_ACCEPTABLE"
    write_json(summary_path, summary)
    write_ultraliser_report(run_root / "report" / "ultraliser_roi003274_report.md", summary)
    return summary
