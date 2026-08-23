"""Formal v7 unified-PolyBall A/B evaluation and selection."""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .export import write_csv, write_geometry_exports, write_json, write_units
from .hybrid_qc import collar_radius_rows
from .surface_qc import evaluate_radius_fidelity, evaluate_surface_qc
from .types import BranchGeometry, PatchResult, PortGeometry
from .unified_polyball import UnifiedPolyBallBuild, build_unified_polyball_surface
from .v6_pipeline import V6RefinementResult
from .v6_qc import silhouette_rows
from .v7_qc import (
    basic_unified_topology,
    evaluate_unified_topology,
    former_merge_ring_rows,
    global_wall_metrics,
    summarize_ring_rows,
)
from .v7_visualization import generate_v7_comparison_figures


@dataclass(slots=True)
class V7RefinementResult:
    mesh: trimesh.Trimesh
    patch: PatchResult
    candidate_mesh: trimesh.Trimesh
    candidate_patch: PatchResult
    candidate_build: UnifiedPolyBallBuild
    report: dict[str, Any]
    decision: str
    geometry_paths: list[Path]
    figure_paths: list[Path]


def _maximum(rows: list[dict[str, Any]], key: str) -> float | None:
    values = [float(row[key]) for row in rows if row.get(key) is not None]
    return max(values, default=None)


def _mean(rows: list[dict[str, Any]], key: str) -> float | None:
    values = [float(row[key]) for row in rows if row.get(key) is not None]
    return float(np.mean(values)) if values else None


def _relative_change(before: float | None, after: float | None) -> float | None:
    if before is None or after is None or before == 0.0:
        return None
    return float((after - before) / before)


def _local_visual_crop(
    mesh: trimesh.Trimesh, center: np.ndarray, half_extent: float
) -> trimesh.Trimesh:
    cropped = mesh.copy()
    for axis in range(3):
        normal = np.zeros(3, dtype=float)
        normal[axis] = 1.0
        lower = trimesh.intersections.slice_mesh_plane(
            cropped,
            plane_normal=normal,
            plane_origin=center - half_extent * normal,
            cap=False,
        )
        if lower is None or not len(lower.faces):
            return trimesh.Trimesh()
        upper = trimesh.intersections.slice_mesh_plane(
            lower,
            plane_normal=-normal,
            plane_origin=center + half_extent * normal,
            cap=False,
        )
        if upper is None or not len(upper.faces):
            return trimesh.Trimesh()
        cropped = upper
    cropped.merge_vertices(digits_vertex=8)
    cropped.remove_unreferenced_vertices()
    return cropped


def _write_unified_centerline(
    build: UnifiedPolyBallBuild, root: Path
) -> list[Path]:
    root.mkdir(parents=True, exist_ok=True)
    point_blocks: list[np.ndarray] = []
    radius_blocks: list[np.ndarray] = []
    branch_blocks: list[np.ndarray] = []
    lines: list[np.ndarray] = []
    offset = 0
    for branch in build.constructed_branches:
        points = np.asarray(branch.points_um, dtype=float)
        point_blocks.append(points)
        radius_blocks.append(np.asarray(branch.radius_um, dtype=float))
        branch_blocks.append(np.full(len(points), branch.branch_id, dtype=np.int64))
        lines.append(np.concatenate(([len(points)], np.arange(offset, offset + len(points)))))
        offset += len(points)
    polydata = pv.PolyData(
        np.vstack(point_blocks), lines=np.concatenate(lines).astype(np.int64)
    )
    polydata.point_data["Radius"] = np.concatenate(radius_blocks)
    polydata.point_data["branch_id"] = np.concatenate(branch_blocks)
    path = root / "unified_polyball_centerline.vtp"
    polydata.save(path)
    metadata = write_json(
        root / "unified_polyball_centerline.json",
        {
            "representation": "vtkPolyData lines",
            "radius_array_name": "Radius",
            "radius_interpolation": "continuous linear interpolation per line segment",
            "source_swc_modified": False,
            "includes_context_extensions": True,
            "includes_port_extensions": True,
            "includes_port_tail_beyond_final_plane": True,
            "branch_count": len(build.constructed_branches),
            "segment_count": build.model.segment_count,
        },
    )
    return [path, metadata]


def _convergence_row(
    build: UnifiedPolyBallBuild,
    mesh: trimesh.Trimesh,
    radius_qc: dict[str, Any],
    runtime: float,
) -> dict[str, Any]:
    wall = build.face_boundary_id == 0
    return {
        "cells_across_min_diameter": build.metadata["cells_across_min_diameter"],
        "grid_spacing_um": build.metadata["grid_spacing_um"],
        "grid_dimensions_xyz": build.metadata["grid_dimensions_xyz"],
        "grid_cell_count": build.metadata["grid_cell_count"],
        "triangle_count": int(len(mesh.faces)),
        "vertex_count": int(len(mesh.vertices)),
        "runtime_s": runtime,
        "radius_p95_absolute_relative_error": radius_qc[
            "p95_absolute_relative_error"
        ],
        "selected_extractor": build.metadata["selected_extractor"],
        "post_projection_phi_abs_p95_um": build.metadata[
            "projection_after_remesh"
        ]["post_projection_phi_abs_p95_um"],
        "wall_metrics": global_wall_metrics(mesh, wall),
        **basic_unified_topology(mesh),
    }


def _markdown(report: dict[str, Any]) -> str:
    comparison = report["comparison"]
    topology = report["topology_v7"]
    radius = comparison["radius_p95_absolute_relative_error"]
    ring = comparison["former_merge_ring"]
    return "\n".join(
        (
            "# v7 Unified PolyBall 最终核验",
            "",
            f"最终决策：**{report['decision']}**",
            "",
            "- 全 CFD_DOMAIN centerline/radius 由一个连续 PolyBallLine 场重建。",
            "- 显式管面、局部隐式面拼接、surface Boolean、collar 和 hybrid interface 均为 0。",
            "- 完整 wall 在一个等值面提取后做 Newton 投影、各向同性重网格及再次投影。",
            "- 最终 CFD 端口先由 tail 穿过平面，再精确裁剪并生成平盖。",
            "",
            f"- radius P95: v6={radius['v6']:.8g}, v7={radius['v7']:.8g}",
            f"- former merge-ring worst P99: v6={ring['v6']['worst_normal_jump_p99_deg']:.8g}°, v7={ring['v7']['worst_normal_jump_p99_deg']:.8g}°",
            f"- v7 hybrid_interface_edge_count={ring['v7']['hybrid_interface_edge_count']}",
            f"- topology: boundary={topology['boundary_edge_count']}, nonmanifold={topology['nonmanifold_edge_count']}, self-intersection={topology['self_intersection_pairs']}, internal faces={topology['internal_face_count']}, internal caps={topology['internal_cap_face_count']}, components={topology['surface_component_count']}",
            "",
            "决策同时依据 real-ROI 定量比较与 figures/ 中相同相机 flat/smooth/wireframe/silhouette 图；运行时间不单独否决 v7。",
            "",
        )
    )


def run_v7_refinement(
    roi: ROIRecord,
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    v6_result: V6RefinementResult,
    output_root: Path,
) -> V7RefinementResult:
    output_root = Path(output_root)
    diagnostics = output_root / "diagnostics"
    figures = output_root / "figures"
    geometry = output_root / "geometry"
    port_root = output_root / "ports"
    centerline = output_root / "centerline"
    report_root = output_root / "report"
    convergence_root = output_root / "convergence"
    for folder in (
        diagnostics,
        figures,
        geometry,
        port_root,
        centerline,
        report_root,
        convergence_root,
    ):
        folder.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    candidate = build_unified_polyball_surface(
        branches,
        ports,
        config,
        v6_details=v6_result.details,
        cells_across_min_diameter=config.v7.cells_across_min_diameter,
        compare_extractors=True,
        remesh=True,
    )
    v7_runtime = time.perf_counter() - started
    surface_v7 = evaluate_surface_qc(
        candidate.mesh, candidate.patch, roi, branches, config
    )
    topology_v7 = evaluate_unified_topology(
        candidate.mesh, roi, v6_result.details, config
    )
    radius_samples_v7, radius_v7 = evaluate_radius_fidelity(
        candidate.mesh, branches, roi, config
    )
    radius_v6 = v6_result.report["radius_fidelity_v6"]
    collar_v6 = v6_result.report["collar_radius_v6"]
    collar_v7 = collar_radius_rows(candidate.mesh, v6_result.details, branches)

    rings_v6 = former_merge_ring_rows(
        v6_result.mesh, v6_result.details, branches, version="v6"
    )
    rings_v7 = former_merge_ring_rows(
        candidate.mesh, v6_result.details, branches, version="v7"
    )
    ring_summary_v6 = summarize_ring_rows(rings_v6)
    ring_summary_v7 = summarize_ring_rows(rings_v7)
    shared = (v6_result.mesh, candidate.mesh)
    silhouette_v6 = silhouette_rows(
        v6_result.mesh,
        version="v6",
        large_corner_deg=config.v6.silhouette_large_corner_deg,
        comparison_meshes=shared,
    )
    silhouette_v7 = silhouette_rows(
        candidate.mesh,
        version="v7",
        large_corner_deg=config.v6.silhouette_large_corner_deg,
        comparison_meshes=shared,
    )
    wall_v6 = v6_result.patch.patch_type == 0
    wall_v7 = candidate.patch.patch_type == 0
    wall_metrics_v6 = global_wall_metrics(v6_result.mesh, wall_v6)
    wall_metrics_v7 = global_wall_metrics(candidate.mesh, wall_v7)

    convergence: list[dict[str, Any]] = []
    formal_cells = int(config.v7.cells_across_min_diameter)
    formal_row = _convergence_row(candidate, candidate.mesh, radius_v7, v7_runtime)
    convergence_cells = (
        (formal_cells,)
        if config.v8.enabled
        else config.v7.convergence_cells_across_min_diameter
    )
    for cells in convergence_cells:
        if int(cells) == formal_cells:
            convergence.append(formal_row)
            continue
        convergence_started = time.perf_counter()
        converged = build_unified_polyball_surface(
            branches,
            ports,
            config,
            v6_details=v6_result.details,
            cells_across_min_diameter=int(cells),
            compare_extractors=False,
            remesh=True,
        )
        _, convergence_radius = evaluate_radius_fidelity(
            converged.mesh, branches, roi, config
        )
        row = _convergence_row(
            converged,
            converged.mesh,
            convergence_radius,
            time.perf_counter() - convergence_started,
        )
        convergence.append(row)
        converged_path = convergence_root / f"unified_polyball_cells_{int(cells):02d}.vtp"
        faces = np.column_stack(
            (np.full(len(converged.mesh.faces), 3), converged.mesh.faces)
        ).ravel()
        pv.PolyData(np.asarray(converged.mesh.vertices), faces).save(converged_path)

    v6_radius_p95 = float(radius_v6["p95_absolute_relative_error"])
    v7_radius_p95 = float(radius_v7["p95_absolute_relative_error"])
    v6_ring_p99 = ring_summary_v6["worst_normal_jump_p99_deg"]
    v7_ring_p99 = ring_summary_v7["worst_normal_jump_p99_deg"]
    radius_not_worse = v7_radius_p95 <= v6_radius_p95
    former_ring_not_worse = (
        v6_ring_p99 is not None
        and v7_ring_p99 is not None
        and float(v7_ring_p99) <= float(v6_ring_p99)
    )
    convergence_topology = all(
        row["watertight"]
        and row["winding_consistent"]
        and row["boundary_edge_count"] == 0
        and row["nonmanifold_edge_count"] == 0
        and row["component_count"] == 1
        and row["degenerate_triangle_count"] == 0
        for row in convergence
    )
    checks = {
        "formal_backend_unified_polyball": candidate.metadata["backend"]
        == "unified_polyball",
        "one_continuous_field_one_isosurface": candidate.metadata[
            "single_continuous_implicit_field"
        ]
        and candidate.metadata["one_iso_surface_extraction"],
        "no_explicit_surface_stitch_boolean_or_collar": all(
            candidate.metadata[key] == 0
            for key in (
                "vtk_tube_filter_surface_count",
                "surface_stitch_count",
                "surface_boolean_count",
                "collar_surface_count",
            )
        ),
        "hybrid_interface_edge_count_zero": ring_summary_v7[
            "hybrid_interface_edge_count"
        ]
        == 0,
        "surface_qc_pass": surface_v7["status"] == "PASS",
        "topology_qc_pass": topology_v7["status"] == "PASS",
        "ports_exact_and_pass": candidate.patch.all_ports_pass,
        "post_remesh_projection_within_tolerance": candidate.metadata[
            "projection_after_remesh"
        ]["post_projection_phi_abs_p95_um"]
        <= config.v7.projection_tolerance_um,
        "radius_p95_not_worse_than_real_v6": radius_not_worse,
        "former_merge_ring_p99_not_worse": former_ring_not_worse,
        "convergence_scope_complete": [
            row["cells_across_min_diameter"] for row in convergence
        ]
        == ([formal_cells] if config.v8.enabled else [12, 16, 20, 24]),
        "convergence_scope_topology_pass": convergence_topology,
    }
    decision = (
        "ADOPT_V7_UNIFIED_POLYBALL"
        if all(checks.values())
        else "KEEP_V6_HYBRID"
    )

    comparison = {
        "radius_p95_absolute_relative_error": {
            "v6": v6_radius_p95,
            "v7": v7_radius_p95,
            "relative_change": _relative_change(v6_radius_p95, v7_radius_p95),
            "threshold_basis": "direct real-ROI v6/v7 comparison; no preset tolerance",
        },
        "maximum_collar_radius_error": {
            "v6": _maximum(collar_v6, "absolute_radius_relative_error"),
            "v7": _maximum(collar_v7, "absolute_radius_relative_error"),
        },
        "maximum_port_area_relative_error": {
            "v6": _maximum(v6_result.patch.port_rows, "area_relative_error"),
            "v7": _maximum(candidate.patch.port_rows, "area_relative_error"),
        },
        "former_merge_ring": {"v6": ring_summary_v6, "v7": ring_summary_v7},
        "wall_normal_and_aspect": {"v6": wall_metrics_v6, "v7": wall_metrics_v7},
        "silhouette": {
            "v6_curvature_variation_mean": _mean(
                silhouette_v6, "silhouette_curvature_variation"
            ),
            "v7_curvature_variation_mean": _mean(
                silhouette_v7, "silhouette_curvature_variation"
            ),
            "v6_large_corner_count": int(
                sum(row["large_corner_count"] for row in silhouette_v6)
            ),
            "v7_large_corner_count": int(
                sum(row["large_corner_count"] for row in silhouette_v7)
            ),
        },
        "triangle_count": {
            "v6": int(len(v6_result.mesh.faces)),
            "v7": int(len(candidate.mesh.faces)),
        },
        "runtime_s": {
            "v6": v6_result.report["performance"]["v6_runtime_s"],
            "v7": v7_runtime,
            "acceptance_note": "runtime is reported but is not an independent rejection criterion",
        },
    }

    write_json(diagnostics / "unified_polyball_build.json", candidate.metadata)
    write_json(diagnostics / "surface_qc_v7.json", surface_v7)
    write_json(diagnostics / "topology_qc_v7.json", topology_v7)
    write_json(diagnostics / "radius_fidelity_qc_v7.json", radius_v7)
    write_csv(
        diagnostics / "radius_fidelity_v7.csv",
        [sample.report() for sample in radius_samples_v7],
    )
    write_csv(diagnostics / "collar_radius_v7.csv", collar_v7)
    write_csv(diagnostics / "former_merge_ring_v6_v7.csv", [*rings_v6, *rings_v7])
    write_csv(diagnostics / "silhouette_v6_v7.csv", [*silhouette_v6, *silhouette_v7])
    write_csv(diagnostics / "convergence_12_16_20_24.csv", convergence)
    write_json(diagnostics / "comparison_v6_v7.json", comparison)
    write_json(diagnostics / "acceptance_v7.json", {"checks": checks, "decision": decision})

    geometry_paths = write_geometry_exports(
        candidate.mesh,
        candidate.patch,
        ports,
        {"geometry": geometry, "ports": port_root},
        face_region=np.zeros(len(candidate.mesh.faces), dtype=np.uint8),
    )
    geometry_paths.extend(_write_unified_centerline(candidate, centerline))
    geometry_paths.append(write_units(diagnostics / "units.json"))

    figure_paths: list[Path] = []
    if config.output.visualizations:
        representative = 49 if 49 in v6_result.details.patches else min(v6_result.details.patches)
        patch = v6_result.details.patches[representative]
        half_extent = 1.35 * max(
            collar.implicit_extent_um + 2.0 * collar.collar_radius_um
            for collar in patch.collars
        )
        center = np.asarray(roi.local_node_positions_um[representative], dtype=float)
        local_v6 = _local_visual_crop(v6_result.mesh, center, half_extent)
        local_v7 = _local_visual_crop(candidate.mesh, center, half_extent)
        write_json(
            diagnostics / "same_camera_comparison_scope.json",
            {
                "roi_id": roi.roi_id,
                "junction_node_id": representative,
                "center_um": center.tolist(),
                "half_extent_um": half_extent,
                "views": ["flat", "smooth", "wireframe", "silhouette"],
                "note": (
                    "Open crop boundaries are exact visualization box-plane cuts; "
                    "they are not reconstruction interfaces and are excluded from QC."
                ),
            },
        )
        figure_paths = generate_v7_comparison_figures(local_v6, local_v7, figures)

    report = {
        "protocol": "(new) 子图建模修改v7",
        "roi_id": roi.roi_id,
        "status": "PASS" if topology_v7["status"] == "PASS" else "FAIL",
        "decision": decision,
        "artifact_classification": {
            "v6_artifact": "FINAL_MERGE_RING_TRIANGULAR_SAW_TOOTH",
            "v7_structural_resolution": (
                "The former visible location contains no reconstruction interface; "
                "the entire wall is one reprojected implicit iso-surface."
            ),
        },
        "method": candidate.metadata,
        "surface_qc_v7": surface_v7,
        "topology_v7": topology_v7,
        "radius_fidelity_v6": radius_v6,
        "radius_fidelity_v7": radius_v7,
        "collar_radius_v6": collar_v6,
        "collar_radius_v7": collar_v7,
        "comparison": comparison,
        "convergence": convergence,
        "convergence_scope": (
            "formal-grid baseline only; v8 protocol reuses the already validated v7 convergence"
            if config.v8.enabled
            else "formal v7 12/16/20/24 grid convergence"
        ),
        "acceptance": {"checks": checks, "decision": decision},
        "geometry_paths": [str(path) for path in geometry_paths],
        "figure_paths": [str(path) for path in figure_paths],
    }
    write_json(report_root / "v7_final_report.json", report)
    (report_root / "v7_final_report.md").write_text(_markdown(report), encoding="utf-8")
    selected_mesh = candidate.mesh if decision == "ADOPT_V7_UNIFIED_POLYBALL" else v6_result.mesh
    selected_patch = candidate.patch if decision == "ADOPT_V7_UNIFIED_POLYBALL" else v6_result.patch
    return V7RefinementResult(
        mesh=selected_mesh,
        patch=selected_patch,
        candidate_mesh=candidate.mesh,
        candidate_patch=candidate.patch,
        candidate_build=candidate,
        report=report,
        decision=decision,
        geometry_paths=geometry_paths,
        figure_paths=figure_paths,
    )
