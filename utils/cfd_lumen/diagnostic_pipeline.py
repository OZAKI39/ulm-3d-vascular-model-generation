"""Strict execution of the `(new) subgraph modeling v2` diagnosis protocol."""

from __future__ import annotations

import logging
import time
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh

from utils.sampling.sampling_types import ROIRecord

from .collision_qc import detect_nonadjacent_collisions
from .config import CFDLumenConfig
from .diagnostic_visualization import generate_diagnostic_figures
from .export import write_csv, write_json
from .geometry_diagnostics import (
    diagnose_junctions,
    diagnose_ports,
    explicit_union_for_diagnostics,
    save_diagnostic_primitives,
)
from .geometry_preprocess import resample_branches, validate_and_extract_branches
from .lumen_builder import build_lumen_primitives, build_lumen_surface
from .mesh_defects import diagnose_mesh_defects
from .port_geometry import construct_port_geometry
from .surface_qc import _section_polygon, evaluate_radius_fidelity
from .synthetic_controls import run_synthetic_controls
from .types import BranchGeometry, CFDRunLayout, PortGeometry


LOGGER = logging.getLogger("cfd_lumen")


PORT_FIELDS = (
    "port_id", "cut_port_id", "exact_cut_x_um", "exact_cut_y_um", "exact_cut_z_um",
    "cut_radius_um", "source_global_edge", "source_edge_start_x_um", "source_edge_start_y_um",
    "source_edge_start_z_um", "source_edge_end_x_um", "source_edge_end_y_um", "source_edge_end_z_um",
    "source_edge_start_radius_um", "source_edge_end_radius_um", "cut_projection_t",
    "radius_recomputed_um", "radius_interpolation_error_um", "branch_id",
    "endpoint_position_error_um", "endpoint_radius_error_um", "endpoint_radius_relative_error",
    "section_offset_epsilon_um", "source_radius_um", "branch_mesh_area_um2",
    "branch_mesh_radius_um", "extension_radius_input_um", "extension_mesh_area_um2",
    "extension_mesh_radius_um", "step_radius_relative_error", "step_area_relative_error",
    "branch_tube_sides", "extension_cylinder_sides", "cross_section_vertex_count_branch",
    "cross_section_vertex_count_extension", "cross_section_polygon_symmetric_difference_relative",
    "expected_overlap_um", "actual_axial_overlap_um", "intersection_volume_um3",
    "intersection_boolean_error", "intersection_boolean_runtime_s", "root_cause_candidates",
)


def _mesh_vtp(mesh: trimesh.Trimesh, path: Path) -> Path:
    faces = np.column_stack((np.full(len(mesh.faces), 3, dtype=np.int64), mesh.faces)).ravel()
    pv.PolyData(np.asarray(mesh.vertices), faces).save(path)
    return path


def _final_port_step(
    mesh: trimesh.Trimesh,
    ports: list[PortGeometry],
) -> float | None:
    errors: list[float] = []
    for port in ports:
        epsilon = max(0.02 * port.radius_um, 1.0e-5)
        before = _section_polygon(
            mesh,
            port.original_position_um - epsilon * port.outward_tangent,
            port.outward_tangent,
        )
        after = _section_polygon(
            mesh,
            port.original_position_um + epsilon * port.outward_tangent,
            port.outward_tangent,
        )
        if before and after:
            errors.append(abs(before[0] - after[0]) / (0.5 * (before[0] + after[0])))
    return max(errors, default=None)


def _backend_row(
    name: str,
    mesh: trimesh.Trimesh,
    runtime_s: float,
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    roi: ROIRecord,
    config: CFDLumenConfig,
    junctions: list[tuple[int, np.ndarray, float]],
    *,
    implicit_grid: dict[str, Any] | None = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    defects, artifacts = diagnose_mesh_defects(mesh, junctions)
    _, radius_qc = evaluate_radius_fidelity(mesh, branches, roi, config)
    row = {
        "backend": name,
        "watertight": defects["watertight"],
        "boundary_edge_count": defects["boundary_edge_count"],
        "non_manifold_edge_count": defects["non_manifold_edge_count"],
        "surface_connected_component_count": defects["surface_connected_component_count"],
        "self_intersection_count": defects["self_intersection_count"],
        "suspected_internal_face_count": defects["suspected_internal_face_count"],
        "is_winding_consistent": defects["is_winding_consistent"],
        "radius_p95_absolute_relative_error": radius_qc["p95_absolute_relative_error"],
        "port_step_area_relative_error_max": _final_port_step(mesh, ports),
        "triangle_count": int(len(mesh.faces)),
        "vertex_count": int(len(mesh.vertices)),
        "surface_area_um2": float(mesh.area),
        "enclosed_volume_um3": float(abs(mesh.volume)) if mesh.is_watertight else None,
        "aspect_ratio_p95": defects["triangle_quality"]["aspect_ratio_p95"],
        "aspect_ratio_max": defects["triangle_quality"]["aspect_ratio_max"],
        "degenerate_triangle_count": defects["triangle_quality"]["degenerate_triangle_count"],
        "runtime_s": runtime_s,
        "implicit_grid": implicit_grid,
    }
    return row, {"defects": defects, "artifacts": artifacts, "radius_qc": radius_qc}


def _root_causes(
    port_rows: list[dict[str, Any]],
    junction_rows: list[dict[str, Any]],
    defects: dict[str, Any],
    collision_qc: dict[str, Any],
) -> tuple[list[str], list[str]]:
    port_causes = {
        cause
        for row in port_rows
        for cause in str(row["root_cause_candidates"]).split(";")
        if cause
    }
    max_step = max((row["step_area_relative_error"] or 0.0 for row in port_rows), default=0.0)
    if not port_causes and max_step <= 0.01:
        port_causes.add("PORT_NORMAL_OR_SHADING_ARTIFACT")
    junction_causes: set[str] = set()
    if defects["boundary_edge_count"]:
        junction_causes.add("BOOLEAN_OPEN_CRACK")
    if defects["non_manifold_edge_count"]:
        junction_causes.add("BOOLEAN_NONMANIFOLD")
    if defects["self_intersection_count"]:
        junction_causes.add("BOOLEAN_SELF_INTERSECTION")
    if defects["suspected_internal_face_count"]:
        junction_causes.add("BOOLEAN_INTERNAL_SURFACE")
    if sum(row["suspected_internal_cap_face_count_after_boolean"] for row in junction_rows):
        junction_causes.add("INTERNAL_CAP_SURVIVED_BOOLEAN")
    if any(
        row["intersection_volume_um3"] is not None and row["intersection_volume_um3"] <= 1.0e-12
        for row in junction_rows
    ):
        junction_causes.add("JUNCTION_INSUFFICIENT_OVERLAP")
    if any(row["junction_to_adjacent_radius_ratio"] < 0.75 for row in junction_rows):
        junction_causes.add("JUNCTION_RADIUS_INCOMPATIBILITY")
    junction_quality = defects.get("junction_triangle_quality", [])
    if any(
        (row.get("aspect_ratio_max") or 0.0) >= 20.0 or (row.get("minimum_angle_deg") or 180.0) <= 1.0
        for row in junction_quality
    ):
        junction_causes.add("BOOLEAN_SLIVER_TRIANGLES")
    if collision_qc["hard_collision_count"]:
        junction_causes.add("SOURCE_GEOMETRY_COLLISION")
    if not junction_causes:
        junction_causes.add("NORMAL_OR_SHADING_ARTIFACT")
    return sorted(port_causes), sorted(junction_causes)


def _recommendations(port_causes: list[str], junction_causes: list[str]) -> list[str]:
    recommendations: list[str] = []
    if "PORT_POLYGON_RESOLUTION_MISMATCH" in port_causes:
        recommendations.append("P0: align branch/extension circumferential frames or continue the terminal tube as one primitive; retain source radius.")
    if "PORT_INSUFFICIENT_OVERLAP" in port_causes:
        recommendations.append("P0: increase only the affected port's positive overlap after verifying the reported intersection volume.")
    if "BOOLEAN_OPEN_CRACK" in junction_causes or "BOOLEAN_NONMANIFOLD" in junction_causes:
        recommendations.append("P0: reject the affected explicit surface from CFD and repair the local Boolean construction before volume meshing.")
    elif "BOOLEAN_SELF_INTERSECTION" in junction_causes or "BOOLEAN_INTERNAL_SURFACE" in junction_causes:
        recommendations.append("P0: reject the affected explicit surface from CFD because internal/self-intersecting faces violate CFD-ready geometry.")
    if "BOOLEAN_SLIVER_TRIANGLES" in junction_causes:
        recommendations.append("P1: replace only the affected junction neighborhood with a local implicit union/remesh; do not smooth or alter SWC radii.")
    if "NORMAL_OR_SHADING_ARTIFACT" in junction_causes or "PORT_NORMAL_OR_SHADING_ARTIFACT" in port_causes:
        recommendations.append("P2: use consistently recomputed display normals; keep the formal geometry unchanged unless topology metrics fail.")
    if not recommendations:
        recommendations.append("P2: no geometry mutation is justified; retain the current reconstruction and archive this diagnosis as baseline evidence.")
    return recommendations


def _report_markdown(
    roi_id: str,
    port_rows: list[dict[str, Any]],
    junction_rows: list[dict[str, Any]],
    defects: dict[str, Any],
    backend_rows: list[dict[str, Any]],
    port_causes: list[str],
    junction_causes: list[str],
    collision_qc: dict[str, Any],
    recommendations: list[str],
    synthetic: dict[str, Any],
) -> str:
    max_position = max((row["endpoint_position_error_um"] for row in port_rows), default=float("nan"))
    max_radius = max((row["endpoint_radius_error_um"] for row in port_rows), default=float("nan"))
    max_area = max((row["step_area_relative_error"] for row in port_rows if row["step_area_relative_error"] is not None), default=float("nan"))
    min_port_overlap = min((row["intersection_volume_um3"] for row in port_rows if row["intersection_volume_um3"] is not None), default=float("nan"))
    port_real = max_position > 1.0e-8 or max_radius > 1.0e-8 or max_area > 0.01
    open_crack = defects["boundary_edge_count"] > 0
    source_issue = collision_qc["hard_collision_count"] > 0
    explicit, implicit = backend_rows
    implicit_removes_defect = (
        implicit["boundary_edge_count"] + implicit["non_manifold_edge_count"] + implicit["self_intersection_count"]
        < explicit["boundary_edge_count"] + explicit["non_manifold_edge_count"] + explicit["self_intersection_count"]
        or implicit["aspect_ratio_p95"] < explicit["aspect_ratio_p95"]
    )
    port_lines = [
        f"- Port {row['port_id']}: source={row['source_radius_um']:.9g} um, branch={row['branch_mesh_radius_um']:.9g} um, "
        f"extension={row['extension_mesh_radius_um']:.9g} um, area difference={row['step_area_relative_error']:.6g}, "
        f"overlap volume={row['intersection_volume_um3']:.9g} um^3"
        for row in port_rows
    ]
    junction_lines = [
        f"- Junction {row['junction_node_id']} / branch {row['branch_id']}: overlap={row['intersection_volume_um3']:.9g} um^3 "
        f"({row['intersection_over_junction_volume']:.6g} of junction solid), A(1D)={row['section_1D_area_um2']:.9g} um^2, "
        f"A(2D)={row['section_2D_area_um2']:.9g} um^2, confirmed internal-cap faces={row['suspected_internal_cap_face_count_after_boolean']}"
        for row in junction_rows
    ]
    scalar_control = next(
        row for row in synthetic["controls"] if row["name"] == "vtk_absolute_radius_scalar_semantics"
    )
    pair_rows = defects["self_intersection_pairs"]
    pair_location_summary = (
        ", ".join(
            f"junction {node}: {sum(pair['nearest_junction_node_id'] == node for pair in pair_rows)} pairs"
            for node in sorted({pair["nearest_junction_node_id"] for pair in pair_rows})
        )
        if pair_rows
        else "none"
    )
    return "\n".join(
        [
            f"# Geometry diagnostic report — `{roi_id}`",
            "",
            "## 1. Port step",
            "",
            f"视觉台阶是否为真实 radius/area discontinuity：**{'YES' if port_real else 'NO'}**。",
            f"最大 endpoint position error = `{max_position:.9g} um`；最大 endpoint radius error = `{max_radius:.9g} um`；最大 area difference = `{max_area:.9g}`。",
            f"最小 positive-overlap intersection volume = `{min_port_overlap:.9g} um^3`。",
            *port_lines,
            f"VTK radius scalar synthetic control：target=`{scalar_control['target_radius_um']:.9g} um`，measured=`{scalar_control['measured_mesh_radius_um']:.9g} um`，status=`{scalar_control['status']}`。",
            f"最可能原因：`{', '.join(port_causes)}`。",
            "",
            "## 2. Junction crack",
            "",
            f"白缝是否是真正 open crack：**{'YES' if open_crack else 'NO'}**。",
            f"boundary edges = `{defects['boundary_edge_count']}`；non-manifold edges = `{defects['non_manifold_edge_count']}`；self intersections = `{defects['self_intersection_count']}`；suspected internal faces = `{defects['suspected_internal_face_count']}`。",
            f"最小 junction/branch overlap volume = `{min((row['intersection_volume_um3'] for row in junction_rows if row['intersection_volume_um3'] is not None), default=float('nan')):.9g} um^3`。",
            *junction_lines,
            f"自交位置归属：`{pair_location_summary}`。",
            f"最可能原因：`{', '.join(junction_causes)}`。",
            "",
            "## 3. Source SWC",
            "",
            f"是否有证据来自 source centerline/radius：**{'YES' if source_issue else 'NO'}**（hard non-adjacent collisions = `{collision_qc['hard_collision_count']}`）。未修改 SWC 或 radius。",
            "",
            "## 4. Reconstruction stage",
            "",
            f"证据定位：port=`{', '.join(port_causes)}`；junction=`{', '.join(junction_causes)}`。",
            "",
            "## 5. Explicit vs implicit",
            "",
            f"定量判断：**{'implicit 消除了 explicit 的局部拓扑/相交缺陷，但并非所有指标全面更优' if implicit_removes_defect else '没有证据表明 implicit 全面优于 explicit'}**。",
            f"Explicit: radius P95=`{explicit['radius_p95_absolute_relative_error']:.9g}`, runtime=`{explicit['runtime_s']:.6g} s`, triangles=`{explicit['triangle_count']}`, boundaries=`{explicit['boundary_edge_count']}`, non-manifold=`{explicit['non_manifold_edge_count']}`, self-intersections=`{explicit['self_intersection_count']}`, degenerate triangles=`{explicit['degenerate_triangle_count']}`。",
            f"Implicit: radius P95=`{implicit['radius_p95_absolute_relative_error']:.9g}`, runtime=`{implicit['runtime_s']:.6g} s`, triangles=`{implicit['triangle_count']}`, boundaries=`{implicit['boundary_edge_count']}`, non-manifold=`{implicit['non_manifold_edge_count']}`, self-intersections=`{implicit['self_intersection_count']}`, degenerate triangles=`{implicit['degenerate_triangle_count']}`。",
            "",
            "## 6. Recommended fixes",
            "",
            *[f"- {item}" for item in recommendations],
            "",
        ]
    )


def diagnose_roi(
    roi: ROIRecord,
    config: CFDLumenConfig,
    run_layout: CFDRunLayout,
    synthetic: dict[str, Any],
) -> dict[str, Any]:
    started = time.perf_counter()
    root = run_layout.run_root / "diagnostics" / roi.roi_id
    primitives_dir = root / "primitives"
    figures_dir = root / "figures"
    logs_dir = root / "logs"
    for folder in (root, primitives_dir, figures_dir, logs_dir):
        folder.mkdir(parents=True, exist_ok=False)
    log_path = logs_dir / "geometry_diagnostics.log"
    log_lines = [f"ROI={roi.roi_id}", "workers=1", "source_geometry_mutated=false"]

    branches, pre_qc = validate_and_extract_branches(roi, config)
    branches = resample_branches(branches, config)
    ports = construct_port_geometry(roi, config)
    collisions, collision_qc = detect_nonadjacent_collisions(branches, config)
    log_lines.append(f"pre_geometry={pre_qc}")
    log_lines.append(f"collision_qc={collision_qc}")

    explicit_started = time.perf_counter()
    primitives = build_lumen_primitives(branches, roi, ports, config)
    explicit_mesh, boolean_runtime = explicit_union_for_diagnostics(primitives)
    explicit_runtime = time.perf_counter() - explicit_started
    save_diagnostic_primitives(primitives, explicit_mesh, primitives_dir)

    port_rows, port_summary = diagnose_ports(roi, branches, ports, primitives, config)
    degree = np.bincount(np.asarray(roi.local_edges).ravel(), minlength=roi.node_count)
    junctions = [
        (int(node_id), np.asarray(roi.local_node_positions_um[node_id]), float(roi.local_node_radius_um[node_id]))
        for node_id in np.flatnonzero(degree >= 3)
    ]
    explicit_row, explicit_detail = _backend_row(
        "explicit_manifold", explicit_mesh, explicit_runtime, branches, ports, roi, config, junctions
    )
    junction_rows, junction_summary = diagnose_junctions(
        roi,
        branches,
        primitives,
        explicit_mesh,
        config,
        boolean_runtime_s=boolean_runtime,
        confirmed_internal_face_ids=set(
            map(int, explicit_detail["artifacts"]["suspected_internal_faces"])
        ),
    )
    implicit_started = time.perf_counter()
    implicit_build = build_lumen_surface(branches, roi, ports, config, backend="implicit")
    implicit_runtime = time.perf_counter() - implicit_started
    implicit_row, implicit_detail = _backend_row(
        "implicit_fallback", implicit_build.mesh, implicit_runtime, branches, ports, roi, config,
        junctions, implicit_grid=implicit_build.implicit_grid,
    )
    backend_rows = [explicit_row, implicit_row]
    _mesh_vtp(implicit_build.mesh, primitives_dir / "implicit_comparison_surface.vtp")

    defects = explicit_detail["defects"]
    artifacts = explicit_detail["artifacts"]
    write_csv(root / "port_diagnostics.csv", port_rows, fieldnames=PORT_FIELDS)
    write_csv(root / "junction_diagnostics.csv", junction_rows)
    write_json(root / "mesh_defects.json", defects)
    write_csv(root / "backend_comparison.csv", backend_rows)
    write_json(root / "synthetic_controls.json", synthetic)
    write_json(root / "collision_qc.json", {**collision_qc, "events": [event.report() for event in collisions]})

    figure_paths = generate_diagnostic_figures(
        explicit_mesh, implicit_build.mesh, branches, ports, primitives, port_rows, junction_rows,
        defects, artifacts, backend_rows, figures_dir,
    )
    port_causes, junction_causes = _root_causes(
        port_rows, junction_rows, defects, collision_qc
    )
    recommendations = _recommendations(port_causes, junction_causes)
    report = _report_markdown(
        roi.roi_id, port_rows, junction_rows, defects, backend_rows, port_causes,
        junction_causes, collision_qc, recommendations,
        synthetic,
    )
    report_path = root / "diagnostic_report.md"
    report_path.write_text(report, encoding="utf-8")
    summary = {
        "roi_id": roi.roi_id,
        "protocol": "(new) 子图建模修改v2",
        "source_geometry_mutated": False,
        "port_issue": {
            "detected": True,
            "affected_ports": [row["port_id"] for row in port_rows],
            "root_causes": port_causes,
            "evidence": port_summary,
        },
        "junction_issue": {
            "detected": True,
            "affected_junctions": sorted({row["junction_node_id"] for row in junction_rows}),
            "root_causes": junction_causes,
            "evidence": {
                **junction_summary,
                "boundary_edges": defects["boundary_edge_count"],
                "nonmanifold_edges": defects["non_manifold_edge_count"],
                "self_intersections": defects["self_intersection_count"],
                "internal_faces": defects["suspected_internal_face_count"],
            },
        },
        "surface": {
            "watertight": defects["watertight"],
            "connected_components": defects["surface_connected_component_count"],
            "nonmanifold_edges": defects["non_manifold_edge_count"],
            "boundary_edges": defects["boundary_edge_count"],
        },
        "synthetic_controls_status": synthetic["status"],
        "backend_comparison": backend_rows,
        "recommended_next_action": recommendations,
        "figure_paths": [str(path) for path in figure_paths],
        "diagnostic_report": str(report_path),
        "runtime_s": time.perf_counter() - started,
    }
    write_json(root / "summary.json", summary)
    log_lines.extend(
        (
            f"port_root_causes={port_causes}",
            f"junction_root_causes={junction_causes}",
            f"runtime_s={summary['runtime_s']}",
            "diagnostic_status=COMPLETE",
        )
    )
    log_path.write_text("\n".join(log_lines) + "\n", encoding="utf-8")
    return summary


def run_geometry_diagnostics(
    rois: list[ROIRecord],
    config: CFDLumenConfig,
    run_layout: CFDRunLayout,
) -> list[dict[str, Any]]:
    """Run v2 diagnosis serially for traceable logs and emit a run-level report."""

    if not rois:
        raise ValueError("No ROIs were provided for geometry diagnosis")
    diagnostics_root = run_layout.run_root / "diagnostics"
    diagnostics_root.mkdir(parents=True, exist_ok=False)
    LOGGER.info("Running strict v2 diagnosis for %d ROI(s) with workers=1", len(rois))
    synthetic = run_synthetic_controls(config)
    write_json(diagnostics_root / "synthetic_controls.json", synthetic)
    summaries = [diagnose_roi(roi, config, run_layout, synthetic) for roi in rois]
    combined = [Path(summary["diagnostic_report"]).read_text(encoding="utf-8") for summary in summaries]
    run_report = diagnostics_root / "diagnostic_report.md"
    run_report.write_text("\n\n---\n\n".join(combined), encoding="utf-8")
    write_json(
        diagnostics_root / "diagnostic_manifest.json",
        {
            "status": "COMPLETE",
            "workers": 1,
            "roi_ids": [roi.roi_id for roi in rois],
            "synthetic_controls_status": synthetic["status"],
            "reports": [summary["diagnostic_report"] for summary in summaries],
        },
    )
    return summaries
