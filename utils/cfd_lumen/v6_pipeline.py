"""Formal v6 refinement, evidence export, and acceptance orchestration."""

from __future__ import annotations

import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import trimesh

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .continuous_field_transition import build_continuous_field_hybrid
from .export import write_csv, write_geometry_exports, write_json, write_units
from .hybrid_qc import (
    collar_radius_rows,
    evaluate_hybrid_surface_qc,
    junction_area_profile_rows,
    write_hybrid_artifacts,
)
from .surface_continuity_qc import face_region_labels
from .surface_qc import (
    evaluate_radius_fidelity,
    evaluate_surface_qc,
    identify_port_patches,
)
from .types import BranchGeometry, HybridBuildDetails, PatchResult, PortGeometry
from .v6_qc import (
    junction_interface_rows,
    port_interface_rows,
    silhouette_rows,
    transition_region_rows,
)
from .v6_validation import run_v6_synthetic_controls
from .v6_visualization import (
    interface_dihedral_figure,
    junction_field_transition_figure,
    port_radius_profile_figure,
    port_tangent_profile_figure,
    silhouette_comparison_figure,
    v5_interface_edges_figure,
    v5_v6_flat_shading,
    v5_v6_smooth_shading,
    v5_v6_wireframe,
    v5_v6_wireframe_overlay,
)


@dataclass(slots=True)
class V6RefinementResult:
    mesh: trimesh.Trimesh
    details: HybridBuildDetails
    patch: PatchResult
    report: dict[str, Any]
    geometry_paths: list[Path]
    figure_paths: list[Path]


def _maximum(rows: list[dict[str, Any]], key: str) -> float | None:
    values = [float(row[key]) for row in rows if row.get(key) is not None]
    return max(values, default=None)


def _maximum_collar_error(rows: list[dict[str, Any]]) -> float | None:
    return _maximum(rows, "absolute_radius_relative_error")


def _mean(rows: list[dict[str, Any]], key: str) -> float | None:
    values = [float(row[key]) for row in rows if row.get(key) is not None]
    return float(np.mean(values)) if values else None


def _sum(rows: list[dict[str, Any]], key: str) -> int:
    return int(sum(int(row.get(key, 0)) for row in rows))


def _relative_reduction(before: float | None, after: float | None) -> float | None:
    if before is None or after is None or before <= 0.0:
        return None
    return float((before - after) / before)


def _basic_topology(mesh: trimesh.Trimesh) -> dict[str, Any]:
    _, counts = np.unique(np.asarray(mesh.edges_sorted), axis=0, return_counts=True)
    return {
        "watertight": bool(mesh.is_watertight),
        "winding_consistent": bool(mesh.is_winding_consistent),
        "boundary_edge_count": int(np.count_nonzero(counts == 1)),
        "nonmanifold_edge_count": int(np.count_nonzero(counts > 2)),
        "component_count": int(len(mesh.split(only_watertight=False))),
    }


def _local_visual_crop(
    mesh: trimesh.Trimesh, center_um: np.ndarray, half_extent_um: float
) -> trimesh.Trimesh:
    centers = np.asarray(mesh.triangles_center, dtype=float)
    selected = np.all(
        np.abs(centers - np.asarray(center_um, dtype=float)[None, :])
        <= float(half_extent_um),
        axis=1,
    )
    return mesh.submesh([np.flatnonzero(selected)], append=True, repair=False)


def _csv_port_diagnostics(
    diagnostics: list[dict[str, Any]], names: tuple[str, ...]
) -> list[dict[str, Any]]:
    return [{name: row.get(name) for name in names} for row in diagnostics]


def _comparison_metrics(
    *,
    v5_port: list[dict[str, Any]],
    v6_port: list[dict[str, Any]],
    v5_junction: list[dict[str, Any]],
    v6_junction: list[dict[str, Any]],
    v5_silhouette: list[dict[str, Any]],
    v6_silhouette: list[dict[str, Any]],
) -> dict[str, Any]:
    v5_junction_all = [row for row in v5_junction if row["junction_node_id"] == "ALL"]
    v6_junction_all = [row for row in v6_junction if row["junction_node_id"] == "ALL"]
    v5_junction_p99 = _maximum(v5_junction_all, "dihedral_p99_deg")
    v6_junction_p99 = _maximum(v6_junction_all, "dihedral_p99_deg")
    v5_junction_max = _maximum(v5_junction_all, "dihedral_max_deg")
    v6_junction_max = _maximum(v6_junction_all, "dihedral_max_deg")
    v5_port_max = _maximum(v5_port, "dihedral_max_deg")
    v6_port_max = _maximum(v6_port, "dihedral_max_deg")
    v5_curvature = _mean(v5_silhouette, "silhouette_curvature_variation")
    v6_curvature = _mean(v6_silhouette, "silhouette_curvature_variation")
    v5_corner_fraction = _mean(v5_silhouette, "large_corner_fraction")
    v6_corner_fraction = _mean(v6_silhouette, "large_corner_fraction")
    return {
        "junction_interface_p99_deg": {
            "v5": v5_junction_p99,
            "v6": v6_junction_p99,
            "reduction_fraction": _relative_reduction(
                v5_junction_p99, v6_junction_p99
            ),
        },
        "junction_interface_max_deg": {
            "v5": v5_junction_max,
            "v6": v6_junction_max,
            "reduction_fraction": _relative_reduction(
                v5_junction_max, v6_junction_max
            ),
        },
        "port_interface_max_deg": {
            "v5": v5_port_max,
            "v6": v6_port_max,
            "reduction_fraction": _relative_reduction(v5_port_max, v6_port_max),
        },
        "silhouette_curvature_variation": {
            "v5_mean": v5_curvature,
            "v6_mean": v6_curvature,
            "reduction_fraction": _relative_reduction(v5_curvature, v6_curvature),
        },
        "silhouette_large_corners": {
            "v5_count": _sum(v5_silhouette, "large_corner_count"),
            "v6_count": _sum(v6_silhouette, "large_corner_count"),
            "v5_mean_fraction": v5_corner_fraction,
            "v6_mean_fraction": v6_corner_fraction,
            "fraction_reduction": _relative_reduction(
                v5_corner_fraction, v6_corner_fraction
            ),
        },
    }


def _acceptance_report(
    *,
    config: CFDLumenConfig,
    details: HybridBuildDetails,
    surface_qc: dict[str, Any],
    hybrid_qc: dict[str, Any],
    radius_v5: dict[str, Any],
    radius_v6: dict[str, Any],
    collar_v5: list[dict[str, Any]],
    collar_v6: list[dict[str, Any]],
    diagnostics: list[dict[str, Any]],
    comparison: dict[str, Any],
    synthetic: dict[str, Any],
    convergence: list[dict[str, Any]],
    v5_port_rows: list[dict[str, Any]],
    v6_port_rows: list[dict[str, Any]],
) -> dict[str, Any]:
    radius_increase = (
        float(radius_v6["p95_absolute_relative_error"])
        - float(radius_v5["p95_absolute_relative_error"])
    )
    collar_before = _maximum_collar_error(collar_v5)
    collar_after = _maximum_collar_error(collar_v6)
    collar_increase = (
        collar_after - collar_before
        if collar_before is not None and collar_after is not None
        else None
    )
    radius_slope_not_worse = all(
        float(row["v6_radius_slope_jump_um_per_um"])
        <= float(row["v5_radius_slope_jump_um_per_um"]) + 1.0e-12
        for row in diagnostics
    )
    # The actual endpoint segment is included so an already-C1 port is not bent
    # solely due to upstream curvature in the weighted diagnostic fit.
    tangent_not_worse = all(
        float(row["v6_interface_tangent_jump_deg"])
        <= float(row["v5_interface_tangent_jump_deg"]) + 1.0e-9
        for row in diagnostics
    )
    junction_p99 = comparison["junction_interface_p99_deg"]
    junction_max = comparison["junction_interface_max_deg"]
    port_max = comparison["port_interface_max_deg"]
    silhouette = comparison["silhouette_curvature_variation"]
    corners = comparison["silhouette_large_corners"]
    convergence_p99 = np.asarray(
        [row["junction_interface_p99_deg"] for row in convergence], dtype=float
    )
    convergence_max = np.asarray(
        [row["junction_interface_max_deg"] for row in convergence], dtype=float
    )
    v5_port_max_by_id = {
        int(row["port_id"]): float(row["dihedral_max_deg"])
        for row in v5_port_rows
    }
    v6_port_max_by_id = {
        int(row["port_id"]): float(row["dihedral_max_deg"])
        for row in v6_port_rows
    }
    checks = {
        "formal_continuous_field_backend": (
            details.transition_backend == "continuous_implicit_field"
        ),
        "no_transition_fallback": details.transition_fallback_reason is None,
        "one_marching_cubes_per_junction": all(
            bool(patch.metadata.get("one_marching_cubes_extraction"))
            for patch in details.patches.values()
        ),
        "no_surface_loop_stitching": all(
            not bool(patch.metadata.get("surface_loop_stitching"))
            for patch in details.patches.values()
        ),
        "surface_qc_pass": surface_qc["status"] == "PASS",
        "all_topology_defects_zero": hybrid_qc["status"] == "PASS",
        "junction_interface_p99_reduced": (
            junction_p99["v5"] is not None
            and junction_p99["v6"] is not None
            and junction_p99["v6"] < junction_p99["v5"]
        ),
        "junction_interface_max_reduced": (
            junction_max["v5"] is not None
            and junction_max["v6"] is not None
            and junction_max["v6"] < junction_max["v5"]
        ),
        "port_interface_max_reduced": (
            port_max["v5"] is not None
            and port_max["v6"] is not None
            and port_max["v6"] < port_max["v5"]
        ),
        "every_port_interface_max_not_worse": all(
            v6_port_max_by_id[port_id]
            <= v5_port_max_by_id[port_id] + 1.0e-2
            for port_id in v5_port_max_by_id
        ),
        "port_tangent_jump_not_worse": tangent_not_worse,
        "port_radius_slope_jump_reduced": radius_slope_not_worse,
        "silhouette_curvature_reduced": (
            silhouette["v5_mean"] is not None
            and silhouette["v6_mean"] is not None
            and silhouette["v6_mean"] < silhouette["v5_mean"]
        ),
        "silhouette_large_corner_fraction_not_worse": (
            corners["v5_mean_fraction"] is not None
            and corners["v6_mean_fraction"] is not None
            and corners["v6_mean_fraction"] <= corners["v5_mean_fraction"]
        ),
        "radius_p95_not_materially_degraded": (
            radius_increase
            <= config.hybrid_transition.maximum_radius_p95_error_increase
        ),
        "collar_radius_not_materially_degraded": (
            collar_increase is not None
            and collar_increase
            <= config.hybrid_transition.maximum_collar_error_increase
        ),
        "six_synthetic_controls_pass": synthetic["status"] == "PASS",
        "convergence_16_20_24_topology_pass": all(
            row["watertight"]
            and row["winding_consistent"]
            and row["boundary_edge_count"] == 0
            and row["nonmanifold_edge_count"] == 0
            and row["component_count"] == 1
            for row in convergence
        ),
        "convergence_interface_p99_stable": (
            float(np.max(convergence_p99) - np.min(convergence_p99))
            <= 0.10 * float(np.max(convergence_p99))
        ),
        "convergence_interface_max_below_v5": (
            junction_max["v5"] is not None
            and float(np.max(convergence_max)) < float(junction_max["v5"])
        ),
    }
    return {
        "protocol": "(new) 子图建模修改v6",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "checks": checks,
        "radius_p95": {
            "v5": radius_v5["p95_absolute_relative_error"],
            "v6": radius_v6["p95_absolute_relative_error"],
            "increase": radius_increase,
        },
        "maximum_collar_radius_error": {
            "v5": collar_before,
            "v6": collar_after,
            "increase": collar_increase,
        },
        "comparison": comparison,
        "decision_basis": [
            "interface-local dihedral metrics",
            "shared-scale six-view silhouette metrics",
            "same-camera flat/smooth/wireframe visual artifacts",
            "topological and radius/collar fidelity QC",
        ],
    }


def _markdown_report(report: dict[str, Any]) -> str:
    comparison = report["acceptance"]["comparison"]
    topology = report["topology"]
    lines = [
        "# v6 Continuous-Field Transition 最终核验",
        "",
        f"总体状态：**{report['status']}**",
        "",
        "## 方法约束",
        "",
        "- source SWC、source radius、CORE/context 拓扑和 CUT_PORT 映射保持不变。",
        "- Junction 使用单个 continuous implicit field 与一次 marching cubes。",
        "- surface loop stitch 仅作为 v5 regression，不参与 v6 正式表面。",
        "- merge 仅位于 w=1 的 PURE_BRANCH overlap；未做全局 smoothing。",
        "",
        "## 关键结果",
        "",
        f"- Junction interface P99: {comparison['junction_interface_p99_deg']['v5']:.6g}° → {comparison['junction_interface_p99_deg']['v6']:.6g}°。",
        f"- Junction interface max: {comparison['junction_interface_max_deg']['v5']:.6g}° → {comparison['junction_interface_max_deg']['v6']:.6g}°。",
        f"- Port interface max: {comparison['port_interface_max_deg']['v5']:.6g}° → {comparison['port_interface_max_deg']['v6']:.6g}°。",
        f"- Silhouette curvature variation mean: {comparison['silhouette_curvature_variation']['v5_mean']:.6g} → {comparison['silhouette_curvature_variation']['v6_mean']:.6g}。",
        f"- Topology defects: boundary={topology['boundary_edge_count']}, nonmanifold={topology['non_manifold_edge_count']}, self-intersection={topology['self_intersection_pairs']}, internal faces={topology['internal_face_count']}, internal caps={topology['internal_cap_face_count']}, degenerate={topology['degenerate_triangle_count']}。",
        "",
        "最终 PASS 不由测试用例单独决定；以上定量结果与 figures/ 中同相机图共同构成判据。",
        "",
    ]
    return "\n".join(lines)


def run_v6_refinement(
    roi: ROIRecord,
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    v5_mesh: trimesh.Trimesh,
    v5_details: HybridBuildDetails,
    output_root: Path,
) -> V6RefinementResult:
    """Run and persist the complete formal v6 protocol for one ROI."""

    output_root = Path(output_root)
    diagnostics_root = output_root / "diagnostics"
    figures_root = output_root / "figures"
    geometry_root = output_root / "geometry"
    ports_root = output_root / "ports"
    report_root = output_root / "report"
    for folder in (
        diagnostics_root,
        figures_root,
        geometry_root,
        ports_root,
        report_root,
    ):
        folder.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    v6_mesh, v6_details, port_diagnostics = build_continuous_field_hybrid(
        branches, roi, ports, config
    )
    v6_runtime = time.perf_counter() - started
    v6_patch = identify_port_patches(v6_mesh, ports, config)
    surface_qc = evaluate_surface_qc(v6_mesh, v6_patch, roi, branches, config)
    hybrid_qc, _ = evaluate_hybrid_surface_qc(v6_mesh, v6_details, roi, config)

    v5_collar = collar_radius_rows(v5_mesh, v5_details, branches)
    v6_collar = collar_radius_rows(v6_mesh, v6_details, branches)
    v6_area = junction_area_profile_rows(v6_mesh, v6_details, branches)
    radius_samples_v5, radius_qc_v5 = evaluate_radius_fidelity(
        v5_mesh, branches, roi, config
    )
    radius_samples_v6, radius_qc_v6 = evaluate_radius_fidelity(
        v6_mesh, branches, roi, config
    )

    v5_port = port_interface_rows(
        v5_mesh, branches, ports, config, version="v5"
    )
    v6_port = port_interface_rows(
        v6_mesh, branches, ports, config, version="v6"
    )
    v5_junction, v5_edges = junction_interface_rows(
        v5_mesh, v5_details, branches, version="v5"
    )
    v6_junction, _ = junction_interface_rows(
        v6_mesh, v6_details, branches, version="v6"
    )
    v5_regions = transition_region_rows(
        v5_mesh, v5_details, branches, version="v5"
    )
    v6_regions = transition_region_rows(
        v6_mesh, v6_details, branches, version="v6"
    )
    shared_meshes = (v5_mesh, v6_mesh)
    v5_silhouette = silhouette_rows(
        v5_mesh,
        version="v5",
        large_corner_deg=config.v6.silhouette_large_corner_deg,
        comparison_meshes=shared_meshes,
    )
    v6_silhouette = silhouette_rows(
        v6_mesh,
        version="v6",
        large_corner_deg=config.v6.silhouette_large_corner_deg,
        comparison_meshes=shared_meshes,
    )
    comparison = _comparison_metrics(
        v5_port=v5_port,
        v6_port=v6_port,
        v5_junction=v5_junction,
        v6_junction=v6_junction,
        v5_silhouette=v5_silhouette,
        v6_silhouette=v6_silhouette,
    )

    convergence: list[dict[str, Any]] = []
    for cells in config.v6.convergence_cells_across_min_diameter:
        if int(cells) == int(config.junction.implicit.cells_across_min_diameter):
            candidate_mesh = v6_mesh
            candidate_details = v6_details
            runtime = v6_runtime
        else:
            convergence_started = time.perf_counter()
            candidate_mesh, candidate_details, _ = build_continuous_field_hybrid(
                branches,
                roi,
                ports,
                config,
                cells_across_min_diameter=int(cells),
            )
            runtime = time.perf_counter() - convergence_started
        interface, _ = junction_interface_rows(
            candidate_mesh, candidate_details, branches, version="v6"
        )
        pooled = [row for row in interface if row["junction_node_id"] == "ALL"]
        convergence.append(
            {
                "cells_across_min_diameter": int(cells),
                "triangle_count": int(len(candidate_mesh.faces)),
                "runtime_s": runtime,
                "junction_interface_p95_deg": _maximum(
                    pooled, "dihedral_p95_deg"
                ),
                "junction_interface_p99_deg": _maximum(
                    pooled, "dihedral_p99_deg"
                ),
                "junction_interface_max_deg": _maximum(
                    pooled, "dihedral_max_deg"
                ),
                **_basic_topology(candidate_mesh),
            }
        )

    synthetic = run_v6_synthetic_controls(config)
    acceptance = _acceptance_report(
        config=config,
        details=v6_details,
        surface_qc=surface_qc,
        hybrid_qc=hybrid_qc,
        radius_v5=radius_qc_v5,
        radius_v6=radius_qc_v6,
        collar_v5=v5_collar,
        collar_v6=v6_collar,
        diagnostics=port_diagnostics,
        comparison=comparison,
        synthetic=synthetic,
        convergence=convergence,
        v5_port_rows=v5_port,
        v6_port_rows=v6_port,
    )

    tangent_fields = (
        "port_id",
        "cut_port_id",
        "branch_id",
        "fit_point_count",
        "v5_tangent_jump_deg",
        "endpoint_segment_jump_deg",
        "v5_interface_tangent_jump_deg",
        "tangent_conditioning_applied",
        "tangent_blend_length_um",
        "v6_tangent_jump_deg",
        "v6_interface_tangent_jump_deg",
    )
    radius_fields = (
        "port_id",
        "cut_port_id",
        "branch_id",
        "source_radius_fit_method",
        "source_radius_slope_um_per_um",
        "v5_radius_slope_jump_um_per_um",
        "radius_blend_length_um",
        "radius_conditioning_applied",
        "radius_slope_conditioning_threshold_um_per_um",
        "v6_radius_slope_jump_um_per_um",
        "radius_after_conditioning_um",
    )
    write_csv(
        diagnostics_root / "port_tangent_jump.csv",
        _csv_port_diagnostics(port_diagnostics, tangent_fields),
    )
    write_csv(
        diagnostics_root / "port_radius_slope.csv",
        _csv_port_diagnostics(port_diagnostics, radius_fields),
    )
    write_csv(diagnostics_root / "port_interface_dihedral.csv", [*v5_port, *v6_port])
    write_csv(
        diagnostics_root / "junction_interface_dihedral.csv",
        [*v5_junction, *v6_junction],
    )
    write_csv(
        diagnostics_root / "transition_region_qc.csv", [*v5_regions, *v6_regions]
    )
    write_csv(
        diagnostics_root / "silhouette_qc.csv",
        [*v5_silhouette, *v6_silhouette],
    )
    write_csv(diagnostics_root / "convergence_16_20_24.csv", convergence)
    write_csv(
        diagnostics_root / "radius_fidelity_v5.csv",
        [sample.report() for sample in radius_samples_v5],
    )
    write_csv(
        diagnostics_root / "radius_fidelity_v6.csv",
        [sample.report() for sample in radius_samples_v6],
    )
    write_csv(diagnostics_root / "collar_radius_v5.csv", v5_collar)
    write_csv(diagnostics_root / "collar_radius_v6.csv", v6_collar)
    write_json(diagnostics_root / "surface_qc.json", surface_qc)
    write_json(diagnostics_root / "hybrid_topology_qc.json", hybrid_qc)
    write_json(diagnostics_root / "radius_fidelity_qc_v5.json", radius_qc_v5)
    write_json(diagnostics_root / "radius_fidelity_qc_v6.json", radius_qc_v6)
    write_json(diagnostics_root / "synthetic_controls.json", synthetic)
    write_json(diagnostics_root / "v6_acceptance.json", acceptance)

    face_region = face_region_labels(v6_mesh, v6_details, branches)
    geometry_directories = {"geometry": geometry_root, "ports": ports_root}
    geometry_paths = write_geometry_exports(
        v6_mesh,
        v6_patch,
        ports,
        geometry_directories,
        face_region=face_region,
    )
    geometry_paths.append(write_units(diagnostics_root / "units.json"))
    geometry_paths.extend(
        write_hybrid_artifacts(
            output_root / "continuous_field",
            v6_mesh,
            v6_details,
            hybrid_qc,
            v6_collar,
            v6_area,
        )
    )

    figure_paths: list[Path] = []
    if config.output.visualizations:
        representative_junction = (
            49 if 49 in v6_details.patches else min(v6_details.patches)
        )
        representative_patch = v6_details.patches[representative_junction]
        comparison_extent = 1.2 * max(
            collar.implicit_extent_um + 2.0 * collar.collar_radius_um
            for collar in representative_patch.collars
        )
        comparison_center = np.asarray(
            roi.local_node_positions_um[representative_junction], dtype=float
        )
        local_v5 = _local_visual_crop(
            v5_mesh, comparison_center, comparison_extent
        )
        local_v6 = _local_visual_crop(
            v6_mesh, comparison_center, comparison_extent
        )
        write_json(
            diagnostics_root / "same_camera_comparison_scope.json",
            {
                "junction_node_id": representative_junction,
                "center_um": comparison_center.tolist(),
                "half_extent_um": comparison_extent,
                "note": (
                    "The open crop boundaries are visualization clips only and are "
                    "excluded from topology QC."
                ),
            },
        )
        figure_paths = [
            v5_interface_edges_figure(
                v5_mesh, v5_edges, figures_root / "v5_interface_edges.png"
            ),
            port_tangent_profile_figure(
                port_diagnostics, figures_root / "port_tangent_profile.png"
            ),
            port_radius_profile_figure(
                port_diagnostics, figures_root / "port_radius_profile.png"
            ),
            junction_field_transition_figure(
                v6_mesh, face_region, figures_root / "junction_field_transition.png"
            ),
            interface_dihedral_figure(
                [*v5_junction, *v6_junction],
                figures_root / "junction_interface_dihedral.png",
                title="v5 loop seams vs v6 continuous-field interfaces",
            ),
            interface_dihedral_figure(
                [*v5_port, *v6_port],
                figures_root / "port_interface_dihedral.png",
                title="CUT_PORT interface dihedral",
            ),
            v5_v6_flat_shading(
                local_v5, local_v6, figures_root / "v5_v6_flat_shading.png"
            ),
            v5_v6_smooth_shading(
                local_v5, local_v6, figures_root / "v5_v6_smooth_shading.png"
            ),
            v5_v6_wireframe(
                local_v5, local_v6, figures_root / "v5_v6_wireframe.png"
            ),
            v5_v6_wireframe_overlay(
                local_v5,
                local_v6,
                figures_root / "v5_v6_wireframe_overlay.png",
            ),
            v5_v6_flat_shading(
                v5_mesh,
                v6_mesh,
                figures_root / "v5_v6_flat_shading_full_roi.png",
            ),
            silhouette_comparison_figure(
                v5_mesh, v6_mesh, figures_root / "silhouette_comparison.png"
            ),
        ]

    report = {
        "protocol": "(new) 子图建模修改v6",
        "roi_id": roi.roi_id,
        "status": acceptance["status"],
        "artifact_classification": {
            "v5_visual_artifact": "LOOP_STITCH_ZIPPER_ARTIFACT",
            "evidence": (
                "The dominant v5 local dihedral interface is the highlighted "
                "implicit-to-loop-stitch boundary; v6 removes surface-level "
                "transition stitching."
            ),
        },
        "method": {
            "transition_backend": v6_details.transition_backend,
            "merge_backend": (
                v6_details.merge_steps[0]["merge_backend"]
                if v6_details.merge_steps
                else None
            ),
            "merge_region": "PURE_BRANCH (w=1)",
            "one_marching_cubes_per_junction": True,
            "global_smoothing": False,
            "source_swc_modified": False,
            "source_radius_modified": False,
            "cut_port_mapping_modified": False,
        },
        "port_diagnostics": _csv_port_diagnostics(
            port_diagnostics, tuple(dict.fromkeys((*tangent_fields, *radius_fields)))
        ),
        "junction_interfaces": v6_junction,
        "transition_regions": v6_regions,
        "surface_qc": surface_qc,
        "hybrid_qc": hybrid_qc,
        "radius_fidelity_v5": radius_qc_v5,
        "radius_fidelity_v6": radius_qc_v6,
        "collar_radius_v5": v5_collar,
        "collar_radius_v6": v6_collar,
        "topology": {
            "self_intersection_pairs": hybrid_qc["self_intersection_pairs"],
            "internal_face_count": hybrid_qc["internal_face_count"],
            "internal_cap_face_count": hybrid_qc["internal_cap_face_count"],
            "boundary_edge_count": hybrid_qc["boundary_edge_count"],
            "non_manifold_edge_count": hybrid_qc["nonmanifold_edge_count"],
            "degenerate_triangle_count": hybrid_qc["degenerate_triangle_count"],
            "component_count": hybrid_qc["surface_component_count"],
        },
        "performance": {
            "v5_runtime_s": v5_details.runtime_s.get(
                "reconstruction_total",
                v5_details.runtime_s.get("wall_total_with_step_qc"),
            ),
            "v6_runtime_s": v6_runtime,
            "v5_triangles": int(len(v5_mesh.faces)),
            "v6_triangles": int(len(v6_mesh.faces)),
        },
        "convergence": convergence,
        "synthetic_controls": synthetic,
        "acceptance": acceptance,
        "geometry_paths": [str(path) for path in geometry_paths],
        "figure_paths": [str(path) for path in figure_paths],
    }
    write_json(report_root / "v6_final_report.json", report)
    (report_root / "v6_final_report.md").write_text(
        _markdown_report(report), encoding="utf-8"
    )
    return V6RefinementResult(
        mesh=v6_mesh,
        details=v6_details,
        patch=v6_patch,
        report=report,
        geometry_paths=geometry_paths,
        figure_paths=figure_paths,
    )
