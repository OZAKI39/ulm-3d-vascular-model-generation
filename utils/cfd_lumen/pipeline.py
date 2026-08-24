"""End-to-end, failure-explicit CFD lumen workflow for one or many saved ROIs."""

from __future__ import annotations

import logging
import os
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import numpy as np

from utils.sampling.sampling_types import ROIRecord

from .branch_local_qc import (
    evaluate_branch_local_cross_section_qc,
    write_branch_local_qc_artifacts,
)
from .collision_qc import detect_nonadjacent_collisions
from .config import CFDLumenConfig
from .context_domain import build_cfd_context_domain
from .export import (
    create_roi_layout,
    verify_volume_mesh,
    write_centerlines,
    write_constructed_centerlines,
    write_csv,
    write_geometry_exports,
    write_json,
    write_source_metadata,
    write_units,
)
from .geometry_preprocess import resample_branches, validate_and_extract_branches
from .hybrid_qc import (
    collar_radius_rows,
    evaluate_hybrid_surface_qc,
    junction_area_profile_rows,
    max_junction_area_ratio,
    write_hybrid_artifacts,
)
from .hybrid_validation import run_v5_synthetic_controls
from .lumen_builder import build_lumen_surface
from .port_geometry import attach_core_port_provenance, construct_port_geometry
from .surface_continuity_qc import evaluate_surface_continuity, face_region_labels
from .surface_qc import evaluate_radius_fidelity, evaluate_surface_qc, identify_port_patches
from .types import CFDRunLayout, GeometryValidationError, ROIProcessResult
from .v6_pipeline import V6RefinementResult, run_v6_refinement
from .v7_pipeline import V7RefinementResult, run_v7_refinement
from .v8_pipeline import V8RefinementResult, run_v8_refinement
from .v9_pipeline import V9RefinementResult, run_v9_refinement
from .visualization import (
    collision_figure,
    cross_section_figure,
    lumen_overlay_figure,
    junction_transition_normals_figure,
    junction_transition_wireframe_figure,
    junction_v4_v5_figure,
    port_transition_profile_figure,
    port_v4_v5_figure,
    ports_figure,
    radius_fidelity_figure,
    resampling_figure,
    run_summary_figure,
    source_geometry_figure,
    wireframe_figure,
)


LOGGER = logging.getLogger("cfd_lumen")


def _maximum_collar_error(rows: list[dict[str, Any]]) -> float | None:
    return max(
        (
            float(row["absolute_radius_relative_error"])
            for row in rows
            if row.get("absolute_radius_relative_error") is not None
        ),
        default=None,
    )


def _reduction(before: float | None, after: float | None) -> float | None:
    if before is None or after is None or before <= 0.0:
        return None
    return float((before - after) / before)


def _v5_comparison_report(
    config: CFDLumenConfig,
    *,
    v4_build: Any,
    v5_build: Any,
    v4_surface_qc: dict[str, Any],
    v5_surface_qc: dict[str, Any],
    v4_hybrid_qc: dict[str, Any],
    v5_hybrid_qc: dict[str, Any],
    v4_radius_qc: dict[str, Any],
    v5_radius_qc: dict[str, Any],
    v4_collar_rows: list[dict[str, Any]],
    v5_collar_rows: list[dict[str, Any]],
    v4_continuity: dict[str, Any],
    v5_continuity: dict[str, Any],
) -> dict[str, Any]:
    transition = config.hybrid_transition
    v4_port_area = v4_continuity["port"][
        "maximum_cut_area_absolute_relative_error"
    ]
    v5_port_area = v5_continuity["port"][
        "maximum_cut_area_absolute_relative_error"
    ]
    v4_port_p95 = v4_continuity["port"]["port_normal_jump"][
        "normal_jump_p95_deg"
    ]
    v5_port_p95 = v5_continuity["port"]["port_normal_jump"][
        "normal_jump_p95_deg"
    ]
    v4_port_p99 = v4_continuity["port"]["port_normal_jump"][
        "normal_jump_p99_deg"
    ]
    v5_port_p99 = v5_continuity["port"]["port_normal_jump"][
        "normal_jump_p99_deg"
    ]
    v4_normal = v4_continuity["normal_jump"]["TRANSITION_COLLAR"]
    v5_normal = v5_continuity["normal_jump"]["TRANSITION_COLLAR"]
    v4_roughness = v4_continuity["roughness"]["regions"]["TRANSITION_COLLAR"]
    v5_roughness = v5_continuity["roughness"]["regions"]["TRANSITION_COLLAR"]
    v4_collar = _maximum_collar_error(v4_collar_rows)
    v5_collar = _maximum_collar_error(v5_collar_rows)
    junction_reduction = _reduction(
        v4_normal["normal_jump_p99_deg"], v5_normal["normal_jump_p99_deg"]
    )
    port_reduction = _reduction(v4_port_p99, v5_port_p99)
    roughness_reduction = _reduction(
        v4_roughness["laplacian_roughness_mean"],
        v5_roughness["laplacian_roughness_mean"],
    )
    topology_checks = dict(v5_hybrid_qc["checks"])
    checks = {
        "default_loop_stitch_used": v5_build.backend_used == "hybrid_loop_stitch",
        "no_transition_fallback": (
            v5_build.hybrid_details is not None
            and v5_build.hybrid_details.transition_fallback_reason is None
        ),
        "no_separate_port_cylinder": (
            v5_build.extension_primitive_count == 0
            and v5_continuity["port"]["separate_cylinder_primitive_count"] == 0
        ),
        "topology_all_pass": all(topology_checks.values()),
        "radius_p95_not_degraded": (
            float(v5_radius_qc["p95_absolute_relative_error"])
            <= float(v4_radius_qc["p95_absolute_relative_error"])
            + transition.maximum_radius_p95_error_increase
        ),
        "collar_radius_not_degraded": (
            v4_collar is not None
            and v5_collar is not None
            and v5_collar <= v4_collar + transition.maximum_collar_error_increase
        ),
        "port_area_not_degraded": (
            v4_port_area is not None
            and v5_port_area is not None
            and v5_port_area
            <= v4_port_area + transition.maximum_port_area_error_increase
        ),
        "port_normal_p99_reduced": (
            port_reduction is not None
            and port_reduction
            >= transition.minimum_port_normal_p99_reduction_fraction
        ),
        "junction_normal_p99_reduced": (
            junction_reduction is not None
            and junction_reduction
            >= transition.minimum_normal_p99_reduction_fraction
        ),
        "transition_roughness_reduced": (
            roughness_reduction is not None
            and roughness_reduction
            >= transition.minimum_transition_roughness_reduction_fraction
        ),
    }
    return {
        "protocol": "(new) 子图建模修改v5",
        "status": "PASS" if all(checks.values()) else "FAIL",
        "threshold_basis": (
            "Frozen after synthetic Y and representative real-ROI sweeps; tolerances "
            "are absolute relative-error increments."
        ),
        "thresholds": {
            "minimum_normal_p99_reduction_fraction": transition.minimum_normal_p99_reduction_fraction,
            "minimum_port_normal_p99_reduction_fraction": transition.minimum_port_normal_p99_reduction_fraction,
            "minimum_transition_roughness_reduction_fraction": transition.minimum_transition_roughness_reduction_fraction,
            "maximum_port_area_error_increase": transition.maximum_port_area_error_increase,
            "maximum_radius_p95_error_increase": transition.maximum_radius_p95_error_increase,
            "maximum_collar_error_increase": transition.maximum_collar_error_increase,
        },
        "checks": checks,
        "topology": {
            "v4": v4_hybrid_qc,
            "v5": v5_hybrid_qc,
        },
        "port": {
            "v4_area_error": v4_port_area,
            "v5_area_error": v5_port_area,
            "v4_normal_jump_p95_deg": v4_port_p95,
            "v5_normal_jump_p95_deg": v5_port_p95,
            "v4_normal_jump_p99_deg": v4_port_p99,
            "v5_normal_jump_p99_deg": v5_port_p99,
            "normal_p99_reduction_fraction": port_reduction,
            "v4_separate_cylinder_count": v4_continuity["port"][
                "separate_cylinder_primitive_count"
            ],
            "v5_separate_cylinder_count": v5_continuity["port"][
                "separate_cylinder_primitive_count"
            ],
        },
        "junction": {
            "v4_normal_jump_p95_deg": v4_normal["normal_jump_p95_deg"],
            "v5_normal_jump_p95_deg": v5_normal["normal_jump_p95_deg"],
            "v4_normal_jump_p99_deg": v4_normal["normal_jump_p99_deg"],
            "v5_normal_jump_p99_deg": v5_normal["normal_jump_p99_deg"],
            "normal_p99_reduction_fraction": junction_reduction,
            "transition_triangle_count": v5_continuity["transition_triangle_count"],
            "v4_transition_roughness_mean": v4_roughness[
                "laplacian_roughness_mean"
            ],
            "v5_transition_roughness_mean": v5_roughness[
                "laplacian_roughness_mean"
            ],
            "v4_transition_roughness_p95": v4_roughness[
                "laplacian_roughness_p95"
            ],
            "v5_transition_roughness_p95": v5_roughness[
                "laplacian_roughness_p95"
            ],
            "roughness_reduction_fraction": roughness_reduction,
            "v4_maximum_collar_radius_error": v4_collar,
            "v5_maximum_collar_radius_error": v5_collar,
            "transition_only_smoothing_enabled": transition.constrained_smoothing,
            "p2_decision": (
                "Transition-only boundary-loop resampling/remeshing is active. "
                "Constrained smoothing remains disabled because the parameter sweep "
                "did not improve absolute transition roughness."
            ),
        },
        "geometry": {
            "v4_radius_p95_error": v4_radius_qc["p95_absolute_relative_error"],
            "v5_radius_p95_error": v5_radius_qc["p95_absolute_relative_error"],
            "v4_surface_volume_um3": v4_surface_qc["enclosed_volume_um3"],
            "v5_surface_volume_um3": v5_surface_qc["enclosed_volume_um3"],
            "v4_triangle_count": v4_surface_qc["triangle_count"],
            "v5_triangle_count": v5_surface_qc["triangle_count"],
            "v4_runtime_s": v4_build.boolean_runtime_s,
            "v5_runtime_s": v5_build.boolean_runtime_s,
        },
    }


def _v4_acceptance_report(
    core_roi: ROIRecord,
    context_report: dict[str, Any],
    branch_local_qc: dict[str, Any] | None,
    controlled_comparison: dict[str, Any],
    hybrid_qc: dict[str, Any] | None,
    radius_qc: dict[str, Any],
    maximum_collar_error: float | None,
) -> dict[str, Any]:
    junction_summaries = {
        int(row["junction_node_id"]): row
        for row in (branch_local_qc or {}).get("junction_summaries", [])
    }
    fidelity_summaries = {
        int(row["junction_node_id"]): row
        for row in (branch_local_qc or {}).get("branch_fidelity_summaries", [])
    }
    j49 = junction_summaries.get(49)
    j49_fidelity = fidelity_summaries.get(49)
    extended = [
        row
        for row in context_report.get("port_mappings", [])
        if float(row.get("added_centerline_length_um", 0.0)) > 0.0
    ]
    topology = hybrid_qc or {}
    answers = [
        {
            "question": "1. J13 是否通过 CFD context extension 解决？",
            "answer": context_report.get("final_JUNCTION_PORT_REGION_CONFLICT_count") == 0,
            "evidence": context_report.get("conflict_history", []),
        },
        {
            "question": "2. 新 CFD boundary 比原 CORE boundary 向外增加多少 um？",
            "answer": [
                {
                    "original_cut_port_id": row["original_cut_port"]["cut_port_id"],
                    "added_centerline_length_um": row["added_centerline_length_um"],
                    "new_cfd_port_ids": [port["cut_port_id"] for port in row["new_cfd_ports"]],
                }
                for row in extended
            ],
        },
        {
            "question": "3. 是否增加了真实 global SWC branches？",
            "answer": bool(context_report.get("added_global_branch_ids")),
            "added_global_branch_ids": context_report.get("added_global_branch_ids", []),
            "source_manifest": context_report.get("source_manifest"),
        },
        {
            "question": "4. CORE ROI 是否保持完全不变？",
            "answer": bool(context_report.get("core_roi_unchanged")),
            "core_roi_id": core_roi.roi_id,
            "signature_before": context_report.get("core_roi_signature_before"),
            "signature_after": context_report.get("core_roi_signature_after"),
        },
        {
            "question": "5. J49 原 11.856 area ratio 中多少是多分支截面污染？",
            "answer": j49,
        },
        {
            "question": "6. J49 branch-local radius / area error 是多少？",
            "answer": j49_fidelity,
        },
        {
            "question": "7. J49 是否存在需进一步修改 geometry 的 artificial bulge？",
            "answer": bool(
                branch_local_qc
                and branch_local_qc.get("controlled_local_implicit_required")
            ),
            "decision_scope": "BRANCH_FIDELITY_REGION only",
        },
        {
            "question": "8. 若不需要，为什么 v3 hybrid 可接受？",
            "answer": {
                "branch_local_qc_status": (branch_local_qc or {}).get("status"),
                "topology_status": topology.get("status"),
                "junction_area_upper_bound_used_as_hard_failure": False,
                "junction_core_minimum_total_plane_area_ratio": (
                    (branch_local_qc or {}).get("junction_core_minimum_total_plane_area_ratio")
                ),
                "junction_core_normal_jump_p99_deg": (
                    (branch_local_qc or {}).get("junction_core_normal_jump_p99_deg")
                ),
            },
        },
        {
            "question": "9. 若需要，controlled local implicit 是否改善？",
            "answer": controlled_comparison,
        },
        {
            "question": "10. 当前完整 representative ROI 是否达到 CFD-ready PASS？",
            "answer": True,
            "radius_p95_absolute_relative_error": radius_qc.get(
                "p95_absolute_relative_error"
            ),
            "maximum_collar_radius_absolute_relative_error": maximum_collar_error,
        },
    ]
    return {
        "status": "PASS",
        "protocol": "(new) 子图建模修改v4",
        "core_roi_id": core_roi.roi_id,
        "answers": answers,
    }


def _identify_patches(build: Any, ports: list[Any], config: CFDLumenConfig) -> Any:
    tolerance = config.ports.plane_tolerance_um
    face_alignment = config.ports.minimum_normal_alignment
    if build.backend_used == "implicit" and build.implicit_grid is not None:
        tolerance = max(tolerance, 0.55 * float(build.implicit_grid["spacing_um"]))
        # Marching cubes creates a one-cell bevel around a mathematically flat
        # cap. Include its outward-facing triangles in the cap-area estimate;
        # the configured normal threshold is still applied to the weighted
        # patch normal by the quantitative QC below.
        face_alignment = 0.0
    return identify_port_patches(
        build.mesh,
        ports,
        config,
        plane_tolerance_um=tolerance,
        face_normal_alignment=face_alignment,
    )


def _convergence_rows(
    roi: ROIRecord,
    branches: list[Any],
    ports: list[Any],
    config: CFDLumenConfig,
    *,
    base_build: Any,
    base_patch: Any,
    base_surface_qc: dict[str, Any],
    base_radius_qc: dict[str, Any],
) -> list[dict[str, Any]]:
    if not config.convergence.enabled:
        return []
    rows: list[dict[str, Any]] = []
    for sides in config.convergence.tube_sides:
        started = time.perf_counter()
        try:
            if int(sides) == config.geometry.tube_sides and base_build.backend_used in {
                "manifold",
                "hybrid_local_implicit",
                "hybrid_controlled_local_implicit",
                "hybrid_loop_stitch",
                "hybrid_manifold_boolean_fallback",
            }:
                build = base_build
                patch = base_patch
                surface_qc = base_surface_qc
                radius_qc = base_radius_qc
            else:
                build = build_lumen_surface(
                    branches,
                    roi,
                    ports,
                    config,
                    tube_sides=int(sides),
                    backend="manifold",
                    controlled_local_implicit=(
                        base_build.backend_used == "hybrid_controlled_local_implicit"
                    ),
                )
                patch = _identify_patches(build, ports, config)
                surface_qc = evaluate_surface_qc(build.mesh, patch, roi, branches, config)
                _, radius_qc = evaluate_radius_fidelity(build.mesh, branches, roi, config)
            port_areas = [float(row["patch_area_um2"]) for row in patch.port_rows]
            rows.append(
                {
                    "tube_sides": int(sides),
                    "status": surface_qc["status"],
                    "backend_used": build.backend_used,
                    "triangle_count": surface_qc["triangle_count"],
                    "surface_area_um2": surface_qc["surface_area_um2"],
                    "enclosed_volume_um3": surface_qc["enclosed_volume_um3"],
                    "mean_port_area_um2": float(np.mean(port_areas)) if port_areas else None,
                    "radius_p95_absolute_relative_error": radius_qc["p95_absolute_relative_error"],
                    "runtime_s": time.perf_counter() - started,
                    "failure_reason": None,
                }
            )
        except Exception as exc:
            rows.append(
                {
                    "tube_sides": int(sides),
                    "status": "FAIL",
                    "backend_used": None,
                    "triangle_count": None,
                    "surface_area_um2": None,
                    "enclosed_volume_um3": None,
                    "mean_port_area_um2": None,
                    "radius_p95_absolute_relative_error": None,
                    "runtime_s": time.perf_counter() - started,
                    "failure_reason": f"{type(exc).__name__}: {exc}",
                }
            )
    return rows


def process_roi(
    roi: ROIRecord,
    config: CFDLumenConfig,
    run_layout: CFDRunLayout,
    sampling_run: Path,
) -> ROIProcessResult:
    """Process one ROI, persist every QC result, and never silently skip a failure."""

    started_total = time.perf_counter()
    core_roi = roi
    directories = create_roi_layout(run_layout, core_roi.roi_id)
    timings: dict[str, float] = {"load": 0.0}
    pre_qc: dict[str, Any] = {"roi_id": core_roi.roi_id, "status": "NOT_RUN"}
    try:
        stage = time.perf_counter()
        context = build_cfd_context_domain(core_roi, sampling_run, config)
        roi = context.cfd_roi
        timings["cfd_context_domain"] = time.perf_counter() - stage
        write_json(directories["qc"] / "cfd_context_domain.json", context.report)
        context_rows: list[dict[str, Any]] = []
        for mapping in context.port_mappings:
            original = mapping["original_cut_port"]
            for replacement in mapping["new_cfd_ports"]:
                context_rows.append(
                    {
                        "original_cut_port_id": original["cut_port_id"],
                        "original_boundary_role": original["boundary_role"],
                        "original_global_edge_id": original["global_edge_id"],
                        "original_x_um": original["intersection_x_um"],
                        "original_y_um": original["intersection_y_um"],
                        "original_z_um": original["intersection_z_um"],
                        "original_radius_um": original["radius_at_cut_um"],
                        "original_distance_to_nearest_junction_um": mapping[
                            "original_distance_to_nearest_junction_um"
                        ],
                        "new_cfd_port_id": replacement["cut_port_id"],
                        "new_boundary_role": replacement["boundary_role"],
                        "new_global_edge_id": replacement["global_edge_id"],
                        "new_x_um": replacement["intersection_x_um"],
                        "new_y_um": replacement["intersection_y_um"],
                        "new_z_um": replacement["intersection_z_um"],
                        "new_radius_um": replacement["radius_at_cut_um"],
                        "added_centerline_length_um": mapping["added_centerline_length_um"],
                        "added_global_edge_ids": ";".join(
                            map(str, mapping["added_global_edge_ids"])
                        ),
                        "added_branch_ids": ";".join(map(str, mapping["added_branch_ids"])),
                        "reason_for_extension": mapping["reason_for_extension"],
                    }
                )
        write_csv(directories["boundary"] / "core_to_cfd_port_mapping.csv", context_rows)
        write_json(
            directories["source"] / "domain_scope.json",
            {
                "CORE_ROI": {
                    "roi_id": core_roi.roi_id,
                    "node_count": core_roi.node_count,
                    "edge_count": core_roi.edge_count,
                    "cut_ports": [port.report() for port in core_roi.cut_ports],
                    "used_for": [
                        "sampling features",
                        "clustering",
                        "representative identity",
                        "default MB/ULM statistics",
                    ],
                },
                "CFD_DOMAIN": {
                    "roi_id": roi.roi_id,
                    "node_count": roi.node_count,
                    "edge_count": roi.edge_count,
                    "cut_ports": [port.report() for port in roi.cut_ports],
                    "used_for": ["CFD geometry", "CFD flow"],
                },
            },
        )

        stage = time.perf_counter()
        branches, pre_qc = validate_and_extract_branches(roi, config)
        pre_qc["core_roi_id"] = core_roi.roi_id
        pre_qc["cfd_domain_id"] = roi.roi_id
        pre_qc["core_roi_unchanged"] = context.report.get("core_roi_unchanged", True)
        pre_qc["status"] = "PASS"
        timings["pre_geometry_qc"] = time.perf_counter() - stage

        stage = time.perf_counter()
        branches = resample_branches(branches, config)
        ports = attach_core_port_provenance(
            construct_port_geometry(roi, config), context.port_mappings
        )
        timings["resample"] = time.perf_counter() - stage
        pre_qc["resampled_point_count"] = int(sum(len(branch.points_um) for branch in branches))
        pre_qc["centerline_smoothing"] = config.geometry.centerline_smoothing
        write_json(directories["qc"] / "pre_geometry_qc.json", pre_qc)
        write_source_metadata(roi, branches, directories, sampling_run=sampling_run)
        if config.output.save_centerlines:
            write_centerlines(branches, directories)

        stage = time.perf_counter()
        collision_events, collision_qc = detect_nonadjacent_collisions(branches, config)
        timings["collision_qc"] = time.perf_counter() - stage
        write_csv(
            directories["qc"] / "collision_report.csv",
            [event.report() for event in collision_events],
            fieldnames=(
                "branch_id_a", "branch_id_b", "segment_index_a", "segment_index_b",
                "distance_um", "radius_a_um", "radius_b_um", "clearance_um", "classification",
                "closest_a_x_um", "closest_a_y_um", "closest_a_z_um",
                "closest_b_x_um", "closest_b_y_um", "closest_b_z_um",
            ),
        )
        write_json(directories["qc"] / "collision_qc.json", collision_qc)
        collision_figure(branches, collision_events, directories["figures"] / "00_collision_qc.png")
        if collision_qc["hard_collision_count"]:
            raise GeometryValidationError(
                f"ROI {roi.roi_id} has {collision_qc['hard_collision_count']} non-adjacent HARD_COLLISION event(s)"
            )

        v4_build: Any | None = None
        if config.hybrid_transition.capture_v4_comparison:
            stage = time.perf_counter()
            v4_build = build_lumen_surface(
                branches,
                roi,
                ports,
                config,
                transition_backend="manifold_boolean",
                continuous_port_extensions=False,
            )
            timings["v4_reference_reconstruction"] = time.perf_counter() - stage
            if (
                v4_build.backend_used
                not in {"hybrid_local_implicit", "hybrid_controlled_local_implicit"}
                or v4_build.hybrid_details is None
            ):
                raise GeometryValidationError(
                    "V5_BASELINE_RECONSTRUCTION_FAILED: v4 manifold-Boolean reference "
                    f"used {v4_build.backend_used}"
                )

        stage = time.perf_counter()
        build = build_lumen_surface(branches, roi, ports, config)
        if build.fallback_reason:
            LOGGER.warning(
                "ROI %s reconstruction fallback: %s",
                roi.roi_id,
                build.fallback_reason,
            )
        timings["v5_surface_reconstruction"] = time.perf_counter() - stage
        timings["primitive_construction_and_boolean"] = timings[
            "v5_surface_reconstruction"
        ]

        stage = time.perf_counter()
        patch = _identify_patches(build, ports, config)
        surface_qc = evaluate_surface_qc(build.mesh, patch, roi, branches, config)
        hybrid_qc: dict[str, Any] | None = None
        hybrid_collar_rows: list[dict[str, Any]] = []
        hybrid_area_rows: list[dict[str, Any]] = []
        branch_local_qc: dict[str, Any] | None = None
        branch_local_rows: list[dict[str, Any]] = []
        branch_ownership_surface: Any | None = None
        controlled_comparison: dict[str, Any] = {
            "triggered": False,
            "reason": "No severe branch-local abnormality outside the junction core",
        }
        hybrid_backends = {
            "hybrid_local_implicit",
            "hybrid_controlled_local_implicit",
            "hybrid_loop_stitch",
            "hybrid_manifold_boolean_fallback",
        }
        if build.backend_used in hybrid_backends and build.hybrid_details is not None:
            hybrid_qc, _ = evaluate_hybrid_surface_qc(
                build.mesh, build.hybrid_details, roi, config
            )
            hybrid_collar_rows = collar_radius_rows(build.mesh, build.hybrid_details, branches)
            hybrid_area_rows = junction_area_profile_rows(
                build.mesh, build.hybrid_details, branches
            )
            if config.branch_local_qc.enabled:
                (
                    branch_ownership_surface,
                    branch_local_rows,
                    branch_local_qc,
                ) = evaluate_branch_local_cross_section_qc(
                    build.mesh, build.hybrid_details, branches, config
                )
        if (
            build.surface_continuity_version != "v5"
            and branch_local_qc
            and branch_local_qc["controlled_local_implicit_required"]
        ):
            v3_snapshot = {
                "backend": build.backend_used,
                "surface_qc": surface_qc,
                "hybrid_qc": hybrid_qc,
                "branch_local_qc": branch_local_qc,
            }
            candidate_root = (
                run_layout.run_root / "hybrid" / core_roi.roi_id / "v3_pre_control"
            )
            if branch_ownership_surface is not None:
                write_branch_local_qc_artifacts(
                    candidate_root / "qc",
                    branch_ownership_surface,
                    branch_local_rows,
                    branch_local_qc,
                )
            controlled_build = build_lumen_surface(
                branches,
                roi,
                ports,
                config,
                controlled_local_implicit=True,
            )
            controlled_patch = _identify_patches(controlled_build, ports, config)
            controlled_surface_qc = evaluate_surface_qc(
                controlled_build.mesh, controlled_patch, roi, branches, config
            )
            if controlled_build.hybrid_details is None:
                raise GeometryValidationError(
                    "CONTROLLED_LOCAL_IMPLICIT did not return hybrid build details"
                )
            controlled_hybrid_qc, _ = evaluate_hybrid_surface_qc(
                controlled_build.mesh, controlled_build.hybrid_details, roi, config
            )
            controlled_collar_rows = collar_radius_rows(
                controlled_build.mesh, controlled_build.hybrid_details, branches
            )
            controlled_area_rows = junction_area_profile_rows(
                controlled_build.mesh, controlled_build.hybrid_details, branches
            )
            (
                controlled_ownership_surface,
                controlled_branch_rows,
                controlled_branch_qc,
            ) = evaluate_branch_local_cross_section_qc(
                controlled_build.mesh,
                controlled_build.hybrid_details,
                branches,
                config,
            )
            controlled_comparison = {
                "triggered": True,
                "trigger_reason": (
                    "Severe branch-local area abnormality outside JUNCTION_CORE_REGION"
                ),
                "v3": v3_snapshot,
                "v4_controlled": {
                    "backend": controlled_build.backend_used,
                    "surface_qc": controlled_surface_qc,
                    "hybrid_qc": controlled_hybrid_qc,
                    "branch_local_qc": controlled_branch_qc,
                },
                "selected_backend": controlled_build.backend_used,
            }
            build = controlled_build
            patch = controlled_patch
            surface_qc = controlled_surface_qc
            hybrid_qc = controlled_hybrid_qc
            hybrid_collar_rows = controlled_collar_rows
            hybrid_area_rows = controlled_area_rows
            branch_ownership_surface = controlled_ownership_surface
            branch_local_rows = controlled_branch_rows
            branch_local_qc = controlled_branch_qc
        if build.backend_used in hybrid_backends and build.hybrid_details is not None:
            hybrid_root = run_layout.run_root / "hybrid" / core_roi.roi_id
            write_hybrid_artifacts(
                hybrid_root,
                build.mesh,
                build.hybrid_details,
                hybrid_qc,
                hybrid_collar_rows,
                hybrid_area_rows,
            )
            if branch_local_qc is not None and branch_ownership_surface is not None:
                write_branch_local_qc_artifacts(
                    hybrid_root / "qc",
                    branch_ownership_surface,
                    branch_local_rows,
                    branch_local_qc,
                )
            write_json(hybrid_root / "qc" / "controlled_implicit_comparison.json", controlled_comparison)
        if (
            (surface_qc["status"] == "FAIL" or (hybrid_qc and hybrid_qc["status"] == "FAIL"))
            and build.backend_used in {"manifold", *hybrid_backends}
            and config.boolean.allow_implicit_fallback
            and config.hybrid_transition.backend != "loop_stitch"
        ):
            explicit_failure = (
                "hybrid-specific QC failed"
                if hybrid_qc and hybrid_qc["status"] == "FAIL"
                else "surface QC failed after explicit/hybrid reconstruction"
            )
            build = build_lumen_surface(branches, roi, ports, config, backend="implicit")
            build.backend_requested = "manifold"
            build.fallback_reason = explicit_failure
            patch = _identify_patches(build, ports, config)
            surface_qc = evaluate_surface_qc(build.mesh, patch, roi, branches, config)
            branch_local_qc = None
        timings["patch_detection_and_surface_qc"] = time.perf_counter() - stage
        write_json(directories["qc"] / "surface_qc.json", surface_qc)
        write_csv(directories["boundary"] / "ports.csv", patch.port_rows)
        write_json(directories["boundary"] / "ports.json", patch.port_rows)

        if (
            build.backend_used != "hybrid_loop_stitch"
            or build.hybrid_details is None
            or build.hybrid_details.transition_fallback_reason is not None
        ):
            raise GeometryValidationError(
                "V5_LOOP_STITCH_REQUIRED: formal v5 geometry may not use manifold-Boolean "
                "or global-implicit fallback"
            )
        if v4_build is None or v4_build.hybrid_details is None:
            raise GeometryValidationError("V5_BASELINE_RECONSTRUCTION_FAILED")
        if hybrid_qc is None:
            raise GeometryValidationError("V5_HYBRID_QC_NOT_RUN")

        stage = time.perf_counter()
        v5_port_continuity_rows, v5_continuity = evaluate_surface_continuity(
            build.mesh, build.hybrid_details, branches, ports, config
        )
        v4_patch = _identify_patches(v4_build, ports, config)
        v4_surface_qc = evaluate_surface_qc(
            v4_build.mesh, v4_patch, roi, branches, config
        )
        v4_hybrid_qc, _ = evaluate_hybrid_surface_qc(
            v4_build.mesh, v4_build.hybrid_details, roi, config
        )
        v4_collar_rows = collar_radius_rows(
            v4_build.mesh, v4_build.hybrid_details, branches
        )
        v4_port_continuity_rows, v4_continuity = evaluate_surface_continuity(
            v4_build.mesh, v4_build.hybrid_details, branches, ports, config
        )
        timings["v4_v5_continuity_qc"] = time.perf_counter() - stage
        write_csv(
            directories["qc"] / "port_transition_profiles_v5.csv",
            v5_port_continuity_rows,
        )
        write_csv(
            directories["qc"] / "port_transition_profiles_v4.csv",
            v4_port_continuity_rows,
        )
        write_json(directories["qc"] / "surface_continuity_v5.json", v5_continuity)
        write_json(directories["qc"] / "surface_continuity_v4.json", v4_continuity)
        write_json(directories["qc"] / "surface_qc_v4_reference.json", v4_surface_qc)
        write_json(directories["qc"] / "hybrid_qc_v4_reference.json", v4_hybrid_qc)
        write_csv(directories["qc"] / "collar_radius_v4_reference.csv", v4_collar_rows)

        stage = time.perf_counter()
        radius_samples, radius_qc = evaluate_radius_fidelity(
            build.mesh, branches, roi, config
        )
        v4_radius_samples, v4_radius_qc = evaluate_radius_fidelity(
            v4_build.mesh, branches, roi, config
        )
        timings["radius_fidelity_v4_v5"] = time.perf_counter() - stage
        timings["radius_fidelity"] = timings["radius_fidelity_v4_v5"]
        write_csv(
            directories["qc"] / "radius_fidelity.csv",
            [sample.report() for sample in radius_samples],
        )
        write_json(directories["qc"] / "radius_fidelity_qc.json", radius_qc)
        write_csv(
            directories["qc"] / "radius_fidelity_v4_reference.csv",
            [sample.report() for sample in v4_radius_samples],
        )
        write_json(
            directories["qc"] / "radius_fidelity_qc_v4_reference.json",
            v4_radius_qc,
        )
        v5_comparison = _v5_comparison_report(
            config,
            v4_build=v4_build,
            v5_build=build,
            v4_surface_qc=v4_surface_qc,
            v5_surface_qc=surface_qc,
            v4_hybrid_qc=v4_hybrid_qc,
            v5_hybrid_qc=hybrid_qc,
            v4_radius_qc=v4_radius_qc,
            v5_radius_qc=radius_qc,
            v4_collar_rows=v4_collar_rows,
            v5_collar_rows=hybrid_collar_rows,
            v4_continuity=v4_continuity,
            v5_continuity=v5_continuity,
        )
        write_json(directories["qc"] / "v4_v5_comparison.json", v5_comparison)
        write_json(directories["qc"] / "v5_acceptance_report.json", v5_comparison)

        if surface_qc["status"] != "PASS":
            raise GeometryValidationError(
                f"ROI {roi.roi_id} final lumen failed surface QC: "
                + ", ".join(
                    name for name, passed in surface_qc["checks"].items() if not passed
                )
            )
        if hybrid_qc["status"] != "PASS":
            raise GeometryValidationError(
                f"ROI {roi.roi_id} final hybrid lumen failed hybrid-specific QC: "
                + ", ".join(
                    name for name, passed in hybrid_qc["checks"].items() if not passed
                )
            )
        if branch_local_qc is not None and branch_local_qc["status"] != "PASS":
            raise GeometryValidationError(
                f"ROI {roi.roi_id} failed BRANCH_LOCAL_CROSS_SECTION_QC: "
                + ", ".join(
                    name
                    for name, passed in branch_local_qc["checks"].items()
                    if not passed
                )
            )
        if radius_qc["status"] != "PASS":
            raise GeometryValidationError(
                f"ROI {roi.roi_id} radius-fidelity P95 exceeds configured threshold"
            )
        if v5_comparison["status"] != "PASS":
            raise GeometryValidationError(
                "ROI "
                f"{roi.roi_id} failed v5 continuity acceptance: "
                + ", ".join(
                    name
                    for name, passed in v5_comparison["checks"].items()
                    if not passed
                )
            )

        v6_result: V6RefinementResult | None = None
        if config.v6.enabled:
            stage = time.perf_counter()
            v6_result = run_v6_refinement(
                roi,
                branches,
                ports,
                config,
                build.mesh,
                build.hybrid_details,
                run_layout.run_root / "v6" / core_roi.roi_id,
            )
            timings["v6_refinement_and_acceptance"] = time.perf_counter() - stage
            if v6_result.report["status"] != "PASS":
                raise GeometryValidationError(
                    "ROI "
                    f"{roi.roi_id} failed v6 combined acceptance: "
                    + ", ".join(
                        name
                        for name, passed in v6_result.report["acceptance"][
                            "checks"
                        ].items()
                        if not passed
                    )
                )

        v7_result: V7RefinementResult | None = None
        if config.v7.enabled:
            if v6_result is None:
                raise GeometryValidationError(
                    "V7 requires the topology-valid v6 baseline for formal A/B selection"
                )
            stage = time.perf_counter()
            v7_result = run_v7_refinement(
                roi,
                branches,
                ports,
                config,
                v6_result,
                run_layout.run_root / "v7" / core_roi.roi_id,
            )
            timings["v7_unified_polyball_and_acceptance"] = (
                time.perf_counter() - stage
            )

        v8_result: V8RefinementResult | None = None
        if config.v8.enabled:
            if v6_result is None or v7_result is None:
                raise GeometryValidationError(
                    "V8 requires the topology-valid v6 baseline and v7 hard-min build"
                )
            stage = time.perf_counter()
            v8_result = run_v8_refinement(
                roi,
                branches,
                ports,
                config,
                v6_result,
                v7_result,
                run_layout.run_root / "v8" / core_roi.roi_id,
            )
            timings["v8_localization_and_smooth_union_sensitivity"] = (
                time.perf_counter() - stage
            )

        v9_result: V9RefinementResult | None = None
        if config.v9.enabled:
            if v6_result is None or v7_result is None or v8_result is None:
                raise GeometryValidationError(
                    "V9 requires v6 topology details, v7 hard-min, and v8 branch-union controls"
                )
            stage = time.perf_counter()
            v9_result = run_v9_refinement(
                roi,
                branches,
                ports,
                config,
                v6_result,
                v7_result,
                v8_result,
                run_layout.run_root / "v9" / core_roi.roi_id,
            )
            timings["v9_segment_crease_spline_and_competition_union"] = (
                time.perf_counter() - stage
            )

        stage = time.perf_counter()
        final_mesh = (
            v9_result.mesh
            if v9_result is not None
            else v8_result.mesh
            if v8_result is not None
            else v7_result.mesh
            if v7_result is not None
            else v6_result.mesh
            if v6_result is not None
            else build.mesh
        )
        final_details = (
            v6_result.details if v6_result is not None else build.hybrid_details
        )
        final_patch = (
            v9_result.patch
            if v9_result is not None
            else v8_result.patch
            if v8_result is not None
            else v7_result.patch
            if v7_result is not None
            else v6_result.patch
            if v6_result is not None
            else patch
        )
        final_face_region = (
            np.zeros(len(final_mesh.faces), dtype=np.uint8)
            if v9_result is not None
            or v8_result is not None
            or v7_result is not None
            and v7_result.decision == "ADOPT_V7_UNIFIED_POLYBALL"
            else face_region_labels(final_mesh, final_details, branches)
        )
        geometry_paths = write_geometry_exports(
            final_mesh,
            final_patch,
            ports,
            directories,
            face_region=final_face_region,
        )
        geometry_paths.extend(
            write_constructed_centerlines(final_details, directories)
        )
        write_units(directories["qc"] / "units.json")
        timings["export"] = time.perf_counter() - stage

        stage = time.perf_counter()
        convergence = _convergence_rows(
            roi,
            branches,
            ports,
            config,
            base_build=build,
            base_patch=patch,
            base_surface_qc=surface_qc,
            base_radius_qc=radius_qc,
        )
        write_csv(directories["qc"] / "geometry_convergence.csv", convergence)
        timings["convergence"] = time.perf_counter() - stage

        surface_stl = directories["geometry"] / "lumen_surface_m.stl"
        volume_mesh = verify_volume_mesh(
            surface_stl,
            directories["mesh"] / "lumen_volume.msh",
            minimum_radius_um=float(roi.local_node_radius_um.min()),
            config=config,
        )
        write_json(directories["mesh"] / "volume_mesh_qc.json", volume_mesh)

        stage = time.perf_counter()
        figure_paths: list[Path] = [directories["figures"] / "00_collision_qc.png"]
        if config.output.visualizations:
            figure_paths.extend(
                (
                    source_geometry_figure(roi, branches, directories["figures"] / "01_source_centerline.png"),
                    resampling_figure(branches, directories["figures"] / "02_resampled_centerline.png"),
                    lumen_overlay_figure(build.mesh, roi, branches, ports, directories["figures"] / "03_lumen_overlay.png"),
                    ports_figure(build.mesh, patch, ports, directories["figures"] / "04_ports_and_normals.png"),
                    wireframe_figure(build.mesh, directories["figures"] / "05_surface_wireframe.png"),
                    radius_fidelity_figure(radius_samples, directories["figures"] / "06_radius_fidelity.png"),
                    cross_section_figure(radius_samples, directories["figures"] / "07_cross_sections.png"),
                )
            )
            figure_paths.extend(
                (
                    port_transition_profile_figure(
                        v4_port_continuity_rows,
                        v5_port_continuity_rows,
                        directories["figures"] / "port_transition_area_profile.png",
                        metric="area",
                    ),
                    port_transition_profile_figure(
                        v4_port_continuity_rows,
                        v5_port_continuity_rows,
                        directories["figures"] / "port_normal_profile.png",
                        metric="normal",
                    ),
                    port_v4_v5_figure(
                        v4_build.mesh,
                        build.mesh,
                        ports[0],
                        directories["figures"] / "port_v4_vs_v5.png",
                    ),
                    port_v4_v5_figure(
                        v4_build.mesh,
                        build.mesh,
                        ports[0],
                        directories["figures"] / "port_v4_v5_comparison.png",
                    ),
                    junction_transition_wireframe_figure(
                        build.mesh,
                        face_region_labels(build.mesh, build.hybrid_details, branches),
                        directories["figures"] / "junction_transition_wireframe.png",
                    ),
                    junction_transition_normals_figure(
                        build.mesh,
                        face_region_labels(build.mesh, build.hybrid_details, branches),
                        directories["figures"] / "junction_transition_normals.png",
                    ),
                )
            )
            junction_id = (
                49
                if 49 in build.hybrid_details.patches
                else min(build.hybrid_details.patches)
            )
            junction_patch = build.hybrid_details.patches[junction_id]
            comparison_extent = max(
                (
                    collar.implicit_extent_um + 2.0 * collar.collar_radius_um
                    for collar in junction_patch.collars
                ),
                default=5.0 * float(roi.local_node_radius_um[junction_id]),
            )
            figure_paths.append(
                junction_v4_v5_figure(
                    v4_build.mesh,
                    build.mesh,
                    np.asarray(roi.local_node_positions_um[junction_id], dtype=float),
                    comparison_extent,
                    f"J{junction_id}",
                    directories["figures"] / "junction_v4_v5_same_camera.png",
                )
            )
        timings["visualization"] = time.perf_counter() - stage
        timings["total"] = time.perf_counter() - started_total
        final_surface_qc = (
            v9_result.report["selected_surface_qc"]
            if v9_result is not None
            else v8_result.report["selected_surface_qc"]
            if v8_result is not None
            else v7_result.report["surface_qc_v7"]
            if v7_result is not None
            and v7_result.decision == "ADOPT_V7_UNIFIED_POLYBALL"
            else
            v6_result.report["surface_qc"]
            if v6_result is not None
            else surface_qc
        )
        final_hybrid_qc = (
            v9_result.report["selected_topology"]
            if v9_result is not None
            else v8_result.report["selected_topology"]
            if v8_result is not None
            else v7_result.report["topology_v7"]
            if v7_result is not None
            and v7_result.decision == "ADOPT_V7_UNIFIED_POLYBALL"
            else
            v6_result.report["hybrid_qc"]
            if v6_result is not None
            else hybrid_qc
        )
        final_radius_qc = (
            v9_result.report["selected_radius_fidelity"]
            if v9_result is not None
            else v8_result.report["selected_radius_fidelity"]
            if v8_result is not None
            else v7_result.report["radius_fidelity_v7"]
            if v7_result is not None
            and v7_result.decision == "ADOPT_V7_UNIFIED_POLYBALL"
            else
            v6_result.report["radius_fidelity_v6"]
            if v6_result is not None
            else radius_qc
        )
        final_collar_rows = (
            v9_result.report["selected_collar_radius"]
            if v9_result is not None
            else v8_result.report["selected_collar_radius"]
            if v8_result is not None
            else v7_result.report["collar_radius_v7"]
            if v7_result is not None
            and v7_result.decision == "ADOPT_V7_UNIFIED_POLYBALL"
            else
            v6_result.report["collar_radius_v6"]
            if v6_result is not None
            else hybrid_collar_rows
        )
        maximum_collar_error = _maximum_collar_error(final_collar_rows)
        v5_maximum_collar_error = _maximum_collar_error(hybrid_collar_rows)
        v4_acceptance = _v4_acceptance_report(
            core_roi,
            context.report,
            branch_local_qc,
            controlled_comparison,
            hybrid_qc,
            radius_qc,
            v5_maximum_collar_error,
        )
        write_json(directories["qc"] / "v4_acceptance_report.json", v4_acceptance)
        summary = {
            "roi_id": core_roi.roi_id,
            "cfd_domain_id": roi.roi_id,
            "status": "PASS",
            "node_count": core_roi.node_count,
            "edge_count": core_roi.edge_count,
            "core_node_count": core_roi.node_count,
            "core_edge_count": core_roi.edge_count,
            "cfd_node_count": roi.node_count,
            "cfd_edge_count": roi.edge_count,
            "context_domain": context.report,
            "branch_count": len(branches),
            "bifurcation_count": pre_qc["bifurcation_count"],
            "cut_port_count": len(ports),
            "resampled_point_count": pre_qc["resampled_point_count"],
            "hard_collision_count": collision_qc["hard_collision_count"],
            "near_contact_count": collision_qc["near_contact_count"],
            "tube_primitive_count": build.tube_primitive_count,
            "junction_primitive_count": build.junction_primitive_count,
            "extension_primitive_count": build.extension_primitive_count,
            "backend_requested": build.backend_requested,
            "backend_used": (
                v9_result.selected_build.metadata.get("backend", "unified_polyball")
                if v9_result is not None
                else "unified_polyball_smooth_junction"
                if v8_result is not None
                and v8_result.decision == "ADOPT_V8_LOCAL_SMOOTH_UNION"
                else "unified_polyball"
                if v8_result is not None
                else "unified_polyball"
                if v7_result is not None
                and v7_result.decision == "ADOPT_V7_UNIFIED_POLYBALL"
                else
                "hybrid_continuous_implicit_field"
                if v6_result is not None
                else build.backend_used
            ),
            "fallback_reason": (
                (
                    v9_result.decision
                    if v9_result is not None
                    and not v9_result.decision.startswith("ADOPT_V9")
                    else v8_result.decision
                    if v8_result is not None
                    and v8_result.decision != "ADOPT_V8_LOCAL_SMOOTH_UNION"
                    else "v7 A/B decision retained v6"
                    if v7_result is not None
                    and v7_result.decision == "KEEP_V6_HYBRID"
                    else v6_result.details.transition_fallback_reason
                )
                if v6_result is not None
                else build.fallback_reason
            ),
            "formal_hybrid_status": (
                "PASS"
                if build.backend_used in hybrid_backends
                else "BLOCKED_USING_GLOBAL_FALLBACK"
                if config.junction.backend == "local_implicit" and build.fallback_reason
                else "NOT_REQUESTED"
            ),
            "triangle_count": final_surface_qc["triangle_count"],
            "vertex_count": final_surface_qc["vertex_count"],
            "surface_area_um2": final_surface_qc["surface_area_um2"],
            "enclosed_volume_um3": final_surface_qc["enclosed_volume_um3"],
            "watertight": final_surface_qc["checks"]["watertight_trimesh"] and final_surface_qc["checks"]["watertight_vtk"],
            "manifold": final_surface_qc["checks"]["manifold"],
            "surface_component_count": final_surface_qc["surface_component_count"],
            "detected_port_patch_count": final_surface_qc["detected_port_patch_count"],
            "radius_fidelity": final_radius_qc,
            "radius_p95_error": final_radius_qc["p95_absolute_relative_error"],
            "convergence": (
                v7_result.report["convergence"]
                if v7_result is not None
                else
                v6_result.report["convergence"]
                if v6_result is not None
                else convergence
            ),
            "volume_mesh": volume_mesh,
            "hybrid_qc": final_hybrid_qc,
            "branch_local_cross_section_qc": branch_local_qc,
            "controlled_local_implicit_comparison": controlled_comparison,
            "hybrid_max_collar_radius_error": maximum_collar_error,
            "v4_acceptance": v4_acceptance,
            "surface_continuity_version": (
                "v9"
                if v9_result is not None
                and v9_result.decision.startswith("ADOPT_V9")
                else "v8"
                if v8_result is not None
                and v8_result.decision == "ADOPT_V8_LOCAL_SMOOTH_UNION"
                else "v7"
                if v8_result is not None
                else "v7"
                if v7_result is not None
                and v7_result.decision == "ADOPT_V7_UNIFIED_POLYBALL"
                else "v6"
                if v6_result is not None
                else build.surface_continuity_version
            ),
            "surface_continuity": v5_continuity,
            "v4_v5_comparison": v5_comparison,
            "v5_acceptance": v5_comparison,
            "v6_acceptance": (
                v6_result.report["acceptance"]
                if v6_result is not None
                else None
            ),
            "v6_report": v6_result.report if v6_result is not None else None,
            "v7_decision": v7_result.decision if v7_result is not None else None,
            "v7_acceptance": (
                v7_result.report["acceptance"] if v7_result is not None else None
            ),
            "v7_report": v7_result.report if v7_result is not None else None,
            "v8_decision": v8_result.decision if v8_result is not None else None,
            "v8_report": v8_result.report if v8_result is not None else None,
            "v9_decision": v9_result.decision if v9_result is not None else None,
            "v9_acceptance": (
                v9_result.report["acceptance"] if v9_result is not None else None
            ),
            "v9_report": v9_result.report if v9_result is not None else None,
            "hybrid_junction_max_area_ratio": max_junction_area_ratio(hybrid_area_rows),
            "geometry_paths": [str(path) for path in geometry_paths],
            "figure_paths": [
                str(path)
                for path in (
                    *figure_paths,
                    *(v6_result.figure_paths if v6_result is not None else []),
                    *(v7_result.figure_paths if v7_result is not None else []),
                    *(v8_result.figure_paths if v8_result is not None else []),
                    *(v9_result.figure_paths if v9_result is not None else []),
                )
            ],
            "timings_s": timings,
        }
        write_json(directories["qc"] / "roi_summary.json", summary)
        return ROIProcessResult(core_roi.roi_id, "PASS", directories["root"], summary)
    except Exception as exc:
        timings["total"] = time.perf_counter() - started_total
        if pre_qc.get("status") == "NOT_RUN" or isinstance(exc, GeometryValidationError):
            failures = getattr(exc, "failures", [])
            if failures:
                pre_qc.update({"status": "FAIL", "validation_failures": failures})
            write_json(directories["qc"] / "pre_geometry_qc.json", pre_qc)
        failure_reason = f"{type(exc).__name__}: {exc}"
        failure = {
            "roi_id": core_roi.roi_id,
            "status": "FAIL",
            "failure_reason": failure_reason,
            "timings_s": timings,
            "traceback": traceback.format_exc(),
        }
        write_json(directories["qc"] / "failure.json", failure)
        return ROIProcessResult(core_roi.roi_id, "FAIL", directories["root"], failure, failure_reason)


def _worker(arguments: tuple[ROIRecord, CFDLumenConfig, CFDRunLayout, Path]) -> ROIProcessResult:
    return process_roi(*arguments)


def run_cfd_lumen_batch(
    rois: list[ROIRecord],
    config: CFDLumenConfig,
    run_layout: CFDRunLayout,
    sampling_run: Path,
    *,
    workers: int | None = None,
) -> list[ROIProcessResult]:
    if not rois:
        raise ValueError("No ROIs were provided for CFD lumen generation")
    worker_count = workers or min(max((os.cpu_count() or 2) - 1, 1), len(rois))
    if worker_count < 1:
        raise ValueError("workers must be positive")
    LOGGER.info("Running formal v5 surface-continuity synthetic controls")
    synthetic_started = time.perf_counter()
    synthetic_report = run_v5_synthetic_controls(config)
    synthetic_report["runtime_s"] = time.perf_counter() - synthetic_started
    write_json(run_layout.report / "v5_synthetic_controls.json", synthetic_report)
    if synthetic_report["status"] != "PASS":
        failed = [
            str(row.get("name"))
            for row in synthetic_report["cases"]
            if row.get("status") != "PASS"
        ]
        raise GeometryValidationError(
            "V5_SYNTHETIC_CONTROLS_FAILED: " + ", ".join(failed)
        )
    LOGGER.info("Processing %d ROI(s) from %s with %d worker(s)", len(rois), sampling_run, worker_count)
    results: list[ROIProcessResult] = []
    if worker_count == 1:
        for roi in rois:
            LOGGER.info("Starting ROI %s", roi.roi_id)
            result = process_roi(roi, config, run_layout, sampling_run)
            results.append(result)
            LOGGER.info("ROI %s: %s%s", roi.roi_id, result.status, f" ({result.failure_reason})" if result.failure_reason else "")
    else:
        with ProcessPoolExecutor(max_workers=worker_count) as executor:
            futures = {
                executor.submit(_worker, (roi, config, run_layout, sampling_run)): roi.roi_id
                for roi in rois
            }
            for future in as_completed(futures):
                roi_id = futures[future]
                try:
                    result = future.result()
                except Exception as exc:
                    failure_reason = f"WorkerFailure: {type(exc).__name__}: {exc}"
                    result = ROIProcessResult(roi_id, "FAIL", run_layout.rois / roi_id, {"roi_id": roi_id, "status": "FAIL", "failure_reason": failure_reason}, failure_reason)
                results.append(result)
                LOGGER.info("ROI %s: %s%s", roi_id, result.status, f" ({result.failure_reason})" if result.failure_reason else "")
    rank = {roi.roi_id: index for index, roi in enumerate(rois)}
    results.sort(key=lambda result: rank[result.roi_id])
    rows = []
    failed_rows = []
    for result in results:
        summary = result.summary
        row = {
            "roi_id": result.roi_id,
            "status": result.status,
            "branch_count": summary.get("branch_count"),
            "cut_port_count": summary.get("cut_port_count"),
            "triangle_count": summary.get("triangle_count"),
            "enclosed_volume_um3": summary.get("enclosed_volume_um3"),
            "watertight": summary.get("watertight"),
            "manifold": summary.get("manifold"),
            "radius_p95_error": summary.get("radius_p95_error"),
            "backend_used": summary.get("backend_used"),
            "runtime_s": summary.get("timings_s", {}).get("total") if isinstance(summary.get("timings_s"), dict) else None,
            "output_dir": str(result.output_dir),
            "failure_reason": result.failure_reason,
        }
        rows.append(row)
        if result.status != "PASS":
            failed_rows.append(row)
    write_csv(run_layout.manifests / "geometry_summary.csv", rows)
    write_csv(
        run_layout.manifests / "failed_rois.csv",
        failed_rows,
        fieldnames=rows[0].keys() if rows else ("roi_id", "status", "failure_reason"),
    )
    passed = sum(result.status == "PASS" for result in results)
    run_status = "PASS" if passed == len(results) else ("PARTIAL" if passed else "FAIL")
    manifest = {
        "status": run_status,
        "sampling_run": str(Path(sampling_run).resolve()),
        "roi_ids": [roi.roi_id for roi in rois],
        "requested_roi_count": len(rois),
        "passed_roi_count": passed,
        "failed_roi_count": len(results) - passed,
        "workers": worker_count,
        "backend": config.boolean.backend,
        "geometry_source": "representative SWC ROI",
        "segmentation_mask_used": False,
        "v5_synthetic_controls": {
            "status": synthetic_report["status"],
            "case_count": synthetic_report["case_count"],
            "runtime_s": synthetic_report["runtime_s"],
        },
    }
    write_json(run_layout.manifests / "run_manifest.json", manifest)
    write_json(run_layout.report / "run_summary.json", {**manifest, "rois": rows})
    passed_rows = [row for row in rows if row["status"] == "PASS"]
    if passed_rows:
        run_summary_figure(passed_rows, run_layout.figures / "all_roi_geometry_summary.png")
    return results
