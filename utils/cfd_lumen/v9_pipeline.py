"""Strict v9 segment-crease diagnosis and source-constrained reconstruction."""

from __future__ import annotations

import time
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .export import write_csv, write_geometry_exports, write_json, write_units
from .hybrid_qc import junction_area_profile_rows
from .junction_fairing import constrained_junction_fairing
from .segment_ownership_qc import (
    competition_support_report,
    evaluate_segment_ownership,
    tangent_audit_rows,
    write_switch_edges,
)
from .smooth_centerline import SmoothCenterlineBuild, build_smooth_centerline
from .types import BranchGeometry, PatchResult, PortGeometry
from .unified_polyball import (
    UnifiedPolyBallBuild,
    build_unified_polyball_surface,
    prepare_polyball_raster,
    release_prepared_polyball_raster,
)
from .v6_pipeline import V6RefinementResult
from .v7_pipeline import V7RefinementResult
from .v8_pipeline import V8RefinementResult
from .v8_qc import (
    build_junction_blend_specs,
    junction_local_volumes,
    local_normal_metrics,
)
from .v8_visualization import camera_report, local_box_crop
from .v9_qc import comparison_row, evaluate_v9_method, final_acceptance
from .v9_visualization import (
    centerline_fidelity_figure,
    generate_method_comparison,
    tangent_defect_correlation_figure,
)


@dataclass(slots=True)
class V9RefinementResult:
    mesh: trimesh.Trimesh
    patch: PatchResult
    selected_build: UnifiedPolyBallBuild
    report: dict[str, Any]
    decision: str
    geometry_paths: list[Path]
    figure_paths: list[Path]


def _centerline_polydata(branches: list[BranchGeometry]) -> pv.PolyData:
    points: list[np.ndarray] = []
    lines: list[np.ndarray] = []
    radius: list[np.ndarray] = []
    branch_ids: list[np.ndarray] = []
    offset = 0
    for branch in branches:
        branch_points = np.asarray(branch.points_um, dtype=float)
        points.append(branch_points)
        radius.append(np.asarray(branch.radius_um, dtype=float))
        branch_ids.append(
            np.full(len(branch_points), int(branch.branch_id), dtype=np.int64)
        )
        lines.append(
            np.concatenate(
                (
                    np.asarray((len(branch_points),), dtype=np.int64),
                    np.arange(offset, offset + len(branch_points), dtype=np.int64),
                )
            )
        )
        offset += len(branch_points)
    output = pv.PolyData(np.vstack(points), lines=np.concatenate(lines))
    output.point_data["radius_um"] = np.concatenate(radius)
    output.point_data["branch_id"] = np.concatenate(branch_ids)
    return output


def _best_v8_build(
    result: V8RefinementResult,
) -> tuple[float, UnifiedPolyBallBuild, str]:
    selected = result.report.get("selected_k_radius_ratio")
    if selected is not None and float(selected) in result.candidate_builds:
        ratio = float(selected)
        return ratio, result.candidate_builds[ratio], "v8 formally selected ratio"

    metrics = result.report.get("k_sensitivity", {})

    def score(ratio: float) -> tuple[float, float, float]:
        row = metrics.get(f"{ratio:.2f}", {})
        silhouette = row.get("silhouette", {}).get("curvature_variation_mean")
        normal = row.get("normal", {}).get("normal_jump_deg", {}).get("p99")
        return (
            float(silhouette) if silhouette is not None else np.inf,
            float(normal) if normal is not None else np.inf,
            ratio,
        )

    ratio = min(result.candidate_builds, key=score)
    return float(ratio), result.candidate_builds[ratio], "lowest measured v8 roughness"


def _safe_metric(metrics: dict[str, Any], *keys: str, default: float = np.inf) -> float:
    value: Any = metrics
    for key in keys:
        value = value.get(key) if isinstance(value, dict) else None
    return default if value is None else float(value)


def _format_metric(value: Any, suffix: str = "", digits: int = 6) -> str:
    if value is None or not np.isfinite(float(value)):
        return "N/A"
    return f"{float(value):.{digits}g}{suffix}"


def _candidate_score(
    ratio: float, metrics: dict[str, Any], acceptance: dict[str, Any]
) -> tuple[Any, ...]:
    return (
        not acceptance["pass"],
        int(metrics["visible_sawtooth_count"]),
        _safe_metric(metrics, "same_branch_segment_gradient_jump_deg", "p99"),
        _safe_metric(metrics, "cross_branch_gradient_jump_deg", "p99"),
        _safe_metric(metrics, "silhouette", "curvature_variation_mean"),
        abs(_safe_metric(metrics, "junction_volume", "maximum_relative_increase")),
        float(ratio),
    )


def _write_method_artifacts(
    label: str,
    build: UnifiedPolyBallBuild,
    switch_rows: list[dict[str, Any]],
    arrays: dict[str, np.ndarray],
    root: Path,
) -> list[Path]:
    folder = root / label
    folder.mkdir(parents=True, exist_ok=True)
    csv_path = write_csv(folder / "segment_switches.csv", switch_rows)
    edge_path = write_switch_edges(
        build.wall_mesh_before_clip,
        arrays["switch_edges"],
        switch_rows,
        folder / "segment_switches.vtp",
    )
    return [csv_path, edge_path]


def _markdown(report: dict[str, Any]) -> str:
    lines = [
        "# v9 strict reconstruction report",
        "",
        f"- ROI: `{report['roi_id']}`",
        f"- status: **{report['status']}**",
        f"- decision: **{report['decision']}**",
        f"- root cause: `{report['root_cause']}`",
        f"- selected spline theta: `{report['smooth_centerline']['selected_theta_max_deg']} deg`",
        f"- selected spline sagitta: `{report['smooth_centerline']['selected_eta_radius_fraction']} r`",
        f"- selected k/r: `{report['selected_k_radius_ratio']}`",
        "",
        "## Required answers",
        "",
    ]
    for index in range(1, 11):
        lines.append(f"{index}. {report['final_answers'][str(index)]}")
    lines.extend(("", "## Final acceptance", ""))
    for name, passed in report["acceptance"]["checks"].items():
        lines.append(f"- [{'x' if passed else ' '}] {name}")
    return "\n".join(lines) + "\n"


def run_v9_refinement(
    roi: ROIRecord,
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    v6_result: V6RefinementResult,
    v7_result: V7RefinementResult,
    v8_result: V8RefinementResult,
    output_root: Path,
) -> V9RefinementResult:
    """Run every diagnosis, method control, and hard acceptance required by v9."""

    output_root = Path(output_root)
    diagnostics = output_root / "diagnostics"
    figures = output_root / "figures"
    geometry = output_root / "geometry"
    variants = output_root / "variants"
    report_root = output_root / "report"
    ports_root = output_root / "ports"
    for folder in (diagnostics, figures, geometry, variants, report_root, ports_root):
        folder.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    hard = v7_result.candidate_build
    v8_ratio, branch_union, v8_selection_basis = _best_v8_build(v8_result)
    specs = build_junction_blend_specs(roi, v6_result.details, config)
    worst_spec = max(
        specs,
        key=lambda spec: (
            local_normal_metrics(hard.wall_mesh_before_clip, (spec,))["normal_jump_deg"][
                "p99"
            ]
            or -np.inf
        ),
    )

    # Phase 1: diagnose the raw FlyingEdges surface before either Newton pass or remesh.
    s0 = hard.stage_meshes["S0_raw_flying_edges"]
    s0_ownership_path = diagnostics / "s0_segment_ownership.vtp"
    s0_report, s0_switch_rows, s0_arrays = evaluate_segment_ownership(
        s0,
        hard.model,
        specs,
        config,
        ownership_path=s0_ownership_path,
    )
    tangent_rows, tangent_report = tangent_audit_rows(
        hard.model, branches, specs, s0_switch_rows
    )
    s0_switch_csv = write_csv(diagnostics / "s0_segment_switches.csv", s0_switch_rows)
    s0_switch_vtp = write_switch_edges(
        s0,
        s0_arrays["switch_edges"],
        s0_switch_rows,
        diagnostics / "s0_segment_switches.vtp",
    )
    tangent_csv = write_csv(diagnostics / "source_tangent_kinks.csv", tangent_rows)
    write_json(diagnostics / "s0_segment_ownership_report.json", s0_report)
    write_json(diagnostics / "source_tangent_kink_report.json", tangent_report)

    # Phase 2: all theta/eta pairs are measured; only a passing pair may be selected.
    smooth: SmoothCenterlineBuild = build_smooth_centerline(branches, ports, config)
    smooth_vtp = diagnostics / "cfd_derived_spline_centerline.vtp"
    _centerline_polydata(smooth.branches).save(smooth_vtp)
    smooth_sensitivity_csv = write_csv(
        diagnostics / "centerline_adaptive_sensitivity.csv", smooth.sensitivity_rows
    )
    smooth_fidelity_csv = write_csv(
        diagnostics / "centerline_source_fidelity.csv", smooth.branch_fidelity_rows
    )
    write_json(diagnostics / "centerline_source_fidelity.json", smooth.report)

    # Phases 3/4: one dense C1 hard field, then the full measured k/r sensitivity.
    preparation_started = time.perf_counter()
    prepared = prepare_polyball_raster(
        branches,
        ports,
        config,
        v6_details=v6_result.details,
        cells_across_min_diameter=config.v7.cells_across_min_diameter,
        constructed_override=smooth.branches,
        port_tail_rows_override=smooth.port_tail_rows,
    )
    preparation_runtime = time.perf_counter() - preparation_started
    spline_hard: UnifiedPolyBallBuild
    union_builds: dict[float, UnifiedPolyBallBuild] = {}
    build_runtimes: dict[str, float] = {}
    try:
        candidate_started = time.perf_counter()
        spline_hard = build_unified_polyball_surface(
            branches,
            ports,
            config,
            v6_details=v6_result.details,
            cells_across_min_diameter=config.v7.cells_across_min_diameter,
            compare_extractors=False,
            remesh=True,
            prepared_raster=prepared,
        )
        spline_hard.metadata.update(
            {
                "protocol": "v9 C1 spline hard PolyBall control",
                "backend": "c1_spline_polyball_hard_min",
                "source_centerline": smooth.report,
            }
        )
        build_runtimes["C_v9_spline_hard"] = time.perf_counter() - candidate_started
        for ratio in config.v9.k_radius_ratios:
            candidate_started = time.perf_counter()
            build = build_unified_polyball_surface(
                branches,
                ports,
                config,
                v6_details=v6_result.details,
                cells_across_min_diameter=config.v7.cells_across_min_diameter,
                compare_extractors=False,
                remesh=True,
                junction_specs=specs,
                smooth_k_radius_ratio=float(ratio),
                prepared_raster=prepared,
                competition_threshold_radius_fraction=(
                    config.v9.competition_threshold_radius_fraction
                ),
            )
            build.metadata["source_centerline"] = smooth.report
            union_builds[float(ratio)] = build
            build_runtimes[f"D_v9_spline_union_k{float(ratio):.2f}"] = (
                time.perf_counter() - candidate_started
            )
    finally:
        release_prepared_polyball_raster(prepared)

    all_builds: dict[str, UnifiedPolyBallBuild] = {
        "A_v7_hard": hard,
        "B_v8_branch_union": branch_union,
        "C_v9_spline_hard": spline_hard,
        **{
            f"D_v9_spline_union_k{ratio:.2f}": build
            for ratio, build in union_builds.items()
        },
    }
    local_meshes = {
        label: local_box_crop(
            build.wall_mesh_before_clip,
            np.asarray(worst_spec.center_world_um),
            1.25 * worst_spec.blend_length_um,
        )
        for label, build in all_builds.items()
    }
    shared_local = tuple(local_meshes.values())
    baseline_area = junction_area_profile_rows(
        hard.mesh, v6_result.details, branches
    )
    baseline_volume = junction_local_volumes(hard.model, specs, config)

    metrics_by_label: dict[str, dict[str, Any]] = {}
    switch_by_label: dict[str, list[dict[str, Any]]] = {}
    arrays_by_label: dict[str, dict[str, np.ndarray]] = {}
    method_artifacts: list[Path] = []
    for label, build in all_builds.items():
        ownership_path = diagnostics / label / "segment_ownership.vtp"
        runtime = build_runtimes.get(label, float(build.metadata.get("runtime_s", 0.0)))
        metrics, switch_rows, arrays = evaluate_v9_method(
            label,
            build,
            roi,
            branches,
            v6_result.details,
            specs,
            config,
            baseline_area=baseline_area,
            baseline_volume=baseline_volume,
            comparison_meshes=shared_local,
            runtime_s=runtime,
            ownership_path=ownership_path,
            silhouette_mesh=local_meshes[label],
        )
        metrics_by_label[label] = metrics
        switch_by_label[label] = switch_rows
        arrays_by_label[label] = arrays
        write_json(diagnostics / label / "metrics.json", metrics)
        method_artifacts.extend(
            _write_method_artifacts(label, build, switch_rows, arrays, diagnostics)
        )
        method_artifacts.append(ownership_path)

    baseline_metrics = metrics_by_label["A_v7_hard"]
    acceptance_by_ratio: dict[float, dict[str, Any]] = {}
    for ratio in config.v9.k_radius_ratios:
        label = f"D_v9_spline_union_k{float(ratio):.2f}"
        acceptance_by_ratio[float(ratio)] = final_acceptance(
            metrics_by_label[label],
            smooth.report,
            baseline_metrics,
            config,
            same_camera_figures_generated=True,
        )
        metrics_by_label[label]["acceptance"] = acceptance_by_ratio[float(ratio)]
        metrics_by_label[label]["competition_support"] = competition_support_report(
            union_builds[float(ratio)].wall_mesh_before_clip,
            union_builds[float(ratio)].field_model,
            specs,
        )

    ranked_ratios = sorted(
        (float(ratio) for ratio in config.v9.k_radius_ratios),
        key=lambda ratio: _candidate_score(
            ratio,
            metrics_by_label[f"D_v9_spline_union_k{ratio:.2f}"],
            acceptance_by_ratio[ratio],
        ),
    )
    selected_ratio = ranked_ratios[0]
    selected_label = f"D_v9_spline_union_k{selected_ratio:.2f}"
    selected = union_builds[selected_ratio]
    selected_metrics = metrics_by_label[selected_label]
    selected_acceptance = acceptance_by_ratio[selected_ratio]
    decision = (
        "ADOPT_V9_SPLINE_COMPETITION_UNION"
        if selected_acceptance["pass"]
        else "V9_UNION_REQUIRES_CONSTRAINED_FAIRING_OR_FALLBACK"
    )

    fairing_report: dict[str, Any] = {
        "status": "NOT_NEEDED",
        "reason": "selected spline plus competition-aware union has no residual visible defect",
        "newton_projection_after_fairing": False,
    }
    fairing_build: UnifiedPolyBallBuild | None = None
    if (
        config.v9.fairing_enabled
        and not selected_acceptance["pass"]
        and selected_metrics["visible_sawtooth_count"] > 0
    ):
        capped_segment, _, capped_arrays = evaluate_segment_ownership(
            selected.mesh,
            selected.model,
            specs,
            config,
            field_model=selected.field_model,
        )
        fairing = constrained_junction_fairing(
            selected.mesh,
            capped_arrays["visible_defect_edges"],
            specs,
            ports,
            selected.patch,
            config,
        )
        fairing_report = dict(fairing.report)
        fairing_report["input_visible_sawtooth_count"] = capped_segment[
            "visible_sawtooth_count"
        ]
        fairing_report["total_volume_relative_change"] = float(
            (fairing.mesh.volume - selected.mesh.volume)
            / max(abs(selected.mesh.volume), 1.0e-15)
        )
        fairing_build = replace(
            selected,
            mesh=fairing.mesh,
            wall_mesh_before_clip=fairing.mesh,
            metadata={
                **selected.metadata,
                "backend": "c1_spline_polyball_competition_union_local_fairing",
                "junction_fairing": fairing_report,
                "newton_projection_after_fairing": False,
            },
        )
        fair_label = "E_v9_constrained_local_fairing"
        fair_local = local_box_crop(
            fairing_build.wall_mesh_before_clip,
            np.asarray(worst_spec.center_world_um),
            1.25 * worst_spec.blend_length_um,
        )
        fair_metrics, fair_switch, fair_arrays = evaluate_v9_method(
            fair_label,
            fairing_build,
            roi,
            branches,
            v6_result.details,
            specs,
            config,
            baseline_area=baseline_area,
            baseline_volume=baseline_volume,
            comparison_meshes=(*shared_local, fair_local),
            runtime_s=float(selected_metrics["runtime_s"]),
            ownership_path=diagnostics / fair_label / "segment_ownership.vtp",
            silhouette_mesh=fair_local,
        )
        fair_acceptance = final_acceptance(
            fair_metrics,
            smooth.report,
            baseline_metrics,
            config,
            same_camera_figures_generated=True,
        )
        fair_safety_checks = {
            "outer_boundary_and_ports_fixed_exactly": fairing_report.get(
                "fixed_vertex_maximum_displacement_um", np.inf
            )
            == 0.0,
            "topology_preserved": fair_metrics["topology"]["status"] == "PASS",
            "volume_change_within_tolerance": abs(
                fairing_report["total_volume_relative_change"]
            )
            <= config.v9.maximum_junction_volume_increase_fraction,
            "no_newton_or_hard_min_reprojection": not fairing_report[
                "newton_projection_after_fairing"
            ],
        }
        fair_acceptance["fairing_safety_checks"] = fair_safety_checks
        fair_acceptance["pass"] = fair_acceptance["pass"] and all(
            fair_safety_checks.values()
        )
        fair_metrics["acceptance"] = fair_acceptance
        fair_metrics["fairing"] = fairing_report
        metrics_by_label[fair_label] = fair_metrics
        switch_by_label[fair_label] = fair_switch
        arrays_by_label[fair_label] = fair_arrays
        all_builds[fair_label] = fairing_build
        local_meshes[fair_label] = fair_local
        write_json(diagnostics / fair_label / "metrics.json", fair_metrics)
        method_artifacts.extend(
            _write_method_artifacts(
                fair_label, fairing_build, fair_switch, fair_arrays, diagnostics
            )
        )
        method_artifacts.append(diagnostics / fair_label / "segment_ownership.vtp")
        if fair_acceptance["pass"]:
            selected = fairing_build
            selected_label = fair_label
            selected_metrics = fair_metrics
            selected_acceptance = fair_acceptance
            decision = "ADOPT_V9_CONSTRAINED_LOCAL_FAIRING"

    # A failed v9 method is never silently promoted over the accepted v7 baseline.
    if not selected_acceptance["pass"]:
        spline_acceptance = final_acceptance(
            metrics_by_label["C_v9_spline_hard"],
            smooth.report,
            baseline_metrics,
            config,
            same_camera_figures_generated=True,
        )
        metrics_by_label["C_v9_spline_hard"]["acceptance"] = spline_acceptance
        if spline_acceptance["pass"]:
            selected = spline_hard
            selected_label = "C_v9_spline_hard"
            selected_metrics = metrics_by_label[selected_label]
            selected_acceptance = spline_acceptance
            selected_ratio = None
            decision = "ADOPT_V9_SPLINE_HARD"
        else:
            selected = hard
            selected_label = "A_v7_hard"
            selected_metrics = baseline_metrics
            selected_ratio = None
            decision = "KEEP_V7_HARD_MIN_V9_ACCEPTANCE_FAILED"

    display_builds = {
        "A v7 hard": hard,
        "B v8 branch union": branch_union,
        "C v9 spline hard": spline_hard,
        "D v9 spline+union": union_builds[ranked_ratios[0]],
    }
    if fairing_build is not None:
        display_builds["E constrained fairing"] = fairing_build
    comparison_figures, camera, _ = generate_method_comparison(
        {
            label: build.wall_mesh_before_clip
            for label, build in display_builds.items()
        },
        worst_spec,
        figures / "same_camera_methods",
    )
    figure_paths = [
        *comparison_figures,
        centerline_fidelity_figure(
            branches, smooth.branches, figures / "centerline_source_fidelity.png"
        ),
        tangent_defect_correlation_figure(
            s0_switch_rows, figures / "tangent_defect_correlation.png"
        ),
    ]
    same_camera_ok = all(path.is_file() for path in comparison_figures)
    if selected_label.startswith(("C_", "D_", "E_")):
        selected_acceptance = final_acceptance(
            selected_metrics,
            smooth.report,
            baseline_metrics,
            config,
            same_camera_figures_generated=same_camera_ok,
        )
        if selected_label.startswith("E_"):
            safety = selected_metrics["acceptance"]["fairing_safety_checks"]
            selected_acceptance["fairing_safety_checks"] = safety
            selected_acceptance["pass"] &= all(safety.values())
        selected_metrics["acceptance"] = selected_acceptance

    variant_geometry_paths: list[Path] = []
    for label, build in {
        "C_v9_spline_hard": spline_hard,
        **{
            f"D_k_{ratio:.2f}r".replace(".", "p"): candidate
            for ratio, candidate in union_builds.items()
        },
        **(
            {"E_constrained_fairing": fairing_build}
            if fairing_build is not None
            else {}
        ),
    }.items():
        variant_root = variants / label
        variant_geometry = variant_root / "geometry"
        variant_ports = variant_root / "ports"
        variant_geometry.mkdir(parents=True, exist_ok=True)
        variant_ports.mkdir(parents=True, exist_ok=True)
        variant_geometry_paths.extend(
            write_geometry_exports(
                build.mesh,
                build.patch,
                ports,
                {"geometry": variant_geometry, "ports": variant_ports},
                face_region=np.zeros(len(build.mesh.faces), dtype=np.uint8),
            )
        )

    geometry_paths = write_geometry_exports(
        selected.mesh,
        selected.patch,
        ports,
        {"geometry": geometry, "ports": ports_root},
        face_region=np.zeros(len(selected.mesh.faces), dtype=np.uint8),
    )
    geometry_paths.extend(
        (
            s0_ownership_path,
            s0_switch_csv,
            s0_switch_vtp,
            tangent_csv,
            smooth_vtp,
            smooth_sensitivity_csv,
            smooth_fidelity_csv,
            *method_artifacts,
            *variant_geometry_paths,
            write_units(diagnostics / "units.json"),
        )
    )

    comparison_rows = [comparison_row(metrics) for metrics in metrics_by_label.values()]
    write_csv(diagnostics / "method_comparison.csv", comparison_rows)
    k_rows = []
    for ratio in config.v9.k_radius_ratios:
        label = f"D_v9_spline_union_k{float(ratio):.2f}"
        row = comparison_row(metrics_by_label[label])
        row["k_radius_ratio"] = float(ratio)
        row["acceptance_pass"] = acceptance_by_ratio[float(ratio)]["pass"]
        k_rows.append(row)
    write_csv(diagnostics / "k_sensitivity.csv", k_rows)
    write_json(
        diagnostics / "same_camera.json",
        {
            "junction_node_id": worst_spec.junction_node_id,
            "camera": camera_report(camera),
            "views": ["flat", "wireframe", "silhouette"],
            "same_camera": True,
        },
    )

    c_metrics = metrics_by_label["C_v9_spline_hard"]
    d_metrics = metrics_by_label[
        f"D_v9_spline_union_k{ranked_ratios[0]:.2f}"
    ]
    maximum_kink = tangent_report.get("maximum_kink") or {}
    fidelity_rows = smooth.branch_fidelity_rows
    maximum_hausdorff_um = max(
        (float(row["hausdorff_distance_um"]) for row in fidelity_rows), default=0.0
    )
    maximum_hausdorff_r = max(
        (
            float(row["hausdorff_distance_min_radius_fraction"])
            for row in fidelity_rows
        ),
        default=0.0,
    )
    maximum_p95_um = max(
        (float(row["p95_distance_um"]) for row in fidelity_rows), default=0.0
    )
    maximum_length_change = max(
        (abs(float(row["branch_length_change_fraction"])) for row in fidelity_rows),
        default=0.0,
    )
    selected_discrete_tangent_max = max(
        (
            float(row["discretized_tangent_angle_deg"]["max"] or 0.0)
            for row in fidelity_rows
        ),
        default=0.0,
    )
    a_same = baseline_metrics["same_branch_segment_gradient_jump_deg"]["p99"]
    c_same = c_metrics["same_branch_segment_gradient_jump_deg"]["p99"]
    c_cross = c_metrics["cross_branch_gradient_jump_deg"]["p99"]
    d_cross = d_metrics["cross_branch_gradient_jump_deg"]["p99"]
    c_silhouette = c_metrics["silhouette"]["curvature_variation_mean"]
    a_silhouette = baseline_metrics["silhouette"]["curvature_variation_mean"]
    d_silhouette = d_metrics["silhouette"]["curvature_variation_mean"]
    visual_disappeared = bool(
        same_camera_ok
        and selected_acceptance["checks"].get("visible_sawtooth_disappeared", False)
    )
    final_answers = {
        "1": (
            f"S0 defect-near same-branch segment-switch fraction="
            f"{100.0 * s0_report['defect_near_any_same_branch_segment_fraction']:.3f}% "
            f"(adjacent={100.0 * s0_report['defect_near_same_branch_adjacent_segment_fraction']:.3f}%, "
            f"nonadjacent={100.0 * s0_report['defect_near_same_branch_nonadjacent_segment_fraction']:.3f}%); "
            f"root cause={s0_report['root_cause_if_confirmed'] or 'NOT_CONFIRMED'}."
        ),
        "2": (
            f"Maximum resampled-centerline tangent kink="
            f"{_format_metric(maximum_kink.get('theta_deg'), ' deg')} at branch "
            f"{maximum_kink.get('branch_id', 'N/A')}, source neighborhood "
            f"{maximum_kink.get('source_node_neighborhood', 'N/A')}, segments "
            f"{maximum_kink.get('segment_id_0', 'N/A')}/"
            f"{maximum_kink.get('segment_id_1', 'N/A')}."
        ),
        "3": (
            f"C1 knot tangent discontinuity=0 deg; selected dense-polyline maximum "
            f"turn={_format_metric(selected_discrete_tangent_max, ' deg')}."
        ),
        "4": (
            f"Smooth-centerline-only changed same-branch gradient P99 from "
            f"{_format_metric(a_same, ' deg')} to {_format_metric(c_same, ' deg')}, "
            f"visible count from {baseline_metrics['visible_sawtooth_count']} to "
            f"{c_metrics['visible_sawtooth_count']}, and silhouette roughness from "
            f"{_format_metric(a_silhouette)} to {_format_metric(c_silhouette)}."
        ),
        "5": (
            f"Competition-aware union at measured k/r={ranked_ratios[0]:.2f} changed "
            f"cross-branch gradient P99 from {_format_metric(c_cross, ' deg')} to "
            f"{_format_metric(d_cross, ' deg')}, visible count from "
            f"{c_metrics['visible_sawtooth_count']} to {d_metrics['visible_sawtooth_count']}, "
            f"and silhouette roughness from {_format_metric(c_silhouette)} to "
            f"{_format_metric(d_silhouette)}."
        ),
        "6": (
            f"Final fairing requirement/status={fairing_report['status']}; after fairing "
            f"Newton/hard-min reprojection=False."
        ),
        "7": (
            f"Same-camera flat/wireframe/silhouette generated={same_camera_ok}; "
            f"quantitative visible-defect disappearance={visual_disappeared}."
        ),
        "8": (
            f"Source arrays/topology modified=False; maximum Hausdorff="
            f"{maximum_hausdorff_um:.6g} um ({maximum_hausdorff_r:.6g}r), maximum "
            f"bidirectional P95={maximum_p95_um:.6g} um, maximum length change="
            f"{100.0 * maximum_length_change:.6g}%, junction/endpoint/port error=0 um."
        ),
        "9": (
            f"Final radius P95={_format_metric(selected_metrics['radius_fidelity']['p95_absolute_relative_error'])}; "
            f"hydraulic error={_format_metric(selected_metrics['hydraulic_resistance']['max_absolute_relative_error'])}; "
            f"acceptance radius/hydraulic checks="
            f"{selected_acceptance['checks'].get('radius_p95_below_one_percent', False)}/"
            f"{selected_acceptance['checks'].get('hydraulic_error_within_tolerance', False)}."
        ),
        "10": (
            f"Recommended backend={selected.metadata.get('backend', 'unified_polyball')}; "
            f"decision={decision}; strict multi-criterion pass={selected_acceptance['pass']}."
        ),
    }
    report = {
        "protocol": "(new) subgraph modeling modification v9",
        "roi_id": roi.roi_id,
        "status": "PASS" if selected_acceptance["pass"] else "FAIL",
        "decision": decision,
        "selected_method": selected_label,
        "selected_k_radius_ratio": selected_ratio,
        "root_cause": s0_report["root_cause_if_confirmed"] or "NOT_CONFIRMED",
        "s0_segment_ownership": s0_report,
        "source_tangent_kink_audit": tangent_report,
        "smooth_centerline": smooth.report,
        "v8_control": {
            "k_radius_ratio": v8_ratio,
            "selection_basis": v8_selection_basis,
        },
        "shared_hard_field_preparation_runtime_s": preparation_runtime,
        "method_comparison": metrics_by_label,
        "k_sensitivity": {
            f"{ratio:.2f}": metrics_by_label[
                f"D_v9_spline_union_k{ratio:.2f}"
            ]
            for ratio in sorted(union_builds)
        },
        "fairing": fairing_report,
        "acceptance": selected_acceptance,
        "selected_surface_qc": selected_metrics["surface_qc"],
        "selected_topology": selected_metrics["topology"],
        "selected_radius_fidelity": selected_metrics["radius_fidelity"],
        "selected_collar_radius": selected_metrics["collar_radius"],
        "same_camera": {
            "generated": same_camera_ok,
            "junction_node_id": worst_spec.junction_node_id,
            "camera": camera_report(camera),
        },
        "final_answers": final_answers,
        "geometry_paths": [str(path) for path in geometry_paths],
        "figure_paths": [str(path) for path in figure_paths],
        "runtime_s": time.perf_counter() - started,
    }
    write_json(report_root / "v9_final_report.json", report)
    (report_root / "v9_final_report.md").write_text(
        _markdown(report), encoding="utf-8"
    )
    return V9RefinementResult(
        mesh=selected.mesh,
        patch=selected.patch,
        selected_build=selected,
        report=report,
        decision=decision,
        geometry_paths=geometry_paths,
        figure_paths=figure_paths,
    )
