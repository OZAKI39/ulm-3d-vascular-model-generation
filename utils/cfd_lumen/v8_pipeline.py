"""Strict v8 localization, local smooth-union sensitivity, and selection."""

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
from .hybrid_qc import collar_radius_rows, junction_area_profile_rows
from .surface_qc import evaluate_radius_fidelity, evaluate_surface_qc
from .types import BranchGeometry, PatchResult, PortGeometry
from .unified_polyball import (
    UnifiedPolyBallBuild,
    build_unified_polyball_surface,
    prepare_polyball_raster,
    release_prepared_polyball_raster,
)
from .v6_pipeline import V6RefinementResult
from .v6_qc import silhouette_rows
from .v7_pipeline import V7RefinementResult
from .v7_qc import evaluate_unified_topology
from .v8_qc import (
    build_junction_blend_specs,
    compare_volume_rows,
    evaluate_ownership_switching,
    field_gradient_angles_on_edges,
    junction_local_volumes,
    local_hydraulic_resistance_from_area_profiles,
    local_normal_metrics,
    stage_localization_rows,
)
from .v8_visualization import (
    camera_report,
    generate_k_sensitivity_figures,
    generate_stage_localization_figures,
    mesh_polydata,
    ownership_and_overlay_figures,
)


@dataclass(slots=True)
class V8RefinementResult:
    mesh: trimesh.Trimesh
    patch: PatchResult
    selected_build: UnifiedPolyBallBuild
    report: dict[str, Any]
    decision: str
    geometry_paths: list[Path]
    figure_paths: list[Path]
    candidate_builds: dict[float, UnifiedPolyBallBuild]


def _maximum(rows: list[dict[str, Any]], key: str) -> float | None:
    values = [float(row[key]) for row in rows if row.get(key) is not None]
    return max(values, default=None)


def _mean(rows: list[dict[str, Any]], key: str) -> float | None:
    values = [float(row[key]) for row in rows if row.get(key) is not None]
    return float(np.mean(values)) if values else None


def _silhouette_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "curvature_variation_mean": _mean(
            rows, "silhouette_curvature_variation"
        ),
        "turning_p95_worst_deg": _maximum(rows, "silhouette_turning_p95_deg"),
        "large_corner_count": int(sum(int(row["large_corner_count"]) for row in rows)),
        "large_corner_fraction_mean": _mean(rows, "large_corner_fraction"),
    }


def _relative_reduction(baseline: float | None, candidate: float | None) -> float | None:
    if baseline is None or candidate is None or baseline <= 0.0:
        return None
    return float((baseline - candidate) / baseline)


def _write_edge_polydata(
    mesh: trimesh.Trimesh,
    edges: np.ndarray,
    rows: list[dict[str, Any]],
    path: Path,
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    edges = np.asarray(edges, dtype=np.int64).reshape((-1, 2))
    if not len(edges):
        pv.PolyData().save(path)
        return path
    points = np.asarray(mesh.vertices, dtype=float)[edges].reshape((-1, 3))
    lines = np.column_stack(
        (
            np.full(len(edges), 2, dtype=np.int64),
            2 * np.arange(len(edges), dtype=np.int64),
            2 * np.arange(len(edges), dtype=np.int64) + 1,
        )
    ).ravel()
    data = pv.PolyData(points, lines=lines)
    for key in (
        "winner_branch_0",
        "winner_branch_1",
        "gradient_jump_angle_deg",
        "surface_dihedral_angle_deg",
        "is_sawtooth_defect",
    ):
        data.cell_data[key] = np.asarray([row[key] for row in rows])
    data.save(path)
    return path


def _ports_unchanged(
    hard: UnifiedPolyBallBuild,
    candidate: UnifiedPolyBallBuild,
) -> bool:
    hard_rows = hard.metadata["port_clip_rows"]
    candidate_rows = candidate.metadata["port_clip_rows"]
    if len(hard_rows) != len(candidate_rows):
        return False
    for first, second in zip(hard_rows, candidate_rows):
        if first["port_id"] != second["port_id"]:
            return False
        if not np.allclose(first["plane_center_um"], second["plane_center_um"], atol=0.0):
            return False
        if not np.allclose(
            first["plane_outward_normal"], second["plane_outward_normal"], atol=0.0
        ):
            return False
    return True


def _per_junction_normal_metrics(
    mesh: trimesh.Trimesh, specs: tuple[Any, ...]
) -> list[dict[str, Any]]:
    return [
        {
            "junction_node_id": spec.junction_node_id,
            **local_normal_metrics(mesh, (spec,)),
        }
        for spec in specs
    ]


def _markdown(report: dict[str, Any]) -> str:
    answers = report["final_answers"]
    lines = [
        "# v8 local smooth-union 最终诊断",
        "",
        f"最终结论：**{report['decision']}**",
        "",
    ]
    for index in range(1, 11):
        lines.append(f"{index}. {answers[str(index)]}")
    lines.extend(
        (
            "",
            "所有 S0–S3 与 k/r 对照图均使用同一相机、flat shading、wireframe；",
            "smooth candidate 的 FlyingEdges、第一次 Newton 和第二次 Newton 始终使用同一 Phi_v8。",
            "端口中心线延伸、端口平面和端口法向未修改。",
            "",
        )
    )
    return "\n".join(lines)


def run_v8_refinement(
    roi: ROIRecord,
    branches: list[BranchGeometry],
    ports: list[PortGeometry],
    config: CFDLumenConfig,
    v6_result: V6RefinementResult,
    v7_result: V7RefinementResult,
    output_root: Path,
) -> V8RefinementResult:
    output_root = Path(output_root)
    diagnostics = output_root / "diagnostics"
    figures = output_root / "figures"
    geometry = output_root / "geometry"
    variants_root = output_root / "variants"
    report_root = output_root / "report"
    stage_root = output_root / "pipeline_stages"
    for folder in (
        diagnostics,
        figures,
        geometry,
        variants_root,
        report_root,
        stage_root,
    ):
        folder.mkdir(parents=True, exist_ok=True)

    started = time.perf_counter()
    hard = v7_result.candidate_build
    specs = build_junction_blend_specs(roi, v6_result.details, config)
    stage_rows, earliest_stage = stage_localization_rows(
        hard.stage_meshes, hard.model, specs, config
    )
    write_json(diagnostics / "pipeline_stage_localization.json", {
        "stages": stage_rows,
        "earliest_sawtooth_stage": earliest_stage,
    })
    write_csv(diagnostics / "pipeline_stage_localization.csv", stage_rows)
    stage_paths: list[Path] = []
    for stage, mesh in hard.stage_meshes.items():
        path = stage_root / f"{stage}.vtp"
        mesh_polydata(mesh).save(path)
        stage_paths.append(path)

    worst_spec = max(
        specs,
        key=lambda spec: (
            local_normal_metrics(hard.stage_meshes["S3_final_projected_before_port_clip"], (spec,))[
                "normal_jump_deg"
            ]["p99"]
            or -np.inf
        ),
    )
    stage_figure_paths, camera, _ = generate_stage_localization_figures(
        hard.stage_meshes, worst_spec, figures / "pipeline_stages"
    )

    ownership_path = diagnostics / "polyball_ownership.vtp"
    ownership, switch_rows, ownership_arrays = evaluate_ownership_switching(
        hard.stage_meshes["S3_final_projected_before_port_clip"],
        hard.model,
        specs,
        config,
        ownership_path=ownership_path,
    )
    switch_path = _write_edge_polydata(
        hard.stage_meshes["S3_final_projected_before_port_clip"],
        ownership_arrays["switch_edges"],
        switch_rows,
        diagnostics / "ownership_switching_boundary.vtp",
    )
    write_json(diagnostics / "polyball_ownership_summary.json", ownership)
    write_csv(diagnostics / "ownership_gradient_jump.csv", switch_rows)
    ownership_figures = ownership_and_overlay_figures(
        hard.stage_meshes["S3_final_projected_before_port_clip"],
        ownership_arrays["winner"],
        ownership_arrays["switch_edges"],
        ownership_arrays["defect_edges"],
        worst_spec,
        figures / "ownership",
        camera,
    )

    hard_surface_qc = evaluate_surface_qc(
        hard.mesh, hard.patch, roi, branches, config
    )
    hard_topology = evaluate_unified_topology(
        hard.mesh, roi, v6_result.details, config
    )
    _, hard_radius = evaluate_radius_fidelity(
        hard.mesh, branches, roi, config
    )
    hard_collar = collar_radius_rows(hard.mesh, v6_result.details, branches)
    hard_area = junction_area_profile_rows(
        hard.mesh, v6_result.details, branches
    )
    hard_volume = junction_local_volumes(hard.model, specs, config)
    hard_normal = local_normal_metrics(hard.wall_mesh_before_clip, specs)
    hard_gradient_on_switch = field_gradient_angles_on_edges(
        hard.stage_meshes["S3_final_projected_before_port_clip"],
        hard.field_model,
        ownership_arrays["switch_edges"],
    )

    candidates: dict[float, UnifiedPolyBallBuild] = {}
    raw_metrics: dict[float, dict[str, Any]] = {}
    prepared_started = time.perf_counter()
    prepared = prepare_polyball_raster(
        branches,
        ports,
        config,
        v6_details=v6_result.details,
        cells_across_min_diameter=config.v7.cells_across_min_diameter,
    )
    prepared_runtime = time.perf_counter() - prepared_started
    try:
        for ratio in config.v8.k_radius_ratios:
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
            )
            candidates[float(ratio)] = build
            surface = evaluate_surface_qc(build.mesh, build.patch, roi, branches, config)
            topology = evaluate_unified_topology(
                build.mesh, roi, v6_result.details, config
            )
            _, radius = evaluate_radius_fidelity(
                build.mesh, branches, roi, config
            )
            collar = collar_radius_rows(build.mesh, v6_result.details, branches)
            area = junction_area_profile_rows(build.mesh, v6_result.details, branches)
            volume = junction_local_volumes(build.field_model, specs, config)
            volume_comparison = compare_volume_rows(hard_volume, volume)
            hydraulic = local_hydraulic_resistance_from_area_profiles(hard_area, area)
            gradient_on_switch = field_gradient_angles_on_edges(
                hard.stage_meshes["S3_final_projected_before_port_clip"],
                build.field_model,
                ownership_arrays["switch_edges"],
            )
            raw_metrics[float(ratio)] = {
                "k_radius_ratio": float(ratio),
                "runtime_s": time.perf_counter() - candidate_started,
                "shared_hard_field_preparation_runtime_s": prepared_runtime,
                "surface_qc": surface,
                "topology": topology,
                "radius_fidelity": radius,
                "collar_radius": collar,
                "branch_local_area_profile": area,
                "junction_volume": volume_comparison,
                "hydraulic_resistance": hydraulic,
                "normal": local_normal_metrics(build.wall_mesh_before_clip, specs),
                "normal_per_junction": _per_junction_normal_metrics(
                    build.wall_mesh_before_clip, specs
                ),
                "gradient_angle_across_hard_switch_edges": gradient_on_switch,
                "ports_unchanged": _ports_unchanged(hard, build),
                "same_phi_v8_all_stages": build.metadata[
                    "same_field_for_flying_edges_and_both_newton_projections"
                ],
            }
    finally:
        release_prepared_polyball_raster(prepared)

    comparison_meshes = {
        "v7_hard_min": hard.wall_mesh_before_clip,
        "v8_k0p10r": candidates[0.10].wall_mesh_before_clip,
        "v8_k0p20r": candidates[0.20].wall_mesh_before_clip,
        "v8_k0p30r": candidates[0.30].wall_mesh_before_clip,
    }
    sensitivity_figures, sensitivity_camera, local_meshes = (
        generate_k_sensitivity_figures(
            comparison_meshes, worst_spec, figures / "k_sensitivity", camera=camera
        )
    )
    shared_local_meshes = tuple(local_meshes.values())
    hard_silhouette_rows = silhouette_rows(
        local_meshes["v7_hard_min"],
        version="v7_hard_min",
        large_corner_deg=config.v6.silhouette_large_corner_deg,
        comparison_meshes=shared_local_meshes,
    )
    hard_silhouette = _silhouette_summary(hard_silhouette_rows)

    sensitivity_rows: list[dict[str, Any]] = []
    all_silhouette_rows = list(hard_silhouette_rows)
    for ratio in config.v8.k_radius_ratios:
        label = f"v8_k0p{int(round(100 * ratio)):02d}r"
        silhouette = silhouette_rows(
            local_meshes[label],
            version=label,
            large_corner_deg=config.v6.silhouette_large_corner_deg,
            comparison_meshes=shared_local_meshes,
        )
        all_silhouette_rows.extend(silhouette)
        silhouette_summary = _silhouette_summary(silhouette)
        metric = raw_metrics[float(ratio)]
        metric["silhouette"] = silhouette_summary
        roughness_reduction = _relative_reduction(
            hard_silhouette["curvature_variation_mean"],
            silhouette_summary["curvature_variation_mean"],
        )
        hard_normal_p99 = hard_normal["normal_jump_deg"]["p99"]
        candidate_normal_p99 = metric["normal"]["normal_jump_deg"]["p99"]
        hard_gradient_p99 = hard_gradient_on_switch["p99"]
        candidate_gradient_p99 = metric[
            "gradient_angle_across_hard_switch_edges"
        ]["p99"]
        collar_increase = (
            (_maximum(metric["collar_radius"], "absolute_radius_relative_error") or 0.0)
            - (_maximum(hard_collar, "absolute_radius_relative_error") or 0.0)
        )
        volume_increase = metric["junction_volume"]["maximum_relative_increase"]
        hydraulic_error = metric["hydraulic_resistance"][
            "max_absolute_relative_error"
        ]
        silhouette_eliminated = bool(
            roughness_reduction is not None
            and roughness_reduction
            >= config.v8.minimum_silhouette_roughness_reduction_fraction
            and silhouette_summary["large_corner_count"]
            < hard_silhouette["large_corner_count"]
        )
        checks = {
            "topology_all_pass": metric["topology"]["status"] == "PASS",
            "ports_unchanged_and_pass": metric["ports_unchanged"]
            and candidates[float(ratio)].patch.all_ports_pass,
            "same_phi_v8_for_extraction_and_both_projections": metric[
                "same_phi_v8_all_stages"
            ],
            "post_projection_within_tolerance": candidates[float(ratio)].metadata[
                "projection_after_remesh"
            ]["post_projection_phi_abs_p95_um"]
            <= config.v7.projection_tolerance_um,
            "silhouette_sawtooth_eliminated": silhouette_eliminated,
            "normal_p99_improved": candidate_normal_p99 is not None
            and hard_normal_p99 is not None
            and candidate_normal_p99 < hard_normal_p99,
            "field_gradient_p99_improved": candidate_gradient_p99 is not None
            and hard_gradient_p99 is not None
            and candidate_gradient_p99 < hard_gradient_p99,
            "junction_volume_not_artificially_bulged": volume_increase is not None
            and volume_increase
            <= config.v8.maximum_junction_volume_increase_fraction,
            "radius_p95_acceptable": metric["radius_fidelity"][
                "p95_absolute_relative_error"
            ]
            <= config.v8.maximum_radius_p95_error,
            "collar_radius_error_not_significantly_worse": collar_increase
            <= config.v8.maximum_collar_radius_error_increase,
            "hydraulic_resistance_not_significantly_changed": hydraulic_error
            is not None
            and hydraulic_error <= config.v8.maximum_hydraulic_resistance_error,
        }
        metric["comparison_to_v7"] = {
            "silhouette_roughness_reduction_fraction": roughness_reduction,
            "normal_p99_change_deg": (
                candidate_normal_p99 - hard_normal_p99
                if candidate_normal_p99 is not None and hard_normal_p99 is not None
                else None
            ),
            "gradient_p99_change_deg": (
                candidate_gradient_p99 - hard_gradient_p99
                if candidate_gradient_p99 is not None and hard_gradient_p99 is not None
                else None
            ),
            "maximum_junction_volume_increase_fraction": volume_increase,
            "maximum_collar_error_increase": collar_increase,
            "maximum_hydraulic_resistance_error": hydraulic_error,
        }
        metric["acceptance"] = {"checks": checks, "pass": all(checks.values())}
        sensitivity_rows.append(
            {
                "k_radius_ratio": float(ratio),
                "normal_p95_deg": metric["normal"]["normal_jump_deg"]["p95"],
                "normal_p99_deg": candidate_normal_p99,
                "normal_max_deg": metric["normal"]["normal_jump_deg"]["max"],
                "gradient_jump_p99_deg": candidate_gradient_p99,
                "silhouette_curvature_variation_mean": silhouette_summary[
                    "curvature_variation_mean"
                ],
                "silhouette_large_corner_count": silhouette_summary[
                    "large_corner_count"
                ],
                "junction_volume_max_relative_increase": volume_increase,
                "radius_p95_absolute_relative_error": metric["radius_fidelity"][
                    "p95_absolute_relative_error"
                ],
                "collar_max_absolute_relative_error": _maximum(
                    metric["collar_radius"], "absolute_radius_relative_error"
                ),
                "hydraulic_resistance_max_absolute_relative_error": hydraulic_error,
                "topology_status": metric["topology"]["status"],
                "acceptance_pass": all(checks.values()),
            }
        )

    eligible = [
        ratio
        for ratio in config.v8.k_radius_ratios
        if raw_metrics[float(ratio)]["acceptance"]["pass"]
    ]
    remeshing_root_cause = (
        not ownership["hard_min_switch_crease_confirmed"]
        and earliest_stage
        in {"S2_pyacvd_before_second_projection", "S3_final_projected_before_port_clip"}
    )
    if ownership["hard_min_switch_crease_confirmed"] and eligible:
        selected_ratio = min(
            eligible,
            key=lambda ratio: (
                raw_metrics[float(ratio)]["silhouette"][
                    "curvature_variation_mean"
                ],
                raw_metrics[float(ratio)]["normal"]["normal_jump_deg"]["p99"],
                abs(
                    raw_metrics[float(ratio)]["junction_volume"][
                        "maximum_relative_increase"
                    ]
                ),
            ),
        )
        decision = "ADOPT_V8_LOCAL_SMOOTH_UNION"
        selected = candidates[float(selected_ratio)]
    elif remeshing_root_cause:
        selected_ratio = None
        decision = "ROOT_CAUSE_IS_REMESHING_NOT_FIELD"
        selected = hard
    else:
        selected_ratio = None
        decision = "KEEP_V7_HARD_MIN"
        selected = hard

    variant_geometry_paths: list[Path] = []
    for ratio, build in candidates.items():
        label = f"k_{ratio:.2f}r".replace(".", "p")
        variant_root = variants_root / label
        (variant_root / "geometry").mkdir(parents=True, exist_ok=True)
        (variant_root / "ports").mkdir(parents=True, exist_ok=True)
        variant_geometry_paths.extend(
            write_geometry_exports(
                build.mesh,
                build.patch,
                ports,
                {
                    "geometry": variant_root / "geometry",
                    "ports": variant_root / "ports",
                },
                face_region=np.zeros(len(build.mesh.faces), dtype=np.uint8),
            )
        )
        write_json(variant_root / "metrics.json", raw_metrics[ratio])

    (output_root / "ports").mkdir(parents=True, exist_ok=True)
    geometry_paths = write_geometry_exports(
        selected.mesh,
        selected.patch,
        ports,
        {"geometry": geometry, "ports": output_root / "ports"},
        face_region=np.zeros(len(selected.mesh.faces), dtype=np.uint8),
    )
    geometry_paths.extend((ownership_path, switch_path, *stage_paths, *variant_geometry_paths))
    geometry_paths.append(write_units(diagnostics / "units.json"))

    write_csv(diagnostics / "k_sensitivity.csv", sensitivity_rows)
    write_csv(diagnostics / "silhouette_v7_v8.csv", all_silhouette_rows)
    write_csv(diagnostics / "junction_area_profile_v7.csv", hard_area)
    write_csv(diagnostics / "collar_radius_v7.csv", hard_collar)
    write_json(diagnostics / "v7_hard_min_baseline.json", {
        "surface_qc": hard_surface_qc,
        "topology": hard_topology,
        "radius_fidelity": hard_radius,
        "normal": hard_normal,
        "gradient_angle_across_switch_edges": hard_gradient_on_switch,
        "silhouette": hard_silhouette,
        "junction_volume": hard_volume,
    })
    write_json(diagnostics / "k_sensitivity_full.json", {
        f"{ratio:.2f}": metric for ratio, metric in raw_metrics.items()
    })
    write_json(diagnostics / "same_camera.json", {
        "junction_node_id": worst_spec.junction_node_id,
        "camera": camera_report(sensitivity_camera),
        "stage_views": ["flat", "wireframe"],
        "sensitivity_views": ["flat", "wireframe", "silhouette"],
    })

    selected_metric = raw_metrics.get(float(selected_ratio)) if selected_ratio else None
    selected_volume = (
        selected_metric["junction_volume"]["maximum_relative_increase"]
        if selected_metric
        else 0.0
    )
    selected_radius = (
        selected_metric["radius_fidelity"]["p95_absolute_relative_error"]
        if selected_metric
        else hard_radius["p95_absolute_relative_error"]
    )
    selected_hydraulic = (
        selected_metric["hydraulic_resistance"]["max_absolute_relative_error"]
        if selected_metric
        else 0.0
    )
    selected_silhouette = (
        selected_metric["acceptance"]["checks"]["silhouette_sawtooth_eliminated"]
        if selected_metric
        else False
    )
    final_answers = {
        "1": f"锯齿最早出现于 {earliest_stage or '未能自动确认的阶段'}。",
        "2": (
            "是；缺陷与 winner-branch switching boundary 的近邻覆盖率为 "
            f"defect→switch {ownership['defect_near_switch_fraction']:.3%}，"
            f"switch→defect {ownership['switch_near_defect_fraction']:.3%}。"
            if ownership["hard_min_switch_crease_confirmed"]
            else "否；未达到 ownership switching overlap 的确认条件。"
        ),
        "3": (
            "hard-min gradient jump P99 = "
            f"{ownership['gradient_jump_angle_deg']['p99']:.6g}°。"
        ),
        "4": "否。" if not remeshing_root_cause else "是。",
        "5": (
            "是。" if selected_silhouette else "否；三档参数均未通过 silhouette 消除条件。"
        ),
        "6": (
            f"k/r = {selected_ratio:.2f}。" if selected_ratio is not None else "无可采纳参数。"
        ),
        "7": f"所选结果最大 junction local volume 增量为 {selected_volume:.3%}。",
        "8": (
            f"radius P95 = {selected_radius:.3%}，"
            + ("可接受。" if selected_radius <= config.v8.maximum_radius_p95_error else "不可接受。")
        ),
        "9": (
            f"最大 branch-local hydraulic resistance 相对变化为 {selected_hydraulic:.3%}，"
            + (
                "未显著改变。"
                if selected_hydraulic <= config.v8.maximum_hydraulic_resistance_error
                else "发生显著改变。"
            )
        ),
        "10": decision,
    }
    figure_paths = [
        *stage_figure_paths,
        *ownership_figures,
        *sensitivity_figures,
    ]
    report = {
        "protocol": "v8 local smooth-union root-cause protocol",
        "roi_id": roi.roi_id,
        "status": "PASS" if selected.patch.all_ports_pass else "FAIL",
        "decision": decision,
        "selected_k_radius_ratio": selected_ratio,
        "root_cause": (
            "HARD_MIN_POLYBALL_SWITCH_CREASE"
            if ownership["hard_min_switch_crease_confirmed"]
            else "REMESHING"
            if remeshing_root_cause
            else "NOT_CONFIRMED"
        ),
        "pipeline_stage_localization": {
            "rows": stage_rows,
            "earliest_sawtooth_stage": earliest_stage,
        },
        "ownership": ownership,
        "hard_min_gradient_jump": ownership["gradient_jump_angle_deg"],
        "hard_min_surface_dihedral_at_switch": ownership[
            "surface_dihedral_at_switch_deg"
        ],
        "selected_surface_qc": (
            selected_metric["surface_qc"] if selected_metric else hard_surface_qc
        ),
        "selected_topology": (
            selected_metric["topology"] if selected_metric else hard_topology
        ),
        "selected_radius_fidelity": (
            selected_metric["radius_fidelity"] if selected_metric else hard_radius
        ),
        "selected_collar_radius": (
            selected_metric["collar_radius"] if selected_metric else hard_collar
        ),
        "v7_hard_min": {
            "surface_qc": hard_surface_qc,
            "topology": hard_topology,
            "radius_fidelity": hard_radius,
            "collar_radius": hard_collar,
            "normal": hard_normal,
            "silhouette": hard_silhouette,
            "junction_volume": hard_volume,
        },
        "k_sensitivity": {
            f"{ratio:.2f}": metric for ratio, metric in raw_metrics.items()
        },
        "final_answers": final_answers,
        "geometry_paths": [str(path) for path in geometry_paths],
        "figure_paths": [str(path) for path in figure_paths],
        "runtime_s": time.perf_counter() - started,
    }
    write_json(report_root / "v8_final_report.json", report)
    (report_root / "v8_final_report.md").write_text(
        _markdown(report), encoding="utf-8"
    )
    return V8RefinementResult(
        mesh=selected.mesh,
        patch=selected.patch,
        selected_build=selected,
        report=report,
        decision=decision,
        geometry_paths=geometry_paths,
        figure_paths=figure_paths,
        candidate_builds=candidates,
    )
