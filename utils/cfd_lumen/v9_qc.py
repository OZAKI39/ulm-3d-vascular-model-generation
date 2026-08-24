"""Combined morphology, fidelity, and acceptance metrics for v9."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import trimesh

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .hybrid_qc import collar_radius_rows, junction_area_profile_rows
from .segment_ownership_qc import (
    CROSS_BRANCH_SWITCH,
    SAME_BRANCH_ADJACENT_SEGMENT_SWITCH,
    SAME_BRANCH_NONADJACENT_SEGMENT_SWITCH,
    evaluate_segment_ownership,
)
from .surface_qc import evaluate_radius_fidelity, evaluate_surface_qc
from .types import BranchGeometry, HybridBuildDetails
from .unified_polyball import JunctionBlendSpec, UnifiedPolyBallBuild
from .v6_qc import silhouette_rows
from .v7_qc import evaluate_unified_topology
from .v8_qc import (
    compare_volume_rows,
    junction_local_volumes,
    local_hydraulic_resistance_from_area_profiles,
    local_normal_metrics,
)


def _summary(values: list[float]) -> dict[str, float | int | None]:
    data = np.asarray(values, dtype=float)
    if not len(data):
        return {"count": 0, "mean": None, "p95": None, "p99": None, "max": None}
    return {
        "count": int(len(data)),
        "mean": float(np.mean(data)),
        "p95": float(np.percentile(data, 95)),
        "p99": float(np.percentile(data, 99)),
        "max": float(np.max(data)),
    }


def _maximum(rows: list[dict[str, Any]], key: str) -> float | None:
    values = [float(row[key]) for row in rows if row.get(key) is not None]
    return max(values, default=None)


def silhouette_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    roughness = [
        float(row["silhouette_curvature_variation"])
        for row in rows
        if row.get("silhouette_curvature_variation") is not None
    ]
    turning = [
        float(row["silhouette_turning_p95_deg"])
        for row in rows
        if row.get("silhouette_turning_p95_deg") is not None
    ]
    return {
        "curvature_variation_mean": float(np.mean(roughness)) if roughness else None,
        "turning_p95_worst_deg": max(turning, default=None),
        "large_corner_count": int(sum(int(row["large_corner_count"]) for row in rows)),
    }


def evaluate_v9_method(
    label: str,
    build: UnifiedPolyBallBuild,
    roi: ROIRecord,
    source_branches: list[BranchGeometry],
    details: HybridBuildDetails,
    specs: tuple[JunctionBlendSpec, ...],
    config: CFDLumenConfig,
    *,
    baseline_area: list[dict[str, Any]],
    baseline_volume: list[dict[str, Any]],
    comparison_meshes: tuple[trimesh.Trimesh, ...],
    runtime_s: float,
    ownership_path: Path | None = None,
    silhouette_mesh: trimesh.Trimesh | None = None,
) -> tuple[dict[str, Any], list[dict[str, Any]], dict[str, np.ndarray]]:
    segment, switch_rows, arrays = evaluate_segment_ownership(
        build.wall_mesh_before_clip,
        build.model,
        specs,
        config,
        ownership_path=ownership_path,
        field_model=build.field_model,
    )
    surface = evaluate_surface_qc(build.mesh, build.patch, roi, source_branches, config)
    topology = evaluate_unified_topology(build.mesh, roi, details, config)
    _, radius = evaluate_radius_fidelity(build.mesh, source_branches, roi, config)
    collar = collar_radius_rows(build.mesh, details, source_branches)
    area = junction_area_profile_rows(build.mesh, details, source_branches)
    volume = compare_volume_rows(
        baseline_volume, junction_local_volumes(build.field_model, specs, config)
    )
    hydraulic = local_hydraulic_resistance_from_area_profiles(baseline_area, area)
    silhouette = silhouette_summary(
        silhouette_rows(
            silhouette_mesh or build.wall_mesh_before_clip,
            version=label,
            large_corner_deg=config.v6.silhouette_large_corner_deg,
            comparison_meshes=comparison_meshes,
        )
    )
    same_gradient = _summary(
        [
            float(row["gradient_jump_angle_deg"])
            for row in switch_rows
            if row["switch_type"]
            in {
                SAME_BRANCH_ADJACENT_SEGMENT_SWITCH,
                SAME_BRANCH_NONADJACENT_SEGMENT_SWITCH,
            }
        ]
    )
    cross_gradient = segment["types"][CROSS_BRANCH_SWITCH][
        "gradient_jump_angle_deg"
    ]
    metrics = {
        "method": label,
        "segment_ownership": segment,
        "cross_branch_gradient_jump_deg": cross_gradient,
        "same_branch_segment_gradient_jump_deg": same_gradient,
        "surface_normal_jump": local_normal_metrics(
            build.wall_mesh_before_clip, specs
        )["normal_jump_deg"],
        "silhouette": silhouette,
        "visible_sawtooth_count": segment["visible_sawtooth_count"],
        "radius_fidelity": radius,
        "collar_radius": collar,
        "hydraulic_resistance": hydraulic,
        "junction_volume": volume,
        "surface_qc": surface,
        "topology": topology,
        "triangle_count": int(len(build.mesh.faces)),
        "runtime_s": float(runtime_s),
        "maximum_collar_absolute_relative_error": _maximum(
            collar, "absolute_radius_relative_error"
        ),
    }
    return metrics, switch_rows, arrays


def final_acceptance(
    metrics: dict[str, Any],
    centerline_report: dict[str, Any],
    baseline_metrics: dict[str, Any],
    config: CFDLumenConfig,
    *,
    same_camera_figures_generated: bool,
) -> dict[str, Any]:
    same_p99 = metrics["same_branch_segment_gradient_jump_deg"]["p99"]
    cross_p99 = metrics["cross_branch_gradient_jump_deg"]["p99"]
    baseline_same_p99 = baseline_metrics["same_branch_segment_gradient_jump_deg"][
        "p99"
    ]
    same_reduction = (
        (baseline_same_p99 - same_p99) / baseline_same_p99
        if baseline_same_p99 is not None
        and same_p99 is not None
        and baseline_same_p99 > 0.0
        else None
    )
    hydraulic = metrics["hydraulic_resistance"]["max_absolute_relative_error"]
    volume = metrics["junction_volume"]["maximum_relative_increase"]
    type_rows = metrics["segment_ownership"]["types"]
    same_switch_count = sum(
        int(type_rows[name]["switch_edge_count"])
        for name in (
            SAME_BRANCH_ADJACENT_SEGMENT_SWITCH,
            SAME_BRANCH_NONADJACENT_SEGMENT_SWITCH,
        )
    )
    cross_switch_count = int(type_rows[CROSS_BRANCH_SWITCH]["switch_edge_count"])
    same_removed = same_switch_count == 0 or (
        same_p99 is not None
        and same_p99 <= config.v9.maximum_same_branch_gradient_p99_deg
        and same_reduction is not None
        and same_reduction >= 0.5
    )
    cross_smoothed = cross_switch_count == 0 or (
        cross_p99 is not None
        and cross_p99 <= config.v9.maximum_cross_branch_gradient_p99_deg
    )
    checks = {
        "same_camera_flat_wireframe_generated": same_camera_figures_generated,
        "visible_sawtooth_disappeared": metrics["visible_sawtooth_count"] == 0,
        "same_branch_segment_switch_defect_significantly_removed": same_removed,
        "cross_branch_switch_smoothed": cross_smoothed,
        "radius_p95_below_one_percent": metrics["radius_fidelity"][
            "p95_absolute_relative_error"
        ]
        < config.v9.maximum_radius_p95_error,
        "hydraulic_error_within_tolerance": hydraulic is not None
        and hydraulic <= config.v9.maximum_hydraulic_resistance_error,
        "junction_volume_within_tolerance": volume is not None
        and volume <= config.v9.maximum_junction_volume_increase_fraction,
        "topology_all_pass": metrics["topology"]["status"] == "PASS",
        "surface_qc_pass": metrics["surface_qc"]["status"] == "PASS",
        "source_centerline_fidelity_pass": centerline_report["status"] == "PASS",
        "junction_position_exact": centerline_report[
            "junction_position_error_um"
        ]
        == 0.0,
        "endpoint_position_exact": centerline_report[
            "endpoint_position_error_um"
        ]
        == 0.0,
        "port_position_exact": centerline_report["port_position_error_um"] == 0.0,
    }
    return {
        "checks": checks,
        "pass": all(checks.values()),
        "same_branch_gradient_p99_reduction_fraction": same_reduction,
    }


def comparison_row(metrics: dict[str, Any]) -> dict[str, Any]:
    return {
        "method": metrics["method"],
        "cross_branch_gradient_p99_deg": metrics["cross_branch_gradient_jump_deg"][
            "p99"
        ],
        "same_branch_segment_gradient_p99_deg": metrics[
            "same_branch_segment_gradient_jump_deg"
        ]["p99"],
        "surface_normal_p99_deg": metrics["surface_normal_jump"]["p99"],
        "silhouette_roughness": metrics["silhouette"][
            "curvature_variation_mean"
        ],
        "visible_sawtooth_count": metrics["visible_sawtooth_count"],
        "radius_p95": metrics["radius_fidelity"][
            "p95_absolute_relative_error"
        ],
        "hydraulic_error": metrics["hydraulic_resistance"][
            "max_absolute_relative_error"
        ],
        "volume_change": metrics["junction_volume"][
            "maximum_relative_increase"
        ],
        "triangle_count": metrics["triangle_count"],
        "runtime_s": metrics["runtime_s"],
        "topology_status": metrics["topology"]["status"],
    }
