"""End-to-end real connected-ROI sampling, comparison, validation, and export."""

from __future__ import annotations

import json
import logging
import time
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np

from ..reporting.logger import configure_logging
from ..rodent_vasculature.swc_io import load_normalized_swc
from .clustering import deterministic_kmeans, exploratory_cluster_scan
from .feature_scaling import apply_group_weights, robust_scale
from .representative_selection import select_representatives
from .roi_extraction import extract_candidate_rois, generate_anchor_ids, global_model_from_swc
from .roi_features import build_feature_matrix, populate_roi_features
from .sampling_config import SamplingConfig
from .sampling_io import (
    SamplingOutputLayout,
    create_sampling_layout,
    experiment_report,
    write_candidate_tables,
    write_experiment,
    write_global_edge_manifest,
    write_json,
    write_roi_library,
    write_rows,
    write_sampling_config,
)
from .sampling_types import (
    GlobalVascularModel,
    ROIRecord,
    SamplingExperiment,
    SamplingRunResult,
)
from .sampling_validation import validate_sampling
from .sampling_visualization import (
    cluster_size_figure,
    feature_mode_comparison_figure,
    global_roi_overview,
    pca_feature_space_figure,
    radius_distribution_figure,
    roi_preview,
    scalar_distribution_figure,
    silhouette_scan_figure,
)


def load_models_from_rodent_run(run_root: Path) -> list[GlobalVascularModel]:
    """Thin adapter from saved rodent preprocessing to the sampling core."""

    models: list[GlobalVascularModel] = []
    sample_root = Path(run_root) / "samples"
    for sample_dir in sorted(path for path in sample_root.iterdir() if path.is_dir()):
        manifest_path = sample_dir / "preprocess_manifest.json"
        if not manifest_path.is_file():
            continue
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        record = manifest["record"]
        spacing = tuple(float(value) for value in manifest["spacing_xyz_um"])
        shape_zyx = tuple(int(value) for value in manifest["image_metadata"]["shape_zyx"])
        swc = load_normalized_swc(
            Path(manifest["normalized_swc_path"]),
            source_path=Path(record["swc_path"]),
            spacing_xyz_um=spacing,
            volume_shape_zyx=shape_zyx,
        )
        maximum = (np.asarray(shape_zyx[::-1], dtype=float) - 1.0) * np.asarray(spacing)
        bounds = (0.0, float(maximum[0]), 0.0, float(maximum[1]), 0.0, float(maximum[2]))
        models.append(
            global_model_from_swc(
                swc,
                source_model_id=str(record["sample_id"]),
                source_mouse_id=str(record.get("parent_group_id") or record["sample_id"]),
                model_bounds_xyz_um=bounds,
            )
        )
    if not models:
        raise FileNotFoundError(f"No completed rodent samples found in {run_root}")
    return models


def _run_experiment(
    rois: list[ROIRecord],
    config: SamplingConfig,
    *,
    feature_mode: str,
    stage_timings: dict[str, float] | None = None,
) -> SamplingExperiment:
    started = time.perf_counter()
    raw, feature_names, groups = build_feature_matrix(rois, config, feature_mode=feature_mode)
    scaled, scaler = robust_scale(raw, feature_names)
    weighted = apply_group_weights(
        scaled,
        groups,
        radius_weight=config.radius_feature_weight,
        structure_weight=config.structure_feature_weight,
    )
    if stage_timings is not None:
        stage_timings["scaling"] = time.perf_counter() - started
    started = time.perf_counter()
    clustering = deterministic_kmeans(
        weighted,
        n_clusters=config.n_clusters,
        feature_names=feature_names,
        seed=config.seed,
        max_iter=config.kmeans_max_iter,
    )
    if stage_timings is not None:
        stage_timings["clustering"] = time.perf_counter() - started
    started = time.perf_counter()
    selection = select_representatives(rois, clustering, config)
    if stage_timings is not None:
        stage_timings["representative_selection"] = time.perf_counter() - started
    started = time.perf_counter()
    validation = validate_sampling(rois, clustering, selection)
    if stage_timings is not None:
        stage_timings["validation"] = time.perf_counter() - started
    return SamplingExperiment(
        feature_mode,
        feature_names,
        scaler,
        clustering,
        selection,
        validation,
    )


def _assign_primary_results(rois: list[ROIRecord], experiment: SamplingExperiment) -> None:
    selected_rank = {index: rank for rank, index in enumerate(experiment.selection.selected_indices, 1)}
    for index, roi in enumerate(rois):
        roi.cluster_id = int(experiment.clustering.assignments[index])
        roi.distance_to_cluster_center = float(experiment.clustering.distances_to_center[index])
        roi.is_representative = index in selected_rank
        roi.selection_rank = selected_rank.get(index, -1)


def _write_figures(
    layout: SamplingOutputLayout,
    models: list[GlobalVascularModel],
    rois: list[ROIRecord],
    experiment: SamplingExperiment,
    comparison: dict[str, SamplingExperiment],
    scan: list[dict[str, float | int | None]],
    config: SamplingConfig,
) -> list[Path]:
    selected = [rois[index] for index in experiment.selection.selected_indices]
    paths = [
        global_roi_overview(models, rois, layout.figures / "global_candidate_overview.png", selected_only=False),
        global_roi_overview(models, selected, layout.figures / "global_selected_overview.png", selected_only=True),
        radius_distribution_figure(rois, selected, layout.figures / "radius_distribution_global_vs_selected.png"),
        scalar_distribution_figure(rois, selected, "branch_count", "Branch-count distribution", layout.figures / "branch_count_global_vs_selected.png"),
        scalar_distribution_figure(rois, selected, "bifurcation_count", "Bifurcation-count distribution", layout.figures / "bifurcation_count_global_vs_selected.png"),
        scalar_distribution_figure(rois, selected, "total_vessel_length_um", "Vessel-length distribution", layout.figures / "vessel_length_global_vs_selected.png"),
        scalar_distribution_figure(rois, selected, "cycle_rank", "Cycle-rank distribution", layout.figures / "cycle_rank_global_vs_selected.png"),
        cluster_size_figure(experiment, layout.figures / "cluster_size_distribution.png"),
        pca_feature_space_figure(experiment, layout.figures / "pca_feature_space.png"),
        silhouette_scan_figure(scan, layout.figures / "silhouette_scan.png"),
    ]
    if len(comparison) >= 2:
        paths.append(
            feature_mode_comparison_figure(
                comparison,
                layout.figures / "radius_only_vs_radius_plus_structure.png",
            )
        )
    for rank, roi in enumerate(selected[: config.max_roi_previews]):
        paths.append(
            roi_preview(
                roi,
                layout.roi_previews / f"cluster_{roi.cluster_id:03d}_rep_{rank:02d}.png",
            )
        )
    return paths


def run_sampling_pipeline(
    models: list[GlobalVascularModel],
    config: SamplingConfig,
    *,
    verbose: bool = False,
) -> SamplingRunResult:
    """Run all Sampling Project phases once over one or more real source models."""

    config.validate()
    total_start = time.perf_counter()
    layout = create_sampling_layout(config)
    logger = configure_logging(layout.log_file, verbose=verbose)
    write_sampling_config(layout, config)
    write_json(
        layout.report / "repository_inspection.json",
        {
            "global_model_source": "existing rodent_vasculature SWCData arrays",
            "global_node_id": "SWC node_id, preserved without renumbering",
            "global_edge_id": (
                "deterministic parent-to-current edge index in SWC row order; disambiguated by "
                "source_model_id and exported in manifests/global_edges.csv"
            ),
            "parent_id": "SWC parent_ids",
            "radius_storage": "SWC radius_raw_um per global node; interpolated only at exact cuts",
            "existing_graphs": ["source_graph", "junction_graph", "branch_graph"],
            "previous_roi_definition": (
                "connected descendant branch subtree for Figure 1(b)-style visualization"
            ),
            "sampling_roi_definition": (
                "axis-aligned physical box clipping followed by anchor-containing connected component"
            ),
            "spatial_query": "cKDTree edge-midpoint query plus exact edge-AABB and segment-box tests",
            "ui_integration": (
                "saved candidate/selected manifests and ROI NPZ -> existing PyVista display; "
                "no sampling or clustering in callbacks"
            ),
        },
    )
    timings: dict[str, float] = {}
    logger.info("Sampling scope: full real vascular network -> representative real ROI library")
    logger.info("Input models: %s", ", ".join(model.source_model_id for model in models))
    logger.info(
        "Global node count=%d; global edge count=%d",
        sum(model.node_count for model in models),
        sum(model.edge_count for model in models),
    )
    logger.info(
        "ROI size um=%s; anchor mode=%s; minimum anchor distance um=%.3f; seed=%d",
        config.roi_size_um,
        config.anchor_mode,
        config.min_anchor_distance_um,
        config.seed,
    )

    candidates: list[ROIRecord] = []
    anchor_count = 0
    rejected: Counter[str] = Counter(
        {
            "empty": 0,
            "disconnected": 0,
            "too_few_branches": 0,
            "invalid_radius": 0,
            "invalid_graph": 0,
            "excessive_overlap": 0,
        }
    )
    anchor_runtime = 0.0
    extraction_runtime = 0.0
    for model in models:
        stage_start = time.perf_counter()
        anchors = generate_anchor_ids(model, config)
        anchor_runtime += time.perf_counter() - stage_start
        stage_start = time.perf_counter()
        batch = extract_candidate_rois(model, config, anchor_ids=anchors)
        extraction_runtime += time.perf_counter() - stage_start
        anchor_count += len(batch.anchor_ids)
        candidates.extend(batch.candidates)
        rejected.update(batch.rejected_reasons)
    timings["anchor_generation"] = anchor_runtime
    timings["roi_extraction"] = extraction_runtime

    stage_start = time.perf_counter()
    valid: list[ROIRecord] = []
    for roi in candidates:
        populate_roi_features(roi)
        if roi.branch_count < config.min_branch_count:
            rejected["too_few_branches"] += 1
            continue
        valid.append(roi)
    candidates = valid
    timings["feature_extraction"] = time.perf_counter() - stage_start
    if not candidates:
        raise ValueError(
            "Sampling produced no valid connected ROIs; lower --sampling-min-branches or "
            "--sampling-min-anchor-distance, or increase the ROI size"
        )
    logger.info(
        "Candidate anchors=%d; valid candidate ROIs=%d; rejected=%d; reasons=%s",
        anchor_count,
        len(candidates),
        sum(rejected.values()),
        dict(sorted(rejected.items())),
    )

    primary = _run_experiment(
        candidates,
        config,
        feature_mode=config.feature_mode,
        stage_timings=timings,
    )
    _assign_primary_results(candidates, primary)
    logger.info("Feature mode=%s; names=%s; dimension=%d", primary.feature_mode, primary.feature_names, len(primary.feature_names))
    logger.info(
        "Clustering method=%s; K=%d; sizes=%s; silhouette=%s",
        primary.clustering.method,
        primary.clustering.n_clusters,
        primary.clustering.cluster_sizes,
        primary.clustering.silhouette_score,
    )
    logger.info(
        "Selection mode=%s (%s); selected=%d/%d; overlap rejections=%d",
        primary.selection.selection_mode,
        primary.selection.scientific_label,
        len(primary.selection.selected_indices),
        primary.selection.requested_count,
        primary.selection.overlap_rejection_count,
    )

    stage_start = time.perf_counter()
    comparison: dict[str, SamplingExperiment] = {config.feature_mode: primary}
    if config.compare_feature_modes:
        for mode in ("radius_only", "radius_plus_structure"):
            if mode not in comparison:
                comparison[mode] = _run_experiment(candidates, config, feature_mode=mode)
    timings["feature_mode_comparison"] = time.perf_counter() - stage_start

    primary_raw, primary_names, primary_groups = build_feature_matrix(
        candidates, config, feature_mode=config.feature_mode
    )
    primary_scaled, _ = robust_scale(primary_raw, primary_names)
    primary_weighted = apply_group_weights(
        primary_scaled,
        primary_groups,
        radius_weight=config.radius_feature_weight,
        structure_weight=config.structure_feature_weight,
    )
    scan = exploratory_cluster_scan(
        primary_weighted,
        feature_names=primary_names,
        candidate_k=config.exploratory_k,
        seed=config.seed,
        max_iter=config.kmeans_max_iter,
    )

    stage_start = time.perf_counter()
    write_global_edge_manifest(layout, models)
    write_candidate_tables(layout, candidates)
    write_roi_library(layout, candidates)
    write_experiment(layout, primary, candidates)
    write_rows(layout.clustering / "silhouette_scan.csv", scan)
    for mode, experiment in comparison.items():
        write_experiment(
            layout,
            experiment,
            candidates,
            subdirectory=layout.comparison / mode,
        )
    timings["export"] = time.perf_counter() - stage_start

    stage_start = time.perf_counter()
    figures = _write_figures(layout, models, candidates, primary, comparison, scan, config)
    timings["visualization"] = time.perf_counter() - stage_start
    timings["total"] = time.perf_counter() - total_start

    integrity = primary.validation["roi_integrity"]
    selected_count = len(primary.selection.selected_indices)
    status = "PASS" if (
        selected_count > 0
        and integrity["all_connected"]
        and integrity["all_globally_traceable"]
        and integrity["all_boundary_semantics_valid"]
    ) else "FAIL"
    summary: dict[str, Any] = {
        "run_id": layout.run_root.name,
        "status": status,
        "seed": config.seed,
        "input_models": [model.source_model_id for model in models],
        "global_node_count": sum(model.node_count for model in models),
        "global_edge_count": sum(model.edge_count for model in models),
        "roi_size_um": list(config.roi_size_um),
        "anchor_mode": config.anchor_mode,
        "minimum_anchor_distance_um": config.min_anchor_distance_um,
        "candidate_anchor_count": anchor_count,
        "candidate_count": len(candidates) + sum(rejected.values()),
        "valid_candidate_count": len(candidates),
        "rejected_candidate_count": sum(rejected.values()),
        "reject_reasons": dict(sorted(rejected.items())),
        "feature_mode": config.feature_mode,
        "features": list(primary.feature_names),
        "cluster_method": primary.clustering.method,
        "n_clusters": primary.clustering.n_clusters,
        "selection_mode": primary.selection.selection_mode,
        "selection_scientific_label": primary.selection.scientific_label,
        "selected_count": selected_count,
        "validation": primary.validation,
        "feature_mode_comparison": {
            mode: experiment_report(experiment) for mode, experiment in comparison.items()
        },
        "exploratory_k": scan,
        "runtime_seconds": timings,
        "output_directory": str(layout.run_root),
        "physiological_inlet_outlet_inference_performed": False,
        "flow_solver_performed": False,
    }
    write_json(layout.summary_file, summary)
    write_json(layout.report / "run_status.json", {"status": status, "run_id": layout.run_root.name})
    logger.info("Validation metrics=%s", primary.validation)
    logger.info("Runtime seconds=%s", timings)
    logger.info("Output directory=%s", layout.run_root)
    logger.info("Sampling run finished: %s", status)
    return SamplingRunResult(
        run_root=layout.run_root,
        status=status,
        candidates=candidates,
        primary_experiment=primary,
        comparison_experiments=comparison,
        summary_path=layout.summary_file,
        log_path=layout.log_file,
        figure_paths=figures,
    )


def run_sampling_from_rodent_run(
    rodent_run_root: Path,
    config: SamplingConfig,
    *,
    verbose: bool = False,
) -> SamplingRunResult:
    load_start = time.perf_counter()
    models = load_models_from_rodent_run(rodent_run_root)
    model_loading_seconds = time.perf_counter() - load_start
    result = run_sampling_pipeline(models, config, verbose=verbose)
    with result.summary_path.open(encoding="utf-8") as stream:
        summary = json.load(stream)
    summary["runtime_seconds"]["model_loading"] = model_loading_seconds
    summary["runtime_seconds"]["total_including_model_loading"] = (
        summary["runtime_seconds"]["total"] + model_loading_seconds
    )
    write_json(result.summary_path, summary)
    logging.getLogger("ulm_3d_vascular").info(
        "Model loading seconds=%.6f; total including model loading seconds=%.6f",
        model_loading_seconds,
        summary["runtime_seconds"]["total_including_model_loading"],
    )
    return result
