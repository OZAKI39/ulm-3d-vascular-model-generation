"""End-to-end NNE2 Step 1-3 equivalent and directed hierarchy pipeline."""

from __future__ import annotations

import gc
import gzip
import hashlib
import json
import logging
import pickle
import platform
import sys
import time
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterator

import matplotlib
import networkx as nx
import nibabel as nib
import numpy as np
import scipy
import skimage

from ..io import write_csv, write_json
from ..reporting.logger import configure_logging
from .catalog import NNE2Record, load_nne2_catalog
from .centerline import (
    build_centerline_graph,
    build_graph_from_centerline,
    extract_centerline,
    graph_config_from_nne2,
)
from .centerline import NNE2CenterlineResult
from .config import NNE2Config
from .export import (
    export_hierarchy,
    export_undirected_graph,
    write_acceptance_html,
    write_stack_nifti,
)
from .hierarchy import build_directed_hierarchy
from .landmarks import match_measurement_anchors
from .manifest import load_source_manifests
from .layout import (
    NNE2OutputLayout,
    create_nne2_output_layout,
    stack_output_dir,
    tree_output_dir,
)
from .segmentation import SegmentationResult, segment_vessels
from .preprocess_pipeline import load_step2_stack, write_step1_step2_artifacts
from .graph_pipeline import export_and_validate_step3_graph
from .stack_io import StackMetadata, inspect_stack, load_stack
from .validation import (
    required_files_status,
    validate_hierarchy,
    validate_preprocess_result,
    validate_stack_result,
)
from .visualization import (
    render_anchor_registration,
    render_component_decisions,
    render_hierarchy,
    render_segmentation_overview,
    render_step1_component_cleanup,
    render_step2_radius_and_centerline,
    render_step3_diagnostics,
    render_tree_map,
    render_undirected_graph,
)


@dataclass(frozen=True, slots=True)
class NNE2PipelineRun:
    run_root: Path
    html_report: Path
    status: str
    summary: dict[str, Any]


@contextmanager
def _timed_stage(
    name: str, logger: logging.Logger, timings: dict[str, float]
) -> Iterator[None]:
    logger.info("Starting: %s", name)
    started = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - started
        timings[name] = timings.get(name, 0.0) + elapsed
        logger.info("Finished: %s (%.3f s)", name, elapsed)


def _status(
    layout: NNE2OutputLayout,
    state: str,
    message: str,
    started_at: str,
    timings: dict[str, float],
    summary: dict[str, Any] | None = None,
) -> None:
    write_json(
        {
            "status": state,
            "message": message,
            "started_at": started_at,
            "updated_at": datetime.now().isoformat(timespec="seconds"),
            "current_stage": "NNE2 Step 1-3 and directed hierarchy",
            "missing_data_policy": "skip_record_without_processing",
            "run_directory": str(layout.run_root.resolve()),
            "timings_seconds": timings,
            "summary": summary or {},
        },
        layout.status_file,
    )


def _environment() -> dict[str, Any]:
    return {
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "scikit_image": skimage.__version__,
        "networkx": nx.__version__,
        "matplotlib": matplotlib.__version__,
        "nibabel": nib.__version__,
    }


def _stack_cache_key(metadata: StackMetadata, config: NNE2Config) -> str:
    source = {
        "cache_schema": "nne2_step1_step2_v2",
        "stack": metadata.stack_name,
        "frame_count": len(metadata.frame_files),
        "total_bytes": sum(item.stat().st_size for item in metadata.frame_files),
        "latest_mtime_ns": max(item.stat().st_mtime_ns for item in metadata.frame_files),
        "processed_spacing_xyz_um": metadata.processed_spacing_xyz_um,
        "segmentation": {
            "gaussian_sigma_um": config.gaussian_sigma_um,
            "background_sigma_um": config.background_sigma_um,
            "foreground_quantile": config.foreground_quantile,
            "min_component_voxels": config.min_component_voxels,
            "closing_iterations": config.closing_iterations,
        },
        "graph": asdict(graph_config_from_nne2(config)),
    }
    return hashlib.sha256(json.dumps(source, sort_keys=True).encode("utf-8")).hexdigest()[:16]


def _load_cache(
    path: Path,
) -> tuple[SegmentationResult, np.ndarray] | None:
    if not path.is_file():
        return None
    try:
        with np.load(path, allow_pickle=False) as data:
            normalized = np.asarray(data["normalized_zyx"], dtype=np.float32)
            candidate = np.asarray(data["candidate_mask_zyx"], dtype=bool)
            mask = np.asarray(data["mask_zyx"], dtype=bool)
            removed_mask = np.asarray(data["removed_mask_zyx"], dtype=bool)
            skeleton = np.asarray(data["skeleton_zyx"], dtype=bool)
            threshold = float(data["threshold"])
            before = int(data["component_count_before"])
            after = int(data["component_count_after"])
            removed = int(data["removed_voxel_count"])
            decisions = json.loads(str(data["component_decisions_json"].item()))
        if not (
            normalized.shape == candidate.shape == mask.shape == removed_mask.shape == skeleton.shape
        ):
            return None
        segmentation = SegmentationResult(
            normalized_zyx=normalized,
            vessel_score_zyx=normalized,
            candidate_mask_zyx=candidate,
            mask_zyx=mask,
            removed_mask_zyx=removed_mask,
            threshold=threshold,
            component_count_before=before,
            component_count_after=after,
            removed_voxel_count=removed,
            component_decisions=decisions,
        )
        return segmentation, skeleton
    except (OSError, ValueError, KeyError):
        return None


def _write_cache(
    path: Path, segmentation: SegmentationResult, skeleton_zyx: np.ndarray
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary,
        normalized_zyx=np.asarray(segmentation.normalized_zyx, dtype=np.float16),
        candidate_mask_zyx=np.asarray(segmentation.candidate_mask_zyx, dtype=np.uint8),
        mask_zyx=np.asarray(segmentation.mask_zyx, dtype=np.uint8),
        removed_mask_zyx=np.asarray(segmentation.removed_mask_zyx, dtype=np.uint8),
        skeleton_zyx=np.asarray(skeleton_zyx, dtype=np.uint8),
        threshold=np.asarray(segmentation.threshold),
        component_count_before=np.asarray(segmentation.component_count_before),
        component_count_after=np.asarray(segmentation.component_count_after),
        removed_voxel_count=np.asarray(segmentation.removed_voxel_count),
        component_decisions_json=np.asarray(json.dumps(segmentation.component_decisions)),
    )
    temporary.replace(path)


def _load_graph_cache(
    path: Path,
    expected_shape_zyx: tuple[int, int, int],
) -> NNE2CenterlineResult | None:
    if not path.is_file():
        return None
    try:
        with gzip.open(path, "rb") as stream:
            result = pickle.load(stream)  # noqa: S301 - trusted, locally generated cache only
        if not isinstance(result, NNE2CenterlineResult):
            return None
        if result.skeleton_zyx.shape != expected_shape_zyx:
            return None
        return result
    except (OSError, EOFError, pickle.UnpicklingError, AttributeError, ValueError):
        return None


def _write_graph_cache(path: Path, result: NNE2CenterlineResult) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with gzip.open(temporary, "wb", compresslevel=1) as stream:
        pickle.dump(result, stream, protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(path)


def _consistent_stack_metadata(records: list[NNE2Record]) -> tuple[float, float]:
    xy_values = {round(float(item.stack_um_per_px), 9) for item in records if item.stack_um_per_px}
    z_values = {round(float(item.stack_step_um), 9) for item in records if item.stack_step_um}
    if len(xy_values) != 1 or len(z_values) != 1:
        raise ValueError("Complete records disagree on stack pixel size or Z step")
    return next(iter(xy_values)), next(iter(z_values))


def _group_tree_records(
    records: list[NNE2Record],
) -> dict[tuple[str, str], list[NNE2Record]]:
    output: dict[tuple[str, str], list[NNE2Record]] = defaultdict(list)
    for item in records:
        assert item.tree_key is not None and item.stack_name is not None
        output[(item.tree_key, item.stack_name)].append(item)
    return dict(output)


def _run_nne2_pipeline_legacy(
    config: NNE2Config, *, verbose: bool = False
) -> NNE2PipelineRun:
    config.validate()
    layout = create_nne2_output_layout(config.output_root, config.input_dir)
    logger = configure_logging(layout.log_file, verbose=verbose)
    started_at = datetime.now().isoformat(timespec="seconds")
    timings: dict[str, float] = {}
    write_json({"stage": "nne2_step1_step3_directed", "nne2": config.report()}, layout.config_file)
    _status(layout, "running", "Pipeline started.", started_at, timings)
    logger.info("NNE2 input: %s", config.input_dir)
    logger.info("Output run: %s", layout.run_root)
    logger.info("Missing critical data policy: skip without processing")
    logger.info("TIFF I/O workers: %d", config.io_workers)

    stack_summaries: list[dict[str, Any]] = []
    tree_summaries: list[dict[str, Any]] = []
    runtime_skips: list[dict[str, Any]] = []
    report_images: list[Path] = []
    required_outputs: list[Path] = []
    try:
        with _timed_stage("catalog_and_complete_record_filter", logger, timings):
            catalog = load_nne2_catalog(config.input_dir)
            write_csv([item.row() for item in catalog.records], layout.inventory / "all_records.csv")
            write_csv(
                [item.row() for item in catalog.complete_records],
                layout.inventory / "complete_records.csv",
            )
            write_csv(
                [item.row() for item in catalog.skipped_records],
                layout.inventory / "skipped_missing_records.csv",
            )
            write_json(catalog.report(), layout.reports / "catalog_report.json")
            selected = list(
                catalog.select(
                    subject_id=config.subject_id,
                    tree_id=config.tree_id,
                    stack_name=config.stack_name,
                )
            )
            logger.info(
                "Catalog: %d total, %d complete, %d skipped missing, %d selected",
                len(catalog.records),
                len(catalog.complete_records),
                len(catalog.skipped_records),
                len(selected),
            )
        if not selected:
            raise ValueError("No complete NNE2 records match the requested filters")

        if config.stage == "inventory":
            summary = {
                **catalog.report(),
                "selected_complete_record_count": len(selected),
                "status": "PASS",
                "processing_note": "Inventory only; image stacks were not segmented.",
            }
            write_json(summary, layout.reports / "acceptance_summary.json")
            write_acceptance_html(layout.html_report, summary, [], layout.run_root)
            _status(layout, "completed", "Inventory completed.", started_at, timings, summary)
            write_json(
                {"run_directory": str(layout.run_root), "status": "completed"},
                layout.latest_file,
            )
            return NNE2PipelineRun(layout.run_root, layout.html_report, "completed", summary)

        by_stack: dict[str, list[NNE2Record]] = defaultdict(list)
        for item in selected:
            assert item.stack_name is not None
            by_stack[item.stack_name].append(item)
        stack_names = sorted(by_stack)
        if config.max_stacks is not None:
            stack_names = stack_names[: config.max_stacks]
        selected = [item for name in stack_names for item in by_stack[name]]
        tree_groups = _group_tree_records(selected)
        eligible_tree_groups: dict[tuple[str, str], list[NNE2Record]] = {}
        for key, records in tree_groups.items():
            if any(item.branching_order == 0 for item in records):
                eligible_tree_groups[key] = records
            else:
                runtime_skips.append(
                    {
                        "item_type": "tree_stack_group",
                        "item_id": f"{key[0]}__{key[1]}",
                        "reason": "no_complete_diving_trunk_branching_order_zero_record",
                    }
                )
        eligible_stacks = {key[1] for key in eligible_tree_groups}
        stack_names = [name for name in stack_names if name in eligible_stacks]
        if not stack_names:
            raise ValueError("No selected tree/stack group has a complete diving-trunk record")

        for stack_number, stack_name in enumerate(stack_names, start=1):
            records = [item for item in selected if item.stack_name == stack_name]
            logger.info("Stack %d/%d: %s", stack_number, len(stack_names), stack_name)
            try:
                xy_spacing, z_spacing = _consistent_stack_metadata(records)
                assert records[0].stack_dir is not None
                metadata = inspect_stack(
                    records[0].stack_dir,
                    xy_spacing_um=xy_spacing,
                    z_spacing_um=z_spacing,
                    target_xy_spacing_um=config.target_xy_spacing_um,
                )
                valid_records = []
                for item in records:
                    assert item.stack_index is not None
                    if item.stack_index > len(metadata.frame_files):
                        runtime_skips.append(
                            {
                                "item_type": "record",
                                "item_id": item.record_id,
                                "reason": "stack_index_outside_available_frames",
                            }
                        )
                    else:
                        valid_records.append(item)
                records = valid_records
                valid_keys = {
                    (item.tree_key, item.stack_name)
                    for item in records
                    if item.branching_order == 0
                }
                if not valid_keys:
                    raise ValueError("No valid Branching Order 0 record remains after frame checks")
                output_dir = stack_output_dir(layout, stack_name)
                cache_file = (
                    layout.cache / stack_name / f"{_stack_cache_key(metadata, config)}.npz"
                )
                graph_cache_file = cache_file.with_suffix(".graph.pkl.gz")
                cached = _load_cache(cache_file) if config.use_cache else None
                if cached is None:
                    with _timed_stage("load_and_downsample_unique_stack", logger, timings):
                        raw = load_stack(metadata, workers=config.io_workers)
                    with _timed_stage("step1_segment_and_clean", logger, timings):
                        segmentation = segment_vessels(
                            raw, metadata.processed_spacing_xyz_um, config, logger
                        )
                    del raw
                    with _timed_stage("step2_centerline_and_step3_undirected_graph", logger, timings):
                        centerline = build_centerline_graph(
                            segmentation.mask_zyx,
                            metadata.processed_spacing_xyz_um,
                            graph_connectivity=config.graph_connectivity,
                            logger=logger,
                        )
                    if config.use_cache:
                        with _timed_stage("write_reusable_stack_cache", logger, timings):
                            _write_cache(cache_file, segmentation, centerline.skeleton_zyx)
                            _write_graph_cache(graph_cache_file, centerline)
                    cache_hit = False
                else:
                    segmentation, cached_skeleton = cached
                    centerline = _load_graph_cache(
                        graph_cache_file, tuple(int(value) for value in segmentation.mask_zyx.shape)
                    )
                    if centerline is None:
                        with _timed_stage("cached_step2_step3_graph_rebuild", logger, timings):
                            centerline = build_centerline_graph(
                                segmentation.mask_zyx,
                                metadata.processed_spacing_xyz_um,
                                cached_skeleton_zyx=cached_skeleton,
                                graph_connectivity=config.graph_connectivity,
                                logger=logger,
                            )
                        if config.use_cache:
                            _write_graph_cache(graph_cache_file, centerline)
                    else:
                        logger.info("Reused branch-graph cache: %s", graph_cache_file)
                    cache_hit = True
                    logger.info("Reused segmentation and skeleton cache: %s", cache_file)

                write_json(metadata.report(), output_dir / "reports" / "stack_metadata.json")
                write_json(segmentation.report(), output_dir / "reports" / "segmentation_report.json")
                write_json(centerline.report(), output_dir / "reports" / "centerline_report.json")
                stack_validation = validate_stack_result(
                    segmentation, centerline.skeleton_zyx, centerline.graph
                )
                write_json(stack_validation, output_dir / "reports" / "stack_acceptance.json")
                exported = export_undirected_graph(centerline.graph, output_dir)
                required_outputs.extend(exported)
                if config.write_nifti:
                    required_outputs.append(
                        write_stack_nifti(
                            segmentation.mask_zyx,
                            output_dir / "volumes" / "cleaned_vessel_mask.nii.gz",
                            metadata.processed_spacing_xyz_um,
                            "NNE2 cleaned segmented vessel mask; stack XYZ coordinates",
                        )
                    )
                    required_outputs.append(
                        write_stack_nifti(
                            centerline.skeleton_zyx,
                            output_dir / "volumes" / "coarse_centerline.nii.gz",
                            metadata.processed_spacing_xyz_um,
                            "NNE2 coarse centerline; stack XYZ coordinates",
                        )
                    )
                if config.visualizations_enabled:
                    segmentation_image = render_segmentation_overview(
                        segmentation,
                        centerline.skeleton_zyx,
                        output_dir / "visualizations" / "segmentation_and_centerline.png",
                    )
                    graph_image = render_undirected_graph(
                        centerline.graph,
                        output_dir / "visualizations" / "undirected_branch_graph.png",
                    )
                    report_images.extend((segmentation_image, graph_image))
                stack_summary = {
                    "stack_name": stack_name,
                    "status": stack_validation["status"],
                    "cache_hit": cache_hit,
                    "input_shape_zyx": metadata.original_shape_zyx,
                    "processed_shape_zyx": metadata.processed_shape_zyx,
                    "foreground_fraction": float(np.mean(segmentation.mask_zyx)),
                    "skeleton_voxel_count": int(np.count_nonzero(centerline.skeleton_zyx)),
                    "branch_count": len(centerline.graph.branches),
                }
                stack_summaries.append(stack_summary)

                for (tree_key, group_stack), group_records in sorted(eligible_tree_groups.items()):
                    if group_stack != stack_name:
                        continue
                    group_records = [item for item in group_records if item in records]
                    if not group_records or not any(item.branching_order == 0 for item in group_records):
                        continue
                    try:
                        with _timed_stage("reference_registration_and_anchor_matching", logger, timings):
                            anchors = match_measurement_anchors(
                                group_records,
                                segmentation.normalized_zyx,
                                centerline.graph,
                                metadata.processed_spacing_xyz_um,
                                config,
                                logger,
                            )
                        with _timed_stage("directed_hierarchy", logger, timings):
                            hierarchy = build_directed_hierarchy(
                                tree_key, stack_name, centerline.graph, anchors, logger
                            )
                        tree_dir = tree_output_dir(layout, tree_key, stack_name)
                        tree_files = export_hierarchy(hierarchy, tree_dir)
                        write_csv(
                            [item.row() for item in group_records],
                            tree_dir / "tables" / "source_measurement_records.csv",
                        )
                        hierarchy_validation = validate_hierarchy(hierarchy)
                        write_json(
                            hierarchy_validation,
                            tree_dir / "reports" / "hierarchy_acceptance.json",
                        )
                        required_outputs.extend(tree_files)
                        if config.visualizations_enabled:
                            map_image = render_tree_map(
                                group_records[0],
                                tree_dir / "visualizations" / "subject_tree_map.png",
                            )
                            registration_image = render_anchor_registration(
                                hierarchy,
                                group_records,
                                segmentation.normalized_zyx,
                                metadata.processed_spacing_xyz_um,
                                tree_dir / "visualizations" / "reference_stack_registration.png",
                            )
                            hierarchy_image = render_hierarchy(
                                hierarchy,
                                centerline.graph,
                                tree_dir / "visualizations" / "directed_hierarchy.png",
                            )
                            component_image = render_component_decisions(
                                hierarchy,
                                centerline.graph,
                                tree_dir / "visualizations" / "hierarchy_component_decisions.png",
                            )
                            report_images.extend(
                                (map_image, registration_image, component_image, hierarchy_image)
                            )
                        tree_summaries.append(
                            {
                                **hierarchy.report(),
                                "status": hierarchy_validation["status"],
                            }
                        )
                    except Exception as exc:
                        logger.exception("Tree hierarchy failed for %s/%s: %s", tree_key, stack_name, exc)
                        runtime_skips.append(
                            {
                                "item_type": "tree_stack_group",
                                "item_id": f"{tree_key}__{stack_name}",
                                "reason": f"processing_error:{type(exc).__name__}:{exc}",
                            }
                        )
            except Exception as exc:
                logger.exception("Stack failed for %s: %s", stack_name, exc)
                stack_summaries.append(
                    {"stack_name": stack_name, "status": "FAIL", "error": f"{type(exc).__name__}: {exc}"}
                )
            finally:
                gc.collect()

        write_csv(stack_summaries, layout.reports / "stack_summary.csv")
        write_csv(tree_summaries, layout.reports / "tree_hierarchy_summary.csv")
        write_csv(runtime_skips, layout.inventory / "runtime_skips.csv")
        write_json(_environment(), layout.reports / "environment.json")
        file_check = required_files_status(required_outputs)
        successful_stacks = sum(item.get("status") in {"PASS", "WARN"} for item in stack_summaries)
        successful_trees = sum(item.get("status") in {"PASS", "WARN"} for item in tree_summaries)
        if successful_stacks == 0 or successful_trees == 0 or file_check["status"] != "PASS":
            overall_status = "FAIL"
        elif any(item.get("status") == "WARN" for item in tree_summaries):
            overall_status = "WARN"
        else:
            overall_status = "PASS"
        summary = {
            "status": overall_status,
            "catalog_total_records": len(catalog.records),
            "catalog_complete_records": len(catalog.complete_records),
            "catalog_skipped_missing_records": len(catalog.skipped_records),
            "selected_stack_count": len(stack_names),
            "successful_stack_count": successful_stacks,
            "failed_stack_count": len(stack_summaries) - successful_stacks,
            "successful_hierarchy_count": successful_trees,
            "failed_or_skipped_hierarchy_count": len(runtime_skips),
            "required_output_status": file_check["status"],
            "cache_enabled": config.use_cache,
            "acceleration": (
                "unique-stack processing, parallel TIFF reads, XY downsampling, no graph smoothing, "
                "persistent segmentation/skeleton/branch-graph cache"
            ),
            "missing_data_policy": "skipped and never sent to image processing",
        }
        write_json(summary, layout.reports / "acceptance_summary.json")
        write_json(file_check, layout.reports / "required_files_check.json")
        write_acceptance_html(
            layout.html_report,
            summary,
            report_images[:30],
            layout.run_root,
        )
        status = "failed" if overall_status == "FAIL" else "completed"
        _status(
            layout,
            status,
            "NNE2 pipeline completed." if status == "completed" else "NNE2 pipeline completed with failed acceptance checks.",
            started_at,
            timings,
            summary,
        )
        write_json(
            {
                "run_directory": str(layout.run_root.resolve()),
                "status": status,
                "acceptance_status": overall_status,
                "updated_at": datetime.now().isoformat(timespec="seconds"),
            },
            layout.latest_file,
        )
        return NNE2PipelineRun(layout.run_root, layout.html_report, status, summary)
    except Exception as exc:
        logger.exception("NNE2 pipeline failed: %s", exc)
        _status(
            layout,
            "failed",
            f"{type(exc).__name__}: {exc}",
            started_at,
            timings,
        )
        raise


def _graph_cache_for_source(
    layout: NNE2OutputLayout,
    stack_name: str,
    source_cache_key: str,
    config: NNE2Config,
) -> Path:
    graph_settings = asdict(graph_config_from_nne2(config))
    digest = hashlib.sha256(
        json.dumps(graph_settings, sort_keys=True).encode("utf-8")
    ).hexdigest()[:12]
    return layout.cache / stack_name / f"{source_cache_key}_{digest}.graph.pkl.gz"


def _process_tree_hierarchies(
    *,
    layout: NNE2OutputLayout,
    stack_name: str,
    records: list[NNE2Record],
    normalized_zyx: np.ndarray,
    spacing_xyz_um: tuple[float, float, float],
    centerline: NNE2CenterlineResult,
    config: NNE2Config,
    logger: logging.Logger,
    timings: dict[str, float],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[Path], list[Path]]:
    summaries: list[dict[str, Any]] = []
    skips: list[dict[str, Any]] = []
    outputs: list[Path] = []
    images: list[Path] = []
    groups = _group_tree_records(records)
    for (tree_key, group_stack), group_records in sorted(groups.items()):
        if group_stack != stack_name:
            continue
        if not any(item.branching_order == 0 for item in group_records):
            skips.append(
                {
                    "item_type": "tree_stack_group",
                    "item_id": f"{tree_key}__{stack_name}",
                    "reason": "no_complete_diving_trunk_branching_order_zero_record",
                }
            )
            continue
        valid_records = [
            item for item in group_records
            if item.stack_index is not None and item.stack_index <= normalized_zyx.shape[0]
        ]
        if not any(item.branching_order == 0 for item in valid_records):
            skips.append(
                {
                    "item_type": "tree_stack_group",
                    "item_id": f"{tree_key}__{stack_name}",
                    "reason": "diving_trunk_stack_index_outside_available_frames",
                }
            )
            continue
        try:
            with _timed_stage("reference_registration_and_anchor_matching", logger, timings):
                anchors = match_measurement_anchors(
                    valid_records,
                    normalized_zyx,
                    centerline.graph,
                    spacing_xyz_um,
                    config,
                    logger,
                )
            with _timed_stage("directed_hierarchy", logger, timings):
                hierarchy = build_directed_hierarchy(
                    tree_key, stack_name, centerline.graph, anchors, logger
                )
            tree_dir = tree_output_dir(layout, tree_key, stack_name)
            tree_files = export_hierarchy(hierarchy, tree_dir)
            source_records = tree_dir / "tables" / "source_measurement_records.csv"
            write_csv([item.row() for item in valid_records], source_records)
            hierarchy_validation = validate_hierarchy(hierarchy)
            validation_path = tree_dir / "reports" / "hierarchy_acceptance.json"
            write_json(hierarchy_validation, validation_path)
            outputs.extend([*tree_files, source_records, validation_path])
            if config.write_nifti:
                retained_xyz = np.zeros_like(centerline.graph.source_skeleton, dtype=bool)
                retained_branch_ids = {item.branch_id for item in hierarchy.branches}
                retained_node_ids = set(hierarchy.directed_graph.nodes)
                for node_id in retained_node_ids:
                    node = centerline.graph.nodes[int(node_id)]
                    indices = np.asarray(node.voxel_indices_xyz, dtype=int)
                    retained_xyz[tuple(indices.T)] = True
                for branch in centerline.graph.branches:
                    if branch.branch_id in retained_branch_ids:
                        indices = np.asarray(branch.voxel_indices_xyz, dtype=int)
                        retained_xyz[tuple(indices.T)] = True
                excluded_xyz = centerline.graph.source_skeleton & ~retained_xyz
                retained_path = write_stack_nifti(
                    np.transpose(retained_xyz, (2, 1, 0)),
                    tree_dir / "volumes" / "retained_root_component_centerline.nii.gz",
                    spacing_xyz_um,
                    "Per-tree centerline retained in the BO0 root component",
                )
                excluded_path = write_stack_nifti(
                    np.transpose(excluded_xyz, (2, 1, 0)),
                    tree_dir / "volumes" / "excluded_island_centerlines.nii.gz",
                    spacing_xyz_um,
                    "Centerline islands outside the BO0 root component",
                )
                component_report_path = tree_dir / "reports" / "root_component_report.json"
                write_json(
                    {
                        "root_node_id": hierarchy.root_node_id,
                        "root_branch_id": hierarchy.root_branch_id,
                        "retained_centerline_voxel_count": int(np.count_nonzero(retained_xyz)),
                        "excluded_centerline_voxel_count": int(np.count_nonzero(excluded_xyz)),
                        "reconstruction_matches_source_skeleton": bool(
                            np.array_equal(retained_xyz | excluded_xyz, centerline.graph.source_skeleton)
                        ),
                        "retained_and_excluded_do_not_overlap": not bool(
                            np.any(retained_xyz & excluded_xyz)
                        ),
                    },
                    component_report_path,
                )
                outputs.extend((retained_path, excluded_path, component_report_path))
            if config.visualizations_enabled:
                tree_images = [
                    render_tree_map(
                        valid_records[0], tree_dir / "visualizations" / "subject_tree_map.png"
                    ),
                    render_anchor_registration(
                        hierarchy,
                        valid_records,
                        normalized_zyx,
                        spacing_xyz_um,
                        tree_dir / "visualizations" / "reference_stack_registration.png",
                    ),
                    render_component_decisions(
                        hierarchy,
                        centerline.graph,
                        tree_dir / "visualizations" / "hierarchy_component_decisions.png",
                    ),
                    render_hierarchy(
                        hierarchy,
                        centerline.graph,
                        tree_dir / "visualizations" / "directed_hierarchy.png",
                    ),
                ]
                images.extend(tree_images)
                outputs.extend(tree_images)
            summaries.append({**hierarchy.report(), "status": hierarchy_validation["status"]})
        except Exception as exc:
            logger.exception("Tree hierarchy failed for %s/%s: %s", tree_key, stack_name, exc)
            skips.append(
                {
                    "item_type": "tree_stack_group",
                    "item_id": f"{tree_key}__{stack_name}",
                    "reason": f"processing_error:{type(exc).__name__}:{exc}",
                }
            )
    return summaries, skips, outputs, images


def run_nne2_pipeline(
    config: NNE2Config, *, verbose: bool = False
) -> NNE2PipelineRun:
    """Run inventory, independent Step 1-2, independent Step 3, or the full workflow."""
    config.validate()
    layout = create_nne2_output_layout(
        config.output_root, config.input_dir, stage=config.stage
    )
    logger = configure_logging(layout.log_file, verbose=verbose)
    started_at = datetime.now().isoformat(timespec="seconds")
    timings: dict[str, float] = {}
    write_json(
        {
            "schema_version": "2.0",
            "stage": config.stage,
            "nne2": config.report(),
            "graph": asdict(graph_config_from_nne2(config)),
        },
        layout.config_file,
    )
    _status(layout, "running", f"NNE2 stage {config.stage} started.", started_at, timings)
    logger.info("NNE2 input: %s", config.input_dir)
    logger.info("Output run: %s", layout.run_root)
    logger.info("Stage: %s; missing critical data policy: skip", config.stage)

    stack_summaries: list[dict[str, Any]] = []
    tree_summaries: list[dict[str, Any]] = []
    runtime_skips: list[dict[str, Any]] = []
    report_images: list[Path] = []
    required_outputs: list[Path] = [layout.config_file, layout.log_file]
    try:
        with _timed_stage("catalog_and_complete_record_filter", logger, timings):
            catalog = load_nne2_catalog(config.input_dir)
            all_records_path = layout.inventory / "all_records.csv"
            complete_records_path = layout.inventory / "complete_records.csv"
            skipped_records_path = layout.inventory / "skipped_missing_records.csv"
            catalog_report_path = layout.reports / "catalog_report.json"
            write_csv([item.row() for item in catalog.records], all_records_path)
            write_csv([item.row() for item in catalog.complete_records], complete_records_path)
            write_csv([item.row() for item in catalog.skipped_records], skipped_records_path)
            write_json(catalog.report(), catalog_report_path)
            required_outputs.extend((all_records_path, complete_records_path, catalog_report_path))
            if catalog.skipped_records:
                required_outputs.append(skipped_records_path)
            selected = list(
                catalog.select(
                    subject_id=config.subject_id,
                    tree_id=config.tree_id,
                    stack_name=config.stack_name,
                )
            )
        if not selected:
            raise ValueError("No complete NNE2 records match the requested filters")

        if config.stage == "inventory":
            summary = {
                **catalog.report(),
                "selected_complete_record_count": len(selected),
                "status": "PASS",
                "processing_note": "Inventory only; image stacks were not segmented.",
            }
            summary_path = layout.reports / "acceptance_summary.json"
            write_json(summary, summary_path)
            write_acceptance_html(layout.html_report, summary, [], layout.run_root)
            _status(layout, "completed", "Inventory completed.", started_at, timings, summary)
            write_json({"run_directory": str(layout.run_root), "status": "completed"}, layout.latest_file)
            return NNE2PipelineRun(layout.run_root, layout.html_report, "completed", summary)

        by_stack: dict[str, list[NNE2Record]] = defaultdict(list)
        for item in selected:
            assert item.stack_name is not None
            by_stack[item.stack_name].append(item)
        stack_names = sorted(by_stack)
        source_manifests: dict[str, dict[str, Any]] = {}
        if config.stage == "hierarchical-graph":
            assert config.source_run is not None
            source_manifests = load_source_manifests(config.source_run)
            stack_names = [name for name in stack_names if name in source_manifests]
            source_path = layout.reports / "source_run.json"
            write_json(
                {
                    "source_run": str(config.source_run),
                    "available_stack_manifests": len(source_manifests),
                    "selected_stack_names": stack_names,
                },
                source_path,
            )
            required_outputs.append(source_path)
        if config.max_stacks is not None:
            stack_names = stack_names[: config.max_stacks]
        if not stack_names:
            raise ValueError("No selected stack is available for this stage")

        graph_config = graph_config_from_nne2(config)
        for stack_number, stack_name in enumerate(stack_names, start=1):
            logger.info("Stack %d/%d: %s", stack_number, len(stack_names), stack_name)
            output_dir = stack_output_dir(layout, stack_name)
            records = by_stack[stack_name]
            try:
                if config.stage == "hierarchical-graph":
                    assert config.source_run is not None
                    with _timed_stage("load_and_verify_step1_step2_artifacts", logger, timings):
                        loaded = load_step2_stack(config.source_run, source_manifests[stack_name])
                    provenance_path = output_dir / "reports" / "source_provenance.json"
                    write_json(
                        {
                            "source_run": str(config.source_run),
                            "manifest_path": source_manifests[stack_name]["manifest_path"],
                            "manifest_schema_version": source_manifests[stack_name]["schema_version"],
                            "verified_artifacts": source_manifests[stack_name]["artifacts"],
                        },
                        provenance_path,
                    )
                    required_outputs.append(provenance_path)
                    normalized = loaded.normalized_zyx
                    mask = loaded.mask_zyx
                    skeleton_zyx = loaded.skeleton_zyx
                    radius_zyx_um = loaded.radius_zyx_um
                    spacing = loaded.spacing_xyz_um
                    skeleton_zyx, computed_radius, skeleton_report = extract_centerline(
                        mask,
                        spacing,
                        cached_skeleton_zyx=skeleton_zyx,
                        logger=logger,
                    )
                    if radius_zyx_um.shape != computed_radius.shape:
                        raise ValueError("Stored Step 2 radius shape is inconsistent")
                    source_cache_key = str(source_manifests[stack_name]["cache_key"])
                    graph_cache = _graph_cache_for_source(
                        layout, stack_name, source_cache_key, config
                    )
                    centerline = _load_graph_cache(graph_cache, mask.shape) if config.use_cache else None
                    if centerline is None:
                        with _timed_stage("step3_branch_graph", logger, timings):
                            graph = build_graph_from_centerline(
                                skeleton_zyx, mask, spacing, graph_config, logger=logger
                            )
                        centerline = NNE2CenterlineResult(
                            skeleton_zyx=skeleton_zyx,
                            radius_zyx_um=radius_zyx_um,
                            skeleton_report=skeleton_report,
                            graph=graph,
                        )
                        if config.use_cache:
                            _write_graph_cache(graph_cache, centerline)
                        cache_hit = False
                    else:
                        logger.info("Reused Step 3 branch-graph cache: %s", graph_cache)
                        cache_hit = True
                else:
                    xy_spacing, z_spacing = _consistent_stack_metadata(records)
                    assert records[0].stack_dir is not None
                    metadata = inspect_stack(
                        records[0].stack_dir,
                        xy_spacing_um=xy_spacing,
                        z_spacing_um=z_spacing,
                        target_xy_spacing_um=config.target_xy_spacing_um,
                    )
                    if int(np.prod(metadata.processed_shape_zyx)) > config.max_voxel_count:
                        raise MemoryError(
                            f"Processed stack would contain {np.prod(metadata.processed_shape_zyx):,} "
                            f"voxels, above max_voxel_count={config.max_voxel_count:,}"
                        )
                    cache_key = _stack_cache_key(metadata, config)
                    cache_file = layout.cache / stack_name / f"{cache_key}.npz"
                    cached = _load_cache(cache_file) if config.use_cache else None
                    if cached is None:
                        with _timed_stage("load_and_downsample_unique_stack", logger, timings):
                            raw = load_stack(metadata, workers=config.io_workers)
                        with _timed_stage("step1_segment_and_component_cleanup", logger, timings):
                            segmentation = segment_vessels(raw, metadata.processed_spacing_xyz_um, config, logger)
                        del raw
                        with _timed_stage("step2_radius_and_coarse_centerline", logger, timings):
                            skeleton_zyx, radius_zyx_um, skeleton_report = extract_centerline(
                                segmentation.mask_zyx,
                                metadata.processed_spacing_xyz_um,
                                logger=logger,
                            )
                        if config.use_cache:
                            _write_cache(cache_file, segmentation, skeleton_zyx)
                        cache_hit = False
                    else:
                        segmentation, cached_skeleton = cached
                        with _timed_stage("reuse_step1_and_validate_step2", logger, timings):
                            skeleton_zyx, radius_zyx_um, skeleton_report = extract_centerline(
                                segmentation.mask_zyx,
                                metadata.processed_spacing_xyz_um,
                                cached_skeleton_zyx=cached_skeleton,
                                logger=logger,
                            )
                        cache_hit = True
                        logger.info("Reused Step 1-2 cache: %s", cache_file)
                    normalized = segmentation.normalized_zyx
                    mask = segmentation.mask_zyx
                    spacing = metadata.processed_spacing_xyz_um
                    preprocess_validation = validate_preprocess_result(
                        segmentation,
                        skeleton_zyx,
                        radius_zyx_um,
                        island_warning_fraction=config.island_warning_fraction,
                        island_fail_fraction=config.island_fail_fraction,
                    )
                    acceptance_path = output_dir / "reports" / "step1_step2_acceptance.json"
                    write_json(preprocess_validation, acceptance_path)
                    artifact_paths, manifest_path = write_step1_step2_artifacts(
                        run_root=layout.run_root,
                        output_dir=output_dir,
                        stack_name=stack_name,
                        segmentation=segmentation,
                        skeleton_zyx=skeleton_zyx,
                        radius_zyx_um=radius_zyx_um,
                        spacing_xyz_um=spacing,
                        stack_metadata_report=metadata.report(),
                        skeleton_report=skeleton_report.report(),
                        cache_key=cache_key,
                        write_nifti=config.write_nifti,
                    )
                    required_outputs.extend([*artifact_paths, acceptance_path])
                    if config.visualizations_enabled:
                        stack_images = [
                            render_segmentation_overview(
                                segmentation,
                                skeleton_zyx,
                                output_dir / "visualizations" / "segmentation_and_centerline.png",
                            ),
                            render_step1_component_cleanup(
                                segmentation,
                                output_dir / "visualizations" / "step1_component_cleanup.png",
                            ),
                            render_step2_radius_and_centerline(
                                segmentation,
                                skeleton_zyx,
                                radius_zyx_um,
                                output_dir / "visualizations" / "step2_radius_and_centerline.png",
                            ),
                        ]
                        report_images.extend(stack_images)
                        required_outputs.extend(stack_images)
                    if config.stage == "preprocess":
                        stack_summaries.append(
                            {
                                "stack_name": stack_name,
                                "status": preprocess_validation["status"],
                                "cache_hit": cache_hit,
                                "processed_shape_zyx": metadata.processed_shape_zyx,
                                "foreground_fraction": float(np.mean(mask)),
                                "skeleton_voxel_count": int(np.count_nonzero(skeleton_zyx)),
                                "manifest": str(manifest_path.relative_to(layout.run_root)),
                            }
                        )
                        continue

                    graph_cache = cache_file.with_suffix(".graph.pkl.gz")
                    centerline = _load_graph_cache(graph_cache, mask.shape) if config.use_cache else None
                    if centerline is None:
                        with _timed_stage("step3_branch_graph", logger, timings):
                            graph = build_graph_from_centerline(
                                skeleton_zyx, mask, spacing, graph_config, logger=logger
                            )
                        centerline = NNE2CenterlineResult(
                            skeleton_zyx=skeleton_zyx,
                            radius_zyx_um=radius_zyx_um,
                            skeleton_report=skeleton_report,
                            graph=graph,
                        )
                        if config.use_cache:
                            _write_graph_cache(graph_cache, centerline)
                    else:
                        logger.info("Reused Step 3 branch-graph cache: %s", graph_cache)

                with _timed_stage("step3_graph_export_and_acceptance", logger, timings):
                    graph_files, graph_acceptance = export_and_validate_step3_graph(
                        centerline.graph, output_dir, graph_config
                    )
                required_outputs.extend(graph_files)
                if config.visualizations_enabled:
                    graph_image = render_undirected_graph(
                        centerline.graph,
                        output_dir / "visualizations" / "undirected_branch_graph.png",
                    )
                    report_images.append(graph_image)
                    required_outputs.append(graph_image)
                    detailed_image = render_step3_diagnostics(
                        centerline.graph,
                        output_dir / "visualizations" / "step3_diagnostics.png",
                        short_branch_warning_um=config.short_branch_warning_um,
                        large_junction_warning_voxels=config.large_junction_warning_voxels,
                        high_degree_warning=config.high_degree_warning,
                    )
                    report_images.append(detailed_image)
                    required_outputs.append(detailed_image)
                hierarchy_summaries, hierarchy_skips, hierarchy_files, hierarchy_images = (
                    _process_tree_hierarchies(
                        layout=layout,
                        stack_name=stack_name,
                        records=records,
                        normalized_zyx=normalized,
                        spacing_xyz_um=spacing,
                        centerline=centerline,
                        config=config,
                        logger=logger,
                        timings=timings,
                    )
                )
                tree_summaries.extend(hierarchy_summaries)
                runtime_skips.extend(hierarchy_skips)
                required_outputs.extend(hierarchy_files)
                report_images.extend(hierarchy_images)
                stack_summaries.append(
                    {
                        "stack_name": stack_name,
                        "status": graph_acceptance["overall_status"].replace("WARNING", "WARN"),
                        "cache_hit": cache_hit,
                        "processed_shape_zyx": tuple(int(value) for value in mask.shape),
                        "foreground_fraction": float(np.mean(mask)),
                        "skeleton_voxel_count": int(np.count_nonzero(skeleton_zyx)),
                        "branch_count": len(centerline.graph.branches),
                        "cycle_count": len(centerline.graph.cycles),
                    }
                )
            except Exception as exc:
                logger.exception("Stack failed for %s: %s", stack_name, exc)
                stack_summaries.append(
                    {"stack_name": stack_name, "status": "FAIL", "error": f"{type(exc).__name__}: {exc}"}
                )
            finally:
                gc.collect()

        stack_summary_path = layout.reports / "stack_summary.csv"
        tree_summary_path = layout.reports / "tree_hierarchy_summary.csv"
        runtime_skips_path = layout.inventory / "runtime_skips.csv"
        environment_path = layout.reports / "environment.json"
        write_csv(stack_summaries, stack_summary_path)
        write_csv(tree_summaries, tree_summary_path)
        write_csv(runtime_skips, runtime_skips_path)
        write_json(_environment(), environment_path)
        required_outputs.extend((stack_summary_path, environment_path))
        if tree_summaries:
            required_outputs.append(tree_summary_path)
        if runtime_skips:
            required_outputs.append(runtime_skips_path)
        file_check = required_files_status(required_outputs)
        successful_stacks = sum(item.get("status") in {"PASS", "WARN"} for item in stack_summaries)
        successful_trees = sum(item.get("status") in {"PASS", "WARN"} for item in tree_summaries)
        needs_trees = config.stage in {"all", "hierarchical-graph"}
        if successful_stacks == 0 or file_check["status"] != "PASS" or (needs_trees and successful_trees == 0):
            overall_status = "FAIL"
        elif any(item.get("status") == "WARN" for item in [*stack_summaries, *tree_summaries]):
            overall_status = "WARN"
        else:
            overall_status = "PASS"
        summary = {
            "status": overall_status,
            "stage": config.stage,
            "source_run": str(config.source_run) if config.source_run else None,
            "catalog_total_records": len(catalog.records),
            "catalog_complete_records": len(catalog.complete_records),
            "catalog_skipped_missing_records": len(catalog.skipped_records),
            "selected_stack_count": len(stack_names),
            "successful_stack_count": successful_stacks,
            "failed_stack_count": len(stack_summaries) - successful_stacks,
            "successful_hierarchy_count": successful_trees,
            "runtime_skipped_item_count": len(runtime_skips),
            "required_output_status": file_check["status"],
            "cache_enabled": config.use_cache,
            "acceleration": (
                "unique-stack processing, parallel TIFF reads, XY downsampling, staged reuse, "
                "persistent Step 1-2 and branch-graph caches"
            ),
            "missing_data_policy": "skipped and never sent to image processing",
        }
        summary_path = layout.reports / "acceptance_summary.json"
        file_check_path = layout.reports / "required_files_check.json"
        write_json(summary, summary_path)
        write_json(file_check, file_check_path)
        write_acceptance_html(layout.html_report, summary, report_images[:30], layout.run_root)
        status = "failed" if overall_status == "FAIL" else "completed"
        _status(
            layout,
            status,
            f"NNE2 stage {config.stage} completed with {overall_status} acceptance.",
            started_at,
            timings,
            summary,
        )
        write_json(
            {
                "run_directory": str(layout.run_root.resolve()),
                "status": status,
                "acceptance_status": overall_status,
                "updated_at": datetime.now().isoformat(timespec="seconds"),
            },
            layout.latest_file,
        )
        return NNE2PipelineRun(layout.run_root, layout.html_report, status, summary)
    except Exception as exc:
        logger.exception("NNE2 pipeline failed: %s", exc)
        _status(layout, "failed", f"{type(exc).__name__}: {exc}", started_at, timings)
        raise
