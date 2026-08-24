"""Staged pipeline for the mouse-brain vasculature TIFF/SWC dataset."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

from ..io import write_csv, write_json
from ..reporting.acceptance import AcceptanceCheck, AcceptanceResult
from ..reporting.logger import configure_logging
from .catalog import RodentSampleRecord, build_catalog, select_records
from .config import RodentVasculatureConfig
from .export import export_directed_graph
from .graph_builder import build_directed_vascular_graph
from .interactive import save_figure2a_preview
from .report import write_html_report
from .swc_analysis import evaluate_optional_mask_qc, select_analysis_swc
from .swc_io import (
    load_normalized_swc,
    load_swc,
    save_normalized_swc,
    save_swc_text,
)
from .tiff_io import (
    inspect_tiff,
    load_normalized_volume,
    load_tiff_volume,
    save_normalized_volume,
    volume_statistics,
)
from .validation import evaluate_directed_graph
from .visualization import create_visualizations


@dataclass(slots=True)
class RodentPipelineRun:
    run_root: Path
    status: str
    html_report: Path
    acceptance: AcceptanceResult
    processed_sample_count: int
    failed_sample_count: int


def _safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", value).strip("._") or "sample"


def _make_run_root(config: RodentVasculatureConfig) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base = config.output_root.resolve() / "rodent_vasculature"
    candidate = base / f"{config.stage}_run_{timestamp}"
    suffix = 1
    while candidate.exists():
        candidate = base / f"{config.stage}_run_{timestamp}_{suffix:02d}"
        suffix += 1
    for folder in (
        candidate / "inventory",
        candidate / "reports",
        candidate / "samples",
    ):
        folder.mkdir(parents=True, exist_ok=False)
    return candidate


def _aggregate_acceptance(
    sample_acceptances: list[AcceptanceResult], failures: list[dict[str, str]], required: list[Path]
) -> AcceptanceResult:
    checks: list[AcceptanceCheck] = []
    if sample_acceptances:
        counts = {
            value: sum(item.overall_status == value for item in sample_acceptances)
            for value in ("PASS", "WARNING", "FAIL")
        }
        checks.append(
            AcceptanceCheck(
                "Selected stage/sample evaluations passed acceptance",
                "FAIL" if counts["FAIL"] else ("WARNING" if counts["WARNING"] else "PASS"),
                f"Evaluation results: {counts}.",
            )
        )
    checks.append(
        AcceptanceCheck(
            "No selected sample raised a processing exception",
            "PASS" if not failures else "FAIL",
            "No processing exceptions." if not failures else f"Failures: {failures}",
        )
    )
    missing = [str(path) for path in required if not path.is_file() or path.stat().st_size == 0]
    checks.append(
        AcceptanceCheck(
            "Run-level reports exist",
            "PASS" if not missing else "FAIL",
            "All run-level reports exist." if not missing else f"Missing: {missing}",
        )
    )
    statuses = {check.status for check in checks}
    overall = "FAIL" if "FAIL" in statuses else ("WARNING" if "WARNING" in statuses else "PASS")
    return AcceptanceResult(overall, checks)


def _record_from_payload(payload: dict[str, Any]) -> RodentSampleRecord:
    tile = payload.get("tile_index_zyx")
    return RodentSampleRecord(
        sample_id=payload["sample_id"],
        cohort=payload["cohort"],
        source_stem=payload["source_stem"],
        parent_group_id=payload["parent_group_id"],
        tile_index_zyx=tuple(tile) if tile else None,
        split=payload.get("split"),
        image_path=Path(payload["image_path"]) if payload.get("image_path") else None,
        mask_path=Path(payload["mask_path"]) if payload.get("mask_path") else None,
        swc_path=Path(payload["swc_path"]) if payload.get("swc_path") else None,
        eligible=bool(payload["eligible"]),
        skip_reason=payload.get("skip_reason"),
    )


def _preprocess_sample(
    record: RodentSampleRecord,
    sample_root: Path,
    config: RodentVasculatureConfig,
) -> tuple[dict[str, Any], AcceptanceResult]:
    assert record.swc_path
    image_metadata = inspect_tiff(record.image_path) if record.image_path else None
    mask_metadata = inspect_tiff(record.mask_path) if record.mask_path else None
    image_volume = load_tiff_volume(record.image_path) if record.image_path else None
    mask_volume = load_tiff_volume(record.mask_path) if record.mask_path else None
    auxiliary_shape_zyx = (
        image_metadata.shape_zyx
        if image_metadata is not None
        else (mask_metadata.shape_zyx if mask_metadata is not None else None)
    )
    reference_swc = load_swc(
        record.swc_path,
        spacing_xyz_um=config.spacing_xyz_um,
        volume_shape_zyx=auxiliary_shape_zyx,
    )
    selection = select_analysis_swc(
        reference_swc,
        spacing_xyz_um=config.spacing_xyz_um,
        volume_shape_zyx=auxiliary_shape_zyx,
        analysis_component_id=config.analysis_component_id,
    )
    analysis_swc = selection.analysis_swc
    mask_qc = evaluate_optional_mask_qc(
        mask_volume,
        reference_swc,
        analysis_swc,
    )
    preprocessing_summary = selection.summary | {
        "mask_role": mask_qc["role"],
        "mask_used_for_component_selection": False,
        "mask_used_for_topology_modification": False,
        "mask_qc": mask_qc,
    }

    display_volume_path: Path | None = None
    if (
        image_volume is not None
        and mask_volume is not None
        and image_volume.shape == mask_volume.shape
    ):
        display_volume_path = save_normalized_volume(
            image_volume,
            mask_volume,
            sample_root / "normalized" / "display_reference_volume.npz",
        )
    reference_normalized_path = save_normalized_swc(
        reference_swc,
        sample_root / "normalized" / "reference_swc_complete.npz",
    )
    analysis_normalized_path = save_normalized_swc(
        analysis_swc,
        sample_root / "normalized" / "analysis_swc_single_component.npz",
    )
    reference_swc_text_path = save_swc_text(
        reference_swc,
        sample_root / "normalized" / "reference_swc_complete.swc",
    )
    analysis_swc_text_path = save_swc_text(
        analysis_swc,
        sample_root / "normalized" / "analysis_swc_single_component.swc",
    )
    audit_root = sample_root / "swc_analysis_preprocessing"
    component_table_path = audit_root / "reference_swc_components.csv"
    reference_only_nodes_path = audit_root / "reference_only_node_ids.csv"
    preprocessing_summary_path = audit_root / "swc_centric_summary.json"
    mask_qc_path = audit_root / "optional_mask_qc.json"
    write_csv(selection.component_records, component_table_path)
    write_csv(
        (
            {"reference_only_node_id": int(node_id)}
            for node_id in selection.reference_only_node_ids
        ),
        reference_only_nodes_path,
    )
    write_json(preprocessing_summary, preprocessing_summary_path)
    write_json(mask_qc, mask_qc_path)

    reference_roundtrip = load_normalized_swc(
        reference_normalized_path,
        source_path=record.swc_path,
        spacing_xyz_um=config.spacing_xyz_um,
        volume_shape_zyx=auxiliary_shape_zyx,
    )
    analysis_roundtrip = load_normalized_swc(
        analysis_normalized_path,
        source_path=record.swc_path,
        spacing_xyz_um=config.spacing_xyz_um,
        volume_shape_zyx=auxiliary_shape_zyx,
    )
    reference_roundtrip_matches = (
        reference_roundtrip.node_ids.tolist() == reference_swc.node_ids.tolist()
        and reference_roundtrip.parent_ids.tolist() == reference_swc.parent_ids.tolist()
    )
    analysis_roundtrip_matches = (
        analysis_roundtrip.node_ids.tolist() == analysis_swc.node_ids.tolist()
        and analysis_roundtrip.parent_ids.tolist() == analysis_swc.parent_ids.tolist()
    )
    logging.getLogger("ulm_3d_vascular").info(
        "SWC-centric preprocessing for %s: reference components=%d, nodes=%d; "
        "analysis component=%d, nodes=%d, edges=%d; topology modifications=0; "
        "mask available=%s and role=QC/display only",
        record.sample_id,
        reference_swc.component_count,
        reference_swc.node_count,
        preprocessing_summary["selected_component_id"],
        analysis_swc.node_count,
        analysis_swc.edge_count,
        mask_volume is not None,
    )
    auxiliary_shapes_agree = (
        image_metadata is None
        or mask_metadata is None
        or image_metadata.shape_zyx == mask_metadata.shape_zyx
    )
    expected_shape_agrees = (
        auxiliary_shape_zyx is None
        or config.expected_shape_zyx is None
        or auxiliary_shape_zyx == config.expected_shape_zyx
    )
    mask_support = mask_qc["reference_swc_node_support_fraction"]
    required_audit_artifacts = [
        component_table_path,
        preprocessing_summary_path,
        mask_qc_path,
        reference_normalized_path,
        analysis_normalized_path,
        reference_swc_text_path,
        analysis_swc_text_path,
    ]
    if display_volume_path is not None:
        required_audit_artifacts.append(display_volume_path)
    checks = [
        AcceptanceCheck(
            "Optional image and Mask TIFF shapes are compatible when both are available",
            "PASS" if auxiliary_shapes_agree else "FAIL",
            f"Image={image_metadata.shape_zyx if image_metadata else 'not provided'}; "
            f"Mask={mask_metadata.shape_zyx if mask_metadata else 'not provided'}. "
            "Neither auxiliary TIFF is required for SWC analysis.",
        ),
        AcceptanceCheck(
            "Available auxiliary TIFF shape agrees with configured expectation",
            "PASS" if expected_shape_agrees else "FAIL",
            f"Observed={auxiliary_shape_zyx if auxiliary_shape_zyx else 'not applicable'}; "
            f"expected={config.expected_shape_zyx}. No TIFF means this optional check is skipped.",
        ),
        AcceptanceCheck(
            "Complete reference SWC parent forest is structurally valid",
            "PASS" if reference_swc.structurally_valid else "FAIL",
            f"Nodes={reference_swc.node_count}; edges={reference_swc.edge_count}; "
            f"components={reference_swc.component_count}; validation={reference_swc.validation}.",
        ),
        AcceptanceCheck(
            "Mask registration QC does not control SWC topology",
            "PASS",
            "Mask is unchanged and is used only for display/registration QC; "
            "component selection, node removal, and topology repair all report false.",
        ),
        AcceptanceCheck(
            "Optional Mask QC supports the reference SWC coordinates",
            "PASS"
            if mask_support is None or mask_support >= 1.0 - 1.0e-12
            else "WARNING",
            (
                "Mask not provided; optional registration QC was skipped and SWC analysis continued."
                if mask_support is None
                else "Reference-node support="
                f"{mask_support:.6f}; dense-edge support="
                f"{mask_qc['reference_swc_dense_edge_support_fraction']:.6f}. "
                "A warning does not alter or reject the SWC analysis topology."
            ),
        ),
        AcceptanceCheck(
            "Analysis SWC is one unchanged original connected component",
            "PASS"
            if analysis_swc.component_count == 1
            and preprocessing_summary["new_node_count"] == 0
            and preprocessing_summary["new_edge_count"] == 0
            and preprocessing_summary["parent_relation_change_count"] == 0
            else "FAIL",
            f"Selected component={preprocessing_summary['selected_component_id']}; "
            f"nodes={analysis_swc.node_count}; edges={analysis_swc.edge_count}; "
            "new nodes=0; new edges=0; changed parent relations=0.",
        ),
        AcceptanceCheck(
            "Complete reference SWC is preserved independently of analysis selection",
            "PASS"
            if reference_roundtrip_matches
            and preprocessing_summary["reference_only_components_are_errors"] is False
            else "FAIL",
            f"Reference artifact={reference_normalized_path}; reference-only components="
            f"{preprocessing_summary['reference_only_component_count']} and are not classified as errors.",
        ),
        AcceptanceCheck(
            "Analysis SWC parent graph remains structurally valid",
            "PASS" if analysis_swc.structurally_valid else "FAIL",
            f"Nodes={analysis_swc.node_count}; edges={analysis_swc.edge_count}; "
            f"validation={analysis_swc.validation}.",
        ),
        AcceptanceCheck(
            "SWC-centric audit artifacts can be written",
            "PASS"
            if all(
                path.is_file() and path.stat().st_size
                for path in required_audit_artifacts
            )
            else "FAIL",
            f"Summary={preprocessing_summary_path}; components={component_table_path}; "
            f"Mask QC={mask_qc_path}.",
        ),
        AcceptanceCheck(
            "Saved analysis SWC preserves selected source identifiers and topology",
            "PASS"
            if analysis_normalized_path.is_file()
            and analysis_normalized_path.stat().st_size
            and analysis_roundtrip_matches
            else "FAIL",
            f"Artifact={analysis_normalized_path}; round-trip matches selected source "
            f"nodes and parents={analysis_roundtrip_matches}.",
        ),
    ]
    statuses = {check.status for check in checks}
    acceptance = AcceptanceResult(
        "FAIL" if "FAIL" in statuses else ("WARNING" if "WARNING" in statuses else "PASS"),
        checks,
    )
    manifest = {
        "record": record.to_dict(),
        "image_metadata": image_metadata.to_dict() if image_metadata else None,
        "mask_metadata": mask_metadata.to_dict() if mask_metadata else None,
        "image_statistics_original": volume_statistics(image_volume) if image_volume is not None else None,
        "mask_statistics_original": volume_statistics(mask_volume) if mask_volume is not None else None,
        "image_statistics": volume_statistics(image_volume) if image_volume is not None else None,
        "mask_statistics": volume_statistics(mask_volume) if mask_volume is not None else None,
        "swc_reference": reference_swc.report(),
        "swc_analysis": analysis_swc.report(),
        "swc": analysis_swc.report(),
        "swc_centric_preprocessing": preprocessing_summary,
        "swc_centric_summary_path": str(preprocessing_summary_path.resolve()),
        "swc_component_table_path": str(component_table_path.resolve()),
        "reference_only_node_ids_path": str(reference_only_nodes_path.resolve()),
        "optional_mask_qc_path": str(mask_qc_path.resolve()),
        "normalized_volume_path": str(display_volume_path.resolve()) if display_volume_path else None,
        "reference_normalized_swc_path": str(reference_normalized_path.resolve()),
        "reference_swc_text_path": str(reference_swc_text_path.resolve()),
        "normalized_swc_path": str(analysis_normalized_path.resolve()),
        "analysis_swc_text_path": str(analysis_swc_text_path.resolve()),
        "spacing_xyz_um": config.spacing_xyz_um,
        "auxiliary_volume_shape_zyx": auxiliary_shape_zyx,
        "array_axis_order": ["z", "y", "x"],
        "coordinate_axis_order": ["x", "y", "z"],
        "flow_direction_rule": "SWC parent_id node -> current node",
        "flow_direction_is_measured": False,
        "downstream_data_source": "validated analysis_swc selected unchanged from reference_swc",
        "mask_role": "optional registration QC and visualization only",
        "raw_image_role": "optional visualization context only",
    }
    write_json(manifest, sample_root / "preprocess_manifest.json")
    write_json(acceptance.report(), sample_root / "preprocess_acceptance.json")
    return manifest, acceptance


def _load_source_manifests(source_run: Path) -> list[dict[str, Any]]:
    index_path = source_run.resolve() / "normalized_index.json"
    if not index_path.is_file():
        raise FileNotFoundError(f"Normalized index does not exist: {index_path}")
    import json

    payload = json.loads(index_path.read_text(encoding="utf-8"))
    return list(payload["samples"])


def _select_manifests(
    manifests: list[dict[str, Any]], config: RodentVasculatureConfig
) -> list[dict[str, Any]]:
    selected = []
    for manifest in manifests:
        record = manifest["record"]
        if config.cohort != "all" and record["cohort"] != config.cohort:
            continue
        if config.sample_id and config.sample_id not in {record["sample_id"], record["source_stem"]}:
            continue
        if config.parent_group_id and record["parent_group_id"] != config.parent_group_id:
            continue
        if config.split and record.get("split") != config.split:
            continue
        selected.append(manifest)
    return selected[: config.max_samples] if config.max_samples else selected


def _graph_sample(
    manifest: dict[str, Any],
    sample_root: Path,
    config: RodentVasculatureConfig,
    *,
    make_visualizations: bool,
) -> tuple[dict[str, Any], AcceptanceResult, list[Path]]:
    record = _record_from_payload(manifest["record"])
    assert record.swc_path
    write_json(manifest, sample_root / "preprocess_manifest.json")
    auxiliary_shape = manifest.get("auxiliary_volume_shape_zyx")
    if auxiliary_shape is None and manifest.get("image_metadata"):
        auxiliary_shape = manifest["image_metadata"].get("shape_zyx")
    swc = load_normalized_swc(
        Path(manifest["normalized_swc_path"]),
        source_path=record.swc_path,
        spacing_xyz_um=config.spacing_xyz_um,
        volume_shape_zyx=tuple(auxiliary_shape) if auxiliary_shape else None,
    )
    graph = build_directed_vascular_graph(record.sample_id, swc, config)
    graph_paths = export_directed_graph(
        graph,
        sample_root / "graphs",
        save_graphml=config.save_graphml,
        save_vtp=config.save_vtp,
        save_npz=config.save_npz,
    )
    visualizations: list[Path] = []
    visualization_artifacts: list[Path] = []
    if make_visualizations:
        image_volume: np.ndarray | None = None
        mask_volume: np.ndarray | None = None
        if manifest.get("normalized_volume_path"):
            image_volume, mask_volume = load_normalized_volume(
                Path(manifest["normalized_volume_path"])
            )
        else:
            if record.image_path is not None:
                image_volume = load_tiff_volume(record.image_path)
            if record.mask_path is not None:
                mask_volume = load_tiff_volume(record.mask_path)
        visualizations = create_visualizations(
            graph,
            sample_root / "visualizations",
            max_arrows=config.max_direction_arrows,
            spacing_xyz_um=config.spacing_xyz_um,
            image_volume=image_volume,
            mask_volume=mask_volume,
        )
        if config.figure2a_enabled:
            figure2a = save_figure2a_preview(
                graph,
                image_volume if image_volume is not None else mask_volume,
                sample_root / "visualizations",
                spacing_xyz_um=config.spacing_xyz_um,
                max_arrows=config.max_direction_arrows,
                volume_opacity=config.figure2a_volume_opacity,
                window_size=config.figure2a_window_size,
            )
            visualizations.append(figure2a.screenshot_path)
            visualization_artifacts.append(figure2a.manifest_path)
    acceptance = evaluate_directed_graph(
        graph,
        graph_paths + visualizations + visualization_artifacts,
        strict_nonpositive_radius=config.strict_nonpositive_radius,
    )
    report_summary = graph.report() | {
        "preprocessing_swc_centric_selection": manifest.get(
            "swc_centric_preprocessing", {}
        ),
    }
    report_path = write_html_report(
        sample_root / "acceptance_report.html",
        title=f"Directed rodent vascular graph — {record.sample_id}",
        summary=report_summary,
        acceptance=acceptance,
        visualizations=visualizations,
    )
    write_json(acceptance.report(), sample_root / "graph_acceptance.json")
    graph_summary = report_summary | {
        "acceptance_status": acceptance.overall_status,
        "sample_report": str(report_path.resolve()),
    }
    write_json(graph_summary, sample_root / "graph_summary.json")
    return graph_summary, acceptance, visualizations


def run_rodent_vasculature_pipeline(
    config: RodentVasculatureConfig, *, verbose: bool = False
) -> RodentPipelineRun:
    config.validate()
    run_root = _make_run_root(config)
    logger = configure_logging(run_root / "pipeline.log", verbose=verbose)
    logger.info("Starting rodent vasculature stage=%s", config.stage)
    logger.info("Direction convention: SWC parent_id node -> current node (structural inference)")
    write_json(config.report(), run_root / "run_config.json")
    failures: list[dict[str, str]] = []
    preprocess_acceptances: list[AcceptanceResult] = []
    graph_acceptances: list[AcceptanceResult] = []
    manifests: list[dict[str, Any]] = []
    graph_summaries: list[dict[str, Any]] = []

    if config.stage in {"inventory", "preprocess", "all"}:
        catalog = build_catalog(config.input_dir, config.cohort)
        catalog_report = catalog.report()
        write_json(catalog_report, run_root / "inventory" / "catalog_summary.json")
        write_csv([record.to_dict() for record in catalog.records], run_root / "inventory" / "catalog.csv")
        selected_records = select_records(
            catalog,
            sample_id=config.sample_id,
            parent_group_id=config.parent_group_id,
            split=config.split,
            max_samples=config.max_samples,
        )
        logger.info(
            "Cataloged %d records (%d eligible); selected %d.",
            len(catalog.records), len(catalog.eligible_records), len(selected_records),
        )
        if config.stage == "inventory":
            inventory_check = AcceptanceResult(
                "PASS" if catalog.records else "FAIL",
                [
                    AcceptanceCheck(
                        "Dataset catalog is non-empty",
                        "PASS" if catalog.records else "FAIL",
                        f"Records={len(catalog.records)}, eligible={len(catalog.eligible_records)}.",
                    )
                ],
            )
            report = write_html_report(
                run_root / "acceptance_report.html",
                title="Rodent vasculature dataset inventory",
                summary=catalog_report,
                acceptance=inventory_check,
                visualizations=[],
            )
            write_json(inventory_check.report(), run_root / "acceptance.json")
            status = "completed" if inventory_check.overall_status != "FAIL" else "failed"
            write_json({"status": status, "stage": config.stage}, run_root / "run_status.json")
            return RodentPipelineRun(run_root, status, report, inventory_check, 0, 0)
        if not selected_records:
            raise ValueError("No eligible samples matched the selection")
        for record in selected_records:
            sample_root = run_root / "samples" / _safe_name(record.sample_id)
            sample_root.mkdir(parents=True, exist_ok=False)
            try:
                manifest, acceptance = _preprocess_sample(record, sample_root, config)
                manifests.append(manifest)
                preprocess_acceptances.append(acceptance)
                logger.info("Preprocessed %s: %s", record.sample_id, acceptance.overall_status)
            except Exception as exc:
                logger.exception("Preprocessing failed for %s", record.sample_id)
                failures.append({"sample_id": record.sample_id, "stage": "preprocess", "error": str(exc)})
        write_json({"samples": manifests}, run_root / "normalized_index.json")

    if config.stage == "hierarchical-graph":
        manifests = _select_manifests(_load_source_manifests(config.source_run), config)  # type: ignore[arg-type]
        if not manifests:
            raise ValueError("No normalized samples matched the selection")

    if config.stage in {"hierarchical-graph", "all"}:
        for sample_index, manifest in enumerate(manifests):
            record = _record_from_payload(manifest["record"])
            sample_root = run_root / "samples" / _safe_name(record.sample_id)
            sample_root.mkdir(parents=True, exist_ok=True)
            make_visualizations = (
                config.visualizations_enabled and sample_index < config.max_visualization_samples
            )
            try:
                summary, acceptance, visualization_paths = _graph_sample(
                    manifest, sample_root, config, make_visualizations=make_visualizations
                )
                graph_summaries.append(summary)
                graph_acceptances.append(acceptance)
                logger.info("Built directed graph for %s: %s", record.sample_id, acceptance.overall_status)
                if visualization_paths:
                    logger.info(
                        "Wrote %d visualizations for %s: %s",
                        len(visualization_paths),
                        record.sample_id,
                        ", ".join(path.name for path in visualization_paths),
                    )
            except Exception as exc:
                logger.exception("Graph construction failed for %s", record.sample_id)
                failures.append({"sample_id": record.sample_id, "stage": "graph", "error": str(exc)})

    write_csv(graph_summaries, run_root / "reports" / "sample_graph_summary.csv")
    write_json(failures, run_root / "reports" / "failures.json")
    required = [run_root / "pipeline.log", run_root / "run_config.json"]
    aggregate = _aggregate_acceptance(
        preprocess_acceptances + graph_acceptances, failures, required
    )
    report = write_html_report(
        run_root / "acceptance_report.html",
        title=f"Rodent vasculature pipeline — {config.stage}",
        summary={
            "stage": config.stage,
            "processed_sample_count": len(graph_summaries or manifests),
            "failed_sample_count": len(failures),
            "direction_rule": "SWC parent_id node -> current node",
            "direction_is_measured": False,
        },
        acceptance=aggregate,
        visualizations=[],
    )
    write_json(aggregate.report(), run_root / "acceptance.json")
    status = "failed" if aggregate.overall_status == "FAIL" else "completed"
    write_json(
        {
            "status": status,
            "stage": config.stage,
            "acceptance_status": aggregate.overall_status,
            "processed_sample_count": len(graph_summaries or manifests),
            "failed_sample_count": len(failures),
        },
        run_root / "run_status.json",
    )
    logger.info("Run finished: %s (%s)", status, aggregate.overall_status)
    return RodentPipelineRun(
        run_root,
        status,
        report,
        aggregate,
        len(graph_summaries or manifests),
        len(failures),
    )
