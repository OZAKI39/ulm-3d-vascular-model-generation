"""Orchestration for the directed Schmid Step 1-3 equivalent workflow."""

from __future__ import annotations

import logging
import platform
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterator

import matplotlib
import networkx as nx
import nibabel as nib
import numpy as np
import scipy
import vtk

from ..io import write_json
from ..reporting.acceptance import AcceptanceResult
from ..reporting.logger import configure_logging
from .cleanup import clean_schmid_input
from .config import SchmidPKLConfig
from .export import (
    export_graph_artifacts,
    write_acceptance_html,
    write_preview_volume,
)
from .graph_builder import build_directed_hierarchical_graph
from .layout import SchmidOutputLayout, create_schmid_output_layout
from .loader import load_schmid_input
from .model import DirectedGraphResult
from .validation import evaluate_schmid_acceptance
from .visualization import render_acceptance_dashboard, render_schmid_visualizations


@dataclass(frozen=True, slots=True)
class SchmidPipelineRun:
    run_root: Path
    html_report: Path
    status: str
    acceptance: AcceptanceResult
    result: DirectedGraphResult


@contextmanager
def _timed_stage(
    name: str, logger: logging.Logger, timings: dict[str, float]
) -> Iterator[None]:
    logger.info("Starting: %s", name)
    started = time.perf_counter()
    try:
        yield
    finally:
        duration = time.perf_counter() - started
        timings[name] = duration
        logger.info("Finished: %s (%.3f s)", name, duration)


def _status_payload(
    layout: SchmidOutputLayout,
    status: str,
    message: str,
    started_at: str,
    timings: dict[str, float],
    acceptance_status: str | None = None,
) -> dict[str, object]:
    return {
        "status": status,
        "acceptance_status": acceptance_status,
        "message": message,
        "started_at": started_at,
        "updated_at": datetime.now().isoformat(timespec="seconds"),
        "source_stage": "Schmid PKL dictionaries",
        "current_stage": "directed Step 1-3 equivalent",
        "run_directory": str(layout.run_root.resolve()),
        "timings_seconds": timings,
    }


def _direction_report(result: DirectedGraphResult) -> dict[str, object]:
    source = result.cleanup.source
    raw_known = [edge for edge in result.cleanup.edges if edge.direction_status == "known"]
    raw_unresolved = [edge for edge in result.cleanup.edges if edge.direction_status != "known"]
    branch_known = [branch for branch in result.branches if branch.direction_status == "known"]
    branch_unresolved = [branch for branch in result.branches if branch.direction_status != "known"]
    return {
        "method": "higher endpoint pressure to lower endpoint pressure",
        "tuple_order_used_as_direction": False,
        "flow_values_interpreted_as": "non-negative magnitude",
        "raw_edge_known_count": len(raw_known),
        "raw_edge_unresolved_count": len(raw_unresolved),
        "raw_edge_unresolved_ids": [edge.edge_id for edge in raw_unresolved],
        "branch_known_count": len(branch_known),
        "branch_unresolved_count": len(branch_unresolved),
        "branch_unresolved_ids": [branch.branch_id for branch in branch_unresolved],
        "pressure_range_mmhg": [
            float(np.min(source.pressure_mmhg)),
            float(np.max(source.pressure_mmhg)),
        ],
        "directed_graph_is_acyclic": result.directed_is_acyclic,
        "important_interpretation": (
            "These are directions from the published simulation fields, not direct in-vivo "
            "measurements. Multiple parents and children are retained."
        ),
    }


def _flow_report(result: DirectedGraphResult) -> dict[str, object]:
    eligible = [
        row
        for row in result.flow_conservation
        if not row["is_pressure_boundary"]
        and row["unresolved_incident_edge_count"] == 0
        and row["incoming_flow_um3_per_ms"] > 0
        and row["outgoing_flow_um3_per_ms"] > 0
    ]
    relative = np.asarray([row["relative_imbalance"] for row in eligible], dtype=float)
    return {
        "eligible_internal_node_count": len(eligible),
        "median_relative_imbalance": float(np.median(relative)) if len(relative) else None,
        "p95_relative_imbalance": float(np.quantile(relative, 0.95)) if len(relative) else None,
        "p99_relative_imbalance": float(np.quantile(relative, 0.99)) if len(relative) else None,
        "maximum_relative_imbalance": float(np.max(relative)) if len(relative) else None,
        "table": "../tables/flow_conservation.csv",
    }


def _write_reports(
    result: DirectedGraphResult,
    layout: SchmidOutputLayout,
    preview_report: dict[str, object],
) -> list[Path]:
    reports = {
        "cleanup_report.json": result.cleanup.report(),
        "directed_hierarchical_graph_report.json": result.report(),
        "direction_report.json": _direction_report(result),
        "flow_conservation_report.json": _flow_report(result),
        "cycle_report.json": {
            "cycle_rank": (
                result.all_connectivity_graph.number_of_edges()
                - result.all_connectivity_graph.number_of_nodes()
                + result.weak_component_count
            ),
            "stored_cycle_basis_count": len(result.cycles),
            "known_direction_graph_is_acyclic": result.directed_is_acyclic,
            "cycles": [item.to_dict() for item in result.cycles],
        },
        "source_provenance.json": {
            "source_directory": str(
                Path(result.cleanup.source.source_files["verticesDict.pkl"]["path"]).parent
            ),
            "source_files": result.cleanup.source.source_files,
            "coordinate_statement": (
                "Coordinates are retained as published source XYZ micrometers. No LPS/RAS "
                "anatomical orientation is asserted."
            ),
        },
        "preview_volume_report.json": preview_report,
        "environment.json": {
            "python": sys.version,
            "python_executable": sys.executable,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "scipy": scipy.__version__,
            "networkx": nx.__version__,
            "matplotlib": matplotlib.__version__,
            "nibabel": nib.__version__,
            "vtk": vtk.vtkVersion.GetVTKVersion(),
        },
    }
    output: list[Path] = []
    for name, payload in reports.items():
        path = layout.reports / name
        write_json(payload, path)
        output.append(path)
    return output


def run_schmid_pkl_pipeline(
    config: SchmidPKLConfig, *, verbose: bool = False
) -> SchmidPipelineRun:
    config.validate()
    layout = create_schmid_output_layout(config.output_root, config.input_dir)
    logger = configure_logging(layout.log_file, verbose=verbose)
    started_at = datetime.now().isoformat(timespec="seconds")
    timings: dict[str, float] = {}
    write_json(
        {"stage": "directed-schmid-step1-step3", "schmid_pkl": config.report()},
        layout.config_file,
    )
    write_json(
        _status_payload(layout, "running", "Pipeline started.", started_at, timings),
        layout.status_file,
    )
    logger.info("Input directory: %s", config.input_dir.resolve())
    logger.info("Output run: %s", layout.run_root.resolve())
    logger.info("Python interpreter: %s", sys.executable)

    try:
        with _timed_stage("load_and_validate_external_pickles", logger, timings):
            source = load_schmid_input(config.input_dir)
            logger.info("Loaded %d vertices and %d edges", source.vertex_count, source.edge_count)
            logger.info("RBC trajectories intentionally not loaded for Step 1-3")

        with _timed_stage("step1_graph_record_cleanup", logger, timings):
            cleanup = clean_schmid_input(source, config, logger)

        with _timed_stage("step2_centerline_standardization", logger, timings):
            # Geometry is measured while branches are assembled. Raw source points are always retained.
            logger.info(
                "Using published centerlines and diameter profiles directly; no STL voxelization is performed"
            )

        with _timed_stage("step3_directed_hierarchical_graph", logger, timings):
            result = build_directed_hierarchical_graph(cleanup, config, logger)

        with _timed_stage("portable_exports", logger, timings):
            exported = export_graph_artifacts(result, layout)

        preview_report: dict[str, object] = {
            "requested": config.write_preview_volume,
            "written": False,
            "used_for_graph_construction": False,
            "represents_original_lumen_surface": False,
        }
        preview_path: Path | None = None
        if config.write_preview_volume:
            with _timed_stage("derived_preview_volume", logger, timings):
                preview_path, preview_report = write_preview_volume(result, config, layout)
                if preview_path is None:
                    logger.warning("Preview volume skipped: %s", preview_report.get("reason"))
                else:
                    logger.info("Wrote derived centerline preview: %s", preview_path)
                    exported.append(preview_path)

        with _timed_stage("reports_before_acceptance", logger, timings):
            report_files = _write_reports(result, layout, preview_report)
            exported.extend(report_files)

        images: list[Path] = []
        if config.visualizations_enabled:
            with _timed_stage("visualizations", logger, timings):
                images = render_schmid_visualizations(result, config, layout.visualizations)

        required = [
            layout.graphs / "directed_hierarchical_vascular_graph.json",
            layout.graphs / "branch_geometry.npz",
            layout.graphs / "directed_junction_branch_graph.graphml",
            layout.graphs / "directed_branch_as_node_graph.graphml",
            layout.graphs / "all_connectivity_graph.graphml",
            layout.graphs / "branch_centerlines.vtp",
            layout.graphs / "vascular_nodes.vtp",
            layout.tables / "nodes.csv",
            layout.tables / "branches.csv",
            layout.tables / "parent_child_relations.csv",
            layout.tables / "cleanup_decisions.csv",
            layout.reports / "direction_report.json",
            layout.reports / "source_provenance.json",
        ]
        with _timed_stage("acceptance", logger, timings):
            acceptance = evaluate_schmid_acceptance(result, config, required)
            acceptance_path = layout.reports / "directed_graph_acceptance_summary.json"
            write_json(acceptance.report(), acceptance_path)
            exported.append(acceptance_path)
            if config.visualizations_enabled:
                images.append(render_acceptance_dashboard(acceptance, layout.visualizations))
            write_acceptance_html(layout, result, acceptance, images, exported + images)

        status = "failed" if acceptance.overall_status == "FAIL" else "completed"
        message = (
            "Directed Schmid graph pipeline completed."
            if status == "completed"
            else "Pipeline completed, but one or more acceptance checks failed."
        )
        write_json(
            _status_payload(
                layout,
                status,
                message,
                started_at,
                timings,
                acceptance.overall_status,
            ),
            layout.status_file,
        )
        write_json(
            {
                "run_directory": str(layout.run_root.resolve()),
                "status": status,
                "acceptance_status": acceptance.overall_status,
                "updated_at": datetime.now().isoformat(timespec="seconds"),
            },
            layout.latest_file,
        )
        logger.info("Acceptance status: %s", acceptance.overall_status)
        logger.info("Acceptance report: %s", layout.html_report)
        return SchmidPipelineRun(
            run_root=layout.run_root,
            html_report=layout.html_report,
            status=status,
            acceptance=acceptance,
            result=result,
        )
    except Exception as exc:
        logger.exception("Directed Schmid pipeline failed: %s", exc)
        write_json(
            _status_payload(
                layout,
                "failed",
                f"{type(exc).__name__}: {exc}",
                started_at,
                timings,
                "FAIL",
            ),
            layout.status_file,
        )
        raise
