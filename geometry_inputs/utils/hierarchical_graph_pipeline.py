"""Independent Step 3 pipeline using artifacts from a completed preprocessing run."""

from __future__ import annotations

import hashlib
import importlib.metadata
import json
import logging
import platform
import sys
import time
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

import nibabel as nib
import networkx as nx
import numpy as np
import psutil
import vtk

from .config import HierarchicalGraphConfig
from .graph.export import export_hierarchical_graph, verify_graph_exports
from .graph.extraction import build_hierarchical_graph
from .graph.validation import evaluate_hierarchical_graph_acceptance
from .graph.visualization import render_hierarchical_graph_visualizations
from .io import write_json, write_nifti_mask
from .reporting.acceptance import AcceptanceResult
from .reporting.hierarchical_graph_report import write_hierarchical_graph_html_report
from .reporting.logger import configure_logging
from .reporting.output_layout import (
    HierarchicalGraphOutputLayout,
    create_hierarchical_graph_output_layout,
)
from .reporting.visualization import render_acceptance_dashboard


@dataclass(frozen=True, slots=True)
class HierarchicalGraphPipelineResult:
    run_root: Path
    html_report: Path
    status: str
    acceptance: AcceptanceResult


@contextmanager
def _timed_stage(
    name: str, logger: logging.Logger, timings: dict[str, float]
) -> Iterator[None]:
    logger.info("Stage started: %s", name)
    started = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - started
        timings[name] = elapsed
        logger.info("Stage finished: %s (%.2f seconds)", name, elapsed)


def _version(distribution: str) -> str:
    try:
        return importlib.metadata.version(distribution)
    except importlib.metadata.PackageNotFoundError:
        return "not installed"


def _environment_report() -> dict[str, object]:
    memory = psutil.virtual_memory()
    return {
        "python_executable": sys.executable,
        "python_version": sys.version,
        "platform": platform.platform(),
        "logical_cpu_count": psutil.cpu_count(logical=True),
        "physical_memory_total_gib": memory.total / (1024**3),
        "physical_memory_available_gib_at_start": memory.available / (1024**3),
        "packages": {
            "numpy": np.__version__,
            "scipy": _version("scipy"),
            "networkx": nx.__version__,
            "nibabel": nib.__version__,
            "vtk": vtk.vtkVersion.GetVTKVersion(),
            "matplotlib": _version("matplotlib"),
        },
    }


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_status(
    layout: HierarchicalGraphOutputLayout,
    *,
    status: str,
    started_at: str,
    message: str,
    timings: dict[str, float],
    acceptance_status: str | None = None,
) -> None:
    payload = {
        "status": status,
        "acceptance_status": acceptance_status,
        "message": message,
        "started_at": started_at,
        "updated_at": datetime.now(timezone.utc).isoformat(),
        "source_stage": "Step 2 coarse skeleton",
        "current_stage": "Step 3 hierarchical vascular representation",
        "run_directory": str(layout.run_root),
        "timings_seconds": timings,
    }
    write_json(payload, layout.status_file)
    write_json(
        {
            "latest_run": str(layout.run_root),
            "latest_run_relative": layout.run_root.relative_to(layout.sample_root).as_posix(),
            "status": status,
            "acceptance_status": acceptance_status,
            "updated_at": payload["updated_at"],
        },
        layout.latest_file,
    )


def _validate_source_run(source_run: Path) -> dict[str, Path]:
    source_run = source_run.expanduser().resolve()
    if not source_run.is_dir():
        raise FileNotFoundError(f"Source run directory does not exist: {source_run}")
    paths = {
        "skeleton": source_run / "volumes" / "coarse_skeleton.nii.gz",
        "mask": source_run / "volumes" / "voxel_mask.nii.gz",
        "spatial_metadata": source_run / "volumes" / "spatial_metadata.json",
        "run_config": source_run / "run_config.json",
        "skeleton_report": source_run / "reports" / "skeleton_report.json",
    }
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Source run is missing required Step 2 files: {missing}")
    return paths


def run_hierarchical_graph_pipeline(
    source_run: Path,
    config: HierarchicalGraphConfig,
    *,
    verbose: bool = False,
) -> HierarchicalGraphPipelineResult:
    config.validate()
    source_run = Path(source_run).expanduser().resolve()
    source_files = _validate_source_run(source_run)
    layout = create_hierarchical_graph_output_layout(source_run)
    logger = configure_logging(layout.log_file, verbose=verbose)
    timings: dict[str, float] = {}
    started_at = datetime.now(timezone.utc).isoformat()
    _write_status(
        layout,
        status="running",
        started_at=started_at,
        message="Step 3 hierarchical graph pipeline is running.",
        timings=timings,
    )
    write_json(
        {
            "stage": "hierarchical-graph",
            "source_run": str(source_run),
            "hierarchical_graph": asdict(config),
        },
        layout.config_file,
    )
    environment_path = layout.reports / "environment.json"
    source_path = layout.reports / "source_provenance.json"
    write_json(_environment_report(), environment_path)
    write_json(
        {
            "source_run": str(source_run),
            "source_files": {
                name: {
                    "path": str(path),
                    "size_bytes": path.stat().st_size,
                    "sha256": _sha256(path),
                }
                for name, path in source_files.items()
            },
        },
        source_path,
    )
    logger.info("Source run: %s", source_run)
    logger.info("Output run: %s", layout.run_root)
    logger.info("Python interpreter: %s", sys.executable)

    try:
        with _timed_stage("load_step2_artifacts", logger, timings):
            skeleton_image = nib.load(str(source_files["skeleton"]))
            mask_image = nib.load(str(source_files["mask"]))
            skeleton = np.asarray(skeleton_image.dataobj) > 0
            mask = np.asarray(mask_image.dataobj) > 0
            metadata = json.loads(
                source_files["spatial_metadata"].read_text(encoding="utf-8")
            )
            origin_lps_um = tuple(float(value) for value in metadata["origin_lps_um"])
            spacing_um = tuple(float(value) for value in metadata["spacing_um"])
            if skeleton.shape != mask.shape:
                raise ValueError(
                    f"Skeleton shape {skeleton.shape} differs from mask shape {mask.shape}"
                )
            if tuple(metadata["dimensions_xyz"]) != skeleton.shape:
                raise ValueError("spatial_metadata dimensions do not match NIfTI arrays")
            logger.info(
                "Loaded skeleton: shape=%s, voxels=%s, spacing=%s um",
                skeleton.shape,
                f"{np.count_nonzero(skeleton):,}",
                spacing_um,
            )

        with _timed_stage("hierarchical_graph_extraction", logger, timings):
            graph_result = build_hierarchical_graph(
                skeleton,
                mask,
                origin_lps_um,
                spacing_um,
                config,
                logger,
            )

        with _timed_stage("graph_export", logger, timings):
            exported = export_hierarchical_graph(
                graph_result,
                layout.graphs,
                layout.tables,
                save_graphml=config.save_graphml,
                save_vtp=config.save_vtp,
                save_npz=config.save_npz,
            )
            reconstructed_path = layout.volumes / "graph_reconstructed_skeleton.nii.gz"
            difference_path = layout.volumes / "skeleton_reconstruction_difference.nii.gz"
            write_nifti_mask(
                graph_result.reconstructed_skeleton,
                reconstructed_path,
                origin_lps_um,
                spacing_um,
                "Skeleton reconstructed losslessly from hierarchical graph",
            )
            write_nifti_mask(
                graph_result.source_skeleton ^ graph_result.reconstructed_skeleton,
                difference_path,
                origin_lps_um,
                spacing_um,
                "Difference between source and graph-reconstructed skeleton",
            )
            export_errors = verify_graph_exports(
                exported, len(graph_result.nodes), len(graph_result.branches)
            )
            if export_errors:
                logger.warning("Export verification issues: %s", export_errors)

        with _timed_stage("reports", logger, timings):
            graph_report = graph_result.report()
            graph_report_path = layout.reports / "hierarchical_graph_report.json"
            topology_report_path = layout.reports / "topology_report.json"
            geometry_report_path = layout.reports / "branch_geometry_report.json"
            cycle_report_path = layout.reports / "cycle_report.json"
            fidelity_report_path = layout.reports / "representation_fidelity_report.json"
            write_json(graph_report, graph_report_path)
            write_json(
                {
                    "flow_direction_known": False,
                    "parent_daughter_fields_available": False,
                    "node_count": len(graph_result.nodes),
                    "branch_count": len(graph_result.branches),
                    "graph_component_count": graph_result.graph_component_count,
                    "nodes": [item.summary() for item in graph_result.nodes],
                },
                topology_report_path,
            )
            write_json(
                {
                    "quality_level": "coarse navigation",
                    "radius_source": graph_result.radius_source,
                    "approved_for_cfd": False,
                    "approved_as_final_geometry_training_truth": False,
                    "branch_length_um": graph_report["branch_length_um"],
                    "branches": [item.summary() for item in graph_result.branches],
                    "junction_local_geometry": [
                        item.to_dict() for item in graph_result.junctions
                    ],
                },
                geometry_report_path,
            )
            write_json(
                {
                    "cycle_rank": graph_result.cycle_rank,
                    "stored_cycle_basis_count": len(graph_result.cycles),
                    "cycles": [item.to_dict() for item in graph_result.cycles],
                },
                cycle_report_path,
            )
            write_json(
                {
                    "source_skeleton_voxel_count": graph_result.skeleton_voxel_count,
                    "represented_voxel_count": graph_result.represented_voxel_count,
                    "missing_voxel_count": graph_result.missing_voxel_count,
                    "extra_voxel_count": graph_result.extra_voxel_count,
                    "duplicate_interior_voxel_count": (
                        graph_result.duplicate_interior_voxel_count
                    ),
                    "source_component_count": graph_result.skeleton_component_count,
                    "graph_component_count": graph_result.graph_component_count,
                    "cycle_rank": graph_result.cycle_rank,
                    "resolution_sensitivity": {
                        "status": "not_evaluated",
                        "reason": "Only one source run was supplied.",
                    },
                },
                fidelity_report_path,
            )

        with _timed_stage("graph_visualization", logger, timings):
            images = render_hierarchical_graph_visualizations(
                graph_result, config, layout.visualizations
            )

        with _timed_stage("graph_acceptance", logger, timings):
            required_files = [
                layout.config_file,
                environment_path,
                source_path,
                reconstructed_path,
                difference_path,
                graph_report_path,
                topology_report_path,
                geometry_report_path,
                cycle_report_path,
                fidelity_report_path,
                *exported,
                *images,
            ]
            acceptance = evaluate_hierarchical_graph_acceptance(
                graph_result, config, required_files, export_errors
            )
            acceptance_path = layout.reports / "graph_acceptance_summary.json"
            write_json(acceptance.report(), acceptance_path)
            dashboard = render_acceptance_dashboard(
                acceptance, layout.visualizations / "graph_acceptance_dashboard.png"
            )
            images.append(dashboard)
            report_files = [
                *exported,
                reconstructed_path,
                difference_path,
                graph_report_path,
                topology_report_path,
                geometry_report_path,
                cycle_report_path,
                fidelity_report_path,
                acceptance_path,
                source_path,
                environment_path,
                layout.config_file,
                layout.log_file,
            ]
            write_hierarchical_graph_html_report(
                layout.html_report,
                source_run=source_run,
                run_root=layout.run_root,
                acceptance=acceptance,
                report=graph_report,
                images=images,
                files=report_files,
            )

        pipeline_status = {
            "PASS": "success",
            "WARNING": "warning",
            "FAIL": "failed",
        }[acceptance.overall_status]
        message = (
            "Step 3 completed with acceptance status "
            f"{acceptance.overall_status}."
        )
        _write_status(
            layout,
            status=pipeline_status,
            acceptance_status=acceptance.overall_status,
            started_at=started_at,
            message=message,
            timings=timings,
        )
        logger.info(message)
        logger.info("Acceptance report: %s", layout.html_report)
        return HierarchicalGraphPipelineResult(
            run_root=layout.run_root,
            html_report=layout.html_report,
            status=pipeline_status,
            acceptance=acceptance,
        )
    except Exception as exc:
        logger.exception("Step 3 pipeline failed: %s", exc)
        _write_status(
            layout,
            status="failed",
            started_at=started_at,
            message=f"Step 3 raised {type(exc).__name__}: {exc}",
            timings=timings,
        )
        raise
