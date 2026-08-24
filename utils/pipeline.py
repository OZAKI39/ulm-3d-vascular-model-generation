"""End-to-end orchestration for mesh cleanup, voxelization, and reporting."""

from __future__ import annotations

import importlib.metadata
import logging
import platform
import sys
import time
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator

import numpy as np
import psutil
import vtk

from .config import PipelineConfig
from .io import read_stl, write_binary_stl, write_csv, write_json, write_nifti_mask
from .mesh.cleanup import cleanup_mesh
from .reporting.acceptance import AcceptanceResult, evaluate_acceptance
from .reporting.html_report import write_html_report
from .reporting.logger import configure_logging
from .reporting.output_layout import OutputLayout, create_output_layout
from .reporting.visualization import (
    render_acceptance_dashboard,
    render_mesh_visualizations,
    render_skeleton_visualizations,
    render_voxel_connectivity_visualization,
    render_voxel_visualizations,
)
from .voxel.skeleton import extract_coarse_skeleton
from .voxel.voxelize import voxelize_surface


@dataclass(frozen=True, slots=True)
class PipelineResult:
    run_root: Path
    html_report: Path
    status: str
    acceptance: AcceptanceResult


@contextmanager
def _timed_stage(
    name: str,
    logger: logging.Logger,
    timings: dict[str, float],
) -> Iterator[None]:
    logger.info("Stage started: %s", name)
    started = time.perf_counter()
    try:
        yield
    finally:
        elapsed = time.perf_counter() - started
        timings[name] = elapsed
        logger.info("Stage finished: %s (%.2f seconds)", name, elapsed)


def _package_version(distribution: str) -> str:
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
            "vtk": vtk.vtkVersion.GetVTKVersion(),
            "pyvista": _package_version("pyvista"),
            "scikit-image": _package_version("scikit-image"),
            "nibabel": _package_version("nibabel"),
            "pymeshfix": _package_version("pymeshfix"),
            "matplotlib": _package_version("matplotlib"),
        },
    }


def _write_status(
    layout: OutputLayout,
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


def run_pipeline(config: PipelineConfig, *, verbose: bool = False) -> PipelineResult:
    config.validate()
    layout = create_output_layout(config.output_root, config.input_stl)
    logger = configure_logging(layout.log_file, verbose=verbose)
    timings: dict[str, float] = {}
    started_at = datetime.now(timezone.utc).isoformat()
    _write_status(
        layout,
        status="running",
        started_at=started_at,
        message="Pipeline is running.",
        timings=timings,
    )
    write_json(config.to_dict(), layout.config_file)
    environment_file = layout.reports / "environment.json"
    write_json(_environment_report(), environment_file)

    logger.info("Input STL: %s", config.input_stl)
    logger.info("Run directory: %s", layout.run_root)
    logger.info("Python interpreter: %s", sys.executable)
    logger.info("Physical units: %s; source coordinate system: %s", config.physical_unit, config.coordinate_system)

    try:
        with _timed_stage("read_input_stl", logger, timings):
            input_mesh = read_stl(config.input_stl)
            logger.info(
                "Input surface loaded: %s points, %s cells, %.2f MiB",
                f"{input_mesh.GetNumberOfPoints():,}",
                f"{input_mesh.GetNumberOfCells():,}",
                config.input_stl.stat().st_size / (1024**2),
            )

        with _timed_stage("mesh_cleanup", logger, timings):
            cleanup = cleanup_mesh(input_mesh, config.mesh, logger)
            cleaned_path = layout.meshes / "cleaned.stl"
            removed_path = layout.meshes / "removed_components.stl"
            removed_fragments_path = layout.meshes / "removed_small_fragments.stl"
            removed_islands_path = layout.meshes / "removed_island_networks.stl"
            write_binary_stl(cleanup.cleaned_mesh, cleaned_path)
            if cleanup.removed_mesh is not None and cleanup.removed_mesh.GetNumberOfCells():
                write_binary_stl(cleanup.removed_mesh, removed_path)
            if (
                cleanup.removed_small_fragments_mesh is not None
                and cleanup.removed_small_fragments_mesh.GetNumberOfCells()
            ):
                write_binary_stl(
                    cleanup.removed_small_fragments_mesh, removed_fragments_path
                )
            if (
                cleanup.removed_island_networks_mesh is not None
                and cleanup.removed_island_networks_mesh.GetNumberOfCells()
            ):
                write_binary_stl(
                    cleanup.removed_island_networks_mesh, removed_islands_path
                )
            cleanup_summary = cleanup.summary()
            cleanup_report_path = layout.reports / "cleanup_report.json"
            main_selection_path = layout.reports / "main_network_selection.json"
            components_path = layout.reports / "mesh_components.csv"
            write_json(cleanup_summary, cleanup_report_path)
            write_json(cleanup.main_network_selection.to_dict(), main_selection_path)
            component_rows = [item.to_dict() for item in sorted(cleanup.components, key=lambda x: x.area_rank)]
            write_csv(component_rows, components_path)
            logger.info(
                "Cleaned main network: %s triangles, %d component; removed "
                "%d small fragments and %d island networks (area %.5f%%)",
                f"{cleanup.cleaned_quality.triangle_count:,}",
                cleanup.cleaned_quality.connected_component_count,
                cleanup_summary["small_fragment_count"],
                cleanup_summary["island_network_count"],
                cleanup_summary["removed_surface_area_fraction"] * 100,
            )
            reopened = read_stl(cleaned_path)
            if reopened.GetNumberOfCells() != cleanup.cleaned_quality.triangle_count:
                raise RuntimeError(
                    "Reopened cleaned STL has a different triangle count: "
                    f"{reopened.GetNumberOfCells()} vs {cleanup.cleaned_quality.triangle_count}"
                )

        with _timed_stage("mesh_visualization", logger, timings):
            images = render_mesh_visualizations(
                input_mesh,
                cleanup.cleaned_mesh,
                cleanup.removed_mesh,
                cleanup.removed_small_fragments_mesh,
                cleanup.removed_island_networks_mesh,
                cleanup.components,
                layout.visualizations,
                config.visualization,
            )

        with _timed_stage("voxelization", logger, timings):
            voxel = voxelize_surface(cleanup.cleaned_mesh, config.voxel, logger)
            mask_path = layout.volumes / "voxel_mask.nii.gz"
            voxel_removed_islands_path = layout.volumes / "voxel_removed_islands.nii.gz"
            mask_affine = write_nifti_mask(
                voxel.mask,
                mask_path,
                voxel.origin_lps_um,
                voxel.spacing_um,
                "Uniform vascular lumen voxel mask",
            )
            if voxel.removed_island_voxel_count:
                write_nifti_mask(
                    voxel.removed_islands_mask,
                    voxel_removed_islands_path,
                    voxel.origin_lps_um,
                    voxel.spacing_um,
                    "Voxel islands removed from the final main vascular mask",
                )
            voxel_report_path = layout.reports / "voxel_report.json"
            voxel_connectivity_report_path = (
                layout.reports / "voxel_connectivity_report.json"
            )
            voxel_report = voxel.report()
            write_json(voxel_report, voxel_report_path)
            write_json(
                {
                    "initial_connected_component_count": (
                        voxel.initial_connected_component_count
                    ),
                    "final_connected_component_count": voxel.connected_component_count,
                    "initial_foreground_voxel_count": (
                        voxel.initial_foreground_voxel_count
                    ),
                    "final_foreground_voxel_count": voxel.foreground_voxel_count,
                    "removed_island_voxel_count": voxel.removed_island_voxel_count,
                    "removed_island_fraction": voxel.removed_island_fraction,
                    "component_voxel_counts_top20": (
                        voxel.component_voxel_counts_top20
                    ),
                },
                voxel_connectivity_report_path,
            )
            spatial_metadata_path = layout.volumes / "spatial_metadata.json"
            write_json(
                {
                    "array_axis_order": ["x", "y", "z"],
                    "source_coordinate_system": "LPS",
                    "nifti_world_coordinate_system": "RAS",
                    "origin_lps_um": voxel.origin_lps_um,
                    "spacing_um": voxel.spacing_um,
                    "dimensions_xyz": voxel.dimensions_xyz,
                    "nifti_affine_ras": mask_affine,
                    "lps_to_ras_rule": "RAS coordinates equal (-LPS_x, -LPS_y, LPS_z)",
                },
                spatial_metadata_path,
            )
            if config.visualization.enabled:
                images.extend(render_voxel_visualizations(voxel.mask, layout.visualizations))
                images.append(
                    render_voxel_connectivity_visualization(
                        voxel.mask,
                        voxel.removed_islands_mask,
                        layout.visualizations,
                    )
                )

        if not config.voxel.skeletonize:
            raise ValueError("Step 2 acceptance requires skeletonize=True")
        with _timed_stage("coarse_skeleton", logger, timings):
            skeleton = extract_coarse_skeleton(voxel.mask, voxel.spacing_um, logger)
            skeleton_path = layout.volumes / "coarse_skeleton.nii.gz"
            write_nifti_mask(
                skeleton.skeleton,
                skeleton_path,
                voxel.origin_lps_um,
                voxel.spacing_um,
                "Coarse navigation skeleton; not final CFD geometry",
            )
            skeleton_report_path = layout.reports / "skeleton_report.json"
            skeleton_report = skeleton.report()
            write_json(skeleton_report, skeleton_report_path)
            if config.visualization.enabled:
                images.extend(
                    render_skeleton_visualizations(
                        voxel.mask,
                        skeleton.skeleton,
                        voxel.origin_lps_um,
                        voxel.spacing_um,
                        layout.visualizations,
                        config.visualization,
                    )
                )

        with _timed_stage("acceptance_report", logger, timings):
            required_files = [
                layout.config_file,
                environment_file,
                cleaned_path,
                mask_path,
                skeleton_path,
                cleanup_report_path,
                main_selection_path,
                components_path,
                voxel_report_path,
                voxel_connectivity_report_path,
                skeleton_report_path,
                spatial_metadata_path,
                *images,
            ]
            if voxel_removed_islands_path.is_file():
                required_files.append(voxel_removed_islands_path)
            acceptance = evaluate_acceptance(
                cleanup,
                voxel,
                skeleton,
                required_files,
                config.voxel,
            )
            acceptance_path = layout.reports / "acceptance_summary.json"
            write_json(acceptance.report(), acceptance_path)
            if config.visualization.enabled:
                dashboard = render_acceptance_dashboard(
                    acceptance, layout.visualizations / "acceptance_dashboard.png"
                )
                images.append(dashboard)
            report_files = [
                cleaned_path,
                removed_path,
                removed_fragments_path,
                removed_islands_path,
                mask_path,
                voxel_removed_islands_path,
                skeleton_path,
                cleanup_report_path,
                main_selection_path,
                components_path,
                voxel_report_path,
                voxel_connectivity_report_path,
                skeleton_report_path,
                acceptance_path,
                spatial_metadata_path,
                environment_file,
                layout.config_file,
                layout.log_file,
            ]
            write_html_report(
                layout.html_report,
                input_stl=config.input_stl,
                run_root=layout.run_root,
                acceptance=acceptance,
                cleanup_summary=cleanup_summary,
                voxel_report=voxel_report,
                skeleton_report=skeleton_report,
                images=images,
                files=report_files,
            )

        pipeline_status = {
            "PASS": "success",
            "WARNING": "warning",
            "FAIL": "failed",
        }[acceptance.overall_status]
        message = f"Pipeline completed with acceptance status {acceptance.overall_status}."
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
        return PipelineResult(
            run_root=layout.run_root,
            html_report=layout.html_report,
            status=pipeline_status,
            acceptance=acceptance,
        )
    except Exception as exc:
        logger.exception("Pipeline failed: %s", exc)
        _write_status(
            layout,
            status="failed",
            started_at=started_at,
            message=f"Pipeline raised {type(exc).__name__}: {exc}",
            timings=timings,
        )
        raise
