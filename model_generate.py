"""CLI orchestration for representative SWC ROI to CFD-ready vascular lumen."""

from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

from utils.cfd_lumen import (
    load_cfd_lumen_config,
    load_sampling_rois,
    resolve_sampling_run,
    run_cfd_lumen_batch,
)
from utils.cfd_lumen.diagnostic_pipeline import run_geometry_diagnostics
from utils.cfd_lumen.export import create_run_layout, write_json, write_resolved_config


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT = PROJECT_ROOT / "outputs"
DEFAULT_CONFIG = PROJECT_ROOT / "cfd_lumen_config.yaml"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct source-traceable, watertight CFD lumen surfaces directly from "
            "representative SWC ROI centerlines and radii."
        )
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--sampling-run",
        type=Path,
        default=None,
        help="Sampling run root; defaults to the latest complete outputs/sampling run.",
    )
    parser.add_argument(
        "--rodent-run",
        type=Path,
        default=None,
        help=(
            "Rodent preprocessing run containing the unchanged global SWC model used for "
            "v4 CFD-domain context extension; otherwise a geometry-verified run is resolved."
        ),
    )
    parser.add_argument("--roi-id", default=None, help="Process exactly one candidate or selected ROI ID.")
    diagnosis = parser.add_mutually_exclusive_group()
    diagnosis.add_argument(
        "--diagnose-roi",
        default=None,
        metavar="ROI_ID",
        help="Run the strict v2 root-cause protocol for exactly one saved ROI.",
    )
    diagnosis.add_argument(
        "--diagnose-all",
        action="store_true",
        help="Run the strict v2 root-cause protocol for every selected ROI, serially.",
    )
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument(
        "--selected-only",
        action="store_true",
        help="Process the selected representative ROI manifest (the default).",
    )
    selection.add_argument(
        "--all-selected",
        action="store_true",
        help="Explicitly process all selected representative ROIs.",
    )
    parser.add_argument("--backend", choices=("manifold", "implicit"), default=None)
    parser.add_argument(
        "--surface-version",
        choices=("v5", "v6", "v7", "v8"),
        default="v8",
        help=(
            "Formal surface protocol; v8 localizes the v7 artifact and evaluates "
            "compact-support smooth union only in true junction cores."
        ),
    )
    parser.add_argument(
        "--junction-backend",
        choices=("local_implicit", "legacy_sphere"),
        default=None,
        help="Junction reconstruction; local_implicit is the v3 formal default.",
    )
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--headless", action="store_true", help="Force off-screen rendering.")
    parser.add_argument("--build-volume-mesh", action="store_true")
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-convergence", action="store_true")
    parser.add_argument("--no-visualizations", action="store_true")
    return parser


def _configure_logging(log_path: Path) -> logging.Logger:
    logger = logging.getLogger("cfd_lumen")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(processName)s | %(levelname)s | %(message)s")
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(formatter)
    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(console)
    logger.addHandler(file_handler)
    return logger


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.workers is not None and args.workers < 1:
        parser.error("--workers must be positive")
    diagnosis_mode = args.diagnose_roi is not None or args.diagnose_all
    if diagnosis_mode and args.roi_id is not None:
        parser.error("--roi-id cannot be combined with --diagnose-roi/--diagnose-all")
    if diagnosis_mode and args.workers not in (None, 1):
        parser.error("v2 diagnosis requires --workers 1 for traceable logs")
    if args.headless:
        os.environ["PYVISTA_OFF_SCREEN"] = "true"
    try:
        config = load_cfd_lumen_config(args.config)
        if args.backend is not None:
            config.boolean.backend = args.backend
        config.v6.enabled = args.surface_version in {"v6", "v7", "v8"}
        config.v7.enabled = args.surface_version in {"v7", "v8"}
        config.v8.enabled = args.surface_version == "v8"
        if args.junction_backend is not None:
            config.junction.backend = args.junction_backend
            config.reconstruction.junction_backend = args.junction_backend
        if args.build_volume_mesh:
            config.volume_mesh.enabled = True
        if args.no_convergence:
            config.convergence.enabled = False
        if args.no_visualizations:
            config.output.visualizations = False
        if args.rodent_run is not None:
            config.context_domain.source_rodent_run = str(args.rodent_run.resolve())
        config.validate()
        sampling_run = resolve_sampling_run(args.sampling_run, project_root=PROJECT_ROOT)
        rois = load_sampling_rois(
            sampling_run,
            roi_id=args.diagnose_roi if diagnosis_mode else args.roi_id,
            selected_only=True,
        )
        run_layout = create_run_layout(args.output_dir, args.run_id)
        logger = _configure_logging(run_layout.logs / "model_generate.log")
        write_resolved_config(run_layout.config / "resolved_config.yaml", config)
        write_json(
            run_layout.report / "repository_inspection.json",
            {
                "sampling_run": str(sampling_run),
                "roi_loader": "utils.sampling.sampling_io.load_sampling_display_rois",
                "roi_geometry_archives": str(sampling_run / "roi_library"),
                "selected_manifest": str(sampling_run / "manifests" / "selected_rois.csv"),
                "cut_port_manifest": str(sampling_run / "manifests" / "cut_ports.csv"),
                "source_coordinate_unit": "um",
                "source_radius_unit": "um",
                "segmentation_mask_used": False,
                "available_roi_fields_verified": [
                    "local_node_ids",
                    "local_node_positions_um",
                    "local_node_radius_um",
                    "local_edges",
                    "local_node_global_ids",
                    "local_edge_global_ids",
                    "cut_port_positions_um",
                    "cut_port_radius_um",
                    "cut_port_global_edge_ids",
                    "cut_port_boundary_faces",
                    "cut_port_boundary_roles",
                    "bbox_min_um",
                    "bbox_max_um",
                ],
            },
        )
        logger.info("Python interpreter: %s", sys.executable)
        logger.info("Input sampling run: %s", sampling_run)
        logger.info("ROI IDs: %s", ", ".join(roi.roi_id for roi in rois))
        logger.info("Backend: %s", config.boolean.backend)
        logger.info("Surface continuity protocol: %s", args.surface_version)
        if diagnosis_mode:
            summaries = run_geometry_diagnostics(rois, config, run_layout)
            logger.info("Run directory: %s", run_layout.run_root)
            logger.info("Strict v2 diagnosis complete for %d ROI(s)", len(summaries))
            return 0
        results = run_cfd_lumen_batch(
            rois,
            config,
            run_layout,
            sampling_run,
            workers=args.workers,
        )
        passed = sum(result.status == "PASS" for result in results)
        logger.info("Run directory: %s", run_layout.run_root)
        logger.info("Run result: %d/%d ROI(s) PASS", passed, len(results))
        return 0 if passed == len(results) else 2
    except Exception as exc:
        logging.getLogger("cfd_lumen").exception("Model generation failed: %s", exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
