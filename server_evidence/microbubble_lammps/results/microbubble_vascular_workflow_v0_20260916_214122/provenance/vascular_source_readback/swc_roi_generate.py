"""YAML-driven SWC preprocessing and representative connected-ROI generation.

With no command-line argument, this entry point reads
``configs/swc_roi_generate.yaml``.  A different YAML file may be supplied as
the sole positional argument.  All scientific and runtime hyperparameters
live in YAML; this file only orchestrates validated processing stages.

Orange arrows encode the SWC parent-to-current relation and do not claim a
measured blood-flow direction.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

from utils.rodent_vasculature import run_rodent_vasculature_pipeline
from utils.rodent_vasculature.interactive import show_saved_run
from utils.sampling.pipeline import run_sampling_from_rodent_run
from utils.swc_roi_yaml_config import SWCROIRunConfig, load_swc_roi_yaml_config


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "swc_roi_generate.yaml"
USAGE = (
    "Usage:\n"
    "  python swc_roi_generate.py\n"
    "  python swc_roi_generate.py <config.yaml>\n\n"
    "All processing and sampling hyperparameters are defined in YAML.\n"
    f"Default configuration: {DEFAULT_CONFIG}"
)


def _configuration_path(argv: list[str]) -> Path | None:
    if not argv:
        return DEFAULT_CONFIG
    if len(argv) == 1 and argv[0] in {"-h", "--help"}:
        print(USAGE)
        return None
    if len(argv) != 1 or argv[0].startswith("-"):
        raise ValueError("Expected no argument or one YAML configuration path")
    candidate = Path(argv[0]).expanduser()
    return candidate.resolve() if candidate.is_absolute() else (Path.cwd() / candidate).resolve()


def _copy_source_configuration(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)
    return destination


def _load_settings(argv: list[str]) -> SWCROIRunConfig | None:
    config_path = _configuration_path(argv)
    if config_path is None:
        return None
    return load_swc_roi_yaml_config(config_path, project_root=PROJECT_ROOT)


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    try:
        settings = _load_settings(arguments)
    except (FileNotFoundError, OSError, TypeError, ValueError) as exc:
        print(f"ERROR: invalid SWC/ROI YAML configuration: {exc}", file=sys.stderr)
        print(USAGE, file=sys.stderr)
        return 1
    if settings is None:
        return 0

    config = settings.rodent
    graph_stage = config.stage in {"all", "hierarchical-graph"}
    try:
        run = run_rodent_vasculature_pipeline(config, verbose=settings.verbose)
        copied_config = _copy_source_configuration(
            settings.source_path,
            run.run_root / "source_swc_roi_generate.yaml",
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    print(f"Configuration: {settings.source_path}")
    print(f"Saved source configuration: {copied_config}")
    print(f"Run directory: {run.run_root}")
    print(f"Acceptance report: {run.html_report}")
    print(f"Acceptance status: {run.acceptance.overall_status}")

    sampling_run = None
    if run.status != "failed" and graph_stage and settings.sampling_enabled:
        try:
            sampling_run = run_sampling_from_rodent_run(
                run.run_root,
                settings.sampling,
                verbose=settings.verbose,
            )
            sampling_source_config = _copy_source_configuration(
                settings.source_path,
                sampling_run.run_root / "config" / "source_swc_roi_generate.yaml",
            )
        except Exception as exc:
            print(f"ERROR: connected-ROI sampling failed: {exc}", file=sys.stderr)
            return 1
        print(f"Sampling run directory: {sampling_run.run_root}")
        print(f"Sampling source configuration: {sampling_source_config}")
        print(f"Sampling summary: {sampling_run.summary_path}")
        print(f"Sampling status: {sampling_run.status}")
        if sampling_run.status == "FAIL":
            return 2

    sampling_gui_preview = (
        sampling_run.run_root / "figures" / "interactive_sampling_layer_preview.png"
        if sampling_run is not None
        else None
    )
    if (
        run.status != "failed"
        and graph_stage
        and config.figure2a_enabled
        and not settings.show_gui
        and sampling_run is not None
    ):
        try:
            show_saved_run(
                run.run_root,
                sample_id=config.sample_id,
                max_arrows=config.max_direction_arrows,
                volume_opacity=config.figure2a_volume_opacity,
                window_size=config.figure2a_window_size,
                sampling_run_root=sampling_run.run_root,
                screenshot_path=sampling_gui_preview,
                show=False,
            )
            print(f"Sampling GUI preview: {sampling_gui_preview}")
        except Exception as exc:
            print(f"ERROR: sampling GUI preview could not be written: {exc}", file=sys.stderr)
            return 1

    if (
        run.status != "failed"
        and graph_stage
        and config.figure2a_enabled
        and settings.show_gui
    ):
        print("Opening interactive Figure 2(a)-style window; close it to finish the program.")
        if sampling_run is not None:
            print(
                "GUI controls: R/S=selected representative ROIs, A=all candidate ROIs, "
                "C=next ROI cluster; left-click a box to inspect it."
            )
        try:
            show_saved_run(
                run.run_root,
                sample_id=config.sample_id,
                max_arrows=config.max_direction_arrows,
                volume_opacity=config.figure2a_volume_opacity,
                window_size=config.figure2a_window_size,
                sampling_run_root=sampling_run.run_root if sampling_run else None,
                screenshot_path=sampling_gui_preview,
            )
        except Exception as exc:
            print(f"ERROR: interactive window could not be opened: {exc}", file=sys.stderr)
            return 1
    return 2 if run.status == "failed" else 0


if __name__ == "__main__":
    raise SystemExit(main())
