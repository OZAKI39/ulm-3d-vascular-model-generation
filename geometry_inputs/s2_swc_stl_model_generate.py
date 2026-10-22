"""YAML-driven orchestration for the validated Ultraliser surface workflow.

With no command-line argument, this entry point reads
``configs/swc_stl_model_generate.yaml``.  A different YAML file may be
supplied as the sole positional argument.  ROI selection, reconstruction,
quality-control, path, and runtime hyperparameters are all defined in YAML.
"""

from __future__ import annotations

import sys
from pathlib import Path

from utils.cfd_lumen import (
    SWCSTLRunConfig,
    load_sampling_rois,
    load_swc_stl_yaml_config,
    resolve_sampling_run,
)
from utils.cfd_lumen.ultraliser_backend import (
    UltraliserBackendError,
    run_ultraliser_reconstruction,
)


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "swc_stl_model_generate.yaml"
USAGE = (
    "Usage:\n"
    "  python swc_stl_model_generate.py\n"
    "  python swc_stl_model_generate.py <config.yaml>\n\n"
    "All ROI-selection, reconstruction, and QC hyperparameters are defined in YAML.\n"
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


def _select_roi(settings: SWCSTLRunConfig, sampling_run: Path):
    if settings.roi_id is not None:
        matches = load_sampling_rois(
            sampling_run,
            roi_id=settings.roi_id,
            selected_only=True,
        )
        selector = settings.roi_id
    else:
        selected = load_sampling_rois(sampling_run, selected_only=True)
        matches = [
            roi
            for roi in selected
            if int(roi.anchor_id) == int(settings.roi_anchor)
        ]
        selector = f"anchor {settings.roi_anchor}"
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one saved ROI for {selector}, found {len(matches)}")
    return matches[0]


def main(argv: list[str] | None = None) -> int:
    arguments = list(sys.argv[1:] if argv is None else argv)
    try:
        config_path = _configuration_path(arguments)
        if config_path is None:
            return 0
        settings = load_swc_stl_yaml_config(config_path, project_root=PROJECT_ROOT)
        sampling_run = resolve_sampling_run(
            settings.sampling_run,
            project_root=PROJECT_ROOT,
        )
        roi = _select_roi(settings, sampling_run)
        summary = run_ultraliser_reconstruction(
            roi,
            settings.lumen,
            output_root=settings.output_root,
            run_id=settings.run_id,
            ultraliser_root=settings.ultraliser_root,
            executable_path=settings.ultraliser_executable,
            source_config_path=settings.source_path,
        )
        print(f"Configuration: {settings.source_path}")
        print(f"Sampling run: {sampling_run}")
        print(f"Selected ROI: {roi.roi_id}")
        print(f"Run directory: {summary['run_root']}")
        print(f"Saved source configuration: {summary['source_configuration']}")
        print(f"Status: {summary['status']}")
        return 0 if summary["status"] == "PASS" else 2
    except (UltraliserBackendError, FileNotFoundError, KeyError, OSError, TypeError, ValueError) as exc:
        print(f"Model generation failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
