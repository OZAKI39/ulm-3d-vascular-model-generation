"""Create a locally extended, tagged CFD surface from saved PASS inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from utils.cfd_surface_prepare.config import load_surface_prepare_config
from utils.cfd_surface_prepare.pipeline import print_result, run_surface_prepare


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "cfd_surface_prepare.yaml"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Locally flatten and extend saved CFD boundaries without modifying the original "
            "Ultraliser surface or creating a volume mesh."
        )
    )
    parser.add_argument(
        "config",
        nargs="?",
        default=DEFAULT_CONFIG,
        type=Path,
        help=f"strict YAML configuration (default: {DEFAULT_CONFIG})",
    )
    return parser


def main() -> int:
    args = _parser().parse_args()
    try:
        config = load_surface_prepare_config(args.config, project_root=PROJECT_ROOT)
        result = run_surface_prepare(config, project_root=PROJECT_ROOT)
        print_result(result)
        return 0 if result.status == "CFD_SURFACE_PREPARE_PASS_PENDING_MANUAL_REVIEW" else 2
    except Exception as error:
        print(f"CFD surface preparation failed: {error}")
        print("Final status: CFD_SURFACE_PREPARE_FAILED")
        print("NEXT: REVIEW CFD SURFACE PREPARATION FAILURE")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
