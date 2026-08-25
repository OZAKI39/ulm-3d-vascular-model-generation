"""Prepare global-1D-derived boundary data for a saved vascular ROI."""

from __future__ import annotations

import argparse
from pathlib import Path

from utils.cfd_preprocess.config import load_cfd_preprocess_config
from utils.cfd_preprocess.pipeline import print_result, run_cfd_preprocess


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "cfd_preprocess.yaml"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Transfer a sparse global 1D baseline solution to saved ROI CUT_PORTs. "
            "This command does not run 3D CFD or regenerate ROI/STL data."
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
        config = load_cfd_preprocess_config(args.config, project_root=PROJECT_ROOT)
        result = run_cfd_preprocess(config, project_root=PROJECT_ROOT)
        print_result(result)
        return 0 if result.status == "CFD_PREPROCESS_BASELINE_PASS" else 2
    except Exception as error:  # top-level CLI status reporting
        print(f"CFD preprocessing failed: {error}")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
