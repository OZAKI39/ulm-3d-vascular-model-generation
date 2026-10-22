"""Create the formal locally extended and tagged CFD surface."""

from __future__ import annotations

import argparse
from pathlib import Path

from utils.cfd_surface_prepare.config import load_surface_prepare_config
from utils.cfd_surface_prepare.io import SurfacePrepareError
from utils.cfd_surface_prepare.vmtk_pipeline import (
    SUCCESS_STATUS,
    print_vmtk_result,
    run_vmtk_surface_prepare,
)


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "cfd_surface_prepare.yaml"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create one official VMTK TPS boundary-normal extension, remesh the "
            "cross-seam active entity while preserving FAR_CORE, and cap the "
            "validated OPEN surface."
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
        result = run_vmtk_surface_prepare(config, project_root=PROJECT_ROOT)
        print_vmtk_result(result)
        return 0 if result.status == SUCCESS_STATUS else 2
    except SurfacePrepareError as error:
        status = str(error).split(":", maxsplit=1)[0]
        print(f"CFD surface preparation failed: {error}")
        print(f"STATUS: {status}")
        print("NEXT: REVIEW CFD SURFACE PREPARATION FAILURE")
        return 1
    except (FileNotFoundError, ValueError) as error:
        print(f"CFD surface preparation failed: {error}")
        print("STATUS: CFD_SURFACE_INPUT_INVALID")
        print("NEXT: REVIEW CFD SURFACE PREPARATION FAILURE")
        return 1
    except Exception as error:  # internal implementation failure remains distinct
        print(f"CFD surface preparation failed: {error}")
        print("STATUS: CFD_SURFACE_PREPARE_INTERNAL_ERROR")
        print("NEXT: REVIEW CFD SURFACE PREPARATION INTERNAL ERROR")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
