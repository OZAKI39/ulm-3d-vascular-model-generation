"""Create a locally extended, tagged CFD surface from saved PASS inputs."""

from __future__ import annotations

import argparse
from pathlib import Path

from utils.cfd_surface_prepare.config import load_surface_prepare_config
from utils.cfd_surface_prepare.vmtk_pipeline import (
    PASS_STATUSES,
    print_vmtk_result,
    run_vmtk_surface_prepare,
)


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "cfd_surface_prepare.yaml"


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create one official VMTK TPS flow-extension candidate without modifying "
            "the original Ultraliser surface or creating a volume mesh."
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
        return 0 if result.status in PASS_STATUSES else 2
    except Exception as error:
        allowed_failures = {
            "VMTK_ENVIRONMENT_BLOCKED",
            "VMTK_TPS_EXTENSION_FAILED",
            "VMTK_EXTENSION_GEOMETRY_FAILED",
            "VMTK_REMESH_CORE_FIDELITY_FAILED",
            "VMTK_SURFACE_QC_FAILED",
            "INVALID_VMTK_EXTENSION_MODE",
            "BOUNDARY_NORMAL_INPUT_PLANE_MISMATCH",
            "VMTK_BOUNDARY_NORMAL_RAW_GEOMETRY_FAILED",
            "VMTK_BOUNDARY_NORMAL_FINAL_SURFACE_FAILED",
            "ORIGINAL_ULTRALISER_GEOMETRY_MODIFIED",
        }
        failure = str(error).split(":", maxsplit=1)[0]
        if failure not in allowed_failures:
            failure = "VMTK_TPS_EXTENSION_FAILED"
        print(f"CFD surface preparation failed: {error}")
        print(f"Final status: {failure}")
        print("NEXT: REVIEW VMTK BOUNDARY-NORMAL TPS FAILURE")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
