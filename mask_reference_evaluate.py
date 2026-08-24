"""Run the v5-1 experimental Mask provenance and A/B/C assessment."""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from utils.cfd_lumen.mask_assessment import AssessmentPaths, run_mask_assessment


PROJECT_ROOT = Path(__file__).resolve().parent
DATASET_ROOT = (
    PROJECT_ROOT / "vessel_model" / "T - A high-resolution dataset of mouse brain vasculature"
)
DEFAULT_SAMPLING_RUN = PROJECT_ROOT / "outputs/sampling/20260822_171206_radius_plus_structure_k5"
DEFAULT_FORMAL_V5_RUN = (
    PROJECT_ROOT / "outputs/model_generate/validation_v5_final_roi003274_20260823"
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Experimental-only v5-1 verification of dataset Mask provenance, SWC dependence, "
            "and local junction utility. The formal SWC-only pipeline is not modified."
        )
    )
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--dataset-root", type=Path, default=DATASET_ROOT)
    parser.add_argument("--sampling-run", type=Path, default=DEFAULT_SAMPLING_RUN)
    parser.add_argument("--formal-v5-run", type=Path, default=DEFAULT_FORMAL_V5_RUN)
    parser.add_argument("--inventory-workers", type=int, default=8)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.inventory_workers < 1:
        raise SystemExit("--inventory-workers must be positive")
    if any(character not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-" for character in args.run_id):
        raise SystemExit("--run-id may contain only letters, digits, dot, underscore, and hyphen")
    os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
    paths = AssessmentPaths(
        project_root=PROJECT_ROOT.resolve(),
        dataset_root=args.dataset_root.resolve(),
        sampling_run=args.sampling_run.resolve(),
        formal_v5_run=args.formal_v5_run.resolve(),
        run_root=(PROJECT_ROOT / "outputs/model_generate" / args.run_id).resolve(),
    )
    print(f"Python interpreter: {sys.executable}")
    print("experimental_only: true")
    print(f"run_root: {paths.run_root}")
    summary = run_mask_assessment(paths, inventory_workers=args.inventory_workers)
    print(f"Mask provenance: {summary['current_roi_mask_provenance']['mask_provenance']}")
    print(f"Global Dice: {summary['global_mask_comparison']['dice']:.6f}")
    print(f"Recommendation: {summary['recommendation']['choice']}. {summary['recommendation']['decision']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
