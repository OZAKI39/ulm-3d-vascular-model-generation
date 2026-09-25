#!/usr/bin/env python3
"""Prepare a derived mouse dataset using the native VascularMD model in physical micrometres."""
from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vascular_processing.mouse_dataset import prepare_mouse_dataset


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=ROOT / "configs/swc_roi_generate_mouse_raw.yaml")
    parser.add_argument("--output-dir", type=Path, required=True, help="New derived dataset directory, outside the original dataset")
    parser.add_argument("--accept-native-merges", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--auto-resample", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--model-mode", choices=("network", "branches"), default="network",
                        help="branches: official per-branch AIC with fixed original endpoints; native merges disabled")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    report = prepare_mouse_dataset(args.config, args.output_dir, project_root=ROOT,
                                   accept_native_merges=args.accept_native_merges, auto_resample=args.auto_resample,
                                   model_mode=args.model_mode)
    print(f"Mouse VascularMD preparation: {report['status']}; {args.output_dir.resolve() / 'vascularmd_mouse_manifest.json'}")
    return 0 if report["status"] == "success" else 1


if __name__ == "__main__":
    raise SystemExit(main())
