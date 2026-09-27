#!/usr/bin/env python3
"""Model BraVa SWC with VascularMD's native AIC penalized splines."""
from __future__ import annotations

import argparse
import json
import logging
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vascular_processing.pipeline import Options, process_file


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--recursive", action="store_true")
    parser.add_argument("--sample-mode", choices=("original-count", "spacing"), default="original-count")
    parser.add_argument("--model-mode", choices=("network", "branches"), default="network", help="branches: fixed original endpoints, native per-branch AIC, no Nfurcation; requires --no-surface")
    parser.add_argument("--spacing-mm", type=float, help="Arc spacing for BraVa mm data; input numbers are never scaled")
    resample = parser.add_mutually_exclusive_group()
    resample.add_argument("--auto-resample", dest="auto_resample", action="store_true")
    resample.add_argument("--no-auto-resample", dest="auto_resample", action="store_false")
    parser.set_defaults(auto_resample=None)
    parser.add_argument("--surface", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--circumferential-n", type=int, default=24)
    parser.add_argument("--longitudinal-density", type=float, default=0.2, help="Native d: section spacing relative to radius, NOT sections per mm")
    parser.add_argument("--qc-plot", action="store_true")
    parser.add_argument("--plot-branches", type=int, default=6)
    parser.add_argument("--blas-threads", type=int, default=1, help="Avoid oversubscribing small native least-squares systems")
    parser.add_argument("--accept-native-merges", action="store_true", help="Accept native contraction of nearby bifurcations; still require every original terminal and a valid outward tree")
    parser.add_argument("--input-units", choices=("mm", "um"), default="mm", help="Physical unit already shared by XYZ and radius; this flag records units, it does not rescale input")
    args = parser.parse_args(argv)
    options = Options(**{name: getattr(args, name) for name in Options.__dataclass_fields__})
    try:
        options.validate()
    except ValueError as exc:
        parser.error(str(exc))
    input_path, output = args.input.resolve(), args.output_dir.resolve()
    if input_path.is_file():
        inputs = [input_path]
    elif input_path.is_dir():
        iterator = input_path.rglob("*") if args.recursive else input_path.iterdir()
        inputs = sorted(p for p in iterator if p.is_file() and p.suffix.lower() == ".swc"
                        and not p.stem.endswith("_vmd_smooth") and not p.is_relative_to(output))
    else:
        parser.error(f"Input does not exist: {input_path}")
    if not inputs:
        parser.error("No source SWC files found")
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
    outcomes = []
    for index, path in enumerate(inputs, 1):
        target = output / path.parent.relative_to(input_path) if input_path.is_dir() else output
        logging.info("[%s/%s] %s", index, len(inputs), path)
        try:
            report = process_file(path, target, options)
            outcomes.append({"input": str(path), "status": report["status"], "warnings": report["warnings"],
                             "failed_components": report["failed_components"], "outputs": report["outputs"]})
        except Exception as exc:
            logging.exception("Input failed before processing; continuing batch: %s", path)
            outcomes.append({"input": str(path), "status": "failed", "error": str(exc)})
    output.mkdir(parents=True, exist_ok=True)
    summary = output / f"vascularmd_batch_{datetime.now():%Y%m%d_%H%M%S_%f}.json"
    summary.write_text(json.dumps(outcomes, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Batch report: {summary}")
    return 0 if all(item["status"] == "success" for item in outcomes) else 1


if __name__ == "__main__":
    raise SystemExit(main())
