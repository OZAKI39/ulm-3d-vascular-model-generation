#!/usr/bin/env python3
"""BraVa landmark-assisted MeVO extraction, with exact geometry and distal coverage."""
import argparse
import json
import logging
from pathlib import Path
import sys

import yaml

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from vascular_processing.mevo_annotation import annotate
from vascular_processing.mevo_export import write_json
from vascular_processing.mevo_pipeline import extract_file, inspect_file
from vascular_processing.roi_landmarks import ROI_DEFINITIONS, sha256


def modeling_args(parser):
    parser.add_argument("--allow-natural-terminal", action="store_true")
    parser.add_argument("--proximal-context", choices=("none", "upstream-edge"), default="none")
    parser.add_argument("--model", action="store_true")
    parser.add_argument("--surface", action="store_true", help="Implies modeling; use the official VascularMD mesh pipeline")
    parser.add_argument("--accept-native-merges", action="store_true")


def batch(directory, landmarks_dir, output, *, recursive=False, **options):
    output = output.resolve()
    if not directory.is_dir() or not landmarks_dir.is_dir():
        raise ValueError("Batch input and landmarks paths must be directories")
    output.mkdir(parents=True, exist_ok=False)
    candidates = {}
    annotation_errors = []
    for path in sorted(landmarks_dir.rglob("*.yaml")):
        try:
            doc = yaml.safe_load(path.read_text())
            if isinstance(doc, dict) and isinstance(doc.get("source"), dict):
                rois = doc.get("rois", {})
                if isinstance(rois, dict) and all(isinstance(roi, dict) and roi.get("enabled") is False for roi in rois.values()):
                    continue  # Blank inspection templates are not competing annotations.
                candidates.setdefault(doc["source"].get("sha256", "").lower(), []).append(path)
        except Exception as exc:
            annotation_errors.append({"file": str(path), "error": str(exc)})
    inputs = sorted(directory.rglob("*.swc") if recursive else directory.glob("*.swc"))
    inputs = [p for p in inputs if not p.resolve().is_relative_to(output)]
    outcomes = []
    for source in inputs:
        matches = candidates.get(sha256(source), [])
        if len(matches) != 1:
            outcomes.append({"source": str(source), "status": "NEEDS_MANUAL_REVIEW",
                "error": "Missing or ambiguous source-hash-matched landmark file", "matches": list(map(str, matches))})
            continue
        target = output / source.relative_to(directory).with_suffix("")
        result = extract_file(source, matches[0], target, **options)
        outcomes.append({"source": str(source), "status": result["status"], "manifest": str(target / result["manifest"])})
    status = "PASS" if outcomes and not annotation_errors and all(x["status"] == "PASS" for x in outcomes) else "PASS_WITH_WARNINGS"
    if not outcomes or annotation_errors or any(x["status"] not in {"PASS", "PASS_WITH_WARNINGS"} for x in outcomes):
        status = "NEEDS_MANUAL_REVIEW"
    report = {"status": status, "matching_rule": "exact source SHA256; filename checked during extraction; never file order",
              "files": outcomes, "annotation_errors": annotation_errors}
    write_json(output / "mevo_batch_manifest.json", report)
    return report


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    inspect = commands.add_parser("inspect")
    inspect.add_argument("input", type=Path)
    inspect.add_argument("--output-dir", type=Path, required=True)
    inspect.add_argument("--input-units", choices=("mm", "um"), default="mm")
    inspect.add_argument("--subject-id")
    annotation = commands.add_parser("annotate")
    annotation.add_argument("input", type=Path)
    annotation.add_argument("--output", type=Path, required=True)
    annotation.add_argument("--roi", choices=tuple(ROI_DEFINITIONS), default="RMCA_M2M3")
    annotation.add_argument("--input-units", choices=("mm", "um"), default="mm")
    extraction = commands.add_parser("extract")
    extraction.add_argument("input", type=Path)
    extraction.add_argument("--landmarks", type=Path, required=True)
    extraction.add_argument("--output-dir", type=Path, required=True)
    extraction.add_argument("--allow-source-mismatch", action="store_true")
    extraction.add_argument("--territory", choices=("MCA", "ACA", "PCA"))
    extraction.add_argument("--side", choices=("L", "R"))
    modeling_args(extraction)
    batches = commands.add_parser("batch")
    batches.add_argument("input", type=Path)
    batches.add_argument("--landmarks-dir", type=Path, required=True)
    batches.add_argument("--output-dir", type=Path, required=True)
    batches.add_argument("--recursive", action="store_true")
    modeling_args(batches)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    try:
        if args.command == "inspect":
            report = inspect_file(args.input, args.output_dir, units=args.input_units, subject_id=args.subject_id)
        elif args.command == "annotate":
            report = annotate(args.input, args.output, roi_name=args.roi, units=args.input_units)
        else:
            options = {key: getattr(args, key) for key in ("allow_natural_terminal", "proximal_context", "model", "surface", "accept_native_merges")}
            if args.command == "extract":
                report = extract_file(args.input, args.landmarks, args.output_dir, **options,
                    allow_source_mismatch=args.allow_source_mismatch, territory=args.territory, side=args.side)
            else:
                report = batch(args.input, args.landmarks_dir, args.output_dir, recursive=args.recursive, **options)
        print(json.dumps({k: report[k] for k in ("status", "needs_manual_review", "error", "message", "source", "topology_node_count", "outputs") if k in report},
                         indent=2, ensure_ascii=False))
        return 0 if report["status"] in {"PASS", "PASS_WITH_WARNINGS"} else 2 if report["status"] == "NEEDS_MANUAL_REVIEW" else 1
    except Exception as exc:
        logging.error("FAILED: %s: %s", type(exc).__name__, exc)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
