"""Verify this curated export without starting a solver or changing scientific files."""

import argparse
import ast
import hashlib
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def checked_path(relative):
    path = Path(relative)
    assert not path.is_absolute() and ".." not in path.parts, relative
    result = ROOT / path
    assert result.is_file() and not result.is_symlink(), relative
    return result


def verify(scientific=False):
    manifest = json.loads((ROOT / "sync/FILE_MANIFEST.json").read_text())
    entries = manifest["files"]
    assert len({e["path"] for e in entries}) == len(entries)
    for entry in entries:
        path = checked_path(entry["path"])
        assert path.stat().st_size == entry["bytes"], entry["path"]
        assert sha(path) == entry["sha256"], entry["path"]
        assert entry["bytes"] < 100_000_000, entry["path"]
    original = json.loads((ROOT / "sync/SOURCE_TRACKED_FILES.json").read_text())
    for name, digest in original["files"].items():
        assert sha(checked_path(name)) == digest, name
    parsed = 0
    for folder in ["src", "scripts", "tests", "benchmarks", "references", "validation", "sync"]:
        for path in (ROOT / folder).rglob("*.py"):
            ast.parse(path.read_text(), filename=str(path.relative_to(ROOT)))
            parsed += 1
    policy = json.loads((ROOT / "sync/EXPORT_POLICY.json").read_text())
    gold_dir = ROOT / "results/stage4/performance_refactor/gold_baseline"
    gold_path = gold_dir / "GOLD_BASELINE_MANIFEST.json"
    assert sha(gold_path) == policy["original_gold_manifest_sha256"]
    gold = json.loads(gold_path.read_text())
    omitted = []
    for name, entry in gold["files"].items():
        path = gold_dir / name
        if not path.is_file():
            omitted.append(name)
        else:
            assert sha(path) == entry["sha256"], name
    assert sorted(omitted) == ["cases/B2/particles/events.jsonl", "cases/B3/particles/events.jsonl"], omitted
    delivery = ROOT / "results/stage4/visualization/B2_B3_clean_enlarged_v7"
    hashes = json.loads((delivery / "DELIVERY_SHA256.json").read_text())
    for name, entry in hashes.items():
        path = delivery / name
        assert path.stat().st_size == entry["bytes"] and sha(path) == entry["sha256"], name
    for relative in ["OPEN_RESULTS.html", "shape_audit/SHAPE_AUDIT.html"]:
        page = delivery / relative
        for link in re.findall(r'(?:href|src|poster)="([^"]+)"', page.read_text()):
            assert (page.parent / link).is_file(), link
    external = json.loads((ROOT / "sync/EXTERNAL_REFERENCES.json").read_text())
    for name in ["RENDER_MANIFEST.json", "shape_audit/SHAPE_AUDIT.json"]:
        data = json.loads((delivery / name).read_text())
        for source, digest in (data.get("protected_inputs", {}) | data.get("renderer_sources", {})).items():
            prefix = policy["local_source"] + "/"
            if not source.startswith(prefix):
                assert source in external, source
                assert external[source]["sha256"] == digest
                assert external[source]["role"] in {"visual_style_reference_only", "rendering_font"}
                continue
            assert sha(checked_path(source.removeprefix(prefix))) == digest, source
    for name in ["particles/trajectories.h5", "particles/events.jsonl"]:
        assert not (ROOT / policy["active_J_binary_data_excluded"] / name).exists()
    report = dict(
        status="PASS",
        files=len(entries),
        bytes=sum(e["bytes"] for e in entries),
        original_tracked_files_identical=len(original["files"]),
        python_files_parsed=parsed,
        latest_delivery_files_checked=len(hashes),
        original_gold_manifest_unchanged=True,
        gold_export="PARTIAL; exactly two large event streams omitted",
        external_visual_references_declared=len(external),
        active_J_binary_data_excluded=True,
        solvers_started=False,
    )
    if scientific:
        import h5py
        import numpy as np

        cases = {}
        for name, case in gold["cases"].items():
            trajectory = ROOT / case["result_id"] / "particles/trajectories.h5"
            assert sha(trajectory) == gold["files"][f"cases/{name}/particles/trajectories.h5"]["sha256"]
            with h5py.File(trajectory, "r") as f:
                assert f.attrs["complete"]
                assert f["positions_m"].shape == (301, 256, 3)
                assert np.allclose(f["time_s"][:], np.arange(301) * 0.02, rtol=0, atol=1e-12)
                assert len(np.unique(f["particle_id"][:])) == 256
                assert np.allclose(f["radii_m"][:], 1.78e-6, rtol=0, atol=1e-18)
                assert np.isfinite(f["positions_m"][:]).all()
            cases[name] = dict(
                complete=True, frames=301, bubbles=256, duration_s=6, branch_counts=case["branch_counts"]
            )
        report["scientific_data_checks"] = cases
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scientific", action="store_true")
    parser.add_argument(
        "--report", type=Path, help="Optional new verification report; no scientific files modified"
    )
    args = parser.parse_args()
    report = verify(args.scientific)
    text = json.dumps(report, indent=2) + "\n"
    if args.report:
        args.report.write_text(text)
    print(text, end="")
