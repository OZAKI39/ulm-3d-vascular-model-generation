"""Export the F10L–F10Q scientific delta onto the published F10K sync.

Read-only with respect to the scientific worktree; never runs a solver.
"""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

PARENT = "15c42e0e422ff5d1521e5cc8c12ede2688b3de59"
BASELINE = "850fd2df3225351f4fe9564f39cbe72a38662b30"
SOURCE = "41ae8ed2b88f2cbbed834fffccaa0b5de6e91589"
BRANCH = "sync/mevo-acoustic-f10q-20261014"
PREFIX = "computation_examples/mevo_acoustic_microbubbles/"
DEFAULT_SOURCE = Path("/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10q_20261014")


def git(root, *args):
    return subprocess.check_output(["git", *args], cwd=root)


def tree(root, ref):
    result = {}
    for line in git(root, "ls-tree", "-rl", ref, PREFIX).decode().splitlines():
        metadata, name = line.split("\t", 1)
        mode, _, oid, size = metadata.split()
        result[name[len(PREFIX):]] = dict(mode=mode, git_blob=oid, bytes=int(size))
    return result


def raw(path):
    return os.readlink(path).encode() if path.is_symlink() else path.read_bytes()


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def exclusion(name):
    p = Path(name)
    if any(part in {".venv", "__pycache__", ".cache", ".ssh"} for part in p.parts):
        return "Environment, cache or private local configuration"
    if p.suffix in {".mp4", ".avi", ".mov", ".whl", ".pyc"}:
        return "Rendered movie or installed binary; source and scientific evidence retained"
    if p.suffix in {".png", ".pdf", ".jpg"} and any(
        part.startswith("review_history") or part in {"pdf_previews", "pdf_review", "previews", "frames"}
        for part in p.parts
    ):
        return "Superseded/rebuildable visual preview; final PNG/PDF and review records retained"
    if name.startswith("results/straight_tube_exp1/f10l/acoustic/basis_") and p.name == "field.npz":
        return "Per-element intermediate field; combined original acoustic/basis.h5 retained"
    if name == "results/straight_tube_exp1/f10m/history/pre_top_port/materials/material_maps.h5":
        return "Superseded pre-port material grid; current grid, configuration and history reports retained"
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    ap.add_argument("--destination-root", type=Path, required=True)
    args = ap.parse_args()
    source, target = args.source_root.resolve(), args.destination_root.resolve()
    assert source != target
    assert git(source, "rev-parse", "HEAD").decode().strip() == SOURCE
    assert not git(source, "status", "--porcelain").strip(), "Scientific source must be clean"
    assert git(target, "rev-parse", "HEAD").decode().strip() == PARENT
    assert git(target, "branch", "--show-current").decode().strip() == BRANCH
    before, current = tree(source, BASELINE), tree(source, SOURCE)
    assert set(before).issubset(current), "Review any historical deletion before exporting"
    delta = {n: e for n, e in current.items() if before.get(n) != e}
    project = target / PREFIX
    out = project / "sync/20261014"
    copied, excluded = {}, {}
    for name, entry in sorted(delta.items()):
        original = source / PREFIX / name
        data = raw(original)
        assert len(data) == entry["bytes"]
        assert hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest() == entry["git_blob"]
        record = dict(**entry, sha256=sha(data), source_path=str(original), source_commit=SOURCE)
        reason = exclusion(name)
        if reason:
            excluded[name] = dict(**record, reason=reason)
            continue
        assert entry["bytes"] < 30_000_000, "Review any new file over 30 MB: " + name
        dest = project / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        assert not dest.is_symlink(), name
        if entry["mode"] == "120000":
            dest.symlink_to(os.readlink(original))
        else:
            shutil.copy2(original, dest)
        copied[name] = record
    write(out / "COPIED_FILES.json", copied)
    write(out / "EXCLUDED_FILES.json", excluded)
    code = {}
    for name, entry in current.items():
        if name.split("/")[0] in {"src", "scripts", "configs", "tests", "benchmarks", "validation"} or (
            "/" not in name and (name.startswith("requirements") or name == "pyproject.toml")
        ):
            assert sha(raw(project / name)) == sha(raw(source / PREFIX / name)), name
            code[name] = dict(**entry, sha256=sha(raw(project / name)))
    write(out / "SCIENTIFIC_CODE_MANIFEST.json", code)
    progress = Path("/home/lzy/projects/microbubble_aggregation.md")
    expected = json.loads((source / PREFIX / "results/straight_tube_exp1/f10q/progress/UPDATE_AUDIT.json").read_text())["after_sha256"]
    assert sha(progress.read_bytes()) == expected
    relative = "docs/MICROBUBBLE_AGGREGATION_PROGRESS_20261014.md"
    shutil.copy2(progress, project / relative)
    write(out / "PROGRESS_PROVENANCE.json", dict(source=str(progress), destination=relative, sha256=expected, original_bytes_unchanged=True))
    write(out / "EXPORT_POLICY.json", dict(
        date="2026-10-14", branch=BRANCH, parent_published_sync=PARENT,
        previous_scientific_source=BASELINE, complete_scientific_delivery=SOURCE,
        source_root=str(source), rule_source="scripts/export_sync_20261014.py",
        operation="Read-only scientific packaging; no solver execution or altered scientific bytes",
        copied_count=len(copied), copied_bytes=sum(e["bytes"] for e in copied.values()),
        excluded_count=len(excluded), excluded_bytes=sum(e["bytes"] for e in excluded.values()),
        source_code_count=len(code), maximum_new_file_bytes=30_000_000,
        essential_larger_files={
            "results/straight_tube_exp1/f10n/cache/fine.h5": "Original 1.5 um field cache required for F10N–Q convergence replay (24.6 MB)",
            "results/straight_tube_exp1/f10q/b_navigation_v1/runs/NAV_Y_DT05/trajectory.h5": "Complete 0.5 us Y-navigation trajectory and contact evidence (13.9 MB); failed shape/ID convergence must remain auditable",
        },
        include=["All current scientific code, configs, tests and dependency declarations", "F10L–Q reports, logs, successes AND failures, original trajectories and analysis arrays", "Four-PZT geometry, APECSS response, coarse/fine water-field caches", "F10Q all-step parameter/force summary", "Final publication PNG/PDF and immutable source-only execution archives", "Unchanged current Chinese progress document"],
        exclude=["New movies, environments, downloaded packages and credentials", "Superseded visual previews and old pre-port material grid", "Three per-element F10L intermediates already represented in original combined basis", "Previously omitted historical bulk remains omitted"],
        inherited_history="Based on published F10K sync; uncurated development commits are not ancestors and are not pushed",
        old_manifests="Historical full-delivery manifests are unchanged and can name deliberately omitted artifacts. This sync manifest defines actual new payload.",
        recovery="Excluded originals remain at each recorded source_path with original SHA-256. Old source commits are provenance identities, not guaranteed GitHub-reachable commits. Immutable execution-source archives are included.",
        scientific_status="EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / training_ready=false",
        scientific_computation_started=False,
    ))
    print(json.dumps(dict(copied_files=len(copied), copied_bytes=sum(e["bytes"] for e in copied.values()), excluded_files=len(excluded), excluded_bytes=sum(e["bytes"] for e in excluded.values()), exact_current_code_files=len(code)), indent=2))


if __name__ == "__main__":
    main()
