"""Curate the completed F10K worktree onto the published F5 sync, without its bulk history."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

PARENT = "c8ecf9deece1bb0e86b0d0f1c574b984a30ebbd0"
SOURCE = "850fd2df3225351f4fe9564f39cbe72a38662b30"
PREFIX = "computation_examples/mevo_acoustic_microbubbles/"
BRANCH = "sync/mevo-acoustic-f10k-20261013"
DEFAULT_SOURCE = Path("/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10k_20261012")
ESSENTIAL_H5 = {
    *(
        f"results/fast_pzt/f10a/extended_{dt}_v1/{mode}/trajectory.h5"
        for dt in ["10", "5"]
        for mode in ["OFF", "S_1000", "O_1000"]
    ),
    *(f"results/fast_pzt/f10d/campaign_v2_diversity/prepare_{seed}/trajectory.h5" for seed in [1101, 1102]),
    "results/fast_pzt/f10h/b2_dt5_v1/O_1000_1102_B2_GEOM_K2_DT5/trajectory.h5",
    "results/fast_pzt/f10h/b2_short_v1/O_1000_1102_B2_GEOM_K2_DT10/trajectory.h5",
    "results/fast_pzt/f10k/topology_dt10_v1/trial/trajectory.h5",
    "results/fast_pzt/f10k/topology_dt5_v1/trial/trajectory.h5",
}


def git(root, *args):
    return subprocess.check_output(["git", *args], cwd=root)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def tree(root, ref):
    result = {}
    for record in git(root, "ls-tree", "-rl", ref, PREFIX).decode().splitlines():
        meta, path = record.split("\t", 1)
        mode, _, oid, size = meta.split()
        result[path[len(PREFIX) :]] = dict(mode=mode, git_blob=oid, bytes=int(size))
    return result


def exclusion(name, entry):
    p = Path(name)
    if not name.startswith("results/"):
        return None
    if any(part in p.parts for part in ["pdf_review", "__pycache__", "frames", "previews"]):
        return "rebuildable preview/cache"
    if "/radial/" in name:
        return "raw radial waveforms; retain response summaries and scientific identities"
    if p.suffix in {".mp4", ".mov", ".avi", ".tiff", ".tif", ".svg", ".zip", ".tar"}:
        return "video, duplicate visual export or bulk archive"
    if p.name in {"SOFT_STEPS.json.gz", "ALL_STEPS.json.gz", "CONTACT_TRACE.json.gz"}:
        return "full per-step solver observer; retained in complete local archive"
    if p.suffix == ".h5" and name not in ESSENTIAL_H5 and entry["bytes"] > 256_000:
        return "bulk/repeated historical trajectory; keep summaries and selected current scientific sources"
    if p.suffix == ".gz" and entry["bytes"] > 256_000:
        return "bulk compressed event/observer stream"
    if p.suffix == ".npz" and entry["bytes"] > 2_000_000:
        return "large derived/repeated trajectory arrays"
    if entry["bytes"] > 10_000_000:
        return "large detailed export; retain gate/report and full-source recovery hash"
    return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--destination-root", type=Path, required=True)
    args = parser.parse_args()
    source, target = args.source_root.resolve(), args.destination_root.resolve()
    assert source != target
    assert git(source, "rev-parse", "HEAD").decode().strip() == SOURCE
    assert not git(source, "status", "--porcelain").strip(), "Source worktree must remain clean."
    assert git(target, "branch", "--show-current").decode().strip() == BRANCH
    original, current = tree(source, PARENT), tree(source, SOURCE)
    assert not (set(original) - set(current)), "Review historical tracked deletions before exporting."
    delta = {
        n: d for n, d in current.items() if n not in original or original[n]["git_blob"] != d["git_blob"]
    }
    dest_project = target / PREFIX
    out = dest_project / "sync/20261013"
    included, excluded = {}, {}
    for name, entry in sorted(delta.items()):
        path = source / PREFIX / name
        data = os.readlink(path).encode() if entry["mode"] == "120000" else path.read_bytes()
        assert len(data) == entry["bytes"]
        assert (
            hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest() == entry["git_blob"]
        )
        evidence = dict(**entry, sha256=digest(data), source_path=str(path), source_commit=SOURCE)
        reason = exclusion(name, entry)
        if reason:
            excluded[name] = dict(**evidence, reason=reason)
            continue
        dest = dest_project / name
        dest.parent.mkdir(parents=True, exist_ok=True)
        if entry["mode"] == "120000":
            if dest.is_symlink():
                dest.unlink()
            dest.symlink_to(os.readlink(path))
        else:
            shutil.copy2(path, dest)
        included[name] = evidence
    write(out / "COPIED_FILES.json", included)
    write(out / "EXCLUDED_FILES.json", excluded)
    write(
        out / "EXPORT_POLICY.json",
        dict(
            date="2026-10-13",
            branch=BRANCH,
            parent_published_sync=PARENT,
            complete_scientific_delivery=SOURCE,
            source_root=str(source),
            operation="Read-only scientific packaging, no numerical work, no force push",
            include=[
                "All current source, tests, configurations and reports",
                "Compact F6-F10K results, plots, execution identities and success/failure logs",
                "Selected original current scientific trajectories",
                "Inherited published inputs and earlier evidence unchanged",
                "F10K immutable source-only execution snapshots",
                "Current research progress document",
            ],
            exclude=[
                "New videos and preview frames",
                "Environments, caches, downloads and credentials",
                "Repeated batch trajectories and full per-step observer output",
                "Large derived arrays and detailed exports beyond registered sync size rules",
                "Unrelated project changes",
            ],
            rule_source="scripts/export_sync_20261013.py",
            essential_h5=sorted(ESSENTIAL_H5),
            inherited_files=len(original),
            inherited_bytes=sum(d["bytes"] for d in original.values()),
            copied_count=len(included),
            copied_bytes=sum(d["bytes"] for d in included.values()),
            excluded_count=len(excluded),
            excluded_bytes=sum(d["bytes"] for d in excluded.values()),
            limits_note="10MB maximum new copied artifact; selected HDF5 are below2.1MB. Already-public33MB acoustic basis is inherited, not duplicated. Existing published historical data is not purged.",
            history_note="A curated commit is based on the published sync. The intervening approximately2GB local result history is not made an ancestor or uploaded.",
            scientific_status="EXPLORATORY / PHYSICS_UNCALIBRATED / training_ready=false",
            old_manifests="Original manifests remain byte-identical and may reference deliberately excluded files; this sync inventory defines actual presence.",
            scientific_computation_started=False,
        ),
    )
    # Recover immutable F10K execution source even though bulky local Git history is not pushed.
    remote_sources = json.loads(
        (source / PREFIX / "results/fast_pzt/f10k/REMOTE_SOURCE_MANIFEST.json").read_text()
    )["sources"]
    snapshots = {}
    for commit, files in sorted(remote_sources.items()):
        archive = out / "execution_sources" / (commit + ".tar.gz")
        archive.parent.mkdir(exist_ok=True)
        archive.write_bytes(git(source, "archive", "--format=tar.gz", commit, *[PREFIX + n for n in files]))
        snapshots[commit] = dict(
            path=str(archive.relative_to(dest_project)),
            bytes=archive.stat().st_size,
            sha256=digest(archive.read_bytes()),
            source_file_count=len(files),
            purpose="Original source only; no bulk numerical outputs; exact file hashes in REMOTE_SOURCE_MANIFEST.json",
        )
    write(out / "EXECUTION_SNAPSHOTS.json", snapshots)
    progress = Path("/home/lzy/projects/microbubble_aggregation.md")
    expected = json.loads(
        (source / PREFIX / "results/fast_pzt/f10k/progress_document/CHECK.json").read_text()
    )["file_sha256"][str(progress)]
    assert digest(progress.read_bytes()) == expected
    exported_progress = dest_project / "docs/MICROBUBBLE_AGGREGATION_PROGRESS_20261013.md"
    shutil.copy2(progress, exported_progress)
    write(
        out / "PROGRESS_PROVENANCE.json",
        dict(
            source=str(progress),
            destination=str(exported_progress.relative_to(dest_project)),
            sha256=expected,
            original_bytes_unchanged=True,
        ),
    )
    print(
        json.dumps(
            dict(
                copied_files=len(included),
                copied_MB=sum(d["bytes"] for d in included.values()) / 1e6,
                excluded_files=len(excluded),
                excluded_GB=sum(d["bytes"] for d in excluded.values()) / 1e9,
                immutable_source_snapshots=len(snapshots),
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
