"""Read-only, size-aware export of local F10T/F10U and device interface evidence.

The published F10S tree is the parent. Development histories/environments are
never grafted into the public branch. Run only in a fresh isolated worktree.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

PARENT = "c7363522ef492d09f7902c0c3382572a70beb68a"
ACOUSTIC = "31a44d3f342e5d0bca55f7f9b83f4d32b8189696"
VISION = "0561769933dfa799dee0778f638a8c26270c3b95"
BRANCH = "sync/mevo-acoustic-vision-f10u-20261017"
PREFIX = "computation_examples/mevo_acoustic_microbubbles/"
EXTERNAL = Path("/home/lzy/projects/external_connect")
SCIENTIFIC = Path("/home/lzy/projects/computation_examples/_archives/mevo_fast_pzt_f10s_20261015")
OUT = Path(__file__).resolve().parent
REPO = OUT.parents[3]
PRUNE = {".git", "__pycache__", ".cache", ".pytest_cache", ".ruff_cache", "third_party", "models", "node_modules"}
TEXT = {".py", ".md", ".txt", ".sh", ".cmd", ".ps1", ".json", ".jsonl", ".csv", ".yaml", ".yml", ".toml", ".ini", ".log", ".patch", ".diff"}


def git(root, *args):
    return subprocess.check_output(["git", *args], cwd=root).decode().strip()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def inventory(root, omitted):
    for parent, dirs, files in os.walk(root):
        for d in list(dirs):
            if d in PRUNE or d.startswith(".venv") or (root.name == "CHRONOS_control" and d == "tools") or (root.name == "Tektronix_control" and d == "drivers"):
                dirs.remove(d)
                p = Path(parent) / d
                omitted.append({"source_path": str(p), "kind": "directory", "reason": "Environment, downloaded dependency/model/binary, Git database or cache; no individual inventory/hash claimed"})
        for f in sorted(files):
            yield Path(parent) / f


def choice(area, name, tracked, annotated):
    p = Path(name)
    if p.suffix in {".pyc", ".pyo", ".pt", ".pth", ".safetensors", ".exe", ".dll", ".whl"} or p.name in {"service.lock", "gpu.lock"}:
        return "Generated lock, compiled binary or model checkpoint"
    if p.name in {".env", "credentials", "id_rsa", "id_ed25519"}:
        return "Local private configuration"
    if area == "vision_tracking_research" and name in tracked:
        return None
    if area == "vision_tracking_research":
        if name == "collective_motion_v1/results/real_observation.json":
            return "Identical content retained losslessly in real_observation.json.gz"
        if name.startswith("smart_online_tracking_v2/results/expert_freezes/"):
            return None  # includes original image/mask TIFF training exports
        if name.startswith("evidence/official/"):
            return "Downloaded reference-page cache; references and license evidence retained"
        if name.startswith("smart_online_tracking_v2/data/candidates/"):
            return None if p.stem in annotated else "Unannotated candidate frame; original manifest and annotated frames retained"
        if p.suffix == ".mp4":
            return None if name.startswith(("collective_motion_v1/videos/", "smart_online_tracking_v2/videos/")) else "Rebuildable intermediate movie"
        if p.suffix in {".png", ".jpg", ".jpeg", ".pdf", ".npy", ".npz"}:
            keep = name.startswith(("figures/", "smart_online_tracking_v2/figures/", "collective_motion_v1/figures/", "collective_motion_v1/results/", "smart_online_tracking_v2/results/expert_freezes/", "smart_online_tracking_v2/data/annotations/", "smart_online_tracking_v2/results/supervised_20261016_v1/predictions/", "smart_online_tracking_v2/results/heldout_49_50_20261016_v1/"))
            return None if keep else "Repeated/raw frame export, temporary proposal or superseded preview; metrics, code and selected evidence retained"
    elif p.suffix in {".png", ".jpg", ".jpeg", ".mp4", ".npy", ".npz"}:
        return "Device capture images/bulk media; interface source, configuration and logs retained"
    if p.suffix in TEXT or p.name.startswith("requirements") or p.name in {".gitignore", ".gitattributes", "LICENSE"}:
        return None
    return "Nonessential generated artifact or dependency"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plan-only", action="store_true")
    args = ap.parse_args()
    assert git(REPO, "rev-parse", "HEAD") == PARENT
    assert git(REPO, "branch", "--show-current") == BRANCH
    vroot = EXTERNAL / "vision_tracking_research"
    assert git(vroot, "rev-parse", "HEAD") == VISION
    assert not git(vroot, "status", "--porcelain", "--untracked-files=no")
    assert git(SCIENTIFIC, "rev-parse", "HEAD") == ACOUSTIC
    assert not git(SCIENTIFIC, "status", "--porcelain")
    tracked = set(git(vroot, "ls-files", "-z").split("\0"))
    annotated = {p.parent.name for p in (vroot / "smart_online_tracking_v2/data/annotations").glob("*/latest.json")}
    copied, excluded, compressed, directories = {}, {}, {}, []
    source_status = git(vroot, "status", "--porcelain", "--untracked-files=all")
    for area in ["vision_tracking_research", "CHRONOS_control", "Tektronix_control"]:
        source = EXTERNAL / area
        for path in inventory(source, directories):
            assert not path.is_symlink(), str(path)
            name = path.relative_to(source).as_posix()
            destname = f"external_connect/{area}/{name}"
            raw = path.read_bytes()
            record = {"source_path": str(path), "bytes": len(raw), "sha256": digest(raw), "source_tracked": name in tracked if area == "vision_tracking_research" else False}
            reason = choice(area, name, tracked, annotated)
            if reason:
                excluded[destname] = dict(record, reason=reason)
                continue
            if len(raw) > 1_000_000 and path.suffix in {".json", ".jsonl", ".log"} and not record["source_tracked"]:
                compressed[destname] = dict(record, destination=destname + ".gz")
                destname += ".gz"
                raw = gzip.compress(raw, compresslevel=9, mtime=0)
                record = dict(record, original_bytes=record["bytes"], original_sha256=record["sha256"], bytes=len(raw), sha256=digest(raw), encoding="gzip; exact original bytes recoverable")
            assert len(raw) < 10_000_000, "Review new file over 10 MB: " + destname
            copied[destname] = record
            if not args.plan_only:
                dest = REPO / destname
                assert not dest.exists(), "Refusing to overwrite " + str(dest)
                dest.parent.mkdir(parents=True, exist_ok=True)
                dest.write_bytes(raw)
                shutil.copymode(path, dest)
    # Raw videos stay local; include their original identities, not their bytes.
    for path in sorted((EXTERNAL / "sample_videos").glob("*")):
        if path.is_file():
            h = hashlib.sha256()
            with path.open("rb") as f:
                for block in iter(lambda: f.read(2**20), b""):
                    h.update(block)
            excluded["external_connect/sample_videos/" + path.name] = {"source_path": str(path), "bytes": path.stat().st_size, "sha256": h.hexdigest(), "reason": "Large original recording; selected labeled images, actual analysis data and short derived videos retained"}
    # Preserve the acoustic code identity already published in F10S.
    inherited = json.loads((REPO / PREFIX / "sync/20261015/SCIENTIFIC_CODE_MANIFEST.json").read_text())
    for name, row in inherited.items():
        expected = row["sha256"]
        assert digest((SCIENTIFIC / PREFIX / name).read_bytes()) == expected, name
        assert digest((REPO / PREFIX / name).read_bytes()) == expected, name
    assert git(vroot, "status", "--porcelain", "--untracked-files=all") == source_status
    progress = Path("/home/lzy/projects/microbubble_aggregation.md")
    progress_name = PREFIX + "docs/MICROBUBBLE_AGGREGATION_PROGRESS_20261017.md"
    raw = progress.read_bytes()
    copied[progress_name] = {"source_path": str(progress), "bytes": len(raw), "sha256": digest(raw), "source_tracked": False}
    if not args.plan_only:
        (REPO / progress_name).write_bytes(raw)
    summary = dict(branch=BRANCH, parent_published_sync=PARENT, acoustic_source=ACOUSTIC, vision_source=VISION,
        source_worktrees_unchanged=True, acoustic_code_files_verified=len(inherited), scientific_computation_started=False,
        included_files=len(copied), included_bytes=sum(v["bytes"] for v in copied.values()),
        excluded_files=len(excluded), excluded_bytes=sum(v["bytes"] for v in excluded.values()),
        excluded_directory_count=len(directories), directory_sizes_not_counted=True,
        losslessly_compressed_files=len(compressed), annotated_candidate_ids=sorted(annotated),
        maximum_new_file_bytes=max(v["bytes"] for v in copied.values()),
        raw_originals_unchanged=True, current_training_ready=False,
        note="Historical absolute paths, source Git identities and original manifests are preserved. Source commits are provenance, not promised reachable public commits. See README for restore instructions.")
    print(json.dumps(summary, indent=2))
    if not args.plan_only:
        write(OUT / "COPIED_FILES.json", copied)
        write(OUT / "EXCLUDED_FILES.json", excluded)
        write(OUT / "EXCLUDED_DIRECTORIES.json", directories)
        write(OUT / "COMPRESSED_FILES.json", compressed)
        write(OUT / "EXPORT_POLICY.json", summary)
        write(OUT / "SOURCE_STATUS_BEFORE.json", {"vision_head": VISION, "vision_status": source_status, "acoustic_head": ACOUSTIC, "acoustic_status": ""})


if __name__ == "__main__":
    main()
