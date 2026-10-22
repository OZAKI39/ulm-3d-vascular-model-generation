"""Verify the curated public package without a solver, camera or training job."""
import argparse
import ast
import gzip
import hashlib
import json
from pathlib import Path
import re
import subprocess

import numpy as np

OUT = Path(__file__).resolve().parent
REPO = OUT.parents[3]
PREFIX = "computation_examples/mevo_acoustic_microbubbles/"
PARENT = "c7363522ef492d09f7902c0c3382572a70beb68a"


def js(path):
    return json.loads(path.read_text())


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def git(*args):
    return subprocess.check_output(["git", *args], cwd=REPO).decode()


def write(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def restore():
    restored = []
    for name, entry in js(OUT / "COMPRESSED_FILES.json").items():
        raw = gzip.decompress((REPO / entry["destination"]).read_bytes())
        assert sha(raw) == entry["sha256"]
        target = REPO / name
        if target.exists():
            assert sha(target.read_bytes()) == entry["sha256"], "Refusing to replace changed file: " + name
        else:
            target.write_bytes(raw)
            restored.append(name)
    print(json.dumps({"restored_files": len(restored), "overwritten_files": 0}))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seal", action="store_true")
    ap.add_argument("--check-sources", action="store_true")
    ap.add_argument("--restore-compressed", action="store_true")
    args = ap.parse_args()
    if args.restore_compressed:
        restore()
        return
    checks = {}
    copied = js(OUT / "COPIED_FILES.json")
    for name, entry in copied.items():
        raw = (REPO / name).read_bytes()
        assert len(raw) == entry["bytes"] and sha(raw) == entry["sha256"], name
        if "original_sha256" in entry:
            assert sha(gzip.decompress(raw)) == entry["original_sha256"], name
        if args.check_sources:
            assert sha(Path(entry["source_path"]).read_bytes()) == entry.get("original_sha256", entry["sha256"]), name
    checks["copied_files_verified"] = len(copied)
    checks["source_bytes_checked"] = args.check_sources
    checks["compressed_originals_verified"] = len(js(OUT / "COMPRESSED_FILES.json"))
    excluded = js(OUT / "EXCLUDED_FILES.json")
    for name in excluded:
        assert not (REPO / name).exists(), "Excluded artifact unexpectedly present: " + name
    checks["excluded_files_absent"] = len(excluded)
    for name, entry in js(REPO / PREFIX / "sync/20261015/SCIENTIFIC_CODE_MANIFEST.json").items():
        assert sha((REPO / PREFIX / name).read_bytes()) == entry["sha256"], name
    checks["inherited_acoustic_code_files_verified"] = 988
    remote = js(OUT / "REMOTE_READONLY_AUDIT.json")["comparison"]
    assert not remote["missing"] and not remote["mismatched"]
    checks["server_original_files_match_local_mirror"] = remote["expected_local_mirror_files"]
    vision = REPO / "external_connect/vision_tracking_research"
    current = vision / "smart_online_tracking_v2"
    freezes = {}
    for p in (current / "results/expert_freezes").glob("*/frozen.json"):
        rows = js(p)["rows"]
        for row in rows:
            assert sha((p.parent / row["frozen_image_file"]).read_bytes()) == row["record"]["image_sha256"]
            assert sha((p.parent / row["frozen_mask_file"]).read_bytes()) == row["annotation"]["mask_sha256"]
        freezes[p.parent.name] = len(rows)
    checks["frozen_annotation_image_and_mask_rows"] = freezes
    annotations = {}
    for p in (current / "data/annotations").glob("*/latest.json"):
        d = js(p)
        assert sha((p.parent / d["mask_file"]).read_bytes()) == d["mask_sha256"]
        annotations[p.parent.name] = dict(status=d.get("status"), review_scope=d.get("review_scope"))
    checks["current_saved_annotations"] = annotations
    gate = js(vision / "collective_motion_v1/FINAL_GATE.json")
    assert gate["training_ready"] is False
    checks["F10U_scientific_gate_preserved"] = True
    duplicate = excluded["external_connect/vision_tracking_research/collective_motion_v1/results/real_observation.json"]
    assert sha(gzip.decompress((vision / "collective_motion_v1/results/real_observation.json.gz").read_bytes())) == duplicate["sha256"]

    seal_omissions = {str((OUT / name).relative_to(REPO)) for name in ["PACKAGE_MANIFEST.json", "VERIFICATION.json"]}
    if args.seal:
        before = {line.split("\t", 1)[1]: line.split("\t", 1)[0].split()[2] for line in git("ls-tree", "-r", "-z", PARENT).split("\0") if line}
        after = {line.split("\t", 1)[1]: line.split("\t", 1)[0].split()[1] for line in git("ls-files", "--stage", "-z").split("\0") if line}
        assert set(before) <= set(after), "Historical files deleted"
        assert all(before[name] == after[name] for name in before if name != "MEVO_SNAPSHOT.md"), "Unexpected historical-file edit"
        changed = set(after) - set(before) | {"MEVO_SNAPSHOT.md"}
        assert set(copied) <= changed, "Exported payload missing from index"
        assert all(n == "MEVO_SNAPSHOT.md" or n.startswith((PREFIX, "external_connect/")) for n in changed)
        for name in changed:
            assert git("hash-object", "--", name).strip() == after[name], "Unstaged change: " + name
        files = {n: {"bytes": (REPO / n).stat().st_size, "sha256": sha((REPO / n).read_bytes())} for n in sorted(changed - seal_omissions)}
    else:
        files = js(OUT / "PACKAGE_MANIFEST.json")["files"]
    assert set(copied) <= set(files), "Incomplete package manifest"
    patterns = {
        "private_key": rb"-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----",
        "github_token": rb"\bgh[pousr]_[A-Za-z0-9]{30,}\b|\bgithub_pat_[A-Za-z0-9_]{50,}\b",
        "hf_token": rb"\bhf_[a-zA-Z0-9]{30,}\b",
        "aws_access": rb"\bAKIA[0-9A-Z]{16}\b",
        "openai_key": rb"\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{40,}\b",
        "credential_url": rb"(?:https?|rtsp)://[^\s/:@]{2,}:[^\s/@]{6,}@",
    }
    counts = dict(json=0, jsonl=0, python=0, npz=0, arrays=0, credential_scanned=0)
    invalid_json = []
    for name, entry in files.items():
        path = REPO / name
        raw = path.read_bytes()
        assert len(raw) == entry["bytes"] and sha(raw) == entry["sha256"], name
        assert len(raw) < 10_000_000, name
        assert not any(x in {".venv", ".git", "__pycache__", "third_party", "models", ".ssh", ".cache"} or x.startswith(".venv") for x in path.relative_to(REPO).parts), name
        assert path.suffix not in {".pt", ".pth", ".safetensors", ".whl", ".pyc", ".dll", ".exe"}
        suffix = path.suffix
        if suffix == ".gz":
            raw = gzip.decompress(raw)
            suffix = Path(name[:-3]).suffix
        for label, pattern in patterns.items():
            assert not re.search(pattern, raw), "Credential pattern " + label + ": " + name
        counts["credential_scanned"] += 1
        if suffix == ".json":
            try:
                json.loads(raw)
                counts["json"] += 1
            except (ValueError, UnicodeError) as e:
                invalid_json.append({"path": name, "error": str(e)})
        if suffix == ".jsonl":
            for line in raw.splitlines():
                if line.strip():
                    json.loads(line)
            counts["jsonl"] += 1
        if suffix == ".py":
            ast.parse(raw, filename=name)
            counts["python"] += 1
        if suffix == ".npz":
            with np.load(path, allow_pickle=False) as arrays:
                for key in arrays.files:
                    _ = arrays[key]
                    counts["arrays"] += 1
            counts["npz"] += 1
    assert not invalid_json, invalid_json
    checks.update(counts)
    checks.update(new_files=len(files), new_bytes=sum(e["bytes"] for e in files.values()), max_new_bytes=max(e["bytes"] for e in files.values()), unrelated_files_unchanged=True, no_deleted_history=True)
    if args.seal:
        write(OUT / "PACKAGE_MANIFEST.json", {"parent": PARENT, "self_exclusions": sorted(seal_omissions), "files": files})
        write(OUT / "VERIFICATION.json", checks)
    print(json.dumps(checks, indent=2))


if __name__ == "__main__":
    main()
