"""Verify the curated Git snapshot without starting scientific computation."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile

import h5py

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT / "sync/20261013"
PREFIX = "computation_examples/mevo_acoustic_microbubbles/"
PARENT = "c8ecf9deece1bb0e86b0d0f1c574b984a30ebbd0"


def js(path):
    return json.loads(path.read_text())


def data(path):
    return os.readlink(path).encode() if path.is_symlink() else path.read_bytes()


def sha(path):
    return hashlib.sha256(data(path)).hexdigest()


def git(*args):
    return subprocess.check_output(["git", *args], cwd=REPO)


def write(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def tree(ref):
    result = {}
    for line in git("ls-tree", "-rl", ref).decode().splitlines():
        metadata, name = line.split("\t", 1)
        mode, _, oid, size = metadata.split()
        result[name] = dict(mode=mode, oid=oid, bytes=int(size))
    return result


def scan_credentials(name, raw):
    # Report only locations and pattern names; never echo a matched value.
    patterns = {
        "private_key": rb"-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----",
        "github_token": rb"\bgh[pousr]_[A-Za-z0-9]{30,}\b|\bgithub_pat_[A-Za-z0-9_]{50,}\b",
        "aws_access_id": rb"\bAKIA[0-9A-Z]{16}\b",
        "openai_key": rb"\bsk-(?:proj-|svcacct-)?[A-Za-z0-9_-]{40,}\b",
        "credential_url": rb"https?://[^\s/:@]{3,}:[^\s/@]{8,}@",
    }
    hits = [
        (name, label, raw[: m.start()].count(b"\n") + 1)
        for label, pattern in patterns.items()
        for m in re.finditer(pattern, raw)
    ]
    assert not hits, "Review credential-pattern locations before publication: " + str(hits)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seal", action="store_true", help="Write package manifest after staging this sync")
    args = ap.parse_args()
    checks = {}
    copied = js(OUT / "COPIED_FILES.json")
    excluded = js(OUT / "EXCLUDED_FILES.json")
    policy = js(OUT / "EXPORT_POLICY.json")
    for name, entry in copied.items():
        assert sha(ROOT / name) == entry["sha256"], name
        assert len(data(ROOT / name)) == entry["bytes"]
    for name in excluded:
        assert not (ROOT / name).exists(), name
    assert len(copied) == policy["copied_count"] and len(excluded) == policy["excluded_count"]
    checks["copied_scientific_bytes_identical"] = len(copied)
    checks["omitted_delta_paths_explicit_and_absent"] = len(excluded)
    p = js(OUT / "PROGRESS_PROVENANCE.json")
    assert sha(ROOT / p["destination"]) == p["sha256"]
    checks["research_progress_copy_unchanged"] = True
    source_manifest = js(ROOT / "results/fast_pzt/f10k/REMOTE_SOURCE_MANIFEST.json")["sources"]
    for commit, entry in js(OUT / "EXECUTION_SNAPSHOTS.json").items():
        path = ROOT / entry["path"]
        assert sha(path) == entry["sha256"]
        with tarfile.open(path, "r:gz") as archive:
            actual = {}
            for member in archive:
                if not member.isfile():
                    assert member.isdir(), "Unexpected archive link/member"
                    continue
                assert member.name.startswith(PREFIX) and ".." not in Path(member.name).parts
                raw = archive.extractfile(member).read()
                actual[member.name[len(PREFIX) :]] = hashlib.sha256(raw).hexdigest()
                scan_credentials(entry["path"] + ":" + member.name, raw)
            assert actual == source_manifest[commit]
            assert len(actual) == entry["source_file_count"]
    checks["immutable_F10K_source_snapshots"] = len(source_manifest)
    datasets = 0
    for name in policy["essential_h5"]:
        assert (ROOT / name).exists(), name
        with h5py.File(ROOT / name) as h:

            def count(_name, obj):
                nonlocal datasets
                if isinstance(obj, h5py.Dataset):
                    _ = obj[()]
                    datasets += 1

            h.visititems(count)
    checks["selected_original_HDF5_fully_read"] = len(policy["essential_h5"])
    checks["HDF5_datasets_read"] = datasets
    figure_sources = set()
    for figure in js(ROOT / "results/fast_pzt/f10k/figures/FIGURE_PROVENANCE.json")["figures"]:
        for name, value in figure["source_sha256"].items():
            assert sha(ROOT / name) == value, name
            figure_sources.add(name)
        for ext, value in figure["files"].items():
            assert sha(ROOT / "results/fast_pzt/f10k/figures" / (figure["figure"] + "." + ext)) == value
    checks["F10K_four_PNG_PDF_and_all_provenance_inputs"] = len(figure_sources)
    gate = js(ROOT / "results/fast_pzt/f10k/analysis/FINAL_GATE.json")
    assert gate["instantaneous_edge_gate"] == "FAIL" and gate["navigation"] == "NOT_ESTABLISHED"
    assert gate["physics"] == "PHYSICS_UNCALIBRATED" and not gate["training_ready"]
    checks["scientific_failures_and_scope_preserved"] = True
    if args.seal:
        old = tree(PARENT)
        index = git("ls-files", "--stage").decode().splitlines()
        changes = []
        for line in index:
            a, name = line.split("\t", 1)
            mode, oid, stage = a.split()
            assert stage == "0"
            if name not in old or old[name]["oid"] != oid:
                assert name == "MEVO_SNAPSHOT.md" or name.startswith(PREFIX)
                changes.append(name)
        assert set(old).issubset({line.split("\t", 1)[1] for line in index})
        omit = {PREFIX + "sync/20261013/PACKAGE_MANIFEST.json", PREFIX + "sync/20261013/VERIFICATION.json"}
        names = sorted(set(changes) - omit)
        files = {n: dict(bytes=len(data(REPO / n)), sha256=sha(REPO / n)) for n in names}
    else:
        files = js(OUT / "PACKAGE_MANIFEST.json")["files"]
    structured = 0
    scanned = 0
    for name, entry in files.items():
        path = REPO / name
        assert sha(path) == entry["sha256"] and len(data(path)) == entry["bytes"], name
        assert entry["bytes"] <= 10_000_000, name
        assert not any(part in {".venv", "__pycache__", ".ssh", ".cache"} for part in Path(name).parts)
        assert path.suffix not in {".mp4", ".avi", ".mov", ".whl", ".pyc"}
        if path.suffix == ".json":
            js(path)
            structured += 1
        if path.suffix in {
            ".py",
            ".json",
            ".jsonl",
            ".txt",
            ".log",
            ".md",
            ".yaml",
            ".yml",
            ".toml",
            ".sh",
            ".diff",
        }:
            scan_credentials(name, data(path))
            scanned += 1
    checks["new_package_files"] = len(files)
    checks["new_package_bytes"] = sum(e["bytes"] for e in files.values())
    checks["JSON_files_parsed"] = structured
    checks["text_files_credential_patterns_scanned"] = scanned
    checks["maximum_new_file_bytes"] = max(e["bytes"] for e in files.values())
    checks["no_unrelated_project_change"] = True
    if args.seal:
        write(
            OUT / "PACKAGE_MANIFEST.json",
            dict(
                parent=PARENT,
                source_delivery=policy["complete_scientific_delivery"],
                scope="New and changed files; inherited parent files unchanged; excludes manifest and verification to avoid circular hashes",
                files=files,
            ),
        )
        write(
            OUT / "VERIFICATION.json",
            dict(
                status="PASS",
                checks=checks,
                date="2026-10-13",
                scientific_computation_started=False,
                manifest_sha256=sha(OUT / "PACKAGE_MANIFEST.json"),
                source_worktree_modified=False,
                note="Packaging verification only. Historical scientific PASS files retain their original full-archive scope; no claim of full raw-data availability.",
            ),
        )
    else:
        assert sha(OUT / "PACKAGE_MANIFEST.json") == js(OUT / "VERIFICATION.json")["manifest_sha256"]
    print(json.dumps(dict(status="PASS", checks=checks), indent=2))


if __name__ == "__main__":
    main()
