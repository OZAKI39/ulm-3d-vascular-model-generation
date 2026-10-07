"""Inventory actual files, strict JSON audit and size-bounded Git selection.

Run inventory on server after experiment completion, then locally after syncing
and rendering. Missing server assets are listed as external, never as local files.
"""

import argparse
import gzip
import json
import subprocess
from scripts.stage5a_r2.seal import ROOT, BASE
from scripts.stage4l.dense_probe import sha, write

EXCLUDED = {
    "ARTIFACT_MANIFEST.json",
    "SERVER_ARTIFACT_MANIFEST.json",
    "GIT_ARTIFACTS.txt",
    "JSON_AUDIT.json",
    "DELIVERY_IDENTITY.json",
}
HEAVY_NAMES = {
    "SOURCE.zip",
    "COMPLETED_TRIALS.json",
    "RESOURCE_SAMPLES.json",
    "ACCEPTED.jsonl",
    "REJECTED.jsonl",
    "TRANSLATION_SOLVES.jsonl",
}


def strict(path):
    def bad(value):
        raise ValueError(f"Nonstandard JSON number {value}")

    if path.suffix == ".jsonl":
        for line in path.read_text().splitlines():
            if line:
                json.loads(line, parse_constant=bad)
    else:
        json.loads(path.read_text(), parse_constant=bad)


def compact_ledgers():
    """Lossless copies, not subsampling. Raw ledgers remain untouched."""
    root = BASE / "final/compact_ledgers"
    root.mkdir(exist_ok=True)
    for folder in [BASE / "runs", BASE / "benchmarks/frozen_reference"]:
        for p in folder.glob("*/*"):
            if p.name not in HEAVY_NAMES - {"SOURCE.zip"}:
                continue
            out = root / p.parent.name / (p.name + ".gz")
            out.parent.mkdir(exist_ok=True)
            raw = p.read_bytes()
            compressed = gzip.compress(raw, mtime=0)
            if not out.exists() or out.read_bytes() != compressed:
                out.write_bytes(compressed)


def main(server=False):
    if not server:
        compact_ledgers()
    rows = []
    invalid = []
    checked = 0
    for p in sorted(BASE.rglob("*")):
        if not p.is_file() or p.name in EXCLUDED or "__pycache__" in p.parts or ".runtime" in p.parts:
            continue
        if p.suffix in [".so", ".pyc"] or any("compiled" in s for s in p.parts):
            continue
        if p.suffix in [".json", ".jsonl"]:
            try:
                strict(p)
                checked += 1
            except Exception as e:
                invalid.append(dict(path=str(p), error=repr(e)))
        rel = p.relative_to(ROOT)
        size = p.stat().st_size
        track = bool(
            size <= 4 * 2**20
            and p.suffix not in [".mp4", ".h5", ".zip"]
            and p.name not in HEAVY_NAMES
            and not p.name.startswith("PREVIEW_")
        )
        rows.append(
            dict(
                relative_path=str(rel),
                actual_path=str(p),
                bytes=size,
                sha256=sha(p),
                git_selected=track,
                reason="Necessary bounded evidence"
                if track
                else "Large/redundant raw asset retained by path and hash",
            )
        )
    write(
        BASE / "final/JSON_AUDIT.json",
        dict(status="PASS" if not invalid else "FAIL", checked_files=checked, failures=invalid),
    )
    assert not invalid, invalid[:3]
    if server:
        write(BASE / "final/SERVER_ARTIFACT_MANIFEST.json", dict(root=str(ROOT), files=rows))
        print("Server files:", len(rows))
        return
    remote_path = BASE / "final/SERVER_ARTIFACT_MANIFEST.json"
    remote = json.loads(remote_path.read_text()) if remote_path.exists() else None
    current = {r["relative_path"]: r for r in rows}
    external = []
    if remote:
        for r in remote["files"]:
            q = current.get(r["relative_path"])
            if q:
                q["server_path"] = r["actual_path"]
                q["server_sha256"] = r["sha256"]
                q["server_copy_identical"] = q["sha256"] == r["sha256"]
            else:
                external.append(dict(r, git_selected=False, local_copy_exists=False))
    write(
        BASE / "final/ARTIFACT_MANIFEST.json",
        dict(
            root=str(ROOT),
            files=rows,
            server_only_files=external,
            excluded_from_self_hash=sorted(EXCLUDED),
            large_asset_policy="Full HDF5/source ZIP/MP4/per-step raw ledgers stay at real paths. Losslessly compressed ledgers, checkpoints, fixtures, sampled actual frames and figures are included when <=4 MiB each. Scientific input fields remain external, with identities in each run.",
            source_of_remote_inventory=str(remote_path) if remote else None,
        ),
    )
    selected = [r["relative_path"] for r in rows if r["git_selected"]]
    selected += [
        str((BASE / "final" / name).relative_to(ROOT))
        for name in ["ARTIFACT_MANIFEST.json", "JSON_AUDIT.json"]
    ]
    (BASE / "final/GIT_ARTIFACTS.txt").write_text("\n".join(selected) + "\n")
    parent = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    write(
        BASE / "final/DELIVERY_IDENTITY.json",
        dict(
            executed_source_commit=(ROOT / "SOURCE_COMMIT_R2.txt").read_text().strip(),
            archive_parent_commit=parent,
            final_archive_commit="Resolve the commit containing this file; no recursive self-commit hash.",
            branch=subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
            base_commit="7c39828d93509fff5551b7a16537ca3aa4ee9c23",
            protocol_sha256=json.loads((BASE / "final/MOTION_SOLVER_R2_GATE.json").read_text())[
                "frozen_protocol_sha256"
            ],
            environment_sha256=sha(BASE / "environment/IDENTITY.json"),
        ),
    )
    print("Local files:", len(rows), "external-only:", len(external), "Git evidence:", len(selected))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--server", action="store_true")
    main(p.parse_args().server)
