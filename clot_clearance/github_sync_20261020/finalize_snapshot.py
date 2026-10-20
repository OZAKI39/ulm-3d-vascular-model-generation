"""Inventory a completed export. Run after all transfers and document edits."""
from pathlib import Path
from collections import Counter
import gzip
import hashlib
import json
import os
import re
import subprocess
from export_files import sha, gzip_write

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "github_sync_20261020"
SOURCE = Path("/home/lzy/projects/temp_storage/clot_github_sync_20261015/clot_clearance")


def dump(name, obj):
    (META / name).write_text(json.dumps(obj, ensure_ascii=False, indent=2) + "\n")


def main():
    tracked = set(subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0"))
    links = []
    for base, dirs, files in os.walk(SOURCE):
        dirs[:] = [d for d in dirs if d not in {".git", "__pycache__"}]
        for name in files:
            p = Path(base) / name
            rel = str(p.relative_to(SOURCE))
            if not p.is_symlink() or rel in tracked:
                continue
            target = os.readlink(p)
            row = dict(path=rel, source_target=target)
            out = ROOT / rel
            resolved = (out.parent / target).resolve()
            if not resolved.is_relative_to(ROOT) or ".so" in name:
                row.update(action="exclude", reason="binary or external link")
            elif out.exists() and not out.is_symlink():
                row["action"] = "regular_file_supplied"
            elif resolved.is_file() or Path(str(resolved) + ".gz").is_file():
                out.parent.mkdir(parents=True, exist_ok=True)
                if not out.is_symlink():
                    out.symlink_to(target)
                assert os.readlink(out) == target
                row.update(action="preserve", requires_text_restore=not resolved.is_file())
            else:
                row.update(action="exclude", reason="target retained only in source/compact archive")
            links.append(row)
    dump("SYMLINKS.json", links)

    rows = []
    for f in sorted(META.glob("EXPORT_RECORDS.*.json.gz")):
        with gzip.open(f, "rt") as stream:
            part = json.load(stream)
        for row in part:
            if row["source"] == "ACTIVE_DEVELOPMENT_ENTRYPOINTS.md":
                row.update(action="documentation_update", export_sha256=sha(ROOT / row["export"]),
                           export_bytes=(ROOT / row["export"]).stat().st_size,
                           reason="Added current snapshot and fixed-camera replay navigation")
        gzip_write(f, json.dumps(part, ensure_ascii=False, indent=2).encode() + b"\n")
        rows.extend(part)
    dump("SYNC_SUMMARY.json", {
        "base_commit": "5036268",
        "branch": "clot-clearance/project-sync-20261020",
        "actions": dict(Counter(x["action"] for x in rows)),
        "exported_payload_bytes": sum(x.get("export_bytes", 0) for x in rows),
        "excluded_source_bytes": sum(x.get("source_bytes", 0) for x in rows if x["action"] == "exclude"),
        "compact_state_count": sum(x["action"] == "compact_state" for x in rows),
        "source_files_deleted": 0, "new_physics_steps": 0,
        "symlinks": dict(Counter(x["action"] for x in links)),
    })
    dump("SOURCE_COVERAGE.json", {
        "local_source": str(SOURCE), "server_source": "/workspace/y_clot/project",
        "local_inventory": "clot_local_inventory.json.gz",
        "server_inventory": "clot_server_inventory.json.gz",
        "source_records": ["EXPORT_RECORDS.local.json.gz", "EXPORT_RECORDS.server.json.gz"],
        "common_file_comparison": "SOURCE_COMMON_FILE_COMPARISON.json",
        "inherited_git_base": "5036268", "previous_remote_snapshot": "450b37f",
        "older_server_exclusions": "../github_sync_20261015/EXCLUDED_FILES.json.gz",
        "scope": "Inherited accepted snapshots plus all new local research files and server-only additions for current research tasks; caches/build/upstream copies and dense redundant payload excluded explicitly.",
        "documentation_choice": "Preserve curated Git README/AGENTS/path adapter; update current navigation separately from scientific evidence.",
    })
    exclusions = {"SNAPSHOT_MANIFEST.json", "SNAPSHOT_MANIFEST.sha256", "VALIDATION.json",
                  "PUSH_RECEIPT.json"}
    files = []
    secret_hits = []
    patterns = [
        re.compile(rb"-----BEGIN (?:RSA |OPENSSH |EC )?PRIVATE KEY-----"),
        re.compile(rb"(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{60,})"),
    ]
    text_suffixes = {".py", ".cpp", ".h", ".md", ".txt", ".json", ".log",
                     ".sh", ".conf", ".yaml", ".yml", ".patch", ".diff"}
    for base, dirs, names in os.walk(ROOT):
        dirs[:] = [d for d in dirs if d not in {".git", "__pycache__", ".pytest_cache"}]
        for name in names:
            p = Path(base) / name
            rel = str(p.relative_to(ROOT))
            if p.parent == META and name in exclusions:
                continue
            if p.is_symlink():
                files.append(dict(path=rel, kind="symlink", target=os.readlink(p)))
            elif p.is_file():
                if rel not in tracked:
                    assert p.stat().st_size < 32 * 1024**2, rel
                    if p.suffix in text_suffixes:
                        b = p.read_bytes()
                        if any(pattern.search(b) for pattern in patterns):
                            secret_hits.append(rel)
                files.append(dict(path=rel, kind="file", bytes=p.stat().st_size, sha256=sha(p)))
    assert not secret_hits, secret_hits
    files.sort(key=lambda x: x["path"])
    dump("SNAPSHOT_MANIFEST.json", dict(base_commit="5036268", files=files,
         self_exclusions=sorted(exclusions), generated_caches_excluded=True))
    (META / "SNAPSHOT_MANIFEST.sha256").write_text(sha(META / "SNAPSHOT_MANIFEST.json") + "  SNAPSHOT_MANIFEST.json\n")
    print(json.dumps({"files": len(files), "bytes": sum(x.get("bytes", 0) for x in files),
                      "new_text_secret_scan": "PASS"}, indent=2))


if __name__ == "__main__":
    main()
