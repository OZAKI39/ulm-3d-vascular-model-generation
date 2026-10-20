"""Verify the Git snapshot; optionally restore whole-file gzip transforms."""
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import shutil
import zipfile

ROOT = Path(__file__).resolve().parents[1]
META = ROOT / "github_sync_20261020"


def sha(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(2**20), b""):
            h.update(b)
    return h.hexdigest()


def records():
    for f in sorted(META.glob("EXPORT_RECORDS.*.json.gz")):
        with gzip.open(f, "rt") as stream:
            yield from json.load(stream)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--restore-text", action="store_true",
                   help="Restore whole-file gzip transforms; never overwrite different data")
    a = p.parse_args()
    manifest = json.loads((META / "SNAPSHOT_MANIFEST.json").read_text())
    inherited_link_limits = {
        x["path"]: x["target"] for x in
        json.loads((META / "INHERITED_LINK_LIMITATIONS.json").read_text())
    }
    for row in manifest["files"]:
        path = ROOT / row["path"]
        assert path.resolve().is_relative_to(ROOT.resolve())
        if row.get("kind") == "symlink":
            import os
            assert path.is_symlink() and os.readlink(path) == row["target"], str(path)
            available = path.exists() or Path(str(path.resolve()) + ".gz").exists()
            assert available or inherited_link_limits.get(row["path"]) == row["target"], str(path)
        else:
            assert path.is_file() and path.stat().st_size == row["bytes"] and sha(path) == row["sha256"], str(path)
    compact = restored = compressed = 0
    for row in records():
        if "export" in row:
            assert sha(ROOT / row["export"]) == row["export_sha256"], row["export"]
        if row["action"] in {"copy", "already_present"}:
            assert row["export_sha256"] == row["source_sha256"], row["source"]
        if row["action"] == "compact_state":
            with zipfile.ZipFile(ROOT / row["export"]) as z:
                for name, want in row["fields_sha256"].items():
                    assert hashlib.sha256(z.read(name)).hexdigest() == want, (row["export"], name)
            compact += 1
        if row["action"] in {"gzip", "represented_by_existing_gzip"}:
            h = hashlib.sha256()
            with gzip.open(ROOT / row["export"], "rb") as stream:
                for block in iter(lambda: stream.read(2**20), b""):
                    h.update(block)
            assert h.hexdigest() == row["source_sha256"], row["source"]
            compressed += 1
            if a.restore_text:
                target = ROOT / row["source"]
                assert target.resolve().is_relative_to(ROOT.resolve())
                if target.exists():
                    assert sha(target) == row["source_sha256"], f"Refusing overwrite: {target}"
                else:
                    target.parent.mkdir(parents=True, exist_ok=True)
                    temporary = target.with_name(target.name + ".restoring")
                    with gzip.open(ROOT / row["export"], "rb") as source, temporary.open("wb") as out:
                        shutil.copyfileobj(source, out)
                    assert sha(temporary) == row["source_sha256"]
                    temporary.rename(target)
                    restored += 1
    print(json.dumps({"status": "PASS", "files": len(manifest["files"]),
                      "bytes": sum(x.get("bytes", 0) for x in manifest["files"]),
                      "compact_states": compact, "whole_file_gzip": compressed,
                      "documented_inherited_unavailable_links": len(inherited_link_limits),
                      "restored": restored}, indent=2))


if __name__ == "__main__":
    main()
