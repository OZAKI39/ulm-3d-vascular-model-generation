"""Seal finished Stage2 evidence without modifying the Stage1 manifest."""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

from mevo_mb.io.artifacts import git_commit


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    root = Path.cwd()
    audit = root / "results/stage2/audit"
    tracked = subprocess.check_output(
        ["git", "-c", f"safe.directory={root}", "ls-files", "-z"], text=True
    )
    source_files = [root / name for name in tracked.split("\0") if name]
    identity = {str(p.relative_to(root)): digest(p) for p in sorted(set(source_files))}
    archive = audit / "FINAL_SOURCE.zip"
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as stream:
        for name in identity:
            stream.write(root / name, name)
    source = {
        "created_UTC": datetime.now(timezone.utc).isoformat(),
        "git_commit": git_commit(),
        "source_sha256": identity,
        "archive_sha256": digest(archive),
        "third_party_commits": {
            "apecss": "a0ee071ca247b5489e997d90e83b5f53c40e4ec7",
            "proteus": "fc24fcbcf9623006e8f79989c0126c45a76bc23b",
        },
        "hardware": "results/stage2/system/SERVER_DOCTOR.json",
        "scientific_status": "PARTIALLY_VALIDATED_STAGE2G_FAILED_H_J_NOT_RUN",
    }
    (audit / "FINAL_SOURCE.json").write_text(json.dumps(source, indent=2) + "\n")
    target = audit / "STAGE2_MANIFEST.json"
    files = {
        str(p.relative_to(root)): {"bytes": p.stat().st_size, "sha256": digest(p)}
        for p in sorted((root / "results/stage2").rglob("*"))
        if p.is_file() and p != target
    }
    manifest = {
        "git_commit": git_commit(),
        "file_count": len(files),
        "total_bytes": sum(v["bytes"] for v in files.values()),
        "files": files,
    }
    target.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps({k: v for k, v in manifest.items() if k != "files"}), flush=True)


if __name__ == "__main__":
    main()
