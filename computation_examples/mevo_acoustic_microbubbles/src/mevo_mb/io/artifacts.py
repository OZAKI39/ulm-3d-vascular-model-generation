"""Immutable run directories and machine-readable provenance."""
from __future__ import annotations

import json
import hashlib
import subprocess
import zipfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import yaml


def write_json(path: Path, data: dict) -> None:
    def convert(value):
        if isinstance(value, np.ndarray):
            return value.tolist()
        if isinstance(value, np.generic):
            return value.item()
        if isinstance(value, Path):
            return str(value)
        raise TypeError(type(value).__name__)
    path.write_text(json.dumps(data, indent=2, default=convert, allow_nan=False) + "\n")


def git_commit() -> str | None:
    r = subprocess.run(["git", "-c", f"safe.directory={Path.cwd()}", "rev-parse", "HEAD"], capture_output=True, text=True, check=False)
    return r.stdout.strip() if r.returncode == 0 else None


def new_run(root: Path, name: str, config: dict) -> Path:
    run = root / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S_%fZ_") + name)
    run.mkdir(parents=True, exist_ok=False)
    for sub in ["logs", "acoustic", "particles", "vtk", "figures"]:
        (run / sub).mkdir()
    (run / "config.yaml").write_text(yaml.safe_dump(config, sort_keys=False))
    project = Path(__file__).resolve().parents[3]
    paths = sorted(p for directory in ["src", "configs", "scripts", "tests"] for p in (project/directory).rglob("*")
                   if p.is_file() and p.suffix in {".py", ".yaml"} and "__pycache__" not in p.parts)
    identities = {str(p.relative_to(project)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    write_json(run/"SOURCE_IDENTITY.json", {"git_commit": git_commit(), "source_sha256": identities,
               "note": "Source archive and hashes define exact runtime code even with uncommitted development changes."})
    with zipfile.ZipFile(run/"logs/source_snapshot.zip", "x", zipfile.ZIP_DEFLATED) as archive:
        for p in paths:
            archive.write(p, p.relative_to(project))
    write_json(run / "metadata.json", {"status": "CREATED", "git_commit": git_commit(), "case": name})
    return run
