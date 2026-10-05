"""Immutable run directories and machine-readable provenance."""
from __future__ import annotations

import json
import subprocess
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
    write_json(run / "metadata.json", {"status": "CREATED", "git_commit": git_commit(), "case": name})
    return run
