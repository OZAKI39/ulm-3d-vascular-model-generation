"""Versioned Step 1-2 manifests used by independent NNE2 Step 3 runs."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from ..io import write_json


SCHEMA_VERSION = "2.0"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_stack_manifest(
    path: Path,
    *,
    run_root: Path,
    stack_name: str,
    shape_zyx: tuple[int, int, int],
    spacing_xyz_um: tuple[float, float, float],
    artifacts: dict[str, Path],
    cache_key: str,
) -> Path:
    payload = {
        "schema_version": SCHEMA_VERSION,
        "stage": "NNE2 Step 1-2 preprocess",
        "stack_name": stack_name,
        "shape_zyx": list(shape_zyx),
        "spacing_xyz_um": list(spacing_xyz_um),
        "cache_key": cache_key,
        "artifacts": {
            name: {
                "relative_path": file.relative_to(run_root).as_posix(),
                "size_bytes": file.stat().st_size,
                "sha256": sha256_file(file),
            }
            for name, file in artifacts.items()
        },
    }
    write_json(payload, path)
    return path


def load_source_manifests(source_run: Path) -> dict[str, dict[str, Any]]:
    source_run = Path(source_run).expanduser().resolve()
    if not source_run.is_dir():
        raise FileNotFoundError(f"Source run directory does not exist: {source_run}")
    status_path = source_run / "run_status.json"
    if not status_path.is_file():
        raise FileNotFoundError(f"Source run is missing run_status.json: {source_run}")
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("status") not in {"completed", "success", "warning"}:
        raise ValueError(f"Source run is not complete: status={status.get('status')!r}")
    manifests: dict[str, dict[str, Any]] = {}
    pattern = "*/reports/step1_step2_manifest.json"
    for path in sorted((source_run / "stacks").glob(pattern)):
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("schema_version") != SCHEMA_VERSION:
            raise ValueError(f"Unsupported NNE2 manifest schema in {path}")
        payload["manifest_path"] = str(path)
        manifests[str(payload["stack_name"])] = payload
    if not manifests:
        raise FileNotFoundError(f"No completed Step 1-2 stack manifest found in: {source_run}")
    return manifests


def resolve_and_verify_artifacts(
    source_run: Path, manifest: dict[str, Any]
) -> dict[str, Path]:
    source_run = Path(source_run).expanduser().resolve()
    required = {"preprocess_arrays", "stack_metadata", "step1_report", "step2_report"}
    missing_keys = sorted(required - set(manifest.get("artifacts", {})))
    if missing_keys:
        raise ValueError(f"Step 1-2 manifest is missing artifact keys: {missing_keys}")
    output: dict[str, Path] = {}
    for name, metadata in manifest["artifacts"].items():
        path = (source_run / metadata["relative_path"]).resolve()
        try:
            path.relative_to(source_run)
        except ValueError as exc:
            raise ValueError(f"Artifact path leaves source run: {path}") from exc
        if not path.is_file() or path.stat().st_size != int(metadata["size_bytes"]):
            raise FileNotFoundError(f"Missing or size-mismatched Step 1-2 artifact: {path}")
        if sha256_file(path) != metadata["sha256"]:
            raise ValueError(f"Step 1-2 artifact checksum mismatch: {path}")
        output[name] = path
    return output
