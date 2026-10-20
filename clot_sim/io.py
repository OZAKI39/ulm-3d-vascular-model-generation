"""Atomic diagnostic/state output and input integrity checks."""

from pathlib import Path
import hashlib
import json
import os
import numpy as np
from .config import ROOT, INPUTS


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def cv(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def dump(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.json")
    tmp.write_text(json.dumps(value, indent=2, default=cv, allow_nan=False) + "\n")
    os.replace(tmp, path)


def load_npz(path):
    with np.load(path, allow_pickle=False) as data:
        return {key: data[key] for key in data.files}


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp.npz")
    np.savez_compressed(tmp, **value)
    os.replace(tmp, path)


def verify_inputs():
    manifest = json.loads((INPUTS / "manifest.json").read_text())
    for item in manifest:
        path = ROOT / item["destination"]
        if sha(path) != item["sha256"]:
            raise ValueError(f"Input checksum mismatch: {path}")
    return manifest


def verify_native_build():
    """Verify the configured library and current C++ sources against its build receipt."""
    from .config import NATIVE, LIB

    receipt = Path(os.environ.get("CLOT_NATIVE_MANIFEST", NATIVE / "build/BUILD.json"))
    report = json.loads(receipt.read_text())
    if sha(LIB) != report["library_sha256"]:
        raise ValueError(
            "Native library differs from build receipt; rebuild or provide CLOT_NATIVE_MANIFEST"
        )
    for name, want in report["source_sha256"].items():
        if sha(NATIVE / "src" / name) != want:
            raise ValueError(f"Native source changed since build: {name}")
    return report
