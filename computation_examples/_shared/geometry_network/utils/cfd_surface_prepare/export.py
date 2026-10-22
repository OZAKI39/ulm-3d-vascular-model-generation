"""Small shared export helpers for CFD surface production."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any, Iterable

from .io import SurfacePrepareError


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> Path:
    records = list(rows)
    if not records:
        raise SurfacePrepareError(f"Cannot write empty CSV: {path}")
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    return path
