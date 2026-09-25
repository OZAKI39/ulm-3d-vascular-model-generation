"""Unit-explicit JSON and CSV writers for reconstruction evidence."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable


def write_json(path: Path, payload: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def write_csv(
    path: Path,
    rows: Iterable[dict[str, Any]],
    *,
    fieldnames: Iterable[str] | None = None,
) -> Path:
    payload = list(rows)
    names = list(fieldnames or (payload[0].keys() if payload else ()))
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names, extrasaction="ignore")
        if names:
            writer.writeheader()
            writer.writerows(payload)
    return path
