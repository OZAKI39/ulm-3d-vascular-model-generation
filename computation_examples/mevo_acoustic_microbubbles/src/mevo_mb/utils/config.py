"""YAML loading with explicit SI units and fail-fast physical constraints."""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import yaml


def merge(base: dict, override: dict) -> dict:
    result = copy.deepcopy(base)
    for key, value in override.items():
        result[key] = merge(result[key], value) if isinstance(value, dict) and isinstance(result.get(key), dict) else value
    return result


def load_config(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    data = yaml.safe_load(path.read_text())
    if not isinstance(data, dict):
        raise ValueError("Configuration must be a mapping")
    if "extends" in data:
        data = merge(load_config(path.parent / data.pop("extends")), data)
    validate_config(data)
    return data


def validate_config(c: dict) -> None:
    if c["acoustics"].get("BonA") is not None:
        raise ValueError("Nonlinear BonA propagation is incompatible with acoustic basis superposition")
    for value in [*c["grid"]["shape"], *c["grid"]["spacing_m"], c["medium"]["sound_speed_m_s"],
                  c["medium"]["density_kg_m3"], c["medium"]["viscosity_Pa_s"], c["acoustics"]["frequency_Hz"],
                  c["vessel"]["radius_m"], c["vessel"]["length_m"]]:
        if not np.isfinite(value) or value <= 0:
            raise ValueError("Grid, medium, frequency and vessel parameters must be positive finite SI values")
    if any(int(n) != n for n in c["grid"]["shape"]) or len(c["grid"]["shape"]) != 3:
        raise ValueError("Grid must have three integer dimensions")
    ids = [p["id"] for p in c["pzt"]]
    if not ids or len(set(ids)) != len(ids):
        raise ValueError("Source IDs must be unique")
    if not 0 < c["acoustics"]["CFL"] <= 0.3:
        raise ValueError("Demo CFL must lie in (0, 0.3]")
    if c["particles"]["dt_s"] <= 0 or c["particles"]["duration_s"] <= 0:
        raise ValueError("Positive particle duration and dt required")


def fingerprint(data: dict) -> str:
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
