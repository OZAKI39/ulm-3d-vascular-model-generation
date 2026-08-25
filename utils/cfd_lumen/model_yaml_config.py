"""Strict YAML configuration for saved-ROI to STL model generation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping

import yaml

from .config import CFDLumenConfig, SurfaceQCConfig, UltraliserConfig


@dataclass(frozen=True, slots=True)
class SWCSTLRunConfig:
    """Validated orchestration and reconstruction settings from one YAML file."""

    source_path: Path
    sampling_run: Path | None
    roi_anchor: int | None
    roi_id: str | None
    surface_backend: str
    output_root: Path
    run_id: str | None
    ultraliser_root: Path
    ultraliser_executable: Path | None
    lumen: CFDLumenConfig


def _mapping(value: Any, label: str) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{label} must be a YAML mapping")
    return dict(value)


def _section(parent: Mapping[str, Any], key: str, allowed: set[str]) -> dict[str, Any]:
    if key not in parent:
        raise ValueError(f"Missing YAML section: {key}")
    values = _mapping(parent[key], key)
    unknown = sorted(set(values) - allowed)
    missing = sorted(allowed - set(values))
    if unknown:
        raise ValueError(f"Unknown keys in {key}: {', '.join(unknown)}")
    if missing:
        raise ValueError(f"Missing keys in {key}: {', '.join(missing)}")
    return values


def _string(value: Any, label: str, *, nullable: bool = False) -> str | None:
    if value is None and nullable:
        return None
    if not isinstance(value, str) or not value.strip():
        suffix = " or null" if nullable else ""
        raise ValueError(f"{label} must be a non-empty string{suffix}")
    return value.strip()


def _integer(value: Any, label: str, *, nullable: bool = False) -> int | None:
    if value is None and nullable:
        return None
    if isinstance(value, bool) or not isinstance(value, int):
        suffix = " or null" if nullable else ""
        raise ValueError(f"{label} must be an integer{suffix}")
    return int(value)


def _resolve_path(
    value: Any,
    label: str,
    project_root: Path,
    *,
    nullable: bool = False,
) -> Path | None:
    text = _string(value, label, nullable=nullable)
    if text is None:
        return None
    path = Path(text).expanduser()
    return path.resolve() if path.is_absolute() else (project_root / path).resolve()


def _typed_dataclass_values(
    values: Mapping[str, Any], defaults: Mapping[str, Any], label: str
) -> dict[str, Any]:
    unknown = sorted(set(values) - set(defaults))
    missing = sorted(set(defaults) - set(values))
    if unknown:
        raise ValueError(f"Unknown keys in {label}: {', '.join(unknown)}")
    if missing:
        raise ValueError(f"Missing keys in {label}: {', '.join(missing)}")
    converted: dict[str, Any] = {}
    for key, default in defaults.items():
        value = values[key]
        field_label = f"{label}.{key}"
        if isinstance(default, bool):
            if not isinstance(value, bool):
                raise ValueError(f"{field_label} must be true or false")
            converted[key] = value
        elif isinstance(default, int):
            if isinstance(value, bool) or not isinstance(value, int):
                raise ValueError(f"{field_label} must be an integer")
            converted[key] = int(value)
        elif isinstance(default, float):
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                raise ValueError(f"{field_label} must be numeric")
            converted[key] = float(value)
        elif isinstance(default, str):
            converted[key] = _string(value, field_label)
        else:
            raise TypeError(f"Unsupported configuration field type for {field_label}")
    return converted


def load_swc_stl_yaml_config(path: Path, *, project_root: Path) -> SWCSTLRunConfig:
    """Load a complete model-generation YAML and reject silent fallbacks."""

    source_path = Path(path).expanduser().resolve()
    if not source_path.is_file():
        raise FileNotFoundError(f"SWC/STL YAML configuration does not exist: {source_path}")
    payload = yaml.safe_load(source_path.read_text(encoding="utf-8"))
    root = _mapping(payload, "YAML root")
    expected_root = {
        "schema_version",
        "paths",
        "selection",
        "reconstruction",
        "ultraliser",
        "surface_qc",
    }
    unknown = sorted(set(root) - expected_root)
    missing = sorted(expected_root - set(root))
    if unknown:
        raise ValueError(f"Unknown top-level YAML keys: {', '.join(unknown)}")
    if missing:
        raise ValueError(f"Missing top-level YAML keys: {', '.join(missing)}")
    schema_version = _integer(root["schema_version"], "schema_version")
    if schema_version != 1:
        raise ValueError(f"Unsupported schema_version: {schema_version}; expected 1")

    paths = _section(
        root,
        "paths",
        {"sampling_run", "output_root", "ultraliser_executable"},
    )
    selection = _section(root, "selection", {"roi_anchor", "roi_id"})
    reconstruction = _section(root, "reconstruction", {"surface_backend", "run_id"})
    ultraliser_values = _mapping(root["ultraliser"], "ultraliser")
    surface_qc_values = _mapping(root["surface_qc"], "surface_qc")

    project_root = Path(project_root).resolve()
    sampling_run = _resolve_path(
        paths["sampling_run"], "paths.sampling_run", project_root, nullable=True
    )
    output_root = _resolve_path(paths["output_root"], "paths.output_root", project_root)
    ultraliser_executable = _resolve_path(
        paths["ultraliser_executable"],
        "paths.ultraliser_executable",
        project_root,
        nullable=True,
    )
    assert output_root is not None

    roi_anchor = _integer(selection["roi_anchor"], "selection.roi_anchor", nullable=True)
    roi_id = _string(selection["roi_id"], "selection.roi_id", nullable=True)
    if (roi_anchor is None) == (roi_id is None):
        raise ValueError("Exactly one of selection.roi_anchor and selection.roi_id must be set")
    if roi_anchor is not None and roi_anchor < 0:
        raise ValueError("selection.roi_anchor must be non-negative")
    surface_backend = str(
        _string(reconstruction["surface_backend"], "reconstruction.surface_backend")
    )
    if surface_backend != "ultraliser":
        raise ValueError("Ultraliser is the only supported surface backend")
    run_id = _string(reconstruction["run_id"], "reconstruction.run_id", nullable=True)

    ultraliser = UltraliserConfig(
        **_typed_dataclass_values(
            ultraliser_values,
            asdict(UltraliserConfig()),
            "ultraliser",
        )
    )
    surface_qc = SurfaceQCConfig(
        **_typed_dataclass_values(
            surface_qc_values,
            asdict(SurfaceQCConfig()),
            "surface_qc",
        )
    )
    lumen = CFDLumenConfig(ultraliser=ultraliser, surface_qc=surface_qc)
    lumen.validate()
    ultraliser_root = _resolve_path(
        ultraliser.ultraliser_root,
        "ultraliser.ultraliser_root",
        project_root,
    )
    assert ultraliser_root is not None
    if not ultraliser_root.is_dir():
        raise FileNotFoundError(f"Ultraliser repository does not exist: {ultraliser_root}")

    return SWCSTLRunConfig(
        source_path=source_path,
        sampling_run=sampling_run,
        roi_anchor=roi_anchor,
        roi_id=roi_id,
        surface_backend=surface_backend,
        output_root=output_root,
        run_id=run_id,
        ultraliser_root=ultraliser_root,
        ultraliser_executable=ultraliser_executable,
        lumen=lumen,
    )
