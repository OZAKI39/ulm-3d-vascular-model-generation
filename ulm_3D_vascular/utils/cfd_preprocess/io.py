"""Resolve and validate saved inputs without regenerating upstream data."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from utils.cfd_lumen.roi_io import load_sampling_rois, resolve_sampling_run
from utils.sampling.pipeline import load_models_from_rodent_run
from utils.sampling.sampling_types import GlobalVascularModel, ROIRecord


class InputValidationError(RuntimeError):
    """An immutable upstream input is absent or incompatible."""


class GeometryReferenceError(InputValidationError):
    """The saved Ultraliser geometry is not the required validated surface."""


@dataclass(frozen=True, slots=True)
class GeometryReference:
    run_root: Path
    run_id: str
    roi_id: str
    radius_scale: float
    radius_p95_absolute_relative_error: float
    surface_um_stl: Path
    surface_um_vtp: Path
    surface_m_stl: Path
    sha256: dict[str, str]

    def report(self) -> dict[str, Any]:
        return {
            "status": "PASS",
            "run_id": self.run_id,
            "run_root": str(self.run_root),
            "roi_id": self.roi_id,
            "radius_scale": self.radius_scale,
            "radius_p95_absolute_relative_error": self.radius_p95_absolute_relative_error,
            "source_geometry_modified": False,
            "files": {
                "lumen_surface_um.stl": str(self.surface_um_stl),
                "lumen_surface_um.vtp": str(self.surface_um_vtp),
                "lumen_surface_m.stl": str(self.surface_m_stl),
            },
            "sha256": self.sha256,
        }


def select_roi(
    sampling_run: Path, *, anchor: int | None, roi_id: str | None
) -> ROIRecord:
    rois = load_sampling_rois(sampling_run, selected_only=True)
    matches = [
        roi
        for roi in rois
        if (anchor is not None and roi.anchor_id == anchor)
        or (roi_id is not None and roi.roi_id == roi_id)
    ]
    if len(matches) != 1:
        raise InputValidationError(
            f"Expected exactly one selected ROI, found {len(matches)} for anchor={anchor}, roi_id={roi_id}"
        )
    return matches[0]


def _compatible_rodent_run(candidate: Path, source_model_id: str) -> bool:
    sample = candidate / "samples" / source_model_id
    return (sample / "preprocess_manifest.json").is_file() and (
        sample / "normalized" / "analysis_swc_single_component.npz"
    ).is_file()


def resolve_rodent_run(
    configured: Path | None,
    *,
    output_root: Path,
    source_model_id: str,
) -> tuple[Path, bool]:
    if configured is not None:
        if not _compatible_rodent_run(configured, source_model_id):
            raise InputValidationError(
                f"Configured rodent run is not compatible with {source_model_id}: {configured}"
            )
        return configured.resolve(), False
    candidates = (
        sorted(
            (
                path
                for path in output_root.iterdir()
                if path.is_dir() and _compatible_rodent_run(path, source_model_id)
            ),
            key=lambda path: (path.stat().st_mtime_ns, path.name),
            reverse=True,
        )
        if output_root.is_dir()
        else []
    )
    if not candidates:
        raise InputValidationError(
            f"No compatible rodent run contains {source_model_id} under {output_root}"
        )
    return candidates[0].resolve(), True


def load_matching_global_model(
    run_root: Path, source_model_id: str
) -> GlobalVascularModel:
    matches = [
        model
        for model in load_models_from_rodent_run(run_root)
        if model.source_model_id == source_model_id
    ]
    if len(matches) != 1:
        raise InputValidationError(
            f"Expected exactly one global source model {source_model_id!r}, found {len(matches)}"
        )
    return matches[0]


def verify_global_edge_manifest(
    model: GlobalVascularModel,
    manifest_path: Path,
) -> None:
    """Require exact source-model and parent/current edge identity agreement."""

    if not manifest_path.is_file():
        raise InputValidationError(
            f"GLOBAL_EDGE_MAPPING_MISMATCH: missing {manifest_path}"
        )
    with manifest_path.open(encoding="utf-8-sig", newline="") as stream:
        rows = [
            row
            for row in csv.DictReader(stream)
            if row.get("source_model_id") == model.source_model_id
        ]
    expected = [
        (
            model.source_model_id,
            edge.edge_id,
            edge.upstream_node_id,
            edge.downstream_node_id,
        )
        for edge in model.edges
    ]

    def node_value(row: dict[str, str], base: str) -> int:
        for key in (base, base.replace("_node_id", "_global_node_id")):
            if row.get(key) not in (None, ""):
                return int(row[key])
        raise KeyError(base)

    try:
        observed = [
            (
                row["source_model_id"],
                int(row["global_edge_id"]),
                node_value(row, "upstream_node_id"),
                node_value(row, "downstream_node_id"),
            )
            for row in rows
        ]
    except (KeyError, TypeError, ValueError) as error:
        raise InputValidationError(
            f"GLOBAL_EDGE_MAPPING_MISMATCH: invalid manifest schema: {error}"
        ) from error
    if observed != expected:
        mismatch = next(
            (
                index
                for index, pair in enumerate(zip(observed, expected))
                if pair[0] != pair[1]
            ),
            min(len(observed), len(expected)),
        )
        raise InputValidationError(
            "GLOBAL_EDGE_MAPPING_MISMATCH: "
            f"observed={len(observed)}, expected={len(expected)}, first_mismatch={mismatch}"
        )


def _load_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise GeometryReferenceError(f"CFD_GEOMETRY_REFERENCE_INVALID: missing {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise GeometryReferenceError(
            f"CFD_GEOMETRY_REFERENCE_INVALID: invalid JSON {path}"
        )
    return value


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_model_run(
    run_root: Path,
    *,
    roi_id: str,
    required_radius_scale: float = 0.91,
    maximum_radius_p95_error: float = 0.05,
) -> GeometryReference:
    root = Path(run_root).resolve()
    summary = _load_json(root / "qc" / "run_summary.json")
    surface_qc = _load_json(root / "qc" / "surface_qc.json")
    radius_qc = _load_json(root / "qc" / "radius_fidelity.json")
    reasons: list[str] = []
    if summary.get("status") != "PASS":
        reasons.append("model run status is not PASS")
    if summary.get("roi_id") != roi_id:
        reasons.append(f"ROI mismatch ({summary.get('roi_id')!r} != {roi_id!r})")
    radius_scale = float(summary.get("radius_scale", float("nan")))
    if abs(radius_scale - required_radius_scale) > 1.0e-12:
        reasons.append(f"radius_scale={radius_scale!r}")
    if surface_qc.get("status") != "PASS":
        reasons.append("surface QC is not PASS")
    p95 = float(radius_qc.get("p95_absolute_relative_error", float("nan")))
    if not (p95 <= maximum_radius_p95_error):
        reasons.append(f"radius P95={p95!r}")
    geometry = root / "geometry"
    files = {
        "lumen_surface_um.stl": geometry / "lumen_surface_um.stl",
        "lumen_surface_um.vtp": geometry / "lumen_surface_um.vtp",
        "lumen_surface_m.stl": geometry / "lumen_surface_m.stl",
    }
    missing = [name for name, path in files.items() if not path.is_file()]
    if missing:
        reasons.append(f"missing geometry: {', '.join(missing)}")
    if reasons:
        raise GeometryReferenceError(
            "CFD_GEOMETRY_REFERENCE_INVALID: " + "; ".join(reasons)
        )
    return GeometryReference(
        run_root=root,
        run_id=root.name,
        roi_id=roi_id,
        radius_scale=radius_scale,
        radius_p95_absolute_relative_error=p95,
        surface_um_stl=files["lumen_surface_um.stl"],
        surface_um_vtp=files["lumen_surface_um.vtp"],
        surface_m_stl=files["lumen_surface_m.stl"],
        sha256={name: _sha256(path) for name, path in files.items()},
    )


def resolve_model_run(
    configured: Path | None,
    *,
    model_output_root: Path,
    roi_id: str,
) -> tuple[GeometryReference, bool]:
    if configured is not None:
        return validate_model_run(configured, roi_id=roi_id), False
    candidates = (
        sorted(
            (path for path in model_output_root.iterdir() if path.is_dir()),
            key=lambda path: (path.stat().st_mtime_ns, path.name),
            reverse=True,
        )
        if model_output_root.is_dir()
        else []
    )
    failures: list[str] = []
    for candidate in candidates:
        try:
            return validate_model_run(candidate, roi_id=roi_id), True
        except GeometryReferenceError as error:
            failures.append(str(error))
    detail = failures[0] if failures else "no model run directories"
    raise GeometryReferenceError(
        f"CFD_GEOMETRY_REFERENCE_INVALID: no matching PASS run under {model_output_root}; {detail}"
    )


__all__ = [
    "GeometryReference",
    "GeometryReferenceError",
    "InputValidationError",
    "load_matching_global_model",
    "resolve_model_run",
    "resolve_rodent_run",
    "resolve_sampling_run",
    "select_roi",
    "validate_model_run",
    "verify_global_edge_manifest",
]
