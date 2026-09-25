"""Anatomical intent is supplied by a person; no branch-order/size classifier."""
from datetime import datetime, timezone
import hashlib
import os
from pathlib import Path

import yaml

ROI_DEFINITIONS = {f"{side}{territory}_{''.join(segments)}": (side, territory, segments)
                   for territory, segments in (("MCA", ["M2", "M3"]), ("ACA", ["A2", "A3"]), ("PCA", ["P2", "P3"]))
                   for side in ("L", "R")}


class LandmarkError(ValueError):
    def __init__(self, message, *, context=None):
        super().__init__(message)
        self.context = context or {}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def landmark_template(source: Path, output: Path, *, units="mm", subject_id=None):
    return {"version": 1, "source": {"path": os.path.relpath(source.resolve(), output.parent.resolve()),
                "filename": source.name, "sha256": sha256(source), "subject_id": subject_id or source.stem,
                "units": units}, "annotation_date": datetime.now(timezone.utc).isoformat(),
            "annotator": "", "rois": {name: {"enabled": False, "side": side, "territory": territory,
                "segments": segments.copy(), "proximal_nodes": [], "distal_nodes": [],
                "exclude_subtree_roots": [], "proximal_boundary_status": "unknown",
                "manually_verified": False, "notes": "", "internal_landmarks": {}}
                for name, (side, territory, segments) in ROI_DEFINITIONS.items()}}


def write_landmarks(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as stream:
        yaml.safe_dump(payload, stream, sort_keys=False, allow_unicode=True)


def node_list(value, label):
    if not isinstance(value, list) or any(type(n) is not int or n < 0 for n in value):
        raise LandmarkError(f"{label} must be a list of original nonnegative integer SWC IDs")
    if len(set(value)) != len(value):
        raise LandmarkError(f"{label} contains duplicate IDs")
    return value


def load_landmarks(path: Path, source: Path, *, allow_source_mismatch=False):
    doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(doc, dict) or type(doc.get("version")) is not int or doc["version"] != 1:
        raise LandmarkError("Landmark YAML version must be 1")
    info = doc.get("source")
    if not isinstance(info, dict) or not isinstance(info.get("sha256"), str) or len(info["sha256"]) != 64:
        raise LandmarkError("source.sha256 (64 hex digits) is required")
    try:
        int(info["sha256"], 16)
    except ValueError as exc:
        raise LandmarkError("source.sha256 must contain hex digits") from exc
    if info.get("units", "mm") not in {"mm", "um"}:
        raise LandmarkError("source.units must be mm or um")
    mismatch = info["sha256"].lower() != sha256(source) or info.get("filename", source.name) != source.name
    warnings = []
    if mismatch:
        message = "SOURCE MISMATCH: landmark hash/filename does not match input; anatomical verification is invalid for this source"
        if not allow_source_mismatch:
            raise LandmarkError(message)
        warnings.append(message)
    rois = doc.get("rois")
    if not isinstance(rois, dict):
        raise LandmarkError("rois must be a mapping")
    for name, roi in rois.items():
        if name not in ROI_DEFINITIONS or not isinstance(roi, dict):
            raise LandmarkError(f"Unsupported ROI definition: {name}")
        side, territory, segments = ROI_DEFINITIONS[name]
        if (roi.get("side"), roi.get("territory"), roi.get("segments")) != (side, territory, segments):
            raise LandmarkError(f"{name}: side/territory/segments must be {side}/{territory}/{segments}")
        for key in ("enabled", "manually_verified"):
            if type(roi.get(key)) is not bool:
                raise LandmarkError(f"{name}.{key} must be true or false")
        for key in ("proximal_nodes", "distal_nodes", "exclude_subtree_roots"):
            node_list(roi.get(key), f"{name}.{key}")
        if roi.get("proximal_boundary_status") not in {"unknown", "confirmed", "truncated_to_available_data"}:
            raise LandmarkError(f"{name}: declare proximal_boundary_status as confirmed, truncated_to_available_data or unknown")
        if not isinstance(roi.get("internal_landmarks", {}), dict):
            raise LandmarkError(f"{name}.internal_landmarks must be a mapping")
        for label, nodes in roi.get("internal_landmarks", {}).items():
            node_list(nodes, f"{name}.internal_landmarks.{label}")
    return doc, warnings
