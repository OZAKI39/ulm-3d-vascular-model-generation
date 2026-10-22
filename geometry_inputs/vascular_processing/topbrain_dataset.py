"""Read-only discovery and physical-coordinate validation of TopBrain NIfTI pairs."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import re

import nibabel as nib
import numpy as np

from .topbrain_labels import LabelMap, load_labelmap
from .topbrain_qc import TopBrainError, sha256

EXPECTED = "imagesTr_topbrain_mr/, labelsTr_topbrain_mr/, imagesTr_topbrain_ct/, labelsTr_topbrain_ct/, itksnap_labelmap_txt/"


@dataclass(frozen=True)
class CasePaths:
    case_id: str
    modality: str
    image: Path
    label: Path
    version: str


def file_identity(path):
    match = re.search(r"(?:topcow|topbrain)_(mr|ct)_(\d+)(?=_|\.)", path.name, re.I)
    if not match:
        raise TopBrainError("CASE_NAME_UNRECOGNIZED", f"Cannot reliably identify case/modality: {path}")
    return match.group(1).lower(), match.group(2).zfill(3)


def discover(root, version="auto"):
    root = Path(root).expanduser().resolve()
    if not root.is_dir():
        raise TopBrainError("TOPBRAIN_DATA_NOT_FOUND", f"{root}; expected {EXPECTED}. Use --topbrain-root; no automatic download.")
    images, labels = {}, {}
    for path in sorted(root.rglob("*.nii*")):
        if not (path.name.endswith(".nii.gz") or path.name.endswith(".nii")):
            continue
        dirs = [p.name.lower() for p in path.parents if p != root.parent]
        image_dirs = [d for d in dirs if d.startswith("imagestr_topbrain_")]
        label_dirs = [d for d in dirs if d.startswith("labelstr_topbrain_")]
        if not image_dirs and not label_dirs:
            continue
        key = file_identity(path)
        for d in image_dirs+label_dirs:
            suffix = re.search(r"_(mr|ct)$", d)
            if suffix and suffix.group(1) != key[0]:
                raise TopBrainError("IMAGE_LABEL_MISMATCH", f"Directory modality disagrees with filename: {path}")
        if image_dirs:
            images.setdefault(key, []).append(path)
        if label_dirs:
            v = "ta36" if any("topaneu36class" in d or "ta36" in d for d in label_dirs) else "2025"
            labels.setdefault((key, v), []).append(path)
    chosen = version
    if chosen == "auto":
        chosen = "ta36" if any(v == "ta36" for _, v in labels) else "2025"
    found = []
    for (key, v), masks in labels.items():
        if v != chosen:
            continue
        scans = images.get(key, [])
        if len(scans) != 1 or len(masks) != 1:
            raise TopBrainError("IMAGE_LABEL_MISMATCH", f"Case {key}: {len(scans)} images, {len(masks)} labels; pairing is not unique")
        found.append(CasePaths(key[1], key[0], scans[0], masks[0], v))
    if not found:
        raise TopBrainError("TOPBRAIN_DATA_NOT_FOUND", f"No paired public training labels in {root}; expected {EXPECTED}")
    return sorted(found, key=lambda item: (item.modality, item.case_id))


def image_metadata(image):
    q, qc = image.get_qform(coded=True)
    s, sc = image.get_sform(coded=True)
    return {"shape": list(image.shape), "dtype": str(image.get_data_dtype()),
            "affine": image.affine.tolist(), "voxel_spacing": list(map(float, image.header.get_zooms()[:3])),
            "orientation": list(nib.aff2axcodes(image.affine)), "spatial_units": image.header.get_xyzt_units()[0],
            "qform": None if q is None else q.tolist(), "qform_code": int(qc),
            "sform": None if s is None else s.tolist(), "sform_code": int(sc)}


@dataclass
class TopBrainCase:
    paths: CasePaths
    image: object
    label_image: object
    labels: np.ndarray
    affine_mm: np.ndarray
    spacing_mm: np.ndarray
    label_map: LabelMap
    provenance: dict

    def mask(self, name):
        return self.labels == self.label_map.resolve(name)


def load_case(root, case_id="001", modality="mr", *, version="auto", labelmap=None, annotation_source="gold_annotation"):
    pairs = discover(root, version)
    selected = [p for p in pairs if p.case_id == str(case_id).zfill(3) and p.modality == modality]
    if len(selected) != 1:
        raise TopBrainError("CASE_NOT_FOUND", f"Requested {modality}/{case_id}; available {[(p.modality,p.case_id) for p in pairs]}")
    p = selected[0]
    scan, mask = nib.load(p.image), nib.load(p.label)
    a, b = image_metadata(scan), image_metadata(mask)
    if (len(scan.shape) != 3 or scan.shape != mask.shape or
        not np.allclose(scan.affine, mask.affine, rtol=0, atol=1e-5) or
        not np.allclose(a["voxel_spacing"], b["voxel_spacing"], rtol=0, atol=1e-6) or
        a["orientation"] != b["orientation"]):
        raise TopBrainError("IMAGE_LABEL_MISMATCH", "Shape/affine/spacing/orientation/units differ; automatic resampling is disabled", image=a, label=b)
    units = a["spatial_units"]
    factors = {"mm": 1., "meter": 1000., "micron": .001}
    unit_warning = None
    if a["spatial_units"] != b["spatial_units"]:
        if b["spatial_units"]=="unknown" and units in factors:
            # Official 2025 labels omit xyzt_units. Their affine/shape/spacing
            # exactly match the paired scan, whose header explicitly states mm.
            unit_warning = "LABEL_UNITS_FROM_PAIRED_IMAGE: segmentation header omits units; using the geometry-matched image header"
        else:
            raise TopBrainError("IMAGE_LABEL_MISMATCH", "Image and label have conflicting physical units")
    if units not in factors:
        raise TopBrainError("NIFTI_UNITS_UNKNOWN", f"Spatial units={units}; do not assume millimetres")
    affine = mask.affine.copy()
    affine[:3] *= factors[units]
    if not np.isfinite(affine).all() or abs(np.linalg.det(affine[:3,:3])) < 1e-15:
        raise TopBrainError("NIFTI_AFFINE_INVALID", "Nonfinite or singular physical transform")
    labels = np.asanyarray(mask.dataobj)
    if not np.isfinite(labels).all() or np.any(labels < 0) or not np.equal(labels, np.floor(labels)).all():
        raise TopBrainError("LABEL_INVALID", "Segmentation must contain finite nonnegative integer labels")
    labels = labels.astype(np.min_scalar_type(int(labels.max())),copy=False)
    labels.flags.writeable = False
    lm = load_labelmap(root, modality, p.version, labelmap)
    if unit_warning:
        lm.warnings.append(unit_warning)
    spacing = np.linalg.norm(affine[:3,:3], axis=0)
    prov = {"dataset": "TopBrain", "dataset_version": p.version, "case_id": p.case_id,
            "modality": p.modality, "annotation_source": annotation_source,
            "image_path": str(p.image), "label_path": str(p.label), "image_sha256": sha256(p.image),
            "label_sha256": sha256(p.label), "image_metadata": a, "label_metadata": b,
            "image_shape": list(scan.shape), "spacing_mm": spacing.tolist(), "affine_mm": affine.tolist(),
            "orientation": a["orientation"], "labelmap": lm.report(), "labelmap_path": str(lm.path) if lm.path else None,
            "effective_spatial_units": units, "unit_resolution": unit_warning or "both NIfTI headers agree",
            "coordinate_transform": "NIfTI affine with header units converted to mm; no flips, resampling or canonicalization",
            "observed_labels": {int(k): int(v) for k,v in zip(*np.unique(labels, return_counts=True))}}
    for parent in p.image.resolve().parents:
        provenance_file = parent/"dataset_provenance.json"
        if provenance_file.is_file():
            archive = json.loads(provenance_file.read_text())
            prov["dataset_release"] = {key:archive.get(key) for key in
                ("dataset","zenodo_record","doi","source_url","archive_md5","licenses")}
            break
    return TopBrainCase(p, scan, mask, labels, affine, spacing, lm, prov)


def save_mask(path, mask, case):
    path = Path(path)
    if path.exists() or path.resolve() in {case.paths.image.resolve(), case.paths.label.resolve()}:
        raise FileExistsError(f"Refusing overwrite: {path}")
    header = case.label_image.header.copy()
    header.set_data_dtype(np.uint8)
    header.set_xyzt_units(case.provenance["effective_spatial_units"])
    result = nib.Nifti1Image(mask.astype(np.uint8), case.label_image.affine, header)
    result.set_qform(case.label_image.get_qform(), int(header["qform_code"]))
    result.set_sform(case.label_image.get_sform(), int(header["sform_code"]))
    nib.save(result, path)
