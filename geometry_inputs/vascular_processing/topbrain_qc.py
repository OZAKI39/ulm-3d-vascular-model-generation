"""Small TopBrain-specific QC utilities; no changes to the existing pipelines."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

ALGORITHM_VERSION = "topbrain-mevo-1.0"
CONNECTIVITY = np.ones((3, 3, 3), dtype=bool)


class TopBrainError(ValueError):
    def __init__(self, code, message, **details):
        self.code, self.details = code, details
        super().__init__(f"{code}: {message}")


def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    def convert(obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.generic):
            return obj.item()
        if isinstance(obj, Path):
            return str(obj)
        raise TypeError(type(obj).__name__)
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, default=convert, allow_nan=False)+"\n")


def components(mask, spacing):
    labels, count = ndi.label(mask, structure=CONNECTIVITY)
    sizes = np.bincount(labels.ravel())[1:]
    return labels, {"component_count": int(count), "component_voxel_counts": sizes.tolist(),
                    "component_volume_mm3": (sizes * np.prod(spacing)).tolist(),
                    "voxel_count": int(mask.sum())}


def crop(mask, margin=1):
    points = np.argwhere(mask)
    if not len(points):
        raise TopBrainError("EMPTY_ROI", "No foreground voxels")
    low = np.maximum(points.min(axis=0)-margin, 0)
    high = np.minimum(points.max(axis=0)+margin+1, mask.shape)
    slices = tuple(slice(int(a), int(b)) for a, b in zip(low, high))
    return slices, low


def physical_points(voxels, affine_mm):
    return np.asarray(voxels, dtype=float) @ affine_mm[:3, :3].T + affine_mm[:3, 3]


def contact(a, b):
    """Foreground voxels in a that touch b in the original 26-neighbour grid."""
    result = np.zeros_like(a,dtype=bool)
    if not b.any():
        return result
    slices,_ = crop(b,margin=1)
    result[slices] = a[slices] & ndi.binary_dilation(b[slices],structure=CONNECTIVITY)
    return result
