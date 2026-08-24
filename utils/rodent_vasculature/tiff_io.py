"""LZW-capable multipage TIFF inspection and loading."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from PIL import Image


@dataclass(frozen=True, slots=True)
class TiffMetadata:
    path: Path
    shape_zyx: tuple[int, int, int]
    dtype: str
    mode: str
    compression: int | str | None

    def to_dict(self) -> dict[str, Any]:
        return {
            "path": str(self.path),
            "shape_zyx": self.shape_zyx,
            "dtype": self.dtype,
            "mode": self.mode,
            "compression": self.compression,
            "array_axis_order": ["z", "y", "x"],
        }


def inspect_tiff(path: Path) -> TiffMetadata:
    with Image.open(path) as image:
        frame = np.asarray(image)
        shape = (int(getattr(image, "n_frames", 1)), int(image.height), int(image.width))
        compression = image.tag_v2.get(259) if hasattr(image, "tag_v2") else None
        return TiffMetadata(path, shape, str(frame.dtype), image.mode, compression)


def load_tiff_volume(path: Path) -> np.ndarray:
    frames: list[np.ndarray] = []
    with Image.open(path) as image:
        for index in range(int(getattr(image, "n_frames", 1))):
            image.seek(index)
            frames.append(np.asarray(image).copy())
    if not frames:
        raise ValueError(f"TIFF contains no frames: {path}")
    volume = np.stack(frames, axis=0)
    if volume.ndim != 3:
        raise ValueError(f"Expected scalar 3D TIFF, found shape {volume.shape}: {path}")
    return volume


def save_tiff_volume(volume_zyx: np.ndarray, path: Path) -> Path:
    """Save a scalar 3-D array as a lossless LZW multipage TIFF."""

    volume = np.asarray(volume_zyx)
    if volume.ndim != 3 or not len(volume):
        raise ValueError("A non-empty scalar 3-D array is required")
    path.parent.mkdir(parents=True, exist_ok=True)
    frames = [Image.fromarray(np.asarray(frame)) for frame in volume]
    frames[0].save(
        path,
        save_all=True,
        append_images=frames[1:],
        compression="tiff_lzw",
    )
    return path


def save_normalized_volume(
    image_volume_zyx: np.ndarray,
    mask_volume_zyx: np.ndarray,
    path: Path,
) -> Path:
    """Save an unchanged auxiliary image/Mask pair for display and optional QC."""

    image = np.asarray(image_volume_zyx)
    mask = np.asarray(mask_volume_zyx)
    if image.ndim != 3 or image.shape != mask.shape:
        raise ValueError("Normalized image and mask must be matching 3-D arrays")
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, image_volume_zyx=image, mask_volume_zyx=mask)
    return path


def load_normalized_volume(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Load an auxiliary image/Mask pair saved by :func:`save_normalized_volume`."""

    with np.load(path, allow_pickle=False) as arrays:
        image = np.asarray(arrays["image_volume_zyx"]).copy()
        mask = np.asarray(arrays["mask_volume_zyx"]).copy()
    if image.ndim != 3 or image.shape != mask.shape:
        raise ValueError(f"Invalid normalized volume artifact: {path}")
    return image, mask


def volume_statistics(volume: np.ndarray) -> dict[str, Any]:
    values = np.asarray(volume)
    return {
        "shape_zyx": values.shape,
        "dtype": str(values.dtype),
        "minimum": float(np.min(values)),
        "maximum": float(np.max(values)),
        "mean": float(np.mean(values)),
        "nonzero_fraction": float(np.count_nonzero(values) / values.size),
    }
