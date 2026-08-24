"""Dataset-Mask loading and SWC-to-Mask reconstruction for v5-1 experiments.

This module is deliberately independent of the formal lumen builder.  The
reconstruction implemented here approximates the paper's one-voxel skeleton
resampling and local-radius-aware Gaussian/tubular annotation procedure; it is
not a replacement for the SWC-only CFD reconstruction.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import cv2
import numpy as np
import tifffile


SPACING_XYZ_UM = (1.0, 1.0, 2.0)
SPACING_ZYX_UM = SPACING_XYZ_UM[::-1]
PAPER_FOREGROUND_THRESHOLD = 0.404
PAPER_FOREGROUND_UINT8 = int(round(255.0 * PAPER_FOREGROUND_THRESHOLD))


@dataclass(frozen=True, slots=True)
class SWCReference:
    """Minimal, source-faithful SWC representation used by the experiment."""

    path: Path
    node_ids: np.ndarray
    points_voxel_xyz: np.ndarray
    radius_um: np.ndarray
    parent_ids: np.ndarray

    @property
    def points_um(self) -> np.ndarray:
        return self.points_voxel_xyz * np.asarray(SPACING_XYZ_UM, dtype=float)


@dataclass(frozen=True, slots=True)
class RadiusCalibration:
    """Held-out annotation-footprint calibration for the approximate generator."""

    radius_scale: float = 0.35
    voxel_floor: float = 0.75
    threshold: float = PAPER_FOREGROUND_THRESHOLD
    coordinate_metric: str = "voxel_xyz"
    calibration_split: str = "train"

    def effective_radius_voxel(self, radius_um: np.ndarray | float) -> np.ndarray:
        return self.radius_scale * np.asarray(radius_um, dtype=float) + self.voxel_floor


def read_tiff_volume(path: Path) -> np.ndarray:
    """Read an LZW multi-page TIFF as a z-y-x NumPy array."""

    success, frames = cv2.imreadmulti(str(Path(path)), flags=cv2.IMREAD_UNCHANGED)
    if not success or not frames:
        raise ValueError(f"Unable to read TIFF stack: {path}")
    array = np.stack(frames, axis=0)
    if array.ndim != 3:
        raise ValueError(f"Expected a 3-D TIFF stack, got {array.shape}: {path}")
    return array


def write_tiff_volume(path: Path, array: np.ndarray) -> Path:
    """Write a z-y-x TIFF without requiring an optional compression codec."""

    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    tifffile.imwrite(
        target,
        np.asarray(array),
        photometric="minisblack",
        metadata={"axes": "ZYX", "spacing_xyz_um": SPACING_XYZ_UM},
    )
    return target


def read_swc_reference(path: Path) -> SWCReference:
    """Read the seven canonical SWC columns without changing topology."""

    rows: list[tuple[int, float, float, float, float, int]] = []
    for raw_line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 7:
            continue
        rows.append(
            (
                int(round(float(fields[0]))),
                float(fields[2]),
                float(fields[3]),
                float(fields[4]),
                float(fields[5]),
                int(round(float(fields[6]))),
            )
        )
    if not rows:
        raise ValueError(f"SWC contains no usable records: {path}")
    array = np.asarray(rows, dtype=float)
    return SWCReference(
        path=Path(path).resolve(),
        node_ids=array[:, 0].astype(np.int64),
        points_voxel_xyz=array[:, 1:4].astype(float),
        radius_um=array[:, 4].astype(float),
        parent_ids=array[:, 5].astype(np.int64),
    )


def resample_swc_one_voxel(swc: SWCReference) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Uniformly sample every SWC edge at no more than one voxel interval.

    Sampling uses the annotation lattice because the paper describes one-voxel
    sampling.  Physical 1x1x2 um spacing is applied later to distance metrics and
    surfaces.  Returned labels identify the nearest source child node.
    """

    index_by_id = {int(node_id): index for index, node_id in enumerate(swc.node_ids)}
    point_blocks: list[np.ndarray] = []
    radius_blocks: list[np.ndarray] = []
    owner_blocks: list[np.ndarray] = []
    for index, parent_id in enumerate(swc.parent_ids):
        if int(parent_id) < 0 or int(parent_id) not in index_by_id:
            point_blocks.append(swc.points_voxel_xyz[index : index + 1])
            radius_blocks.append(swc.radius_um[index : index + 1])
            owner_blocks.append(np.asarray((index,), dtype=np.int64))
            continue
        parent_index = index_by_id[int(parent_id)]
        first = swc.points_voxel_xyz[parent_index]
        second = swc.points_voxel_xyz[index]
        length = float(np.linalg.norm(second - first))
        interval_count = max(1, int(np.ceil(length)))
        fractions = np.linspace(0.0, 1.0, interval_count + 1)[1:]
        point_blocks.append(first[None, :] + fractions[:, None] * (second - first)[None, :])
        radius_blocks.append(
            swc.radius_um[parent_index]
            + fractions * (swc.radius_um[index] - swc.radius_um[parent_index])
        )
        owner_blocks.append(np.full(len(fractions), index, dtype=np.int64))
    return (
        np.concatenate(point_blocks, axis=0),
        np.concatenate(radius_blocks, axis=0),
        np.concatenate(owner_blocks, axis=0),
    )


def _splat_gaussian(
    field: np.ndarray,
    point_xyz: np.ndarray,
    effective_radius_voxel: float,
    *,
    minimum_probability: float,
) -> None:
    # exp(-c) = the paper's 0.404 foreground threshold at d=r_eff.
    coefficient = -float(np.log(PAPER_FOREGROUND_THRESHOLD))
    extent = effective_radius_voxel * np.sqrt(
        -np.log(max(minimum_probability, 1.0e-12)) / coefficient
    )
    shape_xyz = np.asarray(field.shape[::-1], dtype=int)
    lower = np.maximum(np.floor(point_xyz - extent).astype(int), 0)
    upper = np.minimum(np.ceil(point_xyz + extent).astype(int) + 1, shape_xyz)
    if np.any(upper <= lower):
        return
    x = np.arange(lower[0], upper[0], dtype=float)
    y = np.arange(lower[1], upper[1], dtype=float)
    z = np.arange(lower[2], upper[2], dtype=float)
    zz, yy, xx = np.meshgrid(z, y, x, indexing="ij")
    squared = (xx - point_xyz[0]) ** 2 + (yy - point_xyz[1]) ** 2 + (
        zz - point_xyz[2]
    ) ** 2
    probability = np.exp(-coefficient * squared / effective_radius_voxel**2)
    view = field[lower[2] : upper[2], lower[1] : upper[1], lower[0] : upper[0]]
    np.maximum(view, probability.astype(field.dtype, copy=False), out=view)


def reconstruct_mask_from_swc(
    swc: SWCReference,
    shape_zyx: tuple[int, int, int],
    *,
    calibration: RadiusCalibration = RadiusCalibration(),
    minimum_probability: float = 1.0 / 255.0,
) -> tuple[np.ndarray, dict[str, object]]:
    """Approximate the skeleton+radius-derived morphology probability map.

    ``radius_scale`` and ``voxel_floor`` describe the annotation footprint, not
    a revised biological radius.  Their defaults were fixed on three training
    blocks and are evaluated on the current held-out test block.
    """

    points, radii, _ = resample_swc_one_voxel(swc)
    probability = np.zeros(tuple(map(int, shape_zyx)), dtype=np.float32)
    effective = calibration.effective_radius_voxel(radii)
    for point, radius in zip(points, effective):
        _splat_gaussian(
            probability,
            point,
            float(radius),
            minimum_probability=minimum_probability,
        )
    reconstructed = np.rint(np.clip(probability, 0.0, 1.0) * 255.0).astype(np.uint8)
    metadata: dict[str, object] = {
        "method": "approximated_one_voxel_resampling_plus_radius_aware_gaussian_tubes",
        "formal_pipeline": False,
        "experimental_only": True,
        "shape_zyx": list(map(int, shape_zyx)),
        "spacing_xyz_um": list(SPACING_XYZ_UM),
        "array_axis_order": ["z", "y", "x"],
        "coordinate_axis_order": ["x", "y", "z"],
        "resampled_centerline_point_count": int(len(points)),
        "radius_scale": calibration.radius_scale,
        "voxel_floor": calibration.voxel_floor,
        "foreground_threshold_normalized": calibration.threshold,
        "foreground_threshold_uint8": PAPER_FOREGROUND_UINT8,
        "calibration_split": calibration.calibration_split,
        "calibration_note": (
            "Fixed on three train-split blocks; ROI003274 belongs to test split and was not "
            "used to select the defaults. This approximates, rather than reproduces, the "
            "unreleased annotation-map generator."
        ),
    }
    return reconstructed, metadata


def threshold_probability_mask(
    array: np.ndarray,
    threshold: int = PAPER_FOREGROUND_UINT8,
) -> np.ndarray:
    return np.asarray(array) > int(threshold)


def swc_edges(swc: SWCReference) -> Iterable[tuple[int, int]]:
    index_by_id = {int(node_id): index for index, node_id in enumerate(swc.node_ids)}
    for child, parent_id in enumerate(swc.parent_ids):
        parent = index_by_id.get(int(parent_id))
        if parent is not None:
            yield parent, child
