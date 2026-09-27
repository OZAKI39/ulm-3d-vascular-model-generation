"""Fast, deterministic loading of NNE2 TIFF Z-stacks and microscope metadata."""

from __future__ import annotations

import re
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from xml.etree import ElementTree

import numpy as np
from PIL import Image
from scipy import ndimage


_FRAME_RE = re.compile(r"_(\d+)$")


@dataclass(frozen=True, slots=True)
class StackMetadata:
    stack_name: str
    stack_dir: Path
    frame_files: tuple[Path, ...]
    original_shape_zyx: tuple[int, int, int]
    processed_shape_zyx: tuple[int, int, int]
    original_spacing_xyz_um: tuple[float, float, float]
    processed_spacing_xyz_um: tuple[float, float, float]
    xml_file: Path

    def report(self) -> dict[str, object]:
        return {
            "stack_name": self.stack_name,
            "stack_dir": str(self.stack_dir),
            "frame_count": len(self.frame_files),
            "original_shape_zyx": self.original_shape_zyx,
            "processed_shape_zyx": self.processed_shape_zyx,
            "original_spacing_xyz_um": self.original_spacing_xyz_um,
            "processed_spacing_xyz_um": self.processed_spacing_xyz_um,
            "xml_file": str(self.xml_file),
        }


def _frame_number(path: Path) -> int:
    match = _FRAME_RE.search(path.stem)
    if not match:
        raise ValueError(f"Cannot read frame number from TIFF filename: {path.name}")
    return int(match.group(1))


def list_stack_frames(stack_dir: Path) -> tuple[Path, ...]:
    frames = tuple(sorted(stack_dir.glob("*Ch2*.tif"), key=_frame_number))
    if not frames:
        frames = tuple(sorted(stack_dir.glob("*Ch2*.tiff"), key=_frame_number))
    if not frames:
        raise FileNotFoundError(f"No Ch2 TIFF frames found in: {stack_dir}")
    numbers = [_frame_number(item) for item in frames]
    if len(numbers) != len(set(numbers)):
        raise ValueError(f"Duplicate TIFF frame numbers in: {stack_dir}")
    if numbers != list(range(numbers[0], numbers[0] + len(numbers))):
        raise ValueError(f"Missing TIFF frame inside stack: {stack_dir}")
    return frames


def _read_frame(path: Path) -> np.ndarray:
    with Image.open(path) as image:
        array = np.asarray(image)
    if array.ndim != 2:
        raise ValueError(f"Expected one grayscale plane: {path}")
    return np.asarray(array, dtype=np.uint16)


def _xml_spacing(xml_file: Path) -> tuple[float | None, float | None, float | None]:
    try:
        root = ElementTree.parse(xml_file).getroot()
    except ElementTree.ParseError:
        return None, None, None
    values: dict[str, float] = {}
    for element in root.iter("Key"):
        key = element.attrib.get("key", "")
        if key not in {"micronsPerPixel_XAxis", "micronsPerPixel_YAxis", "positionCurrent_ZAxis"}:
            continue
        raw = element.attrib.get("value")
        try:
            value = float(raw) if raw is not None else None
        except ValueError:
            value = None
        if value is not None and key not in values:
            values[key] = value
    return (
        values.get("micronsPerPixel_XAxis"),
        values.get("micronsPerPixel_YAxis"),
        None,
    )


def inspect_stack(
    stack_dir: Path,
    *,
    xy_spacing_um: float,
    z_spacing_um: float,
    target_xy_spacing_um: float,
) -> StackMetadata:
    frames = list_stack_frames(stack_dir)
    with Image.open(frames[0]) as image:
        width, height = image.size
    xml_files = tuple(stack_dir.glob("*.xml"))
    if len(xml_files) != 1:
        raise ValueError(f"Expected exactly one XML metadata file in: {stack_dir}")
    xml_x, xml_y, _ = _xml_spacing(xml_files[0])
    if xml_x is not None and abs(xml_x - xy_spacing_um) > 0.02:
        raise ValueError(
            f"vdb/XML X spacing conflict for {stack_dir.name}: {xy_spacing_um} vs {xml_x}"
        )
    if xml_y is not None and abs(xml_y - xy_spacing_um) > 0.02:
        raise ValueError(
            f"vdb/XML Y spacing conflict for {stack_dir.name}: {xy_spacing_um} vs {xml_y}"
        )
    out_x = max(1, int(round(width * xy_spacing_um / target_xy_spacing_um)))
    out_y = max(1, int(round(height * xy_spacing_um / target_xy_spacing_um)))
    actual_x = width * xy_spacing_um / out_x
    actual_y = height * xy_spacing_um / out_y
    return StackMetadata(
        stack_name=stack_dir.name,
        stack_dir=stack_dir,
        frame_files=frames,
        original_shape_zyx=(len(frames), height, width),
        processed_shape_zyx=(len(frames), out_y, out_x),
        original_spacing_xyz_um=(xy_spacing_um, xy_spacing_um, z_spacing_um),
        processed_spacing_xyz_um=(actual_x, actual_y, z_spacing_um),
        xml_file=xml_files[0],
    )


def load_stack(metadata: StackMetadata, *, workers: int = 1) -> np.ndarray:
    with ThreadPoolExecutor(max_workers=workers) as pool:
        planes = list(pool.map(_read_frame, metadata.frame_files))
    shape = planes[0].shape
    if any(item.shape != shape for item in planes):
        raise ValueError(f"TIFF frame dimensions differ inside: {metadata.stack_dir}")
    volume = np.stack(planes, axis=0)
    target = metadata.processed_shape_zyx
    if volume.shape != target:
        zoom = (1.0, target[1] / volume.shape[1], target[2] / volume.shape[2])
        volume = ndimage.zoom(volume, zoom=zoom, order=1, prefilter=False)
    return np.asarray(volume, dtype=np.float32)


def load_original_frame(metadata: StackMetadata, frame_index_one_based: int) -> np.ndarray:
    if frame_index_one_based < 1 or frame_index_one_based > len(metadata.frame_files):
        raise IndexError(
            f"Frame {frame_index_one_based} outside 1..{len(metadata.frame_files)} for "
            f"{metadata.stack_name}"
        )
    return _read_frame(metadata.frame_files[frame_index_one_based - 1])
