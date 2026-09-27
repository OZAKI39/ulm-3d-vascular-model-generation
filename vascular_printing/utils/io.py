"""File readers and writers used by the preprocessing pipeline."""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping

import nibabel as nib
import numpy as np
import vtk


def read_stl(path: Path) -> vtk.vtkPolyData:
    reader = vtk.vtkSTLReader()
    reader.SetFileName(str(path))
    reader.MergingOn()
    reader.Update()
    mesh = vtk.vtkPolyData()
    mesh.DeepCopy(reader.GetOutput())
    if mesh.GetNumberOfPoints() == 0 or mesh.GetNumberOfCells() == 0:
        raise ValueError(f"STL contains no usable surface: {path}")
    return mesh


def write_binary_stl(mesh: vtk.vtkPolyData, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = vtk.vtkSTLWriter()
    writer.SetFileName(str(path))
    writer.SetFileTypeToBinary()
    writer.SetInputData(mesh)
    if writer.Write() != 1:
        raise OSError(f"Failed to write STL: {path}")


def lps_to_ras_affine(origin_lps: tuple[float, float, float], spacing: tuple[float, float, float]) -> np.ndarray:
    """Return a NIfTI affine while preserving the source mesh's LPS geometry."""

    ox, oy, oz = origin_lps
    sx, sy, sz = spacing
    return np.asarray(
        [
            [-sx, 0.0, 0.0, -ox],
            [0.0, -sy, 0.0, -oy],
            [0.0, 0.0, sz, oz],
            [0.0, 0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )


def write_nifti_mask(
    data_xyz: np.ndarray,
    path: Path,
    origin_lps: tuple[float, float, float],
    spacing: tuple[float, float, float],
    description: str,
) -> np.ndarray:
    path.parent.mkdir(parents=True, exist_ok=True)
    affine = lps_to_ras_affine(origin_lps, spacing)
    image = nib.Nifti1Image(np.asarray(data_xyz, dtype=np.uint8), affine)
    image.set_qform(affine, code=1)
    image.set_sform(affine, code=1)
    image.header.set_xyzt_units("micron")
    image.header["descrip"] = description[:79]
    nib.save(image, str(path))
    return affine


def _json_default(value: Any) -> Any:
    if is_dataclass(value):
        return asdict(value)
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(f"Object is not JSON serializable: {type(value).__name__}")


def write_json(payload: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, default=_json_default) + "\n",
        encoding="utf-8",
    )


def write_csv(rows: Iterable[Mapping[str, Any]], path: Path) -> None:
    materialized = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    if not materialized:
        path.write_text("", encoding="utf-8-sig")
        return
    with path.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(materialized[0].keys()))
        writer.writeheader()
        writer.writerows(materialized)

