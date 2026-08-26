"""Output layout and file writers for the CFD-only derived surface."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pyvista as pv
import trimesh

from .io import BoundaryInput, SurfacePrepareError
from .types import TaggedSurface


@dataclass(frozen=True, slots=True)
class OutputLayout:
    root: Path
    input: Path
    geometry: Path
    boundaries: Path
    bc: Path
    qc: Path
    figures: Path
    report: Path


def create_layout(output_root: Path, run_id: str) -> OutputLayout:
    root = Path(output_root) / run_id
    if root.exists():
        raise SurfacePrepareError(f"Output run already exists: {root}")
    directories = {
        name: root / name
        for name in ("input", "geometry", "boundaries", "bc", "qc", "figures", "report")
    }
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=False)
    return OutputLayout(root=root, **directories)


def _json_ready(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    return value


def write_json(path: Path, payload: Any) -> Path:
    path.write_text(
        json.dumps(_json_ready(payload), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    return path


def write_csv(path: Path, rows: Iterable[dict[str, Any]]) -> Path:
    records = list(rows)
    if not records:
        raise SurfacePrepareError(f"Cannot write empty CSV: {path}")
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    return path


def _submesh(surface: TaggedSurface, mask: np.ndarray) -> trimesh.Trimesh:
    mesh = trimesh.Trimesh(
        vertices=surface.vertices.copy(),
        faces=surface.faces[np.asarray(mask, dtype=bool)].copy(),
        process=False,
    )
    mesh.remove_unreferenced_vertices()
    return mesh


def _export_stl(mesh: trimesh.Trimesh, path: Path) -> Path:
    mesh.export(path, file_type="stl")
    if not path.is_file() or path.stat().st_size == 0:
        raise SurfacePrepareError(f"Failed to export STL: {path}")
    return path


def export_geometry(
    surface: TaggedSurface,
    boundaries: Iterable[BoundaryInput],
    layout: OutputLayout,
    *,
    create_meter_copy: bool,
) -> dict[str, Any]:
    """Write combined/wall/cap STL files and the authoritative tagged VTP."""

    boundary_list = list(boundaries)
    full_mesh = surface.mesh()
    full_um = _export_stl(
        full_mesh, layout.geometry / "cfd_surface_refined_um.stl"
    )
    faces = np.column_stack(
        (np.full(len(surface.faces), 3, dtype=np.int64), surface.faces)
    ).ravel()
    polydata = pv.PolyData(surface.vertices, faces)
    polydata.cell_data["boundary_type_code"] = surface.boundary_type
    polydata.cell_data["boundary_index"] = surface.boundary_index
    polydata.cell_data["boundary_origin_code"] = surface.boundary_origin
    polydata.cell_data["extension_boundary_index"] = surface.extension_index
    polydata.cell_data["extension_band"] = surface.extension_band
    port_ids = np.full(len(surface.faces), "", dtype=f"<U{max(map(len, [item.port_id for item in boundary_list]))}")
    for boundary in boundary_list:
        port_ids[surface.boundary_index == boundary.index] = boundary.port_id
    polydata.cell_data["port_id"] = port_ids
    full_vtp = layout.geometry / "cfd_surface_refined_um.vtp"
    polydata.save(full_vtp, binary=True)

    wall_mask = surface.boundary_type == 0
    wall_path = _export_stl(
        _submesh(surface, wall_mask), layout.geometry / "cfd_wall_refined_um.stl"
    )
    meter_path: Path | None = None
    if create_meter_copy:
        meter_mesh = trimesh.Trimesh(
            vertices=surface.vertices * 1.0e-6,
            faces=surface.faces.copy(),
            process=False,
        )
        meter_path = _export_stl(
            meter_mesh, layout.geometry / "cfd_surface_refined_m.stl"
        )

    manifest: list[dict[str, Any]] = []
    boundary_paths: list[Path] = []
    for boundary in boundary_list:
        role = "inlet" if boundary.role == "ASSUMED_INLET" else "outlet"
        origin = boundary.boundary_origin.lower()
        name = f"boundary_{boundary.index:02d}_{role}_{origin}_um.stl"
        mask = (surface.boundary_type > 0) & (
            surface.boundary_index == boundary.index
        )
        path = _export_stl(_submesh(surface, mask), layout.boundaries / name)
        boundary_paths.append(path)
        manifest.append(
            {
                "boundary_index": boundary.index,
                "port_id": boundary.port_id,
                "role": boundary.role,
                "boundary_origin": boundary.boundary_origin,
                "boundary_type_code": 1 if role == "inlet" else 2,
                "boundary_origin_code": 1
                if boundary.boundary_origin == "CUT_PORT"
                else 2,
                "triangle_count": int(np.count_nonzero(mask)),
                "path": str(path.resolve()),
            }
        )
    manifest_path = write_csv(
        layout.boundaries / "boundary_manifest.csv", manifest
    )
    return {
        "cfd_surface_refined_um_stl": full_um.resolve(),
        "cfd_surface_refined_um_vtp": full_vtp.resolve(),
        "cfd_surface_refined_m_stl": meter_path.resolve() if meter_path else None,
        "cfd_wall_refined_um_stl": wall_path.resolve(),
        "boundary_stl_directory": layout.boundaries.resolve(),
        "boundary_stl_paths": [path.resolve() for path in boundary_paths],
        "boundary_manifest_csv": manifest_path.resolve(),
    }


def meter_scale_qc(um_path: Path, meter_path: Path) -> dict[str, Any]:
    um_mesh = trimesh.load_mesh(um_path, process=False)
    meter_mesh = trimesh.load_mesh(meter_path, process=False)
    um_extent = np.asarray(um_mesh.bounds[1] - um_mesh.bounds[0], dtype=float)
    meter_extent = np.asarray(
        meter_mesh.bounds[1] - meter_mesh.bounds[0], dtype=float
    )
    valid = bool(np.allclose(meter_extent, um_extent * 1.0e-6, rtol=2.0e-6, atol=1.0e-12))
    return {
        "status": "PASS" if valid else "FAIL",
        "scale_factor": 1.0e-6,
        "um_extent": um_extent.tolist(),
        "meter_extent": meter_extent.tolist(),
        "check": valid,
    }
