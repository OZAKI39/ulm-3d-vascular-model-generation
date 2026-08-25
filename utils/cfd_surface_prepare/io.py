"""Read and cross-check immutable PASS preprocessing inputs."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh


class SurfacePrepareError(RuntimeError):
    """A strict input, local surgery, or geometry QC failure."""


@dataclass(frozen=True, slots=True)
class BoundaryInput:
    index: int
    port_id: str
    boundary_origin: str
    role: str
    global_node_id: int | None
    global_edge_id: int
    center_um: np.ndarray
    source_radius_um: float
    pressure_original_pa: float
    expected_flow_m3_s: float
    simulation_tangent: np.ndarray
    outward_normal: np.ndarray
    extension_length_um: float
    extension_end_um: np.ndarray


@dataclass(frozen=True, slots=True)
class SurfaceInputs:
    preprocess_run: Path
    preprocess_summary: dict[str, Any]
    original_boundary_conditions: dict[str, Any]
    geometry_reference: dict[str, Any]
    original_surface_um_stl: Path
    original_surface_um_vtp: Path
    original_surface_m_stl: Path
    boundaries: tuple[BoundaryInput, ...]


def read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise SurfacePrepareError(f"Missing required input: {path}")
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SurfacePrepareError(f"Expected a JSON object: {path}")
    return value


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise SurfacePrepareError(f"Missing required input: {path}")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def _vector(row: dict[str, str], prefix: str) -> np.ndarray:
    return np.asarray([float(row[f"{prefix}_{axis}"]) for axis in "xyz"], dtype=float)


def load_surface_inputs(run_root: Path, *, expected_boundary_count: int) -> SurfaceInputs:
    """Load exactly one frozen PASS run without recalculating any upstream quantity."""

    root = Path(run_root).resolve()
    summary = read_json(root / "qc" / "run_summary.json")
    if summary.get("status") != "CFD_PREPROCESS_BASELINE_PASS":
        raise SurfacePrepareError("CFD_PREPROCESS_INPUT_NOT_PASS")
    original_bc = read_json(root / "roi" / "boundary_conditions.json")
    geometry = read_json(root / "input" / "geometry_reference.json")
    if geometry.get("status") != "PASS":
        raise SurfacePrepareError("CFD_PREPROCESS_GEOMETRY_REFERENCE_NOT_PASS")
    files = geometry.get("files")
    if not isinstance(files, dict):
        raise SurfacePrepareError("Invalid geometry_reference.json files mapping")
    surface_um_stl = Path(str(files.get("lumen_surface_um.stl", ""))).resolve()
    surface_um_vtp = Path(str(files.get("lumen_surface_um.vtp", ""))).resolve()
    surface_m_stl = Path(str(files.get("lumen_surface_m.stl", ""))).resolve()
    for path in (surface_um_stl, surface_um_vtp, surface_m_stl):
        if not path.is_file():
            raise SurfacePrepareError(f"Missing immutable geometry reference: {path}")
    expected_sha = geometry.get("sha256", {}).get("lumen_surface_um.stl")
    if expected_sha != sha256_file(surface_um_stl):
        raise SurfacePrepareError("ORIGINAL_ULTRALISER_SURFACE_HASH_MISMATCH")

    classification = _read_csv(root / "roi" / "port_classification.csv")
    extensions = _read_csv(root / "roi" / "port_extension_plan.csv")
    port_vtp = root / "roi" / "port_planes.vtp"
    if not port_vtp.is_file():
        raise SurfacePrepareError(f"Missing required input: {port_vtp}")
    if len(classification) != expected_boundary_count or len(extensions) != expected_boundary_count:
        raise SurfacePrepareError(
            f"Expected exactly {expected_boundary_count} CFD boundaries"
        )
    extension_by_id = {row["port_id"]: row for row in extensions}
    if len(extension_by_id) != expected_boundary_count:
        raise SurfacePrepareError("Duplicate boundary ID in extension plan")

    inlet = original_bc.get("inlet")
    outlets = original_bc.get("outlets")
    if not isinstance(inlet, dict) or not isinstance(outlets, list):
        raise SurfacePrepareError("Invalid boundary_conditions.json structure")
    bc_by_id = {str(inlet.get("port_id")): inlet}
    for outlet in outlets:
        if not isinstance(outlet, dict):
            raise SurfacePrepareError("Invalid outlet record")
        bc_by_id[str(outlet.get("port_id"))] = outlet
    ids = [row["port_id"] for row in classification]
    if set(ids) != set(extension_by_id) or set(ids) != set(bc_by_id):
        raise SurfacePrepareError("CFD boundary IDs disagree across saved PASS artifacts")

    boundaries: list[BoundaryInput] = []
    for index, row in enumerate(classification):
        port_id = row["port_id"]
        extension = extension_by_id[port_id]
        bc = bc_by_id[port_id]
        role = row["role"]
        origin = row["boundary_origin"]
        if role not in {"ASSUMED_INLET", "ASSUMED_OUTLET"}:
            raise SurfacePrepareError(f"Unsupported saved boundary role: {role}")
        if origin not in {"CUT_PORT", "TRUE_TERMINAL"}:
            raise SurfacePrepareError(f"Unsupported saved boundary origin: {origin}")
        normal = _vector(row, "outward_normal")
        tangent = _vector(row, "simulation_tangent")
        if not np.all(np.isfinite(normal)) or abs(np.linalg.norm(normal) - 1.0) > 1.0e-10:
            raise SurfacePrepareError(f"Invalid outward normal for {port_id}")
        center = np.asarray([float(row[f"{axis}_um"]) for axis in "xyz"], dtype=float)
        extension_end = np.asarray(
            [float(extension[f"extension_end_{axis}_um"]) for axis in "xyz"],
            dtype=float,
        )
        length = float(extension["extension_length_um"])
        planned_end = center + normal * length
        if np.linalg.norm(planned_end - extension_end) > 1.0e-8:
            raise SurfacePrepareError(f"Saved extension plan is inconsistent for {port_id}")
        if role == "ASSUMED_INLET":
            pressure = float(row["P_1D_pa"])
            flow_rate = float(bc["flow_rate_m3_s"])
            if bc.get("type") != "VOLUMETRIC_FLOW_RATE":
                raise SurfacePrepareError(f"Inlet BC mismatch for {port_id}")
        else:
            pressure = float(bc["pressure_pa"])
            flow_rate = float(bc["expected_1d_flow_m3_s"])
            if bc.get("type") != "PRESSURE_DIRICHLET":
                raise SurfacePrepareError(f"Outlet BC mismatch for {port_id}")
        global_node_text = row.get("global_node_id", "").strip()
        boundaries.append(
            BoundaryInput(
                index=index,
                port_id=port_id,
                boundary_origin=origin,
                role=role,
                global_node_id=int(global_node_text) if global_node_text else None,
                global_edge_id=int(row["global_edge_id"]),
                center_um=center,
                source_radius_um=float(row["radius_um"]),
                pressure_original_pa=pressure,
                expected_flow_m3_s=flow_rate,
                simulation_tangent=tangent,
                outward_normal=normal,
                extension_length_um=length,
                extension_end_um=extension_end,
            )
        )
    inlet_count = sum(item.role == "ASSUMED_INLET" for item in boundaries)
    outlet_count = sum(item.role == "ASSUMED_OUTLET" for item in boundaries)
    if inlet_count != 1 or outlet_count != 3:
        raise SurfacePrepareError(
            f"Expected 1 inlet and 3 outlets, got {inlet_count} and {outlet_count}"
        )
    if sum(item.boundary_origin == "TRUE_TERMINAL" for item in boundaries) != 1:
        raise SurfacePrepareError("Expected one saved TRUE_TERMINAL boundary")
    return SurfaceInputs(
        preprocess_run=root,
        preprocess_summary=summary,
        original_boundary_conditions=original_bc,
        geometry_reference=geometry,
        original_surface_um_stl=surface_um_stl,
        original_surface_um_vtp=surface_um_vtp,
        original_surface_m_stl=surface_m_stl,
        boundaries=tuple(boundaries),
    )


def load_original_surface(path: Path) -> trimesh.Trimesh:
    """Load the immutable VTP as a unique-vertex triangular working copy."""

    polydata = pv.read(path).triangulate()
    faces = np.asarray(polydata.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    mesh = trimesh.Trimesh(
        vertices=np.asarray(polydata.points, dtype=float).copy(),
        faces=faces.copy(),
        process=False,
    )
    if len(mesh.faces) == 0 or not mesh.is_watertight:
        raise SurfacePrepareError("Immutable input VTP is not a watertight triangle surface")
    return mesh
