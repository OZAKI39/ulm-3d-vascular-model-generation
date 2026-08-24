"""Restricted pickle loading and schema normalization for Schmid dictionaries."""

from __future__ import annotations

import hashlib
import _codecs
import pickle
from pathlib import Path
from typing import Any

import numpy as np
from numpy._core import multiarray as numpy_multiarray

from .model import SchmidInputData


VERTEX_KEYS = {"pressure", "coords", "pBC"}
EDGE_KEYS = {
    "diameter",
    "tuple",
    "flow",
    "httBC",
    "nkind",
    "length",
    "htt",
    "nRBC",
    "diameters",
    "points",
}


class RestrictedNumpyUnpickler(pickle.Unpickler):
    """Allow only the NumPy constructors used by the published dictionaries."""

    def find_class(self, module: str, name: str) -> Any:
        allowed = {
            ("_codecs", "encode"): _codecs.encode,
            ("numpy", "ndarray"): np.ndarray,
            ("numpy", "dtype"): np.dtype,
            ("numpy.core.multiarray", "_reconstruct"): numpy_multiarray._reconstruct,
            ("numpy.core.multiarray", "scalar"): numpy_multiarray.scalar,
            ("numpy._core.multiarray", "_reconstruct"): numpy_multiarray._reconstruct,
            ("numpy._core.multiarray", "scalar"): numpy_multiarray.scalar,
        }
        target = allowed.get((module, name))
        if target is None:
            raise pickle.UnpicklingError(f"Blocked global in external pickle: {module}.{name}")
        return target


def _load_restricted(path: Path) -> dict[str, Any]:
    with path.open("rb") as stream:
        payload = RestrictedNumpyUnpickler(
            stream, fix_imports=True, encoding="latin1", errors="strict"
        ).load()
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a dictionary in {path.name}, got {type(payload).__name__}")
    return payload


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _optional_float_array(values: list[Any], expected: int, field: str) -> np.ndarray:
    if len(values) != expected:
        raise ValueError(f"{field} has {len(values)} values; expected {expected}")
    output = np.full(expected, np.nan, dtype=np.float64)
    for index, value in enumerate(values):
        if value is not None:
            output[index] = float(value)
    return output


def load_schmid_input(input_dir: Path) -> SchmidInputData:
    vertex_path = input_dir / "verticesDict.pkl"
    edge_path = input_dir / "edgesDict.pkl"
    vertices = _load_restricted(vertex_path)
    edges = _load_restricted(edge_path)
    missing_vertices = VERTEX_KEYS - set(vertices)
    missing_edges = EDGE_KEYS - set(edges)
    if missing_vertices:
        raise ValueError(f"verticesDict.pkl is missing keys: {sorted(missing_vertices)}")
    if missing_edges:
        raise ValueError(f"edgesDict.pkl is missing keys: {sorted(missing_edges)}")

    coordinates = np.asarray(vertices["coords"], dtype=np.float64)
    pressure = np.asarray(vertices["pressure"], dtype=np.float64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(f"coords must have shape (N, 3), got {coordinates.shape}")
    if pressure.shape != (len(coordinates),):
        raise ValueError("pressure length does not match coords")
    vertex_count = len(coordinates)
    pressure_boundary = _optional_float_array(vertices["pBC"], vertex_count, "pBC")

    edge_tuples = np.asarray(edges["tuple"], dtype=np.int64)
    if edge_tuples.ndim != 2 or edge_tuples.shape[1] != 2:
        raise ValueError(f"tuple must have shape (M, 2), got {edge_tuples.shape}")
    edge_count = len(edge_tuples)

    def numeric(name: str, dtype: Any = np.float64) -> np.ndarray:
        array = np.asarray(edges[name], dtype=dtype)
        if array.shape != (edge_count,):
            raise ValueError(f"{name} has shape {array.shape}; expected ({edge_count},)")
        return array

    diameter_profiles: list[np.ndarray] = []
    point_sequences: list[np.ndarray] = []
    if len(edges["diameters"]) != edge_count or len(edges["points"]) != edge_count:
        raise ValueError("diameters/points list length does not match tuple")
    for edge_id, (diameters, points) in enumerate(zip(edges["diameters"], edges["points"], strict=True)):
        diameter_array = np.asarray(diameters, dtype=np.float64).reshape(-1)
        point_array = np.asarray(points, dtype=np.float64)
        if point_array.ndim != 2 or point_array.shape[1] != 3:
            raise ValueError(f"points[{edge_id}] must have shape (K, 3), got {point_array.shape}")
        diameter_profiles.append(diameter_array)
        point_sequences.append(point_array)

    source_files: dict[str, dict[str, Any]] = {}
    for path in (vertex_path, edge_path):
        source_files[path.name] = {
            "path": str(path.resolve()),
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        }
    rbc_path = input_dir / "RBC_trajectories.pkl"
    if rbc_path.is_file():
        source_files[rbc_path.name] = {
            "path": str(rbc_path.resolve()),
            "size_bytes": rbc_path.stat().st_size,
            "sha256": None,
            "loaded": False,
            "reason": "Not required for Step 1-3; intentionally not loaded or hashed.",
        }

    return SchmidInputData(
        coordinates_um=coordinates,
        pressure_mmhg=pressure,
        pressure_boundary_mmhg=pressure_boundary,
        edge_tuples=edge_tuples,
        mean_diameter_um=numeric("diameter"),
        flow_um3_per_ms=numeric("flow"),
        hematocrit_boundary=_optional_float_array(edges["httBC"], edge_count, "httBC"),
        vessel_type_code=numeric("nkind", np.int64),
        source_length_um=numeric("length"),
        hematocrit=numeric("htt"),
        red_blood_cell_count=numeric("nRBC"),
        diameter_profiles_um=diameter_profiles,
        point_sequences_um=point_sequences,
        source_files=source_files,
    )
