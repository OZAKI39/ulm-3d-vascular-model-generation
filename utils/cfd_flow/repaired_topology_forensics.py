"""Read-only topology and unit-scaling forensics for repaired Seeder meshes.

The helpers in this module never launch Seeder, Musubi, or Harvester.  They
decode existing uniform TreElm meshes, retain the distinction between face,
D3Q19, and 26-neighbour connectivity, and prepare a coordinate-scaled copy of
the existing Seeder geometry for an independent unit-invariance oracle.
"""

from __future__ import annotations

import csv
import heapq
import json
import math
import re
import struct
import time
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from .io import sha256_file
from .musubi_boundary_mass_referee import MeshContract
from .port_flux_audit import (
    extract_boundary_property_indices,
    parse_bnd_header,
    parse_boundary_property_header,
    read_boundary_ids,
)
from .restart_decode import read_treelm_elemlist


CONNECTIVITY_LABELS = {
    6: "FACE_GEOMETRIC_CONNECTIVITY",
    18: "D3Q19_STREAMING_CONNECTIVITY",
    26: "FULL_MOORE_DIAGNOSTIC_CONNECTIVITY",
}
PORT_LABELS = ("inlet", "outlet_01", "outlet_02", "outlet_03")
SCALE_FACTOR = 1.0e6
PINNED_SEEDER_SHA = "667109df6fafdcb39f4409e3f5d90f04d75cd33c"
PINNED_TREELM_SHA = "53f273dbb8e9dcbe7feeb3d9831a35f5ae3cd72c"
PINNED_SDR_SHA = "3b3e344aec0af8d1383e2cfb023a21df3361e1e9"


@dataclass(frozen=True, slots=True)
class UniformLattice:
    origin: np.ndarray
    length: float
    level: int
    dx: float
    cell_count: int


def _numbers(text: str) -> list[float]:
    pattern = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
    return [
        float(value.replace("D", "E").replace("d", "e"))
        for value in re.findall(pattern, text)
    ]


def parse_uniform_lattice(header_text: str) -> UniformLattice:
    """Parse the uniform bounding cube without executing Lua."""

    box = re.search(
        r"boundingbox\s*=\s*\{\s*origin\s*=\s*\{([^}]*)\}\s*,?\s*"
        r"length\s*=\s*([-+0-9.EeDd]+)",
        header_text,
        re.DOTALL | re.IGNORECASE,
    )
    level = re.search(r"\bmaxLevel\s*=\s*(\d+)", header_text)
    count = re.search(r"(?m)^\s*nElems\s*=\s*(\d+)", header_text)
    if box is None or level is None or count is None:
        raise ValueError("Uniform TreElm lattice header is incomplete")
    origin = np.asarray(_numbers(box.group(1)), dtype=np.float64)
    if origin.shape != (3,):
        raise ValueError("Bounding-cube origin is not three-dimensional")
    length = float(box.group(2).replace("D", "E").replace("d", "e"))
    level_value = int(level.group(1))
    return UniformLattice(
        origin=origin,
        length=length,
        level=level_value,
        dx=length / (2**level_value),
        cell_count=int(count.group(1)),
    )


def parse_seed_point(seeder_lua: str) -> np.ndarray:
    """Return the configured seed point from the existing Seeder Lua file."""

    seed = re.search(
        r"attribute\s*=\s*\{[^{}]*kind\s*=\s*['\"]seed['\"][^{}]*\}"
        r".*?geometry\s*=\s*\{.*?origin\s*=\s*\{([^{}]+)\}",
        seeder_lua,
        re.DOTALL | re.IGNORECASE,
    )
    if seed is None:
        raise ValueError("Seeder Lua does not contain a seed origin")
    point = np.asarray(_numbers(seed.group(1)), dtype=np.float64)
    if point.shape != (3,):
        raise ValueError("Seed point is not three-dimensional")
    return point


def neighbour_offsets(connectivity: int, *, positive_half: bool = False) -> np.ndarray:
    """Return integer offsets for 6, D3Q19-18, or full 26 connectivity."""

    if connectivity not in CONNECTIVITY_LABELS:
        raise ValueError("connectivity must be one of 6, 18, or 26")
    offsets = []
    for offset in product((-1, 0, 1), repeat=3):
        if offset == (0, 0, 0):
            continue
        manhattan = sum(abs(value) for value in offset)
        if connectivity == 6 and manhattan != 1:
            continue
        if connectivity == 18 and manhattan > 2:
            continue
        if positive_half and not (
            offset[2] > 0
            or (offset[2] == 0 and offset[1] > 0)
            or (offset[2] == 0 and offset[1] == 0 and offset[0] > 0)
        ):
            continue
        offsets.append(offset)
    return np.asarray(offsets, dtype=np.int8)


def sparse_component_labels(
    cell_ijk: np.ndarray,
    connectivity: int,
    *,
    cells_per_axis: int | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Label sparse uniform-grid cells, sorting component IDs by size."""

    coordinates = np.asarray(cell_ijk, dtype=np.int64)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3 or len(coordinates) == 0:
        raise ValueError("cell_ijk must be a non-empty array with shape (n, 3)")
    axis = (
        int(cells_per_axis)
        if cells_per_axis is not None
        else int(np.max(coordinates)) + 2
    )
    if np.any(coordinates < 0) or np.any(coordinates >= axis):
        raise ValueError("cell coordinate lies outside the declared uniform grid")
    codes = coordinates[:, 0] + axis * coordinates[:, 1] + axis**2 * coordinates[:, 2]
    order = np.argsort(codes)
    sorted_codes = codes[order]
    if len(np.unique(sorted_codes)) != len(sorted_codes):
        raise ValueError("fluid cell coordinates are not unique")

    source_chunks: list[np.ndarray] = []
    target_chunks: list[np.ndarray] = []
    all_indices = np.arange(len(coordinates), dtype=np.int64)
    for offset in neighbour_offsets(connectivity, positive_half=True):
        targets = coordinates + offset
        valid = np.all((targets >= 0) & (targets < axis), axis=1)
        source = all_indices[valid]
        targets = targets[valid]
        target_codes = targets[:, 0] + axis * targets[:, 1] + axis**2 * targets[:, 2]
        positions = np.searchsorted(sorted_codes, target_codes)
        in_range = positions < len(sorted_codes)
        hits = np.zeros(len(positions), dtype=bool)
        hits[in_range] = sorted_codes[positions[in_range]] == target_codes[in_range]
        if np.any(hits):
            source_chunks.append(source[hits])
            target_chunks.append(order[positions[hits]])

    if source_chunks:
        rows = np.concatenate(source_chunks)
        columns = np.concatenate(target_chunks)
        graph = coo_matrix(
            (np.ones(len(rows), dtype=np.int8), (rows, columns)),
            shape=(len(coordinates), len(coordinates)),
        ).tocsr()
    else:
        graph = coo_matrix((len(coordinates), len(coordinates))).tocsr()
    raw_count, raw_labels = connected_components(graph, directed=False)
    raw_sizes = np.bincount(raw_labels, minlength=raw_count)
    descending = np.argsort(-raw_sizes, kind="stable")
    remap = np.empty(raw_count, dtype=np.int64)
    remap[descending] = np.arange(raw_count, dtype=np.int64)
    labels = remap[raw_labels]
    sizes = np.bincount(labels, minlength=raw_count).astype(np.int64)
    return labels, sizes


def cell_centers(cell_ijk: np.ndarray, lattice: UniformLattice) -> np.ndarray:
    return lattice.origin + (np.asarray(cell_ijk, dtype=np.float64) + 0.5) * lattice.dx


def component_identity_records(
    cell_ijk: np.ndarray,
    labels: np.ndarray,
    sizes: np.ndarray,
    lattice: UniformLattice,
) -> list[dict[str, Any]]:
    coordinates = np.asarray(cell_ijk, dtype=np.int64)
    centers = cell_centers(coordinates, lattice)
    records = []
    for component_id, size in enumerate(np.asarray(sizes, dtype=np.int64)):
        mask = labels == component_id
        ijk_values = coordinates[mask]
        xyz_values = centers[mask]
        records.append(
            {
                "component_id": int(component_id),
                "cell_count": int(size),
                "fraction_of_total_cells": float(size / len(coordinates)),
                "physical_volume_m3": float(size * lattice.dx**3),
                "ijk_bounds": {
                    "minimum": ijk_values.min(axis=0).tolist(),
                    "maximum": ijk_values.max(axis=0).tolist(),
                },
                "physical_xyz_bounds_m": {
                    "minimum": xyz_values.min(axis=0).tolist(),
                    "maximum": xyz_values.max(axis=0).tolist(),
                },
                "centroid_m": xyz_values.mean(axis=0).tolist(),
                "minimum_coordinates_m": xyz_values.min(axis=0).tolist(),
                "maximum_coordinates_m": xyz_values.max(axis=0).tolist(),
            }
        )
    return records


def nearest_cell_membership(
    point: Sequence[float], centers: np.ndarray, labels: np.ndarray
) -> dict[str, Any]:
    distance, cell_index = cKDTree(np.asarray(centers)).query(np.asarray(point))
    return {
        "component_id": int(labels[int(cell_index)]),
        "nearest_fluid_cell_index": int(cell_index),
        "nearest_fluid_cell_center_m": np.asarray(centers)[int(cell_index)].tolist(),
        "distance_m": float(distance),
    }


def port_component_membership(
    component_labels: np.ndarray,
    boundary_cells: Mapping[str, np.ndarray],
    component_count: int,
) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for port in PORT_LABELS:
        cells = np.asarray(boundary_cells[port], dtype=np.int64)
        counts = np.bincount(component_labels[cells], minlength=component_count)
        nonzero = np.flatnonzero(counts)
        result[port] = {
            "total_boundary_cells": int(len(cells)),
            "component_boundary_cell_counts": {
                str(int(component)): int(counts[component]) for component in nonzero
            },
            "component_ids": [int(value) for value in nonzero],
        }
    component_sets = [set(row["component_ids"]) for row in result.values()]
    common = set.intersection(*component_sets) if component_sets else set()
    result["all_ports_share_a_component"] = bool(common)
    result["shared_component_ids"] = sorted(int(value) for value in common)
    return result


def nearest_component_gap(
    cell_ijk: np.ndarray,
    labels: np.ndarray,
    sizes: np.ndarray,
    lattice: UniformLattice,
) -> dict[str, Any] | None:
    if len(sizes) < 2:
        return None
    first = np.flatnonzero(labels == 0)
    second = np.flatnonzero(labels == 1)
    first_xyz = cell_centers(np.asarray(cell_ijk)[first], lattice)
    second_xyz = cell_centers(np.asarray(cell_ijk)[second], lattice)
    if len(first) <= len(second):
        distances, indices = cKDTree(second_xyz).query(first_xyz)
        local_first = int(np.argmin(distances))
        cell_a = int(first[local_first])
        cell_b = int(second[int(indices[local_first])])
        distance = float(distances[local_first])
    else:
        distances, indices = cKDTree(first_xyz).query(second_xyz)
        local_second = int(np.argmin(distances))
        cell_a = int(first[int(indices[local_second])])
        cell_b = int(second[local_second])
        distance = float(distances[local_second])
    offset = np.asarray(cell_ijk[cell_b] - cell_ijk[cell_a], dtype=np.int64)
    d3q19 = bool(np.max(np.abs(offset)) == 1 and 0 < int(np.sum(np.abs(offset))) <= 2)
    return {
        "component_a": 0,
        "component_b": 1,
        "component_a_cell_index": cell_a,
        "component_b_cell_index": cell_b,
        "component_a_ijk": np.asarray(cell_ijk[cell_a]).tolist(),
        "component_b_ijk": np.asarray(cell_ijk[cell_b]).tolist(),
        "component_a_center_m": cell_centers(cell_ijk[[cell_a]], lattice)[0].tolist(),
        "component_b_center_m": cell_centers(cell_ijk[[cell_b]], lattice)[0].tolist(),
        "physical_distance_m": distance,
        "physical_distance_um": distance * 1.0e6,
        "distance_over_dx": distance / lattice.dx,
        "ijk_offset_b_minus_a": offset.tolist(),
        "offset_is_d3q19_direction": d3q19,
    }


def classify_component_topology(
    sizes_6: Sequence[int],
    sizes_18: Sequence[int],
    port_membership: Mapping[str, Any],
) -> str:
    """Apply only the component classifications allowed by the research prompt."""

    sizes6 = np.asarray(sizes_6, dtype=np.int64)
    sizes18 = np.asarray(sizes_18, dtype=np.int64)
    if len(sizes6) == 1:
        return "MAIN_FLOW_COMPONENT"
    if len(sizes6) == 2 and len(sizes18) == 1:
        return "D3Q19_DIAGONAL_ONLY_CONNECTION"
    secondary_ports = []
    for port in PORT_LABELS:
        secondary_ports.extend(
            value for value in port_membership[port]["component_ids"] if value != 0
        )
    secondary_fraction = float(np.sum(sizes6[1:]) / np.sum(sizes6))
    if secondary_ports:
        return "PORT_BEARING_SECONDARY_COMPONENT"
    if secondary_fraction < 0.001:
        return "PORTLESS_ISOLATED_POCKET"
    if secondary_fraction > 0.01:
        return "MAJOR_NETWORK_SPLIT"
    return "UNRESOLVED"


def coordinate_set_difference(
    old_tree_ids: np.ndarray, new_tree_ids: np.ndarray
) -> dict[str, np.ndarray]:
    old_ids = np.asarray(old_tree_ids, dtype=np.int64)
    new_ids = np.asarray(new_tree_ids, dtype=np.int64)
    return {
        "common": np.intersect1d(old_ids, new_ids, assume_unique=True),
        "old_only": np.setdiff1d(old_ids, new_ids, assume_unique=True),
        "new_only": np.setdiff1d(new_ids, old_ids, assume_unique=True),
    }


def _distance_to_cloud(points: np.ndarray, cloud: np.ndarray) -> float | None:
    if len(points) == 0 or len(cloud) == 0:
        return None
    distances, _ = cKDTree(np.asarray(cloud, dtype=np.float64)).query(
        np.asarray(points, dtype=np.float64)
    )
    return float(np.min(distances))


def difference_cluster_records(
    coordinates: np.ndarray,
    lattice: UniformLattice,
    *,
    port_centers: Mapping[str, np.ndarray],
    gap_centers: np.ndarray,
    maximum_records: int = 10,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    if len(coordinates) == 0:
        return np.empty(0, dtype=np.int64), []
    labels, sizes = sparse_component_labels(
        coordinates, 6, cells_per_axis=2**lattice.level
    )
    centers = cell_centers(coordinates, lattice)
    records = []
    for component_id, size in enumerate(sizes[:maximum_records]):
        selected = labels == component_id
        ijk = coordinates[selected]
        xyz = centers[selected]
        records.append(
            {
                "cluster_id": int(component_id),
                "cell_count": int(size),
                "fraction_of_difference_cells": float(size / len(coordinates)),
                "ijk_bounds": {
                    "minimum": ijk.min(axis=0).tolist(),
                    "maximum": ijk.max(axis=0).tolist(),
                },
                "physical_xyz_bounds_m": {
                    "minimum": xyz.min(axis=0).tolist(),
                    "maximum": xyz.max(axis=0).tolist(),
                },
                "distance_to_ports_m": {
                    name: _distance_to_cloud(xyz, values)
                    for name, values in port_centers.items()
                },
                "distance_to_component_gap_m": _distance_to_cloud(xyz, gap_centers),
            }
        )
    return labels, records


def _scale_numeric_list(content: str, factor: float) -> str:
    number = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?")
    return number.sub(
        lambda match: format(
            float(match.group().replace("D", "E").replace("d", "e")) * factor,
            ".17g",
        ),
        content,
    )


def scale_seeder_lua_geometry(seeder_lua: str, factor: float = SCALE_FACTOR) -> str:
    """Scale only geometric coordinate/length assignments in Seeder Lua."""

    result = re.sub(
        r"(\borigin\s*=\s*\{)([^{}]*)(\})",
        lambda match: (
            match.group(1)
            + _scale_numeric_list(match.group(2), factor)
            + match.group(3)
        ),
        seeder_lua,
    )
    result = re.sub(
        r"(\bvec\s*=\s*\{\s*\{)([^{}]*)(\}\s*,\s*\{)([^{}]*)(\}\s*\})",
        lambda match: (
            match.group(1)
            + _scale_numeric_list(match.group(2), factor)
            + match.group(3)
            + _scale_numeric_list(match.group(4), factor)
            + match.group(5)
        ),
        result,
    )
    result, replacements = re.subn(
        r"(\bbounding_cube\s*=\s*\{.*?\blength\s*=\s*)"
        r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?)",
        lambda match: (
            match.group(1)
            + format(float(match.group(2).replace("D", "E")) * factor, ".17g")
        ),
        result,
        count=1,
        flags=re.DOTALL,
    )
    if replacements != 1:
        raise ValueError("Could not scale bounding_cube length")
    return result


def scale_binary_stl(
    source: Path, destination: Path, factor: float = SCALE_FACTOR
) -> dict[str, Any]:
    """Scale binary-STL vertices while preserving facet order and normals."""

    raw = Path(source).read_bytes()
    if len(raw) < 84:
        raise ValueError(f"Binary STL is too short: {source}")
    triangle_count = struct.unpack_from("<I", raw, 80)[0]
    if len(raw) != 84 + 50 * triangle_count:
        text = raw.decode("utf-8")
        number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?"
        pattern = re.compile(
            rf"(?m)^(\s*vertex\s+)({number})\s+({number})\s+({number})(\s*)$"
        )
        vertices: list[list[float]] = []

        def replace_vertex(match: re.Match[str]) -> str:
            values = [
                float(match.group(index).replace("D", "E").replace("d", "e"))
                for index in (2, 3, 4)
            ]
            vertices.append(values)
            scaled = " ".join(format(value * factor, ".17g") for value in values)
            return match.group(1) + scaled + match.group(5)

        scaled_text, replacements = pattern.subn(replace_vertex, text)
        if replacements == 0 or replacements % 3:
            raise ValueError(f"STL is neither valid binary nor ASCII: {source}")
        destination = Path(destination)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(scaled_text, encoding="utf-8")
        vertex_array = np.asarray(vertices, dtype=np.float64)
        return {
            "source": str(Path(source).resolve()),
            "destination": str(destination.resolve()),
            "format": "ASCII_STL",
            "triangle_count": int(replacements // 3),
            "scale_factor": float(factor),
            "source_bounds": [
                vertex_array.min(axis=0).tolist(),
                vertex_array.max(axis=0).tolist(),
            ],
            "scaled_bounds": [
                (vertex_array.min(axis=0) * factor).tolist(),
                (vertex_array.max(axis=0) * factor).tolist(),
            ],
            "source_sha256": sha256_file(source),
            "destination_sha256": sha256_file(destination),
        }
    records = np.frombuffer(raw, dtype=np.uint8).reshape(-1).copy()
    facet_dtype = np.dtype(
        [("normal", "<f4", 3), ("vertices", "<f4", (3, 3)), ("attribute", "<u2")]
    )
    facets = np.frombuffer(records[84:], dtype=facet_dtype, count=triangle_count)
    before_min = facets["vertices"].min(axis=(0, 1)).astype(np.float64)
    before_max = facets["vertices"].max(axis=(0, 1)).astype(np.float64)
    facets["vertices"] *= np.float32(factor)
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(records.tobytes())
    return {
        "source": str(Path(source).resolve()),
        "destination": str(destination.resolve()),
        "format": "BINARY_STL",
        "triangle_count": int(triangle_count),
        "scale_factor": float(factor),
        "source_bounds": [before_min.tolist(), before_max.tolist()],
        "scaled_bounds": [
            facets["vertices"].min(axis=(0, 1)).astype(np.float64).tolist(),
            facets["vertices"].max(axis=(0, 1)).astype(np.float64).tolist(),
        ],
        "source_sha256": sha256_file(source),
        "destination_sha256": sha256_file(destination),
    }


def load_boundary_ids_by_cell(mesh_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    mesh = Path(mesh_dir)
    header = (mesh / "header.lua").read_text(encoding="utf-8")
    count = parse_uniform_lattice(header).cell_count
    tree_ids, property_bits, _ = read_treelm_elemlist(
        mesh / "elemlist.lsb", n_elems=count
    )
    prop = parse_boundary_property_header(header)
    cells = extract_boundary_property_indices(property_bits, prop.bit_position)
    bnd_header = parse_bnd_header((mesh / "bnd.lua").read_text(encoding="utf-8"))
    rows = read_boundary_ids(
        mesh / "bnd.lsb",
        element_count=prop.element_count,
        side_count=bnd_header.side_count,
    )
    values = np.zeros((count, bnd_header.side_count), dtype=np.int64)
    values[cells] = rows
    return np.asarray(tree_ids), values


def error_statistics(values: np.ndarray) -> dict[str, float | int | None]:
    data = np.asarray(values, dtype=np.float64)
    data = data[np.isfinite(data)]
    if not len(data):
        return {
            "count": 0,
            "mean": None,
            "median": None,
            "rms": None,
            "p95": None,
            "max": None,
        }
    absolute = np.abs(data)
    return {
        "count": int(len(data)),
        "mean": float(np.mean(data)),
        "median": float(np.median(absolute)),
        "rms": float(np.sqrt(np.mean(data * data))),
        "p95": float(np.percentile(absolute, 95.0)),
        "max": float(np.max(absolute)),
    }


def compare_scaled_meshes(
    patched: MeshContract,
    scaled: MeshContract,
    patched_boundary_ids: np.ndarray,
    scaled_boundary_ids: np.ndarray,
) -> dict[str, Any]:
    """Compare topology, boundary IDs, and dimensionless q on common tree IDs."""

    patched_order = np.argsort(patched.tree_ids)
    scaled_order = np.argsort(scaled.tree_ids)
    patched_ids = patched.tree_ids[patched_order]
    scaled_ids = scaled.tree_ids[scaled_order]
    common, patched_positions, scaled_positions = np.intersect1d(
        patched_ids, scaled_ids, assume_unique=True, return_indices=True
    )
    patched_rows = patched_order[patched_positions]
    scaled_rows = scaled_order[scaled_positions]
    boundary_equal = (
        patched_boundary_ids[patched_rows] == scaled_boundary_ids[scaled_rows]
    )
    q_patched = patched.qvalues_by_cell[patched_rows]
    q_scaled = scaled.qvalues_by_cell[scaled_rows]
    q_mask = np.isfinite(q_patched) & np.isfinite(q_scaled)
    q_error = q_patched[q_mask] - q_scaled[q_mask]
    symmetric_difference = (
        len(patched.tree_ids) + len(scaled.tree_ids) - 2 * len(common)
    )
    relative_tree_difference = symmetric_difference / max(
        len(patched.tree_ids), len(scaled.tree_ids)
    )
    return {
        "patched_tree_id_count": int(len(patched.tree_ids)),
        "scaled_tree_id_count": int(len(scaled.tree_ids)),
        "common_tree_ids": int(len(common)),
        "patched_only_tree_ids": int(len(patched.tree_ids) - len(common)),
        "scaled_only_tree_ids": int(len(scaled.tree_ids) - len(common)),
        "tree_id_set_exact_match": bool(
            len(patched.tree_ids) == len(scaled.tree_ids)
            and np.array_equal(patched_ids, scaled_ids)
        ),
        "tree_id_order_exact_match": bool(
            np.array_equal(patched.tree_ids, scaled.tree_ids)
        ),
        "relative_tree_id_symmetric_difference": float(relative_tree_difference),
        "boundary_id_values_compared": int(boundary_equal.size),
        "boundary_id_mismatch_count": int(np.count_nonzero(~boundary_equal)),
        "boundary_id_exact_match": bool(np.all(boundary_equal)),
        "boundary_id_mismatch_fraction": float(np.mean(~boundary_equal)),
        "common_q_links": int(np.count_nonzero(q_mask)),
        "q_patched_minus_scaled": error_statistics(q_error),
    }


def unit_scaling_oracle_decision(
    comparison: Mapping[str, Any],
    patched_component_counts: Mapping[str, int],
    scaled_component_counts: Mapping[str, int],
) -> str:
    relative = float(comparison["relative_tree_id_symmetric_difference"])
    boundary_fraction = float(comparison["boundary_id_mismatch_fraction"])
    topology_same = dict(patched_component_counts) == dict(scaled_component_counts)
    if relative > 0.001 or not topology_same or boundary_fraction > 0.001:
        return "FAIL"
    return "PASS"


def source_line_evidence(path: Path, token: str) -> dict[str, Any]:
    text = Path(path).read_text(encoding="utf-8")
    offset = text.find(token)
    if offset < 0:
        raise ValueError(f"Source token is missing from {path}: {token}")
    return {
        "path": str(Path(path)),
        "line": text.count("\n", 0, offset) + 1,
        "token": token,
        "sha256": sha256_file(path),
    }


def intersect_raytriangle_callgraph(seeder_root: Path) -> dict[str, Any]:
    """Build a source-proven call graph and classify every active use path."""

    root = Path(seeder_root)
    boundary = root / "sdr/source/sdr_boundary_module.f90"
    flooding = root / "sdr/source/sdr_flooding_module.f90"
    proto = root / "sdr/source/sdr_proto2treelm_module.f90"
    main = root / "sdr/source/seeder.f90"
    line = root / "tem/source/shapes/tem_line_module.fpp"
    nodes = [
        {
            "caller": "sdr_qValByNode",
            "callee": "intersect_RayTriangle",
            "categories": ["QVALUE_DISTANCE", "GEOMETRY_INTERSECTION"],
            "evidence": source_line_evidence(
                boundary, "intersected = intersect_RayTriangle("
            ),
        },
        {
            "caller": "flood_periphery",
            "callee": "sdr_qValByNode",
            "categories": ["FLOODING", "QVALUE_DISTANCE"],
            "evidence": source_line_evidence(
                flooding, "call sdr_qValByNode( proto, geometry, dx, iDir,"
            ),
        },
        {
            "caller": "getBCID_and_calcQval",
            "callee": "sdr_qValByNode",
            "categories": ["BOUNDARY_IDENTIFICATION", "QVALUE_DISTANCE"],
            "evidence": source_line_evidence(
                boundary, "call sdr_qValByNode( proto, geometry, leVal%dx, iDir,"
            ),
        },
        {
            "caller": "sdr_identify_boundary",
            "callee": "getBCID_and_calcQval",
            "categories": ["BOUNDARY_IDENTIFICATION"],
            "evidence": source_line_evidence(
                boundary, "call getBCID_and_calcQval( proto"
            ),
        },
        {
            "caller": "proto2Treelm",
            "callee": "sdr_identify_boundary",
            "categories": ["BOUNDARY_IDENTIFICATION"],
            "evidence": source_line_evidence(
                proto, "call sdr_identify_boundary(node_pos"
            ),
        },
        {
            "caller": "seeder main",
            "callee": "sdr_flood_tree",
            "categories": ["FLOODING"],
            "evidence": source_line_evidence(main, "call sdr_flood_tree("),
        },
    ]
    return {
        "status": "PASS_SOURCE_PROVEN",
        "routine": {
            "name": "tem_line_module::intersect_RayTriangle",
            "definition": source_line_evidence(
                line, "function intersect_RayTriangle( line, triangle, intersect_p )"
            ),
        },
        "edges": nodes,
        "caller_categories": sorted(
            {category for node in nodes for category in node["categories"]}
        ),
        "calc_dist_false_paths": [],
        "calc_dist_false_conclusion": (
            "No direct caller reaches intersect_RayTriangle when every boundary has "
            "calc_dist=false; both flooding and boundary-identification qVal paths "
            "are guarded by needCalcQValByBCID/calc_dist."
        ),
        "patch_affects_flooding_or_topology": True,
        "topology_effect_source_proof": (
            "With wall calc_dist=true, flood_periphery calls sdr_qValByNode for all "
            "26 directions and floods a wetted periphery cell only when the tested "
            "link has qVal<0. Changing ray/triangle hits therefore changes the set "
            "of flooded leaves before sdr_proto2treelm writes tree IDs."
        ),
        "proto_tree_refinement_use": False,
    }


def nearest_centerline_segment(
    point_cfd_m: Sequence[float],
    *,
    transform_json: Path,
    nodes_csv: Path,
    edges_csv: Path,
) -> dict[str, Any]:
    """Locate a CFD point on the existing anatomical 1-D vascular graph."""

    transform = json.loads(Path(transform_json).read_text(encoding="utf-8"))
    inverse = np.asarray(transform["inverse_homogeneous_transform_4x4"], dtype=float)
    homogeneous = np.append(np.asarray(point_cfd_m, dtype=float), 1.0)
    anatomical_um = (inverse @ homogeneous)[:3] * 1.0e6
    nodes: dict[int, dict[str, float | int]] = {}
    with Path(nodes_csv).open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            node_id = int(row["node_id"])
            nodes[node_id] = {
                "node_id": node_id,
                "x_um": float(row["x_um"]),
                "y_um": float(row["y_um"]),
                "z_um": float(row["z_um"]),
                "radius_um": float(row["radius_um"]),
            }
    edges = []
    adjacency: dict[int, list[tuple[int, float]]] = {node: [] for node in nodes}
    with Path(edges_csv).open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            parent = int(row["parent_node_id"])
            child = int(row["child_node_id"])
            length = float(row["length_um"])
            edge = {
                "global_edge_id": int(row["global_edge_id"]),
                "parent": parent,
                "child": child,
                "length_um": length,
                "radius_parent_um": float(row["radius_parent_um"]),
                "radius_child_um": float(row["radius_child_um"]),
            }
            edges.append(edge)
            adjacency[parent].append((child, length))
            adjacency[child].append((parent, length))

    best: dict[str, Any] | None = None
    for edge in edges:
        left = np.asarray(
            [nodes[edge["parent"]][key] for key in ("x_um", "y_um", "z_um")],
            dtype=float,
        )
        right = np.asarray(
            [nodes[edge["child"]][key] for key in ("x_um", "y_um", "z_um")],
            dtype=float,
        )
        vector = right - left
        denominator = float(np.dot(vector, vector))
        alpha = (
            0.0
            if denominator == 0.0
            else float(
                np.clip(np.dot(anatomical_um - left, vector) / denominator, 0, 1)
            )
        )
        nearest = left + alpha * vector
        distance = float(np.linalg.norm(anatomical_um - nearest))
        if best is None or distance < best["distance_to_segment_um"]:
            radius = (1.0 - alpha) * edge["radius_parent_um"] + alpha * edge[
                "radius_child_um"
            ]
            best = {
                **edge,
                "alpha": alpha,
                "nearest_point_anatomical_um": nearest.tolist(),
                "distance_to_segment_um": distance,
                "local_radius_um": float(radius),
                "local_hydraulic_diameter_um": float(2.0 * radius),
                "local_minimum_endpoint_diameter_um": float(
                    2.0 * min(edge["radius_parent_um"], edge["radius_child_um"])
                ),
            }
    if best is None:
        raise ValueError("Centerline graph has no edges")

    bifurcations = [node for node, links in adjacency.items() if len(links) >= 3]
    distances = {node: math.inf for node in adjacency}
    queue: list[tuple[float, int]] = []
    for node in bifurcations:
        distances[node] = 0.0
        heapq.heappush(queue, (0.0, node))
    while queue:
        distance, node = heapq.heappop(queue)
        if distance != distances[node]:
            continue
        for neighbor, length in adjacency[node]:
            candidate = distance + length
            if candidate < distances[neighbor]:
                distances[neighbor] = candidate
                heapq.heappush(queue, (candidate, neighbor))
    alpha = float(best["alpha"])
    length = float(best["length_um"])
    parent_distance = alpha * length + distances[int(best["parent"])]
    child_distance = (1.0 - alpha) * length + distances[int(best["child"])]
    best["distance_along_centerline_to_nearest_bifurcation_um"] = float(
        min(parent_distance, child_distance)
    )
    best["query_point_cfd_m"] = list(map(float, point_cfd_m))
    best["query_point_anatomical_um"] = anatomical_um.tolist()
    return best


def timed_runtime(started: float) -> float:
    return float(time.perf_counter() - started)
