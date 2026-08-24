"""Official ``ultraVessMorpho2Mesh`` backend for one saved CORE ROI.

This module deliberately contains no surface reconstruction algorithm.  It
serializes the immutable ROI geometry, invokes the upstream Ultraliser binary,
and records enough provenance to reproduce the invocation.
"""

from __future__ import annotations

import hashlib
import math
import os
import re
import shlex
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import h5py
import networkx as nx
import numpy as np

from utils.sampling.sampling_types import ROIRecord

from .config import CFDLumenConfig
from .export import write_csv, write_json


ULTRALISER_EXPECTED_ANCHOR = 3274
ULTRALISER_MAX_VOXELS = 150_000_000
ULTRALISER_TARGET_CELLS = 12
ULTRALISER_FALLBACK_CELLS = 10
ULTRALISER_EDGE_GAP = 0.05
ULTRALISER_PACKING_ALGORITHM = "polylines-with-spheres"
ULTRALISER_VOXELIZATION_AXIS = "xyz"


class UltraliserBackendError(RuntimeError):
    """Base exception carrying a stable diagnostic code."""

    code = "ULTRALISER_BACKEND_ERROR"

    def __init__(self, message: str) -> None:
        super().__init__(f"{self.code}: {message}")


class UltraliserCycleConflict(UltraliserBackendError):
    code = "ULTRALISER_SWC_SERIALIZATION_CYCLE_CONFLICT"


class UltraliserBuildBlocked(UltraliserBackendError):
    code = "ULTRALISER_BUILD_ENVIRONMENT_BLOCKED"


class UltraliserResolutionBlocked(UltraliserBackendError):
    code = "ULTRALISER_RESOLUTION_MEMORY_LIMIT"


class UltraliserRunFailed(UltraliserBackendError):
    code = "ULTRALISER_EXECUTION_FAILED"


@dataclass(frozen=True, slots=True)
class UltraliserLayout:
    run_root: Path
    input: Path
    build: Path
    smoke: Path
    final: Path
    raw_ultraliser: Path
    qc: Path
    figures: Path
    report: Path


@dataclass(frozen=True, slots=True)
class ROIUltraliserInput:
    swc_path: Path
    h5_path: Path
    metadata_path: Path
    cut_port_mapping_path: Path
    topology: dict[str, Any]
    sections: tuple[tuple[int, ...], ...]


@dataclass(frozen=True, slots=True)
class UltraliserInvocation:
    stage: str
    command: tuple[str, ...]
    command_text: str
    runtime_seconds: float
    return_code: int
    output_directory: Path


def create_ultraliser_layout(
    output_root: Path,
    *,
    roi: ROIRecord,
    run_id: str | None,
) -> UltraliserLayout:
    base = Path(output_root).resolve() / "ultraliser"
    base.mkdir(parents=True, exist_ok=True)
    if run_id is None:
        run_id = f"ultraliser_roi{roi.anchor_id:06d}_{datetime.now():%Y%m%d_%H%M%S}"
    if not re.fullmatch(r"[A-Za-z0-9._-]+", run_id):
        raise ValueError("--run-id may contain only letters, digits, dot, underscore, and hyphen")
    run_root = base / run_id
    if run_root.exists():
        raise FileExistsError(f"Refusing to overwrite an existing Ultraliser run: {run_root}")
    input_directory = run_root / "input"
    build = run_root / "build"
    smoke = run_root / "smoke"
    final = run_root / "final"
    raw = final / "raw_ultraliser"
    qc = run_root / "qc"
    figures = run_root / "figures"
    report = run_root / "report"
    for directory in (input_directory, build, smoke, final, raw, qc, figures, report):
        directory.mkdir(parents=True, exist_ok=False)
    return UltraliserLayout(
        run_root=run_root,
        input=input_directory,
        build=build,
        smoke=smoke,
        final=final,
        raw_ultraliser=raw,
        qc=qc,
        figures=figures,
        report=report,
    )


def _roi_graph(roi: ROIRecord) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(map(int, roi.local_node_ids))
    graph.add_edges_from((int(first), int(second)) for first, second in roi.local_edges)
    return graph


def _source_parent_map(roi: ROIRecord, graph: nx.Graph) -> tuple[dict[int, int], bool]:
    """Use saved parent-to-current edges, with deterministic fallback for serialization only."""

    directed = nx.DiGraph()
    directed.add_nodes_from(graph.nodes)
    directed.add_edges_from((int(first), int(second)) for first, second in roi.local_edges)
    roots = sorted(node for node, degree in directed.in_degree if degree == 0)
    source_direction_valid = (
        nx.is_directed_acyclic_graph(directed)
        and len(roots) == 1
        and all(degree <= 1 for _, degree in directed.in_degree)
        and len(nx.descendants(directed, roots[0])) + 1 == graph.number_of_nodes()
    )
    if source_direction_valid:
        parent = {int(roots[0]): -1}
        parent.update({int(child): int(upstream) for upstream, child in directed.edges})
        return parent, False

    root = int(roi.anchor_id) if int(roi.anchor_id) in graph else min(map(int, graph.nodes))
    parent = {root: -1}
    for upstream, child in nx.bfs_edges(graph, root, sort_neighbors=sorted):
        parent[int(child)] = int(upstream)
    return parent, True


def _ordered_tree_nodes(parent: dict[int, int]) -> list[int]:
    children: dict[int, list[int]] = {node: [] for node in parent}
    root = next(node for node, value in parent.items() if value < 0)
    for node, upstream in parent.items():
        if upstream >= 0:
            children[upstream].append(node)
    for values in children.values():
        values.sort()
    ordered: list[int] = []
    queue = [root]
    while queue:
        current = queue.pop(0)
        ordered.append(current)
        queue.extend(children[current])
    return ordered


def _directed_sections(parent: dict[int, int], graph: nx.Graph) -> tuple[tuple[int, ...], ...]:
    children: dict[int, list[int]] = {node: [] for node in parent}
    root = next(node for node, upstream in parent.items() if upstream < 0)
    for node, upstream in parent.items():
        if upstream >= 0:
            children[upstream].append(node)
    for values in children.values():
        values.sort()

    sections: list[tuple[int, ...]] = []
    starts = [root] + sorted(node for node in graph if node != root and graph.degree[node] != 2)
    seen_edges: set[tuple[int, int]] = set()
    for start in starts:
        for child in children[start]:
            edge = (int(start), int(child))
            if edge in seen_edges:
                continue
            path = [int(start), int(child)]
            seen_edges.add(edge)
            current = int(child)
            while graph.degree[current] == 2 and len(children[current]) == 1:
                following = int(children[current][0])
                seen_edges.add((current, following))
                path.append(following)
                current = following
            sections.append(tuple(path))

    expected = {(int(upstream), int(node)) for node, upstream in parent.items() if upstream >= 0}
    if seen_edges != expected:
        missing = sorted(expected - seen_edges)
        raise UltraliserBackendError(f"section decomposition missed source edges: {missing[:5]}")
    return tuple(sections)


def _write_swc(
    roi: ROIRecord,
    path: Path,
    parent: dict[int, int],
) -> dict[int, int]:
    ordered = _ordered_tree_nodes(parent)
    swc_id = {local_id: index + 1 for index, local_id in enumerate(ordered)}
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float)
    lines = [
        "# Canonical CORE ROI morphology for Ultraliser provenance",
        "# units: x y z radius are micrometres (um)",
        "# columns: node_id type x_um y_um z_um radius_um parent_id",
    ]
    for local_id in ordered:
        position = positions[local_id]
        parent_id = -1 if parent[local_id] < 0 else swc_id[parent[local_id]]
        lines.append(
            f"{swc_id[local_id]} 3 "
            f"{position[0]:.17g} {position[1]:.17g} {position[2]:.17g} "
            f"{radii[local_id]:.17g} {parent_id}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return swc_id


def _write_h5_adapter(
    roi: ROIRecord,
    path: Path,
    sections: tuple[tuple[int, ...], ...],
) -> dict[str, Any]:
    """Write the unmodified geometry in Ultraliser's supported vascular H5 schema.

    The upstream reader interprets the fourth ``points`` value as diameter and
    multiplies it by 0.5.  Storing ``2 * source_radius`` therefore preserves the
    source radius used by the official proxy generator.
    """

    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float)
    point_rows: list[list[float]] = []
    structure: list[tuple[int, int]] = []
    for section in sections:
        first = len(point_rows)
        for local_id in section:
            point_rows.append([*map(float, positions[local_id]), 2.0 * float(radii[local_id])])
        structure.append((first, len(point_rows) - 1))

    connectivity: list[tuple[int, int]] = []
    for parent_index, path_nodes in enumerate(sections):
        terminal = path_nodes[-1]
        for child_index, child_path in enumerate(sections):
            if child_index != parent_index and child_path[0] == terminal:
                connectivity.append((parent_index, child_index))
    points = np.asarray(point_rows, dtype=np.float32)
    structure_array = np.asarray(structure, dtype=np.int64).reshape((-1, 2))
    connectivity_array = np.asarray(connectivity, dtype=np.int64).reshape((-1, 2))
    with h5py.File(path, "w") as stream:
        stream.create_dataset("points", data=points)
        stream.create_dataset("structure", data=structure_array)
        stream.create_dataset("connectivity", data=connectivity_array)
        stream.attrs["canonical_source"] = "roi_core.swc"
        stream.attrs["coordinate_unit"] = "um"
        stream.attrs["points_fourth_column"] = "diameter_um"
        stream.attrs["source_radius_preserved_after_official_reader_times_0_5"] = True

    recovered_positions = points[:, :3].astype(float)
    recovered_radii = 0.5 * points[:, 3].astype(float)
    expanded_positions = np.asarray(
        [positions[local_id] for section in sections for local_id in section], dtype=float
    )
    expanded_radii = np.asarray(
        [radii[local_id] for section in sections for local_id in section], dtype=float
    )
    return {
        "adapter_format": "Ultraliser vascular H5",
        "reason": (
            "local ultraVessMorpho2Mesh commit accepts only H5/VMV; H5 avoids the VMV "
            "reader's fixed 0.8 radius multiplier"
        ),
        "canonical_input": "roi_core.swc",
        "executable_input": path.name,
        "points_fourth_column": "diameter_um",
        "official_reader_operation": "internal_radius_um = points[:,3] * 0.5",
        "section_count": len(sections),
        "connectivity_count": len(connectivity),
        "point_row_count_with_section_endpoint_duplicates": len(points),
        "maximum_float32_coordinate_quantization_um": float(
            np.max(np.abs(recovered_positions - expanded_positions))
        ),
        "maximum_float32_radius_quantization_um": float(
            np.max(np.abs(recovered_radii - expanded_radii))
        ),
        "source_geometry_resampled": False,
        "source_radius_scaled_internally": False,
    }


def export_roi_for_ultraliser(roi: ROIRecord, input_directory: Path) -> ROIUltraliserInput:
    graph = _roi_graph(roi)
    component_count = nx.number_connected_components(graph) if graph else 0
    cycle_rank = graph.number_of_edges() - graph.number_of_nodes() + component_count
    if component_count != 1:
        raise UltraliserBackendError(
            f"ROI must contain exactly one connected component, got {component_count}"
        )
    if cycle_rank > 0:
        raise UltraliserCycleConflict(
            f"ROI {roi.roi_id} has cycle_rank={cycle_rank}; no cyclic edge was removed"
        )
    if int(roi.anchor_id) != ULTRALISER_EXPECTED_ANCHOR:
        raise UltraliserBackendError(
            f"this formal experiment is restricted to anchor_{ULTRALISER_EXPECTED_ANCHOR:06d}"
        )

    parent, serialization_direction_only = _source_parent_map(roi, graph)
    sections = _directed_sections(parent, graph)
    swc_path = input_directory / "roi_core.swc"
    h5_path = input_directory / "roi_core.h5"
    metadata_path = input_directory / "roi_core_metadata.json"
    cut_mapping_path = input_directory / "cut_port_mapping.csv"
    swc_ids = _write_swc(roi, swc_path, parent)
    adapter = _write_h5_adapter(roi, h5_path, sections)

    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float)
    edges = np.asarray(roi.local_edges, dtype=np.int64)
    degree = np.bincount(edges.ravel(), minlength=roi.node_count)
    topology = {
        "connected_components": int(component_count),
        "cycle_rank": int(cycle_rank),
        "source_parent_child_relation_used": not serialization_direction_only,
        "serialization_direction_only": bool(serialization_direction_only),
        "serialization_is_physiological_flow_direction": False,
        "root_local_node_id": int(next(node for node, value in parent.items() if value < 0)),
    }
    metadata = {
        "roi_id": roi.roi_id,
        "source_model_id": roi.source_model_id,
        "anchor_id": int(roi.anchor_id),
        "scope": "CORE_ROI",
        "saved_roi_reused": True,
        "roi_resampled": False,
        "roi_regenerated": False,
        "coordinate_unit": "um",
        "radius_unit": "um",
        "source_node_count": roi.node_count,
        "source_edge_count": roi.edge_count,
        "branch_count": len(sections),
        "bifurcation_count": int(np.count_nonzero(degree >= 3)),
        "r_min_um": float(radii.min()),
        "r_median_um": float(np.median(radii)),
        "r_max_um": float(radii.max()),
        "bbox_min_um": positions.min(axis=0).tolist(),
        "bbox_max_um": positions.max(axis=0).tolist(),
        "bbox_dimensions_um": np.ptp(positions, axis=0).tolist(),
        "saved_roi_bbox_min_um": list(map(float, roi.bbox_min_um)),
        "saved_roi_bbox_max_um": list(map(float, roi.bbox_max_um)),
        "cut_port_count": len(roi.cut_ports),
        "all_cut_ports_preserved": True,
        "topology": topology,
        "swc_node_id_by_local_node_id": {str(key): value for key, value in swc_ids.items()},
        "ultraliser_input_adapter": adapter,
    }
    write_json(metadata_path, metadata)
    write_csv(
        cut_mapping_path,
        [
            {
                "cut_port_id": port.cut_port_id,
                "swc_node_id": swc_ids[int(port.local_node_id)],
                "global_edge_id": int(port.global_edge_id),
                "x_um": float(port.intersection_position_um[0]),
                "y_um": float(port.intersection_position_um[1]),
                "z_um": float(port.intersection_position_um[2]),
                "radius_um": float(port.radius_at_cut_um),
                "boundary_face": port.boundary_face,
            }
            for port in roi.cut_ports
        ],
        fieldnames=(
            "cut_port_id",
            "swc_node_id",
            "global_edge_id",
            "x_um",
            "y_um",
            "z_um",
            "radius_um",
            "boundary_face",
        ),
    )
    return ROIUltraliserInput(
        swc_path=swc_path,
        h5_path=h5_path,
        metadata_path=metadata_path,
        cut_port_mapping_path=cut_mapping_path,
        topology=topology,
        sections=sections,
    )


def _resolution_case(
    positions: np.ndarray,
    radii: np.ndarray,
    *,
    cells_across_minimum_diameter: int,
    voxels_per_micron: float,
) -> dict[str, Any]:
    p_min = np.min(positions - radii[:, None], axis=0)
    p_max = np.max(positions + radii[:, None], axis=0)
    morphology_dimensions = p_max - p_min
    largest = float(morphology_dimensions.max())
    base_resolution = max(1, int(float(voxels_per_micron) * largest))
    expanded_dimensions = morphology_dimensions * (1.0 + 2.0 * ULTRALISER_EDGE_GAP)
    voxel_size = float(expanded_dimensions.max()) / base_resolution
    dimensions = np.maximum(1, np.rint(expanded_dimensions / voxel_size).astype(np.int64))
    voxel_count = int(np.prod(dimensions, dtype=np.int64))
    return {
        "cells_across_minimum_diameter": int(cells_across_minimum_diameter),
        "voxels_per_micron": float(voxels_per_micron),
        "ultraliser_base_resolution": base_resolution,
        "predicted_Nx": int(dimensions[0]),
        "predicted_Ny": int(dimensions[1]),
        "predicted_Nz": int(dimensions[2]),
        "predicted_voxel_count": voxel_count,
        "predicted_bit_volume_bytes": int(math.ceil(voxel_count / 8.0)),
        "predicted_byte_equivalent_bytes": voxel_count,
        "predicted_bit_volume_mib": float(voxel_count / 8.0 / 1024.0**2),
        "predicted_byte_equivalent_mib": float(voxel_count / 1024.0**2),
        "edge_gap_expansion_ratio": ULTRALISER_EDGE_GAP,
        "morphology_radius_inclusive_bbox_dimensions_um": morphology_dimensions.tolist(),
        "expanded_bbox_dimensions_um": expanded_dimensions.tolist(),
        "predicted_voxel_size_um": voxel_size,
        "predicted_effective_voxels_per_micron": 1.0 / voxel_size,
    }


def estimate_resolution(roi: ROIRecord) -> dict[str, Any]:
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float)
    r_min = float(radii.min())
    d_min = 2.0 * r_min
    target_vpm = int(math.ceil(ULTRALISER_TARGET_CELLS / d_min))
    requested = _resolution_case(
        positions,
        radii,
        cells_across_minimum_diameter=ULTRALISER_TARGET_CELLS,
        voxels_per_micron=target_vpm,
    )
    selected = requested
    fallback_applied = False
    if requested["predicted_voxel_count"] > ULTRALISER_MAX_VOXELS:
        target_vpm = int(math.ceil(ULTRALISER_FALLBACK_CELLS / d_min))
        selected = _resolution_case(
            positions,
            radii,
            cells_across_minimum_diameter=ULTRALISER_FALLBACK_CELLS,
            voxels_per_micron=target_vpm,
        )
        fallback_applied = True
    if selected["predicted_voxel_count"] > ULTRALISER_MAX_VOXELS:
        raise UltraliserResolutionBlocked(
            f"10-cell fallback still predicts {selected['predicted_voxel_count']:,} voxels"
        )
    smoke_vpm = max(2.0, 0.5 * float(selected["voxels_per_micron"]))
    smoke = _resolution_case(
        positions,
        radii,
        cells_across_minimum_diameter=0,
        voxels_per_micron=smoke_vpm,
    )
    return {
        "r_min_um": r_min,
        "D_min_um": d_min,
        "formula": "ceil(target_cells_across_minimum_diameter / D_min_um)",
        "maximum_allowed_voxels": ULTRALISER_MAX_VOXELS,
        "requested_12_cell_case": requested,
        "fallback_to_10_cells_applied": fallback_applied,
        "selected_final_case": selected,
        "smoke_case": smoke,
    }


def discover_ultraliser_executable(ultraliser_root: Path) -> Path:
    root = Path(ultraliser_root).resolve()
    candidates = [
        root / "build-wsl" / "bin" / "ultraVessMorpho2Mesh",
        root / "build" / "bin" / "ultraVessMorpho2Mesh",
        root / "build" / "bin" / "ultraVessMorpho2Mesh.exe",
    ]
    candidates.extend(root.glob("**/bin/ultraVessMorpho2Mesh"))
    candidates.extend(root.glob("**/bin/ultraVessMorpho2Mesh.exe"))
    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()
    raise UltraliserBuildBlocked(f"no ultraVessMorpho2Mesh executable found under {root}")


def _wsl_path(path: Path) -> str:
    resolved = Path(path).resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        return resolved.as_posix()
    remainder = resolved.as_posix()[2:].lstrip("/")
    return f"/mnt/{drive}/{remainder}"


def _execution_command(executable: Path, arguments: list[str]) -> list[str]:
    if os.name == "nt" and executable.suffix.lower() != ".exe":
        converted = [
            _wsl_path(Path(value)) if index in {1, 3} else value
            for index, value in enumerate(arguments)
        ]
        return ["wsl", "-d", "Ubuntu", "--", _wsl_path(executable), *converted]
    return [str(executable), *arguments]


def _help_command(executable: Path) -> list[str]:
    if os.name == "nt" and executable.suffix.lower() != ".exe":
        return ["wsl", "-d", "Ubuntu", "--", _wsl_path(executable), "--help"]
    return [str(executable), "--help"]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def capture_ultraliser_build_metadata(
    ultraliser_root: Path,
    executable: Path,
    build_directory: Path,
) -> dict[str, Any]:
    command = _help_command(executable)
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    help_text = completed.stdout + completed.stderr
    (build_directory / "ultraliser_help.txt").write_text(help_text, encoding="utf-8")
    if completed.returncode != 0:
        raise UltraliserRunFailed(f"--help returned {completed.returncode}")
    commit = subprocess.run(
        ["git", "-C", str(Path(ultraliser_root).resolve()), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=True,
    ).stdout.strip()
    dirty_rows = subprocess.run(
        ["git", "-C", str(Path(ultraliser_root).resolve()), "status", "--porcelain"],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=True,
    ).stdout.splitlines()
    payload = {
        "ultraliser_git_commit": commit,
        "ultraliser_source_dirty": bool(dirty_rows),
        "ultraliser_source_status": dirty_rows,
        "executable_path_windows": str(executable),
        "executable_path_execution_environment": _wsl_path(executable)
        if os.name == "nt" and executable.suffix.lower() != ".exe"
        else str(executable),
        "executable_sha256": _sha256(executable),
        "build_mode": "Release",
        "build_environment": "WSL2 Ubuntu",
        "compiler_compatibility_flag": "-include cstdint",
        "ultraliser_cpp_source_modified": False,
        "help_command": shlex.join(command),
        "help_return_code": completed.returncode,
        "vascular_reader_formats_in_local_source": ["h5", "vmv"],
        "vascular_swc_direct_read_supported_by_local_source": False,
        "packing_algorithm": ULTRALISER_PACKING_ALGORITHM,
        "packing_algorithm_basis": (
            "paper vascular workflow and source POLYLINE_SPHERE_PACKING branch; help default is polylines"
        ),
    }
    write_json(build_directory / "ultraliser_version.json", payload)
    return payload


def build_ultraliser_command(
    executable: Path,
    morphology: Path,
    output_directory: Path,
    *,
    prefix: str,
    voxels_per_micron: float,
    threads: int,
) -> tuple[list[str], str]:
    arguments = [
        "--morphology",
        str(morphology.resolve()),
        "--output-directory",
        str(output_directory.resolve()),
        "--prefix",
        prefix,
        "--scaled-resolution",
        "--voxels-per-micron",
        f"{voxels_per_micron:g}",
        "--solid",
        "--voxelization-axis",
        ULTRALISER_VOXELIZATION_AXIS,
        "--packing-algorithm",
        ULTRALISER_PACKING_ALGORITHM,
        "--isosurface-technique",
        "dmc",
        "--adaptive-optimization",
        "--optimization-iterations",
        "5",
        "--smooth-iterations",
        "5",
        "--laplacian-iterations",
        "10",
        "--export-stl-mesh",
        "--stats",
        "--threads",
        str(threads),
    ]
    command = _execution_command(executable, arguments)
    return command, shlex.join(command)


def invoke_ultraliser(
    executable: Path,
    morphology: Path,
    output_directory: Path,
    *,
    stage: str,
    prefix: str,
    voxels_per_micron: float,
    threads: int,
    command_alias_path: Path | None = None,
) -> UltraliserInvocation:
    command, command_text = build_ultraliser_command(
        executable,
        morphology,
        output_directory,
        prefix=prefix,
        voxels_per_micron=voxels_per_micron,
        threads=threads,
    )
    output_directory.mkdir(parents=True, exist_ok=True)
    (output_directory / "command.txt").write_text(command_text + "\n", encoding="utf-8")
    if command_alias_path is not None:
        command_alias_path.write_text(command_text + "\n", encoding="utf-8")
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        check=False,
    )
    runtime = time.perf_counter() - started
    (output_directory / "stdout.log").write_text(completed.stdout, encoding="utf-8")
    (output_directory / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    invocation = UltraliserInvocation(
        stage=stage,
        command=tuple(command),
        command_text=command_text,
        runtime_seconds=float(runtime),
        return_code=int(completed.returncode),
        output_directory=output_directory,
    )
    write_json(
        output_directory / "invocation.json",
        {
            "stage": stage,
            "command": command_text,
            "runtime_seconds": runtime,
            "return_code": completed.returncode,
            "threads": threads,
            "scaled_resolution": True,
            "voxels_per_micron": voxels_per_micron,
            "solid_voxelization": True,
            "voxelization_axis": ULTRALISER_VOXELIZATION_AXIS,
            "packing_algorithm": ULTRALISER_PACKING_ALGORITHM,
        },
    )
    if completed.returncode != 0:
        raise UltraliserRunFailed(
            f"{stage} returned {completed.returncode}; see {output_directory / 'stderr.log'}"
        )
    return invocation


def discover_watertight_stl(output_directory: Path) -> Path:
    candidates = sorted(output_directory.glob("meshes/*-watertight.stl"))
    if not candidates:
        raise UltraliserRunFailed(
            f"official watertight STL was not produced under {output_directory}"
        )
    return candidates[-1]


def run_ultraliser_roi_experiment(
    roi: ROIRecord,
    config: CFDLumenConfig,
    *,
    output_root: Path,
    run_id: str | None,
    ultraliser_root: Path,
    executable_path: Path | None,
    previous_surface: Path | None,
    threads: int | None,
) -> dict[str, Any]:
    """Run one smoke reconstruction and one formal reconstruction for anchor_003274."""

    layout = create_ultraliser_layout(output_root, roi=roi, run_id=run_id)
    exported = export_roi_for_ultraliser(roi, layout.input)
    resolution = estimate_resolution(roi)
    write_json(layout.qc / "resolution_estimate.json", resolution)
    executable = (
        Path(executable_path).resolve()
        if executable_path is not None
        else discover_ultraliser_executable(ultraliser_root)
    )
    build = capture_ultraliser_build_metadata(ultraliser_root, executable, layout.build)
    selected_threads = int(threads or min(8, os.cpu_count() or 1))

    smoke_vpm = float(resolution["smoke_case"]["voxels_per_micron"])
    smoke = invoke_ultraliser(
        executable,
        exported.h5_path,
        layout.smoke,
        stage="smoke",
        prefix="roi003274_smoke",
        voxels_per_micron=smoke_vpm,
        threads=selected_threads,
    )
    smoke_surface = discover_watertight_stl(layout.smoke)

    final_vpm = float(resolution["selected_final_case"]["voxels_per_micron"])
    final = invoke_ultraliser(
        executable,
        exported.h5_path,
        layout.raw_ultraliser,
        stage="final",
        prefix="roi003274_final",
        voxels_per_micron=final_vpm,
        threads=selected_threads,
        command_alias_path=layout.run_root / "ultraliser_command.txt",
    )
    raw_surface = discover_watertight_stl(layout.raw_ultraliser)

    from .ultraliser_qc import finalize_ultraliser_outputs

    qc_result = finalize_ultraliser_outputs(
        roi,
        config,
        layout=layout,
        raw_surface=raw_surface,
        previous_surface=previous_surface,
    )
    summary = {
        "roi_id": roi.roi_id,
        "anchor_id": roi.anchor_id,
        "run_root": str(layout.run_root),
        "canonical_swc": str(exported.swc_path),
        "ultraliser_morphology_input": str(exported.h5_path),
        "ultraliser_input_adapter_required": True,
        "ultraliser_direct_swc_reader_available": False,
        "ultraliser_git_commit": build["ultraliser_git_commit"],
        "ultraliser_executable": build["executable_path_execution_environment"],
        "packing_algorithm": ULTRALISER_PACKING_ALGORITHM,
        "packing_algorithm_basis": build["packing_algorithm_basis"],
        "scaled_resolution": True,
        "voxels_per_micron": final_vpm,
        "solid_voxelization": True,
        "voxelization_axis": ULTRALISER_VOXELIZATION_AXIS,
        "threads": selected_threads,
        "smoke_runtime_seconds": smoke.runtime_seconds,
        "final_runtime_seconds": final.runtime_seconds,
        "smoke_surface": str(smoke_surface),
        "raw_final_surface": str(raw_surface),
        "surface_qc": qc_result["surface_qc"],
        "radius_fidelity": qc_result["radius_fidelity"],
        "visual_assessment": "PENDING_MANUAL_REVIEW",
        "ULTRALISER_RAW_SURFACE_ACCEPTED": None,
        "decision": "PENDING_MANUAL_REVIEW",
        "old_v7_v8_v9_rerun": False,
        "cfd_port_handling_run": False,
        "third_ultraliser_run": False,
    }
    write_json(layout.qc / "run_summary.json", summary)
    from .ultraliser_qc import write_ultraliser_report

    write_ultraliser_report(layout.report / "ultraliser_roi003274_report.md", summary)
    return summary
