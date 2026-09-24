"""Official ``ultraVessMorpho2Mesh`` backend for saved vascular ROIs."""

from __future__ import annotations

import hashlib
import math
import os
import re
import shlex
import shutil
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

from .config import CFDLumenConfig, UltraliserConfig
from .export import write_csv, write_json


ULTRALISER_RADIUS_SCALE = 0.91
ULTRALISER_PACKING_ALGORITHM = "polylines-with-spheres"
ULTRALISER_VOXELIZATION_AXIS = "xyz"
ULTRALISER_ISOSURFACE_TECHNIQUE = "dmc"


class UltraliserBackendError(RuntimeError):
    code = "ULTRALISER_BACKEND_ERROR"

    def __init__(self, message: str) -> None:
        super().__init__(f"{self.code}: {message}")


class UltraliserCycleConflict(UltraliserBackendError):
    code = "ULTRALISER_SWC_SERIALIZATION_CYCLE_CONFLICT"


class UltraliserExecutableNotFound(UltraliserBackendError):
    code = "ULTRALISER_EXECUTABLE_NOT_FOUND"


class UltraliserRunFailed(UltraliserBackendError):
    code = "ULTRALISER_EXECUTION_FAILED"


@dataclass(frozen=True, slots=True)
class UltraliserLayout:
    run_root: Path
    input: Path
    geometry: Path
    qc: Path
    report: Path
    work: Path


@dataclass(frozen=True, slots=True)
class ROIUltraliserInput:
    swc_path: Path
    h5_path: Path
    metadata_path: Path
    radius_feed_mapping_path: Path
    cut_port_mapping_path: Path
    topology: dict[str, Any]
    sections: tuple[tuple[int, ...], ...]


@dataclass(frozen=True, slots=True)
class UltraliserInvocation:
    command: tuple[str, ...]
    command_text: str
    runtime_seconds: float
    return_code: int
    output_directory: Path


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def create_ultraliser_layout(output_root: Path, roi: ROIRecord, run_id: str | None) -> UltraliserLayout:
    root = Path(output_root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    selected_id = run_id or f"ultraliser_anchor{int(roi.anchor_id):06d}_{datetime.now():%Y%m%d_%H%M%S}"
    if not re.fullmatch(r"[A-Za-z0-9._-]+", selected_id):
        raise ValueError("--run-id may contain only letters, digits, dot, underscore, and hyphen")
    run_root = root / selected_id
    if run_root.exists():
        raise FileExistsError(f"Refusing to overwrite reconstruction run: {run_root}")
    directories = {
        name: run_root / name for name in ("input", "geometry", "qc", "report")
    }
    directories["work"] = run_root / ".ultraliser_work"
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=False)
    return UltraliserLayout(
        run_root=run_root,
        input=directories["input"],
        geometry=directories["geometry"],
        qc=directories["qc"],
        report=directories["report"],
        work=directories["work"],
    )


def _roi_graph(roi: ROIRecord) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(map(int, roi.local_node_ids))
    graph.add_edges_from((int(first), int(second)) for first, second in roi.local_edges)
    return graph


def _validate_roi_arrays(roi: ROIRecord) -> None:
    node_ids = np.asarray(roi.local_node_ids, dtype=np.int64)
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float)
    edges = np.asarray(roi.local_edges, dtype=np.int64).reshape((-1, 2))
    if not np.array_equal(node_ids, np.arange(len(node_ids), dtype=np.int64)):
        raise UltraliserBackendError("local node IDs must be contiguous zero-based indices")
    if positions.shape != (len(node_ids), 3) or not np.all(np.isfinite(positions)):
        raise UltraliserBackendError("source coordinates must be finite N x 3 micrometre values")
    if radii.shape != (len(node_ids),) or np.any(~np.isfinite(radii) | (radii <= 0.0)):
        raise UltraliserBackendError("all source radii must be finite and positive")
    if len(edges) == 0 or int(edges.min()) < 0 or int(edges.max()) >= len(node_ids):
        raise UltraliserBackendError("ROI edges are empty or reference an invalid local node")


def _source_parent_map(roi: ROIRecord, graph: nx.Graph) -> tuple[dict[int, int], bool]:
    """Use saved parent-to-current edges, with deterministic serialization if necessary."""

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
    anchor = int(roi.anchor_id)
    root = anchor if anchor in graph else min(map(int, graph.nodes))
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


def _directed_sections(
    parent: dict[int, int],
    graph: nx.Graph,
) -> tuple[tuple[int, ...], ...]:
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


def _write_swc(roi: ROIRecord, path: Path, parent: dict[int, int]) -> dict[int, int]:
    ordered = _ordered_tree_nodes(parent)
    swc_id = {local_id: index + 1 for index, local_id in enumerate(ordered)}
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float)
    lines = [
        "# Canonical saved ROI morphology for Ultraliser provenance",
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
    *,
    radius_scale: float = ULTRALISER_RADIUS_SCALE,
) -> dict[str, Any]:
    """Write H5 diameter rows; the official reader multiplies column four by 0.5."""

    if not math.isfinite(radius_scale) or radius_scale <= 0.0:
        raise ValueError("radius_scale must be finite and positive")
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    source_radii = np.asarray(roi.local_node_radius_um, dtype=float)
    feed_radii = source_radii * float(radius_scale)
    point_rows: list[list[float]] = []
    structure: list[tuple[int, int]] = []
    for section in sections:
        first = len(point_rows)
        for local_id in section:
            point_rows.append(
                [*map(float, positions[local_id]), 2.0 * float(feed_radii[local_id])]
            )
        structure.append((first, len(point_rows) - 1))
    connectivity = [
        (parent_index, child_index)
        for parent_index, parent_path in enumerate(sections)
        for child_index, child_path in enumerate(sections)
        if child_index != parent_index and child_path[0] == parent_path[-1]
    ]
    points = np.asarray(point_rows, dtype=np.float32)
    with h5py.File(path, "w") as stream:
        stream.create_dataset("points", data=points)
        stream.create_dataset("structure", data=np.asarray(structure, dtype=np.int64).reshape((-1, 2)))
        stream.create_dataset(
            "connectivity",
            data=np.asarray(connectivity, dtype=np.int64).reshape((-1, 2)),
        )
        stream.attrs["canonical_source"] = "roi_core.swc"
        stream.attrs["coordinate_unit"] = "um"
        stream.attrs["points_fourth_column"] = "diameter_um"
        stream.attrs["official_reader_operation"] = "radius_um = points[:,3] * 0.5"
        stream.attrs["radius_scale_for_ultraliser"] = float(radius_scale)
        stream.attrs["source_radius_modified"] = False
        stream.attrs["source_swc_modified"] = False
    expanded_positions = np.asarray(
        [positions[node] for section in sections for node in section], dtype=float
    )
    expanded_feed = np.asarray(
        [feed_radii[node] for section in sections for node in section], dtype=float
    )
    return {
        "adapter_format": "Ultraliser vascular H5",
        "executable_input": path.name,
        "points_fourth_column": "diameter_um",
        "official_reader_operation": "internal_radius_um = points[:,3] * 0.5",
        "radius_scale_for_ultraliser": float(radius_scale),
        "source_radius_modified": False,
        "source_swc_modified": False,
        "section_count": len(sections),
        "connectivity_count": len(connectivity),
        "point_row_count_with_section_endpoint_duplicates": len(points),
        "maximum_float32_coordinate_quantization_um": float(
            np.max(np.abs(points[:, :3].astype(float) - expanded_positions))
        ),
        "maximum_float32_radius_quantization_um": float(
            np.max(np.abs(0.5 * points[:, 3].astype(float) - expanded_feed))
        ),
    }


def _write_radius_feed_mapping(
    roi: ROIRecord,
    path: Path,
    radius_scale: float,
) -> Path:
    source = np.asarray(roi.local_node_radius_um, dtype=float)
    feed = source * float(radius_scale)
    if np.any(feed <= 0.0):
        raise UltraliserBackendError("H5 feed radii must all be positive")
    return write_csv(
        path,
        [
            {
                "local_node_id": int(node),
                "source_radius_um": float(source[node]),
                "radius_scale": float(radius_scale),
                "feed_radius_um": float(feed[node]),
                "h5_diameter_um": float(2.0 * feed[node]),
                "absolute_difference_um": float(abs(feed[node] - source[node])),
                "relative_difference": float((feed[node] - source[node]) / source[node]),
            }
            for node in range(roi.node_count)
        ],
    )


def export_roi_for_ultraliser(
    roi: ROIRecord,
    input_directory: Path,
    *,
    radius_scale_for_ultraliser: float = ULTRALISER_RADIUS_SCALE,
    h5_filename: str = "roi_core.h5",
) -> ROIUltraliserInput:
    """Export canonical SWC plus the radius-scaled executable H5 adapter."""

    _validate_roi_arrays(roi)
    graph = _roi_graph(roi)
    component_count = nx.number_connected_components(graph)
    cycle_rank = graph.number_of_edges() - graph.number_of_nodes() + component_count
    if component_count != 1:
        raise UltraliserBackendError(
            f"ROI must contain exactly one connected component, got {component_count}"
        )
    if cycle_rank > 0:
        raise UltraliserCycleConflict(
            f"ROI {roi.roi_id} has cycle_rank={cycle_rank}; cyclic edges are never removed"
        )
    input_directory.mkdir(parents=True, exist_ok=True)
    parent, serialization_direction_only = _source_parent_map(roi, graph)
    sections = _directed_sections(parent, graph)
    swc_path = input_directory / "roi_core.swc"
    h5_path = input_directory / h5_filename
    metadata_path = input_directory / "metadata.json"
    radius_mapping_path = input_directory / "radius_feed_mapping.csv"
    cut_mapping_path = input_directory / "cut_port_mapping.csv"
    swc_ids = _write_swc(roi, swc_path, parent)
    adapter = _write_h5_adapter(
        roi,
        h5_path,
        sections,
        radius_scale=radius_scale_for_ultraliser,
    )
    _write_radius_feed_mapping(roi, radius_mapping_path, radius_scale_for_ultraliser)
    write_csv(
        cut_mapping_path,
        [
            {
                "cut_port_id": port.cut_port_id,
                "local_node_id": int(port.local_node_id),
                "swc_node_id": swc_ids[int(port.local_node_id)],
                "global_edge_id": int(port.global_edge_id),
                "x_um": float(port.intersection_position_um[0]),
                "y_um": float(port.intersection_position_um[1]),
                "z_um": float(port.intersection_position_um[2]),
                "source_radius_um": float(port.radius_at_cut_um),
                "boundary_face": port.boundary_face,
            }
            for port in roi.cut_ports
        ],
    )
    topology = {
        "connected_component_count": int(component_count),
        "cycle_rank": int(cycle_rank),
        "source_parent_child_relation_used": not serialization_direction_only,
        "serialization_direction_only": bool(serialization_direction_only),
        "serialization_is_physiological_flow_direction": False,
        "root_local_node_id": int(next(node for node, value in parent.items() if value < 0)),
    }
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float)
    metadata = {
        "roi_id": roi.roi_id,
        "source_model_id": roi.source_model_id,
        "anchor_id": int(roi.anchor_id),
        "saved_roi_reused": True,
        "source_geometry_modified": False,
        "source_radius_modified": False,
        "source_swc_modified": False,
        "coordinate_unit": "um",
        "radius_unit": "um",
        "source_node_count": roi.node_count,
        "source_edge_count": roi.edge_count,
        "cut_port_count": len(roi.cut_ports),
        "radius_scale_for_ultraliser": float(radius_scale_for_ultraliser),
        "feed_radius_equation": "feed_radius_um = source_radius_um * radius_scale",
        "h5_diameter_equation": "h5_diameter_um = 2 * feed_radius_um",
        "source_radius_min_um": float(radii.min()),
        "source_radius_median_um": float(np.median(radii)),
        "source_radius_max_um": float(radii.max()),
        "bbox_min_um": positions.min(axis=0).tolist(),
        "bbox_max_um": positions.max(axis=0).tolist(),
        "canonical_swc_sha256": _sha256(swc_path),
        "topology": topology,
        "swc_node_id_by_local_node_id": {str(key): value for key, value in swc_ids.items()},
        "ultraliser_input_adapter": adapter,
    }
    write_json(metadata_path, metadata)
    return ROIUltraliserInput(
        swc_path=swc_path,
        h5_path=h5_path,
        metadata_path=metadata_path,
        radius_feed_mapping_path=radius_mapping_path,
        cut_port_mapping_path=cut_mapping_path,
        topology=topology,
        sections=sections,
    )


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
    raise UltraliserExecutableNotFound(
        f"no ultraVessMorpho2Mesh executable found under {root}"
    )


def _wsl_path(path: Path) -> str:
    resolved = Path(path).resolve()
    drive = resolved.drive.rstrip(":").lower()
    if not drive:
        return resolved.as_posix()
    return f"/mnt/{drive}/{resolved.as_posix()[2:].lstrip('/')}"


def _execution_command(executable: Path, arguments: list[str]) -> list[str]:
    if os.name == "nt" and executable.suffix.lower() != ".exe":
        converted = list(arguments)
        for index, value in enumerate(arguments):
            if index > 0 and arguments[index - 1] in {"--morphology", "--output-directory"}:
                converted[index] = _wsl_path(Path(value))
        return ["wsl", "-d", "Ubuntu", "--", _wsl_path(executable), *converted]
    return [str(executable), *arguments]


def build_ultraliser_command(
    executable: Path,
    morphology: Path,
    output_directory: Path,
    *,
    prefix: str,
    voxels_per_micron: float = 6.0,
    threads: int = 8,
    packing_algorithm: str = ULTRALISER_PACKING_ALGORITHM,
    voxelization_axis: str = ULTRALISER_VOXELIZATION_AXIS,
    isosurface_technique: str = ULTRALISER_ISOSURFACE_TECHNIQUE,
    solid_voxelization: bool = True,
    adaptive_optimization: bool = True,
    optimization_iterations: int = 5,
    smooth_iterations: int = 5,
    laplacian_iterations: int = 10,
    export_stl: bool = True,
) -> tuple[list[str], str]:
    if voxels_per_micron <= 0.0 or threads < 1:
        raise ValueError("voxels_per_micron and threads must be positive")
    if min(optimization_iterations, smooth_iterations, laplacian_iterations) < 0:
        raise ValueError("Ultraliser iteration counts cannot be negative")
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
    ]
    if solid_voxelization:
        arguments.append("--solid")
    arguments.extend(
        (
            "--voxelization-axis",
            voxelization_axis,
            "--packing-algorithm",
            packing_algorithm,
            "--isosurface-technique",
            isosurface_technique,
        )
    )
    if adaptive_optimization:
        arguments.append("--adaptive-optimization")
    arguments.extend(
        (
            "--optimization-iterations",
            str(int(optimization_iterations)),
            "--smooth-iterations",
            str(int(smooth_iterations)),
            "--laplacian-iterations",
            str(int(laplacian_iterations)),
        )
    )
    if export_stl:
        arguments.append("--export-stl-mesh")
    arguments.extend(("--stats", "--threads", str(int(threads))))
    command = _execution_command(executable, arguments)
    return command, shlex.join(command)


def invoke_ultraliser(
    executable: Path,
    morphology: Path,
    output_directory: Path,
    *,
    prefix: str,
    settings: UltraliserConfig,
) -> UltraliserInvocation:
    command, command_text = build_ultraliser_command(
        executable,
        morphology,
        output_directory,
        prefix=prefix,
        voxels_per_micron=settings.voxels_per_micron,
        threads=settings.threads,
        packing_algorithm=settings.packing_algorithm,
        voxelization_axis=settings.voxelization_axis,
        isosurface_technique=settings.isosurface_technique,
        solid_voxelization=settings.solid_voxelization,
        adaptive_optimization=settings.adaptive_optimization,
        optimization_iterations=settings.optimization_iterations,
        smooth_iterations=settings.smooth_iterations,
        laplacian_iterations=settings.laplacian_iterations,
        export_stl=settings.export_stl,
    )
    output_directory.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    with (output_directory / "stdout.log").open("wb") as stdout_stream, (
        output_directory / "stderr.log"
    ).open("wb") as stderr_stream:
        completed = subprocess.run(
            command,
            stdout=stdout_stream,
            stderr=stderr_stream,
            check=False,
        )
    runtime = time.perf_counter() - started
    if completed.returncode != 0:
        raise UltraliserRunFailed(
            f"official process returned {completed.returncode}; see {output_directory}"
        )
    return UltraliserInvocation(
        command=tuple(command),
        command_text=command_text,
        runtime_seconds=float(runtime),
        return_code=int(completed.returncode),
        output_directory=output_directory,
    )


def discover_watertight_stl(output_directory: Path) -> Path:
    candidates = sorted(output_directory.glob("meshes/*-watertight.stl"))
    if not candidates:
        raise UltraliserRunFailed(
            f"official watertight STL was not produced under {output_directory}"
        )
    return candidates[-1]


def _remove_successful_work_directory(layout: UltraliserLayout) -> None:
    run_root = layout.run_root.resolve()
    work = layout.work.resolve()
    if work.parent != run_root or work.name != ".ultraliser_work":
        raise UltraliserBackendError(f"refusing to remove unexpected work directory: {work}")
    shutil.rmtree(work)


def _preserve_source_configuration(input_directory: Path, source_path: Path) -> Path:
    source = Path(source_path).resolve()
    if not source.is_file():
        raise FileNotFoundError(f"model-generation source configuration not found: {source}")
    destination = Path(input_directory) / "source_swc_stl_model_generate.yaml"
    shutil.copy2(source, destination)
    return destination


def run_ultraliser_reconstruction(
    roi: ROIRecord,
    config: CFDLumenConfig,
    *,
    output_root: Path,
    run_id: str | None,
    ultraliser_root: Path,
    executable_path: Path | None = None,
    source_config_path: Path | None = None,
) -> dict[str, Any]:
    """Run the only supported reconstruction once; no smoke run or fallback exists."""

    config.validate()
    settings = config.ultraliser
    layout = create_ultraliser_layout(output_root, roi, run_id)
    source_config_copy: Path | None = None
    if source_config_path is not None:
        source_config_copy = _preserve_source_configuration(
            layout.input,
            source_config_path,
        )
    source_positions = np.asarray(roi.local_node_positions_um, dtype=float).copy()
    source_radii = np.asarray(roi.local_node_radius_um, dtype=float).copy()
    exported = export_roi_for_ultraliser(
        roi,
        layout.input,
        radius_scale_for_ultraliser=settings.radius_scale,
    )
    swc_hash_before = _sha256(exported.swc_path)
    executable = (
        Path(executable_path).resolve()
        if executable_path is not None
        else discover_ultraliser_executable(ultraliser_root)
    )
    invocation = invoke_ultraliser(
        executable,
        exported.h5_path,
        layout.work,
        prefix=f"anchor{int(roi.anchor_id):06d}",
        settings=settings,
    )
    raw_surface = discover_watertight_stl(layout.work)
    from .ultraliser_qc import finalize_ultraliser_outputs, write_reconstruction_report

    qc_result = finalize_ultraliser_outputs(
        roi,
        config,
        raw_surface=raw_surface,
        geometry_directory=layout.geometry,
        qc_directory=layout.qc,
    )
    swc_hash_after = _sha256(exported.swc_path)
    source_unchanged = bool(
        swc_hash_before == swc_hash_after
        and np.array_equal(source_positions, np.asarray(roi.local_node_positions_um))
        and np.array_equal(source_radii, np.asarray(roi.local_node_radius_um))
    )
    if not source_unchanged:
        raise UltraliserBackendError("source ROI or canonical SWC changed during reconstruction")
    status = (
        "PASS"
        if qc_result["surface_qc"]["status"] == "PASS"
        and qc_result["radius_fidelity"]["status"] == "PASS"
        else "FAIL"
    )
    summary = {
        "status": status,
        "roi_id": roi.roi_id,
        "anchor_id": int(roi.anchor_id),
        "run_root": str(layout.run_root),
        "source_geometry_modified": False,
        "source_radius_modified": False,
        "source_swc_modified": False,
        "source_swc_unchanged": True,
        "canonical_swc_sha256_before": swc_hash_before,
        "canonical_swc_sha256_after": swc_hash_after,
        "radius_scale": float(settings.radius_scale),
        "feed_radius_equation": "feed_radius_um = source_radius_um * radius_scale",
        "ultraliser_invocation_count": 1,
        "ultraliser_command": invocation.command_text,
        "ultraliser_runtime_seconds": invocation.runtime_seconds,
        "ultraliser_executable": str(executable),
        "source_configuration": str(source_config_copy) if source_config_copy else None,
        "ultraliser_settings": config.report()["ultraliser"],
        "source_qc": qc_result["source_qc"],
        "surface_qc": qc_result["surface_qc"],
        "radius_fidelity": {
            key: value
            for key, value in qc_result["radius_fidelity"].items()
            if key != "samples"
        },
        "outputs": {
            "surface_um_stl": str(layout.geometry / "lumen_surface_um.stl"),
            "surface_um_vtp": str(layout.geometry / "lumen_surface_um.vtp"),
            "surface_m_stl": str(layout.geometry / "lumen_surface_m.stl"),
            "surface_qc": str(layout.qc / "surface_qc.json"),
            "radius_fidelity": str(layout.qc / "radius_fidelity.json"),
            "source_configuration": str(source_config_copy) if source_config_copy else None,
        },
        "fallback_used": False,
        "smoke_run": False,
        "cfd_port_preparation_run": False,
        "volume_mesh_run": False,
        "cfd_run": False,
        "microbubble_simulation_run": False,
    }
    write_json(layout.qc / "run_summary.json", summary)
    write_reconstruction_report(layout.report / "reconstruction_report.md", summary)
    _remove_successful_work_directory(layout)
    return summary
