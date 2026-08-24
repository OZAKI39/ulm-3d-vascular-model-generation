"""Generate document figures explaining representative connected-ROI sampling.

The figures are derived from saved real preprocessing and sampling artifacts.
No vascular geometry, cluster assignment, or representative ROI is synthesized.

Example
-------
python doc_visualize.py
python doc_visualize.py --roi-id 6359 --dpi 220
"""

from __future__ import annotations

import argparse
import csv
import json
import logging
import re
from collections import Counter
from dataclasses import dataclass, replace
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
from mpl_toolkits.mplot3d.art3d import Line3DCollection, Poly3DCollection
from scipy.ndimage import label as label_connected_components
from skimage.measure import marching_cubes

from utils.rodent_vasculature.swc_io import SWCData, load_normalized_swc, load_swc
from utils.rodent_vasculature.tiff_io import load_normalized_volume, load_tiff_volume
from utils.sampling.pipeline import load_models_from_rodent_run
from utils.sampling.roi_boundary import clip_segment_to_box, point_in_box
from utils.sampling.roi_extraction import (
    EdgeSpatialIndex,
    _interior_anchor_candidates,
    generate_anchor_ids,
)
from utils.sampling.sampling_config import SamplingConfig
from utils.sampling.sampling_io import load_sampling_display_rois
from utils.sampling.sampling_types import GlobalVascularModel, ROIRecord


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "outputs"
DEFAULT_FIGURE_ROOT = DEFAULT_OUTPUT_ROOT / "documentation" / "representative_roi_sampling"
FEATURE_LABELS = {
    "r10": "$r_{10}$",
    "r25": "$r_{25}$",
    "r50": "$r_{50}$",
    "r75": "$r_{75}$",
    "r90": "$r_{90}$",
    "branch_count": "分支数",
    "bifurcation_count": "分叉数",
    "total_vessel_length_um": "总长度",
    "cycle_rank": "环路秩",
}
CLUSTER_COLORS = (
    "#2A9D8F",
    "#E76F51",
    "#E9C46A",
    "#457B9D",
    "#9B5DE5",
    "#F15BB5",
)


@dataclass(frozen=True, slots=True)
class RawClipSegment:
    points_um: np.ndarray
    global_edge_id: int
    component_id: int
    is_anchor_component: bool


@dataclass(frozen=True, slots=True)
class MaskSurface:
    """Display mesh extracted from a real binary segmentation mask."""

    vertices_xyz_um: np.ndarray
    faces: np.ndarray


@dataclass(frozen=True, slots=True)
class HierarchyBranch:
    """One real aggregated branch loaded from the directed-graph export."""

    branch_id: int
    upstream_node_id: int
    downstream_node_id: int
    source_node_ids: tuple[int, ...]
    depth: int
    parent_branch_ids: tuple[int, ...]
    daughter_branch_ids: tuple[int, ...]
    strahler_order: int
    horsfield_order: int
    downstream_terminal_count: int
    downstream_branch_count: int
    raw_length_um: float
    mean_radius_um: float


@dataclass(frozen=True, slots=True)
class AnchorStrategyResult:
    """One anchor mode evaluated on the same real interior critical-node pool."""

    mode: str
    anchor_ids: tuple[int, ...]
    closest_pair_ids: tuple[int, int]
    minimum_pairwise_distance_um: float
    median_nearest_neighbor_distance_um: float
    mean_nearest_neighbor_distance_um: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualize the real connected-ROI representative sampling process."
    )
    parser.add_argument(
        "--rodent-run",
        type=Path,
        default=None,
        help="Saved rodent_vasculature run. The latest compatible run is used by default.",
    )
    parser.add_argument(
        "--sampling-run",
        type=Path,
        default=None,
        help="Saved sampling run. The latest complete run is used by default.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_FIGURE_ROOT,
        help="Directory for figures, manifest, and log.",
    )
    parser.add_argument(
        "--roi-id",
        default=None,
        help="Full ROI id or anchor id. Defaults to the first selected representative.",
    )
    parser.add_argument("--dpi", type=int, default=200)
    parser.add_argument("--show", action="store_true", help="Also open the generated figures.")
    return parser.parse_args()


def _configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.sans-serif": ["Microsoft YaHei", "SimHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#4B5563",
            "axes.labelcolor": "#263238",
            "xtick.color": "#455A64",
            "ytick.color": "#455A64",
            "text.color": "#1F2933",
            "axes.titleweight": "semibold",
        }
    )


def _configure_logger(output_dir: Path) -> logging.Logger:
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("doc_visualize")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    stream = logging.StreamHandler()
    stream.setFormatter(formatter)
    file_handler = logging.FileHandler(output_dir / "visualization.log", encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(stream)
    logger.addHandler(file_handler)
    return logger


def _latest_sampling_run(root: Path) -> Path:
    candidates = [
        path
        for path in (root / "sampling").glob("*")
        if path.is_dir()
        and (path / "report" / "sampling_summary.json").is_file()
        and (path / "manifests" / "candidate_rois.csv").is_file()
    ]
    if not candidates:
        raise FileNotFoundError(f"No completed sampling runs found under {root / 'sampling'}")
    return max(candidates, key=lambda path: (path.stat().st_mtime, path.name))


def _compatible_rodent_run(root: Path, source_model_ids: set[str]) -> Path:
    candidates = sorted(
        (path for path in (root / "rodent_vasculature").glob("*") if path.is_dir()),
        key=lambda path: (path.stat().st_mtime, path.name),
        reverse=True,
    )
    for run_root in candidates:
        for manifest_path in (run_root / "samples").glob("*/preprocess_manifest.json"):
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            if str(manifest.get("record", {}).get("sample_id")) in source_model_ids:
                return run_root
    raise FileNotFoundError(
        "No rodent preprocessing run contains the source model used by the sampling run"
    )


def _read_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _simple_yaml_number(path: Path, key: str, default: float) -> float:
    if not path.is_file():
        return default
    match = re.search(
        rf"(?m)^{re.escape(key)}:\s*([-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*$",
        path.read_text(encoding="utf-8"),
    )
    return float(match.group(1)) if match else default


def _choose_roi(rois: list[ROIRecord], requested: str | None) -> ROIRecord:
    if requested:
        token = str(requested).strip()
        for roi in rois:
            if roi.roi_id == token or str(roi.anchor_id) == token:
                return roi
        raise ValueError(f"ROI {requested!r} is not present in the saved candidate library")
    representatives = [roi for roi in rois if roi.is_representative]
    if representatives:
        return min(
            representatives,
            key=lambda roi: roi.selection_rank if roi.selection_rank > 0 else 100000,
        )
    return rois[0]


def _select_model(models: list[GlobalVascularModel], roi: ROIRecord) -> GlobalVascularModel:
    for model in models:
        if model.source_model_id == roi.source_model_id:
            return model
    raise ValueError(f"Source model {roi.source_model_id!r} is absent from the rodent run")


def _load_original_input_data(
    rodent_run: Path,
    source_model_id: str,
) -> tuple[dict[str, object], np.ndarray, SWCData]:
    """Load the original mask and SWC paired with the selected real model."""

    for manifest_path in (rodent_run / "samples").glob("*/preprocess_manifest.json"):
        manifest = _read_json(manifest_path)
        record = manifest.get("record", {})
        if not isinstance(record, dict) or str(record.get("sample_id")) != source_model_id:
            continue
        spacing = tuple(float(value) for value in manifest["spacing_xyz_um"])  # type: ignore[index]
        shape = tuple(int(value) for value in manifest["image_metadata"]["shape_zyx"])  # type: ignore[index]
        mask_path = Path(str(record["mask_path"]))
        swc_path = Path(str(record["swc_path"]))
        mask = load_tiff_volume(mask_path)
        swc = load_swc(
            swc_path,
            spacing_xyz_um=spacing,  # type: ignore[arg-type]
            volume_shape_zyx=shape,  # type: ignore[arg-type]
        )
        return manifest, mask, swc
    raise FileNotFoundError(
        f"No preprocess manifest for source model {source_model_id!r} in {rodent_run}"
    )


def _parse_id_list(value: str | None) -> tuple[int, ...]:
    return tuple(int(token) for token in str(value or "").split(";") if token)


def _load_hierarchy_artifacts(
    rodent_run: Path,
    source_model_id: str,
) -> tuple[Path, list[HierarchyBranch], dict[int, str], list[tuple[int, int]], dict[str, object]]:
    """Load real point/branch hierarchy exports for one processed SWC model."""

    for sample_root in (rodent_run / "samples").iterdir():
        manifest_path = sample_root / "preprocess_manifest.json"
        if not manifest_path.is_file():
            continue
        manifest = _read_json(manifest_path)
        record = manifest.get("record", {})
        if not isinstance(record, dict) or str(record.get("sample_id")) != source_model_id:
            continue
        graph_root = sample_root / "graphs"
        with (graph_root / "directed_branches.csv").open(
            encoding="utf-8-sig", newline=""
        ) as stream:
            branches = [
                HierarchyBranch(
                    branch_id=int(row["branch_id"]),
                    upstream_node_id=int(row["upstream_node_id"]),
                    downstream_node_id=int(row["downstream_node_id"]),
                    source_node_ids=_parse_id_list(row["source_node_ids"]),
                    depth=int(row["depth"]),
                    parent_branch_ids=_parse_id_list(row["parent_branch_ids"]),
                    daughter_branch_ids=_parse_id_list(row["daughter_branch_ids"]),
                    strahler_order=int(row["strahler_order"]),
                    horsfield_order=int(row["horsfield_order"]),
                    downstream_terminal_count=int(row["downstream_terminal_count"]),
                    downstream_branch_count=int(row["downstream_branch_count"]),
                    raw_length_um=float(row["raw_length_um"]),
                    mean_radius_um=float(row["mean_radius_um"]),
                )
                for row in csv.DictReader(stream)
            ]
        with (graph_root / "directed_nodes.csv").open(
            encoding="utf-8-sig", newline=""
        ) as stream:
            node_roles = {
                int(row["node_id"]): str(row["role"])
                for row in csv.DictReader(stream)
            }
        with (graph_root / "source_parent_to_current_edges.csv").open(
            encoding="utf-8-sig", newline=""
        ) as stream:
            source_edges = [
                (int(row["parent_node_id"]), int(row["current_node_id"]))
                for row in csv.DictReader(stream)
            ]
        return (
            graph_root,
            branches,
            node_roles,
            source_edges,
            _read_json(graph_root / "directed_graph_summary.json"),
        )
    raise FileNotFoundError(
        f"No graph exports for source model {source_model_id!r} in {rodent_run}"
    )


def _hierarchy_graph(branches: list[HierarchyBranch]) -> nx.DiGraph:
    graph = nx.DiGraph()
    included_ids = {branch.branch_id for branch in branches}
    graph.add_nodes_from(included_ids)
    graph.add_edges_from(
        (branch.branch_id, daughter_id)
        for branch in branches
        for daughter_id in branch.daughter_branch_ids
        if daughter_id in included_ids
    )
    return graph


def _root_subtree_branch_ids(
    branches: list[HierarchyBranch], *, maximum_relative_depth: int = 3
) -> tuple[int, ...]:
    """Choose the largest real root subtree and keep a readable number of levels."""

    graph = _hierarchy_graph(branches)
    by_id = {branch.branch_id: branch for branch in branches}
    roots = [node for node, indegree in graph.in_degree() if indegree == 0]
    if not roots:
        raise ValueError("The exported branch hierarchy has no structural root")
    root = max(roots, key=lambda branch_id: by_id[branch_id].downstream_branch_count)
    distances = nx.single_source_shortest_path_length(
        graph, root, cutoff=maximum_relative_depth
    )
    return tuple(
        sorted(distances, key=lambda branch_id: (distances[branch_id], branch_id))
    )


def _branch_point_paths(
    branches: Iterable[HierarchyBranch], swc: SWCData
) -> dict[int, np.ndarray]:
    index_by_id = {int(node_id): index for index, node_id in enumerate(swc.node_ids)}
    paths: dict[int, np.ndarray] = {}
    for branch in branches:
        missing = [node_id for node_id in branch.source_node_ids if node_id not in index_by_id]
        if missing:
            raise ValueError(
                f"Branch {branch.branch_id} refers to missing SWC nodes: {missing[:5]}"
            )
        paths[branch.branch_id] = np.asarray(
            [swc.points_um[index_by_id[node_id]] for node_id in branch.source_node_ids],
            dtype=float,
        )
    return paths


def _hierarchy_audit(
    branches: list[HierarchyBranch], source_edges: list[tuple[int, int]]
) -> dict[str, object]:
    represented_edges = [
        edge
        for branch in branches
        for edge in zip(branch.source_node_ids[:-1], branch.source_node_ids[1:])
    ]
    source_counts = Counter(source_edges)
    represented_counts = Counter(represented_edges)
    missing = list((source_counts - represented_counts).elements())
    unexpected = list((represented_counts - source_counts).elements())
    duplicates = sum(max(count - 1, 0) for count in represented_counts.values())
    graph = _hierarchy_graph(branches)
    return {
        "represented_source_edge_count": len(represented_edges),
        "unique_represented_source_edge_count": len(represented_counts),
        "missing_source_edge_count": len(missing),
        "unexpected_or_reversed_edge_count": len(unexpected),
        "duplicate_representation_count": int(duplicates),
        "branch_relation_count": graph.number_of_edges(),
        "branch_hierarchy_is_dag": nx.is_directed_acyclic_graph(graph),
        "branch_root_count": sum(indegree == 0 for _, indegree in graph.in_degree()),
    }


def _mask_surface(
    mask_zyx: np.ndarray,
    spacing_xyz_um: tuple[float, float, float],
    *,
    step_size: int,
    crop_min_xyz_um: np.ndarray | None = None,
    crop_max_xyz_um: np.ndarray | None = None,
) -> MaskSurface:
    """Create a decimated display surface while retaining source-mask coordinates."""

    foreground = np.asarray(mask_zyx) > 0
    shape_xyz = np.asarray(foreground.shape[::-1], dtype=int)
    spacing_xyz = np.asarray(spacing_xyz_um, dtype=float)
    start_xyz = np.zeros(3, dtype=int)
    stop_xyz = shape_xyz.copy()
    if crop_min_xyz_um is not None and crop_max_xyz_um is not None:
        start_xyz = np.floor(np.asarray(crop_min_xyz_um) / spacing_xyz).astype(int)
        stop_xyz = np.ceil(np.asarray(crop_max_xyz_um) / spacing_xyz).astype(int) + 1
        start_xyz = np.clip(start_xyz, 0, shape_xyz - 1)
        stop_xyz = np.clip(stop_xyz, start_xyz + 1, shape_xyz)
    start_zyx = start_xyz[::-1]
    stop_zyx = stop_xyz[::-1]
    cropped = foreground[
        start_zyx[0] : stop_zyx[0],
        start_zyx[1] : stop_zyx[1],
        start_zyx[2] : stop_zyx[2],
    ]
    if not np.any(cropped):
        return MaskSurface(np.empty((0, 3), dtype=float), np.empty((0, 3), dtype=int))
    padded = np.pad(cropped.astype(np.uint8), 1, mode="constant")
    spacing_zyx = spacing_xyz[::-1]
    vertices_zyx, faces, _normals, _values = marching_cubes(
        padded,
        level=0.5,
        spacing=tuple(spacing_zyx),
        step_size=max(1, int(step_size)),
        allow_degenerate=False,
    )
    vertices_zyx += (start_zyx - 1) * spacing_zyx
    return MaskSurface(vertices_zyx[:, ::-1], np.asarray(faces, dtype=int))


def _swc_segments(swc: SWCData) -> np.ndarray:
    index_by_id = {int(node_id): index for index, node_id in enumerate(swc.node_ids)}
    segments = [
        (swc.points_um[index_by_id[int(parent_id)]], swc.points_um[index])
        for index, parent_id in enumerate(swc.parent_ids)
        if int(parent_id) != -1 and int(parent_id) in index_by_id
    ]
    return np.asarray(segments, dtype=float).reshape((-1, 2, 3))


def _swc_component_segments(swc: SWCData) -> tuple[np.ndarray, ...]:
    """Group real SWC edges by weakly connected component in source order."""

    graph = swc.directed_graph()
    components = [set(nodes) for nodes in nx.weakly_connected_components(graph)]
    component_by_node = {
        int(node_id): component_id
        for component_id, nodes in enumerate(components)
        for node_id in nodes
    }
    index_by_id = {int(node_id): index for index, node_id in enumerate(swc.node_ids)}
    grouped: list[list[tuple[np.ndarray, np.ndarray]]] = [[] for _ in components]
    for index, parent_id_raw in enumerate(swc.parent_ids):
        parent_id = int(parent_id_raw)
        node_id = int(swc.node_ids[index])
        if parent_id == -1 or parent_id not in index_by_id:
            continue
        component_id = component_by_node[node_id]
        grouped[component_id].append(
            (swc.points_um[index_by_id[parent_id]], swc.points_um[index])
        )
    return tuple(
        np.asarray(segments, dtype=float).reshape((-1, 2, 3))
        if segments
        else np.empty((0, 2, 3), dtype=float)
        for segments in grouped
    )


def _style_3d_axis(
    axis: plt.Axes,
    bounds_xyz_um: tuple[float, float, float, float, float, float],
    *,
    elev: float = 22.0,
    azim: float = -58.0,
) -> None:
    axis.set_xlim(bounds_xyz_um[0], bounds_xyz_um[1])
    axis.set_ylim(bounds_xyz_um[2], bounds_xyz_um[3])
    axis.set_zlim(bounds_xyz_um[4], bounds_xyz_um[5])
    axis.set_xlabel("x (μm)", labelpad=7)
    axis.set_ylabel("y (μm)", labelpad=7)
    axis.set_zlabel("z (μm)", labelpad=7)
    span = np.maximum(np.asarray(bounds_xyz_um[1::2]) - np.asarray(bounds_xyz_um[0::2]), 1.0)
    axis.set_box_aspect(span)
    axis.view_init(elev=elev, azim=azim)
    axis.tick_params(labelsize=8, pad=1)
    axis.grid(True, linewidth=0.4, alpha=0.45)


def _add_mask_surface_3d(
    axis: plt.Axes,
    surface: MaskSurface,
    *,
    color: str = "#90A4AE",
    alpha: float = 0.62,
) -> None:
    if not len(surface.faces):
        return
    mesh = Poly3DCollection(
        surface.vertices_xyz_um[surface.faces],
        facecolor=color,
        edgecolor="none",
        alpha=alpha,
        linewidth=0.0,
    )
    mesh.set_rasterized(True)
    axis.add_collection3d(mesh)


def _add_segments_3d(
    axis: plt.Axes,
    segments_xyz_um: np.ndarray,
    *,
    color: str,
    linewidth: float | np.ndarray,
    alpha: float,
) -> None:
    if not len(segments_xyz_um):
        return
    collection = Line3DCollection(
        np.asarray(segments_xyz_um, dtype=float),
        colors=color,
        linewidths=linewidth,
        alpha=alpha,
    )
    collection.set_rasterized(True)
    axis.add_collection3d(collection)


def _add_box_3d(
    axis: plt.Axes,
    low_xyz_um: np.ndarray,
    high_xyz_um: np.ndarray,
    *,
    color: str = "#D1495B",
    linewidth: float = 1.8,
    alpha: float = 1.0,
) -> None:
    low = np.asarray(low_xyz_um, dtype=float)
    high = np.asarray(high_xyz_um, dtype=float)
    corners = np.asarray(
        [
            (x, y, z)
            for x in (low[0], high[0])
            for y in (low[1], high[1])
            for z in (low[2], high[2])
        ],
        dtype=float,
    )
    edges = [
        (left, right)
        for left in range(8)
        for right in range(left + 1, 8)
        if np.count_nonzero(corners[left] != corners[right]) == 1
    ]
    _add_segments_3d(
        axis,
        np.asarray([(corners[left], corners[right]) for left, right in edges]),
        color=color,
        linewidth=linewidth,
        alpha=alpha,
    )


def _sampling_config_from_summary(
    summary: dict[str, object], sampling_run: Path
) -> SamplingConfig:
    config_path = sampling_run / "config" / "sampling_config.yaml"
    return SamplingConfig(
        output_root=DEFAULT_OUTPUT_ROOT,
        seed=int(summary["seed"]),
        anchor_mode=str(summary["anchor_mode"]),
        min_anchor_distance_um=float(summary["minimum_anchor_distance_um"]),
        max_candidate_anchors=int(
            _simple_yaml_number(config_path, "max_candidate_anchors", 80)
        ),
        roi_size_um=tuple(float(value) for value in summary["roi_size_um"]),  # type: ignore[arg-type]
        min_branch_count=int(_simple_yaml_number(config_path, "min_branch_count", 2)),
    )


def _anchor_strategy_results(
    model: GlobalVascularModel,
    config: SamplingConfig,
) -> tuple[tuple[int, ...], dict[str, tuple[int, ...]], dict[str, AnchorStrategyResult]]:
    """Re-evaluate all implemented anchor modes on the exact same real candidate pool."""

    candidate_indices = _interior_anchor_candidates(
        model,
        np.asarray(config.roi_size_um, dtype=float),
    )
    candidate_ids = tuple(int(model.node_ids[index]) for index in candidate_indices)
    root_ids = {
        int(node_id)
        for node_id, parent_id in zip(model.node_ids, model.parent_ids)
        if int(parent_id) == -1
    }
    candidate_set = set(candidate_ids)
    role_ids = {
        "root": tuple(sorted(candidate_set & root_ids)),
        "terminal": tuple(
            sorted(
                node_id
                for node_id in candidate_set
                if model.global_degree_by_node[node_id] == 1 and node_id not in root_ids
            )
        ),
        "bifurcation": tuple(
            sorted(
                node_id
                for node_id in candidate_set
                if model.global_degree_by_node[node_id] >= 3 and node_id not in root_ids
            )
        ),
    }
    classified_ids = set().union(*(set(ids) for ids in role_ids.values()))
    role_ids["other_critical"] = tuple(sorted(candidate_set - classified_ids))

    results: dict[str, AnchorStrategyResult] = {}
    for mode in ("random", "farthest_point", "poisson_disk"):
        anchor_ids = generate_anchor_ids(model, replace(config, anchor_mode=mode))
        positions = np.asarray(
            [model.node_positions_um[model.node_index_by_id[node_id]] for node_id in anchor_ids],
            dtype=float,
        )
        if len(positions) < 2:
            closest_pair_ids = (-1, -1)
            minimum_distance = float("nan")
            median_nearest = float("nan")
            mean_nearest = float("nan")
        else:
            distances = np.linalg.norm(
                positions[:, np.newaxis, :] - positions[np.newaxis, :, :], axis=2
            )
            np.fill_diagonal(distances, np.inf)
            left, right = np.unravel_index(int(np.argmin(distances)), distances.shape)
            nearest = np.min(distances, axis=1)
            closest_pair_ids = (int(anchor_ids[left]), int(anchor_ids[right]))
            minimum_distance = float(distances[left, right])
            median_nearest = float(np.median(nearest))
            mean_nearest = float(np.mean(nearest))
        results[mode] = AnchorStrategyResult(
            mode=mode,
            anchor_ids=anchor_ids,
            closest_pair_ids=closest_pair_ids,
            minimum_pairwise_distance_um=minimum_distance,
            median_nearest_neighbor_distance_um=median_nearest,
            mean_nearest_neighbor_distance_um=mean_nearest,
        )
    return candidate_ids, role_ids, results


def _anchor_sampling_process_trace(
    model: GlobalVascularModel,
    config: SamplingConfig,
    candidate_ids: tuple[int, ...],
    strategy_results: dict[str, AnchorStrategyResult],
) -> dict[str, object]:
    """Replay the real decision sequence instead of retaining only final anchors."""

    ordered_candidate_ids = np.asarray(sorted(map(int, candidate_ids)), dtype=np.int64)
    positions = np.asarray(
        [model.node_positions_um[model.node_index_by_id[int(node_id)]] for node_id in ordered_candidate_ids],
        dtype=float,
    )
    rng = np.random.default_rng(config.seed)
    random_order = rng.permutation(len(ordered_candidate_ids)).tolist()
    selected_local: list[int] = []
    random_events: list[dict[str, object]] = []
    for scan_index, candidate_local in enumerate(random_order, start=1):
        candidate_local = int(candidate_local)
        selected_before = tuple(int(ordered_candidate_ids[index]) for index in selected_local)
        if selected_local:
            distances = np.linalg.norm(
                positions[selected_local] - positions[candidate_local], axis=1
            )
            nearest_offset = int(np.argmin(distances))
            nearest_node_id: int | None = int(
                ordered_candidate_ids[selected_local[nearest_offset]]
            )
            nearest_distance_um: float | None = float(distances[nearest_offset])
            accepted = nearest_distance_um >= config.min_anchor_distance_um
        else:
            nearest_node_id = None
            nearest_distance_um = None
            accepted = True
        if accepted:
            selected_local.append(candidate_local)
        random_events.append(
            {
                "scan_index": scan_index,
                "candidate_id": int(ordered_candidate_ids[candidate_local]),
                "accepted": accepted,
                "nearest_selected_id": nearest_node_id,
                "nearest_selected_distance_um": nearest_distance_um,
                "selected_ids_before": selected_before,
                "selected_ids_after": tuple(
                    int(ordered_candidate_ids[index]) for index in selected_local
                ),
            }
        )
        if len(selected_local) >= config.max_candidate_anchors:
            break
    random_accepted = tuple(
        int(event["candidate_id"])
        for event in random_events
        if bool(event["accepted"])
    )
    if random_accepted != strategy_results["random"].anchor_ids:
        raise ValueError("Random anchor process replay disagrees with generate_anchor_ids")

    centroid = np.mean(positions, axis=0)
    first_local = int(np.argmin(np.linalg.norm(positions - centroid, axis=1)))
    selected_local = [first_local]
    farthest_events: list[dict[str, object]] = [
        {
            "selection_index": 1,
            "candidate_id": int(ordered_candidate_ids[first_local]),
            "criterion": "nearest_to_candidate_centroid",
            "criterion_distance_um": float(np.linalg.norm(positions[first_local] - centroid)),
            "nearest_selected_id": None,
            "selected_ids_before": (),
            "selected_ids_after": (int(ordered_candidate_ids[first_local]),),
        }
    ]
    minimum_distance = np.linalg.norm(positions - positions[first_local], axis=1)
    minimum_distance[first_local] = -np.inf
    while len(selected_local) < min(config.max_candidate_anchors, len(ordered_candidate_ids)):
        next_local = int(np.argmax(minimum_distance))
        next_distance = float(minimum_distance[next_local])
        if next_distance < config.min_anchor_distance_um:
            break
        selected_before_local = tuple(selected_local)
        distances_to_selected = np.linalg.norm(
            positions[list(selected_before_local)] - positions[next_local], axis=1
        )
        nearest_offset = int(np.argmin(distances_to_selected))
        nearest_local = selected_before_local[nearest_offset]
        selected_local.append(next_local)
        farthest_events.append(
            {
                "selection_index": len(selected_local),
                "candidate_id": int(ordered_candidate_ids[next_local]),
                "criterion": "maximum_distance_to_selected_set",
                "criterion_distance_um": next_distance,
                "nearest_selected_id": int(ordered_candidate_ids[nearest_local]),
                "selected_ids_before": tuple(
                    int(ordered_candidate_ids[index]) for index in selected_before_local
                ),
                "selected_ids_after": tuple(
                    int(ordered_candidate_ids[index]) for index in selected_local
                ),
            }
        )
        minimum_distance = np.minimum(
            minimum_distance,
            np.linalg.norm(positions - positions[next_local], axis=1),
        )
        minimum_distance[selected_local] = -np.inf
    farthest_selected = tuple(
        int(event["candidate_id"]) for event in farthest_events
    )
    if farthest_selected != strategy_results["farthest_point"].anchor_ids:
        raise ValueError("Farthest-point process replay disagrees with generate_anchor_ids")
    return {
        "candidate_centroid_xyz_um": tuple(map(float, centroid)),
        "random_events": tuple(random_events),
        "farthest_events": tuple(farthest_events),
    }


def _raw_clipped_segments(
    model: GlobalVascularModel,
    roi: ROIRecord,
) -> tuple[list[RawClipSegment], int]:
    """Reproduce the exact pre-connectivity clipping for one saved ROI."""

    bbox_min = np.asarray(roi.bbox_min_um, dtype=float)
    bbox_max = np.asarray(roi.bbox_max_um, dtype=float)
    spatial_index = EdgeSpatialIndex(model)
    raw_edges: list[tuple[int, tuple[object, ...], tuple[object, ...], np.ndarray]] = []
    graph = nx.Graph()
    for edge_index in spatial_index.query(bbox_min, bbox_max).tolist():
        edge = model.edges[int(edge_index)]
        clipped = clip_segment_to_box(
            edge.upstream_position_um,
            edge.downstream_position_um,
            bbox_min,
            bbox_max,
        )
        if clipped is None:
            continue
        start_global = clipped.start_t <= 1.0e-10 and point_in_box(
            edge.upstream_position_um, bbox_min, bbox_max
        )
        end_global = clipped.end_t >= 1.0 - 1.0e-10 and point_in_box(
            edge.downstream_position_um, bbox_min, bbox_max
        )
        start_key: tuple[object, ...] = (
            ("g", edge.upstream_node_id)
            if start_global
            else ("c", edge.edge_id, "entry")
        )
        end_key: tuple[object, ...] = (
            ("g", edge.downstream_node_id)
            if end_global
            else ("c", edge.edge_id, "exit")
        )
        points = np.asarray((clipped.start_um, clipped.end_um), dtype=float)
        raw_edges.append((edge.edge_id, start_key, end_key, points))
        graph.add_edge(start_key, end_key)
    components = list(nx.connected_components(graph))
    anchor_key = ("g", int(roi.anchor_id))
    anchor_component = next(
        (index for index, component in enumerate(components) if anchor_key in component),
        -1,
    )
    component_by_node = {
        node: component_id
        for component_id, component in enumerate(components)
        for node in component
    }
    segments = [
        RawClipSegment(
            points_um=points,
            global_edge_id=edge_id,
            component_id=component_by_node[start_key],
            is_anchor_component=component_by_node[start_key] == anchor_component,
        )
        for edge_id, start_key, _end_key, points in raw_edges
    ]
    return segments, len(components)


def _global_edge_array(model: GlobalVascularModel) -> np.ndarray:
    return np.asarray(
        [
            (edge.upstream_position_um, edge.downstream_position_um)
            for edge in model.edges
        ],
        dtype=float,
    )


def _add_line_collection(
    axis: plt.Axes,
    segments_xyz: np.ndarray | Iterable[np.ndarray],
    *,
    dims: tuple[int, int] = (0, 2),
    color: str | list[str] = "#78909C",
    linewidth: float | np.ndarray = 0.7,
    alpha: float = 1.0,
    zorder: int = 1,
) -> None:
    segments = np.asarray(list(segments_xyz), dtype=float)
    if not len(segments):
        return
    projected = segments[:, :, list(dims)]
    axis.add_collection(
        LineCollection(
            projected,
            colors=color,
            linewidths=linewidth,
            alpha=alpha,
            zorder=zorder,
        )
    )


def _draw_projected_box(
    axis: plt.Axes,
    roi: ROIRecord,
    *,
    dims: tuple[int, int] = (0, 2),
    color: str = "#D1495B",
    linewidth: float = 1.8,
    linestyle: str = "-",
) -> None:
    low = np.asarray(roi.bbox_min_um)
    high = np.asarray(roi.bbox_max_um)
    axis.add_patch(
        Rectangle(
            (low[dims[0]], low[dims[1]]),
            high[dims[0]] - low[dims[0]],
            high[dims[1]] - low[dims[1]],
            fill=False,
            edgecolor=color,
            linewidth=linewidth,
            linestyle=linestyle,
            zorder=8,
        )
    )


def _style_spatial_axis(axis: plt.Axes, title: str) -> None:
    axis.set_title(title, loc="left", fontsize=12, pad=8)
    axis.set_xlabel("x (μm)")
    axis.set_ylabel("z (μm)")
    axis.set_aspect("equal", adjustable="box")
    axis.grid(True, color="#CFD8DC", linewidth=0.45, alpha=0.65)


def _save_figure(figure: plt.Figure, path: Path, dpi: int) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(path, dpi=dpi, bbox_inches="tight", facecolor="white")
    return path


def _figure_input_concepts_3d(
    mask_zyx: np.ndarray,
    mask_surface: MaskSurface,
    swc: SWCData,
    swc_segments: np.ndarray,
    model: GlobalVascularModel,
    roi: ROIRecord,
    spacing_xyz_um: tuple[float, float, float],
    output_path: Path,
    dpi: int,
) -> Path:
    """Show the three opening concepts as real three-dimensional data."""

    shape_xyz = np.asarray(mask_zyx.shape[::-1], dtype=float)
    maximum = np.maximum(shape_xyz - 1.0, 0.0) * np.asarray(spacing_xyz_um)
    global_bounds = (0.0, maximum[0], 0.0, maximum[1], 0.0, maximum[2])
    figure = plt.figure(figsize=(18, 7.2), constrained_layout=True)
    figure.suptitle(
        "输入概念的三维实数据示例：分割掩膜、SWC 中心线与感兴趣区域",
        fontsize=18,
        fontweight="semibold",
    )

    axis = figure.add_subplot(1, 3, 1, projection="3d")
    _add_mask_surface_3d(axis, mask_surface, color="#78909C", alpha=0.72)
    _style_3d_axis(axis, global_bounds)
    axis.set_title(
        "A  分割掩膜（原始输入）\n"
        f"{np.count_nonzero(mask_zyx):,} 个前景体素：表示被标记为血管的三维占据区域",
        loc="left",
        fontsize=11.5,
    )

    axis = figure.add_subplot(1, 3, 2, projection="3d")
    _add_segments_3d(
        axis,
        swc_segments,
        color="#007C91",
        linewidth=0.65,
        alpha=0.78,
    )
    root_index = {int(node_id): index for index, node_id in enumerate(swc.node_ids)}
    root_points = np.asarray(
        [swc.points_um[root_index[node_id]] for node_id in swc.root_ids if node_id in root_index]
    )
    if len(root_points):
        axis.scatter(
            root_points[:, 0],
            root_points[:, 1],
            root_points[:, 2],
            s=16,
            marker="o",
            color="#E76F51",
            depthshade=False,
            label="SWC 结构根",
        )
        axis.legend(loc="lower left", fontsize=8)
    _style_3d_axis(axis, global_bounds)
    axis.set_title(
        "B  SWC 中心线（原始输入）\n"
        f"{swc.node_count:,} 个节点 / {swc.edge_count:,} 条父子连接：表示骨架与拓扑关系",
        loc="left",
        fontsize=11.5,
    )

    axis = figure.add_subplot(1, 3, 3, projection="3d")
    roi_widths = np.clip(np.mean(roi.local_edge_radius_um, axis=1) * 0.75, 0.8, 4.0)
    _add_segments_3d(
        axis,
        roi.local_edge_points_um,
        color="#007C91",
        linewidth=roi_widths,
        alpha=0.96,
    )
    _add_box_3d(axis, roi.bbox_min_um, roi.bbox_max_um, color="#D1495B", linewidth=1.8)
    anchor = model.node_positions_um[model.node_index_by_id[int(roi.anchor_id)]]
    axis.scatter(
        [anchor[0]], [anchor[1]], [anchor[2]],
        marker="*", s=125, color="#F4A261", edgecolor="#7F5539", linewidth=0.8,
        depthshade=False, label="锚点",
    )
    if roi.cut_ports:
        cut_points = np.asarray([port.intersection_position_um for port in roi.cut_ports])
        axis.scatter(
            cut_points[:, 0], cut_points[:, 1], cut_points[:, 2],
            marker="s", s=28, color="#D1495B", depthshade=False, label="CUT_PORT",
        )
    if roi.true_terminal_local_ids:
        terminals = roi.local_node_positions_um[list(roi.true_terminal_local_ids)]
        axis.scatter(
            terminals[:, 0], terminals[:, 1], terminals[:, 2],
            marker="o", s=28, color="#2A9D8F", depthshade=False, label="源图终端",
        )
    roi_bounds = (
        roi.bbox_min_um[0], roi.bbox_max_um[0],
        roi.bbox_min_um[1], roi.bbox_max_um[1],
        roi.bbox_min_um[2], roi.bbox_max_um[2],
    )
    _style_3d_axis(axis, roi_bounds, elev=24.0, azim=-52.0)
    axis.legend(loc="lower left", fontsize=8)
    axis.set_title(
        "C  感兴趣区域 ROI（真实采样结果）\n"
        f"anchor {roi.anchor_id}；{roi.node_count} 节点 / {roi.branch_count} 分支：表示待分析的局部连通血管",
        loc="left",
        fontsize=11.5,
    )
    figure.text(
        0.5,
        0.015,
        "三图来自同一真实样本；数组按 (z,y,x) 存储，图中统一转换为 (x,y,z) 微米坐标。掩膜表面仅为显示而降采样，不改变源体素数据。",
        ha="center",
        fontsize=10,
        color="#455A64",
    )
    return _save_figure(figure, output_path, dpi)


def _figure_input_relationship_3d(
    crop_surface: MaskSurface,
    swc_segments: np.ndarray,
    model: GlobalVascularModel,
    roi: ROIRecord,
    display_min_xyz_um: np.ndarray,
    display_max_xyz_um: np.ndarray,
    output_path: Path,
    dpi: int,
) -> Path:
    """Overlay the actual mask, SWC skeleton, and retained ROI in one local 3-D view."""

    segment_min = np.min(swc_segments, axis=1)
    segment_max = np.max(swc_segments, axis=1)
    nearby = swc_segments[
        np.all(segment_max >= display_min_xyz_um, axis=1)
        & np.all(segment_min <= display_max_xyz_um, axis=1)
    ]
    figure = plt.figure(figsize=(11.5, 9.5), constrained_layout=True)
    axis = figure.add_subplot(1, 1, 1, projection="3d")
    _add_mask_surface_3d(axis, crop_surface, color="#B0BEC5", alpha=0.20)
    _add_segments_3d(
        axis,
        nearby,
        color="#546E7A",
        linewidth=0.65,
        alpha=0.48,
    )
    roi_widths = np.clip(np.mean(roi.local_edge_radius_um, axis=1) * 0.95, 1.2, 5.0)
    _add_segments_3d(
        axis,
        roi.local_edge_points_um,
        color="#00A6B4",
        linewidth=roi_widths,
        alpha=1.0,
    )
    _add_box_3d(axis, roi.bbox_min_um, roi.bbox_max_um, color="#D1495B", linewidth=2.2)
    anchor = model.node_positions_um[model.node_index_by_id[int(roi.anchor_id)]]
    axis.scatter(
        [anchor[0]], [anchor[1]], [anchor[2]],
        marker="*", s=155, color="#F4A261", edgecolor="#7F5539", linewidth=0.9,
        depthshade=False,
    )
    if roi.cut_ports:
        cut_points = np.asarray([port.intersection_position_um for port in roi.cut_ports])
        axis.scatter(
            cut_points[:, 0], cut_points[:, 1], cut_points[:, 2],
            marker="s", s=38, color="#D1495B", depthshade=False,
        )
    display_bounds = (
        display_min_xyz_um[0], display_max_xyz_um[0],
        display_min_xyz_um[1], display_max_xyz_um[1],
        display_min_xyz_um[2], display_max_xyz_um[2],
    )
    _style_3d_axis(axis, display_bounds, elev=25.0, azim=-48.0)
    axis.set_title(
        "同一空间中的三类信息：掩膜体积 → SWC 骨架 → 真实连通 ROI\n"
        f"示例 ROI：anchor {roi.anchor_id}，红色立方体为 {roi.bbox_size_um[0]:.0f} × "
        f"{roi.bbox_size_um[1]:.0f} × {roi.bbox_size_um[2]:.0f} μm³ 的裁切边界",
        loc="left",
        fontsize=14,
        pad=14,
    )
    axis.legend(
        handles=[
            Line2D([0], [0], color="#B0BEC5", linewidth=7, alpha=0.35, label="分割掩膜表面"),
            Line2D([0], [0], color="#546E7A", linewidth=1.2, label="原始 SWC 中心线"),
            Line2D([0], [0], color="#00A6B4", linewidth=3.2, label="保留的连通 ROI"),
            Line2D([0], [0], color="#D1495B", linewidth=2.2, label="ROI 空间边界"),
            Line2D([0], [0], marker="*", color="none", markerfacecolor="#F4A261", markeredgecolor="#7F5539", markersize=11, label="锚点"),
            Line2D([0], [0], marker="s", color="none", markerfacecolor="#D1495B", markersize=7, label="CUT_PORT"),
        ],
        loc="lower left",
        fontsize=9,
        frameon=True,
    )
    return _save_figure(figure, output_path, dpi)


def _figure_connectivity_preprocessing_3d(
    original_mask: np.ndarray,
    mask_surface: MaskSurface,
    reference_swc: SWCData,
    analysis_swc: SWCData,
    component_segments: tuple[np.ndarray, ...],
    component_records: list[dict[str, str]],
    preprocessing: dict[str, object],
    spacing_xyz_um: tuple[float, float, float],
    output_path: Path,
    dpi: int,
) -> Path:
    """Render SWC reference preservation and unchanged analysis-component selection."""

    shape_xyz = np.asarray(original_mask.shape[::-1], dtype=float)
    maximum = np.maximum(shape_xyz - 1.0, 0.0) * np.asarray(spacing_xyz_um)
    bounds = (0.0, maximum[0], 0.0, maximum[1], 0.0, maximum[2])
    selected_component_id = int(preprocessing["selected_component_id"])
    reference_only_ids = [
        int(record["swc_component_id"])
        for record in component_records
        if record["decision"] == "REFERENCE_ONLY"
    ]
    reference_only_segments = (
        np.concatenate([component_segments[index] for index in reference_only_ids], axis=0)
        if reference_only_ids
        else np.empty((0, 2, 3), dtype=float)
    )
    analysis_segments = _swc_segments(analysis_swc)
    reference_segments = _swc_segments(reference_swc)
    mask_qc = dict(preprocessing["mask_qc"])  # type: ignore[arg-type]

    figure = plt.figure(figsize=(18, 14.5), constrained_layout=True)
    figure.suptitle(
        "4.3.3 主血管网络连通性预处理：完整 SWC 参考—分量分析—单网络选择",
        fontsize=18,
        fontweight="semibold",
    )

    axis = figure.add_subplot(2, 3, 1, projection="3d")
    _add_segments_3d(axis, reference_segments, color="#546E7A", linewidth=0.58, alpha=0.76)
    _style_3d_axis(axis, bounds)
    axis.set_title(
        "A  完整人工修订 SWC（reference_swc）\n"
        f"{reference_swc.component_count} 个分量；{reference_swc.node_count:,} 节点 / "
        f"{reference_swc.edge_count:,} 边，全部保留",
        loc="left",
        fontsize=11,
    )

    axis = figure.add_subplot(2, 3, 2, projection="3d")
    for component_id, segments in enumerate(component_segments):
        _add_segments_3d(
            axis,
            segments,
            color=CLUSTER_COLORS[component_id % len(CLUSTER_COLORS)],
            linewidth=0.72,
            alpha=0.84,
        )
    _style_3d_axis(axis, bounds)
    axis.set_title(
        "B  仅按 SWC 的显式父子边识别分量\n"
        "不同颜色代表不同弱连通分量；不依据 Mask 合并或删除",
        loc="left",
        fontsize=11,
    )

    axis = figure.add_subplot(2, 3, 3, projection="3d")
    _add_segments_3d(axis, reference_only_segments, color="#E76F51", linewidth=0.62, alpha=0.46)
    _add_segments_3d(axis, analysis_segments, color="#007C91", linewidth=1.02, alpha=0.96)
    _style_3d_axis(axis, bounds)
    axis.set_title(
        "C  选择当前单网络分析对象\n"
        f"青色 C{selected_component_id} 为最长中心线分量；橙色分量仍保存在 reference_swc",
        loc="left",
        fontsize=11,
    )

    axis = figure.add_subplot(2, 3, 4, projection="3d")
    _add_segments_3d(axis, analysis_segments, color="#007C91", linewidth=0.82, alpha=0.90)
    _style_3d_axis(axis, bounds)
    axis.set_title(
        "D  analysis_swc：进入建图的原始分量\n"
        f"{analysis_swc.node_count:,} 节点 / {analysis_swc.edge_count:,} 边；新增节点=0，新增边=0",
        loc="left",
        fontsize=11,
    )

    axis = figure.add_subplot(2, 3, 5, projection="3d")
    _add_mask_surface_3d(axis, mask_surface, color="#A7C9CE", alpha=0.17)
    _add_segments_3d(axis, reference_segments, color="#546E7A", linewidth=0.52, alpha=0.72)
    _style_3d_axis(axis, bounds)
    axis.set_title(
        "E  可选原始 Mask（本样本存在）仅作为空间配准 QC\n"
        f"Mask 保持 {int(mask_qc['component_count_26'])} 个 26 邻域块；reference_swc 节点支持率 "
        f"{float(mask_qc['reference_swc_node_support_fraction']):.2%}",
        loc="left",
        fontsize=11,
    )

    axis = figure.add_subplot(2, 3, 6, projection="3d")
    _add_mask_surface_3d(axis, mask_surface, color="#A7C9CE", alpha=0.13)
    _add_segments_3d(axis, analysis_segments, color="#007C91", linewidth=0.82, alpha=0.94)
    _style_3d_axis(axis, bounds)
    axis.set_title(
        "F  下游数据角色明确分离\n"
        "青色 analysis_swc 用于 Graph/ROI；可选淡灰 Mask 只用于界面背景与 QC",
        loc="left",
        fontsize=11,
    )
    return _save_figure(figure, output_path, dpi)


def _figure_connectivity_preprocessing_audit(
    component_records: list[dict[str, str]],
    preprocessing: dict[str, object],
    output_path: Path,
    dpi: int,
) -> Path:
    """Explain SWC-only selection while preserving the complete reference forest."""

    figure, axes = plt.subplots(2, 2, figsize=(16, 11.5), constrained_layout=True)
    figure.suptitle(
        "SWC 中心化连通性预处理的定量依据与验收证据",
        fontsize=18,
        fontweight="semibold",
    )

    ordered = sorted(component_records, key=lambda record: int(record["length_rank"]))
    component_ids = np.asarray([int(record["swc_component_id"]) for record in ordered])
    lengths = np.asarray([float(record["total_length_um"]) for record in ordered])
    nodes = np.asarray([int(record["node_count"]) for record in ordered])
    diagonals = np.asarray([float(record["bbox_diagonal_um"]) for record in ordered])
    selected_component_id = int(preprocessing["selected_component_id"])
    selected_index = int(np.flatnonzero(component_ids == selected_component_id)[0])
    mask_qc = dict(preprocessing["mask_qc"])  # type: ignore[arg-type]

    axis = axes[0, 0]
    colors = ["#007C91" if index == selected_index else "#B0BEC5" for index in range(len(ordered))]
    axis.bar(np.arange(1, len(ordered) + 1), lengths, color=colors, width=0.82)
    axis.set_yscale("log")
    axis.set_xlabel("按中心线总长度从大到小排列的 SWC 分量")
    axis.set_ylabel("中心线总长度（μm，对数尺度）")
    axis.set_title(
        "A  analysis_swc 的预定义选择规则\n最长中心线分量用于当前单网络流程；其他分量仍完整保存",
        loc="left",
    )
    axis.annotate(
        f"C{selected_component_id}：{lengths[selected_index]:,.1f} μm",
        (selected_index + 1, lengths[selected_index]),
        xytext=(14, -12),
        textcoords="offset points",
        fontsize=9,
        color="#007C91",
        arrowprops={"arrowstyle": "->", "color": "#007C91"},
    )
    axis.grid(True, axis="y", which="both", color="#CFD8DC", linewidth=0.45)

    axis = axes[0, 1]
    sizes = 30.0 + 170.0 * diagonals / max(float(np.max(diagonals)), 1.0)
    axis.scatter(
        lengths, nodes, s=sizes, color="#B0BEC5", alpha=0.70,
        edgecolor="white", linewidth=0.5, label="REFERENCE_ONLY 或候选分量",
    )
    axis.scatter(
        lengths[selected_index], nodes[selected_index], s=245, marker="*", color="#007C91",
        edgecolor="#004D5A", linewidth=0.9, label=f"ANALYSIS_NETWORK：C{selected_component_id}", zorder=8,
    )
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlabel("中心线总长度（μm，对数尺度）")
    axis.set_ylabel("节点数（对数尺度）")
    axis.set_title(
        "B  选择依据的交叉检查\n横轴为长度，纵轴为节点数，气泡大小为空间跨度",
        loc="left",
    )
    axis.grid(True, which="both", color="#CFD8DC", linewidth=0.45)
    axis.legend(loc="lower right", fontsize=9)

    axis = axes[1, 0]
    labels = ["节点", "边", "中心线长度"]
    reference_values = np.asarray([float(preprocessing["reference_node_count"]), float(preprocessing["reference_edge_count"]), float(preprocessing["reference_total_length_um"])])
    analysis_values = np.asarray([float(preprocessing["analysis_node_count"]), float(preprocessing["analysis_edge_count"]), float(preprocessing["analysis_total_length_um"])])
    fractions = analysis_values / reference_values
    bars = axis.barh(labels, fractions, color="#007C91", alpha=0.88)
    axis.set_xlim(0, 1.05)
    axis.set_xlabel("analysis_swc / reference_swc")
    for bar, fraction, value, total in zip(bars, fractions, analysis_values, reference_values):
        axis.text(min(float(fraction) + 0.02, 0.92), bar.get_y() + bar.get_height() / 2, f"{fraction:.1%}  ({value:,.1f}/{total:,.1f})", va="center", fontsize=9)
    axis.set_title(
        "C  分析对象是参考数据的可追溯子集\n其余 42 个分量没有删除，仍保存在 reference_swc",
        loc="left",
    )
    axis.grid(True, axis="x", color="#CFD8DC", linewidth=0.45)

    axis = axes[1, 1]
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    boxes = [
        (0.02, 0.55, 0.20, 0.28, "reference_swc", f"完整保存\n{int(preprocessing['reference_component_count'])} 个分量"),
        (0.275, 0.55, 0.20, 0.28, "SWC 合法性验证", "编号 / 父节点 / 坐标\n半径 / 环路 / 根节点"),
        (0.53, 0.55, 0.20, 0.28, "分量分析与选择", f"最长中心线规则\n选择 C{selected_component_id}"),
        (0.785, 0.55, 0.20, 0.28, "analysis_swc", "原始单分量\n进入 Graph 与 ROI"),
    ]
    for x, y0, width, height, heading, detail in boxes:
        axis.add_patch(
            FancyBboxPatch(
                (x, y0), width, height,
                boxstyle="round,pad=0.012,rounding_size=0.015",
                facecolor="#E8F1F5" if x < 0.75 else "#E4F4EA",
                edgecolor="#607D8B",
                linewidth=1.2,
            )
        )
        axis.text(x + width / 2, y0 + 0.19, heading, ha="center", va="center", fontsize=10.5, fontweight="semibold")
        axis.text(x + width / 2, y0 + 0.09, detail, ha="center", va="center", fontsize=9.0)
    for start, end in ((0.22, 0.275), (0.475, 0.53), (0.73, 0.785)):
        axis.add_patch(
            FancyArrowPatch((start, 0.69), (end, 0.69), arrowstyle="-|>", mutation_scale=14, linewidth=1.4, color="#455A64")
        )
    axis.text(
        0.50,
        0.27,
        f"原始 Mask：{int(mask_qc['component_count_26'])} 个 26 邻域块，保持不变；"
        f"reference_swc 节点配准支持率 {float(mask_qc['reference_swc_node_support_fraction']):.2%}\n"
        "可选 Mask / raw image → 仅用于 QC 与可视化，不参与分量选择、删边或补边\n"
        f"拓扑改动：新增节点 {int(preprocessing['new_node_count'])} / 新增边 "
        f"{int(preprocessing['new_edge_count'])} / 父子关系变化 "
        f"{int(preprocessing['parent_relation_change_count'])}",
        ha="center",
        va="center",
        fontsize=10.5,
        color="#263238",
    )
    axis.set_title(
        "D  v2 数据职责与严格验收\n人工修订 SWC 决定拓扑；可选 Mask 不驱动任何拓扑修改",
        loc="left",
    )
    return _save_figure(figure, output_path, dpi)


def _path_segments(paths: Iterable[np.ndarray]) -> np.ndarray:
    segments = [
        np.stack((path[:-1], path[1:]), axis=1)
        for path in paths
        if len(path) >= 2
    ]
    return np.concatenate(segments, axis=0) if segments else np.empty((0, 2, 3))


def _paths_bounds(paths: Iterable[np.ndarray], padding_fraction: float = 0.08) -> tuple[float, ...]:
    usable = [np.asarray(path, dtype=float) for path in paths if len(path)]
    points = np.concatenate(usable, axis=0)
    lower = np.min(points, axis=0)
    upper = np.max(points, axis=0)
    padding = np.maximum((upper - lower) * padding_fraction, 1.0)
    lower -= padding
    upper += padding
    return (
        float(lower[0]), float(upper[0]),
        float(lower[1]), float(upper[1]),
        float(lower[2]), float(upper[2]),
    )


def _add_branch_direction_arrows(
    axis: plt.Axes,
    paths: Iterable[np.ndarray],
    *,
    color: str = "#F28E2B",
) -> int:
    count = 0
    for path in paths:
        if len(path) < 2:
            continue
        lengths = np.linalg.norm(np.diff(path, axis=0), axis=1)
        usable = np.flatnonzero(lengths > 1.0e-9)
        if not len(usable):
            continue
        segment_index = int(usable[len(usable) // 2])
        start = path[segment_index]
        vector = path[segment_index + 1] - start
        axis.quiver(
            start[0], start[1], start[2],
            vector[0], vector[1], vector[2],
            color=color,
            linewidth=1.25,
            arrow_length_ratio=0.42,
            length=5.5,
            normalize=True,
        )
        count += 1
    return count


def _figure_hierarchical_graph_construction_3d(
    analysis_swc: SWCData,
    branches: list[HierarchyBranch],
    node_roles: dict[int, str],
    graph_summary: dict[str, object],
    output_path: Path,
    dpi: int,
) -> tuple[Path, tuple[int, ...]]:
    """Show a real SWC root subtree changing from point samples into a branch hierarchy."""

    focus_ids = _root_subtree_branch_ids(branches, maximum_relative_depth=3)
    by_id = {branch.branch_id: branch for branch in branches}
    focus_branches = [by_id[branch_id] for branch_id in focus_ids]
    paths = _branch_point_paths(focus_branches, analysis_swc)
    ordered_paths = [paths[branch_id] for branch_id in focus_ids]
    bounds = _paths_bounds(ordered_paths)
    raw_segments = _path_segments(ordered_paths)
    source_node_ids = tuple(
        dict.fromkeys(
            node_id for branch in focus_branches for node_id in branch.source_node_ids
        )
    )
    index_by_id = {int(node_id): index for index, node_id in enumerate(analysis_swc.node_ids)}
    source_points = np.asarray(
        [analysis_swc.points_um[index_by_id[node_id]] for node_id in source_node_ids]
    )
    focus_parent_edge_count = sum(len(branch.source_node_ids) - 1 for branch in focus_branches)
    branch_colors = {
        branch_id: CLUSTER_COLORS[index % len(CLUSTER_COLORS)]
        for index, branch_id in enumerate(focus_ids)
    }

    figure = plt.figure(figsize=(18, 12.5), constrained_layout=True)
    figure.suptitle(
        "4.3.4 血管层级图构建：点级父子边—关键节点—聚合分支—分支层级",
        fontsize=18,
        fontweight="semibold",
    )

    axis = figure.add_subplot(2, 2, 1, projection="3d")
    _add_segments_3d(axis, raw_segments, color="#90A4AE", linewidth=0.82, alpha=0.84)
    axis.scatter(
        source_points[:, 0], source_points[:, 1], source_points[:, 2],
        s=8, color="#37474F", alpha=0.90, depthshade=False,
    )
    _style_3d_axis(axis, bounds, elev=24, azim=-55)
    axis.set_title(
        "A  点级 SWC：连续采样点与显式父子边\n"
        f"真实根部子树前 4 层：{len(source_node_ids)} 个不同节点 / "
        f"{focus_parent_edge_count} 条父子边",
        loc="left",
        fontsize=11,
    )

    axis = figure.add_subplot(2, 2, 2, projection="3d")
    for branch in focus_branches:
        path = paths[branch.branch_id]
        _add_segments_3d(
            axis,
            _path_segments((path,)),
            color=branch_colors[branch.branch_id],
            linewidth=2.1,
            alpha=0.96,
        )
    critical_ids = tuple(
        dict.fromkeys(
            [focus_branches[0].upstream_node_id]
            + [branch.downstream_node_id for branch in focus_branches]
        )
    )
    critical_points = np.asarray(
        [analysis_swc.points_um[index_by_id[node_id]] for node_id in critical_ids]
    )
    axis.scatter(
        critical_points[:, 0], critical_points[:, 1], critical_points[:, 2],
        s=30, color="#F4C95D", edgecolor="#6D4C41", linewidth=0.55,
        depthshade=False,
    )
    _style_3d_axis(axis, bounds, elev=24, azim=-55)
    axis.set_title(
        "B  以关键节点为边界压缩连续点串\n"
        f"相同父子边被整理为 {len(focus_branches)} 条真实分支；原节点序列仍保留",
        loc="left",
        fontsize=11,
    )

    axis = figure.add_subplot(2, 2, 3, projection="3d")
    _add_segments_3d(axis, raw_segments, color="#00889A", linewidth=1.45, alpha=0.91)
    arrow_count = _add_branch_direction_arrows(axis, ordered_paths)
    role_styles = {
        "inferred_inlet": ("#2CA25F", "结构根（不是已确认入口）", "o"),
        "divergence_junction": ("#F4C95D", "分叉节点", "D"),
        "inferred_outlet": ("#D95F59", "结构叶（不是已确认出口）", "s"),
    }
    for role, (color, _label, marker) in role_styles.items():
        ids = [node_id for node_id in critical_ids if node_roles.get(node_id) == role]
        if not ids:
            continue
        points = np.asarray([analysis_swc.points_um[index_by_id[node_id]] for node_id in ids])
        axis.scatter(
            points[:, 0], points[:, 1], points[:, 2],
            s=48, color=color, marker=marker, edgecolor="#37474F",
            linewidth=0.55, depthshade=False,
        )
    _style_3d_axis(axis, bounds, elev=24, azim=-55)
    axis.legend(
        handles=[
            Line2D([0], [0], marker=marker, color="none", markerfacecolor=color,
                   markeredgecolor="#37474F", label=label, markersize=7)
            for color, label, marker in role_styles.values()
        ]
        + [Line2D([0], [0], color="#F28E2B", linewidth=2, label="SWC parent → current")],
        loc="upper right",
        fontsize=8,
    )
    axis.set_title(
        "C  结构方向与关键节点角色\n"
        f"显示 {arrow_count} 个分支方向箭头；箭头来自 parent_id，不是实测血流",
        loc="left",
        fontsize=11,
    )

    axis = figure.add_subplot(2, 2, 4)
    hierarchy = _hierarchy_graph(focus_branches)
    root_id = next(node for node, indegree in hierarchy.in_degree() if indegree == 0)
    relative_depth = nx.single_source_shortest_path_length(hierarchy, root_id)
    raw_y: dict[int, float] = {}
    next_leaf_y = 0

    def assign_tree_y(branch_id: int) -> float:
        nonlocal next_leaf_y
        daughters = list(hierarchy.successors(branch_id))
        if not daughters:
            raw_y[branch_id] = float(next_leaf_y)
            next_leaf_y += 1
        else:
            raw_y[branch_id] = float(np.mean([assign_tree_y(child) for child in daughters]))
        return raw_y[branch_id]

    assign_tree_y(root_id)
    maximum_leaf_y = max(float(next_leaf_y - 1), 1.0)
    positions = {
        branch_id: (
            0.10 + relative_depth[branch_id] * 0.285,
            0.92 - 0.84 * raw_y[branch_id] / maximum_leaf_y,
        )
        for branch_id in focus_ids
    }
    for parent_id, daughter_id in hierarchy.edges:
        axis.add_patch(
            FancyArrowPatch(
                positions[parent_id], positions[daughter_id],
                arrowstyle="-|>", mutation_scale=13, linewidth=1.15,
                color="#607D8B", connectionstyle="arc3,rad=0.0",
                shrinkA=17, shrinkB=17, zorder=2,
            )
        )
    for branch_id in focus_ids:
        branch = by_id[branch_id]
        x, y = positions[branch_id]
        axis.scatter(
            x, y, s=540, color=branch_colors[branch_id], alpha=0.92,
            edgecolor="#37474F", linewidth=0.75, zorder=3,
        )
        axis.text(
            x, y,
            f"B{branch_id}\nd{branch.depth} · S{branch.strahler_order} · H{branch.horsfield_order}",
            ha="center", va="center", fontsize=7.3, zorder=4,
        )
    axis.set_xlim(0.0, 1.05)
    axis.set_ylim(0.0, 1.0)
    axis.axis("off")
    axis.set_title(
        "D  同一真实子树的分支级有向层级图\n"
        "B=分支编号，d=全局深度，S=Strahler 阶次，H=Horsfield 阶次",
        loc="left",
        fontsize=11,
    )
    figure.text(
        0.5,
        0.006,
        f"完整 analysis_swc：{int(graph_summary['source_node_count']):,} 点级节点 / "
        f"{int(graph_summary['source_edge_count']):,} 父子边 → "
        f"{int(graph_summary['critical_node_count']):,} 关键节点 / "
        f"{int(graph_summary['branch_count']):,} 聚合分支；全部几何和编号来自实际导出数据",
        ha="center",
        fontsize=10.5,
    )
    return _save_figure(figure, output_path, dpi), focus_ids


def _branch_path_midpoint(path: np.ndarray) -> np.ndarray:
    """Return the half-arclength point of one real branch polyline."""

    points = np.asarray(path, dtype=float)
    if len(points) == 1:
        return points[0]
    segment_lengths = np.linalg.norm(np.diff(points, axis=0), axis=1)
    total_length = float(np.sum(segment_lengths))
    if total_length <= 1.0e-12:
        return np.mean(points, axis=0)
    cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
    target = total_length * 0.5
    segment_index = min(int(np.searchsorted(cumulative, target, side="right") - 1), len(points) - 2)
    segment_length = float(segment_lengths[segment_index])
    fraction = 0.0 if segment_length <= 1.0e-12 else (
        target - float(cumulative[segment_index])
    ) / segment_length
    return points[segment_index] + fraction * (points[segment_index + 1] - points[segment_index])


def _add_metric_colored_branches_3d(
    axis: plt.Axes,
    branch_ids: Iterable[int],
    paths: dict[int, np.ndarray],
    branches_by_id: dict[int, HierarchyBranch],
    *,
    metric_name: str,
    cmap_name: str,
    linewidth: float = 2.7,
    alpha: float = 0.96,
) -> tuple[object, Normalize]:
    """Draw real branch paths colored by one exported integer hierarchy metric."""

    ordered_ids = tuple(branch_ids)
    values = [int(getattr(branches_by_id[branch_id], metric_name)) for branch_id in ordered_ids]
    lower = min(values)
    upper = max(values)
    norm = Normalize(vmin=lower - 0.5, vmax=upper + 0.5)
    cmap = plt.get_cmap(cmap_name)
    for branch_id, value in zip(ordered_ids, values):
        _add_segments_3d(
            axis,
            _path_segments((paths[branch_id],)),
            color=cmap(norm(value)),
            linewidth=linewidth,
            alpha=alpha,
        )
    return cmap, norm


def _add_integer_metric_colorbar(
    figure: plt.Figure,
    axis: plt.Axes,
    cmap: object,
    norm: Normalize,
    *,
    label: str,
    ticks: Iterable[int],
) -> None:
    mappable = ScalarMappable(norm=norm, cmap=cmap)  # type: ignore[arg-type]
    mappable.set_array([])
    colorbar = figure.colorbar(mappable, ax=axis, shrink=0.72, pad=0.025)
    colorbar.set_label(label)
    colorbar.set_ticks(tuple(ticks))


def _label_branch_metric_3d(
    axis: plt.Axes,
    branch_ids: Iterable[int],
    paths: dict[int, np.ndarray],
    branches_by_id: dict[int, HierarchyBranch],
    *,
    metric_name: str,
    symbol: str,
) -> None:
    for branch_id in branch_ids:
        midpoint = _branch_path_midpoint(paths[branch_id])
        value = int(getattr(branches_by_id[branch_id], metric_name))
        axis.text(
            midpoint[0], midpoint[1], midpoint[2],
            f"B{branch_id}: {symbol}={value}",
            fontsize=8.2,
            ha="center",
            va="bottom",
            bbox={"boxstyle": "round,pad=0.18", "facecolor": "white", "alpha": 0.80, "edgecolor": "#90A4AE"},
        )


def _deepest_branch_path_ids(branches: list[HierarchyBranch]) -> tuple[int, ...]:
    graph = _hierarchy_graph(branches)
    by_id = {branch.branch_id: branch for branch in branches}
    target = max(branches, key=lambda branch: (branch.depth, branch.raw_length_um, -branch.branch_id))
    roots = [node for node, indegree in graph.in_degree() if indegree == 0]
    usable_roots = [root for root in roots if nx.has_path(graph, root, target.branch_id)]
    if not usable_roots:
        raise ValueError(f"No structural root reaches deepest branch {target.branch_id}")
    root = min(usable_roots)
    path_ids = tuple(nx.shortest_path(graph, root, target.branch_id))
    depths = [by_id[branch_id].depth for branch_id in path_ids]
    if any(child_depth != parent_depth + 1 for parent_depth, child_depth in zip(depths[:-1], depths[1:])):
        raise ValueError("Depth example path does not increase by one at every branch relation")
    return path_ids


def _maximum_horsfield_path_ids(branches: list[HierarchyBranch]) -> tuple[int, ...]:
    graph = _hierarchy_graph(branches)
    by_id = {branch.branch_id: branch for branch in branches}
    roots = [node for node, indegree in graph.in_degree() if indegree == 0]
    if not roots:
        raise ValueError("The exported branch hierarchy has no structural root")
    current = max(roots, key=lambda branch_id: (by_id[branch_id].horsfield_order, -branch_id))
    path_ids = [current]
    while True:
        daughters = list(graph.successors(current))
        if not daughters:
            break
        next_branch = max(
            daughters,
            key=lambda branch_id: (by_id[branch_id].horsfield_order, -branch_id),
        )
        if by_id[next_branch].horsfield_order != by_id[current].horsfield_order - 1:
            raise ValueError("Horsfield example path does not decrease by one toward the structural leaf")
        path_ids.append(next_branch)
        current = next_branch
    return tuple(path_ids)


def _select_strahler_example(
    branches: list[HierarchyBranch],
) -> tuple[int, tuple[int, ...], tuple[int, ...]]:
    graph = _hierarchy_graph(branches)
    by_id = {branch.branch_id: branch for branch in branches}
    candidates: list[HierarchyBranch] = []
    for branch in branches:
        daughter_ids = tuple(graph.successors(branch.branch_id))
        if len(daughter_ids) < 2:
            continue
        daughter_orders = [by_id[branch_id].strahler_order for branch_id in daughter_ids]
        maximum = max(daughter_orders)
        if branch.strahler_order == maximum + 1 and daughter_orders.count(maximum) >= 2:
            candidates.append(branch)
    if not candidates:
        raise ValueError("No real equal-order Strahler increment is present in the branch hierarchy")
    preferred = [branch for branch in candidates if branch.strahler_order == 3] or candidates
    root = min(
        preferred,
        key=lambda branch: (branch.downstream_branch_count, branch.depth, branch.branch_id),
    )
    subtree_ids = tuple(
        sorted(
            {root.branch_id, *nx.descendants(graph, root.branch_id)},
            key=lambda branch_id: (by_id[branch_id].depth, branch_id),
        )
    )
    daughter_ids = tuple(graph.successors(root.branch_id))
    return root.branch_id, subtree_ids, daughter_ids


def _select_horsfield_example(
    branches: list[HierarchyBranch],
) -> tuple[int, tuple[int, ...], tuple[int, ...]]:
    graph = _hierarchy_graph(branches)
    by_id = {branch.branch_id: branch for branch in branches}
    candidates: list[HierarchyBranch] = []
    for branch in branches:
        daughter_ids = tuple(graph.successors(branch.branch_id))
        if len(daughter_ids) < 2:
            continue
        daughter_orders = [by_id[branch_id].horsfield_order for branch_id in daughter_ids]
        if len(set(daughter_orders)) > 1 and branch.horsfield_order == max(daughter_orders) + 1:
            candidates.append(branch)
    if not candidates:
        raise ValueError("No real unequal-daughter Horsfield example is present in the branch hierarchy")
    root = min(
        candidates,
        key=lambda branch: (
            branch.downstream_branch_count,
            branch.horsfield_order,
            branch.depth,
            branch.branch_id,
        ),
    )
    subtree_ids = tuple(
        sorted(
            {root.branch_id, *nx.descendants(graph, root.branch_id)},
            key=lambda branch_id: (by_id[branch_id].depth, branch_id),
        )
    )
    daughter_ids = tuple(graph.successors(root.branch_id))
    return root.branch_id, subtree_ids, daughter_ids


def _figure_branch_depth_3d(
    analysis_swc: SWCData,
    branches: list[HierarchyBranch],
    output_path: Path,
    dpi: int,
) -> tuple[Path, dict[str, object]]:
    """Explain branch depth with a real deepest path and one parent-child calculation."""

    by_id = {branch.branch_id: branch for branch in branches}
    path_ids = _deepest_branch_path_ids(branches)
    path_branches = [by_id[branch_id] for branch_id in path_ids]
    paths = _branch_point_paths(path_branches, analysis_swc)
    ordered_paths = [paths[branch_id] for branch_id in path_ids]
    depths = [by_id[branch_id].depth for branch_id in path_ids]
    calculation_ids = path_ids[: min(3, len(path_ids))]

    figure = plt.figure(figsize=(17.5, 8.3), constrained_layout=True)
    figure.suptitle("分支深度：从结构根沿父子关系逐层向下计数", fontsize=18, fontweight="semibold")

    axis = figure.add_subplot(1, 2, 1, projection="3d")
    cmap, norm = _add_metric_colored_branches_3d(
        axis, path_ids, paths, by_id,
        metric_name="depth", cmap_name="viridis", linewidth=3.0,
    )
    _add_branch_direction_arrows(axis, ordered_paths[::4], color="#E76F51")
    root_point = ordered_paths[0][0]
    leaf_point = ordered_paths[-1][-1]
    axis.scatter(root_point[0], root_point[1], root_point[2], s=70, marker="s", color="#D1495B", depthshade=False)
    axis.scatter(leaf_point[0], leaf_point[1], leaf_point[2], s=70, marker="o", color="#2A9D8F", depthshade=False)
    axis.text(root_point[0], root_point[1], root_point[2], " 结构根：d=0", fontsize=9)
    axis.text(leaf_point[0], leaf_point[1], leaf_point[2], f" 最深分支：d={depths[-1]}", fontsize=9)
    _style_3d_axis(axis, _paths_bounds(ordered_paths), elev=23, azim=-55)
    _add_integer_metric_colorbar(
        figure, axis, cmap, norm,
        label="分支深度 d",
        ticks=range(depths[0], depths[-1] + 1, 2),
    )
    axis.set_title(
        "A  概念：从根到当前分支经过多少级连接\n"
        f"真实最深路径包含 {len(path_ids)} 条分支，d={depths[0]}～{depths[-1]}",
        loc="left", fontsize=11,
    )

    axis = figure.add_subplot(1, 2, 2, projection="3d")
    calc_paths = [paths[branch_id] for branch_id in calculation_ids]
    cmap, norm = _add_metric_colored_branches_3d(
        axis, calculation_ids, paths, by_id,
        metric_name="depth", cmap_name="viridis", linewidth=4.0,
    )
    _add_branch_direction_arrows(axis, calc_paths, color="#E76F51")
    _label_branch_metric_3d(
        axis, calculation_ids, paths, by_id,
        metric_name="depth", symbol="d",
    )
    _style_3d_axis(axis, _paths_bounds(calc_paths, padding_fraction=0.16), elev=24, azim=-52)
    _add_integer_metric_colorbar(
        figure, axis, cmap, norm,
        label="分支深度 d",
        ticks=[by_id[branch_id].depth for branch_id in calculation_ids],
    )
    equations = "；".join(
        f"d(B{child_id})=d(B{parent_id})+1={by_id[child_id].depth}"
        for parent_id, child_id in zip(calculation_ids[:-1], calculation_ids[1:])
    )
    axis.set_title(
        "B  真实计算例子：每经过一次父分支→子分支，深度加 1\n" + equations,
        loc="left", fontsize=11,
    )
    metadata = {
        "source_stage": "unchanged analysis_swc branch hierarchy",
        "coordinate_unit": "um",
        "deepest_path_branch_ids": list(path_ids),
        "deepest_path_depths": depths,
        "calculation_branch_ids": list(calculation_ids),
        "calculation_rule": "depth(child_branch) = depth(parent_branch) + 1",
    }
    return _save_figure(figure, output_path, dpi), metadata


def _figure_strahler_order_3d(
    analysis_swc: SWCData,
    branches: list[HierarchyBranch],
    output_path: Path,
    dpi: int,
) -> tuple[Path, dict[str, object]]:
    """Explain Strahler recursion on one compact real three-order subtree."""

    graph = _hierarchy_graph(branches)
    by_id = {branch.branch_id: branch for branch in branches}
    root_id, subtree_ids, daughter_ids = _select_strahler_example(branches)
    subtree_branches = [by_id[branch_id] for branch_id in subtree_ids]
    paths = _branch_point_paths(subtree_branches, analysis_swc)
    subtree_paths = [paths[branch_id] for branch_id in subtree_ids]
    calculation_ids = (root_id, *daughter_ids)
    child_orders = [by_id[branch_id].strahler_order for branch_id in daughter_ids]
    maximum_order = max(by_id[branch_id].strahler_order for branch_id in subtree_ids)

    figure = plt.figure(figsize=(17.5, 8.3), constrained_layout=True)
    figure.suptitle("Strahler 阶次：相同阶次的下游结构同时出现时，上游阶次才增加", fontsize=18, fontweight="semibold")

    axis = figure.add_subplot(1, 2, 1, projection="3d")
    cmap, norm = _add_metric_colored_branches_3d(
        axis, subtree_ids, paths, by_id,
        metric_name="strahler_order", cmap_name="cividis", linewidth=3.2,
    )
    leaf_ids = [branch_id for branch_id in subtree_ids if graph.out_degree(branch_id) == 0]
    leaf_points = np.asarray([paths[branch_id][-1] for branch_id in leaf_ids])
    axis.scatter(
        leaf_points[:, 0], leaf_points[:, 1], leaf_points[:, 2],
        s=48, marker="o", color="#D95F59", edgecolor="#6D2825", linewidth=0.5,
        depthshade=False, label="结构末端：S=1",
    )
    _label_branch_metric_3d(
        axis, subtree_ids, paths, by_id,
        metric_name="strahler_order", symbol="S",
    )
    _style_3d_axis(axis, _paths_bounds(subtree_paths, padding_fraction=0.12), elev=25, azim=-58)
    _add_integer_metric_colorbar(
        figure, axis, cmap, norm,
        label="Strahler 阶次 S",
        ticks=range(1, maximum_order + 1),
    )
    axis.legend(loc="upper right", fontsize=8)
    axis.set_title(
        "A  概念：从结构末端 S=1 开始，向上游逐级递推\n"
        f"真实子树 B{root_id} 含 {len(subtree_ids)} 条分支，完整呈现 S=1～{maximum_order}",
        loc="left", fontsize=11,
    )

    axis = figure.add_subplot(1, 2, 2, projection="3d")
    calc_paths = [paths[branch_id] for branch_id in calculation_ids]
    cmap, norm = _add_metric_colored_branches_3d(
        axis, calculation_ids, paths, by_id,
        metric_name="strahler_order", cmap_name="cividis", linewidth=4.2,
    )
    _add_branch_direction_arrows(axis, [path[::-1] for path in calc_paths], color="#E76F51")
    _label_branch_metric_3d(
        axis, calculation_ids, paths, by_id,
        metric_name="strahler_order", symbol="S",
    )
    junction = paths[root_id][-1]
    axis.scatter(junction[0], junction[1], junction[2], s=72, marker="D", color="#F4C95D", depthshade=False)
    _style_3d_axis(axis, _paths_bounds(calc_paths, padding_fraction=0.18), elev=25, azim=-58)
    _add_integer_metric_colorbar(
        figure, axis, cmap, norm,
        label="Strahler 阶次 S",
        ticks=range(min(child_orders), by_id[root_id].strahler_order + 1),
    )
    daughter_description = "、".join(
        f"B{branch_id}（S={by_id[branch_id].strahler_order}）" for branch_id in daughter_ids
    )
    axis.set_title(
        "B  真实计算例子：两个最高下游阶次相同，上游阶次加 1\n"
        f"{daughter_description} → B{root_id}：S={max(child_orders)}+1={by_id[root_id].strahler_order}",
        loc="left", fontsize=11,
    )
    metadata = {
        "source_stage": "unchanged analysis_swc branch hierarchy",
        "coordinate_unit": "um",
        "example_root_branch_id": root_id,
        "example_subtree_branch_ids": list(subtree_ids),
        "example_daughter_branch_ids": list(daughter_ids),
        "daughter_strahler_orders": child_orders,
        "root_strahler_order": by_id[root_id].strahler_order,
        "calculation_rule": "increment only when at least two maximum daughter orders are equal",
    }
    return _save_figure(figure, output_path, dpi), metadata


def _figure_horsfield_order_3d(
    analysis_swc: SWCData,
    branches: list[HierarchyBranch],
    output_path: Path,
    dpi: int,
) -> tuple[Path, dict[str, object]]:
    """Explain Horsfield recursion with a real longest path and unequal daughters."""

    by_id = {branch.branch_id: branch for branch in branches}
    longest_path_ids = _maximum_horsfield_path_ids(branches)
    example_root_id, example_subtree_ids, daughter_ids = _select_horsfield_example(branches)
    needed_ids = tuple(dict.fromkeys((*longest_path_ids, *example_subtree_ids)))
    paths = _branch_point_paths([by_id[branch_id] for branch_id in needed_ids], analysis_swc)
    longest_paths = [paths[branch_id] for branch_id in longest_path_ids]
    horsfield_values = [by_id[branch_id].horsfield_order for branch_id in longest_path_ids]
    calculation_ids = (example_root_id, *daughter_ids)
    child_orders = [by_id[branch_id].horsfield_order for branch_id in daughter_ids]

    figure = plt.figure(figsize=(17.5, 8.3), constrained_layout=True)
    figure.suptitle("Horsfield 阶次：由最长的下游层级路径决定", fontsize=18, fontweight="semibold")

    axis = figure.add_subplot(1, 2, 1, projection="3d")
    cmap, norm = _add_metric_colored_branches_3d(
        axis, longest_path_ids, paths, by_id,
        metric_name="horsfield_order", cmap_name="plasma", linewidth=3.0,
    )
    _add_branch_direction_arrows(axis, [path[::-1] for path in longest_paths[::4]], color="#2A9D8F")
    root_point = longest_paths[0][0]
    leaf_point = longest_paths[-1][-1]
    axis.scatter(root_point[0], root_point[1], root_point[2], s=70, marker="s", color="#D1495B", depthshade=False)
    axis.scatter(leaf_point[0], leaf_point[1], leaf_point[2], s=70, marker="o", color="#2A9D8F", depthshade=False)
    axis.text(root_point[0], root_point[1], root_point[2], f" H={horsfield_values[0]}", fontsize=9)
    axis.text(leaf_point[0], leaf_point[1], leaf_point[2], " H=1", fontsize=9)
    _style_3d_axis(axis, _paths_bounds(longest_paths), elev=23, azim=-55)
    _add_integer_metric_colorbar(
        figure, axis, cmap, norm,
        label="Horsfield 阶次 H",
        ticks=range(1, max(horsfield_values) + 1, 2),
    )
    axis.set_title(
        "A  概念：末端 H=1，沿最长路径向上游每级加 1\n"
        f"真实最长路径含 {len(longest_path_ids)} 条分支，H={horsfield_values[0]}～1",
        loc="left", fontsize=11,
    )

    axis = figure.add_subplot(1, 2, 2, projection="3d")
    calc_paths = [paths[branch_id] for branch_id in calculation_ids]
    cmap, norm = _add_metric_colored_branches_3d(
        axis, calculation_ids, paths, by_id,
        metric_name="horsfield_order", cmap_name="plasma", linewidth=4.2,
    )
    _add_branch_direction_arrows(axis, [path[::-1] for path in calc_paths], color="#2A9D8F")
    _label_branch_metric_3d(
        axis, calculation_ids, paths, by_id,
        metric_name="horsfield_order", symbol="H",
    )
    junction = paths[example_root_id][-1]
    axis.scatter(junction[0], junction[1], junction[2], s=72, marker="D", color="#F4C95D", depthshade=False)
    _style_3d_axis(axis, _paths_bounds(calc_paths, padding_fraction=0.18), elev=25, azim=-58)
    _add_integer_metric_colorbar(
        figure, axis, cmap, norm,
        label="Horsfield 阶次 H",
        ticks=range(min(child_orders), by_id[example_root_id].horsfield_order + 1),
    )
    daughter_description = "、".join(
        f"B{branch_id}（H={by_id[branch_id].horsfield_order}）" for branch_id in daughter_ids
    )
    axis.set_title(
        "B  真实计算例子：只取下游最大阶次，再加 1\n"
        f"{daughter_description} → B{example_root_id}：H=max{tuple(child_orders)}+1={by_id[example_root_id].horsfield_order}",
        loc="left", fontsize=11,
    )
    metadata = {
        "source_stage": "unchanged analysis_swc branch hierarchy",
        "coordinate_unit": "um",
        "maximum_horsfield_path_branch_ids": list(longest_path_ids),
        "maximum_horsfield_path_orders": horsfield_values,
        "example_root_branch_id": example_root_id,
        "example_subtree_branch_ids": list(example_subtree_ids),
        "example_daughter_branch_ids": list(daughter_ids),
        "daughter_horsfield_orders": child_orders,
        "root_horsfield_order": by_id[example_root_id].horsfield_order,
        "calculation_rule": "H(parent_branch) = max(H(daughter_branches)) + 1",
    }
    return _save_figure(figure, output_path, dpi), metadata


def _figure_hierarchical_graph_audit(
    branches: list[HierarchyBranch],
    node_roles: dict[int, str],
    output_path: Path,
    dpi: int,
) -> Path:
    """Quantify the real hierarchy metrics and critical-node roles in four panels."""

    figure, axes = plt.subplots(2, 2, figsize=(14.5, 10.8), constrained_layout=True)
    figure.suptitle(
        "完整 analysis_swc 的层级指标与关键节点分布",
        fontsize=18,
        fontweight="semibold",
    )

    def count_values(values: Iterable[int]) -> tuple[np.ndarray, np.ndarray]:
        counts = Counter(values)
        keys = np.asarray(sorted(counts), dtype=int)
        return keys, np.asarray([counts[int(key)] for key in keys], dtype=int)

    axis = axes[0, 0]
    keys, values = count_values(branch.depth for branch in branches)
    axis.bar(keys, values, color="#2A9D8F", width=0.78)
    axis.set_xlabel("分支深度 d")
    axis.set_ylabel("分支数量")
    axis.set_xticks(keys)
    axis.grid(True, axis="y", color="#CFD8DC", linewidth=0.5)
    axis.set_title(f"A  深度分布\n从结构根向下共 {int(np.max(keys)) + 1} 层", loc="left")

    axis = axes[0, 1]
    keys, values = count_values(branch.strahler_order for branch in branches)
    bars = axis.bar(keys, values, color="#457B9D", width=0.72)
    axis.bar_label(bars, padding=2, fontsize=8)
    axis.set_xlabel("Strahler 阶次 S")
    axis.set_ylabel("分支数量")
    axis.set_xticks(keys)
    axis.grid(True, axis="y", color="#CFD8DC", linewidth=0.5)
    axis.set_title("B  Strahler 阶次分布\n同阶下游分支相遇时阶次增加", loc="left")

    axis = axes[1, 0]
    keys, values = count_values(branch.horsfield_order for branch in branches)
    bars = axis.bar(keys, values, color="#9B5DE5", width=0.76)
    axis.bar_label(bars, padding=2, fontsize=7.5)
    axis.set_xlabel("Horsfield 阶次 H")
    axis.set_ylabel("分支数量")
    axis.set_xticks(keys)
    axis.tick_params(axis="x", labelsize=8)
    axis.grid(True, axis="y", color="#CFD8DC", linewidth=0.5)
    axis.set_title("C  Horsfield 阶次分布\n上游阶次由最大下游阶次递推", loc="left")

    axis = axes[1, 1]
    role_order = ("inferred_inlet", "divergence_junction", "inferred_outlet")
    role_labels = ("结构根", "分叉节点", "结构叶")
    role_colors = ("#2CA25F", "#F4C95D", "#D95F59")
    role_counts = Counter(node_roles.values())
    values = [role_counts[role] for role in role_order]
    bars = axis.bar(role_labels, values, color=role_colors, width=0.62)
    axis.bar_label(bars, padding=3, fontsize=9)
    axis.set_ylabel("关键节点数量")
    axis.grid(True, axis="y", color="#CFD8DC", linewidth=0.5)
    axis.set_title(
        "D  关键节点角色\n角色来自图入度/出度，不等于生理入口或出口",
        loc="left",
    )
    return _save_figure(figure, output_path, dpi)


def _figure_anchor_sampling_final_results_3d(
    model: GlobalVascularModel,
    config: SamplingConfig,
    candidate_ids: tuple[int, ...],
    role_ids: dict[str, tuple[int, ...]],
    strategy_results: dict[str, AnchorStrategyResult],
    output_path: Path,
    dpi: int,
) -> Path:
    """Compare all implemented anchor modes on the real 3-D vascular model."""

    figure = plt.figure(figsize=(17.5, 13.0), constrained_layout=True)
    axes = [figure.add_subplot(2, 2, index + 1, projection="3d") for index in range(4)]
    figure.suptitle(
        "真实血管模型中的三维锚点候选与三种采样方式",
        fontsize=18,
        fontweight="semibold",
    )
    global_edges = _global_edge_array(model)
    global_bounds = model.model_bounds_xyz_um
    bounds_low = np.asarray(global_bounds[0::2], dtype=float)
    bounds_high = np.asarray(global_bounds[1::2], dtype=float)
    effective_roi_size = np.minimum(
        np.asarray(config.roi_size_um, dtype=float),
        bounds_high - bounds_low,
    )
    feasible_low = bounds_low + effective_roi_size * 0.5
    feasible_high = bounds_high - effective_roi_size * 0.5
    role_specs = {
        "root": ("结构根", "s", "#D1495B"),
        "terminal": ("源图末端", "o", "#2A9D8F"),
        "bifurcation": ("分叉节点", "D", "#9B5DE5"),
        "other_critical": ("其他关键节点", "P", "#F4A261"),
    }

    def points_for(node_ids: Iterable[int]) -> np.ndarray:
        ids = tuple(int(node_id) for node_id in node_ids)
        if not ids:
            return np.empty((0, 3), dtype=float)
        return np.asarray(
            [model.node_positions_um[model.node_index_by_id[node_id]] for node_id in ids],
            dtype=float,
        )

    axis = axes[0]
    _add_segments_3d(axis, global_edges, color="#90A4AE", linewidth=0.38, alpha=0.24)
    _add_box_3d(
        axis,
        feasible_low,
        feasible_high,
        color="#457B9D",
        linewidth=1.35,
        alpha=0.72,
    )
    role_handles: list[Line2D] = []
    for role, (label, marker, color) in role_specs.items():
        points = points_for(role_ids[role])
        if not len(points):
            continue
        axis.scatter(
            points[:, 0],
            points[:, 1],
            points[:, 2],
            marker=marker,
            s=36 if role == "bifurcation" else 48,
            color=color,
            edgecolor="white",
            linewidth=0.55,
            depthshade=False,
        )
        role_handles.append(
            Line2D(
                [],
                [],
                marker=marker,
                linestyle="none",
                color=color,
                markeredgecolor="white",
                label=f"{label}（{len(points)}）",
            )
        )
    role_handles.append(
        Line2D([], [], color="#457B9D", linewidth=1.5, label="允许放置 ROI 中心的内部范围")
    )
    axis.legend(handles=role_handles, loc="upper left", fontsize=8.0, frameon=True)
    axis.set_title(
        f"A  从真实结构关键节点建立候选池\n"
        f"{len(candidate_ids)} 个内部候选：{len(role_ids['root'])} 根、"
        f"{len(role_ids['terminal'])} 末端、{len(role_ids['bifurcation'])} 分叉",
        loc="left",
        pad=12,
    )
    _style_3d_axis(axis, global_bounds, elev=23, azim=-61)

    panels = (
        (
            axes[1],
            "random",
            "B  随机模式",
            "先按固定种子打乱候选顺序，再按最小间距逐个接纳",
        ),
        (
            axes[2],
            "farthest_point",
            "C  最远点模式（当前正式运行采用）",
            "从几何中心附近开始，每次选择离已选集合最远的候选",
        ),
        (
            axes[3],
            "poisson_disk",
            "D  Poisson-disk 模式",
            "当前代码与随机模式共用同一随机顺序和最小间距筛选",
        ),
    )
    candidate_points = points_for(candidate_ids)
    random_anchor_ids = strategy_results["random"].anchor_ids
    for panel_index, (axis, mode, heading, description) in enumerate(panels):
        result = strategy_results[mode]
        _add_segments_3d(axis, global_edges, color="#B0BEC5", linewidth=0.35, alpha=0.19)
        _add_box_3d(
            axis,
            feasible_low,
            feasible_high,
            color="#457B9D",
            linewidth=0.9,
            alpha=0.35,
        )
        axis.scatter(
            candidate_points[:, 0],
            candidate_points[:, 1],
            candidate_points[:, 2],
            marker="o",
            s=10,
            color="#90A4AE",
            alpha=0.56,
            depthshade=False,
        )
        selected_set = set(result.anchor_ids)
        for role, (_label, marker, color) in role_specs.items():
            selected_role_ids = [
                node_id for node_id in role_ids[role] if node_id in selected_set
            ]
            points = points_for(selected_role_ids)
            if not len(points):
                continue
            axis.scatter(
                points[:, 0],
                points[:, 1],
                points[:, 2],
                marker=marker,
                s=60,
                color=color,
                edgecolor="white",
                linewidth=0.65,
                depthshade=False,
            )
        closest_points = points_for(result.closest_pair_ids)
        if len(closest_points) == 2:
            _add_segments_3d(
                axis,
                np.asarray([[closest_points[0], closest_points[1]]]),
                color="#C62828",
                linewidth=2.0,
                alpha=0.95,
            )
            midpoint = np.mean(closest_points, axis=0)
            axis.text(
                midpoint[0],
                midpoint[1],
                midpoint[2],
                f"最近 {result.minimum_pairwise_distance_um:.1f} μm",
                color="#9B1C1C",
                fontsize=7.5,
            )
        for order, node_id in enumerate(result.anchor_ids[:3], start=1):
            point = model.node_positions_um[model.node_index_by_id[node_id]]
            axis.text(
                point[0],
                point[1],
                point[2] + 2.5,
                str(order),
                fontsize=7.2,
                color="#263238",
                fontweight="semibold",
                bbox={"boxstyle": "circle,pad=0.12", "facecolor": "white", "edgecolor": "#607D8B", "alpha": 0.88},
            )
        identical_note = ""
        if mode == "poisson_disk" and result.anchor_ids == random_anchor_ids:
            identical_note = "；本次与随机模式 13/13 完全相同"
        axis.set_title(
            f"{heading}\n{description}\n"
            f"选出 {len(result.anchor_ids)} 个；最近间距 "
            f"{result.minimum_pairwise_distance_um:.1f} μm，最近邻中位数 "
            f"{result.median_nearest_neighbor_distance_um:.1f} μm{identical_note}",
            loc="left",
            pad=12,
            fontsize=11.2,
        )
        if panel_index == 0:
            axis.legend(
                handles=[
                    Line2D([], [], marker="o", linestyle="none", color="#90A4AE", markersize=5, label="未选候选"),
                    Line2D([], [], marker="D", linestyle="none", color="#9B5DE5", markersize=7, label="已选关键节点（形状表示角色）"),
                    Line2D([], [], color="#C62828", linewidth=2.0, label="本组最近锚点对"),
                ],
                loc="upper left",
                fontsize=7.7,
                frameon=True,
            )
        _style_3d_axis(axis, global_bounds, elev=23, azim=-61)
    return _save_figure(figure, output_path, dpi)


def _figure_anchor_sampling_strategies_3d(
    model: GlobalVascularModel,
    config: SamplingConfig,
    candidate_ids: tuple[int, ...],
    role_ids: dict[str, tuple[int, ...]],
    strategy_results: dict[str, AnchorStrategyResult],
    process_trace: dict[str, object],
    output_path: Path,
    dpi: int,
) -> Path:
    """Explain each anchor algorithm through real intermediate 3-D decisions."""

    del role_ids, strategy_results  # Final counts are audited elsewhere; this figure shows process.
    figure = plt.figure(figsize=(19.5, 19.0), constrained_layout=False)
    figure.subplots_adjust(
        left=0.035,
        right=0.985,
        bottom=0.035,
        top=0.94,
        wspace=0.10,
        hspace=0.34,
    )
    axes = [figure.add_subplot(3, 3, index + 1, projection="3d") for index in range(9)]
    figure.suptitle(
        "真实 analysis_swc 上三种锚点方式的逐步三维采样原理",
        fontsize=18,
        fontweight="semibold",
    )
    global_edges = _global_edge_array(model)
    global_bounds = model.model_bounds_xyz_um
    bounds_low = np.asarray(global_bounds[0::2], dtype=float)
    bounds_high = np.asarray(global_bounds[1::2], dtype=float)
    effective_roi_size = np.minimum(
        np.asarray(config.roi_size_um, dtype=float),
        bounds_high - bounds_low,
    )
    feasible_low = bounds_low + effective_roi_size * 0.5
    feasible_high = bounds_high - effective_roi_size * 0.5
    candidate_ids = tuple(sorted(map(int, candidate_ids)))
    candidate_points = np.asarray(
        [model.node_positions_um[model.node_index_by_id[node_id]] for node_id in candidate_ids],
        dtype=float,
    )
    point_by_id = {
        node_id: model.node_positions_um[model.node_index_by_id[node_id]]
        for node_id in candidate_ids
    }
    random_events = list(process_trace["random_events"])  # type: ignore[arg-type]
    farthest_events = list(process_trace["farthest_events"])  # type: ignore[arg-type]
    centroid = np.asarray(process_trace["candidate_centroid_xyz_um"], dtype=float)

    def points_for(node_ids: Iterable[int]) -> np.ndarray:
        ids = tuple(int(node_id) for node_id in node_ids)
        if not ids:
            return np.empty((0, 3), dtype=float)
        return np.asarray([point_by_id[node_id] for node_id in ids], dtype=float)

    def add_context(axis: plt.Axes) -> None:
        _add_segments_3d(axis, global_edges, color="#B0BEC5", linewidth=0.28, alpha=0.13)
        _add_box_3d(
            axis,
            feasible_low,
            feasible_high,
            color="#457B9D",
            linewidth=0.8,
            alpha=0.28,
        )
        axis.scatter(
            candidate_points[:, 0],
            candidate_points[:, 1],
            candidate_points[:, 2],
            marker="o",
            s=10,
            color="#90A4AE",
            alpha=0.58,
            depthshade=False,
        )

    def scatter_ids(
        axis: plt.Axes,
        node_ids: Iterable[int],
        *,
        marker: str,
        color: str,
        size: float,
        edgecolor: str = "white",
    ) -> None:
        points = points_for(node_ids)
        if not len(points):
            return
        marker_edge_kwargs = (
            {} if marker == "x" else {"edgecolor": edgecolor, "linewidth": 0.65}
        )
        axis.scatter(
            points[:, 0],
            points[:, 1],
            points[:, 2],
            marker=marker,
            s=size,
            color=color,
            depthshade=False,
            **marker_edge_kwargs,
        )

    def connect_nodes(
        axis: plt.Axes,
        left_id: int,
        right_id: int,
        *,
        color: str,
        label: str,
    ) -> None:
        points = np.asarray([[point_by_id[left_id], point_by_id[right_id]]], dtype=float)
        _add_segments_3d(axis, points, color=color, linewidth=2.0, alpha=0.95)
        midpoint = np.mean(points[0], axis=0)
        axis.text(
            midpoint[0], midpoint[1], midpoint[2], label,
            color=color, fontsize=7.2,
            bbox={"boxstyle": "round,pad=0.15", "facecolor": "white", "edgecolor": color, "alpha": 0.84},
        )

    def add_exclusion_sphere(axis: plt.Axes, node_id: int) -> None:
        center = point_by_id[node_id]
        u = np.linspace(0.0, 2.0 * np.pi, 18)
        v = np.linspace(0.0, np.pi, 10)
        radius = float(config.min_anchor_distance_um)
        x = center[0] + radius * np.outer(np.cos(u), np.sin(v))
        y = center[1] + radius * np.outer(np.sin(u), np.sin(v))
        z = center[2] + radius * np.outer(np.ones_like(u), np.cos(v))
        axis.plot_wireframe(
            x, y, z,
            rstride=2,
            cstride=2,
            color="#E9C46A",
            linewidth=0.45,
            alpha=0.28,
        )

    def finish_axis(axis: plt.Axes, title: str) -> None:
        axis.set_title(title, loc="left", pad=9, fontsize=10.2)
        _style_3d_axis(axis, global_bounds, elev=23, azim=-61)

    # Row A: replay the shuffled scan order and the actual accept/reject test.
    first_six = random_events[:6]
    axis = axes[0]
    add_context(axis)
    ordered_ids = [int(event["candidate_id"]) for event in first_six]
    scatter_ids(axis, ordered_ids, marker="o", color="#F4A261", size=48)
    for order, node_id in enumerate(ordered_ids, start=1):
        point = point_by_id[node_id]
        axis.text(
            point[0], point[1], point[2] + 2.5, str(order),
            fontsize=7.2, color="#263238", fontweight="semibold",
            bbox={"boxstyle": "circle,pad=0.12", "facecolor": "white", "edgecolor": "#F4A261", "alpha": 0.9},
        )
    axis.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="none", color="#90A4AE", label="真实候选池"),
            Line2D([], [], marker="o", linestyle="none", color="#F4A261", label="随机扫描顺序 1—6"),
        ],
        loc="upper left", fontsize=7.6, frameon=True,
    )
    finish_axis(
        axis,
        "A1  固定种子打乱 42 个真实候选\n数字表示接下来逐个检查的顺序，并不表示最终都会入选",
    )

    reject_event = random_events[3]
    axis = axes[1]
    add_context(axis)
    scatter_ids(
        axis,
        reject_event["selected_ids_before"],  # type: ignore[arg-type]
        marker="*", color="#2A9D8F", size=72,
    )
    rejected_id = int(reject_event["candidate_id"])
    scatter_ids(axis, [rejected_id], marker="x", color="#C62828", size=76, edgecolor="#C62828")
    nearest_id = int(reject_event["nearest_selected_id"])
    reject_distance = float(reject_event["nearest_selected_distance_um"])
    connect_nodes(
        axis, rejected_id, nearest_id,
        color="#C62828", label=f"{reject_distance:.1f} μm < 45 μm",
    )
    axis.legend(
        handles=[
            Line2D([], [], marker="*", linestyle="none", color="#2A9D8F", markersize=10, label="已经接纳"),
            Line2D([], [], marker="x", linestyle="none", color="#C62828", markersize=8, label="当前候选被拒绝"),
        ],
        loc="upper left", fontsize=7.6, frameon=True,
    )
    finish_axis(
        axis,
        f"A2  扫描到第 4 个候选（节点 {rejected_id}）\n"
        f"它离已选节点 {nearest_id} 只有 {reject_distance:.1f} μm，因此跳过",
    )

    accept_event = random_events[5]
    axis = axes[2]
    add_context(axis)
    accepted_before = tuple(int(node_id) for node_id in accept_event["selected_ids_before"])
    rejected_before = tuple(
        int(event["candidate_id"])
        for event in first_six[:5]
        if not bool(event["accepted"])
    )
    scatter_ids(axis, accepted_before, marker="*", color="#2A9D8F", size=72)
    scatter_ids(axis, rejected_before, marker="x", color="#C62828", size=66, edgecolor="#C62828")
    accepted_id = int(accept_event["candidate_id"])
    nearest_id = int(accept_event["nearest_selected_id"])
    accept_distance = float(accept_event["nearest_selected_distance_um"])
    scatter_ids(axis, [accepted_id], marker="^", color="#F4A261", size=72)
    connect_nodes(
        axis, accepted_id, nearest_id,
        color="#2A9D8F", label=f"{accept_distance:.1f} μm ≥ 45 μm",
    )
    finish_axis(
        axis,
        f"A3  扫描到第 6 个候选（节点 {accepted_id}）\n"
        f"最近距离 {accept_distance:.1f} μm，满足阈值后接纳；随后继续扫描",
    )

    # Row B: nearest to centroid first, then repeatedly maximize distance to the set.
    first_farthest = farthest_events[0]
    axis = axes[3]
    add_context(axis)
    first_id = int(first_farthest["candidate_id"])
    axis.scatter(
        centroid[0], centroid[1], centroid[2],
        marker="P", s=74, color="#263238", edgecolor="white", linewidth=0.6,
        depthshade=False,
    )
    scatter_ids(axis, [first_id], marker="*", color="#2A9D8F", size=82)
    _add_segments_3d(
        axis,
        np.asarray([[centroid, point_by_id[first_id]]]),
        color="#455A64", linewidth=1.7, alpha=0.9,
    )
    axis.legend(
        handles=[
            Line2D([], [], marker="P", linestyle="none", color="#263238", label="42 个候选的几何中心"),
            Line2D([], [], marker="*", linestyle="none", color="#2A9D8F", markersize=10, label="第一个锚点"),
        ],
        loc="upper left", fontsize=7.6, frameon=True,
    )
    finish_axis(
        axis,
        f"B1  选择起点\n先计算候选池几何中心，再取离中心最近的真实节点 {first_id}",
    )

    for panel_offset, event in enumerate(farthest_events[1:3], start=1):
        axis = axes[3 + panel_offset]
        add_context(axis)
        selected_before = tuple(int(node_id) for node_id in event["selected_ids_before"])
        remaining_ids = tuple(node_id for node_id in candidate_ids if node_id not in selected_before)
        selected_points = points_for(selected_before)
        remaining_points = points_for(remaining_ids)
        distances = np.min(
            np.linalg.norm(
                remaining_points[:, np.newaxis, :] - selected_points[np.newaxis, :, :],
                axis=2,
            ),
            axis=1,
        )
        axis.scatter(
            remaining_points[:, 0], remaining_points[:, 1], remaining_points[:, 2],
            c=distances, cmap="viridis", s=24, alpha=0.84, depthshade=False,
        )
        scatter_ids(axis, selected_before, marker="*", color="#2A9D8F", size=76)
        next_id = int(event["candidate_id"])
        nearest_id = int(event["nearest_selected_id"])
        next_distance = float(event["criterion_distance_um"])
        scatter_ids(axis, [next_id], marker="*", color="#F4A261", size=105, edgecolor="#7F5539")
        connect_nodes(
            axis, next_id, nearest_id,
            color="#E76F51", label=f"d_min={next_distance:.1f} μm",
        )
        axis.text2D(
            0.02, 0.04, "暗色：靠近已选集合　亮色：远离已选集合",
            transform=axis.transAxes, fontsize=7.2, color="#455A64",
        )
        selection_index = int(event["selection_index"])
        finish_axis(
            axis,
            f"B{panel_offset + 1}  更新每个候选到已选集合的最近距离\n"
            f"取其中最大者作为第 {selection_index} 个锚点：节点 {next_id}，{next_distance:.1f} μm",
        )

    # Row C: the same current random-order implementation explained as exclusion balls.
    first_random = random_events[0]
    axis = axes[6]
    add_context(axis)
    seed_id = int(first_random["candidate_id"])
    scatter_ids(axis, [seed_id], marker="*", color="#2A9D8F", size=82)
    add_exclusion_sphere(axis, seed_id)
    axis.legend(
        handles=[
            Line2D([], [], marker="*", linestyle="none", color="#2A9D8F", markersize=10, label="已接纳锚点"),
            Line2D([], [], color="#E9C46A", linewidth=1.2, label="45 μm 排斥球"),
        ],
        loc="upper left", fontsize=7.6, frameon=True,
    )
    finish_axis(
        axis,
        f"C1  接纳随机顺序中的第一个节点 {seed_id}\n"
        "以该点为中心建立 45 μm 排斥球，球内不再接纳新锚点",
    )

    second_random = random_events[1]
    axis = axes[7]
    add_context(axis)
    second_id = int(second_random["candidate_id"])
    second_distance = float(second_random["nearest_selected_distance_um"])
    for node_id in (seed_id, second_id):
        add_exclusion_sphere(axis, node_id)
    scatter_ids(axis, [seed_id], marker="*", color="#2A9D8F", size=76)
    scatter_ids(axis, [second_id], marker="^", color="#F4A261", size=78)
    connect_nodes(
        axis, seed_id, second_id,
        color="#2A9D8F", label=f"{second_distance:.1f} μm > 45 μm",
    )
    finish_axis(
        axis,
        f"C2  检查下一节点 {second_id}\n"
        f"它位于已有排斥球之外（{second_distance:.1f} μm），因此接纳并增加新排斥球",
    )

    poisson_reject = random_events[3]
    axis = axes[8]
    add_context(axis)
    poisson_selected = tuple(int(node_id) for node_id in poisson_reject["selected_ids_before"])
    for node_id in poisson_selected:
        add_exclusion_sphere(axis, node_id)
    scatter_ids(axis, poisson_selected, marker="*", color="#2A9D8F", size=72)
    rejected_id = int(poisson_reject["candidate_id"])
    nearest_id = int(poisson_reject["nearest_selected_id"])
    reject_distance = float(poisson_reject["nearest_selected_distance_um"])
    scatter_ids(axis, [rejected_id], marker="x", color="#C62828", size=82, edgecolor="#C62828")
    connect_nodes(
        axis, rejected_id, nearest_id,
        color="#C62828", label=f"{reject_distance:.1f} μm < 45 μm",
    )
    finish_axis(
        axis,
        f"C3  后续候选 {rejected_id} 落入已有排斥球\n"
        f"距锚点 {nearest_id} 仅 {reject_distance:.1f} μm，因此拒绝并继续扫描",
    )
    return _save_figure(figure, output_path, dpi)


def _figure_process_overview(
    summary: dict[str, object],
    roi: ROIRecord,
    raw_component_count: int,
    output_path: Path,
    dpi: int,
) -> Path:
    figure, axis = plt.subplots(figsize=(16, 9))
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1)
    axis.axis("off")
    figure.suptitle(
        "4.3.5 真实连通局部区域的代表性采样：从完整血管到真实代表 ROI",
        fontsize=20,
        fontweight="semibold",
        y=0.97,
    )
    steps = [
        ("1  analysis_swc", f"选定原始单分量\n{int(summary['global_node_count']):,} 节点 / {int(summary['global_edge_count']):,} 边"),
        ("2  选择空间锚点", f"最远点策略，间距 ≥ {float(summary['minimum_anchor_distance_um']):.0f} μm\n实际得到 {int(summary['candidate_anchor_count'])} 个锚点"),
        ("3  建立 ROI 盒", "以锚点为中心\n80 × 80 × 120 μm³"),
        ("4  精确裁切", f"线段与立方体求交\n示例盒内有 {raw_component_count} 个连通片段"),
        ("5  保留锚点分量", f"只保留含 anchor {roi.anchor_id} 的网络\n长度 {roi.retained_component_length_um:.1f} μm"),
        ("6  保存边界与映射", f"local → global 可追溯\nCUT_PORT {roi.cut_port_count} / 源图终端 {roi.true_terminal_count}"),
        ("7  计算 9 维特征", "r10–r90 + 分支数 + 分叉数\n+ 总长度 + 环路秩"),
        ("8  特征分组", f"稳健尺度变换 + K-means\n{int(summary['valid_candidate_count'])} 个有效 ROI → {int(summary['n_clusters'])} 组"),
        ("9  选择真实代表", f"聚类中心只用于算距离\n最终选择 {int(summary['selected_count'])} 个真实 ROI"),
    ]
    centers = [
        (0.16, 0.75), (0.50, 0.75), (0.84, 0.75),
        (0.84, 0.46), (0.50, 0.46), (0.16, 0.46),
        (0.16, 0.17), (0.50, 0.17), (0.84, 0.17),
    ]
    fills = ["#E8F1F5"] * 6 + ["#FFF3D6", "#F2EAFE", "#E4F4EA"]
    width, height = 0.27, 0.19
    for index, ((heading, detail), (cx, cy), fill) in enumerate(
        zip(steps, centers, fills)
    ):
        axis.add_patch(
            FancyBboxPatch(
                (cx - width / 2, cy - height / 2),
                width,
                height,
                boxstyle="round,pad=0.012,rounding_size=0.018",
                facecolor=fill,
                edgecolor="#607D8B",
                linewidth=1.25,
            )
        )
        axis.text(cx, cy + 0.045, heading, ha="center", va="center", fontsize=13, fontweight="semibold")
        axis.text(cx, cy - 0.035, detail, ha="center", va="center", fontsize=10.5, linespacing=1.45)
        if index < len(centers) - 1:
            nx_, ny_ = centers[index + 1]
            dx, dy = nx_ - cx, ny_ - cy
            if abs(dx) > abs(dy):
                start = (cx + np.sign(dx) * width / 2, cy)
                end = (nx_ - np.sign(dx) * width / 2, ny_)
            else:
                start = (cx, cy - height / 2)
                end = (nx_, ny_ + height / 2)
            axis.add_patch(
                FancyArrowPatch(
                    start,
                    end,
                    arrowstyle="-|>",
                    mutation_scale=15,
                    linewidth=1.6,
                    color="#455A64",
                )
            )
    axis.text(
        0.5,
        0.015,
        "关键边界：ROI 来自真实完整模型；当前不推断生理入口/出口，也不计算血流或微泡轨迹。",
        ha="center",
        va="bottom",
        fontsize=10.5,
        color="#5D4037",
    )
    return _save_figure(figure, output_path, dpi)


def _figure_extraction_trace(
    model: GlobalVascularModel,
    rois: list[ROIRecord],
    all_anchor_ids: tuple[int, ...],
    roi: ROIRecord,
    raw_segments: list[RawClipSegment],
    raw_component_count: int,
    output_path: Path,
    dpi: int,
) -> Path:
    figure, axes = plt.subplots(2, 2, figsize=(15.5, 12), constrained_layout=True)
    figure.suptitle(
        f"真实 ROI 裁切与可追溯性示例：anchor {roi.anchor_id}（代表序号 {roi.selection_rank}）",
        fontsize=18,
        fontweight="semibold",
    )
    global_edges = _global_edge_array(model)
    valid_anchor_ids = {candidate.anchor_id for candidate in rois}
    selected_anchor_ids = {candidate.anchor_id for candidate in rois if candidate.is_representative}
    anchor_points = np.asarray(
        [model.node_positions_um[model.node_index_by_id[node_id]] for node_id in all_anchor_ids],
        dtype=float,
    )

    axis = axes[0, 0]
    _add_line_collection(axis, global_edges, color="#90A4AE", linewidth=0.45, alpha=0.34)
    for index, anchor_id in enumerate(all_anchor_ids):
        point = anchor_points[index]
        if anchor_id not in valid_anchor_ids:
            axis.scatter(point[0], point[2], marker="x", s=45, color="#C62828", zorder=9)
        elif anchor_id in selected_anchor_ids:
            axis.scatter(point[0], point[2], marker="*", s=70, color="#2A9D8F", edgecolor="white", linewidth=0.5, zorder=9)
        else:
            axis.scatter(point[0], point[2], marker="o", s=24, facecolor="white", edgecolor="#457B9D", linewidth=1.0, zorder=8)
    chosen = np.asarray(roi.anchor_position_um)
    axis.scatter(chosen[0], chosen[2], marker="*", s=190, color="#F4A261", edgecolor="#7F5539", linewidth=1.2, zorder=10)
    _draw_projected_box(axis, roi, color="#D1495B", linewidth=2.2)
    bounds = model.model_bounds_xyz_um
    axis.set_xlim(bounds[0], bounds[1])
    axis.set_ylim(bounds[4], bounds[5])
    _style_spatial_axis(
        axis,
        f"A  完整模型中选择锚点\n{len(all_anchor_ids)} 个锚点：{len(valid_anchor_ids)} 个有效，红色 × 为被拒绝锚点",
    )

    axis = axes[0, 1]
    other_components = sorted(
        {segment.component_id for segment in raw_segments if not segment.is_anchor_component}
    )
    pale_colors = ["#B0BEC5", "#CFD8DC", "#BCAAA4", "#C5CAE9", "#D7CCC8"]
    for segment in raw_segments:
        color = (
            "#00A6A6"
            if segment.is_anchor_component
            else pale_colors[other_components.index(segment.component_id) % len(pale_colors)]
        )
        _add_line_collection(axis, [segment.points_um], color=color, linewidth=1.35, alpha=0.95)
    axis.scatter(chosen[0], chosen[2], marker="*", s=150, color="#F4A261", edgecolor="#7F5539", linewidth=1.0, zorder=10, label="锚点")
    _draw_projected_box(axis, roi, color="#D1495B")
    axis.set_xlim(roi.bbox_min_um[0] - 3, roi.bbox_max_um[0] + 3)
    axis.set_ylim(roi.bbox_min_um[2] - 3, roi.bbox_max_um[2] + 3)
    _style_spatial_axis(
        axis,
        f"B  精确裁切后的原始片段\n同一立方体内出现 {raw_component_count} 个互不连通的分量；青色分量包含锚点",
    )
    axis.legend(loc="lower right", fontsize=9, frameon=True)

    axis = axes[1, 0]
    local_widths = np.clip(np.mean(roi.local_edge_radius_um, axis=1) * 0.45, 0.65, 3.2)
    _add_line_collection(
        axis,
        roi.local_edge_points_um,
        color="#007C91",
        linewidth=local_widths,
        alpha=0.92,
    )
    axis.scatter(chosen[0], chosen[2], marker="*", s=160, color="#F4A261", edgecolor="#7F5539", linewidth=1.0, zorder=10, label="anchor")
    if roi.cut_ports:
        cut_points = np.asarray([port.intersection_position_um for port in roi.cut_ports])
        axis.scatter(cut_points[:, 0], cut_points[:, 2], marker="s", s=55, color="#D1495B", edgecolor="white", linewidth=0.7, zorder=10, label="CUT_PORT")
    if roi.true_terminal_local_ids:
        terminals = roi.local_node_positions_um[list(roi.true_terminal_local_ids)]
        axis.scatter(terminals[:, 0], terminals[:, 2], marker="o", s=55, color="#2A9D8F", edgecolor="white", linewidth=0.7, zorder=10, label="源图终端")
    _draw_projected_box(axis, roi, color="#D1495B")
    axis.set_xlim(roi.bbox_min_um[0] - 3, roi.bbox_max_um[0] + 3)
    axis.set_ylim(roi.bbox_min_um[2] - 3, roi.bbox_max_um[2] + 3)
    kept_fraction = roi.retained_component_length_um / roi.raw_total_vessel_length_um
    _style_spatial_axis(
        axis,
        f"C  只保留锚点所在连通分量\n{roi.branch_count} 分支 / {roi.bifurcation_count} 分叉；保留原始盒内长度的 {kept_fraction:.0%}",
    )
    axis.legend(loc="lower right", fontsize=8.5, frameon=True)

    axis = axes[1, 1]
    margin = np.asarray((12.0, 12.0, 18.0))
    query_min = np.maximum(
        np.asarray(roi.bbox_min_um) - margin,
        np.asarray(model.model_bounds_xyz_um[0::2]),
    )
    query_max = np.minimum(
        np.asarray(roi.bbox_max_um) + margin,
        np.asarray(model.model_bounds_xyz_um[1::2]),
    )
    nearby_ids = EdgeSpatialIndex(model).query(query_min, query_max)
    nearby = np.asarray(
        [
            (
                model.edges[int(edge_id)].upstream_position_um,
                model.edges[int(edge_id)].downstream_position_um,
            )
            for edge_id in nearby_ids
        ]
    )
    _add_line_collection(axis, nearby, color="#B0BEC5", linewidth=0.7, alpha=0.45)
    _add_line_collection(axis, roi.local_edge_points_um, color="#007C91", linewidth=1.9, alpha=0.9, zorder=4)
    lengths = np.linalg.norm(
        roi.local_edge_points_um[:, 1] - roi.local_edge_points_um[:, 0], axis=1
    )
    highlight_indices = np.argsort(lengths)[-min(4, len(lengths)) :]
    highlight_colors = ("#E76F51", "#9B5DE5", "#F4A261", "#2A9D8F")
    for color, index in zip(highlight_colors, highlight_indices):
        points = roi.local_edge_points_um[index]
        _add_line_collection(axis, [points], color=color, linewidth=3.0, alpha=1.0, zorder=7)
        midpoint = np.mean(points, axis=0)
        axis.annotate(
            f"local e{int(roi.local_edge_ids[index])}\n→ global e{int(roi.local_edge_global_ids[index])}",
            xy=(midpoint[0], midpoint[2]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=7.5,
            color=color,
            bbox={"boxstyle": "round,pad=0.18", "fc": "white", "ec": color, "alpha": 0.85},
        )
    _draw_projected_box(axis, roi, color="#D1495B")
    axis.set_xlim(query_min[0], query_max[0])
    axis.set_ylim(query_min[2], query_max[2])
    _style_spatial_axis(
        axis,
        "D  局部—整体映射\n每条局部边保存源全局边编号；边界外灰线表示完整网络仍在继续",
    )
    return _save_figure(figure, output_path, dpi)


def _figure_representative_roi_sampling_3d(
    model: GlobalVascularModel,
    rois: list[ROIRecord],
    all_anchor_ids: tuple[int, ...],
    roi: ROIRecord,
    raw_segments: list[RawClipSegment],
    raw_component_count: int,
    summary: dict[str, object],
    output_path: Path,
    dpi: int,
) -> Path:
    """Show the saved representative-ROI workflow in its real 3-D coordinates."""

    figure = plt.figure(figsize=(19.5, 13.2), constrained_layout=True)
    axes = [figure.add_subplot(2, 3, index + 1, projection="3d") for index in range(6)]
    figure.suptitle(
        "真实连通局部区域的代表性采样：从完整网络到可追溯的三维代表 ROI",
        fontsize=18,
        fontweight="semibold",
    )

    global_edges = _global_edge_array(model)
    global_bounds = model.model_bounds_xyz_um
    valid_anchor_ids = {candidate.anchor_id for candidate in rois}
    representative_rois = sorted(
        (candidate for candidate in rois if candidate.is_representative),
        key=lambda candidate: candidate.selection_rank,
    )
    representative_anchor_ids = {candidate.anchor_id for candidate in representative_rois}
    chosen = np.asarray(roi.anchor_position_um, dtype=float)
    local_low = np.asarray(roi.bbox_min_um, dtype=float)
    local_high = np.asarray(roi.bbox_max_um, dtype=float)
    local_padding = np.asarray((3.0, 3.0, 3.0), dtype=float)
    local_bounds = tuple(
        value
        for low, high in zip(local_low - local_padding, local_high + local_padding)
        for value in (float(low), float(high))
    )

    # A: anchors are real global SWC nodes; rejected anchors never become ROI records.
    axis = axes[0]
    _add_segments_3d(axis, global_edges, color="#90A4AE", linewidth=0.38, alpha=0.24)
    for anchor_id in all_anchor_ids:
        point = model.node_positions_um[model.node_index_by_id[int(anchor_id)]]
        if anchor_id not in valid_anchor_ids:
            axis.scatter(
                point[0], point[1], point[2], marker="x", s=38,
                color="#C62828", linewidth=1.3, depthshade=False,
            )
        elif anchor_id in representative_anchor_ids:
            candidate = next(item for item in representative_rois if item.anchor_id == anchor_id)
            axis.scatter(
                point[0], point[1], point[2], marker="*", s=62,
                color=CLUSTER_COLORS[candidate.cluster_id % len(CLUSTER_COLORS)],
                edgecolor="white", linewidth=0.45, depthshade=False,
            )
        else:
            axis.scatter(
                point[0], point[1], point[2], marker="o", s=22,
                facecolor="white", edgecolor="#457B9D", linewidth=0.9,
                depthshade=False,
            )
    axis.scatter(
        chosen[0], chosen[1], chosen[2], marker="*", s=180,
        color="#F4A261", edgecolor="#7F5539", linewidth=1.1,
        depthshade=False,
    )
    axis.set_title(
        f"A  在完整 analysis_swc 中布置锚点\n"
        f"{len(all_anchor_ids)} 个锚点：{len(rois)} 个有效，"
        f"{len(all_anchor_ids) - len(rois)} 个因条件不足被拒绝",
        loc="left",
        pad=12,
    )
    axis.legend(
        handles=[
            Line2D([], [], marker="o", linestyle="none", markerfacecolor="white", markeredgecolor="#457B9D", label="有效候选锚点"),
            Line2D([], [], marker="*", linestyle="none", color="#2A9D8F", markersize=10, label="已选代表锚点"),
            Line2D([], [], marker="x", linestyle="none", color="#C62828", label="被拒绝锚点"),
            Line2D([], [], marker="*", linestyle="none", color="#F4A261", markeredgecolor="#7F5539", markersize=12, label=f"当前示例：{roi.anchor_id}"),
        ],
        loc="upper left",
        fontsize=7.7,
        frameon=True,
    )
    _style_3d_axis(axis, global_bounds, elev=23, azim=-61)

    # B: every wireframe is an actual saved candidate ROI box.
    axis = axes[1]
    _add_segments_3d(axis, global_edges, color="#B0BEC5", linewidth=0.34, alpha=0.18)
    for candidate in rois:
        is_chosen = candidate.roi_id == roi.roi_id
        color = (
            "#D1495B"
            if is_chosen
            else CLUSTER_COLORS[candidate.cluster_id % len(CLUSTER_COLORS)]
        )
        _add_box_3d(
            axis,
            np.asarray(candidate.bbox_min_um),
            np.asarray(candidate.bbox_max_um),
            color=color,
            linewidth=2.2 if is_chosen else (1.15 if candidate.is_representative else 0.65),
            alpha=1.0 if is_chosen else (0.82 if candidate.is_representative else 0.34),
        )
        center = np.asarray(candidate.bbox_center_um)
        axis.scatter(
            center[0], center[1], center[2],
            marker="*" if candidate.is_representative else "o",
            s=42 if candidate.is_representative else 10,
            color=color,
            depthshade=False,
        )
    axis.set_title(
        f"B  建立真实三维候选区域\n"
        f"{len(rois)} 个有效 ROI；细框为候选，粗框/星号为已选代表，红框为当前示例",
        loc="left",
        pad=12,
    )
    _style_3d_axis(axis, global_bounds, elev=23, azim=-61)

    # C: exact line-box clipping may reveal several unrelated components in one cube.
    axis = axes[2]
    component_ids = sorted({segment.component_id for segment in raw_segments})
    anchor_component_id = next(
        segment.component_id for segment in raw_segments if segment.is_anchor_component
    )
    component_colors = (
        "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462",
        "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F",
    )
    for component_index, component_id in enumerate(component_ids):
        component_segments = np.asarray(
            [segment.points_um for segment in raw_segments if segment.component_id == component_id],
            dtype=float,
        )
        is_anchor_component = component_id == anchor_component_id
        _add_segments_3d(
            axis,
            component_segments,
            color="#007C91" if is_anchor_component else component_colors[component_index % len(component_colors)],
            linewidth=2.2 if is_anchor_component else 0.95,
            alpha=0.98 if is_anchor_component else 0.72,
        )
    _add_box_3d(axis, local_low, local_high, color="#D1495B", linewidth=1.5, alpha=0.8)
    axis.scatter(
        chosen[0], chosen[1], chosen[2], marker="*", s=145,
        color="#F4A261", edgecolor="#7F5539", linewidth=1.0,
        depthshade=False,
    )
    axis.set_title(
        f"C  按立方体边界精确裁切\n"
        f"{len(raw_segments)} 条盒内线段组成 {raw_component_count} 个互不连通的分量；深青色分量含锚点",
        loc="left",
        pad=12,
    )
    axis.legend(
        handles=[
            Line2D([], [], color="#007C91", linewidth=2.5, label="包含锚点的分量"),
            Line2D([], [], color="#BEBADA", linewidth=2.0, label="同盒内的其他分量"),
            Line2D([], [], color="#D1495B", linewidth=1.5, label="ROI 空间边界"),
        ],
        loc="upper left",
        fontsize=7.7,
        frameon=True,
    )
    _style_3d_axis(axis, local_bounds, elev=20, azim=-53)

    # D: the final candidate is exactly one connected, anchor-containing component.
    axis = axes[3]
    local_widths = np.clip(np.mean(roi.local_edge_radius_um, axis=1) * 0.45, 0.65, 3.2)
    _add_segments_3d(
        axis,
        roi.local_edge_points_um,
        color="#007C91",
        linewidth=local_widths,
        alpha=0.96,
    )
    local_degree = np.zeros(roi.node_count, dtype=int)
    np.add.at(local_degree, roi.local_edges[:, 0], 1)
    np.add.at(local_degree, roi.local_edges[:, 1], 1)
    bifurcation_points = roi.local_node_positions_um[local_degree >= 3]
    if len(bifurcation_points):
        axis.scatter(
            bifurcation_points[:, 0], bifurcation_points[:, 1], bifurcation_points[:, 2],
            marker="D", s=34, color="#9B5DE5", edgecolor="white", linewidth=0.45,
            depthshade=False,
        )
    if roi.true_terminal_local_ids:
        terminal_points = roi.local_node_positions_um[list(roi.true_terminal_local_ids)]
        axis.scatter(
            terminal_points[:, 0], terminal_points[:, 1], terminal_points[:, 2],
            marker="o", s=46, color="#2A9D8F", edgecolor="white", linewidth=0.6,
            depthshade=False,
        )
    if roi.cut_ports:
        cut_points = np.asarray([port.intersection_position_um for port in roi.cut_ports])
        axis.scatter(
            cut_points[:, 0], cut_points[:, 1], cut_points[:, 2],
            marker="s", s=48, color="#D1495B", edgecolor="white", linewidth=0.6,
            depthshade=False,
        )
    axis.scatter(
        chosen[0], chosen[1], chosen[2], marker="*", s=145,
        color="#F4A261", edgecolor="#7F5539", linewidth=1.0,
        depthshade=False,
    )
    _add_box_3d(axis, local_low, local_high, color="#D1495B", linewidth=1.2, alpha=0.58)
    kept_fraction = roi.retained_component_length_um / roi.raw_total_vessel_length_um
    axis.set_title(
        f"D  仅保留锚点所在连通分量\n"
        f"{roi.node_count} 节点 / {roi.edge_count} 边 / {roi.branch_count} 分支；"
        f"保留盒内总长度的 {kept_fraction:.0%}",
        loc="left",
        pad=12,
    )
    axis.legend(
        handles=[
            Line2D([], [], color="#007C91", linewidth=2.5, label="保留的真实血管"),
            Line2D([], [], marker="D", linestyle="none", color="#9B5DE5", label="分叉节点"),
            Line2D([], [], marker="o", linestyle="none", color="#2A9D8F", label="源图终端"),
            Line2D([], [], marker="s", linestyle="none", color="#D1495B", label="CUT_PORT"),
        ],
        loc="upper left",
        fontsize=7.5,
        frameon=True,
    )
    _style_3d_axis(axis, local_bounds, elev=20, azim=-53)

    # E: retain enough of the surrounding global graph to make provenance visible.
    axis = axes[4]
    trace_margin = np.asarray((12.0, 12.0, 18.0), dtype=float)
    global_low = np.asarray(global_bounds[0::2], dtype=float)
    global_high = np.asarray(global_bounds[1::2], dtype=float)
    trace_low = np.maximum(local_low - trace_margin, global_low)
    trace_high = np.minimum(local_high + trace_margin, global_high)
    nearby_ids = EdgeSpatialIndex(model).query(trace_low, trace_high)
    nearby_segments = np.asarray(
        [
            (
                model.edges[int(edge_id)].upstream_position_um,
                model.edges[int(edge_id)].downstream_position_um,
            )
            for edge_id in nearby_ids
        ],
        dtype=float,
    ).reshape((-1, 2, 3))
    _add_segments_3d(axis, nearby_segments, color="#B0BEC5", linewidth=0.55, alpha=0.35)
    _add_segments_3d(axis, roi.local_edge_points_um, color="#007C91", linewidth=1.8, alpha=0.94)
    edge_lengths = np.linalg.norm(
        roi.local_edge_points_um[:, 1] - roi.local_edge_points_um[:, 0], axis=1
    )
    highlight_indices = np.unique(
        np.linspace(0, len(edge_lengths) - 1, min(3, len(edge_lengths)), dtype=int)
    )
    highlight_colors = ("#E76F51", "#9B5DE5", "#F4A261")
    for color, edge_index in zip(highlight_colors, highlight_indices):
        points = roi.local_edge_points_um[int(edge_index)]
        _add_segments_3d(
            axis, np.asarray([points]), color=color, linewidth=3.0, alpha=1.0
        )
        midpoint = np.mean(points, axis=0)
        axis.text(
            midpoint[0], midpoint[1], midpoint[2],
            f"L{int(roi.local_edge_ids[edge_index])}→G{int(roi.local_edge_global_ids[edge_index])}",
            color=color,
            fontsize=7.2,
        )
    _add_box_3d(axis, local_low, local_high, color="#D1495B", linewidth=1.4, alpha=0.75)
    trace_bounds = tuple(
        value
        for low, high in zip(trace_low, trace_high)
        for value in (float(low), float(high))
    )
    axis.set_title(
        "E  保存局部—整体身份映射\n"
        "青色局部边仍指向源全局边；盒外灰线说明血管在完整网络中继续延伸",
        loc="left",
        pad=12,
    )
    _style_3d_axis(axis, trace_bounds, elev=20, azim=-53)

    # F: selected representatives remain real 3-D subgraphs at their original positions.
    axis = axes[5]
    _add_segments_3d(axis, global_edges, color="#B0BEC5", linewidth=0.34, alpha=0.16)
    for candidate in representative_rois:
        color = CLUSTER_COLORS[candidate.cluster_id % len(CLUSTER_COLORS)]
        _add_segments_3d(
            axis,
            candidate.local_edge_points_um,
            color=color,
            linewidth=1.35,
            alpha=0.96,
        )
        _add_box_3d(
            axis,
            np.asarray(candidate.bbox_min_um),
            np.asarray(candidate.bbox_max_um),
            color=color,
            linewidth=1.15,
            alpha=0.76,
        )
        anchor = np.asarray(candidate.anchor_position_um)
        axis.text(
            anchor[0], anchor[1], anchor[2],
            f"R{candidate.selection_rank}",
            color=color,
            fontsize=7.4,
            fontweight="semibold",
        )
    represented_clusters = {candidate.cluster_id for candidate in representative_rois}
    all_clusters = {candidate.cluster_id for candidate in rois}
    missing_clusters = sorted(all_clusters - represented_clusters)
    axis.set_title(
        f"F  从特征分组中选出真实代表区域\n"
        f"{len(rois)} 个候选 → {int(summary['n_clusters'])} 个聚类 → "
        f"{len(representative_rois)} 个代表，实际覆盖 {len(represented_clusters)}/{len(all_clusters)} 个聚类",
        loc="left",
        pad=12,
    )
    if missing_clusters:
        axis.text2D(
            0.42,
            0.91,
            f"未覆盖聚类 {', '.join(map(str, missing_clusters))}：其候选未通过代表选择的空间去冗余约束",
            transform=axis.transAxes,
            fontsize=7.4,
            color="#7F1D1D",
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "edgecolor": "#D7CCC8", "alpha": 0.88},
        )
    axis.legend(
        handles=[
            Line2D(
                [],
                [],
                color=CLUSTER_COLORS[cluster_id % len(CLUSTER_COLORS)],
                linewidth=2.5,
                linestyle="--" if cluster_id in missing_clusters else "-",
                label=(
                    f"聚类 {cluster_id}（本次无代表）"
                    if cluster_id in missing_clusters
                    else f"聚类 {cluster_id}"
                ),
            )
            for cluster_id in sorted(all_clusters)
        ],
        loc="upper left",
        fontsize=7.3,
        ncol=2,
        frameon=True,
    )
    _style_3d_axis(axis, global_bounds, elev=23, azim=-61)

    return _save_figure(figure, output_path, dpi)


def _feature_matrix(
    rois: list[ROIRecord], feature_names: list[str]
) -> np.ndarray:
    rows: list[list[float]] = []
    for roi in rois:
        rows.append(
            [
                float(
                    roi.radius_features[name]
                    if name in roi.radius_features
                    else roi.structural_features[name]
                )
                for name in feature_names
            ]
        )
    return np.asarray(rows, dtype=float)


def _weighted_radius_ecdf(roi: ROIRecord) -> tuple[np.ndarray, np.ndarray]:
    radii = np.mean(roi.local_edge_radius_um, axis=1)
    lengths = np.linalg.norm(
        roi.local_edge_points_um[:, 1] - roi.local_edge_points_um[:, 0], axis=1
    )
    order = np.argsort(radii, kind="stable")
    ordered_radius = radii[order]
    cumulative = np.cumsum(lengths[order]) / np.sum(lengths)
    return ordered_radius, cumulative


def _point_radius_ecdf(roi: ROIRecord) -> tuple[np.ndarray, np.ndarray]:
    radii = np.sort(roi.local_node_radius_um)
    return radii, np.arange(1, len(radii) + 1, dtype=float) / len(radii)


def _read_cluster_centers(path: Path, feature_names: list[str]) -> np.ndarray:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    rows.sort(key=lambda row: int(row["cluster_id"]))
    return np.asarray(
        [[float(row[name]) for name in feature_names] for row in rows], dtype=float
    )


def _figure_features_and_selection(
    rois: list[ROIRecord],
    roi: ROIRecord,
    sampling_run: Path,
    summary: dict[str, object],
    output_path: Path,
    dpi: int,
) -> Path:
    scaler = _read_json(sampling_run / "features" / "scaler.json")
    feature_names = [str(name) for name in scaler["feature_order"]]  # type: ignore[index]
    median = np.asarray(scaler["median"], dtype=float)  # type: ignore[index]
    iqr = np.asarray(scaler["iqr"], dtype=float)  # type: ignore[index]
    raw_features = _feature_matrix(rois, feature_names)
    scaled = (raw_features - median) / np.where(np.abs(iqr) > 1.0e-12, iqr, 1.0)
    config_path = sampling_run / "config" / "sampling_config.yaml"
    radius_weight = _simple_yaml_number(config_path, "radius_feature_weight", 1.0)
    structure_weight = _simple_yaml_number(config_path, "structure_feature_weight", 1.0)
    weights = np.asarray(
        [radius_weight if name.startswith("r") else structure_weight for name in feature_names]
    )
    weighted = scaled * weights
    centers = _read_cluster_centers(
        sampling_run / "clustering" / "cluster_centers.csv", feature_names
    )
    centered = weighted - np.mean(weighted, axis=0)
    _u, singular, vt = np.linalg.svd(centered, full_matrices=False)
    basis = vt[:2].T
    scores = centered @ basis
    center_scores = (centers - np.mean(weighted, axis=0)) @ basis
    variance = singular**2
    variance_ratio = variance / np.sum(variance) if np.sum(variance) > 0 else variance

    figure, axes = plt.subplots(2, 2, figsize=(16, 11.5), constrained_layout=True)
    figure.suptitle(
        "从真实 ROI 特征到聚类与代表选择",
        fontsize=18,
        fontweight="semibold",
    )

    axis = axes[0, 0]
    weighted_r, weighted_f = _weighted_radius_ecdf(roi)
    point_r, point_f = _point_radius_ecdf(roi)
    axis.step(weighted_r, weighted_f, where="post", color="#007C91", linewidth=2.4, label="按血管弧长加权（实际采用）")
    axis.step(point_r, point_f, where="post", color="#90A4AE", linewidth=1.4, linestyle="--", label="按 SWC 点数统计（仅作对照）")
    quantile_colors = ("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51")
    for name, color in zip(("r10", "r25", "r50", "r75", "r90"), quantile_colors):
        value = roi.radius_features[name]
        axis.axvline(value, color=color, linewidth=1.0, alpha=0.75)
        axis.text(value, 0.03, f"{name}={value:.2f}", rotation=90, va="bottom", ha="right", fontsize=8, color=color)
    axis.set_title(
        f"A  弧长加权半径分布（anchor {roi.anchor_id}）\n长血管段贡献更大，避免中心线采样密度造成偏差",
        loc="left",
    )
    axis.set_xlabel("血管半径 r (μm)")
    axis.set_ylabel("累计血管长度比例 F(r)")
    axis.set_ylim(0, 1.03)
    axis.grid(True, color="#CFD8DC", linewidth=0.5)
    axis.legend(loc="lower right", fontsize=9)

    axis = axes[0, 1]
    ordered_indices = sorted(
        range(len(rois)),
        key=lambda index: (
            rois[index].cluster_id,
            0 if rois[index].is_representative else 1,
            rois[index].distance_to_cluster_center,
        ),
    )
    ordered = scaled[ordered_indices]
    vmax = max(1.0, float(np.nanpercentile(np.abs(ordered), 95)))
    image = axis.imshow(
        ordered,
        aspect="auto",
        cmap="coolwarm",
        vmin=-vmax,
        vmax=vmax,
        interpolation="nearest",
    )
    axis.set_xticks(range(len(feature_names)), [FEATURE_LABELS.get(name, name) for name in feature_names], rotation=35, ha="right")
    labels = [
        f"C{rois[index].cluster_id}  anchor {rois[index].anchor_id}"
        + (f"  ★R{rois[index].selection_rank}" if rois[index].is_representative else "")
        for index in ordered_indices
    ]
    axis.set_yticks(range(len(labels)), labels, fontsize=8.5)
    axis.set_title(
        "B  稳健尺度变换后的 9 维描述符\n每行是真实 ROI；★ 表示最终代表",
        loc="left",
    )
    colorbar = figure.colorbar(image, ax=axis, shrink=0.82, pad=0.02)
    colorbar.set_label("相对中位数的稳健尺度值")

    axis = axes[1, 0]
    for cluster_id in sorted({candidate.cluster_id for candidate in rois}):
        indices = [index for index, candidate in enumerate(rois) if candidate.cluster_id == cluster_id]
        axis.scatter(
            scores[indices, 0],
            scores[indices, 1],
            s=58,
            color=CLUSTER_COLORS[cluster_id % len(CLUSTER_COLORS)],
            edgecolor="white",
            linewidth=0.7,
            label=f"聚类 {cluster_id}",
            zorder=4,
        )
    selected_indices = [index for index, candidate in enumerate(rois) if candidate.is_representative]
    axis.scatter(
        scores[selected_indices, 0],
        scores[selected_indices, 1],
        marker="*",
        s=230,
        facecolor="none",
        edgecolor="#1F2933",
        linewidth=1.4,
        label="真实代表 ROI",
        zorder=7,
    )
    axis.scatter(
        center_scores[:, 0],
        center_scores[:, 1],
        marker="X",
        s=115,
        color="#1F2933",
        edgecolor="white",
        linewidth=0.7,
        label="聚类中心（仅特征坐标）",
        zorder=8,
    )
    for index in selected_indices:
        axis.annotate(
            f"R{rois[index].selection_rank}",
            (scores[index, 0], scores[index, 1]),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )
    axis.set_xlabel(f"主成分 1（{variance_ratio[0]:.1%}）")
    axis.set_ylabel(f"主成分 2（{variance_ratio[1]:.1%}）")
    axis.set_title(
        "C  多维特征的二维投影\n颜色表示相似组；中心 × 不会被当作人工血管输出",
        loc="left",
    )
    axis.grid(True, color="#CFD8DC", linewidth=0.5)
    axis.legend(loc="best", fontsize=8, ncol=2)

    axis = axes[1, 1]
    cluster_ids = sorted({candidate.cluster_id for candidate in rois})
    candidate_counts = [sum(candidate.cluster_id == cluster for candidate in rois) for cluster in cluster_ids]
    representative_counts = [
        sum(candidate.cluster_id == cluster and candidate.is_representative for candidate in rois)
        for cluster in cluster_ids
    ]
    x = np.arange(len(cluster_ids))
    width = 0.38
    axis.bar(x - width / 2, candidate_counts, width, color="#90A4AE", label="候选 ROI")
    axis.bar(x + width / 2, representative_counts, width, color="#2A9D8F", label="已选真实代表")
    for xi, value in zip(x - width / 2, candidate_counts):
        axis.text(xi, value + 0.08, str(value), ha="center", fontsize=9)
    for xi, value in zip(x + width / 2, representative_counts):
        axis.text(xi, value + 0.08, str(value), ha="center", fontsize=9)
    axis.set_xticks(x, [f"聚类 {cluster}" for cluster in cluster_ids])
    axis.set_ylabel("ROI 数量")
    axis.set_ylim(0, max(candidate_counts) + 1.2)
    represented_cluster_count = len(
        {candidate.cluster_id for candidate in rois if candidate.is_representative}
    )
    axis.set_title(
        f"D  覆盖均衡的代表选择\n{int(summary['valid_candidate_count'])} 个有效候选 → "
        f"{int(summary['selected_count'])} 个真实代表，实际覆盖 "
        f"{represented_cluster_count}/{int(summary['n_clusters'])} 个聚类",
        loc="left",
    )
    axis.grid(True, axis="y", color="#CFD8DC", linewidth=0.5)
    axis.legend(loc="upper right", fontsize=9)
    return _save_figure(figure, output_path, dpi)


def generate_visualizations(args: argparse.Namespace) -> dict[str, object]:
    _configure_matplotlib()
    output_dir = args.output_dir.resolve()
    logger = _configure_logger(output_dir)
    sampling_run = (
        args.sampling_run.resolve()
        if args.sampling_run is not None
        else _latest_sampling_run(DEFAULT_OUTPUT_ROOT)
    )
    rois = load_sampling_display_rois(sampling_run)
    if not rois:
        raise FileNotFoundError(f"No saved candidate ROI geometry found in {sampling_run}")
    source_model_ids = {roi.source_model_id for roi in rois}
    rodent_run = (
        args.rodent_run.resolve()
        if args.rodent_run is not None
        else _compatible_rodent_run(DEFAULT_OUTPUT_ROOT, source_model_ids)
    )
    summary = _read_json(sampling_run / "report" / "sampling_summary.json")
    roi = _choose_roi(rois, args.roi_id)
    model = _select_model(load_models_from_rodent_run(rodent_run), roi)
    preprocess_manifest, original_mask, original_swc = _load_original_input_data(
        rodent_run,
        model.source_model_id,
    )
    spacing_xyz_um = tuple(
        float(value) for value in preprocess_manifest["spacing_xyz_um"]  # type: ignore[index]
    )
    source_record = preprocess_manifest["record"]
    volume_shape_zyx = tuple(int(value) for value in original_mask.shape)
    _display_image, display_mask = load_normalized_volume(
        Path(str(preprocess_manifest["normalized_volume_path"]))
    )
    analysis_swc = load_normalized_swc(
        Path(str(preprocess_manifest["normalized_swc_path"])),
        source_path=Path(str(source_record["swc_path"])),  # type: ignore[index]
        spacing_xyz_um=spacing_xyz_um,  # type: ignore[arg-type]
        volume_shape_zyx=volume_shape_zyx,  # type: ignore[arg-type]
    )
    graph_root, hierarchy_branches, hierarchy_node_roles, source_parent_edges, graph_summary = (
        _load_hierarchy_artifacts(rodent_run, model.source_model_id)
    )
    hierarchy_audit = _hierarchy_audit(hierarchy_branches, source_parent_edges)
    preprocessing = dict(  # type: ignore[arg-type]
        preprocess_manifest["swc_centric_preprocessing"]
    )
    with Path(str(preprocess_manifest["swc_component_table_path"])).open(
        encoding="utf-8-sig",
        newline="",
    ) as stream:
        component_records = list(csv.DictReader(stream))
    mask_labels, mask_component_count = label_connected_components(
        np.asarray(original_mask) > 0,
        structure=np.ones((3, 3, 3), dtype=np.uint8),
    )
    full_mask_surface = _mask_surface(
        original_mask,
        spacing_xyz_um,  # type: ignore[arg-type]
        step_size=4,
    )
    swc_components = _swc_component_segments(original_swc)
    connectivity_checks = {
        "mask component count": (
            int(mask_component_count),
            int(preprocessing["mask_qc"]["component_count_26"]),  # type: ignore[index]
        ),
        "reference SWC component count": (
            len(swc_components),
            int(preprocessing["reference_component_count"]),
        ),
        "analysis SWC node count": (
            analysis_swc.node_count,
            int(preprocessing["analysis_node_count"]),
        ),
        "analysis SWC edge count": (
            analysis_swc.edge_count,
            int(preprocessing["analysis_edge_count"]),
        ),
        "display mask is unchanged": (
            int(np.array_equal(np.asarray(display_mask) > 0, np.asarray(original_mask) > 0)),
            1,
        ),
        "topology modification count": (
            int(preprocessing["new_node_count"])
            + int(preprocessing["new_edge_count"])
            + int(preprocessing["parent_relation_change_count"]),
            0,
        ),
    }
    mismatches = {
        name: {"reconstructed": actual, "manifest": expected}
        for name, (actual, expected) in connectivity_checks.items()
        if actual != expected
    }
    if mismatches:
        raise ValueError(f"Connectivity audit reconstruction disagrees with manifest: {mismatches}")
    hierarchy_checks = {
        "source node count": (analysis_swc.node_count, int(graph_summary["source_node_count"])),
        "source edge count": (len(source_parent_edges), int(graph_summary["source_edge_count"])),
        "critical node count": (len(hierarchy_node_roles), int(graph_summary["critical_node_count"])),
        "branch count": (len(hierarchy_branches), int(graph_summary["branch_count"])),
        "branch relation count": (int(hierarchy_audit["branch_relation_count"]), int(graph_summary["branch_relation_count"])),
        "missing source edges": (int(hierarchy_audit["missing_source_edge_count"]), 0),
        "unexpected or reversed edges": (int(hierarchy_audit["unexpected_or_reversed_edge_count"]), 0),
        "duplicate edge representations": (int(hierarchy_audit["duplicate_representation_count"]), 0),
        "branch hierarchy DAG": (int(bool(hierarchy_audit["branch_hierarchy_is_dag"])), 1),
    }
    hierarchy_mismatches = {
        name: {"reconstructed": actual, "exported": expected}
        for name, (actual, expected) in hierarchy_checks.items()
        if actual != expected
    }
    if hierarchy_mismatches:
        raise ValueError(
            "Hierarchy reconstruction disagrees with graph exports: "
            f"{hierarchy_mismatches}"
        )
    original_swc_segments = _swc_segments(original_swc)
    global_min = np.asarray(model.model_bounds_xyz_um[0::2], dtype=float)
    global_max = np.asarray(model.model_bounds_xyz_um[1::2], dtype=float)
    relationship_margin = np.asarray((10.0, 10.0, 16.0), dtype=float)
    relationship_min = np.maximum(np.asarray(roi.bbox_min_um) - relationship_margin, global_min)
    relationship_max = np.minimum(np.asarray(roi.bbox_max_um) + relationship_margin, global_max)
    crop_mask_surface = _mask_surface(
        original_mask,
        spacing_xyz_um,  # type: ignore[arg-type]
        step_size=2,
        crop_min_xyz_um=relationship_min,
        crop_max_xyz_um=relationship_max,
    )
    config = _sampling_config_from_summary(summary, sampling_run)
    candidate_anchor_ids, anchor_role_ids, anchor_strategy_results = _anchor_strategy_results(
        model,
        config,
    )
    anchor_process_trace = _anchor_sampling_process_trace(
        model,
        config,
        candidate_anchor_ids,
        anchor_strategy_results,
    )
    all_anchor_ids = anchor_strategy_results["farthest_point"].anchor_ids
    if all_anchor_ids != generate_anchor_ids(model, config):
        raise ValueError("Anchor strategy audit disagrees with the configured sampling mode")
    raw_segments, raw_component_count = _raw_clipped_segments(model, roi)
    if raw_component_count != roi.raw_component_count:
        logger.warning(
            "Reconstructed raw component count differs from manifest: %s != %s",
            raw_component_count,
            roi.raw_component_count,
        )
    logger.info("Rodent run: %s", rodent_run)
    logger.info("Sampling run: %s", sampling_run)
    logger.info(
        "Original input concepts: mask foreground voxels=%s; SWC nodes=%s, edges=%s, components=%s",
        np.count_nonzero(original_mask),
        original_swc.node_count,
        original_swc.edge_count,
        original_swc.component_count,
    )
    logger.info(
        "SWC-centric preprocessing: reference components=%s, nodes=%s; analysis component=%s, "
        "nodes=%s, edges=%s; topology modifications=0; original Mask components=%s and unchanged",
        original_swc.component_count,
        original_swc.node_count,
        preprocessing["selected_component_id"],
        analysis_swc.node_count,
        analysis_swc.edge_count,
        mask_component_count,
    )
    logger.info(
        "Selected example ROI: %s; anchor=%s; cluster=%s; representative_rank=%s",
        roi.roi_id,
        roi.anchor_id,
        roi.cluster_id,
        roi.selection_rank,
    )
    logger.info(
        "Actual counts: global nodes=%s, edges=%s, anchors=%s, valid ROIs=%s, clusters=%s, representatives=%s",
        summary["global_node_count"],
        summary["global_edge_count"],
        len(all_anchor_ids),
        len(rois),
        summary["n_clusters"],
        summary["selected_count"],
    )
    logger.info(
        "3-D anchor strategy audit: candidate pool=%s (root=%s, terminal=%s, bifurcation=%s); "
        "random=%s, farthest_point=%s, poisson_disk=%s; random_equals_poisson=%s",
        len(candidate_anchor_ids),
        len(anchor_role_ids["root"]),
        len(anchor_role_ids["terminal"]),
        len(anchor_role_ids["bifurcation"]),
        len(anchor_strategy_results["random"].anchor_ids),
        len(anchor_strategy_results["farthest_point"].anchor_ids),
        len(anchor_strategy_results["poisson_disk"].anchor_ids),
        anchor_strategy_results["random"].anchor_ids
        == anchor_strategy_results["poisson_disk"].anchor_ids,
    )
    visualized_random_events = anchor_process_trace["random_events"]  # type: ignore[assignment]
    visualized_farthest_events = anchor_process_trace["farthest_events"]  # type: ignore[assignment]
    logger.info(
        "3-D anchor process figure: random events 1..6 visualized (event 4 rejects node %s at %.4f um); "
        "farthest selections 1..3 visualized (%s); Poisson panels replay random events 1, 2, and 4 "
        "with %.1f um exclusion spheres",
        visualized_random_events[3]["candidate_id"],
        visualized_random_events[3]["nearest_selected_distance_um"],
        [event["candidate_id"] for event in visualized_farthest_events[:3]],
        config.min_anchor_distance_um,
    )
    represented_cluster_ids = sorted(
        {candidate.cluster_id for candidate in rois if candidate.is_representative}
    )
    logger.info(
        "3-D representative ROI audit: example raw components=%s, retained nodes=%s, "
        "edges=%s, branches=%s, retained length fraction=%.4f; represented clusters=%s/%s (%s)",
        raw_component_count,
        roi.node_count,
        roi.edge_count,
        roi.branch_count,
        roi.retained_component_length_um / roi.raw_total_vessel_length_um,
        len(represented_cluster_ids),
        summary["n_clusters"],
        represented_cluster_ids,
    )
    logger.info(
        "Hierarchy graph: critical nodes=%s, branches=%s, relations=%s, depth=0..%s, "
        "Strahler=1..%s, Horsfield=1..%s; source-edge audit PASS",
        graph_summary["critical_node_count"],
        graph_summary["branch_count"],
        graph_summary["branch_relation_count"],
        max(branch.depth for branch in hierarchy_branches),
        max(branch.strahler_order for branch in hierarchy_branches),
        max(branch.horsfield_order for branch in hierarchy_branches),
    )

    connectivity_3d = _figure_connectivity_preprocessing_3d(
        original_mask,
        full_mask_surface,
        original_swc,
        analysis_swc,
        swc_components,
        component_records,
        preprocessing,
        spacing_xyz_um,  # type: ignore[arg-type]
        output_dir / "00_connectivity_preprocessing_3d.png",
        args.dpi,
    )
    connectivity_audit = _figure_connectivity_preprocessing_audit(
        component_records,
        preprocessing,
        output_dir / "00_connectivity_preprocessing_audit.png",
        args.dpi,
    )
    input_concepts = _figure_input_concepts_3d(
        original_mask,
        full_mask_surface,
        original_swc,
        original_swc_segments,
        model,
        roi,
        spacing_xyz_um,  # type: ignore[arg-type]
        output_dir / "00_input_concepts_3d.png",
        args.dpi,
    )
    input_relationship = _figure_input_relationship_3d(
        crop_mask_surface,
        original_swc_segments,
        model,
        roi,
        relationship_min,
        relationship_max,
        output_dir / "00_input_spatial_relationship_3d.png",
        args.dpi,
    )
    hierarchy_construction, hierarchy_focus_ids = _figure_hierarchical_graph_construction_3d(
        analysis_swc,
        hierarchy_branches,
        hierarchy_node_roles,
        graph_summary,
        output_dir / "01_hierarchical_graph_construction_3d.png",
        args.dpi,
    )
    depth_3d, depth_3d_metadata = _figure_branch_depth_3d(
        analysis_swc,
        hierarchy_branches,
        output_dir / "01_branch_depth_concept_3d.png",
        args.dpi,
    )
    strahler_3d, strahler_3d_metadata = _figure_strahler_order_3d(
        analysis_swc,
        hierarchy_branches,
        output_dir / "01_strahler_order_concept_3d.png",
        args.dpi,
    )
    horsfield_3d, horsfield_3d_metadata = _figure_horsfield_order_3d(
        analysis_swc,
        hierarchy_branches,
        output_dir / "01_horsfield_order_concept_3d.png",
        args.dpi,
    )
    hierarchy_audit_figure = _figure_hierarchical_graph_audit(
        hierarchy_branches,
        hierarchy_node_roles,
        output_dir / "01_hierarchical_graph_audit.png",
        args.dpi,
    )
    logger.info(
        "3-D hierarchy examples: depth path=%s; Strahler root=B%s daughters=%s; "
        "Horsfield root=B%s daughters=%s",
        depth_3d_metadata["deepest_path_branch_ids"],
        strahler_3d_metadata["example_root_branch_id"],
        strahler_3d_metadata["example_daughter_branch_ids"],
        horsfield_3d_metadata["example_root_branch_id"],
        horsfield_3d_metadata["example_daughter_branch_ids"],
    )
    anchor_sampling_3d = _figure_anchor_sampling_strategies_3d(
        model,
        config,
        candidate_anchor_ids,
        anchor_role_ids,
        anchor_strategy_results,
        anchor_process_trace,
        output_dir / "02_anchor_sampling_strategies_3d.png",
        args.dpi,
    )
    overview = _figure_process_overview(
        summary,
        roi,
        raw_component_count,
        output_dir / "01_representative_roi_sampling_overview.png",
        args.dpi,
    )
    representative_sampling_3d = _figure_representative_roi_sampling_3d(
        model,
        rois,
        all_anchor_ids,
        roi,
        raw_segments,
        raw_component_count,
        summary,
        output_dir / "02_representative_roi_sampling_3d.png",
        args.dpi,
    )
    extraction = _figure_extraction_trace(
        model,
        rois,
        all_anchor_ids,
        roi,
        raw_segments,
        raw_component_count,
        output_dir / "02_real_roi_extraction_trace.png",
        args.dpi,
    )
    features = _figure_features_and_selection(
        rois,
        roi,
        sampling_run,
        summary,
        output_dir / "03_features_clustering_and_representatives.png",
        args.dpi,
    )
    output_files = [
        connectivity_3d,
        connectivity_audit,
        input_concepts,
        input_relationship,
        hierarchy_construction,
        depth_3d,
        strahler_3d,
        horsfield_3d,
        hierarchy_audit_figure,
        anchor_sampling_3d,
        overview,
        representative_sampling_3d,
        extraction,
        features,
    ]
    record = source_record
    manifest: dict[str, object] = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "rodent_run": str(rodent_run),
        "sampling_run": str(sampling_run),
        "source_model_id": model.source_model_id,
        "example_roi_id": roi.roi_id,
        "example_anchor_id": roi.anchor_id,
        "example_cluster_id": roi.cluster_id,
        "example_selection_rank": roi.selection_rank,
        "example_raw_component_count_manifest": roi.raw_component_count,
        "example_raw_component_count_reconstructed": raw_component_count,
        "global_node_count": int(summary["global_node_count"]),
        "global_edge_count": int(summary["global_edge_count"]),
        "candidate_anchor_count": len(all_anchor_ids),
        "valid_candidate_roi_count": len(rois),
        "cluster_count": int(summary["n_clusters"]),
        "representative_roi_count": sum(roi.is_representative for roi in rois),
        "connectivity_preprocessing": preprocessing
        | {
            "visualized_mask_component_count_reconstructed": int(mask_component_count),
            "visualized_reference_edge_count": original_swc.edge_count,
            "visualized_analysis_edge_count": analysis_swc.edge_count,
            "component_audit_source": str(preprocess_manifest["swc_component_table_path"]),
            "reference_only_nodes_source": str(
                preprocess_manifest["reference_only_node_ids_path"]
            ),
            "downstream_source": "unchanged selected analysis_swc component",
        },
        "input_concepts_3d": {
            "mask_path": str(record["mask_path"]),  # type: ignore[index]
            "mask_source_stage": "original segmentation mask; QC and visualization only",
            "mask_shape_zyx": list(map(int, original_mask.shape)),
            "mask_foreground_voxel_count": int(np.count_nonzero(original_mask)),
            "mask_surface_display_step_size": 4,
            "mask_surface_vertex_count": len(full_mask_surface.vertices_xyz_um),
            "mask_surface_face_count": len(full_mask_surface.faces),
            "swc_path": str(record["swc_path"]),  # type: ignore[index]
            "swc_source_stage": "complete human-corrected reference_swc",
            "swc_node_count": original_swc.node_count,
            "swc_edge_count": original_swc.edge_count,
            "swc_component_count": original_swc.component_count,
            "roi_source_stage": "saved connected ROI extracted from unchanged analysis_swc",
            "roi_id": roi.roi_id,
            "roi_bbox_min_xyz_um": list(map(float, roi.bbox_min_um)),
            "roi_bbox_max_xyz_um": list(map(float, roi.bbox_max_um)),
            "coordinate_order": "x, y, z",
            "coordinate_unit": "um",
        },
        "hierarchical_graph_construction": graph_summary
        | hierarchy_audit
        | {
            "graph_artifact_root": str(graph_root),
            "visualized_root_subtree_branch_ids": list(hierarchy_focus_ids),
            "visualized_root_subtree_relative_depth": [0, 3],
            "branch_depth_range": [
                min(branch.depth for branch in hierarchy_branches),
                max(branch.depth for branch in hierarchy_branches),
            ],
            "strahler_order_range": [
                min(branch.strahler_order for branch in hierarchy_branches),
                max(branch.strahler_order for branch in hierarchy_branches),
            ],
            "horsfield_order_range": [
                min(branch.horsfield_order for branch in hierarchy_branches),
                max(branch.horsfield_order for branch in hierarchy_branches),
            ],
            "direction_interpretation": "structural SWC parent-to-current only; not measured flow",
        },
        "hierarchy_metric_examples_3d": {
            "branch_depth": depth_3d_metadata,
            "strahler_order": strahler_3d_metadata,
            "horsfield_order": horsfield_3d_metadata,
        },
        "anchor_sampling_process_3d": {
            "source_stage": "unchanged analysis_swc global nodes and edges",
            "visualization_focus": "real intermediate decisions, not final anchor distributions",
            "candidate_pool_rule": (
                "nodes inside the ROI-center feasible bounds with global undirected degree != 2; "
                "fallbacks are not used for this dataset"
            ),
            "coordinate_order": "x, y, z",
            "coordinate_unit": "um",
            "seed": config.seed,
            "roi_size_xyz_um": list(map(float, config.roi_size_um)),
            "minimum_anchor_distance_um": config.min_anchor_distance_um,
            "maximum_candidate_anchors": config.max_candidate_anchors,
            "candidate_pool_count": len(candidate_anchor_ids),
            "candidate_role_counts": {
                role: len(node_ids) for role, node_ids in anchor_role_ids.items()
            },
            "strategies": {
                mode: {
                    "anchor_count": len(result.anchor_ids),
                    "anchor_ids": list(result.anchor_ids),
                    "closest_pair_ids": list(result.closest_pair_ids),
                    "minimum_pairwise_distance_um": result.minimum_pairwise_distance_um,
                    "median_nearest_neighbor_distance_um": result.median_nearest_neighbor_distance_um,
                    "mean_nearest_neighbor_distance_um": result.mean_nearest_neighbor_distance_um,
                }
                for mode, result in anchor_strategy_results.items()
            },
            "random_equals_poisson_disk_for_current_implementation_and_seed": (
                anchor_strategy_results["random"].anchor_ids
                == anchor_strategy_results["poisson_disk"].anchor_ids
            ),
            "random_farthest_shared_anchor_count": len(
                set(anchor_strategy_results["random"].anchor_ids)
                & set(anchor_strategy_results["farthest_point"].anchor_ids)
            ),
            "visualized_process_trace": {
                "candidate_centroid_xyz_um": list(
                    anchor_process_trace["candidate_centroid_xyz_um"]  # type: ignore[arg-type]
                ),
                "random_first_six_events": list(
                    anchor_process_trace["random_events"][:6]  # type: ignore[index]
                ),
                "farthest_first_three_events": list(
                    anchor_process_trace["farthest_events"][:3]  # type: ignore[index]
                ),
                "poisson_explained_with_random_event_indices": [1, 2, 4],
            },
        },
        "representative_roi_sampling_3d": {
            "source_stage": "saved connected ROI records extracted from unchanged analysis_swc",
            "coordinate_order": "x, y, z",
            "coordinate_unit": "um",
            "candidate_anchor_count": len(all_anchor_ids),
            "valid_candidate_roi_count": len(rois),
            "rejected_anchor_count": len(all_anchor_ids) - len(rois),
            "cluster_count": int(summary["n_clusters"]),
            "represented_cluster_ids": represented_cluster_ids,
            "represented_cluster_count": len(represented_cluster_ids),
            "representative_roi_count": len(
                [candidate for candidate in rois if candidate.is_representative]
            ),
            "example_roi_id": roi.roi_id,
            "example_anchor_id": roi.anchor_id,
            "example_bbox_min_xyz_um": list(map(float, roi.bbox_min_um)),
            "example_bbox_max_xyz_um": list(map(float, roi.bbox_max_um)),
            "example_raw_clipped_segment_count": len(raw_segments),
            "example_raw_component_count": raw_component_count,
            "example_retained_node_count": roi.node_count,
            "example_retained_edge_count": roi.edge_count,
            "example_retained_branch_count": roi.branch_count,
            "example_retained_bifurcation_count": roi.bifurcation_count,
            "example_cut_port_count": roi.cut_port_count,
            "example_true_terminal_count": roi.true_terminal_count,
            "example_raw_total_vessel_length_um": roi.raw_total_vessel_length_um,
            "example_retained_component_length_um": roi.retained_component_length_um,
            "example_retained_length_fraction": (
                roi.retained_component_length_um / roi.raw_total_vessel_length_um
            ),
            "panel_meanings": {
                "A": "real anchors on the complete analysis_swc",
                "B": "all saved valid candidate ROI boxes at their real positions",
                "C": "exact line-box clipping and raw connected components",
                "D": "retained anchor-containing connected component",
                "E": "local-edge to global-edge provenance",
                "F": "saved representative ROI geometries at their original positions",
            },
        },
        "figures": [str(path) for path in output_files],
        "scientific_scope": {
            "roi_geometry_is_real": True,
            "local_to_global_mapping_preserved": True,
            "physiological_inlet_outlet_inference": False,
            "flow_or_microbubble_simulation": False,
        },
    }
    manifest_path = output_dir / "visualization_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    logger.info("Saved figures: %s", ", ".join(str(path) for path in output_files))
    logger.info("Saved manifest: %s", manifest_path)
    if args.show:
        plt.show()
    else:
        plt.close("all")
    return manifest


def main() -> int:
    args = parse_args()
    if args.dpi < 72:
        raise ValueError("--dpi must be at least 72")
    generate_visualizations(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
