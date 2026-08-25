"""Visual explanation of the shared-domain SWC-to-STL reconstruction idea.

The script uses the saved ROI that is also consumed by
``swc_stl_model_generate.py``.  It builds an *explanatory* local occupancy grid
from the real centreline and feed radii, then compares that grid with the
accepted Ultraliser surface.  The local rasterisation is deliberately kept
separate from the production backend: it illustrates the shared-domain union
but is not claimed to be an exported Ultraliser intermediate volume.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from pathlib import Path
from typing import Iterable

import matplotlib.image as mpimg
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.font_manager import FontProperties
from matplotlib.patches import Patch
import networkx as nx
import numpy as np
import pyvista as pv
from scipy import ndimage

from utils.cfd_lumen import load_sampling_rois, load_swc_stl_yaml_config, resolve_sampling_run
from utils.sampling.structural_features import _branch_paths


PROJECT_ROOT = Path(__file__).resolve().parent
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "swc_stl_model_generate.yaml"
DEFAULT_SURFACE = (
    PROJECT_ROOT
    / "outputs"
    / "reference_validation"
    / "roi003274_ultraliser_radius091"
    / "geometry"
    / "lumen_surface_um.vtp"
)
DEFAULT_OUTPUT_DIR = (
    PROJECT_ROOT / "outputs" / "documentation" / "swc_to_stl_shared_voxel_domain"
)
DEFAULT_PROCESSING_OUTPUT_DIR = (
    PROJECT_ROOT / "outputs" / "documentation" / "swc_to_stl_processing_stages"
)
DEFAULT_VOXELIZATION_OUTPUT_DIR = (
    PROJECT_ROOT / "outputs" / "documentation" / "swc_to_stl_voxelization_explained"
)
DEFAULT_RADIUS_VALIDATION_OUTPUT_DIR = (
    PROJECT_ROOT
    / "outputs"
    / "reference_validation"
    / "roi003274_ultraliser_radius091"
    / "figures"
)
DEFAULT_ULTRALISER_ROOT = PROJECT_ROOT / "Ultraliser"

BRANCH_COLORS = (
    "#0072B2",
    "#E69F00",
    "#009E73",
    "#CC79A7",
    "#56B4E9",
    "#D55E00",
)
BACKGROUND = "#F5F5F2"
JUNCTION_COLOR = "#7A0177"
CUT_PORT_COLOR = "#F0C419"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Visualise how saved SWC centreline/radius data form one shared "
            "occupancy domain and one accepted triangle surface."
        )
    )
    parser.add_argument("--config", type=Path, default=DEFAULT_CONFIG)
    parser.add_argument(
        "--sampling-run",
        type=Path,
        default=None,
        help="Saved ROI sampling run; defaults to the latest complete run.",
    )
    parser.add_argument("--roi-anchor", type=int, default=3274)
    parser.add_argument(
        "--surface",
        type=Path,
        default=DEFAULT_SURFACE,
        help="Accepted Ultraliser surface in micrometres.",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--processing-output-dir",
        type=Path,
        default=DEFAULT_PROCESSING_OUTPUT_DIR,
        help="Output directory for the detailed voxel/DMC/mesh-processing explanation.",
    )
    parser.add_argument(
        "--voxelization-output-dir",
        type=Path,
        default=DEFAULT_VOXELIZATION_OUTPUT_DIR,
        help=(
            "Output directory for the focused centreline/sphere/grid/shell/xyz-fill "
            "explanation of section 4.4.5."
        ),
    )
    parser.add_argument(
        "--ultraliser-root",
        type=Path,
        default=DEFAULT_ULTRALISER_ROOT,
        help=(
            "Official Ultraliser checkout. It is invoked once only when the official "
            "DMC/Laplacian/optimized stage meshes have not yet been captured."
        ),
    )
    parser.add_argument(
        "--junction-node",
        type=int,
        default=None,
        help="Local junction node. By default, use the junction furthest from the ROI boundary.",
    )
    parser.add_argument(
        "--crop-half-extent-um",
        type=float,
        default=7.5,
        help="Half-width of the local cubic explanation region.",
    )
    parser.add_argument(
        "--panel-width",
        type=int,
        default=1100,
        help="Width in pixels of each rendered 3-D panel.",
    )
    parser.add_argument(
        "--panel-height",
        type=int,
        default=820,
        help="Height in pixels of each rendered 3-D panel.",
    )
    parser.add_argument(
        "--radius-baseline-run",
        type=Path,
        default=None,
        help="Model run reconstructed with radius_scale=1.0.",
    )
    parser.add_argument(
        "--radius-compensated-run",
        type=Path,
        default=None,
        help="Model run reconstructed with the accepted compensated radius scale.",
    )
    parser.add_argument(
        "--radius-validation-output-dir",
        type=Path,
        default=DEFAULT_RADIUS_VALIDATION_OUTPUT_DIR,
        help="Output directory for the radius-before/after documentation figures.",
    )
    parser.add_argument(
        "--radius-comparison-only",
        action="store_true",
        help="Only regenerate the two radius-validation figures.",
    )
    return parser


def _build_graph(roi) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(int(value) for value in roi.local_node_ids)
    graph.add_edges_from((int(a), int(b)) for a, b in roi.local_edges)
    if not nx.is_connected(graph):
        raise ValueError(f"ROI {roi.roi_id} is not one connected graph")
    return graph


def _select_roi(sampling_run: Path, anchor: int):
    matches = [
        roi
        for roi in load_sampling_rois(sampling_run, selected_only=True)
        if int(roi.anchor_id) == int(anchor)
    ]
    if len(matches) != 1:
        raise ValueError(f"Expected one selected ROI for anchor {anchor}, found {len(matches)}")
    return matches[0]


def _select_junction(roi, graph: nx.Graph, requested: int | None) -> int:
    junctions = [int(node) for node in graph if graph.degree(node) >= 3]
    if not junctions:
        raise ValueError(f"ROI {roi.roi_id} has no degree-three-or-higher junction")
    if requested is not None:
        if requested not in junctions:
            raise ValueError(
                f"Requested node {requested} is not a junction; available junctions: {junctions}"
            )
        return int(requested)

    bbox_min = np.asarray(roi.bbox_min_um, dtype=float)
    bbox_max = np.asarray(roi.bbox_max_um, dtype=float)

    def boundary_clearance(node: int) -> tuple[float, float]:
        point = np.asarray(roi.local_node_positions_um[node], dtype=float)
        clearance = float(np.min(np.r_[point - bbox_min, bbox_max - point]))
        return clearance, float(roi.local_node_radius_um[node])

    return max(junctions, key=boundary_clearance)


def _polyline(points: np.ndarray, radii_um: np.ndarray | None = None) -> pv.PolyData:
    points = np.asarray(points, dtype=float)
    if len(points) < 2:
        raise ValueError("A polyline needs at least two points")
    lines = np.r_[len(points), np.arange(len(points), dtype=np.int64)]
    mesh = pv.PolyData(points, lines=lines)
    if radii_um is not None:
        mesh.point_data["feed_radius_um"] = np.asarray(radii_um, dtype=float)
    return mesh


def _centerline_mesh(roi, graph: nx.Graph, radius_scale: float) -> pv.PolyData:
    edges = np.asarray(roi.local_edges, dtype=np.int64)
    lines = np.column_stack((np.full(len(edges), 2, dtype=np.int64), edges)).ravel()
    mesh = pv.PolyData(np.asarray(roi.local_node_positions_um, dtype=float), lines=lines)
    mesh.point_data["local_node_id"] = np.asarray(roi.local_node_ids, dtype=np.int64)
    mesh.point_data["source_radius_um"] = np.asarray(roi.local_node_radius_um, dtype=float)
    mesh.point_data["feed_radius_um"] = (
        np.asarray(roi.local_node_radius_um, dtype=float) * float(radius_scale)
    )
    mesh.point_data["degree"] = np.asarray(
        [graph.degree(int(node)) for node in roi.local_node_ids], dtype=np.int64
    )
    return mesh


def _crop_bounds(center_um: np.ndarray, half_extent_um: float) -> tuple[float, ...]:
    low = np.asarray(center_um, dtype=float) - float(half_extent_um)
    high = np.asarray(center_um, dtype=float) + float(half_extent_um)
    return (
        float(low[0]),
        float(high[0]),
        float(low[1]),
        float(high[1]),
        float(low[2]),
        float(high[2]),
    )


def _capsule_union_occupancy(
    edge_points_um: np.ndarray,
    edge_radii_um: np.ndarray,
    bounds: tuple[float, ...],
    voxels_per_micron: float,
) -> tuple[pv.ImageData, pv.UnstructuredGrid, pv.PolyData, dict[str, object]]:
    """Rasterise variable-radius edge capsules into one logical union.

    This is a compact teaching implementation.  It shares the production
    input units and spacing, but it does not reproduce or replace Ultraliser's
    C++ ``polylines-with-spheres`` implementation.
    """

    spacing = 1.0 / float(voxels_per_micron)
    requested_low = np.asarray([bounds[0], bounds[2], bounds[4]], dtype=float)
    requested_high = np.asarray([bounds[1], bounds[3], bounds[5]], dtype=float)
    dimensions = np.maximum(
        1, np.rint((requested_high - requested_low) / spacing).astype(np.int64)
    )
    size = dimensions.astype(float) * spacing
    center = 0.5 * (requested_low + requested_high)
    origin = center - 0.5 * size
    actual_high = origin + size
    mask = np.zeros(tuple(int(value) for value in dimensions), dtype=bool)

    used_edges = 0
    for endpoints, radii in zip(edge_points_um, edge_radii_um, strict=True):
        p0, p1 = np.asarray(endpoints, dtype=float)
        r0, r1 = np.asarray(radii, dtype=float)
        maximum_radius = float(max(r0, r1))
        edge_low = np.minimum(p0, p1) - maximum_radius
        edge_high = np.maximum(p0, p1) + maximum_radius
        if np.any(edge_high < origin) or np.any(edge_low > actual_high):
            continue

        index_low = np.maximum(
            0, np.floor((edge_low - origin) / spacing - 0.5).astype(np.int64)
        )
        index_high = np.minimum(
            dimensions - 1,
            np.ceil((edge_high - origin) / spacing - 0.5).astype(np.int64),
        )
        if np.any(index_low > index_high):
            continue
        axes = [
            origin[axis]
            + (np.arange(index_low[axis], index_high[axis] + 1, dtype=float) + 0.5)
            * spacing
            for axis in range(3)
        ]
        xx, yy, zz = np.meshgrid(*axes, indexing="ij")
        samples = np.stack((xx, yy, zz), axis=-1)
        direction = p1 - p0
        squared_length = float(np.dot(direction, direction))
        if squared_length <= np.finfo(float).eps:
            parameter = np.zeros(samples.shape[:-1], dtype=float)
        else:
            parameter = np.clip(
                np.einsum("...i,i->...", samples - p0, direction) / squared_length,
                0.0,
                1.0,
            )
        closest = p0 + parameter[..., None] * direction
        local_radius = r0 + parameter * (r1 - r0)
        inside = np.linalg.norm(samples - closest, axis=-1) <= local_radius
        slices = tuple(
            slice(int(index_low[axis]), int(index_high[axis]) + 1) for axis in range(3)
        )
        mask[slices] |= inside
        used_edges += 1

    grid = pv.ImageData(
        dimensions=tuple(int(value) + 1 for value in dimensions),
        spacing=(spacing, spacing, spacing),
        origin=tuple(float(value) for value in origin),
    )
    grid.cell_data["occupied"] = mask.astype(np.uint8).ravel(order="F")
    occupied = grid.threshold(0.5, scalars="occupied")
    boundary = occupied.extract_surface(algorithm="dataset_surface").clean()
    details: dict[str, object] = {
        "spacing_um": spacing,
        "dimensions_cells": [int(value) for value in dimensions],
        "origin_um": [float(value) for value in origin],
        "bounds_um": [
            [float(value) for value in origin],
            [float(value) for value in actual_high],
        ],
        "total_cell_count": int(mask.size),
        "occupied_cell_count": int(np.count_nonzero(mask)),
        "edge_count_contributing_to_crop": int(used_edges),
    }
    return grid, occupied, boundary, details


def _font() -> FontProperties:
    candidates = (
        Path("C:/Windows/Fonts/msyh.ttc"),
        Path("C:/Windows/Fonts/msyhbd.ttc"),
        Path("C:/Windows/Fonts/simhei.ttf"),
    )
    for path in candidates:
        if path.is_file():
            return FontProperties(fname=str(path))
    return FontProperties(family="DejaVu Sans")


def _camera(center: np.ndarray, half_extent: float) -> tuple[list[float], list[float], list[float]]:
    direction = np.asarray([1.55, -2.0, 1.25], dtype=float)
    direction /= np.linalg.norm(direction)
    position = np.asarray(center, dtype=float) + direction * float(half_extent) * 4.2
    return position.tolist(), np.asarray(center, dtype=float).tolist(), [0.0, 0.0, 1.0]


def _set_local_camera(plotter: pv.Plotter, center: np.ndarray, half_extent: float) -> None:
    plotter.camera_position = _camera(center, half_extent)
    plotter.camera.parallel_projection = True
    plotter.camera.zoom(1.18)


def _add_axes(plotter: pv.Plotter, bounds: tuple[float, ...]) -> None:
    plotter.show_bounds(
        bounds=bounds,
        grid="front",
        location="outer",
        all_edges=True,
        xtitle="x (um)",
        ytitle="y (um)",
        ztitle="z (um)",
        color="#474747",
        font_size=10,
    )


def _new_plotter(window_size: tuple[int, int]) -> pv.Plotter:
    plotter = pv.Plotter(off_screen=True, window_size=window_size)
    plotter.set_background(BACKGROUND)
    plotter.enable_anti_aliasing("ssaa")
    return plotter


def _save_panel(plotter: pv.Plotter, path: Path) -> None:
    plotter.screenshot(str(path), return_img=False)
    plotter.close()


def _render_overview_panel(
    roi,
    graph: nx.Graph,
    paths: list[list[int]],
    centerline: pv.PolyData,
    junction_node: int,
    local_bounds: tuple[float, ...],
    output_path: Path,
    window_size: tuple[int, int],
) -> None:
    plotter = _new_plotter(window_size)
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    for index, path in enumerate(paths):
        plotter.add_mesh(
            _polyline(positions[path]),
            color=BRANCH_COLORS[index % len(BRANCH_COLORS)],
            line_width=7,
            render_lines_as_tubes=True,
        )

    points = pv.PolyData(positions)
    points.point_data["feed_radius_um"] = centerline.point_data["feed_radius_um"]
    sphere = pv.Sphere(radius=1.0, theta_resolution=18, phi_resolution=18)
    glyphs = points.glyph(
        scale="feed_radius_um", orient=False, factor=1.0, geom=sphere
    )
    plotter.add_mesh(
        glyphs,
        color="#78C6E7",
        opacity=0.27,
        smooth_shading=True,
    )
    port_nodes = [int(port.local_node_id) for port in roi.cut_ports]
    if port_nodes:
        plotter.add_points(
            positions[port_nodes],
            color=CUT_PORT_COLOR,
            point_size=15,
            render_points_as_spheres=True,
        )
    junction_position = positions[junction_node]
    plotter.add_points(
        junction_position[None, :],
        color=JUNCTION_COLOR,
        point_size=19,
        render_points_as_spheres=True,
    )
    plotter.add_point_labels(
        junction_position[None, :],
        [f"local node {junction_node}"],
        point_size=0,
        font_size=15,
        text_color=JUNCTION_COLOR,
        shape_color="white",
        shape_opacity=0.82,
    )
    plotter.add_mesh(
        pv.Box(bounds=local_bounds),
        style="wireframe",
        color=JUNCTION_COLOR,
        line_width=4,
    )
    roi_center = 0.5 * (
        np.asarray(roi.bbox_min_um, dtype=float) + np.asarray(roi.bbox_max_um, dtype=float)
    )
    roi_span = float(np.max(np.asarray(roi.bbox_size_um, dtype=float)))
    plotter.camera_position = _camera(roi_center, roi_span * 0.50)
    plotter.camera.parallel_projection = True
    plotter.camera.parallel_scale = roi_span * 0.62
    full_bounds = (
        float(roi.bbox_min_um[0]),
        float(roi.bbox_max_um[0]),
        float(roi.bbox_min_um[1]),
        float(roi.bbox_max_um[1]),
        float(roi.bbox_min_um[2]),
        float(roi.bbox_max_um[2]),
    )
    _add_axes(plotter, full_bounds)
    _save_panel(plotter, output_path)


def _incident_paths(paths: Iterable[list[int]], junction_node: int) -> list[list[int]]:
    return [list(path) for path in paths if int(junction_node) in path]


def _render_separate_branches_panel(
    roi,
    incident_paths: list[list[int]],
    radius_scale: float,
    junction_node: int,
    local_bounds: tuple[float, ...],
    output_path: Path,
    window_size: tuple[int, int],
    half_extent: float,
) -> None:
    plotter = _new_plotter(window_size)
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float) * float(radius_scale)
    crop_box = pv.Box(bounds=local_bounds)
    for index, path in enumerate(incident_paths):
        line = _polyline(positions[path], radii[path])
        tube = line.tube(
            radius=None,
            scalars="feed_radius_um",
            absolute=True,
            capping=True,
            n_sides=28,
        )
        clipped = (
            tube.clip_box(crop_box, invert=False)
            .extract_surface(algorithm="dataset_surface")
            .triangulate()
            .clean()
        )
        plotter.add_mesh(
            clipped,
            color=BRANCH_COLORS[index % len(BRANCH_COLORS)],
            opacity=0.78,
            show_edges=True,
            edge_color="#414141",
            edge_opacity=0.20,
            smooth_shading=False,
            lighting=True,
        )
        plotter.add_mesh(
            line.clip_box(crop_box, invert=False),
            color="#202020",
            line_width=4,
            render_lines_as_tubes=True,
        )
    center = positions[junction_node]
    plotter.add_points(
        center[None, :],
        color=JUNCTION_COLOR,
        point_size=18,
        render_points_as_spheres=True,
    )
    _add_axes(plotter, local_bounds)
    _set_local_camera(plotter, center, half_extent)
    _save_panel(plotter, output_path)


def _render_occupancy_panel(
    boundary: pv.PolyData,
    roi,
    graph: nx.Graph,
    junction_node: int,
    local_bounds: tuple[float, ...],
    output_path: Path,
    window_size: tuple[int, int],
    half_extent: float,
) -> None:
    plotter = _new_plotter(window_size)
    plotter.add_mesh(
        boundary,
        color="#32B6C5",
        opacity=0.86,
        show_edges=True,
        edge_color="#174F56",
        edge_opacity=0.24,
        lighting=True,
    )
    centerline = _centerline_mesh(roi, graph, radius_scale=1.0)
    cropped_line = centerline.clip_box(pv.Box(bounds=local_bounds), invert=False)
    plotter.add_mesh(
        cropped_line,
        color="white",
        line_width=5,
        render_lines_as_tubes=True,
    )
    center = np.asarray(roi.local_node_positions_um[junction_node], dtype=float)
    plotter.add_points(
        center[None, :],
        color=JUNCTION_COLOR,
        point_size=17,
        render_points_as_spheres=True,
    )
    _add_axes(plotter, local_bounds)
    _set_local_camera(plotter, center, half_extent)
    _save_panel(plotter, output_path)


def _render_surface_panel(
    local_surface: pv.PolyData,
    roi,
    graph: nx.Graph,
    junction_node: int,
    local_bounds: tuple[float, ...],
    output_path: Path,
    window_size: tuple[int, int],
    half_extent: float,
) -> None:
    plotter = _new_plotter(window_size)
    plotter.add_mesh(
        local_surface,
        color="#C83C32",
        show_edges=True,
        edge_color="#43201E",
        edge_opacity=0.33,
        smooth_shading=True,
        specular=0.20,
        specular_power=18,
    )
    centerline = _centerline_mesh(roi, graph, radius_scale=1.0)
    cropped_line = centerline.clip_box(pv.Box(bounds=local_bounds), invert=False)
    plotter.add_mesh(
        cropped_line,
        color="#FFE35B",
        line_width=5,
        render_lines_as_tubes=True,
    )
    center = np.asarray(roi.local_node_positions_um[junction_node], dtype=float)
    plotter.add_points(
        center[None, :],
        color=JUNCTION_COLOR,
        point_size=17,
        render_points_as_spheres=True,
    )
    _add_axes(plotter, local_bounds)
    _set_local_camera(plotter, center, half_extent)
    _save_panel(plotter, output_path)


def _compose_figure(
    panel_paths: list[Path],
    output_path: Path,
    *,
    roi,
    branch_count: int,
    junction_node: int,
    radius_scale: float,
    voxels_per_micron: float,
    surface_cells: int,
) -> None:
    font = _font()
    titles = (
        (
            "A  实际保存的 ROI：中心线与局部半径",
            f"{len(roi.local_node_ids)} 个节点，{len(roi.local_edges)} 条边，{branch_count} 条拓扑分支；紫框为局部示例",
        ),
        (
            "B  只把三条分支分别建成管面（用于对比）",
            "颜色代表三个独立表面；它们在同一分叉位置相交，但仍是表面的并列集合",
        ),
        (
            "C  所有分支写入同一个占据域",
            f"实际 feed radius = 源半径 × {radius_scale:g}；体素间距 = {1.0 / voxels_per_micron:.4f} µm",
        ),
        (
            "D  从整体形态提取一个三角形边界",
            f"正式 Ultraliser 结果的局部视图；完整表面共 {surface_cells:,} 个三角形",
        ),
    )
    figure, axes = plt.subplots(2, 2, figsize=(18, 14), dpi=160)
    for axis, path, (title, subtitle) in zip(axes.ravel(), panel_paths, titles, strict=True):
        axis.imshow(mpimg.imread(path))
        axis.set_axis_off()
        axis.text(
            0.5,
            1.075,
            title,
            transform=axis.transAxes,
            ha="center",
            va="bottom",
            fontproperties=font,
            fontsize=17,
            color="#171717",
        )
        axis.text(
            0.5,
            1.025,
            subtitle,
            transform=axis.transAxes,
            ha="center",
            va="bottom",
            fontproperties=font,
            fontsize=10.5,
            color="#505050",
        )
    figure.suptitle(
        "实际 ROI003274：中心线 + 半径 → 共享空间占据 → 统一三角形表面",
        fontproperties=font,
        fontsize=23,
        y=0.985,
        color="#111111",
    )
    legend = (
        Patch(facecolor="#78C6E7", edgecolor="#0072B2", label="中心线采样点的 feed 半径"),
        Patch(facecolor="#32B6C5", edgecolor="#174F56", label="解释性共享体素占据"),
        Patch(facecolor="#C83C32", edgecolor="#43201E", label="正式统一三角形表面"),
        Patch(facecolor=JUNCTION_COLOR, edgecolor=JUNCTION_COLOR, label=f"分叉点 local node {junction_node}"),
    )
    figure.legend(
        handles=legend,
        loc="lower center",
        ncol=4,
        frameon=False,
        prop=font,
        bbox_to_anchor=(0.5, 0.045),
    )
    figure.text(
        0.5,
        0.020,
        (
            "阅读顺序 A→C→D。C 由正式输入数据按相同空间分辨率构造，仅用于解释“共享占据/逻辑并集”；"
            "它不是 Ultraliser 导出的内部中间文件，也不替代正式重建。"
        ),
        ha="center",
        va="bottom",
        fontproperties=font,
        fontsize=10.8,
        color="#3F3F3F",
    )
    figure.subplots_adjust(
        left=0.018,
        right=0.982,
        top=0.890,
        bottom=0.085,
        wspace=0.035,
        hspace=0.235,
    )
    figure.savefig(output_path, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _write_readme(path: Path, metadata: dict[str, object]) -> None:
    occupancy = metadata["explanatory_occupancy"]
    surface = metadata["accepted_surface"]
    text = f"""# SWC 到统一 STL：共享体素域可视化

本目录由 `doc_visualize_v2.py` 基于当前建模代码实际使用的数据生成。

## 如何阅读总览图

1. **A：实际输入。** 保存的 ROI 包含中心线节点、连接边和每个节点的源半径；正式送入 Ultraliser 的半径为源半径乘以 `{metadata['radius_scale']}`。
2. **B：反例式对照。** 三条分支各自生成管面时，虽然都位于正确的三维坐标，它们仍是三个相交或重叠的表面对象；分叉处可能含接缝、内表面或重叠关系。
3. **C：共享占据域。** 每条中心线段及其两端半径决定一组占据单元，所有分支对同一个布尔体素数组执行逻辑并集，因此分叉处只保留“内部/外部”的整体分类。这里的单元间距是 `{occupancy['spacing_um']:.6f} µm`（`{metadata['voxels_per_micron']}` voxels/µm）。
4. **D：统一边界。** 正式 Ultraliser 输出从整体形态提取一个三角形外边界；完整表面有 `{surface['triangle_count']}` 个三角形、`{surface['component_count']}` 个连通分量，并且 watertight=`{surface['watertight']}`。

因此，原文所说的“统一内腔边界”不是指把若干分支 STL 文件摆到一起，而是先在同一个空间域中决定整个血管树哪里属于管腔内部，再只提取这一个整体的外边界。

## 文件

- `swc_to_stl_shared_domain.png`：四阶段总览图。
- `panels/*.png`：四个可单独放大的三维面板。
- `source_centerline_and_radius.vtp`：实际 ROI 中心线，包含 `source_radius_um`、`feed_radius_um`、`degree` 和 `local_node_id` 点数据。
- `explanatory_shared_occupancy.vtu`：局部解释性占据单元，可在 ParaView 中旋转和裁切。
- `explanatory_occupancy_boundary.vtp`：上述解释性占据的外层方格边界。
- `accepted_ultraliser_surface_local_view.vtp`：正式表面在图示立方体内的裁切视图；仅用于显示，不用于 QC。
- `visualization_metadata.json`：数据来源、参数、计数和解释边界。

## 重要边界

`explanatory_shared_occupancy.vtu` 是根据实际中心线、实际 feed 半径和正式体素间距编写的教学性胶囊并集。它说明“全部分支进入一个共享占据域”这一数据结构和几何含义，但不是 Ultraliser C++ 程序导出的内部体素文件，不能用它替代正式 Ultraliser 结果或正式 QC。
"""
    path.write_text(text, encoding="utf-8")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _official_stage_paths(output_dir: Path) -> dict[str, Path]:
    stage_dir = output_dir / "official_stages"
    return {
        "dmc": stage_dir / "01_official_dmc.vtp",
        "laplacian": stage_dir / "02_official_laplacian10.vtp",
        "optimized": stage_dir / "03_official_adaptive_optimized.vtp",
        "watertight": stage_dir / "04_official_watertight.vtp",
    }


def _mesh_summary(mesh: pv.PolyData, path: Path) -> dict[str, object]:
    return {
        "path": str(path),
        "point_count": int(mesh.n_points),
        "triangle_count": int(mesh.n_cells),
        "bounds_um": [float(value) for value in mesh.bounds],
    }


def _load_or_capture_official_stages(
    roi,
    config,
    output_dir: Path,
    ultraliser_root: Path,
    accepted_surface_path: Path,
) -> tuple[dict[str, pv.PolyData], dict[str, object]]:
    """Load retained official stages or capture them once with the formal command."""

    stage_paths = _official_stage_paths(output_dir)
    manifest_path = output_dir / "official_stage_manifest.json"
    if all(path.is_file() for path in stage_paths.values()) and manifest_path.is_file():
        meshes = {
            name: pv.read(path).extract_surface(algorithm="dataset_surface").triangulate().clean()
            for name, path in stage_paths.items()
        }
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        return meshes, manifest

    from utils.cfd_lumen.ultraliser_backend import (
        build_ultraliser_command,
        discover_ultraliser_executable,
        export_roi_for_ultraliser,
        invoke_ultraliser,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    capture_root = output_dir / "_official_capture"
    raw_directory = capture_root / "raw"
    raw_stage_paths = {
        "dmc": raw_directory / "meshes" / "anchor003274_docflow-dmc.stl",
        "laplacian": raw_directory / "meshes" / "anchor003274_docflow-laplacian.stl",
        "optimized": raw_directory / "meshes" / "anchor003274_docflow-optimized.stl",
        "watertight": raw_directory / "meshes" / "anchor003274_docflow-watertight.stl",
    }
    raw_capture_complete = all(path.is_file() for path in raw_stage_paths.values())
    if capture_root.exists() and not raw_capture_complete:
        raise RuntimeError(
            f"Incomplete official stage capture exists at {capture_root}; inspect it before retrying"
        )

    exported = export_roi_for_ultraliser(
        roi,
        capture_root / "input",
        radius_scale_for_ultraliser=float(config.ultraliser.radius_scale),
    )
    executable = discover_ultraliser_executable(ultraliser_root)
    command, command_text = build_ultraliser_command(
        executable,
        exported.h5_path,
        raw_directory,
        prefix="anchor003274_docflow",
        voxels_per_micron=float(config.ultraliser.voxels_per_micron),
        threads=int(config.ultraliser.threads),
        packing_algorithm=config.ultraliser.packing_algorithm,
        voxelization_axis=config.ultraliser.voxelization_axis,
        isosurface_technique=config.ultraliser.isosurface_technique,
        solid_voxelization=bool(config.ultraliser.solid_voxelization),
        adaptive_optimization=bool(config.ultraliser.adaptive_optimization),
        optimization_iterations=int(config.ultraliser.optimization_iterations),
        smooth_iterations=int(config.ultraliser.smooth_iterations),
        laplacian_iterations=int(config.ultraliser.laplacian_iterations),
        export_stl=bool(config.ultraliser.export_stl),
    )
    runtime_seconds: float | None = None
    reused_completed_capture = raw_capture_complete
    if not raw_capture_complete:
        invocation = invoke_ultraliser(
            executable,
            exported.h5_path,
            raw_directory,
            prefix="anchor003274_docflow",
            settings=config.ultraliser,
        )
        command_text = invocation.command_text
        runtime_seconds = float(invocation.runtime_seconds)
    if not all(path.is_file() for path in raw_stage_paths.values()):
        raise RuntimeError("Official Ultraliser run did not produce all four required stage meshes")

    stage_paths["dmc"].parent.mkdir(parents=True, exist_ok=True)
    meshes: dict[str, pv.PolyData] = {}
    raw_hashes: dict[str, str] = {}
    for name, raw_path in raw_stage_paths.items():
        mesh = pv.read(raw_path).extract_surface(algorithm="dataset_surface").triangulate().clean()
        if mesh.n_points == 0 or mesh.n_cells == 0:
            raise RuntimeError(f"Official {name} stage is empty: {raw_path}")
        mesh.save(stage_paths[name], binary=True)
        meshes[name] = mesh
        raw_hashes[name] = _sha256(raw_path)

    accepted_stl_path = accepted_surface_path.with_suffix(".stl")
    accepted_hash = _sha256(accepted_stl_path) if accepted_stl_path.is_file() else None
    manifest: dict[str, object] = {
        "stage_source": "official ultraVessMorpho2Mesh stage exports",
        "explanatory_or_official": {
            "dmc": "official",
            "laplacian": "official",
            "optimized": "official",
            "watertight": "official",
        },
        "roi_id": roi.roi_id,
        "command": command_text,
        "command_argv": command,
        "runtime_seconds": runtime_seconds,
        "reused_completed_temporary_capture": bool(reused_completed_capture),
        "settings": {
            "radius_scale": float(config.ultraliser.radius_scale),
            "voxels_per_micron": float(config.ultraliser.voxels_per_micron),
            "packing_algorithm": config.ultraliser.packing_algorithm,
            "solid_voxelization": bool(config.ultraliser.solid_voxelization),
            "voxelization_axis": config.ultraliser.voxelization_axis,
            "isosurface_technique": config.ultraliser.isosurface_technique,
            "laplacian_iterations": int(config.ultraliser.laplacian_iterations),
            "adaptive_optimization": bool(config.ultraliser.adaptive_optimization),
            "optimization_iterations": int(config.ultraliser.optimization_iterations),
            "smooth_iterations_inside_adaptive_optimization": int(
                config.ultraliser.smooth_iterations
            ),
            "threads": int(config.ultraliser.threads),
        },
        "stages": {
            name: {
                **_mesh_summary(meshes[name], stage_paths[name]),
                "raw_stl_sha256": raw_hashes[name],
            }
            for name in stage_paths
        },
        "accepted_reference_stl": str(accepted_stl_path) if accepted_stl_path.is_file() else None,
        "accepted_reference_stl_sha256": accepted_hash,
        "official_watertight_stl_sha256_matches_accepted": bool(
            accepted_hash is not None and raw_hashes["watertight"] == accepted_hash
        ),
    }
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    (output_dir / "official_invocation.txt").write_text(command_text + "\n", encoding="utf-8")

    resolved_capture = capture_root.resolve()
    if resolved_capture.parent != output_dir.resolve() or resolved_capture.name != "_official_capture":
        raise RuntimeError(f"Refusing to remove unexpected capture path: {resolved_capture}")
    shutil.rmtree(resolved_capture)
    return meshes, manifest


def _selected_cells(grid: pv.ImageData, mask: np.ndarray, name: str) -> pv.UnstructuredGrid:
    selected = grid.copy(deep=True)
    selected.cell_data.clear()
    selected.cell_data[name] = np.asarray(mask, dtype=np.uint8).ravel(order="F")
    return selected.threshold(0.5, scalars=name)


def _split_shell_and_interior(
    grid: pv.ImageData,
) -> tuple[pv.UnstructuredGrid, pv.UnstructuredGrid, dict[str, int]]:
    dimensions = tuple(int(value) - 1 for value in grid.dimensions)
    occupancy = np.asarray(grid.cell_data["occupied"], dtype=bool).reshape(
        dimensions, order="F"
    )
    interior = ndimage.binary_erosion(
        occupancy,
        structure=np.ones((3, 3, 3), dtype=bool),
        iterations=1,
        border_value=0,
    )
    shell = occupancy & ~interior
    return (
        _selected_cells(grid, shell, "surface_shell"),
        _selected_cells(grid, interior, "solid_interior"),
        {
            "surface_shell_cell_count": int(np.count_nonzero(shell)),
            "solid_interior_cell_count": int(np.count_nonzero(interior)),
            "combined_occupied_cell_count": int(np.count_nonzero(occupancy)),
        },
    )


def _render_polyline_spheres_panel(
    roi,
    graph: nx.Graph,
    radius_scale: float,
    junction_node: int,
    local_bounds: tuple[float, ...],
    output_path: Path,
    window_size: tuple[int, int],
    half_extent: float,
) -> None:
    plotter = _new_plotter(window_size)
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float) * float(radius_scale)
    low = np.asarray([local_bounds[0], local_bounds[2], local_bounds[4]], dtype=float)
    high = np.asarray([local_bounds[1], local_bounds[3], local_bounds[5]], dtype=float)
    selected = np.all((positions + radii[:, None]) >= low, axis=1) & np.all(
        (positions - radii[:, None]) <= high, axis=1
    )
    sphere_points = pv.PolyData(positions[selected])
    sphere_points.point_data["feed_radius_um"] = radii[selected]
    spheres = sphere_points.glyph(
        scale="feed_radius_um",
        orient=False,
        factor=1.0,
        geom=pv.Sphere(radius=1.0, theta_resolution=24, phi_resolution=24),
    )
    plotter.add_mesh(
        spheres,
        color="#62B5D1",
        opacity=0.52,
        show_edges=True,
        edge_color="#235767",
        edge_opacity=0.22,
        smooth_shading=True,
    )
    centerline = _centerline_mesh(roi, graph, radius_scale)
    cropped_line = centerline.clip_box(pv.Box(bounds=local_bounds), invert=False)
    plotter.add_mesh(
        cropped_line,
        color="#151515",
        line_width=5,
        render_lines_as_tubes=True,
    )
    plotter.add_points(
        positions[junction_node][None, :],
        color=JUNCTION_COLOR,
        point_size=17,
        render_points_as_spheres=True,
    )
    _add_axes(plotter, local_bounds)
    _set_local_camera(plotter, positions[junction_node], half_extent)
    _save_panel(plotter, output_path)


def _render_shell_fill_panel(
    shell_cells: pv.UnstructuredGrid,
    interior_cells: pv.UnstructuredGrid,
    roi,
    graph: nx.Graph,
    junction_node: int,
    local_bounds: tuple[float, ...],
    output_path: Path,
    window_size: tuple[int, int],
    half_extent: float,
) -> None:
    plotter = _new_plotter(window_size)
    center = np.asarray(roi.local_node_positions_um[junction_node], dtype=float)
    shell_surface = shell_cells.extract_surface(algorithm="dataset_surface")
    interior_cutaway = interior_cells.clip(
        normal=(1.0, -0.65, 0.15), origin=center, invert=False
    )
    plotter.add_mesh(
        shell_surface,
        color="#E69F00",
        opacity=0.26,
        show_edges=True,
        edge_color="#7A4D00",
        edge_opacity=0.22,
    )
    plotter.add_mesh(
        interior_cutaway,
        color="#2A9D8F",
        opacity=0.88,
        show_edges=True,
        edge_color="#174F49",
        edge_opacity=0.22,
    )
    centerline = _centerline_mesh(roi, graph, radius_scale=1.0)
    plotter.add_mesh(
        centerline.clip_box(pv.Box(bounds=local_bounds), invert=False),
        color="white",
        line_width=5,
        render_lines_as_tubes=True,
    )
    plotter.add_points(
        center[None, :],
        color=JUNCTION_COLOR,
        point_size=17,
        render_points_as_spheres=True,
    )
    _add_axes(plotter, local_bounds)
    _set_local_camera(plotter, center, half_extent)
    _save_panel(plotter, output_path)


def _render_official_stage_panel(
    mesh: pv.PolyData,
    color: str,
    roi,
    graph: nx.Graph,
    junction_node: int,
    local_bounds: tuple[float, ...],
    output_path: Path,
    local_mesh_path: Path,
    window_size: tuple[int, int],
    half_extent: float,
) -> pv.PolyData:
    local_mesh = (
        mesh.clip_box(pv.Box(bounds=local_bounds), invert=False)
        .extract_surface(algorithm="dataset_surface")
        .triangulate()
        .clean()
    )
    local_mesh.save(local_mesh_path, binary=True)
    plotter = _new_plotter(window_size)
    plotter.add_mesh(
        local_mesh,
        color=color,
        show_edges=True,
        edge_color="#282828",
        edge_opacity=0.28,
        smooth_shading=False,
        specular=0.12,
        specular_power=14,
    )
    centerline = _centerline_mesh(roi, graph, radius_scale=1.0)
    plotter.add_mesh(
        centerline.clip_box(pv.Box(bounds=local_bounds), invert=False),
        color="#FFE35B",
        line_width=4,
        render_lines_as_tubes=True,
    )
    center = np.asarray(roi.local_node_positions_um[junction_node], dtype=float)
    plotter.add_points(
        center[None, :],
        color=JUNCTION_COLOR,
        point_size=16,
        render_points_as_spheres=True,
    )
    _add_axes(plotter, local_bounds)
    _set_local_camera(plotter, center, half_extent)
    _save_panel(plotter, output_path)
    return local_mesh


def _compose_processing_figure(
    panel_paths: list[Path],
    output_path: Path,
    *,
    config,
    stage_manifest: dict[str, object],
    shell_counts: dict[str, int],
) -> None:
    font = _font()
    stages = stage_manifest["stages"]
    settings = config.ultraliser
    titles = (
        (
            "A  带球端的折线表示",
            "黑线是实际中心线；半透明球的半径为实际 feed radius",
        ),
        (
            "B  体素壳与实体填充",
            (
                f"橙色壳 {shell_counts['surface_shell_cell_count']:,} cells；"
                f"蓝绿色剖切内部 {shell_counts['solid_interior_cell_count']:,} cells"
            ),
        ),
        (
            "C  官方 DMC 初始三角面",
            f"从实体占据域提取；完整表面 {stages['dmc']['triangle_count']:,} triangles",
        ),
        (
            f"D  官方 Laplacian 处理（{settings.laplacian_iterations} 次）",
            (
                "移动相邻顶点以缓和局部起伏；拓扑规模保持为 "
                f"{stages['laplacian']['triangle_count']:,} triangles"
            ),
        ),
        (
            "E  官方自适应优化与表面平滑",
            (
                f"optimization={settings.optimization_iterations}，内部 smooth="
                f"{settings.smooth_iterations}；降至 {stages['optimized']['triangle_count']:,} triangles"
            ),
        ),
        (
            "F  官方 watertight 最终表面",
            (
                f"{stages['watertight']['triangle_count']:,} triangles；与当前验收 STL 的 SHA-256 "
                f"{'完全一致' if stage_manifest['official_watertight_stl_sha256_matches_accepted'] else '不一致'}"
            ),
        ),
    )
    figure, axes = plt.subplots(3, 2, figsize=(18, 20), dpi=160)
    for axis, path, (title, subtitle) in zip(axes.ravel(), panel_paths, titles, strict=True):
        axis.imshow(mpimg.imread(path))
        axis.set_axis_off()
        axis.text(
            0.5,
            1.075,
            title,
            transform=axis.transAxes,
            ha="center",
            va="bottom",
            fontproperties=font,
            fontsize=16.5,
            color="#171717",
        )
        axis.text(
            0.5,
            1.025,
            subtitle,
            transform=axis.transAxes,
            ha="center",
            va="bottom",
            fontproperties=font,
            fontsize=10.3,
            color="#505050",
        )
    figure.suptitle(
        "实际 ROI003274：polylines-with-spheres → solid fill → DMC → 网格处理",
        fontproperties=font,
        fontsize=22,
        y=0.989,
        color="#111111",
    )
    legend = (
        Patch(facecolor="#62B5D1", edgecolor="#235767", label="实际中心线采样球"),
        Patch(facecolor="#E69F00", edgecolor="#7A4D00", label="解释性体素壳"),
        Patch(facecolor="#2A9D8F", edgecolor="#174F49", label="解释性实体内部"),
        Patch(facecolor="#6C4CCF", edgecolor="#282828", label="官方中间/最终三角面"),
    )
    figure.legend(
        handles=legend,
        loc="lower center",
        ncol=4,
        frameon=False,
        prop=font,
        bbox_to_anchor=(0.5, 0.039),
    )
    figure.text(
        0.5,
        0.016,
        (
            "实际代码顺序：surface voxelization → xyz solid fill → DMC → Laplacian(10) → "
            "adaptive optimization〔optimization=5，内部 smooth=5〕→ watertight。"
            "A–B 为基于真实输入的教学性展开；C–F 均为同一次官方 Ultraliser 调用导出的真实阶段网格。"
        ),
        ha="center",
        va="bottom",
        fontproperties=font,
        fontsize=10.6,
        color="#3F3F3F",
    )
    figure.subplots_adjust(
        left=0.018,
        right=0.982,
        top=0.925,
        bottom=0.070,
        wspace=0.035,
        hspace=0.235,
    )
    figure.savefig(output_path, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _write_processing_readme(path: Path, metadata: dict[str, object]) -> None:
    stages = metadata["official_stage_manifest"]["stages"]
    settings = metadata["settings"]
    text = f"""# 4.4.2 整体处理流程：三维可视化

本目录继续使用 anchor 3274 的真实保存 ROI 和当前 `swc_stl_model_generate.py` 配置，解释 `polylines-with-spheres → solid fill → DMC → mesh processing`。

## 六个面板分别表示什么

1. **A：带球端折线。** 折线负责连接实际中心线采样点；每个采样点的球由其 feed radius 决定大小。球和相邻线段共同标记血管应占据的空间，而不等于最终 STL 三角面。
2. **B：壳与实体。** 表面体素化先得到封闭壳；`--solid --voxelization-axis xyz` 再把壳内部标为实体。图中橙色是一层教学性壳单元，蓝绿色是剖切后的内部单元。该体素展开使用真实输入和 `{settings['voxels_per_micron']}` voxels/µm，但不是 Ultraliser 导出的内部 volume。
3. **C：DMC。** 官方 `--isosurface-technique dmc` 从实体域的内外分界产生初始网格，共 `{stages['dmc']['triangle_count']}` 个三角形。
4. **D：Laplacian。** 当前程序先执行 `{settings['laplacian_iterations']}` 次 Laplacian 处理。它主要根据邻域移动顶点，缓和局部高频起伏；此处三角形数仍为 `{stages['laplacian']['triangle_count']}`。
5. **E：自适应优化与平滑。** `--adaptive-optimization --optimization-iterations {settings['optimization_iterations']} --smooth-iterations {settings['smooth_iterations_inside_adaptive_optimization']}` 调用的优化过程包含法向处理、表面平滑、细化以及稠密/平坦区域简化。网格减少到 `{stages['optimized']['triangle_count']}` 个三角形，而不是简单地对所有三角形做同一种移动。
6. **F：watertight。** 最终阶段检查并修复为统一水密表面。捕获的 watertight STL 与当前验收 STL 的 SHA-256 完全一致，因此 C–F 确实来自生成当前验收结果的同一参数组合。

## 文件真实性边界

- `official_stages/*.vtp` 和 `official_local_views/*.vtp`：由官方 `ultraVessMorpho2Mesh` 的 DMC、Laplacian、optimized 和 watertight STL 转换得到，属于真实阶段几何。
- `explanatory_surface_shell.vtu`、`explanatory_solid_interior.vtu`：根据同一实际中心线、feed radius 与体素间距构造的教学性体素拆分；用于说明“壳内填充”，不冒充 Ultraliser 内部 volume。
- `official_invocation.txt`：本次阶段捕获使用的完整命令。
- `processing_stage_metadata.json`：阶段计数、参数、SHA-256 核对和文件路径。
"""
    path.write_text(text, encoding="utf-8")


def _generate_processing_stage_visualization(
    *,
    roi,
    graph: nx.Graph,
    config,
    occupancy_grid: pv.ImageData,
    junction_node: int,
    local_bounds: tuple[float, ...],
    half_extent: float,
    window_size: tuple[int, int],
    output_dir: Path,
    ultraliser_root: Path,
    accepted_surface_path: Path,
) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    panels_dir = output_dir / "panels"
    local_views_dir = output_dir / "official_local_views"
    panels_dir.mkdir(parents=True, exist_ok=True)
    local_views_dir.mkdir(parents=True, exist_ok=True)

    stage_meshes, stage_manifest = _load_or_capture_official_stages(
        roi,
        config,
        output_dir,
        ultraliser_root,
        accepted_surface_path,
    )
    shell_cells, interior_cells, shell_counts = _split_shell_and_interior(occupancy_grid)
    shell_cells.save(output_dir / "explanatory_surface_shell.vtu", binary=True)
    interior_cells.save(output_dir / "explanatory_solid_interior.vtu", binary=True)

    panel_paths = [
        panels_dir / "01_polylines_with_spheres.png",
        panels_dir / "02_surface_shell_and_solid_fill.png",
        panels_dir / "03_official_dmc.png",
        panels_dir / "04_official_laplacian10.png",
        panels_dir / "05_official_adaptive_optimized.png",
        panels_dir / "06_official_watertight.png",
    ]
    _render_polyline_spheres_panel(
        roi,
        graph,
        float(config.ultraliser.radius_scale),
        junction_node,
        local_bounds,
        panel_paths[0],
        window_size,
        half_extent,
    )
    _render_shell_fill_panel(
        shell_cells,
        interior_cells,
        roi,
        graph,
        junction_node,
        local_bounds,
        panel_paths[1],
        window_size,
        half_extent,
    )
    stage_colors = {
        "dmc": "#6C4CCF",
        "laplacian": "#3B78B4",
        "optimized": "#E07A2F",
        "watertight": "#C83C32",
    }
    local_mesh_metadata: dict[str, object] = {}
    for panel_index, stage_name in enumerate(
        ("dmc", "laplacian", "optimized", "watertight"), start=2
    ):
        local_path = local_views_dir / f"{panel_index - 1:02d}_{stage_name}_local.vtp"
        local_mesh = _render_official_stage_panel(
            stage_meshes[stage_name],
            stage_colors[stage_name],
            roi,
            graph,
            junction_node,
            local_bounds,
            panel_paths[panel_index],
            local_path,
            window_size,
            half_extent,
        )
        local_mesh_metadata[stage_name] = _mesh_summary(local_mesh, local_path)

    figure_path = output_dir / "swc_to_stl_processing_stages.png"
    _compose_processing_figure(
        panel_paths,
        figure_path,
        config=config,
        stage_manifest=stage_manifest,
        shell_counts=shell_counts,
    )
    metadata: dict[str, object] = {
        "purpose": "Explain section 4.4.2 with actual ROI003274 data and official stage meshes",
        "roi_id": roi.roi_id,
        "junction_node": int(junction_node),
        "junction_position_um": [
            float(value) for value in roi.local_node_positions_um[junction_node]
        ],
        "stage_order_used_by_code": [
            "polylines-with-spheres surface voxelization",
            "xyz solid voxelization",
            "DMC extraction",
            "Laplacian smoothing",
            "adaptive optimization including surface smoothing",
            "watertightness processing",
        ],
        "settings": stage_manifest["settings"],
        "explanatory_shell_fill": {
            "classification": "teaching reconstruction, not an Ultraliser-exported volume",
            **shell_counts,
        },
        "official_stage_manifest": stage_manifest,
        "official_local_views": local_mesh_metadata,
        "overview_figure": str(figure_path),
        "panel_paths": [str(path) for path in panel_paths],
    }
    (output_dir / "processing_stage_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    _write_processing_readme(output_dir / "README.md", metadata)
    return figure_path


def _junction_star_input(
    roi,
    graph: nx.Graph,
    junction_node: int,
    radius_scale: float,
    voxels_per_micron: float,
) -> dict[str, object]:
    """Extract the three-edge junction star used for the focused explanation."""

    neighbors = sorted(int(value) for value in graph.neighbors(junction_node))
    if len(neighbors) < 3:
        raise ValueError(f"Node {junction_node} does not define a branch junction")

    local_edges = np.asarray(roi.local_edges, dtype=np.int64)
    edge_indices = [
        index
        for index, endpoints in enumerate(local_edges)
        if int(junction_node) in (int(endpoints[0]), int(endpoints[1]))
    ]
    edge_points = np.asarray(roi.local_edge_points_um, dtype=float)[edge_indices]
    edge_source_radii = np.asarray(roi.local_edge_radius_um, dtype=float)[edge_indices]
    edge_feed_radii = edge_source_radii * float(radius_scale)

    node_ids = np.asarray([int(junction_node), *neighbors], dtype=np.int64)
    node_positions = np.asarray(roi.local_node_positions_um, dtype=float)[node_ids]
    node_source_radii = np.asarray(roi.local_node_radius_um, dtype=float)[node_ids]
    node_feed_radii = node_source_radii * float(radius_scale)
    node_lookup = {int(node): index for index, node in enumerate(node_ids)}
    lines = np.asarray(
        [
            [2, node_lookup[int(a)], node_lookup[int(b)]]
            for a, b in local_edges[edge_indices]
        ],
        dtype=np.int64,
    ).ravel()
    centerline = pv.PolyData(node_positions, lines=lines)
    centerline.point_data["local_node_id"] = node_ids
    centerline.point_data["source_radius_um"] = node_source_radii
    centerline.point_data["feed_radius_um"] = node_feed_radii

    padding_um = 3.0 / float(voxels_per_micron)
    lower = np.min(edge_points - edge_feed_radii[..., None], axis=(0, 1)) - padding_um
    upper = np.max(edge_points + edge_feed_radii[..., None], axis=(0, 1)) + padding_um
    bounds = (
        float(lower[0]),
        float(upper[0]),
        float(lower[1]),
        float(upper[1]),
        float(lower[2]),
        float(upper[2]),
    )
    return {
        "junction_node": int(junction_node),
        "neighbor_nodes": neighbors,
        "node_ids": node_ids,
        "node_positions_um": node_positions,
        "node_source_radii_um": node_source_radii,
        "node_feed_radii_um": node_feed_radii,
        "edge_indices": edge_indices,
        "edge_node_ids": local_edges[edge_indices],
        "edge_points_um": edge_points,
        "edge_source_radii_um": edge_source_radii,
        "edge_feed_radii_um": edge_feed_radii,
        "centerline": centerline,
        "bounds": bounds,
        "center_um": np.asarray(roi.local_node_positions_um[junction_node], dtype=float),
        "half_extent_um": float(np.max(upper - lower) / 2.0),
        "padding_um": padding_um,
    }


def _axis_fill_reconstruction(grid: pv.ImageData) -> dict[str, object]:
    """Make a transparent teaching reconstruction of the xyz solid-fill rule.

    The shell is derived from the actual-data capsule occupancy.  For each
    principal axis, closed areas are filled independently on every orthogonal
    two-dimensional slice.  The three candidate volumes are intersected and
    finally united with the original shell, matching the logical structure of
    Ultraliser's xyz solid voxelisation without claiming bitwise identity with
    its internal volume.
    """

    shape = tuple(int(value) - 1 for value in grid.dimensions)
    source_solid = np.asarray(grid.cell_data["occupied"], dtype=bool).reshape(
        shape, order="F"
    )
    eroded = ndimage.binary_erosion(
        source_solid,
        structure=np.ones((3, 3, 3), dtype=bool),
        iterations=1,
        border_value=0,
    )
    shell = source_solid & ~eroded

    candidates: list[np.ndarray] = []
    for axis in range(3):
        candidate = shell.copy()
        for index in range(shape[axis]):
            selector: list[int | slice] = [slice(None), slice(None), slice(None)]
            selector[axis] = index
            section = tuple(selector)
            candidate[section] = ndimage.binary_fill_holes(shell[section])
        candidates.append(candidate)

    xyz_agreement = np.logical_and.reduce(candidates)
    reconstructed_solid = shell | xyz_agreement
    reconstructed_interior = xyz_agreement & ~shell
    classification = np.zeros(shape, dtype=np.uint8)
    classification[shell] = 1
    classification[reconstructed_interior] = 2

    overlap = int(np.count_nonzero(reconstructed_solid & source_solid))
    union = int(np.count_nonzero(reconstructed_solid | source_solid))
    metrics: dict[str, object] = {
        "shell_cell_count": int(np.count_nonzero(shell)),
        "source_solid_cell_count": int(np.count_nonzero(source_solid)),
        "reconstructed_interior_cell_count": int(
            np.count_nonzero(reconstructed_interior)
        ),
        "reconstructed_solid_cell_count": int(np.count_nonzero(reconstructed_solid)),
        "x_candidate_cell_count": int(np.count_nonzero(candidates[0])),
        "y_candidate_cell_count": int(np.count_nonzero(candidates[1])),
        "z_candidate_cell_count": int(np.count_nonzero(candidates[2])),
        "missing_from_reconstruction": int(
            np.count_nonzero(source_solid & ~reconstructed_solid)
        ),
        "extra_in_reconstruction": int(
            np.count_nonzero(reconstructed_solid & ~source_solid)
        ),
        "intersection_over_union": float(overlap / union) if union else 1.0,
    }
    return {
        "source_solid": source_solid,
        "shell": shell,
        "candidates": candidates,
        "xyz_agreement": xyz_agreement,
        "interior": reconstructed_interior,
        "solid": reconstructed_solid,
        "classification": classification,
        "metrics": metrics,
    }


def _add_star_centerline(
    plotter: pv.Plotter,
    star: dict[str, object],
    *,
    colored: bool,
    line_width: float = 7.0,
) -> None:
    for index, endpoints in enumerate(np.asarray(star["edge_points_um"], dtype=float)):
        plotter.add_mesh(
            _polyline(endpoints),
            color=(BRANCH_COLORS[index % len(BRANCH_COLORS)] if colored else "#171717"),
            line_width=line_width,
            render_lines_as_tubes=True,
        )


def _add_star_envelopes(plotter: pv.Plotter, star: dict[str, object]) -> None:
    for index, (endpoints, radii) in enumerate(
        zip(
            np.asarray(star["edge_points_um"], dtype=float),
            np.asarray(star["edge_feed_radii_um"], dtype=float),
            strict=True,
        )
    ):
        tube = _polyline(endpoints, radii).tube(
            radius=None,
            scalars="feed_radius_um",
            absolute=True,
            capping=True,
            n_sides=28,
        )
        plotter.add_mesh(
            tube,
            color=BRANCH_COLORS[index % len(BRANCH_COLORS)],
            opacity=0.52,
            show_edges=True,
            edge_color="#3C3C3C",
            edge_opacity=0.18,
            smooth_shading=True,
        )

    sphere_points = pv.PolyData(np.asarray(star["node_positions_um"], dtype=float))
    sphere_points.point_data["feed_radius_um"] = np.asarray(
        star["node_feed_radii_um"], dtype=float
    )
    spheres = sphere_points.glyph(
        scale="feed_radius_um",
        orient=False,
        factor=1.0,
        geom=pv.Sphere(radius=1.0, theta_resolution=28, phi_resolution=28),
    )
    plotter.add_mesh(
        spheres,
        color="#79C6DB",
        opacity=0.38,
        show_edges=True,
        edge_color="#245B68",
        edge_opacity=0.14,
        smooth_shading=True,
    )


def _grid_plane_wireframe(
    grid: pv.ImageData,
    normal_axis: int,
    coordinate: float,
    stride: int = 2,
) -> pv.PolyData:
    origin = np.asarray(grid.origin, dtype=float)
    spacing = np.asarray(grid.spacing, dtype=float)
    dimensions = np.asarray(grid.dimensions, dtype=np.int64)
    boundaries = [
        origin[axis] + np.arange(dimensions[axis], dtype=float) * spacing[axis]
        for axis in range(3)
    ]
    plane_index = int(
        np.clip(
            np.rint((float(coordinate) - origin[normal_axis]) / spacing[normal_axis]),
            0,
            dimensions[normal_axis] - 1,
        )
    )
    plane_coordinate = float(boundaries[normal_axis][plane_index])
    in_plane = [axis for axis in range(3) if axis != normal_axis]
    first, second = in_plane
    segments: list[tuple[np.ndarray, np.ndarray]] = []
    for value in boundaries[second][::stride]:
        start = origin.copy()
        end = origin.copy()
        start[normal_axis] = plane_coordinate
        end[normal_axis] = plane_coordinate
        start[second] = value
        end[second] = value
        start[first] = boundaries[first][0]
        end[first] = boundaries[first][-1]
        segments.append((start, end))
    for value in boundaries[first][::stride]:
        start = origin.copy()
        end = origin.copy()
        start[normal_axis] = plane_coordinate
        end[normal_axis] = plane_coordinate
        start[first] = value
        end[first] = value
        start[second] = boundaries[second][0]
        end[second] = boundaries[second][-1]
        segments.append((start, end))
    points = np.asarray([point for segment in segments for point in segment], dtype=float)
    lines = np.asarray(
        [[2, 2 * index, 2 * index + 1] for index in range(len(segments))],
        dtype=np.int64,
    ).ravel()
    return pv.PolyData(points, lines=lines)


def _cutaway_mask(grid: pv.ImageData, mask: np.ndarray, center: np.ndarray) -> np.ndarray:
    shape = tuple(int(value) - 1 for value in grid.dimensions)
    axes = [
        float(grid.origin[axis])
        + (np.arange(shape[axis], dtype=float) + 0.5) * float(grid.spacing[axis])
        for axis in range(3)
    ]
    xx, yy, zz = np.meshgrid(*axes, indexing="ij")
    points = np.stack((xx, yy, zz), axis=-1)
    view_direction = np.asarray([1.55, -2.0, 1.25], dtype=float)
    view_direction /= np.linalg.norm(view_direction)
    keep_back_half = np.einsum("...i,i->...", points - center, view_direction) <= 0.0
    return np.asarray(mask, dtype=bool) & keep_back_half


def _best_axis_slice(mask: np.ndarray, axis: int) -> int:
    reduce_axes = tuple(value for value in range(3) if value != axis)
    return int(np.argmax(np.count_nonzero(mask, axis=reduce_axes)))


def _render_voxelization_panels(
    *,
    star: dict[str, object],
    grid: pv.ImageData,
    occupied_cells: pv.UnstructuredGrid,
    occupancy_boundary: pv.PolyData,
    masks: dict[str, object],
    panel_paths: list[Path],
    window_size: tuple[int, int],
) -> None:
    bounds = tuple(float(value) for value in star["bounds"])
    center = np.asarray(star["center_um"], dtype=float)
    half_extent = float(star["half_extent_um"])
    node_positions = np.asarray(star["node_positions_um"], dtype=float)
    node_ids = np.asarray(star["node_ids"], dtype=np.int64)

    plotter = _new_plotter(window_size)
    _add_star_centerline(plotter, star, colored=True)
    plotter.add_points(
        node_positions,
        color="#202020",
        point_size=15,
        render_points_as_spheres=True,
    )
    plotter.add_point_labels(
        node_positions,
        [str(int(value)) for value in node_ids],
        font_size=19,
        text_color="#202020",
        point_size=0,
        shape_opacity=0.55,
    )
    _add_axes(plotter, bounds)
    _set_local_camera(plotter, center, half_extent)
    _save_panel(plotter, panel_paths[0])

    plotter = _new_plotter(window_size)
    _add_star_envelopes(plotter, star)
    _add_star_centerline(plotter, star, colored=False, line_width=5)
    plotter.add_points(
        center[None, :],
        color=JUNCTION_COLOR,
        point_size=18,
        render_points_as_spheres=True,
    )
    _add_axes(plotter, bounds)
    _set_local_camera(plotter, center, half_extent)
    _save_panel(plotter, panel_paths[1])

    plotter = _new_plotter(window_size)
    plotter.add_mesh(
        occupied_cells,
        color="#65BDD1",
        opacity=0.72,
        show_edges=True,
        edge_color="#2F5F69",
        edge_opacity=0.32,
    )
    for axis, color in enumerate(("#D55E00", "#009E73", "#0072B2")):
        plotter.add_mesh(
            _grid_plane_wireframe(grid, axis, center[axis]),
            color=color,
            opacity=0.31,
            line_width=1.0,
        )
    plotter.add_mesh(grid.outline(), color="#252525", line_width=2.0)
    _add_star_centerline(plotter, star, colored=False, line_width=4)
    _add_axes(plotter, bounds)
    _set_local_camera(plotter, center, half_extent)
    _save_panel(plotter, panel_paths[2])

    shell = np.asarray(masks["shell"], dtype=bool)
    cut_shell = _selected_cells(
        grid,
        _cutaway_mask(grid, shell, center),
        "cutaway_surface_shell",
    )
    plotter = _new_plotter(window_size)
    plotter.add_mesh(
        occupancy_boundary,
        color="#E69F00",
        opacity=0.08,
        style="wireframe",
        line_width=1.0,
    )
    plotter.add_mesh(
        cut_shell,
        color="#E69F00",
        opacity=0.94,
        show_edges=True,
        edge_color="#6D4700",
        edge_opacity=0.42,
    )
    _add_star_centerline(plotter, star, colored=False, line_width=4)
    _add_axes(plotter, bounds)
    _set_local_camera(plotter, center, half_extent)
    _save_panel(plotter, panel_paths[3])

    candidates = [np.asarray(value, dtype=bool) for value in masks["candidates"]]
    axis_colors = ("#D55E00", "#009E73", "#0072B2")
    plotter = _new_plotter(window_size)
    plotter.add_mesh(
        _selected_cells(grid, shell, "axis_fill_shell"),
        color="#E69F00",
        opacity=0.13,
        show_edges=True,
        edge_color="#6D4700",
        edge_opacity=0.12,
    )
    shape = shell.shape
    for axis, (candidate, color) in enumerate(zip(candidates, axis_colors, strict=True)):
        candidate_interior = candidate & ~shell
        index = _best_axis_slice(candidate_interior, axis)
        slab = np.zeros(shape, dtype=bool)
        selector: list[int | slice] = [slice(None), slice(None), slice(None)]
        selector[axis] = index
        section = tuple(selector)
        slab[section] = candidate_interior[section]
        plotter.add_mesh(
            _selected_cells(grid, slab, f"axis_{axis}_representative_slice"),
            color=color,
            opacity=0.78,
            show_edges=True,
            edge_color="#303030",
            edge_opacity=0.24,
        )
    low = np.asarray([bounds[0], bounds[2], bounds[4]], dtype=float)
    high = np.asarray([bounds[1], bounds[3], bounds[5]], dtype=float)
    arrow_length = float(np.min(high - low) * 0.22)
    arrow_start = low + (high - low) * np.asarray([0.13, 0.13, 0.16])
    for axis, color in enumerate(axis_colors):
        vector = np.zeros(3, dtype=float)
        vector[axis] = arrow_length
        plotter.add_arrows(arrow_start[None, :], vector[None, :], color=color, mag=1.0)
    _add_axes(plotter, bounds)
    _set_local_camera(plotter, center, half_extent)
    _save_panel(plotter, panel_paths[4])

    plotter = _new_plotter(window_size)
    plotter.add_mesh(
        _selected_cells(grid, shell, "classified_shell"),
        color="#E69F00",
        opacity=0.42,
        show_edges=True,
        edge_color="#6D4700",
        edge_opacity=0.35,
    )
    plotter.add_mesh(
        _selected_cells(grid, np.asarray(masks["interior"], dtype=bool), "classified_interior"),
        color="#2A9D8F",
        opacity=0.94,
        show_edges=True,
        edge_color="#174F49",
        edge_opacity=0.28,
    )
    plotter.add_mesh(grid.outline(), color="#252525", line_width=2.0)
    _add_axes(plotter, bounds)
    _set_local_camera(plotter, center, half_extent)
    _save_panel(plotter, panel_paths[5])


def _compose_voxelization_figure(
    panel_paths: list[Path],
    output_path: Path,
    *,
    details: dict[str, object],
    metrics: dict[str, object],
    star: dict[str, object],
) -> None:
    font = _font()
    spacing = float(details["spacing_um"])
    dimensions = details["dimensions_cells"]
    titles = (
        (
            "A  真实中心线采样点与折线",
            (
                f"分叉节点 {star['junction_node']} 连接真实邻点 "
                f"{', '.join(str(value) for value in star['neighbor_nodes'])}；折线只说明延伸方向"
            ),
        ),
        (
            "B  给采样点加入实际局部半径",
            "半透明球表示局部粗细；相邻球与线段共同扫出连续作用范围",
        ),
        (
            "C  把作用范围放进规则小立方网格",
            (
                f"cell 边长 {spacing:.4f} µm；网格 {dimensions[0]}×{dimensions[1]}×"
                f"{dimensions[2]}，蓝色为被作用范围覆盖的 cell"
            ),
        ),
        (
            "D  表面体素化先留下封闭外壳",
            (
                f"剖去靠近观察者的一半后可见空腔；橙色壳含 "
                f"{metrics['shell_cell_count']:,} 个 cell"
            ),
        ),
        (
            "E  分别进行 x、y、z 三组切片填充",
            "红/绿/蓝为三种方向的代表性截面；只有三者共同认定的内部才被接受",
        ),
        (
            "F  外壳与内部得到最终体素分类",
            (
                f"橙色壳 + 蓝绿色内部 = {metrics['reconstructed_solid_cell_count']:,} cells；"
                "此刻仍是小立方单元，三角形数量为 0"
            ),
        ),
    )
    figure, axes = plt.subplots(3, 2, figsize=(18, 20), dpi=160)
    for axis, path, (title, subtitle) in zip(axes.ravel(), panel_paths, titles, strict=True):
        axis.imshow(mpimg.imread(path))
        axis.set_axis_off()
        axis.text(
            0.5,
            1.075,
            title,
            transform=axis.transAxes,
            ha="center",
            va="bottom",
            fontproperties=font,
            fontsize=16.2,
            color="#171717",
        )
        axis.text(
            0.5,
            1.025,
            subtitle,
            transform=axis.transAxes,
            ha="center",
            va="bottom",
            fontproperties=font,
            fontsize=10.1,
            color="#505050",
            wrap=True,
        )
    figure.suptitle(
        "实际 ROI003274 分叉：中心线与半径如何变成有内外标记的体素实体",
        fontproperties=font,
        fontsize=21.5,
        y=0.989,
        color="#111111",
    )
    legend = (
        Patch(facecolor="#65BDD1", edgecolor="#2F5F69", label="作用范围覆盖的体素"),
        Patch(facecolor="#E69F00", edgecolor="#6D4700", label="封闭体素壳"),
        Patch(facecolor="#2A9D8F", edgecolor="#174F49", label="xyz 共同认定的内部体素"),
    )
    figure.legend(
        handles=legend,
        loc="lower center",
        ncol=3,
        frameon=False,
        prop=font,
        bbox_to_anchor=(0.5, 0.043),
    )
    figure.text(
        0.5,
        0.016,
        (
            "逻辑关系：实体体素 = 原始外壳 ∪（X 方向候选 ∩ Y 方向候选 ∩ Z 方向候选）。"
            f"本局部重建与原始占据范围的 IoU = {metrics['intersection_over_union']:.6f}。"
            "该展开使用真实输入与真实体素间距，是教学性内部重建，不冒充 Ultraliser 导出的 volume。"
        ),
        ha="center",
        va="bottom",
        fontproperties=font,
        fontsize=10.4,
        color="#3F3F3F",
    )
    figure.subplots_adjust(
        left=0.018,
        right=0.982,
        top=0.925,
        bottom=0.073,
        wspace=0.035,
        hspace=0.235,
    )
    figure.savefig(output_path, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _slice_view(
    shell: np.ndarray,
    fill: np.ndarray,
    axis: int,
    index: int,
    grid: pv.ImageData,
) -> tuple[np.ndarray, tuple[float, float, float, float], str, str, float, np.ndarray, np.ndarray]:
    category = np.zeros(shell.shape, dtype=np.uint8)
    category[shell] = 1
    category[fill & ~shell] = 2
    origin = np.asarray(grid.origin, dtype=float)
    spacing = np.asarray(grid.spacing, dtype=float)
    shape = np.asarray(shell.shape, dtype=np.int64)
    boundaries = [
        origin[value] + np.arange(shape[value] + 1, dtype=float) * spacing[value]
        for value in range(3)
    ]
    coordinate = float(origin[axis] + (index + 0.5) * spacing[axis])
    if axis == 0:
        image = category[index, :, :].T
        horizontal, vertical = 1, 2
    elif axis == 1:
        image = category[:, index, :].T
        horizontal, vertical = 0, 2
    else:
        image = category[:, :, index].T
        horizontal, vertical = 0, 1
    extent = (
        float(boundaries[horizontal][0]),
        float(boundaries[horizontal][-1]),
        float(boundaries[vertical][0]),
        float(boundaries[vertical][-1]),
    )
    axis_names = ("x", "y", "z")
    return (
        image,
        extent,
        axis_names[horizontal],
        axis_names[vertical],
        coordinate,
        boundaries[horizontal],
        boundaries[vertical],
    )


def _compose_xyz_slice_figure(
    output_path: Path,
    *,
    grid: pv.ImageData,
    masks: dict[str, object],
) -> None:
    font = _font()
    shell = np.asarray(masks["shell"], dtype=bool)
    candidates = [np.asarray(value, dtype=bool) for value in masks["candidates"]]
    combined = np.asarray(masks["solid"], dtype=bool)
    views: list[tuple[np.ndarray, tuple[float, float, float, float], str, str, float, np.ndarray, np.ndarray]] = []
    fixed_axes = ("x", "y", "z")
    for axis in range(3):
        index = _best_axis_slice(candidates[axis] & ~shell, axis)
        views.append(_slice_view(shell, candidates[axis], axis, index, grid))
    combined_index = _best_axis_slice(np.asarray(masks["interior"], dtype=bool), 2)
    views.append(_slice_view(shell, combined, 2, combined_index, grid))

    colors = ("#D55E00", "#009E73", "#0072B2", "#2A9D8F")
    titles = (
        "A  X 方向：逐个 yz 截面填充",
        "B  Y 方向：逐个 xz 截面填充",
        "C  Z 方向：逐个 xy 截面填充",
        "D  X ∩ Y ∩ Z，并与外壳合并",
    )
    figure, axes = plt.subplots(2, 2, figsize=(16, 13), dpi=180)
    for panel_axis, view, color, title, fixed_axis in zip(
        axes.ravel(), views, colors, titles, (*fixed_axes, "z"), strict=True
    ):
        image, extent, x_label, y_label, coordinate, x_boundaries, y_boundaries = view
        panel_axis.imshow(
            image,
            origin="lower",
            interpolation="nearest",
            extent=extent,
            cmap=ListedColormap(("#F7F7F4", "#E69F00", color)),
            vmin=0,
            vmax=2,
            aspect="equal",
        )
        panel_axis.set_xticks(x_boundaries, minor=True)
        panel_axis.set_yticks(y_boundaries, minor=True)
        panel_axis.grid(which="minor", color="white", linewidth=0.28, alpha=0.58)
        panel_axis.tick_params(which="minor", length=0)
        panel_axis.set_xlabel(f"{x_label} (µm)", fontproperties=font, fontsize=11)
        panel_axis.set_ylabel(f"{y_label} (µm)", fontproperties=font, fontsize=11)
        panel_axis.set_title(title, fontproperties=font, fontsize=15, pad=22)
        panel_axis.text(
            0.5,
            1.015,
            f"固定 {fixed_axis} = {coordinate:.3f} µm；每个小方格是一枚体素截面",
            transform=panel_axis.transAxes,
            ha="center",
            va="bottom",
            fontproperties=font,
            fontsize=9.8,
            color="#505050",
        )
    figure.suptitle(
        "同一真实分叉的正交截面：为什么需要 x、y、z 三次判断",
        fontproperties=font,
        fontsize=21,
        y=0.985,
    )
    figure.legend(
        handles=(
            Patch(facecolor="#F7F7F4", edgecolor="#B5B5B5", label="外部体素（0）"),
            Patch(facecolor="#E69F00", edgecolor="#6D4700", label="外壳体素（1）"),
            Patch(facecolor="#D55E00", edgecolor="#743300", label="X 候选内部"),
            Patch(facecolor="#009E73", edgecolor="#005E45", label="Y 候选内部"),
            Patch(facecolor="#0072B2", edgecolor="#003F63", label="Z 候选内部"),
            Patch(facecolor="#2A9D8F", edgecolor="#174F49", label="最终内部（1）"),
        ),
        loc="lower center",
        ncol=6,
        frameon=False,
        prop=font,
        bbox_to_anchor=(0.5, 0.043),
    )
    figure.text(
        0.5,
        0.014,
        (
            "前三幅图各取一张代表性二维截面；实际运算会把同方向的全部截面重新叠成三维候选体。"
            "最后只保留三组候选的交集，再把原始橙色外壳并回去。"
        ),
        ha="center",
        va="bottom",
        fontproperties=font,
        fontsize=10.5,
        color="#3F3F3F",
    )
    figure.subplots_adjust(
        left=0.065,
        right=0.975,
        top=0.905,
        bottom=0.095,
        wspace=0.16,
        hspace=0.23,
    )
    figure.savefig(output_path, bbox_inches="tight", facecolor="white")
    plt.close(figure)


def _write_voxelization_readme(path: Path, metadata: dict[str, object]) -> None:
    inputs = metadata["actual_input"]
    counts = metadata["voxel_counts"]
    grid = metadata["grid"]
    text = f"""# 4.4.5 体素化、表面提取与网格处理：体素阶段详解

本目录使用当前建模输入中的真实 ROI `{metadata['roi_id']}`。为了把每一枚小立方体看清楚，可视化只截取局部分叉节点 {inputs['junction_node']} 及其三个直接邻点 {', '.join(str(value) for value in inputs['neighbor_nodes'])}；三条线段的坐标与半径均直接来自保存的 ROI，未使用人为 Y 形示意数据。

## 六幅三维图如何阅读

1. `A` 只画中心线采样点与相邻点之间的折线。折线回答血管向哪里延伸，但没有描述粗细。
2. `B` 在真实采样点处加入由局部半径决定的球形范围，并用相邻线段连接。三条分支在节点 {inputs['junction_node']} 处共同占据同一片空间，因此分叉不是三个互不相干的表面。
3. `C` 把上述范围放入规则网格。当前间距为 {grid['spacing_um']:.6f} µm，即每微米 {grid['voxels_per_micron']} 枚体素。网格共 {grid['dimensions_cells'][0]}×{grid['dimensions_cells'][1]}×{grid['dimensions_cells'][2]} 个 cell；蓝色 cell 的中心落在血管作用范围内。
4. `D` 只保留占据域最外侧的一层 cell，得到教学性封闭外壳。剖去靠近观察者的一半后，能够看到里面仍是空的。
5. `E` 分别按 x、y、z 三个方向，把每一组正交二维切片中的封闭区填满；三组切片重新叠成三个三维候选体。当前 Ultraliser 源码随后对三个候选体执行逻辑交集，再与原始外壳执行逻辑并集。
6. `F` 中橙色为外壳，蓝绿色为三组判断共同接受的内部。两类 cell 合计 {counts['reconstructed_solid_cell_count']} 个，但它们仍然只是具有 0/1 内外标记的小立方单元；DMC 尚未开始，因此此阶段没有三角形。

## 正交截面图如何阅读

`xyz_fill_orthogonal_slices.png` 把三维过程切开。白色小格是外部，橙色小格是封闭壳，彩色小格是某一方向在二维截面内填入的区域。前三幅图分别代表 x、y、z 三组切片；第四幅图显示三组结果的交集与外壳合并后的分类。

逻辑关系为：`solid = shell ∪ (candidate_x ∩ candidate_y ∩ candidate_z)`。

本局部教学重建含 {counts['shell_cell_count']} 个外壳 cell 和 {counts['reconstructed_interior_cell_count']} 个内部 cell，与同一真实线段和半径直接计算的占据范围相比，IoU 为 {counts['intersection_over_union']:.6f}，漏填 {counts['missing_from_reconstruction']} 个，误填 {counts['extra_in_reconstruction']} 个。

## 真实性边界

- 中心线、节点坐标、源半径、0.91 半径缩放和 6 voxels/µm 均来自当前真实输入与配置。
- `surface → x/y/z slice fill → AND → OR shell` 的组合顺序依据当前 `Ultraliser/ultraliser/data/volumes/Volume.cpp`。
- 本目录中的壳与填充体由 Python 以同一实际数据重建，作用是把 Ultraliser 内部不便直接观察的逻辑展开；它不是从 Ultraliser 导出的内部 volume，也不用于替代正式建模。

## 输出文件

- `voxelization_mechanism_3d.png`：六阶段三维总图。
- `xyz_fill_orthogonal_slices.png`：x、y、z 三组填充及逻辑合并的正交截面图。
- `panels/*.png`：总图中的六张独立三维面板。
- `actual_junction_star_centerline.vtp`：四个真实节点及三条真实相邻线段。
- `teaching_voxel_domain.vti`：包含原始占据、壳、x/y/z 候选、交集和最终分类的规则网格。
- `teaching_*.vtu`：各类体素单元的可独立查看版本。
- `teaching_occupancy_boundary.vtp`：教学占据域的外边界，仅用于辅助观察。
- `voxelization_metadata.json`：节点坐标、半径、网格规模、分类计数和一致性核对。
"""
    path.write_text(text, encoding="utf-8")


def _generate_voxelization_explanation(
    *,
    roi,
    graph: nx.Graph,
    config,
    junction_node: int,
    output_dir: Path,
    window_size: tuple[int, int],
) -> tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    panels_dir = output_dir / "panels"
    panels_dir.mkdir(parents=True, exist_ok=True)

    radius_scale = float(config.ultraliser.radius_scale)
    voxels_per_micron = float(config.ultraliser.voxels_per_micron)
    star = _junction_star_input(
        roi,
        graph,
        junction_node,
        radius_scale,
        voxels_per_micron,
    )
    grid, occupied_cells, occupancy_boundary, details = _capsule_union_occupancy(
        np.asarray(star["edge_points_um"], dtype=float),
        np.asarray(star["edge_feed_radii_um"], dtype=float),
        tuple(float(value) for value in star["bounds"]),
        voxels_per_micron,
    )
    masks = _axis_fill_reconstruction(grid)
    metrics = masks["metrics"]

    panel_paths = [
        panels_dir / "01_actual_points_and_polylines.png",
        panels_dir / "02_actual_radius_spheres.png",
        panels_dir / "03_regular_voxel_grid.png",
        panels_dir / "04_hollow_surface_shell.png",
        panels_dir / "05_xyz_slice_fill.png",
        panels_dir / "06_classified_voxel_solid.png",
    ]
    _render_voxelization_panels(
        star=star,
        grid=grid,
        occupied_cells=occupied_cells,
        occupancy_boundary=occupancy_boundary,
        masks=masks,
        panel_paths=panel_paths,
        window_size=window_size,
    )

    figure_path = output_dir / "voxelization_mechanism_3d.png"
    _compose_voxelization_figure(
        panel_paths,
        figure_path,
        details=details,
        metrics=metrics,
        star=star,
    )
    slice_figure_path = output_dir / "xyz_fill_orthogonal_slices.png"
    _compose_xyz_slice_figure(slice_figure_path, grid=grid, masks=masks)

    grid.cell_data["source_solid"] = np.asarray(
        masks["source_solid"], dtype=np.uint8
    ).ravel(order="F")
    grid.cell_data["surface_shell"] = np.asarray(masks["shell"], dtype=np.uint8).ravel(
        order="F"
    )
    for axis_name, candidate in zip(
        ("x_fill_candidate", "y_fill_candidate", "z_fill_candidate"),
        masks["candidates"],
        strict=True,
    ):
        grid.cell_data[axis_name] = np.asarray(candidate, dtype=np.uint8).ravel(order="F")
    grid.cell_data["xyz_agreement"] = np.asarray(
        masks["xyz_agreement"], dtype=np.uint8
    ).ravel(order="F")
    grid.cell_data["reconstructed_solid"] = np.asarray(
        masks["solid"], dtype=np.uint8
    ).ravel(order="F")
    grid.cell_data["classification_0outside_1shell_2interior"] = np.asarray(
        masks["classification"], dtype=np.uint8
    ).ravel(order="F")

    star["centerline"].save(output_dir / "actual_junction_star_centerline.vtp", binary=True)
    grid.save(output_dir / "teaching_voxel_domain.vti", binary=True)
    occupied_cells.save(output_dir / "teaching_source_occupied_cells.vtu", binary=True)
    occupancy_boundary.save(output_dir / "teaching_occupancy_boundary.vtp", binary=True)
    _selected_cells(grid, np.asarray(masks["shell"], dtype=bool), "surface_shell").save(
        output_dir / "teaching_surface_shell.vtu", binary=True
    )
    for axis_name, candidate in zip(
        ("x", "y", "z"), masks["candidates"], strict=True
    ):
        _selected_cells(
            grid,
            np.asarray(candidate, dtype=bool),
            f"{axis_name}_fill_candidate",
        ).save(output_dir / f"teaching_{axis_name}_fill_candidate.vtu", binary=True)
    _selected_cells(
        grid,
        np.asarray(masks["interior"], dtype=bool),
        "xyz_agreement_interior",
    ).save(output_dir / "teaching_xyz_agreement_interior.vtu", binary=True)
    _selected_cells(
        grid,
        np.asarray(masks["solid"], dtype=bool),
        "reconstructed_solid",
    ).save(output_dir / "teaching_reconstructed_solid.vtu", binary=True)

    node_records = []
    for node, position, source_radius, feed_radius in zip(
        np.asarray(star["node_ids"], dtype=np.int64),
        np.asarray(star["node_positions_um"], dtype=float),
        np.asarray(star["node_source_radii_um"], dtype=float),
        np.asarray(star["node_feed_radii_um"], dtype=float),
        strict=True,
    ):
        node_records.append(
            {
                "local_node_id": int(node),
                "position_um": [float(value) for value in position],
                "source_radius_um": float(source_radius),
                "feed_radius_um": float(feed_radius),
            }
        )
    edge_records = []
    for edge_index, node_ids, points, source_radii, feed_radii in zip(
        star["edge_indices"],
        np.asarray(star["edge_node_ids"], dtype=np.int64),
        np.asarray(star["edge_points_um"], dtype=float),
        np.asarray(star["edge_source_radii_um"], dtype=float),
        np.asarray(star["edge_feed_radii_um"], dtype=float),
        strict=True,
    ):
        edge_records.append(
            {
                "roi_edge_index": int(edge_index),
                "local_node_ids": [int(value) for value in node_ids],
                "points_um": [[float(value) for value in point] for point in points],
                "source_radii_um": [float(value) for value in source_radii],
                "feed_radii_um": [float(value) for value in feed_radii],
            }
        )

    metadata: dict[str, object] = {
        "purpose": "Explain section 4.4.5 from actual centreline/radius data to classified voxels",
        "production_equivalence_boundary": (
            "Actual saved ROI inputs and production spacing are used. The capsule occupancy, "
            "surface shell, and xyz slice fills are a teaching reconstruction rather than an "
            "Ultraliser-exported internal volume."
        ),
        "roi_id": roi.roi_id,
        "anchor_id": int(roi.anchor_id),
        "actual_input": {
            "junction_node": int(star["junction_node"]),
            "neighbor_nodes": [int(value) for value in star["neighbor_nodes"]],
            "nodes": node_records,
            "edges": edge_records,
            "radius_scale": radius_scale,
        },
        "grid": {
            **details,
            "voxels_per_micron": voxels_per_micron,
            "classification": {
                "0": "outside",
                "1": "surface shell",
                "2": "xyz-agreement interior",
            },
        },
        "xyz_fill_logic": {
            "per_axis_operation": (
                "fill closed areas independently on every 2-D slice orthogonal to that axis"
            ),
            "equation": "solid = shell OR (candidate_x AND candidate_y AND candidate_z)",
            "official_source_basis": (
                "Ultraliser/ultraliser/data/volumes/Volume.cpp: solidVoxelization and "
                "_floodFillAlongAxis"
            ),
        },
        "voxel_counts": metrics,
        "triangles_at_this_stage": 0,
        "outputs": {
            "overview_3d": str(figure_path),
            "orthogonal_slices": str(slice_figure_path),
            "panels": [str(value) for value in panel_paths],
        },
    }
    (output_dir / "voxelization_metadata.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    _write_voxelization_readme(output_dir / "README.md", metadata)
    return figure_path, slice_figure_path


def _read_radius_validation_run(run_root: Path) -> dict[str, object]:
    root = run_root.resolve()
    summary_path = root / "qc" / "run_summary.json"
    radius_path = root / "qc" / "radius_fidelity.json"
    surface_path = root / "geometry" / "lumen_surface_um.vtp"
    missing = [
        path for path in (summary_path, radius_path, surface_path) if not path.is_file()
    ]
    if missing:
        raise FileNotFoundError(
            "Incomplete model run for radius validation: "
            + ", ".join(str(path) for path in missing)
        )
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    radius = json.loads(radius_path.read_text(encoding="utf-8"))
    samples = radius.get("samples", [])
    if not isinstance(samples, list) or not samples:
        raise ValueError(f"Radius validation run has no successful samples: {root}")
    return {
        "root": root,
        "summary_path": summary_path,
        "radius_path": radius_path,
        "surface_path": surface_path,
        "summary": summary,
        "radius": radius,
        "samples": samples,
    }


def _radius_sample_map(run: dict[str, object]) -> dict[tuple[int, int], dict[str, object]]:
    samples = run["samples"]
    assert isinstance(samples, list)
    mapped: dict[tuple[int, int], dict[str, object]] = {}
    for sample in samples:
        if not isinstance(sample, dict):
            raise TypeError("Radius-fidelity samples must be JSON objects")
        key = (int(sample["branch_id"]), int(sample["sample_index"]))
        if key in mapped:
            raise ValueError(f"Duplicate radius-fidelity sample key: {key}")
        mapped[key] = sample
    return mapped


def _metric_percent(report: dict[str, object], name: str) -> float:
    value = report.get(name)
    if value is None:
        raise ValueError(f"Missing radius-fidelity metric: {name}")
    return 100.0 * float(value)


def _compose_junction_before_after(
    *,
    baseline: dict[str, object],
    compensated: dict[str, object],
    roi,
    graph: nx.Graph,
    junction_node: int,
    local_bounds: tuple[float, ...],
    half_extent: float,
    output_dir: Path,
) -> Path:
    panels_dir = output_dir / "panels"
    panels_dir.mkdir(parents=True, exist_ok=True)
    baseline_panel = panels_dir / "junction_radius100.png"
    compensated_panel = panels_dir / "junction_radius091.png"
    window_size = (1100, 820)

    for run, panel in ((baseline, baseline_panel), (compensated, compensated_panel)):
        surface = (
            pv.read(Path(run["surface_path"]))
            .extract_surface(algorithm="dataset_surface")
            .triangulate()
            .clean()
        )
        local_surface = surface.clip_box(
            pv.Box(bounds=local_bounds), invert=False
        ).extract_surface(algorithm="dataset_surface")
        if local_surface.n_cells == 0:
            raise RuntimeError(
                f"Radius-validation surface does not intersect junction crop: {run['surface_path']}"
            )
        _render_surface_panel(
            local_surface,
            roi,
            graph,
            junction_node,
            local_bounds,
            panel,
            window_size,
            half_extent,
        )

    baseline_report = baseline["radius"]
    compensated_report = compensated["radius"]
    assert isinstance(baseline_report, dict) and isinstance(compensated_report, dict)
    font = _font()
    figure, axes = plt.subplots(1, 2, figsize=(10.91, 5.46), dpi=100)
    panels = (
        (
            baseline_panel,
            "A  未补偿：输入比例 1.00",
            baseline_report,
        ),
        (
            compensated_panel,
            "B  半径补偿：输入比例 0.91",
            compensated_report,
        ),
    )
    for axis, (panel, title, report) in zip(axes, panels, strict=True):
        axis.imshow(mpimg.imread(panel))
        axis.set_axis_off()
        axis.set_title(title, fontproperties=font, fontsize=13.5, pad=9, color="#171717")
        axis.text(
            0.5,
            -0.025,
            (
                f"中位有符号误差 {_metric_percent(report, 'median_signed_relative_error'):+.4f}%；"
                f"P95 绝对误差 {_metric_percent(report, 'p95_absolute_relative_error'):.4f}%"
            ),
            transform=axis.transAxes,
            ha="center",
            va="top",
            fontproperties=font,
            fontsize=9.5,
            color="#444444",
        )
    figure.suptitle(
        f"实际 ROI003274 分叉节点 {junction_node}：相同参数下的半径补偿前后对比",
        fontproperties=font,
        fontsize=15.5,
        y=0.975,
        color="#111111",
    )
    figure.subplots_adjust(left=0.015, right=0.985, bottom=0.075, top=0.89, wspace=0.03)
    output_path = output_dir / "junction_before_after.png"
    figure.savefig(output_path, facecolor="white")
    plt.close(figure)
    return output_path


def _compose_radius_before_after(
    *,
    baseline: dict[str, object],
    compensated: dict[str, object],
    output_dir: Path,
) -> tuple[Path, list[tuple[int, int]]]:
    baseline_samples = _radius_sample_map(baseline)
    compensated_samples = _radius_sample_map(compensated)
    keys = sorted(set(baseline_samples) & set(compensated_samples))
    if not keys:
        raise ValueError("The two model runs have no common real cross-section samples")

    baseline_source = np.asarray(
        [float(baseline_samples[key]["source_radius_um"]) for key in keys], dtype=float
    )
    compensated_source = np.asarray(
        [float(compensated_samples[key]["source_radius_um"]) for key in keys], dtype=float
    )
    if not np.allclose(baseline_source, compensated_source, rtol=0.0, atol=1.0e-10):
        raise ValueError("Source radii differ between paired validation samples")
    baseline_reconstructed = np.asarray(
        [float(baseline_samples[key]["reconstructed_radius_um"]) for key in keys],
        dtype=float,
    )
    compensated_reconstructed = np.asarray(
        [float(compensated_samples[key]["reconstructed_radius_um"]) for key in keys],
        dtype=float,
    )
    baseline_error = (baseline_reconstructed - baseline_source) / baseline_source * 100.0
    compensated_error = (
        (compensated_reconstructed - compensated_source) / compensated_source * 100.0
    )

    baseline_report = baseline["radius"]
    compensated_report = compensated["radius"]
    assert isinstance(baseline_report, dict) and isinstance(compensated_report, dict)
    font = _font()
    blue = "#0072B2"
    orange = "#D55E00"
    figure, axes = plt.subplots(2, 2, figsize=(13.5, 11.12), dpi=100)
    all_values = np.concatenate(
        (baseline_source, baseline_reconstructed, compensated_reconstructed)
    )
    low = max(0.0, float(np.min(all_values)) * 0.92)
    high = float(np.max(all_values)) * 1.06

    scatter_specs = (
        (
            axes[0, 0],
            baseline_reconstructed,
            blue,
            "A  未补偿：输入比例 1.00",
            baseline_report,
        ),
        (
            axes[0, 1],
            compensated_reconstructed,
            orange,
            "B  半径补偿：输入比例 0.91",
            compensated_report,
        ),
    )
    for axis, reconstructed, color, title, report in scatter_specs:
        axis.plot([low, high], [low, high], "--", color="#555555", linewidth=1.5)
        axis.scatter(
            baseline_source,
            reconstructed,
            s=38,
            color=color,
            edgecolor="white",
            linewidth=0.6,
            alpha=0.9,
        )
        axis.set_xlim(low, high)
        axis.set_ylim(low, high)
        axis.set_aspect("equal", adjustable="box")
        axis.set_title(title, fontproperties=font, fontsize=14, pad=11)
        axis.set_xlabel("源半径（μm）", fontproperties=font, fontsize=11)
        axis.set_ylabel("重建等效半径（μm）", fontproperties=font, fontsize=11)
        axis.grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.8)
        axis.text(
            0.04,
            0.95,
            (
                f"真实截面 n={len(keys)}\n"
                f"中位有符号误差：{_metric_percent(report, 'median_signed_relative_error'):+.4f}%\n"
                f"P95 绝对误差：{_metric_percent(report, 'p95_absolute_relative_error'):.4f}%"
            ),
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontproperties=font,
            fontsize=10,
            bbox={"boxstyle": "round,pad=0.35", "facecolor": "white", "alpha": 0.88, "edgecolor": "#BBBBBB"},
        )

    order = np.argsort(baseline_source)
    sample_order = np.arange(1, len(keys) + 1)
    axes[1, 0].axhline(0.0, color="#555555", linestyle="--", linewidth=1.2)
    axes[1, 0].plot(
        sample_order,
        baseline_error[order],
        "o-",
        color=blue,
        markersize=4.2,
        linewidth=1.3,
        label="未补偿（1.00）",
    )
    axes[1, 0].plot(
        sample_order,
        compensated_error[order],
        "o-",
        color=orange,
        markersize=4.2,
        linewidth=1.3,
        label="补偿后（0.91）",
    )
    axes[1, 0].set_title(
        "C  同一批真实截面的有符号误差",
        fontproperties=font,
        fontsize=14,
        pad=11,
    )
    axes[1, 0].set_xlabel("按源半径从小到大排列的截面", fontproperties=font, fontsize=11)
    axes[1, 0].set_ylabel("有符号相对误差（%）", fontproperties=font, fontsize=11)
    axes[1, 0].grid(True, color="#D9D9D9", linewidth=0.7, alpha=0.8)
    axes[1, 0].legend(prop=font, frameon=False, loc="best")

    absolute_errors = (np.abs(baseline_error), np.abs(compensated_error))
    boxes = axes[1, 1].boxplot(
        absolute_errors,
        positions=[1, 2],
        widths=0.48,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#222222", "linewidth": 1.6},
    )
    for box, color in zip(boxes["boxes"], (blue, orange), strict=True):
        box.set_facecolor(color)
        box.set_alpha(0.32)
        box.set_edgecolor(color)
    offsets = np.linspace(-0.13, 0.13, len(keys))
    axes[1, 1].scatter(
        1.0 + offsets,
        absolute_errors[0],
        s=22,
        color=blue,
        alpha=0.72,
        edgecolor="white",
        linewidth=0.4,
    )
    axes[1, 1].scatter(
        2.0 + offsets,
        absolute_errors[1],
        s=22,
        color=orange,
        alpha=0.72,
        edgecolor="white",
        linewidth=0.4,
    )
    axes[1, 1].scatter(
        [1, 2],
        [
            _metric_percent(baseline_report, "p95_absolute_relative_error"),
            _metric_percent(compensated_report, "p95_absolute_relative_error"),
        ],
        marker="D",
        s=75,
        color="#111111",
        label="P95",
        zorder=5,
    )
    axes[1, 1].set_xticks([1, 2], ["未补偿\n1.00", "补偿后\n0.91"], fontproperties=font)
    axes[1, 1].set_ylabel("绝对相对误差（%）", fontproperties=font, fontsize=11)
    axes[1, 1].set_title(
        "D  绝对误差分布与第 95 百分位",
        fontproperties=font,
        fontsize=14,
        pad=11,
    )
    axes[1, 1].grid(True, axis="y", color="#D9D9D9", linewidth=0.7, alpha=0.8)
    axes[1, 1].legend(prop=font, frameon=False, loc="upper right")

    figure.suptitle(
        "实际 ROI003274：源半径与重建截面半径的补偿前后对比",
        fontproperties=font,
        fontsize=18,
        y=0.975,
    )
    figure.text(
        0.5,
        0.015,
        "虚线表示重建半径与源半径完全相等；全部点均为两个正式运行共同获得的真实有效截面。",
        ha="center",
        va="bottom",
        fontproperties=font,
        fontsize=10.5,
        color="#444444",
    )
    figure.subplots_adjust(left=0.085, right=0.97, bottom=0.075, top=0.91, wspace=0.23, hspace=0.28)
    output_path = output_dir / "radius_before_after.png"
    figure.savefig(output_path, facecolor="white")
    plt.close(figure)
    return output_path, keys


def _generate_radius_validation_figures(
    *,
    baseline_run: Path,
    compensated_run: Path,
    roi,
    graph: nx.Graph,
    junction_node: int,
    local_bounds: tuple[float, ...],
    half_extent: float,
    output_dir: Path,
) -> tuple[Path, Path, Path]:
    baseline = _read_radius_validation_run(baseline_run)
    compensated = _read_radius_validation_run(compensated_run)
    baseline_summary = baseline["summary"]
    compensated_summary = compensated["summary"]
    assert isinstance(baseline_summary, dict) and isinstance(compensated_summary, dict)
    if baseline_summary.get("roi_id") != compensated_summary.get("roi_id"):
        raise ValueError("Baseline and compensated model runs refer to different ROIs")
    if str(baseline_summary.get("roi_id")) != str(roi.roi_id):
        raise ValueError("Radius-validation runs do not match the selected saved ROI")
    if not np.isclose(float(baseline_summary.get("radius_scale")), 1.0):
        raise ValueError("The baseline model run must use radius_scale=1.0")
    compensated_scale = float(compensated_summary.get("radius_scale"))
    if not (0.0 < compensated_scale < 1.0):
        raise ValueError("The compensated model run must use a radius scale between 0 and 1")

    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    junction_path = _compose_junction_before_after(
        baseline=baseline,
        compensated=compensated,
        roi=roi,
        graph=graph,
        junction_node=junction_node,
        local_bounds=local_bounds,
        half_extent=half_extent,
        output_dir=output_dir,
    )
    radius_path, paired_keys = _compose_radius_before_after(
        baseline=baseline,
        compensated=compensated,
        output_dir=output_dir,
    )

    validation_root = output_dir.parent
    geometry_dir = validation_root / "geometry"
    qc_dir = validation_root / "qc"
    geometry_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)
    baseline_surface_copy = geometry_dir / "lumen_surface_radius100_um.vtp"
    compensated_surface_copy = geometry_dir / "lumen_surface_um.vtp"
    baseline_qc_copy = qc_dir / "radius_fidelity_radius100.json"
    compensated_qc_copy = qc_dir / "radius_fidelity_radius091.json"
    shutil.copy2(Path(baseline["surface_path"]), baseline_surface_copy)
    shutil.copy2(Path(compensated["surface_path"]), compensated_surface_copy)
    shutil.copy2(Path(baseline["radius_path"]), baseline_qc_copy)
    shutil.copy2(Path(compensated["radius_path"]), compensated_qc_copy)

    baseline_report = baseline["radius"]
    compensated_report = compensated["radius"]
    assert isinstance(baseline_report, dict) and isinstance(compensated_report, dict)
    manifest = {
        "purpose": "Recreate thesis figures comparing unscaled and radius-compensated Ultraliser reconstruction",
        "roi_id": str(roi.roi_id),
        "junction_node": int(junction_node),
        "baseline": {
            "run_root": str(baseline["root"]),
            "radius_scale": float(baseline_summary["radius_scale"]),
            "surface_sha256": _sha256(Path(baseline["surface_path"])),
            "successful_sample_count": int(baseline_report["successful_sample_count"]),
            "median_signed_relative_error": float(baseline_report["median_signed_relative_error"]),
            "p95_absolute_relative_error": float(baseline_report["p95_absolute_relative_error"]),
        },
        "compensated": {
            "run_root": str(compensated["root"]),
            "radius_scale": compensated_scale,
            "surface_sha256": _sha256(Path(compensated["surface_path"])),
            "successful_sample_count": int(compensated_report["successful_sample_count"]),
            "median_signed_relative_error": float(compensated_report["median_signed_relative_error"]),
            "p95_absolute_relative_error": float(compensated_report["p95_absolute_relative_error"]),
        },
        "paired_real_cross_section_count": len(paired_keys),
        "paired_sample_keys": [
            {"branch_id": int(branch_id), "sample_index": int(sample_index)}
            for branch_id, sample_index in paired_keys
        ],
        "copied_reference_artifacts": {
            "baseline_surface": str(baseline_surface_copy),
            "compensated_surface": str(compensated_surface_copy),
            "baseline_radius_qc": str(baseline_qc_copy),
            "compensated_radius_qc": str(compensated_qc_copy),
        },
        "figures": {
            "junction_before_after": str(junction_path),
            "radius_before_after": str(radius_path),
        },
    }
    manifest_path = validation_root / "visualization_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    return junction_path, radius_path, manifest_path


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.crop_half_extent_um <= 0.0:
        raise ValueError("--crop-half-extent-um must be positive")
    if args.panel_width < 400 or args.panel_height < 300:
        raise ValueError("Panel dimensions are too small for a readable 3-D figure")

    model_settings = load_swc_stl_yaml_config(args.config, project_root=PROJECT_ROOT)
    config = model_settings.lumen
    sampling_run = resolve_sampling_run(args.sampling_run, project_root=PROJECT_ROOT)
    roi = _select_roi(sampling_run, args.roi_anchor)
    graph = _build_graph(roi)
    paths = _branch_paths(graph)
    junction_node = _select_junction(roi, graph, args.junction_node)
    junction_position = np.asarray(roi.local_node_positions_um[junction_node], dtype=float)
    local_bounds = _crop_bounds(junction_position, args.crop_half_extent_um)

    radius_runs = (args.radius_baseline_run, args.radius_compensated_run)
    if any(path is not None for path in radius_runs) and not all(
        path is not None for path in radius_runs
    ):
        raise ValueError(
            "--radius-baseline-run and --radius-compensated-run must be supplied together"
        )
    if args.radius_comparison_only:
        if not all(path is not None for path in radius_runs):
            raise ValueError(
                "--radius-comparison-only requires both radius-validation model runs"
            )
        junction_path, radius_path, manifest_path = _generate_radius_validation_figures(
            baseline_run=args.radius_baseline_run,
            compensated_run=args.radius_compensated_run,
            roi=roi,
            graph=graph,
            junction_node=junction_node,
            local_bounds=local_bounds,
            half_extent=float(args.crop_half_extent_um),
            output_dir=args.radius_validation_output_dir,
        )
        print(f"ROI: {roi.roi_id}")
        print(f"Junction: local node {junction_node} at {junction_position.tolist()} um")
        print(f"Junction comparison: {junction_path}")
        print(f"Radius comparison: {radius_path}")
        print(f"Validation manifest: {manifest_path}")
        return 0

    surface_path = args.surface.resolve()
    if not surface_path.is_file():
        raise FileNotFoundError(f"Accepted micrometre surface not found: {surface_path}")
    surface = (
        pv.read(surface_path)
        .extract_surface(algorithm="dataset_surface")
        .triangulate()
        .clean()
    )
    local_surface = surface.clip_box(pv.Box(bounds=local_bounds), invert=False).extract_surface(
        algorithm="dataset_surface"
    )

    radius_scale = float(config.ultraliser.radius_scale)
    voxels_per_micron = float(config.ultraliser.voxels_per_micron)
    feed_edge_radii = np.asarray(roi.local_edge_radius_um, dtype=float) * radius_scale
    occupancy_grid, occupied, occupancy_boundary, occupancy_details = _capsule_union_occupancy(
        np.asarray(roi.local_edge_points_um, dtype=float),
        feed_edge_radii,
        local_bounds,
        voxels_per_micron,
    )
    if occupied.n_cells == 0:
        raise RuntimeError("The explanatory local occupancy grid is empty")

    output_dir = args.output_dir.resolve()
    panels_dir = output_dir / "panels"
    panels_dir.mkdir(parents=True, exist_ok=True)
    panel_paths = [
        panels_dir / "01_saved_roi_centerline_radius.png",
        panels_dir / "02_separate_branch_envelopes.png",
        panels_dir / "03_shared_occupancy_grid.png",
        panels_dir / "04_accepted_unified_surface.png",
    ]
    window_size = (int(args.panel_width), int(args.panel_height))
    centerline = _centerline_mesh(roi, graph, radius_scale)
    incident = _incident_paths(paths, junction_node)
    if len(incident) < 3:
        raise RuntimeError(
            f"Junction {junction_node} is represented by fewer than three branch paths"
        )

    _render_overview_panel(
        roi,
        graph,
        paths,
        centerline,
        junction_node,
        local_bounds,
        panel_paths[0],
        window_size,
    )
    _render_separate_branches_panel(
        roi,
        incident,
        radius_scale,
        junction_node,
        local_bounds,
        panel_paths[1],
        window_size,
        args.crop_half_extent_um,
    )
    _render_occupancy_panel(
        occupancy_boundary,
        roi,
        graph,
        junction_node,
        local_bounds,
        panel_paths[2],
        window_size,
        args.crop_half_extent_um,
    )
    _render_surface_panel(
        local_surface,
        roi,
        graph,
        junction_node,
        local_bounds,
        panel_paths[3],
        window_size,
        args.crop_half_extent_um,
    )

    figure_path = output_dir / "swc_to_stl_shared_domain.png"
    _compose_figure(
        panel_paths,
        figure_path,
        roi=roi,
        branch_count=len(paths),
        junction_node=junction_node,
        radius_scale=radius_scale,
        voxels_per_micron=voxels_per_micron,
        surface_cells=int(surface.n_cells),
    )

    centerline.save(output_dir / "source_centerline_and_radius.vtp", binary=True)
    occupied.save(output_dir / "explanatory_shared_occupancy.vtu", binary=True)
    occupancy_boundary.save(output_dir / "explanatory_occupancy_boundary.vtp", binary=True)
    local_surface.save(output_dir / "accepted_ultraliser_surface_local_view.vtp", binary=True)

    junction_feed_radius = float(roi.local_node_radius_um[junction_node]) * radius_scale
    metadata: dict[str, object] = {
        "purpose": "Explain shared-domain occupancy and unified-surface extraction",
        "production_equivalence_boundary": (
            "The shared occupancy is an explanatory variable-radius capsule union built from "
            "the production input data and spacing; it is not an Ultraliser-exported intermediate."
        ),
        "sampling_run": str(sampling_run),
        "roi_id": roi.roi_id,
        "anchor_id": int(roi.anchor_id),
        "node_count": int(len(roi.local_node_ids)),
        "edge_count": int(len(roi.local_edges)),
        "branch_count": int(len(paths)),
        "cut_port_count": int(len(roi.cut_ports)),
        "junction_node": int(junction_node),
        "junction_degree": int(graph.degree(junction_node)),
        "junction_position_um": [float(value) for value in junction_position],
        "junction_source_radius_um": float(roi.local_node_radius_um[junction_node]),
        "junction_feed_radius_um": junction_feed_radius,
        "radius_scale": radius_scale,
        "feed_radius_equation": "feed_radius_um = source_radius_um * radius_scale",
        "voxels_per_micron": voxels_per_micron,
        "explanatory_occupancy_method": (
            "logical union of variable-radius line-segment capsules sampled at voxel centres"
        ),
        "explanatory_occupancy": occupancy_details,
        "accepted_surface": {
            "path": str(surface_path),
            "point_count": int(surface.n_points),
            "triangle_count": int(surface.n_cells),
            "component_count": int(surface.connectivity().cell_data["RegionId"].max() + 1),
            "watertight": bool(surface.is_manifold),
            "local_view_point_count": int(local_surface.n_points),
            "local_view_triangle_count": int(local_surface.n_cells),
        },
        "outputs": {
            "overview_figure": str(figure_path),
            "panels": [str(path) for path in panel_paths],
        },
    }
    metadata_path = output_dir / "visualization_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    _write_readme(output_dir / "README.md", metadata)

    processing_figure_path = _generate_processing_stage_visualization(
        roi=roi,
        graph=graph,
        config=config,
        occupancy_grid=occupancy_grid,
        junction_node=junction_node,
        local_bounds=local_bounds,
        half_extent=float(args.crop_half_extent_um),
        window_size=window_size,
        output_dir=args.processing_output_dir.resolve(),
        ultraliser_root=args.ultraliser_root.resolve(),
        accepted_surface_path=surface_path,
    )
    voxelization_figure_path, xyz_slice_figure_path = _generate_voxelization_explanation(
        roi=roi,
        graph=graph,
        config=config,
        junction_node=junction_node,
        output_dir=args.voxelization_output_dir.resolve(),
        window_size=window_size,
    )

    validation_paths: tuple[Path, Path, Path] | None = None
    if all(path is not None for path in radius_runs):
        validation_paths = _generate_radius_validation_figures(
            baseline_run=args.radius_baseline_run,
            compensated_run=args.radius_compensated_run,
            roi=roi,
            graph=graph,
            junction_node=junction_node,
            local_bounds=local_bounds,
            half_extent=float(args.crop_half_extent_um),
            output_dir=args.radius_validation_output_dir,
        )

    print(f"ROI: {roi.roi_id}")
    print(f"Junction: local node {junction_node} at {junction_position.tolist()} um")
    print(
        "Explanatory occupancy: "
        f"{occupancy_details['occupied_cell_count']} / "
        f"{occupancy_details['total_cell_count']} cells"
    )
    print(f"Figure: {figure_path}")
    print(f"Processing-stage figure: {processing_figure_path}")
    print(f"Voxelization figure: {voxelization_figure_path}")
    print(f"XYZ slice figure: {xyz_slice_figure_path}")
    if validation_paths is not None:
        print(f"Junction comparison: {validation_paths[0]}")
        print(f"Radius comparison: {validation_paths[1]}")
        print(f"Validation manifest: {validation_paths[2]}")
    print(f"Output directory: {output_dir}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, KeyError, RuntimeError, ValueError) as exc:
        print(f"Documentation visualisation failed: {exc}", file=sys.stderr)
        raise SystemExit(1) from exc
