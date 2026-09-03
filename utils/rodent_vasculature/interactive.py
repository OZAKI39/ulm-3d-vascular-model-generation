"""Interactive Figure 2(a)-style rendering for mouse-brain vascular blocks.

The preprocessed main-network volume is rendered as white vessels on black,
matching the visual language of Figure 2(a). SWC centerlines, critical nodes, and arrows
are overlays; arrow direction is the annotation relationship
``parent_id node -> current node`` and is not a measured flow direction.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .model import DirectedVascularGraph
from .tiff_io import load_normalized_volume, load_tiff_volume
from .visualization import sample_direction_arrows
from ..sampling.sampling_io import load_sampling_display_rois
from ..sampling.sampling_types import ROIRecord


UI_FONT_FAMILY = "arial"
COORDINATE_TICK_FONT_SIZE = 11
COORDINATE_TITLE_FONT_SIZE = 12
ORIENTATION_AXIS_FONT_SIZE = 12
LEGEND_FONT_SIZE = 11


@dataclass(frozen=True, slots=True)
class InteractiveSceneGeometry:
    """Small, renderer-independent subset of a directed vascular graph."""

    branches_um: tuple[np.ndarray, ...]
    arrow_points_um: np.ndarray
    arrow_vectors_xyz: np.ndarray
    critical_points_um: np.ndarray
    critical_roles: tuple[str, ...]

    @property
    def branch_count(self) -> int:
        return len(self.branches_um)

    @property
    def arrow_count(self) -> int:
        return len(self.arrow_points_um)

@dataclass(frozen=True, slots=True)
class Figure2aArtifacts:
    screenshot_path: Path
    manifest_path: Path


def _style_text_property(
    text_property: Any,
    *,
    font_size: int,
    bold: bool,
) -> None:
    """Apply the shared readable Arial typography to a VTK text property."""

    text_property.SetFontFamilyToArial()
    text_property.SetFontSize(font_size)
    text_property.SetBold(bold)


def _style_orientation_axes(axes_actor: Any) -> None:
    for caption_getter in (
        axes_actor.GetXAxisCaptionActor2D,
        axes_actor.GetYAxisCaptionActor2D,
        axes_actor.GetZAxisCaptionActor2D,
    ):
        _style_text_property(
            caption_getter().GetCaptionTextProperty(),
            font_size=ORIENTATION_AXIS_FONT_SIZE,
            bold=True,
        )


def _style_bounds_axes(bounds_actor: Any) -> None:
    for axis_index in range(3):
        _style_text_property(
            bounds_actor.GetLabelTextProperty(axis_index),
            font_size=COORDINATE_TICK_FONT_SIZE,
            bold=False,
        )
        _style_text_property(
            bounds_actor.GetTitleTextProperty(axis_index),
            font_size=COORDINATE_TITLE_FONT_SIZE,
            bold=True,
        )


def _style_legend(legend_actor: Any) -> None:
    _style_text_property(
        legend_actor.GetEntryTextProperty(),
        font_size=LEGEND_FONT_SIZE,
        bold=False,
    )


def scene_geometry_from_graph(
    result: DirectedVascularGraph,
    *,
    max_arrows: int,
) -> InteractiveSceneGeometry:
    arrow_points, arrow_vectors, _ = sample_direction_arrows(result, max_arrows)
    ordered_nodes = sorted(result.junction_graph.nodes(data=True), key=lambda item: int(item[0]))
    critical_points = np.asarray(
        [
            (float(data["x_um"]), float(data["y_um"]), float(data["z_um"]))
            for _, data in ordered_nodes
        ],
        dtype=float,
    )
    if not len(critical_points):
        critical_points = np.empty((0, 3), dtype=float)
    branches = tuple(branch.derived_points_um.copy() for branch in result.branches)
    return InteractiveSceneGeometry(
        branches_um=branches,
        arrow_points_um=arrow_points,
        arrow_vectors_xyz=arrow_vectors,
        critical_points_um=critical_points,
        critical_roles=tuple(str(data["role"]) for _, data in ordered_nodes),
    )


def _model_bounds(
    volume_zyx: np.ndarray | None,
    spacing_xyz_um: tuple[float, float, float],
    geometry: InteractiveSceneGeometry | None = None,
) -> tuple[float, float, float, float, float, float]:
    if volume_zyx is not None:
        shape_xyz = np.asarray(volume_zyx.shape[::-1], dtype=float)
        maximum = np.maximum(shape_xyz - 1.0, 0.0) * np.asarray(spacing_xyz_um, dtype=float)
        return (0.0, float(maximum[0]), 0.0, float(maximum[1]), 0.0, float(maximum[2]))
    point_arrays: list[np.ndarray] = []
    if geometry is not None:
        point_arrays.extend(branch for branch in geometry.branches_um if len(branch))
        if len(geometry.critical_points_um):
            point_arrays.append(geometry.critical_points_um)
    if not point_arrays:
        return (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    points = np.concatenate(point_arrays, axis=0)
    lower = np.min(points, axis=0)
    upper = np.max(points, axis=0)
    upper = np.where(upper > lower, upper, lower + 1.0)
    return (
        float(lower[0]), float(upper[0]),
        float(lower[1]), float(upper[1]),
        float(lower[2]), float(upper[2]),
    )


def _normalized_volume(volume_zyx: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    values = np.asarray(volume_zyx, dtype=np.float32)
    finite = values[np.isfinite(values)]
    positive = finite[finite > 0]
    if positive.size:
        lower = float(np.percentile(positive, 12.0))
        upper = float(np.percentile(positive, 99.8))
    elif finite.size:
        lower, upper = float(np.min(finite)), float(np.max(finite))
    else:
        raise ValueError("The TIFF volume contains no finite intensity values")
    if upper <= lower:
        lower = 0.0
        upper = max(float(np.max(finite)), 1.0)
    normalized = np.nan_to_num((values - lower) / (upper - lower), nan=0.0)
    normalized = np.clip(normalized, 0.0, 1.0)
    return normalized, {"source_lower": lower, "source_upper": upper}


def _volume_grid(volume_zyx: np.ndarray, spacing_xyz_um: tuple[float, float, float]) -> Any:
    import pyvista as pv

    z_size, y_size, x_size = volume_zyx.shape
    grid = pv.ImageData(
        dimensions=(x_size, y_size, z_size),
        spacing=spacing_xyz_um,
        origin=(0.0, 0.0, 0.0),
    )
    xyz = np.transpose(volume_zyx, (2, 1, 0))
    grid.point_data["intensity"] = xyz.ravel(order="F")
    return grid


def _branch_mesh(
    branches_um: tuple[np.ndarray, ...], radius_um: tuple[np.ndarray, ...] | None = None
) -> Any | None:
    import pyvista as pv

    usable = [np.asarray(points, dtype=float) for points in branches_um if len(points) >= 2]
    if not usable:
        return None
    points = np.concatenate(usable, axis=0)
    cells: list[np.ndarray] = []
    offset = 0
    for branch in usable:
        indices = np.arange(offset, offset + len(branch), dtype=np.int64)
        cells.append(np.concatenate(([len(branch)], indices)))
        offset += len(branch)
    mesh = pv.PolyData(points)
    mesh.lines = np.concatenate(cells)
    if radius_um is not None:
        mesh.point_data["radius_um"] = np.concatenate(radius_um)
    return mesh


def _add_physical_coordinate_axes(
    plotter: Any,
    bounds: tuple[float, float, float, float, float, float],
    *,
    add_orientation_axes: bool = True,
) -> None:
    """Show direction and numerical bounds in the shared physical coordinate system."""

    import pyvista as pv

    if add_orientation_axes:
        axes_actor = plotter.add_axes(
            xlabel="X (um)",
            ylabel="Y (um)",
            zlabel="Z (um)",
            color="white",
        )
        _style_orientation_axes(axes_actor)
    reference_box = pv.Box(bounds=bounds)
    bounds_actor = plotter.show_bounds(
        mesh=reference_box,
        axes_ranges=bounds,
        xtitle="X (um)",
        ytitle="Y (um)",
        ztitle="Z (um)",
        n_xlabels=4,
        n_ylabels=4,
        n_zlabels=5,
        fmt="%.1f",
        font_size=COORDINATE_TICK_FONT_SIZE,
        font_family=UI_FONT_FAMILY,
        color="#BFC7D5",
        grid="back",
        location="outer",
        ticks="outside",
        all_edges=True,
        use_3d_text=False,
    )
    _style_bounds_axes(bounds_actor)


def _set_full_scene_title(plotter: Any, *, sampling_available: bool) -> None:
    detail = (
        "Representative connected ROIs\n"
        "A: all ROI candidates | R/S: selected | C: next cluster | left-click: inspect"
        if sampling_available
        else "Orange arrows = SWC parent -> current"
    )
    plotter.add_text(
        f"Preprocessed main vascular network\n{detail}",
        position="upper_left",
        font_size=11,
        color="white",
        font=UI_FONT_FAMILY,
        name="full_scene_title",
    )


def _add_full_scene(
    plotter: Any,
    volume_zyx: np.ndarray | None,
    geometry: InteractiveSceneGeometry,
    *,
    spacing_xyz_um: tuple[float, float, float],
    volume_opacity: float,
    sample_id: str,
    sampling_available: bool = False,
) -> dict[str, Any]:
    import pyvista as pv

    plotter.set_background("black")
    intensity_window: dict[str, float] | None = None
    if volume_zyx is not None:
        normalized, intensity_window = _normalized_volume(volume_zyx)
        grid = _volume_grid(normalized, spacing_xyz_um)
        opacity = np.asarray([0.0, 0.0, 0.005, 0.02, 0.08, 0.22, 0.58, 1.0])
        opacity *= float(volume_opacity)
        plotter.add_volume(
            grid,
            scalars="intensity",
            cmap="gray",
            opacity=opacity,
            clim=(0.0, 1.0),
            shade=True,
            ambient=0.20,
            diffuse=0.85,
            specular=0.25,
            show_scalar_bar=False,
            pickable=False,
        )

    branch_mesh = _branch_mesh(geometry.branches_um)
    if branch_mesh is not None:
        plotter.add_mesh(
            branch_mesh,
            color="#35D6E3",
            line_width=2.4,
            opacity=0.86,
            label="SWC centerline",
            pickable=False,
        )

    if geometry.arrow_count:
        arrows = pv.PolyData(geometry.arrow_points_um)
        arrows.point_data["parent_to_current"] = geometry.arrow_vectors_xyz
        bounds = _model_bounds(volume_zyx, spacing_xyz_um, geometry)
        diagonal = float(
            np.linalg.norm(
                np.asarray((bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]))
            )
        )
        glyphs = arrows.glyph(
            orient="parent_to_current",
            scale=False,
            factor=max(2.0, diagonal * 0.022),
        )
        plotter.add_mesh(
            glyphs,
            color="#FF9F1C",
            label="parent -> current",
            lighting=True,
            pickable=False,
        )

    role_styles = {
        "inferred_inlet": ("#4ADE80", 12.0, "structural root (not flow inlet)"),
        "inferred_outlet": ("#F87171", 10.0, "structural leaf (not flow outlet)"),
        "divergence_junction": ("#FACC15", 12.0, "divergence junction"),
        "convergence_junction": ("#C084FC", 12.0, "convergence junction"),
    }
    roles = np.asarray(geometry.critical_roles, dtype=object)
    for role, (color, size, display_label) in role_styles.items():
        selected = geometry.critical_points_um[roles == role]
        if len(selected):
            plotter.add_points(
                selected,
                color=color,
                point_size=size,
                render_points_as_spheres=True,
                label=display_label,
                pickable=False,
            )

    _set_full_scene_title(plotter, sampling_available=sampling_available)
    model_bounds = _model_bounds(volume_zyx, spacing_xyz_um, geometry)
    legend = plotter.add_legend(
        bcolor="#101010",
        face=None,
        size=(0.255, 0.235),
        font_family=UI_FONT_FAMILY,
    )
    _style_legend(legend)
    plotter.view_isometric()
    plotter.camera.zoom(1.12)
    _add_physical_coordinate_axes(plotter, model_bounds)
    return {
        "sample_id": sample_id,
        "volume_shape_zyx": list(map(int, volume_zyx.shape)) if volume_zyx is not None else None,
        "optional_background_volume_available": volume_zyx is not None,
        "spacing_xyz_um": list(map(float, spacing_xyz_um)),
        "branch_count": geometry.branch_count,
        "arrow_count": geometry.arrow_count,
        "critical_node_count": len(geometry.critical_points_um),
        "global_coordinate_bounds_xyz_um": list(model_bounds),
        "coordinate_units": "um",
        "interface_style": {
            "font_family": "Arial",
            "coordinate_tick_font_size": COORDINATE_TICK_FONT_SIZE,
            "coordinate_title_font_size": COORDINATE_TITLE_FONT_SIZE,
            "legend_font_size": LEGEND_FONT_SIZE,
            "overlay_emphasis": "enhanced in both viewports",
        },
        "intensity_window": intensity_window,
        "direction_rule": "SWC parent_id node -> current node",
        "direction_is_measured_flow": False,
        "interaction": {
            "selection": "left-click a sampled ROI cube in the full-model viewport"
            if sampling_available
            else "rotate, pan, and zoom the preprocessed full model",
            "result": "the right viewport is replaced by the selected connected ROI"
            if sampling_available
            else "the full-model view remains interactive",
        },
    }


def _sampling_color(cluster_id: int) -> str:
    palette = (
        "#4DD0E1", "#FF8A65", "#FFD54F", "#AB47BC", "#66BB6A",
        "#42A5F5", "#EC407A", "#9CCC65", "#FFA726", "#7E57C2",
    )
    return palette[int(cluster_id) % len(palette)] if cluster_id >= 0 else "#B0BEC5"


def _add_sampling_roi_scene(
    plotter: Any,
    roi: ROIRecord,
    *,
    add_orientation_axes: bool = True,
) -> None:
    """Render one saved connected ROI without running sampling or clustering."""

    import pyvista as pv

    plotter.set_background("black")
    bounds = (
        roi.bbox_min_um[0], roi.bbox_max_um[0],
        roi.bbox_min_um[1], roi.bbox_max_um[1],
        roi.bbox_min_um[2], roi.bbox_max_um[2],
    )
    box = pv.Box(bounds=bounds)
    plotter.add_mesh(box, color="#FF5C5C", opacity=0.025, pickable=False)
    plotter.add_mesh(
        box,
        style="wireframe",
        color="#FF5C5C",
        line_width=4.0,
        label="spatial ROI boundary",
        pickable=False,
    )
    branches = tuple(np.asarray(edge, dtype=float) for edge in roi.local_edge_points_um)
    radii = tuple(np.asarray(radius, dtype=float) for radius in roi.local_edge_radius_um)
    mesh = _branch_mesh(branches, radii)
    if mesh is not None:
        tubes = mesh.tube(
            scalars="radius_um",
            absolute=True,
            radius_factor=1.0,
            n_sides=12,
            capping=True,
        )
        plotter.add_mesh(
            tubes,
            color="#F5F7FA",
            smooth_shading=True,
            opacity=0.88,
            label="source SWC radius",
            pickable=False,
        )
        plotter.add_mesh(
            mesh,
            color="#51E5FF",
            line_width=3.0,
            label="connected ROI centerline",
            pickable=False,
        )
    if len(roi.local_edge_points_um):
        starts = roi.local_edge_points_um[:, 0]
        ends = roi.local_edge_points_um[:, 1]
        vectors = ends - starts
        lengths = np.linalg.norm(vectors, axis=1)
        valid = lengths > 1.0e-12
        arrow_points = (starts[valid] + ends[valid]) * 0.5
        directions = vectors[valid] / lengths[valid, None]
        if len(arrow_points):
            arrows = pv.PolyData(arrow_points)
            arrows.point_data["parent_to_current"] = directions
            diagonal = float(np.linalg.norm(np.asarray(roi.bbox_size_um)))
            glyphs = arrows.glyph(
                orient="parent_to_current",
                scale=False,
                factor=max(2.0, diagonal * 0.055),
            )
            plotter.add_mesh(
                glyphs,
                color="#FF9F1C",
                label="global parent -> current",
                pickable=False,
            )
    if roi.true_terminal_local_ids:
        points = roi.local_node_positions_um[list(roi.true_terminal_local_ids)]
        plotter.add_points(
            points,
            color="#4ADE80",
            point_size=19,
            render_points_as_spheres=True,
            label="TRUE_TERMINAL",
            pickable=False,
        )
    if roi.cut_ports:
        points = np.asarray([port.intersection_position_um for port in roi.cut_ports], dtype=float)
        plotter.add_points(
            points,
            color="#F87171",
            point_size=21,
            render_points_as_spheres=True,
            label="CUT_PORT",
            pickable=False,
        )
    radius = roi.radius_features
    structure = roi.structural_features
    plotter.add_text(
        f"{roi.roi_id} | cluster {roi.cluster_id}\n"
        f"nodes {roi.node_count} | branches {roi.branch_count} | bifurcations {roi.bifurcation_count}\n"
        f"radius P10/P25/P50/P75/P90: "
        f"{radius.get('r10', float('nan')):.2f} / {radius.get('r25', float('nan')):.2f} / "
        f"{radius.get('r50', float('nan')):.2f} / {radius.get('r75', float('nan')):.2f} / "
        f"{radius.get('r90', float('nan')):.2f} um\n"
        f"length {structure.get('total_vessel_length_um', roi.retained_component_length_um):.1f} um | "
        f"cycle rank {structure.get('cycle_rank', float('nan')):.0f} | "
        f"true terminals {roi.true_terminal_count} | cut ports {roi.cut_port_count}",
        position="upper_left",
        font_size=10,
        color="white",
        font=UI_FONT_FAMILY,
        name="sampling_roi_information",
    )
    legend = plotter.add_legend(
        bcolor="#101010",
        face=None,
        size=(0.285, 0.22),
        loc="lower right",
        font_family=UI_FONT_FAMILY,
    )
    _style_legend(legend)
    plotter.view_isometric()
    plotter.reset_camera(bounds=bounds)
    _add_physical_coordinate_axes(
        plotter,
        bounds,
        add_orientation_axes=add_orientation_axes,
    )


def _add_sampling_boxes(
    plotter: Any,
    rois: tuple[ROIRecord, ...],
    indices: list[int],
    *,
    mode_label: str,
) -> tuple[list[Any], dict[str, int]]:
    import pyvista as pv

    actors: list[Any] = []
    actor_to_index: dict[str, int] = {}
    label_points: list[np.ndarray] = []
    labels: list[str] = []
    for index in indices:
        roi = rois[index]
        bounds = (
            roi.bbox_min_um[0], roi.bbox_max_um[0],
            roi.bbox_min_um[1], roi.bbox_max_um[1],
            roi.bbox_min_um[2], roi.bbox_max_um[2],
        )
        color = _sampling_color(roi.cluster_id)
        actor = plotter.add_mesh(
            pv.Box(bounds=bounds),
            color=color,
            opacity=0.16 if roi.is_representative else 0.05,
            show_edges=True,
            edge_color=color,
            line_width=4.0 if roi.is_representative else 1.6,
            pickable=True,
            reset_camera=False,
            name=f"sampling_roi_pick_{index}",
        )
        actors.append(actor)
        actor_to_index[actor.memory_address] = index
        if roi.is_representative or mode_label.startswith("cluster"):
            label_points.append(np.asarray(roi.bbox_center_um, dtype=float))
            labels.append(f"C{roi.cluster_id}\nR{roi.selection_rank if roi.selection_rank > 0 else '-'}")
    if label_points:
        actors.append(
            plotter.add_point_labels(
                np.asarray(label_points),
                labels,
                font_size=9,
                font_family=UI_FONT_FAMILY,
                text_color="white",
                shape_color="#202020",
                shape_opacity=0.70,
                always_visible=True,
                pickable=False,
                reset_camera=False,
                name="sampling_roi_labels",
            )
        )
    actors.append(
        plotter.add_text(
            f"Sampling layer: {mode_label} ({len(indices)} boxes)",
            position="lower_left",
            font_size=10,
            color="white",
            font=UI_FONT_FAMILY,
            name="sampling_layer_mode",
        )
    )
    return actors, actor_to_index


def _add_sampling_active_outline(plotter: Any, roi: ROIRecord) -> Any:
    import pyvista as pv

    bounds = (
        roi.bbox_min_um[0], roi.bbox_max_um[0],
        roi.bbox_min_um[1], roi.bbox_max_um[1],
        roi.bbox_min_um[2], roi.bbox_max_um[2],
    )
    return plotter.add_mesh(
        pv.Box(bounds=bounds),
        style="wireframe",
        color="white",
        line_width=5.0,
        pickable=False,
        reset_camera=False,
        name="active_sampling_roi_outline",
    )


def _install_sampling_layer(plotter: Any, rois: tuple[ROIRecord, ...]) -> dict[str, Any] | None:
    """Install display-only mode switching and ROI inspection callbacks."""

    if not rois:
        return None
    selected_indices = [index for index, roi in enumerate(rois) if roi.is_representative]
    default_index = min(
        selected_indices or list(range(len(rois))),
        key=lambda index: rois[index].selection_rank if rois[index].selection_rank > 0 else index + 100000,
    )
    cluster_ids = sorted({roi.cluster_id for roi in rois if roi.cluster_id >= 0})
    overlay_actors: list[Any] = []
    actor_to_index: dict[str, int] = {}
    active_outline: list[Any | None] = [None]
    cluster_position = [0]
    interaction_busy = [False]

    def redraw(indices: list[int], mode_label: str) -> None:
        if interaction_busy[0]:
            return
        interaction_busy[0] = True
        try:
            plotter.subplot(0, 0)
            for actor in overlay_actors:
                plotter.remove_actor(actor, reset_camera=False, render=False)
            overlay_actors.clear()
            new_actors, mapping = _add_sampling_boxes(
                plotter,
                rois,
                indices,
                mode_label=mode_label,
            )
            overlay_actors.extend(new_actors)
            actor_to_index.clear()
            actor_to_index.update(mapping)
            plotter.render()
        finally:
            plotter.subplot(0, 0)
            interaction_busy[0] = False

    def select_roi(actor: Any) -> None:
        if interaction_busy[0]:
            return
        index = actor_to_index.get(getattr(actor, "memory_address", ""))
        if index is None:
            return
        interaction_busy[0] = True
        try:
            plotter.subplot(0, 0)
            if active_outline[0] is not None:
                plotter.remove_actor(active_outline[0], reset_camera=False, render=False)
            active_outline[0] = _add_sampling_active_outline(plotter, rois[index])
            plotter.subplot(0, 1)
            plotter.renderer.clear_actors()
            _add_sampling_roi_scene(
                plotter,
                rois[index],
                add_orientation_axes=False,
            )
            plotter.subplot(0, 0)
            plotter.render()
        finally:
            plotter.subplot(0, 0)
            interaction_busy[0] = False

    def show_all() -> None:
        redraw(list(range(len(rois))), "all candidates")

    def show_selected() -> None:
        redraw(selected_indices or list(range(len(rois))), "selected representatives")

    def show_next_cluster() -> None:
        if not cluster_ids:
            return
        cluster_id = cluster_ids[cluster_position[0] % len(cluster_ids)]
        cluster_position[0] += 1
        redraw(
            [index for index, roi in enumerate(rois) if roi.cluster_id == cluster_id],
            f"cluster {cluster_id}",
        )

    plotter.subplot(0, 0)
    redraw(selected_indices or list(range(len(rois))), "selected representatives")
    active_outline[0] = _add_sampling_active_outline(plotter, rois[default_index])
    key_actions = {
        "a": show_all,
        "r": show_selected,
        "s": show_selected,
        "c": show_next_cluster,
    }
    for key, callback in key_actions.items():
        for variant in (key, key.upper()):
            plotter.clear_events_for_key(variant)
            plotter.add_key_event(variant, callback)
    plotter.enable_mesh_picking(
        callback=select_roi,
        show=False,
        show_message=False,
        left_clicking=True,
        use_actor=True,
    )
    return {
        "show_all": show_all,
        "show_selected": show_selected,
        "show_next_cluster": show_next_cluster,
        "select_roi": select_roi,
        "actor_to_index": actor_to_index,
    }


def render_figure2a_scene(
    volume_zyx: np.ndarray | None,
    geometry: InteractiveSceneGeometry,
    *,
    spacing_xyz_um: tuple[float, float, float],
    sample_id: str,
    volume_opacity: float = 0.32,
    window_size: tuple[int, int] = (1800, 900),
    screenshot_path: Path | None = None,
    show: bool = False,
    sampling_rois: tuple[ROIRecord, ...] = (),
) -> dict[str, Any]:
    """Render a static preview or open the native interactive PyVista window."""

    import pyvista as pv

    if sampling_rois:
        plotter = pv.Plotter(
            shape=(1, 2),
            border=True,
            border_color="#606060",
            off_screen=not show,
            window_size=window_size,
        )
        plotter.subplot(0, 0)
    else:
        plotter = pv.Plotter(off_screen=not show, window_size=window_size)
    metadata = _add_full_scene(
        plotter,
        volume_zyx,
        geometry,
        spacing_xyz_um=spacing_xyz_um,
        volume_opacity=volume_opacity,
        sample_id=sample_id,
        sampling_available=bool(sampling_rois),
    )
    if sampling_rois:
        selected = [roi for roi in sampling_rois if roi.is_representative]
        default_roi = min(
            selected or list(sampling_rois),
            key=lambda roi: roi.selection_rank if roi.selection_rank > 0 else 100000,
        )
        plotter.subplot(0, 1)
        _add_sampling_roi_scene(plotter, default_roi)
        plotter.subplot(0, 0)
        metadata["sampling_layer"] = {
            "candidate_count": len(sampling_rois),
            "selected_count": len(selected),
            "cluster_ids": sorted({roi.cluster_id for roi in sampling_rois}),
            "default_roi_id": default_roi.roi_id,
            "display_modes": [
                "all ROI candidates",
                "selected representatives",
                "cluster X",
                "ROI X",
            ],
            "core_recomputed_in_ui": False,
        }
        metadata["available_display_layers"] = ["representative connected ROIs"]
        metadata["default_display_layer"] = "representative connected ROIs"
        if show:
            _install_sampling_layer(plotter, sampling_rois)
            metadata["display_switch_keys"] = {
                "R": "selected representative ROIs",
                "A": "all candidate ROIs",
                "S": "selected representative ROIs",
                "C": "next ROI cluster",
            }
            metadata["interaction"] = {
                "selection": "choose an ROI display mode, then left-click a colored box",
                "result": "the right viewport shows the selected connected ROI",
            }
        else:
            indices = [index for index, roi in enumerate(sampling_rois) if roi.is_representative]
            _add_sampling_boxes(
                plotter,
                sampling_rois,
                indices or list(range(len(sampling_rois))),
                mode_label="selected representatives",
            )
            _add_sampling_active_outline(plotter, default_roi)
    if screenshot_path is not None:
        screenshot_path.parent.mkdir(parents=True, exist_ok=True)
    if show:
        plotter.show(
            title=f"Mouse brain vasculature - {sample_id}",
            screenshot=str(screenshot_path) if screenshot_path else None,
            auto_close=True,
        )
    else:
        plotter.show(
            screenshot=str(screenshot_path) if screenshot_path else None,
            auto_close=True,
        )
    return metadata


def save_figure2a_preview(
    result: DirectedVascularGraph,
    image_volume: np.ndarray | None,
    output_dir: Path,
    *,
    spacing_xyz_um: tuple[float, float, float],
    max_arrows: int,
    volume_opacity: float,
    window_size: tuple[int, int],
) -> Figure2aArtifacts:
    """Write acceptance evidence without opening a GUI."""

    screenshot = output_dir / "figure2a_interactive_preview.png"
    manifest = output_dir / "figure2a_scene_manifest.json"
    geometry = scene_geometry_from_graph(
        result,
        max_arrows=max_arrows,
    )
    metadata = render_figure2a_scene(
        image_volume,
        geometry,
        spacing_xyz_um=spacing_xyz_um,
        sample_id=result.sample_id,
        volume_opacity=volume_opacity,
        window_size=window_size,
        screenshot_path=screenshot,
        show=False,
    )
    manifest.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return Figure2aArtifacts(screenshot, manifest)


def _geometry_from_exports(
    graph_dir: Path,
    max_arrows: int,
) -> InteractiveSceneGeometry:
    archive = np.load(graph_dir / "directed_branch_geometry.npz")
    points = np.asarray(archive["points_um"], dtype=float)
    vectors = np.asarray(archive["direction_parent_to_current_xyz"], dtype=float)
    offsets = np.asarray(archive["branch_offsets"], dtype=np.int64)
    branches = tuple(points[start:end] for start, end in zip(offsets[:-1], offsets[1:]))
    candidates: list[int] = []
    for start, end in zip(offsets[:-1], offsets[1:]):
        if end - start < 2:
            continue
        local = np.linspace(start, end - 1, min(4, end - start), dtype=int)
        candidates.extend(int(value) for value in np.unique(local))
    candidates = [index for index in candidates if np.linalg.norm(vectors[index]) > 0]
    if len(candidates) > max_arrows:
        keep = np.linspace(0, len(candidates) - 1, max_arrows, dtype=int)
        candidates = [candidates[index] for index in keep]

    node_points: list[tuple[float, float, float]] = []
    node_roles: list[str] = []
    with (graph_dir / "directed_nodes.csv").open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            node_points.append((float(row["x_um"]), float(row["y_um"]), float(row["z_um"])))
            node_roles.append(row["role"])
    return InteractiveSceneGeometry(
        branches_um=branches,
        arrow_points_um=points[candidates] if candidates else np.empty((0, 3)),
        arrow_vectors_xyz=vectors[candidates] if candidates else np.empty((0, 3)),
        critical_points_um=np.asarray(node_points, dtype=float).reshape((-1, 3)),
        critical_roles=tuple(node_roles),
    )


def show_saved_run(
    run_root: Path,
    *,
    sample_id: str | None = None,
    max_arrows: int = 600,
    volume_opacity: float = 0.32,
    window_size: tuple[int, int] = (1800, 900),
    sampling_run_root: Path | None = None,
    screenshot_path: Path | None = None,
    show: bool = True,
) -> Path:
    """Open the first (or requested) saved sample in a zoomable native window."""

    sample_dirs = sorted(path for path in (run_root / "samples").iterdir() if path.is_dir())
    if sample_id:
        sample_dirs = [path for path in sample_dirs if sample_id in {path.name, path.name.split("__", 1)[-1]}]
    if not sample_dirs:
        raise FileNotFoundError(f"No processed sample found in {run_root}")
    sample_root = sample_dirs[0]
    manifest_path = sample_root / "preprocess_manifest.json"
    graph_dir = sample_root / "graphs"
    if not manifest_path.is_file() or not graph_dir.is_dir():
        raise FileNotFoundError(f"Interactive artifacts are incomplete in {sample_root}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    spacing = tuple(float(value) for value in manifest["spacing_xyz_um"])
    if manifest.get("normalized_volume_path"):
        image_volume, _ = load_normalized_volume(Path(manifest["normalized_volume_path"]))
    elif manifest["record"].get("image_path"):
        image_volume = load_tiff_volume(Path(manifest["record"]["image_path"]))
    elif manifest["record"].get("mask_path"):
        image_volume = load_tiff_volume(Path(manifest["record"]["mask_path"]))
    else:
        image_volume = None
    geometry = _geometry_from_exports(
        graph_dir,
        max_arrows,
    )
    sampling_rois = tuple(
        load_sampling_display_rois(sampling_run_root)
        if sampling_run_root is not None
        else ()
    )
    render_figure2a_scene(
        image_volume,
        geometry,
        spacing_xyz_um=spacing,  # type: ignore[arg-type]
        sample_id=str(manifest["record"]["sample_id"]),
        volume_opacity=volume_opacity,
        window_size=window_size,
        screenshot_path=screenshot_path,
        show=show,
        sampling_rois=sampling_rois,
    )
    return sample_root
