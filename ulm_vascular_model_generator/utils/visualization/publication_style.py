"""Publication-oriented geometry and scene refinements for vascular trees."""

import math

import numpy as np


def add_smooth_junctions(
    pv,
    vessel_mesh,
    tree,
    radius_scale,
    sides,
):
    """
    Cover internal cylinder seams with smooth, radius-matched junctions.

    SWC coordinates and radii remain unchanged. The spheres affect only the
    display mesh and carry the original physical radius as scalar data.
    """
    degree = {node_id: 0 for node_id in tree.nodes}
    for parent_id, child_id in tree.edges:
        degree[parent_id] += 1
        degree[child_id] += 1

    junctions = []
    theta_resolution = max(int(sides), 24)
    phi_resolution = max(theta_resolution // 2, 16)

    for node_id, node in tree.nodes.items():
        if degree[node_id] < 2:
            continue

        junction = pv.Sphere(
            radius=node.radius_um * radius_scale * 0.98,
            center=node.xyz,
            theta_resolution=theta_resolution,
            phi_resolution=phi_resolution,
        )
        junction.cell_data["radius_um"] = np.full(
            junction.n_cells,
            node.radius_um,
        )
        junction.cell_data["edge_id"] = np.full(
            junction.n_cells,
            -1,
            dtype=np.int64,
        )
        junctions.append(junction)

    if not junctions:
        return vessel_mesh

    combined = pv.merge(
        [vessel_mesh, *junctions],
        merge_points=False,
    )
    try:
        combined = combined.compute_normals(
            cell_normals=False,
            point_normals=True,
            split_vertices=False,
            auto_orient_normals=True,
        )
    except Exception:
        pass
    return combined


def _nice_scale_length(span):
    target = max(0.18 * span, 1e-12)
    magnitude = 10.0 ** math.floor(math.log10(target))
    for multiplier in (5.0, 2.0, 1.0):
        candidate = multiplier * magnitude
        if candidate <= target:
            return candidate
    return magnitude


def add_scale_bar(plotter, pv, mesh, view, font_size):
    """Add a data-space scale bar for the orthographic publication views."""
    bounds = mesh.bounds
    span = max(
        bounds[1] - bounds[0],
        bounds[3] - bounds[2],
        bounds[5] - bounds[4],
        1.0,
    )
    length = _nice_scale_length(span)
    margin = 0.055 * span
    tick = 0.012 * span

    if view == "xz":
        depth = bounds[2] - 0.02 * span
        start = np.array([bounds[0] + margin, depth, bounds[4] + margin])
        end = start + np.array([length, 0.0, 0.0])
        tick_axis = np.array([0.0, 0.0, tick])
        label_point = 0.5 * (start + end) + 1.8 * tick_axis
    elif view == "xy":
        depth = bounds[5] + 0.02 * span
        start = np.array([bounds[0] + margin, bounds[2] + margin, depth])
        end = start + np.array([length, 0.0, 0.0])
        tick_axis = np.array([0.0, tick, 0.0])
        label_point = 0.5 * (start + end) + 1.8 * tick_axis
    elif view == "yz":
        depth = bounds[0] - 0.02 * span
        start = np.array([depth, bounds[2] + margin, bounds[4] + margin])
        end = start + np.array([0.0, length, 0.0])
        tick_axis = np.array([0.0, 0.0, tick])
        label_point = 0.5 * (start + end) + 1.8 * tick_axis
    else:
        return

    lines = pv.merge(
        [
            pv.Line(start, end),
            pv.Line(start - tick_axis, start + tick_axis),
            pv.Line(end - tick_axis, end + tick_axis),
        ],
        merge_points=False,
    )
    plotter.add_mesh(
        lines,
        color="#20262d",
        line_width=4,
        lighting=False,
        render_lines_as_tubes=True,
    )
    plotter.add_point_labels(
        [label_point],
        [f"{length:g} \u00b5m"],
        show_points=False,
        shape=None,
        bold=False,
        font_size=max(int(round(0.72 * font_size)), 12),
        text_color="#20262d",
        always_visible=True,
    )


def _add_planar_axes(
    plotter,
    pv,
    horizontal_range,
    vertical_range,
    make_point,
    horizontal_title,
    vertical_title,
    font_size,
):
    """Draw uncluttered data-space axes for one orthographic plane."""
    horizontal_min, horizontal_max = horizontal_range
    vertical_min, vertical_max = vertical_range
    span = max(
        horizontal_max - horizontal_min,
        vertical_max - vertical_min,
        1.0,
    )
    axis_offset = 0.035 * span
    tick_length = 0.012 * span
    label_offset = 0.045 * span
    title_offset = 0.085 * span

    horizontal_axis_value = vertical_min - axis_offset
    vertical_axis_value = horizontal_min - axis_offset
    horizontal_ticks = np.linspace(horizontal_min, horizontal_max, 5)
    vertical_ticks = np.linspace(vertical_min, vertical_max, 5)

    segments = [
        pv.Line(
            make_point(horizontal_min, horizontal_axis_value),
            make_point(horizontal_max, horizontal_axis_value),
        ),
        pv.Line(
            make_point(vertical_axis_value, vertical_min),
            make_point(vertical_axis_value, vertical_max),
        ),
    ]
    for value in horizontal_ticks:
        segments.append(
            pv.Line(
                make_point(value, horizontal_axis_value - tick_length),
                make_point(value, horizontal_axis_value + tick_length),
            )
        )
    for value in vertical_ticks:
        segments.append(
            pv.Line(
                make_point(vertical_axis_value - tick_length, value),
                make_point(vertical_axis_value + tick_length, value),
            )
        )

    axes_mesh = pv.merge(segments, merge_points=False)
    plotter.add_mesh(
        axes_mesh,
        color="#20262d",
        line_width=3,
        lighting=False,
        render_lines_as_tubes=True,
    )

    tick_font_size = max(int(round(0.72 * font_size)), 14)
    plotter.add_point_labels(
        [
            make_point(value, horizontal_axis_value - label_offset)
            for value in horizontal_ticks
        ],
        [f"{value:.0f}" for value in horizontal_ticks],
        show_points=False,
        shape=None,
        bold=False,
        font_size=tick_font_size,
        text_color="#20262d",
        always_visible=True,
    )
    plotter.add_point_labels(
        [
            make_point(vertical_axis_value - label_offset, value)
            for value in vertical_ticks
        ],
        [f"{value:.0f}" for value in vertical_ticks],
        show_points=False,
        shape=None,
        bold=False,
        font_size=tick_font_size,
        text_color="#20262d",
        always_visible=True,
    )
    plotter.add_point_labels(
        [
            make_point(
                0.5 * (horizontal_min + horizontal_max),
                horizontal_axis_value - title_offset,
            ),
            make_point(
                vertical_axis_value - label_offset,
                vertical_max + 0.55 * title_offset,
            ),
        ],
        [horizontal_title, vertical_title],
        show_points=False,
        shape=None,
        bold=False,
        font_size=font_size,
        text_color="#20262d",
        always_visible=True,
    )


def add_publication_axes(plotter, pv, mesh, view, font_size):
    """Draw quantitative axes with micrometre units for the selected view."""
    bounds = mesh.bounds
    span = max(
        bounds[1] - bounds[0],
        bounds[3] - bounds[2],
        bounds[5] - bounds[4],
        1.0,
    )

    if view == "xz":
        depth = bounds[2] - 0.02 * span
        _add_planar_axes(
            plotter,
            pv,
            (bounds[0], bounds[1]),
            (bounds[4], bounds[5]),
            lambda horizontal, vertical: np.array(
                [horizontal, depth, vertical]
            ),
            "X (\u00b5m)",
            "Z (\u00b5m)",
            font_size,
        )
    elif view == "xy":
        depth = bounds[5] + 0.02 * span
        _add_planar_axes(
            plotter,
            pv,
            (bounds[0], bounds[1]),
            (bounds[2], bounds[3]),
            lambda horizontal, vertical: np.array(
                [horizontal, vertical, depth]
            ),
            "X (\u00b5m)",
            "Y (\u00b5m)",
            font_size,
        )
    elif view == "yz":
        depth = bounds[0] - 0.02 * span
        _add_planar_axes(
            plotter,
            pv,
            (bounds[2], bounds[3]),
            (bounds[4], bounds[5]),
            lambda horizontal, vertical: np.array(
                [depth, horizontal, vertical]
            ),
            "Y (\u00b5m)",
            "Z (\u00b5m)",
            font_size,
        )
    else:
        plotter.show_bounds(
            mesh=mesh,
            xtitle="X (\u00b5m)",
            ytitle="Y (\u00b5m)",
            ztitle="Z (\u00b5m)",
            n_xlabels=4,
            n_ylabels=4,
            n_zlabels=4,
            font_size=font_size,
            font_family="arial",
            color="#20262d",
            fmt="%.0f",
            location="outer",
        )


def add_publication_effects(
    plotter,
    pv,
    mesh,
    tree,
    view,
    font_size,
    show_scale_bar,
):
    """Apply restrained depth cues and quantitative annotations."""
    try:
        plotter.add_silhouette(
            mesh,
            color="#26323b",
            line_width=1.1,
            opacity=0.24,
            feature_angle=30.0,
        )
    except Exception:
        pass

    try:
        median_radius = float(np.median(tree.radii_um))
        plotter.enable_ssao(
            radius=max(2.5 * median_radius, 1.0),
            bias=max(0.03 * median_radius, 0.01),
            kernel_size=256,
            blur=True,
        )
    except Exception:
        pass

    if show_scale_bar:
        add_scale_bar(plotter, pv, mesh, view, font_size)
