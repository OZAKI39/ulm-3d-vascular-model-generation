"""Required old-custom versus official-VMTK manual review figures."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pyvista as pv

from .io import BoundaryInput
from .local_cut import orthogonal_basis


def _parts(data: pv.PolyData, boundaries: list[BoundaryInput]) -> np.ndarray:
    faces = np.asarray(data.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    centers = np.asarray(data.points)[faces].mean(axis=1)
    parts = np.zeros(len(faces), dtype=np.uint8)
    for boundary in boundaries:
        relative = centers - boundary.center_um
        axial = relative @ boundary.outward_normal
        radial = np.linalg.norm(relative - np.outer(axial, boundary.outward_normal), axis=1)
        mask = (
            (axial >= -0.15 * boundary.source_radius_um)
            & (axial <= 1.08 * boundary.extension_length_um)
            & (radial <= 2.5 * boundary.source_radius_um)
        )
        parts[mask] = 1 if boundary.role == "ASSUMED_INLET" else 2
    return parts


def _camera(plotter: pv.Plotter, boundary: BoundaryInput) -> None:
    first, second = orthogonal_basis(boundary.outward_normal)
    focal = boundary.center_um + 1.3 * boundary.source_radius_um * boundary.outward_normal
    distance = max(16.0 * boundary.source_radius_um, boundary.extension_length_um)
    plotter.camera.focal_point = focal
    plotter.camera.position = focal + first * distance + 0.3 * boundary.outward_normal * distance
    plotter.camera.up = second
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = 3.8 * boundary.source_radius_um


def _add_surface(plotter: pv.Plotter, data: pv.PolyData, *, wireframe: bool = False) -> None:
    if "review_part" not in data.cell_data:
        raise ValueError("review_part cell array is required")
    if wireframe:
        plotter.add_mesh(
            data,
            scalars="review_part",
            preference="cell",
            categories=True,
            cmap=["#dce6eb", "#2878b5", "#c51f1f"],
            show_scalar_bar=False,
            show_edges=True,
            edge_color="#13232d",
            line_width=0.8,
        )
    else:
        plotter.add_mesh(
            data,
            scalars="review_part",
            preference="cell",
            categories=True,
            cmap=["#b8c4ca", "#2878b5", "#c51f1f"],
            show_scalar_bar=False,
            smooth_shading=True,
        )


def _local(data: pv.PolyData, boundary: BoundaryInput) -> pv.PolyData:
    faces = np.asarray(data.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    centers = np.asarray(data.points)[faces].mean(axis=1)
    radius = max(4.5 * boundary.source_radius_um, 0.4 * boundary.extension_length_um)
    mask = np.linalg.norm(centers - boundary.center_um, axis=1) <= radius
    return data.extract_cells(np.flatnonzero(mask)).extract_surface().triangulate()


def _extension_local(data: pv.PolyData, boundary: BoundaryInput) -> pv.PolyData:
    faces = np.asarray(data.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    centers = np.asarray(data.points)[faces].mean(axis=1)
    relative = centers - boundary.center_um
    axial = relative @ boundary.outward_normal
    radial = np.linalg.norm(
        relative - np.outer(axial, boundary.outward_normal), axis=1
    )
    mask = (
        (axial >= -2.0 * boundary.source_radius_um)
        & (axial <= 1.25 * boundary.extension_length_um)
        & (radial <= 7.0 * boundary.source_radius_um)
    )
    return data.extract_cells(np.flatnonzero(mask)).extract_surface().triangulate()


def save_vmtk_review_figures(
    *,
    old_custom_vtp: Path,
    centerline_raw_vtp: Path,
    boundarynormal_raw_vtp: Path,
    boundarynormal_remeshed_open_vtp: Path,
    boundarynormal_final_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    output_directory: Path,
) -> tuple[Path, ...]:
    """Create the five mandatory manual-review comparisons with fixed cameras."""

    output_directory.mkdir(parents=True, exist_ok=True)
    boundary_list = list(boundaries)
    old = pv.read(old_custom_vtp).triangulate()
    centerline = pv.read(centerline_raw_vtp).triangulate()
    raw = pv.read(boundarynormal_raw_vtp).triangulate()
    remeshed = pv.read(boundarynormal_remeshed_open_vtp).triangulate()
    final = pv.read(boundarynormal_final_vtp).triangulate()
    for data in (old, centerline, raw, remeshed, final):
        data.cell_data["review_part"] = _parts(data, boundary_list)

    surface_path = (
        output_directory / "old_custom_vs_centerline_vs_boundarynormal_surface.png"
    )
    plotter = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(2400, 800))
    plotter.set_background("white")
    surfaces = (
        ("OLD CUSTOM REFINED", old),
        ("VMTK CENTERLINE RAW", centerline),
        ("VMTK BOUNDARYNORMAL FINAL", final),
    )
    for column, (title, data) in enumerate(surfaces):
        plotter.subplot(0, column)
        plotter.add_text(title, color="black", font_size=14)
        _add_surface(plotter, data)
        plotter.view_isometric()
    plotter.link_views()
    plotter.show(screenshot=surface_path, auto_close=True)

    interface_path = output_directory / "three_way_interface_closeups.png"
    interfaces = pv.Plotter(shape=(4, 3), off_screen=True, window_size=(2400, 2400))
    interfaces.set_background("white")
    wireframe_path = output_directory / "three_way_interface_wireframe.png"
    wires = pv.Plotter(shape=(4, 3), off_screen=True, window_size=(2400, 2400))
    wires.set_background("white")
    for row, boundary in enumerate(boundary_list):
        comparisons = (
            ("OLD CUSTOM", old),
            ("VMTK CENTERLINE", centerline),
            ("VMTK BOUNDARYNORMAL", raw),
        )
        for column, (label, data) in enumerate(comparisons):
            local = _local(data, boundary)
            short = boundary.port_id.rsplit("__", 1)[-1]
            interfaces.subplot(row, column)
            interfaces.add_text(f"{short} | {label}", color="black", font_size=11)
            _add_surface(interfaces, local)
            _camera(interfaces, boundary)
            wires.subplot(row, column)
            wires.add_text(f"{short} | {label}", color="black", font_size=11)
            _add_surface(wires, local, wireframe=True)
            _camera(wires, boundary)
    interfaces.show(screenshot=interface_path, auto_close=True)
    wires.show(screenshot=wireframe_path, auto_close=True)

    raw_remeshed_path = output_directory / "boundarynormal_raw_vs_remeshed.png"
    comparison = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(1800, 850))
    comparison.set_background("white")
    raw_remeshed = (
        ("VMTK BOUNDARYNORMAL RAW", raw),
        ("VMTK BOUNDARYNORMAL REMESHED", remeshed),
    )
    for column, (title, data) in enumerate(raw_remeshed):
        comparison.subplot(0, column)
        comparison.add_text(title, color="black", font_size=14)
        _add_surface(comparison, data, wireframe=True)
        comparison.view_isometric()
    comparison.link_views()
    comparison.show(screenshot=raw_remeshed_path, auto_close=True)

    cut002 = next(
        boundary for boundary in boundary_list if boundary.port_id.endswith("__cut_002")
    )
    cut002_path = output_directory / "cut002_centerline_vs_boundarynormal.png"
    cut_plot = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(1800, 900))
    cut_plot.set_background("white")
    for column, (title, data) in enumerate(
        (
            ("CUT_002 CENTERLINE_DIRECTION", centerline),
            ("CUT_002 BOUNDARY_NORMAL", raw),
        )
    ):
        cut_plot.subplot(0, column)
        cut_plot.add_text(title, color="black", font_size=14)
        _add_surface(cut_plot, _extension_local(data, cut002), wireframe=True)
        _camera(cut_plot, cut002)
    cut_plot.show(screenshot=cut002_path, auto_close=True)

    paths = (
        surface_path,
        interface_path,
        wireframe_path,
        raw_remeshed_path,
        cut002_path,
    )
    if not all(path.is_file() and path.stat().st_size > 0 for path in paths):
        raise RuntimeError("Required VMTK manual review figures were not created")
    return tuple(path.resolve() for path in paths)


def _hotspot_local(data: pv.PolyData, center: np.ndarray, radius: float) -> pv.PolyData:
    faces = np.asarray(data.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    centers = np.asarray(data.points)[faces].mean(axis=1)
    mask = np.linalg.norm(centers - center, axis=1) <= radius
    return data.extract_cells(np.flatnonzero(mask)).extract_surface().triangulate()


def _hotspot_camera(plotter: pv.Plotter, center: np.ndarray, radius: float) -> None:
    view = np.asarray((1.0, 0.8, 0.65), dtype=float)
    view /= np.linalg.norm(view)
    plotter.camera.focal_point = center
    plotter.camera.position = center + 8.0 * radius * view
    plotter.camera.up = (0.0, 0.0, 1.0)
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = 1.15 * radius


def _add_plain_surface(
    plotter: pv.Plotter,
    data: pv.PolyData,
    *,
    wireframe: bool = False,
    color: str = "#b8c4ca",
) -> None:
    plotter.add_mesh(
        data,
        color=color,
        show_edges=wireframe,
        edge_color="#13232d",
        line_width=0.7,
        smooth_shading=not wireframe,
    )


def save_caponly_review_figures(
    *,
    original_vtp: Path,
    raw_vtp: Path,
    previous_global_final_vtp: Path,
    caponly_final_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    hotspots: list[dict[str, float]],
    output_directory: Path,
) -> tuple[Path, ...]:
    """Create quantitative, automatically centered cap-only review figures."""

    output_directory.mkdir(parents=True, exist_ok=True)
    boundary_list = list(boundaries)
    original = pv.read(original_vtp).triangulate()
    raw = pv.read(raw_vtp).triangulate()
    previous = pv.read(previous_global_final_vtp).triangulate()
    caponly = pv.read(caponly_final_vtp).triangulate()
    for data in (raw, previous, caponly):
        data.cell_data["review_part"] = _parts(data, boundary_list)

    whole_path = output_directory / "global_remesh_vs_caponly_whole_surface.png"
    whole = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(2000, 900))
    whole.set_background("white")
    for column, (title, data) in enumerate(
        (("PREVIOUS GLOBAL-REMESHED FINAL", previous), ("NEW CAP-ONLY FINAL", caponly))
    ):
        whole.subplot(0, column)
        whole.add_text(title, color="black", font_size=14)
        _add_surface(whole, data)
        whole.view_isometric()
    whole.link_views()
    whole.show(screenshot=whole_path, auto_close=True)

    selected = hotspots[: max(5, min(8, len(hotspots)))]
    if len(selected) < 5:
        raise RuntimeError("At least five automatic core hotspots are required")
    radius = max(3.0 * float(np.median([item.source_radius_um for item in boundary_list])), 2.0)
    bifurcation_path = output_directory / "bifurcation_artifact_comparison.png"
    bifurcation = pv.Plotter(
        shape=(len(selected), 4), off_screen=True, window_size=(3000, 650 * len(selected))
    )
    bifurcation.set_background("white")
    three_path = output_directory / "previous_global_remesh_artifact_hotspots.png"
    three = pv.Plotter(
        shape=(len(selected), 3), off_screen=True, window_size=(2400, 650 * len(selected))
    )
    three.set_background("white")
    for row, hotspot in enumerate(selected):
        center = np.asarray(
            (hotspot["center_x_um"], hotspot["center_y_um"], hotspot["center_z_um"]),
            dtype=float,
        )
        four_columns = (
            ("ORIGINAL ULTRALISER", original),
            ("BOUNDARYNORMAL RAW", raw),
            ("CAP-ONLY FINAL", caponly),
            ("PREVIOUS GLOBAL-REMESHED", previous),
        )
        for column, (label, data) in enumerate(four_columns):
            bifurcation.subplot(row, column)
            bifurcation.add_text(
                f"H{hotspot['hotspot_id']} | {label}", color="black", font_size=10
            )
            _add_plain_surface(bifurcation, _hotspot_local(data, center, radius))
            _hotspot_camera(bifurcation, center, radius)
        for column, (label, data) in enumerate(
            (
                ("ORIGINAL ULTRALISER", original),
                ("BOUNDARYNORMAL RAW", raw),
                ("PREVIOUS GLOBAL-REMESHED", previous),
            )
        ):
            three.subplot(row, column)
            three.add_text(
                f"H{hotspot['hotspot_id']} | {label}", color="black", font_size=10
            )
            _add_plain_surface(three, _hotspot_local(data, center, radius))
            _hotspot_camera(three, center, radius)
    bifurcation.show(screenshot=bifurcation_path, auto_close=True)
    three.show(screenshot=three_path, auto_close=True)

    closeup_path = output_directory / "extension_caponly_closeups.png"
    closeups = pv.Plotter(shape=(4, 2), off_screen=True, window_size=(1800, 2400))
    closeups.set_background("white")
    wireframe_path = output_directory / "extension_caponly_wireframe.png"
    wires = pv.Plotter(shape=(4, 2), off_screen=True, window_size=(1800, 2400))
    wires.set_background("white")
    for row, boundary in enumerate(boundary_list):
        for column, (label, data) in enumerate((("RAW", raw), ("CAP-ONLY", caponly))):
            local_data = _extension_local(data, boundary)
            short = boundary.port_id.rsplit("__", 1)[-1]
            closeups.subplot(row, column)
            closeups.add_text(f"{short} | {label}", color="black", font_size=11)
            _add_surface(closeups, local_data)
            _camera(closeups, boundary)
            wires.subplot(row, column)
            wires.add_text(f"{short} | {label}", color="black", font_size=11)
            _add_surface(wires, local_data, wireframe=True)
            _camera(wires, boundary)
    closeups.show(screenshot=closeup_path, auto_close=True)
    wires.show(screenshot=wireframe_path, auto_close=True)

    normal_path = output_directory / "core_normal_deviation_hotspots.png"
    normals = pv.Plotter(
        shape=(len(selected), 3), off_screen=True, window_size=(2400, 650 * len(selected))
    )
    normals.set_background("white")
    for row, hotspot in enumerate(selected):
        center = np.asarray(
            (hotspot["center_x_um"], hotspot["center_y_um"], hotspot["center_z_um"]),
            dtype=float,
        )
        for column, (label, data, color) in enumerate(
            (
                ("ORIGINAL NORMALS", original, "#a8b7bf"),
                ("GLOBAL-REMESH DEVIATION", previous, "#c94f45"),
                ("CAP-ONLY DEVIATION", caponly, "#2c7fb8"),
            )
        ):
            normals.subplot(row, column)
            normals.add_text(
                f"H{hotspot['hotspot_id']} | {label}", color="black", font_size=10
            )
            _add_plain_surface(
                normals,
                _hotspot_local(data, center, radius),
                wireframe=True,
                color=color,
            )
            _hotspot_camera(normals, center, radius)
    normals.show(screenshot=normal_path, auto_close=True)

    paths = (whole_path, bifurcation_path, closeup_path, wireframe_path, normal_path)
    required = (*paths, three_path)
    if not all(path.is_file() and path.stat().st_size > 0 for path in required):
        raise RuntimeError("Required cap-only manual review figures were not created")
    return tuple(path.resolve() for path in paths)
