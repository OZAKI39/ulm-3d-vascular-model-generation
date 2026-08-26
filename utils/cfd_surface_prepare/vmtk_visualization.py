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


def save_vmtk_review_figures(
    *,
    old_custom_vtp: Path,
    raw_vtp: Path,
    final_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    output_directory: Path,
) -> tuple[Path, ...]:
    output_directory.mkdir(parents=True, exist_ok=True)
    boundary_list = list(boundaries)
    old = pv.read(old_custom_vtp).triangulate()
    raw = pv.read(raw_vtp).triangulate()
    final = pv.read(final_vtp).triangulate()
    for data in (old, raw, final):
        data.cell_data["review_part"] = _parts(data, boundary_list)

    surface_path = output_directory / "old_custom_vs_vmtk_tps_surface.png"
    plotter = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(1800, 850))
    plotter.set_background("white")
    for column, (title, data) in enumerate((("OLD CUSTOM", old), ("VMTK TPS", final))):
        plotter.subplot(0, column)
        plotter.add_text(title, color="black", font_size=14)
        _add_surface(plotter, data)
        plotter.view_isometric()
    plotter.link_views()
    plotter.show(screenshot=surface_path, auto_close=True)

    interface_path = output_directory / "old_custom_vs_vmtk_tps_interfaces.png"
    interfaces = pv.Plotter(shape=(4, 2), off_screen=True, window_size=(1800, 2400))
    interfaces.set_background("white")
    wireframe_path = output_directory / "old_custom_vs_vmtk_tps_wireframe.png"
    wires = pv.Plotter(shape=(4, 2), off_screen=True, window_size=(1800, 2400))
    wires.set_background("white")
    for row, boundary in enumerate(boundary_list):
        for column, (label, data) in enumerate((("OLD CUSTOM", old), ("VMTK TPS", raw))):
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

    final_path = output_directory / "vmtk_tps_final_surface.png"
    final_plot = pv.Plotter(off_screen=True, window_size=(1400, 1100))
    final_plot.set_background("white")
    final_plot.add_text("VMTK TPS REMESHED + SIMPLE CAPS", color="black", font_size=14)
    _add_surface(final_plot, final)
    final_plot.view_isometric()
    final_plot.show(screenshot=final_path, auto_close=True)

    raw_remeshed_path = output_directory / "vmtk_raw_vs_remeshed.png"
    comparison = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(1800, 850))
    comparison.set_background("white")
    for column, (title, data) in enumerate((("VMTK RAW", raw), ("VMTK REMESHED", final))):
        comparison.subplot(0, column)
        comparison.add_text(title, color="black", font_size=14)
        _add_surface(comparison, data, wireframe=True)
        comparison.view_isometric()
    comparison.link_views()
    comparison.show(screenshot=raw_remeshed_path, auto_close=True)
    paths = (surface_path, interface_path, wireframe_path, final_path, raw_remeshed_path)
    if not all(path.is_file() and path.stat().st_size > 0 for path in paths):
        raise RuntimeError("Required VMTK manual review figures were not created")
    return tuple(path.resolve() for path in paths)
