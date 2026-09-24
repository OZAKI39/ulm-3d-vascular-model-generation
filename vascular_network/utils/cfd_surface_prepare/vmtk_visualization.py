"""The two manual-review figures emitted by the formal surface pipeline."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pyvista as pv

from .io import BoundaryInput


def _camera(plotter: pv.Plotter, boundary: BoundaryInput) -> None:
    normal = np.asarray(boundary.outward_normal, dtype=float)
    normal /= np.linalg.norm(normal)
    up = np.asarray((0.0, 0.0, 1.0))
    if abs(float(np.dot(normal, up))) > 0.9:
        up = np.asarray((0.0, 1.0, 0.0))
    center = np.asarray(boundary.center_um, dtype=float)
    plotter.camera.focal_point = center
    plotter.camera.position = center + 8.0 * boundary.source_radius_um * normal
    plotter.camera.up = up
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = 2.8 * boundary.source_radius_um


def _local(data: pv.PolyData, boundary: BoundaryInput) -> pv.PolyData:
    center = np.asarray(boundary.center_um, dtype=float)
    radius = 4.0 * boundary.source_radius_um
    points = np.asarray(data.points, dtype=float)
    mask = np.linalg.norm(points - center, axis=1) <= radius
    return data.extract_points(mask, adjacent_cells=True).extract_surface()


def _add_surface(
    plotter: pv.Plotter, data: pv.PolyData, *, wireframe: bool = False
) -> None:
    plotter.add_mesh(
        data,
        color="#b8c4ca",
        show_edges=wireframe,
        edge_color="#13232d",
        line_width=0.7,
        smooth_shading=not wireframe,
    )


def _original_cut_seam_lines(raw: pv.PolyData) -> pv.PolyData:
    faces = np.asarray(raw.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    regions = np.asarray(raw.cell_data["SurfaceRegionId"], dtype=np.uint8)
    edge_faces: dict[tuple[int, int], list[int]] = {}
    for face_id, face in enumerate(faces):
        for first, second in zip(face, np.roll(face, -1)):
            edge_faces.setdefault(
                tuple(sorted((int(first), int(second)))), []
            ).append(face_id)
    seam_edges = [
        edge
        for edge, linked in edge_faces.items()
        if len(linked) == 2
        and {int(regions[linked[0]]), int(regions[linked[1]])} == {0, 1}
    ]
    if not seam_edges:
        return pv.PolyData()
    points = np.asarray(raw.points)[np.asarray(seam_edges, dtype=np.int64)].reshape(
        (-1, 3)
    )
    lines = np.column_stack(
        (
            np.full(len(seam_edges), 2, dtype=np.int64),
            np.arange(0, 2 * len(seam_edges), 2, dtype=np.int64),
            np.arange(1, 2 * len(seam_edges), 2, dtype=np.int64),
        )
    ).ravel()
    return pv.PolyData(points, lines=lines)


def save_production_review_figures(
    *,
    raw_vtp: Path,
    remeshed_open_vtp: Path,
    final_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    output_directory: Path,
) -> tuple[Path, Path]:
    """Write the interface close-ups and the final capped-surface review."""

    output_directory.mkdir(parents=True, exist_ok=True)
    boundary_list = list(boundaries)
    raw = pv.read(raw_vtp).triangulate()
    remeshed = pv.read(remeshed_open_vtp).triangulate()
    final = pv.read(final_vtp).triangulate()
    seam_lines = _original_cut_seam_lines(raw)

    closeup_path = output_directory / "crossseam_interface_closeups.png"
    closeups = pv.Plotter(shape=(len(boundary_list), 1), off_screen=True, window_size=(1100, 2400))
    closeups.set_background("white")
    for row, boundary in enumerate(boundary_list):
        closeups.subplot(row, 0)
        closeups.add_text(
            f"{boundary.port_id.rsplit('__', 1)[-1]} | cross-seam OPEN | red = original seam",
            color="black",
            font_size=10,
        )
        _add_surface(closeups, _local(remeshed, boundary), wireframe=True)
        closeups.add_mesh(seam_lines, color="#e6194b", line_width=6)
        _camera(closeups, boundary)
    closeups.show(screenshot=closeup_path, auto_close=True)

    final_path = output_directory / "final_surface_review.png"
    final_plot = pv.Plotter(off_screen=True, window_size=(1600, 1200))
    final_plot.set_background("white")
    final_plot.add_text(
        "FINAL CAPPED CROSS-SEAM CFD SURFACE | MANUAL REVIEW REQUIRED",
        color="black",
        font_size=13,
    )
    _add_surface(final_plot, final)
    final_plot.view_isometric()
    final_plot.enable_parallel_projection()
    final_plot.show(screenshot=final_path, auto_close=True)

    paths = (closeup_path, final_path)
    if not all(path.is_file() and path.stat().st_size > 0 for path in paths):
        raise RuntimeError("Required production review figures were not created")
    return tuple(path.resolve() for path in paths)
