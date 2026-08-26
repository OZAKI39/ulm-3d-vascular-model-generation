"""Before/after and wireframe evidence for manual extension-mesh review."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pyvista as pv

from .io import BoundaryInput
from .local_cut import orthogonal_basis
from .types import TaggedSurface


COLORS = ["#a9b6bf", "#2878b5", "#c51f1f"]


def _polydata(
    vertices: np.ndarray, faces: np.ndarray, parts: np.ndarray | None = None
) -> pv.PolyData:
    vtk_faces = np.column_stack(
        (np.full(len(faces), 3, dtype=np.int64), faces)
    ).ravel()
    data = pv.PolyData(vertices, vtk_faces)
    if parts is not None:
        data.cell_data["part"] = parts
    return data


def _camera(plotter: pv.Plotter, boundary: BoundaryInput) -> None:
    first, second = orthogonal_basis(boundary.outward_normal)
    focal = (
        boundary.center_um
        + 1.5 * boundary.source_radius_um * boundary.outward_normal
    )
    distance = max(
        18.0 * boundary.source_radius_um, 1.5 * boundary.extension_length_um
    )
    plotter.camera.focal_point = focal
    plotter.camera.position = (
        focal + first * distance + 0.35 * boundary.outward_normal * distance
    )
    plotter.camera.up = second
    plotter.enable_parallel_projection()
    plotter.camera.parallel_scale = 3.5 * boundary.source_radius_um


def _previous_side_mask(
    centers: np.ndarray, boundary: BoundaryInput, boundary_type: np.ndarray
) -> np.ndarray:
    relative = centers - boundary.center_um
    axial = relative @ boundary.outward_normal
    radial = np.linalg.norm(
        relative - np.outer(axial, boundary.outward_normal), axis=1
    )
    return (
        (boundary_type == 0)
        & (axial > 0.0)
        & (axial < boundary.extension_length_um + 1.0e-8)
        & (radial < 2.0 * boundary.source_radius_um)
    )


def save_refinement_review_figures(
    previous: TaggedSurface,
    refined: TaggedSurface,
    boundaries: Iterable[BoundaryInput],
    output_directory: Path,
) -> tuple[Path, Path, Path]:
    output_directory.mkdir(parents=True, exist_ok=True)
    boundary_list = list(boundaries)
    previous_data = _polydata(
        previous.vertices, previous.faces, previous.face_kind
    )
    refined_data = _polydata(refined.vertices, refined.faces, refined.face_kind)

    whole = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(1800, 850))
    whole.set_background("white")
    pairs = (
        ("Previous direct extrusion", previous_data),
        ("Refined multi-ring extension", refined_data),
    )
    for column, (title, data) in enumerate(pairs):
        whole.subplot(0, column)
        whole.add_text(title, color="black", font_size=13)
        whole.add_mesh(
            data,
            scalars="part",
            preference="cell",
            categories=True,
            cmap=COLORS,
            show_scalar_bar=False,
        )
        whole.view_isometric()
    whole.link_views()
    whole_path = output_directory / "previous_vs_refined_surface.png"
    whole.show(screenshot=whole_path, auto_close=True)

    previous_centers = previous.vertices[previous.faces].mean(axis=1)
    refined_centers = refined.vertices[refined.faces].mean(axis=1)
    closeups = pv.Plotter(shape=(4, 2), off_screen=True, window_size=(1800, 2400))
    closeups.set_background("white")
    wireframe = pv.Plotter(shape=(4, 2), off_screen=True, window_size=(1800, 2400))
    wireframe.set_background("white")
    for row, boundary in enumerate(boundary_list):
        radius = max(
            4.5 * boundary.source_radius_um, 0.45 * boundary.extension_length_um
        )
        old_mask = (
            np.linalg.norm(previous_centers - boundary.center_um, axis=1) <= radius
        )
        new_mask = (
            np.linalg.norm(refined_centers - boundary.center_um, axis=1) <= radius
        )
        old_local = _polydata(
            previous.vertices,
            previous.faces[old_mask],
            previous.face_kind[old_mask],
        )
        new_local = _polydata(
            refined.vertices,
            refined.faces[new_mask],
            refined.face_kind[new_mask],
        )
        short_id = boundary.port_id.rsplit("__", 1)[-1]
        for column, (label, data) in enumerate(
            (("BEFORE", old_local), ("AFTER", new_local))
        ):
            closeups.subplot(row, column)
            closeups.add_text(
                f"{short_id} | {label}", color="black", font_size=11
            )
            closeups.add_mesh(
                data,
                scalars="part",
                preference="cell",
                categories=True,
                cmap=COLORS,
                show_scalar_bar=False,
            )
            _camera(closeups, boundary)
            wireframe.subplot(row, column)
            wireframe.add_text(
                f"{short_id} | {label} WIREFRAME", color="black", font_size=11
            )
            wireframe.add_mesh(
                data,
                style="surface",
                color="#dce6eb",
                opacity=0.35,
                show_edges=True,
                edge_color="#15232d",
                line_width=1.0,
            )
            _camera(wireframe, boundary)
    closeup_path = output_directory / "extension_interface_before_after.png"
    closeups.show(screenshot=closeup_path, auto_close=True)
    wireframe_path = output_directory / "extension_wireframe_before_after.png"
    wireframe.show(screenshot=wireframe_path, auto_close=True)
    return whole_path.resolve(), closeup_path.resolve(), wireframe_path.resolve()


def load_previous_tagged_surface(
    path: Path, boundaries: Iterable[BoundaryInput]
) -> TaggedSurface:
    """Load the immutable old VTP and classify only its direct extensions."""

    data = pv.read(path).triangulate()
    vertices = np.asarray(data.points, dtype=float)
    faces = np.asarray(data.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    boundary_type = np.asarray(
        data.cell_data["boundary_type_code"], dtype=np.uint8
    )
    boundary_index = np.asarray(data.cell_data["boundary_index"], dtype=np.int32)
    boundary_origin = np.asarray(
        data.cell_data["boundary_origin_code"], dtype=np.uint8
    )
    centers = vertices[faces].mean(axis=1)
    face_kind = np.zeros(len(faces), dtype=np.uint8)
    extension_index = np.full(len(faces), -1, dtype=np.int32)
    for boundary in boundaries:
        side = _previous_side_mask(centers, boundary, boundary_type)
        face_kind[side] = 1
        extension_index[side] = boundary.index
    face_kind[boundary_type > 0] = 2
    return TaggedSurface(
        vertices=vertices,
        faces=faces,
        boundary_type=boundary_type,
        boundary_index=boundary_index,
        boundary_origin=boundary_origin,
        face_kind=face_kind,
        extension_index=extension_index,
        extension_band=np.full(len(faces), -1, dtype=np.int32),
        source_vertex_index=np.full(len(vertices), -1, dtype=np.int64),
    )
