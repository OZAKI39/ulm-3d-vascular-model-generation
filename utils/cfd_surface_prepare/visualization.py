"""Off-screen visual evidence for the required manual surface review."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pyvista as pv
import trimesh

from .io import BoundaryInput
from .local_cut import orthogonal_basis
from .types import TaggedSurface


def _polydata(vertices: np.ndarray, faces: np.ndarray) -> pv.PolyData:
    vtk_faces = np.column_stack(
        (np.full(len(faces), 3, dtype=np.int64), faces)
    ).ravel()
    return pv.PolyData(vertices, vtk_faces)


def save_review_figures(
    original: trimesh.Trimesh,
    final: TaggedSurface,
    boundaries: Iterable[BoundaryInput],
    output_directory: Path,
) -> tuple[Path, Path]:
    """Render one whole-model comparison and one four-boundary close-up sheet."""

    output_directory.mkdir(parents=True, exist_ok=True)
    original_data = _polydata(
        np.asarray(original.vertices), np.asarray(original.faces)
    )
    final_data = _polydata(final.vertices, final.faces)
    final_data.cell_data["part"] = final.face_kind
    comparison = pv.Plotter(shape=(1, 2), off_screen=True, window_size=(1800, 850))
    comparison.set_background("white")
    comparison.subplot(0, 0)
    comparison.add_text("Original validated Ultraliser surface", font_size=12, color="black")
    comparison.add_mesh(original_data, color="#b8c4cf", smooth_shading=True)
    comparison.view_isometric()
    comparison.subplot(0, 1)
    comparison.add_text("Derived CFD surface: 4 straight extensions", font_size=12, color="black")
    comparison.add_mesh(
        final_data,
        scalars="part",
        preference="cell",
        categories=True,
        cmap=["#b8c4cf", "#2f77b4", "#d62728"],
        show_scalar_bar=False,
        smooth_shading=False,
    )
    comparison.view_isometric()
    comparison.link_views()
    comparison_path = output_directory / "original_vs_cfd_surface.png"
    comparison.show(screenshot=comparison_path, auto_close=True)

    closeups = pv.Plotter(shape=(2, 2), off_screen=True, window_size=(1800, 1500))
    closeups.set_background("white")
    face_centers_original = np.asarray(original.triangles_center)
    face_centers_final = np.mean(final.vertices[final.faces], axis=1)
    for boundary in boundaries:
        row, column = divmod(boundary.index, 2)
        closeups.subplot(row, column)
        radius = max(4.0 * boundary.source_radius_um, 1.25 * boundary.extension_length_um)
        original_mask = np.linalg.norm(
            face_centers_original - boundary.center_um, axis=1
        ) <= radius
        final_mask = np.linalg.norm(
            face_centers_final - boundary.center_um, axis=1
        ) <= radius * 1.5
        if np.any(original_mask):
            closeups.add_mesh(
                _polydata(
                    np.asarray(original.vertices),
                    np.asarray(original.faces)[original_mask],
                ),
                color="#b0b0b0",
                opacity=0.35,
            )
        local_final = _polydata(final.vertices, final.faces[final_mask])
        local_final.cell_data["part"] = final.face_kind[final_mask]
        closeups.add_mesh(
            local_final,
            scalars="part",
            preference="cell",
            categories=True,
            cmap=["#7aa6c2", "#2878b5", "#d62728"],
            show_scalar_bar=False,
        )
        label_role = "INLET" if boundary.role == "ASSUMED_INLET" else "OUTLET"
        closeups.add_text(
            f"{label_role} | {boundary.boundary_origin}\nL={boundary.extension_length_um:.3f} um",
            font_size=11,
            color="black",
        )
        closeups.add_axes()
        first, second = orthogonal_basis(boundary.outward_normal)
        focal_point = (
            boundary.center_um
            + 0.5 * boundary.extension_length_um * boundary.outward_normal
        )
        camera_distance = max(
            4.0 * boundary.extension_length_um,
            16.0 * boundary.source_radius_um,
        )
        closeups.camera.focal_point = focal_point
        closeups.camera.position = (
            focal_point
            + first * camera_distance
            + boundary.outward_normal * 0.45 * camera_distance
        )
        closeups.camera.up = second
        closeups.enable_parallel_projection()
        closeups.camera.parallel_scale = max(
            0.75 * boundary.extension_length_um,
            4.0 * boundary.source_radius_um,
        )
    closeup_path = output_directory / "boundary_closeups.png"
    closeups.show(screenshot=closeup_path, auto_close=True)
    return comparison_path.resolve(), closeup_path.resolve()
