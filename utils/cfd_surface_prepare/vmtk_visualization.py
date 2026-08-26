"""Required old-custom versus official-VMTK manual review figures."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pyvista as pv

from .io import BoundaryInput
from .local_cut import orthogonal_basis
from .mesh_quality import triangle_metrics


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


def save_entityremesh_review_figures(
    *,
    original_vtp: Path,
    raw_vtp: Path,
    previous_global_final_vtp: Path,
    entity_remeshed_open_vtp: Path,
    entity_final_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    output_directory: Path,
) -> tuple[Path, ...]:
    """Create the five required entity-aware remesh manual-review figures."""

    output_directory.mkdir(parents=True, exist_ok=True)
    boundary_list = list(boundaries)
    original = pv.read(original_vtp).triangulate()
    raw = pv.read(raw_vtp).triangulate()
    previous = pv.read(previous_global_final_vtp).triangulate()
    remeshed = pv.read(entity_remeshed_open_vtp).triangulate()
    final = pv.read(entity_final_vtp).triangulate()
    for data in (raw, previous, remeshed, final):
        data.cell_data["review_part"] = _parts(data, boundary_list)

    closeup_path = output_directory / "extension_three_way_closeups.png"
    closeups = pv.Plotter(shape=(4, 3), off_screen=True, window_size=(2400, 2400))
    closeups.set_background("white")
    wire_path = output_directory / "extension_three_way_wireframe.png"
    wires = pv.Plotter(shape=(4, 3), off_screen=True, window_size=(2400, 2400))
    wires.set_background("white")
    comparisons = (
        ("CAP-ONLY RAW / NO REMESH", raw),
        ("PREVIOUS GLOBAL REMESH", previous),
        ("NEW ENTITY-AWARE REMESH", final),
    )
    for row, boundary in enumerate(boundary_list):
        for column, (label, data) in enumerate(comparisons):
            local = _extension_local(data, boundary)
            short = boundary.port_id.rsplit("__", 1)[-1]
            closeups.subplot(row, column)
            closeups.add_text(f"{short} | {label}", color="black", font_size=10)
            _add_surface(closeups, local)
            _camera(closeups, boundary)
            wires.subplot(row, column)
            wires.add_text(f"{short} | {label}", color="black", font_size=10)
            _add_surface(wires, local, wireframe=True)
            _camera(wires, boundary)
    closeups.show(screenshot=closeup_path, auto_close=True)
    wires.show(screenshot=wire_path, auto_close=True)

    core_path = output_directory / "core_original_vs_global_vs_entityremesh.png"
    core = pv.Plotter(shape=(1, 3), off_screen=True, window_size=(2700, 900))
    core.set_background("white")
    for column, (label, data) in enumerate(
        (
            ("ORIGINAL ULTRALISER", original),
            ("PREVIOUS GLOBAL REMESH", previous),
            ("NEW ENTITY-AWARE REMESH", final),
        )
    ):
        core.subplot(0, column)
        core.add_text(label, color="black", font_size=13)
        _add_plain_surface(core, data, wireframe=True)
        core.view_isometric()
    core.link_views()
    core.show(screenshot=core_path, auto_close=True)

    interface_path = output_directory / "entity_interface_closeups.png"
    interfaces = pv.Plotter(shape=(4, 2), off_screen=True, window_size=(1800, 2400))
    interfaces.set_background("white")
    for row, boundary in enumerate(boundary_list):
        for column, (label, data) in enumerate(
            (("RAW ENTITY INTERFACE", raw), ("REMESHED ENTITY INTERFACE", remeshed))
        ):
            local = _local(data, boundary)
            short = boundary.port_id.rsplit("__", 1)[-1]
            interfaces.subplot(row, column)
            interfaces.add_text(f"{short} | {label}", color="black", font_size=10)
            _add_surface(interfaces, local, wireframe=True)
            _camera(interfaces, boundary)
    interfaces.show(screenshot=interface_path, auto_close=True)

    tail_path = output_directory / "extension_mesh_tail_hotspots.png"
    tails = pv.Plotter(shape=(4, 2), off_screen=True, window_size=(1800, 2400))
    tails.set_background("white")
    remesh_faces = np.asarray(remeshed.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    regions = np.asarray(remeshed.cell_data["SurfaceRegionId"], dtype=np.uint8)
    extension_ids = np.flatnonzero(regions == 1)
    extension_centers = np.asarray(remeshed.points)[remesh_faces[extension_ids]].mean(axis=1)
    scores = []
    for boundary in boundary_list:
        relative = extension_centers - boundary.center_um
        axial = relative @ boundary.outward_normal
        radial = np.linalg.norm(
            relative - np.outer(axial, boundary.outward_normal), axis=1
        )
        before = np.maximum(-axial, 0.0)
        after = np.maximum(axial - boundary.extension_length_um, 0.0)
        scores.append(
            (radial / boundary.source_radius_um) ** 2
            + ((before + after) / boundary.source_radius_um) ** 2
        )
    assignment = np.argmin(np.column_stack(scores), axis=1)
    for row, boundary in enumerate(boundary_list):
        selected = extension_ids[assignment == row]
        metrics = triangle_metrics(np.asarray(remeshed.points), remesh_faces[selected])
        worst = (
            ("LOWEST MINIMUM ANGLE", selected[int(np.argmin(metrics.minimum_angles_deg))]),
            ("HIGHEST ASPECT RATIO", selected[int(np.argmax(metrics.aspect_ratios))]),
        )
        for column, (label, face_id) in enumerate(worst):
            center = np.asarray(remeshed.points)[remesh_faces[face_id]].mean(axis=0)
            tails.subplot(row, column)
            tails.add_text(
                f"{boundary.port_id.rsplit('__', 1)[-1]} | {label}",
                color="black",
                font_size=10,
            )
            _add_surface(tails, _extension_local(remeshed, boundary), wireframe=True)
            tails.add_mesh(
                pv.Sphere(
                    radius=max(0.12 * boundary.source_radius_um, 0.03),
                    center=center,
                ),
                color="#ffcc00",
            )
            _camera(tails, boundary)
    tails.show(screenshot=tail_path, auto_close=True)

    paths = (closeup_path, wire_path, core_path, interface_path, tail_path)
    if not all(path.is_file() and path.stat().st_size > 0 for path in paths):
        raise RuntimeError("Required entity-remesh manual review figures were not created")
    return tuple(path.resolve() for path in paths)


def _add_guarded_surface(
    plotter: pv.Plotter, data: pv.PolyData, *, wireframe: bool = True
) -> None:
    plotter.add_mesh(
        data,
        scalars="RemeshEntityId",
        preference="cell",
        categories=True,
        cmap=["#8da0ae", "#f2a93b", "#2878b5"],
        clim=(1, 3),
        show_scalar_bar=False,
        show_edges=wireframe,
        edge_color="#17242b",
        line_width=0.7,
        smooth_shading=not wireframe,
    )


def save_previous_collision_closeup(
    *,
    remeshed_vtp: Path,
    diagnosis: dict[str, object],
    boundary: BoundaryInput,
    output_path: Path,
) -> Path:
    """Show the exact saved CORE/BODY penetration with both triangles highlighted."""

    output_path.parent.mkdir(parents=True, exist_ok=True)
    data = pv.read(remeshed_vtp).triangulate()
    faces = np.asarray(data.faces, dtype=np.int64).reshape((-1, 4))[:, 1:]
    core_id = int(diagnosis["core_face_id"])
    extension_id = int(diagnosis["extension_face_id"])
    center = np.asarray(data.points)[faces[[core_id, extension_id]]].mean(axis=(0, 1))
    face_centers = np.asarray(data.points)[faces].mean(axis=1)
    local_ids = np.flatnonzero(
        np.linalg.norm(face_centers - center, axis=1)
        <= 2.5 * boundary.source_radius_um
    )
    local = data.extract_cells(local_ids).extract_surface().triangulate()
    plotter = pv.Plotter(off_screen=True, window_size=(1500, 1200))
    plotter.set_background("white")
    plotter.add_text(
        "PREVIOUS ENTITY-REMESH TRUE CORE/EXTENSION PENETRATION | cut_002",
        color="black",
        font_size=13,
    )
    plotter.add_mesh(
        local,
        scalars="RemeshEntityId",
        preference="cell",
        categories=True,
        cmap=["#8da0ae", "#2878b5"],
        show_scalar_bar=False,
        show_edges=True,
        edge_color="#1b252a",
    )
    plotter.add_mesh(
        data.extract_cells([core_id]).extract_surface(),
        color="#ffe119",
        show_edges=True,
        edge_color="#111111",
        line_width=5,
    )
    plotter.add_mesh(
        data.extract_cells([extension_id]).extract_surface(),
        color="#f032e6",
        show_edges=True,
        edge_color="#111111",
        line_width=5,
    )
    start = diagnosis.get("intersection_segment_start_xyz")
    end = diagnosis.get("intersection_segment_end_xyz")
    if start is not None and end is not None:
        plotter.add_mesh(pv.Line(start, end), color="#e6194b", line_width=10)
    _camera(plotter, boundary)
    plotter.camera.focal_point = center
    plotter.camera.parallel_scale = 2.2 * boundary.source_radius_um
    plotter.show(screenshot=output_path, auto_close=True)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError("Previous collision close-up was not created")
    return output_path.resolve()


def save_guarded_open_figures(
    *,
    raw_vtp: Path,
    remeshed_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    tail_records: list[dict[str, object]],
    intersections: list[dict[str, object]],
    output_directory: Path,
) -> tuple[Path, ...]:
    """Create guarded-region, interface, wireframe, tail, and collision evidence."""

    output_directory.mkdir(parents=True, exist_ok=True)
    boundary_list = list(boundaries)
    raw = pv.read(raw_vtp).triangulate()
    remeshed = pv.read(remeshed_vtp).triangulate()
    guard_path = output_directory / "guard_region_visualization.png"
    guard_plot = pv.Plotter(shape=(4, 1), off_screen=True, window_size=(1000, 2400))
    guard_plot.set_background("white")
    for row, boundary in enumerate(boundary_list):
        guard_plot.subplot(row, 0)
        guard_plot.add_text(
            f"{boundary.port_id.rsplit('__', 1)[-1]} | CORE / GUARD / BODY",
            color="black",
            font_size=11,
        )
        _add_guarded_surface(guard_plot, _extension_local(raw, boundary))
        _camera(guard_plot, boundary)
    guard_plot.show(screenshot=guard_path, auto_close=True)

    interface_path = output_directory / "guarded_entityremesh_interface_closeups.png"
    wire_path = output_directory / "guarded_entityremesh_wireframe.png"
    interfaces = pv.Plotter(shape=(4, 2), off_screen=True, window_size=(1800, 2400))
    wires = pv.Plotter(shape=(4, 2), off_screen=True, window_size=(1800, 2400))
    interfaces.set_background("white")
    wires.set_background("white")
    for row, boundary in enumerate(boundary_list):
        for column, (label, data) in enumerate(
            (("RAW GUARDED ENTITIES", raw), ("GUARDED REMESH", remeshed))
        ):
            local = _local(data, boundary)
            title = f"{boundary.port_id.rsplit('__', 1)[-1]} | {label}"
            interfaces.subplot(row, column)
            interfaces.add_text(title, color="black", font_size=10)
            _add_guarded_surface(interfaces, local, wireframe=False)
            _camera(interfaces, boundary)
            wires.subplot(row, column)
            wires.add_text(title, color="black", font_size=10)
            _add_guarded_surface(wires, local, wireframe=True)
            _camera(wires, boundary)
    interfaces.show(screenshot=interface_path, auto_close=True)
    wires.show(screenshot=wire_path, auto_close=True)

    tail_path = output_directory / "mesh_tail_hotspots.png"
    tails = pv.Plotter(shape=(4, 2), off_screen=True, window_size=(1800, 2400))
    tails.set_background("white")
    remeshed_faces = np.asarray(remeshed.faces).reshape((-1, 4))[:, 1:]
    for row, boundary in enumerate(boundary_list):
        port_records = [
            record for record in tail_records if record["port_id"] == boundary.port_id
        ]
        for column, entity in enumerate(("PROXIMAL_GUARD", "EXTENSION_BODY")):
            records = [record for record in port_records if record["entity"] == entity]
            tails.subplot(row, column)
            tails.add_text(
                f"{boundary.port_id.rsplit('__', 1)[-1]} | {entity}",
                color="black",
                font_size=10,
            )
            _add_guarded_surface(tails, _extension_local(remeshed, boundary))
            if records:
                worst = max(records, key=lambda item: float(item["aspect_ratio"]))
                face_id = int(worst["triangle_id"])
                center = np.asarray(remeshed.points)[remeshed_faces[face_id]].mean(axis=0)
                tails.add_mesh(
                    pv.Sphere(0.12 * boundary.source_radius_um, center=center),
                    color="#e6194b",
                )
            _camera(tails, boundary)
    tails.show(screenshot=tail_path, auto_close=True)

    collision_path = output_directory / "guarded_entityremesh_collision_closeups.png"
    collision = pv.Plotter(off_screen=True, window_size=(1500, 1200))
    collision.set_background("white")
    if intersections:
        first = intersections[0]
        first_id = int(first["first_face_id"])
        second_id = int(first["second_face_id"])
        center = np.asarray(remeshed.points)[
            remeshed_faces[[first_id, second_id]]
        ].mean(axis=(0, 1))
        boundary = boundary_list[
            int(
                np.argmin(
                    [np.linalg.norm(center - item.center_um) for item in boundary_list]
                )
            )
        ]
        collision.add_text("TRUE GUARDED-REMESH INTERSECTION", color="black")
        _add_guarded_surface(collision, _local(remeshed, boundary))
        for face_id, color in ((first_id, "#ffe119"), (second_id, "#f032e6")):
            collision.add_mesh(
                remeshed.extract_cells([face_id]).extract_surface(),
                color=color,
                show_edges=True,
                edge_color="black",
                line_width=5,
            )
        _camera(collision, boundary)
        collision.camera.focal_point = center
    else:
        collision.add_text(
            "NO TRUE SELF-INTERSECTIONS DETECTED",
            color="#137333",
            font_size=18,
        )
        _add_guarded_surface(collision, remeshed)
        collision.view_isometric()
    collision.show(screenshot=collision_path, auto_close=True)
    paths = (guard_path, interface_path, wire_path, tail_path, collision_path)
    if not all(path.is_file() and path.stat().st_size > 0 for path in paths):
        raise RuntimeError("Required guarded-open figures were not created")
    return tuple(path.resolve() for path in paths)


def save_guarded_three_way_comparison(
    *,
    caponly_vtp: Path,
    previous_entityremesh_vtp: Path,
    guarded_final_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    output_path: Path,
) -> Path:
    """Compare cap-only, previous entity-remesh, and new guarded final surfaces."""

    boundary_list = list(boundaries)
    surfaces = (
        ("CAP-ONLY", pv.read(caponly_vtp).triangulate()),
        ("PREVIOUS ENTITY REMESH", pv.read(previous_entityremesh_vtp).triangulate()),
        ("NEW GUARDED ENTITY REMESH", pv.read(guarded_final_vtp).triangulate()),
    )
    plotter = pv.Plotter(shape=(4, 3), off_screen=True, window_size=(2400, 2400))
    plotter.set_background("white")
    for row, boundary in enumerate(boundary_list):
        for column, (label, data) in enumerate(surfaces):
            plotter.subplot(row, column)
            plotter.add_text(
                f"{boundary.port_id.rsplit('__', 1)[-1]} | {label}",
                color="black",
                font_size=9,
            )
            _add_plain_surface(
                plotter, _extension_local(data, boundary), wireframe=True
            )
            _camera(plotter, boundary)
    plotter.show(screenshot=output_path, auto_close=True)
    if not output_path.is_file() or output_path.stat().st_size == 0:
        raise RuntimeError("Three-way guarded comparison was not created")
    return output_path.resolve()


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


def _add_crossseam_assignment(
    plotter: pv.Plotter, data: pv.PolyData, *, wireframe: bool = True
) -> None:
    display = data.copy()
    entities = np.asarray(display.cell_data["RemeshEntityId"], dtype=np.int32)
    regions = np.asarray(display.cell_data["SurfaceRegionId"], dtype=np.uint8)
    display.cell_data["CrossSeamDisplayRegion"] = np.where(
        entities == 1, 0, np.where(regions == 0, 1, 2)
    ).astype(np.uint8)
    plotter.add_mesh(
        display,
        scalars="CrossSeamDisplayRegion",
        preference="cell",
        categories=True,
        cmap=["#9aa9b2", "#f2a93b", "#2878b5"],
        clim=(0, 2),
        show_scalar_bar=False,
        show_edges=wireframe,
        edge_color="#17242b",
        line_width=0.7,
        smooth_shading=not wireframe,
    )


def save_crossseam_review_figures(
    *,
    raw_vtp: Path,
    global_remesh_vtp: Path,
    guarded_remesh_vtp: Path,
    crossseam_final_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    output_directory: Path,
) -> tuple[Path, ...]:
    """Create exactly four cross-seam manual-review figures."""

    output_directory.mkdir(parents=True, exist_ok=True)
    boundary_list = list(boundaries)
    raw = pv.read(raw_vtp).triangulate()
    global_remesh = pv.read(global_remesh_vtp).triangulate()
    guarded = pv.read(guarded_remesh_vtp).triangulate()
    crossseam = pv.read(crossseam_final_vtp).triangulate()
    seam_lines = _original_cut_seam_lines(raw)

    active_path = output_directory / "active_region_visualization.png"
    active_plot = pv.Plotter(
        shape=(4, 1), off_screen=True, window_size=(1100, 2400)
    )
    active_plot.set_background("white")
    for row, boundary in enumerate(boundary_list):
        active_plot.subplot(row, 0)
        active_plot.add_text(
            f"{boundary.port_id.rsplit('__', 1)[-1]} | FAR_CORE / ACTIVE CORE COLLAR / EXTENSION",
            color="black",
            font_size=10,
        )
        _add_crossseam_assignment(
            active_plot, _extension_local(raw, boundary), wireframe=True
        )
        active_plot.add_mesh(seam_lines, color="#e6194b", line_width=6)
        _camera(active_plot, boundary)
    active_plot.show(screenshot=active_path, auto_close=True)

    three_way_path = output_directory / "seam_wireframe_three_way.png"
    three_way = pv.Plotter(
        shape=(4, 3), off_screen=True, window_size=(2400, 2400)
    )
    three_way.set_background("white")
    surfaces = (
        ("GLOBAL REMESH", global_remesh),
        ("GUARDED ENTITY REMESH", guarded),
        ("NEW CROSS-SEAM REMESH", crossseam),
    )
    for row, boundary in enumerate(boundary_list):
        for column, (label, data) in enumerate(surfaces):
            three_way.subplot(row, column)
            three_way.add_text(
                f"{boundary.port_id.rsplit('__', 1)[-1]} | {label}",
                color="black",
                font_size=9,
            )
            _add_plain_surface(
                three_way, _extension_local(data, boundary), wireframe=True
            )
            _camera(three_way, boundary)
    three_way.show(screenshot=three_way_path, auto_close=True)

    closeup_path = output_directory / "crossseam_interface_closeups.png"
    closeups = pv.Plotter(
        shape=(4, 1), off_screen=True, window_size=(1100, 2400)
    )
    closeups.set_background("white")
    for row, boundary in enumerate(boundary_list):
        closeups.subplot(row, 0)
        closeups.add_text(
            f"{boundary.port_id.rsplit('__', 1)[-1]} | CROSS-SEAM REMESH | red = original seam",
            color="black",
            font_size=10,
        )
        _add_plain_surface(
            closeups, _local(crossseam, boundary), wireframe=True
        )
        closeups.add_mesh(seam_lines, color="#e6194b", line_width=6)
        _camera(closeups, boundary)
        closeups.camera.parallel_scale = 2.8 * boundary.source_radius_um
    closeups.show(screenshot=closeup_path, auto_close=True)

    final_path = output_directory / "final_surface_review.png"
    final_plot = pv.Plotter(off_screen=True, window_size=(1600, 1200))
    final_plot.set_background("white")
    final_plot.add_text(
        "FINAL CROSS-SEAM CFD SURFACE | VISUAL ACCEPTANCE REQUIRES MANUAL REVIEW",
        color="black",
        font_size=13,
    )
    _add_plain_surface(final_plot, crossseam, wireframe=False)
    final_plot.view_isometric()
    final_plot.enable_parallel_projection()
    final_plot.show(screenshot=final_path, auto_close=True)

    paths = (active_path, three_way_path, closeup_path, final_path)
    if not all(path.is_file() and path.stat().st_size > 0 for path in paths):
        raise RuntimeError("Required cross-seam figures were not created")
    return tuple(path.resolve() for path in paths)


def save_crossseam_resume_figures(
    *,
    raw_vtp: Path,
    global_remesh_vtp: Path,
    guarded_remesh_vtp: Path,
    frozen_open_vtp: Path,
    final_vtp: Path,
    boundaries: Iterable[BoundaryInput],
    hotspots: list[dict[str, object]],
    output_directory: Path,
) -> tuple[Path, ...]:
    """Create exactly four figures for frozen-open post-QC continuation."""

    output_directory.mkdir(parents=True, exist_ok=True)
    boundary_list = list(boundaries)
    global_remesh = pv.read(global_remesh_vtp).triangulate()
    guarded = pv.read(guarded_remesh_vtp).triangulate()
    frozen = pv.read(frozen_open_vtp).triangulate()
    final = pv.read(final_vtp).triangulate()
    seam_lines = _original_cut_seam_lines(pv.read(raw_vtp).triangulate())

    three_way_path = output_directory / "seam_wireframe_three_way.png"
    three_way = pv.Plotter(
        shape=(4, 3), off_screen=True, window_size=(2400, 2400)
    )
    three_way.set_background("white")
    surfaces = (
        ("GLOBAL REMESH", global_remesh),
        ("GUARDED ENTITY REMESH", guarded),
        ("FROZEN CROSS-SEAM REMESH", frozen),
    )
    for row, boundary in enumerate(boundary_list):
        for column, (label, data) in enumerate(surfaces):
            three_way.subplot(row, column)
            three_way.add_text(
                f"{boundary.port_id.rsplit('__', 1)[-1]} | {label}",
                color="black",
                font_size=9,
            )
            _add_plain_surface(
                three_way, _extension_local(data, boundary), wireframe=True
            )
            _camera(three_way, boundary)
    three_way.show(screenshot=three_way_path, auto_close=True)

    closeup_path = output_directory / "crossseam_interface_closeups.png"
    closeups = pv.Plotter(
        shape=(4, 1), off_screen=True, window_size=(1100, 2400)
    )
    closeups.set_background("white")
    for row, boundary in enumerate(boundary_list):
        closeups.subplot(row, 0)
        closeups.add_text(
            f"{boundary.port_id.rsplit('__', 1)[-1]} | frozen OPEN | red = original seam",
            color="black",
            font_size=10,
        )
        _add_plain_surface(closeups, _local(frozen, boundary), wireframe=True)
        closeups.add_mesh(seam_lines, color="#e6194b", line_width=6)
        _camera(closeups, boundary)
        closeups.camera.parallel_scale = 2.8 * boundary.source_radius_um
    closeups.show(screenshot=closeup_path, auto_close=True)

    hotspot_path = output_directory / "active_collar_distance_hotspots.png"
    displayed = hotspots[:10]
    rows = 2
    columns = 5
    hotspot_plot = pv.Plotter(
        shape=(rows, columns), off_screen=True, window_size=(2500, 1000)
    )
    hotspot_plot.set_background("white")
    for index in range(rows * columns):
        hotspot_plot.subplot(index // columns, index % columns)
        if index >= len(displayed):
            hotspot_plot.add_text("No hotspot", color="black", font_size=10)
            continue
        record = displayed[index]
        sample = np.asarray(
            (
                record["sample_x_um"],
                record["sample_y_um"],
                record["sample_z_um"],
            ),
            dtype=float,
        )
        closest = np.asarray(
            (
                record["closest_x_um"],
                record["closest_y_um"],
                record["closest_z_um"],
            ),
            dtype=float,
        )
        boundary = next(
            item for item in boundary_list if item.port_id == record["port_id"]
        )
        radius = max(1.5 * boundary.source_radius_um, 0.5)
        _add_plain_surface(
            hotspot_plot, _hotspot_local(frozen, sample, radius), wireframe=True
        )
        hotspot_plot.add_mesh(
            pv.Sphere(max(0.035 * boundary.source_radius_um, 0.015), center=sample),
            color="#e6194b",
        )
        hotspot_plot.add_mesh(
            pv.Sphere(
                max(0.035 * boundary.source_radius_um, 0.015), center=closest
            ),
            color="#3cb44b",
        )
        line = pv.Line(sample, closest)
        hotspot_plot.add_mesh(line, color="#ffe119", line_width=5)
        hotspot_plot.add_text(
            f"#{record['rank']} | {float(record['distance_um']):.4g} um",
            color="black",
            font_size=9,
        )
        _hotspot_camera(hotspot_plot, sample, radius)
    hotspot_plot.show(screenshot=hotspot_path, auto_close=True)

    final_path = output_directory / "final_surface_review.png"
    final_plot = pv.Plotter(off_screen=True, window_size=(1600, 1200))
    final_plot.set_background("white")
    final_plot.add_text(
        "FINAL CAPPED CROSS-SEAM CFD SURFACE | MANUAL REVIEW REQUIRED",
        color="black",
        font_size=13,
    )
    _add_plain_surface(final_plot, final, wireframe=False)
    final_plot.view_isometric()
    final_plot.enable_parallel_projection()
    final_plot.show(screenshot=final_path, auto_close=True)

    paths = (three_way_path, closeup_path, hotspot_path, final_path)
    if not all(path.is_file() and path.stat().st_size > 0 for path in paths):
        raise RuntimeError("Required frozen cross-seam figures were not created")
    return tuple(path.resolve() for path in paths)
