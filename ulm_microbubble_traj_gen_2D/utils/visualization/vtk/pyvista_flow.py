import csv
import importlib.util
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from ...core.types import FlowField, GridDomain, RasterizedVessels
from .vtk_flow_grid import LIC_ARRAY, PRESSURE_ARRAY, SPEED_ARRAY, VELOCITY_ARRAY, build_vtk_stage_grid
from .vtk_streamlines import StreamlineContinuityError, trace_root_to_outlets, write_streamline_continuity_csv


@dataclass(frozen=True)
class CfdStageArtifacts:
    html_path: Path
    preview_path: Path
    field_vti_path: Path
    formal_streamlines_vtp_path: Path
    diagnostic_streamlines_vtp_path: Path
    continuity_csv_path: Path


def validate_cfd_flow_dependencies() -> None:
    """
    Check that the packages needed to draw and export the flow are available.
    """
    required = {
        "pyvista"       : "PyVista",
        "vtk"           : "VTK",
        "trame"         : "trame",
        "trame_vtk"     : "trame-vtk",
        "trame_vuetify" : "trame-vuetify",
    }

    missing = [display for module, display in required.items() if importlib.util.find_spec(module) is None]
    if missing:
        raise ImportError(
            "PyVista/VTK CFD visualization is missing: "
            + ", ".join(missing)
            + '. Install it with: python -m pip install "pyvista[jupyter]>=0.48"'
        )
    
    import pyvista as pv
    import vtk 
    from packaging.version import Version
    from trame_vtk.tools.vtksz2html import write_html  
    from vtkmodules.vtkRenderingLICOpenGL2 import vtkImageDataLIC2D  
    
    if Version(pv.__version__) < Version("0.48"):
        raise ImportError(f"PyVista >= 0.48 is required, found {pv.__version__}.")
    if not hasattr(pv.Plotter, "export_html") or not hasattr(pv.PolyData, "decimate_polyline"):
        raise ImportError("The installed PyVista version is missing APIs required by the CFD renderer.")


def render_cfd_flow_fields(
    domain: GridDomain,
    raster: RasterizedVessels,
    flow: FlowField,
    output_dir: Path,
    *,
    vessels,
    effective_thickness_um: float,
    boundary_depth_cells: float,
) -> tuple[Path, Path]:
    """
    Draw and save both the initial and final blood-flow views.
    """

    validate_cfd_flow_dependencies()
    from ...flow.flow_boundaries import build_flux_boundaries

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Use the same inlet and outlet openings as the blood-flow calculation.
    open_boundaries = build_flux_boundaries(
        domain,
        raster,
        vessels,
        effective_thickness_um=float(effective_thickness_um),
        depth_cells=float(boundary_depth_cells),
    )
    initial = render_cfd_flow_field(
        domain,
        raster,
        flow,
        output_dir,
        stage="initial",
        open_boundaries=open_boundaries,
        vessels=vessels,
    )
    final = render_cfd_flow_field(
        domain,
        raster,
        flow,
        output_dir,
        stage="final",
        open_boundaries=open_boundaries,
        vessels=vessels,
    )
    return initial.html_path, final.html_path


def render_cfd_flow_field(
    domain: GridDomain,
    raster: RasterizedVessels,
    flow: FlowField,
    output_dir: Path,
    *,
    stage: str,
    open_boundaries=None,
    vessels=None,
) -> CfdStageArtifacts:
    """Save one blood-flow view and the data used to check its flow paths."""

    validate_cfd_flow_dependencies()
    if stage not in {"initial", "final"}:
        raise ValueError("stage must be either 'initial' or 'final'.")
    stage_name = stage
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    artifacts = CfdStageArtifacts(
        html_path=output_dir / f"{stage_name}_flow_field_cfd.html",
        preview_path=output_dir / f"{stage_name}_flow_field_cfd.png",
        field_vti_path=output_dir / f"{stage_name}_flow_field.vti",
        formal_streamlines_vtp_path=output_dir / f"{stage_name}_streamlines.vtp",
        diagnostic_streamlines_vtp_path=(
            output_dir / f"{stage_name}_streamlines_diagnostic.vtp"
        ),
        continuity_csv_path=(
            output_dir / f"{stage_name}_root_to_outlet_continuity.csv"
        ),
    )

    stage_grid = build_vtk_stage_grid(domain, raster, flow, stage=stage_name, include_lic=True)
    # Save the full field even when no complete flow path can be drawn.
    stage_grid.image_grid.save(artifacts.field_vti_path)
    try:
        trace = trace_root_to_outlets(
            domain,
            raster,
            flow,
            stage_grid,
            open_boundaries=open_boundaries,
            vessels=vessels,
        )
    except StreamlineContinuityError as exc:
        # Keep the failed paths for inspection, but never draw them as valid.
        _write_native_artifacts(stage_grid.image_grid, exc.result, artifacts)
        raise
    except ValueError as exc:
        _write_unavailable_trace_artifacts(artifacts, stage_name, str(exc))
        raise

    _write_native_artifacts(stage_grid.image_grid, trace, artifacts)
    plotter = _new_cfd_plotter()
    try:
        _populate_cfd_plotter(plotter, domain, stage_grid, trace)
        # Create a preview image before exporting the view used in a browser.
        plotter.show(
            auto_close=False,
            interactive=False,
            screenshot=str(artifacts.preview_path),
        )
        plotter.export_html(artifacts.html_path)
        _customize_exported_html(
            artifacts.html_path,
            "Initial microvascular CFD field" if stage_name == "initial" else "Converged microvascular CFD field",
        )
    finally:
        plotter.close()
    return artifacts


def _new_cfd_plotter():
    import pyvista as pv

    return pv.Plotter(
        shape=(1, 2),
        off_screen=True,
        border=True,
        border_color="#475569",
        window_size=(2000, 900),
    )


def _populate_cfd_plotter(plotter, domain, stage_grid, trace):
    spacing = float(domain.spacing_um)
    plotter.set_background("white")

    fluid = stage_grid.fluid_grid
    outline = fluid.extract_feature_edges(
        boundary_edges=True,
        feature_edges=False,
        manifold_edges=False,
        non_manifold_edges=False,
    )
    outline = outline.copy(deep=True)
    outline.points[:, 2] = 0.04 * spacing

    render_lines = trace.render_lines.copy(deep=True)
    # Simplify only the lines used in the picture. The saved paths keep every
    # calculated point.
    decimated_lines = render_lines.decimate_polyline(
        reduction=0.60,
        maximum_error=0.20 * spacing,
    )
    if _polyline_segments_stay_in_fluid(decimated_lines, fluid, spacing_um=spacing):
        render_lines = decimated_lines
    render_lines.points[:, 2] = 0.10 * spacing
    tubes = render_lines.tube(radius=0.30 * spacing, n_sides=8, capping=True)
    glyphs = _build_sparse_direction_glyphs(render_lines, spacing_um=spacing)

    finite_speed = np.asarray(stage_grid.speed_um_s, dtype=float)
    finite_speed = finite_speed[np.isfinite(finite_speed)]
    maximum_speed = (
        float(np.percentile(finite_speed, 99.5)) if finite_speed.size else 0.0
    )
    if not np.isfinite(maximum_speed) or maximum_speed <= np.finfo(float).eps:
        maximum_speed = 1.0
    speed_clim = (0.0, maximum_speed)

    finite_pressure = np.abs(
        np.asarray(stage_grid.pressure_projection, dtype=float)
    )
    finite_pressure = finite_pressure[np.isfinite(finite_pressure)]
    pressure_extent = (
        float(np.percentile(finite_pressure, 99.5))
        if finite_pressure.size
        else 0.0
    )
    if (
        not np.isfinite(pressure_extent)
        or pressure_extent <= np.finfo(float).eps
    ):
        pressure_extent = 1.0
    pressure_clim = (-pressure_extent, pressure_extent)
    stage_title = "Initial pre-projection field" if stage_grid.stage == "initial" else "Converged projected field"

    # The left view shows speed and complete paths from the inlet to the outlets.
    plotter.subplot(0, 0)
    plotter.add_mesh(
        fluid,
        scalars=SPEED_ARRAY,
        preference="point",
        cmap="turbo",
        clim=speed_clim,
        opacity=0.90,
        lighting=False,
        show_scalar_bar=False,
        nan_opacity=0.0,
        name="speed_background",
    )
    if LIC_ARRAY in fluid.point_data:
        lic_layer = fluid.copy(deep=True)
        lic_layer.points[:, 2] = 0.025 * spacing
        plotter.add_mesh(
            lic_layer,
            scalars=LIC_ARRAY,
            preference="point",
            cmap="gray",
            clim=(0.0, 1.0),
            opacity=0.20,
            lighting=False,
            show_scalar_bar=False,
            name="lic_overlay",
        )
    plotter.add_mesh(outline, color="#334155", line_width=1.0, lighting=False, name="lumen_outline")
    plotter.add_mesh(
        tubes,
        scalars=SPEED_ARRAY,
        preference="point",
        cmap="turbo",
        clim=speed_clim,
        smooth_shading=True,
        ambient=0.30,
        diffuse=0.75,
        specular=0.28,
        specular_power=24.0,
        scalar_bar_args=_scalar_bar_args("Velocity magnitude [um/s]"),
        name="root_to_outlet_stream_tubes",
    )
    plotter.add_title(f"{stage_title}\nVelocity + LIC + continuous root-to-outlet tubes", font_size=13, color="#111827")
    _configure_planar_renderer(plotter, domain)

    # The right view shows pressure and a small number of direction arrows.
    plotter.subplot(0, 1)
    plotter.add_mesh(
        fluid,
        scalars=PRESSURE_ARRAY,
        preference="point",
        cmap="coolwarm",
        clim=pressure_clim,
        opacity=0.96,
        lighting=False,
        nan_opacity=0.0,
        scalar_bar_args=_scalar_bar_args(
            "Pressure reference" if stage_grid.stage == "initial" else "Pressure [mmHg]"
        ),
        name="pressure_background",
    )
    plotter.add_mesh(
        render_lines,
        color="white",
        line_width=1.0,
        opacity=0.38,
        lighting=False,
        render_lines_as_tubes=True,
        show_scalar_bar=False,
        name="streamline_context",
    )
    if glyphs.n_cells:
        plotter.add_mesh(
            glyphs,
            color="#111827",
            smooth_shading=True,
            ambient=0.35,
            diffuse=0.70,
            specular=0.15,
            show_scalar_bar=False,
            name="sparse_velocity_glyphs",
        )
    plotter.add_mesh(outline, color="#334155", line_width=1.0, lighting=False, name="lumen_outline")
    pressure_title = "zero pressure reference" if stage_grid.stage == "initial" else "pressure distribution"
    plotter.add_title(f"{stage_title}\n{pressure_title} + sparse velocity glyphs", font_size=13, color="#111827")
    _configure_planar_renderer(plotter, domain)

    # Keep both views at the same scale and prevent the vessel shape from stretching.
    plotter.link_views()
    plotter.subplot(0, 0)
    plotter.reset_camera()
    plotter.camera.zoom(1.03)


def _build_sparse_direction_glyphs(
    lines,
    *,
    spacing_um,
    maximum_glyphs=600,
):
    import pyvista as pv

    if lines.n_lines == 0 or VELOCITY_ARRAY not in lines.point_data:
        return pv.PolyData()
    connectivity = np.asarray(lines.lines, dtype=np.int64)
    velocity = np.asarray(lines.point_data[VELOCITY_ARRAY], dtype=float)
    selected = []
    occupied = set()
    cursor = 0
    physical_interval = 35.0 * float(spacing_um)
    bin_size = 18.0 * float(spacing_um)

    for _ in range(lines.n_lines):
        count = int(connectivity[cursor])
        ids = connectivity[cursor + 1 : cursor + 1 + count]
        cursor += count + 1
        if ids.size < 3:
            continue
        points = np.asarray(lines.points[ids, :2], dtype=float)
        cumulative = np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(points, axis=0), axis=1))))
        if cumulative[-1] <= physical_interval:
            candidates = np.asarray([ids[ids.size // 2]], dtype=int)
        else:
            targets = np.arange(0.15 * cumulative[-1], 0.92 * cumulative[-1], physical_interval)
            local = np.clip(np.searchsorted(cumulative, targets), 0, ids.size - 1)
            candidates = ids[local]
        for point_id in candidates:
            point = lines.points[int(point_id)]
            key = (int(np.floor(point[0] / bin_size)), int(np.floor(point[1] / bin_size)))
            if key in occupied:
                continue
            direction = velocity[int(point_id), :2]
            if not np.all(np.isfinite(direction)) or np.linalg.norm(direction) <= np.finfo(float).eps:
                continue
            occupied.add(key)
            selected.append(int(point_id))
            if len(selected) >= int(maximum_glyphs):
                break
        if len(selected) >= int(maximum_glyphs):
            break

    if not selected:
        return pv.PolyData()
    point_ids = np.asarray(selected, dtype=int)
    points = np.asarray(lines.points[point_ids], dtype=float).copy()
    points[:, 2] = 0.18 * float(spacing_um)
    direction = np.asarray(velocity[point_ids], dtype=float).copy()
    norm = np.linalg.norm(direction, axis=1, keepdims=True)
    direction = np.divide(direction, np.maximum(norm, np.finfo(float).eps))
    cloud = pv.PolyData(points)
    cloud.point_data["direction"] = direction
    arrow = pv.Arrow(
        tip_length=0.32,
        tip_radius=0.12,
        shaft_radius=0.035,
        tip_resolution=8,
        shaft_resolution=8,
    )
    return cloud.glyph(orient="direction", scale=False, factor=7.0 * float(spacing_um), geom=arrow)


def _polyline_segments_stay_in_fluid(lines, fluid_grid, *, spacing_um):
    """Return false when a shortened flow path cuts through a vessel wall."""

    if lines.n_lines == 0:
        return False
    connectivity = np.asarray(lines.lines, dtype=np.int64)
    cursor = 0
    samples = []
    sample_spacing = 0.25 * float(spacing_um)
    for _ in range(lines.n_lines):
        count = int(connectivity[cursor])
        ids = connectivity[cursor + 1 : cursor + 1 + count]
        cursor += count + 1
        points = np.asarray(lines.points[ids], dtype=float)
        for start, end in zip(points[:-1], points[1:]):
            length = float(np.linalg.norm(end - start))
            divisions = max(2, int(np.ceil(length / sample_spacing)) + 1)
            fraction = np.linspace(0.0, 1.0, divisions, dtype=float)[1:-1]
            if fraction.size:
                samples.append(start[None, :] + fraction[:, None] * (end - start)[None, :])
    if not samples:
        return False
    containing = np.asarray(fluid_grid.find_containing_cell(np.vstack(samples)), dtype=int)
    return bool(np.all(containing >= 0))


def _configure_planar_renderer(plotter, domain):
    plotter.view_xy()
    plotter.enable_parallel_projection()
    plotter.show_bounds(
        show_zaxis=False,
        show_zlabels=False,
        xtitle="X [um]",
        ytitle="Z [um]",
        use_2d=True,
        color="#334155",
        font_size=9,
        n_xlabels=5,
        n_ylabels=5,
        ticks="outside",
        all_edges=True,
        padding=0.01,
    )
    plotter.reset_camera()


def _write_native_artifacts(image_grid, trace, artifacts):
    image_grid.save(artifacts.field_vti_path)
    trace.formal_lines.save(artifacts.formal_streamlines_vtp_path)
    trace.diagnostic_lines.save(artifacts.diagnostic_streamlines_vtp_path)
    write_streamline_continuity_csv(artifacts.continuity_csv_path, trace)


def _write_unavailable_trace_artifacts(artifacts, stage, message):
    import pyvista as pv

    for path, note in (
        (artifacts.formal_streamlines_vtp_path, "No formal streamline was available."),
        (artifacts.diagnostic_streamlines_vtp_path, message),
    ):
        empty = pv.PolyData()
        empty.field_data["note"] = np.asarray([note])
        empty.save(path)
    with artifacts.continuity_csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["stage", "is_complete_root_to_outlet", "semantic_termination"])
        writer.writeheader()
        writer.writerow(
            {
                "stage": stage,
                "is_complete_root_to_outlet": 0,
                "semantic_termination": message,
            }
        )


def _scalar_bar_args(title):
    return {
        "title": title,
        "vertical": True,
        "position_x": 0.84,
        "position_y": 0.12,
        "width": 0.08,
        "height": 0.70,
        "title_font_size": 11,
        "label_font_size": 9,
        "color": "#111827",
        "font_family": "arial",
        "fmt": "%.3g",
    }


def _customize_exported_html(path, title):
    """Replace the default page title and add a title above the flow views."""

    html = Path(path).read_text(encoding="utf-8")
    html = html.replace("<title>VTK.js | Example - OfflineLocalView</title>", f"<title>{title}</title>", 1)
    banner = (
        '<div id="cfd-global-title" style="position:fixed;top:4px;left:50%;transform:translateX(-50%);'
        'z-index:10;pointer-events:none;padding:4px 12px;border-radius:6px;background:rgba(255,255,255,.82);'
        'color:#111827;font:600 15px Arial,sans-serif;">'
        f"{title}</div>"
    )
    html = html.replace('<div id="vtk-root"', banner + '\n    <div id="vtk-root"', 1)
    Path(path).write_text(html, encoding="utf-8")
