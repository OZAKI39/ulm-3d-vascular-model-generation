"""Matplotlib rendering for microbubble flow in field-based result folders."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from scipy import ndimage

from .bubble_style import (
    BubbleGlowStyle,
    bubble_highlight_geometry,
    build_bubble_glow_style,
    glow_layer_diameters_points,
)
from .molecular_target_overlay import target_display_masks
from .result_loader import FieldVisualizationData

ColorMode = Literal[
    "speed",
    "wall_shear",
    "local_shear",
    "vessel_id",
    "diameter",
    "clearance",
    "active",
]
FlowStage = Literal["initial", "final"]

_NEAR_WALL_RING_RGBA = np.asarray([1.0, 0.65, 0.0, 1.0])
_TARGET_ADHESION_RING_RGBA = np.asarray([1.0, 0.1, 0.65, 1.0])
_INVALID_GAP_RING_RGBA = np.asarray([1.0, 0.0, 0.0, 1.0])


@dataclass
class _BubbleArtists:
    """Keep the real disk and its purely visual glow layers synchronized in animations."""

    core: object
    halos: tuple[object, ...]
    contact_ring: object | None
    highlight: object | None
    halo_points_per_um: float

    @property
    def all(self) -> tuple[object, ...]:
        artists = (*self.halos, self.core)
        if self.contact_ring is not None:
            artists = (*artists, self.contact_ring)
        return artists if self.highlight is None else (*artists, self.highlight)


def render_snapshot(
    data: FieldVisualizationData,
    output_path: Path,
    *,
    frame_index: int = -1,
    max_bubbles: int = 0,
    tail_length: int = 25,
    color_mode: ColorMode = "speed",
    dpi: int = 180,
    glow_enabled: bool = True,
    glow_scale: float = 3.0,
) -> Path:
    """Save a single X-Z frame showing vessel lumen, bubble positions, and tails."""

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.collections import EllipseCollection, LineCollection

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    frames = data.positions_um.shape[0]
    frame = _normalize_frame_index(frame_index, frames)
    selected = _select_bubbles(data, max_bubbles)
    frame_selected = _active_frame_selection(data, selected, frame)

    fig, ax = plt.subplots(figsize=(10.5, 8.0), constrained_layout=True)
    color_limits = _color_limits(data, selected, color_mode)
    _draw_background(
        ax,
        data,
        color_mode=color_mode,
        color_limits=color_limits,
    )
    _finish_axes(ax, data, _particle_frame_title(data, frame, frames, displayed=frame_selected.size))
    fig.canvas.draw()
    _draw_tails(ax, data, frame_selected, frame, tail_length, LineCollection)
    bubble_artists, label = _draw_bubbles(
        ax,
        data,
        frame_selected,
        frame,
        color_mode,
        EllipseCollection,
        color_limit_selection=selected,
        glow_style=build_bubble_glow_style(enabled=glow_enabled, outer_diameter_scale=glow_scale),
        color_limits=color_limits,
    )
    fig.colorbar(bubble_artists.core, ax=ax, label=label, fraction=0.046, pad=0.03)
    _draw_bubble_state_legend(ax)
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return output_path


def render_animation(
    data: FieldVisualizationData,
    output_path: Path,
    *,
    max_bubbles: int = 0,
    tail_length: int = 20,
    stride: int = 3,
    fps: int = 20,
    color_mode: ColorMode = "speed",
    dpi: int = 130,
    show_progress: bool = True,
    glow_enabled: bool = True,
    glow_scale: float = 3.0,
) -> Path:
    """Save an animated GIF or MP4 of bubbles moving through the vessel field."""

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.collections import EllipseCollection, LineCollection

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    selected = _select_bubbles(data, max_bubbles)
    frame_indices = np.arange(0, data.positions_um.shape[0], max(1, int(stride)), dtype=int)
    if frame_indices[-1] != data.positions_um.shape[0] - 1:
        frame_indices = np.append(frame_indices, data.positions_um.shape[0] - 1)

    fig, ax = plt.subplots(
        figsize=(10.5, 8.0),
        dpi=max(1, int(dpi)),
        constrained_layout=True,
    )
    color_limits = _color_limits(data, selected, color_mode)
    _draw_background(
        ax,
        data,
        color_mode=color_mode,
        color_limits=color_limits,
    )
    _finish_axes(ax, data, "")
    fig.canvas.draw()
    tail_collection = LineCollection([], colors="#ffffff", linewidths=0.85, alpha=0.55, zorder=2.5)
    ax.add_collection(tail_collection)
    first_frame = int(frame_indices[0])
    first_selected = _active_frame_selection(data, selected, first_frame)
    glow_style = build_bubble_glow_style(enabled=glow_enabled, outer_diameter_scale=glow_scale)
    bubble_artists, label = _draw_bubbles(
        ax,
        data,
        first_selected,
        first_frame,
        color_mode,
        EllipseCollection,
        color_limit_selection=selected,
        glow_style=glow_style,
        color_limits=color_limits,
    )
    title = ax.set_title(
        _particle_frame_title(
            data,
            first_frame,
            data.positions_um.shape[0],
            displayed=first_selected.size,
        ),
        fontsize=11,
    )
    fig.colorbar(bubble_artists.core, ax=ax, label=label, fraction=0.046, pad=0.03)
    _draw_bubble_state_legend(ax)

    def update(frame: int):
        frame_selected = _active_frame_selection(data, selected, frame)
        segments = _tail_segments(data, frame_selected, frame, tail_length)
        tail_collection.set_segments(segments)
        xy = _bubble_xy(data, frame_selected, frame)
        values = _bubble_color_values(data, frame_selected, frame, color_mode)
        diameters = _bubble_physical_diameters(data, frame_selected, frame)
        _update_bubble_artists(
            bubble_artists,
            xy,
            values,
            diameters,
            _bubble_edge_colors(data, frame_selected, frame),
            _bubble_contact_ring_colors(data, frame_selected, frame),
            _bubble_rotation_angles(data, frame_selected, frame),
            glow_style,
        )
        title.set_text(
            _particle_frame_title(data, frame, data.positions_um.shape[0], displayed=frame_selected.size)
        )
        return tail_collection, *bubble_artists.all, title

    dynamic_artists = (tail_collection, *bubble_artists.all, title)
    for artist in dynamic_artists:
        artist.set_animated(True)
    fig.canvas.draw()
    fig.set_layout_engine("none")
    background = fig.canvas.copy_from_bbox(fig.bbox)

    try:
        _save_blitted_animation(
            fig,
            frame_indices,
            update,
            dynamic_artists,
            background,
            output_path,
            fps=max(1, int(fps)),
            show_progress=show_progress,
        )
    finally:
        plt.close(fig)
    return output_path


def render_hole_map(
    data: FieldVisualizationData,
    output_path: Path,
    *,
    dpi: int = 240,
) -> Path:
    """Save a high-contrast static map of filled holes in the lumen mask."""

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    extent = _image_extent(data)
    fig, ax = plt.subplots(figsize=(12.5, 9.5), constrained_layout=True)
    _draw_binary_lumen_context(ax, data, extent)
    hole_count = int(np.count_nonzero(data.hole_mask))
    if hole_count > 0:
        glow = ndimage.binary_dilation(data.hole_mask, iterations=_visibility_iterations(data))
        _imshow_boolean(ax, glow, extent, "#ff1744", alpha=0.32)
        _imshow_boolean(ax, data.hole_mask, extent, "#ffff00", alpha=1.0)
        ax.contour(
            data.x_coordinates_um,
            data.z_coordinates_um,
            data.hole_mask.T.astype(float),
            levels=[0.5],
            colors="#ffffff",
            linewidths=1.2,
            alpha=1.0,
        )
        _mark_components(ax, data.hole_mask, data, prefix="H")
    _finish_quality_axes(ax, data, f"Lumen holes | cells={hole_count}")
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return output_path


def render_narrow_cell_map(
    data: FieldVisualizationData,
    output_path: Path,
    *,
    dpi: int = 240,
) -> Path:
    """Save a high-contrast static map of narrow lumen cells."""

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    extent = _image_extent(data)
    fig, ax = plt.subplots(figsize=(12.5, 9.5), constrained_layout=True)
    _draw_binary_lumen_context(ax, data, extent)
    narrow_count = int(np.count_nonzero(data.narrow_lumen_mask))
    if narrow_count > 0:
        glow = ndimage.binary_dilation(data.narrow_lumen_mask, iterations=_visibility_iterations(data)) & data.lumen_mask
        _imshow_boolean(ax, glow, extent, "#fff200", alpha=0.3)
        _imshow_boolean(ax, data.narrow_lumen_mask, extent, "#fff200", alpha=1.0)
        ax.contour(
            data.x_coordinates_um,
            data.z_coordinates_um,
            data.narrow_lumen_mask.T.astype(float),
            levels=[0.5],
            colors="#ff3b30",
            linewidths=0.12,
            alpha=0.75,
        )
        _mark_components(ax, data.narrow_lumen_mask, data, prefix="N")
    fraction = narrow_count / max(int(np.count_nonzero(data.lumen_mask)), 1)
    title = (
        f"Narrow lumen cells | cells={narrow_count} ({100.0 * fraction:.2f}%) | "
        f"effective diameter < 8 px ({8.0 * data.grid_spacing_um:.3g} um)"
    )
    _finish_quality_axes(ax, data, title)
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return output_path


def render_flow_field_map(
    data: FieldVisualizationData,
    output_path: Path,
    *,
    stage: FlowStage,
    dpi: int = 220,
) -> Path:
    """Save a static velocity/pressure visualization for one flow-field stage."""

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    velocity, speed, pressure, label = _flow_stage_arrays(data, stage)
    fig, axes = plt.subplots(1, 2, figsize=(15.0, 6.8), constrained_layout=True)
    speed_artist = _draw_speed_stream_panel(
        axes[0],
        data,
        speed,
        velocity,
        f"{label} speed and streamlines",
    )
    pressure_artist = _draw_pressure_quiver_panel(
        axes[1],
        data,
        pressure,
        velocity,
        _pressure_panel_title(stage, label),
    )
    fig.colorbar(speed_artist, ax=axes[0], label="Velocity [um/s]", fraction=0.046, pad=0.03)
    fig.colorbar(
        pressure_artist,
        ax=axes[1],
        label="Pressure [mmHg]" if stage == "final" else "Pressure reference",
        fraction=0.046,
        pad=0.03,
    )
    fig.suptitle(f"{data.result_dir.name} | {label} flow field", fontsize=12)
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return output_path


def render_wall_shear_stress_map(
    data: FieldVisualizationData,
    output_path: Path,
    *,
    dpi: int = 220,
) -> Path:
    """Save WSS only on closed solid-wall cells."""

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lower_pa, upper_pa = _wall_shear_color_limits(data)
    values = _wall_shear_values(data)
    median_pa = float(np.nanmedian(values)) if values.size else 0.0
    percentile_95_pa = float(np.nanpercentile(values, 95.0)) if values.size else 0.0

    fig, ax = plt.subplots(figsize=(10.5, 8.0), constrained_layout=True)
    artist = _draw_wall_shear_panel(
        ax,
        data,
        lower_pa=lower_pa,
        upper_pa=upper_pa,
        title=(
            "Final wall shear stress and flow direction\n"
            f"median={median_pa:.3g} Pa | P95={percentile_95_pa:.3g} Pa | "
            f"display capped at P99.5={upper_pa:.3g} Pa"
        ),
    )
    fig.colorbar(
        artist,
        ax=ax,
        label="Wall shear stress [Pa]",
        fraction=0.046,
        pad=0.03,
    )
    fig.suptitle(
        f"{data.result_dir.name} | closed-wall WSS distribution",
        fontsize=12,
    )
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return output_path


def render_local_shear_stress_map(
    data: FieldVisualizationData,
    output_path: Path,
    *,
    dpi: int = 220,
) -> Path:
    """Save the Newtonian local viscous shear magnitude throughout the lumen."""

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    lower_pa, upper_pa = _local_shear_color_limits(data)
    values = _local_shear_values(data)
    median_pa = float(np.nanmedian(values)) if values.size else 0.0
    percentile_95_pa = float(np.nanpercentile(values, 95.0)) if values.size else 0.0

    fig, ax = plt.subplots(figsize=(10.5, 8.0), constrained_layout=True)
    artist = _draw_local_shear_panel(
        ax,
        data,
        lower_pa=lower_pa,
        upper_pa=upper_pa,
        title=(
            "Final whole-lumen local viscous shear stress\n"
            f"median={median_pa:.3g} Pa | P95={percentile_95_pa:.3g} Pa | "
            f"display capped at P99.5={upper_pa:.3g} Pa"
        ),
    )
    fig.colorbar(
        artist,
        ax=ax,
        label="Local viscous shear stress [Pa]",
        fraction=0.046,
        pad=0.03,
    )
    fig.suptitle(
        f"{data.result_dir.name} | whole-lumen local shear distribution",
        fontsize=12,
    )
    fig.savefig(output_path, dpi=dpi)
    plt.close(fig)
    return output_path


def _save_blitted_animation(
    fig,
    frame_indices: np.ndarray,
    update,
    dynamic_artists: tuple[object, ...],
    background,
    output_path: Path,
    *,
    fps: int,
    show_progress: bool,
) -> None:
    """Draw only moving artists and stream encoded frames directly to disk."""

    writer = _animation_stream_writer(output_path, fps)
    progress = _progress_bar(frame_indices.size, show_progress)
    canvas = fig.canvas
    try:
        with writer:
            for frame in frame_indices:
                canvas.restore_region(background)
                update(int(frame))
                for artist in dynamic_artists:
                    fig.draw_artist(artist)
                canvas.blit(fig.bbox)
                writer.append(np.asarray(canvas.buffer_rgba())[..., :3])
                progress.update(1)
    finally:
        progress.close()


def _animation_stream_writer(output_path: Path, fps: int):
    """Create a bounded-memory GIF or MP4 frame sink."""

    if output_path.suffix.lower() == ".mp4":
        return _OpenCvMp4Writer(output_path, fps)
    return _PillowGifWriter(output_path, fps)


class _PillowGifWriter:
    """Stream globally-paletted delta frames instead of retaining an animation in RAM."""

    def __init__(self, output_path: Path, fps: int) -> None:
        self.output_path = Path(output_path)
        self.duration_ms = max(10, int(round(1000.0 / max(1, int(fps)))))
        self._file = None
        self._palette_frame = None
        self._previous_frame = None

    def __enter__(self):
        self._file = self.output_path.open("wb")
        return self

    def append(self, rgb: np.ndarray) -> None:
        from PIL import GifImagePlugin, Image, ImageChops

        image = Image.fromarray(np.asarray(rgb, dtype=np.uint8))
        if self._palette_frame is None:
            frame = image.quantize(colors=256, method=Image.Quantize.MEDIANCUT)
            self._palette_frame = frame.copy()
            self._file.writelines(
                GifImagePlugin._get_global_header(
                    frame,
                    {"duration": self.duration_ms, "loop": 0},
                )
            )
            offset = (0, 0)
        else:
            full_frame = image.quantize(
                palette=self._palette_frame,
                dither=Image.Dither.NONE,
            )
            changed = ImageChops.difference(self._previous_frame, full_frame).getbbox()
            if changed is None:
                changed = (0, 0, 1, 1)
            offset = (int(changed[0]), int(changed[1]))
            frame = full_frame.crop(changed)
        GifImagePlugin._write_frame_data(
            self._file,
            frame,
            offset,
            {"duration": self.duration_ms},
        )
        self._previous_frame = self._palette_frame.copy() if self._previous_frame is None else full_frame

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        if self._file is None:
            return
        if exc_type is None:
            self._file.write(b";")
        self._file.close()


class _OpenCvMp4Writer:
    """Encode MP4 frames incrementally with OpenCV's native video backend."""

    def __init__(self, output_path: Path, fps: int) -> None:
        self.output_path = Path(output_path)
        self.fps = max(1, int(fps))
        self._cv2 = None
        self._writer = None

    def __enter__(self):
        return self

    def append(self, rgb: np.ndarray) -> None:
        if self._cv2 is None:
            try:
                import cv2
            except ImportError as exc:  # pragma: no cover - depends on optional environment package.
                raise RuntimeError(
                    "MP4 output requires OpenCV (opencv-python); use a .gif output path otherwise."
                ) from exc
            self._cv2 = cv2

        frame = np.asarray(rgb, dtype=np.uint8)
        if self._writer is None:
            height, width = frame.shape[:2]
            codec = self._cv2.VideoWriter_fourcc(*"mp4v")
            self._writer = self._cv2.VideoWriter(
                str(self.output_path),
                codec,
                float(self.fps),
                (int(width), int(height)),
            )
            if not self._writer.isOpened():
                raise RuntimeError(f"OpenCV could not open the MP4 output: {self.output_path}")
        self._writer.write(self._cv2.cvtColor(frame, self._cv2.COLOR_RGB2BGR))

    def __exit__(self, exc_type, exc_value, traceback) -> None:
        if self._writer is not None:
            self._writer.release()


def _progress_bar(total: int, enabled: bool):
    """Return a tqdm progress bar, or a tiny no-op object when tqdm is unavailable."""

    try:
        import sys

        from tqdm.auto import tqdm

        return tqdm(total=int(total), desc="Rendering microbubble flow", unit="frame", disable=not bool(enabled), file=sys.stdout)
    except ImportError:
        return _NoProgress()


class _NoProgress:
    """Minimal fallback used when tqdm is not installed."""

    def update(self, _amount: int) -> None:
        return None

    def close(self) -> None:
        return None


def _flow_stage_arrays(
    data: FieldVisualizationData,
    stage: FlowStage,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    if stage == "initial":
        pressure = np.zeros_like(data.pressure, dtype=float)
        pressure[~data.lumen_mask] = np.nan
        return data.initial_velocity_xz_um_s, data.initial_speed_um_s, pressure, "Initial"
    return data.final_velocity_xz_um_s, data.speed_um_s, data.pressure, "Final"


def _pressure_panel_title(stage: FlowStage, label: str) -> str:
    if stage == "initial":
        return f"{label} pressure reference and velocity vectors"
    return f"{label} pressure-like scalar and velocity vectors"


def _draw_speed_stream_panel(
    ax,
    data: FieldVisualizationData,
    speed_um_s: np.ndarray,
    velocity_xz_um_s: np.ndarray,
    title: str,
):
    import matplotlib.pyplot as plt

    x, z, speed, lumen = _subsample_scalar_grid(data, speed_um_s, max_axis_points=850)
    levels = _contour_levels(speed_um_s, data.lumen_mask, count=48, zero_floor=True)
    masked_speed = np.ma.array(speed.T, mask=~lumen.T)
    artist = ax.contourf(x, z, masked_speed, levels=levels, cmap=plt.get_cmap("turbo"), extend="max")
    _draw_lumen_outline(ax, data, line_color="#111827", line_width=0.35)
    _draw_streamlines(ax, data, velocity_xz_um_s)
    _finish_flow_axes(ax, data, title)
    return artist


def _draw_wall_shear_panel(
    ax,
    data: FieldVisualizationData,
    *,
    lower_pa: float,
    upper_pa: float,
    title: str,
):
    """Draw WSS on a visually thickened closed-wall band."""

    import matplotlib.pyplot as plt

    shear = np.asarray(data.wall_shear_grid_pa, dtype=float)
    wall_mask = _wall_shear_mask(data)
    display_mask = (
        ndimage.binary_dilation(wall_mask, iterations=1)
        & np.asarray(data.lumen_mask, dtype=bool)
    )
    levels = np.linspace(float(lower_pa), float(upper_pa), 64)
    masked_shear = np.ma.array(
        shear.T,
        mask=(~display_mask.T) | (~np.isfinite(shear.T)) | (shear.T < 0.0),
    )
    artist = ax.imshow(
        masked_shear,
        origin="lower",
        extent=_image_extent(data),
        cmap=plt.get_cmap("turbo"),
        interpolation="nearest",
        aspect="equal",
        vmin=float(levels[0]),
        vmax=float(levels[-1]),
    )
    _draw_lumen_outline(ax, data, line_color="#111827", line_width=0.35)
    _draw_streamlines(ax, data, data.final_velocity_xz_um_s)
    _finish_flow_axes(ax, data, title)
    return artist


def _draw_local_shear_panel(
    ax,
    data: FieldVisualizationData,
    *,
    lower_pa: float,
    upper_pa: float,
    title: str,
):
    """Draw local viscous shear throughout the lumen."""

    import matplotlib.pyplot as plt

    x, z, shear, lumen = _subsample_scalar_grid(
        data,
        data.local_shear_grid_pa,
        max_axis_points=850,
    )
    levels = np.linspace(float(lower_pa), float(upper_pa), 64)
    shear_t = shear.T
    masked_shear = np.ma.array(
        shear_t,
        mask=(~lumen.T) | (~np.isfinite(shear_t)) | (shear_t < 0.0),
    )
    artist = ax.contourf(
        x,
        z,
        masked_shear,
        levels=levels,
        cmap=plt.get_cmap("turbo"),
        extend="max",
    )
    _draw_lumen_outline(ax, data, line_color="#111827", line_width=0.35)
    _draw_streamlines(ax, data, data.final_velocity_xz_um_s)
    _finish_flow_axes(ax, data, title)
    return artist


def _draw_pressure_quiver_panel(
    ax,
    data: FieldVisualizationData,
    pressure: np.ndarray,
    velocity_xz_um_s: np.ndarray,
    title: str,
):
    import matplotlib.pyplot as plt

    x, z, scalar, lumen = _subsample_scalar_grid(data, pressure, max_axis_points=850)
    levels = _contour_levels(pressure, data.lumen_mask, count=48, zero_floor=False)
    masked_pressure = np.ma.array(scalar.T, mask=~lumen.T)
    artist = ax.contourf(x, z, masked_pressure, levels=levels, cmap=plt.get_cmap("coolwarm"), extend="both")
    _draw_lumen_outline(ax, data, line_color="#111827", line_width=0.35)
    _draw_quiver(ax, data, velocity_xz_um_s)
    _finish_flow_axes(ax, data, title)
    return artist


def _draw_streamlines(ax, data: FieldVisualizationData, velocity_xz_um_s: np.ndarray) -> None:
    x, z, velocity, lumen = _subsample_vector_grid(data, velocity_xz_um_s, max_axis_points=220)
    speed = np.linalg.norm(velocity, axis=-1)
    if np.count_nonzero(lumen & np.isfinite(speed) & (speed > 0.0)) < 2:
        return
    u = np.ma.array(velocity[:, :, 0].T, mask=~lumen.T)
    v = np.ma.array(velocity[:, :, 1].T, mask=~lumen.T)
    for color, linewidth, zorder in (("#111827", 1.05, 5), ("#ffffff", 0.48, 6)):
        try:
            ax.streamplot(
                x,
                z,
                u,
                v,
                color=color,
                linewidth=linewidth,
                density=1.35,
                arrowsize=0.65,
                minlength=0.08,
                zorder=zorder,
            )
        except ValueError:
            return


def _draw_quiver(ax, data: FieldVisualizationData, velocity_xz_um_s: np.ndarray) -> None:
    x, z, velocity, lumen = _subsample_vector_grid(data, velocity_xz_um_s, max_axis_points=46)
    speed = np.linalg.norm(velocity, axis=-1)
    valid = lumen & np.isfinite(speed) & (speed > 0.0)
    if np.count_nonzero(valid) == 0:
        return
    grid_x, grid_z = np.meshgrid(x, z, indexing="xy")
    u = velocity[:, :, 0].T
    v = velocity[:, :, 1].T
    valid_t = valid.T
    ax.quiver(
        grid_x[valid_t],
        grid_z[valid_t],
        u[valid_t],
        v[valid_t],
        color="#111827",
        angles="xy",
        scale_units="xy",
        scale=_quiver_scale(data, velocity_xz_um_s),
        width=0.0024,
        alpha=0.85,
        zorder=6,
    )


def _draw_lumen_outline(ax, data: FieldVisualizationData, *, line_color: str, line_width: float) -> None:
    ax.contour(
        data.x_coordinates_um,
        data.z_coordinates_um,
        data.lumen_mask.T.astype(float),
        levels=[0.5],
        colors=line_color,
        linewidths=float(line_width),
        alpha=0.95,
        zorder=7,
    )


def _finish_flow_axes(ax, data: FieldVisualizationData, title: str) -> None:
    x_min, x_max, z_min, z_max = _image_extent(data)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(z_min, z_max)
    ax.set_xlabel("X [um]")
    ax.set_ylabel("Z [um]")
    ax.set_title(title, fontsize=11)
    ax.set_facecolor("#808080")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    for spine in ax.spines.values():
        spine.set_visible(False)


def _subsample_scalar_grid(
    data: FieldVisualizationData,
    scalar: np.ndarray,
    *,
    max_axis_points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    stride = _grid_stride(data, max_axis_points)
    return (
        data.x_coordinates_um[::stride],
        data.z_coordinates_um[::stride],
        np.asarray(scalar, dtype=float)[::stride, ::stride],
        data.lumen_mask[::stride, ::stride],
    )


def _subsample_vector_grid(
    data: FieldVisualizationData,
    velocity_xz_um_s: np.ndarray,
    *,
    max_axis_points: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    stride = _grid_stride(data, max_axis_points)
    return (
        data.x_coordinates_um[::stride],
        data.z_coordinates_um[::stride],
        np.asarray(velocity_xz_um_s, dtype=float)[::stride, ::stride, :],
        data.lumen_mask[::stride, ::stride],
    )


def _grid_stride(data: FieldVisualizationData, max_axis_points: int) -> int:
    longest_axis = max(int(data.x_coordinates_um.size), int(data.z_coordinates_um.size))
    return max(1, int(np.ceil(longest_axis / max(1, int(max_axis_points)))))


def _contour_levels(scalar: np.ndarray, mask: np.ndarray, *, count: int, zero_floor: bool) -> np.ndarray:
    values = np.asarray(scalar, dtype=float)[np.asarray(mask, dtype=bool)]
    values = values[np.isfinite(values)]
    if values.size == 0:
        return np.linspace(0.0, 1.0, int(count))
    lo = 0.0 if zero_floor else float(np.nanpercentile(values, 1.0))
    hi = float(np.nanpercentile(values, 99.0))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        center = float(np.nanmedian(values)) if values.size else 0.0
        span = max(abs(center) * 0.1, 1.0)
        lo = 0.0 if zero_floor else center - span
        hi = max(center + span, lo + span)
    return np.linspace(float(lo), float(hi), int(count))


def _wall_shear_values(data: FieldVisualizationData) -> np.ndarray:
    """Return finite, non-negative WSS values on closed solid-wall cells."""

    shear = np.asarray(data.wall_shear_grid_pa, dtype=float)
    values = shear[_wall_shear_mask(data)]
    return values[np.isfinite(values) & (values >= 0.0)]


def _wall_shear_mask(data: FieldVisualizationData) -> np.ndarray:
    """Return the closed-wall display mask, with a legacy test-data fallback."""

    lumen = np.asarray(data.lumen_mask, dtype=bool)
    explicit = getattr(data, "wall_shear_display_mask", None)
    if explicit is not None:
        mask = np.asarray(explicit, dtype=bool)
        if mask.shape == lumen.shape:
            return mask & lumen
    wall = getattr(data, "wall_mask", None)
    if wall is not None:
        mask = np.asarray(wall, dtype=bool)
        if mask.shape == lumen.shape:
            return mask & lumen
    return lumen


def _local_shear_values(data: FieldVisualizationData) -> np.ndarray:
    """Return finite, non-negative local shear values throughout the lumen."""

    shear = np.asarray(data.local_shear_grid_pa, dtype=float)
    values = shear[np.asarray(data.lumen_mask, dtype=bool)]
    return values[np.isfinite(values) & (values >= 0.0)]


def _wall_shear_color_limits(
    data: FieldVisualizationData,
) -> tuple[float, float]:
    """Use one robust WSS scale for static maps, backgrounds, and bubbles."""

    values = _wall_shear_values(data)
    if values.size == 0:
        return 0.0, 1.0
    upper = float(np.nanpercentile(values, 99.5))
    if not np.isfinite(upper) or upper <= np.finfo(float).eps:
        upper = max(float(np.nanmax(values, initial=0.0)), 1.0)
    return 0.0, upper


def _local_shear_color_limits(
    data: FieldVisualizationData,
) -> tuple[float, float]:
    """Use one robust scale for the whole-lumen local shear field."""

    values = _local_shear_values(data)
    if values.size == 0:
        return 0.0, 1.0
    upper = float(np.nanpercentile(values, 99.5))
    if not np.isfinite(upper) or upper <= np.finfo(float).eps:
        upper = max(float(np.nanmax(values, initial=0.0)), 1.0)
    return 0.0, upper


def _quiver_scale(data: FieldVisualizationData, velocity_xz_um_s: np.ndarray) -> float:
    speed = np.linalg.norm(np.asarray(velocity_xz_um_s, dtype=float), axis=-1)
    values = speed[np.asarray(data.lumen_mask, dtype=bool) & np.isfinite(speed)]
    if values.size == 0:
        return 1.0
    speed_scale = float(np.nanpercentile(values, 95.0))
    if speed_scale <= np.finfo(float).eps:
        return 1.0
    x_min, x_max, z_min, z_max = _image_extent(data)
    target_arrow_um = max(x_max - x_min, z_max - z_min) / 42.0
    return speed_scale / max(target_arrow_um, np.finfo(float).eps)


def _draw_background(
    ax,
    data: FieldVisualizationData,
    *,
    color_mode: ColorMode = "speed",
    color_limits: tuple[float, float] | None = None,
) -> None:
    """Draw the spatial field corresponding to the selected bubble scalar."""

    wall_shear_mode = color_mode == "wall_shear"
    local_shear_mode = color_mode == "local_shear"
    if wall_shear_mode:
        scalar = np.asarray(data.wall_shear_grid_pa, dtype=float)
        spatial_mask = _wall_shear_mask(data)
    elif local_shear_mode:
        scalar = np.asarray(data.local_shear_grid_pa, dtype=float)
        spatial_mask = np.asarray(data.lumen_mask, dtype=bool)
    else:
        scalar = np.asarray(data.speed_um_s, dtype=float)
        spatial_mask = np.asarray(data.lumen_mask, dtype=bool)
    invalid_scalar = ~np.isfinite(scalar)
    if wall_shear_mode or local_shear_mode:
        invalid_scalar |= scalar < 0.0
    masked_scalar = np.ma.array(
        scalar,
        mask=(~spatial_mask) | invalid_scalar,
    )
    extent = _image_extent(data)
    image = ax.imshow(
        masked_scalar.T,
        origin="lower",
        extent=extent,
        cmap="turbo" if (wall_shear_mode or local_shear_mode) else "magma",
        alpha=0.88 if (wall_shear_mode or local_shear_mode) else 0.62,
        interpolation="nearest" if wall_shear_mode else "bilinear",
        aspect="equal",
    )
    if color_limits is None:
        if wall_shear_mode:
            color_limits = _wall_shear_color_limits(data)
        elif local_shear_mode:
            color_limits = _local_shear_color_limits(data)
        else:
            color_limits = (
                0.0,
                max(
                    float(
                        np.nanpercentile(
                            data.speed_um_s[data.lumen_mask],
                            99.0,
                        )
                    ),
                    1.0,
                ),
            )
    image.set_clim(*color_limits)
    _draw_molecular_target_overlay(ax, data, extent)
    ax.contour(
        data.x_coordinates_um,
        data.z_coordinates_um,
        data.lumen_mask.T.astype(float),
        levels=[0.5],
        colors="#d9f2ff",
        linewidths=0.45,
        alpha=0.9,
    )


def _draw_molecular_target_overlay(
    ax,
    data: FieldVisualizationData,
    extent: tuple[float, float, float, float],
) -> None:
    """Draw exact target-positive wall cells with a separate subtle visibility halo."""

    target = data.molecular_target
    if not target.visible:
        return

    from matplotlib.patches import Patch

    exact_mask, halo_mask = target_display_masks(target.target_wall_mask, halo_cells=2)
    _imshow_boolean(ax, halo_mask, extent, "#65ffd4", alpha=0.24, zorder=1.25)
    _imshow_boolean(ax, exact_mask, extent, "#00f5ad", alpha=0.96, zorder=1.5)
    ax.contour(
        data.x_coordinates_um,
        data.z_coordinates_um,
        exact_mask.T.astype(float),
        levels=[0.5],
        colors="#eafff8",
        linewidths=0.75,
        alpha=1.0,
        zorder=2.2,
    )

    density_um2 = float(target.target_density_molecules_per_m2) * 1.0e-12
    label = (
        "Molecular target wall (geometry only; density=0)"
        if density_um2 <= 0.0
        else f"Molecular target wall ({density_um2:.4g} molecules/um^2)"
    )
    ax.legend(
        handles=[Patch(facecolor="#00f5ad", edgecolor="#eafff8", label=label)],
        loc="upper right",
        fontsize=8,
        framealpha=0.82,
    )


def _draw_binary_lumen_context(ax, data: FieldVisualizationData, extent: tuple[float, float, float, float]) -> None:
    """Draw a dark, simple lumen background for mask-quality maps."""

    _imshow_boolean(ax, data.lumen_mask, extent, "#1f2937", alpha=1.0)
    ax.contour(
        data.x_coordinates_um,
        data.z_coordinates_um,
        data.lumen_mask.T.astype(float),
        levels=[0.5],
        colors="#f8fafc",
        linewidths=0.45,
        alpha=0.95,
    )


def _imshow_boolean(
    ax,
    mask: np.ndarray,
    extent: tuple[float, float, float, float],
    color: str,
    *,
    alpha: float,
    zorder: float = 0.0,
) -> None:
    """Draw a boolean mask as one solid color."""

    from matplotlib.colors import ListedColormap

    cmap = ListedColormap([color])
    cmap.set_bad(alpha=0.0)
    layer = np.ma.masked_where(~np.asarray(mask, dtype=bool).T, np.ones_like(mask.T, dtype=float))
    ax.imshow(
        layer,
        origin="lower",
        extent=extent,
        cmap=cmap,
        alpha=float(alpha),
        interpolation="nearest",
        aspect="equal",
        zorder=float(zorder),
    )


def _mark_components(ax, mask: np.ndarray, data: FieldVisualizationData, *, prefix: str) -> None:
    """Draw large component markers so small defects are visible in full-network views."""

    labels, count = ndimage.label(np.asarray(mask, dtype=bool), structure=np.ones((3, 3), dtype=bool))
    if count <= 0:
        return
    centers = ndimage.center_of_mass(mask, labels, index=np.arange(1, count + 1))
    xs = np.asarray([np.interp(center[0], np.arange(data.x_coordinates_um.size), data.x_coordinates_um) for center in centers])
    zs = np.asarray([np.interp(center[1], np.arange(data.z_coordinates_um.size), data.z_coordinates_um) for center in centers])
    ax.scatter(xs, zs, s=280, facecolors="none", edgecolors="#ffffff", linewidths=2.2, zorder=8)
    ax.scatter(xs, zs, s=80, c="#ff1744", edgecolors="#000000", linewidths=0.8, zorder=9)
    for idx, (x, z) in enumerate(zip(xs, zs, strict=False), start=1):
        ax.text(
            float(x),
            float(z),
            f"{prefix}{idx}",
            color="#ffffff",
            fontsize=8,
            fontweight="bold",
            ha="left",
            va="bottom",
            zorder=10,
        )


def _visibility_iterations(data: FieldVisualizationData) -> int:
    """Return a dilation radius that keeps tiny hole cells visible on full maps."""

    return max(3, int(round(14.0 / max(float(data.grid_spacing_um), np.finfo(float).eps))))


def _finish_quality_axes(ax, data: FieldVisualizationData, title: str) -> None:
    """Apply high-contrast axis styling for static mask-quality maps."""

    x_min, x_max, z_min, z_max = _image_extent(data)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(z_min, z_max)
    ax.set_xlabel("X (um)")
    ax.set_ylabel("Z (um)")
    ax.set_title(title, fontsize=12, color="#111827")
    ax.set_facecolor("#05070b")
    ax.grid(False)


def _draw_tails(ax, data: FieldVisualizationData, selected: np.ndarray, frame: int, tail_length: int, line_collection_cls) -> None:
    """Add short trajectory tails behind the current bubble positions."""

    collection = line_collection_cls(
        _tail_segments(data, selected, frame, tail_length),
        colors="#ffffff",
        linewidths=0.85,
        alpha=0.55,
        zorder=2.5,
    )
    ax.add_collection(collection)


def _draw_bubbles(
    ax,
    data: FieldVisualizationData,
    selected: np.ndarray,
    frame: int,
    color_mode: ColorMode,
    ellipse_collection_cls,
    *,
    color_limit_selection: np.ndarray | None = None,
    glow_style: BubbleGlowStyle | None = None,
    color_limits: tuple[float, float] | None = None,
):
    """Draw true-size disks plus a non-physical, gradually fading visibility halo."""

    xy = _bubble_xy(data, selected, frame)
    values = _bubble_color_values(data, selected, frame, color_mode)
    label = _colorbar_label(color_mode)
    limit_selection = selected if color_limit_selection is None else color_limit_selection
    vmin, vmax = (
        _color_limits(data, limit_selection, color_mode)
        if color_limits is None
        else color_limits
    )
    diameters = _bubble_physical_diameters(data, selected, frame)
    rotation_angles = _bubble_rotation_angles(data, selected, frame)
    style = glow_style or BubbleGlowStyle()

    # 光晕位于实体圆盘之后；它只增强辨识度，可以越过壁线，但不代表微泡几何尺寸。
    halo_artists: list[object] = []
    halo_points_per_um = _points_per_data_um(ax)
    layer_diameters = glow_layer_diameters_points(diameters, halo_points_per_um, style)
    layer_alphas = style.layer_alphas if style.enabled else ()
    for layer_index, (halo_diameters, alpha) in enumerate(
        zip(layer_diameters, layer_alphas, strict=True)
    ):
        halo = ellipse_collection_cls(
            widths=halo_diameters,
            heights=halo_diameters,
            angles=0.0,
            units="points",
            offsets=xy,
            offset_transform=ax.transData,
            facecolors=(*style.glow_rgb, 1.0),
            edgecolors="none",
            linewidths=0.0,
            alpha=float(alpha),
            zorder=3.0 + 0.1 * layer_index,
        )
        ax.add_collection(halo)
        halo_artists.append(halo)

    # Near-wall transport and target adhesion use a separate state ring outside
    # the visibility halo, leaving the scalar-filled physical core unchanged.
    contact_ring = ellipse_collection_cls(
        widths=_contact_ring_diameters_points(diameters, halo_points_per_um, style),
        heights=_contact_ring_diameters_points(diameters, halo_points_per_um, style),
        angles=0.0,
        units="points",
        offsets=xy,
        offset_transform=ax.transData,
        facecolors="none",
        edgecolors=_bubble_contact_ring_colors(data, selected, frame),
        linewidths=1.15,
        zorder=3.8,
    )
    ax.add_collection(contact_ring)

    core = ellipse_collection_cls(
        widths=diameters,
        heights=diameters,
        angles=0.0,
        units="xy",
        offsets=xy,
        offset_transform=ax.transData,
        array=values,
        cmap=(
            "turbo"
            if color_mode in ("wall_shear", "local_shear")
            else "viridis"
        ),
        edgecolors=_bubble_edge_colors(data, selected, frame),
        linewidths=0.32,
        alpha=0.98,
        zorder=4.0,
    )
    core.set_clim(vmin, vmax)
    ax.add_collection(core)

    # 高光完全位于真实圆盘内部，使圆形读起来更接近发光球体，而不会扩大实体边界。
    highlight_offsets, highlight_diameters = bubble_highlight_geometry(
        xy,
        diameters,
        style,
        rotation_angles,
    )
    highlight = ellipse_collection_cls(
        widths=highlight_diameters,
        heights=highlight_diameters,
        angles=0.0,
        units="xy",
        offsets=highlight_offsets,
        offset_transform=ax.transData,
        facecolors=(1.0, 1.0, 1.0, float(style.highlight_alpha)),
        edgecolors="none",
        linewidths=0.0,
        zorder=4.2,
    )
    ax.add_collection(highlight)
    return _BubbleArtists(
        core=core,
        halos=tuple(halo_artists),
        contact_ring=contact_ring,
        highlight=highlight,
        halo_points_per_um=halo_points_per_um,
    ), label


def _update_bubble_artists(
    artists: _BubbleArtists,
    xy: np.ndarray,
    values: np.ndarray,
    diameters: np.ndarray,
    edgecolors: np.ndarray,
    contact_ring_colors: np.ndarray,
    rotation_angles: np.ndarray,
    style: BubbleGlowStyle,
) -> None:
    """Update every glow/material layer for one animation frame."""

    for halo, layer_diameters in zip(
        artists.halos,
        glow_layer_diameters_points(diameters, artists.halo_points_per_um, style),
        strict=True,
    ):
        halo.set_offsets(xy)
        halo.set_widths(layer_diameters)
        halo.set_heights(layer_diameters)

    artists.core.set_offsets(xy)
    artists.core.set_array(values)
    artists.core.set_widths(diameters)
    artists.core.set_heights(diameters)
    artists.core.set_edgecolors(edgecolors)

    if artists.contact_ring is not None:
        contact_diameters = _contact_ring_diameters_points(
            diameters,
            artists.halo_points_per_um,
            style,
        )
        artists.contact_ring.set_offsets(xy)
        artists.contact_ring.set_widths(contact_diameters)
        artists.contact_ring.set_heights(contact_diameters)
        artists.contact_ring.set_edgecolors(contact_ring_colors)

    if artists.highlight is not None:
        highlight_offsets, highlight_diameters = bubble_highlight_geometry(
            xy,
            diameters,
            style,
            rotation_angles,
        )
        artists.highlight.set_offsets(highlight_offsets)
        artists.highlight.set_widths(highlight_diameters)
        artists.highlight.set_heights(highlight_diameters)


def _points_per_data_um(ax) -> float:
    """Measure the current equal-aspect transform so a halo can have a modest screen-size floor."""

    origin_px = np.asarray(ax.transData.transform((0.0, 0.0)), dtype=float)
    x_unit_px = np.linalg.norm(np.asarray(ax.transData.transform((1.0, 0.0))) - origin_px)
    z_unit_px = np.linalg.norm(np.asarray(ax.transData.transform((0.0, 1.0))) - origin_px)
    pixels_per_um = min(float(x_unit_px), float(z_unit_px))
    return max(pixels_per_um * 72.0 / float(ax.figure.dpi), np.finfo(float).eps)


def _finish_axes(ax, data: FieldVisualizationData, title: str) -> None:
    """Apply axis labels, limits, and title consistently."""

    x_min, x_max, z_min, z_max = _image_extent(data)
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(z_min, z_max)
    ax.set_xlabel("X (um)")
    ax.set_ylabel("Z (um)")
    ax.set_title(title, fontsize=11)
    ax.set_facecolor("#101722")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)


def _select_bubbles(data: FieldVisualizationData, max_bubbles: int) -> np.ndarray:
    """Select all current lanes by default or an evenly distributed display subset."""

    n_bubbles = data.positions_um.shape[1]
    requested = int(max_bubbles)
    if requested <= 0 or requested >= n_bubbles:
        return np.arange(n_bubbles, dtype=int)
    return np.unique(np.rint(np.linspace(0, n_bubbles - 1, requested)).astype(int))


def _active_frame_selection(data: FieldVisualizationData, selected: np.ndarray, frame: int) -> np.ndarray:
    """Return only real, active microbubble observations from one stored frame."""

    indices = np.asarray(selected, dtype=int)
    if indices.size == 0:
        return indices
    frame_index = int(frame)
    positions = data.positions_um[frame_index, indices, :]
    valid = (
        data.active[frame_index, indices]
        & (data.bubble_id[frame_index, indices] >= 0)
        & np.all(np.isfinite(positions), axis=1)
    )
    return indices[valid]


def _tail_segments(data: FieldVisualizationData, selected: np.ndarray, frame: int, tail_length: int) -> list[np.ndarray]:
    """Build tails without connecting two permanent bubble IDs that used one display lane."""

    start = max(0, int(frame) - max(1, int(tail_length)))
    positions = data.positions_um[start : int(frame) + 1, selected, :]
    bubble_ids = data.bubble_id[start : int(frame) + 1, selected]
    active = data.active[start : int(frame) + 1, selected]
    segments: list[np.ndarray] = []
    for bubble_idx in range(positions.shape[1]):
        if not bool(active[-1, bubble_idx]):
            continue
        points = positions[:, bubble_idx, :][:, [0, 2]]
        current_id = int(bubble_ids[-1, bubble_idx])
        same_physical_bubble = bubble_ids[:, bubble_idx] == current_id
        finite = np.all(np.isfinite(points), axis=1) & same_physical_bubble & active[:, bubble_idx]
        mismatches = np.flatnonzero(~same_physical_bubble)
        if mismatches.size:
            finite[: int(mismatches[-1]) + 1] = False
        if np.count_nonzero(finite) >= 2:
            segments.append(points[finite])
    return segments


def _particle_frame_title(
    data: FieldVisualizationData,
    frame: int,
    frames: int,
    *,
    displayed: int | None = None,
) -> str:
    """Return population-aware title text for snapshots and animations."""

    active = int(data.active_count_per_frame[int(frame)])
    injected = int(data.injected_count_per_frame[int(frame)])
    terminated = int(data.terminated_count_per_frame[int(frame)])
    title = (
        f"{data.result_dir.name} | frame {int(frame)}/{int(frames) - 1} | "
        f"active {active} (peak {int(data.peak_active_bubbles)}) | entered +{injected} / exited -{terminated}"
    )
    frame_gaps = np.asarray(data.bubble_wall_gap_um[int(frame)], dtype=float)
    finite_gaps = frame_gaps[np.isfinite(frame_gaps) & np.asarray(data.active[int(frame)], dtype=bool)]
    if finite_gaps.size:
        title += f" | min gap {float(np.min(finite_gaps)):.3g} um"
    contact_count = int(
        np.count_nonzero(
            np.asarray(data.bubble_wall_contact[int(frame)], dtype=bool)
            & np.asarray(data.active[int(frame)], dtype=bool)
        )
    )
    if contact_count:
        title += f" | contact {contact_count}"
    near_wall_count = int(
        np.count_nonzero(
            _bubble_near_wall_mask(data, np.flatnonzero(data.active[int(frame)]), frame)
        )
    )
    adhesion_count = int(
        np.count_nonzero(
            _bubble_target_adhesion_mask(
                data,
                np.flatnonzero(data.active[int(frame)]),
                frame,
            )
        )
    )
    if near_wall_count:
        title += f" | near wall {near_wall_count}"
    if adhesion_count:
        title += f" | target adhesion {adhesion_count}"
    if displayed is not None and int(displayed) != active:
        title += f" | displayed {int(displayed)}"
    return title


def _bubble_xy(data: FieldVisualizationData, selected: np.ndarray, frame: int) -> np.ndarray:
    """Return current X-Z coordinates for selected bubbles."""

    return data.positions_um[int(frame), selected, :][:, [0, 2]]


def _bubble_color_values(data: FieldVisualizationData, selected: np.ndarray, frame: int, color_mode: ColorMode) -> np.ndarray:
    """Return per-bubble scalar values for the chosen color mode."""

    frame = int(frame)
    if color_mode == "wall_shear":
        return data.bubble_wall_shear_pa[frame, selected].astype(float)
    if color_mode == "local_shear":
        return _sample_scalar_grid_at_bubbles(
            data,
            data.local_shear_grid_pa,
            selected,
            frame,
        )
    if color_mode == "vessel_id":
        return data.bubble_vessel_id[frame, selected].astype(float)
    if color_mode == "diameter":
        return data.bubble_diameter_um[frame, selected].astype(float)
    if color_mode == "clearance":
        return data.bubble_wall_gap_um[frame, selected].astype(float)
    if color_mode == "active":
        return data.active[frame, selected].astype(float)
    velocity = _display_velocity_array(data)
    return np.linalg.norm(velocity[frame, selected, :][:, [0, 2]], axis=1)


def _color_limits(data: FieldVisualizationData, selected: np.ndarray, color_mode: ColorMode) -> tuple[float, float]:
    """Compute stable color limits across selected bubbles and all frames."""

    if color_mode == "wall_shear":
        return _wall_shear_color_limits(data)
    elif color_mode == "local_shear":
        return _local_shear_color_limits(data)
    elif color_mode == "vessel_id":
        values = data.bubble_vessel_id[:, selected].astype(float).ravel()
        values = values[values > 0]
    elif color_mode == "diameter":
        values = data.bubble_diameter_um[:, selected].astype(float).ravel()
    elif color_mode == "clearance":
        values = data.bubble_wall_gap_um[:, selected].astype(float).ravel()
    elif color_mode == "active":
        return 0.0, 1.0
    else:
        velocity = _display_velocity_array(data)
        values = np.linalg.norm(velocity[:, selected, :][..., [0, 2]], axis=-1).ravel()

    active = data.active[:, selected].ravel()
    if values.size == active.size:
        values = values[active]
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 0.0, 1.0
    lo = float(np.nanpercentile(values, 1.0))
    hi = float(np.nanpercentile(values, 99.0))
    if hi <= lo:
        hi = lo + 1.0
    return lo, hi


def _colorbar_label(color_mode: ColorMode) -> str:
    """Map color modes to human-readable colorbar labels."""

    if color_mode == "wall_shear":
        return "sampled wall shear stress (Pa)"
    if color_mode == "local_shear":
        return "local viscous shear stress (Pa)"
    if color_mode == "vessel_id":
        return "sampled vessel ID"
    if color_mode == "diameter":
        return "bubble diameter (um)"
    if color_mode == "clearance":
        return "bubble-to-wall gap (um)"
    if color_mode == "active":
        return "active"
    return "realized center speed (um/s)"


def _sample_scalar_grid_at_bubbles(
    data: FieldVisualizationData,
    scalar_grid: np.ndarray,
    selected: np.ndarray,
    frame: int,
) -> np.ndarray:
    """Bilinearly sample a regular X-Z scalar grid at bubble centres."""

    positions = _bubble_xy(data, selected, int(frame))
    x = np.asarray(data.x_coordinates_um, dtype=float)
    z = np.asarray(data.z_coordinates_um, dtype=float)
    scalar = np.asarray(scalar_grid, dtype=float)
    if scalar.shape != (x.size, z.size):
        raise ValueError("Scalar grid shape does not match visualization coordinates.")
    if x.size < 2 or z.size < 2:
        return np.full(positions.shape[0], np.nan, dtype=float)

    result = np.full(positions.shape[0], np.nan, dtype=float)
    valid = np.all(np.isfinite(positions), axis=1)
    if not np.any(valid):
        return result
    ix = np.interp(positions[valid, 0], x, np.arange(x.size, dtype=float))
    iz = np.interp(positions[valid, 1], z, np.arange(z.size, dtype=float))
    i0 = np.clip(np.floor(ix).astype(np.int64), 0, x.size - 2)
    j0 = np.clip(np.floor(iz).astype(np.int64), 0, z.size - 2)
    fx = np.clip(ix - i0, 0.0, 1.0)
    fz = np.clip(iz - j0, 0.0, 1.0)
    weights = np.column_stack(
        (
            (1.0 - fx) * (1.0 - fz),
            fx * (1.0 - fz),
            (1.0 - fx) * fz,
            fx * fz,
        )
    )
    samples = np.column_stack(
        (
            scalar[i0, j0],
            scalar[i0 + 1, j0],
            scalar[i0, j0 + 1],
            scalar[i0 + 1, j0 + 1],
        )
    )
    finite = np.isfinite(samples)
    effective_weights = np.where(finite, weights, 0.0)
    weight_sum = np.sum(effective_weights, axis=1)
    sampled = np.divide(
        np.sum(np.where(finite, samples, 0.0) * effective_weights, axis=1),
        weight_sum,
        out=np.full(weight_sum.shape, np.nan, dtype=float),
        where=weight_sum > 0.0,
    )
    result[valid] = sampled
    return result


def _display_velocity_array(data: FieldVisualizationData) -> np.ndarray:
    """Prefer post-constraint, same-ID displacement velocity for visual speed."""

    realized = getattr(data, "realized_velocities_um_s", None)
    if realized is not None:
        values = np.asarray(realized, dtype=float)
        if values.shape == np.asarray(data.positions_um).shape:
            return values
    return np.asarray(data.velocities_um_s, dtype=float)


def _bubble_physical_diameters(data: FieldVisualizationData, selected: np.ndarray, frame: int) -> np.ndarray:
    """Return physical disk diameters in the X-Z data-coordinate unit, micrometers."""

    diameters = np.asarray(data.bubble_diameter_um[int(frame), selected], dtype=float)
    return np.where(np.isfinite(diameters) & (diameters > 0.0), diameters, 2.0)


def _bubble_rotation_angles(data: FieldVisualizationData, selected: np.ndarray, frame: int) -> np.ndarray:
    """Return saved +Y rotation angles, or a stationary marker for legacy results."""

    rotation = getattr(data, "bubble_rotation_angle_rad", None)
    if rotation is None:
        return np.zeros(np.asarray(selected).size, dtype=float)
    values = np.asarray(rotation, dtype=float)[int(frame), selected]
    return np.where(np.isfinite(values), values, 0.0)


def _contact_ring_diameters_points(
    diameters_um: np.ndarray,
    points_per_um: float,
    style: BubbleGlowStyle,
) -> np.ndarray:
    """Place the diagnostic contact ring outside the core and optional halo."""

    physical_points = np.asarray(diameters_um, dtype=float) * float(points_per_um)
    if style.enabled:
        halo_layers = glow_layer_diameters_points(diameters_um, points_per_um, style)
        base = halo_layers[0] if halo_layers else physical_points
    else:
        base = physical_points
    return np.maximum(np.asarray(base, dtype=float) * 1.14, 6.8)


def _bubble_edge_colors(data: FieldVisualizationData, selected: np.ndarray, frame: int) -> np.ndarray:
    """Keep the physical core neutral, reserving red for numerical penetration."""

    gaps = np.asarray(data.bubble_wall_gap_um[int(frame), selected], dtype=float)
    colors = np.tile(np.asarray([0.0, 0.0, 0.0, 1.0]), (selected.size, 1))
    invalid_below_um = float(data.wall_gap_invalid_below_um)
    colors[np.isfinite(gaps) & (gaps < invalid_below_um)] = np.asarray(
        [1.0, 0.0, 0.0, 1.0]
    )
    return colors


def _bubble_contact_ring_colors(data: FieldVisualizationData, selected: np.ndarray, frame: int) -> np.ndarray:
    """Color state rings: orange near-wall, pink target adhesion, red invalid gap."""

    gaps = np.asarray(data.bubble_wall_gap_um[int(frame), selected], dtype=float)
    colors = np.zeros((selected.size, 4), dtype=float)
    colors[_bubble_near_wall_mask(data, selected, frame)] = _NEAR_WALL_RING_RGBA
    colors[_bubble_target_adhesion_mask(data, selected, frame)] = (
        _TARGET_ADHESION_RING_RGBA
    )
    invalid_below_um = float(data.wall_gap_invalid_below_um)
    colors[np.isfinite(gaps) & (gaps < invalid_below_um)] = _INVALID_GAP_RING_RGBA
    return colors


def _bubble_near_wall_mask(
    data: FieldVisualizationData,
    selected: np.ndarray,
    frame: int,
) -> np.ndarray:
    """Use saved near-wall weights, falling back to the configured gap ratio."""

    selected = np.asarray(selected, dtype=np.int64)
    weights = getattr(data, "bubble_near_wall_weight", None)
    if weights is None:
        frame_weights = np.full(selected.size, np.nan, dtype=float)
    else:
        frame_weights = np.asarray(weights, dtype=float)[int(frame), selected]

    near_wall = np.isfinite(frame_weights) & (frame_weights > 0.0)
    missing = ~np.isfinite(frame_weights)
    if np.any(missing):
        gaps = np.asarray(data.bubble_wall_gap_um[int(frame), selected], dtype=float)
        contacts = getattr(data, "bubble_wall_contact", None)
        if contacts is not None:
            near_wall[missing] = np.asarray(contacts, dtype=bool)[int(frame), selected][
                missing
            ]
        diameters = getattr(data, "bubble_diameter_um", None)
        if diameters is not None:
            frame_diameters = np.asarray(diameters, dtype=float)[int(frame), selected]
            xi_far = max(float(getattr(data, "near_wall_xi_far", 1.0)), 0.0)
            fallback_limit_um = 0.5 * frame_diameters * xi_far
            near_wall[missing] |= (
                np.isfinite(gaps[missing])
                & np.isfinite(fallback_limit_um[missing])
                & (gaps[missing] <= fallback_limit_um[missing])
            )
    return near_wall


def _bubble_target_adhesion_mask(
    data: FieldVisualizationData,
    selected: np.ndarray,
    frame: int,
) -> np.ndarray:
    """Treat a positive expected molecular bond count as target adhesion."""

    selected = np.asarray(selected, dtype=np.int64)
    bond_counts = getattr(data, "bubble_bond_count_expected", None)
    if bond_counts is None:
        return np.zeros(selected.size, dtype=bool)
    values = np.asarray(bond_counts, dtype=float)[int(frame), selected]
    return np.isfinite(values) & (values > 0.0)


def _draw_bubble_state_legend(ax) -> None:
    """Add ring-color keys while preserving the molecular-target wall key."""

    from matplotlib.lines import Line2D

    handles: list[object] = []
    labels: list[str] = []
    existing = ax.get_legend()
    if existing is not None:
        handles.extend(getattr(existing, "legend_handles", ()))
        labels.extend(text.get_text() for text in existing.get_texts())

    for color, label in (
        (_NEAR_WALL_RING_RGBA, "Near-wall bubble"),
        (_TARGET_ADHESION_RING_RGBA, "Target-adhered bubble"),
    ):
        handles.append(
            Line2D(
                [],
                [],
                linestyle="none",
                marker="o",
                markerfacecolor="none",
                markeredgecolor=color,
                markeredgewidth=1.8,
                markersize=7.0,
            )
        )
        labels.append(label)
    ax.legend(
        handles=handles,
        labels=labels,
        loc="upper right",
        fontsize=8,
        framealpha=0.82,
    )


def _image_extent(data: FieldVisualizationData) -> tuple[float, float, float, float]:
    """Return Matplotlib imshow extent in world X-Z coordinates."""

    return (
        float(data.x_coordinates_um[0]),
        float(data.x_coordinates_um[-1]),
        float(data.z_coordinates_um[0]),
        float(data.z_coordinates_um[-1]),
    )


def _normalize_frame_index(frame_index: int, frame_count: int) -> int:
    """Convert negative frame indices and clamp them into valid range."""

    frame = int(frame_index)
    if frame < 0:
        frame = frame_count + frame
    return int(np.clip(frame, 0, frame_count - 1))
