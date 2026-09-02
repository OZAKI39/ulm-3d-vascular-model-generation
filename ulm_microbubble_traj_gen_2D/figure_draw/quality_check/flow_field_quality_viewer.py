"""Interactive quality-control viewer for saved CFD flow-field results.

The viewer reads ``velocity_and_wall_shear_field.npz`` written by
``utils.workflows.runner`` and keeps the inspection focused on numerical flow
results.  Scalar fields can be switched without rebuilding the window, while
the continuous wall, anatomical openings, Cartesian cells, and sparse velocity
arrows remain available as independent overlays.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Mapping

import numpy as np


PACKAGE_DIR = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = PACKAGE_DIR / "results"
FIELD_FILENAME = "velocity_and_wall_shear_field.npz"
_QUIVER_CAPACITY = 900


@dataclass(frozen=True)
class FlowFieldQualityData:
    """Validated Cartesian flow fields and continuous vessel geometry."""

    source_path: Path
    x_coordinates_um: np.ndarray
    z_coordinates_um: np.ndarray
    spacing_um: float
    lumen_mask: np.ndarray
    boundary_mask: np.ndarray
    scalar_fields: Mapping[str, np.ndarray]
    vector_fields: Mapping[str, np.ndarray]
    continuous_wall_start_xz_um: np.ndarray
    continuous_wall_end_xz_um: np.ndarray
    open_section_point_xz_um: np.ndarray
    open_section_tangent_xz: np.ndarray
    open_section_half_width_um: np.ndarray
    open_section_kind: np.ndarray
    pressure_unit: str
    pressure_semantics: str

    @property
    def shape(self) -> tuple[int, int]:
        return (
            int(self.x_coordinates_um.size),
            int(self.z_coordinates_um.size),
        )

    @property
    def extent(self) -> tuple[float, float, float, float]:
        half = 0.5 * self.spacing_um
        return (
            float(self.x_coordinates_um[0] - half),
            float(self.x_coordinates_um[-1] + half),
            float(self.z_coordinates_um[0] - half),
            float(self.z_coordinates_um[-1] + half),
        )


@dataclass(frozen=True)
class _Layer:
    key: str
    label: str
    unit: str
    values: np.ndarray
    cmap: str
    vector_key: str = "final_velocity"
    mask_key: str = "lumen"
    symmetric: bool = False
    zero_floor: bool = False


def _resolve_field_archive(path: str | Path | None) -> Path:
    if path is not None:
        candidate = Path(path).expanduser().resolve()
        if candidate.is_dir():
            candidate = candidate / FIELD_FILENAME
        if not candidate.is_file():
            raise FileNotFoundError(f"Field archive does not exist: {candidate}")
        return candidate

    candidates = sorted(
        DEFAULT_RESULTS_DIR.glob(f"*/{FIELD_FILENAME}"),
        key=lambda item: item.stat().st_mtime_ns,
        reverse=True,
    )
    if not candidates:
        raise FileNotFoundError(
            f"No {FIELD_FILENAME} was found below {DEFAULT_RESULTS_DIR}."
        )
    return candidates[0].resolve()


def _required_array(
    archive: np.lib.npyio.NpzFile,
    key: str,
    *,
    ndim: int | None = None,
) -> np.ndarray:
    if key not in archive.files:
        raise ValueError(f"Field archive is missing required array {key!r}.")
    value = np.asarray(archive[key])
    if ndim is not None and value.ndim != ndim:
        raise ValueError(
            f"Field array {key!r} must have {ndim} dimensions, "
            f"found {value.shape}."
        )
    return value


def _required_grid(
    archive: np.lib.npyio.NpzFile,
    key: str,
    shape: tuple[int, int],
) -> np.ndarray:
    value = _required_array(archive, key)
    if value.shape != shape:
        raise ValueError(
            f"Field array {key!r} has shape {value.shape}; expected {shape}."
        )
    return np.ascontiguousarray(value)


def _required_vector_grid(
    archive: np.lib.npyio.NpzFile,
    key: str,
    shape: tuple[int, int],
) -> np.ndarray:
    value = _required_array(archive, key)
    expected = (*shape, 2)
    if value.shape != expected:
        raise ValueError(
            f"Field array {key!r} has shape {value.shape}; expected {expected}."
        )
    return np.ascontiguousarray(value)


def _optional_grid(
    archive: np.lib.npyio.NpzFile,
    key: str,
    shape: tuple[int, int],
) -> np.ndarray | None:
    if key not in archive.files:
        return None
    return _required_grid(archive, key, shape)


def _optional_vector_grid(
    archive: np.lib.npyio.NpzFile,
    key: str,
    shape: tuple[int, int],
) -> np.ndarray | None:
    if key not in archive.files:
        return None
    return _required_vector_grid(archive, key, shape)


def _optional_scalar_string(
    archive: np.lib.npyio.NpzFile,
    key: str,
    default: str,
) -> str:
    if key not in archive.files:
        return default
    values = np.asarray(archive[key]).reshape(-1)
    if values.size != 1:
        raise ValueError(f"Field array {key!r} must contain exactly one value.")
    return str(values[0])


def _optional_wall_geometry(
    archive: np.lib.npyio.NpzFile,
) -> tuple[np.ndarray, np.ndarray]:
    start_key = "continuous_wall_start_xz_um"
    end_key = "continuous_wall_end_xz_um"
    if start_key not in archive.files and end_key not in archive.files:
        empty = np.empty((0, 2), dtype=np.float64)
        return empty, empty.copy()
    starts = np.asarray(_required_array(archive, start_key, ndim=2), dtype=float)
    ends = np.asarray(_required_array(archive, end_key, ndim=2), dtype=float)
    if starts.shape != ends.shape or starts.shape[1:] != (2,):
        raise ValueError(
            "Continuous-wall start/end arrays must both have shape (N, 2)."
        )
    return np.ascontiguousarray(starts), np.ascontiguousarray(ends)


def _optional_open_sections(
    archive: np.lib.npyio.NpzFile,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    keys = (
        "continuous_open_section_point_xz_um",
        "continuous_open_section_tangent_xz",
        "continuous_open_section_half_width_um",
        "continuous_open_section_kind",
    )
    if not any(key in archive.files for key in keys):
        return (
            np.empty((0, 2), dtype=np.float64),
            np.empty((0, 2), dtype=np.float64),
            np.empty(0, dtype=np.float64),
            np.empty(0, dtype=np.int8),
        )
    if not all(key in archive.files for key in keys):
        missing = ", ".join(key for key in keys if key not in archive.files)
        raise ValueError(
            "Continuous open-section geometry is incomplete; missing "
            f"{missing}."
        )
    points = np.asarray(_required_array(archive, keys[0], ndim=2), dtype=float)
    tangents = np.asarray(_required_array(archive, keys[1], ndim=2), dtype=float)
    half_width = np.asarray(_required_array(archive, keys[2], ndim=1), dtype=float)
    kind = np.asarray(_required_array(archive, keys[3], ndim=1), dtype=np.int8)
    count = points.shape[0]
    if (
        points.shape[1:] != (2,)
        or tangents.shape != points.shape
        or half_width.shape != (count,)
        or kind.shape != (count,)
    ):
        raise ValueError(
            "Continuous open-section points/tangents must have shape (N, 2), "
            "and half-width/kind arrays must have shape (N,)."
        )
    return (
        np.ascontiguousarray(points),
        np.ascontiguousarray(tangents),
        np.ascontiguousarray(half_width),
        np.ascontiguousarray(kind),
    )


def load_flow_field_quality_data(
    path: str | Path | None = None,
) -> FlowFieldQualityData:
    """Load flow fields saved by :func:`utils.io.field_io.save_field_npz`."""

    field_path = _resolve_field_archive(path)
    with np.load(field_path, allow_pickle=False) as archive:
        x = np.asarray(
            _required_array(archive, "x_coordinates_um", ndim=1),
            dtype=np.float64,
        )
        z = np.asarray(
            _required_array(archive, "z_coordinates_um", ndim=1),
            dtype=np.float64,
        )
        if x.size < 2 or z.size < 2:
            raise ValueError("Flow viewing requires at least a 2 x 2 grid.")
        if (
            not np.all(np.isfinite(x))
            or not np.all(np.isfinite(z))
            or np.any(np.diff(x) <= 0.0)
            or np.any(np.diff(z) <= 0.0)
        ):
            raise ValueError(
                "Grid coordinates must be finite and strictly increasing."
            )
        shape = (int(x.size), int(z.size))

        spacing_values = np.asarray(
            _required_array(archive, "spacing_um"), dtype=np.float64
        ).reshape(-1)
        if spacing_values.size != 1:
            raise ValueError("spacing_um must contain exactly one value.")
        spacing = float(spacing_values[0])
        if not math.isfinite(spacing) or spacing <= 0.0:
            raise ValueError("spacing_um must be finite and positive.")
        tolerance = 1.0e-7 * max(spacing, 1.0)
        if not np.allclose(np.diff(x), spacing, rtol=0.0, atol=tolerance):
            raise ValueError("X coordinates are inconsistent with spacing_um.")
        if not np.allclose(np.diff(z), spacing, rtol=0.0, atol=tolerance):
            raise ValueError("Z coordinates are inconsistent with spacing_um.")

        lumen = np.asarray(
            _required_grid(archive, "lumen_mask", shape), dtype=bool
        )
        if not np.any(lumen):
            raise ValueError("lumen_mask does not contain any fluid cells.")

        final_velocity = np.asarray(
            _required_vector_grid(archive, "velocity_xz_um_s", shape),
            dtype=np.float32,
        )
        speed = np.asarray(
            _required_grid(archive, "speed_um_s", shape), dtype=np.float32
        )
        wall_shear = np.asarray(
            _required_grid(archive, "wall_shear_stress_pa", shape),
            dtype=np.float32,
        )

        scalar_fields: dict[str, np.ndarray] = {
            "speed_um_s": speed,
            "velocity_x_um_s": final_velocity[..., 0],
            "velocity_z_um_s": final_velocity[..., 1],
            "wall_shear_stress_pa": wall_shear,
        }
        vector_fields: dict[str, np.ndarray] = {
            "final_velocity": final_velocity,
        }

        optional_scalars = (
            "pressure",
            "divergence_s_inv",
            "wall_penetration_um_s",
            "open_boundary_flux_um2_s",
            "boundary_weight",
        )
        for key in optional_scalars:
            value = _optional_grid(archive, key, shape)
            if value is not None:
                scalar_fields[key] = np.ascontiguousarray(value)

        initial_velocity = _optional_vector_grid(
            archive, "initial_velocity_xz_um_s", shape
        )
        initial_speed = _optional_grid(archive, "initial_speed_um_s", shape)
        if (initial_velocity is None) != (initial_speed is None):
            raise ValueError(
                "initial_velocity_xz_um_s and initial_speed_um_s must either "
                "both be present or both be absent."
            )
        if initial_velocity is not None:
            initial_velocity = np.asarray(initial_velocity, dtype=np.float32)
            vector_fields["initial_velocity"] = initial_velocity
            scalar_fields["initial_velocity_x_um_s"] = initial_velocity[..., 0]
            scalar_fields["initial_velocity_z_um_s"] = initial_velocity[..., 1]
        if initial_speed is not None:
            initial_speed = np.asarray(initial_speed, dtype=np.float32)
            scalar_fields["initial_speed_um_s"] = initial_speed
            scalar_fields["speed_change_um_s"] = np.ascontiguousarray(
                speed - initial_speed
            )

        boundary_velocity = _optional_vector_grid(
            archive, "boundary_velocity_xz_um_s", shape
        )
        if boundary_velocity is not None:
            boundary_velocity = np.asarray(boundary_velocity, dtype=np.float32)
            vector_fields["boundary_velocity"] = boundary_velocity
            scalar_fields["boundary_speed_um_s"] = np.ascontiguousarray(
                np.linalg.norm(boundary_velocity, axis=2)
            )

        boundary_weight = scalar_fields.get("boundary_weight")
        if boundary_weight is None:
            boundary_mask = np.zeros(shape, dtype=bool)
            open_flux = scalar_fields.get("open_boundary_flux_um2_s")
            if open_flux is not None:
                boundary_mask = np.isfinite(open_flux) & (open_flux != 0.0)
            if boundary_velocity is not None:
                boundary_mask |= np.linalg.norm(boundary_velocity, axis=2) > 0.0
        else:
            boundary_mask = np.isfinite(boundary_weight) & (boundary_weight > 0.0)

        wall_starts, wall_ends = _optional_wall_geometry(archive)
        open_points, open_tangents, open_half_width, open_kind = (
            _optional_open_sections(archive)
        )
        pressure_unit = _optional_scalar_string(
            archive, "solver_sampled_pressure_unit", "saved units"
        )
        pressure_semantics = _optional_scalar_string(
            archive, "solver_sampled_pressure_semantics", "unspecified"
        )

    return FlowFieldQualityData(
        source_path=field_path,
        x_coordinates_um=np.ascontiguousarray(x),
        z_coordinates_um=np.ascontiguousarray(z),
        spacing_um=spacing,
        lumen_mask=np.ascontiguousarray(lumen),
        boundary_mask=np.ascontiguousarray(boundary_mask),
        scalar_fields=scalar_fields,
        vector_fields=vector_fields,
        continuous_wall_start_xz_um=wall_starts,
        continuous_wall_end_xz_um=wall_ends,
        open_section_point_xz_um=open_points,
        open_section_tangent_xz=open_tangents,
        open_section_half_width_um=open_half_width,
        open_section_kind=open_kind,
        pressure_unit=pressure_unit,
        pressure_semantics=pressure_semantics,
    )


def _build_layers(data: FlowFieldQualityData) -> list[_Layer]:
    fields = data.scalar_fields
    layers = [
        _Layer(
            "speed_um_s",
            "[Final] Speed",
            "um/s",
            fields["speed_um_s"],
            "viridis",
            zero_floor=True,
        ),
        _Layer(
            "velocity_x_um_s",
            "[Final] X velocity",
            "um/s",
            fields["velocity_x_um_s"],
            "coolwarm",
            symmetric=True,
        ),
        _Layer(
            "velocity_z_um_s",
            "[Final] Z velocity",
            "um/s",
            fields["velocity_z_um_s"],
            "coolwarm",
            symmetric=True,
        ),
    ]
    if "pressure" in fields:
        layers.append(
            _Layer(
                "pressure",
                "[Final] Pressure",
                data.pressure_unit,
                fields["pressure"],
                "coolwarm",
            )
        )
    layers.append(
        _Layer(
            "wall_shear_stress_pa",
            "[Final] Wall shear stress",
            "Pa",
            fields["wall_shear_stress_pa"],
            "magma",
            zero_floor=True,
        )
    )
    if "divergence_s_inv" in fields:
        layers.append(
            _Layer(
                "divergence_s_inv",
                "[QC] Divergence",
                "1/s",
                fields["divergence_s_inv"],
                "coolwarm",
                symmetric=True,
            )
        )
    if "wall_penetration_um_s" in fields:
        layers.append(
            _Layer(
                "wall_penetration_um_s",
                "[QC] Wall penetration",
                "um/s",
                fields["wall_penetration_um_s"],
                "coolwarm",
                symmetric=True,
            )
        )
    if "speed_change_um_s" in fields:
        layers.append(
            _Layer(
                "speed_change_um_s",
                "[QC] Final - initial speed",
                "um/s",
                fields["speed_change_um_s"],
                "coolwarm",
                symmetric=True,
            )
        )
    if "initial_speed_um_s" in fields:
        layers.extend(
            [
                _Layer(
                    "initial_speed_um_s",
                    "[Initial] Speed",
                    "um/s",
                    fields["initial_speed_um_s"],
                    "viridis",
                    vector_key="initial_velocity",
                    zero_floor=True,
                ),
                _Layer(
                    "initial_velocity_x_um_s",
                    "[Initial] X velocity",
                    "um/s",
                    fields["initial_velocity_x_um_s"],
                    "coolwarm",
                    vector_key="initial_velocity",
                    symmetric=True,
                ),
                _Layer(
                    "initial_velocity_z_um_s",
                    "[Initial] Z velocity",
                    "um/s",
                    fields["initial_velocity_z_um_s"],
                    "coolwarm",
                    vector_key="initial_velocity",
                    symmetric=True,
                ),
            ]
        )
    if "boundary_speed_um_s" in fields:
        layers.append(
            _Layer(
                "boundary_speed_um_s",
                "[Boundary] Prescribed speed",
                "um/s",
                fields["boundary_speed_um_s"],
                "viridis",
                vector_key="boundary_velocity",
                mask_key="boundary",
                zero_floor=True,
            )
        )
    if "open_boundary_flux_um2_s" in fields:
        layers.append(
            _Layer(
                "open_boundary_flux_um2_s",
                "[Boundary] Open-face flux",
                "um^2/s",
                fields["open_boundary_flux_um2_s"],
                "coolwarm",
                mask_key="boundary",
                symmetric=True,
            )
        )
    return layers


class FlowFieldQualityViewer:
    """Matplotlib viewer for saved flow fields and solver diagnostics."""

    def __init__(
        self,
        data: FlowFieldQualityData,
        *,
        robust_color_limits: bool = True,
    ):
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        from matplotlib.widgets import CheckButtons, RadioButtons

        self.data = data
        self.layers = _build_layers(data)
        self.layer_by_label = {layer.label: layer for layer in self.layers}
        self.current_layer = self.layers[0]
        self._show_wall = data.continuous_wall_start_xz_um.shape[0] > 0
        self._show_open_sections = data.open_section_point_xz_um.shape[0] > 0
        self._show_velocity_arrows = True
        self._show_cell_grid = False
        self._robust_color_limits = bool(robust_color_limits)
        self._updating_limits = False

        self.figure, self.axes = plt.subplots(figsize=(15.0, 9.0))
        self.figure.subplots_adjust(
            left=0.065,
            right=0.755,
            bottom=0.09,
            top=0.94,
        )
        self.radio_axes = self.figure.add_axes((0.78, 0.405, 0.20, 0.53))
        self.check_axes = self.figure.add_axes((0.78, 0.055, 0.20, 0.225))
        self.colorbar_axes = self.figure.add_axes((0.725, 0.17, 0.016, 0.62))

        self.image = self.axes.imshow(
            self._display_values(self.current_layer).T,
            origin="lower",
            extent=data.extent,
            interpolation="nearest",
            aspect="equal",
            cmap=self.current_layer.cmap,
            zorder=1,
        )
        self._apply_layer_limits()
        self.colorbar = self.figure.colorbar(self.image, cax=self.colorbar_axes)
        self._configure_colorbar()

        wall_segments = np.stack(
            (
                data.continuous_wall_start_xz_um,
                data.continuous_wall_end_xz_um,
            ),
            axis=1,
        )
        self.wall_collection = LineCollection(
            wall_segments,
            colors="#111111",
            linewidths=0.75,
            alpha=0.9,
            zorder=6,
        )
        self.wall_collection.set_visible(self._show_wall)
        self.axes.add_collection(self.wall_collection)

        open_segments, open_colors = self._open_section_segments_and_colors()
        self.open_section_collection = LineCollection(
            open_segments,
            colors=open_colors,
            linewidths=2.0,
            alpha=0.95,
            zorder=7,
        )
        self.open_section_collection.set_visible(self._show_open_sections)
        self.axes.add_collection(self.open_section_collection)

        self.grid_collection = LineCollection(
            [],
            colors="#555555",
            linewidths=0.35,
            alpha=0.48,
            zorder=5,
        )
        self.axes.add_collection(self.grid_collection)

        self.velocity_quiver = self.axes.quiver(
            np.zeros(_QUIVER_CAPACITY),
            np.zeros(_QUIVER_CAPACITY),
            np.ma.masked_all(_QUIVER_CAPACITY),
            np.ma.masked_all(_QUIVER_CAPACITY),
            angles="xy",
            scale_units="xy",
            scale=1.0,
            pivot="middle",
            color="#151515",
            width=0.0022,
            headwidth=3.4,
            headlength=4.2,
            alpha=0.78,
            zorder=8,
        )

        self.radio = RadioButtons(
            self.radio_axes,
            [layer.label for layer in self.layers],
            active=0,
            activecolor="#1f77b4",
        )
        for label in self.radio.labels:
            label.set_fontsize(7.8)
        self.radio.on_clicked(self._on_layer_selected)
        self.radio_axes.set_title("Flow-field result", fontsize=10, loc="left")

        self.checks = CheckButtons(
            self.check_axes,
            (
                "Continuous wall",
                "Anatomical openings",
                "Velocity arrows",
                "Cartesian cell edges",
                "Robust color limits",
            ),
            (
                self._show_wall,
                self._show_open_sections,
                self._show_velocity_arrows,
                self._show_cell_grid,
                self._robust_color_limits,
            ),
        )
        for label in self.checks.labels:
            label.set_fontsize(8.3)
        self.checks.on_clicked(self._on_option_toggled)

        self.summary_text = self.figure.text(
            0.78,
            0.375,
            "",
            ha="left",
            va="top",
            fontsize=8.0,
            family="monospace",
        )
        self.hover_text = self.figure.text(
            0.065,
            0.025,
            "Move over a lumen cell for exact flow values.",
            ha="left",
            va="center",
            fontsize=8.8,
        )

        self.axes.set_xlabel("X (um)")
        self.axes.set_ylabel("Z (um)")
        self.axes.set_xlim(data.extent[:2])
        self.axes.set_ylim(data.extent[2:])
        self.axes.callbacks.connect("xlim_changed", self._on_limits_changed)
        self.axes.callbacks.connect("ylim_changed", self._on_limits_changed)
        self.figure.canvas.mpl_connect("scroll_event", self._on_scroll)
        self.figure.canvas.mpl_connect("motion_notify_event", self._on_motion)
        self.figure.canvas.mpl_connect("key_press_event", self._on_key)

        self._update_title_and_summary()
        self._refresh_dynamic_overlays()

    def _mask_for_layer(self, layer: _Layer) -> np.ndarray:
        if layer.mask_key == "boundary":
            return self.data.boundary_mask
        return self.data.lumen_mask

    def _display_values(self, layer: _Layer) -> np.ma.MaskedArray:
        values = np.asarray(layer.values)
        valid = self._mask_for_layer(layer) & np.isfinite(values)
        return np.ma.array(values, mask=~valid, copy=False)

    def _finite_values(self, layer: _Layer) -> np.ndarray:
        return self._display_values(layer).compressed()

    def _color_limits(self, layer: _Layer) -> tuple[float, float]:
        values = self._finite_values(layer)
        if values.size == 0:
            return 0.0, 1.0
        finite = np.asarray(values, dtype=np.float64)
        if layer.symmetric:
            absolute = np.abs(finite)
            limit = float(
                np.percentile(absolute, 99.5)
                if self._robust_color_limits
                else np.max(absolute)
            )
            if not math.isfinite(limit) or limit <= np.finfo(float).eps:
                limit = max(float(np.max(absolute)), 1.0e-12)
            return -limit, limit

        if self._robust_color_limits and finite.size > 20:
            minimum, maximum = np.percentile(finite, (0.5, 99.5))
        else:
            minimum = float(np.min(finite))
            maximum = float(np.max(finite))
        minimum = float(minimum)
        maximum = float(maximum)
        if layer.zero_floor and minimum >= 0.0:
            minimum = 0.0
        if minimum == maximum:
            padding = max(abs(minimum) * 0.05, 1.0e-12)
            minimum -= padding
            maximum += padding
        return minimum, maximum

    def _apply_layer_limits(self) -> None:
        self.image.set_clim(*self._color_limits(self.current_layer))

    def _configure_colorbar(self) -> None:
        layer = self.current_layer
        self.colorbar.ax.yaxis.set_label_position("left")
        self.colorbar.set_label(
            layer.label if not layer.unit else f"{layer.label} ({layer.unit})"
        )
        from matplotlib.ticker import AutoLocator

        self.colorbar.locator = AutoLocator()
        self.colorbar.update_ticks()

    def _open_section_segments_and_colors(
        self,
    ) -> tuple[np.ndarray, list[str]]:
        points = self.data.open_section_point_xz_um
        if points.shape[0] == 0:
            return np.empty((0, 2, 2), dtype=float), []
        offsets = (
            self.data.open_section_tangent_xz
            * self.data.open_section_half_width_um[:, None]
        )
        segments = np.stack((points - offsets, points + offsets), axis=1)
        colors = [
            "#2364aa" if kind < 0 else "#f28e2b" if kind > 0 else "#666666"
            for kind in self.data.open_section_kind
        ]
        return segments, colors

    def _on_layer_selected(self, label: str) -> None:
        self.current_layer = self.layer_by_label[str(label)]
        self.image.set_data(self._display_values(self.current_layer).T)
        self.image.set_cmap(self.current_layer.cmap)
        self._apply_layer_limits()
        self.colorbar.update_normal(self.image)
        self._configure_colorbar()
        self._update_title_and_summary()
        self._refresh_dynamic_overlays()
        self.figure.canvas.draw_idle()

    def select_layer(self, key: str) -> None:
        """Select a scalar field by its stable key."""

        for index, layer in enumerate(self.layers):
            if layer.key == key:
                self.radio.set_active(index)
                return
        available = ", ".join(layer.key for layer in self.layers)
        raise KeyError(f"Unknown flow layer {key!r}; available: {available}")

    def _on_option_toggled(self, label: str) -> None:
        if label == "Continuous wall":
            self._show_wall = not self._show_wall
            self.wall_collection.set_visible(self._show_wall)
        elif label == "Anatomical openings":
            self._show_open_sections = not self._show_open_sections
            self.open_section_collection.set_visible(self._show_open_sections)
        elif label == "Velocity arrows":
            self._show_velocity_arrows = not self._show_velocity_arrows
        elif label == "Cartesian cell edges":
            self._show_cell_grid = not self._show_cell_grid
        elif label == "Robust color limits":
            self._robust_color_limits = not self._robust_color_limits
            self._apply_layer_limits()
            self._update_title_and_summary()
        self._refresh_dynamic_overlays()
        self.figure.canvas.draw_idle()

    def _on_limits_changed(self, _axes: object) -> None:
        if not self._updating_limits:
            self._refresh_dynamic_overlays()

    def _visible_index_bounds(self) -> tuple[int, int, int, int]:
        xmin, xmax = sorted(self.axes.get_xlim())
        zmin, zmax = sorted(self.axes.get_ylim())
        x = self.data.x_coordinates_um
        z = self.data.z_coordinates_um
        half = 0.5 * self.data.spacing_um
        ix0 = max(0, int(np.searchsorted(x + half, xmin, side="left")))
        ix1 = min(x.size, int(np.searchsorted(x - half, xmax, side="right")))
        iz0 = max(0, int(np.searchsorted(z + half, zmin, side="left")))
        iz1 = min(z.size, int(np.searchsorted(z - half, zmax, side="right")))
        return ix0, max(ix0, ix1), iz0, max(iz0, iz1)

    def _refresh_dynamic_overlays(self) -> None:
        ix0, ix1, iz0, iz1 = self._visible_index_bounds()
        nx = ix1 - ix0
        nz = iz1 - iz0
        self._refresh_cell_grid(ix0, ix1, iz0, iz1)
        self._refresh_velocity_arrows(ix0, ix1, iz0, iz1, nx, nz)

    def _refresh_cell_grid(
        self,
        ix0: int,
        ix1: int,
        iz0: int,
        iz1: int,
    ) -> None:
        nx = ix1 - ix0
        nz = iz1 - iz0
        segments = []
        if self._show_cell_grid and nx > 0 and nz > 0 and nx * nz <= 12_000:
            h = self.data.spacing_um
            x_edges = self.data.x_coordinates_um[ix0:ix1] - 0.5 * h
            z_edges = self.data.z_coordinates_um[iz0:iz1] - 0.5 * h
            x_edges = np.append(
                x_edges, self.data.x_coordinates_um[ix1 - 1] + 0.5 * h
            )
            z_edges = np.append(
                z_edges, self.data.z_coordinates_um[iz1 - 1] + 0.5 * h
            )
            segments.extend(
                ((float(value), float(z_edges[0])), (float(value), float(z_edges[-1])))
                for value in x_edges
            )
            segments.extend(
                ((float(x_edges[0]), float(value)), (float(x_edges[-1]), float(value)))
                for value in z_edges
            )
        self.grid_collection.set_segments(segments)

    def _refresh_velocity_arrows(
        self,
        ix0: int,
        ix1: int,
        iz0: int,
        iz1: int,
        nx: int,
        nz: int,
    ) -> None:
        if not self._show_velocity_arrows or nx <= 0 or nz <= 0:
            self._set_empty_quiver()
            return
        vector = self.data.vector_fields.get(self.current_layer.vector_key)
        if vector is None:
            self._set_empty_quiver()
            return

        step_x = max(1, int(math.ceil(nx / 30.0)))
        step_z = max(1, int(math.ceil(nz / 30.0)))
        selected = vector[ix0:ix1:step_x, iz0:iz1:step_z]
        mask = self._mask_for_layer(self.current_layer)[
            ix0:ix1:step_x, iz0:iz1:step_z
        ]
        speed = np.linalg.norm(selected, axis=2)
        valid = mask & np.all(np.isfinite(selected), axis=2) & (speed > 0.0)
        if not np.any(valid):
            self._set_empty_quiver()
            return

        x = self.data.x_coordinates_um[ix0:ix1:step_x]
        z = self.data.z_coordinates_um[iz0:iz1:step_z]
        x_grid, z_grid = np.meshgrid(x, z, indexing="ij")
        target_length = (
            0.62 * self.data.spacing_um * max(step_x, step_z)
        )
        unit = np.divide(
            selected,
            speed[..., None],
            out=np.zeros_like(selected, dtype=np.float64),
            where=speed[..., None] > 0.0,
        )
        offsets = np.column_stack((x_grid[valid], z_grid[valid]))
        count = offsets.shape[0]
        padded_offsets = np.zeros((_QUIVER_CAPACITY, 2), dtype=np.float64)
        padded_offsets[:count] = offsets
        u_values = np.zeros(_QUIVER_CAPACITY, dtype=np.float64)
        v_values = np.zeros(_QUIVER_CAPACITY, dtype=np.float64)
        u_values[:count] = unit[..., 0][valid] * target_length
        v_values[:count] = unit[..., 1][valid] * target_length
        hidden = np.ones(_QUIVER_CAPACITY, dtype=bool)
        hidden[:count] = False
        self.velocity_quiver.set_offsets(padded_offsets)
        self.velocity_quiver.set_UVC(
            np.ma.array(u_values, mask=hidden),
            np.ma.array(v_values, mask=hidden),
        )
        self.velocity_quiver.set_visible(True)

    def _set_empty_quiver(self) -> None:
        self.velocity_quiver.set_UVC(
            np.ma.masked_all(_QUIVER_CAPACITY),
            np.ma.masked_all(_QUIVER_CAPACITY),
        )
        self.velocity_quiver.set_visible(False)

    def _on_scroll(self, event: object) -> None:
        if getattr(event, "inaxes", None) is not self.axes:
            return
        xdata = getattr(event, "xdata", None)
        zdata = getattr(event, "ydata", None)
        if xdata is None or zdata is None:
            return
        factor = 0.80 if getattr(event, "button", None) == "up" else 1.25
        xmin, xmax = self.axes.get_xlim()
        zmin, zmax = self.axes.get_ylim()
        width = (xmax - xmin) * factor
        height = (zmax - zmin) * factor
        x_fraction = (float(xdata) - xmin) / max(
            xmax - xmin, np.finfo(float).eps
        )
        z_fraction = (float(zdata) - zmin) / max(
            zmax - zmin, np.finfo(float).eps
        )
        self._updating_limits = True
        self.axes.set_xlim(
            float(xdata) - x_fraction * width,
            float(xdata) + (1.0 - x_fraction) * width,
        )
        self.axes.set_ylim(
            float(zdata) - z_fraction * height,
            float(zdata) + (1.0 - z_fraction) * height,
        )
        self._updating_limits = False
        self._refresh_dynamic_overlays()
        self.figure.canvas.draw_idle()

    def _on_motion(self, event: object) -> None:
        if getattr(event, "inaxes", None) is not self.axes:
            return
        x_value = getattr(event, "xdata", None)
        z_value = getattr(event, "ydata", None)
        if x_value is None or z_value is None:
            return
        i = int(
            np.floor((float(x_value) - self.data.extent[0]) / self.data.spacing_um)
        )
        j = int(
            np.floor((float(z_value) - self.data.extent[2]) / self.data.spacing_um)
        )
        if not (0 <= i < self.data.shape[0] and 0 <= j < self.data.shape[1]):
            return
        if not self.data.lumen_mask[i, j]:
            self.hover_text.set_text(
                f"cell=({i}, {j}) | centre="
                f"({self.data.x_coordinates_um[i]:.6g}, "
                f"{self.data.z_coordinates_um[j]:.6g}) um | outside lumen"
            )
            self.figure.canvas.draw_idle()
            return

        fields = self.data.scalar_fields
        velocity = self.data.vector_fields["final_velocity"][i, j]
        details = [
            f"cell=({i}, {j})",
            (
                f"centre=({self.data.x_coordinates_um[i]:.6g}, "
                f"{self.data.z_coordinates_um[j]:.6g}) um"
            ),
            (
                f"current={float(self.current_layer.values[i, j]):.6g} "
                f"{self.current_layer.unit}"
            ).rstrip(),
            f"u=({velocity[0]:.6g}, {velocity[1]:.6g}) um/s",
            f"speed={float(fields['speed_um_s'][i, j]):.6g} um/s",
            (
                f"WSS={float(fields['wall_shear_stress_pa'][i, j]):.6g} Pa"
            ),
        ]
        if "pressure" in fields:
            details.append(
                f"p={float(fields['pressure'][i, j]):.6g} "
                f"{self.data.pressure_unit}"
            )
        if "divergence_s_inv" in fields:
            details.append(
                f"div={float(fields['divergence_s_inv'][i, j]):.6g} 1/s"
            )
        self.hover_text.set_text(" | ".join(details))
        self.figure.canvas.draw_idle()

    def _on_key(self, event: object) -> None:
        key = str(getattr(event, "key", "")).lower()
        if key in {"r", "home"}:
            self.axes.set_xlim(self.data.extent[:2])
            self.axes.set_ylim(self.data.extent[2:])
            self.figure.canvas.draw_idle()
        elif key == "v":
            self.checks.set_active(2)
        elif key == "g":
            self.checks.set_active(3)
        elif key == "w":
            self.checks.set_active(0)
        elif key == "o":
            self.checks.set_active(1)
        elif key == "c":
            self.checks.set_active(4)

    def _update_title_and_summary(self) -> None:
        mode = "robust" if self._robust_color_limits else "full"
        self.axes.set_title(
            f"{self.current_layer.label} - {self._source_label()} "
            f"({mode} color range)",
            loc="left",
        )
        values = self._finite_values(self.current_layer)
        if values.size:
            value_summary = (
                f"min/median/max: {np.min(values):.5g} / "
                f"{np.median(values):.5g} / {np.max(values):.5g}"
            )
        else:
            value_summary = "min/median/max: no finite values"
        lower, upper = self.image.get_clim()
        self.summary_text.set_text(
            f"grid: {self.data.shape[0]} x {self.data.shape[1]}\n"
            f"spacing: {self.data.spacing_um:.6g} um\n"
            f"{value_summary}\n"
            f"color: {lower:.5g} .. {upper:.5g}"
        )

    def _source_label(self) -> str:
        if self.data.source_path.name == FIELD_FILENAME:
            return self.data.source_path.parent.name
        return self.data.source_path.stem

    def save_snapshot(self, path: str | Path, *, dpi: int = 180) -> Path:
        """Save the current viewport as PNG, PDF, or SVG."""

        output_path = Path(path).expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        self.figure.savefig(output_path, dpi=int(dpi), bbox_inches="tight")
        return output_path

    def show(self) -> None:
        import matplotlib.pyplot as plt

        plt.show()


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Inspect velocity, pressure, wall shear, and solver diagnostics "
            "saved by the main flow workflow."
        )
    )
    parser.add_argument(
        "--field",
        type=Path,
        default=None,
        help=(
            "A velocity_and_wall_shear_field.npz or its result directory. "
            "The newest result is used when omitted."
        ),
    )
    parser.add_argument(
        "--layer",
        default="speed_um_s",
        help="Initial stable flow-layer key.",
    )
    parser.add_argument(
        "--full-color-range",
        action="store_true",
        help="Use exact extrema instead of robust percentile color limits.",
    )
    parser.add_argument(
        "--snapshot",
        type=Path,
        default=None,
        help="Optionally save the initial/current view to PNG, PDF, or SVG.",
    )
    parser.add_argument(
        "--no-show",
        action="store_true",
        help="Build the viewer without opening a GUI.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    data = load_flow_field_quality_data(args.field)
    viewer = FlowFieldQualityViewer(
        data,
        robust_color_limits=not bool(args.full_color_range),
    )
    viewer.select_layer(str(args.layer))
    if args.snapshot is not None:
        saved = viewer.save_snapshot(args.snapshot)
        print(f"Saved flow-field snapshot: {saved}")
    print(f"Flow field archive: {data.source_path}")
    print(
        "Controls: mouse wheel=zoom, toolbar=pan/home, "
        "V=velocity arrows, G=Cartesian edges, W=wall, "
        "O=openings, C=color range."
    )
    if not args.no_show:
        viewer.show()


if __name__ == "__main__":
    main()
