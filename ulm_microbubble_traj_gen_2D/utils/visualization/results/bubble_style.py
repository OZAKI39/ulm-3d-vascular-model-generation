"""Pure visual-style helpers for physically sized microbubble markers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


DEFAULT_GLOW_DIAMETER_SCALE = 3.0
DEFAULT_GLOW_OUTER_DIAMETER_POINTS = 6.0


@dataclass(frozen=True)
class BubbleGlowStyle:
    """Describe a restrained three-layer halo around the physical bubble disk."""

    enabled: bool = True
    outer_diameter_scale: float = DEFAULT_GLOW_DIAMETER_SCALE
    # 三层透明圆盘由外向内叠加，得到近似径向渐隐的光晕，而不是新的实体边界。
    layer_positions: tuple[float, float, float] = (1.0, 0.60, 0.25)
    layer_alphas: tuple[float, float, float] = (0.08, 0.14, 0.25)
    # 全景图中真实微泡可能不足一个像素；该下限只约束光晕，不会放大实体圆盘。
    minimum_outer_diameter_points: float = DEFAULT_GLOW_OUTER_DIAMETER_POINTS
    glow_rgb: tuple[float, float, float] = (0.35, 0.94, 1.0)
    highlight_diameter_scale: float = 0.28
    highlight_offset_fraction: float = 0.15
    highlight_alpha: float = 0.46

    def __post_init__(self) -> None:
        scale = float(self.outer_diameter_scale)
        if not np.isfinite(scale) or scale < 1.0:
            raise ValueError("glow diameter scale must be finite and at least 1.0")
        if len(self.layer_positions) != len(self.layer_alphas):
            raise ValueError("glow layer positions and alphas must have equal length")
        if any(not 0.0 <= float(value) <= 1.0 for value in self.layer_positions):
            raise ValueError("glow layer positions must lie within [0, 1]")
        if any(not 0.0 <= float(value) <= 1.0 for value in self.layer_alphas):
            raise ValueError("glow layer alphas must lie within [0, 1]")
        if not np.isfinite(self.minimum_outer_diameter_points) or self.minimum_outer_diameter_points <= 0.0:
            raise ValueError("minimum outer glow diameter must be finite and positive")
        if len(self.glow_rgb) != 3 or any(not 0.0 <= float(value) <= 1.0 for value in self.glow_rgb):
            raise ValueError("glow RGB components must lie within [0, 1]")

    @property
    def layer_diameter_scales(self) -> tuple[float, ...]:
        """Return outer-to-inner halo scales while keeping every layer outside the real disk."""

        excess = float(self.outer_diameter_scale) - 1.0
        return tuple(1.0 + excess * float(position) for position in self.layer_positions)


def build_bubble_glow_style(*, enabled: bool, outer_diameter_scale: float) -> BubbleGlowStyle:
    """Build a style whose physical and overview-screen glow sizes respond to one scale option."""

    scale = float(outer_diameter_scale)
    screen_size = DEFAULT_GLOW_OUTER_DIAMETER_POINTS * scale / DEFAULT_GLOW_DIAMETER_SCALE
    return BubbleGlowStyle(
        enabled=bool(enabled),
        outer_diameter_scale=scale,
        minimum_outer_diameter_points=screen_size,
    )


def glow_layer_diameters(diameters_um: np.ndarray, style: BubbleGlowStyle) -> tuple[np.ndarray, ...]:
    """Return halo diameters without altering the real physical bubble diameters."""

    physical = np.asarray(diameters_um, dtype=float)
    if not style.enabled:
        return ()
    return tuple(physical * scale for scale in style.layer_diameter_scales)


def glow_layer_diameters_points(
    diameters_um: np.ndarray,
    points_per_um: float,
    style: BubbleGlowStyle,
) -> tuple[np.ndarray, ...]:
    """Return screen-visible halo sizes while preserving the requested physical scale when larger."""

    if not style.enabled:
        return ()
    physical_layers_um = glow_layer_diameters(diameters_um, style)
    scales = style.layer_diameter_scales
    outer_scale = max(float(scales[0]), np.finfo(float).eps)
    result: list[np.ndarray] = []
    for physical_um, scale in zip(physical_layers_um, scales, strict=True):
        visible_floor = float(style.minimum_outer_diameter_points) * float(scale) / outer_scale
        result.append(np.maximum(physical_um * float(points_per_um), visible_floor))
    return tuple(result)


def bubble_highlight_geometry(
    xy_um: np.ndarray,
    diameters_um: np.ndarray,
    style: BubbleGlowStyle,
    rotation_angle_rad: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Place a rotating diagnostic highlight fully inside each physical bubble disk."""

    xy = np.asarray(xy_um, dtype=float)
    physical = np.asarray(diameters_um, dtype=float)
    offsets = np.array(xy, copy=True)
    if offsets.size:
        shift = physical * float(style.highlight_offset_fraction)
        angle = (
            np.zeros(physical.shape, dtype=float)
            if rotation_angle_rad is None
            else np.asarray(rotation_angle_rad, dtype=float)
        )
        if angle.shape != physical.shape:
            raise ValueError("rotation_angle_rad must match diameters_um.shape")
        angle = np.where(np.isfinite(angle), angle, 0.0)
        base_x = -shift
        base_z = shift
        cosine = np.cos(angle)
        sine = np.sin(angle)
        # Positive +Y rotation is clockwise when the X-Z plane is viewed from +Y.
        offsets[:, 0] += cosine * base_x + sine * base_z
        offsets[:, 1] += -sine * base_x + cosine * base_z
    return offsets, physical * float(style.highlight_diameter_scale)
