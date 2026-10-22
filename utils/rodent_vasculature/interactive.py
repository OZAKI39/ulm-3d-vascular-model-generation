"""Interactive Figure 2(a)-style rendering for mouse-brain vascular blocks.

The preprocessed main-network volume is rendered as white vessels on black,
matching the visual language of Figure 2(a). SWC centerlines, critical nodes, and arrows
are overlays; arrow direction is the annotation relationship
``parent_id node -> current node`` and is not a measured flow direction.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

import numpy as np

from .model import DirectedVascularGraph
from .roi_surface import build_roi_display_tubes
from .tiff_io import load_normalized_volume, load_tiff_volume
from .visualization import sample_direction_arrows
from ..sampling.sampling_io import load_sampling_display_rois
from ..sampling.sampling_types import ROIRecord


UI_FONT_FAMILY = "arial"
# 主坐标轴数字（数值刻度）的字号。
COORDINATE_TICK_FONT_SIZE = 11
# 主坐标轴标题（X/Y/Z (um)）的字号。
COORDINATE_TITLE_FONT_SIZE = 12
# 左下角方向指示器中 X/Y/Z 标题的字号。
ORIENTATION_AXIS_FONT_SIZE = 12
LEGEND_FONT_SIZE = 11

# 底部科研图例：左下角 0~0.2 的归一化视口留给方向坐标系，图例从其右侧开始。
LEGEND_LEFT_FRACTION = 0.205
LEGEND_RIGHT_MARGIN_FRACTION = 0.015
LEGEND_BOTTOM_FRACTION = 0.012
# 图例整体高度和相邻两行中心之间的距离。减小 ROW_SPACING 即可单独压缩行距；
# HEIGHT 只控制底部深色背景区域的总高度。
LEGEND_HEIGHT_FRACTION = 0.105
LEGEND_ROW_SPACING_FRACTION = 0.032
LEGEND_COLUMN_COUNT = 3
LEGEND_SWATCH_WIDTH_FRACTION = 0.035
LEGEND_SWATCH_HEIGHT_FRACTION = 0.010
LEGEND_LINE_WIDTH_PX = 4.0
LEGEND_ARROW_HEAD_LENGTH_FRACTION = 0.010
LEGEND_ARROW_HEAD_HALF_HEIGHT_FRACTION = 0.006
LEGEND_CIRCLE_RADIUS_FRACTION = 0.008
LEGEND_CIRCLE_RESOLUTION = 40
LEGEND_TEXT_GAP_FRACTION = 0.010
LEGEND_BACKGROUND_COLOR = "#101010"
LEGEND_BACKGROUND_OPACITY = 0.78
LEGEND_TEXT_COLOR = "#F4F4F4"

# 右侧局部 ROI 血管使用直径着色。色条坐标均相对于右侧子图自身，因而窗口
# 缩放后仍固定在右侧；色彩范围在所有候选 ROI 间统一，切换 ROI 时不会跳变。
DIAMETER_SCALAR_NAME = "diameter_um"
DIAMETER_COLORMAP = "viridis"
# False：色条跟随当前右图 ROI，能充分展示该 ROI 内部的直径差异；
# True：所有 ROI 共用全局最小/最大值，适合切换 ROI 时进行绝对颜色比较。
DIAMETER_COLORBAR_USE_GLOBAL_RANGE = False
DIAMETER_COLORBAR_TITLE = "Diameter (um)"
DIAMETER_COLORBAR_POSITION_X = 0.895
DIAMETER_COLORBAR_POSITION_Y = 0.245
DIAMETER_COLORBAR_WIDTH_FRACTION = 0.075
DIAMETER_COLORBAR_HEIGHT_FRACTION = 0.48
DIAMETER_COLORBAR_LABEL_COUNT = 5
DIAMETER_COLORBAR_TITLE_FONT_SIZE = 12
DIAMETER_COLORBAR_LABEL_FONT_SIZE = 10
DIAMETER_COLORBAR_FORMAT = "%.1f"
DIAMETER_COLORBAR_TITLE_SEPARATION_PX = 8

# 两个视口内三维标记的独立尺寸。左图的四类关键节点、右图的真实终端和
# 截断端口均为圆球，因此集中放在这里，便于继续统一缩放。
STRUCTURAL_ROOT_POINT_SIZE_PX = 9.0
STRUCTURAL_LEAF_POINT_SIZE_PX = 7.0
DIVERGENCE_JUNCTION_POINT_SIZE_PX = 9.0
CONVERGENCE_JUNCTION_POINT_SIZE_PX = 9.0
RIGHT_TRUE_TERMINAL_POINT_SIZE_PX = 14.0
RIGHT_CUT_PORT_POINT_SIZE_PX = 15.0

# 左、右方向箭头的尺度必须相互独立。右图保持原尺寸，只限制数量并按空间
# 均匀抽样；左图提高到包围盒对角线的 3%，使全局方向更醒目。
LEFT_PARENT_ARROW_MIN_LENGTH_UM = 2.0
LEFT_PARENT_ARROW_SCALE_FRACTION = 0.030
RIGHT_PARENT_ARROW_MIN_LENGTH_UM = 2.0
RIGHT_PARENT_ARROW_SCALE_FRACTION = 0.055
RIGHT_PARENT_ARROW_MAX_COUNT = 24

# 左右视口共用同一个 ROI 边界红色，保证活动框和局部 ROI 一一对应。
ROI_BOUNDARY_COLOR = "#FF5C5C"
# 左图当前活动 ROI 的灰色半透明填充；其他候选 ROI 继续使用聚类颜色。
ACTIVE_ROI_FILL_COLOR = "#A7A7A7"
ACTIVE_ROI_FILL_OPACITY = 0.34
ACTIVE_ROI_EDGE_COLOR = ROI_BOUNDARY_COLOR
ACTIVE_ROI_EDGE_WIDTH_PX = 5.0
CANDIDATE_ROI_LEGEND_COLOR = "#4DD0E1"

# 主三维坐标框的屏幕布局参数。刻度线、数字与标题均贴在三维轴的投影位置，
# 刻度方向严格垂直于轴并朝外；候选边之间仍使用连续淡入淡出，不会瞬间换边。
COORDINATE_AXIS_COLOR = "#BFC7D5"
# 左图坐标框线宽；框线与刻度共享完全相同的 bounds 几何。
COORDINATE_FRAME_LINE_WIDTH_PX = 1.2
# 左图六个框面的辅助网格样式；网格同样由显式 bounds 几何生成。
COORDINATE_GRID_LINE_WIDTH_PX = 1.0
COORDINATE_GRID_OPACITY = 0.28
# 主刻度短线的长度和粗细，单位都是屏幕像素。
COORDINATE_TICK_LENGTH_PX = 7.0
COORDINATE_TICK_LINE_WIDTH_PX = 1.5
# 数字锚点到框线、标题锚点到数字锚点的外向距离（屏幕像素）。
COORDINATE_TICK_LABEL_OFFSET_PX = 10.0
COORDINATE_TITLE_OFFSET_PX = 30.0
# X/Y/Z 每根活动轴上显示的主刻度数量（包含首尾刻度）。
COORDINATE_X_TICK_COUNT = 4
COORDINATE_Y_TICK_COUNT = 4
COORDINATE_Z_TICK_COUNT = 5
COORDINATE_EDGE_CROSSFADE_PX = 36.0
COORDINATE_TEXT_VISIBILITY_EPSILON = 0.01

# 交互窗口中的自动水平旋转参数。定时器只在 show=True 的实时 GUI 中安装，
# 因此离线验收截图不会在保存期间发生转动。
AUTO_HORIZONTAL_ROTATION_ENABLED = True
HORIZONTAL_ROTATION_INTERVAL_MS = 40
HORIZONTAL_ROTATION_DEGREES_PER_STEP = 0.40
HORIZONTAL_ROTATION_VIEW_UP = (0.0, 0.0, 1.0)

# PyVista 的 timer event 必须给出最大步数。这里约等于连续运行 994 天，实际
# 定时器会在用户关闭窗口时随 VTK interactor 一起结束。
HORIZONTAL_ROTATION_MAX_STEPS = 2_147_483_647


LegendSymbol = Literal["line", "arrow", "circle", "patch"]


@dataclass(frozen=True, slots=True)
class _ScientificLegendEntry:
    """One legend entry with a symbol matching the rendered object type."""

    label: str
    color: str
    symbol: LegendSymbol


@dataclass(frozen=True, slots=True)
class InteractiveSceneGeometry:
    """Small, renderer-independent subset of a directed vascular graph."""

    branches_um: tuple[np.ndarray, ...]
    arrow_points_um: np.ndarray
    arrow_vectors_xyz: np.ndarray
    critical_points_um: np.ndarray
    critical_roles: tuple[str, ...]

    @property
    def branch_count(self) -> int:
        return len(self.branches_um)

    @property
    def arrow_count(self) -> int:
        return len(self.arrow_points_um)


@dataclass(frozen=True, slots=True)
class Figure2aArtifacts:
    screenshot_path: Path
    manifest_path: Path


def _continuous_edge_opacities(
    scores: np.ndarray,
    *,
    prefer_minimum: bool,
) -> np.ndarray:
    """Return smooth edge visibility weights without discrete winner switching."""

    values = np.asarray(scores, dtype=float)
    if values.ndim != 1 or values.size == 0 or not np.all(np.isfinite(values)):
        raise ValueError("Projected edge scores must be a non-empty finite vector")
    if COORDINATE_EDGE_CROSSFADE_PX <= 0:
        raise ValueError("COORDINATE_EDGE_CROSSFADE_PX must be positive")
    extreme = float(np.min(values) if prefer_minimum else np.max(values))
    distance = values - extreme if prefer_minimum else extreme - values
    proximity = np.clip(
        1.0 - distance / COORDINATE_EDGE_CROSSFADE_PX,
        0.0,
        1.0,
    )
    # Smoothstep gives zero slope at both ends, so an annotation neither pops in
    # nor visibly changes speed when it enters or leaves the transition interval.
    return proximity * proximity * (3.0 - 2.0 * proximity)


@dataclass(slots=True)
class _ProjectedAxisAnnotation:
    """A fixed world-space axis whose labels are rendered in screen space."""

    axis_name: str
    point1_xyz: tuple[float, float, float]
    point2_xyz: tuple[float, float, float]
    tick_line_points: Any
    tick_line_actor: Any
    tick_actors: tuple[Any, ...]
    title_actor: Any
    last_outward_xy: tuple[float, float] | None = None
    tick_label_offset_px: float = COORDINATE_TICK_LABEL_OFFSET_PX
    title_offset_px: float = COORDINATE_TITLE_OFFSET_PX
    vertical_axis_name: str = "Z"
    endpoint_label_inset_px: float = 0.0

    def _set_outside_text_alignment(
        self,
        actor: Any,
        outward: np.ndarray,
    ) -> None:
        """Keep horizontal text wholly on the outside half of its axis."""

        text_property = actor.GetTextProperty()
        # 标题和数字始终保持屏幕正向，彻底避免旋转半周后出现倒字。
        text_property.SetOrientation(0.0)
        if self.axis_name != self.vertical_axis_name:
            # 底部活动边以锚点作为文字上边缘，文字只向下展开。
            text_property.SetJustificationToCentered()
            text_property.SetVerticalJustificationToTop()
        else:
            # 竖直轴在左右两侧显示，以靠轴的一侧对齐，文字向画框外展开。
            if float(outward[0]) < 0.0:
                text_property.SetJustificationToRight()
            else:
                text_property.SetJustificationToLeft()
            text_property.SetVerticalJustificationToCentered()

    def update_from_projection(
        self,
        point1_xy: np.ndarray,
        point2_xy: np.ndarray,
        box_center_xy: np.ndarray,
        opacity: float,
    ) -> None:
        """Place every annotation next to its axis and away from the box centre."""

        axis_vector = np.asarray(point2_xy, dtype=float) - np.asarray(point1_xy, dtype=float)
        axis_length = float(np.linalg.norm(axis_vector))
        clamped_opacity = float(np.clip(opacity, 0.0, 1.0))
        visible = (
            axis_length > 1.0e-9
            and clamped_opacity > COORDINATE_TEXT_VISIBILITY_EPSILON
        )
        for actor in (*self.tick_actors, self.title_actor):
            actor.SetVisibility(visible)
            actor.GetTextProperty().SetOpacity(clamped_opacity)
        self.tick_line_actor.SetVisibility(visible)
        self.tick_line_actor.GetProperty().SetOpacity(clamped_opacity)
        if axis_length <= 1.0e-9:
            return

        midpoint = (np.asarray(point1_xy, dtype=float) + np.asarray(point2_xy, dtype=float)) * 0.5
        axis_direction = axis_vector / axis_length
        outward = np.asarray((-axis_direction[1], axis_direction[0]), dtype=float)
        centre_to_axis = midpoint - np.asarray(box_center_xy, dtype=float)
        side_score = float(np.dot(outward, centre_to_axis))
        if side_score < 0.0:
            outward *= -1.0
        elif abs(side_score) <= 1.0e-9 and self.last_outward_xy is not None:
            # 投影恰好退化到中心线时延续上一帧方向，避免数值噪声令刻度来回翻转。
            previous = np.asarray(self.last_outward_xy, dtype=float)
            if float(np.dot(outward, previous)) < 0.0:
                outward *= -1.0
        self.last_outward_xy = (float(outward[0]), float(outward[1]))

        fractions = np.linspace(0.0, 1.0, len(self.tick_actors))
        for tick_index, (fraction, actor) in enumerate(zip(fractions, self.tick_actors)):
            axis_point = np.asarray(point1_xy, dtype=float) + fraction * axis_vector
            tick_end = axis_point + outward * COORDINATE_TICK_LENGTH_PX
            self.tick_line_points.SetPoint(
                tick_index * 2,
                float(axis_point[0]),
                float(axis_point[1]),
                0.0,
            )
            self.tick_line_points.SetPoint(
                tick_index * 2 + 1,
                float(tick_end[0]),
                float(tick_end[1]),
                0.0,
            )
            label_position = axis_point + outward * self.tick_label_offset_px
            if tick_index == 0:
                label_position += axis_direction * self.endpoint_label_inset_px
            elif tick_index == len(self.tick_actors) - 1:
                label_position -= axis_direction * self.endpoint_label_inset_px
            actor.SetPosition(float(label_position[0]), float(label_position[1]))
            self._set_outside_text_alignment(actor, outward)
        self.tick_line_points.Modified()

        title_position = midpoint + outward * (
            self.tick_label_offset_px + self.title_offset_px
        )
        self.title_actor.SetPosition(
            float(title_position[0]),
            float(title_position[1]),
        )
        self._set_outside_text_alignment(self.title_actor, outward)


@dataclass(slots=True)
class _ScreenAnchoredCoordinateAxes:
    """Continuously keep axis annotations attached and facing outside the box."""

    renderer: Any
    camera: Any
    bounds: tuple[float, float, float, float, float, float]
    x_annotations: tuple[_ProjectedAxisAnnotation, ...]
    y_annotations: tuple[_ProjectedAxisAnnotation, ...]
    z_annotations: tuple[_ProjectedAxisAnnotation, ...]
    y_vertical: bool = False
    camera_observer_id: int | None = None
    x_opacities: tuple[float, ...] = ()
    y_opacities: tuple[float, ...] = ()
    z_opacities: tuple[float, ...] = ()
    single_edge_per_axis: bool = False
    _updating: bool = False

    def _project_to_display(self, point_xyz: tuple[float, float, float]) -> np.ndarray:
        self.renderer.SetWorldPoint(
            float(point_xyz[0]),
            float(point_xyz[1]),
            float(point_xyz[2]),
            1.0,
        )
        self.renderer.WorldToDisplay()
        display_point = self.renderer.GetDisplayPoint()
        return np.asarray(
            (
                float(display_point[0]),
                float(display_point[1]),
            ),
            dtype=float,
        )

    def _project_annotations(
        self,
        annotations: tuple[_ProjectedAxisAnnotation, ...],
    ) -> list[tuple[np.ndarray, np.ndarray]]:
        return [
            (
                self._project_to_display(annotation.point1_xyz),
                self._project_to_display(annotation.point2_xyz),
            )
            for annotation in annotations
        ]

    @staticmethod
    def _midpoint_scores(
        projected_axes: list[tuple[np.ndarray, np.ndarray]],
        coordinate: int,
    ) -> np.ndarray:
        return np.asarray(
            [
                float((point1[coordinate] + point2[coordinate]) * 0.5)
                for point1, point2 in projected_axes
            ],
            dtype=float,
        )

    @staticmethod
    def _update_group(
        annotations: tuple[_ProjectedAxisAnnotation, ...],
        projected_axes: list[tuple[np.ndarray, np.ndarray]],
        box_center_xy: np.ndarray,
        opacities: np.ndarray,
    ) -> None:
        for annotation, (point1, point2), opacity in zip(
            annotations,
            projected_axes,
            opacities,
        ):
            annotation.update_from_projection(
                point1,
                point2,
                box_center_xy,
                float(opacity),
            )

    def update(self, *_event_args: Any) -> None:
        """Update positions and edge visibility after automatic or mouse rotation."""

        if self._updating:
            return
        self._updating = True
        try:
            xmin, xmax, ymin, ymax, zmin, zmax = self.bounds
            box_center_xy = self._project_to_display(
                (
                    (xmin + xmax) * 0.5,
                    (ymin + ymax) * 0.5,
                    (zmin + zmax) * 0.5,
                )
            )
            projected_x = self._project_annotations(self.x_annotations)
            projected_y = self._project_annotations(self.y_annotations)
            projected_z = self._project_annotations(self.z_annotations)
            projected_values = [
                box_center_xy,
                *(value for axis in (projected_x, projected_y, projected_z) for value in axis),
            ]
            if not all(np.all(np.isfinite(value)) for value in projected_values):
                return

            x_opacities = _continuous_edge_opacities(
                self._midpoint_scores(projected_x, coordinate=1),
                prefer_minimum=True,
            )
            y_opacities = _continuous_edge_opacities(
                self._midpoint_scores(projected_y, coordinate=1),
                prefer_minimum=True,
            )
            z_x_scores = self._midpoint_scores(projected_z, coordinate=0)
            z_opacities = np.maximum(
                _continuous_edge_opacities(z_x_scores, prefer_minimum=True),
                _continuous_edge_opacities(z_x_scores, prefer_minimum=False),
            )
            if self.y_vertical:
                y_x_scores = self._midpoint_scores(projected_y, coordinate=0)
                y_opacities = np.maximum(
                    _continuous_edge_opacities(y_x_scores, prefer_minimum=True),
                    _continuous_edge_opacities(y_x_scores, prefer_minimum=False),
                )
                z_opacities = _continuous_edge_opacities(
                    self._midpoint_scores(projected_z, coordinate=1),
                    prefer_minimum=True,
                )
            if self.single_edge_per_axis:
                # Pick one bottom edge for each horizontal axis and one left
                # edge for the vertical axis. No cross-fade: even at an edge
                # transition, two copies of the same ticks must never coexist.
                # Degenerate (end-on) edges cannot carry a readable scale.
                unique_opacities = []
                for axis_index, projected in enumerate((projected_x, projected_y, projected_z)):
                    vertical = axis_index == (1 if self.y_vertical else 2)
                    scores = self._midpoint_scores(projected, coordinate=0 if vertical else 1)
                    valid = np.asarray([np.linalg.norm(p2 - p1) > 1.0e-9 for p1, p2 in projected])
                    opacity = np.zeros(len(projected), dtype=float)
                    if np.any(valid):
                        opacity[int(np.argmin(np.where(valid, scores, np.inf)))] = 1.0
                    unique_opacities.append(opacity)
                x_opacities, y_opacities, z_opacities = unique_opacities
            self.x_opacities = tuple(float(value) for value in x_opacities)
            self.y_opacities = tuple(float(value) for value in y_opacities)
            self.z_opacities = tuple(float(value) for value in z_opacities)
            self._update_group(
                self.x_annotations,
                projected_x,
                box_center_xy,
                x_opacities,
            )
            self._update_group(
                self.y_annotations,
                projected_y,
                box_center_xy,
                y_opacities,
            )
            self._update_group(
                self.z_annotations,
                projected_z,
                box_center_xy,
                z_opacities,
            )
        finally:
            self._updating = False

    def install_camera_observer(self) -> None:
        """Update for automatic rotation and for direct mouse camera movement."""

        self.update()
        if self.camera_observer_id is None:
            self.camera_observer_id = int(
                self.camera.AddObserver("ModifiedEvent", self.update)
            )

    def dispose(self) -> None:
        """Detach the camera callback before a viewport is rebuilt or closed."""

        if self.camera_observer_id is not None:
            self.camera.RemoveObserver(self.camera_observer_id)
            self.camera_observer_id = None


def _install_synchronized_horizontal_rotation(
    plotter: Any,
    cameras: tuple[Any, ...],
    *,
    left_view_up: tuple[float, float, float] | None = None,
    right_view_up: tuple[float, float, float] | None = None,
) -> Any | None:
    """为所有给定视口安装同速、同方向的实时水平旋转。

    每个子图保留自己的焦点、缩放和裁剪范围，只共享每个定时步增加的方位角，
    所以左侧全局网络和右侧局部 ROI 不会因为强行共享同一台相机而错位。
    返回回调主要用于单元测试；PyVista 本身会通过 timer 保存并调用它。
    """

    if not AUTO_HORIZONTAL_ROTATION_ENABLED:
        return None
    if HORIZONTAL_ROTATION_INTERVAL_MS <= 0:
        raise ValueError("HORIZONTAL_ROTATION_INTERVAL_MS must be positive")
    if not np.isfinite(HORIZONTAL_ROTATION_DEGREES_PER_STEP):
        raise ValueError("HORIZONTAL_ROTATION_DEGREES_PER_STEP must be finite")
    if HORIZONTAL_ROTATION_MAX_STEPS <= 0:
        raise ValueError("HORIZONTAL_ROTATION_MAX_STEPS must be positive")
    view_up = np.asarray(HORIZONTAL_ROTATION_VIEW_UP if right_view_up is None else right_view_up, dtype=float)
    if view_up.shape != (3,) or not np.all(np.isfinite(view_up)):
        raise ValueError("Horizontal rotation view-up must contain three finite values")
    if float(np.linalg.norm(view_up)) <= 1.0e-12:
        raise ValueError("Horizontal rotation view-up must be non-zero")
    view_up /= np.linalg.norm(view_up)
    left_up = view_up if left_view_up is None else np.asarray(left_view_up, dtype=float)
    if left_up.shape != (3,) or not np.all(np.isfinite(left_up)) or np.linalg.norm(left_up) <= 1.0e-12:
        raise ValueError("left_view_up must contain three finite values and be non-zero")
    left_up = left_up / np.linalg.norm(left_up)

    # 多视口正常情况下各有一台相机；去重可防止调用者以后启用相机链接时，
    # 同一台相机在单个 timer step 中被重复旋转。
    unique_cameras: list[Any] = []
    seen_camera_ids: set[int] = set()
    for camera in cameras:
        identity = id(camera)
        if identity not in seen_camera_ids:
            unique_cameras.append(camera)
            seen_camera_ids.add(identity)
    if not unique_cameras:
        return None

    def rotate_viewports(_step: int) -> None:
        for index, camera in enumerate(unique_cameras):
            # 默认保留物理 Z 轴；调用者可分别指定左右视窗的竖直轴。
            camera.SetViewUp(*(left_up if index == 0 else view_up))
            camera.Azimuth(HORIZONTAL_ROTATION_DEGREES_PER_STEP)

    # PyVista 0.48.x 允许在 show() 前调用 add_timer_event，但此时若 VTK
    # interactor 尚未 Initialize，CreateRepeatingTimer 会静默返回 0；这个无效
    # timer 之后不会被 show() 自动重建。必须先初始化，再注册定时器。
    interactor = getattr(plotter, "iren", None)
    if interactor is not None and not bool(interactor.initialized):
        interactor.initialize()
    plotter.add_timer_event(
        max_steps=HORIZONTAL_ROTATION_MAX_STEPS,
        duration=HORIZONTAL_ROTATION_INTERVAL_MS,
        callback=rotate_viewports,
    )
    if interactor is not None:
        timer_id = getattr(getattr(interactor, "_timer", None), "id", None)
        if not timer_id:
            raise RuntimeError(
                "VTK horizontal-rotation timer could not be created after "
                "interactor initialization"
            )
    return rotate_viewports


def _style_text_property(
    text_property: Any,
    *,
    font_size: int,
    bold: bool,
) -> None:
    """Apply the shared readable Arial typography to a VTK text property."""

    text_property.SetFontFamilyToArial()
    text_property.SetFontSize(font_size)
    text_property.SetBold(bold)


def _style_orientation_axes(axes_actor: Any) -> None:
    for caption_getter in (
        axes_actor.GetXAxisCaptionActor2D,
        axes_actor.GetYAxisCaptionActor2D,
        axes_actor.GetZAxisCaptionActor2D,
    ):
        _style_text_property(
            caption_getter().GetCaptionTextProperty(),
            font_size=ORIENTATION_AXIS_FONT_SIZE,
            bold=True,
        )


def _style_bounds_axes(bounds_actor: Any) -> None:
    for axis_index in range(3):
        _style_text_property(
            bounds_actor.GetLabelTextProperty(axis_index),
            font_size=COORDINATE_TICK_FONT_SIZE,
            bold=False,
        )
        _style_text_property(
            bounds_actor.GetTitleTextProperty(axis_index),
            font_size=COORDINATE_TITLE_FONT_SIZE,
            bold=True,
        )


def _add_coordinate_text_actor(
    plotter: Any,
    *,
    text: str,
    font_size: int,
    bold: bool,
    name: str,
) -> Any:
    """Create a centred 2-D text actor in the current renderer's pixel space."""

    import pyvista as pv
    from vtkmodules.vtkRenderingCore import vtkTextActor

    actor = vtkTextActor()
    actor.SetInput(text)
    actor.GetPositionCoordinate().SetCoordinateSystemToDisplay()
    color = pv.Color(COORDINATE_AXIS_COLOR).float_rgb
    text_property = actor.GetTextProperty()
    text_property.SetColor(color)
    text_property.SetFontFamilyToArial()
    text_property.SetFontSize(font_size)
    text_property.SetBold(bold)
    text_property.SetJustificationToCentered()
    text_property.SetVerticalJustificationToCentered()
    plotter.add_actor(
        actor,
        name=name,
        reset_camera=False,
        pickable=False,
        render=False,
    )
    return actor


def _add_coordinate_tick_actor(
    plotter: Any,
    *,
    tick_count: int,
    name: str,
) -> tuple[Any, Any]:
    """Create one screen-space line actor containing all major tick marks."""

    import pyvista as pv
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
    from vtkmodules.vtkRenderingCore import vtkActor2D, vtkCoordinate, vtkPolyDataMapper2D

    points = vtkPoints()
    points.SetNumberOfPoints(tick_count * 2)
    lines = vtkCellArray()
    for tick_index in range(tick_count):
        lines.InsertNextCell(2)
        lines.InsertCellPoint(tick_index * 2)
        lines.InsertCellPoint(tick_index * 2 + 1)

    polydata = vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    display_coordinate = vtkCoordinate()
    display_coordinate.SetCoordinateSystemToDisplay()
    mapper = vtkPolyDataMapper2D()
    mapper.SetInputData(polydata)
    mapper.SetTransformCoordinate(display_coordinate)
    actor = vtkActor2D()
    actor.SetMapper(mapper)
    color = pv.Color(COORDINATE_AXIS_COLOR).float_rgb
    actor.GetProperty().SetColor(color)
    actor.GetProperty().SetLineWidth(COORDINATE_TICK_LINE_WIDTH_PX)
    actor.GetProperty().SetDisplayLocationToForeground()
    plotter.add_actor(
        actor,
        name=name,
        reset_camera=False,
        pickable=False,
        render=False,
    )
    return points, actor


def _add_projected_axis_annotation(
    plotter: Any,
    *,
    axis_name: str,
    edge_index: int,
    point1_xyz: tuple[float, float, float],
    point2_xyz: tuple[float, float, float],
    minimum_value: float,
    maximum_value: float,
    label_count: int,
) -> _ProjectedAxisAnnotation:
    """Create all screen-facing text actors attached to one fixed 3-D edge."""

    values = np.linspace(minimum_value, maximum_value, label_count)
    tick_line_points, tick_line_actor = _add_coordinate_tick_actor(
        plotter,
        tick_count=label_count,
        name=f"coordinate_{axis_name.lower()}_{edge_index}_tick_lines",
    )
    tick_actors = tuple(
        _add_coordinate_text_actor(
            plotter,
            text=f"{value:.1f}",
            font_size=COORDINATE_TICK_FONT_SIZE,
            bold=False,
            name=f"coordinate_{axis_name.lower()}_{edge_index}_tick_{tick_index}",
        )
        for tick_index, value in enumerate(values)
    )
    title_actor = _add_coordinate_text_actor(
        plotter,
        text=f"{axis_name} (um)",
        font_size=COORDINATE_TITLE_FONT_SIZE,
        bold=True,
        name=f"coordinate_{axis_name.lower()}_{edge_index}_title",
    )
    return _ProjectedAxisAnnotation(
        axis_name=axis_name,
        point1_xyz=point1_xyz,
        point2_xyz=point2_xyz,
        tick_line_points=tick_line_points,
        tick_line_actor=tick_line_actor,
        tick_actors=tick_actors,
        title_actor=title_actor,
    )


def _style_legend(legend_actor: Any) -> None:
    _style_text_property(
        legend_actor.GetEntryTextProperty(),
        font_size=LEGEND_FONT_SIZE,
        bold=False,
    )


def _add_viewport_rectangle_actor(
    plotter: Any,
    *,
    rectangle_xyxy: tuple[float, float, float, float],
    color: str,
    opacity: float,
    name: str,
) -> Any:
    """Add a filled rectangle using normalized coordinates of the active viewport."""

    import pyvista as pv
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
    from vtkmodules.vtkRenderingCore import vtkActor2D, vtkCoordinate, vtkPolyDataMapper2D

    xmin, ymin, xmax, ymax = rectangle_xyxy
    points = vtkPoints()
    for point in (
        (xmin, ymin, 0.0),
        (xmax, ymin, 0.0),
        (xmax, ymax, 0.0),
        (xmin, ymax, 0.0),
    ):
        points.InsertNextPoint(*point)
    polygon = vtkCellArray()
    polygon.InsertNextCell(4)
    for point_index in range(4):
        polygon.InsertCellPoint(point_index)
    polydata = vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(polygon)
    coordinate = vtkCoordinate()
    coordinate.SetCoordinateSystemToNormalizedViewport()
    mapper = vtkPolyDataMapper2D()
    mapper.SetInputData(polydata)
    mapper.SetTransformCoordinate(coordinate)
    actor = vtkActor2D()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(pv.Color(color).float_rgb)
    actor.GetProperty().SetOpacity(float(opacity))
    actor.GetProperty().SetDisplayLocationToForeground()
    plotter.add_actor(
        actor,
        name=name,
        reset_camera=False,
        pickable=False,
        render=False,
    )
    return actor


def _add_viewport_line_actor(
    plotter: Any,
    *,
    point1_xy: tuple[float, float],
    point2_xy: tuple[float, float],
    color: str,
    opacity: float,
    name: str,
) -> Any:
    """Add a true line segment in normalized coordinates of the active viewport."""

    import pyvista as pv
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
    from vtkmodules.vtkRenderingCore import vtkActor2D, vtkCoordinate, vtkPolyDataMapper2D

    points = vtkPoints()
    points.InsertNextPoint(float(point1_xy[0]), float(point1_xy[1]), 0.0)
    points.InsertNextPoint(float(point2_xy[0]), float(point2_xy[1]), 0.0)
    lines = vtkCellArray()
    lines.InsertNextCell(2)
    lines.InsertCellPoint(0)
    lines.InsertCellPoint(1)
    polydata = vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    coordinate = vtkCoordinate()
    coordinate.SetCoordinateSystemToNormalizedViewport()
    mapper = vtkPolyDataMapper2D()
    mapper.SetInputData(polydata)
    mapper.SetTransformCoordinate(coordinate)
    actor = vtkActor2D()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(pv.Color(color).float_rgb)
    actor.GetProperty().SetOpacity(float(opacity))
    actor.GetProperty().SetLineWidth(LEGEND_LINE_WIDTH_PX)
    actor.GetProperty().SetDisplayLocationToForeground()
    plotter.add_actor(
        actor,
        name=name,
        reset_camera=False,
        pickable=False,
        render=False,
    )
    return actor


def _add_viewport_arrow_actor(
    plotter: Any,
    *,
    point1_xy: tuple[float, float],
    point2_xy: tuple[float, float],
    color: str,
    opacity: float,
    name: str,
) -> Any:
    """Add an open line-arrow symbol without a filled triangular glyph."""

    import pyvista as pv
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
    from vtkmodules.vtkRenderingCore import vtkActor2D, vtkCoordinate, vtkPolyDataMapper2D

    start = np.asarray(point1_xy, dtype=float)
    tip = np.asarray(point2_xy, dtype=float)
    arrow_vector = tip - start
    arrow_length = float(np.linalg.norm(arrow_vector))
    if arrow_length <= 1.0e-12:
        raise ValueError("Legend arrow endpoints must be different")
    direction = arrow_vector / arrow_length
    perpendicular = np.asarray((-direction[1], direction[0]), dtype=float)
    head_base = tip - direction * min(
        LEGEND_ARROW_HEAD_LENGTH_FRACTION,
        arrow_length * 0.45,
    )
    head_upper = head_base + perpendicular * LEGEND_ARROW_HEAD_HALF_HEIGHT_FRACTION
    head_lower = head_base - perpendicular * LEGEND_ARROW_HEAD_HALF_HEIGHT_FRACTION

    points = vtkPoints()
    for point in (start, tip, head_upper, head_lower):
        points.InsertNextPoint(float(point[0]), float(point[1]), 0.0)
    lines = vtkCellArray()
    for point_indices in ((0, 1), (1, 2), (1, 3)):
        lines.InsertNextCell(2)
        lines.InsertCellPoint(point_indices[0])
        lines.InsertCellPoint(point_indices[1])
    polydata = vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(lines)
    coordinate = vtkCoordinate()
    coordinate.SetCoordinateSystemToNormalizedViewport()
    mapper = vtkPolyDataMapper2D()
    mapper.SetInputData(polydata)
    mapper.SetTransformCoordinate(coordinate)
    actor = vtkActor2D()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(pv.Color(color).float_rgb)
    actor.GetProperty().SetOpacity(float(opacity))
    actor.GetProperty().SetLineWidth(LEGEND_LINE_WIDTH_PX)
    actor.GetProperty().SetDisplayLocationToForeground()
    plotter.add_actor(
        actor,
        name=name,
        reset_camera=False,
        pickable=False,
        render=False,
    )
    return actor


def _add_viewport_circle_actor(
    plotter: Any,
    *,
    center_xy: tuple[float, float],
    radius_fraction: float,
    color: str,
    opacity: float,
    name: str,
) -> Any:
    """Add a filled circular marker in normalized coordinates of the viewport."""

    import pyvista as pv
    from vtkmodules.vtkCommonCore import vtkPoints
    from vtkmodules.vtkCommonDataModel import vtkCellArray, vtkPolyData
    from vtkmodules.vtkRenderingCore import vtkActor2D, vtkCoordinate, vtkPolyDataMapper2D

    if LEGEND_CIRCLE_RESOLUTION < 8:
        raise ValueError("LEGEND_CIRCLE_RESOLUTION must be at least 8")
    # 归一化视口的 X/Y 比例在非正方形子图中不同。按当前 renderer 的像素
    # 长宽比修正 X 半径，确保屏幕上看到的始终是圆，而不是椭圆。
    viewport_width, viewport_height = plotter.renderer.GetSize()
    x_radius = float(radius_fraction)
    if viewport_width > 0 and viewport_height > 0:
        x_radius *= float(viewport_height) / float(viewport_width)
    angles = np.linspace(0.0, 2.0 * np.pi, LEGEND_CIRCLE_RESOLUTION, endpoint=False)
    points = vtkPoints()
    for angle in angles:
        points.InsertNextPoint(
            float(center_xy[0] + x_radius * np.cos(angle)),
            float(center_xy[1] + radius_fraction * np.sin(angle)),
            0.0,
        )
    polygon = vtkCellArray()
    polygon.InsertNextCell(LEGEND_CIRCLE_RESOLUTION)
    for point_index in range(LEGEND_CIRCLE_RESOLUTION):
        polygon.InsertCellPoint(point_index)
    polydata = vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetPolys(polygon)
    coordinate = vtkCoordinate()
    coordinate.SetCoordinateSystemToNormalizedViewport()
    mapper = vtkPolyDataMapper2D()
    mapper.SetInputData(polydata)
    mapper.SetTransformCoordinate(coordinate)
    actor = vtkActor2D()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(pv.Color(color).float_rgb)
    actor.GetProperty().SetOpacity(float(opacity))
    actor.GetProperty().SetDisplayLocationToForeground()
    plotter.add_actor(
        actor,
        name=name,
        reset_camera=False,
        pickable=False,
        render=False,
    )
    return actor


def _add_scientific_bottom_legend(
    plotter: Any,
    entries: list[_ScientificLegendEntry],
    *,
    name_prefix: str,
) -> tuple[Any, ...]:
    """Draw a compact row-major legend with line, circle, and patch symbols."""

    if not entries:
        return ()
    if LEGEND_COLUMN_COUNT <= 0:
        raise ValueError("LEGEND_COLUMN_COUNT must be positive")

    import pyvista as pv
    from vtkmodules.vtkRenderingCore import vtkTextActor

    column_count = min(LEGEND_COLUMN_COUNT, len(entries))
    row_count = int(np.ceil(len(entries) / column_count))
    left = LEGEND_LEFT_FRACTION
    right = 1.0 - LEGEND_RIGHT_MARGIN_FRACTION
    bottom = LEGEND_BOTTOM_FRACTION
    height = LEGEND_HEIGHT_FRACTION
    width = right - left
    cell_width = width / column_count
    required_content_height = (row_count - 1) * LEGEND_ROW_SPACING_FRACTION
    if required_content_height >= height:
        raise ValueError(
            "Legend rows exceed LEGEND_HEIGHT_FRACTION; reduce "
            "LEGEND_ROW_SPACING_FRACTION or increase the legend height"
        )
    first_row_y = bottom + height * 0.5 + required_content_height * 0.5
    actors: list[Any] = []
    actors.append(
        _add_viewport_rectangle_actor(
            plotter,
            rectangle_xyxy=(left, bottom, right, bottom + height),
            color=LEGEND_BACKGROUND_COLOR,
            opacity=LEGEND_BACKGROUND_OPACITY,
            name=f"{name_prefix}_background",
        )
    )
    for entry_index, entry in enumerate(entries):
        column_index = entry_index % column_count
        row_index = entry_index // column_count
        cell_left = left + column_index * cell_width
        centre_y = first_row_y - row_index * LEGEND_ROW_SPACING_FRACTION
        swatch_left = cell_left + 0.010
        swatch_right = swatch_left + LEGEND_SWATCH_WIDTH_FRACTION
        symbol_name = f"{name_prefix}_symbol_{entry_index}"
        if entry.symbol == "line":
            symbol_actor = _add_viewport_line_actor(
                plotter,
                point1_xy=(swatch_left, centre_y),
                point2_xy=(swatch_right, centre_y),
                color=entry.color,
                opacity=1.0,
                name=symbol_name,
            )
        elif entry.symbol == "arrow":
            symbol_actor = _add_viewport_arrow_actor(
                plotter,
                point1_xy=(swatch_left, centre_y),
                point2_xy=(swatch_right, centre_y),
                color=entry.color,
                opacity=1.0,
                name=symbol_name,
            )
        elif entry.symbol == "circle":
            symbol_actor = _add_viewport_circle_actor(
                plotter,
                center_xy=((swatch_left + swatch_right) * 0.5, centre_y),
                radius_fraction=LEGEND_CIRCLE_RADIUS_FRACTION,
                color=entry.color,
                opacity=1.0,
                name=symbol_name,
            )
        elif entry.symbol == "patch":
            symbol_actor = _add_viewport_rectangle_actor(
                plotter,
                rectangle_xyxy=(
                    swatch_left,
                    centre_y - LEGEND_SWATCH_HEIGHT_FRACTION * 0.5,
                    swatch_right,
                    centre_y + LEGEND_SWATCH_HEIGHT_FRACTION * 0.5,
                ),
                color=entry.color,
                opacity=1.0,
                name=symbol_name,
            )
        else:  # pragma: no cover - Literal keeps normal callers out of this branch.
            raise ValueError(f"Unsupported legend symbol: {entry.symbol}")
        actors.append(symbol_actor)
        text_actor = vtkTextActor()
        text_actor.SetInput(entry.label)
        text_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        text_actor.SetPosition(
            swatch_right + LEGEND_TEXT_GAP_FRACTION,
            centre_y,
        )
        text_property = text_actor.GetTextProperty()
        text_property.SetColor(pv.Color(LEGEND_TEXT_COLOR).float_rgb)
        text_property.SetFontFamilyToArial()
        text_property.SetFontSize(LEGEND_FONT_SIZE)
        text_property.SetBold(False)
        text_property.SetJustificationToLeft()
        text_property.SetVerticalJustificationToCentered()
        plotter.add_actor(
            text_actor,
            name=f"{name_prefix}_text_{entry_index}",
            reset_camera=False,
            pickable=False,
            render=False,
        )
        actors.append(text_actor)
    return tuple(actors)


def scene_geometry_from_graph(
    result: DirectedVascularGraph,
    *,
    max_arrows: int,
) -> InteractiveSceneGeometry:
    arrow_points, arrow_vectors, _ = sample_direction_arrows(result, max_arrows)
    ordered_nodes = sorted(result.junction_graph.nodes(data=True), key=lambda item: int(item[0]))
    critical_points = np.asarray(
        [
            (float(data["x_um"]), float(data["y_um"]), float(data["z_um"]))
            for _, data in ordered_nodes
        ],
        dtype=float,
    )
    if not len(critical_points):
        critical_points = np.empty((0, 3), dtype=float)
    branches = tuple(branch.derived_points_um.copy() for branch in result.branches)
    return InteractiveSceneGeometry(
        branches_um=branches,
        arrow_points_um=arrow_points,
        arrow_vectors_xyz=arrow_vectors,
        critical_points_um=critical_points,
        critical_roles=tuple(str(data["role"]) for _, data in ordered_nodes),
    )


def _model_bounds(
    volume_zyx: np.ndarray | None,
    spacing_xyz_um: tuple[float, float, float],
    geometry: InteractiveSceneGeometry | None = None,
) -> tuple[float, float, float, float, float, float]:
    if volume_zyx is not None:
        shape_xyz = np.asarray(volume_zyx.shape[::-1], dtype=float)
        maximum = np.maximum(shape_xyz - 1.0, 0.0) * np.asarray(spacing_xyz_um, dtype=float)
        return (0.0, float(maximum[0]), 0.0, float(maximum[1]), 0.0, float(maximum[2]))
    point_arrays: list[np.ndarray] = []
    if geometry is not None:
        point_arrays.extend(branch for branch in geometry.branches_um if len(branch))
        if len(geometry.critical_points_um):
            point_arrays.append(geometry.critical_points_um)
    if not point_arrays:
        return (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    points = np.concatenate(point_arrays, axis=0)
    lower = np.min(points, axis=0)
    upper = np.max(points, axis=0)
    upper = np.where(upper > lower, upper, lower + 1.0)
    return (
        float(lower[0]), float(upper[0]),
        float(lower[1]), float(upper[1]),
        float(lower[2]), float(upper[2]),
    )


def _normalized_volume(volume_zyx: np.ndarray) -> tuple[np.ndarray, dict[str, float]]:
    values = np.asarray(volume_zyx, dtype=np.float32)
    finite = values[np.isfinite(values)]
    positive = finite[finite > 0]
    if positive.size:
        lower = float(np.percentile(positive, 12.0))
        upper = float(np.percentile(positive, 99.8))
    elif finite.size:
        lower, upper = float(np.min(finite)), float(np.max(finite))
    else:
        raise ValueError("The TIFF volume contains no finite intensity values")
    if upper <= lower:
        lower = 0.0
        upper = max(float(np.max(finite)), 1.0)
    normalized = np.nan_to_num((values - lower) / (upper - lower), nan=0.0)
    normalized = np.clip(normalized, 0.0, 1.0)
    return normalized, {"source_lower": lower, "source_upper": upper}


def _volume_grid(volume_zyx: np.ndarray, spacing_xyz_um: tuple[float, float, float]) -> Any:
    import pyvista as pv

    z_size, y_size, x_size = volume_zyx.shape
    grid = pv.ImageData(
        dimensions=(x_size, y_size, z_size),
        spacing=spacing_xyz_um,
        origin=(0.0, 0.0, 0.0),
    )
    xyz = np.transpose(volume_zyx, (2, 1, 0))
    grid.point_data["intensity"] = xyz.ravel(order="F")
    return grid


def _branch_mesh(
    branches_um: tuple[np.ndarray, ...], radius_um: tuple[np.ndarray, ...] | None = None
) -> Any | None:
    import pyvista as pv

    usable = [np.asarray(points, dtype=float) for points in branches_um if len(points) >= 2]
    if not usable:
        return None
    points = np.concatenate(usable, axis=0)
    cells: list[np.ndarray] = []
    offset = 0
    for branch in usable:
        indices = np.arange(offset, offset + len(branch), dtype=np.int64)
        cells.append(np.concatenate(([len(branch)], indices)))
        offset += len(branch)
    mesh = pv.PolyData(points)
    mesh.lines = np.concatenate(cells)
    if radius_um is not None:
        mesh.point_data["radius_um"] = np.concatenate(radius_um)
    return mesh


def _add_coordinate_reference_geometry(
    plotter: Any,
    bounds: tuple[float, float, float, float, float, float],
    *,
    draw_frame: bool,
    draw_grid: bool,
) -> None:
    """Draw a frame and face grid from the exact geometry used by tick projection."""

    import pyvista as pv

    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    if draw_frame:
        reference_box = pv.Box(bounds=bounds)
        plotter.add_mesh(
            reference_box,
            style="wireframe",
            color=COORDINATE_AXIS_COLOR,
            line_width=COORDINATE_FRAME_LINE_WIDTH_PX,
            reset_camera=False,
            pickable=False,
            render=False,
        )
    if not draw_grid:
        return

    # 六个框面都使用较淡的显式网格。框线、网格和刻度由完全相同的物理坐标
    # 构造，因此不会出现 CubeAxesActor 框线与屏幕刻度投影不重合的问题。
    segments: list[tuple[tuple[float, float, float], tuple[float, float, float]]] = []
    x_values = np.linspace(xmin, xmax, COORDINATE_X_TICK_COUNT)[1:-1]
    y_values = np.linspace(ymin, ymax, COORDINATE_Y_TICK_COUNT)[1:-1]
    z_values = np.linspace(zmin, zmax, COORDINATE_Z_TICK_COUNT)[1:-1]
    for x_value in x_values:
        x = float(x_value)
        segments.extend(
            (
                ((x, ymin, zmin), (x, ymax, zmin)),
                ((x, ymin, zmax), (x, ymax, zmax)),
                ((x, ymin, zmin), (x, ymin, zmax)),
                ((x, ymax, zmin), (x, ymax, zmax)),
            )
        )
    for y_value in y_values:
        y = float(y_value)
        segments.extend(
            (
                ((xmin, y, zmin), (xmax, y, zmin)),
                ((xmin, y, zmax), (xmax, y, zmax)),
                ((xmin, y, zmin), (xmin, y, zmax)),
                ((xmax, y, zmin), (xmax, y, zmax)),
            )
        )
    for z_value in z_values:
        z = float(z_value)
        segments.extend(
            (
                ((xmin, ymin, z), (xmax, ymin, z)),
                ((xmin, ymax, z), (xmax, ymax, z)),
                ((xmin, ymin, z), (xmin, ymax, z)),
                ((xmax, ymin, z), (xmax, ymax, z)),
            )
        )
    if not segments:
        return
    grid_points = np.asarray(
        [point for segment in segments for point in segment],
        dtype=float,
    )
    grid_cells = np.asarray(
        [
            (2, segment_index * 2, segment_index * 2 + 1)
            for segment_index in range(len(segments))
        ],
        dtype=np.int64,
    )
    grid_mesh = pv.PolyData(grid_points)
    grid_mesh.lines = grid_cells.ravel()
    plotter.add_mesh(
        grid_mesh,
        color=COORDINATE_AXIS_COLOR,
        line_width=COORDINATE_GRID_LINE_WIDTH_PX,
        opacity=COORDINATE_GRID_OPACITY,
        lighting=False,
        reset_camera=False,
        pickable=False,
        render=False,
    )


def _add_physical_coordinate_axes(
    plotter: Any,
    bounds: tuple[float, float, float, float, float, float],
    *,
    add_orientation_axes: bool = True,
    show_reference_frame: bool = True,
    show_reference_grid: bool | None = None,
    y_vertical: bool = False,
) -> _ScreenAnchoredCoordinateAxes:
    """Show camera-aware physical axes, optionally with a white reference grid."""

    if add_orientation_axes:
        axes_actor = plotter.add_axes(
            xlabel="X (um)",
            ylabel="Y (um)",
            zlabel="Z (um)",
            color="white",
        )
        _style_orientation_axes(axes_actor)
    if show_reference_grid is None:
        show_reference_grid = show_reference_frame
    if show_reference_frame or show_reference_grid:
        _add_coordinate_reference_geometry(
            plotter,
            bounds,
            draw_frame=show_reference_frame,
            draw_grid=show_reference_grid,
        )
    xmin, xmax, ymin, ymax, zmin, zmax = bounds
    x_edges = (
        ((xmin, ymin, zmin), (xmax, ymin, zmin)),
        ((xmin, ymax, zmin), (xmax, ymax, zmin)),
    )
    y_edges = (
        ((xmin, ymin, zmin), (xmin, ymax, zmin)),
        ((xmax, ymin, zmin), (xmax, ymax, zmin)),
    )
    z_edges = tuple(
        ((x_value, y_value, zmin), (x_value, y_value, zmax))
        for x_value in (xmin, xmax)
        for y_value in (ymin, ymax)
    )
    if y_vertical:
        # Same bottom/side annotation layout, with Y taking the vertical role.
        x_edges = tuple(
            ((xmin, ymin, z_value), (xmax, ymin, z_value)) for z_value in (zmin, zmax)
        )
        z_edges = tuple(
            ((x_value, ymin, zmin), (x_value, ymin, zmax)) for x_value in (xmin, xmax)
        )
        y_edges = tuple(
            ((x_value, ymin, z_value), (x_value, ymax, z_value))
            for x_value in (xmin, xmax) for z_value in (zmin, zmax)
        )
    x_annotations = tuple(
        _add_projected_axis_annotation(
            plotter,
            axis_name="X",
            edge_index=index,
            point1_xyz=point1,
            point2_xyz=point2,
            minimum_value=xmin,
            maximum_value=xmax,
            label_count=COORDINATE_X_TICK_COUNT,
        )
        for index, (point1, point2) in enumerate(x_edges)
    )
    y_annotations = tuple(
        _add_projected_axis_annotation(
            plotter,
            axis_name="Y",
            edge_index=index,
            point1_xyz=point1,
            point2_xyz=point2,
            minimum_value=ymin,
            maximum_value=ymax,
            label_count=COORDINATE_Y_TICK_COUNT,
        )
        for index, (point1, point2) in enumerate(y_edges)
    )
    z_annotations = tuple(
        _add_projected_axis_annotation(
            plotter,
            axis_name="Z",
            edge_index=index,
            point1_xyz=point1,
            point2_xyz=point2,
            minimum_value=zmin,
            maximum_value=zmax,
            label_count=COORDINATE_Z_TICK_COUNT,
        )
        for index, (point1, point2) in enumerate(z_edges)
    )
    controller = _ScreenAnchoredCoordinateAxes(
        renderer=plotter.renderer,
        camera=plotter.camera,
        bounds=bounds,
        x_annotations=x_annotations,
        y_annotations=y_annotations,
        z_annotations=z_annotations,
        y_vertical=y_vertical,
    )
    controller.install_camera_observer()
    return controller


def _set_full_scene_title(plotter: Any, *, sampling_available: bool) -> None:
    detail = (
        "Representative connected ROIs\n"
        "A: all ROI candidates | R/S: selected | C: next cluster | left-click: inspect"
        if sampling_available
        else "Orange arrows = SWC parent -> current"
    )
    plotter.add_text(
        f"Preprocessed main vascular network\n{detail}",
        position="upper_left",
        font_size=11,
        color="white",
        font=UI_FONT_FAMILY,
        name="full_scene_title",
    )


def _add_full_scene(
    plotter: Any,
    volume_zyx: np.ndarray | None,
    geometry: InteractiveSceneGeometry,
    *,
    spacing_xyz_um: tuple[float, float, float],
    volume_opacity: float,
    sample_id: str,
    sampling_available: bool = False,
    left_view_up: tuple[float, float, float] | None = None,
) -> tuple[dict[str, Any], _ScreenAnchoredCoordinateAxes]:
    import pyvista as pv

    plotter.set_background("black")
    legend_entries: list[_ScientificLegendEntry] = []
    intensity_window: dict[str, float] | None = None
    if volume_zyx is not None:
        normalized, intensity_window = _normalized_volume(volume_zyx)
        grid = _volume_grid(normalized, spacing_xyz_um)
        opacity = np.asarray([0.0, 0.0, 0.005, 0.02, 0.08, 0.22, 0.58, 1.0])
        opacity *= float(volume_opacity)
        plotter.add_volume(
            grid,
            scalars="intensity",
            cmap="gray",
            opacity=opacity,
            clim=(0.0, 1.0),
            shade=True,
            ambient=0.20,
            diffuse=0.85,
            specular=0.25,
            show_scalar_bar=False,
            pickable=False,
        )

    branch_mesh = _branch_mesh(geometry.branches_um)
    if branch_mesh is not None:
        plotter.add_mesh(
            branch_mesh,
            color="#35D6E3",
            line_width=2.4,
            opacity=0.86,
            label="SWC centerline",
            pickable=False,
        )
        legend_entries.append(
            _ScientificLegendEntry("SWC centerline", "#35D6E3", "line")
        )

    if geometry.arrow_count:
        arrows = pv.PolyData(geometry.arrow_points_um)
        arrows.point_data["parent_to_current"] = geometry.arrow_vectors_xyz
        bounds = _model_bounds(volume_zyx, spacing_xyz_um, geometry)
        diagonal = float(
            np.linalg.norm(
                np.asarray((bounds[1] - bounds[0], bounds[3] - bounds[2], bounds[5] - bounds[4]))
            )
        )
        glyphs = arrows.glyph(
            orient="parent_to_current",
            scale=False,
            factor=max(
                LEFT_PARENT_ARROW_MIN_LENGTH_UM,
                diagonal * LEFT_PARENT_ARROW_SCALE_FRACTION,
            ),
        )
        plotter.add_mesh(
            glyphs,
            color="#FF9F1C",
            label="parent -> current",
            lighting=True,
            pickable=False,
        )
        legend_entries.append(
            _ScientificLegendEntry("parent -> current", "#FF9F1C", "arrow")
        )

    role_styles = {
        "inferred_inlet": (
            "#4ADE80",
            STRUCTURAL_ROOT_POINT_SIZE_PX,
            "structural root (not flow inlet)",
            "structural root",
        ),
        "inferred_outlet": (
            "#F87171",
            STRUCTURAL_LEAF_POINT_SIZE_PX,
            "structural leaf (not flow outlet)",
            "structural leaf",
        ),
        "divergence_junction": (
            "#FACC15",
            DIVERGENCE_JUNCTION_POINT_SIZE_PX,
            "divergence junction",
            "divergence junction",
        ),
        "convergence_junction": (
            "#C084FC",
            CONVERGENCE_JUNCTION_POINT_SIZE_PX,
            "convergence junction",
            "convergence junction",
        ),
    }
    roles = np.asarray(geometry.critical_roles, dtype=object)
    for role, (color, size, display_label, legend_label) in role_styles.items():
        selected = geometry.critical_points_um[roles == role]
        if len(selected):
            plotter.add_points(
                selected,
                color=color,
                point_size=size,
                render_points_as_spheres=True,
                label=display_label,
                pickable=False,
            )
            legend_entries.append(
                _ScientificLegendEntry(legend_label, color, "circle")
            )

    _set_full_scene_title(plotter, sampling_available=sampling_available)
    model_bounds = _model_bounds(volume_zyx, spacing_xyz_um, geometry)
    if sampling_available:
        legend_entries.extend(
            (
                _ScientificLegendEntry(
                    "candidate ROI",
                    CANDIDATE_ROI_LEGEND_COLOR,
                    "patch",
                ),
                _ScientificLegendEntry(
                    "active ROI",
                    ACTIVE_ROI_FILL_COLOR,
                    "patch",
                ),
            )
        )
    _add_scientific_bottom_legend(
        plotter,
        legend_entries,
        name_prefix="full_scene_legend",
    )
    plotter.view_isometric()
    if left_view_up is not None:
        plotter.camera.SetViewUp(*left_view_up)
    # Leave enough screen margin for the outward-facing coordinate annotations.
    plotter.camera.zoom(1.02)
    coordinate_axes = _add_physical_coordinate_axes(
        plotter, model_bounds, y_vertical=left_view_up == (0.0, 1.0, 0.0),
    )
    metadata = {
        "sample_id": sample_id,
        "volume_shape_zyx": list(map(int, volume_zyx.shape)) if volume_zyx is not None else None,
        "optional_background_volume_available": volume_zyx is not None,
        "spacing_xyz_um": list(map(float, spacing_xyz_um)),
        "branch_count": geometry.branch_count,
        "arrow_count": geometry.arrow_count,
        "critical_node_count": len(geometry.critical_points_um),
        "global_coordinate_bounds_xyz_um": list(model_bounds),
        "coordinate_units": "um",
        "interface_style": {
            "font_family": "Arial",
            "coordinate_tick_font_size": COORDINATE_TICK_FONT_SIZE,
            "coordinate_title_font_size": COORDINATE_TITLE_FONT_SIZE,
            "coordinate_frame_line_width_px": COORDINATE_FRAME_LINE_WIDTH_PX,
            "coordinate_grid_line_width_px": COORDINATE_GRID_LINE_WIDTH_PX,
            "coordinate_grid_opacity": COORDINATE_GRID_OPACITY,
            "coordinate_tick_length_px": COORDINATE_TICK_LENGTH_PX,
            "coordinate_tick_line_width_px": COORDINATE_TICK_LINE_WIDTH_PX,
            "coordinate_tick_counts_xyz": [
                COORDINATE_X_TICK_COUNT,
                COORDINATE_Y_TICK_COUNT,
                COORDINATE_Z_TICK_COUNT,
            ],
            "legend_font_size": LEGEND_FONT_SIZE,
            "font_scope": "Arial for all numeric and English text",
            "critical_point_sizes_px": {
                "structural_root": STRUCTURAL_ROOT_POINT_SIZE_PX,
                "structural_leaf": STRUCTURAL_LEAF_POINT_SIZE_PX,
                "divergence_junction": DIVERGENCE_JUNCTION_POINT_SIZE_PX,
                "convergence_junction": CONVERGENCE_JUNCTION_POINT_SIZE_PX,
            },
            "left_parent_arrow_scale_fraction": LEFT_PARENT_ARROW_SCALE_FRACTION,
            "overlay_emphasis": "enhanced in both viewports",
        },
        "intensity_window": intensity_window,
        "direction_rule": "SWC parent_id node -> current node",
        "direction_is_measured_flow": False,
        "interaction": {
            "selection": "left-click a sampled ROI cube in the full-model viewport"
            if sampling_available
            else "rotate, pan, and zoom the preprocessed full model",
            "result": "the right viewport is replaced by the selected connected ROI"
            if sampling_available
            else "the full-model view remains interactive",
        },
    }
    return metadata, coordinate_axes


def _sampling_color(cluster_id: int) -> str:
    palette = (
        "#4DD0E1", "#FF8A65", "#FFD54F", "#AB47BC", "#66BB6A",
        "#42A5F5", "#EC407A", "#9CCC65", "#FFA726", "#7E57C2",
    )
    return palette[int(cluster_id) % len(palette)] if cluster_id >= 0 else "#B0BEC5"


def _sampling_diameter_limits_um(
    rois: tuple[ROIRecord, ...],
) -> tuple[float, float] | None:
    """Return one finite positive diameter range shared by all displayed ROIs."""

    diameter_chunks: list[np.ndarray] = []
    for roi in rois:
        for edge_radius_um in getattr(roi, "local_edge_radius_um", ()):
            diameter_um = 2.0 * np.asarray(edge_radius_um, dtype=float).ravel()
            usable = diameter_um[np.isfinite(diameter_um) & (diameter_um > 0.0)]
            if usable.size:
                diameter_chunks.append(usable)
    if not diameter_chunks:
        return None
    diameter_values = np.concatenate(diameter_chunks)
    lower = float(np.min(diameter_values))
    upper = float(np.max(diameter_values))
    if upper <= lower:
        padding = max(abs(lower) * 0.05, 0.1)
        lower = max(0.0, lower - padding)
        upper += padding
    return lower, upper


def _add_diameter_colorbar(plotter: Any, mapper: Any) -> Any:
    """Add the right-view scientific color bar for physical vessel diameter."""

    scalar_bar = plotter.add_scalar_bar(
        title=DIAMETER_COLORBAR_TITLE,
        mapper=mapper,
        n_labels=DIAMETER_COLORBAR_LABEL_COUNT,
        title_font_size=DIAMETER_COLORBAR_TITLE_FONT_SIZE,
        label_font_size=DIAMETER_COLORBAR_LABEL_FONT_SIZE,
        color=LEGEND_TEXT_COLOR,
        font_family=UI_FONT_FAMILY,
        width=DIAMETER_COLORBAR_WIDTH_FRACTION,
        height=DIAMETER_COLORBAR_HEIGHT_FRACTION,
        position_x=DIAMETER_COLORBAR_POSITION_X,
        position_y=DIAMETER_COLORBAR_POSITION_Y,
        vertical=True,
        interactive=False,
        fmt=DIAMETER_COLORBAR_FORMAT,
        use_opacity=False,
        outline=False,
        n_colors=256,
        unconstrained_font_size=True,
        render=False,
    )
    # PyVista 的字体参数负责初始样式；这里再直接约束 VTK 文本属性，确保标题
    # 和所有数字刻度始终为 Arial。正 separation 将标题向色条上方移开。
    scalar_bar.GetTitleTextProperty().SetFontFamilyToArial()
    scalar_bar.GetLabelTextProperty().SetFontFamilyToArial()
    scalar_bar.SetVerticalTitleSeparation(DIAMETER_COLORBAR_TITLE_SEPARATION_PX)
    return scalar_bar


def _remove_diameter_colorbar(plotter: Any) -> None:
    """Remove the previous ROI diameter bar before rebuilding the right view."""

    scalar_bars = getattr(plotter, "scalar_bars", None)
    if scalar_bars is not None and DIAMETER_COLORBAR_TITLE in scalar_bars:
        plotter.remove_scalar_bar(DIAMETER_COLORBAR_TITLE, render=False)


def _spatially_spread_indices(
    points_xyz: np.ndarray,
    *,
    maximum_count: int,
) -> np.ndarray:
    """Select deterministic, spatially distributed points for sparse arrows."""

    points = np.asarray(points_xyz, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("Arrow candidate points must have shape (N, 3)")
    if maximum_count <= 0:
        raise ValueError("maximum_count must be positive")
    if len(points) <= maximum_count:
        return np.arange(len(points), dtype=np.int64)

    # 从最接近整体中心的点开始，再反复选择距当前集合最远的点。与简单按数组
    # 下标抽取相比，这能覆盖不同空间分支并避免箭头重新挤在局部血管段上。
    centre = np.mean(points, axis=0)
    first_index = int(np.argmin(np.sum((points - centre) ** 2, axis=1)))
    selected = [first_index]
    selected_mask = np.zeros(len(points), dtype=bool)
    selected_mask[first_index] = True
    minimum_squared_distance = np.sum(
        (points - points[first_index]) ** 2,
        axis=1,
    )
    while len(selected) < maximum_count:
        candidate_scores = np.where(
            selected_mask,
            -np.inf,
            minimum_squared_distance,
        )
        next_index = int(np.argmax(candidate_scores))
        selected.append(next_index)
        selected_mask[next_index] = True
        squared_distance = np.sum((points - points[next_index]) ** 2, axis=1)
        minimum_squared_distance = np.minimum(
            minimum_squared_distance,
            squared_distance,
        )
    return np.asarray(selected, dtype=np.int64)


def _add_sampling_roi_scene(
    plotter: Any,
    roi: ROIRecord,
    *,
    add_orientation_axes: bool = True,
    diameter_clim_um: tuple[float, float] | None = None,
    view_up: tuple[float, float, float] | None = None,
) -> _ScreenAnchoredCoordinateAxes:
    """Render one saved connected ROI without running sampling or clustering."""

    import pyvista as pv

    plotter.set_background("black")
    legend_entries = [
        _ScientificLegendEntry("ROI boundary", ROI_BOUNDARY_COLOR, "line")
    ]
    bounds = (
        roi.bbox_min_um[0], roi.bbox_max_um[0],
        roi.bbox_min_um[1], roi.bbox_max_um[1],
        roi.bbox_min_um[2], roi.bbox_max_um[2],
    )
    box = pv.Box(bounds=bounds)
    plotter.add_mesh(box, color=ROI_BOUNDARY_COLOR, opacity=0.025, pickable=False)
    plotter.add_mesh(
        box,
        style="wireframe",
        color=ROI_BOUNDARY_COLOR,
        line_width=4.0,
        label="spatial ROI boundary",
        pickable=False,
    )
    display_geometry = build_roi_display_tubes(roi)
    if display_geometry is not None:
        mesh, tubes = display_geometry
        if diameter_clim_um is None:
            diameter_clim_um = _sampling_diameter_limits_um((roi,))
        tube_actor = plotter.add_mesh(
            tubes,
            scalars=DIAMETER_SCALAR_NAME,
            cmap=DIAMETER_COLORMAP,
            clim=diameter_clim_um,
            smooth_shading=True,
            opacity=0.88,
            show_scalar_bar=False,
            pickable=False,
        )
        _add_diameter_colorbar(plotter, tube_actor.mapper)
        plotter.add_mesh(
            mesh,
            color="#51E5FF",
            line_width=3.0,
            label="connected ROI centerline",
            name="sampling_roi_centerline",
            pickable=False,
        )
        legend_entries.extend(
            (
                _ScientificLegendEntry("ROI centerline", "#51E5FF", "line"),
            )
        )
    if len(roi.local_edge_points_um):
        starts = roi.local_edge_points_um[:, 0]
        ends = roi.local_edge_points_um[:, 1]
        vectors = ends - starts
        lengths = np.linalg.norm(vectors, axis=1)
        valid = lengths > 1.0e-12
        candidate_points = (starts[valid] + ends[valid]) * 0.5
        candidate_directions = vectors[valid] / lengths[valid, None]
        arrow_indices = _spatially_spread_indices(
            candidate_points,
            maximum_count=RIGHT_PARENT_ARROW_MAX_COUNT,
        )
        arrow_points = candidate_points[arrow_indices]
        directions = candidate_directions[arrow_indices]
        if len(arrow_points):
            arrows = pv.PolyData(arrow_points)
            arrows.point_data["parent_to_current"] = directions
            diagonal = float(np.linalg.norm(np.asarray(roi.bbox_size_um)))
            glyphs = arrows.glyph(
                orient="parent_to_current",
                scale=False,
                factor=max(
                    RIGHT_PARENT_ARROW_MIN_LENGTH_UM,
                    diagonal * RIGHT_PARENT_ARROW_SCALE_FRACTION,
                ),
            )
            plotter.add_mesh(
                glyphs,
                color="#FF9F1C",
                label="global parent -> current",
                pickable=False,
            )
            legend_entries.append(
                _ScientificLegendEntry("parent -> current", "#FF9F1C", "arrow")
            )
    if roi.true_terminal_local_ids:
        points = roi.local_node_positions_um[list(roi.true_terminal_local_ids)]
        plotter.add_points(
            points,
            color="#4ADE80",
            point_size=RIGHT_TRUE_TERMINAL_POINT_SIZE_PX,
            render_points_as_spheres=True,
            label="TRUE_TERMINAL",
            pickable=False,
        )
        legend_entries.append(
            _ScientificLegendEntry("true terminal", "#4ADE80", "circle")
        )
    if roi.cut_ports:
        points = np.asarray([port.intersection_position_um for port in roi.cut_ports], dtype=float)
        plotter.add_points(
            points,
            color="#F87171",
            point_size=RIGHT_CUT_PORT_POINT_SIZE_PX,
            render_points_as_spheres=True,
            label="CUT_PORT",
            pickable=False,
        )
        legend_entries.append(
            _ScientificLegendEntry("cut port", "#F87171", "circle")
        )
    radius = roi.radius_features
    structure = roi.structural_features
    plotter.add_text(
        f"{roi.roi_id} | cluster {roi.cluster_id}\n"
        f"nodes {roi.node_count} | branches {roi.branch_count} | bifurcations {roi.bifurcation_count}\n"
        f"radius P10/P25/P50/P75/P90: "
        f"{radius.get('r10', float('nan')):.2f} / {radius.get('r25', float('nan')):.2f} / "
        f"{radius.get('r50', float('nan')):.2f} / {radius.get('r75', float('nan')):.2f} / "
        f"{radius.get('r90', float('nan')):.2f} um\n"
        f"length {structure.get('total_vessel_length_um', roi.retained_component_length_um):.1f} um | "
        f"cycle rank {structure.get('cycle_rank', float('nan')):.0f} | "
        f"true terminals {roi.true_terminal_count} | cut ports {roi.cut_port_count}",
        position="upper_left",
        font_size=10,
        color="white",
        font=UI_FONT_FAMILY,
        name="sampling_roi_information",
    )
    _add_scientific_bottom_legend(
        plotter,
        legend_entries,
        name_prefix="sampling_roi_legend",
    )
    plotter.view_isometric()
    if view_up is not None:
        plotter.camera.SetViewUp(*view_up)
    plotter.reset_camera(bounds=bounds)
    return _add_physical_coordinate_axes(
        plotter,
        bounds,
        add_orientation_axes=add_orientation_axes,
        # 右图已有红色 ROI 边界框；不再叠加外层白色 CubeAxes 框。
        # 但保留由同一 bounds 生成的灰色辅助网格。
        show_reference_frame=False,
        show_reference_grid=True,
        y_vertical=view_up == (0.0, 1.0, 0.0),
    )


def _add_sampling_boxes(
    plotter: Any,
    rois: tuple[ROIRecord, ...],
    indices: list[int],
    *,
    mode_label: str,
) -> tuple[list[Any], dict[str, int]]:
    import pyvista as pv

    actors: list[Any] = []
    actor_to_index: dict[str, int] = {}
    label_points: list[np.ndarray] = []
    labels: list[str] = []
    for index in indices:
        roi = rois[index]
        bounds = (
            roi.bbox_min_um[0], roi.bbox_max_um[0],
            roi.bbox_min_um[1], roi.bbox_max_um[1],
            roi.bbox_min_um[2], roi.bbox_max_um[2],
        )
        color = _sampling_color(roi.cluster_id)
        actor = plotter.add_mesh(
            pv.Box(bounds=bounds),
            color=color,
            opacity=0.16 if roi.is_representative else 0.05,
            show_edges=True,
            edge_color=color,
            line_width=4.0 if roi.is_representative else 1.6,
            pickable=True,
            reset_camera=False,
            name=f"sampling_roi_pick_{index}",
        )
        actors.append(actor)
        actor_to_index[actor.memory_address] = index
        if roi.is_representative or mode_label.startswith("cluster"):
            label_points.append(np.asarray(roi.bbox_center_um, dtype=float))
            labels.append(f"C{roi.cluster_id}\nR{roi.selection_rank if roi.selection_rank > 0 else '-'}")
    if label_points:
        actors.append(
            plotter.add_point_labels(
                np.asarray(label_points),
                labels,
                font_size=9,
                font_family=UI_FONT_FAMILY,
                text_color="white",
                shape_color="#202020",
                shape_opacity=0.70,
                always_visible=True,
                pickable=False,
                reset_camera=False,
                name="sampling_roi_labels",
            )
        )
    actors.append(
        plotter.add_text(
            f"Sampling layer: {mode_label} ({len(indices)} boxes)",
            # 底部区域专门留给方向坐标系和科研图例，模式提示放到右上角，
            # 避免与图例、X/Y 轴数字及轴标题相互遮挡。
            position="upper_right",
            font_size=10,
            color="white",
            font=UI_FONT_FAMILY,
            name="sampling_layer_mode",
        )
    )
    return actors, actor_to_index


def _add_sampling_active_outline(plotter: Any, roi: ROIRecord) -> Any:
    """Highlight the active ROI with a translucent grey fill and white outline."""

    import pyvista as pv

    bounds = (
        roi.bbox_min_um[0], roi.bbox_max_um[0],
        roi.bbox_min_um[1], roi.bbox_max_um[1],
        roi.bbox_min_um[2], roi.bbox_max_um[2],
    )
    return plotter.add_mesh(
        pv.Box(bounds=bounds),
        style="surface",
        color=ACTIVE_ROI_FILL_COLOR,
        opacity=ACTIVE_ROI_FILL_OPACITY,
        lighting=False,
        show_edges=True,
        edge_color=ACTIVE_ROI_EDGE_COLOR,
        edge_opacity=1.0,
        line_width=ACTIVE_ROI_EDGE_WIDTH_PX,
        pickable=False,
        reset_camera=False,
        name="active_sampling_roi_outline",
    )


def _install_sampling_layer(
    plotter: Any,
    rois: tuple[ROIRecord, ...],
    *,
    coordinate_axes_controllers: list[_ScreenAnchoredCoordinateAxes | None] | None = None,
    diameter_clim_um: tuple[float, float] | None = None,
) -> dict[str, Any] | None:
    """Install display-only mode switching and ROI inspection callbacks."""

    if not rois:
        return None
    selected_indices = [index for index, roi in enumerate(rois) if roi.is_representative]
    default_index = min(
        selected_indices or list(range(len(rois))),
        key=lambda index: rois[index].selection_rank if rois[index].selection_rank > 0 else index + 100000,
    )
    cluster_ids = sorted({roi.cluster_id for roi in rois if roi.cluster_id >= 0})
    overlay_actors: list[Any] = []
    actor_to_index: dict[str, int] = {}
    active_outline: list[Any | None] = [None]
    cluster_position = [0]
    interaction_busy = [False]

    def redraw(indices: list[int], mode_label: str) -> None:
        if interaction_busy[0]:
            return
        interaction_busy[0] = True
        try:
            plotter.subplot(0, 0)
            for actor in overlay_actors:
                plotter.remove_actor(actor, reset_camera=False, render=False)
            overlay_actors.clear()
            new_actors, mapping = _add_sampling_boxes(
                plotter,
                rois,
                indices,
                mode_label=mode_label,
            )
            overlay_actors.extend(new_actors)
            actor_to_index.clear()
            actor_to_index.update(mapping)
            plotter.render()
        finally:
            plotter.subplot(0, 0)
            interaction_busy[0] = False

    def select_roi(actor: Any) -> None:
        if interaction_busy[0]:
            return
        index = actor_to_index.get(getattr(actor, "memory_address", ""))
        if index is None:
            return
        interaction_busy[0] = True
        try:
            plotter.subplot(0, 0)
            if active_outline[0] is not None:
                plotter.remove_actor(active_outline[0], reset_camera=False, render=False)
            active_outline[0] = _add_sampling_active_outline(plotter, rois[index])
            plotter.subplot(0, 1)
            if (
                coordinate_axes_controllers is not None
                and len(coordinate_axes_controllers) >= 2
                and coordinate_axes_controllers[1] is not None
            ):
                coordinate_axes_controllers[1].dispose()
                coordinate_axes_controllers[1] = None
            _remove_diameter_colorbar(plotter)
            plotter.renderer.clear_actors()
            replacement_axes = _add_sampling_roi_scene(
                plotter,
                rois[index],
                add_orientation_axes=False,
                diameter_clim_um=diameter_clim_um,
            )
            if (
                coordinate_axes_controllers is not None
                and len(coordinate_axes_controllers) >= 2
            ):
                coordinate_axes_controllers[1] = replacement_axes
            plotter.subplot(0, 0)
            plotter.render()
        finally:
            plotter.subplot(0, 0)
            interaction_busy[0] = False

    def show_all() -> None:
        redraw(list(range(len(rois))), "all candidates")

    def show_selected() -> None:
        redraw(selected_indices or list(range(len(rois))), "selected representatives")

    def show_next_cluster() -> None:
        if not cluster_ids:
            return
        cluster_id = cluster_ids[cluster_position[0] % len(cluster_ids)]
        cluster_position[0] += 1
        redraw(
            [index for index, roi in enumerate(rois) if roi.cluster_id == cluster_id],
            f"cluster {cluster_id}",
        )

    plotter.subplot(0, 0)
    redraw(selected_indices or list(range(len(rois))), "selected representatives")
    active_outline[0] = _add_sampling_active_outline(plotter, rois[default_index])
    key_actions = {
        "a": show_all,
        "r": show_selected,
        "s": show_selected,
        "c": show_next_cluster,
    }
    for key, callback in key_actions.items():
        for variant in (key, key.upper()):
            plotter.clear_events_for_key(variant)
            plotter.add_key_event(variant, callback)
    plotter.enable_mesh_picking(
        callback=select_roi,
        show=False,
        show_message=False,
        left_clicking=True,
        use_actor=True,
    )
    return {
        "show_all": show_all,
        "show_selected": show_selected,
        "show_next_cluster": show_next_cluster,
        "select_roi": select_roi,
        "actor_to_index": actor_to_index,
    }


def render_figure2a_scene(
    volume_zyx: np.ndarray | None,
    geometry: InteractiveSceneGeometry,
    *,
    spacing_xyz_um: tuple[float, float, float],
    sample_id: str,
    volume_opacity: float = 0.32,
    window_size: tuple[int, int] = (1800, 900),
    screenshot_path: Path | None = None,
    show: bool = False,
    sampling_rois: tuple[ROIRecord, ...] = (),
    left_view_up: tuple[float, float, float] | None = None,
) -> dict[str, Any]:
    """Render a static preview or open the native interactive PyVista window."""

    import pyvista as pv

    if sampling_rois:
        plotter = pv.Plotter(
            shape=(1, 2),
            border=True,
            border_color="#606060",
            off_screen=not show,
            window_size=window_size,
        )
        plotter.subplot(0, 0)
    else:
        plotter = pv.Plotter(off_screen=not show, window_size=window_size)
    # 为所有后续由 PyVista 自动创建的数字和英文文本设置 Arial；自定义 VTK
    # 坐标轴、图例、点标签和色条同时在各自创建处再次显式指定 Arial。
    plotter.theme.font.family = UI_FONT_FAMILY
    metadata, left_coordinate_axes = _add_full_scene(
        plotter,
        volume_zyx,
        geometry,
        spacing_xyz_um=spacing_xyz_um,
        volume_opacity=volume_opacity,
        sample_id=sample_id,
        sampling_available=bool(sampling_rois),
        left_view_up=left_view_up,
    )
    # 保存左侧全局场景的独立相机引用。后续切换 subplot 只改变当前 renderer，
    # 不会让这个引用指向右侧相机。
    rotation_cameras: list[Any] = [plotter.camera]
    coordinate_axes_controllers: list[_ScreenAnchoredCoordinateAxes | None] = [
        left_coordinate_axes
    ]
    if sampling_rois:
        global_diameter_clim_um = _sampling_diameter_limits_um(sampling_rois)
        diameter_clim_um = (
            global_diameter_clim_um
            if DIAMETER_COLORBAR_USE_GLOBAL_RANGE
            else None
        )
        selected = [roi for roi in sampling_rois if roi.is_representative]
        default_roi = min(
            selected or list(sampling_rois),
            key=lambda roi: roi.selection_rank if roi.selection_rank > 0 else 100000,
        )
        initial_diameter_clim_um = diameter_clim_um or _sampling_diameter_limits_um(
            (default_roi,)
        )
        plotter.subplot(0, 1)
        right_coordinate_axes = _add_sampling_roi_scene(
            plotter,
            default_roi,
            diameter_clim_um=diameter_clim_um,
        )
        coordinate_axes_controllers.append(right_coordinate_axes)
        # 右侧 ROI 场景拥有自己的焦点和缩放，旋转时只同步方位角增量。
        rotation_cameras.append(plotter.camera)
        plotter.subplot(0, 0)
        metadata["sampling_layer"] = {
            "candidate_count": len(sampling_rois),
            "selected_count": len(selected),
            "cluster_ids": sorted({roi.cluster_id for roi in sampling_rois}),
            "default_roi_id": default_roi.roi_id,
            "diameter_colorbar": {
                "scalar": DIAMETER_SCALAR_NAME,
                "unit": "um",
                "initial_range_um": (
                    list(initial_diameter_clim_um)
                    if initial_diameter_clim_um
                    else None
                ),
                "global_range_um": (
                    list(global_diameter_clim_um)
                    if global_diameter_clim_um
                    else None
                ),
                "colormap": DIAMETER_COLORMAP,
                "outline": False,
                "title_separation_px": DIAMETER_COLORBAR_TITLE_SEPARATION_PX,
                "range_mode": (
                    "all_rois"
                    if DIAMETER_COLORBAR_USE_GLOBAL_RANGE
                    else "current_roi"
                ),
            },
            "display_modes": [
                "all ROI candidates",
                "selected representatives",
                "cluster X",
                "ROI X",
            ],
            "core_recomputed_in_ui": False,
            "right_parent_arrow_max_count": RIGHT_PARENT_ARROW_MAX_COUNT,
            "right_parent_arrow_scale_fraction": RIGHT_PARENT_ARROW_SCALE_FRACTION,
            "right_point_sizes_px": {
                "true_terminal": RIGHT_TRUE_TERMINAL_POINT_SIZE_PX,
                "cut_port": RIGHT_CUT_PORT_POINT_SIZE_PX,
            },
        }
        metadata["available_display_layers"] = ["representative connected ROIs"]
        metadata["default_display_layer"] = "representative connected ROIs"
        if show:
            _install_sampling_layer(
                plotter,
                sampling_rois,
                coordinate_axes_controllers=coordinate_axes_controllers,
                diameter_clim_um=diameter_clim_um,
            )
            metadata["display_switch_keys"] = {
                "R": "selected representative ROIs",
                "A": "all candidate ROIs",
                "S": "selected representative ROIs",
                "C": "next ROI cluster",
            }
            metadata["interaction"] = {
                "selection": "choose an ROI display mode, then left-click a colored box",
                "result": "the right viewport shows the selected connected ROI",
            }
        else:
            indices = [index for index, roi in enumerate(sampling_rois) if roi.is_representative]
            _add_sampling_boxes(
                plotter,
                sampling_rois,
                indices or list(range(len(sampling_rois))),
                mode_label="selected representatives",
            )
            _add_sampling_active_outline(plotter, default_roi)
    rotation_callback = (
        _install_synchronized_horizontal_rotation(
            plotter,
            tuple(rotation_cameras),
            left_view_up=left_view_up,
        )
        if show
        else None
    )
    metadata["automatic_horizontal_rotation"] = {
        "enabled": rotation_callback is not None,
        "viewport_count": len(rotation_cameras),
        "interval_ms": HORIZONTAL_ROTATION_INTERVAL_MS,
        "degrees_per_step": HORIZONTAL_ROTATION_DEGREES_PER_STEP,
        "degrees_per_second": (
            HORIZONTAL_ROTATION_DEGREES_PER_STEP
            * 1000.0
            / HORIZONTAL_ROTATION_INTERVAL_MS
        ),
        "view_up": list(HORIZONTAL_ROTATION_VIEW_UP),
    }
    if left_view_up is not None:
        metadata["automatic_horizontal_rotation"]["left_view_up"] = list(left_view_up)
    metadata["coordinate_axis_layout"] = {
        "x_y_labels": "attached to bottom 3-D axes with continuous edge cross-fading",
        "z_labels": "attached to both outside vertical sides with continuous cross-fading",
        "major_ticks": "screen-projected, perpendicular to each axis, and always outward",
        "numbers": "screen-upright and anchored beyond the outward tick tips",
        "titles": "screen-upright and anchored outside the numerical labels",
        "tick_length_px": COORDINATE_TICK_LENGTH_PX,
        "tick_label_offset_px": COORDINATE_TICK_LABEL_OFFSET_PX,
        "title_offset_px": COORDINATE_TITLE_OFFSET_PX,
        "left_viewport_frame": "explicit bounds wireframe shared with tick projection",
        "left_viewport_grid": "explicit six-face grid generated from the same bounds",
        "right_viewport_frame": "red ROI boundary only; white reference frame disabled",
        "right_viewport_grid": "explicit six-face grid generated from the ROI bounds",
        "legend_layout": (
            "bottom compact row-major line/open-arrow/circle/patch symbols aligned beside "
            "orientation axes"
        ),
        "active_roi_fill": {
            "color": ACTIVE_ROI_FILL_COLOR,
            "opacity": ACTIVE_ROI_FILL_OPACITY,
        },
        "updates_during_mouse_rotation": True,
        "edge_crossfade_px": COORDINATE_EDGE_CROSSFADE_PX,
    }
    if screenshot_path is not None:
        screenshot_path.parent.mkdir(parents=True, exist_ok=True)
    try:
        if show:
            plotter.show(
                title=f"Mouse brain vasculature - {sample_id}",
                screenshot=str(screenshot_path) if screenshot_path else None,
                auto_close=True,
            )
        else:
            plotter.show(
                screenshot=str(screenshot_path) if screenshot_path else None,
                auto_close=True,
            )
    finally:
        for controller in coordinate_axes_controllers:
            if controller is not None:
                controller.dispose()
    return metadata


def save_figure2a_preview(
    result: DirectedVascularGraph,
    image_volume: np.ndarray | None,
    output_dir: Path,
    *,
    spacing_xyz_um: tuple[float, float, float],
    max_arrows: int,
    volume_opacity: float,
    window_size: tuple[int, int],
) -> Figure2aArtifacts:
    """Write acceptance evidence without opening a GUI."""

    screenshot = output_dir / "figure2a_interactive_preview.png"
    manifest = output_dir / "figure2a_scene_manifest.json"
    geometry = scene_geometry_from_graph(
        result,
        max_arrows=max_arrows,
    )
    metadata = render_figure2a_scene(
        image_volume,
        geometry,
        spacing_xyz_um=spacing_xyz_um,
        sample_id=result.sample_id,
        volume_opacity=volume_opacity,
        window_size=window_size,
        screenshot_path=screenshot,
        show=False,
    )
    manifest.write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return Figure2aArtifacts(screenshot, manifest)


def _geometry_from_exports(
    graph_dir: Path,
    max_arrows: int,
) -> InteractiveSceneGeometry:
    archive = np.load(graph_dir / "directed_branch_geometry.npz")
    points = np.asarray(archive["points_um"], dtype=float)
    vectors = np.asarray(archive["direction_parent_to_current_xyz"], dtype=float)
    offsets = np.asarray(archive["branch_offsets"], dtype=np.int64)
    branches = tuple(points[start:end] for start, end in zip(offsets[:-1], offsets[1:]))
    candidates: list[int] = []
    for start, end in zip(offsets[:-1], offsets[1:]):
        if end - start < 2:
            continue
        local = np.linspace(start, end - 1, min(4, end - start), dtype=int)
        candidates.extend(int(value) for value in np.unique(local))
    candidates = [index for index in candidates if np.linalg.norm(vectors[index]) > 0]
    if len(candidates) > max_arrows:
        keep = np.linspace(0, len(candidates) - 1, max_arrows, dtype=int)
        candidates = [candidates[index] for index in keep]

    node_points: list[tuple[float, float, float]] = []
    node_roles: list[str] = []
    with (graph_dir / "directed_nodes.csv").open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            node_points.append((float(row["x_um"]), float(row["y_um"]), float(row["z_um"])))
            node_roles.append(row["role"])
    return InteractiveSceneGeometry(
        branches_um=branches,
        arrow_points_um=points[candidates] if candidates else np.empty((0, 3)),
        arrow_vectors_xyz=vectors[candidates] if candidates else np.empty((0, 3)),
        critical_points_um=np.asarray(node_points, dtype=float).reshape((-1, 3)),
        critical_roles=tuple(node_roles),
    )


def show_saved_run(
    run_root: Path,
    *,
    sample_id: str | None = None,
    max_arrows: int = 600,
    volume_opacity: float = 0.32,
    window_size: tuple[int, int] = (1800, 900),
    sampling_run_root: Path | None = None,
    screenshot_path: Path | None = None,
    show: bool = True,
    left_view_up: tuple[float, float, float] | None = None,
) -> Path:
    """Open the first (or requested) saved sample in a zoomable native window."""

    sample_dirs = sorted(path for path in (run_root / "samples").iterdir() if path.is_dir())
    if sample_id:
        sample_dirs = [path for path in sample_dirs if sample_id in {path.name, path.name.split("__", 1)[-1]}]
    if not sample_dirs:
        raise FileNotFoundError(f"No processed sample found in {run_root}")
    sample_root = sample_dirs[0]
    manifest_path = sample_root / "preprocess_manifest.json"
    graph_dir = sample_root / "graphs"
    if not manifest_path.is_file() or not graph_dir.is_dir():
        raise FileNotFoundError(f"Interactive artifacts are incomplete in {sample_root}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    spacing = tuple(float(value) for value in manifest["spacing_xyz_um"])
    if manifest.get("normalized_volume_path"):
        image_volume, _ = load_normalized_volume(Path(manifest["normalized_volume_path"]))
    elif manifest["record"].get("image_path"):
        image_volume = load_tiff_volume(Path(manifest["record"]["image_path"]))
    elif manifest["record"].get("mask_path"):
        image_volume = load_tiff_volume(Path(manifest["record"]["mask_path"]))
    else:
        image_volume = None
    geometry = _geometry_from_exports(
        graph_dir,
        max_arrows,
    )
    sampling_rois = tuple(
        load_sampling_display_rois(sampling_run_root)
        if sampling_run_root is not None
        else ()
    )
    render_figure2a_scene(
        image_volume,
        geometry,
        spacing_xyz_um=spacing,  # type: ignore[arg-type]
        sample_id=str(manifest["record"]["sample_id"]),
        volume_opacity=volume_opacity,
        window_size=window_size,
        screenshot_path=screenshot_path,
        show=show,
        sampling_rois=sampling_rois,
        left_view_up=left_view_up,
    )
    return sample_root
