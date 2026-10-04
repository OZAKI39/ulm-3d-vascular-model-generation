"""Display-only axes extracted from the requested human-spatial right panel.

Source functions retain the reference frame/grid/tick placement. The wrapper
uses already-transformed candidate coordinates in mm (no /1000 conversion).
The reference project's files are never imported or executed on the server.
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Any
import numpy as np

def apply_font(prop):
    prop.SetFontFamilyToArial()

def add_spatial_reference_axes(plotter, bounds, *, frame_color='#FF5C5C'):
    """Right-panel grid and outward ticks, with an optional frame colour."""
    import pyvista as pv
    box = pv.Box(bounds=bounds)
    plotter.add_mesh(box, color=frame_color, opacity=.025, pickable=False, reset_camera=False)
    plotter.add_mesh(box, color=frame_color, style='wireframe', line_width=4.,
                     pickable=False, reset_camera=False)
    controller = _add_physical_coordinate_axes(plotter, tuple(bounds),
        show_reference_frame=False, show_reference_grid=True, y_vertical=False)
    controller.single_edge_per_axis = True
    for index, group in enumerate((controller.x_annotations, controller.y_annotations, controller.z_annotations)):
        for annotation in group:
            values = np.linspace(annotation.point1_xyz[index], annotation.point2_xyz[index], len(annotation.tick_actors))
            for value, actor in zip(values, annotation.tick_actors):
                actor.SetInput(f'{value:.2f}')
                actor.SetTextScaleModeToNone()
                actor.GetTextProperty().SetFontSize(14)
                apply_font(actor.GetTextProperty())
            annotation.title_actor.SetInput(f'{annotation.axis_name} (mm)')
            annotation.title_actor.GetTextProperty().SetFontSize(16)
            apply_font(annotation.title_actor.GetTextProperty())
            annotation.tick_label_offset_px = 20.
            # The preserved camera is Z-up; leave a small gap between its
            # vertical title and a middle tick number at the same height.
            annotation.title_offset_px = 50. if annotation.axis_name == 'Z' else 32.
            annotation.vertical_axis_name = 'Z'
            annotation.endpoint_label_inset_px = 28.
    controller.update()
    labels = StableOrientationLabels(plotter.renderer.axes_widget)
    labels.axes.SetTotalLength(1.4, .8, 1.)
    labels.update()
    return controller, labels

def axis_annotations(controller):
    """Actual visible text extents in top-left image coordinates for safe fit."""
    from vtkmodules.vtkRenderingCore import vtkTextRenderer
    height = controller.renderer.GetSize()[1]
    dpi = controller.renderer.GetRenderWindow().GetDPI()
    result = []
    for group in (controller.x_annotations, controller.y_annotations, controller.z_annotations):
        for annotation in group:
            for actor in (*annotation.tick_actors, annotation.title_actor):
                if not actor.GetVisibility():
                    continue
                bbox = [0,0,0,0]
                vtkTextRenderer.GetInstance().GetBoundingBox(actor.GetTextProperty(), actor.GetInput(), bbox, dpi)
                x,y = actor.GetPosition()
                result.append(dict(text=actor.GetInput(), axis=annotation.axis_name,
                    bbox_px=[x+bbox[0],height-(y+bbox[3]),x+bbox[1],height-(y+bbox[2])],
                    edge_start_mm=list(annotation.point1_xyz),edge_end_mm=list(annotation.point2_xyz)))
    return result

COORDINATE_TICK_FONT_SIZE = 11

COORDINATE_TITLE_FONT_SIZE = 12

ORIENTATION_AXIS_FONT_SIZE = 12

COORDINATE_AXIS_COLOR = "#BFC7D5"

COORDINATE_FRAME_LINE_WIDTH_PX = 1.2

COORDINATE_GRID_LINE_WIDTH_PX = 1.0

COORDINATE_GRID_OPACITY = 0.28

COORDINATE_TICK_LENGTH_PX = 7.0

COORDINATE_TICK_LINE_WIDTH_PX = 1.5

COORDINATE_TICK_LABEL_OFFSET_PX = 10.0

COORDINATE_TITLE_OFFSET_PX = 30.0

COORDINATE_X_TICK_COUNT = 4

COORDINATE_Y_TICK_COUNT = 4

COORDINATE_Z_TICK_COUNT = 5

COORDINATE_EDGE_CROSSFADE_PX = 36.0

COORDINATE_TEXT_VISIBILITY_EPSILON = 0.01

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

class StableOrientationLabels:
    """Fixed-size, centre-anchored letters moving continuously with the triad.

    The original widget and its coloured 3-D axes are retained. vtkCaptionActor2D
    auto-fits and reanchors text as the axes rotate; replacing only those captions
    with projected text removes the resulting size/alignment jumps.
    """
    def __init__(self, widget):
        from vtkmodules.vtkRenderingCore import vtkTextActor

        self.widget = widget
        self.renderer = widget.GetRenderer()
        self.axes = widget.GetOrientationMarker()
        self.camera = self.renderer.GetActiveCamera()
        self.axes.AxisLabelsOff()
        self.actors = []
        for letter in 'XYZ':
            actor = vtkTextActor()
            actor.SetInput(f'{letter} (mm)')
            actor.SetTextScaleModeToNone()
            actor.GetPositionCoordinate().SetCoordinateSystemToDisplay()
            prop = actor.GetTextProperty()
            apply_font(prop)
            prop.SetFontSize(ORIENTATION_AXIS_FONT_SIZE)
            prop.SetBold(True)
            prop.SetColor(1., 1., 1.)
            prop.SetJustificationToCentered()
            prop.SetVerticalJustificationToCentered()
            # VTK 9.7 removed AddActor2D; AddViewProp retains this 2-D actor.
            self.renderer.AddViewProp(actor)
            self.actors.append(actor)
        self.observers = [(self.camera, self.camera.AddObserver('ModifiedEvent', self.update)),
                          (self.renderer, self.renderer.AddObserver('StartEvent', self.update))]
        self.update()

    def update(self, *_):
        lengths = self.axes.GetTotalLength()
        for index, actor in enumerate(self.actors):
            endpoint = [0., 0., 0.]
            endpoint[index] = lengths[index] * 1.20
            self.renderer.SetWorldPoint(*endpoint, 1.)
            self.renderer.WorldToDisplay()
            x, y, _ = self.renderer.GetDisplayPoint()
            if np.isfinite(x) and np.isfinite(y):
                actor.SetPosition(float(x), float(y))

    def dispose(self):
        for source, observer in self.observers:
            source.RemoveObserver(observer)
        self.observers.clear()
        for actor in self.actors:
            self.renderer.RemoveViewProp(actor)
