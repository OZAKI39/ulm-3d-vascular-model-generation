from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from PIL import Image

from utils.rodent_vasculature import interactive as interactive_module
from utils.rodent_vasculature.config import RodentVasculatureConfig
from utils.rodent_vasculature.graph_builder import build_directed_vascular_graph
from utils.rodent_vasculature.pipeline import run_rodent_vasculature_pipeline
from utils.rodent_vasculature.swc_analysis import (
    evaluate_optional_mask_qc,
    select_analysis_swc,
)
from utils.rodent_vasculature.swc_io import load_normalized_swc, load_swc, save_normalized_swc
from utils.rodent_vasculature.tiff_io import load_tiff_volume
from utils.rodent_vasculature.validation import evaluate_directed_graph


SWC_TEXT = """# id type x y z radius parent
1.0 2 0 1 1 2 -1
2.0 2 1 1 1 2 1.0
3.0 2 2 1 1 2 2.0
4.0 2 2 2 1 1 3.0
5.0 2 2 0 1 1 3.0
6.0 2 3 2 1 1 4.0
"""


def _config(tmp_path: Path, **overrides: object) -> RodentVasculatureConfig:
    values = {
        "input_dir": tmp_path,
        "output_root": tmp_path / "outputs",
        "stage": "all",
        "expected_shape_zyx": (4, 4, 4),
        "spacing_xyz_um": (1.0, 1.0, 2.0),
        "save_vtp": False,
    }
    values.update(overrides)
    return RodentVasculatureConfig(**values)  # type: ignore[arg-type]


def test_parent_is_upstream_and_branch_order_is_preserved(tmp_path: Path) -> None:
    path = tmp_path / "tree.swc"
    path.write_text(SWC_TEXT, encoding="utf-8")
    swc = load_swc(path, spacing_xyz_um=(1, 1, 2), volume_shape_zyx=(4, 4, 4))
    result = build_directed_vascular_graph("synthetic", swc, _config(tmp_path))

    assert set(result.source_graph.edges) == {(1, 2), (2, 3), (3, 4), (3, 5), (4, 6)}
    sequences = {tuple(branch.source_node_ids) for branch in result.branches}
    assert sequences == {(1, 2, 3), (3, 4, 6), (3, 5)}
    root_branch = next(branch for branch in result.branches if branch.source_node_ids == [1, 2, 3])
    assert root_branch.downstream_terminal_count == 2
    assert root_branch.strahler_order == 2
    assert len(root_branch.daughter_branch_ids) == 2
    acceptance = evaluate_directed_graph(result, [], strict_nonpositive_radius=False)
    assert acceptance.overall_status == "PASS"


def test_nonpositive_radius_is_preserved_but_warned(tmp_path: Path) -> None:
    path = tmp_path / "radius.swc"
    path.write_text(SWC_TEXT.replace("3.0 2 2 1 1 2", "3.0 2 2 1 1 0"), encoding="utf-8")
    swc = load_swc(path, spacing_xyz_um=(1, 1, 2), volume_shape_zyx=(4, 4, 4))
    result = build_directed_vascular_graph("synthetic", swc, _config(tmp_path))
    assert swc.radius_raw_um[2] == 0
    assert any("raw values were preserved" in warning for warning in result.warnings)
    acceptance = evaluate_directed_graph(result, [], strict_nonpositive_radius=False)
    assert acceptance.overall_status == "WARNING"


def test_swc_centric_selection_preserves_reference_and_never_repairs_from_mask(
    tmp_path: Path,
) -> None:
    path = tmp_path / "two_components.swc"
    path.write_text(
        "1 0 1 1 1 1 -1\n"
        "2 0 2 1 1 1 1\n"
        "3 0 4 1 1 1 -1\n"
        "4 0 5 1 1 1 3\n"
        "5 0 6 1 1 1 4\n",
        encoding="utf-8",
    )
    mask = np.zeros((3, 3, 8), dtype=np.uint8)
    mask[1, 1, 1:7] = 255
    original = load_swc(path, spacing_xyz_um=(1, 1, 1), volume_shape_zyx=mask.shape)
    result = select_analysis_swc(
        original,
        spacing_xyz_um=(1, 1, 1),
        volume_shape_zyx=mask.shape,
    )
    mask_qc = evaluate_optional_mask_qc(mask, original, result.analysis_swc)

    assert original.component_count == 2
    assert result.reference_swc.node_ids.tolist() == [1, 2, 3, 4, 5]
    assert result.analysis_swc.node_ids.tolist() == [3, 4, 5]
    assert result.analysis_swc.parent_ids.tolist() == [-1, 3, 4]
    assert result.reference_only_node_ids.tolist() == [1, 2]
    assert result.summary["new_node_count"] == 0
    assert result.summary["new_edge_count"] == 0
    assert result.summary["parent_relation_change_count"] == 0
    assert result.summary["reference_only_components_are_errors"] is False
    assert mask_qc["used_for_component_selection"] is False
    assert mask_qc["used_for_topology_repair"] is False

    normalized = save_normalized_swc(result.analysis_swc, tmp_path / "analysis.npz")
    reloaded = load_normalized_swc(
        normalized,
        source_path=path,
        spacing_xyz_um=(1, 1, 1),
        volume_shape_zyx=mask.shape,
    )
    assert reloaded.component_count == 1
    assert reloaded.node_ids.tolist() == result.analysis_swc.node_ids.tolist()


def _write_multipage_tiff(path: Path, volume: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frames = [Image.fromarray(frame.astype(np.uint8)) for frame in volume]
    frames[0].save(path, save_all=True, append_images=frames[1:], compression="tiff_lzw")


def _make_dataset(tmp_path: Path) -> None:
    root = tmp_path / "raw_data" / "analysis_data" / "analysis_data"
    stem = "mouse_0_0_0"
    volume = np.zeros((4, 4, 4), dtype=np.uint8)
    for x, y, z in ((0, 1, 1), (1, 1, 1), (2, 1, 1), (2, 2, 1), (2, 0, 1), (3, 2, 1)):
        volume[z, y, x] = 200
    _write_multipage_tiff(root / "images" / f"{stem}.tif", volume)
    _write_multipage_tiff(root / "mask" / f"{stem}.tif", volume)
    (root / "swc").mkdir(parents=True)
    (root / "swc" / f"{stem}.swc").write_text(SWC_TEXT, encoding="utf-8")


def test_all_stage_writes_arrow_visualizations_and_reports(tmp_path: Path) -> None:
    _make_dataset(tmp_path)

    run = run_rodent_vasculature_pipeline(
        _config(
            tmp_path,
            input_dir=tmp_path,
            cohort="raw-analysis",
            max_samples=1,
            figure2a_enabled=True,
        ),
        verbose=False,
    )
    assert run.status == "completed"
    sample_root = next((run.run_root / "samples").iterdir())
    assert (sample_root / "visualizations" / "direction_parent_to_current_3d.png").is_file()
    assert (sample_root / "visualizations" / "direction_parent_to_current_orthogonal.png").is_file()
    assert (sample_root / "visualizations" / "directed_branch_topology_xy.png").is_file()
    assert (
        sample_root / "visualizations" / "figure2a_interactive_preview.png"
    ).is_file()
    scene_manifest = sample_root / "visualizations" / "figure2a_scene_manifest.json"
    assert scene_manifest.is_file()
    scene = json.loads(scene_manifest.read_text(encoding="utf-8"))
    assert scene["direction_rule"] == "SWC parent_id node -> current node"
    assert scene["direction_is_measured_flow"] is False
    assert scene["global_coordinate_bounds_xyz_um"] == [0.0, 3.0, 0.0, 3.0, 0.0, 6.0]
    assert scene["coordinate_units"] == "um"
    assert scene["interaction"]["selection"].startswith("rotate")
    graph_acceptance = json.loads(
        (sample_root / "graph_acceptance.json").read_text(encoding="utf-8")
    )
    check_names = {check["name"] for check in graph_acceptance["checks"]}
    assert all("tree" not in name.lower() for name in check_names)
    assert graph_acceptance["overall_status"] == "PASS"
    assert (sample_root / "graphs" / "source_parent_to_current_edges.csv").is_file()
    report = (sample_root / "acceptance_report.html").read_text(encoding="utf-8")
    assert "parent_id node → current node" in report
    assert "figure2a_interactive_preview.png" in report


def test_swc_only_sample_runs_without_optional_image_or_mask(tmp_path: Path) -> None:
    root = tmp_path / "raw_data" / "analysis_data" / "analysis_data" / "swc"
    root.mkdir(parents=True)
    (root / "mouse_0_0_0.swc").write_text(SWC_TEXT, encoding="utf-8")

    run = run_rodent_vasculature_pipeline(
        _config(
            tmp_path,
            input_dir=tmp_path,
            cohort="raw-analysis",
            max_samples=1,
            figure2a_enabled=True,
        )
    )

    assert run.status == "completed"
    sample_root = next((run.run_root / "samples").iterdir())
    manifest = json.loads((sample_root / "preprocess_manifest.json").read_text("utf-8"))
    assert manifest["record"]["eligible"] is True
    assert manifest["record"]["image_path"] is None
    assert manifest["record"]["mask_path"] is None
    assert manifest["normalized_volume_path"] is None
    assert manifest["swc_centric_preprocessing"]["mask_qc"]["available"] is False
    assert manifest["swc_centric_preprocessing"]["new_node_count"] == 0
    assert (sample_root / "graphs" / "source_parent_to_current_edges.csv").is_file()
    scene_path = sample_root / "visualizations" / "figure2a_scene_manifest.json"
    assert scene_path.is_file()
    scene = json.loads(scene_path.read_text("utf-8"))
    assert scene["optional_background_volume_available"] is False


def test_preprocess_and_hierarchical_graph_can_run_separately(tmp_path: Path) -> None:
    _make_dataset(tmp_path)
    preprocess = run_rodent_vasculature_pipeline(
        _config(
            tmp_path,
            input_dir=tmp_path,
            stage="preprocess",
            cohort="raw-analysis",
            max_samples=1,
            visualizations_enabled=False,
        )
    )
    assert preprocess.status == "completed"
    graph = run_rodent_vasculature_pipeline(
        _config(
            tmp_path,
            stage="hierarchical-graph",
            source_run=preprocess.run_root,
            cohort="raw-analysis",
            max_samples=1,
            visualizations_enabled=False,
        )
    )
    assert graph.status == "completed"
    sample_root = next((graph.run_root / "samples").iterdir())
    assert (sample_root / "graphs" / "branch_hierarchy_parent_to_current.graphml").is_file()


def test_pipeline_preserves_reference_and_selects_analysis_component_without_mask_cleanup(
    tmp_path: Path,
) -> None:
    _make_dataset(tmp_path)
    source_swc = next(tmp_path.rglob("*.swc"))
    source_swc.write_text(
        SWC_TEXT + "7 2 3 3 3 1 -1\n8 2 3.1 3 3 1 7\n",
        encoding="utf-8",
    )
    mask_path = next((tmp_path / "raw_data").rglob("mask/*.tif"))
    image_path = next((tmp_path / "raw_data").rglob("images/*.tif"))
    mask = load_tiff_volume(mask_path)
    image = load_tiff_volume(image_path)
    mask[3, 3, 3] = 200
    image[3, 3, 3] = 200
    _write_multipage_tiff(mask_path, mask)
    _write_multipage_tiff(image_path, image)
    run = run_rodent_vasculature_pipeline(
        _config(
            tmp_path,
            input_dir=tmp_path,
            cohort="raw-analysis",
            max_samples=1,
            visualizations_enabled=False,
        )
    )

    assert run.status == "completed"
    sample_root = next((run.run_root / "samples").iterdir())
    manifest = json.loads((sample_root / "preprocess_manifest.json").read_text("utf-8"))
    graph_summary = json.loads((sample_root / "graph_summary.json").read_text("utf-8"))
    summary = manifest["swc_centric_preprocessing"]
    assert manifest["swc_reference"]["component_count"] == 2
    assert manifest["swc_analysis"]["component_count"] == 1
    assert summary["reference_only_node_count"] == 2
    assert summary["reference_only_components_are_errors"] is False
    assert summary["new_node_count"] == 0
    assert summary["new_edge_count"] == 0
    assert summary["mask_qc"]["component_count_26"] == 2
    assert summary["mask_qc"]["used_for_component_selection"] is False
    assert graph_summary["source_node_count"] == 6
    assert graph_summary["component_count"] == 1
    assert (
        sample_root
        / "swc_analysis_preprocessing"
        / "reference_swc_components.csv"
    ).is_file()
    assert (
        sample_root
        / "swc_analysis_preprocessing"
        / "reference_only_node_ids.csv"
    ).is_file()


class _FakeRenderer:
    def __init__(self) -> None:
        self.clear_count = 0

    def clear_actors(self) -> None:
        self.clear_count += 1


class _FakeButtonRepresentation:
    def __init__(self, state: bool) -> None:
        self.state = int(state)

    def SetState(self, state: int) -> None:
        self.state = int(state)

    def GetState(self) -> int:
        return self.state


class _FakeButtonWidget:
    def __init__(self, state: bool) -> None:
        self.representation = _FakeButtonRepresentation(state)

    def GetRepresentation(self) -> _FakeButtonRepresentation:
        return self.representation


class _FakePyVistaInteractor:
    def __init__(self) -> None:
        self.initialized = 0
        self.initialize_count = 0
        self._timer = SimpleNamespace(id=None)

    def initialize(self) -> None:
        self.initialized = 1
        self.initialize_count += 1


class _FakeInteractivePlotter:
    def __init__(self) -> None:
        self.renderer = _FakeRenderer()
        self.iren = _FakePyVistaInteractor()
        self.render_count = 0
        self.key_events: dict[str, object] = {}
        self.picking_callbacks: list[object] = []
        self.radio_buttons: list[_FakeButtonWidget] = []
        self.window_size = (1800, 900)
        self.text_entries: dict[str, tuple[str, object]] = {}
        self.timer_events: list[dict[str, object]] = []

    def subplot(self, _row: int, _column: int) -> None:
        return None

    def remove_actor(self, *_args: object, **_kwargs: object) -> None:
        return None

    def render(self) -> None:
        self.render_count += 1

    def add_key_event(self, key: str, callback: object) -> None:
        self.key_events[key] = callback

    def clear_events_for_key(self, key: str) -> None:
        self.key_events.pop(key, None)

    def add_radio_button_widget(
        self,
        _callback: object,
        _group: str,
        *,
        value: bool = False,
        **_kwargs: object,
    ) -> _FakeButtonWidget:
        widget = _FakeButtonWidget(value)
        self.radio_buttons.append(widget)
        return widget

    def add_text(
        self,
        value: str,
        *_args: object,
        name: str | None = None,
        color: object = None,
        **_kwargs: object,
    ) -> object:
        if name is not None:
            self.text_entries[name] = (value, color)
        return object()

    def enable_mesh_picking(self, callback: object, **_kwargs: object) -> None:
        self.picking_callbacks.append(callback)

    def add_timer_event(
        self,
        max_steps: int,
        duration: int,
        callback: object,
    ) -> None:
        # VTK cannot create a valid repeating timer until the interactor has
        # been initialized. Model that ordering constraint in the unit test.
        assert self.iren.initialized
        self.iren._timer = SimpleNamespace(id=2)
        self.timer_events.append(
            {
                "max_steps": max_steps,
                "duration": duration,
                "callback": callback,
            }
        )


class _FakeCamera:
    def __init__(self) -> None:
        self.view_up_calls: list[tuple[float, float, float]] = []
        self.azimuth_calls: list[float] = []

    def SetViewUp(self, x: float, y: float, z: float) -> None:
        self.view_up_calls.append((float(x), float(y), float(z)))

    def Azimuth(self, angle: float) -> None:
        self.azimuth_calls.append(float(angle))


class _FakeObservableCamera:
    def __init__(self) -> None:
        self.callback: object | None = None
        self.removed_observer_ids: list[int] = []

    def AddObserver(self, _event: str, callback: object) -> int:
        self.callback = callback
        return 17

    def RemoveObserver(self, observer_id: int) -> None:
        self.removed_observer_ids.append(observer_id)


class _FakeProjectionRenderer:
    def __init__(
        self,
        projected_xy_by_world_point: dict[tuple[float, float, float], tuple[float, float]],
    ) -> None:
        self.projected_xy_by_world_point = projected_xy_by_world_point
        self.world_point = (0.0, 0.0, 0.0)

    def SetWorldPoint(self, x: float, y: float, z: float, _w: float) -> None:
        self.world_point = (float(x), float(y), float(z))

    def WorldToDisplay(self) -> None:
        return None

    def GetDisplayPoint(self) -> tuple[float, float, float]:
        x, y = self.projected_xy_by_world_point[self.world_point]
        return x, y, 0.0

    def GetOrigin(self) -> tuple[int, int]:
        return 0, 0


class _FakeProjectedAxisAnnotation:
    def __init__(
        self,
        point1_xyz: tuple[float, float, float],
        point2_xyz: tuple[float, float, float],
    ) -> None:
        self.point1_xyz = point1_xyz
        self.point2_xyz = point2_xyz
        self.updates: list[tuple[np.ndarray, np.ndarray, np.ndarray, float]] = []

    def update_from_projection(
        self,
        point1_xy: np.ndarray,
        point2_xy: np.ndarray,
        box_center_xy: np.ndarray,
        opacity: float,
    ) -> None:
        self.updates.append((point1_xy, point2_xy, box_center_xy, float(opacity)))


class _FakeTextProperty:
    def __init__(self) -> None:
        self.family = ""
        self.font_size = 0
        self.bold = False
        self.opacity = 1.0
        self.orientation = 0.0
        self.justification = "centered"
        self.vertical_justification = "centered"

    def SetFontFamilyToArial(self) -> None:
        self.family = "Arial"

    def SetFontSize(self, font_size: int) -> None:
        self.font_size = int(font_size)

    def SetBold(self, bold: bool) -> None:
        self.bold = bool(bold)

    def SetOpacity(self, opacity: float) -> None:
        self.opacity = float(opacity)

    def SetOrientation(self, orientation: float) -> None:
        self.orientation = float(orientation)

    def SetJustificationToCentered(self) -> None:
        self.justification = "centered"

    def SetJustificationToLeft(self) -> None:
        self.justification = "left"

    def SetJustificationToRight(self) -> None:
        self.justification = "right"

    def SetVerticalJustificationToCentered(self) -> None:
        self.vertical_justification = "centered"

    def SetVerticalJustificationToTop(self) -> None:
        self.vertical_justification = "top"


class _FakeOverlayTextActor:
    def __init__(self) -> None:
        self.text_property = _FakeTextProperty()
        self.position = (0.0, 0.0)
        self.visible = False

    def SetPosition(self, x: float, y: float) -> None:
        self.position = (float(x), float(y))

    def SetVisibility(self, visible: bool) -> None:
        self.visible = bool(visible)

    def GetTextProperty(self) -> _FakeTextProperty:
        return self.text_property


class _FakeTickLinePoints:
    def __init__(self, count: int) -> None:
        self.values = [(0.0, 0.0, 0.0)] * (count * 2)
        self.modified_count = 0

    def SetPoint(self, index: int, x: float, y: float, z: float) -> None:
        self.values[index] = (float(x), float(y), float(z))

    def Modified(self) -> None:
        self.modified_count += 1


class _FakeTickLineProperty:
    def __init__(self) -> None:
        self.opacity = 1.0

    def SetOpacity(self, opacity: float) -> None:
        self.opacity = float(opacity)


class _FakeTickLineActor:
    def __init__(self) -> None:
        self.property = _FakeTickLineProperty()
        self.visible = False

    def SetVisibility(self, visible: bool) -> None:
        self.visible = bool(visible)

    def GetProperty(self) -> _FakeTickLineProperty:
        return self.property


class _FakeCaptionActor:
    def __init__(self, text_property: _FakeTextProperty) -> None:
        self.text_property = text_property

    def GetCaptionTextProperty(self) -> _FakeTextProperty:
        return self.text_property


class _FakeAxesActor:
    def __init__(self) -> None:
        self.properties = [_FakeTextProperty() for _ in range(3)]
        self.captions = [_FakeCaptionActor(prop) for prop in self.properties]

    def GetXAxisCaptionActor2D(self) -> _FakeCaptionActor:
        return self.captions[0]

    def GetYAxisCaptionActor2D(self) -> _FakeCaptionActor:
        return self.captions[1]

    def GetZAxisCaptionActor2D(self) -> _FakeCaptionActor:
        return self.captions[2]


class _FakeBoundsActor:
    def __init__(self) -> None:
        self.labels = [_FakeTextProperty() for _ in range(3)]
        self.titles = [_FakeTextProperty() for _ in range(3)]

    def GetLabelTextProperty(self, index: int) -> _FakeTextProperty:
        return self.labels[index]

    def GetTitleTextProperty(self, index: int) -> _FakeTextProperty:
        return self.titles[index]


class _FakeLegendActor:
    def __init__(self) -> None:
        self.text_property = _FakeTextProperty()

    def GetEntryTextProperty(self) -> _FakeTextProperty:
        return self.text_property


def test_interactive_typography_uses_arial_with_rebalanced_sizes() -> None:
    axes = _FakeAxesActor()
    bounds = _FakeBoundsActor()
    legend = _FakeLegendActor()

    interactive_module._style_orientation_axes(axes)
    interactive_module._style_bounds_axes(bounds)
    interactive_module._style_legend(legend)

    assert interactive_module.UI_FONT_FAMILY == "arial"
    assert all(prop.family == "Arial" for prop in axes.properties)
    assert all(
        prop.font_size == interactive_module.ORIENTATION_AXIS_FONT_SIZE
        for prop in axes.properties
    )
    assert all(prop.family == "Arial" for prop in bounds.labels + bounds.titles)
    assert all(
        prop.font_size == interactive_module.COORDINATE_TICK_FONT_SIZE
        for prop in bounds.labels
    )
    assert interactive_module.COORDINATE_TICK_FONT_SIZE > 9
    assert legend.text_property.family == "Arial"
    assert legend.text_property.font_size == interactive_module.LEGEND_FONT_SIZE
    assert interactive_module.LEGEND_FONT_SIZE <= interactive_module.COORDINATE_TITLE_FONT_SIZE


def test_bottom_legend_uses_line_circle_and_patch_geometry_with_compact_rows() -> None:
    import pyvista as pv

    plotter = pv.Plotter(off_screen=True, window_size=(600, 600))
    try:
        entries = [
            interactive_module._ScientificLegendEntry("line", "#00FFFF", "line"),
            interactive_module._ScientificLegendEntry("point", "#00FF00", "circle"),
            interactive_module._ScientificLegendEntry("region", "#AAAAAA", "patch"),
            interactive_module._ScientificLegendEntry("direction", "#FFFFFF", "arrow"),
        ]
        actors = interactive_module._add_scientific_bottom_legend(
            plotter,
            entries,
            name_prefix="test_scientific_legend",
        )

        line_data = actors[1].GetMapper().GetInput()
        circle_data = actors[3].GetMapper().GetInput()
        patch_data = actors[5].GetMapper().GetInput()
        arrow_data = actors[7].GetMapper().GetInput()
        assert line_data.GetNumberOfLines() == 1
        assert line_data.GetNumberOfPolys() == 0
        assert circle_data.GetNumberOfPolys() == 1
        assert circle_data.GetNumberOfPoints() == interactive_module.LEGEND_CIRCLE_RESOLUTION
        assert patch_data.GetNumberOfPolys() == 1
        assert patch_data.GetNumberOfPoints() == 4
        assert arrow_data.GetNumberOfLines() == 3
        assert arrow_data.GetNumberOfPoints() == 4
        first_row_y = actors[2].GetPosition()[1]
        second_row_y = actors[8].GetPosition()[1]
        assert np.isclose(
            first_row_y - second_row_y,
            interactive_module.LEGEND_ROW_SPACING_FRACTION,
        )
    finally:
        plotter.close()


def test_sampling_diameter_limits_are_twice_the_positive_radius_values() -> None:
    rois = (
        SimpleNamespace(
            local_edge_radius_um=(
                np.asarray((1.0, 2.5)),
                np.asarray((np.nan, -3.0)),
            )
        ),
        SimpleNamespace(local_edge_radius_um=(np.asarray((0.75, 4.0)),)),
    )

    assert interactive_module._sampling_diameter_limits_um(rois) == (1.5, 8.0)


def test_right_roi_arrow_subsampling_is_sparse_and_spatially_distributed() -> None:
    points = np.column_stack(
        (
            np.linspace(0.0, 99.0, 100),
            np.zeros(100),
            np.zeros(100),
        )
    )

    selected = interactive_module._spatially_spread_indices(
        points,
        maximum_count=interactive_module.RIGHT_PARENT_ARROW_MAX_COUNT,
    )

    assert len(selected) == interactive_module.RIGHT_PARENT_ARROW_MAX_COUNT
    assert len(np.unique(selected)) == len(selected)
    assert 0 in selected
    assert 99 in selected


def test_horizontal_rotation_timer_updates_both_viewport_cameras() -> None:
    """The left and right cameras receive the same horizontal rotation step."""

    plotter = _FakeInteractivePlotter()
    left_camera = _FakeCamera()
    right_camera = _FakeCamera()

    callback = interactive_module._install_synchronized_horizontal_rotation(
        plotter,
        (left_camera, right_camera),
    )

    assert callback is not None
    assert plotter.iren.initialize_count == 1
    assert plotter.iren._timer.id == 2
    assert len(plotter.timer_events) == 1
    assert plotter.timer_events[0]["duration"] == (
        interactive_module.HORIZONTAL_ROTATION_INTERVAL_MS
    )
    assert plotter.timer_events[0]["max_steps"] == (
        interactive_module.HORIZONTAL_ROTATION_MAX_STEPS
    )

    callback(0)
    expected_view_up = tuple(
        float(value) for value in interactive_module.HORIZONTAL_ROTATION_VIEW_UP
    )
    expected_step = interactive_module.HORIZONTAL_ROTATION_DEGREES_PER_STEP
    assert left_camera.view_up_calls == [expected_view_up]
    assert right_camera.view_up_calls == [expected_view_up]
    assert left_camera.azimuth_calls == [expected_step]
    assert right_camera.azimuth_calls == [expected_step]


def test_projected_axis_ticks_and_text_stay_outside_and_upright() -> None:
    tick_actors = tuple(_FakeOverlayTextActor() for _ in range(3))
    tick_line_points = _FakeTickLinePoints(len(tick_actors))
    tick_line_actor = _FakeTickLineActor()
    title_actor = _FakeOverlayTextActor()
    annotation = interactive_module._ProjectedAxisAnnotation(
        axis_name="X",
        point1_xyz=(0.0, 0.0, 0.0),
        point2_xyz=(1.0, 0.0, 0.0),
        tick_line_points=tick_line_points,
        tick_line_actor=tick_line_actor,
        tick_actors=tick_actors,
        title_actor=title_actor,
    )

    annotation.update_from_projection(
        np.asarray((10.0, 50.0)),
        np.asarray((110.0, 50.0)),
        np.asarray((60.0, 100.0)),
        0.75,
    )
    expected_tick_y = 50.0 - interactive_module.COORDINATE_TICK_LABEL_OFFSET_PX
    expected_title_y = 50.0 - (
        interactive_module.COORDINATE_TICK_LABEL_OFFSET_PX
        + interactive_module.COORDINATE_TITLE_OFFSET_PX
    )
    assert all(np.isclose(actor.position[1], expected_tick_y) for actor in tick_actors)
    assert np.isclose(title_actor.position[1], expected_title_y)
    assert all(actor.visible for actor in (*tick_actors, title_actor))
    assert all(
        np.isclose(actor.text_property.opacity, 0.75)
        for actor in (*tick_actors, title_actor)
    )
    assert tick_line_actor.visible
    assert np.isclose(tick_line_actor.property.opacity, 0.75)
    assert tick_line_points.modified_count == 1
    expected_tick_ends = (
        (10.0, 50.0 - interactive_module.COORDINATE_TICK_LENGTH_PX, 0.0),
        (60.0, 50.0 - interactive_module.COORDINATE_TICK_LENGTH_PX, 0.0),
        (110.0, 50.0 - interactive_module.COORDINATE_TICK_LENGTH_PX, 0.0),
    )
    expected_tick_starts = (
        (10.0, 50.0, 0.0),
        (60.0, 50.0, 0.0),
        (110.0, 50.0, 0.0),
    )
    assert all(
        np.allclose(tick_line_points.values[index * 2], expected)
        for index, expected in enumerate(expected_tick_starts)
    )
    assert all(
        np.allclose(tick_line_points.values[index * 2 + 1], expected)
        for index, expected in enumerate(expected_tick_ends)
    )
    projected_axis = np.asarray((100.0, 0.0))
    for tick_index in range(len(tick_actors)):
        start = np.asarray(tick_line_points.values[tick_index * 2][:2])
        end = np.asarray(tick_line_points.values[tick_index * 2 + 1][:2])
        assert np.isclose(np.dot(end - start, projected_axis), 0.0)
    assert all(actor.text_property.orientation == 0.0 for actor in tick_actors)
    assert title_actor.text_property.orientation == 0.0
    assert all(
        actor.text_property.vertical_justification == "top"
        for actor in (*tick_actors, title_actor)
    )


def test_screen_anchored_axes_follow_bottom_and_both_side_edges() -> None:
    """X/Y favour bottom edges and Z remains visible on both outside edges."""

    bounds = (0.0, 1.0, 0.0, 1.0, 0.0, 2.0)
    projected = {
        (0.0, 0.0, 0.0): (0.0, 20.0),
        (0.0, 1.0, 0.0): (80.0, 100.0),
        (1.0, 0.0, 0.0): (100.0, 120.0),
        (1.0, 1.0, 0.0): (180.0, 200.0),
        (0.0, 0.0, 2.0): (0.0, 60.0),
        (0.0, 1.0, 2.0): (80.0, 140.0),
        (1.0, 0.0, 2.0): (100.0, 160.0),
        (1.0, 1.0, 2.0): (180.0, 240.0),
        (0.5, 0.5, 1.0): (90.0, 130.0),
    }
    renderer = _FakeProjectionRenderer(projected)
    camera = _FakeObservableCamera()
    x_annotations = (
        _FakeProjectedAxisAnnotation((0.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
        _FakeProjectedAxisAnnotation((0.0, 1.0, 0.0), (1.0, 1.0, 0.0)),
    )
    y_annotations = (
        _FakeProjectedAxisAnnotation((0.0, 0.0, 0.0), (0.0, 1.0, 0.0)),
        _FakeProjectedAxisAnnotation((1.0, 0.0, 0.0), (1.0, 1.0, 0.0)),
    )
    z_annotations = (
        _FakeProjectedAxisAnnotation((0.0, 0.0, 0.0), (0.0, 0.0, 2.0)),
        _FakeProjectedAxisAnnotation((0.0, 1.0, 0.0), (0.0, 1.0, 2.0)),
        _FakeProjectedAxisAnnotation((1.0, 0.0, 0.0), (1.0, 0.0, 2.0)),
        _FakeProjectedAxisAnnotation((1.0, 1.0, 0.0), (1.0, 1.0, 2.0)),
    )
    controller = interactive_module._ScreenAnchoredCoordinateAxes(
        renderer=renderer,
        camera=camera,
        bounds=bounds,
        x_annotations=x_annotations,  # type: ignore[arg-type]
        y_annotations=y_annotations,  # type: ignore[arg-type]
        z_annotations=z_annotations,  # type: ignore[arg-type]
    )

    controller.install_camera_observer()
    assert np.allclose(controller.x_opacities, (1.0, 0.0))
    assert np.allclose(controller.y_opacities, (1.0, 0.0))
    assert np.allclose(controller.z_opacities, (1.0, 0.0, 0.0, 1.0))
    assert all(len(annotation.updates) == 1 for annotation in (*x_annotations, *y_annotations, *z_annotations))

    # Emulate a later horizontal angle and invoke the actual camera observer.
    renderer.projected_xy_by_world_point = {
        (0.0, 0.0, 0.0): (0.0, 20.0),
        (0.0, 1.0, 0.0): (80.0, 100.0),
        (1.0, 0.0, 0.0): (-100.0, -80.0),
        (1.0, 1.0, 0.0): (-20.0, 0.0),
        (0.0, 0.0, 2.0): (0.0, 60.0),
        (0.0, 1.0, 2.0): (80.0, 140.0),
        (1.0, 0.0, 2.0): (-100.0, -40.0),
        (1.0, 1.0, 2.0): (-20.0, 40.0),
        (0.5, 0.5, 1.0): (-10.0, 30.0),
    }
    assert callable(camera.callback)
    camera.callback(camera, "ModifiedEvent")  # type: ignore[operator]
    assert np.allclose(controller.x_opacities, (1.0, 0.0))
    assert np.allclose(controller.y_opacities, (0.0, 1.0))
    assert np.allclose(controller.z_opacities, (0.0, 1.0, 1.0, 0.0))
    assert all(len(annotation.updates) == 2 for annotation in (*x_annotations, *y_annotations, *z_annotations))

    controller.dispose()
    assert camera.removed_observer_ids == [17]


def test_screen_edge_crossfade_is_continuous_and_has_no_winner_jump() -> None:
    width = interactive_module.COORDINATE_EDGE_CROSSFADE_PX
    scores = np.asarray((0.0, width * 0.5, width))
    assert np.allclose(
        interactive_module._continuous_edge_opacities(
            scores,
            prefer_minimum=True,
        ),
        (1.0, 0.5, 0.0),
    )
    assert np.allclose(
        interactive_module._continuous_edge_opacities(
            scores,
            prefer_minimum=False,
        ),
        (0.0, 0.5, 1.0),
    )
    assert np.allclose(
        interactive_module._continuous_edge_opacities(
            np.asarray((10.0, 10.0)),
            prefer_minimum=True,
        ),
        (1.0, 1.0),
    )


def test_roi_switch_callbacks_preserve_existing_orientation_axes(monkeypatch) -> None:
    """Replacing local actors must not recreate VTK orientation widgets."""

    plotter = _FakeInteractivePlotter()
    actor = SimpleNamespace(memory_address="sampling-actor")
    roi = SimpleNamespace(is_representative=True, selection_rank=1, cluster_id=0)
    sampling_axis_flags: list[bool] = []
    disposed_axes: list[object] = []
    previous_axes = SimpleNamespace(dispose=lambda: disposed_axes.append(previous_axes))
    replacement_axes = object()
    coordinate_axes_controllers = [None, previous_axes]
    monkeypatch.setattr(
        interactive_module,
        "_add_sampling_boxes",
        lambda *_args, **_kwargs: ([actor], {actor.memory_address: 0}),
    )
    monkeypatch.setattr(
        interactive_module,
        "_add_sampling_active_outline",
        lambda *_args, **_kwargs: object(),
    )

    def record_sampling_scene(
        _plotter: object,
        _roi: object,
        *,
        add_orientation_axes: bool = True,
        diameter_clim_um: tuple[float, float] | None = None,
    ) -> object:
        assert diameter_clim_um is None
        sampling_axis_flags.append(add_orientation_axes)
        return replacement_axes

    monkeypatch.setattr(interactive_module, "_add_sampling_roi_scene", record_sampling_scene)
    callbacks = interactive_module._install_sampling_layer(
        plotter,
        (roi,),
        coordinate_axes_controllers=coordinate_axes_controllers,
    )
    assert callbacks is not None
    assert set(plotter.key_events) == {"a", "A", "r", "R", "s", "S", "c", "C"}
    assert "t" not in plotter.key_events and "T" not in plotter.key_events
    assert len(plotter.picking_callbacks) == 1
    callbacks["select_roi"](actor)
    assert sampling_axis_flags == [False]
    assert plotter.renderer.clear_count == 1
    assert disposed_axes == [previous_axes]
    assert coordinate_axes_controllers[1] is replacement_axes
