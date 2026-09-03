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


class _FakeInteractivePlotter:
    def __init__(self) -> None:
        self.renderer = _FakeRenderer()
        self.render_count = 0
        self.key_events: dict[str, object] = {}
        self.picking_callbacks: list[object] = []
        self.radio_buttons: list[_FakeButtonWidget] = []
        self.window_size = (1800, 900)
        self.text_entries: dict[str, tuple[str, object]] = {}

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


class _FakeTextProperty:
    def __init__(self) -> None:
        self.family = ""
        self.font_size = 0
        self.bold = False

    def SetFontFamilyToArial(self) -> None:
        self.family = "Arial"

    def SetFontSize(self, font_size: int) -> None:
        self.font_size = int(font_size)

    def SetBold(self, bold: bool) -> None:
        self.bold = bool(bold)


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


def test_roi_switch_callbacks_preserve_existing_orientation_axes(monkeypatch) -> None:
    """Replacing local actors must not recreate VTK orientation widgets."""

    plotter = _FakeInteractivePlotter()
    actor = SimpleNamespace(memory_address="sampling-actor")
    roi = SimpleNamespace(is_representative=True, selection_rank=1, cluster_id=0)
    sampling_axis_flags: list[bool] = []
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
    ) -> None:
        sampling_axis_flags.append(add_orientation_axes)

    monkeypatch.setattr(interactive_module, "_add_sampling_roi_scene", record_sampling_scene)
    callbacks = interactive_module._install_sampling_layer(plotter, (roi,))
    assert callbacks is not None
    assert set(plotter.key_events) == {"a", "A", "r", "R", "s", "S", "c", "C"}
    assert "t" not in plotter.key_events and "T" not in plotter.key_events
    assert len(plotter.picking_callbacks) == 1
    callbacks["select_roi"](actor)
    assert sampling_axis_flags == [False]
    assert plotter.renderer.clear_count == 1
