from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pyvista as pv
import pytest

import doc_visualize_v3 as visualizer


def _write_json(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def _make_candidate(root: Path, name: str, completed: str) -> Path:
    run = root / "outputs" / "cfd_flow" / name
    vtu = run / "flow" / "production_steady_flow_field.vtu"
    vtu.parent.mkdir(parents=True)
    vtu.write_bytes(name.encode("ascii"))
    _write_json(
        vtu.parent / "production_steady_flow_field_manifest.json",
        {
            "status": "PASS",
            "sha256": visualizer.sha256_file(vtu),
            "cell_count": 1,
        },
    )
    _write_json(
        run / "qc" / "run_summary.json",
        {"status": visualizer.EXPECTED_RUN_STATUS, "completed_at": completed},
    )
    _write_json(run / "qc" / "production_steady_qc.json", {"status": "PASS"})
    _write_json(run / "qc" / "production_primary_metrics.json", {})
    _write_json(run / "steady_replay" / "physical_port_flux.json", {})
    return run


def _small_grid() -> pv.UnstructuredGrid:
    grid = pv.ImageData(dimensions=(3, 2, 2)).cast_to_unstructured_grid()
    velocity = np.array(((1.0, 0.0, 0.0), (0.0, 2.0, 0.0)))
    speed = np.linalg.norm(velocity, axis=1)
    grid.cell_data["velocity_phy"] = velocity
    grid.cell_data["velocity_magnitude_m_s"] = speed
    grid.cell_data["velocity_magnitude_mm_s"] = speed * 1.0e3
    grid.cell_data["pressure_gauge_pa"] = np.array((-1.0, 2.0))
    grid.cell_data["pressure_absolute_solver_pa"] = np.array((9.0, 12.0))
    grid.cell_data["rho_lattice"] = np.array((0.999, 1.001))
    return grid


def _plane_contract(contract_hash: str) -> dict:
    plane = {
        "origin_m": [1.0e-6, 2.0e-6, 3.0e-6],
        "unit_normal": [1.0, 0.0, 0.0],
        "basis_u": [0.0, 1.0, 0.0],
        "basis_v": [0.0, 0.0, 1.0],
        "physical_aperture_area_m2": 1.0e-12,
        "local_hydraulic_diameter_m": 1.0e-6,
        "physical_aperture_contour_uv_m": [
            [1.0e-6 * np.cos(angle), 1.0e-6 * np.sin(angle)]
            for angle in np.linspace(0.0, 2.0 * np.pi, 12, endpoint=False)
        ],
    }
    return {
        "status": "PASS",
        "contract_sha256": contract_hash,
        "ports": {
            label: {"planes": {"central": plane}} for label in visualizer.PORT_ORDER
        },
    }


def test_cli_defaults_and_required_options_are_available():
    parser = visualizer.build_parser()
    args = parser.parse_args([])
    assert args.scalar == "velocity"
    assert (args.window_width, args.window_height) == (1800, 1100)
    assert args.streamline_seeds == 24
    explicit = parser.parse_args(
        [
            "--run-dir",
            "run",
            "--vtu",
            "field.vtu",
            "--scalar",
            "pressure",
            "--no-streamlines",
            "--no-ports",
            "--full-range",
            "--publication-screenshot",
            "figure.png",
            "--publication-suite",
            "figures",
            "--projection",
            "perspective",
            "--ui-mode",
            "analysis",
            "--debug-cells",
            "--self-test",
        ]
    )
    assert explicit.vtu == Path("field.vtu")
    assert explicit.publication_screenshot == Path("figure.png")
    assert explicit.publication_suite == Path("figures")
    assert explicit.projection == "perspective"
    assert explicit.ui_mode == "analysis"


def test_explicit_vtu_has_priority_over_explicit_run(tmp_path: Path):
    vtu = tmp_path / "chosen" / "flow" / "chosen.vtu"
    run, selected = visualizer.resolve_run_and_vtu(
        project_root=tmp_path,
        explicit_run_dir=tmp_path / "ignored",
        explicit_vtu=vtu,
    )
    assert selected == vtu.resolve()
    assert run == vtu.resolve().parent.parent


def test_latest_complete_pass_run_is_discovered(tmp_path: Path):
    older = _make_candidate(
        tmp_path,
        "production_tau1_base_promotion_anchor003274_20260101_000000",
        "2026-01-01T00:00:00",
    )
    newer = _make_candidate(
        tmp_path,
        "production_tau1_base_promotion_anchor003274_20260102_000000",
        "2026-01-02T00:00:00",
    )
    assert visualizer.locate_latest_run(tmp_path) == newer.resolve()
    (newer / "qc" / "production_steady_qc.json").unlink()
    assert visualizer.locate_latest_run(tmp_path) == older.resolve()


def test_field_range_reports_raw_and_robust_limits():
    values = np.r_[np.linspace(1.0, 2.0, 1000), 1000.0]
    result = visualizer.calculate_field_range(values)
    assert result.raw_min == 1.0
    assert result.raw_max == 1000.0
    assert result.percentile_min > 1.0
    assert result.percentile_max < 1000.0
    assert result.selected(False) != result.selected(True)


def test_field_range_rejects_nonfinite_values():
    with pytest.raises(visualizer.VisualizerInputError, match="non-finite"):
        visualizer.calculate_field_range(np.array((1.0, np.nan)))


def test_required_vtu_field_contract_is_strict():
    grid = _small_grid()
    visualizer.validate_required_fields(grid, {"status": "PASS", "cell_count": 2})
    del grid.cell_data["pressure_gauge_pa"]
    with pytest.raises(visualizer.VisualizerInputError, match="missing required"):
        visualizer.validate_required_fields(
            grid, {"status": "PASS", "cell_count": 2}
        )


def test_plane_contract_is_loaded_from_production_snapshot(tmp_path: Path):
    contract_hash = "a" * 64
    contract_path = tmp_path / "contracts" / "planes.json"
    _write_json(contract_path, _plane_contract(contract_hash))
    run = tmp_path / "outputs" / "cfd_flow" / "production_run"
    snapshot = run / "input" / "cfd_flow_promotion_regression.yaml"
    snapshot.parent.mkdir(parents=True)
    snapshot.write_text(
        "paths:\n  physical_plane_contract: contracts/planes.json\n",
        encoding="utf-8",
    )
    loaded = visualizer.load_plane_contract(
        run,
        {"plane_contract_sha256": contract_hash},
        project_root=tmp_path,
    )
    assert loaded["contract_sha256"] == contract_hash
    assert len(loaded["ports"]) == 4


def test_visualizer_has_no_solver_launch_path():
    source = (visualizer.PROJECT_ROOT / "doc_visualize_v3.py").read_text(
        encoding="utf-8"
    )
    for forbidden in ("subprocess", "run_wsl_tool", "mpirun", "generate_musubi_lua"):
        assert forbidden not in source


def _viewer_plane_contract() -> dict:
    angles = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    origins_um = {
        "inlet": (0.0, 0.5, 0.5),
        "outlet_01": (2.0, 0.2, 0.2),
        "outlet_02": (2.0, 0.5, 0.5),
        "outlet_03": (2.0, 0.8, 0.8),
    }
    ports = {}
    for label, origin_um in origins_um.items():
        ports[label] = {
            "planes": {
                "central": {
                    "origin_m": (np.asarray(origin_um) * 1.0e-6).tolist(),
                    "unit_normal": [1.0, 0.0, 0.0],
                    "basis_u": [0.0, 1.0, 0.0],
                    "basis_v": [0.0, 0.0, 1.0],
                    "local_hydraulic_diameter_m": 0.3e-6,
                    "physical_aperture_contour_uv_m": [
                        [0.12e-6 * np.cos(angle), 0.12e-6 * np.sin(angle)]
                        for angle in angles
                    ],
                }
            }
        }
    return {"status": "PASS", "ports": ports}


def _small_visual_data(tmp_path: Path) -> visualizer.VisualData:
    grid = _small_grid()
    centers = np.asarray(grid.cell_centers().points)
    original_hashes = {
        name: visualizer._array_sha256(np.asarray(grid.cell_data[name]))
        for name in visualizer.REQUIRED_ARRAYS
    }
    display = visualizer._build_display_grid(grid)
    surface = visualizer._build_overview_surface(display)
    fields = {
        "velocity": visualizer.FieldSpec(
            "velocity_magnitude_mm_s", "Velocity magnitude", "mm s⁻¹", "viridis"
        ),
        "pressure": visualizer.FieldSpec(
            "pressure_gauge_pa", "Gauge pressure", "Pa", "cividis"
        ),
        "rho": visualizer.FieldSpec("rho_lattice", "Lattice density", "–", "cividis"),
    }
    ranges = {
        key: visualizer.calculate_field_range(grid.cell_data[field.array])
        for key, field in fields.items()
    }
    line = pv.lines_from_points(np.array(((0.0, 0.5, 0.5), (2.0, 0.5, 0.5))))
    line.point_data["velocity_magnitude_mm_s"] = np.array((1.0, 2.0))
    vtu = tmp_path / "flow" / "production_steady_flow_field.vtu"
    return visualizer.VisualData(
        run_dir=tmp_path,
        vtu_path=vtu,
        vtu_sha256="f" * 64,
        manifest={"status": "PASS", "cell_count": 2},
        metrics={
            "iteration": 10,
            "rho_mean": 1.0,
            "Qin_m3_s": 1.0e-15,
            "physical_volume_closure": 1.0e-6,
        },
        steady_qc={"status": "PASS"},
        run_summary={
            "steady_solution_source": "VALIDATED_RESEARCH_BASE_ACCEPTED_RESTART",
            "numerical_contract": {
                "dx_m": 0.2e-6,
                "tau": 1.0,
                "target_volume_flow_m3_s": 1.0e-15,
            },
        },
        physical_flux={"status": "PASS"},
        plane_contract=_viewer_plane_contract(),
        grid_um=grid,
        display_grid_um=display,
        surface_um=surface,
        centers_um=centers,
        center_tree_um=visualizer.cKDTree(centers),
        fields=fields,
        ranges=ranges,
        streamlines_um=line,
        valid_streamline_count=1,
        original_cell_array_sha256=original_hashes,
    )


def test_display_interpolation_is_point_only_and_original_cells_are_unchanged():
    grid = _small_grid()
    original = np.asarray(grid.cell_data["pressure_gauge_pa"]).copy()
    display = visualizer._build_display_grid(grid)
    assert "pressure_gauge_pa" in display.point_data
    assert "pressure_gauge_pa" in display.cell_data
    assert np.array_equal(grid.cell_data["pressure_gauge_pa"], original)
    assert visualizer.RENDER_INTERPOLATION == "CELL_TO_POINT_DISPLAY_ONLY"


def test_academic_style_meets_scalarbar_typography_and_port_contracts():
    style = visualizer.AcademicStyle()
    field = visualizer.FieldSpec("v", "Velocity magnitude", "mm s⁻¹", "viridis")
    scalar_bar = visualizer.AcademicLayout().scalar_bar_args(field, style)
    assert style.background == "#FBFBFA"
    assert style.title_font_size <= 14
    assert 0.020 <= scalar_bar["width"] <= 0.032
    assert 0.32 <= scalar_bar["height"] <= 0.42
    assert scalar_bar["n_labels"] <= 5
    assert scalar_bar["title"] == "Velocity magnitude\nmm s⁻¹"
    assert visualizer.PORT_COLORS == {
        "inlet": "#0072B2",
        "outlet_01": "#009E73",
        "outlet_02": "#E69F00",
        "outlet_03": "#CC79A7",
    }


def test_camera_is_deterministic_and_fitted_to_vessel_only():
    points = np.asarray(
        _small_grid().extract_surface(algorithm="dataset_surface").points
    )
    first = visualizer.academic_camera_parameters(points, 1.5)
    second = visualizer.academic_camera_parameters(points, 1.5)
    assert np.allclose(first["position_um"], second["position_um"])
    assert np.allclose(first["focal_point_um"], np.mean(points, axis=0))
    assert first["projected_height_fraction"] >= 0.55
    assert first["projected_width_fraction"] >= 0.35


def test_clean_overview_defaults_widget_lifecycle_and_original_cell_pick(tmp_path: Path):
    data = _small_visual_data(tmp_path)
    config = visualizer.VisualConfig(
        width=1200, height=800, build_streamlines=True, show_ports=True
    )
    viewer = visualizer.AcademicCFDViewer(data, config, off_screen=True)
    try:
        assert viewer.visual_mode == "overview"
        assert viewer.projection == "parallel"
        assert viewer.visible_helper_actors() == {
            "plane_widget": False,
            "help": False,
            "info": False,
            "picked_marker": False,
            "vectors": False,
            "streamlines": False,
            "bounding_box": False,
            "port_normals": False,
            "ports": True,
        }
        assert "vessel_context" not in viewer.plotter.actors
        assert all(f"port_outline_{label}" in viewer.plotter.actors for label in visualizer.PORT_ORDER)
        assert not any(name.startswith("port_normal_") for name in viewer.plotter.actors)
        viewer.show_plane_widget("clip")
        assert viewer.visual_mode == "clip"
        assert viewer.plane_widget_visible
        viewer.hide_plane_widget()
        assert viewer.visual_mode == "clip"
        assert not viewer.plane_widget_visible
        viewer.set_visual_mode("slice")
        assert viewer.context_actor is not None
        viewer.set_visual_mode("overview")
        viewer._pick_callback(data.centers_um[1])
        assert viewer.picked_cell_id == 1
        assert viewer.picked_marker_visible
        assert float(data.grid_um.cell_data["pressure_gauge_pa"][1]) == 2.0
    finally:
        viewer.plotter.close()


def test_publication_suite_contract_covers_five_clean_scenes():
    assert visualizer.PUBLICATION_SCENES == (
        ("01_after_velocity_overview.png", "velocity", "overview", False),
        ("02_after_pressure_overview.png", "pressure", "overview", False),
        ("03_after_velocity_clip.png", "velocity", "clip", False),
        ("04_after_pressure_slice.png", "pressure", "slice", False),
        ("05_after_streamlines.png", "velocity", "clip", True),
    )
