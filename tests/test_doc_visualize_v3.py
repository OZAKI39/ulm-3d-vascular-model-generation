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
            "--debug-cells",
            "--self-test",
        ]
    )
    assert explicit.vtu == Path("field.vtu")
    assert explicit.publication_screenshot == Path("figure.png")


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
