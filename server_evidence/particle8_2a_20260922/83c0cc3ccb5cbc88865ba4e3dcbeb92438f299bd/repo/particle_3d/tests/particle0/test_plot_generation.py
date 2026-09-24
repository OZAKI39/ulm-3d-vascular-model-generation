"""Regression for explicit x tick formatter plus scalar y offset formatting."""
import importlib.util
from pathlib import Path


def test_shared_face_plot_renders_custom_distance_ticks(real_field, tmp_path, monkeypatch):
    script = Path(__file__).resolve().parents[2] / "scripts/generate_particle0_report.py"
    spec = importlib.util.spec_from_file_location("particle0_review_plot_regression", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "REPORT", tmp_path)
    monkeypatch.setattr(module, "DATA", tmp_path)
    monkeypatch.setattr(module, "FIGURES", tmp_path)
    metrics = module.stage3(real_field)
    assert metrics["passed"]
    assert (tmp_path / "03_shared_face_continuity.png").is_file()
    assert (tmp_path / "03_shared_face_continuity.csv").is_file()
