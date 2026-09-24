import importlib.util
from pathlib import Path


def test_synthetic_plot_functions_rerun_without_real_fem(tmp_path,monkeypatch):
    script=Path(__file__).resolve().parents[2]/'scripts/generate_particle1_report.py'
    spec=importlib.util.spec_from_file_location('p1_plot_regression',script)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    for name in ['REPORT','DATA','FIGURES']: monkeypatch.setattr(module,name,tmp_path)
    module.stage2(); module.stage3(); module.stage4()
    assert len(list(tmp_path.glob('*.png')))==3
    assert len(list(tmp_path.glob('*.md')))==3
