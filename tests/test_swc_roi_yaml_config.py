from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

import swc_roi_generate as entrypoint
from utils.swc_roi_yaml_config import load_swc_roi_yaml_config


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "swc_roi_generate.yaml"


def _payload() -> dict[str, object]:
    value = yaml.safe_load(DEFAULT_CONFIG.read_text(encoding="utf-8"))
    assert isinstance(value, dict)
    return value


def test_shipped_yaml_reproduces_current_swc_roi_hyperparameters() -> None:
    settings = load_swc_roi_yaml_config(DEFAULT_CONFIG, project_root=PROJECT_ROOT)

    assert settings.rodent.stage == "all"
    assert settings.rodent.sample_id == "fMOST_0_5_6_0_0_6_0001_02_01"
    assert settings.rodent.spacing_xyz_um == (1.0, 1.0, 2.0)
    assert settings.rodent.expected_shape_zyx == (192, 192, 192)
    assert settings.rodent.smoothing_enabled is False
    assert settings.rodent.figure2a_enabled is True
    assert settings.show_gui is False
    assert settings.sampling_enabled is True
    assert settings.sampling.roi_size_um == (80.0, 80.0, 120.0)
    assert settings.sampling.feature_mode == "radius_plus_structure"
    assert settings.sampling.n_clusters == 5
    assert settings.sampling.selection_mode == "coverage_balanced"


def test_yaml_loader_rejects_unknown_hyperparameter(tmp_path: Path) -> None:
    payload = deepcopy(_payload())
    sampling = payload["sampling"]
    assert isinstance(sampling, dict)
    roi = sampling["roi"]
    assert isinstance(roi, dict)
    roi["misspelled_size"] = [1, 1, 1]
    path = tmp_path / "unknown.yaml"
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")

    with pytest.raises(ValueError, match="Unknown keys in roi"):
        load_swc_roi_yaml_config(path, project_root=PROJECT_ROOT)


def test_yaml_loader_rejects_invalid_vector_length(tmp_path: Path) -> None:
    payload = deepcopy(_payload())
    coordinates = payload["coordinates"]
    assert isinstance(coordinates, dict)
    coordinates["spacing_xyz_um"] = [1.0, 1.0]
    path = tmp_path / "invalid_vector.yaml"
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")

    with pytest.raises(ValueError, match="must contain exactly 3 values"):
        load_swc_roi_yaml_config(path, project_root=PROJECT_ROOT)


def test_yaml_entrypoint_dispatches_configs_and_preserves_source_yaml(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    payload = deepcopy(_payload())
    input_dir = tmp_path / "input"
    input_dir.mkdir()
    paths = payload["paths"]
    assert isinstance(paths, dict)
    paths["input_dir"] = str(input_dir)
    paths["output_dir"] = str(tmp_path / "outputs")
    visualization = payload["visualization"]
    assert isinstance(visualization, dict)
    visualization["enabled"] = False
    source = tmp_path / "run.yaml"
    source.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")

    rodent_root = tmp_path / "outputs" / "rodent_vasculature" / "all_run_test"
    sampling_root = tmp_path / "outputs" / "sampling" / "sampling_test"
    rodent_root.mkdir(parents=True)
    (sampling_root / "config").mkdir(parents=True)
    observed: dict[str, object] = {}

    def fake_rodent(config, *, verbose: bool):
        observed["rodent"] = config
        observed["rodent_verbose"] = verbose
        return SimpleNamespace(
            run_root=rodent_root,
            status="completed",
            html_report=rodent_root / "acceptance_report.html",
            acceptance=SimpleNamespace(overall_status="PASS"),
        )

    def fake_sampling(run_root, config, *, verbose: bool):
        observed["sampling_source"] = run_root
        observed["sampling"] = config
        observed["sampling_verbose"] = verbose
        return SimpleNamespace(
            run_root=sampling_root,
            status="PASS",
            summary_path=sampling_root / "report" / "sampling_summary.json",
        )

    monkeypatch.setattr(entrypoint, "run_rodent_vasculature_pipeline", fake_rodent)
    monkeypatch.setattr(entrypoint, "run_sampling_from_rodent_run", fake_sampling)

    assert entrypoint.main([str(source)]) == 0
    assert observed["sampling_source"] == rodent_root
    assert observed["rodent_verbose"] is True
    assert observed["sampling_verbose"] is True
    assert (rodent_root / "source_swc_roi_generate.yaml").read_bytes() == source.read_bytes()
    assert (
        sampling_root / "config" / "source_swc_roi_generate.yaml"
    ).read_bytes() == source.read_bytes()
