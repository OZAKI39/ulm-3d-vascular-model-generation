from __future__ import annotations

from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace

import pytest
import yaml

import swc_stl_model_generate as entrypoint
from utils.cfd_lumen.model_yaml_config import load_swc_stl_yaml_config
from utils.cfd_lumen.ultraliser_backend import _preserve_source_configuration


PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONFIG = PROJECT_ROOT / "configs" / "swc_stl_model_generate.yaml"


def _payload() -> dict[str, object]:
    value = yaml.safe_load(DEFAULT_CONFIG.read_text(encoding="utf-8"))
    assert isinstance(value, dict)
    return value


def test_shipped_yaml_reproduces_validated_ultraliser_settings() -> None:
    settings = load_swc_stl_yaml_config(DEFAULT_CONFIG, project_root=PROJECT_ROOT)

    assert settings.roi_anchor == 3274
    assert settings.roi_id is None
    assert settings.surface_backend == "ultraliser"
    assert settings.sampling_run == (
        PROJECT_ROOT / "outputs" / "sampling" / "20260825_133201_radius_plus_structure_k5"
    ).resolve()
    assert settings.lumen.ultraliser.radius_scale == 0.91
    assert settings.lumen.ultraliser.voxels_per_micron == 6.0
    assert settings.lumen.ultraliser.threads == 8
    assert settings.lumen.ultraliser.laplacian_iterations == 10
    assert settings.lumen.surface_qc.max_radius_p95_error == 0.05
    assert settings.lumen.surface_qc.require_watertight is True


def test_yaml_loader_rejects_unknown_reconstruction_key(tmp_path: Path) -> None:
    payload = deepcopy(_payload())
    reconstruction = payload["reconstruction"]
    assert isinstance(reconstruction, dict)
    reconstruction["fallback_backend"] = "custom"
    path = tmp_path / "unknown.yaml"
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")

    with pytest.raises(ValueError, match="Unknown keys in reconstruction"):
        load_swc_stl_yaml_config(path, project_root=PROJECT_ROOT)


def test_yaml_loader_requires_exactly_one_roi_selector(tmp_path: Path) -> None:
    payload = deepcopy(_payload())
    selection = payload["selection"]
    assert isinstance(selection, dict)
    selection["roi_id"] = "another-roi"
    path = tmp_path / "ambiguous.yaml"
    path.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")

    with pytest.raises(ValueError, match="Exactly one"):
        load_swc_stl_yaml_config(path, project_root=PROJECT_ROOT)


def test_yaml_entrypoint_selects_roi_and_passes_source_config(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    roi = SimpleNamespace(
        roi_id="roi-anchor-3274",
        anchor_id=3274,
    )
    observed: dict[str, object] = {}

    monkeypatch.setattr(
        entrypoint,
        "resolve_sampling_run",
        lambda path, *, project_root: Path(path),
    )
    monkeypatch.setattr(
        entrypoint,
        "load_sampling_rois",
        lambda sampling_run, **kwargs: [roi],
    )

    def fake_reconstruction(selected_roi, config, **kwargs):
        observed["roi"] = selected_roi
        observed["config"] = config
        observed.update(kwargs)
        return {
            "run_root": str(PROJECT_ROOT / "outputs" / "model_generate" / "test"),
            "source_configuration": str(kwargs["source_config_path"]),
            "status": "PASS",
        }

    monkeypatch.setattr(entrypoint, "run_ultraliser_reconstruction", fake_reconstruction)

    assert entrypoint.main([str(DEFAULT_CONFIG)]) == 0
    assert observed["roi"] is roi
    assert observed["source_config_path"] == DEFAULT_CONFIG.resolve()
    assert observed["output_root"] == (PROJECT_ROOT / "outputs" / "model_generate").resolve()
    assert observed["ultraliser_root"] == (PROJECT_ROOT / "Ultraliser").resolve()


def test_source_yaml_is_preserved_byte_for_byte_in_model_input(tmp_path: Path) -> None:
    source = tmp_path / "experiment.yaml"
    source.write_bytes(DEFAULT_CONFIG.read_bytes())
    input_directory = tmp_path / "run" / "input"
    input_directory.mkdir(parents=True)

    copied = _preserve_source_configuration(input_directory, source)

    assert copied == input_directory / "source_swc_stl_model_generate.yaml"
    assert copied.read_bytes() == source.read_bytes()
