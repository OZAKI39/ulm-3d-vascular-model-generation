"""Exercise native mouse fitting, voxel round trips and untouched forest data."""
from pathlib import Path
import json

import numpy as np
import pytest
import yaml

from vascular_processing.mouse_dataset import prepare_mouse_dataset
from vascular_processing.pipeline import Options
from utils.rodent_vasculature.swc_io import load_swc

ROOT = Path(__file__).resolve().parents[1]


def mouse_input(tmp_path):
    dataset = tmp_path / "original"
    folder = dataset / "raw_data/analysis_data/analysis_data"
    source = folder / "swc/test_mouse_0001_02_01.swc"
    source.parent.mkdir(parents=True)
    # A reference-only line and a longer analysis component; anisotropic Z.
    rows = np.array([[1, 2, 100, 100, 40, 0.9, -1], [2, 2, 101, 100, 40, 0.8, 1],
                     [10, 1, 0, 5, 1, 1.2, -1], [11, 3, 5, 5, 3, 1.25, 10],
                     [12, 3, 10, 5, 5, 0.65, 11], [13, 3, 15, 5, 7, 1.22, 12],
                     [14, 3, 20, 5, 9, 1.18, 13]], dtype=float)
    np.savetxt(source, rows, fmt="%.17g")
    # Preparation links auxiliary files without decoding or rewriting them.
    for name in ("images", "mask"):
        (folder / name).mkdir()
        (folder / name / (source.stem + ".tif")).write_bytes(b"unchanged auxiliary bytes")
    config = yaml.safe_load((ROOT / "configs/swc_roi_generate_mouse_raw.yaml").read_text())
    config["paths"]["input_dir"] = str(dataset)
    config["pipeline"]["sample_id"] = source.stem
    path = tmp_path / "raw.yaml"
    path.write_text(yaml.safe_dump(config))
    return path, source, rows


@pytest.mark.parametrize("model_mode", ["network", "branches"])
def test_native_mouse_model_preserves_scale_reference_forest_and_auxiliary_files(tmp_path, model_mode):
    config, source, original = mouse_input(tmp_path)
    original_bytes = source.read_bytes()
    report = prepare_mouse_dataset(config, tmp_path / "derived", project_root=ROOT, model_mode=model_mode)
    assert report["status"] == "success", report
    sample, = report["samples"]
    assert sample["original_selection"]["selected_component_id"] == 1
    assert sample["reference_component_count_preserved"] == 2
    assert sample["reference_only_nodes_preserved"] == 2
    output = load_swc(Path(sample["output_swc"]), spacing_xyz_um=(1, 1, 2), volume_shape_zyx=None)
    np.testing.assert_array_equal(output.points_voxel_xyz[:2], original[:2, 2:5])
    np.testing.assert_array_equal(output.radius_raw_um[:2], original[:2, 5])
    np.testing.assert_array_equal(output.parent_ids[:2], [-1, 1])
    assert set(output.node_ids[2:]).isdisjoint(original[:, 0])
    native = json.loads(Path(sample["native_qc"]).read_text())
    physical_input = np.loadtxt(native["source"])
    np.testing.assert_array_equal(physical_input[:, 2:5], original[2:, 2:5] * [1, 1, 2])
    np.testing.assert_array_equal(physical_input[:, 5], original[2:, 5])
    physical_output = np.loadtxt(native["outputs"]["swc"])
    np.testing.assert_allclose(output.points_um[2:], physical_output[:, 2:5], rtol=0, atol=1e-12)
    np.testing.assert_array_equal(output.radius_raw_um[2:], physical_output[:, 5])
    assert not np.array_equal(output.radius_raw_um[2:], original[2:, 5])
    assert native["options"]["input_units"] == "um"
    assert native["preprocessing"]["density_unit"] == "points/um"
    assert source.read_bytes() == original_bytes
    assert Path(sample["original_backup"]["path"]).read_bytes() == original_bytes
    assert not Path(sample["original_backup"]["path"]).is_symlink()
    assert Path(report["config_backup"]["path"]).read_bytes() == config.read_bytes()
    for item in sample["auxiliary"].values():
        link = Path(item["link"])
        assert link.is_symlink() and link.resolve() == Path(item["source"])
        assert link.read_bytes() == b"unchanged auxiliary bytes"


def test_native_failure_never_publishes_a_replacement_mouse_swc(tmp_path, monkeypatch):
    config, source, _ = mouse_input(tmp_path)
    monkeypatch.setattr("vascular_processing.mouse_dataset.process_file", lambda *a, **k: {
        "status": "failed", "failed_components": [{"stage": "native_model", "message": "test failure"}]})
    output = tmp_path / "derived"
    report = prepare_mouse_dataset(config, output, project_root=ROOT)
    assert report["status"] == "failed"
    assert report["samples"][0]["failed_components"][0]["stage"] == "native_model"
    assert not (output / "raw_data/analysis_data/analysis_data/swc" / source.name).exists()
    assert Path(report["samples"][0]["original_backup"]["path"]).read_bytes() == source.read_bytes()


def test_preparation_cannot_write_inside_original_mouse_dataset(tmp_path):
    config, source, _ = mouse_input(tmp_path)
    with pytest.raises(ValueError, match="outside the original"):
        prepare_mouse_dataset(config, source.parent / "derived", project_root=ROOT)


def test_millimetre_spacing_cannot_silently_control_micrometre_sampling():
    with pytest.raises(ValueError, match="--spacing-mm requires"):
        Options(input_units="um", sample_mode="spacing", spacing_mm=0.5).validate()
