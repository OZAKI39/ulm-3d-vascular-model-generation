from __future__ import annotations

from pathlib import Path
import json
import runpy

import numpy as np
import pytest
import yaml

from utils.rodent_vasculature.catalog import build_catalog, select_records
from utils.rodent_vasculature.swc_io import load_normalized_swc, load_swc, save_normalized_swc
from utils.sampling.pipeline import load_models_from_rodent_run
from utils.swc_roi_yaml_config import load_swc_roi_yaml_config


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def test_brava_catalog_reads_flat_files_and_groups_variants(tmp_path: Path) -> None:
    for name in ("BG001.CNG.swc", "BG001_ColorCoded.CNG.swc", "BH0002.CNG.swc", "BG001.CNG_vmd_smooth.swc"):
        (tmp_path / name).write_text("1 1 0 0 0 1 -1\n")
    (tmp_path / "empty.swc").touch()
    (tmp_path / "BG001.CNG.swc:Zone.Identifier").write_text("not SWC data")
    catalog = build_catalog(tmp_path, "brava")
    assert len(catalog.records) == 5
    assert len(catalog.eligible_records) == 4
    selected = select_records(
        catalog, sample_id="BG001.CNG", parent_group_id=None, split=None, max_samples=1,
    )
    assert [record.swc_path for record in selected] == [tmp_path / "BG001.CNG.swc"]
    assert {record.parent_group_id for record in catalog.eligible_records} == {"BG001", "BH0002"}
    assert all(record.image_path is None and record.mask_path is None for record in catalog.records)


def test_physical_swc_scales_coordinates_and_radius_once(tmp_path: Path) -> None:
    source = tmp_path / "human.swc"
    original = "1 1 1 -2 3 0.62 -1\n2 3 2 -2 3 0.31 1\n"
    source.write_text(original)
    swc = load_swc(source, spacing_xyz_um=(1, 1, 1), volume_shape_zyx=None, swc_units="mm")
    np.testing.assert_allclose(swc.points_um, [[1000, -2000, 3000], [2000, -2000, 3000]])
    np.testing.assert_allclose(swc.radius_raw_um, [620, 310])
    assert swc.parent_ids.tolist() == [-1, 1]
    assert swc.structurally_valid
    saved = save_normalized_swc(swc, tmp_path / "normalized.npz")
    restored = load_normalized_swc(
        saved, source_path=source, spacing_xyz_um=(1, 1, 1), volume_shape_zyx=None,
    )
    np.testing.assert_array_equal(restored.points_um, swc.points_um)
    np.testing.assert_array_equal(restored.radius_raw_um, swc.radius_raw_um)
    assert source.read_text() == original
    mouse = load_swc(source, spacing_xyz_um=(1, 1, 2), volume_shape_zyx=None)
    np.testing.assert_allclose(mouse.points_um[0], [1, -2, 6])
    np.testing.assert_allclose(mouse.radius_raw_um, [0.62, 0.31])


def test_physical_swc_rejects_second_coordinate_scaling(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="Physical SWCs require"):
        load_swc(tmp_path / "unused.swc", spacing_xyz_um=(620, 620, 620), volume_shape_zyx=None, swc_units="mm")


@pytest.mark.parametrize("auxiliary", [None, "image_metadata", "mask_metadata"])
def test_sampling_uses_swc_bounds_without_tiff(tmp_path: Path, auxiliary: str | None) -> None:
    source = tmp_path / "human.swc"
    source.write_text("1 1 1 -2 3 0.62 -1\n2 3 2 -1 4 0.31 1\n")
    swc = load_swc(source, spacing_xyz_um=(1, 1, 1), volume_shape_zyx=None, swc_units="mm")
    sample = tmp_path / "samples" / "brava__human"
    saved = save_normalized_swc(swc, sample / "analysis.npz")
    manifest = {
        "record": {"swc_path": str(source), "sample_id": "brava__human", "parent_group_id": "human"},
        "spacing_xyz_um": [1, 1, 1],
        "image_metadata": None,
        "mask_metadata": None,
        "normalized_swc_path": str(saved),
    }
    if auxiliary:
        manifest[auxiliary] = {"shape_zyx": [5001, 5001, 5001]}
    (sample / "preprocess_manifest.json").write_text(json.dumps(manifest))
    model, = load_models_from_rodent_run(tmp_path)
    expected = (0, 5000, 0, 5000, 0, 5000) if auxiliary else (1000, 2000, -2000, -1000, 3000, 4000)
    assert model.model_bounds_xyz_um == expected
    np.testing.assert_array_equal(model.node_positions_um, swc.points_um)
    np.testing.assert_array_equal(model.node_radius_um, [620, 310])


def test_human_defaults_and_mouse_visualization_are_preserved(tmp_path: Path) -> None:
    path = PROJECT_ROOT / "configs/swc_roi_generate_human.yaml"
    human = load_swc_roi_yaml_config(path, project_root=PROJECT_ROOT)
    mouse = load_swc_roi_yaml_config(PROJECT_ROOT / "configs/swc_roi_generate.yaml", project_root=PROJECT_ROOT)
    assert human.rodent.input_dir == PROJECT_ROOT / "vessel_model/T - Brava/swc_files_vmd"
    assert human.rodent.sample_id == "BG001.CNG_vmd_smooth"
    assert human.rodent.output_root == PROJECT_ROOT / "outputs/human_brava_vmd"
    assert human.rodent.max_samples == 1
    assert human.rodent.swc_units == "mm"
    assert human.rodent.expected_shape_zyx is None
    assert mouse.rodent.swc_units == "voxel"
    assert mouse.rodent.spacing_xyz_um == (1, 1, 2)
    human_yaml = yaml.safe_load(path.read_text())
    mouse_yaml = yaml.safe_load(mouse.source_path.read_text())
    assert human_yaml["visualization"] == mouse_yaml["visualization"]
    assert human_yaml["sampling"]["visualization"] == mouse_yaml["sampling"]["visualization"]
    human_yaml["pipeline"]["sample_id"] = "BG0003.CNG"
    alternate = tmp_path / "alternate.yaml"
    alternate.write_text(yaml.safe_dump(human_yaml))
    assert load_swc_roi_yaml_config(alternate, project_root=PROJECT_ROOT).rodent.sample_id == "BG0003.CNG"
    entrypoint = runpy.run_path(str(PROJECT_ROOT / "s1-2_swc_roi_generate_human.py"))
    assert entrypoint["DEFAULT_CONFIG"] == path
