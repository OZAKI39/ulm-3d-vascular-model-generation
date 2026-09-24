"""MeVO data integration uses the unchanged dual-viewport renderer."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest
import yaml

from utils.rodent_vasculature.swc_io import load_swc
from utils.sampling.roi_extraction import global_model_from_swc
from utils.sampling.sampling_io import load_sampling_display_rois
from utils.swc_roi_yaml_config import load_swc_roi_yaml_config
from vascular_processing.mevo_display import display_record, preflight_mevo
from vascular_processing.mevo_graph import load_exact_graph
from vascular_processing.mevo_roi import extract_graph
from vascular_processing.roi_landmarks import LandmarkError, landmark_template, write_landmarks

ROOT = Path(__file__).resolve().parents[1]


def fixture_configuration(tmp_path, *, enabled=True, verified=True):
    data = tmp_path / "data"
    data.mkdir()
    source = data / "synthetic.swc"
    source.write_bytes((ROOT / "tests/data/vascularmd_bifurcation.swc").read_bytes())
    exact = load_exact_graph(source)
    bif, = [n for n in exact.graph if exact.graph.out_degree(n) > 1]
    cuts = [n for n in exact.graph if exact.graph.out_degree(n) == 0]
    landmarks = tmp_path / "landmarks.yaml"
    doc = landmark_template(source, landmarks)
    doc["rois"]["RMCA_M2M3"].update(enabled=enabled, proximal_nodes=[bif] if enabled else [],
        distal_nodes=cuts if enabled else [], proximal_boundary_status="confirmed" if enabled else "unknown",
        manually_verified=verified if enabled else False, notes="SYNTHETIC display test, not a real BraVa anatomical annotation")
    write_landmarks(landmarks, doc)
    cfg = yaml.safe_load((ROOT / "configs/swc_roi_generate_human.yaml").read_text())
    cfg["paths"].update(input_dir=str(data), output_dir=str(tmp_path / "output"))
    cfg["pipeline"]["sample_id"] = "synthetic"
    cfg["visualization"]["figure2a"]["show_gui"] = False
    cfg["mevo"] = {"enabled": True, "landmarks": str(landmarks), "allow_natural_terminal": False}
    configuration = tmp_path / "config.yaml"
    configuration.write_text(yaml.safe_dump(cfg))
    return source, landmarks, doc, configuration


def test_missing_annotations_fail_before_any_spatial_sampling(tmp_path, monkeypatch, capsys):
    _, _, _, path = fixture_configuration(tmp_path, enabled=False)
    spec = importlib.util.spec_from_file_location("mevo_human_entry", ROOT / "s1-2_swc_roi_generate_human.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "run_rodent_vasculature_pipeline", lambda *a, **k: pytest.fail("Must stop before processing"))
    monkeypatch.setattr(module, "run_sampling_from_rodent_run", lambda *a, **k: pytest.fail("No spatial fallback allowed"))
    assert module.main([str(path)]) == 2
    assert "NEEDS_MANUAL_REVIEW" in capsys.readouterr().err
    assert not (tmp_path / "output").exists()


def test_display_nodes_radii_and_edges_are_exact_after_one_mm_conversion(tmp_path):
    source, _, doc, _ = fixture_configuration(tmp_path)
    exact = load_exact_graph(source)
    swc = load_swc(source, spacing_xyz_um=(1, 1, 1), volume_shape_zyx=None, swc_units="mm")
    model = global_model_from_swc(swc, source_model_id="synthetic", source_mouse_id="synthetic")
    definition = doc["rois"]["RMCA_M2M3"]
    extracted = extract_graph(exact, definition)
    roi = display_record(model, exact, definition, extracted, "RMCA_M2M3",
        anatomical_status="UNVERIFIED", units="mm", rank=1)
    assert roi.roi_id.endswith("__RMCA_M2M3__UNVERIFIED")
    for i, n in enumerate(roi.local_node_global_ids):
        np.testing.assert_array_equal(roi.local_node_positions_um[i], exact.graph.nodes[n]["coords"][:3]*1000)
        assert roi.local_node_radius_um[i] == exact.graph.nodes[n]["coords"][3]*1000
    assert set(map(tuple, roi.local_node_global_ids[roi.local_edges])) == set(extracted.graph.edges)
    assert {p.boundary_role for p in roi.cut_ports} == {"ANATOMICAL_PROXIMAL", "ANATOMICAL_DISTAL_AT_SOURCE_TERMINAL"}
    assert len(roi.local_node_ids) == len(extracted.graph)
    assert roi.branch_count == 2 and roi.bifurcation_count == 1


def test_wrong_source_geometry_cannot_be_attached_to_existing_global_view(tmp_path):
    source, _, doc, _ = fixture_configuration(tmp_path)
    exact = load_exact_graph(source)
    swc = load_swc(source, spacing_xyz_um=(1, 1, 1), volume_shape_zyx=None, swc_units="mm")
    model = global_model_from_swc(swc, source_model_id="synthetic", source_mouse_id="synthetic")
    model.node_radius_um *= 2
    definition = doc["rois"]["RMCA_M2M3"]
    with pytest.raises(AssertionError):
        display_record(model, exact, definition, extract_graph(exact, definition), "RMCA_M2M3",
                       anatomical_status="MANUALLY_VERIFIED", units="mm", rank=1)


def test_current_human_preserves_visual_config_and_mouse_default_uses_existing_sampling():
    current = load_swc_roi_yaml_config(ROOT / "configs/swc_roi_generate_human.yaml", project_root=ROOT)
    baseline = yaml.safe_load((ROOT / "configs/swc_roi_generate_human_spatial.yaml").read_text())
    actual = yaml.safe_load((ROOT / "configs/swc_roi_generate_human.yaml").read_text())
    assert {k: v for k, v in actual.items() if k != "mevo"} == baseline
    assert current.mevo_enabled and current.mevo_landmarks.name == "BG001.CNG_vmd_smooth.landmarks.yaml"
    mouse = load_swc_roi_yaml_config(ROOT / "configs/swc_roi_generate.yaml", project_root=ROOT)
    assert not mouse.mevo_enabled


def test_mevo_cannot_silently_fall_back_when_sampling_is_disabled(tmp_path):
    _, _, _, path = fixture_configuration(tmp_path)
    config = yaml.safe_load(path.read_text())
    config["sampling"]["enabled"] = False
    path.write_text(yaml.safe_dump(config))
    with pytest.raises(ValueError, match="mevo.enabled requires sampling.enabled"):
        load_swc_roi_yaml_config(path, project_root=ROOT)


def test_full_human_entry_generates_mevo_roi_library_and_existing_preview(tmp_path):
    _, _, _, path = fixture_configuration(tmp_path)
    process = subprocess.run([sys.executable, str(ROOT / "s1-2_swc_roi_generate_human.py"), str(path)],
                             cwd=ROOT, capture_output=True, text=True, timeout=120)
    assert process.returncode == 0, process.stdout + process.stderr
    run, = (tmp_path / "output/sampling").glob("*_mevo_landmarks")
    summary = json.loads((run / "report/sampling_summary.json").read_text())
    assert summary["status"] == "PASS"
    assert summary["mode"] == "anatomical_landmark_mevo"
    assert not summary["spatial_sampling_performed"] and not summary["clustering_performed"]
    rois = load_sampling_display_rois(run)
    assert len(rois) == 1 and "RMCA_M2M3" in rois[0].roi_id
    assert (run / "figures/interactive_sampling_layer_preview.png").stat().st_size > 1000
