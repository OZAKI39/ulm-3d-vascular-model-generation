"""Anatomical labels in these fixtures are synthetic test annotations only."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

import networkx as nx
import numpy as np
import pyvista as pv
import pytest
import yaml

from vascular_processing.mevo_annotation import annotate
from vascular_processing.mevo_graph import load_exact_graph
from vascular_processing.mevo_pipeline import extract_file, inspect_file
from vascular_processing.mevo_roi import extract_graph
from vascular_processing.roi_landmarks import LandmarkError, landmark_template, load_landmarks, write_landmarks
from vascular_processing.swc_export import read_source

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def network(tmp_path):
    parents = {10: -1, 20: 10, 30: 20, 40: 30, 50: 40, 51: 50, 60: 40, 61: 60,
               70: 30, 80: 70, 81: 80, 90: 70, 91: 90}
    rows = [[n, n % 8, i*2.0, (i % 4)*1.5, i % 3, 0.4+(i % 3)*0.15, p] for i, (n, p) in enumerate(parents.items())]
    source = tmp_path / "synthetic.swc"
    np.savetxt(source, rows, fmt="%.17g")
    landmark = tmp_path / "landmarks.yaml"
    doc = landmark_template(source, landmark)
    roi = doc["rois"]["RMCA_M2M3"]
    roi.update(enabled=True, proximal_nodes=[20], distal_nodes=[50, 60, 80, 90],
               proximal_boundary_status="confirmed", manually_verified=True, notes="SYNTHETIC test only")
    return source, landmark, doc


def run_extract(network, output, **kwargs):
    source, landmarks, document = network
    write_landmarks(landmarks, document)
    return extract_file(source, landmarks, output, **kwargs)


def test_four_distal_cuts_preserve_exact_geometry_types_ids_and_true_connections(network, tmp_path):
    source, _, _ = network
    original = source.read_bytes()
    report = run_extract(network, tmp_path / "result")
    item = report["rois"]["RMCA_M2M3"]
    assert report["status"] == "PASS", report
    assert item["coverage"]["covered_path_count"] == 4
    assert item["coverage"]["unbounded_path_count"] == 0
    assert item["node_count"] == 8 and item["branch_count"] == 7
    assert item["bifurcation_count"] == 3 and item["terminal_count"] == 4
    mapping_path = tmp_path / "result" / item["outputs"]["strict"]["node_mapping"]
    rows = list(csv.DictReader(mapping_path.open()))
    mapping = {int(r["original_swc_id"]): int(r["roi_node_id"]) for r in rows}
    assert set(mapping) == {20, 30, 40, 50, 60, 70, 80, 90}
    saved = read_source(tmp_path / "result" / item["outputs"]["strict"]["swc"])
    raw = read_source(source)
    assert sorted(saved.graph) == list(range(1, 9))
    for old, new in mapping.items():
        np.testing.assert_array_equal(raw.graph.nodes[old]["coords"], saved.graph.nodes[new]["coords"])
        assert raw.graph.nodes[old]["swc_type"] == saved.graph.nodes[new]["swc_type"]
    assert set(saved.graph.edges) == {(mapping[a], mapping[b]) for a, b in raw.graph.edges if a in mapping and b in mapping}
    assert source.read_bytes() == original
    assert item["outputs"]["strict"]["max_xyz_difference_for_mapped_nodes"] == 0
    assert item["outputs"]["strict"]["max_radius_difference_for_mapped_nodes"] == 0


def test_missing_cut_reports_exact_unbounded_path_and_does_not_export(network, tmp_path):
    network[2]["rois"]["RMCA_M2M3"]["distal_nodes"].remove(90)
    report = run_extract(network, tmp_path / "result")
    item = report["rois"]["RMCA_M2M3"]
    assert report["status"] == "NEEDS_MANUAL_REVIEW"
    coverage = item["failure_context"]["coverage"]
    assert coverage["unbounded_path_count"] == 1
    assert coverage["unbounded_paths"][0]["node_ids"] == [20, 30, 70, 90, 91]
    assert coverage["unbounded_paths"][0]["terminal_id"] == 91
    assert coverage["unbounded_paths"][0]["topology_edges"]
    assert not list((tmp_path / "result").glob("*_roi.swc"))


def test_explicit_natural_terminal_opt_in_is_recorded(network, tmp_path):
    network[2]["rois"]["RMCA_M2M3"]["distal_nodes"].remove(90)
    report = run_extract(network, tmp_path / "result", allow_natural_terminal=True)
    item = report["rois"]["RMCA_M2M3"]
    assert report["status"] == "PASS_WITH_WARNINGS"
    assert item["coverage"]["unbounded_path_count"] == 0
    assert item["coverage"]["allowed_natural_terminal_ids"] == [91]
    assert item["warnings"]


def test_exclusion_removes_root_and_all_descendants(network, tmp_path):
    roi = network[2]["rois"]["RMCA_M2M3"]
    roi.update(distal_nodes=[50, 60], exclude_subtree_roots=[70])
    report = run_extract(network, tmp_path / "result")
    item = report["rois"]["RMCA_M2M3"]
    assert item["status"] == "PASS"
    assert item["coverage"]["excluded_path_count"] == 1
    rows = list(csv.DictReader((tmp_path / "result" / item["outputs"]["strict"]["node_mapping"]).open()))
    assert {int(r["original_swc_id"]) for r in rows} == {20, 30, 40, 50, 60}


@pytest.mark.parametrize("change,message", [
    ({"proximal_nodes": [999]}, "do not exist"),
    ({"distal_nodes": [999]}, "do not exist"),
    ({"proximal_nodes": [50], "distal_nodes": [20]}, "not upstream"),
    ({"exclude_subtree_roots": [10]}, "not downstream"),
    ({"distal_nodes": [40, 50, 60, 80, 90]}, "not visited"),
    ({"exclude_subtree_roots": [91]}, "not visited"),
    ({"proximal_nodes": [40, 70]}, "one connected"),
    ({"proximal_nodes": [20, 91]}, "bypass"),
    ({"proximal_boundary_status": "unknown"}, "unknown"),
    ({"distal_nodes": []}, "boundaries are required"),
    ({"internal_landmarks": {"m2_m3_transition": [51]}}, "outside"),
])
def test_invalid_boundaries_fail_explicitly(network, change, message):
    roi = network[2]["rois"]["RMCA_M2M3"]
    roi.update(change)
    with pytest.raises(LandmarkError, match=message):
        extract_graph(load_exact_graph(network[0]), roi)


def test_no_pairwise_path_enumeration_in_extraction(network, monkeypatch):
    exact = load_exact_graph(network[0])
    def forbidden(*args, **kwargs):
        pytest.fail("Main extraction must use linear directed traversal")
    monkeypatch.setattr(nx, "all_simple_paths", forbidden)
    result = extract_graph(exact, network[2]["rois"]["RMCA_M2M3"])
    assert len(result.graph) == 8


def test_multiple_connected_proximals_and_optional_internal_landmarks(network):
    roi = network[2]["rois"]["RMCA_M2M3"]
    roi.update(proximal_nodes=[20, 30], internal_landmarks={"m2_m3_transition": [40, 70]})
    result = extract_graph(load_exact_graph(network[0]), roi)
    assert len(result.graph) == 8 and nx.is_arborescence(result.graph)
    assert result.coverage["unbounded_path_count"] == 0


def test_overlapping_named_rois_keep_all_vtp_memberships(network, tmp_path):
    doc = network[2]
    first = doc["rois"]["RMCA_M2M3"]
    second = doc["rois"]["LACA_A2A3"]
    second.update(enabled=True, proximal_nodes=first["proximal_nodes"].copy(), distal_nodes=first["distal_nodes"].copy(),
                  manually_verified=False, proximal_boundary_status="confirmed", notes="Synthetic overlap QC only")
    report = run_extract(network, tmp_path / "overlap")
    mesh = pv.read(tmp_path / "overlap" / report["outputs"]["qc_vtp"])
    selected = mesh.point_data["RMCA_M2M3_is_roi"].astype(bool)
    assert np.all(mesh.point_data["roi_id"][selected] == -1)
    assert np.all(mesh.point_data["LACA_A2A3_is_roi"][selected] == 1)
    assert np.all(mesh.point_data["roi_membership_bits"][selected] == 6)


def test_truncated_and_unverified_are_never_promoted(network, tmp_path):
    roi = network[2]["rois"]["RMCA_M2M3"]
    roi.update(proximal_boundary_status="truncated_to_available_data", manually_verified=False)
    report = run_extract(network, tmp_path / "result")
    item = report["rois"]["RMCA_M2M3"]
    assert item["anatomical_completeness"] == "partial"
    assert item["anatomical_status"] == "UNVERIFIED" and item["needs_manual_review"]
    assert not roi["manually_verified"]
    assert "ANATOMICAL_STATUS = UNVERIFIED" in (tmp_path / "result" / item["outputs"]["strict"]["swc"]).read_text()


def test_topology_mapping_inspection_and_point_cell_metadata(network, tmp_path):
    report = inspect_file(network[0], tmp_path / "inspect")
    assert report["source"]["node_count"] == 13 and report["topology_node_count"] == 8
    doc = yaml.safe_load((tmp_path / "inspect" / report["outputs"]["landmark_template"]).read_text())
    assert len(doc["rois"]) == 6 and not any(r["enabled"] or r["manually_verified"] for r in doc["rois"].values())
    exact = load_exact_graph(network[0])
    assert exact.locations[10]["topo_node_id"] != 10
    assert exact.locations[20]["topo_node_id"] == -1
    topo = exact.topo
    for a, b, data in topo.edges(data=True):
        for i, original in enumerate(data["full_id"]):
            location = exact.locations[original]
            assert (location["topo_edge_start"], location["topo_edge_end"], location["topo_edge_offset"]) == (a, b, i)
    extracted = run_extract(network, tmp_path / "result")
    mesh = pv.read(tmp_path / "result" / extracted["outputs"]["qc_vtp"])
    for name in ["original_swc_id", "radius", "diameter", "roi_id", "roi_name_code", "is_roi", "is_boundary", "boundary_type", "is_context"]:
        assert name in mesh.point_data and name in mesh.cell_data
    assert mesh.n_points == 13 and mesh.n_lines == 12
    np.testing.assert_array_equal(mesh.point_data["diameter"], mesh.point_data["radius"] * 2)
    assert mesh.point_data["is_roi"].sum() == 8


def test_source_hash_mismatch_requires_override_and_cannot_keep_verified_status(network, tmp_path):
    source, landmarks, doc = network
    doc["source"]["sha256"] = "0" * 64
    write_landmarks(landmarks, doc)
    with pytest.raises(LandmarkError, match="SOURCE MISMATCH"):
        load_landmarks(landmarks, source)
    report = extract_file(source, landmarks, tmp_path / "override", allow_source_mismatch=True)
    item = report["rois"]["RMCA_M2M3"]
    assert item["anatomical_status"] == "UNVERIFIED" and item["needs_manual_review"]
    assert report["warnings"]


@pytest.mark.parametrize("field,value", [("manually_verified", "false"), ("proximal_nodes", [True]),
                                        ("proximal_nodes", [20.0]), ("distal_nodes", [50, 50]),
                                        ("segments", ["M1", "M2"]), ("side", "L")])
def test_schema_rejects_ambiguous_boolean_id_and_anatomical_fields(network, field, value):
    source, landmarks, doc = network
    doc["rois"]["RMCA_M2M3"][field] = value
    write_landmarks(landmarks, doc)
    with pytest.raises(LandmarkError):
        load_landmarks(landmarks, source)


def test_trifurcation_strict_and_context_remain_separate(tmp_path):
    source = tmp_path / "trifurcation.swc"
    np.savetxt(source, [[10, 1, 0, 0, 0, 1, -1], [20, 7, 0, 5, 0, 1, 10],
        [31, 3, -4, 10, 0, .8, 20], [50, 4, 0, 10, 1, .7, 20], [81, 5, 4, 10, 0, .6, 20]], fmt="%.17g")
    path = tmp_path / "annotation.yaml"
    doc = landmark_template(source, path)
    doc["rois"]["RMCA_M2M3"].update(enabled=True, proximal_nodes=[20], distal_nodes=[31, 50, 81],
        manually_verified=True, proximal_boundary_status="confirmed")
    write_landmarks(path, doc)
    report = extract_file(source, path, tmp_path / "result", proximal_context="upstream-edge")
    item = report["rois"]["RMCA_M2M3"]
    assert item["status"] == "PASS", item
    assert item["strict_roi_node_count"] == 4 and item["context_node_count"] == 1
    assert not item["outputs"]["strict"]["vascularmd_reload"]["compatible_with_modeling_adapter"]
    assert item["outputs"]["modelable"]["vascularmd_reload"]["compatible_with_modeling_adapter"]
    context = list(csv.DictReader((tmp_path / "result" / item["outputs"]["modelable"]["node_mapping"]).open()))
    assert [r["original_swc_id"] for r in context if r["is_context"] == "True"] == ["10"]
    assert context[0]["is_anatomical_roi"] == "False"


def test_forest_components_are_supported_without_guessing_territory(network, tmp_path):
    source, landmarks, doc = network
    with source.open("a") as stream:
        stream.write("100 1 99 99 99 1 -1\n110 2 100 99 99 1 100\n")
    from vascular_processing.roi_landmarks import sha256
    doc["source"]["sha256"] = sha256(source)
    report = run_extract(network, tmp_path / "result")
    assert report["status"] == "PASS"
    assert report["source"]["component_count"] == 2
    assert report["rois"]["RMCA_M2M3"]["node_count"] == 8


@pytest.mark.parametrize("text", ["1 1 0 0 0 1 2\n2 3 1 0 0 1 1\n", "1 1 0 0 0 0 -1\n2 3 1 0 0 1 1\n",
                                 "1 1 nan 0 0 1 -1\n2 3 1 0 0 1 1\n", "1 1 0 0 0 1 -1\n2 3 1 0 0 1 99\n"])
def test_invalid_source_never_repaired_or_reoriented(tmp_path, text):
    source = tmp_path / "invalid.swc"
    source.write_text(text)
    with pytest.raises(ValueError):
        load_exact_graph(source)
    assert source.read_text() == text


def test_headless_annotation_saves_template_and_explains_fallback(network, tmp_path, monkeypatch):
    monkeypatch.delenv("DISPLAY", raising=False)
    monkeypatch.delenv("WAYLAND_DISPLAY", raising=False)
    result = annotate(network[0], tmp_path / "headless.yaml")
    assert result["status"] == "NEEDS_MANUAL_REVIEW"
    assert "ParaView" in result["message"]
    assert not any(r["enabled"] for r in yaml.safe_load((tmp_path / "headless.yaml").read_text())["rois"].values())


def test_optional_surface_uses_existing_official_pipeline_and_keeps_strict_roi(tmp_path):
    source = ROOT / "tests/data/vascularmd_bifurcation.swc"
    exact = load_exact_graph(source)
    bif, = [n for n in exact.graph if exact.graph.out_degree(n) > 1]
    tips = [n for n in exact.graph if exact.graph.out_degree(n) == 0]
    path = tmp_path / "model.yaml"
    doc = landmark_template(source, path)
    doc["rois"]["RMCA_M2M3"].update(enabled=True, proximal_nodes=[bif], distal_nodes=tips,
        manually_verified=True, proximal_boundary_status="confirmed", notes="Synthetic vascular modeling regression")
    write_landmarks(path, doc)
    report = extract_file(source, path, tmp_path / "model", proximal_context="upstream-edge", model=True, surface=True)
    item = report["rois"]["RMCA_M2M3"]
    assert item["status"] == "PASS", item
    assert item["model_status"] == "success"
    assert item["outputs"]["strict"]["max_radius_difference_for_mapped_nodes"] == 0
    assert pv.read(tmp_path / "model" / item["model_outputs"]["surface_vtk"]).n_cells > 0
    assert pv.read(tmp_path / "model" / item["model_outputs"]["surface_stl"]).n_cells > 0


def test_model_failure_retains_strict_roi_and_reports_no_fallback(network, tmp_path, monkeypatch):
    monkeypatch.setattr("vascular_processing.mevo_pipeline.process_file", lambda *a, **k: {
        "status": "failed", "outputs": {}, "warnings": [], "failed_components": [{"message": "native test failure"}]})
    report = run_extract(network, tmp_path / "result", model=True)
    assert report["status"] == "FAILED"
    item = report["rois"]["RMCA_M2M3"]
    assert (tmp_path / "result" / item["outputs"]["strict"]["swc"]).is_file()
    assert item["needs_manual_review"]


def test_batch_matches_hash_not_sort_order(network, tmp_path):
    from tools.vascularmd_extract_mevo import batch
    from vascular_processing.roi_landmarks import sha256
    directory, annotations = tmp_path / "input", tmp_path / "annotations"
    directory.mkdir(); annotations.mkdir()
    for filename, annotation_name, offset in [("a.swc", "z.yaml", 0.1), ("b.swc", "a.yaml", 0.2)]:
        rows = np.loadtxt(network[0]); rows[:, 2] += offset
        source = directory / filename; np.savetxt(source, rows, fmt="%.17g")
        doc = yaml.safe_load(yaml.safe_dump(network[2]));doc["source"].update(filename=filename,sha256=sha256(source))
        write_landmarks(annotations / annotation_name, doc)
        blank = annotations / (filename + ".template.yaml")
        write_landmarks(blank, landmark_template(source, blank))
    result = batch(directory, annotations, tmp_path / "batch")
    assert result["status"] == "PASS" and len(result["files"]) == 2
    for row in result["files"]:
        report = json.loads(Path(row["manifest"]).read_text())
        assert report["source_sha256"] == sha256(Path(row["source"]))


def test_real_brava_inspect_smoke(tmp_path):
    source = ROOT / "vessel_model/T - Brava/swc_files/BG001.CNG.swc"
    if not source.exists():
        pytest.skip("BraVa source not installed")
    report = inspect_file(source, tmp_path / "real_inspect")
    assert report["status"] == "PASS_WITH_WARNINGS", report
    assert report["source"]["node_count"] == 2810
    assert report["original_file_preserved"]
    assert report["topology_branch_count"] == 197
