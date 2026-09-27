"""Integration tests exercise the real pinned upstream algorithm, never a filter substitute."""
from __future__ import annotations

import csv
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyvista as pv
import pytest

from third_party.vascularmd.Spline import Spline
from vascular_processing.pipeline import Options, process_file
from vascular_processing.swc_export import read_source
from vascular_processing.vascularmd_adapter import load_tree

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests/data"


def run_native(path, output, **kwargs):
    report = process_file(path, output, Options(**kwargs))
    assert report["status"] == "success", report["failed_components"]
    assert report["original_file_preserved"]
    assert report["round_trip_valid"]
    assert report["topology_preserved"]
    return report


def test_official_example_exports_valid_swc_surface_and_diagnostics(tmp_path):
    path = ROOT / "third_party/vascularmd/Data/example_centerline_ICA.swc"
    report = run_native(path, tmp_path, qc_plot=True)
    raw, smooth = (report["comparison"][name] for name in ("raw", "smooth"))
    for key in ("node_count", "branch_count", "terminal_count", "bifurcation_count"):
        assert raw[key] == smooth[key]
    assert smooth["log_radius_jump"]["p95"] < raw["log_radius_jump"]["p95"]
    mesh = pv.read(report["outputs"]["surface_vtk"])
    assert mesh.n_cells > 100
    assert np.isfinite(mesh.points).all()
    assert Path(report["outputs"]["qc_plot"]).stat().st_size > 1000


def test_five_point_spike_uses_native_separate_aic_models(tmp_path, monkeypatch):
    calls = []
    original = Spline.approximation

    def traced(self, *args, **kwargs):
        result = original(self, *args, **kwargs)
        if kwargs.get("criterion") == "AIC" and kwargs.get("radius_model", True):
            calls.append(self.get_lbd())
        return result

    monkeypatch.setattr(Spline, "approximation", traced)
    report = run_native(DATA / "vascularmd_spike.swc", tmp_path, surface=False)
    assert calls and all(len(lambdas) == 2 and all(value > 0 for value in lambdas) for lambdas in calls)
    assert report["preprocessing"]["original_node_count"] == 5
    assert report["preprocessing"]["vascularmd_preprocessing_node_count"] != 5
    assert report["comparison"]["smooth"]["node_count"] == 5
    raw, smooth = (report["comparison"][name]["log_radius_jump"] for name in ("raw", "smooth"))
    assert smooth["maximum"] < raw["maximum"] * 0.5
    # A genuine change to the data is expected; this is regularization, not exact interpolation.
    output = read_source(Path(report["outputs"]["swc"]))
    assert not np.array_equal(output.rows[:, 5], [1.2, 1.25, 0.65, 1.22, 1.18])


def test_bifurcation_has_one_shared_node_and_distinct_daughter_calibers(tmp_path):
    report = run_native(DATA / "vascularmd_bifurcation.swc", tmp_path)
    graph = read_source(Path(report["outputs"]["swc"])).graph
    forks = [node for node in graph if graph.out_degree(node) > 1]
    assert len(forks) == 1
    assert graph.out_degree(forks[0]) == 2 and graph.in_degree(forks[0]) == 1
    assert graph.nodes[forks[0]]["swc_type"] == 7
    radii = {graph.nodes[n]["swc_type"]: graph.nodes[n]["coords"][3] for n in graph if graph.out_degree(n) == 0}
    assert set(radii) == {5, 6}
    assert abs(radii[5] - radii[6]) > 0.2
    assert report["comparison"]["smooth"]["terminal_count"] == 2
    assert all(b["raw_count"] == b["smooth_count"] for b in report["comparison"]["branches"])


def test_spacing_is_uniform_native_arc_length_and_no_auto_resample(tmp_path):
    report = run_native(DATA / "vascularmd_spike.swc", tmp_path, sample_mode="spacing", spacing_mm=0.7, auto_resample=False, surface=False)
    assert report["preprocessing"]["vascularmd_preprocessing_node_count"] == 5
    with Path(report["outputs"]["radius_csv"]).open() as handle:
        rows = [row for row in csv.DictReader(handle) if row["series"] == "smooth"]
    steps = np.diff([float(row["arc_length"]) for row in rows])
    assert 0 < steps.min() <= steps.max() <= 0.7 + 1e-10
    np.testing.assert_allclose(steps, steps[0], rtol=1e-12)


def test_native_loader_does_not_convert_radius_diameter_or_units():
    source = read_source(DATA / "vascularmd_bifurcation.swc")
    tree, _, _ = load_tree(source, False)
    for node in source.graph:
        np.testing.assert_array_equal(tree.get_full_graph().nodes[node]["coords"], source.graph.nodes[node]["coords"])


def test_model_failure_reports_component_without_silent_fallback(tmp_path, monkeypatch):
    from vascular_processing.vascularmd_adapter import TracedTree

    def fail(self, **kwargs):
        self.active_component = {"stage": "model_bifurcation", "source_swc_node": 21}
        raise RuntimeError("Injected native modeling failure")

    monkeypatch.setattr(TracedTree, "model_network", fail)
    path = DATA / "vascularmd_bifurcation.swc"
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    report = process_file(path, tmp_path, Options())
    assert report["status"] == "failed" and report["needs_manual_review"]
    assert report["failed_components"][0]["component"]["source_swc_node"] == 21
    assert not list(tmp_path.glob("*_smooth.swc"))
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest


@pytest.mark.parametrize("accept_native_merges", [False, True])
def test_native_topology_changes_are_rejected(tmp_path, monkeypatch, accept_native_merges):
    from vascular_processing.vascularmd_adapter import TracedTree
    original = TracedTree.model_network

    def remove_terminal(self, **kwargs):
        original(self, **kwargs)
        graph = self.get_model_graph()
        graph.remove_node(next(node for node in graph if graph.out_degree(node) == 0))

    monkeypatch.setattr(TracedTree, "model_network", remove_terminal)
    report = process_file(DATA / "vascularmd_spike.swc", tmp_path,
                          Options(surface=False, accept_native_merges=accept_native_merges))
    assert report["status"] == "failed" and report["needs_manual_review"]
    assert not list(tmp_path.glob("*_smooth.swc"))


def test_exhausted_native_apex_search_stops_instead_of_refitting_forever(tmp_path, monkeypatch):
    # Simulate an upstream intersection search that reaches the endpoint forever.
    monkeypatch.setattr(Spline, "first_intersectionv2", lambda *args: (np.zeros(3), [1.0, 1.0]))
    source = DATA / "vascularmd_bifurcation.swc"
    report = process_file(source, tmp_path / "output", Options(surface=False, accept_native_merges=True))
    assert report["status"] == "failed"
    failure, = report["failed_components"]
    assert failure["error_type"] == "NativeApexNotFound"
    assert failure["component"]["stage"] == "native_apex_search"
    assert report["original_file_preserved"]
    assert "swc" not in report["outputs"]


def test_scoped_native_point_cache_is_numerically_identical_and_restores_getter():
    from vascular_processing.vascularmd_adapter import cache_native_distance_points
    spline = Spline(np.array([[0, 0, 0, 1], [2, 1, 1, 1.2], [5, 3, 0, 1.1], [8, 4, 2, 0.8]]), order=4)
    data = np.array([[1, 1, 0, 1], [3, 1, 2, 1], [7, 2, 1, 1]], dtype=float)
    expected = spline.distance(data)
    original_distance = Spline.distance
    with cache_native_distance_points():
        actual = spline.distance(data)
        for a, b in zip(actual, expected):
            np.testing.assert_array_equal(a, b)
        assert "get_points" not in spline.__dict__
    assert Spline.distance is original_distance


@pytest.mark.parametrize("accept_native_merges", [False, True])
def test_real_nearby_bifurcations_require_opt_in_and_preserve_all_tips(tmp_path, accept_native_merges):
    path = DATA / "vascularmd_nearby_bifurcations.swc"
    report = process_file(path, tmp_path, Options(surface=False, accept_native_merges=accept_native_merges))
    assert report["original_file_preserved"]
    if not accept_native_merges:
        assert report["status"] == "failed"
        assert report["failed_components"][0]["error_type"] == "TopologyChanged"
        assert not list(tmp_path.glob("*_smooth.swc"))
        return
    assert report["status"] == "success", report["failed_components"]
    assert report["round_trip_valid"] and report["native_model_topology_preserved"]
    assert not report["topology_preserved"]
    assert report["topology"]["all_original_terminals_preserved"]
    assert report["topology"]["merged_bifurcation_swc_ids"] == [10, 12, 489]
    raw, smooth = (report["comparison"][name] for name in ("raw", "smooth"))
    assert raw["terminal_count"] == smooth["terminal_count"] == 11
    assert raw["bifurcation_count"] == 10 and smooth["bifurcation_count"] == 7
    # Raw global QC must count source edges only once, including shared merged trunks.
    graph = read_source(path).graph
    expected_pairs = sum(graph.out_degree(a) <= 1 and graph.out_degree(b) <= 1 for a, b in graph.edges)
    assert raw["log_radius_jump"]["number_of_evaluated_pairs"] == expected_pairs
    assert smooth["log_radius_jump"]["p95"] < raw["log_radius_jump"]["p95"]
    exported = Path(report["outputs"]["swc"]).read_text()
    assert "accept-native-merges" in exported


def test_recursive_batch_continues_after_bad_input_and_excludes_own_outputs(tmp_path):
    inputs = tmp_path / "inputs"
    (inputs / "nested").mkdir(parents=True)
    (inputs / "a_bad.swc").write_text("1 1 0 0 0 1 -1\n2 3 10 0 0 0 1\n")
    (inputs / "nested/good.swc").write_bytes((DATA / "vascularmd_spike.swc").read_bytes())
    output = inputs / "results"
    result = subprocess.run([sys.executable, str(ROOT / "tools/vascularmd_smooth_swc.py"), str(inputs),
                             "--recursive", "--output-dir", str(output), "--no-surface"], capture_output=True, text=True)
    assert result.returncode == 1
    assert (output / "a_bad_vmd_qc.json").exists()
    assert (output / "nested/good_vmd_smooth.swc").exists()
    summary = json.loads(next(output.glob("vascularmd_batch_*.json")).read_text())
    assert len(summary) == 2
    assert [record["status"] for record in summary] == ["failed", "success"]


def test_topology_failure_can_keep_native_surface_only(tmp_path, monkeypatch):
    import vascular_processing.pipeline as pipeline
    from vascular_processing.vascularmd_adapter import TopologyChanged

    def changed(*args, **kwargs):
        raise TopologyChanged({"stage": "topology_validation", "affected_source_swc_nodes_and_ancestors": [21]})

    monkeypatch.setattr(pipeline, "sample_model", changed)
    report = process_file(DATA / "vascularmd_bifurcation.swc", tmp_path, Options())
    assert report["status"] == "partial" and report["needs_manual_review"]
    assert not report["topology_preserved"]
    assert report["surface"]["diagnostic_only"]
    assert Path(report["outputs"]["surface_vtk"]).exists()
    assert not list(tmp_path.glob("*_smooth.swc"))


@pytest.mark.parametrize("kwargs", [{"circumferential_n": 10}, {"circumferential_n": 4},
                                    {"sample_mode": "spacing"}, {"sample_mode": "spacing", "spacing_mm": float("nan")},
                                    {"longitudinal_density": 0}, {"spacing_mm": 1}])
def test_invalid_settings_fail_before_processing(kwargs):
    with pytest.raises(ValueError):
        Options(**kwargs).validate()
