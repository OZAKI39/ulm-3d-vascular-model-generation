"""Native branch fits retain even close junctions, all tips and source endpoints."""
import numpy as np
import pytest

from vascular_processing.pipeline import Options, process_file
from vascular_processing.swc_export import read_source
from third_party.vascularmd.ArterialTree import ArterialTree
from third_party.vascularmd.Spline import Spline
from pathlib import Path

DATA = Path(__file__).parent / "data"


@pytest.mark.parametrize("name", ["vascularmd_bifurcation.swc", "vascularmd_nearby_bifurcations.swc"])
def test_branch_aic_keeps_connections_types_and_all_endpoints(tmp_path, monkeypatch, name):
    calls = []
    original = Spline.approximation

    def observed(self, data, constraint, values, derivatives, **kwargs):
        calls.append((constraint.copy(), kwargs.copy()))
        return original(self, data, constraint, values, derivatives, **kwargs)

    def forbidden(*args, **kwargs):
        pytest.fail("branch mode must not call Nfurcation network reconstruction")

    monkeypatch.setattr(Spline, "approximation", observed)
    monkeypatch.setattr(ArterialTree, "model_network", forbidden)
    source = read_source(DATA / name)
    report = process_file(source.path, tmp_path, Options(model_mode="branches", surface=False))
    assert report["status"] == "success", report["failed_components"]
    assert report["topology_preserved"]
    assert report["topology"]["merged_bifurcation_swc_ids"] == []
    assert len(calls) == len(report["branch_fits"]) == report["comparison"]["raw"]["branch_count"]
    for constraints, keywords in calls:
        assert constraints == [True, False, False, True]
        assert keywords == dict(radius_model=True, criterion="AIC", akaike=False, max_distance=6)
    output = read_source(Path(report["outputs"]["swc"]))
    mapping = report["source_to_output_node_ids"]
    assert len(output.graph) == len(source.graph)
    assert set(output.graph.edges) == {(mapping[a], mapping[b]) for a, b in source.graph.edges}
    for node in source.graph:
        raw, fitted = source.graph.nodes[node], output.graph.nodes[mapping[node]]
        assert raw["swc_type"] == fitted["swc_type"]
        if source.graph.in_degree(node) != 1 or source.graph.out_degree(node) != 1:
            np.testing.assert_array_equal(raw["coords"], fitted["coords"])
    for fit in report["branch_fits"]:
        assert fit["output_count"] == fit["original_count"]
        assert fit["fit_count"] >= 6
        assert fit["dense_radius_min"] > 0
        assert np.isfinite([fit["spatial_lambda"], fit["radius_lambda"]]).all()


def test_branch_mode_does_not_claim_a_native_surface_or_merge():
    for kwargs in [dict(surface=True), dict(accept_native_merges=True),
                   dict(sample_mode="spacing", spacing_mm=0.5)]:
        with pytest.raises(ValueError, match="branches mode requires"):
            Options(**(dict(model_mode="branches", surface=False) | kwargs)).validate()


def test_invalid_native_branch_result_fails_without_radius_clipping(tmp_path, monkeypatch):
    original = Spline.approximation

    def invalid(self, *args, **kwargs):
        original(self, *args, **kwargs)
        samples = self.get_points().copy()
        samples[2, 3] = -1
        self.get_points = lambda: samples

    monkeypatch.setattr(Spline, "approximation", invalid)
    report = process_file(DATA / "vascularmd_spike.swc", tmp_path, Options(model_mode="branches", surface=False))
    assert report["status"] == "failed"
    assert "swc" not in report["outputs"]
    assert not list(tmp_path.glob("*_smooth.swc"))
    assert "Nonpositive" in report["failed_components"][0]["message"]
    assert report["failed_components"][0]["component"]["source_branch"]
