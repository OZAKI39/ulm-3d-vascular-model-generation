from __future__ import annotations

import json

import networkx as nx
import numpy as np

from schmid_test_data import write_synthetic_schmid
from utils.schmid_pkl.config import SchmidPKLConfig
from utils.schmid_pkl.pipeline import run_schmid_pkl_pipeline


def test_small_directed_pipeline_writes_reopenable_acceptance_outputs(tmp_path) -> None:
    input_dir = write_synthetic_schmid(tmp_path / "NW1_results")
    config = SchmidPKLConfig(
        input_dir=input_dir,
        output_root=tmp_path / "outputs",
        visualizations_enabled=False,
        smoothing_enabled=False,
    )
    run = run_schmid_pkl_pipeline(config)

    assert run.status == "completed"
    assert run.acceptance.overall_status == "PASS"
    assert run.html_report.is_file()
    graph_json = run.run_root / "graphs" / "directed_hierarchical_vascular_graph.json"
    payload = json.loads(graph_json.read_text(encoding="utf-8"))
    assert payload["representation"]["directed"] is True
    assert payload["summary"]["branch_count"] == 3
    graph = nx.read_graphml(run.run_root / "graphs" / "directed_junction_branch_graph.graphml")
    assert graph.is_directed()
    with np.load(run.run_root / "graphs" / "branch_geometry.npz") as arrays:
        assert arrays["branch_ids"].shape == (3,)
