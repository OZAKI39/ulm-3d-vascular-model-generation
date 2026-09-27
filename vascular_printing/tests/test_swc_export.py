from pathlib import Path

import networkx as nx
import numpy as np
import pytest

from vascular_processing.swc_export import read_source, write_swc


def test_round_trip_nonconsecutive_ids_types_parent_order_and_full_precision(tmp_path):
    path = tmp_path / "original.swc"
    path.write_text("20 7 3.123456789012345 0 0 1.33 4\n"
                    "4 1 0 0 0 2 -1\n30 5 6 3 0 0.7 20\n90 6 6 -3 0 0.9 20\n")
    source = read_source(path)
    before = path.read_bytes()
    output = tmp_path / "smooth.swc"
    ids = write_swc(output, source.graph, path, "test", "original-count")
    result = read_source(output)
    assert result.rows[:, 0].tolist() == [1, 2, 3, 4]
    assert np.all((result.rows[:, 6] == -1) | (result.rows[:, 6] < result.rows[:, 0]))
    assert nx.is_arborescence(result.graph)
    for node, new in ids.items():
        np.testing.assert_array_equal(result.graph.nodes[new]["coords"], source.graph.nodes[node]["coords"])
        assert result.graph.nodes[new]["swc_type"] == source.graph.nodes[node]["swc_type"]
    assert path.read_bytes() == before
    with pytest.raises(FileExistsError):
        write_swc(output, source.graph, path, "test", "original-count")
    with pytest.raises(ValueError, match="overwrite"):
        write_swc(path, source.graph, path, "test", "original-count")


@pytest.mark.parametrize("text", [
    "1 1 0 0 0 1 -1\n1 3 1 0 0 1 1\n",
    "1 1 0 0 0 1 -1\n2 3 1 0 0 1 99\n",
    "1 1 0 0 0 1 2\n2 3 1 0 0 1 1\n",
    "1 1 0 0 0 1 -1\n2 3 1 0 0 1 -1\n",
    "1 1 0 0 0 1 -1\n2 3 1 0 0 1 2\n",
    "1 1 0 0 0 1 -1\n2 3 1 0 0 0 1\n",
    "1 1 0 0 0 1 -1\n2 3 1 nan 0 1 1\n",
])
def test_invalid_input_rejected_without_repair(tmp_path, text):
    path = tmp_path / "bad.swc"
    path.write_text(text)
    with pytest.raises(ValueError):
        read_source(path)
    assert path.read_text() == text
