from __future__ import annotations

import pickle

import pytest

from schmid_test_data import write_synthetic_schmid
from utils.schmid_pkl.loader import RestrictedNumpyUnpickler, load_schmid_input


def test_restricted_loader_reads_expected_schmid_schema(tmp_path) -> None:
    input_dir = write_synthetic_schmid(tmp_path / "NW1_results")
    source = load_schmid_input(input_dir)

    assert source.vertex_count == 5
    assert source.edge_count == 5
    assert source.point_sequences_um[0].shape == (2, 3)
    assert source.source_files["verticesDict.pkl"]["sha256"]


def test_restricted_loader_blocks_unexpected_globals(tmp_path) -> None:
    path = tmp_path / "unsafe.pkl"
    with path.open("wb") as stream:
        pickle.dump(len, stream, protocol=2)
    with path.open("rb") as stream:
        with pytest.raises(pickle.UnpicklingError, match="Blocked global"):
            RestrictedNumpyUnpickler(stream).load()
