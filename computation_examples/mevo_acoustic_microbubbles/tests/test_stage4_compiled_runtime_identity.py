"""Passing short trajectories must not authorize a different compiled binary."""
import json
from pathlib import Path

import pytest


def test_compiled_backend_rejects_binary_drift(tmp_path,monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]))
    from scripts.stage4i5 import common
    from scripts.stage4i5.variants import validate_runtime_identity
    monkeypatch.setattr(common,"ROOT",tmp_path)
    runtime=dict(backend="P1_ordered",ordered_contact_sweep=dict(library_sha256="verified",numpy_norm=dict(library_sha256="blas")))
    rows=[]
    for i in range(4):
        folder=tmp_path/str(i)
        folder.mkdir()
        (folder/"SUMMARY.json").write_text(json.dumps(dict(backend_info=runtime)))
        rows.append(dict(candidate=str(i)))
    decision=dict(candidates=[dict(backend="P1_ordered",strict_correctness=True,regression=rows)])
    validate_runtime_identity(runtime,decision)
    different=dict(runtime,ordered_contact_sweep=dict(library_sha256="changed",numpy_norm=dict(library_sha256="blas")))
    with pytest.raises(ValueError,match="identity differs"):
        validate_runtime_identity(different,decision)
