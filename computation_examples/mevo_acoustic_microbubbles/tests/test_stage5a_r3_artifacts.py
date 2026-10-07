"""The lightweight archive must retain the evidence needed to reproduce its figures."""

from scripts.stage5a_r3 import artifacts


def test_review_keeps_lossless_diagnostics_and_local_acceptance_logs(tmp_path, monkeypatch):
    monkeypatch.setattr(artifacts, "BASE", tmp_path)
    for relative in [
        "review/compact_dynamics/DYNAMIC_INSTRUMENTED_initial.json.gz",
        "review/repaired_B2_terminal/nominal/ACCEPTED.jsonl",
        "review/repaired_B2_terminal/nominal/REJECTED.jsonl",
    ]:
        path = tmp_path / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"{}\n")
        assert artifacts.selected(path)
    raw = tmp_path / "diagnosis/DYNAMIC_INSTRUMENTED_initial.json"
    raw.parent.mkdir()
    raw.write_text("{}\n")
    assert not artifacts.selected(raw)


def test_original_comparison_is_omitted_only_when_review_input_exists(tmp_path, monkeypatch):
    monkeypatch.setattr(artifacts, "BASE", tmp_path)
    original = tmp_path / "comparisons/B2_repeat1_all_active.npz"
    original.parent.mkdir()
    original.write_bytes(b"archive selection does not parse array contents")
    assert artifacts.selected(original)
    compact = tmp_path / "review/all_active_inputs" / original.name
    compact.parent.mkdir(parents=True)
    compact.write_bytes(b"separate reconstruction audit verifies the actual arrays")
    assert artifacts.selected(compact)
    assert not artifacts.selected(original)
