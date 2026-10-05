"""A failed optional backend must never be displayed as an active replay."""
import json
from pathlib import Path


def test_failed_and_completed_replays_not_presented_as_pending(tmp_path,monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]))
    from scripts.stage4 import status
    monkeypatch.setattr(status,"ROOT",tmp_path)
    base=tmp_path/"results/stage4"
    root=base/"performance_refactor/p1_vtk_contact"
    for name in ["failed","completed","pending","prefix"]:
        folder=root/name
        folder.mkdir(parents=True)
        (folder/"PROGRESS.json").write_text(json.dumps(dict(backend="P2_warp",segment="B2_S2",time_s=5.94)))
    (root/"failed/EXECUTION_FAILURE.json").write_text(json.dumps(dict(status="FAIL",error="projection failure")))
    (root/"completed/SUMMARY.json").write_text(json.dumps(dict(status="COMPLETE",backend="P1_static",
        segment=dict(id="B2_S2"),accepted_steps=3162)))
    (root/"prefix/SUMMARY.json").write_text(json.dumps(dict(status="DIAGNOSTIC_PREFIX",backend="P0",
        segment=dict(id="B2_S2"),accepted_steps=5)))
    (base/"CURRENT_EXECUTION_STATE.json").write_text(json.dumps(dict(status="RUNNING",stage="I.5/ordered_replay")))
    result=status.snapshot()["performance_refactor"]
    assert len(result["candidate_failures"])==1
    assert len(result["completed_short_replays"])==1
    assert [Path(r["run"]).name for r in result["unfinished_saved_replays"]]==["pending"]
