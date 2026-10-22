"""Finalize legacy analysis metadata only from existing closed evidence, preserving originals."""
import hashlib
import json
from pathlib import Path

from mevo_mb.io.artifacts import write_json


def finish(run, status, evidence):
    path = run/"metadata.json"
    previous = json.loads(path.read_text())
    if previous["status"] != "CREATED":
        return
    (run/"logs/metadata_before_finalization.json").write_bytes(path.read_bytes())
    write_json(path, {**previous, "status": status, "status_finalized_from_evidence": str(evidence),
                      "evidence_sha256": hashlib.sha256(evidence.read_bytes()).hexdigest(),
                      "note": "Bookkeeping reconciliation after completed outputs; scientific arrays unchanged. Original metadata retained in logs."})


def main():
    root = Path("results")
    for p in root.glob("grid*_STATE.json"):
        state = json.loads(p.read_text())
        if state["status"] == "COMPLETE":
            assert Path(state["refined_field"]).is_file()
            finish(Path(state["run"]), "GRID_FIELD_COMPLETE", p)
    for p in root.glob("*_grid_comparison/acoustic/GRID_CONVERGENCE.json"):
        assert (p.parents[1]/"figures/grid_convergence.pdf").is_file()
        finish(p.parents[1], "DIAGNOSTIC_COMPLETE", p)
    for p in root.glob("*_convergence/CONVERGENCE.json"):
        assert (p.parent/"figures/timestep_convergence.pdf").is_file()
        finish(p.parent, "ANALYSIS_COMPLETE", p)


if __name__ == "__main__":
    main()
