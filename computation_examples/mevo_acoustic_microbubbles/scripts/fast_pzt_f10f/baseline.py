"""Read-only local inventory; captures historical hashes before F10F experiments."""

import gzip
import h5py
import numpy as np
import subprocess
from .common import BASE, HISTORY, ROOT, BASELINE, js, sha, write


def main():
    manifest = js(HISTORY / "ARTIFACT_MANIFEST.json")
    print("manifest schema:", type(manifest), list(manifest)[:8])
    files = {str(f.relative_to(ROOT)): sha(f) for f in HISTORY.rglob("*") if f.is_file()}
    trials = []
    for folder in ["contact_v1/prepare_1104_C0", "contact_v1/O_1000_1102_C0"]:
        path = HISTORY / folder
        with gzip.open(path / "CONTACT_TRACE.json.gz", "rt") as s:
            data = __import__("json").load(s)
        last = data["traces"][-1]
        with h5py.File(path / "trajectory.h5") as h:
            proposal = h["rejected_interval/proposal_m"][0]
            assert np.array_equal(proposal, np.array(last["proposal_positions_m"]))
            trials.append(
                dict(
                    path=str(path.relative_to(ROOT)),
                    sha256=sha(path / "trajectory.h5"),
                    first_rejection=js(path / "OBSERVER.json")["first_rejection_step"],
                    last_saved_s=float(h["time_s"][-1]),
                    sample_id=js(path / "SUMMARY.json")["sample_id"],
                    rejected_proposal_matches_trace=True,
                    source=js(path / "SUMMARY.json")["source"],
                )
            )
    gate = js(HISTORY / "analysis/FINAL_GATE.json")
    assert gate["remote_hashes"] == "PASS" and gate["original_quality_budgets_preserved"]
    write(
        BASE / "BASELINE_AUDIT.json",
        dict(
            baseline_commit=BASELINE,
            f10e_gate=gate["status"],
            baseline_trials=trials,
            historical_files_sha256=files,
            worktrees=subprocess.check_output(["git", "worktree", "list"], cwd=ROOT, text=True),
            applicable_local_agents=[],
            server_guide_sha256="1d264a8e24863e481852969f4bedd2834abfed79c68dab66d13e9d61a2a7360b",
            original_f10e_clean=True,
        ),
    )
    print("PASS", len(files), "historical files")


if __name__ == "__main__":
    main()
