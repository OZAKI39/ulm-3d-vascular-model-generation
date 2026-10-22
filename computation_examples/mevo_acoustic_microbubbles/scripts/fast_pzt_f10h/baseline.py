"""Audit complete F10G delivery and original scientific inputs without re-running it."""

import subprocess
import h5py
import numpy as np
from .common import ROOT, BASE, BASELINE, CASES, state, js, sha, write


def main():
    BASE.mkdir(parents=True, exist_ok=True)
    f = ROOT / "results/fast_pzt/f10g"
    manifest = js(f / "ARTIFACT_MANIFEST.json")["files"]
    for name, item in manifest.items():
        assert sha(f / name) == item["sha256"], name
    inputs = js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json")
    for name, item in inputs.items():
        assert sha(ROOT / "results/fast_pzt/f5" / name) == item["sha256"], name
    history = []
    for case in CASES:
        _, t, path = state(case)
        assert not t["accepted"] and np.array_equal(t["initial_positions_m"], t["positions_after_step_m"])
        hp = f / f"integration_v1/{case}_B2_K2/trajectory.h5"
        with h5py.File(hp) as h:
            assert h["positions_m"].shape[2:] == (32, 3)
            assert not h.attrs["training_ready"]
            history.append(
                dict(
                    case=case,
                    failure_trace_sha256=sha(path),
                    hdf5_sha256=sha(hp),
                    last_saved_time_s=float(h["time_s"][-1]),
                )
            )
    source = js(f / "REMOTE_VERIFICATION.json")
    for commit in source["source_files_against_git_archive"]:
        subprocess.run(["git", "cat-file", "-e", commit + "^{commit}"], check=True, cwd=ROOT)
    selection = js(ROOT / "results/fast_pzt/f10e/SELECTION.json")
    for s in [selection["branch"], selection["preparation"]]:
        assert sha(ROOT / s["relative_path"] / "trajectory.h5") == s["sha256"]
    write(
        BASE / "BASELINE_AUDIT.json",
        dict(
            status="PASS",
            baseline=BASELINE,
            f10g_verified_files=len(manifest),
            inputs=inputs,
            history=history,
            execution_commits=list(source["source_files_against_git_archive"]),
            original_selection=selection,
            worktrees=subprocess.check_output(["git", "worktree", "list"], cwd=ROOT, text=True),
            applicable_local_agents=[],
            remote_guide_sha256="1d264a8e24863e481852969f4bedd2834abfed79c68dab66d13e9d61a2a7360b",
        ),
    )
    print("F10G audit PASS", len(manifest), "artifacts;", len(inputs), "inputs")


if __name__ == "__main__":
    main()
