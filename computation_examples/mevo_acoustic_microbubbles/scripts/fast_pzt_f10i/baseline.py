"""Verify the entire delivered F10H archive before new F10I analysis."""

import subprocess
import h5py
import numpy as np
from .common import BASE, H, ROOT, BASELINE, TRIALS, js, sha, write


def main():
    BASE.mkdir(exist_ok=True)
    manifest = js(H / "ARTIFACT_MANIFEST.json")["files"]
    for name, v in manifest.items():
        assert sha(H / name) == v["sha256"], name
    inputs = js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json")
    for name, v in inputs.items():
        assert sha(ROOT / "results/fast_pzt/f5" / name) == v["sha256"], name
    rows = []
    for mode, p in TRIALS.items():
        with h5py.File(p / "trajectory.h5") as h:
            assert not h.attrs["training_ready"]
            np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
            assert np.all(h["radii_m"][:] == 1.78e-6)
            rows.append(
                dict(
                    mode=mode,
                    path=str(p.relative_to(ROOT)),
                    files={f.name: sha(f) for f in sorted(p.iterdir()) if f.is_file()},
                    saved_frames=len(h["time_s"]),
                    accepted_steps=len(h["substeps/accepted_end_s"]),
                    raw_edge_events=len(h["edge_events/time_s"]),
                )
            )
    gate = js(H / "analysis/FINAL_GATE.json")
    assert gate["B2_TIMESTEP_CONVERGENCE"] == "FAIL" and gate["PRODUCTION_NOT_ADMITTED"]
    write(
        BASE / "BASELINE_AUDIT.json",
        dict(
            status="PASS",
            baseline_commit=BASELINE,
            F10H_artifact_files=len(manifest),
            F10H_manifest_sha256=sha(H / "ARTIFACT_MANIFEST.json"),
            inputs=inputs,
            trials=rows,
            historical_gate={k: v for k, v in gate.items() if k != "source"},
            worktrees=subprocess.check_output(["git", "worktree", "list"], cwd=ROOT, text=True),
            local_agents="No applicable local AGENTS.md in project or ancestors",
            server_guide_sha256="1d264a8e24863e481852969f4bedd2834abfed79c68dab66d13e9d61a2a7360b",
        ),
    )
    print(
        "F10H archive and actual inputs PASS",
        len(manifest),
        len(inputs),
        rows[0]["raw_edge_events"],
        rows[1]["raw_edge_events"],
    )


if __name__ == "__main__":
    main()
