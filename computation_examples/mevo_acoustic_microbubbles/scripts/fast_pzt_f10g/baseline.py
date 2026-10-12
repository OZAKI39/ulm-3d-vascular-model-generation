"""Verify existing F10F artifacts and inputs before starting any new experiment."""

import subprocess
import numpy as np
from .common import ROOT, BASE, BASELINE, CASES, state, js, sha, write


def main():
    BASE.mkdir(parents=True, exist_ok=True)
    f = ROOT / "results/fast_pzt/f10f"
    manifest = js(f / "ARTIFACT_MANIFEST.json")["files"]
    for name, item in manifest.items():
        assert sha(f / name) == item["sha256"], name
    inputs = js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json")
    for name, item in inputs.items():
        assert sha(ROOT / "results/fast_pzt/f5" / name) == item["sha256"], name
    rows = []
    for case in CASES:
        data, t, p = state(case)
        assert not t["accepted"] and t["dt_s"] == 1e-5
        assert np.array_equal(t["initial_positions_m"], t["positions_after_step_m"])
        rows.append(
            dict(
                case=case,
                path=str(p.relative_to(ROOT)),
                sha256=sha(p),
                start_s=t["start_s"],
                end_s=t["attempted_end_s"],
                physics=data["physics"],
                quality=t["quality"],
            )
        )
    write(
        BASE / "BASELINE_AUDIT.json",
        dict(
            baseline=BASELINE,
            f10f_verified_files=len(manifest),
            inputs_verified=len(inputs),
            inputs=inputs,
            historical_sha256={
                str(p.relative_to(ROOT)): sha(p)
                for s in ["f10e", "f10f"]
                for p in (ROOT / "results/fast_pzt" / s).rglob("*")
                if p.is_file()
            },
            cases=rows,
            worktrees=subprocess.check_output(["git", "worktree", "list"], cwd=ROOT, text=True),
            applicable_local_agents=[],
            remote_guide_sha256="1d264a8e24863e481852969f4bedd2834abfed79c68dab66d13e9d61a2a7360b",
        ),
    )
    print("BASELINE PASS", len(manifest), len(inputs))


if __name__ == "__main__":
    main()
