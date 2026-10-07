"""R3 identities and read-only R2 inputs; never redirects R2 output roots."""

import json
import os
from pathlib import Path
import numpy as np
from scripts.stage4l.dense_probe import load_case, sha, write as write

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage5/motion_solver_r3"
R2ROOT = Path(os.environ.get("MEVO_R2_ROOT", "/workspace/mevo_acoustic_microbubbles_stage5a_r2"))
if not R2ROOT.exists():
    R2ROOT = Path(
        "/home/lzy/projects/computation_examples/_archives/mevo_stage5a_r2_20261006/computation_examples/mevo_acoustic_microbubbles"
    )
R2 = R2ROOT / "results/stage5/motion_solver_r2"
B2_RUNS = {
    "nominal": "20261007T004943915783Z_B2_S1_0.005s_validated_v4a_0.005_nominal_repeat0",
    "refined": "20261007T011553819251Z_B2_S1_0.005s_validated_v4a_0.005_refined_repeat0",
}


def read(path):
    p = Path(path)
    return (
        [json.loads(s) for s in p.read_text().splitlines() if s]
        if p.suffix == ".jsonl"
        else json.loads(p.read_text())
    )


def load(case, pair="numpy", components=True, workers=1):
    from mevo_mb.bubbles.tracker import ParticleTracker
    from mevo_mb.bubbles.motion_sweep_r2 import WallSweep

    cfg, field, flow, gold, inputs = load_case(Path("/workspace/mevo_acoustic_microbubbles"), case)
    selection = read(R2ROOT / "results/stage5/motion_solver/audit_v2/BASELINE_SELECTION.json")
    snapshot = Path(next(r["run"] for r in selection if r["key"] == case + "_1ms")) / "INPUT_SNAPSHOT.npz"
    inputs[str(snapshot)] = sha(snapshot)
    with np.load(snapshot) as f:
        state = {k: f[k] for k in f.files}
    cfg["particles"].update(
        translation_backend="constrained_lubrication_r2",
        pair_backend=pair,
        response_backend="numpy",
        lubrication=dict(activation=0.1, floor=0.001, residual_tolerance=1e-11),
        constrained_solver=dict(solver="prox", components=components, component_workers=workers),
    )
    tracker = ParticleTracker(cfg, field, flow)
    wall = Path(
        "/workspace/mevo_acoustic_microbubbles/results/stage4/simvascular/mesh/20261002T141441884393Z/SV_MESH/mesh-surfaces/WALL.vtp"
    )
    inputs[str(wall)] = sha(wall)
    sweep = WallSweep(tracker.cylinder, wall)
    return tracker, sweep, cfg, state, inputs


def identity():
    return dict(
        source_commit=(ROOT / "SOURCE_COMMIT_R3.txt").read_text().strip()
        if (ROOT / "SOURCE_COMMIT_R3.txt").exists()
        else "UNCOMMITTED_DEVELOPMENT",
        source_sha256={
            str(p.relative_to(ROOT)): sha(p)
            for folder in ["src", "scripts/stage5a_r3"]
            for p in (ROOT / folder).rglob("*.py")
        },
        environment_sha256=sha(BASE / "environment/IDENTITY.json"),
        plan_sha256=sha(ROOT / "docs/STAGE5A_R3_PLAN.md"),
    )
