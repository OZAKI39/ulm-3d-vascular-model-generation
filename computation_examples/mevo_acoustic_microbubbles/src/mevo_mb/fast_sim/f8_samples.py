"""F8 task identity with causal control history and explicitly separated run audit."""

import json
import h5py
from .ai_samples import read_group_sample


def read_navigation_sample(path, frame, history_frames=3):
    sample = read_group_sample(path, frame, history_frames)
    now = float(sample["trajectory"]["time_s"][-1])
    run_audit = sample.pop("quality")
    # The inherited API exposes all future controls for audit. F8 features do not.
    sample["control"] = [e for e in sample["control"] if e["time_s"] <= now + 1e-13]
    with h5py.File(path) as h:
        cfg = json.loads(h.attrs["config_json"])
        env = cfg["environment"]
        first = max(0, frame - history_frames + 1)
        gate = json.loads(h["trajectory"].attrs.get("gate_crossing_json", "null"))
        sample["task"] = dict(
            task_type=h.attrs["task_type"],
            target_m=env["target_m"],
            x_gate_m=env.get("x_gate_m"),
            target_tolerance_m=env.get("target_tolerance_m"),
            prediction_mobility=h.attrs["prediction_mobility"],
            production_mobility=h.attrs["production_mobility"],
            contact_model=h.attrs["contact_model"],
            numerical_dt_s=cfg["dt_s"],
            gate_observed=gate if gate and gate.get("time_s", float("inf")) <= now else None,
        )
        sample["quality"] = dict(
            observed_raw=h["quality/raw_F1_diagonal"][first : frame + 1],
            observed_wall=h["wall_stats"][first : frame + 1],
            run_audit_only=run_audit,
            restriction="run_audit_only contains end-of-run labels/quality; never use as causal input",
        )
    sample["training_ready"] = False
    sample["split_restriction"] = (
        "Keep every F5/F6/F7/F8 paired initial-state family together; no independent third family or validation set."
    )
    return sample
