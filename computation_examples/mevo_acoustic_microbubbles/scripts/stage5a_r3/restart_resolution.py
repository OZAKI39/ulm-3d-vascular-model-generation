"""Audit numerical restart against the preregistered resolution, retain bitwise FAIL."""

from pathlib import Path
import h5py
import numpy as np
from scripts.stage5a_r3.common import BASE, read, write
from scripts.stage5a_r3.lineage import summarize
from scripts.stage5a_r3.run import restore


def main():
    old = read(BASE / "diagnosis/RESTART_V2.json")
    sources = [Path(old[k]) for k in ["full", "part", "resumed"]]
    arrays = []
    for p in sources:
        with h5py.File(p / "trajectories.h5") as h:
            arrays.append((h["time_s"][:], h["positions_m"][:], h["active"][0], h["radii_m"][:]))
    t = np.r_[arrays[1][0], arrays[2][0][1:]]
    x = np.concatenate([arrays[1][1], arrays[2][1][1:]])
    mask = arrays[0][2]
    diff = np.linalg.norm(x[:, mask] - arrays[0][1][:, mask], axis=2)
    r = float(arrays[0][3][mask][0])
    budget = dict(
        rms_m=0.01 * r, max_m=0.05 * r, source="STAGE5A_R3_PLAN.md; fixed before restart observations"
    )
    errors = summarize(diff)
    same_time = np.array_equal(t, arrays[0][0])
    state, meta, cache = restore(sources[1] / "CHECKPOINT.npz")
    resumed = read(sources[2] / "IDENTITY.json")
    start_config = read(sources[1] / "IDENTITY.json")
    data = read(sources[2] / "SUMMARY.json")
    safety = (
        data["accepted_kkt_max"] <= 1e-9
        and data["accepted_linear_residual_max"] <= 1e-11
        and data["accepted_pair_swept_gap_min_m"] >= -1e-11
        and data["accepted_wall_swept_gap_bound_min_m"] >= -1e-11
        and data["projection_max_m"] == 0
    )
    mapped = (
        resumed["resume"] == str(sources[1] / "CHECKPOINT.npz")
        and resumed["config"] == start_config["config"]
    )
    status = (
        "PASS"
        if (
            same_time
            and safety
            and mapped
            and errors["rms"] <= budget["rms_m"]
            and errors["max"] <= budget["max_m"]
        )
        else "FAIL"
    )
    report = dict(
        status=status,
        bitwise_gate=old,
        errors_m=errors,
        budget=budget,
        identical_accepted_times=same_time,
        numerical_safety=safety,
        checkpoint_input_mapping=mapped,
        native_workspace_serialized=False,
        limitation="ProxQP numeric workspaces are rebuilt; physical state and canonical accepted cache are restored. Bitwise replay is not guaranteed and the original FAIL is retained. This local restart resolution test is not an independent global accuracy proof.",
    )
    write(BASE / "diagnosis/RESTART_V2_RESOLUTION_AUDIT.json", report)
    print(report, flush=True)


if __name__ == "__main__":
    main()
