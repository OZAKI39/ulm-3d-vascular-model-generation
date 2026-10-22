"""Check saved-state alignment and finite pair-norm diagnostics before sealing."""

import json
from pathlib import Path

import h5py
import numpy as np

from scripts.stage4m.common import BASE, write


def validate():
    comparison = json.loads((BASE / "comparisons/DENSE_COMPARISON.json").read_text())
    rows = []
    for case, policies in comparison.items():
        for policy, metrics in policies.items():
            run = Path(metrics["summary"]["absolute_run"])
            diagnostics = (
                (BASE / "comparisons/dense" / case / "A0") if policy == "A0" else run
            ) / "PAIR_FORCE_DIAGNOSTICS.jsonl"
            data = [json.loads(line) for line in diagnostics.read_text().splitlines()]
            with h5py.File(run / "trajectories.h5", "r") as h:
                t, active, radii = h["time_s"][:], h["active"][:], h["radii_m"][:]
            assert len(t) == len(data) == metrics["summary"]["accepted_steps"] + 1
            np.testing.assert_allclose([d["time_s"] for d in data], t, rtol=0, atol=1e-13)
            assert np.all(np.diff(t) > 0)
            assert np.all(radii == radii[0])
            for d, mask in zip(data, active):
                n = int(mask.sum())
                count = n * (n - 1) // 2
                assert d["unordered_pair_count"] == count
                assert (
                    sum(
                        d[k] for k in ["full_weight_pair_count", "taper_zone_pair_count", "zeroed_pair_count"]
                    )
                    == count
                )
                classical, effective = d["S_classical_N"], d["S_effective_N"]
                assert np.isfinite([classical, effective]).all()
                assert 0 <= effective <= classical * (1 + 2e-14)
                if classical == 0:
                    assert (
                        d["effective_classical_ratio"] is None
                        and d["removed_classical_force_norm_fraction"] is None
                    )
                else:
                    assert abs(d["effective_classical_ratio"] - effective / classical) <= 2e-14
                    assert (
                        abs(d["removed_classical_force_norm_fraction"] - (1 - effective / classical)) <= 2e-14
                    )
                if policy == "A0":
                    assert d["full_weight_pair_count"] == count
                if policy == "A1":
                    assert d["taper_zone_pair_count"] == 0
            rows.append(
                dict(
                    case=case,
                    policy=policy,
                    checked_saved_states=len(t),
                    force_diagnostics_match_saved_time_and_active_pair_count=True,
                    force_norms_finite=True,
                    zero_denominator_semantics_verified=True,
                )
            )
    result = dict(
        status="PASS",
        rows=rows,
        scope="Read-only output validation; no trajectory integration or coordinate changes.",
    )
    write(BASE / "final_decision/ARTIFACT_VALIDATION.json", result)
    print(json.dumps(result), flush=True)
    return result


if __name__ == "__main__":
    validate()
