"""Check causal prefixes and keep incomplete physical validation out of AI labels."""

import numpy as np
from mevo_mb.fast_sim.f9_samples import read_audit_sample
from .common import BASE, ROOT, js, write


def main():
    result = js(BASE / "analysis/RESULTS.json")
    checks, rows = {}, {}
    for key, row in result["cases"].items():
        path = ROOT / row["path"] / "trajectory.h5"
        for frame in [0, 150, 600]:
            sample = read_audit_sample(path, frame)
            now = sample["trajectory"]["time_s"][-1]
            assert all(c["time_s"] <= now + 1e-13 for c in sample["control"])
            ledger = sample["motion_ledger"]
            if ledger is not None:
                m, e = ledger["step_meta"], ledger["wall_contact_events"]
                assert np.all(m[:, 0] + m[:, 1] <= now + 1e-13)
                assert np.all(e[:, 1] + e[:, 2] <= now + 1e-13)
                if frame == 0:
                    assert not len(m) and not len(e)
                    assert np.count_nonzero(ledger["observed_cumulative_m"]) == 0
            assert not sample["training_ready"]
        checks[key] = "PASS"
        rows[key] = dict(
            **result["provenance"][key],
            paired_family="original_Y_F5_300ms",
            training_ready=False,
            numerical_label="NUMERICAL_PASS",
            navigation_label="NAVIGATION_FAIL" if key.startswith("Y0") else "COHESION_ONLY",
            physical_label="PHYSICS_UNCALIBRATED",
            signed_ledger_available="ledger" in row,
            diagnostic_quality_flags=dict(
                STEP_SENSITIVE=result["timestep"][key.split("_")[0]]["STEP_SENSITIVE"],
                WALL_CORRECTION_SIGNIFICANT=row.get("WALL_CORRECTION_SIGNIFICANT"),
                CONTACT_SENSITIVE_direct_group_component=row.get("CONTACT_SENSITIVE_direct_group_component"),
                indirect_contact_effect="NOT_ISOLATED",
                TRANSVERSE_EFFECT_RESOLVED=all(
                    e["TRANSVERSE_EFFECT_RESOLVED"] for e in result["effects"].values()
                )
                if key.startswith("Y2")
                else None,
            ),
        )
    write(
        BASE / "analysis/AI_INTERFACE_CHECKS.json",
        dict(
            status="PASS",
            cases=checks,
            frames_checked=[0, 150, 600],
            scope="Causal ledger/control prefixes; quality.run_audit_only and dataset-level labels must never become causal features. No model trained.",
        ),
    )
    write(
        BASE / "analysis/AI_DATA_MANIFEST.json",
        dict(
            schema="F9_passive_displacement_audit",
            cases=rows,
            independent_new_initial_state_families=0,
            restriction="No training approval; audit reproductions are not independent examples; null flags are unmeasured, not zero.",
        ),
    )


if __name__ == "__main__":
    main()
