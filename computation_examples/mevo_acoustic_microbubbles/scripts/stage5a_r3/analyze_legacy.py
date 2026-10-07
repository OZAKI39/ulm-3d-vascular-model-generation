"""Attribute source, clock and device differences using controlled actual histories."""

from pathlib import Path
import numpy as np
from scripts.stage5a_r3.common import BASE, R2, read, write
from scripts.stage5a_r3.lineage import compare_all


def main():
    rows = []
    for case in ["B2", "B3"]:
        group = np.load(R2 / "comparisons" / f"{case}_0.001_v4a_legacy_reference.npz")["particle_ids"]
        folder = BASE / "diagnosis/legacy_attribution_v1"
        pairs = [
            ("source_only", "archived_500us_auto", "fresh_500us_auto"),
            ("clock_only", "fresh_500us_auto", "fresh_10us_auto"),
            ("device_setting_only", "fresh_10us_auto", "fresh_10us_cupy"),
            ("legacy_precision_01", "precision_L0", "precision_L1"),
            ("legacy_precision_12", "precision_L1", "precision_L2"),
        ]
        for name, left, right in pairs:
            a = folder / (case + "_" + left) / "trajectories.h5"
            b = folder / (case + "_" + right) / "trajectories.h5"
            row = compare_all(case, "attribution_" + name, a, b, group)
            row["controlled_factor"] = name
            rows.append(row)
        for label in ["B2_initial" if case == "B2" else "B3_initial"]:
            runs = read(BASE / "diagnosis" / ("LOCAL_REFINEMENT_" + label + ".json"))
            for r in runs:
                if r["level"] != 2:
                    continue
                a = folder / (case + "_precision_L2") / "trajectories.h5"
                b = Path(r["run"]) / "trajectories.h5"
                rows.append(compare_all(case, "legacy_fine_vs_R3_fine", a, b, group))
    existing = read(BASE / "final/REFERENCE_LINEAGE_AUDIT.json")
    backup = BASE / "final/REFERENCE_LINEAGE_AUDIT_PRE_R3.json"
    if not backup.exists():
        write(backup, existing)
    existing["controlled_attribution"] = rows
    existing["common_state_replay"] = read(BASE / "diagnosis/COMMON_STATE_REPLAY_V1.json")
    existing["status"] = "BASELINE_DISAGREEMENT"
    existing["promotion"] = "NOT_PERMITTED: independent global per-bubble error bound not established"
    write(BASE / "final/REFERENCE_LINEAGE_AUDIT.json", existing)
    print(
        [(r["case"], r["name"], r["groups"]["all_active"]["errors_m"]["raw"]["endpoint"]) for r in rows],
        flush=True,
    )


if __name__ == "__main__":
    main()
