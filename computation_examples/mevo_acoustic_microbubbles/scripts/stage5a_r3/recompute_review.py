"""Reproduce published all-ID error statistics from compact actual coordinates."""

import numpy as np
from scripts.stage5a_r3.common import ROOT, BASE, read, write, sha
from scripts.stage5a_r3.lineage import metrics, summarize


def main():
    index = read(BASE / "review/ALL_ACTIVE_COMPACT_INDEX.json")
    rows = []
    for row in index["rows"]:
        path = ROOT / row["compact"]
        assert sha(path) == row["compact_sha256"]
        a = np.load(path)
        masks = {
            "all_active": np.ones(len(a["particle_id"]), bool),
            "original_group": np.isin(a["particle_id"], a["original_group_particle_id"]),
        }
        masks["complement"] = ~masks["original_group"]
        maximum_difference = 0.0
        for group, mask in masks.items():
            errors = metrics(a["left_positions_m"][:, mask], a["right_positions_m"][:, mask])
            for mode, values in errors.items():
                observed = summarize(values[-1])
                expected = row["expected_groups"][group]["errors_m"][mode]["endpoint"]
                for key in observed:
                    delta = abs(observed[key] - expected[key])
                    maximum_difference = max(maximum_difference, delta)
                    # Reproduction rounding allowance, NOT a physical accuracy
                    # or solver tolerance. Raw arrays and IDs are unmodified.
                    assert np.isclose(observed[key], expected[key], rtol=1e-8, atol=1e-15), (
                        row["name"],
                        group,
                        mode,
                        key,
                    )
        rows.append(
            dict(
                case=row["case"],
                comparison=row["name"],
                status="PASS",
                maximum_statistic_reproduction_difference_m=maximum_difference,
                active_count=len(a["particle_id"]),
                common_frames=len(a["clock_s"]),
            )
        )
    write(
        BASE / "review/COMPACT_DATA_REPRODUCTION.json",
        dict(
            status="PASS",
            compact_index_sha256=sha(BASE / "review/ALL_ACTIVE_COMPACT_INDEX.json"),
            analysis_script_sha256=sha(__file__),
            rows=rows,
            scope="Data/statistics reproduction only. No claim of an independent physical trajectory error bound.",
        ),
    )
    print(dict(compact_comparisons=len(rows), status="PASS"))


if __name__ == "__main__":
    main()
