"""Three matched mature-state runs per numerical mode; no new physics variants."""

import argparse
import itertools
import numpy as np
from mevo_mb.fast_sim.morphodynamics import contrast
from scripts.fast_pzt_f2.analyze import read, adjacency
from .common import BASE, F2, PROTOCOL, write, js, identity
from .experiment import run_case


def pair_comparison(a, b):
    assert np.array_equal(a["time_s"], b["time_s"])
    rows = [contrast(x, y, np.arange(32)) for x, y in zip(a["positions_m"][:, 0], b["positions_m"][:, 0])]
    result = {k: np.asarray([r[k] for r in rows]) for k in rows[0]}
    result["neighbor_disagreement"] = np.sum(adjacency(a)["linked"] != adjacency(b)["linked"], axis=1)
    result["max_ID_position_difference_m"] = np.max(
        np.linalg.norm(a["positions_m"][:, 0] - b["positions_m"][:, 0], axis=-1), axis=1
    )
    return result


def main(label, checks, driving):
    assert js(BASE / checks / "RESULT.json")["status"] == "PASS"
    assert js(BASE / driving / "RESULT.json")["complete"]
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    base = js(F2 / "mature_v1/unchanged/CONFIG.json")
    base["capture_initial_force"] = True
    rows = {}
    modes = {}
    for mode in ["fast_atomic", "stable_reduction"]:
        paths = []
        for i in range(PROTOCOL["repeat_count_per_mode"]):
            name = f"{mode}_{i}"
            cfg = dict(base, execution_mode=mode)
            rows[name] = run_case(out / name, cfg)
            paths.append(out / name)
            write(out / "PARTIAL.json", dict(completed=list(rows), rows=rows))
        data = [read(p) for p in paths]
        comparisons = []
        for i, j in itertools.combinations(range(3), 2):
            c = pair_comparison(data[i], data[j])
            comparisons.append(c)
            np.savez_compressed(out / f"{mode}_pair_{i}_{j}.npz", time_s=data[0]["time_s"], **c)
        envelope = {k: np.max([c[k] for c in comparisons], axis=0) for k in comparisons[0]}
        np.savez_compressed(out / f"{mode}_envelope.npz", time_s=data[0]["time_s"], **envelope)
        modes[mode] = dict(
            position_bitwise_equal=bool(
                all(np.array_equal(data[0]["positions_m"], d["positions_m"]) for d in data[1:])
            ),
            raw_sorted_edge_events_identical=bool(
                all(
                    np.array_equal(data[0]["edge_ids"], d["edge_ids"])
                    and np.array_equal(data[0]["edge_times"], d["edge_times"])
                    for d in data[1:]
                )
            ),
            maximum_over_time={k: float(v.max()) for k, v in envelope.items()},
            endpoint={k: np.asarray(v[-1]).tolist() for k, v in envelope.items()},
            walls_s=[rows[p.name]["wall_s"] for p in paths],
            quality=[rows[p.name]["numerical_quality"] for p in paths],
        )
    qualified = modes["stable_reduction"]["position_bitwise_equal"] and all(
        v == "PASS" for v in modes["stable_reduction"]["quality"]
    )
    result = dict(
        complete=True,
        status="PASS_SAME_DEVICE_BITWISE" if qualified else "NOT_QUALIFIED",
        identity=identity(),
        modes=modes,
        rows=rows,
        selected_execution_mode="stable_reduction" if qualified else "fast_atomic",
        repeat_scope="Three same-device runs; maximum over three repeat-pair contrasts, not a confidence interval; frozen physics",
    )
    write(out / "RESULT.json", result)
    print("REPEATABILITY", result["status"], {n: v["endpoint"] for n, v in modes.items()}, flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--checks", required=True)
    ap.add_argument("--driving", required=True)
    a = ap.parse_args()
    main(a.label, a.checks, a.driving)
