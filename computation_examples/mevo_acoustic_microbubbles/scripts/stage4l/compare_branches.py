"""Branch and actual-step comparisons. Legacy sensitivity is never inherited."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np

from scripts.stage4l.dense_probe import BASE, write


def effective_steps(run):
    summary = json.loads((run / "particles/SUMMARY.json").read_text())
    # RK2 emits initial, midpoint, accepted endpoint. Exit-removal reevaluations
    # duplicate an endpoint time, and are removed before selecting accepted times.
    times = []
    for line in (run / "TRANSLATION_SOLVES.jsonl").open():
        row = json.loads(line)
        if "solver" in row:
            t = row["time_s"]
            if not times or t != times[-1]:
                times.append(t)
    if len(times) != 2 * summary["accepted_steps"] + 1:
        raise ValueError("Cannot reconstruct actual RK2 steps from solve log; do not guess")
    dt = np.diff(np.asarray(times)[::2])
    return dict(
        min_s=float(dt.min()),
        median_s=float(np.median(dt)),
        max_s=float(dt.max()),
        mean_s=float(dt.mean()),
        accepted_steps=len(dt),
        source="Accepted endpoints reconstructed from RK2 solve log; count checked",
    )


def compare(reference, candidate):
    a = json.loads((reference / "BRANCH_METRICS.json").read_text())
    b = json.loads((candidate / "BRANCH_METRICS.json").read_text())
    ia = json.loads((reference / "RUN_IDENTITY.json").read_text())
    ib = json.loads((candidate / "RUN_IDENTITY.json").read_text())
    assert ia["initial_population_sha256"] == ib["initial_population_sha256"]
    if set(a["assignment"]) != set(b["assignment"]):
        raise ValueError("Branch comparison requires identical particle IDs")
    changed = [int(pid) for pid in a["assignment"] if a["assignment"][pid] != b["assignment"][pid]]
    outcomes = ["BRANCH_A", "BRANCH_B", "UNASSIGNED", "WALL_LOSS"]
    switch = {
        origin: {
            dest: sum(
                a["assignment"][pid] == origin and b["assignment"][pid] == dest for pid in a["assignment"]
            )
            for dest in outcomes
        }
        for origin in outcomes
    }
    ea = {e["particle_id"]: e for e in a["crossings"]}
    eb = {e["particle_id"]: e for e in b["crossings"]}
    shared = sorted(set(ea) & set(eb))
    delta = np.array([eb[pid]["crossing_time_s"] - ea[pid]["crossing_time_s"] for pid in shared])
    early = {}
    ma = json.loads((reference / "metrics.json").read_text())["rows"]
    mb = json.loads((candidate / "metrics.json").read_text())["rows"]
    grid = np.linspace(0.0, 0.1, 6)
    for key in [
        "mean_nearest_neighbor_m",
        "near_contact_fraction",
        "largest_cluster_fraction",
        "cluster_count",
    ]:
        va = [(v["time_s"], v[key]) for v in ma if v.get(key) is not None]
        vb = [(v["time_s"], v[key]) for v in mb if v.get(key) is not None]
        ta, ya = zip(*va)
        tb, yb = zip(*vb)
        ya, yb = np.interp(grid, ta, ya), np.interp(grid, tb, yb)
        early[key] = dict(
            relative_L2=float(np.linalg.norm(yb - ya) / max(np.linalg.norm(ya), 1e-30)),
            max_absolute_difference=float(np.max(np.abs(yb - ya))),
        )
    return dict(
        reference=str(reference),
        candidate=str(candidate),
        assignment_agreement=1 - len(changed) / len(a["assignment"]),
        fraction_semantics="f_A and f_B condition on assigned exits; initial_population_fractions include unassigned IDs",
        initial_population_fractions={
            name: {
                key: row[key] / len(row["assignment"])
                for key in ["N_A", "N_B", "N_unassigned", "N_wall_loss"]
            }
            for name, row in [("reference", a), ("candidate", b)]
        },
        reference_backend=ia.get("translation_backend", "legacy_stokes"),
        candidate_backend=ib.get("translation_backend", "legacy_stokes"),
        reference_dt_s=ia["config"]["particles"]["dt_s"],
        candidate_dt_s=ib["config"]["particles"]["dt_s"],
        changed_particle_ids=changed,
        branch_switch_matrix=switch,
        reference_counts={k: a[k] for k in ["N_A", "N_B", "N_unassigned", "N_wall_loss", "f_A", "f_B"]},
        candidate_counts={k: b[k] for k in ["N_A", "N_B", "N_unassigned", "N_wall_loss", "f_A", "f_B"]},
        shared_crossings=len(shared),
        mean_abs_crossing_time_difference_s=float(np.mean(abs(delta))) if len(delta) else None,
        max_abs_crossing_time_difference_s=float(np.max(abs(delta))) if len(delta) else None,
        early_cluster_common_times_s=grid.tolist(),
        early_cluster_metrics=early,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy-root", type=Path, required=True)
    args = parser.parse_args()
    out = BASE / "comparison_legacy" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ_branches")
    out.mkdir(exist_ok=False, parents=True)
    old = json.loads((args.legacy_root / "results/stage4/primary_only/STATE.json").read_text())["rows"]
    nominal = {}
    halves = {}
    for kind, target in [("primary_baseline", nominal), ("timestep_sensitivity", halves)]:
        for path in sorted((BASE / kind).glob("*/SUMMARY.json")):
            d = json.loads(path.read_text())
            if d.get("status") == "COMPLETE" and "legacy_reference" in d:
                case = json.loads((path.parent / "RUN_IDENTITY.json").read_text())["case"]
                target[case] = path.parent
    if set(nominal) != set(["B0", "B1", "B2", "B3"]) or set(halves) != set(["B2", "B3"]):
        raise ValueError("New branch campaign is incomplete")
    legacy_comparison = {
        case: compare(
            args.legacy_root / next(r["run"] for r in old if r["case"] == case and r["seed"] == 20261002), run
        )
        for case, run in nominal.items()
    }
    sensitivity = {
        case: dict(
            comparison=compare(nominal[case], halves[case]),
            nominal_steps=effective_steps(nominal[case]),
            half_steps=effective_steps(halves[case]),
        )
        for case in ["B2", "B3"]
    }
    write(
        out / "SUMMARY.json",
        dict(
            primary_comparison=legacy_comparison,
            timestep_sensitivity=sensitivity,
            gate_status="REQUIRES_SCIENTIFIC_REVIEW",
            overall_stage4_acceptance="NOT_EVALUATED",
        ),
    )
    print(out, flush=True)


if __name__ == "__main__":
    main()
