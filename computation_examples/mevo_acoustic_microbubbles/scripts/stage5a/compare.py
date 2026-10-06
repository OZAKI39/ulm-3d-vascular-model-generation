"""Actual accepted-state comparison. Never interpolate particle coordinates."""

from pathlib import Path
import argparse
import json
import numpy as np
import h5py
from scipy.spatial.distance import pdist
from scripts.stage4l.analyze_dense import largest, edges
from scripts.stage4l.dense_probe import write, sha

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage5/motion_solver"


def sample(path, clock, members=None):
    with h5py.File(path) as h:
        t = h["time_s"][:]
        initial = h["positions_m"][0]
        mask = h["active"][0]
        r = float(h["radii_m"][0])
        ids = h["particle_id"][:]
        if members is None:
            members = ids[np.flatnonzero(mask)[largest(initial[mask], 2.1 * r)]]
        index = np.array([int(np.flatnonzero(ids == i)[0]) for i in members])
        ix = np.maximum(np.searchsorted(t, clock, side="right") - 1, 0)
        unique, inverse = np.unique(ix, return_inverse=True)
        frames = h["positions_m"][unique][:, index][inverse]
        center = frames.mean(axis=1)
        centered = frames - center[:, None]
        ref = centered[0]
        rigid = []
        for frame in centered:
            u, _, v = np.linalg.svd(ref.T @ frame)
            rotation = u @ np.diag([1, 1, np.linalg.det(u @ v)]) @ v
            rigid.append(np.sqrt(np.mean(np.sum((ref @ rotation - frame) ** 2, axis=1))))
        graph = [edges(x, 2.1 * r) for x in frames]
        hard = [edges(x, 2 * r + 1e-11) for x in frames]
        turn = [len(a ^ b) / max(1, len(a | b)) for a, b in zip(graph[:-1], graph[1:])]
        distances = np.array([np.quantile(pdist(x), [0, 0.25, 0.5, 0.75, 1]) for x in frames])
        return dict(
            time_s=clock,
            actual_time_s=t[ix],
            sample_age_s=clock - t[ix],
            particle_ids=members,
            positions_m=frames,
            centroid_m=center,
            rg_m=np.sqrt(np.mean(np.sum(centered**2, axis=2), axis=1)),
            rigid_fit_m=np.array(rigid),
            pair_quantiles_m=distances,
            neighbor_turnover=np.array(turn),
            hard_edge_count=np.array(list(map(len, hard))),
            near_edge_count=np.array(list(map(len, graph))),
            minimum_pair_gap_m=distances[:, 0] - 2 * r,
            radius_m=np.array(r),
        )


def compare(a_path, b_path, out, label):
    out.mkdir(parents=True, exist_ok=True)
    with h5py.File(a_path) as a, h5py.File(b_path) as b:
        t0 = max(float(a["time_s"][0]), float(b["time_s"][0]))
        end = min(float(a["time_s"][-1]), float(b["time_s"][-1]))
    clock = t0 + np.arange(int(np.floor((end - t0) / 1e-5 + 1e-7)) + 1) * 1e-5
    if len(clock) == 1:
        clock = np.array([t0])
    if abs(clock[-1] - end) < 1e-12:
        clock[-1] = end
    a = sample(a_path, clock)
    b = sample(b_path, clock, a["particle_ids"])
    radius = float(a["radius_m"])
    rg_endpoint = float(abs(b["rg_m"][-1] / a["rg_m"][-1] - 1))
    rg_l2 = float(np.linalg.norm(b["rg_m"] - a["rg_m"]) / np.linalg.norm(a["rg_m"]))
    quantile_rel = abs(b["pair_quantiles_m"][-1, 1:4] / a["pair_quantiles_m"][-1, 1:4] - 1)
    travel = float(np.sum(np.linalg.norm(np.diff(a["centroid_m"], axis=0), axis=1)))
    centroid = float(np.linalg.norm(b["centroid_m"][-1] - a["centroid_m"][-1]))
    cb = max(0.1 * radius, 0.02 * travel)
    assessed_duration = float(clock[-1] - clock[0])
    meaningful = len(clock) > 1 and assessed_duration > 0
    passed = (
        meaningful
        and rg_endpoint <= 0.02
        and rg_l2 <= 0.02
        and np.max(quantile_rel) <= 0.02
        and centroid <= cb
    )
    for key, data in [("reference", a), ("candidate", b)]:
        np.savez_compressed(out / (label + "_" + key + ".npz"), **data)
    result = dict(
        label=label,
        reference=str(a_path),
        candidate=str(b_path),
        reference_sha256=sha(a_path),
        candidate_sha256=sha(b_path),
        available_common_duration_s=end - t0,
        compared_duration_s=assessed_duration,
        comparison_status="EVALUATED_PREFIX" if meaningful else "NOT_EVALUABLE_NO_NONZERO_COMMON_CLOCK",
        sampling="actual last accepted state; no coordinate interpolation",
        reference_max_sample_age_s=float(a["sample_age_s"].max()),
        candidate_max_sample_age_s=float(b["sample_age_s"].max()),
        same_ID_particles=a["particle_ids"].tolist(),
        Rg_endpoint_relative=rg_endpoint if meaningful else None,
        Rg_series_relative_L2=rg_l2 if meaningful else None,
        pair_quartiles_relative=quantile_rel.tolist() if meaningful else None,
        centroid_difference_m=centroid if meaningful else None,
        centroid_budget_m=cb,
        rigid_fit_reference_end_m=float(a["rigid_fit_m"][-1]),
        rigid_fit_candidate_end_m=float(b["rigid_fit_m"][-1]),
        reference_neighbor_turnover_mean=float(np.mean(a["neighbor_turnover"]))
        if len(a["neighbor_turnover"])
        else None,
        candidate_neighbor_turnover_mean=float(np.mean(b["neighbor_turnover"]))
        if len(b["neighbor_turnover"])
        else None,
        pass_scientific=bool(passed),
        hard_contact_chatter_is_sole_failure_criterion=False,
    )
    write(out / (label + ".json"), result)
    return result


def gate(duration, qp_solver="sparse_active_set"):
    out = BASE / "comparisons"
    out.mkdir(parents=True, exist_ok=True)
    cases = {}
    for case in ["B2", "B3"]:
        selected = {}
        for p in sorted((BASE / "runs").glob("*/SUMMARY.json")):
            s = json.load(open(p))
            if (
                s["case"] == case
                and s["duration_s"] == duration
                and s.get("qp_solver", "sparse_active_set") == qp_solver
            ):
                selected[s["macro_dt_s"]] = (p.parent, s)
        if not all(dt in selected for dt in [20e-6, 10e-6]):
            cases[case] = dict(pass_gate=False, reason="RUNS_NOT_COMPLETE")
            continue
        nominal, ns = selected[20e-6]
        half, hs = selected[10e-6]
        comparison = compare(
            Path(ns["baseline_run"]) / "trajectories.h5",
            nominal / "trajectories.h5",
            out,
            f"{qp_solver}_{case}_{duration * 1000:g}ms_old_new",
        )
        sensitivity = compare(
            nominal / "trajectories.h5",
            half / "trajectories.h5",
            out,
            f"{qp_solver}_{case}_{duration * 1000:g}ms_dt_half",
        )
        numerical = all(
            s["status"] == "COMPLETE"
            and s["solver_failure"] is None
            and s["source_inputs_unchanged"]
            and s["kkt_max"] is not None
            and s["kkt_max"] <= 1e-9
            and s["linear_residual_max"] <= 1e-11
            and min(s["minimum_pair_gap_m"], s["minimum_wall_gap_m"]) >= -1e-11
            and s["multi_radius_corrections"] == 0
            for s in [ns, hs]
        )
        original = json.load(open(Path(ns["baseline_run"]) / "SUMMARY.json"))
        cases[case] = dict(
            pass_gate=bool(numerical and comparison["pass_scientific"] and sensitivity["pass_scientific"]),
            numerical_pass=numerical,
            scientific_comparison=comparison,
            dt_half_comparison=sensitivity,
            nominal=ns,
            half=hs,
            old_wall_s=original["wall_s"],
            wall_speedup=(original["wall_s"] / ns["wall_s"] if ns["status"] == "COMPLETE" else None),
            speedup_scope="Completed nominal runs only; incomplete attempts have no completion speedup",
        )
    result = dict(
        duration_s=duration,
        qp_solver=qp_solver,
        pass_gate=all(v["pass_gate"] for v in cases.values()),
        cases=cases,
    )
    prefix = "OSQP_" if qp_solver == "osqp" else ""
    write(out / f"{prefix}GATE_{duration * 1000:g}ms.json", result)
    return result


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--duration", type=float, choices=[0.001, 0.005], required=True)
    p.add_argument("--qp-solver", choices=["sparse_active_set", "osqp"], default="sparse_active_set")
    args = p.parse_args()
    print(json.dumps(gate(args.duration, args.qp_solver)), flush=True)
