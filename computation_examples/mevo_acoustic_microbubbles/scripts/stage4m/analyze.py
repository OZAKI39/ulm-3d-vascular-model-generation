"""Read-only scientific comparisons and engineering gates for Stage 4M."""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.spatial.distance import pdist

from scripts.stage4l.analyze_dense import edges, largest
from scripts.stage4l.native_graph import graph_statistics, hard_threshold_motion
from scripts.stage4m.common import BASE, refs, write, sha


def analyze_dense(run, out):
    """Never write into the supplied run (including existing A0)."""
    out.mkdir(parents=True, exist_ok=True)
    summary = json.loads((run / "SUMMARY.json").read_text())
    with h5py.File(run / "trajectories.h5", "r") as h:
        t, ids, r = h["time_s"][:], h["particle_id"][:], h["radii_m"][:]
        active = h["active"][0]
        member = np.flatnonzero(active)[largest(h["positions_m"][0][active], 2.1 * r[0])]
        x = h["positions_m"][:][:, member]
        all_finite = bool(np.isfinite(h["positions_m"][:]).all())
    if len(t) < 2:
        result = dict(status=summary["status"], summary=summary, analysis_status="NO_ACCEPTED_STEPS")
        write(out / "DENSE_METRICS.json", result)
        return result
    centered = x - x.mean(axis=1, keepdims=True)
    native_rg = np.sqrt(np.mean(np.sum(centered**2, axis=-1), axis=1))
    # Inherited common 10 us diagnostic clock, using actual accepted states only.
    samples = t[0] + np.arange(int(np.floor((t[-1] - t[0]) / 1e-5 + 1e-8)) + 1) * 1e-5
    if abs(samples[-1] - t[-1]) < 1e-12:
        samples[-1] = t[-1]
    indices = np.maximum(np.searchsorted(t, samples, side="right") - 1, 0)
    selected = centered[indices]
    reference = centered[0]
    rigid, quantiles, component, graph = [], [], [], []
    for frame in selected:
        u, _, vt = np.linalg.svd(reference.T @ frame)
        q = u @ np.diag([1.0, 1.0, np.linalg.det(u @ vt)]) @ vt
        rigid.append(float(np.sqrt(np.mean(np.sum((reference @ q - frame) ** 2, axis=1)))))
        quantiles.append(np.quantile(pdist(frame), [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]))
        component.append(len(largest(frame, 2.1 * r[0])))
        graph.append(edges(frame, 2.1 * r[0]))
    turnover = [len(a ^ b) / max(1, len(a | b)) for a, b in zip(graph[:-1], graph[1:])]
    rms = np.sqrt(np.mean(np.sum((selected - reference) ** 2, axis=-1), axis=1))
    np.savez_compressed(
        out / "CLUSTER_SERIES.npz",
        time_s=samples,
        actual_time_s=t[indices],
        same_ID_particle_ids=ids[member],
        rigid_fit_residual_m=rigid,
        centered_rms_m=rms,
        radius_of_gyration_m=native_rg[indices],
        largest_component_size=component,
        turnover_fraction=turnover,
        unordered_pair_distance_quantiles_m=quantiles,
    )
    native = dict(
        near_contact=graph_statistics(x, t, 2.1 * r[0]),
        hard_contact=graph_statistics(x, t, 2 * r[0] + 1e-11),
        hard_threshold_motion=hard_threshold_motion(x, r[0]),
    )
    translation = json.loads((run / "TRANSLATION_DIAGNOSTICS.json").read_text())
    result = dict(
        status=summary["status"],
        summary=summary,
        translation=translation,
        same_ID_particle_ids=ids[member].tolist(),
        initial_cluster_size=len(member),
        initial_rg_m=float(native_rg[0]),
        final_rg_m=float(native_rg[-1]),
        max_rg_ratio=float(np.max(native_rg) / native_rg[0]),
        centered_rms_end_m=float(rms[-1]),
        rigid_fit_residual_end_m=rigid[-1],
        neighbor_turnover_mean=float(np.mean(turnover)) if turnover else None,
        final_largest_same_ID_component=component[-1],
        unordered_pair_distance_quantiles_initial_m=quantiles[0].tolist(),
        unordered_pair_distance_quantiles_final_m=quantiles[-1].tolist(),
        native_graph=native,
        all_positions_finite=all_finite,
        all_steps_positive=bool(np.all(np.diff(t) > 0)),
        max_common_clock_sample_age_s=float(np.max(samples - t[indices])),
        trajectory_sha256=sha(run / "trajectories.h5"),
        semantics="All-accepted-state graph rates and censored lifetimes; no position interpolation. Common-clock traces use preceding accepted states; sampling age is explicit. Dense initial clusters were formed by L0, not A2.",
    )
    write(out / "DENSE_METRICS.json", result)
    return result


def dense_gate(rows):
    decisions = {}
    for case in ("B2", "B3"):
        a0, a2 = rows[case]["A0"], rows[case]["A2"]
        old, new = a0["summary"], a2["summary"]
        step_ratio = old["accepted_steps"] / new["accepted_steps"] if new["accepted_steps"] else None
        wall_ratio = old["wall_s"] / new["wall_s"]
        complete = new["status"] == "COMPLETE" and abs(new["completed_duration_s"] - 0.005) < 1e-10
        stable = (
            complete
            and a2.get("all_positions_finite", False)
            and a2.get("all_steps_positive", False)
            and new["error"] is None
            and new["projection"]["multi_radius_corrections"] == 0
            and a2.get("max_rg_ratio", float("inf")) <= 2
            and not a2.get("translation", {}).get("failure")
            and new["source_inputs_unchanged"]
        )
        improvement = bool((step_ratio is not None and step_ratio >= 5) or wall_ratio >= 5)
        decisions[case] = dict(
            pass_gate=bool(stable and improvement),
            completed=complete,
            numerical_and_envelope_checks=bool(stable),
            accepted_step_reduction=step_ratio,
            wall_speedup=wall_ratio,
            improvement_gate=improvement,
            same_ID_max_rg_ratio=a2.get("max_rg_ratio"),
        )
    passed = all(row["pass_gate"] for row in decisions.values())
    result = dict(
        status="PASS_ENGINEERING_CONTINUATION_GATE" if passed else "A2_RUNTIME_GATE_FAILED",
        continue_formation=passed,
        cases=decisions,
        scope="Fivefold is an engineering continuation criterion, not physical accuracy. Historical A0 predates final contact QR fallbacks; wall ratios are not isolated identical-source benchmarks.",
    )
    write(BASE / "comparisons/M0_GATE.json", result)
    return result


def dense_figures(rows, locations):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    colors = dict(A0="#6b4c9a", A1="#969696", A2="#237b9a")
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    out = BASE / "comparisons"

    def save(fig, name):
        for suffix in ("png", "pdf"):
            fig.savefig(out / f"{name}.{suffix}", dpi=300, bbox_inches="tight")
        plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)
    for ax, metric, label in zip(axs, ["wall_s", "accepted_steps"], ["Wall time (s)", "Accepted steps"]):
        for k, policy in enumerate(["A0", "A1", "A2"]):
            bars = ax.bar(
                np.arange(2) + (k - 1) * 0.25,
                [rows[c][policy]["summary"][metric] for c in ["B2", "B3"]],
                0.23,
                color=colors[policy],
                label=policy,
            )
            for bar, case in zip(bars, ["B2", "B3"]):
                if rows[case][policy]["status"] != "COMPLETE":
                    bar.set_hatch("//")
                    ax.text(
                        bar.get_x() + bar.get_width() / 2,
                        bar.get_height(),
                        "Incomplete",
                        rotation=90,
                        va="bottom",
                        ha="center",
                        fontsize=8,
                    )
        ax.set_xticks([0, 1], ["B2", "B3"])
        ax.set_ylabel(label)
        ax.margins(y=0.3)
    axs[0].legend(frameon=False, ncol=3, loc="upper right")
    fig.suptitle("Dense 5 ms diagnostic: runtime and accepted steps")
    save(fig, "Figure_M2_dense_performance")
    fig, axs = plt.subplots(2, 3, figsize=(10, 6), constrained_layout=True)
    for j, case in enumerate(["B2", "B3"]):
        for policy in ["A0", "A2"]:
            file = locations[case][policy] / "CLUSTER_SERIES.npz"
            if not file.exists():
                continue
            with np.load(file) as d:
                elapsed = (d["time_s"] - d["time_s"][0]) * 1e3
                for col, key, factor, ts in [
                    (0, "rigid_fit_residual_m", 1e6, elapsed),
                    (1, "turnover_fraction", 100, elapsed[1:]),
                    (2, "radius_of_gyration_m", 1e6, elapsed),
                ]:
                    axs[j, col].plot(ts, d[key] * factor, color=colors[policy], label=policy, lw=1.1)
        for col, ylabel in enumerate(
            ["Rigid-fit residual (µm)", "Neighbor turnover (%)", "Same-ID gyration radius (µm)"]
        ):
            axs[j, col].set(xlabel="Probe time (ms)", ylabel=ylabel, title=case)
    axs[0, 0].legend(frameon=False)
    save(fig, "Figure_M3_dense_internal_dynamics")
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)
    for ax, case in zip(axs, ["B2", "B3"]):
        for policy in ["A0", "A1", "A2"]:
            file = (
                locations[case][policy] / "PAIR_FORCE_DIAGNOSTICS.jsonl"
                if policy == "A0"
                else Path(rows[case][policy]["summary"]["absolute_run"]) / "PAIR_FORCE_DIAGNOSTICS.jsonl"
            )
            diag = [json.loads(line) for line in file.read_text().splitlines()]
            ax.plot(
                [(d["time_s"] - 3) * 1e3 for d in diag],
                [d["effective_classical_ratio"] for d in diag],
                label=policy,
                color=colors[policy],
                lw=1,
            )
        ax.set(
            xlabel="Probe time (ms)",
            ylabel="Effective / classical pair-norm sum",
            title=case,
            ylim=(-0.02, 1.05),
        )
    axs[0].legend(frameon=False)
    save(fig, "Figure_M5_dense_secondary_removal")


def analyze_all(legacy, stage4l):
    rows, locations = {}, {}
    for case in ("B2", "B3"):
        rows[case], locations[case] = {}, {}
        runs = {"A0": refs(legacy, stage4l, case)["dense_A0"]}
        for policy in ("A1", "A2"):
            matches = sorted((BASE / "dense_probe" / policy).glob(f"*_{case}_{policy}_5ms/SUMMARY.json"))
            if len(matches) != 1:
                raise ValueError(f"Expected one finished {policy}/{case} probe, got {matches}")
            runs[policy] = matches[0].parent
        for policy, run in runs.items():
            out = BASE / "comparisons/dense" / case / policy
            result = analyze_dense(run, out)
            result["summary"]["absolute_run"] = str(run)
            rows[case][policy] = result
            locations[case][policy] = out
    gate = dense_gate(rows)
    write(BASE / "comparisons/DENSE_COMPARISON.json", rows)
    dense_figures(rows, locations)
    print(json.dumps(gate), flush=True)
    return rows, gate


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--legacy-root", type=Path, required=True)
    p.add_argument("--stage4l-root", type=Path, required=True)
    a = p.parse_args()
    analyze_all(a.legacy_root, a.stage4l_root)
