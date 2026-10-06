"""Evidence report and eight figures; unavailable gated results are explicit."""

from pathlib import Path
from collections import Counter
import json
import xml.etree.ElementTree as ET
import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from scripts.stage5a.compare import sample
from scripts.stage4l.dense_probe import write, sha

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage5/motion_solver"
COLORS = ["#64748b", "#0078b5"]
NAMES = ["Frozen Stage4L", "Stage5A experimental"]


def rows(path):
    return [json.loads(line) for line in path.open()] if path.exists() else []


def main():
    final = BASE / "final"
    final.mkdir(parents=True, exist_ok=True)
    figdir = BASE / "figures"
    figdir.mkdir(exist_ok=True)
    old = json.load(open(BASE / "audit_v2/HISTORICAL_AUDIT.json"))
    profiles = json.load(open(BASE / "audit_v2/FIXED_STATE_PROFILES.json"))
    bench = json.load(open(BASE / "benchmarks/OSQP_BENCHMARK_SUMMARY.json"))
    new = {}
    for p in sorted((BASE / "runs").glob("*/SUMMARY.json")):
        s = json.load(open(p))
        if s.get("qp_solver") != "osqp":
            continue
        new[(s["case"], s["duration_s"], s["macro_dt_s"])] = (p.parent, s)
    gates = {}
    for p in (BASE / "comparisons").glob("OSQP_GATE_*ms.json"):
        gates[p.stem.removeprefix("OSQP_")] = json.load(open(p))
    numerical_failure = any(
        s["solver_failure"] is not None or s["status"] == "FAILED" for _, s in new.values()
    )
    scientific_failure = any(
        c.get("numerical_pass")
        and (
            not c["scientific_comparison"]["pass_scientific"]
            or not c["dt_half_comparison"]["pass_scientific"]
        )
        for g in gates.values()
        for c in g["cases"].values()
    )
    promoted = (
        "GATE_5ms" in gates
        and gates["GATE_5ms"]["pass_gate"]
        and all(c["wall_speedup"] >= 5 for c in gates["GATE_5ms"]["cases"].values())
    )
    status = (
        "MOTION_SOLVER_PROMOTED_WITH_SCOPE"
        if promoted
        else (
            "MOTION_SOLVER_NUMERICAL_GATE_FAILED"
            if numerical_failure or scientific_failure
            else "MOTION_SOLVER_RUNTIME_GATE_FAILED"
        )
    )
    gate = dict(
        status=status,
        completion_label="PROMOTED_WITH_SCOPE" if promoted else "FAILED_WITH_EVIDENCE",
        promotion=promoted,
        numerical_failure=numerical_failure,
        scientific_failure=scientific_failure,
        stages=gates,
        scope="Frozen Stage4L normal lubrication/A0 physics; translation numerics only",
        five_ms_not_launched_reason=None if "GATE_5ms" in gates else "Joint 1 ms gate not passed",
        new_long_transport_validated=False,
        numerical_method="STAGE5A_OSQP_PC_V2",
        nominal_runs={case: s for (case, duration, dt), (_, s) in new.items() if dt == 20e-6},
        stopped_previous_campaign=json.load(open(BASE / "STOPPED_BY_USER_20261006T203352Z.json")),
        scientific_scope="Only available prefixes can be compared; failed/incomplete runs do not validate the full target duration",
        primary_regression=json.load(open(BASE / "audit_v2/PRIMARY_TRANSPORT_REGRESSION.json")),
    )
    write(final / "MOTION_SOLVER_GATE.json", gate)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "pdf.fonttype": 42,
            "savefig.dpi": 320,
        }
    )
    figures = []

    def save(fig, name):
        for suffix in ["png", "pdf"]:
            fig.savefig(figdir / (name + "." + suffix), facecolor="white", bbox_inches="tight")
        plt.close(fig)
        figures.append(name)

    def pick(case):
        duration = 0.005 if (case, 0.005, 20e-6) in new else 0.001
        return duration, new[(case, duration, 20e-6)]

    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.4), layout="constrained")
    for ax, case in zip(axes, ["B2", "B3"]):
        for j, duration in enumerate([0.001, 0.005]):
            val = old[f"{case}_{duration * 1000:g}ms"]["summary"]["wall_s"]
            ax.bar(j - 0.18, val, 0.34, color=COLORS[0], label=NAMES[0] if j == 0 else None)
            item = new.get((case, duration, 20e-6))
            if item:
                s = item[1]
                ax.bar(
                    j + 0.18,
                    s["wall_s"],
                    0.34,
                    color=COLORS[1],
                    hatch="//" if s["status"] != "COMPLETE" else "",
                    label=NAMES[1] if j == 0 else None,
                )
                ax.annotate(
                    f"{s['wall_s']:.1f} s\n" + ("complete" if s["status"] == "COMPLETE" else "INCOMPLETE"),
                    (j + 0.18, s["wall_s"]),
                    xytext=(0, 4),
                    textcoords="offset points",
                    ha="center",
                    fontsize=7,
                )
            else:
                ax.text(
                    j + 0.18, val * 0.15, "NOT RUN", ha="center", rotation=90, color=COLORS[1], fontsize=8
                )
        ax.set(
            xticks=[0, 1],
            xticklabels=["1 ms target", "5 ms target"],
            ylabel="Total wall time (s)",
            title=case,
            xlim=(-0.55, 1.55),
        )
        ax.margins(y=0.25)
    axes[0].legend(frameon=False, fontsize=7)
    fig.suptitle("Incomplete runs do not define completion speedup", fontsize=10)
    save(fig, "FIG_M1_RUNTIME")
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.2), layout="constrained")
    native = {}
    for ax, case in zip(axes, ["B2", "B3"]):
        duration, (p, s) = pick(case)
        accepted = rows(p / "ACCEPTED.jsonl")
        native[case] = (p, s, accepted)
        historical = rows(Path(old[f"{case}_{duration * 1000:g}ms"]["run"]) / "steps.jsonl")
        if not accepted:
            ax.text(0.05, 0.38, "New: no accepted timestep", transform=ax.transAxes, color=COLORS[1])
        for data, color, label in [(historical, COLORS[0], NAMES[0]), (accepted, COLORS[1], NAMES[1])]:
            dt = np.sort([r["dt_s"] * 1e6 for r in data])
            if len(dt):
                q = np.quantile(dt, [0.05, 0.5, 0.95])
                ax.plot(
                    dt,
                    np.arange(1, len(dt) + 1) / len(dt),
                    color=color,
                    label=f"{label}\nmedian {q[1]:.3g} µs",
                )
                ax.scatter(q, [0.05, 0.5, 0.95], s=13, color=color)
        ax.set(
            xscale="log",
            xlabel="Accepted dt (µs)",
            ylabel="Empirical cumulative fraction",
            title=f"{case}: {duration * 1e3:g} ms target; new {s['status']}",
        )
        ax.legend(frameon=False, fontsize=7)
    save(fig, "FIG_M2_TIMESTEP")
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.6), layout="constrained")
    labels = ["Response", "Pair force", "Resistance", "Contact/QP", "Wall geometry", "Other"]
    palette = ["#8da0cb", "#66c2a5", "#e5c494", "#fc8d62", "#a6d854", "#d5d5d5"]
    breakdown = {}
    for ax, case in zip(axes, ["B2", "B3"]):
        pp = [r for r in profiles if r["case"] == case]
        den = sum(r["wall_s"] for r in pp)
        before = [
            sum(r["tracker_timings"]["response_queries_s"] for r in pp),
            sum(r["tracker_timings"]["pair_force_s"] for r in pp),
            sum(r["calls"]["mobility"]["seconds"] for r in pp),
            sum(r["calls"]["contact"]["seconds"] for r in pp),
            sum(
                sum(r["calls"][k]["seconds"] for k in ["wall_distance", "wall_gap", "wall_project"])
                for r in pp
            ),
        ]
        before.append(max(0, den - sum(before)))
        _, s, _ = native[case]
        d = s["timings"]
        tr = s["tracker_timings"]
        after = [
            tr["response_queries_s"],
            tr["pair_force_s"],
            d["assembly_s"],
            d["qp_qp_seconds"],
            d["wall_s"] + tr.get("wall_distance_s", 0),
        ]
        after.append(max(0, s["wall_s"] - sum(after)))
        breakdown[case] = dict(
            old_replay_seconds=before,
            new_run_seconds=after,
            old_replay_total=den,
            new_run_total=s["wall_s"],
            note="Old: 9 fixed-state replays; new: full attempted run, including error rejections. Missing early failed trial solve time stays in Other.",
        )
        stack = np.array([before, after])
        stack = stack / stack.sum(axis=1)[:, None] * 100
        base = np.zeros(2)
        for j, label in enumerate(labels):
            ax.bar([0, 1], stack[:, j], bottom=base, color=palette[j], label=label)
            base += stack[:, j]
        ax.set(
            xticks=[0, 1],
            xticklabels=["Old state replays", "New attempted run"],
            ylabel="Measured wall-time share (%)",
            title=case,
        )
    axes[1].legend(frameon=False, bbox_to_anchor=(1, 1), loc="upper left", fontsize=7)
    save(fig, "FIG_M3_SOLVER_BREAKDOWN")
    write(final / "BREAKDOWN_DATA.json", breakdown)
    datasets = {}
    for case in ["B2", "B3"]:
        duration, (p, s) = pick(case)
        baseline = Path(old[f"{case}_{duration * 1000:g}ms"]["run"])
        with h5py.File(baseline / "trajectories.h5") as h:
            t0 = float(h["time_s"][0])
            end = float(h["time_s"][-1])
        clock = np.linspace(t0, end, 501)
        a = sample(baseline / "trajectories.h5", clock)
        with h5py.File(p / "trajectories.h5") as h:
            last = float(h["time_s"][-1])
        # Short partial runs get a separate actual-time clock and explicit extent.
        bclock = (
            np.linspace(t0, last, min(501, max(2, s["accepted_steps"] + 1))) if last > t0 else np.array([t0])
        )
        b = sample(p / "trajectories.h5", bclock, a["particle_ids"])
        datasets[case] = (a, b, s)
        common = sample(baseline / "trajectories.h5", bclock, a["particle_ids"])
        np.savez_compressed(final / f"{case}_PLOT_COMMON_REFERENCE.npz", **common)
        np.savez_compressed(final / f"{case}_PLOT_REFERENCE.npz", **a)
        np.savez_compressed(final / f"{case}_PLOT_CANDIDATE.npz", **b)
    fig, axes = plt.subplots(1, 2, figsize=(7.5, 3.2), layout="constrained")
    for ax, case in zip(axes, ["B2", "B3"]):
        a, b, s = datasets[case]
        for data, color, label in [(a, COLORS[0], NAMES[0]), (b, COLORS[1], NAMES[1])]:
            ax.plot(
                (data["time_s"] - data["time_s"][0]) * 1e3,
                data["rg_m"] * 1e6,
                color=color,
                label=label,
                lw=1.3,
                marker="o" if len(data["time_s"]) == 1 else None,
            )
        ax.set(
            xlabel="Time from checkpoint (ms)",
            ylabel="Fixed-ID radius of gyration (µm)",
            title=f"{case}: {s['duration_s'] * 1e3:g} ms target",
        )
        if s["duration_s"] != 0.005:
            ax.text(0.03, 0.06, "5 ms new comparison gated off", transform=ax.transAxes, fontsize=8)
        if s["status"] != "COMPLETE":
            ax.text(
                0.03,
                0.15,
                f"New ends at {s['completed_duration_s'] * 1e6:.3g} µs",
                transform=ax.transAxes,
                fontsize=8,
            )
    axes[0].legend(frameon=False, fontsize=7)
    save(fig, "FIG_M4_CLUSTER_ENVELOPE")
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5.2), layout="constrained")
    for col, case in enumerate(["B2", "B3"]):
        a, b, s = datasets[case]
        if s["accepted_steps"]:
            a = dict(np.load(final / f"{case}_PLOT_COMMON_REFERENCE.npz"))
        for d, color, label in [(a, COLORS[0], NAMES[0]), (b, COLORS[1], NAMES[1])]:
            t = (d["time_s"] - d["time_s"][0]) * 1e6
            axes[0, col].plot(
                t, d["rigid_fit_m"] * 1e6, color=color, label=label, marker="o" if len(t) == 1 else None
            )
            axes[1, col].plot(t[1:], d["neighbor_turnover"], color=color)
        axes[0, col].set(title=f"{case}: shared partial clock", ylabel="Rigid-fit residual (µm)")
        axes[1, col].set(ylabel="Neighbor edge turnover", xlabel="Time from checkpoint (µs)")
        if s["accepted_steps"] == 0:
            for ax in axes[:, col]:
                ax.text(0.05, 0.3, "No accepted motion\nComparison unavailable", transform=ax.transAxes)
                ax.set_xlim(0, s["duration_s"] * 1e6)
            axes[0, col].set_title(f"{case}: frozen reference; new unavailable")
    axes[0, 0].legend(frameon=False, fontsize=7)
    save(fig, "FIG_M5_INTERNAL_REARRANGEMENT")
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 5), layout="constrained")
    for col, case in enumerate(["B2", "B3"]):
        p, s, accepted = native[case]
        oldsteps = rows(Path(s["baseline_run"]) / "steps.jsonl")
        with h5py.File(Path(s["baseline_run"]) / "trajectories.h5") as h:
            oldtimes = h["time_s"][:]
            extent = oldtimes[0] + (s["completed_duration_s"] or s["duration_s"])
            oldindices = np.flatnonzero(oldtimes <= extent)
            oldgaps = []
            from scipy.spatial.distance import pdist

            radius = float(h["radii_m"][0])
            for index in oldindices:
                points = h["positions_m"][index][h["active"][index]]
                oldgaps.append(float(np.min(pdist(points)) - 2 * radius))
        axes[0, col].plot(
            (oldtimes[oldindices] - oldtimes[0]) * 1e6,
            np.array(oldgaps) / radius,
            color=COLORS[0],
            label=NAMES[0],
        )
        axes[0, col].plot(
            [(d["time_s"] - 3) * 1e6 for d in accepted],
            [d["minimum_pair_gap_m"] / 1.78e-6 for d in accepted],
            color=COLORS[1],
            label=NAMES[1],
        )
        for data, color in [(oldsteps, COLORS[0]), (accepted, COLORS[1])]:
            axes[1, col].plot(
                [(d["time_s"] - 3) * 1e6 for d in data],
                [d.get("projection_max_m", 0) / 1.78e-6 for d in data],
                color=color,
            )
        axes[0, col].axhline(
            -1e-11 / 1.78e-6, color="#b91c1c", ls=":", lw=0.8, label="Existing tolerance bound"
        )
        axes[0, col].set(
            title=f"{case}: actual saved states",
            ylabel="Minimum all-particle gap / R",
            xlim=(0, (s["completed_duration_s"] or s["duration_s"]) * 1e6),
        )
        axes[1, col].set(
            xlabel="Time from checkpoint (µs)",
            ylabel="Maximum projection / R",
            xlim=(0, (s["completed_duration_s"] or s["duration_s"]) * 1e6),
        )
        if not accepted:
            for ax in axes[:, col]:
                ax.text(
                    0.05,
                    0.2,
                    "New: no accepted motion\nNo new gap/projection curve",
                    transform=ax.transAxes,
                    fontsize=8,
                )
    axes[0, 0].legend(frameon=False, fontsize=7)
    save(fig, "FIG_M6_GAP_SAFETY")
    fig, axes = plt.subplots(3, 2, figsize=(7.5, 7.3), layout="constrained")
    for col, case in enumerate(["B2", "B3"]):
        p, s, accepted = native[case]
        observed = accepted or [
            d["diagnostics"] for d in rows(p / "REJECTED.jsonl") if "corrector" in d["diagnostics"]
        ]
        histogram = Counter(n for d in observed for n in d["corrector"]["component_sizes"])
        axes[2, col].bar(list(histogram), list(histogram.values()), color=COLORS[1], width=1.5)
        axes[2, col].set(
            xlabel="Component size (particles)",
            ylabel="Count over " + ("accepted states" if accepted else "rejected trials"),
            yscale="linear",
        )
        t = (
            np.array([d["time_s"] - 3 for d in observed]) * 1e6
            if accepted
            else np.arange(1, len(observed) + 1)
        )
        for key, label in [
            ("component_count", "Count"),
            ("largest_component", "Largest size"),
            ("median_component", "Median size"),
        ]:
            axes[0, col].plot(t, [d["corrector"][key] for d in observed], label=label, lw=1)
        for key, label in [
            ("isolated_fraction", "Isolated particle fraction"),
            ("largest_component_fraction", "Largest component fraction"),
        ]:
            axes[1, col].plot(t, [d["corrector"][key] for d in observed], label=label, lw=1)
        axes[0, col].set(title=case, ylabel="Particles / components")
        axes[1, col].set(
            xlabel="Accepted physical time (µs)" if accepted else "Rejected trial index (not physical time)",
            ylabel="Fraction",
            ylim=(-0.02, 1.02),
        )
    axes[0, 0].legend(frameon=False, fontsize=7)
    axes[1, 0].legend(frameon=False, fontsize=7)
    save(fig, "FIG_M7_COMPONENT_STRUCTURE")
    fig, axes = plt.subplots(2, 2, figsize=(7.5, 6), layout="constrained")
    overlay = {}
    for row, case in enumerate(["B2", "B3"]):
        p, s, _ = native[case]
        a, b, _ = datasets[case]
        clock = np.array([3.0, 3.0 + s["completed_duration_s"]])
        oa = sample(Path(s["baseline_run"]) / "trajectories.h5", clock, a["particle_ids"])
        nb = sample(p / "trajectories.h5", clock, a["particle_ids"])
        overlay[case] = dict(
            clock_s=clock.tolist(),
            reference_actual_s=oa["actual_time_s"].tolist(),
            candidate_actual_s=nb["actual_time_s"].tolist(),
            reference_sample_age_s=oa["sample_age_s"].tolist(),
            candidate_sample_age_s=nb["sample_age_s"].tolist(),
            particle_ids=a["particle_ids"].tolist(),
            coordinate_transform="Subtract common initial reference centroid; SI to micrometres only; no magnification",
        )
        origin = oa["centroid_m"][0]
        for col, dims in enumerate([(0, 1), (0, 2)]):
            ax = axes[row, col]
            for d, color, label in [(oa, COLORS[0], NAMES[0]), (nb, COLORS[1], NAMES[1])]:
                xy = (d["positions_m"][-1] - origin)[:, dims] * 1e6
                ax.scatter(
                    xy[:, 0], xy[:, 1], s=10, facecolors="none", edgecolors=color, linewidths=0.8, label=label
                )
                hull = ConvexHull(xy)
                polygon = xy[np.r_[hull.vertices, hull.vertices[0]]]
                ax.plot(polygon[:, 0], polygon[:, 1], color=color, lw=0.8)
            ax.set_aspect("equal", adjustable="datalim")
            ax.set(
                xlabel="x − initial centroid (µm)",
                ylabel=("y" if col == 0 else "z") + " − initial centroid (µm)",
                title=(
                    f"{case}, clock +{s['completed_duration_s'] * 1e6:.3g} µs\n"
                    f"Sample age old/new: {oa['sample_age_s'][-1] * 1e6:.3g}/{nb['sample_age_s'][-1] * 1e6:.3g} µs"
                ),
            )
    axes[0, 0].legend(frameon=False, fontsize=7)
    fig.suptitle("Checkpoint t₀ = 3 s; no accepted new motion, initial state only", fontsize=10)
    save(fig, "FIG_M8_TRAJECTORY_OVERLAY")
    write(final / "OVERLAY_SAMPLE_AGES.json", overlay)
    write(
        final / "FIGURE_MANIFEST.json",
        dict(
            figures=[
                dict(path=str(figdir / (name + "." + ext)), sha256=sha(figdir / (name + "." + ext)))
                for name in figures
                for ext in ["png", "pdf"]
            ],
            status=status,
            coordinate_interpolation=False,
        ),
    )
    testroot = ET.parse(BASE / "environment/osqp_20261006/combined_tests.xml").getroot()
    counts = {
        key: sum(int(s.attrib.get(key, 0)) for s in testroot.iter("testsuite"))
        for key in ["tests", "failures", "errors", "skipped"]
    }
    tests = counts["tests"] - counts["failures"] - counts["errors"] - counts["skipped"]
    gate["test_counts"] = dict(
        counts,
        passed=tests,
        expected_failure_reason="Historical active-set B2 QP remains xfail; the same fixture passes with OSQP. New trajectory failures are separately recorded.",
    )
    write(final / "MOTION_SOLVER_GATE.json", gate)
    lines = [
        "# Stage5A motion solver development",
        "",
        f"**{gate['completion_label']} — {status}**",
        "",
        "Branch `dev/stage5a-motion-constrained-20261006`, base `40765a4f6252506b9ffdfb35453f62ca537c0797`. Preregistered plan commit `b24fbc0fd`; original active-set implementation `ce1808fd1`; OSQP scientific implementation `7dbc71f19f83687a8ff7e2646bc8402f1d5f634a`. All source/config identities are retained per run. The historical interrupted campaign is preserved separately and was not resumed.",
        "",
        "## Bottleneck and mathematical formulation",
        "",
        "Historical B2/B3 5 ms half-gap limiting fractions were 0.999425/0.999463. Contact consumed 1092.34/1825.85 s; sparse resistance assembly+solve consumed 804.50/899.01 s. Each old RK2 step uses two mobility/evaluator calls, four contact calls and two projections. Historical L-BFGS iteration counts were unavailable; nine independent saved-state replays per case measured them without overwriting old evidence.",
        "",
        "The experimental minimization is 0.5 Uᵀ A U − bᵀ U subject to G U ≥ lower, with unchanged Stage4L A=D+Lpp+Lpw and finite-horizon normal bounds. Pair rows contain at most six nonzeros; wall rows at most three. Vectorized normal outer-product assembly and sparse CSC KKT systems avoid full dense Jacobians/dual matrices. Following explicit user authorization, OSQP 1.1.3 was installed in an isolated dependency overlay. The original NumPy/SciPy/CuPy environment is unchanged, verified by before/after package inventories. OSQP uses builtin CPU sparse algebra; one sparse saddle polish and one cold retry can recover an inaccurate warm candidate only if all original inequality/KKT audits pass. The earlier matrix-free active-set implementation remains separately selectable and its failures remain preserved.",
        "",
        "KKT diagnostics use velocities divided by the maximum absolute target velocity and resistance divided by median Stokes drag. They audit scaled stationarity, inequality feasibility, nonnegative multipliers, active slack and complementarity. The separate resistance residual is relative to the physical force right-hand side. Summary KKT/linear maxima include completed solves in rejected local-error trials as well as accepted states; failed solves remain separate failure records.",
        "",
        "## Accepted-state reuse, components and integration",
        "",
        "Canonical particle/contact IDs carry primal velocity, dimensional dual multipliers, active contacts and sparsity identity only after acceptance. OSQP CSC workspaces are updated when the sparsity pattern matches; every attempted solve explicitly starts from accepted cache values or a cold zero state, not an unaccepted prior trial. Component decomposition includes both lubrication and horizon constraints; walls stay in their owning component. It is exact only for the current local resistance model. It does not generalize to full many-body mobility.",
        "",
        "One constrained predictor and one corrected solve refresh forces, lubrication coefficients and exact wall gaps/normals. Both start and predictor tangent inequalities are retained in the corrected displacement problem. Position/gap errors cause binary subdivision of that physical interval only; future nominal intervals reset. Final exact wall/pair checks and all-pair swept-segment checks guard acceptance. The OSQP attempts accepted no new state: zero projection calls are not successful motion-safety evidence. No force caps, damping, contact-tolerance relaxation, ablation or coordinate smoothing was introduced.",
        "",
        "## CPU/GPU and tests",
        "",
        f"{tests} focused tests passed, including T1–T15, 25 OSQP tests and legacy Stage4L/4M regression; {counts['skipped']} historical active-set real-QP regression remains xfail. That same saved QP passes with OSQP. New dynamic trajectory QPs nevertheless failed, so this test result does not promote the solver. Frozen legacy class-source hashes and force/ablation bytes are unchanged. Archived primary nominal/dt-half transport comparisons were recomputed with input SHA checks; that test does not validate new long transport.",
        "",
        f"Benchmark: `{json.dumps(bench['timings_median_s'])}`. CuPy selected for classical pair forces; CPU sparse QP. Serial components selected because the measured parallel overhead was greater. Identical-state warm-start speedup is {bench['warm_speedup']:.3f}x; this is not a claim about continuous-trajectory speedup. Global warm QP was faster than decomposed warm QP in these small cases. Component decomposition remains explicit and equivalent.",
        "",
        "## Dense probe outcomes",
        "",
        "| Case | Target (ms) | Macro dt (µs) | Status | Reached (µs) | Wall (s) | Accepted | Rejected | Median dt (µs) | Max KKT |",
        "|---|---:|---:|---|---:|---:|---:|---:|---:|---:|",
    ]
    for (case, duration, dt), (p, s) in sorted(new.items()):
        lines.append(
            f"| {case} | {duration * 1e3:g} | {dt * 1e6:g} | {s['status']} | {s['completed_duration_s'] * 1e6:.6g} | {s['wall_s']:.3f} | {s['accepted_steps']} | {s['rejected_intervals']} | {format(s['dt_quantiles_s'][2] * 1e6, '.6g') if s['dt_quantiles_s'] else 'N/A'} | {s['kkt_max']} (completed trial solves only) |"
        )
    lines += [
        "",
        "Complete frozen reference wall times: B2/B3 1 ms 400.37/491.86 s; 5 ms 2216.01/3080.30 s. Ratios involving incomplete new runs are not completion speedups. If the joint 1 ms gate fails, new 5 ms results and their speedups are unavailable, and no longer formation or transport run is launched.",
        "",
        "Machine-readable scientific differences (Rg endpoint/series, centroid, pair quartiles, rigid fit, neighbor turnover and sample ages) are in `results/stage5/motion_solver/comparisons/`. Partial-prefix comparisons cannot certify a complete 1 ms or 5 ms gate. Native hard-contact counts remain descriptive; no chatter filtering changes trajectories.",
        "",
        "## Runtime interpretation and limitations",
        "",
        "This is a coefficient-lagged predictor/corrector with implicit constrained velocities, not a fully implicit nonlinear position solve. The continuous dense predictor/corrector problem is harder than a repeated frozen-state QP. Both OSQP runs rejected ten local-error intervals and failed before accepting a step. Nearly parallel start/predictor tangent constraints challenge both ADMM and active-set polishing. An excellent repeated-state warm-start benchmark is not trajectory acceleration. Accepted dt, new-motion Rg/centroid/quantile differences, and completed-run speedups are unavailable; initial-state identity is not a convergence result.",
        "",
        "Tolerance values were frozen before scientific runs and were not loosened after poor performance. Runtime exhaustion is incomplete evidence, not solver convergence or completed science. The 1,800 s budget is checked after an accepted step, so one in-progress bounded solver call can cause an overrun; it is not a hard real-time deadline. Classical near-contact secondary-force validity, normal-only immobile-interface lubrication, absent tangential/history/deformation physics and long-transport accuracy remain outside this numerical development.",
        "",
        "## Delivered evidence",
        "",
        "Eight publication-size PNG/PDF figures are under `results/stage5/motion_solver/figures/`. Gated-off 5 ms comparisons are explicitly labelled; available shorter/partial data are not relabelled as 5 ms. M3 compares measured original fixed-state replay timing with measured new attempted-run timing; it does not invent unrecorded historical stage timing. M8 uses actual last accepted coordinates at the same clock and records sample age, with a common-origin translation and metre-to-micrometre conversion only.",
        "",
        "Authoritative gate: `results/stage5/motion_solver/final/MOTION_SOLVER_GATE.json`. Full per-run source/config/input identities, accepted/rejected ledgers, actual trajectories, last accepted state, resource observations and failures remain in their exclusive run directories. New JSON is strict and finite. No original Stage4L/M result was overwritten.",
    ]
    attempt_path = final / "ATTEMPT_ANALYSIS.json"
    if attempt_path.exists():
        lines += [
            "",
            "## Attempt cost decomposition",
            "",
            "Completed QP calls below include accepted and local-error-rejected trials. Time inside the terminal failed QP is retained in the failure/campaign evidence and not invented as a completed solve timing.",
            "",
            "| Method / case / macro dt (µs) | Rejection reasons | Predictor QP (s) | Corrector QP (s) | Corrector polish (s) | Corrector primal fallback QPs |",
            "|---|---|---:|---:|---:|---:|",
        ]
        for attempt in json.loads(attempt_path.read_text()):
            corrector = attempt["solvers"]["corrector"]
            lines += [
                f"| {'OSQP' if attempt['run'].endswith('_osqp') else 'Earlier active-set'} / {attempt['case']} / {attempt['macro_dt_s'] * 1e6:g} | {attempt['rejection_reasons']} | {attempt['solvers']['predictor']['time_s']:.3f} | {corrector['time_s']:.3f} | {corrector['polish_s']:.3f} | {corrector['primal_fallback_qps']} |"
            ]
    lines += [
        "",
        "## Scientific comparisons on the available prefix",
        "",
        "The preregistered comparison clock is 10 µs. A prefix with no nonzero common clock is not evaluable. Passing a short prefix does not pass the requested 1 ms gate.",
        "",
        "| Case / comparison | Compared time (µs) | Rg endpoint relative | Rg series L2 | Pair Q1/median/Q3 relative | Centroid difference (m) | Prefix assessment |",
        "|---|---:|---:|---:|---|---:|---|",
    ]
    for stage in gates.values():
        for case, casegate in stage["cases"].items():
            for kind in ["scientific_comparison", "dt_half_comparison"]:
                comparison = casegate.get(kind)
                if not comparison:
                    continue
                lines += [
                    f"| {case} / {kind} | {comparison['compared_duration_s'] * 1e6:.6g} | {comparison['Rg_endpoint_relative']} | {comparison['Rg_series_relative_L2']} | {comparison['pair_quartiles_relative']} | {comparison['centroid_difference_m']} | {comparison['comparison_status']} |"
                ]
    lines += ["", "## Observed numerical failure and blockers", ""]
    for path in sorted((BASE / "diagnostics").glob("*/DIAGNOSTIC.json")):
        diagnostic = json.load(open(path))
        audit = diagnostic.get("matrix_audit", {})
        lines += [
            f"Bounded {diagnostic.get('case', 'B2')} single-interval reconstruction: **{diagnostic['status']}**. "
            + (
                "The earlier active-set cache is reconstructed from one preceding accepted interval. "
                if audit
                else "The OSQP replay starts at the unchanged initial state with an empty accepted cache. "
            )
            + f"No new trajectory is continued. Diagnostic: `{path.relative_to(ROOT)}`."
        ]
        small = diagnostic.get("small_component_audit")
        if small:
            lines += [
                f"OSQP diagnostic constraint matrix dimensions: {small['dimensions']}; numerical rank {small['numerical_rank']}; {small['rows_with_near_parallel_partner']} rows have a nearly parallel partner. Dense rank/cosine checks are diagnostic-only; all production QP matrices remain sparse."
            ]
        if audit:
            lines += [
                f"Failed component: {audit['primal_dimension']} primal unknowns, {audit['constraint_count']} inequalities; active saddle condition {audit['active_saddle_condition']:.6g}. {audit['rows_with_parallel_partner_cos_above_1_minus_1e_8']} rows have a nearly parallel partner. The feasible fallback exhausted {audit['primal_fallback_iterations']} iterations. Its final working multiplier minimum is {audit.get('actual_last_working_min_multiplier')}; this does not satisfy dual feasibility.",
                "The original failure-summary stationarity value was recomputed with zero multipliers after fallback exhaustion; it is not an accepted-state residual. The saved matrix audit additionally reports actual last-working-set stationarity, negative multipliers and active slack. Dense SVD was used only on this diagnostic component, never in production motion.",
            ]
    lines += [
        "",
        "Blocking issues: robust treatment of nearly dependent start/predictor constraints, continuous-state QP convergence and cost, and stiff integration-error control. Passing isolated or frozen-state tests does not resolve these. No tolerance, activation/floor, secondary force, acoustic field or CFD value was altered to make the gate pass.",
        "",
        "## Accepted-state safety and resources",
        "",
    ]
    for (case, duration, dt), (_, result) in sorted(new.items()):
        lines += [
            f"{case}, target {duration * 1e3:g} ms, macro {dt * 1e6:g} µs: accepted pair/wall gaps are unavailable (zero accepted steps); projection calls {result['projection_calls']}. Maximum linear relative residual across completed trial solves: {result['linear_residual_max']}. The unchanged geometric tolerance is 1e-11 m. Peak RAM {result['peak_ram_bytes'] / 2**30:.3f} GiB, GPU allocation {result['peak_gpu_memory_MiB']:.1f} MiB, mean process CPU {result['cpu_utilization_mean_percent']:.1f}% (100%=one core), sampled GPU utilization {result['gpu_utilization_mean_percent']:.2f}%."
        ]
    lines += [
        "",
        "The earlier FORMAT_AUDIT.json applies to the original active-set revision only. OSQP run identities record source commit 7dbc71f19 and full file hashes. Later analysis scripts are separately hashed and do not rewrite executed provenance. The frozen legacy class/force hashes continue to pass.",
    ]
    lines += [
        "",
        "## OSQP installation and dynamic failure details",
        "",
        "Installation lock: requirements-stage5a-osqp.lock.txt. Invocation: `bash scripts/stage5a/python_osqp.sh ...`. Installation, package hashes, smoke tests, failed diagnostic attempts and combined JUnit report are retained in `results/stage5/motion_solver/environment/osqp_20261006/`. No scientific tolerance was changed after the OSQP trajectory attempts.",
        "",
    ]
    for (case, duration, dt), (_, result) in sorted(new.items()):
        failure = result.get("solver_failure") or {}
        attempts = failure.get("attempts", [])
        if attempts:
            last = attempts[-1]
            lines += [
                f"{case}: failed QP has {failure['matrix_dimension']} unknowns and {failure['constraint_count']} constraints; OSQP returned `{last['status']}` after {last['iterations']} iterations, with independent candidate KKT {last['kkt']:.6g}. Sparse polishing did not produce a valid full KKT certificate. This is a rejected trial, not an accepted-state residual."
            ]
    lines += [
        "",
        "B3's ADMM candidate KKT is below 1e-9, but its `solved inaccurate` status and unsuccessful sparse polish do not satisfy the preregistered adapter acceptance rule. B2 also exceeds the independent KKT bound. Neither failed candidate was used to advance physical time. Installing OSQP resolves the dependency gap and a saved QP failure; it does not by itself resolve dynamic near-dependent constraint handling.",
    ]
    (ROOT / "docs/STAGE5A_MOTION_SOLVER_DEVELOPMENT.md").write_text("\n".join(lines) + "\n")
    print(status, flush=True)


if __name__ == "__main__":
    main()
