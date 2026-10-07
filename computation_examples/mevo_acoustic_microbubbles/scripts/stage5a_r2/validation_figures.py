"""R1/R2/R5/R6/R8 from recorded trajectories, ledgers and process observations."""

import json
import textwrap
from pathlib import Path
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts.stage5a_r2.seal import ROOT, BASE
from scripts.stage4l.dense_probe import sha, write

BLUE = "#315f89"
GRAY = "#777a80"
RED = "#b65a53"
PURPLE = "#795390"
plt.rcParams.update(
    {
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    }
)
SOURCES = {}
FIGURE_CONTEXT = {}


def local(path):
    path = str(path)
    marker = "/workspace/mevo_acoustic_microbubbles_stage5a_r2/"
    return ROOT / path.removeprefix(marker) if path.startswith(marker) else Path(path)


def read(path, lines=False):
    p = local(path)
    SOURCES[str(p.relative_to(ROOT))] = sha(p)
    return [json.loads(s) for s in p.read_text().splitlines() if s] if lines else json.loads(p.read_text())


def save(fig, name, description):
    fig.suptitle(name.replace("_", " "), fontweight="bold")
    commit = (ROOT / "SOURCE_COMMIT_R2.txt").read_text().strip()
    kernel = json.loads((BASE / "final/VALIDATION_PROTOCOL_V3.json").read_text())["source_commit"]
    fig.text(
        0.01,
        -0.055,
        "\n".join(
            textwrap.wrap(
                f"R2 S1/S0 V6 | kernel freeze {kernel[:12]} | "
                + "; ".join(f"{k}: {v}" for k, v in FIGURE_CONTEXT.items()),
                width=170,
            )
        ),
        fontsize=6.5,
    )
    for suffix in ["png", "pdf"]:
        fig.savefig(BASE / "figures" / f"{name}.{suffix}", dpi=240, bbox_inches="tight")
    plt.close(fig)
    write(
        BASE / "figures" / f"{name}_SOURCES.json",
        dict(
            definition=description,
            files=SOURCES.copy(),
            script_sha256=sha(Path(__file__)),
            orchestration_checkout_at_render_time=commit,
            numerical_kernel_commit=kernel,
            selected_runs=FIGURE_CONTEXT.copy(),
        ),
    )


def select():
    campaign = read(sorted((BASE / "final").glob("CAMPAIGN_V*/RESULT.json"))[-1])
    entries = [e for e in campaign["ledger"] if e.get("summary")]
    # Longest attempted nominal, retaining failures. Performance comparisons
    # separately filter COMPLETE and equal durations; no partial extrapolation.
    chosen = {}
    for case in ["B2", "B3"]:
        items = [e for e in entries if e["case"] == case and e["kind"] == "nominal" and e["repeat"] == 0]
        chosen[case] = max(items, key=lambda e: e["duration"])
        FIGURE_CONTEXT[case] = Path(chosen[case]["run"]).name
    return campaign, entries, chosen


def progress(entries, chosen):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.3), layout="constrained")
    for ax, (case, e) in zip(axes, chosen.items()):
        stopped = []
        for kind, color in [("nominal", BLUE), ("refined", PURPLE), ("reference", GRAY)]:
            row = next(
                (
                    r
                    for r in entries
                    if r["case"] == case
                    and r["duration"] == e["duration"]
                    and r["kind"] == kind
                    and r["repeat"] == 0
                ),
                None,
            )
            if not row:
                continue
            summary = row["summary"]
            if summary["status"] != "COMPLETE":
                failure = summary.get("failure")
                reason = failure.get("reason") if isinstance(failure, dict) else failure
                reason = reason or summary.get("error") or "See terminal log"
                stopped.append(f"{kind}: {str(reason)[:85]}")
            folder = local(row["run"])
            a = read(folder / "ACCEPTED.jsonl", True)
            ax.plot(
                [0] + [r["wall_s"] for r in a],
                [0] + [r["elapsed_s"] * 1e3 for r in a],
                color=color,
                label=f"{kind}: {row['summary']['status']}",
            )
            if kind != "reference":
                rejected = read(folder / "REJECTED.jsonl", True)
                ax.scatter(
                    [r["wall_s"] for r in rejected],
                    [r["elapsed_s"] * 1e3 for r in rejected],
                    s=5,
                    color=color,
                    marker="x",
                    alpha=0.5,
                )
                recovery = [r for r in a if r.get("substep_method") == "S0"]
                ax.scatter(
                    [r["wall_s"] for r in recovery],
                    [r["elapsed_s"] * 1e3 for r in recovery],
                    s=8,
                    color="#b77e38",
                    label="S0 recovery" if kind == "nominal" else None,
                )
                fb = [r for r in a if r["fallbacks"]]
                ax.scatter(
                    [r["wall_s"] for r in fb], [r["elapsed_s"] * 1e3 for r in fb], s=25, color=RED, marker="^"
                )
        ax.axhline(e["duration"] * 1e3, ls="--", lw=1, color="black", label="Requested duration")
        ax.set(xlabel="Elapsed wall time (s)", ylabel="Accepted physical duration (ms)", title=case)
        ax.legend(fontsize=8)
        if stopped:
            ax.text(
                0.02,
                0.02,
                "\n".join(stopped),
                transform=ax.transAxes,
                fontsize=7,
                va="bottom",
                color=RED,
                bbox=dict(facecolor="white", alpha=0.9, edgecolor="none"),
            )
    save(
        fig,
        "R1_PROGRESS",
        "Actual accepted duration vs wall time. Crosses: rejected intervals; red triangles: fallbacks in accepted substeps. Endpoint status is explicit; no extrapolation.",
    )


def error_and_safety(chosen):
    fig, axes = plt.subplots(3, 2, figsize=(12, 9), layout="constrained")
    safe, saxes = plt.subplots(3, 2, figsize=(12, 9), layout="constrained")
    rejfig, raxes = plt.subplots(2, 2, figsize=(12, 6.5), layout="constrained")
    timing_audit = {}
    for column, (case, e) in enumerate(chosen.items()):
        folder = local(e["run"])
        a = read(folder / "ACCEPTED.jsonl", True)
        rejected = read(folder / "REJECTED.jsonl", True)
        t = np.array([r["elapsed_s"] for r in a]) * 1e3
        axes[0, column].plot(t, [r["dt_s"] * 1e6 for r in a], color=BLUE, lw=0.7, label="Accepted half-step")
        axes[0, column].scatter(
            [r["elapsed_s"] * 1e3 for r in rejected],
            [r["dt_s"] * 1e6 for r in rejected],
            color=RED,
            s=8,
            marker="x",
            label="Rejected full trial",
        )
        axes[0, column].set(yscale="log", ylabel="Actual dt (µs)", title=f"{case} | {e['summary']['status']}")
        axes[0, column].legend(fontsize=8)
        for key, color in zip(["position", "relative", "pair", "wall"], [BLUE, GRAY, PURPLE, RED]):
            axes[1, column].plot(t, [r["error_terms"][key] for r in a], lw=0.6, label=key, color=color)
        axes[1, column].axhline(1, color="black", ls="--", lw=1)
        axes[1, column].set(ylabel="Accepted trial error / budget", yscale="log")
        axes[1, column].legend(fontsize=8, ncol=2)
        reasons = sorted({r["reason"] for r in rejected})
        for key, color in zip(["position", "relative", "pair", "wall"], [BLUE, GRAY, PURPLE, RED]):
            rr = [
                (j, r["diagnostics"]["error_terms"][key])
                for j, r in enumerate(rejected)
                if "error_terms" in r["diagnostics"]
            ]
            raxes[0, column].scatter([j for j, v in rr], [v for j, v in rr], s=8, color=color, label=key)
        raxes[0, column].axhline(1, color="black", ls="--", lw=1)
        raxes[0, column].set(title=case, yscale="log", ylabel="Completed rejected error / budget")
        raxes[0, column].legend(fontsize=7, ncol=2)
        for key, label, color in [
            ("original_equilibrium_kkt", "Completed trial: original equilibrium", BLUE),
            ("terminal_original_nonlinear_kkt", "Failed nonlinear final candidate", RED),
        ]:
            rr = [(j, r["diagnostics"][key]) for j, r in enumerate(rejected) if key in r["diagnostics"]]
            raxes[1, column].scatter([j for j, v in rr], [v for j, v in rr], s=10, color=color, label=label)
        raxes[1, column].axhline(1e-9, color="black", ls="--", lw=1)
        raxes[1, column].set(yscale="log", ylabel="Original-equation KKT")
        raxes[1, column].legend(fontsize=7)
        for ax in raxes[:, column]:
            ax.set_xlabel("Rejected-trial index (not physical time)")
        for i, reason in enumerate(reasons):
            rr = [(j, r) for j, r in enumerate(rejected) if r["reason"] == reason]
            axes[2, column].scatter([j for j, _ in rr], [r["dt_s"] * 1e6 for _, r in rr], s=12, label=reason)
        axes[2, column].set(
            xlabel="Rejected-trial index (not physical time)", ylabel="Rejected trial dt (µs)", yscale="log"
        )
        if reasons:
            axes[2, column].legend(fontsize=7)
        for ax in axes[:2, column]:
            ax.set_xlabel("Accepted physical time since 3 s (ms)")
        for key, label, color in [
            ("minimum_pair_swept_gap_m", "Pair: continuous minimum", BLUE),
            ("minimum_wall_swept_gap_bound_m", "WALL: certified lower bound", RED),
            ("minimum_wall_endpoint_gap_m", "Wall endpoint", GRAY),
        ]:
            saxes[0, column].plot(t, [r["safety"][key] * 1e9 for r in a], lw=0.7, label=label, color=color)
        saxes[0, column].axhline(-0.01, ls="--", color="black", label="Original −0.01 nm tolerance")
        saxes[0, column].set(ylabel="Minimum signed gap / bound (nm)", title=case)
        saxes[0, column].set_yscale("symlog", linthresh=0.01)
        saxes[0, column].legend(fontsize=7)
        for key, label, color in [
            ("kkt", "Numerical QP KKT (S0: start matrix)", BLUE),
            ("nonlinear_kkt", "Original endpoint KKT (S1 only)", PURPLE),
            ("linear_residual", "Resistance residual", GRAY),
        ]:
            saxes[1, column].plot(t, [r[key] for r in a], lw=0.6, color=color, label=label)
        saxes[1, column].axhline(1e-9, color="black", ls="--")
        saxes[1, column].axhline(1e-11, color=GRAY, ls=":")
        saxes[1, column].set(ylabel="Accepted-state residual", yscale="log")
        saxes[1, column].legend(fontsize=8)
        saxes[2, column].plot(
            t, [r["safety"]["projection_max_m"] * 1e9 for r in a], color=BLUE, label="Residual projection"
        )
        saxes[2, column].plot(
            t,
            [r["dt_geometry_error_bound_m"] * 1e9 for r in a],
            color=GRAY,
            label="dt × constraint velocity error",
        )
        saxes[2, column].set(ylabel="Geometric displacement (nm)", ylim=(-0.001, 0.011))
        saxes[2, column].legend(fontsize=8)
        for ax in saxes[:, column]:
            ax.set_xlabel("Accepted physical time since 3 s (ms)")
        s = e["summary"]
        q = s["timings_exclusive_s"]
        partition = {key: q[key] for key in ["forcing", "assembly", "qp"]}
        partition["I/O"] = s["io_seconds"]
        partition["Geometry, control and other overhead"] = e["wall_s_with_process_start"] - sum(
            partition.values()
        )
        assert partition["Geometry, control and other overhead"] >= 0
        timing_audit[case] = dict(
            run=e["run"],
            exclusive_partition_s=partition,
            correction="Raw timings_exclusive_s.error includes three calls also counted in wall. Do not stack error and wall. These overlapping raw counters are excluded from this partition; their union and monitoring/control/init are in the measured remainder.",
            monitoring_inclusive_s=s["monitor_inclusive_seconds"],
            process_start_and_outer_wrapper_s=e["wall_s_with_process_start"] - s["wall_s"],
        )
    save(
        fig,
        "R2_TIMESTEP_ERROR",
        "Actual dt and step-doubling error components. Separate axes for accepted physical time and rejected-trial index; rejected full trial vs accepted half-step are explicitly distinguished.",
    )
    save(
        safe,
        "R5_GEOMETRY_SAFETY",
        "Certified minima along every adopted straight half-step; true endpoint gap and conservative sweep bound are distinguished. Symmetric logarithmic gap axis is linear within the unchanged tolerance band. Complete values are in ACCEPTED.jsonl.",
    )
    save(
        rejfig,
        "R2b_REJECTED_TRIALS",
        "Completed rejected LOCAL_ERROR trials retain four error components; failed nonlinear candidates have separate final original-equation residuals. Missing uncompleted error estimates are not plotted as zero.",
    )
    write(BASE / "figures/TIMING_AUDIT.json", timing_audit)


def metrics():
    original_context = FIGURE_CONTEXT.copy()
    for p in sorted((BASE / "comparisons").glob("*_legacy.json")):
        if "refined_" in p.stem:
            continue
        report = read(p)
        FIGURE_CONTEXT.clear()
        FIGURE_CONTEXT.update({k: Path(report[k]).parent.name for k in ["reference", "candidate"]})
        label = p.stem
        files = [p.with_name(label + "_" + key + ".npz") for key in ["reference", "candidate"]]
        arrays = [dict(np.load(f)) for f in files]
        for f in files:
            SOURCES[str(f.relative_to(ROOT))] = sha(f)
        fig, axes = plt.subplots(2, 3, figsize=(13, 7.2), layout="constrained")
        origin = arrays[0]["centroid_m"][0]
        for d, name, color in zip(arrays, ["Frozen Stage4L", "R2 S1/S0 V6"], [GRAY, BLUE]):
            t = (d["requested_time_s"] - d["requested_time_s"][0]) * 1e3
            axes[0, 0].plot(t, d["rg_m"] * 1e6, label=name, color=color)
            axes[0, 1].plot(
                t, np.linalg.norm(d["centroid_m"] - origin, axis=1) * 1e6, label=name, color=color
            )
            for j, ls in zip([1, 2, 3], [":", "-", "--"]):
                axes[0, 2].plot(
                    t, d["pair_quantiles_m"][:, j] * 1e6, ls=ls, color=color, label=f"{name}: Q{j}"
                )
            axes[1, 0].plot(t, d["rigid_fit_m"] * 1e6, label=name, color=color)
            axes[1, 1].plot(t[1:], d["neighbor_turnover"], label=name, color=color)
            axes[1, 2].plot(t, d["sample_age_s"] * 1e6, label=name, color=color)
        for ax, y in zip(
            axes.flat,
            [
                "Radius of gyration (µm)",
                "Centroid travel from common origin (µm)",
                "Pair-distance quartiles (µm)",
                "Rigid-fit RMS residual (µm)",
                "Neighbor graph turnover",
                "Sample age (µs)",
            ],
        ):
            ax.set(xlabel="Physical time since 3 s (ms)", ylabel=y)
            ax.legend(fontsize=7)
        title_prefix = (
            "Partial-prefix metrics: "
            if report.get("eligible_for_target_gate") is False
            else "Registered metric comparison: "
        )
        axes[0, 0].set_title(title_prefix + report["status"], fontsize=10)
        if report.get("eligible_for_target_gate") is False:
            FIGURE_CONTEXT["scope"] = "Partial overlap; complete 5 ms target NOT MET"
        max_age = max(float(np.max(abs(d["sample_age_s"]))) for d in arrays)
        if max_age < 1e-15:
            axes[1, 2].set_ylim(-0.01, 0.01)
            axes[1, 2].text(
                0.03,
                0.95,
                f"Max |age| = {max_age:.2g} s (roundoff)",
                transform=axes[1, 2].transAxes,
                va="top",
                fontsize=8,
            )
        save(
            fig,
            f"R6_MOTION_METRICS_{label}",
            f"Fixed original IDs {report['same_ID_particles']}. Actual sampled positions, no coordinate interpolation. Comparison status: {report['status']}; compared duration {report['compared_duration_s']} s.",
        )

    FIGURE_CONTEXT.clear()
    FIGURE_CONTEXT.update(original_context)


def performance(entries, chosen):
    fig, axes = plt.subplots(2, 3, figsize=(14, 8.2), layout="constrained")
    paired = [
        duration
        for duration in {e["duration"] for e in entries}
        if all(
            any(
                e["duration"] == duration
                and e["case"] == case
                and e["kind"] == kind
                and e["summary"]["status"] == "COMPLETE"
                for e in entries
            )
            for case in ["B2", "B3"]
            for kind in ["nominal", "reference"]
        )
    ]
    duration = max(paired, default=max(e["duration"] for e in entries))
    for j, case in enumerate(["B2", "B3"]):
        for k, kind in enumerate(["reference", "nominal"]):
            samples = [
                e["wall_s_with_process_start"]
                for e in entries
                if e["case"] == case
                and e["duration"] == duration
                and e["kind"] == kind
                and e["summary"]["status"] == "COMPLETE"
            ]
            x = j * 3 + k
            if samples:
                axes[0, 0].scatter(
                    np.full(len(samples), x), samples, color=GRAY if kind == "reference" else BLUE, s=30
                )
                axes[0, 0].plot([x - 0.25, x + 0.25], [np.median(samples)] * 2, color="black")
                axes[0, 0].text(x, max(samples) * 1.02, f"n={len(samples)}", ha="center", fontsize=8)
            else:
                axes[0, 0].text(
                    x,
                    0.02,
                    "Incomplete",
                    rotation=90,
                    ha="center",
                    transform=axes[0, 0].get_xaxis_transform(),
                    fontsize=8,
                )
    axes[0, 0].set_xticks([0, 1, 3, 4], ["B2 old", "B2 R2", "B3 old", "B3 R2"])
    axes[0, 0].set(ylabel="Complete-run wall time (s)", title=f"Equal-duration task: {duration * 1e3:g} ms")
    audit = read(BASE / "figures/TIMING_AUDIT.json")
    bottom = np.zeros(2)
    for i, key in enumerate(next(iter(audit.values()))["exclusive_partition_s"]):
        values = [audit[c]["exclusive_partition_s"][key] for c in ["B2", "B3"]]
        axes[0, 1].bar(range(2), values, bottom=bottom, label=key)
        bottom += values
    axes[0, 1].set_xticks(range(2), ["B2 R2", "B3 R2"])
    axes[0, 1].set_ylabel("Mutually exclusive elapsed cost (s)")
    axes[0, 1].set_title(
        "Selected run costs: " + ", ".join(f"{c} {e['duration'] * 1e3:g} ms" for c, e in chosen.items()),
        fontsize=10,
    )
    axes[0, 1].legend(fontsize=7)
    for case, color in [("B2", BLUE), ("B3", PURPLE)]:
        e = chosen[case]
        folder = local(e["run"])
        resources = read(folder / "RESOURCE_SAMPLES.json")
        a = read(folder / "ACCEPTED.jsonl", True)
        t = np.array([r["monotonic_s"] for r in resources])
        t -= t[0]
        axes[0, 2].plot(
            t, [r["process_cpu_percent"] for r in resources], label=f"{case} process", color=color
        )
        axes[0, 2].plot(
            t,
            [r["host_cpu_percent"] for r in resources],
            ls=":",
            alpha=0.7,
            label=f"{case} host mean",
            color=color,
        )
        axes[1, 0].plot(t, [r.get("gpu_percent", np.nan) for r in resources], label=case, color=color)
        axes[1, 1].plot(
            t, [r["rss_bytes"] / 2**30 for r in resources], label=f"{case} process RAM", color=color
        )
        axes[1, 1].plot(
            t,
            [r.get("gpu_memory_MiB", np.nan) / 1024 for r in resources],
            ls="--",
            label=f"{case} device VRAM",
            color=color,
        )
        axes[1, 2].plot(
            [r["elapsed_s"] * 1e3 for r in a],
            [max(r["component_sizes"]) for r in a],
            color=color,
            label=f"{case}: largest graph component",
        )
    for ax, y in [
        (axes[0, 2], "CPU usage (%)"),
        (axes[1, 0], "Device utilization (%)"),
        (axes[1, 1], "Memory (GiB)"),
    ]:
        ax.set(xlabel="Elapsed wall time (s)", ylabel=y)
        ax.legend(fontsize=7)
    axes[1, 2].set(xlabel="Accepted physical time since 3 s (ms)", ylabel="Particles in connected component")
    axes[1, 2].legend(fontsize=7)
    save(
        fig,
        "R8_PERFORMANCE_RESOURCES",
        "Only complete equal-duration runs enter performance samples; all retries and setup included. One scientific worker, BLAS=1. Actual final solve layout is component serial; component graph is diagnostic. Host/GPU readings are device-wide. Cost partition excludes overlapping wall/error subcounters; see TIMING_AUDIT.json.",
    )
    # Full per-core evidence is plotted separately, not replaced by host mean.
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout="constrained")
    for ax, (case, e) in zip(axes, chosen.items()):
        r = read(local(e["run"]) / "RESOURCE_SAMPLES.json")
        array = np.array([q["per_core_percent"] for q in r]).T
        im = ax.imshow(
            array,
            origin="lower",
            aspect="auto",
            vmin=0,
            vmax=100,
            cmap="viridis",
            extent=[0, r[-1]["monotonic_s"] - r[0]["monotonic_s"], -0.5, len(array) - 0.5],
        )
        ax.set(
            title=f"{case}: 1 scientific worker",
            xlabel="Elapsed wall time (s)",
            ylabel="Host logical CPU index",
        )
    fig.colorbar(im, ax=list(axes), label="Per-core CPU usage (%)")
    save(
        fig,
        "R8b_PER_CORE_RESOURCES",
        "Userspace psutil per-core observations, device-wide sampling. Affinity/quota recorded separately; not a claim every core belongs to this task.",
    )


if __name__ == "__main__":
    campaign, entries, chosen = select()
    progress(entries, chosen)
    error_and_safety(chosen)
    metrics()
    performance(entries, chosen)
