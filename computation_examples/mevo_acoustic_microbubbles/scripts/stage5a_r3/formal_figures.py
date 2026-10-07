"""G4/G6/G7/G8: actual formal data, explicit incomplete cases, no synthetic motion."""

import numpy as np
import matplotlib.pyplot as plt
from scripts.stage5a_r3.common import BASE
from scripts.stage5a_r3.paths import local
from scripts.stage5a_r3.current import campaign as current_campaign
from scripts.stage5a_r3.figures import save, source, history, BLUE, RED, PURPLE, GRAY, GOLD


def formal():
    return current_campaign()["ledger"]


def nominal(case):
    return next(
        r
        for r in formal()
        if r["case"] == case
        and r["kind"] == "nominal"
        and r["duration_s"] == 0.005
        and r["repeat"] == 0
        and r.get("run")
    )


def g4():
    rows = source(BASE / "final/FORMAL_COMPARISONS.json")["rows"]
    fig = plt.figure(figsize=(15, 8), layout="constrained")
    arrays = {}
    for row, case in enumerate(["B2", "B3"]):
        r = next(
            x for x in rows if x["case"] == case and x["name"].endswith("_0.005_refined") and x.get("samples")
        )
        a = np.load(local(r["samples"]))
        errors = a["all_active_raw_errors_m"] * 1e9
        t = (a["clock_s"] - a["clock_s"][0]) * 1e3
        arrays[case + "_errors_nm"] = errors
        arrays[case + "_particle_id"] = a["particle_id"]
        arrays[case + "_time_ms"] = t
        ax = fig.add_subplot(2, 3, 3 * row + 1)
        for mode, color, label in [
            ("raw", BLUE, "Raw"),
            ("translation_removed", RED, "Translation removed"),
            ("rigid_fit", PURPLE, "Rigid fit"),
        ]:
            v = a[f"all_active_{mode}_errors_m"][-1] * 1e9
            ax.hist(v, bins=24, histtype="step", lw=1.5, color=color, label=label)
        ax.axvline(89, color=GRAY, ls="--", label="0.05 R target")
        ax.set(
            xlabel="Endpoint same-ID error (nm)", ylabel="Bubble count", title=f"{case}: nominal vs refined"
        )
        ax.legend(fontsize=8)
        ax = fig.add_subplot(2, 3, 3 * row + 2)
        art = ax.pcolormesh(t, np.arange(errors.shape[1]), errors.T, shading="auto", cmap="magma")
        ax.set(
            xlabel="Time after checkpoint (ms)",
            ylabel="Fixed active-ID index",
            title=f"All {errors.shape[1]} active particles",
        )
        fig.colorbar(art, ax=ax, label="Raw error (nm)", shrink=0.85)
        ax = fig.add_subplot(2, 3, 3 * row + 3, projection="3d")
        j, i = np.unravel_index(errors.argmax(), errors.shape)
        nb = np.load(local(r["worst_neighborhood"]))
        main = int(np.flatnonzero(nb["particle_id"] == nb["worst_particle_id"])[0])
        origin = nb["left_positions_m"][0, main]
        trace_keys = ["left_all_accepted_positions_m", "right_all_accepted_positions_m"]
        union = np.concatenate([nb[k].reshape(-1, 3) for k in trace_keys])
        extent = (union - origin) * 1e6
        center = (extent.min(0) + extent.max(0)) / 2
        half = max(np.ptp(extent, axis=0)) * 0.55
        arrays[case + "_neighbor_particle_id"] = nb["particle_id"]
        arrays[case + "_neighbor_origin_m"] = origin
        for key, color, label in [
            (trace_keys[0], RED, "Refined"),
            (trace_keys[1], BLUE, "Nominal"),
        ]:
            x = (nb[key] - origin) * 1e6
            arrays[case + "_neighbor_" + key] = nb[key]
            for k in range(x.shape[1]):
                ax.plot(
                    *x[:, k].T,
                    color=color,
                    alpha=1 if k == main else 0.24,
                    lw=1.5 if k == main else 0.5,
                    ls="--" if label == "Refined" else "-",
                    label=label if k == main else None,
                )
        ax.set(
            xlabel="X offset (µm)",
            ylabel="Y offset (µm)",
            zlabel="Z offset (µm)",
            title=f"Worst ID {int(a['particle_id'][i])} and real neighbors",
        )
        ax.legend(fontsize=8)
        ax.set_box_aspect((1, 1, 1))
        ax.set_xlim(center[0] - half, center[0] + half)
        ax.set_ylim(center[1] - half, center[1] + half)
        ax.set_zlim(center[2] - half, center[2] + half)
    save(
        fig,
        "G4_all_active_formal_precision",
        "All active IDs preserved, no rematching. Error statistics use actual stored frames at common clock, with sampling certification retained. Worst-neighbor paths use every stored accepted state at original coordinates, equal physical axis scales, and no coordinate interpolation. Independent global error bound is NOT_EVALUABLE.",
        arrays,
    )


def g6():
    fig, axs = plt.subplots(2, 3, figsize=(14.5, 8), layout="constrained")
    arrays = {}
    for row, case in enumerate(["B2", "B3"]):
        entry = nominal(case)
        folder = local(entry["run"])
        d = history(folder / "trajectories.h5")
        accepted = source(folder / "ACCEPTED.jsonl")
        active = d["active"][0]
        x = d["positions_m"][:, active]
        dt = (d["time_s"] - d["time_s"][0]) * 1e3
        displacement = np.linalg.norm(x - x[0], axis=2) * 1e6
        axs[row, 0].plot(dt, np.sqrt(np.mean(displacement**2, axis=1)), color=BLUE, label="RMS")
        axs[row, 0].plot(dt, np.quantile(displacement, 0.95, axis=1), color=PURPLE, label="p95")
        axs[row, 0].plot(dt, displacement.max(1), color=GRAY, lw=1, label="Maximum")
        axs[row, 0].set(
            xlabel="Time after checkpoint (ms)",
            ylabel="Actual center displacement (µm)",
            title=case + " | " + entry["status"],
        )
        axs[row, 0].legend()
        ts = np.array([r["elapsed_s"] for r in accepted]) * 1e3
        pair = np.array([r["safety"]["minimum_pair_swept_gap_m"] for r in accepted])
        wall = np.array([r["safety"]["minimum_wall_swept_gap_bound_m"] for r in accepted])
        for values, color, label in [
            (pair, BLUE, "Pair certified lower bound"),
            (wall, RED, "WALL certified lower bound"),
        ]:
            axs[row, 1].plot(ts, values * 1e12, color=color, lw=0.8, label=label)
        axs[row, 1].axhline(-10, color=GRAY, ls="--", label="Unchanged tolerance")
        axs[row, 1].set(
            yscale="symlog",
            xlabel="Time after checkpoint (ms)",
            ylabel="Continuous gap / bound (pm)",
            title="Every adopted half step",
        )
        axs[row, 1].legend(fontsize=8)
        for key, color, label in [
            ("kkt", BLUE, "QP KKT"),
            ("nonlinear_kkt", PURPLE, "S1 endpoint KKT"),
            ("linear_residual", GOLD, "Resistance residual"),
        ]:
            values = np.array([r[key] if r[key] is not None else np.nan for r in accepted])
            arrays[case + "_accepted_" + key] = values
            axs[row, 2].semilogy(ts, np.maximum(values, 1e-18), color=color, lw=0.7, label=label)
        axs[row, 2].axhline(1e-9, color=GRAY, ls="--")
        axs[row, 2].axhline(1e-11, color=GRAY, ls=":")
        axs[row, 2].set(
            xlabel="Time after checkpoint (ms)",
            ylabel="Dimensionless residual",
            title="Accepted states only | zero projection",
        )
        axs[row, 2].legend(fontsize=8)
        arrays.update(
            {
                case + "_actual_time_s": d["time_s"],
                case + "_displacement_um": displacement,
                case + "_pair_swept_gap_m": pair,
                case + "_wall_swept_lower_bound_m": wall,
            }
        )
    save(
        fig,
        "G6_full_motion_continuous_safety",
        "Real all-active displacements and every accepted continuous certificate. Negative tolerance-band gaps retained. S0 endpoint nonlinear KKT is absent, not zero. Rejected candidates are excluded from accepted residuals.",
        arrays,
    )


def g4b():
    rows = source(BASE / "final/FORMAL_COMPARISONS.json")["rows"]
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    arrays = {}
    modes = [
        ("raw", BLUE, "Raw"),
        ("translation_removed", RED, "Translation removed"),
        ("rigid_fit", PURPLE, "Rigid fit"),
    ]
    names = ["all_active", "original_group", "complement"]
    for col, case in enumerate(["B2", "B3"]):
        r = next(
            x for x in rows if x["case"] == case and x["name"].endswith("_0.005_refined") and x.get("samples")
        )
        x = np.arange(3)
        labels = [
            label + f"\nn = {r['groups'][name]['count']}"
            for name, label in zip(names, ["All active", "Original group", "Complement"])
        ]
        for row, metric, target in [(0, "rms", 17.8), (1, "max", 89.0)]:
            ax = axs[row, col]
            for k, (mode, color, label) in enumerate(modes):
                values = (
                    np.array([r["groups"][name]["errors_m"][mode]["endpoint"][metric] for name in names])
                    * 1e9
                )
                ax.bar(x + (k - 1) * 0.24, values, width=0.22, color=color, label=label)
                arrays[f"{case}_{metric}_{mode}_nm"] = values
            ax.axhline(target, color=GRAY, ls="--", label="Frozen all-active target")
            ax.set(
                xticks=x,
                xticklabels=labels,
                yscale="log",
                ylabel=f"Endpoint {metric.upper()} difference (nm)",
                title=f"{case}: nominal / refined | +5 ms",
            )
            if row == 0:
                ax.legend(fontsize=8)
        arrays[case + "_groups"] = np.array(names)
    save(
        fig,
        "G4b_original_group_and_complement",
        "Every active particle is included exactly once in the original group or its complement. Raw/translation-removed/rigid-fit endpoint errors are descriptive fixed-ID comparisons. The frozen all-active resolution targets are shown as guides; small original-group or centroid errors do not establish all-active accuracy. No independent global error bound is available.",
        arrays,
    )


def g5c():
    rows = source(BASE / "final/FORMAL_COMPARISONS.json")["rows"]
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True, layout="constrained")
    data = {}
    for ax, case in zip(axs, ["B2", "B3"]):
        chosen, labels = [], []
        for suffix, label in [
            ("refined", "R3 / refined"),
            ("reference", "R3 / fresh"),
            ("archived", "R3 / archived"),
            ("archived_vs_fresh", "Archived / fresh"),
            ("R2", "R3 / R2"),
        ]:
            r = next(
                x
                for x in rows
                if x["case"] == case and x["name"].endswith("_0.005_" + suffix) and x.get("samples")
            )
            chosen.append(r)
            labels.append(label + f"\n{r['duration_s'] * 1000:.3f} ms common")
        x = np.arange(len(chosen))
        for n, mode, color, label in [
            (0, "raw", BLUE, "Raw"),
            (1, "translation_removed", RED, "Translation removed"),
            (2, "rigid_fit", PURPLE, "Rigid fit"),
        ]:
            values = (
                np.array([r["groups"]["all_active"]["errors_m"][mode]["endpoint"]["rms"] for r in chosen])
                * 1e9
            )
            ax.bar(x + (n - 1) * 0.25, values, width=0.23, color=color, label=label)
            data[case + "_" + mode + "_endpoint_RMS_nm"] = values
        ax.axhline(17.8, color=GRAY, ls="--", label="0.01 R resolution target")
        ax.set(
            xticks=x,
            xticklabels=labels,
            yscale="log",
            ylabel="Endpoint all-active RMS difference (nm)",
            title=case + " | named numerical lineages",
        )
        ax.tick_params(axis="x", labelsize=8)
        data[case + "_names"] = np.array([r["name"] for r in chosen])
        data[case + "_duration_s"] = np.array([r["duration_s"] for r in chosen])
        data[case + "_sampling_bound_certified"] = np.array([r["sampling_bound_certified"] for r in chosen])
    handles, labels = axs[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside upper center", ncol=4, fontsize=9)
    save(
        fig,
        "G5c_full_reference_lineage",
        "Original IDs and all active bubbles. Each bar is a descriptive difference against its named numerical lineage, not a physical accuracy estimate. Shared actual duration is labeled; the incomplete historical B2 R2 path is not extrapolated. Sampling certification and separate raw/translation/rigid errors are retained in the source reports. BASELINE_DISAGREEMENT and NOT_EVALUABLE global accuracy remain.",
        data,
    )


def g7():
    windows = ["initial", "stable", "S0", "event"]
    rows = [source(BASE / "diagnosis" / f"DYNAMIC_INSTRUMENTED_{w}.json") for w in windows]
    fig, axs = plt.subplots(2, 2, figsize=(13, 8), layout="constrained")
    categories = [
        "qp_native_candidate",
        "qp_adapter_and_fallback",
        "qp_partition_and_global_certificate",
        "resistance_assembly",
        "force_response_flow",
        "wall_query",
        "finite_feature_constraints",
        "continuous_geometry",
        "nonlinear_and_substep_other",
        "step_doubling_error_and_other",
        "checkpoint_IO",
        "driver_initialization_IO_monitoring_other",
    ]
    colors = plt.cm.tab20(np.linspace(0, 1, len(categories)))
    bottom = np.zeros(4)
    cost = np.array([[r["exclusive_s"].get(k, 0) for r in rows] for k in categories])
    for k, c, v in zip(categories, colors, cost):
        axs[0, 0].bar(windows, v, bottom=bottom, color=c, label=k.replace("_", " "))
        bottom += v
    axs[0, 0].set(ylabel="Exclusive instrumented wall time (s)", title="Mutually exclusive dynamic cost")
    axs[0, 0].legend(fontsize=6, ncol=2)
    x = np.arange(4)
    stages = ["full", "half1", "half2"]
    stage_counts = []
    for n, stage in enumerate(stages):
        accepted = []
        rejected = []
        for r in rows:
            accepted.append(
                sum(a["count"] for a in r["stage_totals"] if a["stage"] == stage and a["trial_accepted"])
            )
            rejected.append(
                sum(a["count"] for a in r["stage_totals"] if a["stage"] == stage and not a["trial_accepted"])
            )
        axs[0, 1].bar(x + (n - 1) * 0.24, accepted, width=0.22, color=[BLUE, PURPLE, GOLD][n], label=stage)
        axs[0, 1].bar(
            x + (n - 1) * 0.24,
            rejected,
            bottom=accepted,
            width=0.22,
            color=[BLUE, PURPLE, GOLD][n],
            hatch="///",
            alpha=0.5,
        )
        stage_counts.append([accepted, rejected])
    axs[0, 1].set(
        xticks=x,
        xticklabels=windows,
        ylabel="Executed advance calls",
        title="Step doubling | hatched: rejected trial",
    )
    axs[0, 1].legend()
    qps = np.array([len(r["qp_calls"]) for r in rows])
    candidates = np.array([len(r["native_candidates"]) for r in rows])
    axs[1, 0].bar(x - 0.18, qps, 0.35, color=BLUE, label="QP adapter calls")
    axs[1, 0].bar(x + 0.18, candidates, 0.35, color=PURPLE, label="Native candidates (incl. fallback)")
    axs[1, 0].set(xticks=x, xticklabels=windows, ylabel="Call count", title="Changing-state QP work")
    axs[1, 0].legend(fontsize=8)
    bottom = np.zeros(4)
    native_phases = []
    for key, color in [("setup_s", GOLD), ("update_s", PURPLE), ("solve_s", BLUE), ("polish_s", GRAY)]:
        values = np.array(
            [sum((q["diagnostics"] or {}).get(key, 0) for q in r["native_candidates"]) for r in rows]
        )
        axs[1, 1].bar(windows, values, bottom=bottom, color=color, label=key)
        bottom += values
        native_phases.append(values)
    axs[1, 1].set(ylabel="Native recorded time (s)", title="Native solver phases | separate timing view")
    axs[1, 1].legend(fontsize=8)
    save(
        fig,
        "G7_dynamic_QP_and_exclusive_cost",
        "All four windows use real changing full-system states. Full/half1/half2 and rejected trials counted separately. Diagnostic instrumentation overhead included; these times are never used for end-to-end speedup. Native phase view overlaps the exclusive partition.",
        dict(
            exclusive_cost_s=cost,
            windows=np.array(windows),
            categories=np.array(categories),
            qp_calls=qps,
            native_candidates=candidates,
            stages=np.array(stages),
            stage_accept_reject_counts=np.array(stage_counts),
            native_phase_names=np.array(["setup_s", "update_s", "solve_s", "polish_s"]),
            native_phase_time_s=np.array(native_phases),
        ),
    )


def g7b():
    windows = ["B3_global_costly", "B3_serial_costly"]
    rows = [source(BASE / "diagnosis" / f"DYNAMIC_INSTRUMENTED_{w}.json") for w in windows]
    labels = ["Global QP", "Component serial"]
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.8), layout="constrained")
    categories = sorted(set().union(*(r["exclusive_s"] for r in rows)))
    cost = np.array([[r["exclusive_s"].get(k, 0) for r in rows] for k in categories])
    bottom = np.zeros(2)
    for key, color, values in zip(categories, plt.cm.tab20(np.linspace(0, 1, len(categories))), cost):
        axs[0].bar(labels, values, bottom=bottom, color=color, label=key.replace("_", " "))
        bottom += values
    axs[0].set(ylabel="Instrumented wall time (s)", title="B3 +1.5 ms | exclusive cost")
    axs[0].legend(fontsize=5.5, ncol=1, loc="upper left")
    stage_counts, stage_cost, nonlinear = [], [], []
    for stage in ["full", "half1", "half2"]:
        stage_counts.append([sum(s["count"] for s in r["stage_totals"] if s["stage"] == stage) for r in rows])
        stage_cost.append([sum(s["wall_s"] for s in r["stage_totals"] if s["stage"] == stage) for r in rows])
        nonlinear.append(
            [sum(s["nonlinear_iterations"] for s in r["stage_totals"] if s["stage"] == stage) for r in rows]
        )
    bottom = np.zeros(2)
    for stage, color, values in zip(["full", "half1", "half2"], [BLUE, PURPLE, GOLD], stage_cost):
        axs[1].bar(labels, values, bottom=bottom, color=color, label=stage)
        bottom += values
    axs[1].set(ylabel="Advance-call inclusive wall time (s)", title="Step-doubling work | separate view")
    axs[1].legend(fontsize=8)
    qps = np.array([len(r["qp_calls"]) for r in rows])
    native = np.array([len(r["native_candidates"]) for r in rows])
    axs[2].bar(np.arange(2) - 0.18, qps, 0.35, color=BLUE, label="QP adapter")
    axs[2].bar(np.arange(2) + 0.18, native, 0.35, color=PURPLE, label="Native candidate")
    axs[2].set(xticks=[0, 1], xticklabels=labels, ylabel="Actual call count", title="Partitioned QP calls")
    axs[2].legend(fontsize=8)
    save(
        fig,
        "G7b_B3_global_and_component_dynamic_cost",
        "Same actual all-active initial state and requested +15 microsecond interval, with explicit cold cache. Adaptive work and component call counts can differ. Instrumented timings include diagnostic overhead and are not formal completion speedups. Exclusive cost, inclusive stage times and native counts are separate views, not additive totals.",
        dict(
            windows=np.array(windows),
            categories=np.array(categories),
            exclusive_cost_s=cost,
            stages=np.array(["full", "half1", "half2"]),
            stage_counts=np.array(stage_counts),
            stage_wall_s=np.array(stage_cost),
            nonlinear_iterations=np.array(nonlinear),
            qp_calls=qps,
            native_candidates=native,
            actual_elapsed_s=np.array([r["summary"]["segment_completed_duration_s"] for r in rows]),
        ),
    )


def g8():
    gate = source(BASE / "final/MOTION_SOLVER_R3_GATE.json")
    resource = source(BASE / "final/RESOURCE_LAYOUT_BENCHMARK.json")
    fig, axs = plt.subplots(2, 2, figsize=(13, 8), layout="constrained")
    arrays = {}
    for i, row in enumerate(gate["performance"]):
        new = np.array(row["new_wall_s"])
        old = np.array(row["reference_wall_s"])
        for values, offset, color, label in [(new, -0.14, BLUE, "R3"), (old, 0.14, RED, "Frozen Stage4L")]:
            spread = np.linspace(-0.05, 0.05, len(values)) if len(values) > 1 else np.zeros(len(values))
            axs[0, 0].scatter(
                i + offset + spread,
                values / 60,
                color=color,
                label=label if i == 0 else None,
                s=45,
            )
            if len(values):
                axs[0, 0].plot(
                    [i + offset - 0.08, i + offset + 0.08], [np.median(values) / 60] * 2, color=color, lw=2
                )
        ratios = [r["speedup"] for r in row["paired"] if r["speedup"] is not None]
        spread = np.linspace(-0.045, 0.045, len(ratios)) if len(ratios) > 1 else np.zeros(len(ratios))
        axs[0, 1].scatter(i + spread, ratios, color=BLUE, s=45)
        axs[0, 1].text(
            i, 0.04, row["status"], ha="center", transform=axs[0, 1].get_xaxis_transform(), fontsize=9
        )
        arrays[row["case"] + "_new_wall_s"] = new
        arrays[row["case"] + "_reference_wall_s"] = old
    axs[0, 0].set(
        xticks=[0, 1],
        xticklabels=["B2", "B3"],
        ylabel="Complete +5 ms process wall time (min)",
        title="Individual full-task repeats",
    )
    axs[0, 0].legend()
    axs[0, 1].axhline(5, color=GRAY, ls="--", label="5× target")
    axs[0, 1].set(
        xticks=[0, 1],
        xticklabels=["B2", "B3"],
        ylabel="Complete reference / R3 wall time",
        title="No incomplete-reference extrapolation",
    )
    axs[0, 1].legend()
    axs[0, 1].set_xlim(-0.5, 1.5)
    axs[0, 1].set_ylim(0, max(5.5, axs[0, 1].get_ylim()[1]))
    names = sorted(set((r["case"], r["source_time_s"]) for r in resource["full_trial_layouts"]))
    x = np.arange(len(names))
    for n, (layout, color) in enumerate(
        [("global", BLUE), ("serial", PURPLE), ("parallel2", GOLD), ("parallel4", GRAY)]
    ):
        values = [
            next(
                r["warm_median_wall_s"]
                for r in resource["full_trial_layouts"]
                if r["case"] == c and r["source_time_s"] == w and r["layout"] == layout
            )
            for c, w in names
        ]
        axs[1, 0].bar(x + (n - 1.5) * 0.19, np.array(values) * 1e3, 0.18, color=color, label=layout)
        arrays["layout_" + layout + "_warm_wall_s"] = np.array(values)
    arrays["layout_cases"] = np.array([c for c, _ in names])
    arrays["layout_source_time_s"] = np.array([w for _, w in names])
    axs[1, 0].set(
        xticks=x,
        xticklabels=[f"{c}\n+{(w - 3) * 1e3:.1f} ms" for c, w in names],
        ylabel="Complete 5 µs trial wall time (ms)",
        title="Changing full/half1/half2 work",
    )
    axs[1, 0].legend(fontsize=8)
    for n, key, color, label in [
        (0, "CPU_s", RED, "CPU"),
        (1, "GPU_end_to_end_s", BLUE, "GPU incl. transfers/sync"),
    ]:
        v = [r["median"][key] * 1e3 for r in resource["pair_force"]]
        axs[1, 1].bar(np.arange(2) + (n - 0.5) * 0.3, v, 0.28, color=color, label=label)
        arrays["pair_force_" + key] = np.array(v) * 1e-3
    axs[1, 1].set(
        xticks=[0, 1],
        xticklabels=["B2: 186 bubbles", "B3: 189 bubbles"],
        ylabel="Pair-force evaluation (ms)",
        title="Force-only timing, not trajectory speedup",
    )
    axs[1, 1].legend(fontsize=8)
    axs[1, 1].set_ylim(0, 1.25 * max(r["median"]["CPU_s"] * 1e3 for r in resource["pair_force"]))
    save(
        fig,
        "G8_completion_and_resource_performance",
        "Completion ratios require identical initial states and complete equal physical duration. At least three completed repeats per method required for final performance PASS. Deterministic horizontal spreading distinguishes overlapping repeat markers without changing measured times. Layout and pair-force timings are separate microbenchmarks.",
        arrays,
    )


def g8b():
    gate = source(BASE / "final/MOTION_SOLVER_R3_GATE.json")
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 5), layout="constrained")
    arrays = {}
    for i, r in enumerate(gate["protected_R2_repeat_performance"]):
        for offset, key, color, label in [
            (-0.12, "new_wall_s", BLUE, "R3"),
            (0.12, "protected_R2_wall_s", PURPLE, "Protected R2"),
        ]:
            v = np.array(r[key])
            spread = np.linspace(-0.04, 0.04, len(v)) if len(v) > 1 else np.zeros(len(v))
            axs[0].scatter(i + offset + spread, v, color=color, label=label if i == 1 else None, s=42)
            if len(v):
                axs[0].plot([i + offset - 0.06, i + offset + 0.06], [np.median(v)] * 2, color=color, lw=2)
            arrays[r["case"] + "_5ms_" + key] = v
        if not r["protected_R2_wall_s"]:
            axs[0].text(
                i,
                0.9,
                "No complete R2 baseline",
                transform=axs[0].get_xaxis_transform(),
                ha="center",
                fontsize=8,
            )
    axs[0].set(
        xticks=[0, 1],
        xticklabels=["B2", "B3"],
        ylabel="Complete +5 ms process wall time (s)",
        title="R3 / R2 | every complete timing repeat",
    )
    axs[0].legend(fontsize=8)
    for i, case in enumerate(["B2", "B3"]):
        r = next(
            r
            for r in gate["rows"]
            if r["case"] == case and r["duration_s"] == 0.001 and r["method"] == "nominal"
        )
        values = [
            r["charged_wall_s"],
            r["performance_vs_frozen"].get("reference_wall_s"),
            r["performance_vs_R2"].get("reference_wall_s"),
        ]
        for k, (v, color, label) in enumerate(
            zip(values, [BLUE, RED, PURPLE], ["R3", "Fresh Stage4L", "R2"])
        ):
            if v is not None:
                axs[1].bar(i + (k - 1) * 0.24, v, width=0.22, color=color, label=label if i == 0 else None)
        arrays[case + "_1ms_wall_s"] = np.array([np.nan if v is None else v for v in values])
    axs[1].set(
        xticks=[0, 1],
        xticklabels=["B2", "B3"],
        ylabel="Complete +1 ms process wall time (s)",
        title="One completed sample per method",
    )
    axs[1].legend(fontsize=8)
    save(
        fig,
        "G8b_R2_protection_and_1ms_regression",
        "Only complete equal-duration same-initial-state runs are shown. R2 B2 +5 ms remains incomplete and is not extrapolated. At 1 ms, B2 R2 uses the preserved historical execution, while B3 R2 is rerun in the protected campaign; single-sample timing variation is not established. NaN in the standalone plot-data NPZ denotes unavailable timing, not zero; strict audit JSON uses null.",
        arrays,
    )


if __name__ == "__main__":
    g4()
    g4b()
    g5c()
    g6()
    g7()
    g7b()
    g8()
    g8b()
