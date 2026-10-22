"""Four reproducible publication figures from archived measurements; no synthetic science."""

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import BoundaryNorm, ListedColormap
from .common import BASE, ROOT, F5, js, write, sha, identity
from .postprocess_reference import process as reference
from mevo_mb.fast_sim.f10d_metrics import process as optimized
from .validate import compare_results

OUT = BASE / "figures"
PURPLE, BLUE, GRAY, RED = "#6345A0", "#3885AD", "#A7AFB9", "#B75060"
COLORS = ["#bec8d3", "#8da8b9", "#4d849d", "#3e9f9b", "#695199", "#bb8eaf", "#c9b68a", "#d8d8d8"]
SOURCES = {}


def read(rel):
    path = BASE / rel
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return js(path)


def title(ax, text):
    ax.set_title(text, loc="left", fontweight="bold", pad=12)
    ax.spines[["top", "right"]].set_visible(False)


def save(fig, name, note):
    fig.text(0.035, 0.016, note, fontsize=9, color="#333333")
    fig.savefig(OUT / (name + ".png"), dpi=220, facecolor="white")
    fig.savefig(OUT / (name + ".pdf"), facecolor="white")
    plt.close(fig)


def paired(ax, old, new, labels, ylabel):
    xx = np.arange(len(labels))
    for shift, data, color, name in [
        (-0.18, old, GRAY, "Reference pipeline"),
        (0.18, new, PURPLE, "F10D pipeline"),
    ]:
        bars = ax.bar(xx + shift, data, width=0.34, color=color, label=name)
        ax.bar_label(bars, fmt="%.2f", fontsize=9, padding=3)
    ax.set_xticks(xx, labels)
    ax.set_ylabel(ylabel)
    ax.set_ylim(0, max(max(old), max(new)) * 1.24)


def fig1(val, bench):
    fig, axs = plt.subplots(2, 2, figsize=(13.4, 9.5))
    fig.subplots_adjust(left=0.09, right=0.97, bottom=0.12, top=0.94, wspace=0.3, hspace=0.55)
    row = val["rows"][1]
    old, new = row["reference_timing"], row["optimized_timing"]
    metric_keys = [
        "centroid_rg_s",
        "global_tensor_eigen_s",
        "pair_distances_components_s",
        "component_shape_pair_statistics_s",
        "sampled_lineage_s",
    ]
    shape_keys = ["fixed_tensor_eigen_s", "kabsch_s", "aspect_axis_dynamic_s"]
    paired(
        axs[0, 0],
        [
            old["metrics_and_sampled_lineage_s"],
            old["fixed_dynamic_shape_kabsch_s"],
            old["accepted_step_topology_s"],
        ],
        [sum(new[k] for k in metric_keys), sum(new[k] for k in shape_keys), new["accepted_step_topology_s"]],
        ["Frame metrics\n+ sampled lineage", "Fixed/dynamic\nshape + Kabsch", "Accepted-step\ntopology"],
        "Offline wall time (s)",
    )
    title(axs[0, 0], "a  Same saved B64 arrays: CPU postprocessing")
    axs[0, 0].legend(frameon=False, fontsize=9)
    labels = [
        "Centroid + Rg",
        "Tensor / eigen /\nshape statistics",
        "Kabsch",
        "Pairs + components",
        "Accepted-step topology",
    ]
    data = [
        new["centroid_rg_s"],
        sum(
            new[k]
            for k in [
                "global_tensor_eigen_s",
                "fixed_tensor_eigen_s",
                "component_shape_pair_statistics_s",
                "aspect_axis_dynamic_s",
            ]
        ),
        new["kabsch_s"],
        new["pair_distances_components_s"],
        new["accepted_step_topology_s"],
    ]
    ax = axs[0, 1]
    bars = ax.barh(np.arange(5), data, color=[BLUE, BLUE, BLUE, PURPLE, PURPLE])
    ax.bar_label(bars, fmt="%.3f", padding=4, fontsize=9)
    ax.set_yticks(np.arange(5), labels, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.2)
    ax.set_xlabel("Optimized wall time (s)")
    title(ax, "b  Remaining CPU cost")

    ax = axs[1, 0]
    cold = next(r for r in bench["rows"] if r["B"] == 64 and r["startup"] == "cold")
    keys = [
        "static_prepare_s",
        "model_prepare_s",
        "graph_capture_s",
        "integration_quality_wall_s",
        "sync_readback_wall_s",
        "morphology_quality_output_prepare_s",
        "hdf5_output_s",
    ]
    names = [
        "Static inputs",
        "Model + module",
        "CUDA Graph",
        "Integration + quality",
        "Readback",
        "Morphology + diagnostics",
        "HDF5",
    ]
    bottom = np.zeros(2)
    for key, name, color in zip(keys, names, COLORS):
        values = [cold["modes"][m][key] for m in ["legacy", "optimized"]]
        ax.bar([0, 1], values, bottom=bottom, width=0.55, color=color, label=name)
        bottom += values
    total = np.array([cold["modes"][m]["full_process_wall_s"] for m in ["legacy", "optimized"]])
    ax.bar([0, 1], total - bottom, bottom=bottom, width=0.55, color=COLORS[-1], label="Import / JSON / other")
    for i, t in enumerate(total):
        ax.text(i, t + 0.12, f"{t:.2f} s", ha="center", fontsize=10)
    ax.set_xticks([0, 1], ["Reference", "F10D"])
    ax.set_ylabel("Complete fresh-process time (s)")
    ax.set_ylim(0, 10)
    ax.legend(frameon=False, fontsize=8.3, loc="upper right", bbox_to_anchor=(1.02, 0.94))
    title(ax, "c  Full B64 workflow, matched GPU window")
    storage = read("storage_v1/RESULTS.json")
    ax = axs[1, 1]
    paired(
        ax,
        [storage["rows"][0]["hdf5_write_s"], storage["rows"][0]["read_verify_s"]],
        [storage["rows"][1]["hdf5_write_s"], storage["rows"][1]["read_verify_s"]],
        ["Write", "Read + verify"],
        "Wall time (s)",
    )
    handles, _ = ax.get_legend_handles_labels()
    ax.legend(
        handles,
        ["Original gzip 1", "Shuffle + gzip 1 (rejected)"],
        frameon=False,
        fontsize=9,
        loc="upper left",
    )
    ax.set_ylim(0, 0.85)
    ax.text(0.04, 0.81, "26.13 vs 27.57 MB; all 60 datasets exact", transform=ax.transAxes, fontsize=9)
    title(ax, "d  Lossless output candidate: no benefit")
    save(
        fig,
        "Fig_F10D_1_postprocessing_breakdown",
        "B64 = 64 x 32 bubbles; 10 ms, 201 frames. Reference centroid/Rg subcosts were not isolated.\nCPU timings exclude imports; cold total includes them. Final pipeline retains original gzip 1. No quality checks removed.",
    )


def fig2(bench, family):
    fig, axs = plt.subplots(2, 2, figsize=(12.8, 9.2))
    fig.subplots_adjust(left=0.085, right=0.97, bottom=0.12, top=0.94, wspace=0.3, hspace=0.48)
    rows = bench["rows"]
    old, new, labels = [], [], []
    for r in rows:
        key = "full_process_wall_s" if r["startup"] == "cold" else "complete_run_wall_s"
        old.append(r["modes"]["legacy"][key])
        new.append(r["modes"]["optimized"][key])
        labels.append(f"B{r['B']}\n{r['startup']}")
    paired(axs[0, 0], old, new, labels, "Complete run time (s)")
    axs[0, 0].legend(frameon=False, fontsize=9)
    title(axs[0, 0], "a  Same-window pipeline comparison")
    for i, r in enumerate(rows):
        axs[0, 0].text(
            i, max(old[i], new[i]) + 0.75, f"{r['speedup']:.2f}x", ha="center", fontsize=10, color=PURPLE
        )
    paired(
        axs[0, 1],
        [3600 * r["B"] / s / 1000 for r, s in zip(rows, old)],
        [3600 * r["B"] / s / 1000 for r, s in zip(rows, new)],
        labels,
        "Numerical-pass trajectories / hour (thousands)",
    )
    title(axs[0, 1], "b  All benchmark members passed")
    ax = axs[1, 0]
    vals = [family["serial_total_s"], family["grouped_execution_s"] + family["grouping_s"]]
    bars = ax.bar([0, 1], vals, color=[GRAY, PURPLE], width=0.55)
    ax.bar_label(bars, fmt="%.2f s", padding=3)
    ax.set_xticks([0, 1], ["20 individual runs", "3 grouped runs"])
    ax.set_ylabel("Resident complete-output wall time (s)")
    ax.set_ylim(0, 17)
    ax.text(
        0.95,
        0.83,
        f"{family['actual_grouping_speedup_vs_serial']:.2f}x\n20 / 20 serial-batch checks passed",
        transform=ax.transAxes,
        ha="right",
        fontsize=10,
    )
    title(ax, "c  Mixed OFF / S / O: actual serial comparison")
    ax = axs[1, 1]
    cold = [r for r in rows if r["startup"] == "cold"]
    paired(
        ax,
        [r["modes"]["legacy"]["nvml_sampled_peak_MiB"] for r in cold],
        [r["modes"]["optimized"]["nvml_sampled_peak_MiB"] for r in cold],
        ["B16", "B64"],
        "Sampled device GPU memory (MiB)",
    )
    ax.set_ylim(0, 1260)
    ax.text(
        0.03,
        0.91,
        "NVML: 10 ms sampling; allocator peaks: 93 / 117 MiB\n105 resource samples; no conflicting process observed",
        transform=ax.transAxes,
        fontsize=9,
        va="top",
    )
    title(ax, "d  Resource evidence")
    save(
        fig,
        "Fig_F10D_2_throughput_comparison",
        "All runs: full-particle M0, 10 µs steps, 201 frames, unchanged quality budgets. B64 cold speedup: 2.15x (target >= 1.25x).\nOne matched pair per size/startup; BENCHMARK_ONLY. Hourly rates are short-run extrapolations, not independent scientific samples.",
    )


def fig3(family, bench, diversity):
    fig, axs = plt.subplots(2, 2, figsize=(13, 9.4))
    fig.subplots_adjust(left=0.08, right=0.97, bottom=0.12, top=0.94, wspace=0.3, hspace=0.53)
    rows = sorted(family["rows"], key=lambda r: ["OFF", "S_1000", "O_1000"].index(r["control"]))
    names = ["OFF", "S", "O"]
    ax = axs[0, 0]
    ax.bar(np.arange(3), [16] * 3, color="#eeeeef", width=0.55, label="Capacity 16")
    bars = ax.bar(
        np.arange(3),
        [r["valid"] for r in rows],
        color=[GRAY, BLUE, PURPLE],
        width=0.55,
        label="Numerical pass",
    )
    ax.bar_label(
        bars,
        labels=[f"{r['valid']}/{r['batch_size']} passed\n{r['fill_fraction']:.1%} full" for r in rows],
        padding=5,
        fontsize=10,
    )
    ax.set_xticks(np.arange(3), names)
    ax.set_ylim(0, 19)
    ax.set_ylabel("Trajectories / batch")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    title(ax, "a  Partial batches; fixed input identities")
    ax = axs[0, 1]
    bottom = np.zeros(3)
    for keys, name, color in [
        (["model_prepare_s", "graph_capture_s"], "Model + graph", GRAY),
        (["integration_quality_wall_s", "sync_readback_wall_s"], "Integration / quality / readback", BLUE),
        (["morphology_quality_output_prepare_s"], "Morphology / diagnostics", PURPLE),
        (["hdf5_output_s", "serialization_and_tail_s"], "HDF5 / JSON", "#ba92ba"),
    ]:
        data = [sum(r["timing"][k] for k in keys) for r in rows]
        ax.bar(np.arange(3), data, bottom=bottom, color=color, label=name, width=0.55)
        bottom += data
    ax.set_xticks(np.arange(3), names)
    ax.set_ylim(0, 1.65)
    ax.set_ylabel("Complete run stages (s)")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    title(ax, "b  Shared static inputs loaded once")
    ax = axs[1, 0]
    ideal = next(r for r in bench["rows"] if r["startup"] == "resident" and r["B"] == 64)
    values = [
        family["total_valid_samples_per_gpu_hour"],
        family["resident_valid_samples_per_gpu_hour"],
        64 * 3600 / ideal["modes"]["optimized"]["complete_run_wall_s"],
    ]
    bars = ax.bar(np.arange(3), np.array(values) / 1000, color=[BLUE, PURPLE, GRAY], width=0.55)
    ax.bar_label(bars, fmt="%.1f", padding=4)
    ax.set_xticks(
        np.arange(3), ["Mixed 20\n+ static inputs", "Mixed 20\nresident", "Single-family B64\nresident"]
    )
    ax.set_ylim(0, 111)
    ax.set_ylabel("Numerical-pass trajectories / hour (thousands)")
    ax.text(
        0.02,
        0.93,
        f"Grouping: {family['grouping_s'] * 1000:.3f} ms; static inputs: {family['static_prepare_s']:.3f} s",
        transform=ax.transAxes,
        fontsize=9,
    )
    title(ax, "c  Mixed-family yield versus ideal B64")
    ax = axs[1, 1]
    branches = {r["control"]: r for r in diversity["rows"]}
    valid = [diversity["prepared_valid"]] + [
        branches[c]["timing"]["valid_sample_count"] for c in ["OFF", "S_1000", "O_1000"]
    ]
    total = [diversity["requested_seeds"]] + [
        len(branches[c]["sample_ids"]) for c in ["OFF", "S_1000", "O_1000"]
    ]
    ax.bar(np.arange(4), valid, color=BLUE, label="Numerical pass")
    ax.bar(np.arange(4), np.array(total) - valid, bottom=valid, color=RED, label="Failed / batch-aborted")
    ax.set_xticks(np.arange(4), ["Preparation", "OFF", "S", "O"])
    ax.set_ylim(0, 6.5)
    ax.set_yticks([0, 1, 2, 3, 4])
    ax.set_ylabel("Independent-seed attempts / branches")
    ax.legend(frameon=False, fontsize=9, loc="upper left")
    ax.text(
        0.98,
        0.73,
        "4 seeds requested\n2 reached a connected 32-bubble state\nNo seed completed all 3 controls",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=9,
    )
    title(ax, "d  Diversity test: failures retained")
    save(
        fig,
        "Fig_F10D_3_control_family_yield",
        "Panels a-c: BENCHMARK_ONLY (translated / repeated initial states); timings exclude import and regression comparisons.\nPanel d: distinct numerical arrangements under one environment. PHYSICS_UNCALIBRATED; training_ready = false.",
    )


def fig4(val, diversity):
    from matplotlib.gridspec import GridSpec

    fig = plt.figure(figsize=(14, 11.8))
    gs = GridSpec(
        3,
        12,
        figure=fig,
        left=0.075,
        right=0.97,
        bottom=0.115,
        top=0.89,
        hspace=0.53,
        wspace=0.85,
        height_ratios=[1, 1, 1.04],
    )
    initial = BASE / "campaign_v2_diversity/INITIAL_STATES.npz"
    SOURCES[str(initial.relative_to(ROOT))] = sha(initial)
    states = np.load(initial)
    source = F5 / "science_v1/main/trajectory.h5"
    SOURCES[str(source.relative_to(ROOT))] = sha(source)
    with h5py.File(source) as h:
        times = h["time_s"][:]
        x = h["positions_m"][:]
        ei = h["edge_events/batch_i_j_linked"][:]
        et = h["edge_events/time_s"][:]
    old = reference(x, 1.78e-6, times, ei, et)
    local_new = optimized(x, 1.78e-6, times, ei, et)
    comparison = compare_results(old, local_new)
    assert comparison["status"] == "PASS"
    data_path = BASE / "validate_v1/F5_RELEASE_METRICS.npz"
    SOURCES[str(data_path.relative_to(ROOT))] = sha(data_path)
    remote = np.load(data_path)
    new = dict(
        rg_m=local_new[0]["global_rg_m"],
        aspect=local_new[2]["fixed_aspect_ratio"],
        labels=local_new[0]["labels"],
        largest_size=local_new[0]["largest_size"],
    )
    np.testing.assert_array_equal(new["labels"], remote["labels"])
    np.testing.assert_array_equal(new["largest_size"], remote["largest_size"])
    np.testing.assert_allclose(
        new["aspect"],
        remote["aspect"],
        **js(ROOT / "configs/fast_pzt/f10d/protocol.json")["float_gates"]["dimensionless"],
    )
    np.testing.assert_allclose(
        new["rg_m"],
        remote["rg_m"],
        **js(ROOT / "configs/fast_pzt/f10d/protocol.json")["float_gates"]["length_m"],
    )
    write(
        OUT / "WSL_PLOT_RECOMPUTATION.json",
        dict(
            status="PASS",
            numpy=np.__version__,
            comparison=comparison,
            remote_vs_WSL_aspect_max_abs=float(np.max(np.abs(new["aspect"] - remote["aspect"]))),
            scope="Figure overlay recomputed with both implementations in the same WSL environment. Cross-environment eigensolver roundoff is compared under the preregistered tolerance; remote scientific tests retain their own results.",
        ),
    )
    np.savez_compressed(
        OUT / "F5_PLOT_VALUES.npz",
        time_s=times,
        old_rg_m=old[0]["global_rg_m"],
        new_rg_m=new["rg_m"],
        old_aspect=old[2]["fixed_aspect_ratio"],
        new_aspect=new["aspect"],
    )
    # All four clouds share the original placement center, not individual recentering.
    controls = ROOT / "results/fast_pzt/f10a/extended_10_v1/O_1000/CONFIG.json"
    SOURCES[str(controls.relative_to(ROOT))] = sha(controls)
    center = np.mean(js(controls)["environment"]["positions_m"], axis=0)
    cmap = matplotlib.colormaps["viridis"]
    for i, (seed, xx) in enumerate(zip(states["seeds"], states["positions_m"])):
        ax = fig.add_subplot(gs[0, 3 * i : 3 * i + 3])
        xx = (xx - center) * 1e6
        for point in sorted(xx, key=lambda q: q[2]):
            ax.add_patch(
                Circle(
                    point[:2],
                    1.78,
                    facecolor=cmap((point[2] + 12) / 24),
                    edgecolor="#263742",
                    lw=0.45,
                    alpha=0.87,
                )
            )
        ax.set(xlim=(-15, 15), ylim=(-15, 15), aspect="equal", xlabel="x offset (µm)", ylabel="y offset (µm)")
        ax.set_xticks([-10, 0, 10])
        ax.set_yticks([-10, 0, 10])
        title(ax, f"{chr(97 + i)}  Seed {seed}")
    fig.text(
        0.075,
        0.968,
        "Independent initial arrangements: XY projection, physical radius 1.78 µm, common placement origin",
        fontsize=11,
    )
    fig.text(
        0.075,
        0.945,
        f"Origin (mm): ({center[0] * 1e3:.3f}, {center[1] * 1e3:.3f}, {center[2] * 1e3:.3f}); color encodes z from -12 to +12 µm (purple to yellow).",
        fontsize=9,
    )
    for idx, (oldv, newv, ylabel, caption) in enumerate(
        [
            (
                old[0]["global_rg_m"][:, 0] * 1e6,
                new["rg_m"][:, 0] * 1e6,
                "Fixed-ID Rg (µm)",
                "e  F5 aggregation and release",
            ),
            (
                old[2]["fixed_aspect_ratio"][:, 0],
                new["aspect"][:, 0],
                "Fixed-ID aspect ratio",
                "f  Same tensor / eigen definition",
            ),
            (
                old[0]["largest_size"][:, 0],
                new["largest_size"][:, 0],
                "Largest component (bubbles)",
                "g  Accepted-step edge events",
            ),
        ]
    ):
        ax = fig.add_subplot(gs[1, 4 * idx : 4 * idx + 4])
        ax.plot(times * 1000, oldv, color=GRAY, lw=3, label="Reference")
        ax.plot(times * 1000, newv, color=PURPLE, lw=1.2, ls="--", label="F10D")
        ax.axvline(300, color=RED, ls=":", lw=1.2, label="Sound OFF")
        ax.set(xlabel="Physical time (ms)", ylabel=ylabel, xlim=(0, 600))
        title(ax, caption)
        if idx == 0:
            ax.legend(frameon=False, fontsize=8, loc="upper left")
        if idx == 2:
            ax.plot(et * 1000, np.full(len(et), -1.8), "|", color=BLUE, ms=5, alpha=0.4)
            ax.set_ylim(-4, 36)
            ax.text(
                0.97,
                0.16,
                "674 raw changes\n10 µs accepted-step timestamps",
                ha="right",
                transform=ax.transAxes,
                fontsize=8,
            )
    ax = fig.add_subplot(gs[2, :8])
    colors = ListedColormap(
        list(matplotlib.colormaps["tab20"].colors) + list(matplotlib.colormaps["tab20b"].colors[:12])
    )
    label_image = ax.imshow(
        new["labels"][:, 0, :].T,
        origin="lower",
        interpolation="nearest",
        aspect="auto",
        extent=[0, 600, -0.5, 31.5],
        cmap=colors,
        norm=BoundaryNorm(np.arange(-0.5, 32.5), 32),
    )
    fig.colorbar(
        label_image, ax=ax, pad=0.015, fraction=0.025, ticks=[0, 8, 16, 24, 31], label="Component label"
    )
    ax.axvline(300, color="black", ls=":", lw=1.3)
    ax.set(xlabel="Physical time (ms)", ylabel="Fixed particle ID", yticks=[0, 8, 16, 24, 31])
    title(ax, "h  Original and optimized connectivity labels are exactly equal")
    ax = fig.add_subplot(gs[2, 8:])
    ax.axis("off")
    top = read("validate_v1/F5_TOPOLOGY.json")[0]
    summary = read("analysis/SUMMARY.json")
    text = (
        "i  Regression evidence\n\n"
        "F5: 1,201 frames over 600 ms\nSaved frames every 0.5 ms\n32 -> 7 bubbles in largest group\n(after sound OFF at 300 ms)\n\n"
        f"{top['counts']['merge_persistent']} persistent merges\n{top['counts']['split_persistent']} persistent splits; 1 right-censored\n"
        "All labels / masks / events exact\nSame-server Rg / tensor / Kabsch exact\n"
        f"Max axis-angle error: {max(summary['non_bitwise_maxima'].values()):.2e} deg\n"
        "Preregistered tolerance: 1e-5 deg"
    )
    ax.text(0, 1, text, va="top", fontsize=10, linespacing=1.4)
    save(
        fig,
        "Fig_F10D_4_morphology_topology_validation",
        "Top row: actual generated states before preparation, no interpolation or geometric rescaling. F5 overlays use identical saved input.\nColors in panel h denote frame-local component labels, not persistent cluster identities. Graph split/merge does not mean bubble rupture/fusion.\nEXPLORATORY / PHYSICS_UNCALIBRATED; independent-seed failures retained; training_ready = false.",
    )


def main():
    OUT.mkdir(exist_ok=True)
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "axes.labelcolor": "#222222",
            "xtick.color": "#222222",
            "ytick.color": "#222222",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    val = read("validate_v1/RESULTS.json")
    bench = read("campaign_v2_benchmark/RESULTS.json")
    family = read("campaign_v2_families/RESULTS.json")
    diversity = read("campaign_v2_diversity/RESULTS.json")
    read("campaign_v2/RESOURCE_AUDIT.json")
    fig1(val, bench)
    fig2(bench, family)
    fig3(family, bench, diversity)
    fig4(val, diversity)
    write(
        OUT / "FIGURE_SOURCES.json",
        dict(
            identity=identity(),
            script_sha256=sha(__file__),
            inputs=SOURCES,
            outputs={p.name: sha(p) for p in sorted(OUT.iterdir()) if p.suffix in [".png", ".pdf", ".npz"]},
            scope="English scientific plots of measured arrays and recorded wall times; no invented data; no repeated-run confidence intervals.",
        ),
    )
    print("Generated four PNG/PDF figure pairs")


if __name__ == "__main__":
    main()
