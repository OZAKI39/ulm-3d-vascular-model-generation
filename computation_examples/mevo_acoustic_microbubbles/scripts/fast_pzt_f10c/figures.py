"""Publication figures from actual saved full-particle states and measured timings."""

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from .common import BASE, ROOT, js, write, sha, identity

OUT = BASE / "figures"
PALETTE = ["#41366d", "#3679a6", "#af6385", "#8b8995"]


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.titleweight": "bold",
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "axes.edgecolor": "#303039",
            "text.color": "#20202a",
            "xtick.color": "#303039",
            "ytick.color": "#303039",
            "savefig.facecolor": "white",
            "figure.facecolor": "white",
            "grid.alpha": 0.2,
        }
    )


def finish(fig, name, title, note):
    fig.suptitle(title, fontsize=15, fontweight="bold", y=0.99)
    fig.text(0.02, 0.018, note, fontsize=9, color="#42424b", va="bottom")
    fig.subplots_adjust(top=0.90, bottom=0.12, wspace=0.36, hspace=0.55)
    if name == "Fig_F10C_1_performance_breakdown":
        fig.subplots_adjust(left=0.16, right=0.98, top=0.86, bottom=0.23, wspace=0.95)
    fig.savefig(OUT / (name + ".png"), dpi=190)
    fig.savefig(OUT / (name + ".pdf"))
    plt.close(fig)


def figure1(result):
    profile = js(BASE / "profile_v1/RESULTS.json")
    row = result["exclusive_performance"]["batch"][0]["serial"][0]
    fig, ax = plt.subplots(1, 3, figsize=(16, 5.5), gridspec_kw={"width_ratios": [1, 1.25, 1]})
    values = [
        result["exclusive_performance"]["batch"][0]["shared_static_prepare_s"],
        row["model_prepare_s"],
        row["graph_capture_s"],
        row["integration_quality_wall_s"],
        row["sync_readback_wall_s"],
        row["morphology_quality_output_prepare_s"],
        row["hdf5_output_s"],
    ]
    names = [
        "Static inputs",
        "Model construction",
        "Graph capture",
        "Integration + quality",
        "State readback",
        "Morphology / quality",
        "HDF5 output",
    ]
    ax[0].barh(names[::-1], values[::-1], color=PALETTE[0])
    ax[0].set_xlabel("Wall time (s)")
    ax[0].set_title("A  Full32 output cost (B = 1)")
    grouped = {
        "Field + response": 0.0,
        "Primary force": 0.0,
        "Secondary + reduction": 0.0,
        "CFD queries": 0.0,
        "Lubrication + balance": 0.0,
        "PBD contact": 0.0,
        "Wall / boundaries": 0.0,
        "Acceptance / clock": 0.0,
    }
    for key, value in profile["activities"].items():
        if "field_response" in key:
            group = "Field + response"
        elif "consistent_primary" in key:
            group = "Primary force"
        elif "pair_buffer" in key or "pair_sum" in key:
            group = "Secondary + reduction"
        elif "tetra_query" in key:
            group = "CFD queries"
        elif "contact_" in key:
            group = "PBD contact"
        elif any(s in key for s in ["local_mobility", "component_", "predict_", "flow_free"]):
            group = "Lubrication + balance"
        elif any(s in key for s in ["wall_", "boundaries", "geometry", "commit_events"]):
            group = "Wall / boundaries"
        else:
            group = "Acceptance / clock"
        grouped[group] += value["event_span_ms"]
    labels = list(grouped)
    vals = list(grouped.values())
    ax[1].barh(labels[::-1], vals[::-1], color=PALETTE[1])
    ax[1].set_xlabel("Activity spans (ms / 20 steps)")
    ax[1].set_title("B  Instrumented 20-step profile")
    x = np.arange(2)
    r = profile["comparison"]
    w = 0.32
    ax[2].bar(x - w / 2, [v["wall_s"] * 1000 for v in r], w, color=PALETTE[0], label="Wall clock")
    ax[2].bar(
        x + w / 2, [v["cuda_event_span_s"] * 1000 for v in r], w, color=PALETTE[2], label="CUDA Event span"
    )
    ax[2].set_xticks(x, ["Graph", "No Graph"])
    ax[2].set_ylabel("Time for 100 steps (ms)")
    ax[2].set_title("C  Existing CUDA Graph")
    ax[2].set_ylim(0, 80)
    ax[2].legend(frameon=False)
    ax[2].text(
        0.5,
        0.5,
        "Positions: bitwise equal",
        ha="center",
        transform=ax[2].transAxes,
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9),
    )
    finish(
        fig,
        "Fig_F10C_1_performance_breakdown",
        "F10C | Measured execution costs of the full-particle M0 solver",
        "RTX 4090; N = 32; dt = 10 us. Panel A: exclusive replay. B/C: original bounded profile and Graph comparison.\nPanel B is instrumented CUDA-Event activity timing, not pure kernel duration. Pure GPU kernel duration was not measured.",
    )


def figure2(result):
    rows = result["exclusive_performance"]["batch"]
    cold = result["exclusive_performance"]["cold"]
    x = np.arange(len(rows))
    w = 0.25
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    sizes = [r["B"] for r in rows]
    ax[0, 0].bar(x - w, [r["serial_total_s"] for r in rows], w, label="Serial, resident", color=PALETTE[3])
    ax[0, 0].bar(x, [r["batch_total_s"] for r in rows], w, label="Batch, resident", color=PALETTE[0])
    ax[0, 0].bar(
        x + w, [r["full_subprocess_wall_s"] for r in cold], w, label="Batch, fresh process", color=PALETTE[1]
    )
    ax[0, 0].set_yscale("log")
    ax[0, 0].set_ylabel("Measured elapsed time (s)")
    ax[0, 0].set_title("A  Identical full-particle output")
    ax[0, 0].legend(frameon=False, fontsize=9)
    ax[0, 1].plot(
        x, [r["valid_samples_per_gpu_hour"] for r in rows], "o-", color=PALETTE[0], label="Resident batch"
    )
    ax[0, 1].plot(
        x,
        [r["valid_samples_per_gpu_hour"] for r in cold],
        "s--",
        color=PALETTE[1],
        label="Fresh-process batch",
    )
    ax[0, 1].set_ylabel("Numerically valid trajectories / GPU hour")
    ax[0, 1].set_title("B  Exclusive GPU throughput")
    ax[0, 1].legend(frameon=False)
    ax[1, 0].plot(
        x,
        [r["warp_allocator_high_water_MiB"] for r in cold],
        "o-",
        color=PALETTE[0],
        label="F10C Warp allocator high-water",
    )
    right = ax[1, 0].twinx()
    right.plot(
        x,
        [r["peak_observed_gpu_MiB"] for r in cold],
        "s--",
        color=PALETTE[2],
        label="Whole-device sampled NVML peak",
    )
    ax[1, 0].set_ylabel("F10C allocator memory (MiB)")
    right.set_ylabel("Whole-device memory (MiB)", color=PALETTE[2])
    ax[1, 0].set_title("C  Allocations and whole-device memory")
    ax[1, 0].legend(loc="upper left", frameon=False, fontsize=8)
    right.legend(loc="center left", frameon=False, fontsize=8)
    bars = ax[1, 1].bar(x, [r["resident_speedup"] for r in rows], color=PALETTE[0])
    ax[1, 1].axhline(2, color=PALETTE[2], ls="--", label="2x reference threshold")
    for i, (bar, r) in enumerate(zip(bars, rows)):
        ax[1, 1].text(
            i,
            bar.get_height() + 0.15,
            f"{r['resident_speedup']:.2f}x\n{r['B']}/{r['B']} bitwise",
            ha="center",
            fontsize=9,
        )
    ax[1, 1].set_ylim(0, max(r["resident_speedup"] for r in rows) * 1.3)
    ax[1, 1].set_ylabel("Observed serial / batch time ratio")
    ax[1, 1].set_title("D  Bitwise equality and measured speedup")
    ax[1, 1].legend(frameon=False)
    for a in ax.flat:
        a.set_xticks(x, sizes)
        a.set_xlabel("Batch size B")
        a.grid(axis="y")
    finish(
        fig,
        "Fig_F10C_2_batch_throughput",
        "F10C | Full-particle batching: exclusive GPU performance replay",
        "Each trajectory: 32 IDs, 10 ms, 201 saved states, identical morphology/quality output. One timing observation per size; sampled resource audit.\nBENCHMARK_ONLY replicated/translated states; training_ready = false. Earlier shared-load timings are retained separately.",
    )


def figure3(result):
    rows = result["multirate"]
    modes = ["Q0", "Q1", "Q2", "Q3"]
    controls = ["OFF", "S_1000", "O_1000"]
    x = np.arange(4)
    w = 0.24
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    for c, key in enumerate(controls):
        rr = [next(r for r in rows if r["mode"] == mode and r["control"] == key) for mode in modes]
        ax[0, 0].bar(
            x + (c - 1) * w,
            [r["errors_vs_Q0"]["max_centroid_distance_m"] * 1e9 for r in rr],
            w,
            label=key.replace("_1000", ""),
            color=PALETTE[c],
        )
        ax[1, 1].bar(x + (c - 1) * w, [r["timing"]["resident_body_s"] for r in rr], w, color=PALETTE[c])
    ax[0, 0].axhline(10, color="k", ls="--", lw=1, label="10 nm group gate")
    ax[0, 0].set_ylabel("Maximum centroid difference (nm)")
    ax[0, 0].set_title("A  Group motion vs Q0")
    ax[0, 0].legend(frameon=False, fontsize=8)
    for i, other in enumerate(["OFF", "S_1000"]):
        vals = [result["modes"][mode]["paired"][other]["delta_z_m"] for mode in modes]
        ax[0, 1].bar(
            x + (i - 0.5) * 0.32,
            [np.nan if v is None else v * 1e6 for v in vals],
            0.32,
            color=PALETTE[i],
            label="O - " + other.replace("_1000", ""),
        )
    ax[0, 1].text(2, 0.16, "S rejected:\nno 10 ms pair", ha="center", color="#9d3756", fontsize=8)
    ax[0, 1].set_ylim(0, 0.26)
    ax[0, 1].set_ylabel("Additional Z displacement (µm)")
    ax[0, 1].set_title("B  Matched 10 ms endpoints")
    ax[0, 1].legend(frameon=False, fontsize=8)
    for i, (metric, label, gate) in enumerate(
        [("max_rg_relative", "Rg", 0.01), ("max_aspect_relative", "Aspect", 0.02)]
    ):
        vals = [max(r["errors_vs_Q0"][metric] for r in rows if r["mode"] == mode) * 100 for mode in modes]
        ax[0, 2].bar(x + (i - 0.5) * 0.32, vals, 0.32, color=PALETTE[i], label=label)
        ax[0, 2].axhline(gate * 100, color=PALETTE[i], ls="--", lw=1)
    ax[0, 2].set_ylabel("Maximum relative difference (%)")
    ax[0, 2].set_title("C  Shape metrics and registered gates")
    ax[0, 2].legend(frameon=False, fontsize=8)
    for i, (metric, label, gate) in enumerate(
        [("max_particle_error_m", "Same-ID position", 0.178), ("max_aligned_rmsd_m", "Kabsch RMSD", 0.089)]
    ):
        vals = [max(r["errors_vs_Q0"][metric] for r in rows if r["mode"] == mode) * 1e6 for mode in modes]
        ax[1, 0].bar(x + (i - 0.5) * 0.32, vals, 0.32, color=PALETTE[i], label=label)
        ax[1, 0].axhline(gate, color=PALETTE[i], ls="--", lw=1)
    ax[1, 0].set_ylabel("Maximum particle difference (µm)")
    ax[1, 0].set_title("D  Internal rearrangement differs")
    ax[1, 0].set_ylim(0, 4.2)
    ax[1, 0].legend(frameon=False, fontsize=8)
    ax[1, 1].set_ylabel("Accuracy-run body time (s)")
    ax[1, 1].set_title("E  Recorded costs; no admitted speed claim")
    ax[1, 2].axis("off")
    ax[1, 2].set_title("F  Separate numerical admission")
    cells = []
    for mode in modes:
        m = result["modes"][mode]
        cells.append(
            [
                mode,
                "PASS" if m["group_pass"] else "FAIL",
                "PASS" if m["morphology_pass"] else "FAIL",
                "BASELINE" if mode == "Q0" else "REJECT",
            ]
        )
    table = ax[1, 2].table(
        cellText=cells, colLabels=["Mode", "Group", "Morph.", "Full model"], loc="center", cellLoc="center"
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.9)
    ax[1, 2].text(
        0.5,
        0.08,
        "Q2-S: first rejection at 9.16 ms\nOverlap 17.865 nm > 17.8 nm budget\nNo q = 4 or batch + cache production run",
        transform=ax[1, 2].transAxes,
        ha="center",
        fontsize=9,
    )
    for a in list(ax.flat)[:-1]:
        a.set_xticks(x, modes)
        a.grid(axis="y")
    finish(
        fig,
        "Fig_F10C_3_multirate_accuracy_speed",
        "F10C | External-query caching preserves centroid motion better than morphology",
        "Q0: (acoustic, flow) stride (1,1); Q1: (2,1); Q2: (1,2); Q3: (2,2). Mechanics and CFD proposal validity checked every 10 us.\nQ2-S error bars use its accepted prefix through 9.15 ms. Exact edge-event timing is required for morphology; all ON-cache modes fail this gate.",
    )


def figure4():
    paths = [BASE / "batch_v1/B1_serial_000", BASE / "batch64_v1/B64_batch", BASE / "multirate_v1/Q1_O_1000"]
    titles = ["Serial M0 reference", "GPU batch B = 64, member 0", "Q1 acoustic cache (NOT ADMITTED)"]
    fig = plt.figure(figsize=(15, 9))
    norm = colors.Normalize(0, 31)
    cmap = matplotlib.colormaps["viridis"]
    az = np.linspace(0, 2 * np.pi, 18)
    el = np.linspace(0, np.pi, 10)
    sx = 1.78 * np.outer(np.cos(az), np.sin(el))
    sy = 1.78 * np.outer(np.sin(az), np.sin(el))
    sz = 1.78 * np.outer(np.ones_like(az), np.cos(el))
    axes = []
    for col, path in enumerate(paths):
        with h5py.File(path / "trajectory.h5") as h:
            for row, frame in enumerate([0, -1]):
                x = h["positions_m"][frame, 0]
                xc = x.mean(0)
                rel = (x - xc) * 1e6
                rg = float(h["metrics/global_rg_m"][frame, 0]) * 1e6
                e = h["metrics/global_eigenvalues_m2"][frame, 0]
                aspect = np.sqrt(e[0] / e[-1])
                largest = int(h["metrics/largest_size"][frame, 0])
                ax = fig.add_subplot(2, 3, row * 3 + col + 1, projection="3d")
                axes.append(ax)
                for i, point in enumerate(rel):
                    ax.plot_surface(
                        sx + point[0],
                        sy + point[1],
                        sz + point[2],
                        color=cmap(norm(i)),
                        linewidth=0,
                        shade=True,
                        alpha=0.93,
                        rasterized=True,
                    )
                distances = np.linalg.norm(rel[:, None] - rel[None, :], axis=-1)
                for i, j in zip(*np.where(np.triu(distances <= 2.1 * 1.78, 1))):
                    ax.plot(*rel[[i, j]].T, color="#43414b", lw=0.6, alpha=0.45)
                ax.scatter([0], [0], [0], marker="+", color="black", s=60, depthshade=False)
                ax.set(
                    xlim=(-8, 8),
                    ylim=(-8, 8),
                    zlim=(-8, 8),
                    xticks=[-5, 0, 5],
                    yticks=[-5, 0, 5],
                    zticks=[-5, 0, 5],
                )
                ax.set_box_aspect((1, 1, 1))
                ax.view_init(23, -58)
                ax.set_xlabel("x - Xc (µm)", labelpad=0)
                ax.set_ylabel("y - Yc (µm)", labelpad=0)
                ax.set_zlabel("z - Zc (µm)", labelpad=0)
                ax.tick_params(labelsize=8, pad=0)
                ax.set_title(
                    f"{titles[col]}\nt = {0 if frame == 0 else 10} ms | Rg {rg:.3f} um | aspect {aspect:.3f}\nLargest component: {largest}/32",
                    fontsize=10,
                    pad=0,
                )
                for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
                    axis.pane.fill = False
    fig.subplots_adjust(top=0.89, bottom=0.12, left=0.01, right=0.94, hspace=0.35, wspace=0.04)
    cbax = fig.add_axes([0.955, 0.24, 0.012, 0.51])
    cb = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbax)
    cb.set_label("Fixed bubble ID")
    cb.set_ticks([0, 8, 16, 24, 31])
    fig.suptitle(
        "F10C | Real full-particle morphology, common scale and fixed IDs",
        fontsize=15,
        fontweight="bold",
        y=0.99,
    )
    fig.text(
        0.025,
        0.022,
        "Original radius R = 1.78 um; each frame is centered on its actual centroid (no shape rescaling). Serial and batch positions are bitwise identical.\nAll 32 particles remain connected; no component split or merge was observed. Q1 changes individual arrangements despite a small centroid error.",
        fontsize=10,
    )
    for ext in ["png", "pdf"]:
        fig.savefig(OUT / ("Fig_F10C_4_particle_morphology_comparison." + ext), dpi=190)
    plt.close(fig)


def main():
    OUT.mkdir(exist_ok=True)
    style()
    result = js(BASE / "analysis/RESULTS.json")
    figure1(result)
    figure2(result)
    figure3(result)
    figure4()
    sources = js(BASE / "analysis/SOURCES.json")
    for path in [
        BASE / "profile_v1/RESULTS.json",
        BASE / "cold_v1/RESULTS.json",
        BASE / "analysis/RESULTS.json",
        BASE / "batch64_v1/B64_batch/trajectory.h5",
        BASE / "failure_v1/REJECTION.json",
        BASE / "exclusive_v1/COMPLETE.json",
        BASE / "exclusive_v1/RESOURCE_AUDIT.json",
        BASE / "exclusive_v1_batch/RESULTS.json",
        BASE / "exclusive_v1_cold/RESULTS.json",
    ]:
        sources[str(path.relative_to(ROOT))] = sha(path)
    write(
        OUT / "FIGURE_SOURCES.json",
        dict(
            execution_commit=identity()["execution_commit"],
            script_sha256=sha(__file__),
            data_sha256=sources,
            scientific_execution_commits={
                p.parent.name: js(p)["execution_commit"] for p in BASE.glob("*/IDENTITY.json")
            },
            spatial_data="Actual saved full32 coordinates; identical sphere radius and axis limits; no interpolation, no G1 template",
            performance_scope="Figures 1A/2: exclusive_v1 sampled-resource replay; earlier shared-load results retained separately. One timing observation per size.",
        ),
    )


if __name__ == "__main__":
    main()
