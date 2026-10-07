"""Reconstruct accepted-step topology and plot only completed recorded data."""

import json
import h5py
import numpy as np
from scipy.sparse.csgraph import connected_components
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .common import BASE, ROOT, write

FIG = BASE / "figures"
ON = BASE / "balanced_steps_v1/dt_20us"
OFF = BASE / "balanced_controls_v1/formation_off"
C = BASE / "balanced_morphodynamics_v1/controlled"
FIXED = BASE / "balanced_morphodynamics_v1/fixed_from_switch"


def read(run):
    with h5py.File(run / "trajectory.h5") as h:
        data = {k: h[k][:] for k in ["time_s", "positions_m", "radii_m", "control_index"]}
        data.update({k: h["metrics/" + k][:] for k in h["metrics"]})
        data["controls"] = json.loads(h.attrs["controls_json"])
    return data


def partition(labels):
    return {frozenset(np.flatnonzero(labels == i)) for i in np.unique(labels)}


def topology(run):
    """Simultaneous edge changes are applied together, without artificial ordering events."""
    with h5py.File(run / "trajectory.h5") as h:
        x = h["positions_m"][0]
        radius = h["radii_m"][0]
        edge_t = h["edge_events/time_s"][:]
        edge_id = h["edge_events/batch_i_j_linked"][:]
        frame_t = h["time_s"][:]
        saved_labels = h["metrics/labels"][:]
    graph = np.linalg.norm(x[:, :, None] - x[:, None, :], axis=-1) <= 2.1 * radius
    labels = [connected_components(g, directed=False)[1] for g in graph]
    events, cursor, frame_count = [], 0, 0
    for f, frame_time in enumerate(frame_t):
        while cursor < len(edge_t) and edge_t[cursor] <= frame_time + 1e-14:
            t = edge_t[cursor]
            end = int(np.searchsorted(edge_t, t, side="right"))
            changed = np.unique(edge_id[cursor:end, 0])
            for b, i, j, linked in edge_id[cursor:end]:
                graph[b, i, j] = graph[b, j, i] = bool(linked)
            for b in changed:
                old = labels[b]
                new = connected_components(graph[b], directed=False)[1]
                for kind, left, right in [("merge", old, new), ("split", new, old)]:
                    for c in np.unique(right):
                        ids = np.flatnonzero(right == c)
                        other = np.unique(left[ids])
                        if len(other) > 1:
                            events.append(
                                dict(
                                    kind=kind,
                                    time_s=float(t),
                                    batch=int(b),
                                    particle_ids=ids.tolist(),
                                    participant_groups=[np.flatnonzero(left == a).tolist() for a in other],
                                )
                            )
                labels[b] = new
            cursor = end
        for b, lab in enumerate(labels):
            assert partition(lab) == partition(saved_labels[f, b]), (str(run), f, b)
            frame_count += 1
    assert cursor == len(edge_t)
    result = dict(
        scope="Accepted-substep topology times, not continuous-time crossing estimates; simultaneous changes grouped",
        source="trajectory.h5:edge_events",
        radius_m=float(radius),
        cluster_gap_over_radius=0.1,
        edge_change_count=len(edge_t),
        checked_sample_partitions=frame_count,
        merge_count=sum(e["kind"] == "merge" for e in events),
        split_count=sum(e["kind"] == "split" for e in events),
        events=events,
    )
    write(run / "TOPOLOGY.json", result)
    return result


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "savefig.dpi": 250,
            "pdf.fonttype": 42,
            "axes.grid": True,
            "grid.alpha": 0.16,
            "legend.frameon": False,
        }
    )


def finish(fig, stem):
    fig.savefig(FIG / (stem + ".png"), facecolor="white")
    fig.savefig(FIG / (stem + ".pdf"), facecolor="white")
    plt.close(fig)


def figure2(on, off):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), layout="constrained")
    for ax, key, title, factor in zip(
        axes,
        ["mean_nearest_neighbor_m", "cluster_count", "largest_size"],
        ["Mean nearest-neighbor distance (µm)", "Number of components", "Largest component (bubbles)"],
        [1e6, 1, 1],
    ):
        for data, color, label in [(on, "#175db1", "Secondary ON"), (off, "#a3602d", "Secondary OFF")]:
            ax.plot(data["time_s"] * 1e3, data[key][:, 0] * factor, color=color, lw=2, label=label)
        ax.set(xlabel="Physical time (ms)", ylabel=title)
    axes[0].axhline(2 * on["radii_m"][0] * 1e6, ls=":", color="0.35", label="Bubble diameter")
    axes[0].legend(fontsize=8)
    fig.suptitle("Secondary interactions drive aggregation | matched initial positions, N = 32", fontsize=12)
    finish(fig, "Fig2_Secondary_ON_OFF")
    np.savez_compressed(
        FIG / "Fig2_data.npz",
        time_s=on["time_s"],
        **{
            f"{k}_{side}": d[k]
            for side, d in [("on", on), ("off", off)]
            for k in ["mean_nearest_neighbor_m", "cluster_count", "largest_size"]
        },
    )


def figure3(d, fixed, events):
    fig, ax = plt.subplots(2, 3, figsize=(13, 7.5), layout="constrained")
    t = d["time_s"] * 1e3
    tf = 80 + fixed["time_s"] * 1e3
    colors = ["#294c9b", "#008472", "#b56a24", "#855fa8"]
    control = d["controls"]
    for k, name in enumerate(["Y+", "Y−", "Z+", "Z−"]):
        a = [e["amplitude_Pa"][k] * bool(e["enabled"][k]) for e in control]
        ph = [e["phase_rad"][k] * 180 / np.pi for e in control]
        ax[0, 0].step([0, 80, 200], [*a, a[-1]], where="post", color=colors[k], label=name)
        ax[0, 1].step([0, 80, 200], [*ph, ph[-1]], where="post", color=colors[k], label=name)
    ax[0, 0].set(ylabel="Source amplitude (Pa)", title="a  PZT amplitudes / element enable")
    ax[0, 1].set(ylabel="Source phase (degrees)", title="b  PZT phases (2.5 MHz)")
    ax[0, 0].legend(ncol=4, fontsize=8)
    c0 = d["global_centroid_m"][0, 0]
    for k, name in enumerate(["x", "y", "z"]):
        ax[0, 2].plot(t, (d["global_centroid_m"][:, 0, k] - c0[k]) * 1e6, color=colors[k], label=name)
        ax[0, 2].plot(tf, (fixed["global_centroid_m"][:, 0, k] - c0[k]) * 1e6, color=colors[k], ls="--", lw=1)
    ax[0, 2].set(ylabel="Centroid displacement (µm)", title="c  Whole-cloud centroid")
    ax[0, 2].legend(ncol=3, fontsize=8)
    ax[1, 0].plot(t, d["global_rg_m"][:, 0] * 1e6, color=colors[0], label="All bubbles")
    ax[1, 0].plot(t, d["largest_rg_m"][:, 0] * 1e6, color=colors[1], label="Largest component")
    ax[1, 0].plot(tf, fixed["largest_rg_m"][:, 0] * 1e6, color=colors[1], ls="--")
    ax[1, 0].set(ylabel="Radius of gyration (µm)", title="d  Radius of gyration")
    ax[1, 0].legend(fontsize=8)
    for k in range(3):
        ax[1, 1].plot(t, d["largest_eigenvalues_m2"][:, 0, k] * 1e12, color=colors[k], label=f"λ{k + 1}")
        ax[1, 1].plot(tf, fixed["largest_eigenvalues_m2"][:, 0, k] * 1e12, color=colors[k], ls="--", lw=1)
    ax[1, 1].set(ylabel="Shape eigenvalues (µm²)", title="e  Largest component shape")
    ax[1, 1].legend(ncol=3, fontsize=8)
    ax[1, 2].plot(t, d["cluster_count"][:, 0], color=colors[0], label="Components")
    ax[1, 2].plot(t, d["largest_size"][:, 0], color=colors[1], label="Largest size")
    for kind, marker, y, color in [("merge", "|", -1, colors[2]), ("split", "|", -3, colors[3])]:
        when = [e["time_s"] * 1e3 for e in events["events"] if e["kind"] == kind]
        if when:
            ax[1, 2].scatter(when, [y] * len(when), marker=marker, color=color, s=18, label=kind.title())
    ax[1, 2].set(ylabel="Number of components / bubbles", title="f  Membership and topology")
    ax[1, 2].legend(ncol=2, fontsize=8)
    for a in ax.flat:
        a.axvline(80, color="0.45", ls=":", lw=1)
        a.set(xlabel="Physical time (ms)", xlim=(0, 200))
    fig.suptitle(
        "Real PZT basis: control switch at 80 ms | dashed: unchanged-control continuation\n"
        "EXPLORATORY · identical switch state · component membership can change",
        fontsize=12,
    )
    finish(fig, "Fig3_PZT_morphodynamics")
    keys = [
        "time_s",
        "global_centroid_m",
        "largest_rg_m",
        "global_rg_m",
        "largest_eigenvalues_m2",
        "global_shape_tensor_m2",
        "principal_axes",
        "axis_identifiable",
        "cluster_count",
        "largest_size",
    ]
    np.savez_compressed(
        FIG / "Fig3_data.npz",
        **{f"{k}_{side}": a[k] for side, a in [("controlled", d), ("fixed", fixed)] for k in keys},
    )


def figure4():
    old = json.loads((BASE / "time_step_v1/RESULT.json").read_text())["rows"]
    final = json.loads((BASE / "balanced_steps_v1/RESULT.json").read_text())["rows"]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), layout="constrained")
    labels = [
        "Accepted substeps (50 ms run)",
        "Simulation wall time (s)",
        "Rg relative RMS difference (%)",
        "Maximum residual overlap (nm)",
        "Maximum contact-pass correction (nm)",
        "Sum |correction| / sum |physical step|",
    ]
    for rows, label, color, fmt in [
        (old, "v1 screening, 4 passes (superseded)", "0.62", "o--"),
        (final, "v2 balanced mobility, 12 passes", "#126c64", "s-"),
    ]:
        x = np.array([r["nominal_dt_s"] * 1e6 for r in rows])
        ys = [
            [r["actual_substep_count"] for r in rows],
            [r["wall_s"] for r in rows],
            [r["comparison"]["global_rg_m_relative_rms"] * 100 for r in rows],
            [r["quality"]["max_accepted_penetration_m"] * 1e9 for r in rows],
            [r["quality"]["max_contact_pass_correction_m"] * 1e9 for r in rows],
            [r["quality"]["correction_fraction"] for r in rows],
        ]
        for ax, y in zip(axes.flat, ys):
            ax.plot(x, y, fmt, color=color, label=label, ms=5)
    for i, (ax, label) in enumerate(zip(axes.flat, labels)):
        ax.set(xscale="log", xlabel="Macro step (µs)", ylabel=label, xticks=[1, 2, 5, 10, 20, 50])
        ax.set_xticklabels([1, 2, 5, 10, 20, 50])
        ax.axvline(20, color="#126c64", alpha=0.3, lw=1)
        if i in [0, 1, 3, 4]:
            ax.set_yscale("log")
    axes[1, 0].axhline(17.8, color="#a14b31", ls=":", label="Predeclared limit: 0.01 R")
    axes[1, 0].legend(fontsize=8, loc="upper left")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle(
        "Step selection: 20 µs macro / one physical substep / 12 contact passes\n"
        "Each version uses its own 1 µs reference; v1 is retained diagnostic evidence",
        fontsize=12,
    )
    finish(fig, "Fig4_step_performance")
    write(FIG / "Fig4_data.json", dict(preliminary_superseded=old, final_balanced=final))


def main():
    FIG.mkdir(exist_ok=True)
    style()
    events = {}
    for run in [
        ON,
        OFF,
        C,
        FIXED,
        BASE / "balanced_throughput_v1/N32_1s",
        BASE / "balanced_throughput_v1/batch8",
    ]:
        events[str(run.relative_to(BASE))] = topology(run)
    on, off, c, fixed = [read(p) for p in [ON, OFF, C, FIXED]]
    np.testing.assert_array_equal(on["positions_m"][0], off["positions_m"][0])
    off_drift = float(np.max(np.linalg.norm(off["positions_m"] - off["positions_m"][0], axis=-1)))
    assert off_drift == 0
    drift = []
    for name in ["basic_v1/N32_secondary_True", "time_step_v1/dt_1us", "contact12_v1/dt_20us"]:
        d = read(BASE / name)
        drift.append(
            dict(
                run=name,
                maximum_uniform_field_centroid_drift_m=float(
                    np.linalg.norm(d["global_centroid_m"][:, 0] - d["global_centroid_m"][0, 0], axis=-1).max()
                ),
            )
        )
    write(
        BASE / "PROTOTYPE_RETIREMENT.json",
        dict(
            status="SUPERSEDED_MOBILITY_DRIFT",
            reason="Truncated local mobility lacked component net-force balance. Final v2 adds a common component velocity correction without changing relative velocities.",
            affected_prefixes=[
                "smoke",
                "basic_v1",
                "time_step_v1",
                "contact8_v1",
                "contact12_v1",
                "short_performance",
                "morphodynamics_v1",
                "core_layout_v2",
            ],
            retained=True,
            measurements=drift,
            final_model="fast_pzt_f1_v2_balanced_mobility",
        ),
    )
    write(
        BASE / "ANALYSIS_CHECKS.json",
        dict(
            status="PASS",
            matched_OFF_max_displacement_m=off_drift,
            topology_checks={k: {a: v for a, v in r.items() if a != "events"} for k, r in events.items()},
            final_uniform_centroid_drift_m=float(
                np.linalg.norm(on["global_centroid_m"][:, 0] - on["global_centroid_m"][0, 0], axis=-1).max()
            ),
        ),
    )
    (ROOT / "configs/fast_pzt/scene_c.json").write_text((C / "CONFIG.json").read_text())
    figure2(on, off)
    figure3(c, fixed, events[str(C.relative_to(BASE))])
    figure4()
    print("Topology partitions verified; Figures 2–4 written", flush=True)


if __name__ == "__main__":
    main()
