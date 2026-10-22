"""F3 scientific panels from unmodified force/trajectory archives."""

import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from scripts.fast_pzt_f2.analyze import read
from .common import BASE, js
from .analyze import BRANCHES

FIG = BASE / "figures"
COLORS = ["#2869ac", "#9b553d", "#7262a6"]


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.labelcolor": "#242932",
            "text.color": "#242932",
            "axes.edgecolor": "#353c46",
            "xtick.color": "#353c46",
            "ytick.color": "#353c46",
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def finish(fig, name):
    for ax in fig.axes:
        if ax.axison:
            ax.grid(alpha=0.13)
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"{name}.{ext}", dpi=225, bbox_inches="tight")
    plt.close(fig)


def driving():
    d = js(BASE / "CONTROL_FORCE_DIAGNOSTICS.json")
    names = ["unchanged", d["selected"]["phase"], d["selected"]["amplitude"]]
    labels = ["Unchanged", "Phase only", "Amp/enable"]
    scales = js(FIG / "VECTOR_SCALES.json")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 9), layout="constrained")
    ax = axes.ravel()
    for k, (n, label) in enumerate(zip(names, labels)):
        ax[k].imshow(Image.open(FIG / f"force_{n}.png"))
        ax[k].axis("off")
        kind = "Primary force" if k == 0 else "Centered control-induced Primary"
        ax[k].set_title(
            f"{chr(97 + k)}  {label}: {kind}\nArrow scale: 1 pN = {scales[n]['arrow_micrometers_per_pN']:.1f} µm",
            fontsize=10,
        )
    rows = [d["candidate_rows"][n] for n in names]
    idx = np.arange(3)
    for j, (key, label) in enumerate([("D_P_N", "Differential Primary, DP"), ("D_S_N", "Secondary, DS")]):
        ax[3].bar(
            idx + (j - 0.5) * 0.3,
            [r["force"][key] * 1e12 for r in rows],
            width=0.28,
            label=label,
            color=[COLORS[0], COLORS[1]][j],
        )
    ax[3].set(
        xticks=idx,
        xticklabels=labels,
        yscale="log",
        ylabel="RMS force (pN)",
        title="d  Direct differential drive is small",
    )
    ax[3].legend(fontsize=8)
    for j, (key, label) in enumerate(
        [("free_velocity_m_s", "Free Stokes"), ("mobility_velocity_m_s", "After local mobility")]
    ):
        ax[4].bar(
            idx + (j - 0.5) * 0.3,
            [r["kinematics"][key]["relative_velocity_RMS_m_s"] * 1e3 for r in rows],
            width=0.28,
            label=label,
            color=[COLORS[0], COLORS[2]][j],
        )
    ax[4].set(
        xticks=idx,
        xticklabels=labels,
        ylabel="Centered predicted velocity RMS (mm/s)",
        title="e  Pre-contact motion; not observed strain",
    )
    ax[4].legend(fontsize=8)
    ax[5].axis("off")
    ax[5].text(
        0,
        1,
        "f  Resolution and interpretation\n\nGrid spacing: 100 µm; occupied cells: 2\n32 actual positions; R = 1.78 µm\n\nDP / DS ≈ 3.6–5.0 × 10⁻⁵\nStored-gradient vs pressure-derivative\nL2 mismatch: 24–38% (eligible controls)\n\nAmplitude control: 72.3% of Secondary\nforce-change energy is global scaling.\n\nLocal mobility also leaves a torque-balance\nresidual; it is not an exact coupled solve.\nFrozen physics retained; EXPLORATORY.",
        va="top",
        fontsize=10,
        linespacing=1.35,
    )
    fig.suptitle(
        "F3-1 | Differential acoustic driving at the same mature state\nActual spheres; arrows alone are scaled · surface-anchored vectors · no particle displacement",
        fontsize=14,
    )
    finish(fig, "Fig_F3_1_differential_driving")
    with h5py.File(BASE / "driving_v3/FORCE_FIELD_DATA.h5") as h:
        arrays = {
            "positions_m": h["positions_m"][:],
            "particle_id": h["particle_id"][:],
            "radii_m": h["radii_m"][:],
        }
        for n in names:
            for k in [
                "primary_N",
                "secondary_N",
                "total_N",
                "control_delta_centered_primary_N",
                "free_velocity_m_s",
                "mobility_velocity_m_s",
            ]:
                arrays[n + "_" + k] = h[n + "/" + k][:]
    np.savez_compressed(FIG / "Fig_F3_1_data.npz", **arrays)


def reproducibility():
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), layout="constrained")
    ax = axes.ravel()
    export = {}
    report = js(BASE / "NUMERICAL_REPRODUCIBILITY.json")
    for j, mode in enumerate(["fast_atomic", "stable_reduction"]):
        color = COLORS[j]
        label = ["F2 atomic", "Stable reduction"][j]
        for i in range(3):
            d = read(BASE / "repeats_v1" / f"{mode}_{i}")
            t = d["time_s"] * 1e3
            ax[0].plot(
                t,
                d["global_rg_m"][:, 0] * 1e6,
                color=color,
                lw=1 if j == 0 else 1.4,
                alpha=0.7,
                label=label if i == 0 else None,
            )
            for e, ls in enumerate(["-", "--", ":"]):
                ax[1].plot(
                    t, d["global_eigenvalues_m2"][:, 0, e] * 1e12, color=color, ls=ls, lw=0.8, alpha=0.7
                )
            export[f"{mode}_{i}_rg_m"] = d["global_rg_m"]
            export[f"{mode}_{i}_eigenvalues_m2"] = d["global_eigenvalues_m2"]
        e = np.load(BASE / "repeats_v1" / f"{mode}_envelope.npz")
        ax[2].plot(e["time_s"] * 1e3, e["aligned_internal_rmsd_m"] * 1e6, color=color, label=label)
        ax[3].plot(e["time_s"] * 1e3, e["neighbor_disagreement"], color=color, label=label)
        walls = report["modes"][mode]["walls_s"]
        ax[4].bar(j, np.mean(walls), color=color, alpha=0.7)
        ax[4].scatter(np.full(3, j), walls, c="#20252c", s=16, zorder=3)
        export[mode + "_aligned_RMSD_m"] = e["aligned_internal_rmsd_m"]
        export[mode + "_neighbor_disagreement"] = e["neighbor_disagreement"]
    ax[0].set(
        xlabel="Time since mature state (ms)",
        ylabel="Global radius of gyration (µm)",
        title="a  Three repeats per mode",
    )
    ax[0].legend(fontsize=8)
    ax[1].set(
        xlabel="Time since mature state (ms)",
        ylabel="Eigenvalues (µm²)",
        title="b  λ1 / λ2 / λ3: solid / dashed / dotted",
    )
    ax[2].set(
        xlabel="Time since mature state (ms)",
        ylabel="Maximum repeat-pair aligned RMSD (µm)",
        title="c  Same IDs; remove best-fit rigid motion",
    )
    ax[2].legend(fontsize=8)
    ax[3].set(
        xlabel="Time since mature state (ms)",
        ylabel="Maximum disagreeing neighbor edges",
        title="d  Repeated graph consistency",
    )
    ax[4].set(
        xticks=[0, 1],
        xticklabels=["F2 atomic", "Stable reduction"],
        ylabel="Simulation wall time (s)",
        title="e  32 bubbles, 120 ms",
    )
    ax[5].axis("off")
    ax[5].text(
        0,
        1,
        "f  What the checks resolve\n\nEight fixed-state repetitions:\nAtomic pair force varies by ≤7.76 × 10⁻²⁶ N.\nStable pair force is bitwise identical.\nComponent-only stabilization retains variation.\n\nThree full stable trajectories:\nAll positions and sorted edge events identical.\nGraph and direct short replay also agree.\n\nSeparate 1 pm / 1 nm initial perturbations\nstill yield ≈3.42 / 3.43 µm internal differences.\nRepeatability does not imply physical accuracy.",
        va="top",
        fontsize=10,
        linespacing=1.4,
    )
    fig.suptitle(
        "F3-2 | Numerical reproducibility with frozen physics\nR = 1.78 µm · Δt = 10 µs · same saved 32-bubble state · same GPU",
        fontsize=14,
    )
    finish(fig, "Fig_F3_2_numerical_reproducibility")
    np.savez_compressed(FIG / "Fig_F3_2_data.npz", time_s=d["time_s"], **export)


def response():
    fig, axes = plt.subplots(3, 3, figsize=(14, 11), layout="constrained")
    ax = axes.ravel()
    report = js(BASE / "MORPHODYNAMICS_COMPARISON.json")
    unc = np.load(BASE / "SHAPE_UNCERTAINTY.npz")
    export = {}
    names = ["unchanged", "phase", "amplitude"]
    labels = ["Unchanged", "Phase only", "Amplitude / enable"]
    for j, (name, label) in enumerate(zip(names, labels)):
        m = np.load(BRANCHES[name] / "MORPHOLOGY.npz")
        c = np.load(BRANCHES[name] / "CONTRAST.npz")
        t = m["time_s"] * 1e3
        color = COLORS[j]
        ax[0].plot(t, c["centroid_distance_m"] * 1e6, c=color, label=label)
        ax[1].plot(t, m["fixed_rg_m"] * 1e6, c=color, label=label)
        ax[2].plot(t, c["eigenvalues_relative"] * 100, c=color, label=label)
        ax[3].plot(t, c["aligned_internal_rmsd_m"] * 1e6, c=color, label=label)
        ax[4].plot(t, m["principal_axis_angle_deg"], c=color, lw=0.8)
        ax[7].plot(t, m["fixed_aspect_ratio"], c=color, label=label)
        export.update({name + "_" + k: m[k] for k in m.files})
    ax[0].set(
        xlabel="Time since switch (ms)",
        ylabel="Centroid distance vs unchanged (µm)",
        title="a  Resolved translation",
    )
    ax[0].legend(fontsize=8)
    ax[1].set(
        xlabel="Time since switch (ms)",
        ylabel="Fixed-ID radius of gyration (µm)",
        title="b  Cluster size; stable repeats overlap",
    )
    for panel, key, factor in [(2, "eigenvalues_relative", 100), (3, "aligned_internal_rmsd_m", 1e6)]:
        ax[panel].fill_between(
            unc["time_s"] * 1e3,
            0,
            unc[key] * factor,
            color=".45",
            alpha=0.2,
            label="Historical model / initial-perturbation envelope",
        )
        ax[panel].plot(
            unc["time_s"] * 1e3,
            3 * unc["repeat_" + key] * factor,
            c="#20252c",
            ls=":",
            label="3 × stable-repeat difference",
        )
    ax[2].set(
        xlabel="Time since switch (ms)",
        ylabel="Relative eigenvalue difference (%)",
        title="c  Shape contrast and sensitivity envelope",
    )
    ax[2].legend(fontsize=7)
    ax[3].set(
        xlabel="Time since switch (ms)",
        ylabel="Rigid-aligned internal RMSD (µm)",
        title="d  Rearrangement is not robust to sensitivity",
    )
    ax[4].set(
        xlabel="Time since switch (ms)",
        ylabel="Principal-axis angle from mature state (deg)",
        title="e  Ambiguous axes masked; not rigid rotation",
        ylim=(0, 90),
    )
    for j, name in enumerate(["phase", "amplitude"]):
        intervals = report["effects"][name]["metrics"]["eigenvalues_relative"][
            "above_repeat_model_initial_windows"
        ]["intervals"]
        ax[5].broken_barh(
            [(v["start_s"] * 1e3, max(0.2, (v["end_s"] - v["start_s"]) * 1e3)) for v in intervals],
            (j - 0.22, 0.44),
            facecolors=COLORS[j + 1],
        )
    ax[5].set(
        xlabel="Time since switch (ms)",
        yticks=[0, 1],
        yticklabels=["Phase only", "Amp/enable"],
        xlim=(0, 120),
        ylim=(-0.6, 1.6),
        title="f  Eigenvalue contrast above all envelopes",
    )
    for j in range(3):
        ax[6].broken_barh([(-10, 10)], (j - 0.22, 0.44), facecolors="#bbc1c9")
        ax[6].broken_barh([(0, 120)], (j - 0.22, 0.44), facecolors=COLORS[j])
    ax[6].axvline(0, c="#303640", ls=":")
    ax[6].set(
        xlabel="Time relative to PZT switch (ms)",
        yticks=[0, 1, 2],
        yticklabels=labels,
        xlim=(-10, 120),
        title="g  Switch at mature absolute time 102.2 ms",
    )
    ax[7].set(
        xlabel="Time since switch (ms)",
        ylabel="Aspect ratio √(λ1 / λ3)",
        title="h  No established lasting elongation",
    )
    ax[8].axis("off")
    ax[8].text(
        0,
        1,
        "i  Interpretation\n\nLongest eigenvalue windows above envelopes:\nphase 2.6 ms; amplitude 8.8 ms\nBoth shorter than the 10 ms sustained rule.\n\nTRANSIENT_SHAPE_RESPONSE\nNO_ROBUST_ENDPOINT_SHAPE_RESPONSE\n\nAll 32 members remain connected.\nPersistent split / merge: 0 / 0\nNo membership exchange.\n\nSame-control stable-repeat differences are zero\nfor raw states; SVD residuals are roundoff only.\nEXPLORATORY; field and mobility limits remain.",
        va="top",
        fontsize=9.5,
        linespacing=1.3,
    )
    fig.suptitle(
        "F3-3 | PZT control versus the same-state unchanged trajectory\nStable reduction · 120 ms continuation · fixed criteria and separately measured sensitivities",
        fontsize=14,
    )
    finish(fig, "Fig_F3_3_PZT_control_response")
    np.savez_compressed(FIG / "Fig_F3_3_data.npz", **export)


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    style()
    driving()
    reproducibility()
    response()
