"""Four-panel sets grounded in archived, unfiltered per-bubble F2 data."""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from .common import BASE
from .analyze import read, js, BRANCHES

FIG = BASE / "figures"
COLORS = ["#2c64ad", "#aa503b", "#7161a9"]


def style():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.labelsize": 10,
            "axes.titlesize": 11,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#30343b",
            "text.color": "#20252c",
            "axes.labelcolor": "#20252c",
            "xtick.color": "#30343b",
            "ytick.color": "#30343b",
            "grid.alpha": 0.15,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "legend.frameon": False,
        }
    )


def finish(fig, name):
    for ax in fig.axes:
        ax.grid(alpha=0.15)
    for ext in ["png", "pdf"]:
        fig.savefig(FIG / f"{name}.{ext}", dpi=230, bbox_inches="tight")
    plt.close(fig)


def near():
    fig, axes = plt.subplots(2, 3, figsize=(13, 7.4), layout="constrained")
    ax = axes.ravel()
    curves = np.load(BASE / "core_v1/PAIR_FORCE_CURVES.npz")
    d = curves["separation_m"] / 3.56e-6
    for j, (model, label) in enumerate([("f1", "F1 smooth"), ("f2", "F2 feedback")]):
        for k, ph in enumerate(["in_phase", "counterphase"]):
            ax[0].plot(
                d,
                curves[f"{ph}_{model}_N"] * 1e12,
                color=COLORS[j],
                ls="-" if k == 0 else "--",
                label=f"{label}, " + ("in phase" if k == 0 else "opposite phase"),
            )
    ax[0].axhline(0, c=".5", lw=0.5)
    ax[0].axvline(1, c=".5", lw=0.6)
    ax[0].set(
        xlabel=r"Separation / $(R_i+R_j)$",
        ylabel="Signed pair force (pN)",
        title="a  Pair interaction; positive = attraction",
    )
    ax[0].legend(fontsize=8)
    export = {k: curves[k] for k in curves.files}
    for j, (model, label) in enumerate(
        [("f1_smooth_near_field", "F1 smooth"), ("f2_near_field_candidate", "F2 feedback")]
    ):
        a = read(BASE / "core_v1" / model)
        x = a["positions_m"][:, 0]
        sep = np.linalg.norm(x[:, 1] - x[:, 0], axis=-1)
        t = a["time_s"]
        ax[1].plot(t * 1e3, sep * 1e6, c=COLORS[j], label=label)
        ax[2].plot((t[1:] + t[:-1]) * 500, -np.diff(sep) / np.diff(t) * 1e3, c=COLORS[j], label=label)
        export[f"{model}_pair_t_s"] = t
        export[f"{model}_pair_separation_m"] = sep
    ax[1].axhline(3.56, c=".5", ls=":", label="Geometric contact")
    ax[1].legend(fontsize=8)
    ax[1].set(xlabel="Time (ms)", ylabel="Separation (µm)", title="b  Actual two-bubble approach")
    ax[2].set(xlabel="Time (ms)", ylabel="Approach speed (mm/s)", title="c  Saved-interval relative speed")
    for j, (name, label) in enumerate([("f1_10us", "F1 near field"), ("near_10us", "F2 near field")]):
        a = read(BASE / "sensitivity_v1" / name)
        t = a["time_s"] * 1e3
        for panel, key, ylabel in [
            (3, "cluster_count", "Connected components"),
            (4, "largest_size", "Largest component (bubbles)"),
            (5, "global_rg_m", "Global radius of gyration (µm)"),
        ]:
            y = a[key][:, 0] * (1e6 if panel == 5 else 1)
            ax[panel].plot(t, y, c=COLORS[j], label=label)
            ax[panel].set(xlabel="Time (ms)", ylabel=ylabel)
            export[f"{name}_{key}"] = a[key]
        export[f"{name}_time_s"] = a["time_s"]
    ax[3].set_title("d  Same 32-bubble initial positions")
    ax[4].set_title("e  Same contact and 10 µs step")
    ax[5].set_title("f  Size is less sensitive than packing")
    ax[5].legend(fontsize=9)
    fig.suptitle(
        "F2-1 | Near-field model comparison\nR = 1.78 µm · 2.5 MHz · uniform 1 kPa · EXPLORATORY", fontsize=14
    )
    finish(fig, "Fig_F2_1_near_field")
    np.savez_compressed(FIG / "Fig_F2_1_data.npz", **export)


def contact():
    fig, axes = plt.subplots(2, 3, figsize=(13, 8), layout="constrained")
    ax = axes.ravel()
    export = {}
    report = js(BASE / "CONTACT_SENSITIVITY.json")
    for j, (name, label) in enumerate(
        [("near_5us", "F1 central PBD"), ("contact_5us", "F2 mobility candidate")]
    ):
        d = read(BASE / "sensitivity_v1" / name)
        t = d["time_s"] * 1e3
        q = d["quality"]
        ax[0].plot(t, q[:, :, 3].max(axis=1) * 1e9, c=COLORS[j], label=label)
        ax[1].plot(t, q[:, :, 4].sum(axis=1) * 1e3, c=COLORS[j], label=label)
        ax[2].plot(t, d["global_rg_m"][:, 0] * 1e6, c=COLORS[j], label=label)
        a = np.load(BASE / "sensitivity_v1" / name / "ADJACENCY.npz")
        # Final panels use actual mesh sphere geometry, with the same camera and limits.
        img = Image.open(FIG / ("contact_F1_PBD.png" if j == 0 else "contact_F2_candidate.png"))
        ax[3 + j].imshow(img)
        ax[3 + j].axis("off")
        ax[3 + j].set_title(("d  " if j == 0 else "e  ") + label + " · 50 ms", fontsize=11)
        export[f"{name}_time_s"] = d["time_s"]
        export[f"{name}_quality"] = q
        export[f"{name}_positions_m"] = d["positions_m"]
        export[f"{name}_neighbors"] = a["linked"]
    ax[0].axhline(17.8, c="#5c626c", ls=":", label="Unchanged 0.01R budget")
    ax[0].legend(fontsize=8)
    ax[0].set(
        xlabel="Time (ms)",
        ylabel="Maximum accepted overlap (nm)",
        title="a  Cumulative maximum, matched 5 µs",
    )
    ax[1].set(
        xlabel="Time (ms)", ylabel="Summed contact correction (mm)", title="b  Accumulated over all particles"
    )
    ax[2].set(
        xlabel="Time (ms)", ylabel="Global radius of gyration (µm)", title="c  Same force; distinct packing"
    )
    ax[2].legend(fontsize=8)
    ax[5].axis("off")
    s = report["comparison"]
    text = (
        "f  Matched comparison and disposition\n\n"
        "Time to Nmax ≥ 16: 10.5 / 10.5 ms\nFinal components: 1 / 1 (32 bubbles)\n"
        f"Final Rg difference: {s['final_rg_relative_difference'] * 100:.2f}%\nEigenvalue difference: {s['final_eigenvalues_relative_difference'] * 100:.2f}%\n"
        f"ID-aligned RMSD: {s['final_ID_aligned_rmsd_m'] * 1e6:.2f} µm\n\n"
        "CONTACT_MODEL_SENSITIVE\nCandidate rejected as default:\nnoncentral corrections induce rotation.\n\nSecondary OFF: no motion or aggregation.\nFinal packing remains unvalidated."
    )
    ax[5].text(0.03, 0.98, text, transform=ax[5].transAxes, va="top", fontsize=10, linespacing=1.45)
    fig.suptitle(
        "F2-2 | Contact and packing sensitivity\nTrue radius and identical 3D view in d–e · force model fixed · EXPLORATORY",
        fontsize=14,
    )
    finish(fig, "Fig_F2_2_contact_sensitivity")
    np.savez_compressed(FIG / "Fig_F2_2_data.npz", **export)


def morphodynamics():
    fig, axes = plt.subplots(3, 3, figsize=(14, 11), layout="constrained")
    ax = axes.ravel()
    report = js(BASE / "PZT_MORPHODYNAMICS.json")
    switch = report["maturity"]["absolute_time_s"] * 1e3
    names = ["unchanged", "phase", "amplitude"]
    labels = ["Unchanged", "Phase only", "Amplitude / enable"]
    a0 = report["controls"]["unchanged"]["amplitude_Pa"]
    a1 = report["controls"]["amplitude"]["event"]["amplitude_Pa"]
    ph = np.degrees(report["controls"]["phase"]["event"]["phase_rad"])
    for j in range(4):
        ax[0].step(
            [0, switch, switch + 120],
            [a0[j], a1[j], a1[j]],
            where="post",
            color=["#325b93", "#ad6840", "#40765e", "#8758a0"][j],
            ls=["-", "--", ":", "-."][j],
            label=["Y+", "Y−", "Z+", "Z−"][j],
        )
    ax[0].axvline(switch, c=".4", lw=0.8)
    ax[0].set(
        xlabel="Absolute physical time (ms)",
        ylabel="Source amplitude (Pa)",
        title="a  Amplitude / enable control",
    )
    ax[0].legend(ncol=2, fontsize=8)
    ax[0].text(
        0.03, 0.96, f"Switch after maturity: {switch:.1f} ms", va="top", transform=ax[0].transAxes, fontsize=8
    )
    ax[1].axis("off")
    ax[1].text(
        0,
        1,
        "b  Phase-only control\n\nY+ / Y− / Z+ / Z−\n"
        + f"0 / 0 / 0 / 0 → {ph[0]:.1f} / {ph[1]:.2f} / 0 / 0 deg\n\n"
        + "Amplitudes unchanged: 151.34 Pa each\nAll four elements remain enabled\nROI pressure RMS matched: 843.77 Pa\n27 fixed diagnostic points; all controls\n\nSame mature state: 32 / 32 bubbles\n10 ms stable membership before switch",
        va="top",
        fontsize=10,
        linespacing=1.3,
    )
    export = {}
    for j, (name, label) in enumerate(zip(names, labels)):
        m = np.load(BRANCHES[name] / "MORPHOLOGY.npz")
        c = np.load(BRANCHES[name] / "CONTRAST.npz")
        t = m["time_s"] * 1e3
        color = COLORS[j]
        ax[2].plot(
            t,
            np.linalg.norm(m["fixed_centroid_m"] - m["fixed_centroid_m"][0], axis=1) * 1e6,
            c=color,
            label=label,
        )
        ax[3].plot(t, m["fixed_rg_m"] * 1e6, c=color, label=label)
        for k, ls in enumerate(["-", "--", ":"]):
            ax[4].plot(t, m["fixed_eigenvalues_m2"][:, k] * 1e12, c=color, ls=ls, lw=0.85)
        ax[5].plot(t, m["principal_axis_angle_deg"], c=color, label=label, lw=0.8)
        ax[6].plot(t, c["aligned_internal_rmsd_m"] * 1e6, c=color, label=label)
        a = np.load(BRANCHES[name] / "ADJACENCY.npz")
        ax[7].plot(t, a["jaccard_with_initial"], c=color, label=label, lw=0.9)
        export.update({f"{name}_{k}": m[k] for k in m.files})
        export[f"{name}_adjacency_jaccard"] = a["jaccard_with_initial"]
    repeat = np.load(BRANCHES["unchanged_repeat"] / "CONTRAST.npz")
    ax[6].fill_between(
        repeat["time_s"] * 1e3,
        0,
        3 * repeat["aligned_internal_rmsd_m"] * 1e6,
        color=".5",
        alpha=0.17,
        label="3 × same-control repeat difference",
    )
    ax[2].set(
        xlabel="Time since mature-state switch (ms)",
        ylabel="Centroid displacement (µm)",
        title="c  Translation from mature state",
    )
    ax[2].legend(fontsize=8)
    ax[3].set(
        xlabel="Time since switch (ms)", ylabel="Fixed-ID radius of gyration (µm)", title="d  Overall size"
    )
    ax[4].set(
        xlabel="Time since switch (ms)",
        ylabel="Tensor eigenvalue (µm²)",
        title="e  Shape eigenvalues: λ1 / λ2 / λ3",
    )
    ax[5].set(
        xlabel="Time since switch (ms)",
        ylabel="Principal-axis angle (deg)",
        title="f  Relative to mature axis; ambiguous masked",
        ylim=(0, 90),
    )
    ax[6].set(
        xlabel="Time since switch (ms)",
        ylabel="Aligned RMSD vs unchanged (µm)",
        title="g  Internal change after rigid alignment",
    )
    ax[6].legend(fontsize=7, loc="upper left")
    ax[7].set(
        xlabel="Time since switch (ms)",
        ylabel="Neighbor Jaccard index",
        title="h  Fixed-ID neighbors vs mature state",
        ylim=(0, 1.05),
    )
    ax[8].axis("off")
    ax[8].text(
        0,
        1,
        "i  Topology and inference\n\nAll branches: one 32-bubble component\nPersistent merge / split: 0 / 0\nTransient component merge / split: 0 / 0\nMembership exchange: none\n\nRaw pair-edge changes:\nunchanged 1479; phase 633; amplitude 979\nThese do not change component membership.\n\nNO_SIGNIFICANT_SHAPE_RESPONSE\nEndpoint shape effects do not exceed\nthe registered repeat envelope.\nTransient eigenvalue windows: up to\n5.4 ms (amplitude), 5.6 ms (phase).\nOne repeat; descriptive evidence only.",
        va="top",
        fontsize=10,
        linespacing=1.25,
    )
    fig.suptitle(
        "F2-3 | PZT control of an already mature cluster\nSame 102.2 ms state · 120 ms continuation · actual basis field · EXPLORATORY",
        fontsize=14,
    )
    finish(fig, "Fig_F2_3_PZT_morphodynamics")
    np.savez_compressed(FIG / "Fig_F2_3_data.npz", **export)


if __name__ == "__main__":
    FIG.mkdir(exist_ok=True)
    style()
    near()
    contact()
    morphodynamics()
