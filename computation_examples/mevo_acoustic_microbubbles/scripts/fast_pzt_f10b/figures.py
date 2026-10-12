"""Publication plots of saved observations; no trajectory generation/interpolation."""

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from .common import ROOT, BASE, F10A, js, write, sha, identity

COLORS = ["#3366A2", "#25958B", "#D18532", "#8159A4"]
CONTROL = {"OFF": "#64748B", "S_1000": "#AA5965", "O_1000": "#256FA2"}
STYLE = {"reference": ("-", 1.8), "G0": ("--", 1.4), "G1": (":", 2.2)}
NAMES = [
    "Fig_F10B_1_contact_failure",
    "Fig_F10B_2_centroid_surrogate",
    "Fig_F10B_3_accuracy_speed",
    "Fig_F10B_4_surrogate_scope",
]


def main():
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#38424C",
            "text.color": "#202A32",
            "xtick.color": "#38424C",
            "ytick.color": "#38424C",
            "legend.frameon": False,
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )
    out = BASE / "figures"
    out.mkdir(exist_ok=True)
    result = js(BASE / "analysis/RESULTS.json")
    checks = js(BASE / "checks_v2/RESULTS.json")
    curves = np.load(BASE / "analysis/curves.npz")
    times = curves["time_s"] * 1e3
    ident = identity()
    sources = {}

    def save(fig, name):
        fig.text(
            0.012,
            0.014,
            "M0 engineering model; physics uncalibrated; no morphology prediction.   Data: F10A / F10B | plot "
            + ident["execution_commit"][:9],
            fontsize=7,
            color="#56616C",
        )
        for ext in ["png", "pdf"]:
            fig.savefig(out / (name + "." + ext), dpi=200)
        sources[name] = dict(
            plot_commit=ident["execution_commit"],
            script_sha256=sha(__file__),
            analysis_sha256=sha(BASE / "analysis/RESULTS.json"),
            curve_sha256=sha(BASE / "analysis/curves.npz"),
            checks_sha256=sha(BASE / "checks_v2/RESULTS.json"),
            science_identity=js(BASE / "rollout_v1/IDENTITY.json")["execution_commit"],
            contact_commits=[
                js(BASE / f / "IDENTITY.json")["execution_commit"] for f in ["contact_v1", "contact_v2"]
            ],
        )
        plt.close(fig)

    def title(ax, s):
        ax.set_title(s, loc="left", fontweight="bold", pad=9)

    def grid(ax):
        ax.grid(alpha=0.17, linewidth=0.6)

    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    fig.subplots_adjust(left=0.075, right=0.97, top=0.88, bottom=0.12, hspace=0.43, wspace=0.37)
    fig.suptitle("3 kPa contact diagnosis: time step and contact iterations", fontsize=16, y=0.965)
    fig.text(
        0.5,
        0.918,
        "Original force / wall models and 17.8 nm residual / 178 nm correction limits retained",
        ha="center",
        fontsize=10,
    )
    contact = []
    for row, color in zip(result["contacts"], COLORS):
        d = js(BASE / row["source"])
        ts = d["traces"]
        t = np.array([v["attempted_end_s"] for v in ts]) * 1e6
        label = f"{row['case']}: {row['dt_s'] * 1e6:g} µs / {row['rounds']} rounds"
        contact.append((row, ts, t, color, label))
        ax[0, 0].plot(
            t, [v["prediction"]["max_overlap_m"] * 1e9 for v in ts], "-o", ms=4, color=color, label=label
        )
        q = ts[0]["rounds"]
        rr = np.arange(1, len(q) + 1)
        ax[0, 1].plot(
            rr, [v["max_overlap_m"] * 1e9 for v in q], "-o", ms=3, color=color, label=label, alpha=0.85
        )
        ax[0, 2].plot(
            t, [max(v["max_correction_m"] for v in s["rounds"]) * 1e9 for s in ts], "-o", ms=4, color=color
        )
        ax[1, 0].plot(t, [v["final"]["max_overlap_m"] * 1e9 for v in ts], "-o", ms=4, color=color)
        if row["first_rejection_s"]:
            for a, vals in [
                (ax[0, 0], [s["prediction"]["max_overlap_m"] * 1e9 for s in ts]),
                (ax[0, 2], [max(v["max_correction_m"] for v in s["rounds"]) * 1e9 for s in ts]),
                (ax[1, 0], [s["final"]["max_overlap_m"] * 1e9 for s in ts]),
            ]:
                a.plot(t[-1], vals[-1], "x", ms=10, mew=2, color=color)
        fp = [np.sqrt(np.mean(np.sum(np.array(s["primary_force_N"]) ** 2, axis=1))) * 1e12 for s in ts]
        fs = [np.sqrt(np.mean(np.sum(np.array(s["secondary_force_N"]) ** 2, axis=1))) * 1e12 for s in ts]
        ax[1, 1].semilogy(t, fs, "-o", color=color, ms=3)
        ax[1, 1].semilogy(t, fp, "--", color=color)
    title(ax[0, 0], "a  Predicted overlap (before PBD)")
    ax[0, 0].set(xlabel="Attempted step end (µs)", ylabel="Maximum overlap (nm)")
    ax[0, 0].legend(fontsize=8, loc="upper left")
    title(ax[0, 1], "b  First step: convergence per round")
    ax[0, 1].set(xlabel="PBD round", ylabel="Residual overlap (nm)", yscale="log")
    ax[0, 1].axhline(17.8, color="#444", ls="--", lw=1)
    ax[0, 1].text(
        0.04, 0.04, "A0/A2 and A1/A3 share initial prefixes", transform=ax[0, 1].transAxes, fontsize=7.8
    )
    title(ax[0, 2], "c  Largest single-round correction")
    ax[0, 2].set(xlabel="Attempted step end (µs)", ylabel="Correction (nm)")
    ax[0, 2].axhline(178, color="#444", ls="--", lw=1, label="178 nm budget")
    ax[0, 2].legend(fontsize=8)
    title(ax[1, 0], "d  Residual after all PBD rounds")
    ax[1, 0].set(xlabel="Attempted step end (µs)", ylabel="Residual overlap (nm)")
    ax[1, 0].axhline(17.8, color="#444", ls="--", lw=1)
    ax[1, 0].text(0.95, 0.85, "×  First rejection", ha="right", transform=ax[1, 0].transAxes, fontsize=9)
    title(ax[1, 1], "e  Actual per-bubble RMS force")
    ax[1, 1].set(xlabel="Attempted step end (µs)", ylabel="Force (pN)")
    ax[1, 1].legend(
        handles=[
            Line2D([0], [0], color="#444", label="Secondary"),
            Line2D([0], [0], color="#444", ls="--", label="Primary"),
        ],
        fontsize=8,
        loc="center right",
    )
    title(ax[1, 2], "f  Accepted prefix and first rejection")
    ax[1, 2].set(
        xlabel="Physical time (µs)", yticks=range(4), yticklabels=["A0", "A1", "A2", "A3"], xlim=(0, 57)
    )
    for i, (row, ts, t, color, label) in enumerate(contact):
        accepted = row["accepted_duration_s"] * 1e6
        ax[1, 2].barh(i, accepted, color=color, height=0.45, alpha=0.7)
        if row["first_rejection_s"]:
            ax[1, 2].plot(t[-1], i, "x", color=color, ms=9, mew=2)
            ax[1, 2].text(t[-1] + 1.3, i, ("residual" if i < 2 else "correction"), va="center", fontsize=8)
        else:
            ax[1, 2].text(accepted + 0.8, i, "pass", va="center", fontsize=8)
    ax[1, 2].invert_yaxis()
    for a in ax.flat:
        grid(a)
    save(fig, NAMES[0])

    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    fig.subplots_adjust(left=0.072, right=0.98, top=0.81, bottom=0.12, hspace=0.45, wspace=0.33)
    fig.suptitle("Independent 10 ms centroid prediction in the original CFD field", fontsize=16, y=0.965)
    handles = [Line2D([0], [0], color=c, label=k.replace("_1000", "")) for k, c in CONTROL.items()] + [
        Line2D(
            [0],
            [0],
            color="#30343B",
            ls=style[0],
            lw=style[1],
            label=("Full M0 reference (5 µs)" if name == "reference" else name + " prediction (10 µs)"),
        )
        for name, style in STYLE.items()
    ]
    fig.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.925),
        ncol=3,
        fontsize=9,
        columnspacing=2.5,
    )
    for axis in range(3):
        a = ax[0, axis]
        for key, color in CONTROL.items():
            for mode, (ls, lw) in STYLE.items():
                c = curves[mode + "_" + key + "_centroid_m"]
                a.plot(times, (c[:, axis] - c[0, axis]) * 1e6, color=color, ls=ls, lw=lw)
        title(a, chr(97 + axis) + "  " + ["Axial X", "Transverse Y", "Target Z"][axis] + " displacement")
        a.set(xlabel="Time (ms)", ylabel="Displacement (µm)")
        grid(a)
    for i, comp in enumerate(["OFF", "S_1000"]):
        a = ax[1, i]
        for mode, (ls, lw) in STYLE.items():
            z = curves[mode + "_O_1000_centroid_m"][:, 2] - curves[mode + "_" + comp + "_centroid_m"][:, 2]
            a.plot(times, (z - z[0]) * 1e6, ls=ls, lw=lw, color="#256FA2")
        title(a, chr(100 + i) + "  O − " + comp.replace("_1000", "") + ": acoustic Z increment")
        a.set(xlabel="Time (ms)", ylabel="Additional Z displacement (µm)")
        grid(a)
    a = ax[1, 2]
    title(a, "f  Endpoint 3D prediction error")
    for j, mode in enumerate(["G0", "G1"]):
        a.bar(
            np.arange(3) + (j - 0.5) * 0.34,
            [result["models"][mode]["cases"][k]["endpoint_3D_error_m"] * 1e9 for k in CONTROL],
            width=0.32,
            color=COLORS[j],
            label=mode,
        )
    a.axhline(50, color="#444", ls="--", lw=1, label="50 nm limit")
    a.set(xticks=range(3), xticklabels=["OFF", "S", "O"], ylabel="Endpoint error (nm)")
    a.legend(fontsize=8)
    grid(a)
    save(fig, NAMES[1])

    fig, ax = plt.subplots(2, 3, figsize=(14, 8))
    fig.subplots_adjust(left=0.075, right=0.98, top=0.86, bottom=0.12, hspace=0.45, wspace=0.4)
    fig.suptitle("Accuracy passes for G1; standalone preparation limits acceleration", fontsize=15, y=0.96)
    fig.text(
        0.5,
        0.91,
        "One timing sample per control / mode; same RTX 4090; all input preparation and comparable centroid output included",
        ha="center",
        fontsize=9,
    )
    for i, mode in enumerate(["G0", "G1"]):
        a = ax[0, i]
        title(a, chr(97 + i) + f"  {mode}: signed Z prediction error")
        for key, color in CONTROL.items():
            a.plot(
                times,
                curves[mode + "_" + key + "_error_m"][:, 2] * 1e9,
                color=color,
                label=key.replace("_1000", ""),
            )
        for threshold, ls in [(20, ":"), (30, "--")]:
            a.axhline(
                threshold,
                color="#555",
                ls=ls,
                lw=0.9,
                label=f"{threshold} nm " + ("endpoint" if threshold == 20 else "path") + " limit",
            )
            a.axhline(-threshold, color="#555", ls=ls, lw=0.9)
        a.set(xlabel="Time (ms)", ylabel="Z error (nm)")
        a.legend(fontsize=7.4, ncol=2, frameon=True, facecolor="white", edgecolor="none", framealpha=0.95)
        grid(a)
    a = ax[0, 2]
    title(a, "c  Endpoint paired-increment error")
    for j, mode in enumerate(["G0", "G1"]):
        a.bar(
            np.arange(2) + (j - 0.5) * 0.34,
            [
                result["models"][mode]["increments"][k]["error_m"] * 1e9
                for k in ["O_minus_OFF", "O_minus_S_1000"]
            ],
            width=0.32,
            color=COLORS[j],
            label=mode,
        )
    a.set(xticks=range(2), xticklabels=["O − OFF", "O − S"], ylabel="Prediction − reference (nm)")
    a.axhline(0, color="#444", lw=0.8)
    a.legend(fontsize=8)
    grid(a)
    rows = [
        next(v for v in result["timing"]["rows"] if v["mode"] == mode and v["control"] == "O_1000")
        for mode in ["FULL", "G0", "G1"]
    ]
    a = ax[1, 0]
    title(a, "d  Standalone program body (O control)")
    bottom = np.zeros(3)
    for key, label, color in [
        ("preparation_s", "Input / model / graph preparation", "#9DAFC0"),
        ("integration_s", "Integration + required validity checks", "#3976A8"),
        ("sync_readback_s", "Readback", "#55A599"),
        ("output_s", "Matched output", "#D6A158"),
    ]:
        v = np.array([r[key] for r in rows])
        a.bar(range(3), v, bottom=bottom, color=color, label=label)
        bottom += v
    a.set(xticks=range(3), xticklabels=["Full M0", "G0", "G1"], ylabel="Wall-clock time (s)", ylim=(0, 1.72))
    a.legend(
        fontsize=7.5, loc="upper right", frameon=True, facecolor="white", edgecolor="none", framealpha=0.95
    )
    grid(a)
    a = ax[1, 1]
    title(a, "e  Measured speedup (same O control)")
    for j, (key, label, color) in enumerate(
        [
            ("integration_speedup", "Integration only", "#9DAFC0"),
            ("body_speedup", "Complete body", "#3976A8"),
            ("subprocess_speedup", "Including process startup", "#55A599"),
        ]
    ):
        a.bar(
            np.arange(2) + (j - 1) * 0.23,
            [result["models"][m]["performance"][key] for m in ["G0", "G1"]],
            width=0.22,
            color=color,
            label=label,
        )
    a.axhline(2, color="#555", ls="--", lw=1, label="2× body threshold")
    a.set(xticks=range(2), xticklabels=["G0", "G1"], ylabel="Full M0 / surrogate")
    a.legend(fontsize=7.5)
    grid(a)
    a = ax[1, 2]
    title(a, "f  Separate evidence gates")
    a.axis("off")
    tab = a.table(
        cellText=[["G0", "FAIL", "1.98×", "FAIL"], ["G1", "PASS", "1.82×", "FAIL"]],
        colLabels=["Model", "Accuracy", "Body speed", "≥2× + accuracy"],
        cellLoc="center",
        loc="upper center",
        colWidths=[0.18, 0.24, 0.25, 0.34],
    )
    tab.auto_set_font_size(False)
    tab.set_fontsize(9)
    tab.scale(1, 1.6)
    a.text(
        0.02,
        0.54,
        "G1 maximum Z error: 3.60 nm\nG1 largest endpoint 3D error: 3.88 nm\n\nNo empirical gain fitted.\nNo additional batch campaign.\ntraining_ready = false",
        transform=a.transAxes,
        fontsize=10,
        va="top",
        linespacing=1.6,
    )
    save(fig, NAMES[2])

    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(
        2, 2, left=0.07, right=0.96, bottom=0.13, top=0.88, hspace=0.53, wspace=0.24, height_ratios=[1, 1]
    )
    a = fig.add_subplot(gs[0, :])
    a.axis("off")
    fig.suptitle("Explicit surrogate scope: admission is not a morphology guarantee", fontsize=15, y=0.965)
    title(a, "a  Actual saved states: static admission / rejection (no rejected-case proxy rollout)")
    data = []
    for row in checks["scope"]:
        data.append(
            [
                row["case"].replace("F10A_", "F10A "),
                f"{row['minimum_wall_gap_m'] * 1e6:.3f}",
                f"{row['maximum_shear_s_inv']:.3f}",
                f"{row['maximum_pressure_Pa'] / 1000:.3f}",
                f"{row['initial_rg_m'] * 1e6:.2f}",
                f"{row['largest_fraction']:.2f}",
                "ADMIT" if not row["reasons"] else "REJECT",
            ]
        )
    tab = a.table(
        cellText=data,
        colLabels=[
            "Saved state",
            "Min wall gap (µm)",
            "Max shear (s⁻¹)",
            "Peak |P| (kPa)",
            "Rg (µm)",
            "Largest fraction",
            "Decision",
        ],
        cellLoc="center",
        bbox=[0, -0.04, 1, 0.94],
        colWidths=[0.21, 0.17, 0.16, 0.14, 0.1, 0.13, 0.1],
    )
    tab.auto_set_font_size(False)
    tab.set_fontsize(9)
    for (i, j), cell in tab.get_celld().items():
        cell.set_edgecolor("#D6DFE5")
        if i == 0:
            cell.set_facecolor("#E8EEF3")
            cell.set_text_props(weight="bold")
        elif j == 6:
            cell.set_facecolor("#E2F0EA" if i <= 3 else "#F4E7E5")
    a.text(
        0.0,
        -0.18,
        "RPY row uses its actual flow-OFF setting (zero applied shear); the current proxy requires the validated CFD-ON setting.",
        fontsize=8.5,
        transform=a.transAxes,
    )
    a = fig.add_subplot(gs[1, 0])
    title(a, "b  F5 actual OFF release: a frozen template would be wrong")
    with h5py.File(ROOT / "results/fast_pzt/f5/science_v1/main/trajectory.h5") as h:
        t = h["time_s"][:]
        sel = t >= 0.3 - 1e-12
        rg = h["metrics/global_rg_m"][:, 0][sel]
        fraction = h["metrics/largest_size"][:, 0][sel] / 32
    a.plot((t[sel] - 0.3) * 1e3, rg / rg[0], color=COLORS[0], label="Actual Rg / initial Rg")
    a.set(xlabel="Time after sound OFF (ms)", ylabel="Rg / Rg at switch", ylim=(0.95, 2.8))
    twin = a.twinx()
    twin.plot((t[sel] - 0.3) * 1e3, fraction, color=COLORS[2], label="Actual largest fraction")
    twin.set(ylabel="Largest connected fraction", ylim=(0, 1.17))
    a.legend(loc="upper left", fontsize=8)
    twin.legend(loc="upper right", fontsize=8)
    grid(a)
    a = fig.add_subplot(gs[1, 1])
    title(a, "c  F10A actual internal evolution, not proxy output")
    for key, color in CONTROL.items():
        with h5py.File(F10A / "extended_5_v1" / key / "trajectory.h5") as h:
            rg = h["metrics/global_rg_m"][:, 0]
        a.plot(times, 100 * (rg / rg[0] - 1), color=color, label=key.replace("_1000", ""))
    a.set(xlabel="Time in validated interior window (ms)", ylabel="Full-reference Rg change (%)")
    a.legend(fontsize=8)
    grid(a)
    fig.text(
        0.5,
        0.063,
        "Required: M0, 32 equal radii, registered mature template + controls, ≤10 ms, gap ≥100 µm, shear ≤2 s⁻¹, |P| ≤1.5 kPa; no exits / morphology request.",
        ha="center",
        fontsize=8.5,
    )
    save(fig, NAMES[3])
    write(out / "FIGURE_SOURCES.json", sources)
    print("Four PNG/PDF figure pairs saved")


if __name__ == "__main__":
    main()
