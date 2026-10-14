"""Data-only English publication figures. No manufactured trajectories."""

import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.colors import ListedColormap
from .common import BASE, ROOT, CONFIG, read, write, sha, identity

PURPLE = "#70559c"
TEAL = "#258e95"
RED = "#b45a60"
GREY = "#697580"
COLORS = [GREY, TEAL, PURPLE]
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelcolor": "#262b32",
        "text.color": "#262b32",
        "xtick.color": "#262b32",
        "ytick.color": "#262b32",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
        "axes.titlesize": 11,
    }
)
OUT = BASE / "figures"
AN = BASE / "analysis"
outputs = []


def finish(fig, name, title, footer):
    footer += "\nEXPLORATORY · HOMOGENEOUS_WATER_CONTROL_PROXY · PHYSICS_UNCALIBRATED · training_ready=false"
    fig.suptitle(title, fontsize=15, y=0.99)
    fig.text(0.5, 0.015, footer, ha="center", va="bottom", fontsize=8, color="#444950")
    fig.tight_layout(rect=(0, 0.055, 1, 0.965), h_pad=2.7, w_pad=2.5)
    for ext in ["png", "pdf"]:
        p = OUT / (name + "." + ext)
        fig.savefig(p, dpi=220)
        outputs.append(p)
    plt.close(fig)


def rows(path):
    with h5py.File(path / "trajectory.h5") as h:
        return h["step_metrics"][:]


def fig1():
    d = read(AN / "TRUE_GAP_FAILURE_DIAGNOSIS.json")
    w = d["worst_pair"]
    f = read(BASE / "local/FOUNDATIONS_v2.json")["frozen_failure"]
    p = read(CONFIG)
    fig, ax = plt.subplots(2, 2, figsize=(11.4, 8.6))
    a, b, c, e = ax.ravel()
    x = np.array(d["all_positions_before_m"])
    y = np.array(d["all_positions_proposal_m"])
    ii, jj = w["pair"]
    n = np.array(w["old_normal"])
    dr = (y[jj] - x[jj]) - (y[ii] - x[ii])
    tan = dr - n * (n @ dr)
    tan /= np.linalg.norm(tan)
    basis = np.array([n, tan])
    mid = (x[ii] + x[jj]) / 2
    for k in [ii, jj]:
        old = (x[k] - mid) @ basis.T * 1e6
        new = (y[k] - mid) @ basis.T * 1e6
        a.add_patch(Circle(old, p["radius_m"] * 1e6, facecolor="#ebedf0", edgecolor=GREY, lw=1.5))
        a.add_patch(Circle(new, p["radius_m"] * 1e6, fill=False, edgecolor=PURPLE, lw=1))
        a.text(*old, f"ID {k}", ha="center", va="center")
    a.set(
        xlim=(-4, 4),
        ylim=(-2.1, 2.1),
        xlabel="Old normal coordinate (µm)",
        ylabel="Tangential coordinate (µm)",
        title="a  Actual rejected pair at 0.464 ms",
        aspect="equal",
    )
    a.text(
        0.5, 0.02, "True radii; displacement is not magnified", transform=a.transAxes, ha="center", fontsize=8
    )
    values = [w["linear_gap_m"], w["actual_gap_m"], w["second_order_m"], w["weighted_true_comp_m"]]
    b.bar(
        ["Linear\ngap", "True\ngap", "Second-order\ncurvature", "Weighted\ncomplementarity"],
        np.array(values) * 1e9,
        color=[GREY, PURPLE, TEAL, RED],
    )
    b.axhline(0.01, color=RED, ls="--", lw=1, label="Weighted threshold: 0.01 nm")
    b.set(ylabel="Gap / equivalent gap (nm)", title="b  Curvature quantitatively explains the gap")
    b.legend(fontsize=8, loc="upper left")
    b.set_ylim(-0.0004, 0.018)
    c.axis("off")
    c.set_title("c  Actual kinematics and reaction", loc="left")
    c.text(
        0,
        0.95,
        f"Tangential relative displacement: {w['relative_tangential_displacement_m'] * 1e9:.6f} nm\nNormal relative displacement: {w['relative_normal_displacement_m'] * 1e9:.6f} nm\nNormal rotation: {w['normal_rotation_deg']:.6f}°\nPair reaction: {w['reaction_N'] * 1e9:.6f} nN\nNormalization force: {d['force_scale_N'] * 1e9:.6f} nN\n\nExact curvature identity error < {d['max_exact_identity_error_m'] * 1e9:.2g} nm\nPair mapping and SI-unit residual reproduced.\nNo new 2.1R graph edge at this rejected step.",
        va="top",
        linespacing=1.7,
    )
    hist = f["iterations"]
    q = np.arange(1, len(hist) + 1)
    force_tol = p["force_residual_atol_N"] + p["force_residual_rtol"] * f["force_scale_N"]
    e.semilogy(
        q,
        [v["true_force_residual_N"] / force_tol for v in hist],
        "o-",
        color=RED,
        label="True-normal force / limit",
    )
    e.semilogy(
        q,
        [v["true_complementarity_m"] / p["actual_geometry_complementarity_m"] for v in hist],
        "s-",
        color=PURPLE,
        label="True complementarity / limit",
    )
    e.axhline(1, color=GREY, ls="--")
    e.set(
        xticks=q,
        xlabel="Local QP within the same physical step",
        ylabel="Residual / unchanged limit",
        title="d  Frozen proposal: all gates pass at QP 4",
    )
    e.legend(fontsize=8)
    finish(
        fig,
        "Fig_F10Q_1_true_gap_mechanism",
        "F10Q | Finite-sphere geometry and contact consistency",
        "Actual F10P rejected proposal and F10Q frozen-state replay. Frozen F/A; repeated QPs do not advance physical time.",
    )


def fig2():
    data = read(AN / "RUN_ANALYSIS.json")
    new = data["STRONG_DT1"]
    oldp = ROOT / "results/straight_tube_exp1/f10p/validation_v1/runs/STRONG_DT1"
    oldb = ROOT / "results/straight_tube_exp1/f10o/validation_v1/runs/STRONG_DT1"
    newp = BASE / "a_strong_v1/runs/STRONG_DT1"
    pr = read(oldp / "RESULT.json")
    br = read(oldb / "RESULT.json")
    fig, ax = plt.subplots(2, 3, figsize=(14, 8.8))
    a, b, c, d, e, f = ax.ravel()
    names = ["B2", "Local hard", "Sequential QP"]
    vals = [br["completed_s"], pr["completed_s"], new["completed_s"]]
    a.barh(names, np.array(vals) * 1e3, color=[GREY, TEAL, PURPLE])
    a.axvline(2, color=RED, ls="--")
    a.set(
        xlim=(0, 2.2),
        xlabel="Accepted physical time (ms)",
        title="a  Same strong control; 2 ms not completed",
    )
    for k, v in enumerate(vals):
        a.text(v * 1e3 + 0.04, k, f"{v * 1e3:.3f}", va="center")
    for folder, col, label in [(oldp, TEAL, "Local hard"), (newp, PURPLE, "Sequential QP")]:
        rr = rows(folder)
        tt = (np.arange(len(rr)) + 1) * 1e-3
        legal = rr[:, 0].astype(bool)
        b.semilogy(tt[legal], np.maximum(rr[legal, 17] * 1e9, 1e-15), color=col, label=label)
        if (~legal).any():
            b.plot(tt[~legal], rr[~legal, 17] * 1e9, "x", color=col, ms=8)
    b.axhline(0.01, color=RED, ls="--")
    b.set(
        xlabel="Attempted physical time (ms)",
        ylabel="True weighted complementarity (nm)",
        title="b  Accepted prefix and rejected proposal",
    )
    b.legend(fontsize=8)
    with h5py.File(newp / "trajectory.h5") as h:
        rr = h["step_metrics"][:]
        extra = h["extra_metrics"][:]
        reaction = abs(h["pair_force_N"][:]).max(1) * 1e9
    tt = (np.arange(len(rr)) + 1) * 1e-3
    legal = rr[:, 0].astype(bool)
    c.plot(tt[legal], reaction[legal], color=PURPLE, label="Maximum pair reaction")
    c.set(
        xlabel="Physical time (ms)",
        ylabel="Pair reaction (nN)",
        title="c  Real reactions and active contacts",
    )
    c2 = c.twinx()
    c2.plot(tt[legal], extra[legal, 0], color=TEAL, lw=1)
    c2.set_ylabel("Active contacts", color=TEAL)
    c2.spines["top"].set_visible(False)
    for folder, col, label in [
        (oldb, GREY, "B2"),
        (oldp, TEAL, "Local hard"),
        (newp, PURPLE, "Sequential QP"),
    ]:
        rr = rows(folder)
        legal = rr[:, 0].astype(bool)
        t = (np.arange(len(rr)) + 1) * 1e-3
        d.plot(t[legal], rr[legal, 2] * 1e9, color=col, label=label)
    d.axhline(17.8, color=RED, ls="--")
    d.set(
        xlabel="Physical time (ms)",
        ylabel="Geometric overlap (nm)",
        title="d  Unchanged 17.8 nm overlap budget",
    )
    d.legend(fontsize=8)
    e.bar(
        names,
        np.array([br["mean_step_s"], pr["mean_step_s"], new["mean_step_s"]]) * 1e3,
        color=[GREY, TEAL, PURPLE],
    )
    e.set(ylabel="Measured mean step cost (ms)", title="e  Actual CPU/GPU end-to-end step cost")
    e.tick_params(axis="x", rotation=15)
    f.axis("off")
    f.set_title("f  Limits and rejected states", loc="left")
    f.text(
        0,
        0.97,
        f"Sequential QP: 1–{new['max_local_qps']} local solves per step\nStrong rejection: active-set cycle at 0.884 ms\nNo trajectory point beyond 0.883 ms accepted.\n\nNew strong dt = 0.5 µs: NOT_RUN\nReason: full dt = 1 µs run failed.\nOld local hard dt = 0.5 µs: 0.5 ms prefix only.\n\nSame-model navigation dt/cache checks: Fig. 4\nB2 complementarity: NOT_APPLICABLE\nNo spring or PBD after the hard QP.",
        va="top",
        fontsize=9,
        linespacing=1.65,
    )
    finish(
        fig,
        "Fig_F10Q_2_sequential_qp_admission",
        "F10Q | Strong-control contact admission and stopping evidence",
        "EXPLORATORY · No tolerance relaxation. Historical B2/local hard data retained; sequential trajectory starts from the original initial state.",
    )


def fig3():
    data = read(AN / "THREE_MODEL_COMPARISON.json")
    ts = np.array(data["common_saved_time_s"]) * 1e3
    fig, ax = plt.subplots(2, 3, figsize=(14, 8.5))
    axs = ax.ravel()
    for axis in range(3):
        a = axs[axis]
        for (name, m), col, marker in zip(data["models"].items(), COLORS, ["o", "s", "^"]):
            a.plot(
                ts,
                np.array(m["COM_m"])[:, axis] * 1e6,
                marker=marker,
                ls="-",
                mfc="none",
                color=col,
                label=name,
                ms=8 - axis,
            )
        a.set(
            xlabel="Exact saved time (ms)",
            ylabel=f"COM {'XYZ'[axis]} (µm)",
            title=f"{'abc'[axis]}  Fixed-ID centroid: {'XYZ'[axis]}",
        )
        a.ticklabel_format(axis="y", useOffset=False)
    axs[0].legend(fontsize=8)
    axs[2].text(0.04, 0.08, "Round-off scale; no resolved Z motion", transform=axs[2].transAxes, fontsize=8)
    for (name, m), col in zip(data["models"].items(), COLORS):
        axs[3].plot(ts, np.array(m["Rg_m"]) * 1e6, "o-", color=col, label=name)
    axs[3].set(
        xlabel="Exact saved time (ms)",
        ylabel="Radius of gyration (µm)",
        title="d  Shape differs despite similar COM",
    )
    for j, ((name, m), col) in enumerate(zip(data["models"].items(), COLORS)):
        q = np.array(m["Q_m2"])
        dq = np.diag(q[-1] - q[0]) * 1e12
        axs[4].bar(np.arange(3) + (0.22 * (j - 1)), dq, width=0.22, label=name, color=col)
    axs[4].set(
        xticks=np.arange(3),
        xticklabels=["ΔQxx", "ΔQyy", "ΔQzz"],
        ylabel="Shape-tensor change (µm²)",
        title="e  Centered shape tensor at 0.1 ms",
    )
    axs[4].legend(fontsize=8)
    names = list(data["contrasts"])
    xx = np.arange(3)
    axs[5].bar(
        xx - 0.18,
        [data["contrasts"][k]["COM_m"] * 1e9 for k in names],
        width=0.36,
        label="COM difference",
        color=TEAL,
    )
    axs[5].bar(
        xx + 0.18,
        [data["contrasts"][k]["fixed_id_m"] * 1e9 for k in names],
        width=0.36,
        label="Maximum fixed-ID difference",
        color=PURPLE,
    )
    axs[5].set(
        xticks=xx,
        xticklabels=[k.replace("_VS_", " / ") for k in names],
        yscale="log",
        ylabel="Cross-model difference (nm)",
        title="f  Model contrast, not dt error",
    )
    axs[5].set_ylim(1e-6, 1e3)
    axs[5].legend(fontsize=8)
    finish(
        fig,
        "Fig_F10Q_3_group_motion_model_comparison",
        "F10Q | Historical common-time group motion: M0/PBD, B2 and hard spheres",
        "Different contact/mobility approximations. Only 0 and 0.1 ms are jointly stored; lines connect samples. No interpolation or longer-window claim.",
    )


def fig4():
    data = read(AN / "RUN_ANALYSIS.json")
    gate = read(AN / "FINAL_GATE.json")
    p = read(CONFIG)
    fig, ax = plt.subplots(2, 3, figsize=(14.5, 9))
    a, b, c, d, e, f = ax.ravel()
    for name, r in data.items():
        if "_P" in name:
            m = np.load(AN / (name + ".npz"))
            a.plot(m["centroid_m"][:, 0] * 1e6, m["centroid_m"][:, 1] * 1e6, color="#ccd1d9", lw=0.7)
    for axis, col in [("X", PURPLE), ("Y", TEAL)]:
        tag = "NAV_" + axis
        key = tag + "_DT1"
        target = np.array(p["XY_targets_m"][axis]) * 1e6
        a.add_patch(Circle(target[:2], 0.075, fill=False, ls="--", edgecolor=col, alpha=0.6))
        a.plot(*target[:2], marker="*", color=col, ms=12, label=axis + " target")
        if key in data:
            m = np.load(AN / (key + ".npz"))
            rr = data[key]
            cent = m["centroid_m"] * 1e6
            a.plot(
                cent[:, 0], cent[:, 1], color=col, label=f"{axis} control ({rr['completed_s'] * 1e3:.3f} ms)"
            )
            a.plot(*cent[-1, :2], "o", color=col)
            error = np.linalg.norm(cent - target, axis=1)
            b.plot(m["time_s"] * 1e3, error, color=col, label=axis + " error")
            c.plot(m["time_s"] * 1e3, (m["rg_m"] - m["rg_m"][0]) * 1e9, color=col, label=axis + " ΔRg")
    a.plot(0, 0, "s", color=GREY, label="OFF")
    a.set(
        xlim=(-0.09, 0.3),
        ylim=(-0.09, 0.3),
        xlabel="COM X (µm)",
        ylabel="COM Y (µm)",
        title="a  Original XY targets and actual trajectories",
        aspect="equal",
    )
    a.legend(fontsize=7, loc="upper right")
    b.axhline(0.075, color=RED, ls="--", label="Target tolerance")
    b.set(
        xlabel="Accepted physical time (ms)",
        ylabel="Distance to assigned target (µm)",
        title="b  Full-window target error",
    )
    b.legend(fontsize=8)
    for suffix, style in [("DT05", "--"), ("FINE", ":")]:
        path = AN / ("NAV_Y_" + suffix + ".npz")
        if path.exists():
            m = np.load(path)
            c.plot(
                m["time_s"] * 1e3,
                (m["rg_m"] - m["rg_m"][0]) * 1e9,
                ls=style,
                color=TEAL,
                alpha=0.7,
                label="Y " + ("half dt" if suffix == "DT05" else "fine cache"),
            )
    c.set(
        xlabel="Accepted physical time (ms)",
        ylabel="Centered ΔRg (nm)",
        title="c  Shape evolution, separate from translation",
    )
    c.legend(fontsize=8)
    names = ["OFF", "NAV_X", "NAV_Y"]
    keys = [
        "BASIC_NUMERICS_PASS",
        "GROUP_POSITION_NUMERICALLY_ADMITTED",
        "GROUP_SHAPE_NUMERICALLY_ADMITTED",
        "GROUP_TOPOLOGY_NUMERICALLY_ADMITTED",
        "PARTICLE_ID_NUMERICALLY_ADMITTED",
    ]
    matrix = np.full((3, 6), 2)
    texts = []
    for i, name in enumerate(names):
        g = gate["groups"].get(name, {})
        rr = []
        for j, key in enumerate(keys):
            status = g.get(key, "NOT_RUN")
            matrix[i, j] = 0 if status == "PASS" else (1 if status == "FAIL" else 2)
            rr.append("PASS" if status == "PASS" else ("FAIL" if status == "FAIL" else "N/R"))
        reached = g.get("TARGET_REACHED", {}).get(name[-1], False) if name != "OFF" else False
        matrix[i, 5] = 0 if reached else 1
        rr.append("YES" if reached else "NO")
        texts.append(rr)
    d.imshow(matrix, cmap=ListedColormap(["#e0eeea", "#f2dfe0", "#eceef1"]), vmin=0, vmax=2, aspect="auto")
    d.set(
        xticks=np.arange(6),
        xticklabels=["Basic", "COM", "Shape", "Graph", "ID", "Target"],
        yticks=np.arange(3),
        yticklabels=["OFF", "X", "Y"],
        title="d  Independent layered admission (2 ms)",
    )
    for i in range(3):
        for j in range(6):
            d.text(j, i, texts[i][j], ha="center", va="center", fontsize=8)
    e.set_title("e  Same-model numerical differences", loc="left")
    e.axis("off")
    lines = []
    comps = read(AN / "COMPARISONS.json")
    for axis in ["X", "Y"]:
        for suf in ["DT05", "FINE"]:
            v = comps.get("NAV_" + axis + "_" + suf)
            lines.append(
                f"{axis} {'dt' if suf == 'DT05' else 'cache'}: "
                + (
                    f"COM {v['COM_m'] * 1e9:.3g} nm; ID {v['fixed_id_m'] * 1e9:.3g} nm\n          ΔRg {v['Rg_m'] * 1e9:.3g} nm; ΔQ/Q {v['Q_relative']:.3g}"
                    if v
                    else "NOT_RUN"
                )
            )
    e.text(0, 0.96, "\n\n".join(lines), va="top", fontsize=9)
    f.set_title("f  Scope, failures and actual computation", loc="left")
    f.axis("off")
    cost = read(AN / "COSTS.json")
    lines = [
        "8 / 8 new controls: legal 0.5 ms screenings",
        f"Route A: {cost['routes_wall_s']['A']:.2f} s; route B: {cost['routes_wall_s']['B']:.2f} s",
    ]
    for axis in ["X", "Y"]:
        r = data.get("NAV_" + axis + "_DT1", {})
        g = gate["groups"].get("NAV_" + axis, {})
        err = g.get("target_error_m", {}).get(axis)
        lines.append(
            f"{axis}: {r.get('status', 'NOT_RUN')}; error " + (f"{err * 1e6:.6f} µm" if err else "NOT_RUN")
        )
        if r.get("reasons"):
            lines.append("  " + ", ".join(r["reasons"]))
    lines += [
        "5 / 10 / 30 ms: NOT_RUN",
        "Exploratory group data; contact-model dependent",
        "Experimental calibration: NOT_VERIFIED",
        "TRAINING_READY = false",
    ]
    f.text(0, 0.96, "\n\n".join(lines), va="top", fontsize=8.5)
    finish(
        fig,
        "Fig_F10Q_4_xy_navigation",
        "F10Q | Independent XY navigation with four PZTs and zero flow",
        "2.5 MHz · Homogeneous-water control proxy · PHYSICS_UNCALIBRATED · Shape/ID gates remain independent of COM admission.",
    )


def main():
    OUT.mkdir(exist_ok=True)
    fig1()
    fig2()
    fig3()
    fig4()
    inputs = (
        list(AN.glob("*.json"))
        + list(AN.glob("*.npz"))
        + [BASE / "local/FOUNDATIONS_v2.json", BASE / "CANDIDATES.json"]
    )
    write(
        OUT / "FIGURE_PROVENANCE.json",
        dict(
            identity=identity(),
            script_sha256=sha(__file__),
            input_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs},
            figures_sha256={p.name: sha(p) for p in outputs},
            style="White publication backgrounds, English text, original scientific metrics; no AI imagery or synthetic trajectories.",
        ),
    )


if __name__ == "__main__":
    main()
