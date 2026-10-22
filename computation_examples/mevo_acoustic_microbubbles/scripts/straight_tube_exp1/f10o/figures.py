"""Four publication figures from immutable actual F10N/F10O trajectories."""

import numpy as np
import textwrap
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from .common import BASE, OLD, ROOT, read, write, sha, identity
from .analyze import load, metrics

OUT = BASE / "figures"
PURPLE = "#665298"
BLUE = "#2878a5"
TEAL = "#247f76"
RED = "#aa3949"
GREY = "#737984"
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "axes.titleweight": "bold",
        "axes.edgecolor": "#343b44",
        "text.color": "#222832",
        "axes.labelcolor": "#222832",
        "xtick.color": "#343b44",
        "ytick.color": "#343b44",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "legend.frameon": False,
        "savefig.facecolor": "white",
    }
)
SOURCES = {}


def used(path):
    SOURCES[str(path.relative_to(ROOT))] = sha(path)
    return path


def trajectory(name, old=False):
    path = (OLD / "campaign_v1/runs" if old else BASE / "validation_v1/runs") / name
    used(path / "trajectory.h5")
    return load(path)


def panel(ax, label, title):
    cols = ax.get_subplotspec().get_gridspec().ncols
    ax.set_title(label + "  " + textwrap.fill(title, width=34 if cols == 3 else 47), loc="left", pad=12)
    ax.grid(alpha=0.15)


def footer(fig):
    fig.text(
        0.5,
        0.015,
        "EXPLORATORY  |  HOMOGENEOUS_WATER_CONTROL_PROXY  |  PHYSICS_UNCALIBRATED  |  training_ready=false",
        ha="center",
        fontsize=9,
        color=GREY,
    )


def save(fig, name, title, records):
    fig.suptitle(title, fontsize=17, fontweight="bold", y=0.99)
    footer(fig)
    fig.subplots_adjust(left=0.135, right=0.96, bottom=0.13, top=0.85, hspace=0.70, wspace=0.36)
    for ext in ["png", "pdf"]:
        fig.savefig(OUT / (name + "." + ext), dpi=260, bbox_inches="tight", pad_inches=0.18)
    plt.close(fig)
    records.append(
        dict(
            figure=name,
            title=title,
            file_sha256={ext: sha(OUT / (name + "." + ext)) for ext in ["png", "pdf"]},
        )
    )


def main():
    OUT.mkdir(exist_ok=True)
    g = read(used(BASE / "analysis/FINAL_GATE.json"))
    runs = read(used(BASE / "analysis/RUN_ANALYSIS.json"))
    checks = read(used(BASE / "pilot_v1/CHECKS.json"))
    comparisons = read(used(BASE / "analysis/COMPARISONS.json"))
    used(BASE / "PREREGISTRATION.json")
    used(BASE / "analysis/COSTS.json")
    old_quality = read(used(OLD / "campaign_v1/runs/OPTIMIZED_REPLAY/RESULT.json"))["quality"]
    old_sensitivity = read(used(OLD / "analysis/NUMERICAL_SENSITIVITY.json"))["navigation_dt"]
    records = []
    fig, ax = plt.subplots(2, 2, figsize=(12.8, 8.8))
    a = ax[0, 0]
    panel(a, "a", "Strong control: actual initial force input")
    ev = checks["evidence"]["OPTIMIZED_REPLAY"]
    keys = [
        "primary_force_N",
        "secondary_force_N",
        "no_feedback_regularized_secondary_N",
        "nearfield_feedback_increment_N",
    ]
    vals = [np.linalg.norm(ev[k], axis=1).max() * 1e12 for k in keys]
    a.bar(range(4), vals, color=[BLUE, PURPLE, GREY, TEAL], width=0.62)
    a.set_xticks(range(4), ["Primary", "F2 secondary", "No-feedback\nsecondary", "Feedback\nincrement"])
    a.set_yscale("log")
    a.set_ylim(top=max(vals) * 4)
    a.set_ylabel("Maximum particle force magnitude (pN)")
    a.text(
        0.02,
        0.96,
        "Same fixed IDs and full-precision control",
        transform=a.transAxes,
        va="top",
        fontsize=9,
    )
    a = ax[0, 1]
    panel(a, "b", "Strong B2: geometric overlap rejection")
    for tag, c, ls in [("STRONG_DT1", PURPLE, "-"), ("STRONG_DT05", TEAL, "--"), ("STRONG_FINE", BLUE, ":")]:
        path = used(BASE / "validation_v1/runs" / tag / "trajectory.h5")
        with h5py.File(path) as h:
            m = h["step_metrics"][:]
        r = runs["validation_v1/" + tag]
        tt = np.arange(1, len(m) + 1) * r["dt_s"] * 1e3
        a.plot(tt, m[:, 2] * 1e9, label=tag.replace("STRONG_", ""), c=c, ls=ls, lw=2)
        a.scatter(tt[-1], m[-1, 2] * 1e9, c=RED, marker="x", s=65, zorder=6)
    a.axhline(17.8, color=RED, ls="--", lw=1.2, label="Unchanged 0.01R limit")
    a.set(xlabel="Physical time (ms)", ylabel="Maximum geometric overlap (nm)")
    a.legend(fontsize=9, loc="lower right")
    a.text(0.03, 0.95, "All three proposals rejected at 0.175 ms", transform=a.transAxes, va="top", color=RED)
    a = ax[1, 0]
    panel(a, "c", "Quality checks are independent")
    f = g["strong_failures"][0]
    values = [f["residual_m"] / 1e-13, f["overlap_m"] / 1.78e-8, f["max_iteration_m"] / 1.78e-7]
    a.barh(["Newton residual", "Geometric overlap", "Newton correction"], values, color=[TEAL, RED, TEAL])
    a.axvline(1, c=GREY, ls="--")
    a.set(xlim=(0, 1.45), xlabel="Measured quantity / unchanged budget")
    for j, v in enumerate(values):
        a.text(v + 0.025, j, f"{v:.5g}", va="center", fontsize=10)
    a.text(
        0.40,
        0.89,
        "Wall correction = 0\nSealed-X / top-scope checks active",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.95),
        transform=a.transAxes,
        fontsize=9,
    )
    a = ax[1, 1]
    a.axis("off")
    panel(a, "d", "Reused models, different contact approximations")
    text = (
        "M0 / PBD (historical F10N)\n"
        "  Two local lubrication iterations; 12 geometric projection passes.\n"
        f"  Strong 2 ms same-ID dt difference: {old_sensitivity['max_fixed_id_position_difference_m'] * 1e6:.3f} µm.\n"
        f"  Cumulative PBD correction / predicted path: {old_quality['correction_fraction']:.3f}.\n\n"
        "CPU B2 (unchanged F10H / F10K)\n"
        "  Solve the same pair-normal resistance matrix.\n"
        "  Implicit unilateral spring on actual final geometry.\n"
        "  k = 0.04641319518378658 N/m; residual tolerance 10⁻¹³ m.\n"
        "  No subsequent PBD; any wall correction rejects this interior test.\n\n"
        "A converged soft-contact equation can still exceed the overlap limit.\n"
        "M0–B2 trajectory differences are model contrasts,\n"
        "not pure time-integration error."
    )
    a.text(
        0,
        0.98,
        "\n".join(textwrap.fill(line, 64) if line else "" for line in text.splitlines()),
        va="top",
        fontsize=9.5,
        linespacing=1.25,
    )
    save(fig, "Fig_F10O_1_contact_mechanism", "F10O | Strong-control contact audit", records)

    fig, ax = plt.subplots(2, 3, figsize=(15, 8.8))
    a = ax[0, 0]
    panel(a, "a", "Same-ID differences: strong valid prefix")
    for key, c, label in [("STRONG_DT05", PURPLE, "1 vs 0.5 µs"), ("STRONG_FINE", TEAL, "3 vs 1.5 µm cache")]:
        data = np.load(used(BASE / "analysis" / f"{key}.npz"))
        a.semilogy(data["time_s"][1:] * 1e3, data["fixed_id_m"][1:] * 1e6, c=c, label=label)
    a.axhline(0.05, c=RED, ls="--", label="0.05 µm gate")
    a.set(xlabel="Time (ms)", ylabel="Maximum same-ID difference (µm)")
    a.legend(fontsize=8, loc="center right")
    a.text(
        0.04,
        0.96,
        "2 ms admission: FAIL (trajectory rejected)",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.9),
        transform=a.transAxes,
        va="top",
        color=RED,
        fontsize=9,
    )
    a = ax[0, 1]
    panel(a, "b", "Weak control: complete 2 ms sensitivity")
    for key, c, label in [("WEAK_DT05", PURPLE, "Time step"), ("WEAK_FINE", TEAL, "Field cache")]:
        data = np.load(used(BASE / "analysis" / f"{key}.npz"))
        a.plot(data["time_s"] * 1e3, data["fixed_id_m"] * 1e9, c=c, label=label)
    a.set(xlabel="Time (ms)", ylabel="Maximum same-ID difference (nm)", ylim=(-0.0001, 0.0062))
    a.legend(fontsize=9)
    a.text(0.03, 0.95, "Both pass the original 50 nm tolerance", transform=a.transAxes, va="top", fontsize=9)
    a = ax[0, 2]
    panel(a, "c", "COM and internal-shape sensitivity")
    keys = ["fixed_id_m", "COM_m", "kabsch_m", "Rg_m"]
    for i, (case, c) in enumerate([("STRONG_DT05", PURPLE), ("WEAK_DT05", TEAL)]):
        values = [max(comparisons[case][k] * 1e9, 1e-10) for k in keys]
        a.bar(
            np.arange(4) + (i - 0.5) * 0.34,
            values,
            width=0.32,
            color=c,
            label="Strong prefix" if i == 0 else "Weak 2 ms",
        )
    a.set_xticks(range(4), ["Same ID", "COM", "Kabsch", "Rg"])
    a.set_yscale("log")
    a.set_ylabel("Time-step difference (nm)")
    a.set_ylim(top=1)
    a.legend(fontsize=8)
    a.text(
        0.01,
        -0.24,
        "Tensor difference (µm²):\n"
        + f"strong {comparisons['STRONG_DT05']['tensor_m2'] * 1e12:.3g}; weak {comparisons['WEAK_DT05']['tensor_m2'] * 1e12:.3g}",
        transform=a.transAxes,
        fontsize=9,
    )
    a = ax[1, 0]
    panel(a, "d", "Weak control: Rg changes sign across models")
    for name, old, c, label in [
        ("SHAPE_REPLAY", True, GREY, "M0 / PBD"),
        ("WEAK_DT1", False, PURPLE, "CPU B2"),
    ]:
        t, x = trajectory(name, old)
        m = metrics(t, x)
        a.plot(t * 1e3, (m["rg_m"] - m["rg_m"][0]) * 1e9, c=c, label=label)
    a.axhline(0, c="#999999", lw=0.6)
    a.set(xlabel="Time (ms)", ylabel="Change in Rg (nm)")
    a.legend()
    a = ax[1, 1]
    panel(a, "e", "Weak control: small X/Y shape response")
    for name, old, c, label in [
        ("SHAPE_REPLAY", True, GREY, "M0 / PBD"),
        ("WEAK_DT1", False, PURPLE, "CPU B2"),
    ]:
        t, x = trajectory(name, old)
        m = metrics(t, x)
        a.plot(t * 1e3, (m["directional_xy"] - 1) * 100, c=c, label=label)
    a.set(xlabel="Time (ms)", ylabel="Change in X/Y ratio (%)")
    a.legend()
    a.text(
        0.03,
        0.96,
        "Registered shape target: +25% (unmet)",
        transform=a.transAxes,
        va="top",
        fontsize=9,
        color=RED,
    )
    a = ax[1, 2]
    panel(a, "f", "Fixed-ID weak B2 geometry, 2 ms")
    _, x = trajectory("WEAK_DT1")
    xx = x[-1] * 1e6
    # Oblique orthographic projection; unchanged coordinates, no visual deformation.
    u = np.array([1, -1, 0]) / np.sqrt(2)
    v = np.array([0.35, 0.35, 1])
    v /= np.linalg.norm(v)
    for layer, c in [(xx[:, 2] < 0, BLUE), (xx[:, 2] >= 0, PURPLE)]:
        a.scatter(xx[layer] @ u, xx[layer] @ v, s=38, c=c, edgecolor="white", linewidth=0.4)
    a.set_aspect("equal")
    a.set(xlabel="Oblique coordinate u (µm)", ylabel="Oblique coordinate v (µm)", ylim=(-5.7, 7.4))
    a.text(
        0.02,
        0.98,
        "32 centres; one component\nNo new graph event",
        transform=a.transAxes,
        va="top",
        fontsize=9,
    )
    save(
        fig,
        "Fig_F10O_2_particle_convergence",
        "F10O | Fixed-ID convergence and model-dependent morphology",
        records,
    )

    fig, ax = plt.subplots(2, 2, figsize=(12.8, 8.6))
    a = ax[0, 0]
    panel(a, "a", "XY targets locked before any new trajectory")
    for label, target, c in [("X", g["XY_targets_m"]["X"], PURPLE), ("Y", g["XY_targets_m"]["Y"], TEAL)]:
        q = np.array(target) * 1e6
        a.scatter(*q[:2], marker="*", s=150, c=c, zorder=5, label=label + " target")
        a.add_patch(Circle(q[:2], g["XY_tolerance_m"] * 1e6, fill=False, ec=c, ls="--"))
    a.scatter(0, 0, c="#222832", s=45, label="Initial COM / OFF")
    a.set(xlabel="X (µm)", ylabel="Y (µm)", xlim=(-0.09, 0.31), ylim=(-0.09, 0.31))
    a.set_aspect("equal")
    a.legend(loc="upper right", fontsize=8)
    a = ax[0, 1]
    panel(a, "b", "Observed COM: numerical controls only")
    for name, c, label in [
        ("OFF_DT1", GREY, "OFF, 2 ms"),
        ("WEAK_DT1", TEAL, "Weak, 2 ms"),
        ("STRONG_DT1", PURPLE, "Strong, valid 0.174 ms"),
    ]:
        t, x = trajectory(name)
        center = x.mean(1) * 1e6
        a.plot(center[:, 0], center[:, 1], c=c, lw=2, label=label)
        a.scatter(*center[-1, :2], c=c, s=30)
    a.set(xlabel="X (µm)", ylabel="Y (µm)")
    a.legend(fontsize=9)
    a.text(
        0.02,
        0.94,
        "No baseline or optimized XY policy was run",
        transform=a.transAxes,
        va="top",
        color=RED,
        fontsize=9,
    )
    a = ax[1, 0]
    panel(a, "c", "Actual distance to registered X target")
    for name, c, label in [
        ("OFF_DT1", GREY, "OFF"),
        ("WEAK_DT1", TEAL, "Weak"),
        ("STRONG_DT1", PURPLE, "Rejected strong audit"),
    ]:
        t, x = trajectory(name)
        err = np.linalg.norm(x.mean(1) - np.array(g["XY_targets_m"]["X"]), axis=1) * 1e6
        a.plot(t * 1e3, err, c=c, label=label)
    a.axhline(0.075, c=RED, ls="--", label="Target tolerance")
    a.set(xlabel="Actual physical time (ms)", ylabel="COM target error (µm)", ylim=(0.06, 0.21))
    a.legend(fontsize=9)
    a = ax[1, 1]
    a.axis("off")
    panel(a, "d", "Navigation and extension admission")
    a.text(
        0,
        0.96,
        "OFF numerical control: completed 2 ms\nBASELINE XY policy: NOT_RUN\nOPTIMIZED X / Y policies: NOT_RUN\n\nStrong B2 failed the unchanged overlap budget.\nNo independent X/Y navigation has been demonstrated.\n\n5 ms: NOT_RUN     10 ms: NOT_RUN     30 ms: NOT_RUN\nNo extrapolated trajectories are shown.\n\nNear-zero Z COM is specific to this symmetric state;\nit is not a theorem of global Z uncontrollability.",
        va="top",
        fontsize=11,
        linespacing=1.5,
    )
    save(
        fig, "Fig_F10O_3_XY_navigation", "F10O | Preregistered XY task stopped at the numerical gate", records
    )

    fig, ax = plt.subplots(2, 2, figsize=(12.8, 8.6))
    a = ax[0, 0]
    a.axis("off")
    panel(a, "a", "Admission levels remain separate")
    names = [
        "Dead-code cleanup",
        "B2 engineering bridge",
        "B2 mathematical residual",
        "Strong 2 ms overlap quality",
        "Strong 2 ms fixed-ID admission",
        "Weak 2 ms fixed-ID sensitivity",
        "New XY target navigation",
        "Experimental calibration",
        "AI training admission",
    ]
    states = [
        "PASS",
        "PASS",
        "PASS",
        "FAIL",
        "NOT ADMITTED",
        "PASS",
        "NOT_RUN",
        "NOT CALIBRATED",
        "training_ready=false",
    ]
    for j, (name, status) in enumerate(zip(names, states)):
        y = 0.94 - j * 0.102
        a.text(0, y, name, fontsize=10)
        a.text(0.69, y, status, fontsize=9.5, color=TEAL if status == "PASS" else RED)
    a = ax[0, 1]
    panel(a, "b", "Measured end-to-end step cost")
    tags = ["OFF", "WEAK", "STRONG"]
    vals = [runs["pilot_v1/" + tag + "_PILOT"]["mean_step_s"] * 1e3 for tag in tags]
    a.bar(tags, vals, color=[GREY, TEAL, PURPLE])
    a.set_ylabel("Wall time / step (ms)")
    for i, v in enumerate(vals):
        a.text(i, v + 0.12, f"{v:.2f}", ha="center")
    a.set_ylim(0, max(vals) * 1.3)
    a.text(
        0.03,
        0.95,
        "150 accepted steps per pilot\nGPU fields / walls + original CPU B2",
        va="top",
        transform=a.transAxes,
        fontsize=9,
    )
    a = ax[1, 0]
    panel(a, "c", "Actual physical time completed")
    labels = ["OFF", "Weak", "Strong dt=1 µs", "Strong dt=0.5 µs", "Strong fine cache"]
    keys = ["OFF_DT1", "WEAK_DT1", "STRONG_DT1", "STRONG_DT05", "STRONG_FINE"]
    vals = [runs["validation_v1/" + k]["completed_s"] * 1e3 for k in keys]
    a.barh(labels, vals, color=[GREY, TEAL, RED, RED, RED])
    a.set(xlabel="Accepted physical time (ms)", xlim=(0, 2.4))
    a.invert_yaxis()
    for j, v in enumerate(vals):
        a.text(v + 0.03, j, f"{v:.4g}", va="center", fontsize=9)
    a = ax[1, 1]
    a.axis("off")
    panel(a, "d", "Cost, stopping evidence and scientific scope")
    text = (
        f"Numerical processes: {g['numerical_process_total_s']:.2f} s\n"
        f"Integration and step evidence: {g['integration_s']:.2f} s\n"
        f"CPU B2 solver calls: {g['CPU_solver_sum_s']:.2f} s; resistance: {g['CPU_resistance_sum_s']:.2f} s\n"
        f"Peak process RSS: {g['peak_process_rss_MiB']:.1f} MiB\n\n"
        "12 runs: 9 complete short windows; 3 strong overlap rejections.\n"
        "Deleted one unused Warp helper (5 lines); history preserved.\n"
        "No new contact solver, stiffness change or relaxed threshold.\n\n"
        "Uniform-water acoustics, near-field and contact laws remain\n"
        "uncalibrated. Graph connectivity is not physical fusion.\n"
        "Times exclude deployment, imports, plotting and documentation.\n"
        "GPU kernel time was not profiled separately."
    )
    a.text(
        0,
        0.98,
        "\n".join(textwrap.fill(line, 63) if line else "" for line in text.splitlines()),
        va="top",
        fontsize=9.5,
        linespacing=1.25,
    )
    save(
        fig,
        "Fig_F10O_4_admission_cost",
        "F10O | Bounded computation and explicit admission decisions",
        records,
    )
    write(
        OUT / "FIGURE_PROVENANCE.json",
        dict(
            identity=identity(),
            source_script_sha256=sha(__file__),
            figures=records,
            input_sha256=SOURCES,
            scientific_execution_commits=[
                read(BASE / p / "CAMPAIGN.json")["identity"]["execution_commit"]
                for p in ["pilot_v1", "validation_v1"]
            ],
            no_synthetic_trajectories=True,
        ),
    )
    print("Saved four PNG/PDF figure pairs")


if __name__ == "__main__":
    main()
