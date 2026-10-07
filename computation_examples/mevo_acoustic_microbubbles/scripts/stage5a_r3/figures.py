"""G1–G8 from actual states and retained diagnostics; PNG/PDF plus source arrays."""

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scripts.stage5a_r3.common import BASE, R2, B2_RUNS, read, write, sha
from scripts.stage5a_r3.paths import local
from scripts.stage5a_r2.motion_visuals import wall_path, glyph, decorate

BLUE = "#356fa7"
RED = "#b85c5c"
PURPLE = "#80649c"
GRAY = "#686d76"
GOLD = "#ba8c47"
plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelcolor": "#25282d",
        "text.color": "#25282d",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
        "axes.grid": False,
    }
)
OUT = BASE / "figures"
DATA = OUT / "data"
SOURCES = {}


def source(path, lines=False):
    path = local(path)
    SOURCES[str(path)] = sha(path)
    return read(path)


def history(path):
    path = local(path)
    SOURCES[str(path)] = sha(path)
    with h5py.File(path) as h:
        return {k: h[k][:] for k in ["time_s", "positions_m", "particle_id", "radii_m", "active"]}


def save(fig, name, caption, arrays=None):
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    if arrays:
        np.savez_compressed(DATA / (name + ".npz"), **arrays)
    for suffix in ["png", "pdf"]:
        fig.savefig(OUT / (name + "." + suffix), dpi=260, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    write(
        OUT / (name + "_SOURCES.json"),
        dict(
            caption=caption,
            sources=SOURCES.copy(),
            script_sha256=sha(__file__),
            array_file=str(DATA / (name + ".npz")) if arrays else None,
            units="SI in stored trajectories; axes state display unit conversions",
        ),
    )


def event_data():
    p = next((BASE / "runs").glob("*local_v1_B2_event_L0/trajectories.h5"))
    new = history(p)
    old = history(R2 / "runs" / B2_RUNS["nominal"] / "trajectories.h5")
    begin = new["time_s"][0]
    take = old["time_s"] >= begin - 1e-15
    old["positions_m"], old["time_s"], old["active"] = (
        old["positions_m"][take],
        old["time_s"][take],
        old["active"][take],
    )
    return old, new, p


def g1():
    import pyvista as pv

    old, new, path = event_data()
    j = int(np.flatnonzero(new["particle_id"] == 152)[0])
    origin = new["positions_m"][0, j]
    radius = float(new["radii_m"][j])
    plot = pv.Plotter(off_screen=True, window_size=(1350, 1050))
    wall = pv.read(wall_path())
    bounds = [-5.0, 5.0, -5.0, 5.0, -5.0, 5.0]
    decorate(plot, wall, bounds, origin, False, "Actual finite WALL and bubble 152")
    positions = new["positions_m"][0][new["active"][0]]
    glyph(plot, (positions - origin) * 1e6, radius * 1e6, BLUE)
    # Highlight the two actual source faces; no reconstructed smooth wall.
    tri = wall.triangulate()
    for cell, color in [(3178, "#bb7372"), (13018, "#ceae9c")]:
        mesh = tri.extract_cells([cell]).extract_surface()
        mesh.points = (mesh.points - origin) * 1e6
        for axis in range(3):
            n = np.eye(3)[axis]
            lo = np.zeros(3)
            hi = np.zeros(3)
            lo[axis] = bounds[2 * axis]
            hi[axis] = bounds[2 * axis + 1]
            mesh = mesh.clip(normal=n, origin=lo, invert=False).clip(normal=n, origin=hi, invert=True)
        plot.add_mesh(mesh, color=color, opacity=0.25, show_edges=True, edge_color=GRAY, line_width=1)
    image = plot.screenshot(return_img=True)
    plot.close()
    fig = plt.figure(figsize=(12, 5.2), layout="constrained")
    left = fig.add_subplot(1, 2, 1)
    left.imshow(image)
    left.axis("off")
    left.set_title("a  Actual radius: 1.78 µm", loc="left", fontweight="bold")
    ax = fig.add_subplot(1, 2, 2, projection="3d")
    for data, color, label in [(old, RED, "R2 accepted path (stalled)"), (new, BLUE, "R3 accepted path")]:
        index = int(np.flatnonzero(data["particle_id"] == 152)[0])
        x = (data["positions_m"][:, index] - origin) * 1e9
        ax.plot(*x.T, color=color, lw=1.8, label=label)
        ax.scatter(*x[-1], color=color, s=24)
    ax.set(xlabel="Δx (nm)", ylabel="Δy (nm)", zlabel="Δz (nm)", title="b  Bubble-center displacement")
    ax.view_init(elev=24, azim=-55)
    ax.legend(loc="upper left", fontsize=8)
    ax.set_box_aspect((1, 1, 1))
    save(
        fig,
        "G1_B2_wall_event",
        "Same full 186-particle system. Left: actual sphere radius and real finite WALL faces; right: center displacement in nanometres. Curves join only recorded accepted endpoints; no interpolated states.",
        dict(
            old_time_s=old["time_s"],
            old_positions_m=old["positions_m"][:, j],
            new_time_s=new["time_s"],
            new_positions_m=new["positions_m"][:, j],
            origin_m=origin,
            radius_m=radius,
        ),
    )


def g2():
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), layout="constrained")
    arrays = {}
    for ax, case in zip(axes, ["nominal", "refined"]):
        report = source(BASE / "diagnosis/reproduction_v1" / case / "RESULT.json")
        rows = report["trials"]
        dt = np.array([r["dt_s"] for r in rows])
        old = np.array([r["old_safety"]["minimum_wall_swept_gap_bound_m"] for r in rows])
        exact = np.array(
            [
                next(z for z in r["exact_geometry"] if z["particle_id"] == 152)["minimum"]["gap_m"]
                for r in rows
            ]
        )
        ax.plot(dt * 1e9, (old + 1e-11) * 1e12, "x--", color=RED, label="R2 lower bound")
        ax.plot(dt * 1e9, (exact + 1e-11) * 1e12, "o-", color=BLUE, label="Exact continuous minimum")
        ax.axhline(0, color="black", lw=0.8)
        ax.set(
            xscale="log",
            yscale="symlog",
            xlabel="Requested step (ns)",
            ylabel="Safety margin: min gap + tolerance (pm)",
            title=f"{case.capitalize()} terminal state",
        )
        ax.set_yscale("symlog", linthresh=1e-4)
        ax.annotate(
            "Safe path rejected by R2",
            (dt[0] * 1e9, (exact[0] + 1e-11) * 1e12),
            xytext=(0.04, 0.86),
            textcoords="axes fraction",
            fontsize=8,
            arrowprops={"arrowstyle": "->", "color": GRAY},
        )
        ax.legend(fontsize=8, loc="lower left")
        arrays[case + "_dt_s"] = dt
        arrays[case + "_old_bound_m"] = old
        arrays[case + "_exact_gap_m"] = exact
    save(
        fig,
        "G2_continuous_wall_certificate",
        "Exact rational finite-triangle minimum on each actual full-system candidate segment. Positive margin is safe; negative margin exceeds the original 1e-11 m tolerance. Rejected candidates never become accepted trajectories.",
        arrays,
    )


def g3():
    old, new, path = event_data()
    terminal = float(old["time_s"][-1])
    r2 = source(R2 / "runs" / B2_RUNS["nominal"] / "ACCEPTED.jsonl")
    r3 = source(path.parent / "ACCEPTED.jsonl")
    reject = source(R2 / "runs" / B2_RUNS["nominal"] / "REJECTED.jsonl")
    fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True, layout="constrained")
    arrays = {}
    for name, rows, color in [("R2", r2, RED), ("R3", r3, BLUE)]:
        rows = [r for r in rows if r["time_s"] >= new["time_s"][0] - 1e-15]
        t = np.array([r["time_s"] - terminal for r in rows]) * 1e6
        dt = np.array([r["dt_s"] for r in rows])
        s0 = np.array([r["substep_method"] == "S0" for r in rows])
        ax[0].plot(t, dt * 1e9, ".-", color=color, label=name + " accepted")
        ax[0].scatter(
            t[s0], dt[s0] * 1e9, s=27, marker="s", color=GOLD, label=name + " S0" if s0.any() else None
        )
        ax[1].plot(t, [r["kkt"] for r in rows], color=color, label=name + " QP KKT")
        ax[2].plot(
            t, [r["safety"]["minimum_wall_swept_gap_bound_m"] * 1e12 for r in rows], color=color, label=name
        )
        arrays[name + "_time_relative_us"] = t
        arrays[name + "_dt_s"] = dt
        arrays[name + "_S0"] = s0
    rejects = [r for r in reject if r["time_s"] >= new["time_s"][0] - 1e-15]
    ax[0].scatter(
        [(r["time_s"] - terminal) * 1e6 for r in rejects],
        [r["dt_s"] * 1e9 for r in rejects],
        color=GRAY,
        marker="x",
        s=22,
        label="R2 rejected trial",
    )
    ax[0].set(yscale="log", ylabel="Physical step (ns)")
    ax[0].legend(ncol=3, fontsize=8)
    ax[1].set(yscale="log", ylabel="Original-QP KKT residual")
    ax[1].axhline(1e-9, color="black", ls="--", lw=0.8, label="Acceptance limit")
    ax[1].legend(fontsize=8)
    ax[2].axhline(-10, color="black", ls="--", lw=0.8)
    ax[2].set(ylabel="Continuous WALL lower bound (pm)", xlabel="Time relative to R2 terminal state (µs)")
    for a in ax:
        a.axvline(0, color=GRAY, ls=":", lw=0.8)
    fig.suptitle("Local continuation through the original WALL feature switch", fontweight="bold")
    save(
        fig,
        "G3_S0_S1_event_diagnostics",
        "Accepted half-step sizes, bounded S0 usage, rejected trials and independent QP residuals. The feature switch occurs at the original R2 terminal time (vertical dotted line); rejection points are diagnostic candidates only.",
        arrays,
    )


def g4():
    report = source(BASE / "final/ALL_ACTIVE_TRAJECTORY_AUDIT.json")
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    arrays = {}
    for col, case in enumerate(["B2", "B3"]):
        row = next(r for r in report["rows"] if r["case"] == case and r["name"] == "r2_refinement")
        data = dict(np.load(local(row["samples"])))
        err = data["all_active_raw_errors_m"] * 1e9
        axes[0, col].hist(err[-1], bins=24, color=BLUE, alpha=0.8)
        axes[0, col].axvline(89, color=RED, ls="--", label="89 nm resolution target")
        axes[0, col].set(
            title=f"{case}: R2 nominal / refined",
            xlabel="Endpoint same-ID displacement (nm)",
            ylabel="Active bubbles",
        )
        axes[0, col].legend(fontsize=8)
        image = axes[1, col].imshow(
            err.T,
            origin="lower",
            aspect="auto",
            extent=[0, (data["clock_s"][-1] - 3) * 1e3, 0, len(data["particle_id"])],
            cmap="magma",
            vmin=0,
            vmax=max(89, np.quantile(err, 0.995)),
        )
        axes[1, col].set(xlabel="Time since 3 s (ms)", ylabel="Original active-ID order")
        fig.colorbar(image, ax=axes[1, col], label="Same-ID error (nm)")
        arrays[case + "_errors_m"] = err * 1e-9
        arrays[case + "_particle_id"] = data["particle_id"]
        arrays[case + "_clock_s"] = data["clock_s"]
    fig.suptitle("All-active lineage audit: original IDs retained", fontweight="bold")
    save(
        fig,
        "G4_all_active_error_distribution",
        "Pre-R3 baseline audit: all 186/189 active particles. Error distribution and heat map use last actual stored states, with sampling uncertainty retained in linked audits. These plots are descriptive, not an independent accuracy PASS.",
        arrays,
    )


def g5():
    report = source(BASE / "final/REFERENCE_LINEAGE_AUDIT.json")
    rows = report.get("controlled_attribution", [])
    if not rows:
        raise RuntimeError("Controlled attribution results not synchronized")
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout="constrained")
    names = [
        "attribution_source_only",
        "attribution_clock_only",
        "attribution_device_setting_only",
        "legacy_fine_vs_R3_fine",
    ]
    labels = [
        "Contact-solver\nsource only",
        "500 / 10 µs\ncutoff clock",
        "CPU / GPU\nsetting",
        "Refined legacy / R3\n(first 10 µs)",
    ]
    arrays = {}
    for ax, case in zip(axes, ["B2", "B3"]):
        chosen = [next(r for r in rows if r["case"] == case and r["name"] == name) for name in names]
        rms = (
            np.array([r["groups"]["all_active"]["errors_m"]["raw"]["endpoint"]["rms"] for r in chosen]) * 1e9
        )
        maximum = (
            np.array([r["groups"]["all_active"]["errors_m"]["raw"]["endpoint"]["max"] for r in chosen]) * 1e9
        )
        xx = np.arange(4)
        ax.bar(xx - 0.16, rms, 0.32, label="RMS", color=BLUE)
        ax.bar(xx + 0.16, maximum, 0.32, label="Maximum", color=PURPLE)
        ax.set_yscale("symlog", linthresh=0.01)
        ax.set_xticks(xx, labels, fontsize=8)
        ax.set(ylabel="Endpoint same-ID difference (nm)", title=case)
        ax.legend(fontsize=8)
        for j, v in enumerate(maximum):
            if v == 0:
                ax.text(j, 0, "0", ha="center", va="bottom", fontsize=8)
        arrays[case + "_rms_m"] = rms * 1e-9
        arrays[case + "_max_m"] = maximum * 1e-9
    save(
        fig,
        "G5_reference_lineage_attribution",
        "Controlled short replays isolate source, cutoff-clock and device factors (40 µs). The separate refined cross-method comparison covers 10 µs. Identical QP common-state velocities agree at about 1e-12 m/s. No legacy source is promoted to truth.",
        arrays,
    )


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    for function in [g1, g2, g3, g4, g5]:
        function()


if __name__ == "__main__":
    main()
