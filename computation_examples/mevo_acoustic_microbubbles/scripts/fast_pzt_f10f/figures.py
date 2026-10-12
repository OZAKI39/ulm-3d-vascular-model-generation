"""Actual accepted coordinates and archived rejected diagnostics; no synthetic frames."""

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .common import BASE, ROOT, js, sha, write, identity
from .analyze import gz

COLORS = {"B0": "#34465e", "SS2": "#d08739", "SS4": "#7d5ba6", "C1": "#2b8a78"}
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9,
        "axes.labelcolor": "#20252c",
        "text.color": "#20252c",
        "axes.edgecolor": "#50555c",
        "xtick.color": "#30353c",
        "ytick.color": "#30353c",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    }
)


def title(ax, text):
    ax.set_title(text, loc="left", fontweight="bold", fontsize=10, pad=10)


def grid(ax):
    ax.grid(alpha=0.18)


def finish(fig, name, caption):
    fig.get_layout_engine().set(rect=(0, 0.045, 1, 0.955))
    fig.text(0.02, 0.008, caption, fontsize=8, color="#454a50")
    fig.savefig(OUT / (name + ".png"), dpi=220, bbox_inches="tight")
    fig.savefig(OUT / (name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def contact():
    fig, aa = plt.subplots(2, 3, figsize=(14, 8), layout="constrained")
    for ax, case, label in zip(
        aa[0, :2], ["O_1000_1102", "prepare_1104"], ["a  O / seed 1102", "b  Preparation / seed 1104"]
    ):
        for r in rows:
            if r["case"] != case:
                continue
            v = np.r_[r["final_predicted_nm"], r["final_rounds_nm"]]
            ax.plot(
                np.arange(len(v)),
                v,
                "o-",
                ms=3,
                color=COLORS[r["mode"]],
                label=f"{r['mode']}: rejected at {r['first_rejection_s'] * 1e3:g} ms",
            )
        ax.axhline(17.8, color="#aa3040", ls="--", lw=1, label="Residual budget: 17.8 nm")
        ax.set(xlabel="PBD round (0 = prediction)", ylabel="Maximum overlap (nm)")
        title(ax, label)
        ax.legend(fontsize=8)
        grid(ax)
    ax = aa[0, 2]
    title(ax, "c  First rejected substep")
    vals = [
        [
            next(r for r in rows if r["case"] == case and r["mode"] == m)["first_rejection_s"] * 1e3
            for m in ["B0", "SS2", "SS4"]
        ]
        for case in ["O_1000_1102", "prepare_1104"]
    ]
    for j, (val, label) in enumerate(zip(vals, ["O / 1102", "Preparation / 1104"])):
        for i, v in enumerate(val):
            ax.scatter(
                i + j * 0.15 - 0.075,
                v,
                s=40,
                color=COLORS[["B0", "SS2", "SS4"][i]],
                marker="o" if j == 0 else "s",
            )
            ax.annotate(
                f"{v:g}",
                (i + j * 0.15 - 0.075, v),
                xytext=(6, 1 if j == 0 else 8),
                textcoords="offset points",
                fontsize=8,
            )
    ax.set(
        xticks=range(3),
        xticklabels=["B0", "SS2", "SS4"],
        yscale="log",
        ylabel="Attempted end time (ms)",
        ylim=(0.001, 10),
    )
    ax.plot([], [], "o", color="#50555c", label="O / 1102")
    ax.plot([], [], "s", color="#50555c", label="Preparation / 1104")
    ax.legend()
    grid(ax)
    ax = aa[1, 0]
    title(ax, "d  Same B0 prediction: local Jacobi only")
    for v in con["same_proposal"]:
        ax.plot(
            range(13),
            np.r_[v["initial_prediction_nm"], v["residual_nm"]],
            "o-",
            ms=3,
            label="O / 1102" if v["case"].startswith("O") else "Preparation / 1104",
        )
    ax.axhline(17.8, color="#aa3040", ls="--")
    ax.set(xlabel="PBD round; no physical time advancement", ylabel="Maximum overlap (nm)")
    ax.legend()
    grid(ax)
    ax = aa[1, 1]
    title(ax, "e  Largest single-pass correction")
    labels = [r["mode"] + ("\nO" if r["case"].startswith("O") else "\nPrep.") for r in rows]
    ax.bar(range(6), [r["maximum_pass_nm"] for r in rows], color=[COLORS[r["mode"]] for r in rows])
    ax.axhline(178, color="#aa3040", ls="--", label="Pass budget: 178 nm")
    ax.set(xticks=range(6), xticklabels=labels, ylabel="Maximum over attempted history (nm)", ylim=(0, 205))
    ax.legend()
    grid(ax)
    ax = aa[1, 2]
    title(ax, "f  Actual accepted substep displacement")
    trace = gz(BASE / "short_v1/prepare_1104_SS4/CONTACT_TRACE.json.gz")["traces"]
    last = trace[-1]
    xx = np.array(last["proposal_positions_m"])
    d = np.linalg.norm(xx[:, None] - xx[None, :], axis=-1)
    d[np.diag_indices(32)] = np.inf
    pair = np.unravel_index(d.argmin(), d.shape)
    accepted = [t for t in trace if t["accepted"]][-5:]
    states = np.array([t["positions_after_step_m"] for t in accepted])
    t = np.array([v["attempted_end_s"] for v in accepted]) * 1e3
    for i in pair:
        for axis, ls in enumerate(["-", "--", ":"]):
            ax.plot(
                t,
                (states[:, i, axis] - states[0, i, axis]) * 1e9,
                ls,
                marker=".",
                ms=4,
                color=plt.cm.viridis(i / 31),
                label=f"ID {i}, {'XYZ'[axis]}",
            )
    ax.set(xlabel="Accepted physical time (ms)", ylabel="Change from first displayed state (nm)")
    ax.legend(ncol=2, fontsize=7)
    grid(ax)
    fig.suptitle("F10F | Smaller predictions did not pass the original residual budget", fontsize=14)
    finish(
        fig,
        "Fig_F10F_1_smallsteps_contact",
        "Panels a–c: independent histories and different rejection times. Panel d: identical saved B0 proposals. Panel f: accepted states only. Radius = 1.78 µm.",
    )


def spheres(ax, x, origin, limits):
    theta = np.linspace(0, 2 * np.pi, 33)
    phi = np.linspace(0, np.pi, 21)
    sx = 1.78 * np.outer(np.cos(theta), np.sin(phi))
    sy = 1.78 * np.outer(np.sin(theta), np.sin(phi))
    sz = 1.78 * np.outer(np.ones_like(theta), np.cos(phi))
    for i, point in enumerate((x - origin) * 1e6):
        ax.plot_surface(
            sx + point[0],
            sy + point[1],
            sz + point[2],
            color=plt.cm.viridis(i / 31),
            linewidth=0,
            antialiased=True,
            shade=True,
            rasterized=True,
        )
    for fun, lim in zip([ax.set_xlim, ax.set_ylim, ax.set_zlim], limits):
        fun(lim)
    ax.set_box_aspect([v[1] - v[0] for v in limits])
    ax.view_init(elev=22, azim=-58)
    ax.set(xlabel="X (µm)", ylabel="Y (µm)", zlabel="Z (µm)")
    ax.tick_params(labelsize=7, pad=0)
    for axis in [ax.xaxis, ax.yaxis, ax.zaxis]:
        axis.pane.set_facecolor((1, 1, 1, 0))
        axis._axinfo["grid"]["color"] = (0.4, 0.45, 0.5, 0.2)


def morphology():
    fig = plt.figure(figsize=(14, 12), layout="constrained")
    fig.get_layout_engine().set(hspace=0.13)
    gs = fig.add_gridspec(3, 3, height_ratios=[1.35, 1, 1])
    arrays = []
    for m in ["B0", "SS2", "SS4"]:
        with h5py.File(BASE / "short_v1" / ("prepare_1104_" + m) / "trajectory.h5") as h:
            k = np.argmin(np.abs(h["time_s"][:] - 0.0017))
            assert abs(h["time_s"][k] - 0.0017) < 1e-14
            arrays.append(h["positions_m"][k, 0])
            origin = h["positions_m"][0, 0].mean(0)
    cloud = (np.concatenate(arrays) - origin) * 1e6
    lo = cloud.min(0) - 2.2
    hi = cloud.max(0) + 2.2
    width = max(hi - lo)
    center = (lo + hi) / 2
    limits = list(zip(center - width / 2, center + width / 2))
    for j, (x, m) in enumerate(zip(arrays, ["B0", "SS2", "SS4"])):
        ax = fig.add_subplot(gs[0, j], projection="3d")
        spheres(ax, x, origin, limits)
        title(ax, f"{'abc'[j]}  {m} | accepted 1.700 ms")
    axes = [fig.add_subplot(gs[i, j]) for i in [1, 2] for j in range(3)]
    for mi, label in [(0, "SS2 vs B0"), (1, "SS4 vs B0")]:
        r = morph["comparisons"][mi]
        s = r["series"]
        time = np.array(s["time_s"]) * 1e3
        c = COLORS["SS2" if mi == 0 else "SS4"]
        for ax, key, scale in [
            (axes[0], "fixed_ID_max_m", 1e9),
            (axes[1], "Kabsch_RMSD_m", 1e9),
            (axes[2], "centroid_difference_m", 1e9),
        ]:
            ax.plot(time, np.array(s[key]) * scale, color=c, label=label)
        axes[3].plot(
            time, (np.array(s["Rg_candidate_m"]) / s["Rg_reference_m"] - 1) * 100, color=c, label=label
        )
        axes[4].plot(
            time, (np.array(s["aspect_candidate"]) / s["aspect_reference"] - 1) * 100, color=c, label=label
        )
    for ax, ttl, ylabel in zip(
        axes[:5],
        [
            "d  Fixed-ID position sensitivity",
            "e  Kabsch-aligned shape sensitivity",
            "f  Centroid sensitivity",
            "g  Radius of gyration",
            "h  Aspect ratio",
        ],
        [
            "Maximum position difference (nm)",
            "RMSD (nm)",
            "3D centroid difference (nm)",
            "Relative difference (%)",
            "Relative difference (%)",
        ],
    ):
        title(ax, ttl)
        ax.set(xlabel="Common accepted physical time (ms)", ylabel=ylabel)
        grid(ax)
    axes[0].axhline(178, ls="--", color="#aa3040", label="178 nm gate")
    axes[1].axhline(89, ls="--", color="#aa3040", label="89 nm gate")
    axes[0].legend(fontsize=8)
    axes[1].legend(fontsize=8)
    ax = axes[5]
    title(ax, "i  Edge events and evidence limits")
    counts = [
        morph["comparisons"][0]["values"]["event_count"],
        morph["comparisons"][1]["values"]["event_count"],
    ]
    ax.bar(
        [0, 1, 3, 4],
        np.array(counts).ravel(),
        color=[COLORS["B0"], COLORS["SS2"], COLORS["B0"], COLORS["SS4"]],
    )
    ax.set(
        xticks=[0, 1, 3, 4],
        xticklabels=["B0\n2.10 ms", "SS2\n2.10 ms", "B0\n1.70 ms", "SS4\n1.70 ms"],
        ylabel="Accepted pair-edge changes",
        ylim=(0, 160),
    )
    ax.text(
        0.04,
        0.96,
        "O / 1102: no accepted SS step\nPersistent events: mostly right-censored\nMORPHOLOGY_NOT_CONVERGED",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
    )
    fig.suptitle("F10F | Fixed-ID geometry remains numerically sensitive", fontsize=14)
    finish(
        fig,
        "Fig_F10F_2_morphology_convergence",
        "Preparation / seed 1104. True radius = 1.78 µm; identical origin, 3D limits and camera. Fixed ID colors: viridis, 0 → 31. No rejected state is drawn as accepted.",
    )


def cost():
    fig, aa = plt.subplots(2, 3, figsize=(14, 8), layout="constrained")
    ax = aa[0, 0]
    title(ax, "a  Diagnostic integration cost")
    for offset, case in [(-0.18, "O_1000_1102"), (0.18, "prepare_1104")]:
        vals = [
            next(r for r in rows if r["case"] == case and r["mode"] == m)["integration_ms_per_10us"]
            for m in ["B0", "SS2", "SS4"]
        ]
        ax.bar(
            np.arange(3) + offset,
            vals,
            0.34,
            label="O / 1102" if case.startswith("O") else "Preparation / 1104",
        )
    ax.set(xticks=range(3), xticklabels=["B0", "SS2", "SS4"], ylabel="Wall ms / attempted 10 µs interval")
    ax.legend()
    grid(ax)
    ax = aa[0, 1]
    title(ax, "b  Fresh work per 10 µs interval")
    for j, (name, val) in enumerate(
        [
            ("Acoustic / pair", [1, 2, 4]),
            ("CFD", [2, 4, 8]),
            ("Lubrication", [2, 4, 8]),
            ("PBD", [12, 12, 12]),
        ]
    ):
        ax.bar(np.arange(3) + (j - 1.5) * 0.18, val, 0.18, label=name)
    ax.set(xticks=range(3), xticklabels=["B0", "SS2", "SS4"], ylabel="Calls (original full-step schedule)")
    ax.legend(fontsize=8)
    grid(ax)
    ax = aa[0, 2]
    title(ax, "c  Actual attempted physics work")
    ax.bar(range(6), [r["attempted_steps"] for r in rows], color=[COLORS[r["mode"]] for r in rows])
    ax.set(
        xticks=range(6),
        xticklabels=[r["mode"] + ("\nO" if r["case"].startswith("O") else "\nPrep.") for r in rows],
        ylabel="Physical substeps including rejection",
    )
    for i, r in enumerate(rows):
        ax.text(i, r["attempted_steps"] + 10, str(r["attempted_steps"]), ha="center", fontsize=8)
    grid(ax)
    ax = aa[1, 0]
    title(ax, "d  All scientific replay costs")
    bottom = np.zeros(6)
    for key, label in [
        ("model_prepare_s", "Model preparation"),
        ("graph_capture_s", "Graph build"),
        ("integration_quality_wall_s", "Integration + diagnostics"),
        ("sync_readback_wall_s", "Saved-state readback"),
        ("morphology_quality_output_prepare_s", "Morphology + quality"),
        ("hdf5_output_s", "HDF5"),
        ("observer_output_s", "Trace export"),
    ]:
        val = np.array([r["timing"][key] for r in rows])
        ax.bar(range(6), val, bottom=bottom, label=label)
        bottom += val
    ax.set(
        xticks=range(6),
        xticklabels=[r["mode"] + ("\nO" if r["case"].startswith("O") else "\nPrep.") for r in rows],
        ylabel="Measured wall time (s)",
    )
    ax.legend(fontsize=6.5, ncol=2)
    grid(ax)
    ax = aa[1, 1]
    title(ax, "e  Effective full-duration yield")
    ax.axis("off")
    ax.text(
        0.02,
        0.92,
        f"SS attempts passing original failure windows: 0 / 4\n\nNew complete 10 ms control trajectories: 0\nNew OFF / S / O triplets: 0\nReused independent seed families: 2\n\nAll cold processes: {costs['cold_all_process_s']:.3f} s\nIncludes engineering checks and failed attempts\n\nNew complete trajectories / reserved GPU hour: 0\nCost per new qualified sample: undefined (zero yield)",
        va="top",
        fontsize=9,
    )
    ax = aa[1, 2]
    title(ax, "f  Performance admission")
    ax.axis("off")
    ax.text(
        0.02,
        0.92,
        "Formal production benchmark: NOT ADMITTED\n\nAll Small Steps candidates failed contact gate.\nDiagnostic CPU readbacks are included in a,d.\nCUDA span includes host gaps; not pure kernel time.\n\nHistorical C1 (5 µs / 12 rounds):\n1 O / 1102 trajectory completed 10 ms;\nno matched SS full-duration timing is available.\n\nNo SS12, full-duration SS or batch expansion.\nNo demonstrated production speedup.",
        va="top",
        fontsize=9,
    )
    fig.suptitle("F10F | Equal PBD rounds do not imply equal total cost", fontsize=14)
    finish(
        fig,
        "Fig_F10F_3_cost_and_yield",
        "All timing bars are diagnostic runs with different stopping times; they are not a matched production speed comparison. GPU-hour denotes exclusive card reservation.",
    )


def admission():
    fig, aa = plt.subplots(2, 2, figsize=(13, 8), layout="constrained")
    ax = aa[0, 0]
    title(ax, "a  Contact outcome by original case")
    ax.axis("off")
    cell = []
    for r in rows:
        cell.append(
            [
                ("O / 1102" if r["case"].startswith("O") else "Prep. / 1104"),
                r["mode"],
                f"{r['first_rejection_s'] * 1e3:g}",
                f"{r['final_residual_nm']:.3f}",
                "FAIL",
            ]
        )
    tb = ax.table(
        cellText=cell,
        colLabels=["Case", "Mode", "Rejected (ms)", "Residual (nm)", "17.8 nm gate"],
        loc="center",
        cellLoc="center",
    )
    tb.auto_set_font_size(False)
    tb.set_fontsize(8)
    tb.scale(1, 1.8)
    ax = aa[0, 1]
    title(ax, "b  Strict common-prefix morphology gates")
    keys = [
        "centroid",
        "Rg",
        "fixed_ID",
        "Kabsch",
        "aspect",
        "tensor",
        "contact_path",
        "component_membership",
        "edge_events",
        "persistent_events",
    ]
    mm = morph["comparisons"][:3]
    arr = np.array(
        [[1 if r["gates"][k] is True else 0 if r["gates"][k] is False else 0.5 for k in keys] for r in mm]
    )
    from matplotlib.colors import ListedColormap

    ax.imshow(arr, vmin=0, vmax=1, cmap=ListedColormap(["#c37479", "#d6d1bd", "#75a99d"]), aspect="auto")
    ax.set(
        yticks=range(3),
        yticklabels=["B0–SS2\n0–2.10 ms", "B0–SS4\n0–1.70 ms", "SS2–SS4\n0–1.70 ms"],
        xticks=range(len(keys)),
        xticklabels=[
            "Centroid",
            "Rg",
            "Fixed ID",
            "Kabsch",
            "Aspect",
            "Tensor",
            "Contact path",
            "Members",
            "Edges",
            "Persistence",
        ],
    )
    plt.setp(ax.get_xticklabels(), rotation=50, ha="right", fontsize=8)
    for i in range(3):
        for j in range(len(keys)):
            ax.text(
                j,
                i,
                "P" if arr[i, j] == 1 else "F" if arr[i, j] == 0 else "N/E",
                ha="center",
                va="center",
                fontsize=8,
            )
    ax = aa[1, 0]
    title(ax, "c  Observed failure mechanism")
    ax.axis("off")
    ax.text(
        0.02,
        0.95,
        "Same O / 1102 starting coordinates:\nPrediction overlap: B0 85.544 nm → SS2 50.535 nm\n                                        → SS4 32.978 nm\n\nSmaller step reduces the predicted penetration.\nFewer Jacobi passes leave 20.920 / 20.361 nm residual.\nThe original 17.8 nm limit rejects both first substeps.\n\nOn the same B0 rejected proposal, one round leaves\n66.191 nm (O) / 55.276 nm (preparation);\n12 rounds still leave 17.831 nm.\nThis local test is not a new physical trajectory.",
        va="top",
        fontsize=9,
    )
    ax = aa[1, 1]
    title(ax, "d  Method admission and scope")
    ax.axis("off")
    ax.text(
        0.02,
        0.95,
        "Engineering checks: PASS\nOriginal historical defaults and results: preserved\nSmall Steps contact windows: FAIL (4 / 4)\nGroup centroid agreement: small on accepted prep prefixes\nFull morphology admission: NOT PASSED\nProduction speed / small batch: NOT ADMITTED\n\nSupported evidence: two fixed original cases, 32 bubbles.\nNot supported: dense long-time production, independent\ngeneralization, morphology training truth or calibrated shells.\n\nNext: local lubrication convergence and near-field audit;\nassess contact prediction only with unchanged acceptance.\nEXPLORATORY / PHYSICS_UNCALIBRATED / training_ready=false",
        va="top",
        fontsize=9,
    )
    fig.suptitle("F10F | Small Steps is not admitted as a production numerical mode", fontsize=14)
    finish(
        fig,
        "Fig_F10F_4_method_admission",
        "P = passed only the displayed common prefix; F = failed; N/E = right-censored. O / 1102 SS morphology is not evaluated: no accepted nonzero substep.",
    )


def main():
    global OUT, con, rows, morph, costs
    OUT = BASE / "figures"
    OUT.mkdir(exist_ok=True)
    con = js(BASE / "analysis/CONTACT.json")
    rows = con["rows"]
    morph = js(BASE / "analysis/MORPHOLOGY.json")
    costs = js(BASE / "analysis/COST_AND_YIELD.json")
    for fn in [contact, morphology, cost, admission]:
        fn()
    write(
        OUT / "SOURCE_MANIFEST.json",
        dict(
            render_source=identity(),
            source_data={
                str(p.relative_to(ROOT)): sha(p)
                for folder in [BASE / "analysis", BASE / "short_v1"]
                for p in folder.rglob("*")
                if p.is_file()
            },
            outputs={p.name: sha(p) for p in OUT.iterdir() if p.suffix in [".png", ".pdf"]},
            radius_m=1.78e-6,
            synthetic_frames=False,
        ),
    )
    print("4 PNG + 4 PDF", flush=True)


if __name__ == "__main__":
    main()
