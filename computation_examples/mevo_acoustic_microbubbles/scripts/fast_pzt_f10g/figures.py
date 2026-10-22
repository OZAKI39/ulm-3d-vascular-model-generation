"""English publication figures from saved F10G numerical evidence only."""

import numpy as np
import h5py
from matplotlib.colors import ListedColormap
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from .common import BASE, ROOT, CASES, state, js, sha, write, identity
from .analyze import gz
from mevo_mb.fast_sim.f10g_contact import shape, overlaps

OUT = BASE / "figures"
OUT.mkdir(exist_ok=True)
C = ["#405a77", "#91629c", "#2d8b83", "#b47932"]
LABEL = ["O / 1102", "Preparation / 1104"]
plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelcolor": "#222733",
        "text.color": "#222733",
        "axes.edgecolor": "#434953",
        "xtick.color": "#30343b",
        "ytick.color": "#30343b",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    }
)
L = js(BASE / "lubrication_v1/RESULTS.json")["cases"]
P = js(BASE / "prototype_v1/DENSE_STATIC.json")["rows"]
S = js(BASE / "prototype_v1/SMALL_SYSTEMS.json")
M = js(BASE / "analysis/MORPHOLOGY.json")["comparisons"]
T = js(BASE / "analysis/TRAJECTORY_AUDIT.json")["rows"]
K = js(BASE / "performance_v2/RESULTS.json")["rows"]
COST = js(BASE / "analysis/COST_AND_YIELD.json")


def title(ax, s):
    ax.set_title(s, loc="left", fontweight="bold", fontsize=11, pad=10)


def grid(ax):
    ax.grid(alpha=0.18)


def finish(fig, name, caption):
    fig.get_layout_engine().set(rect=(0, 0.05, 1, 0.955))
    fig.text(0.025, 0.008, caption, fontsize=9, color="#454953")
    fig.savefig(OUT / (name + ".png"), dpi=300, bbox_inches="tight")
    fig.savefig(OUT / (name + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def lubrication():
    fig, aa = plt.subplots(2, 3, figsize=(15, 9), layout="constrained")
    for j, c in enumerate(L):
        rr = c["rows"]
        xx = np.arange(4)
        aa[0, 0].semilogy(
            xx, [r["residuals"]["balanced"]["max_residual_m_s"] for r in rr], "o-", color=C[j], label=LABEL[j]
        )
        aa[0, 1].plot(
            xx,
            [100 * r["residuals"]["balanced"]["normalized_preconditioned"] for r in rr],
            "o-",
            color=C[j],
            label=LABEL[j],
        )
        aa[0, 2].plot(
            xx, [-min(r["relative_normal_velocity_m_s"]) * 1e3 for r in rr], "o-", color=C[j], label=LABEL[j]
        )
        aa[1, 0].plot(
            xx, [r["predicted_overlap_m"] * 1e9 for r in rr], "o-", color=C[j], label=LABEL[j] + " prediction"
        )
        aa[1, 0].plot(
            xx, [r["final_overlap_m"] * 1e9 for r in rr], "s--", color=C[j], label=LABEL[j] + " after PBD"
        )
        aa[1, 2].plot(
            xx, [r["timing"]["median_s"] * 1e6 for r in rr], "o-", color=C[j], label=LABEL[j] + " CPU solve"
        )
    for ax in [*aa[0], aa[1, 0], aa[1, 2]]:
        ax.set(
            xticks=range(4),
            xticklabels=["2", "4", "8", "Reference"],
            xlabel="Jacobi iterations / direct reference",
        )
        grid(ax)
    title(aa[0, 0], "a  Residual of the same resistance system")
    aa[0, 0].set_ylabel("Max. velocity-unit residual (m/s)")
    aa[0, 0].legend()
    title(aa[0, 1], "b  Block-preconditioned residual")
    aa[0, 1].set_ylabel("Max. residual / max. free speed (%)")
    aa[0, 1].legend()
    title(aa[0, 2], "c  Fastest normal approach")
    aa[0, 2].set_ylabel("Approach speed over neighbor pairs (mm/s)")
    aa[0, 2].legend()
    title(aa[1, 0], "d  Same state, same force, same 10 µs")
    aa[1, 0].axhline(17.8, color="#a84246", ls=":", label="Final-overlap budget")
    aa[1, 0].set_ylabel("Maximum overlap (nm)")
    aa[1, 0].legend(fontsize=8)
    title(aa[1, 1], "e  Original 12 Jacobi PBD passes")
    for i, r in enumerate(L[0]["rows"]):
        aa[1, 1].plot(
            range(13),
            np.r_[r["predicted_overlap_m"], [v["residual_m"] for v in r["pbd_rounds"]]] * 1e9,
            "o-",
            ms=3,
            color=C[i],
            label=r["mode"],
        )
    aa[1, 1].axhline(17.8, color="#a84246", ls=":")
    aa[1, 1].set(xlabel="PBD round (0 = prediction), O / 1102", ylabel="Maximum overlap (nm)")
    aa[1, 1].legend(title="Lubrication")
    grid(aa[1, 1])
    title(aa[1, 2], "f  CPU solve cost, preassembled operator")
    aa[1, 2].set_ylabel("Median wall time (µs), 25 calls")
    aa[1, 2].legend(fontsize=8)
    finish(
        fig,
        "Fig_F10G_1_lubrication_convergence",
        "Frozen F10E first-failure inputs; original gap floor, forces, CFD and 10 µs step. Mean balance is included in panels a–d.\nAlgebraic convergence is not experimental hydrodynamic accuracy. CPU matrix assembly excluded in panel f; no new trajectory in this figure.",
    )


def implicit():
    fig, aa = plt.subplots(2, 3, figsize=(15, 9), layout="constrained")
    for j, mult in enumerate([0.5, 1, 2]):
        rows = [
            r
            for r in S["two"]
            if r["name"] == "overlapping" and r["mode"] == "B1" and r["multiplier"] == mult
        ]
        dt = [r["dt_s"] * 1e6 for r in rows]
        aa[0, 0].plot(
            [r["analytic_overlap_m"] * 1e9 for r in rows],
            [r["result"]["pair_force_N"][0] * 1e12 for r in rows],
            "o-",
            color=C[j],
            label=f"{mult:g} × kc",
        )
        aa[0, 1].plot(
            dt,
            [r["analytic_overlap_m"] * 1e9 for r in rows],
            "o-",
            color=C[j],
            label=f"Implicit, {mult:g} × kc",
        )
        aa[0, 1].plot(
            dt,
            [r["explicit_signed_overlap_m"] * 1e9 for r in rows],
            "x--",
            color=C[j],
            label=f"Explicit, {mult:g} × kc",
        )
    for j, mode in enumerate(["B1", "B2"]):
        rows = [
            r for r in S["two"] if r["name"] == "approaching" and r["mode"] == mode and r["multiplier"] == 1
        ]
        aa[0, 2].plot(
            [r["dt_s"] * 1e6 for r in rows],
            [r["result"]["relative_normal_velocity_m_s"][0] * 1e3 for r in rows],
            "o-",
            color=C[j],
            label=mode,
        )
    title(aa[0, 0], "a  Candidate spring: force and overlap")
    aa[0, 0].set(xlabel="Implicit overlap (nm)", ylabel="Repulsive force (pN)")
    aa[0, 0].legend()
    title(aa[0, 1], "b  Explicit overshoot versus implicit decay")
    aa[0, 1].set(xlabel="Time step (µs)", ylabel="Signed linearized overlap (nm)")
    aa[0, 1].axhline(0, color="0.5", lw=0.8)
    aa[0, 1].legend(fontsize=8, ncol=2)
    title(aa[0, 2], "c  Normal velocity under compression")
    aa[0, 2].set(xlabel="Time step (µs)", ylabel="Separation velocity (mm/s)")
    aa[0, 2].legend()
    for j, mode in enumerate(["B1", "B2"]):
        for index in range(3):
            rows = [r for r in S["four"] if r["mode"] == mode and r["configuration"] == index]
            aa[1, 0].plot(
                [r["multiplier"] for r in rows],
                [max(r["result"]["geometric_overlap_m"]) * 1e9 for r in rows],
                marker=["o", "s", "^"][index],
                ls="-" if j == 0 else "--",
                color=C[index],
                label=f"Config. {index + 1}, {mode}",
            )
    title(aa[1, 0], "d  Four bubbles: stiffness sensitivity")
    aa[1, 0].set(xlabel="Stiffness / kc", ylabel="Maximum geometric overlap (nm)", xticks=[0.5, 1, 2])
    aa[1, 0].legend(fontsize=7, ncol=2)
    r = next(r for r in S["four"] if r["configuration"] == 0 and r["mode"] == "B1" and r["multiplier"] == 1)
    for j, (x, label) in enumerate(
        [
            (r["result"]["prediction_m"], "Prediction"),
            (r["pbd_positions_m"], "PBD12"),
            (r["result"]["positions_m"], "Implicit"),
        ]
    ):
        aa[1, 1].plot(range(6), overlaps(np.array(x), 1.78e-6) * 1e9, "o-", color=C[j], label=label)
    title(aa[1, 1], "e  Shared-node contact redistribution")
    aa[1, 1].set(
        xticks=range(6),
        xticklabels=["0–1", "0–2", "0–3", "1–2", "1–3", "2–3"],
        xlabel="Bubble ID pair, configuration 1 (B1)",
        ylabel="Geometric overlap (nm)",
    )
    aa[1, 1].legend()
    errors = [
        np.max(np.abs(np.array(r["result"]["positions_m"]) - r["reference"]["positions_m"])) * 1e12
        for r in S["four"]
    ]
    aa[1, 2].bar(range(len(errors)), errors, color=C[2])
    title(aa[1, 2], "f  Independent four-bubble force reference")
    aa[1, 2].set(
        xlabel="Three configurations × two mobilities × three kc",
        ylabel="Max. position difference (pm)",
        xticks=[0, 5, 11, 17],
        xticklabels=["1", "6", "12", "18"],
    )
    for ax in aa.flat:
        grid(ax)
    k = js(ROOT / "configs/fast_pzt/f10g/stiffness.json")["kc_N_m"]
    finish(
        fig,
        "Fig_F10G_2_implicit_contact_reference",
        f"kc = {k:.7f} N/m, selected from saved compression forces / 8.9 nm before contact trials; not a measured SonoVue stiffness.\n90 two-bubble cases plus ID swaps and hard-limit check; 18 four-bubble cases. B1: free Stokes; B2: inverse of the original lubrication resistance.",
    )


def spheres(ax, x, origin, lim):
    x = (np.array(x) - origin) * 1e6
    radius = 1.78
    u = np.linspace(0, 2 * np.pi, 18)
    v = np.linspace(0, np.pi, 13)
    sx = np.outer(np.cos(u), np.sin(v))
    sy = np.outer(np.sin(u), np.sin(v))
    sz = np.outer(np.ones_like(u), np.cos(v))
    for i, j in shape(x * 1e-6, 1.78e-6)["edges"]:
        ax.plot(*x[[i, j]].T, color="#485566", alpha=0.55, lw=0.7)
    for i, p in enumerate(x):
        ax.plot_surface(
            p[0] + radius * sx,
            p[1] + radius * sy,
            p[2] + radius * sz,
            color=plt.cm.viridis(i / 31),
            linewidth=0,
            antialiased=True,
            shade=True,
            alpha=0.95,
        )
    # Fixed identity colors for every sphere, explicit IDs on a reproducible subset.
    for i in [0, 7, 15, 23, 28, 31]:
        ax.text(*(x[i] + [0, 0, 1.9]), str(i), fontsize=7, color="#202733")
    ax.set(
        xlim=(-lim, lim),
        ylim=(-lim, lim),
        zlim=(-lim, lim),
        xlabel="x (µm)",
        ylabel="y (µm)",
        zlabel="z (µm)",
    )
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=23, azim=-53)
    ax.tick_params(labelsize=8, pad=0)


def dense():
    fig = plt.figure(figsize=(15.5, 15), layout="constrained")
    gs = fig.add_gridspec(3, 3, height_ratios=[1.15, 1.15, 1])
    fig.get_layout_engine().set(hspace=0.10)
    for j, case in enumerate(CASES):
        _, t, _ = state(case)
        x = np.array(t["initial_positions_m"])
        origin = x.mean(0)
        soft = next(r for r in P if r["case"] == case and r["mode"] == "B2" and r["multiplier"] == 2)
        data = [x, np.array(t["proposal_positions_m"]), np.array(soft["result"]["positions_m"])]
        lim = max(np.abs(p - origin).max() for p in data) * 1e6 + 2.5
        for i, (xx, label) in enumerate(
            zip(
                data,
                ["Saved accepted input", "Original PBD: rejected proposal", "B2 / 2kc: static budget pass"],
            )
        ):
            ax = fig.add_subplot(gs[j, i], projection="3d")
            spheres(ax, xx, origin, lim)
            sh = shape(xx, 1.78e-6)
            title(
                ax,
                f"{chr(97 + j * 3 + i)}  {LABEL[j]} — {label}\nRg {sh['Rg_m'] * 1e6:.3f} µm; aspect {sh['aspect_ratio']:.3f}; overlap {overlaps(xx, 1.78e-6).max() * 1e9:.3f} nm",
            )
    ax = fig.add_subplot(gs[2, 0])
    title(ax, "g  Common accepted prefix: internal positions")
    for j, m in enumerate(M):
        s = m["series"]
        t = np.array(s["time_s"]) * 1e3
        ax.plot(t, np.array(s["fixed_ID_max_m"]) * 1e6, color=C[j], label=LABEL[j] + " max ID")
        ax.plot(t, np.array(s["Kabsch_RMSD_m"]) * 1e6, "--", color=C[j], label=LABEL[j] + " RMSD")
    ax.set(xlabel="Time (ms)", ylabel="Difference from original PBD (µm)")
    ax.legend(fontsize=7)
    grid(ax)
    ax = fig.add_subplot(gs[2, 1])
    title(ax, "h  Rg and aspect: same accepted prefix")
    for j, m in enumerate(M):
        s = m["series"]
        t = np.array(s["time_s"]) * 1e3
        for key, ls in [("Rg", "-"), ("aspect", "--")]:
            ref = np.array(s[key + "_reference_m" if key == "Rg" else key + "_reference"])
            val = np.array(s[key + "_candidate_m" if key == "Rg" else key + "_candidate"])
            ax.plot(t, 100 * (val / ref - 1), ls, color=C[j], label=LABEL[j] + " " + key)
    ax.set(xlabel="Time (ms)", ylabel="Relative change from PBD (%)")
    ax.legend(fontsize=7)
    grid(ax)
    ax = fig.add_subplot(gs[2, 2])
    title(ax, "i  New accepted short windows")
    for j, case in enumerate(CASES):
        v = gz(BASE / f"integration_v1/{case}_B2_K2/SUBSTEP_LEDGER.json.gz")["steps"]
        ax.plot(
            [r["attempted_end_s"] * 1e3 for r in v],
            [r["final_overlap_m"] * 1e9 for r in v],
            color=C[j],
            label=LABEL[j],
        )
    ax.axhline(17.8, color="#a84246", ls="--", label="17.8 nm budget")
    ax.set(xlabel="Time (ms)", ylabel="Maximum accepted overlap (nm)", ylim=(0, 20))
    ax.legend(fontsize=8)
    grid(ax)
    finish(
        fig,
        "Fig_F10G_3_dense_cluster_contact",
        "True radius = 1.78 µm; fixed IDs and viridis colors (0 = purple, 31 = yellow). Equal XYZ scales; lines link centers at d ≤ 2.1R.\nUpper panels: same saved initial state and force, one-step proposals; these are not a time sequence. Lower panels: actual saved trajectories.\nOld PBD has no accepted reference beyond 4.85 ms (O) / 2.20 ms (preparation). Different-model agreement is not physical truth or new-model time convergence.",
    )


def admission():
    fig, aa = plt.subplots(2, 3, figsize=(15, 9), layout="constrained")
    for j, case in enumerate(CASES):
        v = [r for r in K if r["case"] == case and r["scope"] == "lubrication_balance_prediction_PBD12"]
        aa[0, 0].plot(
            [r["iterations"] for r in v],
            [r["gpu_median_s"] * 1e6 for r in v],
            "o-",
            color=C[j],
            label=LABEL[j],
        )
        vals = [
            next(r for r in P if r["case"] == case and r["mode"] == mode and r["multiplier"] == 2)["timing"][
                "median_s"
            ]
            * 1e3
            for mode in ["B1", "B2"]
        ]
        aa[0, 1].bar(np.arange(2) + (j - 0.5) * 0.32, vals, width=0.32, color=C[j], label=LABEL[j])
    title(aa[0, 0], "a  Original GPU kernels, fixed states")
    aa[0, 0].set(xlabel="Jacobi passes + balance + PBD12", ylabel="CUDA event span (µs)", xticks=[2, 4, 8])
    aa[0, 0].legend(fontsize=8)
    title(aa[0, 1], "b  CPU implicit initialization + contact")
    aa[0, 1].set(xticks=[0, 1], xticklabels=["B1 / 2kc", "B2 / 2kc"], ylabel="Median time (ms), A supplied")
    aa[0, 1].legend(fontsize=8)
    ax = aa[0, 2]
    title(ax, "c  Short-window diagnostic cost")
    labels = ["O: 5 ms", "Prep.: 2.35 ms"]
    bottom = np.zeros(2)
    for key, color, label in [
        ("integration_quality_wall_s", C[0], "Integration + checks"),
        ("observer_output_s", C[1], "Detailed evidence output"),
        ("model_prepare_s", C[2], "Preparation"),
    ]:
        vals = [r["timing"][key] for r in COST["trajectory_costs"]["rows"]]
        ax.bar(range(2), vals, bottom=bottom, color=color, label=label)
        bottom += vals
    ax.set(xticks=range(2), xticklabels=labels, ylabel="Wall time (s)")
    ax.legend(fontsize=8)
    ax = aa[1, 0]
    title(ax, "d  Two fixed-state contact budgets")
    labels = ["B1 / 0.5kc", "B1 / kc", "B1 / 2kc", "B2 / 0.5kc", "B2 / kc", "B2 / 2kc"]
    for j, case in enumerate(CASES):
        vals = [max(r["result"]["geometric_overlap_m"]) * 1e9 for r in P if r["case"] == case]
        ax.plot(range(6), vals, "o-", color=C[j], label=LABEL[j])
    ax.axhline(17.8, color="#a84246", ls="--")
    ax.set(xticks=range(6), xticklabels=labels, ylabel="Maximum geometric overlap (nm)")
    ax.tick_params(axis="x", labelrotation=35, labelsize=8)
    ax.legend(fontsize=8)
    ax = aa[1, 1]
    ax.axis("off")
    title(ax, "e  Admission evidence")
    entries = [
        ("Discrete lubrication reference", "PASS"),
        ("Two / four bubble contact", "PASS"),
        ("B2 / 2kc: both 32-state checks", "PASS"),
        ("B2 / 2kc: two short windows", "PASS"),
        ("Same-model dt convergence", "NOT RUN"),
        ("Production GPU soft-contact cost", "NOT RUN"),
        ("Exact geometric spring law", "NO"),
        ("Calibrated shell-contact law", "NO"),
        ("Production / AI admission", "NO"),
    ]
    table = ax.table(
        cellText=entries, colLabels=["Check", "Status"], cellLoc="left", loc="center", colWidths=[0.76, 0.24]
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.65)
    for (i, j), cell in table.get_celld().items():
        cell.set_edgecolor("#d9dde2")
        cell.set_facecolor("#ecf1f4" if i == 0 else "white")
    ax = aa[1, 2]
    ax.axis("off")
    title(ax, "f  Scope and next admission gate")
    text = "Two short windows completed\n735 accepted steps; no contact rejection\n\nOld / new morphology differs:\nmax. fixed-ID difference 3.657 / 2.515 µm\ncentroid difference only 1.566 / 0.259 nm\n\nNo new 10 ms or OFF/S/O triplet\nNo production default change\nNo claim of experimental validation\n\nNext: contact-normal / law refinement,\nsame-model dt convergence,\nthen bounded matrix-free GPU feasibility."
    ax.text(0.02, 0.97, text, va="top", fontsize=10, linespacing=1.5)
    for ax in [*aa[0], aa[1, 0]]:
        grid(ax)
    finish(
        fig,
        "Fig_F10G_4_accuracy_cost_admission",
        "Fixed-state kernel spans and CPU reference costs have different scopes: neither is a full production trajectory benchmark.\nDetailed observer output is charged. All candidates remain EXPLORATORY / PHYSICS_UNCALIBRATED; training_ready = false.",
    )


def network():
    fig, aa = plt.subplots(2, 2, figsize=(12, 10), layout="constrained")
    colors = ListedColormap(["#ffffff", "#9aabba", "#526ba2", "#b56a70"])
    for j, case in enumerate(CASES):
        with (
            h5py.File(ROOT / f"results/fast_pzt/f10e/contact_v1/{case}_C0/trajectory.h5") as old,
            h5py.File(BASE / f"integration_v1/{case}_B2_K2/trajectory.h5") as new,
        ):
            nt = min(len(old["time_s"]), len(new["time_s"]))
            time = float(old["time_s"][nt - 1])
            x = old["positions_m"][nt - 1, 0]
            y = new["positions_m"][nt - 1, 0]

        def adjacency(x):
            out = np.linalg.norm(x[:, None] - x[None, :], axis=-1) <= 2.1 * 1.78e-6
            np.fill_diagonal(out, False)
            return out

        a, b = adjacency(x), adjacency(y)
        aa[j, 0].imshow(
            a,
            cmap=ListedColormap(["white", "#526ba2"]),
            vmin=0,
            vmax=1,
            interpolation="nearest",
            origin="lower",
        )
        changes = np.zeros_like(a, dtype=int)
        changes[a & b] = 1
        changes[a & ~b] = 2
        changes[b & ~a] = 3
        im = aa[j, 1].imshow(
            changes, cmap=colors, vmin=-0.5, vmax=3.5, interpolation="nearest", origin="lower"
        )
        title(aa[j, 0], f"{LABEL[j]}: old PBD at {time * 1e3:g} ms\n{a.sum() // 2} undirected edges")
        title(
            aa[j, 1],
            f"Same time: PBD / B2 edge identities\nShared {(a & b).sum() // 2}; PBD only {(a & ~b).sum() // 2}; B2 only {(b & ~a).sum() // 2}",
        )
        for ax in aa[j]:
            ax.set(
                xlabel="Fixed bubble ID",
                ylabel="Fixed bubble ID",
                xticks=[0, 7, 15, 23, 31],
                yticks=[0, 7, 15, 23, 31],
            )
    bar = fig.colorbar(im, ax=aa[:, 1], ticks=[0, 1, 2, 3], shrink=0.75, pad=0.025)
    bar.ax.set_yticklabels(["No edge", "Shared", "PBD only", "B2 only"])
    finish(
        fig,
        "Fig_F10G_3_dense_cluster_contact_network",
        "Supplement to Figure 3. Exact adjacency at the common accepted cutoff; d ≤ 2.1R, R = 1.78 µm, diagonal excluded.\nMatrices show bubble identities, not resized bubble glyphs. An edge change is geometric connectivity, not physical coalescence or breakup.",
    )


def main():
    lubrication()
    implicit()
    dense()
    network()
    admission()
    inputs = [
        BASE / "lubrication_v1/RESULTS.json",
        BASE / "prototype_v1/SMALL_SYSTEMS.json",
        BASE / "prototype_v1/DENSE_STATIC.json",
        BASE / "performance_v2/RESULTS.json",
        BASE / "analysis/MORPHOLOGY.json",
        BASE / "analysis/TRAJECTORY_AUDIT.json",
        BASE / "analysis/COST_AND_YIELD.json",
    ]
    inputs += [
        BASE / f"integration_v1/{case}_B2_K2/{name}"
        for case in CASES
        for name in ["SUBSTEP_LEDGER.json.gz", "trajectory.h5"]
    ]
    inputs += [ROOT / f"results/fast_pzt/f10e/contact_v1/{case}_C0/trajectory.h5" for case in CASES]
    write(
        OUT / "FIGURE_PROVENANCE.json",
        dict(
            source=identity(),
            data_sha256={str(p.relative_to(ROOT)): sha(p) for p in inputs},
            figures={p.name: sha(p) for p in OUT.glob("Fig_F10G*")},
            radius_m=1.78e-6,
            no_synthetic_frames=True,
        ),
    )
    print("Four figure groups saved")


if __name__ == "__main__":
    main()
