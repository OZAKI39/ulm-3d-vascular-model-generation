"""Summarize completed evidence and benchmark batch queries without new physics runs."""

from collections import Counter
import argparse
import hashlib
import h5py
import json
from pathlib import Path
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
import numpy as np

from mevo_mb.bubbles.adaptive_table import AdaptiveResponseTableBubbleModel
from mevo_mb.validation.regime_maps import PairPhysicsRegimeMap, REGIMES
from common import new_run, finish, atomic_json, watermark


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--figures-only", action="store_true")
    args = parser.parse_args()
    out = new_run(
        "audit",
        {"task": "final evidence analysis", "seed": 250025, "figures_only": args.figures_only},
        ["Existing Stage2.5 benchmark data; no new reference ODE/PDE"],
    )
    root = Path("results/stage2_5")
    table = json.loads((root / "response_table/STATE.json").read_text())
    table_run = Path(table["run"])
    tree = json.loads((table_run / "CHECKPOINT.json").read_text())
    hold = json.loads((table_run / "HOLDOUTS.json").read_text())["independent_uniform_queries"]
    good = [r for r in hold if r["status"] == "PASS"]
    error_max = {k: max(r["errors"][k] for r in good) for k in good[0]["errors"]}
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    axes[0].hist([r["depth"] for r in tree["leaves"]], bins=np.arange(22) - 0.5)
    axes[0].set(
        xlabel="Binary subdivision depth",
        ylabel="Valid leaf count",
        title=f"{len(tree['leaves'])} leaves, full declared domain",
    )
    for k in error_max:
        axes[1].plot([r["errors"][k] for r in good], ".", ms=3, label=k)
    axes[1].axhline(0.05, c="k", ls="--", lw=1)
    axes[1].set(xlabel="Independent holdout query", ylabel="Relative error (phase: radians)", ylim=(0, 0.055))
    axes[1].legend(fontsize=6)
    watermark(fig, "PARTIALLY VALIDATED — NUMERICAL INTERPOLATION, NOT EXPERIMENTAL CALIBRATION")
    fig.savefig(out / "FINAL_ADAPTIVE_VALIDATION.png", dpi=180)
    plt.close(fig)
    pair = json.loads((root / "pair_asymptotics/STATE.json").read_text())
    model = PairPhysicsRegimeMap(root / "pair_asymptotics/STATE.json")
    physical = [r for r in pair["rows"] if r["model"] == "zhang" and r["sound_speed_m_s"] == 1482.0]
    records = []
    for r in physical:
        status = model.classify(
            r["radius1_m"], r["radius2_m"], r["separation_m"], r["frequency_Hz"], r["pressure_Pa"], "uncoated"
        )
        records.append(
            dict(
                frequency_Hz=r["frequency_Hz"],
                pressure_Pa=r["pressure_Pa"],
                ratio=r["separation_m"] / r["radius2_m"],
                status=status,
            )
        )
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), constrained_layout=True)
    freqs = sorted({r["frequency_Hz"] for r in records})
    ratios = [10, 20, 40, 80, 160]
    palette = ListedColormap(["#55aa66", "#eecc66", "#dd7755", "#775599", "#dddddd"])
    norm = BoundaryNorm(np.arange(6) - 0.5, 5)
    for ax, p in zip(axes, [100, 500]):
        grid = np.array(
            [
                [
                    REGIMES.index(
                        next(
                            r["status"]
                            for r in records
                            if r["frequency_Hz"] == f and r["pressure_Pa"] == p and np.isclose(r["ratio"], q)
                        )
                    )
                    for q in ratios
                ]
                for f in freqs
            ]
        )
        im = ax.imshow(grid, origin="lower", aspect="auto", cmap=palette, norm=norm)
        ax.set_xticks(range(5), ratios)
        ax.set_yticks(range(len(freqs)), [f / 1e6 for f in freqs])
        ax.set(xlabel="d / max(R1,R2)", ylabel="Frequency (MHz)", title=f"Uncoated R1=1.5, R2=2 um; {p} Pa")
    cb = fig.colorbar(im, ax=axes, ticks=range(5), shrink=0.7)
    cb.ax.set_yticklabels(
        ["Safe", "Transition", "Reference required", "Contact", "Outside sampled validity"], fontsize=8
    )
    watermark(fig, "PARTIALLY VALIDATED — SAMPLED UNCOATED PAIRS; NO SONOVUE CERTIFICATE")
    fig.savefig(out / "PAIR_REGIME_MAP.png", dpi=180)
    plt.close(fig)
    atomic_json(
        out / "PAIR_REGIME_MAP.json",
        dict(
            records=records,
            counts=dict(Counter(r["status"] for r in records)),
            scope="Exact sampled nodes only; kd<=0.3, radial excursion<=2%, separation change<=1%",
        ),
    )
    # Replot saved failed-case arrays only. Original HDF5/figures remain intact.
    zhang = json.loads((root / "zhang_forensic/STATE.json").read_text())
    source = Path(zhang["run"]) / "ZHANG_FIG3_FAILED_CASE.h5"
    with h5py.File(source) as h:
        data = {k: h[k][:] for k in h.keys()}
        parameters = json.loads(h.attrs["parameters_json"])
    t = data["t"]
    f = parameters["f"]
    last = t * f >= 253
    tx = (t[last] - t[last][0]) * 1e6
    fig, ax = plt.subplots(4, 2, figsize=(12, 12), constrained_layout=True)
    r = data["R"] / np.array([parameters["R10"], parameters["R20"]])
    for j in range(2):
        ax[0, j].plot(tx, r[last, j])
        ax[0, j].set(xlabel="Time in final three cycles (us)", ylabel=f"R{j + 1}/R{j + 1}0")
    ax[1, 0].plot(tx, data["drive_Pa"][last] / 1e3)
    ax[1, 0].set(xlabel="Time in final three cycles (us)", ylabel="Drive (kPa)")
    ax[1, 1].plot(r[last, 0] - 1, r[last, 1] - 1)
    ax[1, 1].set(xlabel="R1/R10 - 1", ylabel="R2/R20 - 1", title="Radial phase relation")
    for j, key, label in [
        (0, "Vdot_product", "Vdot1 Vdot2 (m6/s2)"),
        (1, "force_integrand_N", "Force integrand (N)"),
    ]:
        ax[2, j].plot(tx, data[key][last])
        ax[2, j].set(xlabel="Time in final three cycles (us)", ylabel=label)
    cf = data["cycle_average_coefficient_N_m2"]
    ax[3, 0].plot(np.arange(len(cf)) + 0.5, cf / 1e-18, label="Per cycle")
    ax[3, 0].plot(t * f, data["accumulated_average_N_m2"] / 1e-18, label="From initial rest")
    ax[3, 0].axhline(6.6, color="k", ls=":", label="Approximate Fig.3 digitization")
    ax[3, 0].legend(fontsize=8)
    ax[3, 0].set(xlabel="Elapsed cycles", ylabel="fB / 1e-18 N m²")
    ax[3, 1].plot(t * 1e6, data["separation"] * 1e6)
    ax[3, 1].set(xlabel="Time (us)", ylabel="Separation (um)")
    fig.suptitle("Zhang Fig.3: R1=3 um, R2=6 um, D0=120 um; 100 kPa, 2.25 MHz")
    watermark(fig, "PARTIALLY VALIDATED — FIGURE 3 SIGN UNRESOLVED")
    fig.savefig(out / "ZHANG_FIG3_DEBUG.png", dpi=170)
    plt.close(fig)
    atomic_json(
        out / "FIGURE_SOURCES.json",
        dict(
            zhang_h5=str(source),
            zhang_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            table_state=str(root / "response_table/STATE.json"),
            pair_state=str(root / "pair_asymptotics/STATE.json"),
        ),
    )
    if args.figures_only:
        finish(
            out,
            dict(
                status="FIGURES_REFRESHED",
                holdout_error_max=error_max,
                source_zhang=str(source),
                scientific_data_modified=False,
            ),
        )
        return
    # Both processors query exactly the same points; warm-up is excluded.
    import cupy as cp

    bounds = np.array(tree["config"]["bounds"])
    rng = np.random.default_rng(250025)
    cpu = AdaptiveResponseTableBubbleModel(table["table"])
    gpu = AdaptiveResponseTableBubbleModel(table["table"], backend="cupy")
    timing = []
    for n in [128, 1000, 10000]:
        q = rng.uniform(bounds[:, 0], bounds[:, 1], (n, 3))
        device = cp.asarray(q)
        cpu.query(*q.T)
        gpu.query(*device.T)
        cp.cuda.Stream.null.synchronize()
        values = {}
        times = {}
        for label, m, points in [("numpy", cpu, q), ("cupy", gpu, device)]:
            samples = []
            for _ in range(3):
                cp.cuda.Stream.null.synchronize()
                start = time.perf_counter()
                v = m.query(*points.T)
                cp.cuda.Stream.null.synchronize()
                samples.append(time.perf_counter() - start)
            values[label] = cp.asnumpy(v["harmonics_m3"]) if label == "cupy" else v["harmonics_m3"]
            times[label] = float(np.median(samples))
        error = float(np.linalg.norm(values["numpy"] - values["cupy"]) / np.linalg.norm(values["numpy"]))
        timing.append(
            dict(
                N=n,
                median_seconds=times,
                relative_L2=error,
                gpu_pool_active_bytes=cp.get_default_memory_pool().used_bytes(),
                scope="Warm query kernels with inputs resident; excludes table load and host transfer",
            )
        )
    atomic_json(
        out / "BATCH_INTERPOLATION_BENCHMARK.json",
        dict(rows=timing, status="PASS" if all(r["relative_L2"] < 1e-12 for r in timing) else "FAIL"),
    )
    finish(
        out,
        dict(
            status="ANALYSIS_COMPLETE",
            adaptive_table_run=table["run"],
            holdout_error_max=error_max,
            pair_regime_counts=dict(Counter(r["status"] for r in records)),
            batch_interpolation=timing,
        ),
    )


if __name__ == "__main__":
    main()
