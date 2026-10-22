"""Stage 2G: published N=2 reference trends; never fit source parameters."""

from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
import multiprocessing
from pathlib import Path
import subprocess
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mevo_mb.validation.zhang2011 import (
    Zhang2011Parameters,
    Zhang2011EncapsulatedPairReference,
    secondary_coefficient,
)
from mevo_mb.validation.doinikov_bouakaz_pair import (
    DoinikovBouakazParameters,
    DoinikovBouakazPairReferenceModel,
)


def calculate(job):
    kind, key, parameters, settings, directory = job
    out = Path(directory) / key
    out.mkdir()
    start = time.perf_counter()
    meta = dict(
        parameters=asdict(parameters),
        settings=settings,
        status="FAILED",
        physics_status="partially validated",
    )
    try:
        cfg = settings.copy()
        if kind == "zhang":
            model = Zhang2011EncapsulatedPairReference(parameters)
        else:
            model = DoinikovBouakazPairReferenceModel(parameters, cfg.pop("order"))
        z = model.solve(**cfg)
        coefficient = secondary_coefficient(z, parameters.density_kg_m3, frequency=parameters.frequency_Hz)
        sample = settings.get("samples_per_cycle", 64 if kind == "zhang" else 48)
        last = z["distance_m"][-sample:]
        radial_fraction = np.max(
            np.abs(z["radius_m"] / np.array([parameters.radius1_m, parameters.radius2_m]) - 1)
        )
        meta.update(
            status="CONTACT_STOP" if z["contact_reached"] else "COMPLETE",
            contact_time_s=z["contact_time_s"],
            nfev=z["nfev"],
            metrics=dict(
                secondary_coefficient_N_m2=coefficient,
                final_cycle_distance_m=float(np.mean(last)),
                distance_change_m=float(np.mean(last) - parameters.distance_m),
                min_gap_m=float(np.min(z["distance_m"] - z["radius_m"].sum(axis=1))),
                max_radial_fraction=float(radial_fraction),
                min_distance_over_sum_radii=float(np.min(z["distance_m"] / z["radius_m"].sum(axis=1))),
            ),
        )
        np.savez_compressed(
            out / "trajectory.npz", **{k: v for k, v in z.items() if isinstance(v, np.ndarray)}
        )
    except Exception as error:
        meta["error"] = repr(error)
    meta["runtime_s"] = time.perf_counter() - start
    (out / "RESULT.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(kind, key, meta["status"], meta.get("metrics", {}).get("distance_change_m"), flush=True)
    return kind, key, meta


def main():
    root = Path.cwd()
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    dirs = {
        k: root / "results/stage2" / v / stamp
        for k, v in (("zhang", "pair_zhang2011"), ("doinikov", "pair_doinikov2015"))
    }
    for out in dirs.values():
        out.mkdir(parents=True)
    jobs = []
    for distance in (30.0, 40.0, 50.0, 70.0, 100.0, 120.0):
        p = replace(Zhang2011Parameters(), distance_m=distance * 1e-6)
        jobs.append(("zhang", f"distance_{distance:g}", p, {"cycles": 64}, dirs["zhang"]))
    for r in (1.0, 1.5, 2.0, 2.5, 2.7, 2.8, 2.9, 3.0, 3.1, 3.2, 3.5, 4.0):
        p = replace(Zhang2011Parameters(), radius1_m=r * 1e-6, radius2_m=1.5e-6, distance_m=40e-6)
        jobs.append(("zhang", f"radius_{r:g}", p, {"cycles": 64}, dirs["zhang"]))
    for r in (3.0, 4.0):
        p = replace(Zhang2011Parameters(), radius1_m=r * 1e-6)
        jobs.append(("zhang", f"translation_{r:g}", p, {"cycles": 450}, dirs["zhang"]))
    jobs.append(
        (
            "zhang",
            "distance_50_fine",
            Zhang2011Parameters(),
            {"cycles": 64, "samples_per_cycle": 128, "max_steps_per_cycle": 96, "rtol": 2e-10},
            dirs["zhang"],
        )
    )
    for f in (2.15, 2.95):
        p = replace(DoinikovBouakazParameters(), frequency_Hz=f * 1e6)
        for order in (4, 8):
            jobs.append(("doinikov", f"weak_f{f:g}_m{order}", p, {"order": order}, dirs["doinikov"]))
    p = replace(
        DoinikovBouakazParameters(), distance_m=10.5e-6, frequency_Hz=1.7e6, pressure_Pa=1e5, pulse_cycles=30
    )
    for order in (2, 4, 8, 12):
        jobs.append(("doinikov", f"near_m{order}", p, {"order": order}, dirs["doinikov"]))
    start = time.perf_counter()
    results = {"zhang": {}, "doinikov": {}}
    with ProcessPoolExecutor(max_workers=6, mp_context=multiprocessing.get_context("spawn")) as pool:
        for kind, key, meta in pool.map(calculate, jobs):
            results[kind][key] = meta
    for kind, out in dirs.items():
        summary = dict(
            run=str(out.relative_to(root)),
            cases=results[kind],
            elapsed_s=time.perf_counter() - start,
            git_commit=subprocess.check_output(
                ["git", "-c", f"safe.directory={root}", "rev-parse", "HEAD"], text=True
            ).strip(),
            hardware="results/stage2/system/SERVER_DOCTOR.json",
            workers=6,
            random_seed=None,
            source_sha256={
                p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (root / "src/mevo_mb/validation").glob("*.py")
            },
            physics_status="partially validated",
        )
        if kind == "zhang":
            a = np.load(out / "distance_50/trajectory.npz")
            b = np.load(out / "distance_50_fine/trajectory.npz")
            error = float(np.max(np.abs(a["distance_m"] - b["distance_m"][::2])))
            r = results[kind]
            signs = (
                r["distance_50"]["metrics"]["secondary_coefficient_N_m2"]
                < 0
                < r["distance_120"]["metrics"]["secondary_coefficient_N_m2"]
            )
            translation = (
                r["translation_4"]["metrics"]["distance_change_m"]
                < 0
                < r["translation_3"]["metrics"]["distance_change_m"]
            )
            peak = max(
                (k for k in r if k.startswith("radius_")),
                key=lambda k: r[k]["metrics"]["secondary_coefficient_N_m2"],
            )
            summary.update(
                status="PASS_REPRESENTATIVE_TRENDS"
                if signs and translation and error < 1e-8
                else "FAILED_TREND_GATE",
                timestep_max_distance_error_m=error,
                fig3_sign_inversion=signs,
                fig4_initial_translation_signs=translation,
                fig5_peak_sample=peak,
                limitations=[
                    "64-cycle coefficient uses last eight cycles; paper averaging window unspecified",
                    "Fig4 only first 0.2 ms, no claim of final stable-pair reproduction",
                    "Chi=1 N/m; elastic shell without buckling/rupture; not SonoVue calibration",
                ],
            )
        else:
            a = np.load(out / "near_m8/trajectory.npz")
            b = np.load(out / "near_m12/trajectory.npz")
            n = min(len(a["time_s"]), len(b["time_s"]))
            error = float(np.max(np.abs(a["distance_m"][:n] - b["distance_m"][:n])))
            signs = (
                results[kind]["weak_f2.15_m8"]["metrics"]["distance_change_m"]
                < 0
                < results[kind]["weak_f2.95_m8"]["metrics"]["distance_change_m"]
            )
            summary.update(
                status="PASS_LIMITED_REFERENCE" if signs and error < 1e-8 else "FAILED_REFERENCE_GATE",
                multipole_m8_m12_max_distance_error_m=error,
                weak_field_signs=signs,
                limitations=[
                    "Author Mathematica supplement inaccessible; independently evaluated finite Legendre boundary system",
                    "Gaussian full duration=6 sigma assumed; Fig3 shortened to 60 cycles, not quantitative curve reproduction",
                    "Uncoated gas bubbles, spherical, potential flow; no boundary-layer overlap or contact physics",
                    "Multipole convergence is numerical evidence, not experimental validation",
                ],
            )
        (out / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
        (out.parent / "STATE.json").write_text(json.dumps(summary, indent=2) + "\n")
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    r = results["zhang"]
    out = dirs["zhang"]
    ds = [30, 40, 50, 70, 100, 120]
    ax[0, 0].plot(ds, [r[f"distance_{d}"]["metrics"]["secondary_coefficient_N_m2"] for d in ds], "o-")
    ax[0, 0].axhline(0, color="grey", lw=0.8)
    ax[0, 0].set(xlabel="Initial distance (um)", ylabel="f_B (N m²)", title="Zhang Fig3: interaction sign")
    rs = [1, 1.5, 2, 2.5, 2.7, 2.8, 2.9, 3, 3.1, 3.2, 3.5, 4]
    ax[0, 1].plot(rs, [r[f"radius_{v:g}"]["metrics"]["secondary_coefficient_N_m2"] for v in rs], "o-")
    ax[0, 1].axhline(0, color="grey", lw=0.8)
    ax[0, 1].set(xlabel="Bubble 1 radius (um)", ylabel="f_B (N m²)", title="Zhang Fig5: resonance region")
    for v in (3, 4):
        z = np.load(out / f"translation_{v}/trajectory.npz")
        ax[1, 0].plot(z["time_s"] * 1e6, z["distance_m"] * 1e6, label=f"R1={v} um, R2=6 um")
    ax[1, 0].set(xlabel="Time (us)", ylabel="Center distance (um)", title="Zhang Fig4: initial segment only")
    ax[1, 0].legend()
    for order in (2, 4, 8, 12):
        z = np.load(dirs["doinikov"] / f"near_m{order}/trajectory.npz")
        ax[1, 1].plot(z["time_s"] * 1e6, z["distance_m"] * 1e6, label=f"Legendre M={order}")
    ax[1, 1].set(
        xlabel="Time (us)", ylabel="Center distance (um)", title="Doinikov 2015: assumed Gaussian width"
    )
    ax[1, 1].legend()
    fig.savefig(out / "Figure_C_pair.png", dpi=220)
    plt.close(fig)
    print(
        json.dumps(
            {
                k: {"path": str(v), "status": json.loads((v / "SUMMARY.json").read_text())["status"]}
                for k, v in dirs.items()
            }
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
