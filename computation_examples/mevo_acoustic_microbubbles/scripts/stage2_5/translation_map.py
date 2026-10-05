"""Matched-force translation approximation map; thick shell is NOT SonoVue."""

from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import json
from pathlib import Path
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mevo_mb.validation.doinikov_dayton import DoinikovDaytonParameters
from common import resumable_run, finish, atomic_json, watermark
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "stage2"))
from translation_validate import calculate


def resumed_calculate(job):
    key, parameters, directory, fine = job
    result_path = Path(directory) / key / "RESULT.json"
    if result_path.exists():
        result = json.loads(result_path.read_text())
        if result["parameters"] != asdict(parameters):
            raise ValueError("Translation resume parameter mismatch")
        return key, result
    if result_path.parent.exists():
        attempts = Path(directory) / "resume_attempts"
        attempts.mkdir(exist_ok=True)
        folder = attempts / f"{key}_{len(list(attempts.glob(key + '_*'))):04d}"
        folder.mkdir()
        name, result = calculate((key, parameters, folder, fine))
        result["reused_from"] = str(folder / key)
        atomic_json(result_path, result)
        return name, result
    return calculate(job)


def main():
    start = time.perf_counter()
    old = json.loads(Path("results/stage2/translation/STATE.json").read_text())
    out = resumable_run(
        "translation_map",
        dict(
            radii_um=[1.5, 2, 3],
            pressures_kPa=[5, 40, 150],
            frequencies_MHz=[1, 2.5, 4],
            cycles=20,
            shell="DoinikovDayton2006 15nm solid shell",
            safe_error=0.05,
            approximate_error=0.20,
        ),
        ["10.1121/1.2215228 Eqs6,31-33,36,40,63"],
    )
    jobs = []
    rows = []
    for r in [1.5, 2.0, 3.0]:
        for p in [5.0, 40.0, 150.0]:
            for f in [1.0, 2.5, 4.0]:
                params = DoinikovDaytonParameters(
                    outer_radius_m=r * 1e-6, pressure_Pa=p * 1e3, frequency_Hz=f * 1e6
                )
                name = f"r{r:g}_p{p:g}_f{f:g}"
                cached = next(
                    ((k, v) for k, v in old["cases"].items() if v["parameters"] == asdict(params)), None
                )
                if cached:
                    k, v = cached
                    rows.append((name, dict(v, reused_from=old["run"] + "/" + k)))
                else:
                    jobs.append((name, params, out, False))
    with ProcessPoolExecutor(max_workers=6) as pool:
        rows.extend(pool.map(resumed_calculate, jobs))
    data = []
    for name, result in rows:
        if result["status"] != "COMPLETE":
            data.append(dict(case=name, status="FAILED"))
            continue
        m = result["metrics"]
        p = result["parameters"]
        den = max(abs(m["displacement_m"]), 1e-14)
        e = abs(m["approximation_error_m"]) / den
        status = (
            "FAST_MODEL_SAFE"
            if e <= 0.05
            else ("FAST_MODEL_APPROXIMATE" if e <= 0.2 else "REFERENCE_REQUIRED")
        )
        data.append(
            dict(
                case=name,
                radius_m=p["outer_radius_m"],
                pressure_Pa=p["pressure_Pa"],
                frequency_Hz=p["frequency_Hz"],
                shell_model="doinikov_dayton_2006_solid_15nm",
                cycles=20,
                error=e,
                status=status,
                metrics=m,
                source=result.get("reused_from", str(out / name)),
            )
        )
    atomic_json(
        out / "TRANSLATION_VALIDITY_MAP.json",
        dict(
            scope="matched cycle force;20-cycle sinusoid;no postpulse;not SonoVue",
            thresholds={"safe": 0.05, "approximate": 0.2},
            rows=data,
            interpolation="exact sampled nodes only; no untested interpolation",
        ),
    )
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.5), constrained_layout=True)
    for ax, r in zip(axes, [1.5, 2.0, 3.0]):
        arr = np.array(
            [
                [
                    next(
                        x["error"]
                        for x in data
                        if np.isclose(x["radius_m"], r * 1e-6, atol=1e-15)
                        and x["pressure_Pa"] == p * 1e3
                        and x["frequency_Hz"] == f * 1e6
                    )
                    for f in [1, 2.5, 4]
                ]
                for p in [5, 40, 150]
            ]
        )
        im = ax.imshow(arr * 100, origin="lower", aspect="auto", vmin=0, vmax=max(20, float(arr.max() * 100)))
        for (i, j), value in np.ndenumerate(arr):
            ax.text(j, i, f"{value * 100:.1f}%", ha="center", va="center", color="white")
        ax.set_xticks(range(3), [1, 2.5, 4])
        ax.set_yticks(range(3), [5, 40, 150])
        ax.set(title=f"R0={r:g} um", xlabel="Frequency (MHz)", ylabel="Pressure (kPa)")
        fig.colorbar(im, ax=ax, label="Displacement error (%)")
    watermark(fig, "PARTIALLY VALIDATED — SOLID-SHELL REFERENCE, NOT SONOVUE")
    fig.savefig(out / "TRANSLATION_ERROR_MAP.png", dpi=180)
    plt.close(fig)
    finish(
        out,
        dict(
            status="MAP_COMPLETE",
            rows=data,
            new_cases=len(jobs),
            reused_cases=27 - len(jobs),
            runtime_s=time.perf_counter() - start,
        ),
    )


if __name__ == "__main__":
    main()
