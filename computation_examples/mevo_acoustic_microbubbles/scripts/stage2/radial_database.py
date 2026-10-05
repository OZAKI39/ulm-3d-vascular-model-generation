"""Stage 2E: parallel coarse APECSS table, convergence and withheld-cell checks."""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import resource
import subprocess
import threading
import time

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import psutil

from mevo_mb.bubbles.apecss import ApecssBubbleResponseModel, RadialCase, APECSS_COMMIT
from mevo_mb.bubbles.parameters import SonoVueParameterDistribution, load_parameter_cases
from mevo_mb.bubbles.response_table import ResponseTableBubbleModel


def worker_count():
    cores = min(psutil.cpu_count(logical=False) or 1, len(os.sched_getaffinity(0)))
    path = Path("/sys/fs/cgroup/cpu.max")
    if path.exists():
        q, period = path.read_text().split()
        if q != "max":
            cores = min(cores, float(q) / float(period))
    return max(1, int(cores * 0.85))


def task(entry):
    key, case, out, exe = entry
    try:
        meta, harmonics, _ = ApecssBubbleResponseModel(exe).solve(case, out)
        return key, meta, harmonics
    except Exception as e:
        return key, {"status": "FAILED", "error": repr(e)}, None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", default="auto")
    parser.add_argument("--profile", choices=("coarse", "focused"), default="coarse")
    args = parser.parse_args()
    count = worker_count() if args.workers == "auto" else int(args.workers)
    if count < 1:
        raise ValueError("workers must be positive")
    root = Path.cwd()
    out = root / "results/stage2/single_bubble" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out.mkdir(parents=True)
    exe = root / "third_party/build/apecss/ultrasound_apecss"
    parameter_path = root / "configs/stage2/sonovue_parameters.json"
    params = SonoVueParameterDistribution(parameter_path)
    radii = np.array([1.0, 1.78, 2.5, 3.0]) * 1e-6
    freqs = np.array([1.0, 2.5, 4.0]) * 1e6
    pressures = np.array([5.0, 40.0, 100.0, 150.0, 200.0]) * 1e3
    if args.profile == "focused":
        radii = np.array([1.6, 1.7, 1.78, 1.86, 1.96]) * 1e-6
        freqs = np.array([2.3, 2.4, 2.5, 2.6, 2.7]) * 1e6
        pressures = np.array([5.0, 10.0, 20.0, 30.0, 40.0]) * 1e3
    axes = [radii, freqs, pressures]
    shape = tuple(len(x) for x in axes)

    def make_case(r, f, p):
        ps = params.parameter_set(float(r))
        return RadialCase(
            float(r),
            float(f),
            float(p),
            ps.shell_elasticity_N_m.value,
            ps.shell_viscosity_kg_s.value,
            0.036,
            provenance=json.dumps(ps.to_dict(), sort_keys=True) + "; sigma_init=.036 N/m assumed",
        )

    jobs = []
    cases = {}
    for idx in np.ndindex(shape):
        key = "node_" + "_".join(map(str, idx))
        case = make_case(*[a[i] for a, i in zip(axes, idx)])
        cases[key] = case
        jobs.append((key, case, out / key, exe))
    cell_shape = tuple(s - 1 for s in shape)
    for idx in np.ndindex(cell_shape):
        key = "mid_" + "_".join(map(str, idx))
        case = make_case(*[(a[i] + a[i + 1]) / 2 for a, i in zip(axes, idx)])
        cases[key] = case
        jobs.append((key, case, out / key, exe))
    # Low/mid/high forcing plus shell-state sensitivity and exact published constants.
    for key in ("node_1_1_0", "node_1_1_2", "node_1_1_4", "node_2_0_4"):
        base = cases[key]
        for suffix, fine in [
            ("fine", replace(base, samples_per_cycle=1024, tolerance=1e-11)),
            ("long", replace(base, cycles=96)),
        ]:
            name = key + "_" + suffix
            cases[name] = fine
            jobs.append((name, fine, out / name, exe))
    fixed = load_parameter_cases(parameter_path)["sonovue_table2"]
    for sigma in (0.0, 0.036, 0.072):
        key = f"fixed_table2_sigma{sigma:g}"
        case = RadialCase(
            fixed.radius_m.value,
            2.5e6,
            150e3,
            0.3,
            3.2e-9,
            sigma,
            provenance=json.dumps(fixed.to_dict()) + "; sigma_init explicit hypothesis",
        )
        cases[key] = case
        jobs.append((key, case, out / key, exe))
    official = RadialCase(
        0.975e-6,
        2.9e6,
        130e3,
        1.0,
        15e-9,
        0.0,
        sound_speed_m_s=1480.0,
        viscosity_Pa_s=0.001,
        sigma_water_N_m=0.073,
        polytropic_exponent=1.095,
        ambient_pressure_Pa=1e5,
        provenance="Unmodified official APECSS ultrasound/lipidcoated_simple parameters",
    )
    jobs.append(("official_lipid", official, out / "official_lipid", exe))
    cases["official_lipid"] = official
    hashes = {
        str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
        for folder in ("src", "scripts/stage2", "configs/stage2")
        for p in (root / folder).rglob("*")
        if p.is_file() and "__pycache__" not in p.parts
    }
    provenance = {
        "git_commit": subprocess.check_output(
            ["git", "-c", f"safe.directory={root}", "rev-parse", "HEAD"], text=True
        ).strip(),
        "working_tree_source_hashes": hashes,
        "apecss_commit": APECSS_COMMIT,
        "random_seed": None,
        "profile": args.profile,
        "workers": count,
        "threads_per_worker": 1,
        "hardware": "results/stage2/system/SERVER_DOCTOR.json",
        "cases": {k: asdict(c) for k, c in cases.items()},
    }
    (out / "CONFIG.json").write_text(json.dumps(provenance, indent=2) + "\n")
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[name] = "1"
    peak = [0]
    stop = threading.Event()

    def monitor():
        while not stop.wait(0.1):
            total = 0
            for proc in [psutil.Process()] + psutil.Process().children(recursive=True):
                try:
                    total += proc.memory_info().rss
                except psutil.Error:
                    pass
            peak[0] = max(peak[0], total)

    thread = threading.Thread(target=monitor, daemon=True)
    thread.start()
    start = time.perf_counter()
    cpu0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    results = {}
    with ProcessPoolExecutor(max_workers=count, mp_context=multiprocessing.get_context("spawn")) as pool:
        for future in as_completed([pool.submit(task, job) for job in jobs]):
            key, meta, harmonics = future.result()
            results[key] = (meta, harmonics)
            print(f"{len(results)}/{len(jobs)} {key}: {meta['status']}", flush=True)
    elapsed = time.perf_counter() - start
    stop.set()
    thread.join()
    cpu1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    cpu_s = cpu1.ru_utime + cpu1.ru_stime - cpu0.ru_utime - cpu0.ru_stime
    values = np.zeros(shape + (8,), dtype=complex)
    node_valid = np.zeros(shape, bool)
    resolved = np.zeros(shape, bool)
    rmax, rmin = np.zeros(shape), np.zeros(shape)
    for idx in np.ndindex(shape):
        meta, h = results["node_" + "_".join(map(str, idx))]
        if h is not None:
            values[idx] = h
            m = meta["metrics"]
            node_valid[idx], resolved[idx] = m["valid"], m["harmonic_energy_resolved"]
            rmax[idx], rmin[idx] = m["max_R_over_R0"], m["min_R_over_R0"]
    table_path = out / "bubble_response_database.h5"
    metadata = {
        "radial_model": "APECSS RPAR + reversible Marmottant",
        "physics_status": "partially validated",
        "apecss_commit": APECSS_COMMIT,
        "parameter_file_sha256": hashlib.sha256(parameter_path.read_bytes()).hexdigest(),
        "parameter_model": "Fig4 digitization; sigma_init=0.036 N/m assumption",
        "parameter_calibration_validated": False,
        "experimental_condition": "2.5 MHz, 150 kPa pulses",
        "medium": {"density_kg_m3": 1000.0, "sound_speed_m_s": 1500.0, "viscosity_Pa_s": 0.002},
        "phase_convention": "Re[P exp(iwt)], V_h for real positive fundamental pressure",
        "interpolation_evidence": "midpoint per cell; 5% complex fundamental and weighted harmonic tolerance, not a certified bound",
    }
    with h5py.File(table_path, "x") as h5:
        h5.attrs.update(schema_version=1, complete=True, metadata_json=json.dumps(metadata))
        for name, a in zip(("radius_m", "frequency_Hz", "pressure_Pa"), axes):
            h5.create_dataset(name, data=a)
        for name, a in dict(
            harmonics_m3=values,
            node_valid=node_valid,
            harmonic_resolved=resolved,
            max_R_over_R0=rmax,
            min_R_over_R0=rmin,
            cell_valid=np.zeros(cell_shape, bool),
        ).items():
            h5.create_dataset(name, data=a)
    lookup = ResponseTableBubbleModel(table_path, require_validated_cells=False)
    cell_valid = np.zeros(cell_shape, bool)
    comparisons = []
    for idx in np.ndindex(cell_shape):
        key = "mid_" + "_".join(map(str, idx))
        meta, actual = results[key]
        row = dict(index=idx, status="FAILED")
        if actual is not None:
            try:
                c = cases[key]
                pred = lookup.query(c.radius_m, c.frequency_Hz, c.pressure_Pa)["harmonics_m3"]
                e1 = float(abs(pred[0] - actual[0]) / max(abs(actual[0]), 1e-30))
                h = np.arange(1, 9)
                eh = float(np.linalg.norm(h * (pred - actual)) / max(np.linalg.norm(h * actual), 1e-30))
                cell_valid[idx] = max(e1, eh) < 0.05 and meta["metrics"]["valid"]
                row.update(
                    status="PASS" if cell_valid[idx] else "NEEDS_REFINEMENT",
                    fundamental_relative_error=e1,
                    derivative_weighted_harmonic_relative_error=eh,
                )
            except ValueError as e:
                row["error"] = str(e)
        comparisons.append(row)
    with h5py.File(table_path, "r+") as h5:
        h5["cell_valid"][:] = cell_valid
    convergence = []
    for key in ("node_1_1_0", "node_1_1_2", "node_1_1_4", "node_2_0_4"):
        for suffix in ("fine", "long"):
            a, b = results[key][1], results[key + "_" + suffix][1]
            error = None if a is None or b is None else float(np.linalg.norm(a - b) / np.linalg.norm(b))
            convergence.append(
                dict(
                    case=key,
                    comparison=suffix,
                    relative_complex_harmonic_L2=error,
                    passed=error is not None and error < 0.002,
                )
            )
    official_run = json.loads((root / "results/stage2/apecss_validation/STATE.json").read_text())["run"]
    ref = np.loadtxt(next((root / official_run / "lipidcoated_simple").glob("RPAR_*.txt")))
    raw = np.loadtxt(next((out / "official_lipid").glob("RPAR_*.txt")))
    from scipy.interpolate import PchipInterpolator

    compared = PchipInterpolator(raw[:, 1], raw[:, 3])(ref[:, 1])
    software_error = float(np.linalg.norm(compared - ref[:, 3]) / np.linalg.norm(ref[:, 3]))
    completed = sum(m[0]["status"] == "COMPLETE" for m in results.values())
    gate = completed == len(jobs) and all(c["passed"] for c in convergence) and software_error < 1e-4
    summary = dict(
        status="PASS_SOFTWARE_RADIAL" if gate else "FAILED_RADIAL_GATE",
        run=str(out.relative_to(root)),
        physics_status="partially validated",
        profile=args.profile,
        completed=completed,
        cases=len(jobs),
        workers=count,
        elapsed_s=elapsed,
        cases_per_second=completed / elapsed,
        child_cpu_s=cpu_s,
        average_cpu_cores=cpu_s / elapsed,
        cpu_worker_utilization=cpu_s / (elapsed * count),
        peak_aggregate_RSS_bytes=peak[0],
        RSS_note="sampled 0.1s; sums RSS including shared pages",
        official_lipid_waveform_relative_L2=software_error,
        convergence=convergence,
        periodic_valid_nodes=int(node_valid.sum()),
        total_nodes=int(node_valid.size),
        validated_cells=int(cell_valid.sum()),
        total_cells=int(cell_valid.size),
        cell_comparisons=comparisons,
        table=str(table_path.relative_to(root)),
        limitations=[
            "Software validation is not an experimental SonoVue calibration",
            "Coarse table cannot be used in cells failing midpoint checks",
            "Initial surface tension is a declared assumption",
        ],
    )
    (out / "SUMMARY.json").write_text(json.dumps(summary, indent=2) + "\n")
    (out.parent / "STATE.json").write_text(json.dumps(summary, indent=2) + "\n")
    fig, ax = plt.subplots(2, 2, figsize=(11, 8), constrained_layout=True)
    for pidx in (0, 2, 4):
        key = f"node_1_1_{pidx}"
        z = np.load(out / key / "cycles.npz")
        t = z["time_s"]
        select = t >= t[-1] - 2 / freqs[1]
        ax[0, 0].plot(
            (t[select] - t[select][0]) * 1e6,
            z["radius_m"][select] / radii[1],
            label=f"{pressures[pidx] / 1e3:g} kPa",
        )
    for ridx, r in enumerate(radii):
        ax[0, 1].plot(pressures / 1e3, rmax[ridx, 1], "o-", label=f"R0={r * 1e6:g} um")
    ax[1, 0].plot(ref[:, 1] * 1e6, ref[:, 3] / official.radius_m, label="Official example")
    ax[1, 0].plot(ref[:, 1] * 1e6, compared / official.radius_m, "--", label="Wrapper, refined step")
    for sigma in (0.0, 0.036, 0.072):
        z = np.load(out / f"fixed_table2_sigma{sigma:g}" / "cycles.npz")
        t = z["time_s"]
        select = t >= t[-1] - 2 / 2.5e6
        ax[1, 1].plot(
            (t[select] - t[select][0]) * 1e6, z["radius_m"][select] / 1.78e-6, label=f"sigma0={sigma:g} N/m"
        )
    for a, xlabel, ylabel in zip(
        ax.flat,
        ["Time (us)", "Pressure (kPa)", "Time (us)", "Time (us)"],
        ["R/R0", "Maximum R/R0", "R/R0", "R/R0"],
    ):
        a.set(xlabel=xlabel, ylabel=ylabel)
        a.legend(fontsize=8)
        a.grid(alpha=0.25)
    fig.suptitle("Figure A — APECSS nonlinear radial response; SonoVue mapping is a hypothesis")
    fig.savefig(out / "Figure_A_radial.png", dpi=180)
    plt.close(fig)
    print(
        json.dumps(
            {k: v for k, v in summary.items() if k not in ("cell_comparisons", "convergence")}, indent=2
        )
    )


if __name__ == "__main__":
    main()
