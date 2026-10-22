#!/usr/bin/env python3
"""Bounded sampler-only validation. No simulation engine or transport physics."""
from pathlib import Path
import hashlib
import json
import platform
import subprocess
import sys
import tempfile
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sonovue_sampler import SonoVueDistribution, SAMPLER_VERSION, sha256_file, write_population

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "validation"
PLAN = json.loads((ROOT / "contracts/VALIDATION_PLAN_V0.json").read_text())


def save(path, data):
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def digest_array(array):
    return hashlib.sha256(np.asarray(array, dtype="<f8").tobytes()).hexdigest()


def run():
    if (OUT / "SONOVUE_SAMPLER_VALIDATION.json").exists():
        raise FileExistsError("Frozen validation exists; use a new isolated copy for another audit")
    d = SonoVueDistribution()
    assert d.histogram_sha256 == PLAN["histogram_sha256"]
    checks = {"INPUT_HISTOGRAM_CHECK": True}
    widths = d.high - d.low
    knots = np.unique(np.r_[d.low, d.high])
    knot_cdf = d.cdf(knots)
    static = np.r_[d.cdf_low, d.cdf_high, knot_cdf]
    checks["CDF_MONOTONIC"] = bool(np.all(np.diff(knot_cdf) >= 0))
    checks["CDF_ENDPOINTS"] = bool(d.cdf(d.support_min_um) == 0 and d.cdf(d.support_max_um) == 1)
    checks["CDF_STATIC_CHECK"] = bool(np.isfinite(static).all() and np.all((static >= 0) & (static <= 1))
                                      and np.all(widths > 0) and np.all(d.pdf_per_um >= 0)
                                      and abs(d.normalized_probability_sum - 1) <= 32*np.finfo(float).eps)
    zero = d.probability == 0
    checks["ZERO_BINS_PRESERVED"] = bool(np.all(d.cdf_high[zero] == d.cdf_low[zero])
                                         and np.all(d.pdf_per_um[zero] == 0))
    d.export_cdf(OUT / "SONOVUE_EMPIRICAL_CDF.csv")

    inv_plan = PLAN["inverse_consistency"]
    u = np.random.default_rng(inv_plan["seed"]).random(inv_plan["N"])
    diameter = d.inverse_cdf(u)
    inverse_error = float(np.max(np.abs(d.cdf(diameter) - u)))
    # Independent direct integral over all original bins; does not use the CDF search code.
    direct_u = np.sum(np.clip((diameter[:, None] - d.low) / widths, 0, 1) * d.probability, axis=1)
    independent_error = float(np.max(np.abs(direct_u - u)))
    boundary_u = np.unique(np.r_[0., d.cdf_low, d.cdf_high, 1.])
    boundary_error = float(np.max(np.abs(d.cdf(d.inverse_cdf(boundary_u)) - boundary_u)))
    checks["INVERSE_CDF_CHECK"] = max(inverse_error, independent_error, boundary_error) <= inv_plan["max_absolute_error"]

    sample_plan = PLAN["sampling_smoke"]
    sample_file = OUT / "VALIDATION_SONOVUE_POPULATION_100000.csv"
    population_meta = write_population(d, sample_plan["N"], sample_plan["seed"], sample_file, sample_plan["role"])
    sample = np.genfromtxt(sample_file, delimiter=",", names=True)
    x = sample["diameter_um"]
    n = len(x)
    assert n == sample_plan["N"]
    assert digest_array(x) == population_meta["diameter_sequence_sha256_float64_le"]
    checks["SAMPLE_FINITE"] = bool(np.isfinite(x).all())
    checks["SAMPLE_SUPPORT"] = bool(np.all((x >= d.support_min_um) & (x <= d.support_max_um)))
    # Original input has contiguous edges. No zero bins are removed or merged.
    assert np.array_equal(d.low[1:], d.high[:-1])
    edges = np.r_[d.low, d.high[-1]]
    counts, _ = np.histogram(x, bins=edges)
    frequencies = counts / n
    bin_errors = frequencies - d.probability
    max_bin_error = float(np.max(np.abs(bin_errors)))
    rmse = float(np.sqrt(np.mean(bin_errors**2)))
    checks["HISTOGRAM_REPRODUCTION_SMOKE"] = max_bin_error <= PLAN["max_absolute_bin_probability_error"]
    checks["ZERO_BIN_SAMPLES"] = bool(np.all(counts[zero] == 0))
    source_idx = np.searchsorted(d.low, x, side="right") - 1
    checks["SAMPLE_IN_POSITIVE_SOURCE_INTERVAL"] = bool(np.all(d.probability[source_idx] > 0)
                                                        and np.all(x < d.high[source_idx])
                                                        and np.array_equal(sample["source_bin_low_um"], d.low[source_idx])
                                                        and np.array_equal(sample["source_bin_high_um"], d.high[source_idx]))
    sorted_x = np.sort(x)
    target_at_samples = d.cdf(sorted_x)
    d_plus = float(np.max(np.arange(1, n+1)/n - target_at_samples))
    d_minus = float(np.max(target_at_samples - np.arange(n)/n))
    ks_deviation = max(d_plus, d_minus)
    checks["CDF_REPRODUCTION_SMOKE"] = ks_deviation <= PLAN["max_CDF_deviation"]
    target_mean = d.target_mean_um()
    sample_mean = float(np.mean(x))
    mean_error = abs(sample_mean - target_mean)
    checks["MEAN_SANITY"] = mean_error <= PLAN["max_absolute_mean_difference_um"]
    probabilities = np.array([.1, .5, .9])
    target_quantiles = d.inverse_cdf(probabilities)
    sample_quantiles = np.quantile(x, probabilities, method="linear")
    rows = np.column_stack((np.arange(len(d.low)), d.low, d.high, d.probability, counts, frequencies, bin_errors))
    np.savetxt(OUT / "BIN_REPRODUCTION.csv", rows, delimiter=",", fmt=["%d", "%.17g", "%.17g", "%.17g", "%d", "%.17g", "%.17g"],
               header="bin_index,bin_low_um,bin_high_um,target_probability,sample_count,sample_probability,error", comments="")

    repro = PLAN["reproducibility"]
    repro_results = []
    with tempfile.TemporaryDirectory(prefix="repro_check_", dir=OUT) as tmp:
        paths = [Path(tmp) / "replicate_a.csv", Path(tmp) / "replicate_b.csv"]
        for path in paths:  # Two independent process invocations, not a per-sample loop.
            cmd = [sys.executable, "-B", str(ROOT / "src/sonovue_sampler.py"), "--histogram", str(d.histogram_path),
                   "--n", str(repro["N"]), "--seed", str(repro["seed"]), "--output", str(path), "--role", "SAMPLER_VALIDATION_ONLY"]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            arr = np.genfromtxt(path, delimiter=",", names=True)["diameter_um"]
            meta = json.loads(path.with_suffix(".metadata.json").read_text())
            assert meta["histogram_sha256"] == d.histogram_sha256
            repro_results.append({"returncode": result.returncode, "N": len(arr), "seed": meta["seed"],
                                  "csv_sha256": sha256_file(path), "diameters_float64_le_sha256": digest_array(arr),
                                  "numpy_version": meta["numpy_version"], "sampler_version": meta["sampler_version"]})
        a = np.genfromtxt(paths[0], delimiter=",", names=True)["diameter_um"]
        b = np.genfromtxt(paths[1], delimiter=",", names=True)["diameter_um"]
        exact_elements = bool(np.array_equal(a, b))
        exact_array_bytes = a.astype("<f8").tobytes() == b.astype("<f8").tobytes()
        exact_csv_bytes = paths[0].read_bytes() == paths[1].read_bytes()
    checks["REPRODUCIBILITY_CHECK"] = exact_elements and exact_array_bytes and exact_csv_bytes
    save(OUT / "REPRODUCIBILITY_CHECK.json", {"status": "PASS" if checks["REPRODUCIBILITY_CHECK"] else "FAIL",
         "independent_processes": repro_results, "exact_elements": exact_elements, "exact_array_bytes": exact_array_bytes,
         "exact_csv_bytes": exact_csv_bytes, "temporary_files_removed_after_comparison": True,
         "scope": "Two independent Python processes with the same NumPy 1.26.4 environment; not a promise for untested NumPy versions or architectures"})

    demo = PLAN["demo"]
    write_population(d, demo["N"], demo["seed"], OUT / "DEMO_SONOVUE_POPULATION_1000.csv", demo["role"])
    perf = PLAN["performance"]
    start = time.perf_counter()
    million = d.sample_diameters(perf["N"], perf["seed"])
    elapsed = time.perf_counter() - start
    perf_digest = digest_array(million)
    assert len(million) == perf["N"] and np.isfinite(million).all()
    del million

    # Plain QA plots, using the original bin widths and all 100000 empirical CDF steps.
    with plt.rc_context({"font.size": 10, "figure.dpi": 140}):
        fig, ax = plt.subplots(figsize=(9, 4.8))
        ax.bar(d.low, d.probability*100, width=widths, align="edge", color="#4477aa", alpha=.45,
               edgecolor="white", linewidth=.3, label="Frozen histogram (number weighted)")
        ax.stairs(frequencies*100, edges, color="#cc6633", linewidth=1.4, label="100,000 inverse-CDF samples")
        ax.set(xlabel="Diameter (um)", ylabel="Probability per original bin (%)", title="SonoVue sampler V0: histogram reproduction")
        ax.set_xlim(d.support_min_um, d.support_max_um)
        ax.legend(frameon=False);ax.grid(axis="y", alpha=.2)
        fig.text(.12, .01, "Sampler validation only; seed=20260915; zero bins and source-defined support retained.", fontsize=8)
        fig.tight_layout(rect=(0, .03, 1, 1));fig.savefig(OUT / "SONOVUE_SAMPLER_QA.png");plt.close(fig)
        fig, ax = plt.subplots(figsize=(9, 4.8))
        ax.step(np.r_[d.support_min_um, sorted_x, d.support_max_um], np.r_[0, np.arange(1, n+1)/n, 1],
                where="post", color="#cc6633", linewidth=1.1, label="Sample ECDF (100,000 samples)")
        ax.plot(knots, knot_cdf, color="#225588", linewidth=1.5, linestyle="--", label="Target continuous piecewise-linear CDF")
        ax.set(xlabel="Diameter (um)", ylabel="Cumulative probability", title="SonoVue sampler V0: CDF reproduction",
               xlim=(d.support_min_um, d.support_max_um), ylim=(0, 1.02))
        ax.legend(frameon=False);ax.grid(alpha=.2)
        fig.text(.12, .01, f"KS-like deviation={ks_deviation:.6g}; descriptive smoke gate <=0.02, not a publication-level test.", fontsize=8)
        fig.tight_layout(rect=(0, .03, 1, 1));fig.savefig(OUT / "SONOVUE_CDF_QA.png");plt.close(fig)

    original = Path('/home/lzy/projects/FROZEN_SONOVUE_HISTOGRAM.csv')
    checks["ORIGINAL_INPUT_UNCHANGED"] = sha256_file(original) == d.histogram_sha256
    status = "PASS" if all(checks.values()) else "FAIL"
    record = {
        "status": status, "stage": "SONOVUE_CONTINUOUS_SAMPLER_V0", "checks": {k: "PASS" if v else "FAIL" for k,v in checks.items()},
        "histogram_sha256": d.histogram_sha256, "sampler_version": SAMPLER_VERSION,
        "sampler_source_sha256": sha256_file(ROOT / "src/sonovue_sampler.py"),
        "validation_source_sha256": sha256_file(__file__), "validation_plan_sha256": sha256_file(ROOT / "contracts/VALIDATION_PLAN_V0.json"),
        "bin_count": len(d.low), "zero_probability_bins": int(zero.sum()), "zero_probability_bin_sample_count": int(counts[zero].sum()),
        "bin_width_um": d.bin_width_um, "bin_width_min_um": float(widths.min()), "bin_width_max_um": float(widths.max()),
        "support_min_um": d.support_min_um, "support_max_um": d.support_max_um,
        "INPUT_PROBABILITY_SUM": d.input_probability_sum, "NORMALIZED_PROBABILITY_SUM": d.normalized_probability_sum,
        "cdf_mass_vs_normalized_probability_max_abs_roundoff": d.cdf_mass_roundoff_max,
        "validation_N": n, "validation_seed": sample_plan["seed"],
        "inverse_N": inv_plan["N"], "inverse_seed": inv_plan["seed"],
        "INVERSE_CDF_MAX_ERROR": inverse_error, "independent_integral_inverse_max_error": independent_error,
        "inverse_boundary_max_error": boundary_error, "MAX_BIN_PROBABILITY_ERROR": max_bin_error,
        "RMSE_BIN_PROBABILITY": rmse, "MAX_CDF_DEVIATION": ks_deviation, "KS_like_D_plus": d_plus, "KS_like_D_minus": d_minus,
        "KS_like_definition": "max_i(max(i/N-F(x_i), F(x_i)-(i-1)/N)) over all sorted samples; no p-value or fitted-model claim",
        "TARGET_MEAN_DIAMETER_UM": target_mean, "VALIDATION_SAMPLE_MEAN_UM": sample_mean, "mean_absolute_difference_um": mean_error,
        "REFERENCE_MEAN_DIAMETER_UM": PLAN["reference_mean_um"], "REFERENCE_MEAN_ROLE": "INFORMATIONAL_ONLY",
        "target_minus_reference_mean_um": target_mean - PLAN["reference_mean_um"],
        "target_quantiles_um": dict(zip(["d10", "d50", "d90"], target_quantiles.tolist())),
        "validation_quantiles_um": dict(zip(["d10", "d50", "d90"], sample_quantiles.tolist())),
        "sample_minus_target_quantiles_um": dict(zip(["d10", "d50", "d90"], (sample_quantiles-target_quantiles).tolist())),
        "sample_quantile_method": "numpy.quantile(method='linear'); diagnostic only",
        "performance": {"N": perf["N"], "seed": perf["seed"], "wall_seconds": elapsed,
                        "samples_per_second": perf["N"]/elapsed, "sequence_sha256_float64_le": perf_digest,
                        "scope": "One in-memory draw including default_rng construction and vectorized inverse; excludes CSV parsing, file I/O and QA plots",
                        "role": perf["role"], "performance_tuning": "OUT_OF_SCOPE", "million_sample_file_saved": False},
        "environment": {"python": sys.version, "executable": sys.executable, "numpy": np.__version__, "matplotlib": matplotlib.__version__,
                        "platform": platform.platform(), "dtype": "float64", "bit_generator": type(np.random.default_rng(0).bit_generator).__name__},
        "NO_PARAMETRIC_FIT": "YES", "NO_HISTOGRAM_REWEIGHTING": "YES", "NO_SMOOTHING": "YES",
        "relative_probabilities_changed": False, "tail_truncation": False, "tail_extension": False,
        "formal_simulation_population_generated": False, "automatic_next_stage": False,
    }
    save(OUT / "SONOVUE_SAMPLER_VALIDATION.json", record)
    print(json.dumps(record, indent=2, ensure_ascii=False))
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    run()
