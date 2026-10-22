"""Finite direct checks, unmodified pinned APECSS library and ultrasound driver."""

import argparse
import subprocess
import time
import gzip
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from mevo_mb.bubbles.apecss import RadialCase, ApecssBubbleResponseModel, APECSS_COMMIT
from mevo_mb.bubbles.parameters import SonoVueParameterDistribution
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from .common import ROOT, BASE, js, write, sha, identity


def solve_one(args):
    exe, out, name, chi, visc, sigma, pressure, samples = args
    case = RadialCase(
        radius_m=1.78e-6,
        frequency_Hz=2.5e6,
        pressure_Pa=pressure,
        elasticity_N_m=chi,
        shell_viscosity_kg_s=visc,
        sigma_initial_N_m=sigma,
        samples_per_cycle=samples,
    )
    path = Path(out) / name
    try:
        meta, h, c = ApecssBubbleResponseModel(exe).solve(case, path)
        rb = case.radius_m / np.sqrt(1 + sigma / chi)
        rr = rb * np.sqrt(1 + case.sigma_water_N_m / chi)
        radius = c["radius_m"]
        table = ResponseSection(ROOT / "results/fast_pzt/f5/inputs/f4/response.h5", F2Physics()).query(
            np.array([pressure + 0j])
        )[0]
        primary_error = abs(table[0] - h[0]) / abs(h[0])

        def energy(z):
            return float(np.sum(np.arange(1, 9) ** 2 * abs(z) ** 2))

        row = dict(
            name=name,
            pressure_Pa=pressure,
            shell_family=name.split("_")[0],
            sigma_initial_N_m=sigma,
            elasticity_N_m=chi,
            shell_viscosity_kg_s=visc,
            Rb_m=rb,
            Rr_m=rr,
            Rmin_m=radius.min(),
            Rmax_m=radius.max(),
            contraction_fraction=radius.min() / case.radius_m - 1,
            expansion_fraction=radius.max() / case.radius_m - 1,
            branch_fractions=dict(
                buckled=np.mean(radius < rb),
                elastic=np.mean((radius >= rb) & (radius <= rr)),
                ruptured=np.mean(radius > rr),
            ),
            fundamental_abs_m3=abs(h[0]),
            fundamental_phase_rad=np.angle(h[0]),
            harmonics_real_m3=h.real,
            harmonics_imag_m3=h.imag,
            high_harmonic_dotV_energy_fraction=1 - abs(h[0]) ** 2 / energy(h),
            table_comparison_applicable=name.startswith("production_") and sigma == 0.036,
            table_fundamental_complex_relative_error=primary_error,
            table_dotV_energy_relative_error=abs(energy(table) / energy(h) - 1),
            metrics=meta["metrics"],
            elapsed_s=meta["elapsed_s"],
            samples_per_cycle=samples,
            status=meta["status"],
        )
        row["RADIAL_NUMERICAL_PASS"] = bool(
            meta["metrics"]["valid"] and meta["metrics"]["harmonic_energy_resolved"]
        )
        row["SHELL_BRANCH_CROSSED"] = bool(
            row["branch_fractions"]["buckled"] + row["branch_fractions"]["ruptured"] > 0
        )
        write(path / "ANALYSIS.json", row)
        for raw in path.glob("RPAR_*.txt"):
            with raw.open("rb") as source, gzip.open(str(raw) + ".gz", "wb") as dst:
                dst.write(source.read())
            raw.unlink()  # byte-exact gzip archive, no scientific samples removed
        return row
    except Exception as e:
        return dict(name=name, status="FAILED", error=repr(e), RADIAL_NUMERICAL_PASS=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    write(out / "IDENTITY.json", identity())
    proto = js(ROOT / "configs/fast_pzt/f10a/protocol.json")
    write(out / "REGISTERED_PROTOCOL.json", proto)
    deps = Path("/workspace/mevo_acoustic_microbubbles/f6/results/fast_pzt/f6/dependencies")
    source = deps / "apecss_v1_8"
    lib = deps / "build_99fcbfc0bd7c/libapecss.a"
    driver = source / "examples/ultrasound/src/ultrasound_apecss.c"
    assert (
        subprocess.check_output(["git", "-C", str(source), "rev-parse", "HEAD"], text=True).strip()
        == APECSS_COMMIT
    )
    exe = out / "apecss_ultrasound"
    cmd = [
        "gcc",
        "-O3",
        "-march=native",
        "-DNDEBUG",
        "-I",
        str(source / "include"),
        str(driver),
        str(lib),
        "-lm",
        "-o",
        str(exe),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    (out / "build.log").write_text(result.stdout + result.stderr)
    write(
        out / "DEPENDENCY.json",
        dict(
            commit=APECSS_COMMIT,
            command=cmd,
            library_sha256=sha(lib),
            driver_sha256=sha(driver),
            executable_sha256=sha(exe),
            unmodified_official_driver=True,
        ),
    )
    params = SonoVueParameterDistribution(ROOT / "configs/stage2/sonovue_parameters.json").parameter_set(
        1.78e-6
    )
    # Production parameters are read from the actual manifold, never silently relabelled Table 2.
    chi = params.shell_elasticity_N_m.value
    visc = params.shell_viscosity_kg_s.value
    jobs = []
    for family, ch, vi in [("production", chi, visc), ("TuTable2", 0.3, 3.2e-9)]:
        for pressure in proto["pressures_Pa"]:
            jobs.append((str(exe), str(out), f"{family}_p{pressure}_sigma036", ch, vi, 0.036, pressure, 512))
    for sigma in [0.0, 0.072]:
        for pressure in [3000, 10000]:
            jobs.append(
                (str(exe), str(out), f"TuTable2_p{pressure}_sigma{sigma}", 0.3, 3.2e-9, sigma, pressure, 512)
            )
    for pressure in [3000, 10000]:
        jobs.append((str(exe), str(out), f"production_p{pressure}_fine", chi, visc, 0.036, pressure, 1024))
    with ProcessPoolExecutor(max_workers=4) as pool:
        rows = list(pool.map(solve_one, jobs))
    response = ResponseSection(ROOT / "results/fast_pzt/f5/inputs/f4/response.h5", F2Physics())
    zero = response.query(np.array([0j]))
    assert not np.any(zero)
    rejected = {}
    for pressure in [0.01, 40001]:
        try:
            response.query(np.array([pressure + 0j]))
            rejected[str(pressure)] = "UNEXPECTED_ACCEPT"
        except Exception as e:
            rejected[str(pressure)] = repr(e)
    gate = {}
    for row in rows:
        if row.get("table_comparison_applicable") and row["samples_per_cycle"] == 512:
            gate[str(row["pressure_Pa"])] = bool(
                row["RADIAL_NUMERICAL_PASS"]
                and row["table_fundamental_complex_relative_error"]
                <= proto["radial_table_relative_tolerance"]
                and row["table_dotV_energy_relative_error"] <= proto["radial_table_relative_tolerance"]
            )
    convergence = {}
    for pressure in [3000, 10000]:
        coarse = next(r for r in rows if r["name"] == f"production_p{pressure}_sigma036")
        fine = next(r for r in rows if r["name"] == f"production_p{pressure}_fine")
        convergence[str(pressure)] = dict(
            fundamental_relative=abs(fine["fundamental_abs_m3"] / coarse["fundamental_abs_m3"] - 1),
            Rmax_relative=abs(fine["Rmax_m"] / coarse["Rmax_m"] - 1),
        )
    write(
        out / "RESULTS.json",
        dict(
            rows=rows,
            pressure_admission=gate,
            zero_exact=True,
            out_of_bounds=rejected,
            convergence=convergence,
            process_wall_s=time.perf_counter() - start,
            apecss_sum_solver_s=sum(r.get("elapsed_s", 0) for r in rows),
        ),
    )
    print("radial admission", gate, flush=True)


if __name__ == "__main__":
    main()
