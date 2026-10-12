"""Pinned APECSS low-pressure reference; no fitted production multiplier."""

import argparse
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict
import os
import subprocess
import time
import gzip
import h5py
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import cumulative_trapezoid
from mevo_mb.bubbles.apecss import RadialCase, ApecssBubbleResponseModel, APECSS_COMMIT
from mevo_mb.bubbles.parameters import SonoVueParameterDistribution
from validation.history_force_sensitivity.model import solve
from .common import ROOT, BASE, sha, write, js, identity


def build(out):
    deps = (BASE / "dependencies").resolve()
    deps.mkdir(exist_ok=True)
    src = deps / "apecss_v1_8"
    bld = deps / ("build_" + identity()["execution_commit"][:12])
    bld.mkdir(exist_ok=True)
    commands = []

    def run(cmd):
        with (out / "build.log").open("a") as log:
            subprocess.run(
                list(map(str, cmd)),
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                timeout=180,
                env=dict(os.environ, APECSS_DIR=str(src), USRLIB_DIR="/usr/lib/x86_64-linux-gnu"),
            )
        commands.append(list(map(str, cmd)))

    if not src.exists():
        run(
            ["git", "clone", "--branch", "v1.8", "--depth", "1", "https://github.com/polycfd/apecss.git", src]
        )
    assert (
        subprocess.check_output(["git", "-C", str(src), "rev-parse", "HEAD"], text=True).strip()
        == APECSS_COMMIT
    )
    run(
        ["cmake", "-S", src / "lib", "-B", bld, "-DCMAKE_BUILD_TYPE=Release", "-DCMAKE_C_FLAGS=-march=native"]
    )
    run(["cmake", "--build", bld, "--parallel", "4"])
    source = src / "examples/ultrasound/src/ultrasound_apecss.c"
    # Only the example's applied excitation callbacks change. Radial library untouched.
    pulse = source.read_text()
    callbacks = """
static APECSS_FLOAT f6_pressure(APECSS_FLOAT t, struct APECSS_Bubble *b) {
 double cycle=t*b->Excitation->f;
 if (cycle>=48.0 && cycle<72.0) return b->p0;
 return apecss_bubble_pressure_infinity_sinexcitation(t,b);
}
static APECSS_FLOAT f6_derivative(APECSS_FLOAT t, struct APECSS_Bubble *b) {
 double cycle=t*b->Excitation->f;
 if (cycle>=48.0 && cycle<72.0) return 0.0;
 return apecss_bubble_pressurederivative_infinity_sinexcitation(t,b);
}
"""
    pulse = pulse.replace("int main(", callbacks + "\nint main(", 1)
    pulse = pulse.replace(
        "  apecss_bubble_processoptions(Bubble);",
        "  apecss_bubble_processoptions(Bubble);\n  Bubble->get_pressure_infinity=f6_pressure;\n  Bubble->get_pressurederivative_infinity=f6_derivative;",
        1,
    )
    driver = out / "ultrasound_switch_driver.c"
    driver.write_text(pulse)
    run(
        [
            "gcc",
            "-O3",
            "-march=native",
            "-DNDEBUG",
            "-I",
            src / "include",
            driver,
            bld / "libapecss.a",
            "-lm",
            "-o",
            bld / "f6_switch",
        ]
    )
    write(
        out / "DEPENDENCY.json",
        dict(
            version="v1.8",
            commit=APECSS_COMMIT,
            source_url="https://github.com/polycfd/apecss/tree/v1.8",
            commands=commands,
            upstream_driver_sha256=sha(source),
            switch_driver_sha256=sha(driver),
            executable_sha256=sha(bld / "f6_switch"),
            modification_scope="Example excitation callbacks only; ON 0-48 cycles, OFF48-72, ON72-96; official radial ODE unchanged",
        ),
    )
    return bld / "f6_switch"


def one_case(task):
    out, exe, pressure = task
    out = Path(out)
    out.mkdir()
    par = SonoVueParameterDistribution(ROOT / "configs/stage2/sonovue_parameters.json").parameter_set(1.78e-6)
    case = RadialCase(
        radius_m=1.78e-6,
        frequency_Hz=2.5e6,
        pressure_Pa=pressure,
        elasticity_N_m=par.shell_elasticity_N_m.value,
        shell_viscosity_kg_s=par.shell_viscosity_kg_s.value,
        sigma_initial_N_m=0.036,
        cycles=96,
        samples_per_cycle=512,
        provenance="Existing Tu2009 Fig4 manifold; 0.036 N/m initial tension assumption; low-pressure F6 local reference",
    )
    meta, harmonics, cycles = ApecssBubbleResponseModel(exe).solve(case, out / "radial")
    raw_path = next((out / "radial").glob("RPAR_*.txt"))
    raw = np.loadtxt(raw_path)
    rad = PchipInterpolator(raw[:, 1], raw[:, 3])
    speed = PchipInterpolator(raw[:, 1], raw[:, 4])
    period = 1 / case.frequency_Hz
    convergence = []
    previous = None
    final = {}
    for spc in [64, 128, 256]:
        t = np.arange(96 * spc + 1) * period / spc
        r = rad(t)
        rd = speed(t)
        rdd = speed.derivative()(t)
        volume = 4 * np.pi * r**3 / 3
        ddv = 8 * np.pi * r * rd**2 + 4 * np.pi * r * r * rdd
        on = (t < 48 * period) | (t >= 72 * period)
        # Fixed-position standing wave dp/dx=-k P sin(wt); actual P,f scale.
        primary = volume * (2 * np.pi / (1500 * period)) * pressure * np.sin(2 * np.pi * t / period) * on
        pair = -1000 / (4 * np.pi * (10 * case.radius_m) ** 2) * volume * ddv
        # Secondary residual radial oscillations persist during OFF in reference.
        records = {}
        for label, force in [("single_primary", primary), ("pair_secondary", pair)]:
            steady_slice = slice(40 * spc, 48 * spc + 1)
            mean = np.trapezoid(force[steady_slice], t[steady_slice]) / (8 * period)
            vp = mean / (6 * np.pi * 0.001 * case.radius_m) * on
            vr = force / (6 * np.pi * 0.001 * r)
            models = dict(
                Aprod=dict(displacement_m=cumulative_trapezoid(vp, t, initial=0)),
                Aresolved=dict(displacement_m=cumulative_trapezoid(vr, t, initial=0)),
                B=solve(t, r, force, history=False),
                C=solve(t, r, force, history=True),
            )
            endpoints = {
                window: {model: float(v["displacement_m"][idx * spc]) for model, v in models.items()}
                for window, idx in [("before_off", 48), ("end_off", 72), ("after_reenable", 96)]
            }
            records[label] = endpoints
            if spc == 256:
                with h5py.File(out / (label + ".h5"), "x") as h:
                    h["time_s"] = t
                    h["radius_m"] = r
                    h["force_N"] = force
                    for k, m in models.items():
                        for key, val in m.items():
                            h.create_dataset(k + "/" + key, data=val, compression="gzip")
        if previous is not None:
            convergence.append(
                dict(
                    samples_per_cycle=spc,
                    relative_change={
                        label: {
                            w: {m: abs(v - previous[label][w][m]) / max(abs(v), 1e-15) for m, v in mm.items()}
                            for w, mm in ends.items()
                        }
                        for label, ends in records.items()
                    },
                )
            )
        previous = records
        final = records
    final_change = convergence[-1]["relative_change"]
    passed = all(v <= 0.02 for label in final_change.values() for w in label.values() for v in w.values())
    bias = {
        label: {
            w: dict(
                signed_C_over_Aprod=mm["C"] / mm["Aprod"] if abs(mm["Aprod"]) > 1e-15 else None,
                relative_Aprod_minus_C=(mm["Aprod"] - mm["C"]) / max(abs(mm["C"]), 1e-15),
                same_direction=bool(mm["C"] * mm["Aprod"] > 0),
            )
            for w, mm in rows.items()
        }
        for label, rows in final.items()
    }
    report = dict(
        status="PASS" if passed else "FAIL_CONVERGENCE",
        case=asdict(case),
        endpoints_m=final,
        convergence=convergence,
        bias=bias,
        calibration="UNCALIBRATED",
        radial_sha256=sha(raw_path),
        scope="Fixed-position, imposed identical radii in far-field pair at10R; no mutual radial feedback/walls or long-time .3-1 s history. Two-bubble per-bubble signed pair displacement, counterpart opposite. Switching at48/72 cycles resolves local transients, not complete navigation-time history. No correction fitted.",
    )
    write(out / "SUMMARY.json", report)
    with raw_path.open("rb") as src, gzip.open(str(raw_path) + ".gz", "wb") as dst:
        dst.write(src.read())
    raw_path.unlink()  # only newly generated duplicate uncompressed reference output
    return report


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--protocol", default="preflight_v2")
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    start = time.perf_counter()
    write(out / "IDENTITY.json", identity())
    exe = build(out)
    static = js(BASE / a.protocol / "STATIC_PREFLIGHT.json")
    actual = max(max(v["pressure_Pa"]) for d in static for v in d["controls"].values())
    pressures = [305.0, float(actual)]
    write(
        out / "REFERENCE_PROTOCOL.json",
        dict(
            pressures_Pa=pressures,
            frequency_Hz=2.5e6,
            cycles=96,
            off_cycle=48,
            reenable_cycle=72,
            samples_per_cycle=[64, 128, 256],
            convergence_limit=0.02,
            correction="none",
            max_workers=2,
        ),
    )
    with ProcessPoolExecutor(max_workers=2) as pool:
        rows = list(
            pool.map(one_case, [(out / ("pressure_" + str(i)), exe, p) for i, p in enumerate(pressures)])
        )
    write(
        out / "REFERENCE.json",
        dict(
            status="PASS" if all(r["status"] == "PASS" for r in rows) else "FAIL_CONVERGENCE",
            rows=rows,
            wall_s=time.perf_counter() - start,
            production_correction="none",
            calibration="UNCALIBRATED",
        ),
    )
    print("LOCAL REFERENCE COMPLETE", flush=True)


if __name__ == "__main__":
    main()
