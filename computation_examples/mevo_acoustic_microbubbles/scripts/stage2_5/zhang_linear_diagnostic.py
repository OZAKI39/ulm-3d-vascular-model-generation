"""Diagnostic linearization, not a replacement of Zhang's nonlinear equations."""

from dataclasses import replace
from pathlib import Path
import sys
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from references.zhang2011_minimal import Parameters, residual
from mevo_mb.validation.zhang2011 import Zhang2011EncapsulatedPairReference
from zhang_forensic import production_parameters
from common import new_run, finish, watermark


def linearized(p, relative_step=1e-5):
    p0 = replace(p, Pa=0)
    y = np.array([p.R10, 0, 0, 0, p.R20, 0, p.D0, 0])
    zero = np.zeros(4)
    M = np.column_stack([residual(0, y, np.eye(4)[j], p0) - residual(0, y, zero, p0) for j in range(4)])
    columns = []
    for indices in ([0, 4, 2, 6], [1, 5, 3, 7]):
        mat = []
        for index in indices:
            step = relative_step * (
                min(p.R10, p.R20) if index % 2 == 0 else 2 * np.pi * p.f * min(p.R10, p.R20)
            )
            a, b = y.copy(), y.copy()
            a[index] += step
            b[index] -= step
            mat.append((residual(0, a, zero, p0) - residual(0, b, zero, p0)) / (2 * step))
        columns.append(np.column_stack(mat))
    K, C = columns
    w = 2 * np.pi * p.f
    forcing = np.array([-1j * p.Pa / p.rho, -1j * p.Pa / p.rho, 0, 0])
    displacement = np.linalg.solve(K + 1j * w * C - w * w * M, forcing)
    V = 4 * np.pi * np.array([p.R10, p.R20]) ** 2 * displacement[:2]
    return float(p.rho / (8 * np.pi) * w * w * np.real(V[0] * V[1].conjugate()))


def main():
    start = time.perf_counter()
    out = new_run(
        "zhang_reproduction",
        dict(
            purpose="test whether reported Figure3 resembles small-signal response; never refit Pa or chi",
            distances_um=[30, 40, 50, 60, 70, 80, 100, 120],
        ),
        ["Zhang2011 Eqs3-11, Fig3"],
    )
    rows = []
    for d in [30, 40, 50, 60, 70, 80, 100, 120]:
        p = Parameters(D0=d * 1e-6)
        linear = linearized(p)
        small = replace(p, Pa=100.0)
        z = Zhang2011EncapsulatedPairReference(production_parameters(small)).solve(
            cycles=128, samples_per_cycle=128, rtol=2e-10, max_steps_per_cycle=96
        )
        t = z["time_s"]
        use = t * small.f >= 120 - 1e-7
        dv = 4 * np.pi * z["radius_m"] ** 2 * z["radial_speed_m_s"]
        weak = (
            small.rho
            / (4 * np.pi)
            * np.trapezoid(np.prod(dv[use], axis=1), t[use])
            / (t[use][-1] - t[use][0])
            * (p.Pa / small.Pa) ** 2
        )
        rows.append(
            dict(
                distance_um=d,
                linearized_coefficient_N_m2=linear,
                weak_100Pa_scaled_to_100kPa_N_m2=float(weak),
                derivative_refinement_relative_error=abs(linearized(p, 1e-6) - linear)
                / max(abs(linear), 1e-30),
            )
        )
    import json

    digit = json.loads(Path("validation/zhang2011_forensic/FIG3_DIGITIZATION.json").read_text())
    fig, ax = plt.subplots(figsize=(7, 5), constrained_layout=True)
    ax.errorbar(
        [v["distance_m"] * 1e6 for v in digit["points"]],
        [v["coefficient_N_m2"] / 1e-18 for v in digit["points"]],
        yerr=[v["uncertainty_N_m2"] / 1e-18 for v in digit["points"]],
        fmt="o",
        label="Approximate paper digitization",
    )
    ax.plot(
        [r["distance_um"] for r in rows],
        [r["linearized_coefficient_N_m2"] / 1e-18 for r in rows],
        label="Independent residual linearization",
    )
    ax.plot(
        [r["distance_um"] for r in rows],
        [r["weak_100Pa_scaled_to_100kPa_N_m2"] / 1e-18 for r in rows],
        "--",
        label="100Pa nonlinear solve, scaled P² (diagnostic)",
    )
    old = json.loads(Path("results/stage2/pair_zhang2011/STATE.json").read_text())
    ds = [30, 40, 50, 70, 100, 120]
    ax.plot(
        ds,
        [old["cases"][f"distance_{d}"]["metrics"]["secondary_coefficient_N_m2"] / 1e-18 for d in ds],
        label="Actual100kPa nonlinear equations (Stage2)",
    )
    ax.axhline(0, color="grey", lw=0.7)
    ax.legend(fontsize=8)
    ax.set(
        xlabel="Initial separation (um)",
        ylabel="fB / 1e-18 N m²",
        title="Figure3: distinguish nonlinear equations from small-signal limit",
    )
    watermark(fig, "UNRESOLVED — LINEARIZATION IS A DIAGNOSTIC, NOT A PAPER CORRECTION")
    fig.savefig(out / "LINEAR_VS_NONLINEAR.png", dpi=180)
    finish(
        out,
        dict(
            status="DIAGNOSTIC_COMPLETE_NOT_VALIDATION_PASS", rows=rows, runtime_s=time.perf_counter() - start
        ),
    )


if __name__ == "__main__":
    main()
