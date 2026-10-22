"""Additional comparison and explicit gated Figure D; no clustering simulation."""

from dataclasses import replace
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mevo_mb.validation.zhang2011 import (
    Zhang2011Parameters,
    Zhang2011EncapsulatedPairReference,
    secondary_coefficient,
)


def main():
    root = Path.cwd()
    s = json.loads((root / "results/stage2/pair_zhang2011/STATE.json").read_text())
    out = root / s["run"]
    p = replace(Zhang2011Parameters(), distance_m=1.0)
    z = Zhang2011EncapsulatedPairReference(p).solve(cycles=64)
    coefficient = secondary_coefficient(z, p.density_kg_m3)
    np.savez_compressed(
        out / "uncoupled_radial_limit.npz", **{k: v for k, v in z.items() if isinstance(v, np.ndarray)}
    )
    ds = np.array([30.0, 40.0, 50.0, 70.0, 100.0, 120.0])
    coeff = np.array([s["cases"][f"distance_{d:g}"]["metrics"]["secondary_coefficient_N_m2"] for d in ds])
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    axs[0].plot(ds, coeff / 1e-18, "o-", label="Coupled Zhang equations")
    axs[0].axhline(coefficient / 1e-18, ls="--", label="Independent radial limit")
    axs[0].axhline(0, color="grey", lw=0.6)
    axs[0].set(
        xlabel="Initial distance (um)", ylabel="f_B (10^-18 N m²)", title="Fig3 sign reversal NOT reproduced"
    )
    axs[0].legend()
    axs[1].plot(ds, coeff / (ds * 1e-6) ** 2 * 1e9, "o-", label="Coupled-response far-field estimate")
    axs[1].plot(ds, coefficient / (ds * 1e-6) ** 2 * 1e9, "--", label="Independent-response far-field")
    axs[1].set(
        xlabel="Initial distance (um)",
        ylabel="Secondary force estimate (nN)",
        title="Coefficient comparison, not near-field validation",
    )
    axs[1].legend(fontsize=8)
    fig.savefig(out / "Figure_C_reference_vs_farfield.png", dpi=220)
    plt.close(fig)
    (out / "FARFIELD_COMPARISON.json").write_text(
        json.dumps(
            dict(
                independent_coefficient_N_m2=coefficient,
                method="Zhang radial law at 1m numerical decoupling limit; not acoustic physical propagation at 1m",
                pair_case="3um and6um,100kPa,2.25MHz",
                status="DIAGNOSTIC_ONLY",
            ),
            indent=2,
        )
        + "\n"
    )
    fig, axs = plt.subplots(1, 3, figsize=(10, 3), constrained_layout=True)
    for ax, pressure in zip(axs, (50, 75, 100)):
        ax.text(0.5, 0.55, "NOT RUN", ha="center", va="center", fontsize=17, color="darkred")
        ax.text(
            0.5,
            0.25,
            "Stage 2G validation gate\nNo clustering data generated",
            ha="center",
            va="center",
            fontsize=10,
        )
        ax.set_title(f"1 MHz / {pressure} kPa / 20 ms")
        ax.set_xticks([])
        ax.set_yticks([])
    fig.savefig(root / "results/stage2/lazarus2017/Figure_D_gate_status.png", dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    main()
