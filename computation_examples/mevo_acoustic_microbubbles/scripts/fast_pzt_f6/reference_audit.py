"""Compare the local reference against the actual unchanged production table.

The reference's Aprod averages its newly resolved radial cycles. At low pressure
this is close to, but is not identical to, interpolation of the stored table.
Keep both baselines explicit without editing the original reference files.
"""

import argparse
import numpy as np
import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mevo_mb.fast_sim.models import ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics
from .common import BASE, INPUT, js, write, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="reference_v2")
    a = ap.parse_args()
    folder = BASE / a.label
    out = folder / "table_audit"
    out.mkdir(exist_ok=False)
    ref = js(folder / "REFERENCE.json")
    p = F2Physics()
    response = ResponseSection(INPUT / "f4/response.h5", p)
    omega = 2 * np.pi * p.frequency_Hz
    zeta = 6 * np.pi * p.viscosity_Pa_s * p.radius_m
    rows = []
    fig, axs = plt.subplots(2, 2, figsize=(9, 6), layout="constrained")
    for index, row in enumerate(ref["rows"]):
        pressure = row["case"]["pressure_Pa"]
        v = response.query(np.array([1j * pressure]))[0]
        gradient = 1j * pressure * omega / 1500
        forces = dict(
            single_primary=-0.5 * np.real(v[0] * np.conj(gradient)),
            pair_secondary=p.density_kg_m3
            * omega**2
            / (8 * np.pi * (10 * p.radius_m) ** 2)
            * np.sum(np.arange(1, 9) ** 2 * abs(v) ** 2),
        )
        for col, (label, force) in enumerate(forces.items()):
            endpoints = {}
            for window, on_cycles in [("before_off", 48), ("end_off", 48), ("after_reenable", 72)]:
                displacement = force / zeta * on_cycles / p.frequency_Hz
                old = row["endpoints_m"][label][window]
                endpoints[window] = dict(
                    actual_table_fast_displacement_m=float(displacement),
                    resolved_cycle_mean_displacement_m=old["Aprod"],
                    instantaneous_stokes_m=old["Aresolved"],
                    added_mass_m=old["B"],
                    added_mass_and_history_m=old["C"],
                    cycle_mean_vs_table_relative=float((old["Aprod"] - displacement) / displacement),
                    signed_history_over_table=float(old["C"] / displacement),
                )
            rows.append(
                dict(pressure_Pa=pressure, force_model=label, force_N=float(force), endpoints=endpoints)
            )
            with h5py.File(folder / f"pressure_{index}" / (label + ".h5")) as h:
                t = h["time_s"][:]
                cycles = t * p.frequency_Hz
                on_time = np.minimum(cycles, 48) + np.maximum(cycles - 72, 0)
                ax = axs[index, col]
                ax.plot(
                    cycles,
                    force / zeta * on_time / p.frequency_Hz * 1e9,
                    label="Stored-table fast model",
                    color="#34265c",
                )
                for model, name, color in [
                    ("Aresolved", "Instantaneous Stokes", "#be849f"),
                    ("B", "Added mass", "#7387a5"),
                    ("C", "Added mass + history", "#333333"),
                ]:
                    ax.plot(cycles, h[model + "/displacement_m"][:] * 1e9, label=name, color=color, lw=0.8)
                ax.axvspan(48, 72, color="#999999", alpha=0.15)
                ax.set(
                    title=f"{label.replace('_', ' ').title()}, {pressure:.0f} Pa",
                    xlabel="Acoustic cycles",
                    ylabel="Signed displacement (nm)",
                )
                ax.spines[["top", "right"]].set_visible(False)
    axs[0, 0].legend(fontsize=7)
    fig.savefig(out / "Fig_F6_local_translation.png", dpi=220, facecolor="white")
    fig.savefig(out / "Fig_F6_local_translation.pdf", facecolor="white")
    plt.close(fig)
    write(
        out / "TABLE_BASELINE_AUDIT.json",
        dict(
            identity=identity(),
            rows=rows,
            calibration="UNCALIBRATED",
            scope="Exact original table interpolation and original eight-harmonic far-field force at 10R. Signed local reference only, no full-navigation history correction. Shaded plot interval is sound OFF. No new radial or trajectory integrations.",
        ),
    )


if __name__ == "__main__":
    main()
