"""Early same-ID divergence: actual common frames and three local precision levels."""

import numpy as np
import matplotlib.pyplot as plt
from scripts.stage5a_r3.common import BASE
from scripts.stage5a_r3.paths import local
from scripts.stage5a_r3.figures import source, save, BLUE, RED, PURPLE, GRAY


def main():
    audit = source(BASE / "final/B2_EARLY_PRECISION_AUDIT.json")
    reports = {r["name"]: r for r in audit["comparisons"]}
    data = {}
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 8), layout="constrained")
    formal = source(BASE / "comparisons/B2_R3_V1_0.005_refined_REPORT.json")
    a = np.load(local(formal["samples"]))
    take = a["clock_s"] <= 3.00015 + 1e-14
    t = (a["clock_s"][take] - 3) * 1e6
    errors = a["all_active_raw_errors_m"][take] * 1e9
    for values, color, label in [
        (np.sqrt(np.mean(errors**2, axis=1)), BLUE, "All-active RMS"),
        (errors.max(1), RED, "Maximum same-ID error"),
    ]:
        axs[0, 0].plot(t, values, ".-", color=color, label=label)
    axs[0, 0].axhline(89, ls="--", color=GRAY, label="0.05 R resolution target")
    axs[0, 0].axvspan(30, 40, color=PURPLE, alpha=0.12, label="First target-crossing bracket")
    axs[0, 0].set(
        xlabel="Time after original checkpoint (µs)",
        ylabel="Difference (nm)",
        title="B2 formal nominal / refined | all 186 IDs",
    )
    axs[0, 0].legend(fontsize=8)
    data.update(formal_time_us=t, formal_error_nm=errors)
    x = np.arange(2)
    for kind, color in [("nominal", BLUE), ("refined", RED)]:
        values = [
            reports[f"{kind}_L{k}_L{k + 1}"]["groups"]["all_active"]["errors_m"]["raw"]["endpoint"]
            for k in range(2)
        ]
        for key, style in [("rms", "o-"), ("max", "s--")]:
            y = np.array([v[key] for v in values]) * 1e9
            axs[0, 1].plot(x, y, style, color=color, label=f"{kind.capitalize()} start: {key.upper()}")
            data[f"{kind}_adjacent_{key}_nm"] = y
    axs[0, 1].set(
        xticks=x,
        xticklabels=["L0 vs L1", "L1 vs L2"],
        ylabel="Endpoint difference (nm)",
        title="Same real start at +30 µs; local +20 µs",
    )
    axs[0, 1].legend(fontsize=8)
    b = np.load(local(reports["cross_initial_L2"]["samples"]))
    tt = (b["clock_s"] - 3) * 1e6
    for key, color, label in [
        ("raw", BLUE, "Raw"),
        ("translation_removed", RED, "Translation removed"),
        ("rigid_fit", PURPLE, "Rigid fit"),
    ]:
        e = b[f"all_active_{key}_errors_m"] * 1e9
        values = np.sqrt(np.mean(e**2, axis=1))
        axs[1, 0].plot(tt, values, "o-", color=color, label=label)
        data[f"cross_initial_L2_{key}_RMS_nm"] = values
    axs[1, 0].set(
        xlabel="Time after original checkpoint (µs)",
        ylabel="All-active RMS difference (nm)",
        title="Finest local runs | distinct actual initial states",
    )
    axs[1, 0].legend(fontsize=8)
    data["cross_initial_L2_clock_s"] = b["clock_s"]
    for index, kind in enumerate(("nominal", "refined")):
        c = np.load(local(audit["formal_event_correlations"][index]["arrays"]))
        time = (c["time_s"] - 3) * 1e6
        count = c["id_82_contact_band_count"]
        axs[1, 1].step(time, count, where="post", color=[BLUE, RED][index], lw=0.9, label=kind.capitalize())
        data[kind + "_contact_clock_s"] = c["time_s"]
        data[kind + "_id82_contact_band_count"] = count
    axs[1, 1].set(
        xlabel="Time after original checkpoint (µs)",
        ylabel="Neighbors within h ≤ 10 pm",
        title="ID 82 | actual pair proximity-band changes",
        xlim=(30, 50),
    )
    axs[1, 1].legend(fontsize=8)
    save(
        fig,
        "G5b_B2_early_contact_refinement",
        "Real stored states only, original IDs retained. Adjacent three-level differences decrease locally, but are not an independent global error bound. The two finest cross-lineage runs start from different actual states. Contact-band transitions are geometric correlation, not proof of multiplier activity or a unique cause. No S0 or fallback occurs in this formal +30..80 us interval; the two selected IDs have no WALL face switch.",
        data,
    )


if __name__ == "__main__":
    main()
