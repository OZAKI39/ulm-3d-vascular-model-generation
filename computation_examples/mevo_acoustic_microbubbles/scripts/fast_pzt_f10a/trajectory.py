"""Single-location short comparisons, gated by preregistered static screening."""

import argparse
import time
import copy
import numpy as np
from .common import ROOT, BASE, js, write, identity, objects, clean
from .runner import run_case


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--phase", choices=["pilot", "half", "extend", "extended_half"], required=True)
    ap.add_argument("--key", default="O_1000")
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    proto = js(ROOT / "configs/fast_pzt/f10a/protocol.json")
    screen = js(BASE / "static_v1/RESULTS.json")
    eligible = [
        s
        for s in ["middle", "interior"]
        if all(dict(state=s, pressure_Pa=v) in screen["admitted"] for v in [1000, 3000])
    ]
    assert eligible, "No single location admits both1k and3k"
    state = eligible[0]
    env = copy.deepcopy(js(ROOT / "configs/fast_pzt/f8/protocol.json")["environment_Y"])
    x = np.array(screen["states"][state]["positions_m"])
    center = x.mean(0)
    env.update(
        positions_m=x,
        initial_centroid_m=center,
        target_m=center + np.array([0, 0, 1e-6]),
        task_type="short_positive_Z_increment",
        duration_s=0.002,
        split_group="F10A_translated_original_Y32",
        partition="exploratory_only",
    )
    env.pop("x_gate_m", None)
    env = clean(env)
    base = js(ROOT / "configs/fast_pzt/f8/protocol.json")
    base["sample_interval_s"] = 5e-5
    controls = {"OFF": env["events"]["off"]}
    for r in screen["rows"]:
        if r["state"] == state and r["pressure_Pa"] in [300, 1000, 3000]:
            controls[f"{r['mode']}_{r['pressure_Pa']}"] = r["event"]
    write(out / "IDENTITY.json", identity())
    write(
        out / "REGISTERED_PROTOCOL.json",
        dict(
            protocol=proto,
            phase=a.phase,
            location=state,
            positions_m=x,
            controls=controls,
            sample_interval_s=base["sample_interval_s"],
        ),
    )
    obj = objects()
    runs = {}

    def run(key, dt, duration):
        runs[key] = run_case(out / key, env, controls[key], base, obj, dt=dt, duration=duration, audit=True)
        write(
            out / "CAMPAIGN.json",
            dict(runs=runs, process_wall_s=time.perf_counter() - started, status="IN_PROGRESS"),
        )
        return runs[key]["numerical_quality"] == "PASS"

    if a.phase == "pilot":
        for key in ["OFF", "S_300", "O_300", "S_1000", "O_1000"]:
            run(key, 1e-5, 0.002)
        # Preserve fair paired failure evidence, but never increase pressure after failure.
        if all(s["numerical_quality"] == "PASS" for s in runs.values()):
            for key in ["S_3000", "O_3000"]:
                run(key, 1e-5, 0.002)
        else:
            write(
                out / "STOP_HIGHER_PRESSURE.json",
                dict(
                    reason="Pilot numerical quality failed; no3k trajectory admitted after contact/wall failure",
                    failed=[k for k, v in runs.items() if v["numerical_quality"] != "PASS"],
                ),
            )
    elif a.phase == "half":
        # A failing10us proposal is retained before checking this smaller step.
        run(a.key, 5e-6, 0.002)
    else:
        extension = js(ROOT / "configs/fast_pzt/f10a/extension_admission.json")
        assert a.key == extension["key"]
        pilot = js(BASE / "pilot_v3/CAMPAIGN.json")
        assert pilot["runs"][a.key]["numerical_quality"] == "PASS"
        for key in ["OFF", "S_" + a.key.split("_")[1], a.key]:
            run(key, 5e-6 if a.phase == "extended_half" else 1e-5, 0.01)
    write(
        out / "CAMPAIGN.json",
        dict(runs=runs, process_wall_s=time.perf_counter() - started, status="COMPLETE", state=state),
    )


if __name__ == "__main__":
    main()
