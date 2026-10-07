"""Final-model two-bubble and matched Secondary-OFF formation controls."""

import argparse
import json
from .common import BASE, write
from .experiment import run_case


def main(label, performance):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    step = json.loads((BASE / performance / "SELECTED_STEP.json").read_text())
    cfg = dict(
        count=32,
        seed=20261007,
        cloud_radius_m=20e-6,
        duration_s=0.05,
        dt_s=step["dt_s"],
        sample_interval_s=0.0005,
        physics={"contact_iterations": step["contact_iterations"], "secondary": False},
        events=[{"time_s": 0, "pressure_Pa": 1000}],
    )
    rows = {"formation_off": run_case(out / "formation_off", cfg)}
    for secondary in [False, True]:
        name = "pair_on" if secondary else "pair_off"
        pair = dict(
            cfg,
            count=2,
            positions_m=[[-4e-6, 0, 0], [4e-6, 0, 0]],
            physics=dict(cfg["physics"], secondary=secondary),
        )
        rows[name] = run_case(out / name, pair)
    write(out / "RESULT.json", dict(complete=True, rows=rows))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--performance", required=True)
    args = ap.parse_args()
    main(args.label, args.performance)
