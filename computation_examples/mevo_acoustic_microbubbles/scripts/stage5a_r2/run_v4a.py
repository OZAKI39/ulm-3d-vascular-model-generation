"""Versioned formal entry: V4 contract with unchanged V3/V6 integration loop."""

import argparse
from scripts.stage5a_r2.validation_v4a import verify
from scripts.stage5a_r2.run import run as frozen_run


def parser():
    p = argparse.ArgumentParser()
    p.add_argument("--case", choices=["B2", "B3"], required=True)
    p.add_argument("--duration", type=float, required=True)
    p.add_argument("--max-dt", type=float, required=True)
    p.add_argument("--tol-factor", type=float, required=True)
    p.add_argument("--pair", required=True)
    p.add_argument("--nonsmooth-recovery", action="store_true")
    p.add_argument("--nonlinear-strategy", required=True)
    p.add_argument("--label", required=True)
    p.add_argument("--budget", type=float, required=True)
    return p


def execute(a):
    protocol = verify()
    allowed = [(5e-6, 0.25), (2.5e-6, 0.0625)]
    valid = (
        a.pair == "cupy"
        and a.nonsmooth_recovery
        and a.nonlinear_strategy == "quasi_newton"
        and (a.max_dt, a.tol_factor) in allowed
        and 0 < a.budget <= protocol["resources"]["per_run_s"]
        and a.duration in [0.001, 0.005]
        and a.label.startswith("formal_v4a_")
    )
    if not valid:
        raise ValueError("CONFIG_ERROR: V4A formal parameter/resource contract")
    # The legacy "formal*" label routes through the immutable V3 whitelist.
    # This distinct entry validates the complete stricter V4A contract above,
    # then invokes exactly the same integration loop with an explicit V4A ID.
    return frozen_run(
        a.case,
        a.duration,
        method="S1",
        max_dt=a.max_dt,
        tol_factor=a.tol_factor,
        pair=a.pair,
        budget=a.budget,
        label="validated_v4a_" + a.label.removeprefix("formal_v4a_"),
        initial_dt=1e-8,
        components=True,
        nonlinear_strategy="quasi_newton",
        nonsmooth_recovery=True,
    )


if __name__ == "__main__":
    execute(parser().parse_args())
