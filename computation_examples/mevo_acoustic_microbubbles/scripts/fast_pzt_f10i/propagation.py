"""The single justified 2ms tighter-tolerance control from original initial state."""

import argparse
from mevo_mb.fast_sim.f10i_diagnostic import TightContactDiagnostic
from scripts.fast_pzt_f9.common import objects
from scripts.fast_pzt_f10h.trials import record
from .common import ROOT, BASE, js, write, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    decision = js(BASE / "analysis/DEVELOPMENT_DECISIONS.json")
    assert decision["short_tolerance_trajectory"]["admitted"]
    out = BASE / a.label
    out.mkdir(exist_ok=False, parents=True)
    source = identity()
    write(out / "IDENTITY.json", source)
    row = record(
        out / "O_1102_B2_GEOM_K2_TOL14_DT10",
        js(ROOT / "results/fast_pzt/f10e/SELECTION.json")["branch"],
        objects(),
        "B2_GEOM_K2_TOL14_DT10",
        0.002,
        source,
        simulator=TightContactDiagnostic,
        kc_N_m=js(ROOT / "configs/fast_pzt/f10i/protocol.json")["kc_N_m"],
    )
    write(
        out / "RESULTS.json",
        dict(
            rows=[row],
            tolerance_m=1e-14,
            complete_10ms=False,
            training_ready=False,
            scope="Single tolerance-propagation control, same 10us integrator and original complete controls; not a new independent initial state",
        ),
    )
    print(row["status"], row["accepted_steps"], row["first_rejection"], flush=True)


if __name__ == "__main__":
    main()
