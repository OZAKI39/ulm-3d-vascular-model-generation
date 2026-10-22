"""Matched 0.2ms diagnostic cost, equal outputs and checks, no production-rate claim."""

import argparse
from mevo_mb.fast_sim.f10h_lubrication import LubricationDiagnostic
from scripts.fast_pzt_f9.common import objects
from .common import ROOT, BASE, js, write, identity
from .trials import record
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    source = identity()
    write(out / "IDENTITY.json", source)
    obj = objects()
    p, b, r, g, f = obj
    c = js(ROOT / "results/fast_pzt/f10e/SELECTION.json")["branch"]
    cfg = js(ROOT / c["relative_path"] / "CONFIG.json")
    modes = ["M0_LUB2_PBD12_DT10", "M0_LUB4_PBD12_DT10"]
    for mode in modes:
        s = LubricationDiagnostic(
            np.array(cfg["initial_positions_m"]), p, r, b, vessel=g, flow=f, numerical_mode=mode
        )
        s.configure_schedule(cfg["events"], 1e-5)
        s.control(cfg["events"][0])
        s.advance(1e-5, 2, graphs=False)
    rows = []
    for rep in range(3):
        for mode in modes if rep % 2 == 0 else modes[::-1]:
            result = record(out / f"{mode}_rep{rep}", c, obj, mode, 0.0002, source)
            assert result["status"] == "COMPLETE"
            rows.append(result | dict(repetition=rep))
            write(
                out / "RESULTS.json",
                dict(
                    rows=rows,
                    scope="Matched 0.2ms instrumented B=1 runs; includes same physical kernels, quality, 50us HDF5 and step diagnostic output. CUDA span includes CPU/synchronization gaps. No production throughput claim.",
                    warmup_steps_per_mode=2,
                    training_ready=False,
                ),
            )
    print("Matched diagnostic timings complete", flush=True)


if __name__ == "__main__":
    main()
