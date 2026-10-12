"""Separate pure decomposition and bounded OFF/S/O q=2 experiments."""

import argparse
from time import perf_counter
import numpy as np
from .common import BASE, ROOT, objects, case, identity, write, js
from .runner import run
from .batch import compare


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    obj = objects()
    cfg = case("O_1000")
    x = np.array(cfg["environment"]["positions_m"])
    timing = run(out / "decomposition_O", x, cfg["events"], obj, decomposition=True)
    equality = compare(BASE / "batch_v1/B1_serial_000", out / "decomposition_O", 0)
    assert equality["bitwise"], equality
    write(
        out / "DECOMPOSITION.json",
        dict(
            status="PASS",
            comparison=equality,
            timing=timing,
            reconstruction=js(out / "decomposition_O/SUMMARY.json")["decomposition_max_residual"],
            scope="Every physical step diagnostic Xc+r and Ubar+Urel; original predictor consumes unchanged balanced U. No separate long-term state, no rigid template, no acceleration claim.",
        ),
    )
    rows = []
    protocol = js(ROOT / "configs/fast_pzt/f10c/protocol.json")
    for mode, (qa, qf) in protocol["modes"].items():
        for key in ["OFF", "S_1000", "O_1000"]:
            if mode == "Q0" and key == "O_1000":
                rows.append(
                    dict(mode=mode, control=key, path="batch_v1/B1_batch", source="already measured baseline")
                )
                continue
            cfg = case(key)
            path = out / (mode + "_" + key)
            start = perf_counter()
            timing = run(path, x, cfg["events"], obj, field_stride=qa, flow_stride=qf)
            summary = js(path / "SUMMARY.json")
            row = dict(
                mode=mode,
                control=key,
                path=str(path.relative_to(BASE)),
                timing=timing,
                outer_wall_s=perf_counter() - start,
                status=summary["status"],
                reason=summary["reason"],
                cache_statistics=summary.get("cache_statistics"),
            )
            rows.append(row)
            write(
                out / "RESULTS.json",
                dict(
                    rows=rows,
                    scope="Accuracy experiments; recorded timings do not qualify a rejected mode as an acceleration.",
                    training_ready=False,
                ),
            )
            print(
                dict(
                    mode=mode,
                    control=key,
                    status=summary["status"],
                    wall_s=row["outer_wall_s"],
                    cache=summary.get("cache_statistics"),
                ),
                flush=True,
            )
    write(
        out / "COMPLETE.json",
        dict(status="COMPLETE", tested=["Q0", "Q1", "Q2", "Q3"], q4="PENDING_Q2_NUMERICAL_ADMISSION"),
    )


if __name__ == "__main__":
    main()
