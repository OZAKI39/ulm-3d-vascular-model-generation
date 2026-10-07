"""Bind a fresh manufactured recheck to the final frozen kernel, after timings."""

from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from functools import partial
import json
from pathlib import Path
import time
import numpy as np
from scripts.stage5a_r2.manufactured import experiment
from scripts.stage5a_r2.seal import ROOT, BASE
from scripts.stage5a_r2.freeze import critical_hashes
from scripts.stage5a_r2.validation_v4a import verify
from scripts.stage4l.dense_probe import write, sha


def main():
    protocol = verify()
    assert (BASE / "final/CAMPAIGN_V4A/RESULT.json").exists(), "Run after isolated timing campaign"
    out = BASE / "diagnosis/manufactured_final"
    out.mkdir(exist_ok=False)
    identity = dict(
        time_utc=datetime.now(timezone.utc).isoformat(),
        source_commit=protocol["source_commit"],
        source_commit_scope="Frozen numerical protocol identity; postprocessing drivers are identified separately by their actual content hashes.",
        kernel_sha256=critical_hashes(),
        driver_sha256={
            str(p.relative_to(ROOT)): sha(p)
            for p in [Path(__file__), ROOT / "scripts/stage5a_r2/manufactured.py"]
        },
        environment_sha256=sha(BASE / "environment/IDENTITY.json"),
        method="S1 endpoint resistance, explicit constant forcing, start contact tangents; no S0 needed for smooth manufactured geometry",
        workers=3,
        blas_threads=1,
        role="Post-campaign verification only; not tuning or performance evidence",
    )
    write(out / "IDENTITY.json", identity)
    started = time.perf_counter()
    tasks = [
        (wall, f, g)
        for wall in [False, True]
        for f, g in [(1.0, 1.0), (0.25, 1.0), (0.0625, 1.0), (1.0, 0.5), (1.0, 0.25)]
    ]
    with ProcessPoolExecutor(max_workers=3) as pool:
        rows = list(pool.map(partial(experiment, output=out), tasks))
    write(out / "SUMMARY.json", rows)
    orders = {}
    for case in ["pair", "wall"]:
        r = sorted(
            [x for x in rows if x["case"] == case and x["tolerance_factor"] == 1],
            key=lambda x: -x["max_dt_s"],
        )
        err = np.array([x["endpoint_error_m"] for x in r])
        orders[case] = dict(
            actual_steps=[x["accepted_steps"] for x in r],
            observed_orders=np.log2(err[:-1] / err[1:]).tolist(),
            error_decreases=bool(np.all(np.diff(err) < 0)),
        )
    passed = all(
        r["completed_s"] >= 0.001 - 1e-14
        and r["kkt_max"] <= 1e-9
        and r["nonlinear_kkt_max"] <= 1e-9
        and r["min_gap_m"] >= -1e-11
        and r["projection_max_m"] == 0
        for r in rows
    )
    passed &= all(r["error_decreases"] and min(r["observed_orders"]) > 0.8 for r in orders.values())
    assert critical_hashes() == identity["kernel_sha256"]
    write(
        out / "GATE.json",
        dict(
            status="PASS" if passed else "FAIL",
            wall_s=time.perf_counter() - started,
            orders=orders,
            source_identity_sha256=sha(out / "IDENTITY.json"),
            scope="Independent scalar ODE accuracy and actual grid refinement; development v2 evidence remains untouched.",
        ),
    )
    print(json.dumps({"status": "PASS" if passed else "FAIL", "orders": orders}), flush=True)
    assert passed


if __name__ == "__main__":
    main()
