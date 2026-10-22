"""Compare every saved old/fresh reference frame, with no temporal interpolation."""

import argparse
import h5py
import numpy as np
from scripts.stage5a_r3.common import BASE, R2, read, write, sha
from scripts.stage5a_r3.current import campaign
from scripts.stage5a_r3.paths import local


def main(case):
    old = next(
        r
        for r in read(R2 / "final/CAMPAIGN_V4A/RESULT.json")["ledger"]
        if r["case"] == case and r["duration"] == 0.005 and r["kind"] == "reference" and r.get("repeat") == 0
    )
    new = next(
        r
        for r in campaign()["ledger"]
        if r["case"] == case and r["duration_s"] == 0.005 and r["kind"] == "reference" and r["repeat"] == 0
    )
    assert new["summary"]["status"] == "COMPLETE"
    a_path, b_path = [local(r["run"]) / "trajectories.h5" for r in [old, new]]
    with h5py.File(a_path) as a, h5py.File(b_path) as b:
        n = len(a["time_s"])
        checks = {
            k: bool(np.array_equal(a[k][:], b[k][:n], equal_nan=True))
            for k in ["time_s", "positions_m", "active"]
        }
        checks.update(
            {k: bool(np.array_equal(a[k][:], b[k][:], equal_nan=True)) for k in ["particle_id", "radii_m"]}
        )
        result = dict(
            case=case,
            status="PASS" if all(checks.values()) else "FAIL",
            checks=checks,
            frames_checked=n,
            actual_interval_s=[float(a["time_s"][0]), float(a["time_s"][-1])],
            old_path=str(a_path),
            new_path=str(b_path),
            old_sha256=sha(a_path),
            new_sha256=sha(b_path),
            analysis_script_sha256=sha(__file__),
            scope="Every actually stored old-reference frame compared with the prefix of the new complete reference. No interpolation and no extrapolation beyond the incomplete old reference. This tests replay consistency, not physical trajectory accuracy.",
        )
    write(BASE / "final" / (case + "_FULL_FRESH_PREFIX_REPLAY_AUDIT.json"), result)
    print(result, flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--case", choices=["B2", "B3"], required=True)
    main(p.parse_args().case)
