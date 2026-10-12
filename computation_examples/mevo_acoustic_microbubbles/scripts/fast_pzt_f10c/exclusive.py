"""Bounded performance replay with sampled GPU exclusivity evidence."""

import argparse
from datetime import datetime, timezone
import json
import subprocess
import sys
import threading
import time
from .common import BASE, ROOT, identity, write, js, sha


def snapshot():
    gpu = (
        subprocess.check_output(
            ["nvidia-smi", "--query-compute-apps=pid", "--format=csv,noheader"], text=True
        )
        .strip()
        .splitlines()
    )
    status = subprocess.run(["supervisorctl", "status"], capture_output=True, text=True).stdout
    return dict(
        time_UTC=datetime.now(timezone.utc).isoformat(),
        gpu_pids=sorted(set(gpu)),
        science_running=[
            line
            for line in status.splitlines()
            if "RUNNING" in line and line.startswith(("mevo_fast_", "y_clot_"))
        ],
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    own = "mevo_fast_pzt_f10c_" + a.label
    before = snapshot()
    assert not before["gpu_pids"], before
    assert all(line.split()[0] == own for line in before["science_running"]), before
    write(out / "BEFORE.json", before)
    protocol = ROOT / "configs/fast_pzt/f10c/exclusive_performance.json"
    records, errors = [], []
    stop = threading.Event()
    phase = ["initial"]

    def monitor():
        with (out / "RESOURCE_SAMPLES.jsonl").open("w") as log:
            while not stop.is_set():
                try:
                    row = snapshot()
                    row["phase"] = phase[0]
                    row["exclusive_sample"] = len(row["gpu_pids"]) <= 1 and all(
                        line.split()[0] == own for line in row["science_running"]
                    )
                    records.append(row)
                    log.write(json.dumps(row) + "\n")
                    log.flush()
                    if not row["exclusive_sample"]:
                        errors.append(row)
                except Exception as exc:
                    errors.append(dict(error=repr(exc)))
                stop.wait(0.5)

    thread = threading.Thread(target=monitor, daemon=True)
    thread.start()
    jobs = []
    try:
        for module, extra in [
            ("batch", ["--max-batch", "64"]),
            ("cold", ["--reference-label", a.label + "_batch"]),
        ]:
            phase[0] = module
            label = a.label + "_" + module
            start = time.perf_counter()
            with (out / (module + ".log")).open("w") as log:
                process = subprocess.Popen(
                    [sys.executable, "-m", "scripts.fast_pzt_f10c." + module, "--label", label, *extra],
                    stdout=log,
                    stderr=subprocess.STDOUT,
                )
                while process.poll() is None:
                    if errors:
                        process.terminate()
                        process.wait(timeout=20)
                        raise RuntimeError("Resource contention or monitor failure; retest rejected")
                    time.sleep(0.2)
                assert process.returncode == 0, (module, process.returncode)
            jobs.append(dict(module=module, label=label, wall_s=time.perf_counter() - start))
            print(jobs[-1], flush=True)
    finally:
        stop.set()
        thread.join(timeout=15)
        after = snapshot()
        write(
            out / "RESOURCE_AUDIT.json",
            dict(
                status="PASS" if records and not errors and not after["gpu_pids"] else "FAIL",
                before=before,
                after=after,
                sample_count=len(records),
                violations=errors,
                poll_wait_s=0.5,
                monitor_scope="Sampled CUDA process cardinality <=1 and no other MEVO/Y-clot Supervisor science jobs; this serial campaign creates at most one CUDA context process at a time. Sampling cannot exclude sub-interval external work.",
                max_gpu_processes=max((len(row["gpu_pids"]) for row in records), default=0),
            ),
        )
    assert js(out / "RESOURCE_AUDIT.json")["status"] == "PASS"
    rows = js(BASE / (a.label + "_batch/RESULTS.json"))["rows"]
    assert [r["B"] for r in rows] == [1, 4, 16, 64]
    assert all(r["status"] == "BATCH_NUMERICAL_PASS" for r in rows)
    assert js(BASE / (a.label + "_cold/RESULTS.json"))["status"] == "PASS"
    write(
        out / "COMPLETE.json",
        dict(
            status="PASS",
            jobs=jobs,
            protocol_sha256=sha(protocol),
            resource_scope="EXCLUSIVE_GPU_SAMPLED",
            comparisons=85,
            engineering_2x_pass=all(r["resident_speedup"] >= 2 for r in rows if r["B"] >= 16),
            training_ready=False,
        ),
    )


if __name__ == "__main__":
    main()
