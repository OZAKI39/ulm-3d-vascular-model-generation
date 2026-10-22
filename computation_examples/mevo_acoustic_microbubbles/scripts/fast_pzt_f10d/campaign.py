"""Finite sequential campaign with continuous sampled resource admission."""

import argparse
import json
import subprocess
import sys
import threading
from time import perf_counter
from scripts.fast_pzt_f10c.exclusive import snapshot
from .common import BASE, F10C, identity, write, js


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    assert js(BASE / "validate_v1/RESULTS.json")["status"] == "PASS"
    assert js(BASE / "storage_v1/RESULTS.json")["status"] == "PASS"
    own = "mevo_fast_pzt_f10d_" + a.label
    before = snapshot()
    assert not before["gpu_pids"]
    assert all(line.split()[0] == own for line in before["science_running"])
    samples = []
    violations = []
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
                    samples.append(row)
                    log.write(json.dumps(row) + "\n")
                    log.flush()
                    if not row["exclusive_sample"]:
                        violations.append(row)
                except Exception as exc:
                    violations.append(dict(error=repr(exc)))
                stop.wait(0.5)

    thread = threading.Thread(target=monitor, daemon=True)
    thread.start()
    jobs = []

    def child(module, label, *extra):
        phase[0] = module + " " + " ".join(extra)
        start = perf_counter()
        with (out / (str(len(jobs)) + "_" + module + ".log")).open("w") as log:
            proc = subprocess.Popen(
                [sys.executable, "-m", "scripts.fast_pzt_f10d." + module, "--label", label, *extra],
                stdout=log,
                stderr=subprocess.STDOUT,
            )
            while proc.poll() is None:
                if violations:
                    proc.terminate()
                    proc.wait(timeout=20)
                    raise RuntimeError("Resource contention: timing rejected")
                try:
                    proc.wait(timeout=0.01)
                except subprocess.TimeoutExpired:
                    pass
            wall = perf_counter() - start
            jobs.append(
                dict(
                    module=module,
                    label=label,
                    extra=extra,
                    full_process_wall_s=wall,
                    exit_code=proc.returncode,
                )
            )
            write(out / "JOBS.json", jobs)
            assert proc.returncode == 0, jobs[-1]
        print(jobs[-1], flush=True)
        return wall

    try:
        label = a.label + "_benchmark"
        child("benchmark", label)
        from .pipeline import compare_hdf

        resident = js(BASE / label / "RESIDENT_COMPLETE.json")
        rows = resident["rows"]
        for row in rows:
            size = row["B"]
            reference = compare_hdf(
                F10C / f"exclusive_v1_batch/B{size}_batch", BASE / label / f"resident_B{size}_legacy"
            )
            assert reference["status"] == "PASS"
            row["F10C_archive_comparison"] = reference
        for size in [16, 64]:
            modes = {}
            for mode in ["legacy", "optimized"]:
                wall = child("benchmark", label, "--size", str(size), "--mode", mode)
                value = js(BASE / label / f"cold_B{size}_{mode}/COLD.json")
                value["full_process_wall_s"] = wall
                value["cold_valid_samples_per_gpu_hour"] = 3600 * value["valid_sample_count"] / wall
                modes[mode] = value
            check = compare_hdf(
                BASE / label / f"cold_B{size}_legacy", BASE / label / f"cold_B{size}_optimized"
            )
            assert check["status"] == "PASS"
            rows.append(
                dict(
                    B=size,
                    startup="cold",
                    modes=modes,
                    comparison=check,
                    speedup=modes["legacy"]["full_process_wall_s"]
                    / modes["optimized"]["full_process_wall_s"],
                )
            )
        write(BASE / label / "RESULTS.json", dict(status="PASS", rows=rows, training_ready=False))
        child("families", a.label + "_families")
        child("diversity", a.label + "_diversity")
    finally:
        stop.set()
        thread.join(timeout=15)
        after = snapshot()
        write(
            out / "RESOURCE_AUDIT.json",
            dict(
                status="PASS" if samples and not violations and not after["gpu_pids"] else "FAIL",
                before=before,
                after=after,
                sample_count=len(samples),
                violations=violations,
                poll_wait_s=0.5,
                scope="At most one CUDA child at any time; sampled compute PID cardinality and MEVO/Y-clot Supervisor jobs. No exclusion claim for sub-interval external work; telemetry overhead is included.",
            ),
        )
    assert js(out / "RESOURCE_AUDIT.json")["status"] == "PASS"
    cold64 = next(r for r in rows if r["B"] == 64 and r["startup"] == "cold")
    write(
        out / "COMPLETE.json",
        dict(
            status="PASS",
            jobs=jobs,
            cold_B64_speedup=cold64["speedup"],
            engineering_1p25_pass=cold64["speedup"] >= 1.25,
            physics_scope="EXPLORATORY / PHYSICS_UNCALIBRATED",
            training_ready=False,
        ),
    )


if __name__ == "__main__":
    main()
