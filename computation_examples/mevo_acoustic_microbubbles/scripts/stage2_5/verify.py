"""Persist test and lint results; nonzero exit retains the failure log."""

from pathlib import Path
import subprocess
import time
from common import new_run, finish


def main():
    out = new_run(
        "audit", {"task": "full test suite and lint"}, ["Existing and Stage2.5 software validation tests"]
    )
    results = []
    for name, command in [
        ("pytest", [".venv/bin/python", "-m", "pytest", "-q"]),
        ("ruff", [".venv/bin/ruff", "check", "src", "tests", "scripts", "benchmarks", "references"]),
    ]:
        start = time.perf_counter()
        with (out / (name + ".log")).open("w") as stream:
            p = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
        results.append(
            dict(
                name=name,
                command=command,
                returncode=p.returncode,
                runtime_s=time.perf_counter() - start,
                log=str(out / (name + ".log")),
            )
        )
        print(name, p.returncode, flush=True)
    result = dict(status="PASS" if all(r["returncode"] == 0 for r in results) else "FAIL", checks=results)
    finish(out, result)
    # Dedicated pointer cannot be confused with later audit analyses.
    from common import atomic_json

    atomic_json(Path("results/stage2_5/audit/SOFTWARE_CHECKS.json"), dict(result, run=str(out)))
    if result["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
