"""Standard-library process envelope with sampled GPU ownership, no physics imports."""

import argparse
import json
from pathlib import Path
import subprocess
import sys
import threading
import time
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--module", required=True)
    a, extra = ap.parse_known_args()
    root = Path("results/fast_pzt/f10h")
    log = root / (a.label + "_PROCESS.json")
    assert not log.exists()
    records = []
    stop = threading.Event()

    def monitor():
        while not stop.is_set():
            result = subprocess.run(
                ["nvidia-smi", "--query-compute-apps=pid,process_name,used_memory", "--format=csv,noheader"],
                capture_output=True,
                text=True,
            )
            records.append(
                dict(monotonic_s=time.monotonic(), processes=result.stdout.strip(), error=result.stderr)
            )
            stop.wait(0.5)

    thread = threading.Thread(target=monitor, daemon=True)
    thread.start()
    start = time.perf_counter()
    proc = subprocess.run(
        [sys.executable, "-m", "scripts.fast_pzt_f10h." + a.module, "--label", a.label, *extra]
    )
    elapsed = time.perf_counter() - start
    stop.set()
    thread.join(timeout=10)
    log.write_text(
        json.dumps(
            dict(
                full_process_wall_s=elapsed,
                exit_code=proc.returncode,
                wrapper_pid=os.getpid(),
                resource_samples=records,
            ),
            indent=2,
        )
        + "\n"
    )
    raise SystemExit(proc.returncode)


if __name__ == "__main__":
    main()
