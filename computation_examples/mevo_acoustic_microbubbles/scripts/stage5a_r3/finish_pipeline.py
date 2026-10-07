"""Managed bounded follow-through after formal timing; never overlaps science jobs."""

import subprocess
from time import monotonic, sleep
from scripts.stage5a_r3.common import ROOT, BASE, read, write, sha


def main():
    path = BASE / "final/CAMPAIGN_R3_V2/RESULT.json"
    deadline = monotonic() + 50400
    while not (path.exists() and read(path)["complete"]):
        if monotonic() >= deadline:
            raise TimeoutError("Finite campaign wait expired; no follow-up jobs launched")
        status = subprocess.run(
            ["supervisorctl", "status", "mevo_stage5a_r3_campaign_v2"], capture_output=True, text=True
        ).stdout
        if "RUNNING" not in status and "STARTING" not in status:
            raise RuntimeError("Campaign stopped before completion; follow-through is not a restart")
        sleep(15)
    # Formal campaign is closed. The separate supplement charges the same R2
    # group allocation; every following diagnostic has its own native watchdog.
    tasks = [("protected_repeats", [], 4950)] + [
        ("dynamic_audit", ["--window", w], 330)
        for w in ["initial", "stable", "S0", "event", "B3_global_costly", "B3_serial_costly"]
    ]
    folder = BASE / "diagnosis/finish_pipeline"
    folder.mkdir(exist_ok=False)
    rows = []
    for module, args, cap in tasks:
        label = module + ("_" + args[-1] if args else "")
        log = folder / (label + ".log")
        command = [
            "/usr/bin/timeout",
            "--signal=INT",
            "--kill-after=20",
            str(cap),
            "/bin/bash",
            "scripts/stage5a_r3/python.sh",
            "-m",
            "scripts.stage5a_r3." + module,
            *args,
        ]
        row = dict(module=module, args=args, command=command, log=str(log), status="RUNNING")
        rows.append(row)
        write(folder / "RESULT.json", dict(complete=False, rows=rows, source_sha256=sha(__file__)))
        start = monotonic()
        with log.open("x") as f:
            result = subprocess.run(command, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
        row.update(
            returncode=result.returncode,
            process_wall_s=monotonic() - start,
            status="COMPLETE" if result.returncode == 0 else "FAILED",
        )
        write(folder / "RESULT.json", dict(complete=False, rows=rows, source_sha256=sha(__file__)))
        print(row, flush=True)
        if module == "protected_repeats" and result.returncode != 0:
            raise RuntimeError("Protected supplement did not close; no diagnostic overlap allowed")
    write(
        folder / "RESULT.json",
        dict(
            complete=True,
            rows=rows,
            source_sha256=sha(__file__),
            scope="Follow-through completion only; numerical and accuracy gates are assessed separately",
        ),
    )


if __name__ == "__main__":
    main()
