"""Local-only copy of explicitly closed campaign runs; no growing HDF5 is read."""

import argparse
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import tempfile
from scripts.stage5a_r3.common import BASE, read, write

REMOTE = "/workspace/mevo_acoustic_microbubbles_stage5a_r3/results/stage5/motion_solver_r3/"


def main(host="vast4090"):
    assert str(BASE).startswith("/home/"), "Run this transfer tool on the local workstation"
    folder = BASE / "final/CAMPAIGN_R3_V2"
    folder.mkdir(exist_ok=True)
    cmd = ["rsync", "-az", "--bwlimit=512", "--rsync-path=nice -n 19 ionice -c 3 rsync"]
    subprocess.run(
        cmd + [host + ":" + REMOTE + "final/CAMPAIGN_R3_V2/RESULT.json", str(folder) + "/"], check=True
    )
    ledger = read(folder / "RESULT.json")["ledger"]
    selected = []
    for row in ledger:
        if not row.get("run") or not row.get("summary") or row["status"] == "RUNNING":
            continue
        relative = row["run"].split("results/stage5/motion_solver_r3/", 1)[1]
        assert ".." not in Path(relative).parts
        names = [
            "SUMMARY.json",
            "IDENTITY.json",
            "trajectories.h5",
            "ACCEPTED.jsonl",
            "RESOURCE_SAMPLES.json",
        ]
        if row["kind"] != "reference":
            names += ["REJECTED.jsonl", "CHECKPOINT.npz", "COMPLETED_TRIALS.json", "INPUT_SNAPSHOT.npz"]
        for name in names:
            # Already completed immutable data is not transferred repeatedly.
            if not (BASE / relative / name).is_file():
                selected.append(relative + "/" + name)
    if not selected:
        print("All currently closed run files are already local")
        return
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    record = dict(
        start_utc=stamp,
        remote_root=REMOTE,
        files=selected,
        policy="Read-only completed files, nice=19, ionice=idle, bandwidth<=512 KiB/s. No scientific analysis runs on the server during formal timing.",
        timing_note="This background artifact transfer is disclosed; it can add I/O overhead. It is not used to subtract any measured formal time.",
    )
    output = BASE / "diagnosis" / ("COMPLETED_SYNC_" + stamp + ".json")
    write(output, record)
    with tempfile.NamedTemporaryFile("w", prefix="mevo-r3-completed-", suffix=".txt") as f:
        f.write("\n".join(selected) + "\n")
        f.flush()
        result = subprocess.run(cmd + ["--files-from=" + f.name, host + ":" + REMOTE, str(BASE) + "/"])
    record.update(
        end_utc=datetime.now(timezone.utc).isoformat(),
        returncode=result.returncode,
        missing_local=[x for x in selected if not (BASE / x).is_file()],
    )
    write(output, record)
    assert result.returncode == 0, record
    print(dict(transferred_files=len(selected), log=str(output)))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--host", default="vast4090")
    main(p.parse_args().host)
