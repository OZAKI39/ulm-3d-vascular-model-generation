"""Two extra protected B3 R2 timing repeats, within the original 5400 s group cap.

This presealed validation supplement changes no candidate, benchmark method,
resource setting or numerical threshold. Original campaign ledgers stay intact.
"""

import argparse
from datetime import datetime, timezone
from pathlib import Path
import subprocess
from time import perf_counter
from scripts.stage5a_r3.common import BASE, read, write, sha
from scripts.stage5a_r3.campaign import numerical

PLAN = BASE / "final/PROTECTED_REPEAT_PLAN.json"
EXECUTABLE = Path(
    "/workspace/mevo_acoustic_microbubbles_stage5a_r3/.execution/b5079ad19be81ffe9e3d4f5b7c548717b3d43d67"
)


def freeze():
    assert not PLAN.exists()
    write(
        PLAN,
        dict(
            status="FROZEN",
            created_utc=datetime.now(timezone.utc).isoformat(),
            case="B3",
            duration_s=0.005,
            repeats=[1, 2],
            original_protocol_sha256=sha(BASE / "final/VALIDATION_PROTOCOL_R3_V2.json"),
            wrapper_sha256=sha(__file__),
            original_group_cap_s=5400,
            per_run_cap_s=2400,
            source_commit="b5079ad19be81ffe9e3d4f5b7c548717b3d43d67",
            purpose="Complete three protected R2 B3 timings for R3-vs-R2 repeatability, using only the already allocated R2 group budget",
            prerequisites="Original V2 campaign closed; protected B3 repeat 0 complete and numerical PASS; no concurrent scientific job",
            scope="No new B2 R2 failed-interval restart, no changed method, no budget extension or reset",
        ),
    )


def main():
    plan = read(PLAN)
    assert plan["wrapper_sha256"] == sha(__file__)
    assert plan["original_protocol_sha256"] == sha(BASE / "final/VALIDATION_PROTOCOL_R3_V2.json")
    original = read(BASE / "final/CAMPAIGN_R3_V2/RESULT.json")
    assert original["complete"], "Wait until all original formal timings finish"
    first = next(
        r for r in original["ledger"] if r["kind"] == "r2" and r["case"] == "B3" and r["duration_s"] == 0.005
    )
    out = BASE / "final/CAMPAIGN_R2_PROTECTED_REPEATS"
    out.mkdir(exist_ok=False)
    rows, spent = [], original["spent_s"]["r2"]

    def save(complete=False):
        write(
            out / "RESULT.json",
            dict(
                complete=complete,
                ledger=rows,
                initial_spent_r2_s=original["spent_s"]["r2"],
                total_spent_r2_s=spent,
                original_campaign_sha256=sha(BASE / "final/CAMPAIGN_R3_V2/RESULT.json"),
                plan_sha256=sha(PLAN),
            ),
        )

    save()
    for repeat in plan["repeats"]:
        row = dict(kind="r2", case="B3", duration_s=0.005, repeat=repeat, status="NOT_RUN")
        rows.append(row)
        if not numerical(first.get("summary")):
            row["reason"] = "PROTECTED_REPEAT0_NOT_NUMERICALLY_COMPLETE"
            save()
            continue
        if spent + plan["per_run_cap_s"] + 50 > plan["original_group_cap_s"]:
            row["reason"] = "ORIGINAL_R2_GROUP_BUDGET"
            save()
            continue
        subprocess.run(
            ["/bin/bash", "scripts/stage5a_r3/python.sh", "-m", "scripts.stage5a_r3.freeze_v2", "--verify"],
            cwd=EXECUTABLE,
            check=True,
        )
        run_root = BASE / "benchmarks/r2_protected"
        before = set(run_root.glob("*/SUMMARY.json"))
        command = [
            "/usr/bin/timeout",
            "--signal=INT",
            "--kill-after=20",
            "2430",
            "/bin/bash",
            "scripts/stage5a_r3/python.sh",
            "-m",
            "scripts.stage5a_r3.r2_reference",
            "--case",
            "B3",
            "--duration",
            "0.005",
            "--max-dt",
            "5e-6",
            "--tol-factor",
            "0.25",
            "--initial-dt",
            "1e-8",
            "--pair",
            "cupy",
            "--budget",
            "2400",
            "--nonlinear-strategy",
            "quasi_newton",
            "--nonsmooth-recovery",
            "--label",
            f"R3_protected_0.005_repeat{repeat}",
        ]
        log = out / f"B3_0.005_r2_repeat{repeat}.log"
        row.update(status="RUNNING", command=command, log=str(log))
        save()
        start = perf_counter()
        with log.open("x") as f:
            result = subprocess.run(command, cwd=EXECUTABLE, stdout=f, stderr=subprocess.STDOUT)
        wall = perf_counter() - start
        spent += wall
        created = set(run_root.glob("*/SUMMARY.json")) - before
        path = next(iter(created)) if len(created) == 1 else None
        summary = read(path) if path else None
        row.update(
            status=summary["status"] if summary else "INCOMPLETE",
            summary=summary,
            run=str(path.parent) if path else None,
            charged_wall_s=wall,
            returncode=result.returncode,
            numerical_gate="PASS" if numerical(summary) else "FAIL",
        )
        save()
        print({k: v for k, v in row.items() if k not in ["summary", "command"]}, flush=True)
    save(True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--freeze", action="store_true")
    if parser.parse_args().freeze:
        freeze()
    else:
        main()
