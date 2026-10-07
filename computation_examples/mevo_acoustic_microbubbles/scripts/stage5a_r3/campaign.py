"""Sequential finite R3 campaign; every subprocess is independently watchdog-bound."""

from datetime import datetime, timezone
import subprocess
import sys
from time import perf_counter
from scripts.stage5a_r3.common import ROOT, BASE, read, write, sha
from scripts.stage5a_r3.freeze import verify


def numerical(summary):
    return bool(
        summary
        and summary["status"] == "COMPLETE"
        and summary["accepted_steps"] > 0
        and summary["accepted_kkt_max"] <= 1e-9
        and summary["accepted_linear_residual_max"] <= 1e-11
        and summary["accepted_nonlinear_kkt_max"] <= 1e-9
        and summary["accepted_pair_swept_gap_min_m"] >= -1e-11
        and summary["accepted_wall_swept_gap_bound_min_m"] >= -1e-11
        and summary["projection_max_m"] == 0
        and summary["source_inputs_unchanged"]
    )


def main():
    protocol = verify()
    assert read(BASE / "final/HOLDOUT_R3_V1.json")["status"] == "PASS"
    out = BASE / "final/CAMPAIGN_R3_V1"
    out.mkdir(exist_ok=False)
    ledger = []
    spent = dict(new=0.0, reference=0.0, r2=0.0)
    initial = perf_counter()

    def save(complete=False):
        write(
            out / "RESULT.json",
            dict(
                complete=complete,
                ledger=ledger,
                spent_s=spent,
                wall_s=perf_counter() - initial,
                protocol_sha256=sha(BASE / "final/VALIDATION_PROTOCOL_R3_V1.json"),
                current_resource_policy="One trajectory subprocess; BLAS=1, QP worker=1, GPU pair force, no concurrent diagnostic jobs",
                promotion="BLOCKED_PENDING_ALL_ACTIVE_INDEPENDENT_ACCURACY",
            ),
        )

    def execute(kind, case, duration, repeat=0):
        verify()
        group = "reference" if kind == "reference" else "r2" if kind == "r2" else "new"
        cap = 4500 if group == "reference" else 2400 if duration == 0.005 else 1200
        available = (
            protocol["budget"][
                {"new": "new_total_s", "reference": "reference_total_s", "r2": "r2_protected_total_s"}[group]
            ]
            - spent[group]
        )
        if available < cap:
            ledger.append(
                dict(
                    kind=kind,
                    case=case,
                    duration_s=duration,
                    repeat=repeat,
                    status="NOT_RUN",
                    reason="PREREGISTERED_GROUP_BUDGET",
                )
            )
            save()
            return None
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        label = f"formal_R3_V1_{duration:g}_{kind}_repeat{repeat}"
        if group == "reference":
            module = "scripts.stage5a_r3.reference"
            args = ["--case", case, "--duration", str(duration), "--repeat", str(repeat), "--pair", "cupy"]
            folder = BASE / "benchmarks/frozen_reference"
        else:
            refined = kind == "refined"
            options = protocol["method"]["refined" if refined else "nominal"]
            module = "scripts.stage5a_r3.run" if group == "new" else "scripts.stage5a_r3.r2_reference"
            args = [
                "--case",
                case,
                "--duration",
                str(duration),
                "--max-dt",
                str(options["max_dt_s"]),
                "--tol-factor",
                str(options["tolerance_factor"]),
                "--initial-dt",
                "1e-8",
                "--pair",
                "cupy",
                "--budget",
                str(cap),
                "--nonlinear-strategy",
                "quasi_newton",
                "--nonsmooth-recovery",
                "--label",
                label if group == "new" else f"R3_protected_{duration:g}_repeat{repeat}",
            ]
            if group == "new":
                args += ["--optimized", "--global-qp", "--s0-dt", str(options["S0_cap_s"])]
            folder = BASE / "runs" if group == "new" else BASE / "benchmarks/r2_protected"
        existing = set(folder.glob("*/SUMMARY.json")) if folder.exists() else set()
        command = [
            "/usr/bin/timeout",
            "--signal=INT",
            "--kill-after=20",
            str(cap + 30),
            sys.executable,
            "-m",
            module,
            *args,
        ]
        log = out / f"{stamp}_{case}_{duration:g}_{kind}_{repeat}.log"
        entry = dict(
            kind=kind,
            case=case,
            duration_s=duration,
            repeat=repeat,
            status="RUNNING",
            command=command,
            log=str(log),
        )
        ledger.append(entry)
        save()
        begin = perf_counter()
        with log.open("x") as f:
            result = subprocess.run(command, cwd=ROOT, stdout=f, stderr=subprocess.STDOUT)
        charged = perf_counter() - begin
        spent[group] += charged
        created = set(folder.glob("*/SUMMARY.json")) - existing
        summary_path = next(iter(created)) if len(created) == 1 else None
        summary = read(summary_path) if summary_path else None
        entry.update(
            returncode=result.returncode,
            charged_wall_s=charged,
            run=str(summary_path.parent) if summary_path else None,
            summary=summary,
            status=summary["status"] if summary else "INCOMPLETE",
            numerical_gate="PASS"
            if group != "reference" and numerical(summary)
            else "NOT_EVALUABLE"
            if group == "reference"
            else "FAIL",
        )
        save()
        print({k: v for k, v in entry.items() if k not in ["summary", "command"]}, flush=True)
        return summary

    save()
    # D precedes E. A failed original-start numerical gate closes its longer runs.
    allowed = {case: numerical(execute("nominal", case, 0.001)) for case in ["B2", "B3"]}
    complete5 = {}
    for case in ["B2", "B3"]:
        complete5[case] = numerical(execute("nominal", case, 0.005)) if allowed[case] else False
    for case in ["B2", "B3"]:
        if complete5[case]:
            execute("refined", case, 0.005)
    # Full frozen references were given a larger finite budget before execution.
    for case in ["B2", "B3"]:
        if complete5[case]:
            execute("reference", case, 0.005, 0)
    for repeat in [1, 2]:
        for case in ["B2", "B3"]:
            if complete5[case]:
                execute("nominal", case, 0.005, repeat)
        for case in ["B2", "B3"]:
            if complete5[case]:
                execute("reference", case, 0.005, repeat)
    # Equal-duration 1 ms reference check and protected R2 B3 regressions.
    for case in ["B2", "B3"]:
        if allowed[case]:
            execute("reference", case, 0.001, 0)
    if allowed["B3"]:
        execute("r2", "B3", 0.001, 0)
        execute("r2", "B3", 0.005, 0)
    save(complete=True)


if __name__ == "__main__":
    main()
