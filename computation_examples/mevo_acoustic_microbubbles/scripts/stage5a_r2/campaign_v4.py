"""Same V6 kernel, stricter error configuration; all prior costs are retained."""

from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import time
from scripts.stage5a_r2.validation_v4 import verify
from scripts.stage5a_r2.seal import BASE, ROOT
from scripts.stage5a_r2.compare import compare
from scripts.stage4l.dense_probe import write, sha


def numerical(s):
    return bool(
        s
        and s["status"] == "COMPLETE"
        and s["accepted_steps"] > 0
        and s["accepted_kkt_max"] <= 1e-9
        and s["accepted_linear_residual_max"] <= 1e-11
        and (
            s["accepted_nonlinear_kkt_max"] <= 1e-9
            if s["accepted_nonlinear_kkt_max"] is not None
            else s.get("accepted_method_counts", {}).get("S1", 0) == 0
        )
        and s.get("explicit_recovery_bursts", 0) <= 1024
        and (s.get("explicit_recovery_max_dt_s") is None or s["explicit_recovery_max_dt_s"] <= 1e-7)
        and s["accepted_pair_swept_gap_min_m"] >= -1e-11
        and s["accepted_wall_swept_gap_bound_min_m"] >= -1e-11
        and s["projection_max_m"] == 0
        and s["source_inputs_unchanged"]
    )


def main():
    verify()
    holdout = json.loads((BASE / "final/HOLDOUT_V4.json").read_text())
    assert holdout["status"] == "PASS"
    assert json.loads((BASE / "final/V3_PREVIOUS_FAILURE_REGRESSIONS.json").read_text())["status"] == "PASS"
    assert json.loads((BASE / "final/GEOMETRY_INPUT_AUDIT.json").read_text())["status"] == "PASS"
    # A failed *motion* holdout remains a promotion blocker. It does not erase
    # the separately passed original-QP gate or prohibit scoped original-3s
    # validation. No solver parameters are changed after opening these states.
    assert all(
        max(r["relative_velocity_difference"]) <= 1e-7 and max(r["relative_reaction_difference"]) <= 1e-7
        for r in holdout["rows"]
    )
    out = BASE / "final/CAMPAIGN_V4"
    out.mkdir(exist_ok=False)
    ledger = []
    stages = {}
    prior = json.loads((BASE / "final/CAMPAIGN_V3/RESULT.json").read_text())
    assert prior["complete"]
    spent = prior["spent_s"].copy()
    write(
        out / "PRIOR_BUDGET.json",
        dict(spent_s=spent.copy(), source=str(BASE / "final/CAMPAIGN_V3/RESULT.json")),
    )

    def execute(case, duration, kind, repeat=0):
        verify()
        if duration == 0.001 and repeat == 0 and kind in ["nominal", "reference"]:
            source_kind = "refined" if kind == "nominal" else "reference"
            entry = next(
                e
                for e in prior["ledger"]
                if e["case"] == case
                and e["kind"] == source_kind
                and e["duration"] == duration
                and e["repeat"] == 0
            )
            assert entry["summary"]["status"] == "COMPLETE"
            identity = json.loads((Path(entry["run"]) / "IDENTITY.json").read_text())
            if kind == "nominal":
                assert identity["tolerance_factor"] == 0.25 and identity["max_dt_s"] == 5e-6
                assert numerical(entry["summary"])
            ledger.append(
                dict(
                    entry,
                    kind=kind,
                    reused_from_campaign="V3",
                    original_kind=source_kind,
                    runtime_already_charged=True,
                    reuse_reason="Identical V6 kernel, environment and configuration; the measured run is reused without claiming reexecution.",
                )
            )
            write(out / "EXECUTIONS.json", dict(ledger=ledger, spent_s=spent, stages=stages))
            print(json.dumps(dict(event="REUSE", case=case, kind=kind, run=entry["run"])), flush=True)
            return entry["summary"]
        bucket = "reference" if kind == "reference" else "new"
        remaining = 10800 - spent[bucket]
        label = f"formal_v4_{duration:g}_{kind}_repeat{repeat}"
        if remaining < 60:
            ledger.append(
                dict(case=case, duration=duration, kind=kind, status="NOT_RUN", reason="AGGREGATE_BUDGET")
            )
            return None
        if kind == "reference":
            cmd = [
                sys.executable,
                "-u",
                "-m",
                "scripts.stage5a_r2.frozen_reference",
                "--case",
                case,
                "--duration",
                str(duration),
                "--pair",
                "cupy",
                "--repeat",
                str(repeat),
            ]
            folder = BASE / "benchmarks/frozen_reference"
        else:
            cmd = [
                sys.executable,
                "-u",
                "-m",
                "scripts.stage5a_r2.run",
                "--case",
                case,
                "--duration",
                str(duration),
                "--max-dt",
                str(2.5e-6 if kind == "refined" else 5e-6),
                "--tol-factor",
                str(0.0625 if kind == "refined" else 0.25),
                "--pair",
                "cupy",
                "--nonsmooth-recovery",
                "--nonlinear-strategy",
                "quasi_newton",
                "--label",
                label,
                "--budget",
                str(min(1800, remaining)),
            ]
            folder = BASE / "runs"
        before = set(folder.glob("*"))
        begin = time.perf_counter()
        log = out / f"{case}_{label}.log"
        print(
            json.dumps(
                dict(event="START", case=case, kind=kind, duration_s=duration, repeat=repeat, log=str(log))
            ),
            flush=True,
        )
        with log.open("x") as stream:
            process = subprocess.Popen(
                cmd, cwd=ROOT, stdout=stream, stderr=subprocess.STDOUT, start_new_session=True
            )
            timed_out = False
            try:
                code = process.wait(timeout=min(1830, remaining + 15))
            except subprocess.TimeoutExpired:
                import os
                import signal

                timed_out = True
                os.killpg(process.pid, signal.SIGINT)
                try:
                    code = process.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGKILL)
                    code = process.wait()
        wall = time.perf_counter() - begin
        spent[bucket] += wall
        new = [p for p in set(folder.glob("*")) - before if p.is_dir()]
        assert len(new) <= 1, "unexpected concurrent experiment"
        summary = (
            json.loads((new[0] / "SUMMARY.json").read_text())
            if new and (new[0] / "SUMMARY.json").exists()
            else None
        )
        entry = dict(
            case=case,
            duration=duration,
            kind=kind,
            repeat=repeat,
            wall_s_with_process_start=wall,
            returncode=code,
            watchdog=timed_out,
            run=str(new[0]) if new else None,
            summary=summary,
            command=cmd,
            log=str(log),
            log_sha256=sha(log),
        )
        ledger.append(entry)
        write(out / "EXECUTIONS.json", dict(ledger=ledger, spent_s=spent, stages=stages))
        print(
            json.dumps(
                dict(
                    event="END",
                    case=case,
                    kind=kind,
                    duration_s=duration,
                    status=summary["status"] if summary else "MISSING_SUMMARY",
                    wall_s=wall,
                )
            ),
            flush=True,
        )
        return summary

    for duration in [0.001, 0.005]:
        runs = {case: {} for case in ["B2", "B3"]}
        # Complete both nominal trajectories first; diagnostics/reference have
        # their own processes and never compete with this timed process.
        for kind in ["nominal", "refined", "reference"]:
            for case in runs:
                runs[case][kind] = execute(case, duration, kind)
        stage = dict(duration_s=duration, cases={}, status="INCOMPLETE")
        for case, items in runs.items():
            checks = dict(
                nominal_numerical=numerical(items["nominal"]), refined_numerical=numerical(items["refined"])
            )
            a, b, ref = (items[k] for k in ["nominal", "refined", "reference"])
            if checks["nominal_numerical"] and checks["refined_numerical"]:
                checks["actual_grid_refined"] = b["accepted_steps"] > a["accepted_steps"]
                checks["self_convergence"] = compare(
                    Path(b["run"]) / "trajectories.h5",
                    Path(a["run"]) / "trajectories.h5",
                    BASE / "comparisons",
                    f"{case}_{duration:g}_v4_self",
                )
            if checks["nominal_numerical"] and ref and ref["status"] == "COMPLETE":
                checks["legacy_comparison"] = compare(
                    Path(ref["run"]) / "trajectories.h5",
                    Path(a["run"]) / "trajectories.h5",
                    BASE / "comparisons",
                    f"{case}_{duration:g}_v4_legacy",
                )
            if checks["refined_numerical"] and ref and ref["status"] == "COMPLETE":
                checks["refined_legacy_comparison"] = compare(
                    Path(ref["run"]) / "trajectories.h5",
                    Path(b["run"]) / "trajectories.h5",
                    BASE / "comparisons",
                    f"{case}_{duration:g}_v4_refined_legacy",
                )
            passed = all(
                [
                    checks["nominal_numerical"],
                    checks["refined_numerical"],
                    checks.get("actual_grid_refined", False),
                    checks.get("self_convergence", {}).get("pass_scientific", False),
                    checks.get("legacy_comparison", {}).get("pass_scientific", False),
                    checks.get("refined_legacy_comparison", {}).get("pass_scientific", False),
                ]
            )
            checks["status"] = (
                "PASS"
                if passed
                else "BASELINE_DISAGREEMENT"
                if checks.get("self_convergence", {}).get("pass_scientific", False)
                and checks.get("legacy_comparison", {}).get("status") == "FAIL"
                else "INCOMPLETE_OR_FAIL"
            )
            stage["cases"][case] = checks
        stage["status"] = (
            "PASS" if all(r["status"] == "PASS" for r in stage["cases"].values()) else "INCOMPLETE_OR_FAIL"
        )
        stages[str(duration)] = stage
        write(out / "EXECUTIONS.json", dict(ledger=ledger, spent_s=spent, stages=stages))
        print(json.dumps(stage), flush=True)
        if stage["status"] != "PASS":
            break
    # Repeats use the longest fully compared duration, only when the remaining
    # fixed category budget can plausibly complete another run (20% margin).
    longest = max((float(d) for d, s in stages.items() if s["status"] == "PASS"), default=None)
    if longest:
        for repeat in [1, 2]:
            for kind in ["nominal", "reference"]:
                bucket = "reference" if kind == "reference" else "new"
                for case in ["B2", "B3"]:
                    matching = [
                        e
                        for e in ledger
                        if e.get("summary")
                        and e["case"] == case
                        and e["kind"] == kind
                        and e["duration"] == longest
                        and e["summary"]["status"] == "COMPLETE"
                    ]
                    expected = max((e["wall_s_with_process_start"] for e in matching), default=1800)
                    if 10800 - spent[bucket] >= 1.2 * expected:
                        execute(case, longest, kind, repeat)
    write(
        out / "RESULT.json",
        dict(
            time_utc=datetime.now(timezone.utc).isoformat(),
            ledger=ledger,
            spent_s=spent,
            stages=stages,
            complete=True,
            protocol_sha256=sha(BASE / "final/VALIDATION_PROTOCOL_V4.json"),
        ),
    )


if __name__ == "__main__":
    main()
