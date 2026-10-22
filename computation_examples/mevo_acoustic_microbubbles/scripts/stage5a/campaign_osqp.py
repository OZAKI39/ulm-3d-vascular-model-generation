"""New OSQP identities only; never resume the stopped active-set campaign."""

from scripts.stage5a.run_probe import run, BASE
from scripts.stage5a.compare import gate
from scripts.stage4l.dense_probe import write


def main():
    rows = []
    for duration in [0.001, 0.005]:
        complete = True
        for case in ["B2", "B3"]:
            folder, summary = run(case, duration, 20e-6, qp_solver="osqp")
            rows.append(
                dict(
                    run=str(folder),
                    status=summary["status"],
                    case=case,
                    duration_s=duration,
                    macro_dt_s=20e-6,
                )
            )
            write(BASE / "OSQP_CAMPAIGN_ROWS.json", rows)
            complete &= summary["status"] == "COMPLETE"
            if summary["status"] == "STOPPED_BY_USER":
                return
        # Refining an already incomplete nominal trajectory cannot satisfy the
        # joint completion gate. Do not spend another wall budget on that arm.
        if complete:
            for case in ["B2", "B3"]:
                folder, summary = run(case, duration, 10e-6, qp_solver="osqp")
                rows.append(
                    dict(
                        run=str(folder),
                        status=summary["status"],
                        case=case,
                        duration_s=duration,
                        macro_dt_s=10e-6,
                    )
                )
                write(BASE / "OSQP_CAMPAIGN_ROWS.json", rows)
                if summary["status"] == "STOPPED_BY_USER":
                    return
        decision = gate(duration, qp_solver="osqp")
        write(
            BASE / "OSQP_CAMPAIGN_STATE.json",
            dict(
                status="STAGE_COMPLETE" if decision["pass_gate"] else "STOPPED_AT_GATE",
                last_duration_s=duration,
                pass_gate=decision["pass_gate"],
                nominal_complete=complete,
                rows=rows,
                dt_half_skipped_reason=None if complete else "Nominal completion prerequisite failed",
            ),
        )
        if not decision["pass_gate"]:
            break


if __name__ == "__main__":
    main()
