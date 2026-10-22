"""Sequential gated campaign. Never launches 5 ms unless both 1 ms gates pass."""

from scripts.stage5a.run_probe import run, BASE
from scripts.stage5a.compare import gate
from scripts.stage4l.dense_probe import write


def main():
    for duration in [0.001, 0.005]:
        for case in ["B2", "B3"]:
            for dt in [20e-6, 10e-6]:
                run(case, duration, dt)
        decision = gate(duration)
        write(
            BASE / "CAMPAIGN_STATE.json",
            dict(
                last_duration_s=duration,
                pass_gate=decision["pass_gate"],
                status="STOPPED_AT_GATE" if not decision["pass_gate"] else "STAGE_COMPLETE",
            ),
        )
        if not decision["pass_gate"]:
            break


if __name__ == "__main__":
    main()
