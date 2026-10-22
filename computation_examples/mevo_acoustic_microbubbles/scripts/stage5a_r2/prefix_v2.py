"""Development prefixes only. Continue from atomic accepted checkpoints."""

from scripts.stage5a_r2.run import run
from scripts.stage5a_r2.seal import BASE
from scripts.stage4l.dense_probe import write

if __name__ == "__main__":
    rows = []
    for case in ["B2", "B3"]:
        p, s = run(case, 1e-6, method="S0", label="S0_diagnostic_v2", budget=180)
        rows.append(s)
        write(BASE / "diagnosis/PREFIX_V2_SUMMARY.json", rows)
        prior = None
        for end in [1e-6, 1e-5, 1e-4]:
            p, s = run(case, end, method="S1", label="prefix_v2", budget=600, resume=prior)
            rows.append(s)
            write(BASE / "diagnosis/PREFIX_V2_SUMMARY.json", rows)
            if s["status"] != "COMPLETE":
                break
            prior = p / "CHECKPOINT.npz"
