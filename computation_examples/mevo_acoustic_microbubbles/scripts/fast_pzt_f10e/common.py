from pathlib import Path
import subprocess
from scripts.fast_pzt_f9.common import clean, write, js, sha, objects  # noqa: F401

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/fast_pzt/f10e"
HISTORY = ROOT / "results/fast_pzt/f10d"
CAMPAIGN = HISTORY / "campaign_v2_diversity"
BASELINE = "f5c397a15498b14cf281d6e6fa875649fa9e43ca"


def identity():
    marker = ROOT / "SOURCE_F10E_COMMIT.txt"
    commit = (
        marker.read_text().strip()
        if marker.exists()
        else subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    )
    return dict(
        execution_commit=commit,
        baseline_commit=BASELINE,
        protocol_sha256=sha(ROOT / "configs/fast_pzt/f10e/protocol.json"),
        inputs=js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json"),
        source_sha256={
            str(p.relative_to(ROOT)): sha(p)
            for folder in [
                "src/mevo_mb/fast_sim",
                "scripts/fast_pzt_f10c",
                "scripts/fast_pzt_f10d",
                "scripts/fast_pzt_f10e",
            ]
            for p in sorted((ROOT / folder).glob("*.py"))
        },
    )


def original_records(control):
    cfg = js(CAMPAIGN / control / "CONFIG.json")
    summary = js(CAMPAIGN / control / "SUMMARY.json")
    return [
        dict(
            sample_id=sid,
            positions_m=x,
            events=cfg["events"],
            split_group=family,
            initial_state_family=family,
            sample_class="EXPLORATORY_NUMERICAL_SAMPLE",
            initial_state_source=dict(
                parent_run=control,
                parent_member=i,
                parent_config_sha256=sha(CAMPAIGN / control / "CONFIG.json"),
                parent_hdf5_sha256=sha(CAMPAIGN / control / "trajectory.h5"),
            ),
        )
        for i, (sid, x, family) in enumerate(
            zip(summary["sample_id"], cfg["initial_positions_m"], summary["initial_state_family"])
        )
    ]
