"""Publish only completed scoped evidence; history report is required before cylinder admission."""

import hashlib
import json
from pathlib import Path
from scripts.stage3.common import atomic_json


def evidence(path):
    path = Path(path)
    return dict(path=str(path), sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def main():
    root = Path("results/stage3")
    states = {
        k: json.loads((root / k / "STATE.json").read_text())
        for k in ["acoustic", "secondary", "clustering", "flow", "history"]
    }
    assert states["acoustic"]["validation"]["status"] == "PASS"
    for k in ["secondary", "clustering", "flow", "history"]:
        assert states[k]["status"] == "PASS", k
    radial = Path("results/stage2_5/response_table/20261001T220715485503Z")
    assert json.loads((radial / "SUMMARY.json").read_text())["status"] == "PASS_EMPIRICAL_DOMAIN"
    report = Path("docs/HISTORY_FORCE_SENSITIVITY.md")
    assert report.is_file()
    mapping = {
        "ACOUSTIC_PASS": [root / "acoustic/STATE.json"],
        "BUBBLE_RESPONSE_PASS": [
            radial / "SUMMARY.json",
            radial / "HOLDOUTS.json",
            radial / "adaptive_response.h5",
        ],
        "PRIMARY_PASS": [root / "secondary/STATE.json"],
        "SECONDARY_BASIC_PASS": [root / "secondary/STATE.json"],
        "CLUSTERING_TREND_PASS": [root / "clustering/STATE.json"],
        "FLOW_COUPLING_PASS": [root / "flow/STATE.json"],
    }
    result = dict(
        stage="Stage 3 minimal simulator",
        scope="Water, 2.5 MHz, fixed Tu2009-derived radial manifold. Admission does not imply experimental or near-contact validation.",
        history_prerequisite=[evidence(root / "history/STATE.json"), evidence(report)],
        gates={
            name: dict(status="PASS", scopes=["sonovue_cylinder"], evidence=[evidence(p) for p in paths])
            for name, paths in mapping.items()
        },
    )
    atomic_json(root / "audit/PHYSICS_GATES.json", result)
    print("Stage 3 cylinder prerequisites PASS")


if __name__ == "__main__":
    main()
