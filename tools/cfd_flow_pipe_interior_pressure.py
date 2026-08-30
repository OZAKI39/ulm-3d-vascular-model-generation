"""Fail-closed entry point for the downstream interior-pressure stage."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utils.cfd_flow.port_grid_sensitivity import RESEARCH_RUN  # noqa: E402


if __name__ == "__main__":
    decision_path = (
        ROOT / "outputs" / "cfd_flow" / RESEARCH_RUN / "qc" / "isolation_root_cause_decision.json"
    )
    decision = json.loads(decision_path.read_text(encoding="utf-8"))
    print(
        json.dumps(
            {
                "status": "NOT_RUN_WALL_GATE_NOT_REACHED",
                "root_cause": decision["root_cause_final"],
                "musubi_calls": 0,
                "next": decision["next"],
            },
            indent=2,
        )
    )
    raise SystemExit(2)
