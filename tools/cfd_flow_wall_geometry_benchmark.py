"""Fail-closed entry point for the wall-only periodic stage."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utils.cfd_flow.port_grid_sensitivity import RESEARCH_RUN  # noqa: E402


if __name__ == "__main__":
    qvalue_path = ROOT / "outputs" / "cfd_flow" / RESEARCH_RUN / "qc" / "wall_qvalue_oracle.json"
    qvalue = json.loads(qvalue_path.read_text(encoding="utf-8"))
    if qvalue["status"] != "PASS_NO_QVALUE_CONTRACT_BUG_IDENTIFIED":
        print(
            json.dumps(
                {
                    "status": "NOT_RUN_QVALUE_GATE_FAILED",
                    "seeder_calls": 0,
                    "musubi_calls": 0,
                    "next": "FIX QVALUE CONTRACT BEFORE ANY CFD",
                },
                indent=2,
            )
        )
        raise SystemExit(2)
    raise SystemExit("q-value gate passed; periodic executor is intentionally not auto-launched")
