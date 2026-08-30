"""Run the isolated, zero-CFD port grid-sensitivity forensic."""

from __future__ import annotations

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.port_grid_sensitivity import run_port_grid_sensitivity  # noqa: E402


if __name__ == "__main__":
    report = run_port_grid_sensitivity(PROJECT_ROOT)
    print(json.dumps(report, indent=2, ensure_ascii=False))
