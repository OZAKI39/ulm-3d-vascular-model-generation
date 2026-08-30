"""Run the isolated small Musubi pressure-boundary benchmark."""

from __future__ import annotations

import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.musubi_pressure_bc_benchmark import (  # noqa: E402
    run_pressure_bc_benchmark,
)


if __name__ == "__main__":
    report = run_pressure_bc_benchmark(PROJECT_ROOT)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    raise SystemExit(0 if report["status"] == "PASS" else 2)
