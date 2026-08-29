"""Run the source-proven boundary mass referee on the frozen healthy case."""

from __future__ import annotations

import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.musubi_boundary_mass_referee import run_referee  # noqa: E402


def main() -> int:
    result = run_referee(PROJECT_ROOT, run_one_step=True)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result["status"].endswith("PASS") else 1


if __name__ == "__main__":
    raise SystemExit(main())

