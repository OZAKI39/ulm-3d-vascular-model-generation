"""Run the one-rotation, at-most-one-Seeder inlet mesh preflight."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.axis_aligned_inlet import (  # noqa: E402
    SUCCESS_STATUS,
    run_axis_aligned_inlet_preflight,
)


def main() -> int:
    result = run_axis_aligned_inlet_preflight(PROJECT_ROOT)
    print(f"Run root: {result['run_root']}")
    print(f"Seeder calls: {result['seeder_run_count']}")
    print(f"Seeder return code: {result['seeder_return_code']}")
    print(f"Musubi calls: {result['musubi_run_count']}")
    print(f"Harvester calls: {result['harvester_run_count']}")
    print(f"STATUS: {result['status']}")
    print(f"NEXT: {result['next']}")
    return 0 if result["status"] == SUCCESS_STATUS else 2


if __name__ == "__main__":
    raise SystemExit(main())
