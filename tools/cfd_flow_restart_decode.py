"""Decode the frozen Musubi restart directly, without any APES execution."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.restart_decode import (  # noqa: E402
    PORT_REVIEW_STATUS,
    SUCCESS_STATUS,
    run_direct_restart_decode,
)


def main() -> int:
    summary = run_direct_restart_decode(PROJECT_ROOT)
    print(f"Run root: {summary['run_root']}")
    print(f"External APES calls: {summary['external_apes_executable_calls']}")
    print(f"Seeder calls: {summary['seeder_run_count']}")
    print(f"Musubi calls: {summary['musubi_run_count']}")
    print(f"Harvester calls: {summary['harvester_run_count']}")
    print(f"STATUS: {summary['status']}")
    print(f"NEXT: {summary['next']}")
    return 0 if summary["status"] in {SUCCESS_STATUS, PORT_REVIEW_STATUS} else 2


if __name__ == "__main__":
    raise SystemExit(main())
