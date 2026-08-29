"""Attach a zero-solver-call exact checkpoint watchdog to one active fast run."""

from __future__ import annotations

import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.healthy_capillary_fast_steady import watch_existing_fast_steady_run  # noqa: E402


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: cfd_flow_fast_steady_checkpoint_watchdog.py RUN_ROOT")
    result = watch_existing_fast_steady_run(PROJECT_ROOT, Path(sys.argv[1]))
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result.get("status") in {
        "EXACT_GATE_PASS_STOP_REQUESTED", "NO_EXACT_GATE_PASS_BEFORE_SOLVER_EXIT",
    } else 1


if __name__ == "__main__":
    raise SystemExit(main())
