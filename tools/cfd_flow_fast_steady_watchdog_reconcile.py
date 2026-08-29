"""Reconcile one completed zero-call fast-steady checkpoint watchdog."""

from __future__ import annotations

import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.healthy_capillary_fast_steady import (  # noqa: E402
    reconcile_existing_fast_steady_watchdog,
)


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: cfd_flow_fast_steady_watchdog_reconcile.py RUN_ROOT")
    result = reconcile_existing_fast_steady_watchdog(PROJECT_ROOT, Path(sys.argv[1]))
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
