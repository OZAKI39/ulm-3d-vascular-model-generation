"""Run one isolated 8-rank adaptive-flux steady Musubi baseline."""

from __future__ import annotations

import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.adaptive_flux_steady import (  # noqa: E402
    STEADY_PENDING_AUDIT,
    run_adaptive_flux_steady,
)


def main() -> int:
    result = run_adaptive_flux_steady(PROJECT_ROOT)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result.get("status") == STEADY_PENDING_AUDIT else 1


if __name__ == "__main__":
    raise SystemExit(main())
