"""Run the isolated MPI-selected, exact-checkpoint healthy steady continuation."""

from __future__ import annotations

import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.healthy_capillary_fast_steady import (  # noqa: E402
    STATUS_PASS,
    STATUS_WALLCLOCK,
    run_healthy_capillary_fast_steady,
)


def main() -> int:
    result = run_healthy_capillary_fast_steady(PROJECT_ROOT)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result.get("status") in {STATUS_PASS, STATUS_WALLCLOCK} else 1


if __name__ == "__main__":
    raise SystemExit(main())
