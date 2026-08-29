"""Run the gated isolated adaptive-flux solver validations."""

from __future__ import annotations

import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.adaptive_flux_validation import (  # noqa: E402
    SMOKE_PASS,
    run_adaptive_flux_validation,
)


def main() -> int:
    result = run_adaptive_flux_validation(PROJECT_ROOT)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result.get("status") == SMOKE_PASS else 1


if __name__ == "__main__":
    raise SystemExit(main())

