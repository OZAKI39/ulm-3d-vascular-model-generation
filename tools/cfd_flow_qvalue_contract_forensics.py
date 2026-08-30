"""Run zero-executable qVal contract forensics on existing meshes."""

from __future__ import annotations

import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.qvalue_contract_forensics import (  # noqa: E402
    run_qvalue_contract_forensics,
)


def main() -> int:
    result = run_qvalue_contract_forensics(PROJECT_ROOT)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0 if result["decision_case"] in {"CASE_A", "CASE_B", "CASE_C"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
