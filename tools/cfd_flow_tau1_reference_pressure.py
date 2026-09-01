#!/usr/bin/env python3
"""Validate the research-only Tau1 numerical reference-pressure scaling."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.tau1_reference_pressure import (  # noqa: E402
    audit_reference_scaled_smoke,
    finalize_reference_scaled_smoke_blocked,
    prepare_reference_pressure_zero_run,
    run_reference_scaled_smoke,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "action",
        choices=("zero-run", "run-smoke", "audit-smoke", "diagnose-blocked"),
    )
    args = parser.parse_args()
    if args.action == "zero-run":
        result = prepare_reference_pressure_zero_run(PROJECT_ROOT)
    elif args.action == "run-smoke":
        result = run_reference_scaled_smoke(PROJECT_ROOT)
    elif args.action == "audit-smoke":
        result = audit_reference_scaled_smoke(PROJECT_ROOT)
    else:
        result = finalize_reference_scaled_smoke_blocked(PROJECT_ROOT)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0 if result.get("status") not in {"FAIL"} else 1


if __name__ == "__main__":
    raise SystemExit(main())
