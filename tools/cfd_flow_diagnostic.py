"""Run the authorized low-cost diagnostic from the frozen Musubi restart."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.diagnostics import DIAGNOSTIC_STATUS, run_diagnostic  # noqa: E402
from utils.cfd_flow.io import FlowError  # noqa: E402


def main() -> int:
    try:
        summary = run_diagnostic(PROJECT_ROOT)
    except (FlowError, FileNotFoundError, ValueError) as error:
        status = error.status if isinstance(error, FlowError) else "CFD_FLOW_DIAGNOSTIC_INPUT_INVALID"
        print(f"CFD FLOW diagnostic failed: {error}")
        print(f"STATUS: {status}")
        return 1
    print(f"Run root: {summary['run_root']}")
    print(f"Restart: {summary['restart_validation']['status']}")
    print(f"Harvester runs: {summary['diagnostic_harvester_run_count']}")
    print(f"Additional iterations: {summary['additional_musubi_iterations']}")
    print(f"Classification: {summary['convergence_classification']}")
    print(f"STATUS: {summary['status']}")
    print(f"NEXT_RECOMMENDATION: {summary['next_recommendation']}")
    return 0 if summary["status"] == DIAGNOSTIC_STATUS else 2


if __name__ == "__main__":
    raise SystemExit(main())

