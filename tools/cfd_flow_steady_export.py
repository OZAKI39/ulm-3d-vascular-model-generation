"""Export and validate the frozen project-steady Musubi field once."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.io import FlowError  # noqa: E402
from utils.cfd_flow.steady_export import (  # noqa: E402
    SUCCESS_STATUS,
    run_steady_field_export,
)


def main() -> int:
    try:
        summary = run_steady_field_export(PROJECT_ROOT)
    except (FlowError, FileNotFoundError, ValueError) as error:
        status = error.status if isinstance(error, FlowError) else "CFD_FLOW_STEADY_FIELD_EXPORT_INVALID"
        print(f"CFD FLOW steady field export failed: {error}")
        print(f"STATUS: {status}")
        print("NEXT: REVIEW STEADY FIELD EXPORT EVIDENCE WITHOUT RERUN")
        return 1
    print(f"Run root: {summary['run_root']}")
    print(f"Seeder calls: {summary['seeder_run_count']}")
    print(f"Musubi calls: {summary['musubi_run_count']}")
    print(f"Harvester calls: {summary['harvester_run_count']}")
    print(f"STATUS: {summary['status']}")
    print(f"NEXT: {summary['next']}")
    return 0 if summary["status"] == SUCCESS_STATUS else 2


if __name__ == "__main__":
    raise SystemExit(main())
