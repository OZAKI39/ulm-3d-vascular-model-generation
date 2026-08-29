"""Run the zero-APES frozen-lattice port flux audit."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.port_flux_audit import (  # noqa: E402
    AUDIT_UNRESOLVED,
    BACKFLOW_CONFIRMED,
    INTEGRATION_ARTIFACT,
    run_port_flux_audit,
)


def main() -> int:
    result = run_port_flux_audit(PROJECT_ROOT)
    print(f"Run root: {result['run_root']}")
    print(f"External APES calls: {result['external_apes_executable_calls']}")
    print(f"Seeder calls: {result['seeder_run_count']}")
    print(f"Musubi calls: {result['musubi_run_count']}")
    print(f"Harvester calls: {result['harvester_run_count']}")
    print(f"STATUS: {result['status']}")
    print(f"NEXT: {result['next']}")
    completed = {INTEGRATION_ARTIFACT, BACKFLOW_CONFIRMED, AUDIT_UNRESOLVED}
    return 0 if result["status"] in completed else 2


if __name__ == "__main__":
    raise SystemExit(main())
