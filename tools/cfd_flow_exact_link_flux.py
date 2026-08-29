"""Run the zero-executable exact Musubi boundary-link flux audit."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.exact_link_flux import (  # noqa: E402
    MFR_TARGET_MISMATCH,
    MFR_TARGET_PASS,
    run_exact_link_flux_audit,
)


def main() -> int:
    result = run_exact_link_flux_audit(PROJECT_ROOT)
    print(f"Run root: {result['run_root']}")
    print(f"Seeder calls: {result['seeder_run_count']}")
    print(f"Musubi calls: {result['musubi_run_count']}")
    print(f"Harvester calls: {result['harvester_run_count']}")
    print(f"STATUS: {result['status']}")
    print(f"NEXT: {result['next']}")
    return 0 if result["status"] in {MFR_TARGET_PASS, MFR_TARGET_MISMATCH} else 2


if __name__ == "__main__":
    raise SystemExit(main())
