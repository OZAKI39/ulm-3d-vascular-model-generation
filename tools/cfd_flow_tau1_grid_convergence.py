#!/usr/bin/env python3
"""CLI for the isolated repaired Tau=1 three-grid research workflow."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.physical_port_flux import (  # noqa: E402
    build_interior_plane_contract,
    run_physical_port_flux_validation,
)
from utils.cfd_flow.tau1_grid_convergence import (  # noqa: E402
    build_physical_port_plane_contract as build_legacy_physical_port_plane_contract,
    run_base_physical_flux_preflight as run_legacy_base_physical_flux_preflight,
    write_grid_design,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "action",
        choices=(
            "grid-design",
            "physical-plane-contract",
            "base-preflight",
            "legacy-physical-plane-contract",
            "legacy-base-preflight",
        ),
    )
    args = parser.parse_args()
    if args.action == "grid-design":
        result = write_grid_design(PROJECT_ROOT)
    elif args.action == "physical-plane-contract":
        result = build_interior_plane_contract(PROJECT_ROOT)
    elif args.action == "base-preflight":
        result = run_physical_port_flux_validation(PROJECT_ROOT)
    elif args.action == "legacy-physical-plane-contract":
        result = build_legacy_physical_port_plane_contract(PROJECT_ROOT)
    else:
        result = run_legacy_base_physical_flux_preflight(PROJECT_ROOT)
    print(json.dumps(result, indent=2))
    return 0 if result.get("status") == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
