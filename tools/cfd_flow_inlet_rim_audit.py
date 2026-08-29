"""Run the read-only inlet D3Q19 rim-localization audit."""

from __future__ import annotations

import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from utils.cfd_flow.inlet_rim_audit import run_inlet_rim_audit  # noqa: E402


def main() -> int:
    result = run_inlet_rim_audit(PROJECT_ROOT)
    stats = result.get("diagonal_statistics", {})
    print(f"Mesh source: {result['mesh_source_path']}")
    print(f"Total/cardinal/diagonal: {result['total_inlet_d3q19_cells']}/"
          f"{result['negative_z_cardinal_count']}/{result['diagonal_count']}")
    if stats:
        print(f"Diagonal <=2dx: {stats['distance_le_2dx_count']}/"
              f"{stats['distance_le_2dx_fraction']:.6f}")
        print(f"Diagonal >3dx: {stats['distance_gt_3dx_count']}")
    print(f"RIM_LOCALIZED: {result['rim_localized']}")
    print(f"STATUS: {result['status']}")
    print(f"NEXT: {result['next']}")
    return 0 if result["source_frozen_files_unchanged"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
