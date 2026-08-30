from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utils.cfd_flow.final_base_geometry_validation import (  # noqa: E402
    validate_final_base,
    write_compact_stage_files,
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run full final-BASE continuous-STL and D3Q19 geometry gates."
    )
    parser.add_argument("base_run", type=Path)
    parser.add_argument("continuous_surface", type=Path)
    args = parser.parse_args()
    result = validate_final_base(
        base_run=args.base_run, continuous_surface=args.continuous_surface
    )
    write_compact_stage_files(args.base_run, result)
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
