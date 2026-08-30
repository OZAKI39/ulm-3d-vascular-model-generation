from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utils.cfd_flow.io import write_json  # noqa: E402
from utils.cfd_flow.musubi_wall_force_diagnostics import (  # noqa: E402
    zero_run_baseline_audit,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Zero-run audit of historical Pipe Force evidence")
    parser.add_argument("case_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = zero_run_baseline_audit(args.case_root.resolve())
    write_json(args.output.resolve(), result)
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
