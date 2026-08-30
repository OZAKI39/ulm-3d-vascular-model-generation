from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utils.cfd_flow.periodic_pipe_force import (  # noqa: E402
    CASES,
    audit_case,
    final_decision,
    prepare_case,
)


def main() -> int:
    parser = argparse.ArgumentParser(description="Manage frozen-kernel Periodic Pipe Force cases")
    parser.add_argument("root", type=Path)
    parser.add_argument("action", choices=("prepare", "audit", "decide"))
    parser.add_argument("--case", choices=tuple(CASES))
    args = parser.parse_args()
    root = args.root.resolve()
    if args.action in {"prepare", "audit"} and args.case is None:
        parser.error("--case is required for prepare/audit")
    if args.action == "prepare":
        result = prepare_case(root, args.case)
    elif args.action == "audit":
        result = audit_case(root / "cases" / args.case)
    else:
        result = final_decision(root)
    print(json.dumps(result, indent=2))
    if args.action == "audit":
        return 0 if result["status"] == "PASS" else 1
    if args.action == "decide" and result["status"].endswith("FAILED"):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
