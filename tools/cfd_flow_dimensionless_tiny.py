from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utils.cfd_flow.dimensionless_geometry_kernel import (  # noqa: E402
    semantic_files_success,
)
from utils.cfd_flow.io import write_json  # noqa: E402
from utils.cfd_flow.qvalue_repair_validation import tiny_cylinder_gate  # noqa: E402


REQUIRED_MESH_FILES = (
    "header.lua",
    "elemlist.lsb",
    "bnd.lua",
    "bnd.lsb",
    "qval.lua",
    "qval.lsb",
)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate a final dimensionless-kernel tiny cylinder mesh."
    )
    parser.add_argument("case", type=Path)
    args = parser.parse_args()
    case = args.case.resolve()
    semantic = semantic_files_success(case / "mesh", REQUIRED_MESH_FILES)
    gate = tiny_cylinder_gate(case / "mesh", case / "source_mesh_summary.json")
    summary = json.loads(
        (case / "source_mesh_summary.json").read_text(encoding="utf-8")
    )
    stdout = (case / "seeder_stdout.log").read_text(
        encoding="utf-8", errors="replace"
    )
    runtime_match = re.search(r"Done with Seeder in\s+([0-9.]+)\s+s", stdout)
    runtime_seconds = float(runtime_match.group(1)) if runtime_match else None
    result = {
        "stage": "tiny_analytic_cylinder",
        "status": "PASS"
        if semantic["semantic_success"] and gate["status"] == "PASS"
        else "FAIL",
        "semantic_output": semantic,
        "analytic_gate": gate,
        "input_hashes": {
            "wall_stl_sha256": summary["wall"]["sha256"],
            "seeder_base_sha": summary["seeder_base_sha"],
            "treelm_base_sha": summary["treelm_base_sha"],
            "patch_sha256": summary["patch_sha256"],
        },
        "binary_hashes": {
            "seeder_sha256": "d7be681ca90da706559a4fd7e8f769fdb8f4303b8508f751077205f8e00cc7ed"
        },
        "solver_calls": {
            "seeder_semantic_calls": 1,
            "small_musubi_semantic_calls": 0,
            "vascular_musubi_semantic_calls": 0,
            "launch_failures": 0,
            "preflight_failures": 0,
        },
        "runtime_seconds": runtime_seconds,
        "first_failure": None if gate["status"] == "PASS" else "ANALYTIC_GATE",
        "recovery_attempts": 0,
        "next": "final_base_seeder" if gate["status"] == "PASS" else "STOP",
    }
    write_json(case / "tiny_cylinder_status.json", result)
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
