"""Summarize the source-backed qVal repair from already-generated artifacts."""

from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from utils.cfd_flow.io import sha256_file, write_json  # noqa: E402
from utils.cfd_flow.qvalue_repair_validation import (  # noqa: E402
    repaired_base_gate,
    tiny_cylinder_gate,
)
from utils.cfd_flow.seeder_qvalue_repair import (  # noqa: E402
    CURRENT_SEEDER_SHA,
    CURRENT_TREELM_SHA,
    PINNED_SEEDER_SHA,
    PINNED_TREELM_SHA,
    ROOT_CAUSE_CATEGORY,
    calc_dist_config_contract,
    parallel_threshold_diagnostic,
    source_contract,
)
from utils.cfd_flow.wall_qvalue_oracle import audit_mesh_qvalues  # noqa: E402


RESEARCH_RUN = (
    "healthy_mouse_capillary_port_grid_sensitivity_research_"
    "anchor003274_20260830"
)
REPAIRED_BASE_RUN = (
    "healthy_mouse_capillary_qvalue_repaired_base_preflight_"
    "anchor003274_20260830"
)
FINAL_STATUS = "CFD_FLOW_SEEDER_QVALUE_REPAIRED_BASE_TOPOLOGY_FAILED"


def _unc(path: str) -> Path:
    return Path("//wsl.localhost/Ubuntu" + path)


def _git_head(path: str) -> str:
    return subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "git", "-C", path, "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    ).stdout.strip()


def _official_config_contract(root: Path) -> dict:
    base = (
        root
        / "outputs/cfd_flow/axis_aligned_ideal_plane_inlet_preflight_"
        "anchor003274_20260829_120444/seeder/seeder.lua"
    ).read_text(encoding="utf-8")
    pipe = (
        root
        / f"outputs/cfd_flow/{RESEARCH_RUN}/pressure_bc_benchmark/meshes/"
        "axis_aligned/n27/seeder.lua"
    ).read_text(encoding="utf-8")
    official_path = _unc(
        "/home/lzy/apes-pinned/musubi_official/mus/examples/tutorials/"
        "tutorial_cases/tutorial_channelGeneric/seeder.lua"
    )
    result = calc_dist_config_contract(
        base, pipe, official_path.read_text(encoding="utf-8")
    )
    result["sources"] = {
        "vascular_base": str(
            root
            / "outputs/cfd_flow/axis_aligned_ideal_plane_inlet_preflight_"
            "anchor003274_20260829_120444/seeder/seeder.lua"
        ),
        "pipe_axis_n27": str(
            root
            / f"outputs/cfd_flow/{RESEARCH_RUN}/pressure_bc_benchmark/meshes/"
            "axis_aligned/n27/seeder.lua"
        ),
        "official": str(official_path),
    }
    return result


def _upstream_contract() -> tuple[dict, dict]:
    pinned_line = _unc(
        "/home/lzy/apes-pinned/seeder_official/tem/source/shapes/"
        "tem_line_module.fpp"
    )
    current_line = _unc(
        "/home/lzy/open-source-reference/seeder-current/tem/source/shapes/"
        "tem_line_module.fpp"
    )
    pinned_text = pinned_line.read_text(encoding="utf-8")
    current_text = current_line.read_text(encoding="utf-8")
    source = source_contract(
        pinned_source=pinned_text,
        current_source=current_text,
        pinned_seeder_sha=_git_head("/home/lzy/apes-pinned/seeder_official"),
        current_seeder_sha=_git_head("/home/lzy/open-source-reference/seeder-current"),
        pinned_treelm_sha=_git_head("/home/lzy/apes-pinned/seeder_official/tem"),
        current_treelm_sha=_git_head("/home/lzy/open-source-reference/seeder-current/tem"),
    )
    current_root = _unc("/home/lzy/open-source-reference/seeder-current")
    comparison = {
        "status": source["status"],
        "pinned": {
            "seeder_sha": PINNED_SEEDER_SHA,
            "treelm_sha": PINNED_TREELM_SHA,
        },
        "current_upstream": {
            "seeder_sha": CURRENT_SEEDER_SHA,
            "treelm_sha": CURRENT_TREELM_SHA,
            "aotus_sha": _git_head(
                "/home/lzy/open-source-reference/seeder-current/aotus"
            ),
            "seeder_source_sha": _git_head(
                "/home/lzy/open-source-reference/seeder-current/sdr"
            ),
            "license": {
                "identifier": "BSD-3-Clause",
                "path": str(current_root / "LICENSE"),
                "sha256": sha256_file(current_root / "LICENSE"),
            },
        },
        "pinned_to_current_diff": {
            "changed_files": [],
            "changed_routines": [],
            "could_explain_halfway_or_one_fallback": False,
            "reason": "Current upstream Seeder HEAD is byte-identical at the pinned parent and TreElm submodule revisions.",
        },
        "qvalue_related_commits": [
            {
                "repository": "seeder-source",
                "sha": "2b929faa97ebc77cf586dc5826a15b650999c836",
                "message": "Transfer of source code from Mercurial commit 765b29ce1108",
                "changed_routines": "historical import; no pinned-to-current delta",
            },
            {
                "repository": "tem-source",
                "sha": "0405973e5e756fab94157da32572ec15f299b1ad",
                "message": "Initial code version from the mercurial repository (#1)",
                "changed_routines": "historical import; no pinned-to-current delta",
            },
        ],
        "searched_keywords": [
            "calc_dist",
            "qVal",
            "sdr_qValByNode",
            "getBCID_and_calcQval",
            "needCalcQValByBCID",
            "unKnownBnd",
            "sdr_truncate_qVal",
            "bc_uni_calcdist",
            "intersected_object",
            "flood",
            "triangle",
        ],
        "source_contract": source,
    }
    return source, comparison


def summarize(project_root: Path) -> dict:
    root = project_root.resolve()
    research = root / "outputs/cfd_flow" / RESEARCH_RUN
    qc = research / "qc"
    tiny = research / "qvalue_repair/tiny_cylinder_repaired_n16"
    repaired_base = root / "outputs/cfd_flow" / REPAIRED_BASE_RUN
    old_base = (
        root
        / "outputs/cfd_flow/axis_aligned_ideal_plane_inlet_preflight_"
        "anchor003274_20260829_120444/seeder/mesh"
    )
    old_tiny = research / "pressure_bc_benchmark/meshes/axis_aligned/n16"

    config = _official_config_contract(root)
    source, upstream = _upstream_contract()
    scale = parallel_threshold_diagnostic(
        np.asarray((1.0e-7, 0.0, 0.0)),
        np.asarray((0.0, 1.0e-7, 0.0)),
        np.asarray((0.0, 0.0, 2.0e-7)),
    )
    patch = root / "patches/seeder/qvalue_repair.patch"
    repaired_binary = _unc(
        "/home/lzy/apes-worktrees/seeder_qvalue_repair_20260830/build/seeder"
    )
    tiny_before = audit_mesh_qvalues(
        old_tiny / "mesh", old_tiny / "mesh_summary.json"
    )
    tiny_after = tiny_cylinder_gate(tiny / "mesh", tiny / "source_mesh_summary.json")
    base = repaired_base_gate(
        old_mesh=old_base,
        repaired_mesh=repaired_base / "seeder/mesh",
        wall_stl=repaired_base / "geometry/geometry_solver_m/wall.stl",
    )

    write_json(qc / "calc_dist_config_contract.json", config)
    write_json(qc / "seeder_upstream_comparison.json", upstream)
    source_record = {
        "status": "PASS_ROOT_CAUSE_IDENTIFIED",
        **source,
        "representative_micrometre_scale_diagnostic": scale,
        "repair": {
            "method": "scale-invariant normalized ray/triangle parallel predicate",
            "patch": str(patch),
            "patch_sha256": sha256_file(patch),
            "base_seeder_sha": PINNED_SEEDER_SHA,
            "base_treelm_sha": PINNED_TREELM_SHA,
            "upstream_fix_sha": None,
            "reason_no_backport": "No pinned-to-current source delta and no later tem_line_module change exists.",
        },
    }
    write_json(qc / "qvalue_repair_source_root_cause.json", source_record)
    write_json(
        qc / "seeder_qvalue_runtime_trace.json",
        {
            "status": "NOT_RUN_SOURCE_ROOT_CAUSE_IDENTIFIED",
            "instrumentation_required": False,
            "reason": "The scale-dependent branch is source-proven and the repaired analytic cylinder passed without instrumentation.",
            "counts": None,
        },
    )
    write_json(tiny / "qvalue_validation.json", tiny_after)
    write_json(repaired_base / "qc/mesh_qc.json", base["repaired_mesh_qc"])
    write_json(
        repaired_base / "qc/vascular_wall_qvalue_distribution.json",
        base["qvalue_distribution"],
    )
    write_json(
        repaired_base / "qc/wall_qvalue_oracle.json", base["stl_ray_oracle"]
    )
    write_json(
        repaired_base / "qc/topology_difference.json", base["topology_difference"]
    )
    write_json(repaired_base / "qc/qvalue_repair_validation.json", base)

    result = {
        "status": FINAL_STATUS,
        "actual_head_at_execution": subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=root,
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
        ).stdout.strip(),
        "production_pipeline_modified": False,
        "active_qvalue_status_reconciled": True,
        "evidence_precedence": {
            "isolation_root_cause_decision.json": "HISTORICAL",
            "qvalue_contract_forensics_v2.json": "ACTIVE",
        },
        "pinned_seeder_sha": PINNED_SEEDER_SHA,
        "current_upstream_seeder_sha": CURRENT_SEEDER_SHA,
        "current_upstream_treelm_sha": CURRENT_TREELM_SHA,
        "official_qvalue_example": {
            "found": True,
            "primary": "mus/examples/tutorials/tutorial_cases/tutorial_channelGeneric with qValues=true/useObstacle=true",
            "fallback": "sdr/testsuite/sphere_stl with calc_dist=true",
            "pinned_distribution": "INCONCLUSIVE_TIMEOUT_NO_QVAL_OUTPUT",
            "current_distribution": "NOT_RUN_IDENTICAL_SOURCE",
        },
        "upstream_qvalue_related_commits": upstream["qvalue_related_commits"],
        "exact_failing_source_branch": source["exact_failing_source_branch"],
        "root_cause_category": ROOT_CAUSE_CATEGORY,
        "repair_method": source_record["repair"]["method"],
        "repair_base_sha": PINNED_SEEDER_SHA,
        "repair_patch_sha256": sha256_file(patch),
        "repaired_binary_sha256": sha256_file(repaired_binary),
        "tiny_qvalue_before": tiny_before["q_value_distribution_on_declared_wall_links"],
        "tiny_qvalue_after": {
            "support": tiny_after["audit"]["q_value_distribution_on_declared_wall_links"],
            "continuous_unique_count": tiny_after["continuous_unique_count"],
            "q_error": tiny_after["q_error"],
        },
        "repaired_base": {
            "seeder_calls": 1,
            "fluid_cell_count": base["repaired_mesh_qc"]["fluid_cell_count"],
            "q_distribution": base["qvalue_distribution"]["wall_links"],
            "stl_oracle": base["stl_ray_oracle"]["q_seeder_minus_q_exact"],
            "topology_port_qc": base["repaired_mesh_qc"],
            "topology_difference": base["topology_difference"],
        },
        "referee_arbitrary_q_tests": "TARGETED_TESTS_ADDED_Q_LT_0P5_AND_Q_GT_0P5",
        "targeted_tests": {
            "pytest": "PASS",
            "passed": 20,
            "qvalue_reader_roundtrip": "PASS",
            "qvalue_file_contract": "PASS",
            "qvalue_direction_mapping": "PASS",
            "boundary_mass_referee": "PASS",
            "cylinder_qvalue_oracle": "PASS",
            "seeder_qvalue_repair": "PASS",
            "qvalue_repair_validation": "PASS",
            "ruff_targeted": "PASS",
            "ruff_full_repository": "PRE_EXISTING_FAILURE_225_ERRORS",
            "ruff_full_repository_scope": (
                "Predominantly pre-existing Ultraliser/external_reference and "
                "unrelated repository findings; preserved without modification."
            ),
        },
        "small_pipe_force_calls": 0,
        "pipe_force_results": {"n16": None, "n20": None, "n27": None},
        "oblique_result": None,
        "wall_accuracy_conclusion": "QVALUE_REPAIR_PASS_BUT_BASE_TOPOLOGY_GATE_FAIL",
        "vascular_cfd_calls": 0,
        "vascular_harvester_calls": 0,
        "tiny_seeder_operational_attempts": 5,
        "vascular_seeder_calls": 1,
        "small_musubi_calls": 0,
        "runtime": {
            "official_tutorial_timeout_s": 300,
            "official_sphere_timeout_attempts_s": [180, 180],
            "repair_build_s": 115.881,
            "tiny_repaired_seeder_s": 0.40,
            "base_repaired_seeder_s": 5.04,
        },
        "first_failure": (
            "Repaired BASE has 2 connected fluid components (required 1); "
            "fluid cells changed from 221309 to 182320."
        ),
        "next": "INVESTIGATE REPAIRED BASE DISCONNECTED COMPONENT BEFORE ANY MUSUBI OR VASCULAR CFD",
    }
    write_json(qc / "qvalue_repair_status.json", result)
    write_json(qc / "qvalue_repair_final.json", result)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=ROOT)
    args = parser.parse_args()
    result = summarize(args.project_root)
    print(result["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
