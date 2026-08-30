from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pytest

from utils.cfd_flow.seeder_qvalue_repair import (
    CURRENT_SEEDER_SHA,
    CURRENT_TREELM_SHA,
    PINNED_SEEDER_SHA,
    PINNED_TREELM_SHA,
    ROOT_CAUSE_CATEGORY,
    TRACE_KEYS,
    calc_dist_config_contract,
    forbidden_production_paths_modified,
    parallel_threshold_diagnostic,
    parse_runtime_trace,
    source_contract,
)


def _root() -> Path:
    return Path(__file__).resolve().parents[1]


def _wsl_path(relative: str) -> Path:
    return Path(r"\\wsl.localhost\Ubuntu\home\lzy") / relative.replace("/", "\\")


def test_micrometre_geometry_exposes_dimensional_parallel_bug() -> None:
    diagnostic = parallel_threshold_diagnostic(
        np.asarray((1.0e-7, 0.0, 0.0)),
        np.asarray((0.0, 1.0e-7, 0.0)),
        np.asarray((0.0, 0.0, 2.0e-7)),
    )

    assert diagnostic["determinant_b"] == pytest.approx(2.0e-21)
    assert diagnostic["pinned_absolute_predicate_parallel"] is True
    assert diagnostic["repaired_scale_invariant_predicate_parallel"] is False
    assert np.isclose(diagnostic["normalized_abs_dot"], 1.0)


def test_calc_dist_config_matches_official_stl_obstacle_contract() -> None:
    root = _root()
    base = (
        root
        / "outputs/cfd_flow/axis_aligned_ideal_plane_inlet_preflight_"
        "anchor003274_20260829_120444/seeder/seeder.lua"
    ).read_text(encoding="utf-8")
    pipe = (
        root
        / "outputs/cfd_flow/healthy_mouse_capillary_port_grid_sensitivity_"
        "research_anchor003274_20260830/pressure_bc_benchmark/meshes/"
        "axis_aligned/n27/seeder.lua"
    ).read_text(encoding="utf-8")
    official = (
        _wsl_path(
            "apes-pinned/musubi_official/mus/examples/tutorials/"
            "tutorial_cases/tutorial_channelGeneric/seeder.lua"
        )
    ).read_text(encoding="utf-8")

    contract = calc_dist_config_contract(base, pipe, official)

    assert contract["status"] == "PASS"
    assert contract["selected"]["vascular_base"]["calc_dist"] == "true"
    assert contract["selected"]["pipe_axis_n27"]["calc_dist"] == "true"
    assert (
        contract["selected"]["official_tutorial_channelGeneric"]["calc_dist"]
        == "qValues"
    )


def test_pinned_and_current_upstream_share_exact_failing_source_branch() -> None:
    pinned_path = _wsl_path(
        "apes-pinned/seeder_official/tem/source/shapes/tem_line_module.fpp"
    )
    current_path = _wsl_path(
        "open-source-reference/seeder-current/tem/source/shapes/tem_line_module.fpp"
    )
    contract = source_contract(
        pinned_source=pinned_path.read_text(encoding="utf-8"),
        current_source=current_path.read_text(encoding="utf-8"),
        pinned_seeder_sha=PINNED_SEEDER_SHA,
        current_seeder_sha=CURRENT_SEEDER_SHA,
        pinned_treelm_sha=PINNED_TREELM_SHA,
        current_treelm_sha=CURRENT_TREELM_SHA,
    )

    assert contract["status"] == "PASS"
    assert contract["pinned_equals_current"] is True
    assert contract["upstream_fix_available"] is False
    assert contract["root_cause_category"] == ROOT_CAUSE_CATEGORY


def test_qvalue_runtime_trace_parser_accepts_aggregate_counters() -> None:
    text = "\n".join(f"{key}: {index}" for index, key in enumerate(TRACE_KEYS, 1))

    parsed = parse_runtime_trace(text)

    assert parsed["wall_bcid_count"] == 1
    assert parsed["missing_intersected_object_count"] == len(TRACE_KEYS)


def test_patch_and_worktree_do_not_modify_production_pipeline() -> None:
    root = _root()
    patch = (root / "patches/seeder/qvalue_repair.patch").read_text(encoding="utf-8")
    assert "source/shapes/tem_line_module.fpp" in patch
    assert "sdr_boundary_module" not in patch
    status = subprocess.run(
        ["git", "status", "--porcelain"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
        encoding="utf-8",
    ).stdout
    changed = [line[3:].strip().replace("\\", "/") for line in status.splitlines()]
    assert forbidden_production_paths_modified(changed) == []
