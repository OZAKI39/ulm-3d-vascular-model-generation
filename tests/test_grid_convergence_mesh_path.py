from __future__ import annotations

from pathlib import Path

import pytest
import numpy as np

from utils.cfd_flow.exact_link_flux import reconstruct_boundary
from utils.cfd_flow.grid_convergence import (
    MESH_QVAL_FILES,
    MESH_REQUIRED_FILES,
    can_reuse_completed_seeder_output,
    generate_grid_musubi_lua,
    grid_lua_contract,
    grid_specs,
    mesh_file_contract,
    musubi_mesh_assignment,
    require_complete_mesh,
    resolve_mesh_directory,
    three_grid_scalar_analysis,
)
from utils.cfd_flow.io import FlowError


def _complete_fake_mesh(mesh_dir: Path) -> None:
    mesh_dir.mkdir(parents=True)
    for name in MESH_REQUIRED_FILES + MESH_QVAL_FILES:
        (mesh_dir / name).write_bytes(b"mesh-contract-test\n")


def test_seeder_output_resolution_adds_mesh_exactly_once(tmp_path: Path) -> None:
    run_root = tmp_path / "coarse" / "seeder"
    run_root.mkdir(parents=True)
    first = resolve_mesh_directory(run_root)
    second = resolve_mesh_directory(first.mesh_dir)
    assert first.mesh_dir == run_root / "mesh"
    assert second.mesh_dir == first.mesh_dir
    assert second.mesh_dir.name == "mesh"
    assert second.mesh_dir.parent.name == "seeder"


def test_windows_to_wsl_conversion_targets_the_same_mesh_dir(tmp_path: Path) -> None:
    run_root = tmp_path / "coarse" / "seeder"
    run_root.mkdir(parents=True)
    contract = resolve_mesh_directory(run_root)
    assert contract.mesh_windows == run_root / "mesh"
    assert contract.mesh_wsl.endswith("/coarse/seeder/mesh")
    assert "\\" not in contract.mesh_wsl


def test_missing_bnd_lsb_hard_blocks_musubi(tmp_path: Path) -> None:
    mesh_dir = tmp_path / "seeder" / "mesh"
    _complete_fake_mesh(mesh_dir)
    (mesh_dir / "bnd.lsb").unlink()
    with pytest.raises(FlowError, match="bnd.lsb"):
        require_complete_mesh(mesh_dir)
    with pytest.raises(FlowError, match="bnd.lsb"):
        musubi_mesh_assignment(mesh_dir)


def test_complete_mesh_musubi_assignment_uses_real_mesh_dir(tmp_path: Path) -> None:
    mesh_dir = tmp_path / "coarse" / "seeder" / "mesh"
    _complete_fake_mesh(mesh_dir)
    assignment = musubi_mesh_assignment(mesh_dir)
    assert assignment.startswith("mesh = '")
    assert assignment.endswith("/coarse/seeder/mesh/'")
    assert "/mesh/mesh/" not in assignment


def test_completed_output_is_reusable_but_failed_or_incomplete_is_not(
    tmp_path: Path,
) -> None:
    mesh_dir = tmp_path / "seeder" / "mesh"
    _complete_fake_mesh(mesh_dir)
    passed = {"status": "PASS", "seeder_returncode": 0}
    failed = {"status": "FAIL", "seeder_returncode": 2}
    assert mesh_file_contract(mesh_dir)["status"] == "PASS"
    assert can_reuse_completed_seeder_output(mesh_dir, passed)
    assert not can_reuse_completed_seeder_output(mesh_dir, failed)
    (mesh_dir / "elemlist.lsb").unlink()
    assert not can_reuse_completed_seeder_output(mesh_dir, passed)


def test_grid_qc_can_record_zero_aggregate_normal_without_inventing_one() -> None:
    boundary_ids = np.zeros((1, 26), dtype=np.int64)
    boundary_ids[0, 0] = 1
    boundary_ids[0, 3] = 1
    with pytest.raises(ValueError, match="cannot be zero"):
        reconstruct_boundary(
            boundary_ids,
            np.asarray((0,)),
            label="inlet",
            boundary_id=1,
        )
    tolerant = reconstruct_boundary(
        boundary_ids,
        np.asarray((0,)),
        label="inlet",
        boundary_id=1,
        allow_zero_normals=True,
    )
    assert tolerant.normal_indices.tolist() == [-1]


def test_grid_musubi_lua_freezes_physics_and_uses_canonical_mesh(
    tmp_path: Path,
) -> None:
    mesh_dir = tmp_path / "coarse" / "seeder" / "mesh"
    _complete_fake_mesh(mesh_dir)
    spec = grid_specs()["coarse"]
    lua = generate_grid_musubi_lua(
        spec,
        mesh_dir=mesh_dir,
        maximum_iterations=2_000,
        simulation_name="gc_c_test",
        write_restarts=False,
    )
    contract = grid_lua_contract(lua, spec, mesh_dir)
    assert contract["status"] == "PASS"
    assert "restart =" not in lua
    assert "read =" not in lua
    steady = generate_grid_musubi_lua(
        spec,
        mesh_dir=mesh_dir,
        maximum_iterations=1_000_000,
        simulation_name="gc_c_steady_test",
        write_restarts=True,
    )
    assert grid_lua_contract(steady, spec, mesh_dir)["status"] == "PASS"
    assert "write = 'restart/'" in steady
    assert "interval = { iter = 5000 }" in steady
    resumed = generate_grid_musubi_lua(
        spec,
        mesh_dir=mesh_dir,
        maximum_iterations=1_000_000,
        simulation_name="gc_c_resume_test",
        write_restarts=True,
        resume_header="restart/checkpoint_header.lua",
        first_restart_iteration=15_000,
    )
    assert grid_lua_contract(
        resumed,
        spec,
        mesh_dir,
        resume_header="restart/checkpoint_header.lua",
    )["status"] == "PASS"
    assert "min = { iter = 15000 }" in resumed


def test_three_grid_analysis_reports_valid_and_unavailable_cases() -> None:
    # Exact second-order sequence around an extrapolated value of 1.0.
    r = 1.3
    fine = 1.0 + 0.01
    base = 1.0 + 0.01 * r**2
    coarse = 1.0 + 0.01 * r**4
    valid = three_grid_scalar_analysis(coarse, base, fine, refinement_ratio=r)
    assert valid["status"] == "AVAILABLE"
    assert valid["observed_order_p"] == pytest.approx(2.0)
    assert valid["richardson_extrapolation"] == pytest.approx(1.0)
    oscillatory = three_grid_scalar_analysis(1.1, 0.9, 1.0, refinement_ratio=r)
    assert oscillatory["status"] == "UNAVAILABLE"
    assert oscillatory["observed_order_p"] is None
