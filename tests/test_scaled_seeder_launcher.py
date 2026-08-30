from __future__ import annotations

from pathlib import Path

from utils.cfd_flow.repaired_topology_forensics import (
    REQUIRED_SEEDER_MESH_FILES,
    evaluate_scaled_seeder_preflight,
    scaled_seeder_command,
    scaled_seeder_preflight_command,
    seeder_mesh_semantic_success,
)


def _preflight_stdout(workdir: str) -> str:
    return "\n".join(
        (
            f"PREFLIGHT_CWD={workdir}",
            f"SEEDER_LUA_REALPATH={workdir}/seeder.lua",
            "SEEDER_LUA_SHA256=lua-sha",
            "WALL_STL_SHA256=wall-sha",
            "SEEDER_BINARY_SHA256=binary-sha",
            "SEEDER_SOURCE_HEAD=source-sha",
        )
    )


def test_launcher_commands_set_exact_working_directory_and_visible_lua() -> None:
    workdir = "/mnt/e/research/unpatched scaled/seeder"
    preflight = scaled_seeder_preflight_command(workdir, binary_wsl="/bin/seeder")
    launch = scaled_seeder_command(workdir, binary_wsl="/bin/seeder")

    assert "Start-Process" not in preflight + launch
    assert "\n" not in preflight + launch
    assert "cd '/mnt/e/research/unpatched scaled/seeder'" in preflight
    assert "test -f seeder.lua" in preflight
    assert "exec /bin/seeder seeder.lua" in launch


def test_preflight_requires_exact_working_directory_and_hashes() -> None:
    workdir = "/mnt/e/research/scaled/seeder"
    result = evaluate_scaled_seeder_preflight(
        returncode=0,
        stdout=_preflight_stdout(workdir),
        expected_workdir_wsl=workdir,
        expected_lua_sha256="lua-sha",
        expected_wall_sha256="wall-sha",
        expected_binary_sha256="binary-sha",
        expected_source_head="source-sha",
    )
    wrong_directory = evaluate_scaled_seeder_preflight(
        returncode=0,
        stdout=_preflight_stdout("/wrong/directory"),
        expected_workdir_wsl=workdir,
        expected_lua_sha256="lua-sha",
        expected_wall_sha256="wall-sha",
        expected_binary_sha256="binary-sha",
        expected_source_head="source-sha",
    )

    assert result["status"] == "PASS"
    assert result["checks"]["seeder_lua_visible"] is True
    assert wrong_directory["status"] == "FAIL"


def test_semantic_success_requires_every_nonempty_mesh_file(tmp_path: Path) -> None:
    mesh = tmp_path / "mesh"
    mesh.mkdir()
    for name in REQUIRED_SEEDER_MESH_FILES:
        (mesh / name).write_bytes(b"evidence")

    result = seeder_mesh_semantic_success(mesh, returncode=0, stdout="Seeder done")

    assert result["status"] == "PASS"
    assert result["semantic_success"] is True


def test_returncode_zero_without_mesh_is_semantic_failure(tmp_path: Path) -> None:
    result = seeder_mesh_semantic_success(
        tmp_path / "missing-mesh",
        returncode=0,
        stdout="Cannot load configuration file: cannot open seeder.lua",
    )

    assert result["status"] == "FAIL"
    assert result["semantic_success"] is False
    assert result["checks"]["returncode_zero"] is True
    assert result["checks"]["configuration_loaded"] is False
    assert result["checks"]["all_required_mesh_files_nonempty"] is False
