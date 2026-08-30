from __future__ import annotations

from pathlib import Path

from utils.cfd_flow.dimensionless_geometry_kernel import semantic_files_success


def test_semantic_success_requires_every_nonempty_artifact(tmp_path: Path) -> None:
    (tmp_path / "header.lua").write_text("header", encoding="utf-8")
    (tmp_path / "mesh.lsb").write_bytes(b"")

    failed = semantic_files_success(tmp_path, ("header.lua", "mesh.lsb"))
    (tmp_path / "mesh.lsb").write_bytes(b"mesh")
    passed = semantic_files_success(tmp_path, ("header.lua", "mesh.lsb"))

    assert failed["semantic_success"] is False
    assert passed["semantic_success"] is True
