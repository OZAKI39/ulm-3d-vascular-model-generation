from __future__ import annotations

import struct
from pathlib import Path

import numpy as np

from utils.cfd_flow.repaired_topology_forensics import (
    scale_binary_stl,
    scale_seeder_lua_geometry,
)


def test_lua_geometry_scaling_preserves_level_and_scales_coordinates() -> None:
    source = """minlevel = 9
bounding_cube = { origin = { 1e-6, 2e-6, 3e-6 }, length = 4e-6 }
spatial_object = {{
 attribute = { kind = 'boundary', label = 'outlet_01', level = minlevel },
 geometry = { kind = 'canoND', object = { origin = { 5e-6, 6e-6, 7e-6 }, vec = { { 8e-6, 0, 0 }, { 0, 9e-6, 0 } } } }
}}
"""
    scaled = scale_seeder_lua_geometry(source)
    assert "minlevel = 9" in scaled
    assert "outlet_01" in scaled
    assert "origin = { 1, 2, 3 }" in scaled
    assert "length = 4" in scaled
    assert "origin = { 5, 6, 7 }" in scaled
    assert "vec = { { 8, 0, 0 }, { 0, 9, 0 } }" in scaled


def test_binary_stl_scaling_changes_vertices_not_normals(
    tmp_path: Path,
) -> None:
    source = tmp_path / "source.stl"
    destination = tmp_path / "scaled.stl"
    header = b"scale-test".ljust(80, b"\0") + struct.pack("<I", 1)
    facet = struct.pack(
        "<12fH",
        0.0,
        0.0,
        1.0,
        1.0e-6,
        0.0,
        0.0,
        0.0,
        2.0e-6,
        0.0,
        0.0,
        0.0,
        3.0e-6,
        0,
    )
    source.write_bytes(header + facet)
    scale_binary_stl(source, destination)
    values = struct.unpack_from("<12fH", destination.read_bytes(), 84)
    assert np.allclose(values[:3], (0.0, 0.0, 1.0))
    assert np.allclose(values[3:12], (1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0))


def test_ascii_stl_scaling_changes_only_vertex_coordinates(tmp_path: Path) -> None:
    source = tmp_path / "source_ascii.stl"
    destination = tmp_path / "scaled_ascii.stl"
    source.write_text(
        """solid test
  facet normal 0 0 1
    outer loop
      vertex 1e-6 0 0
      vertex 0 2e-6 0
      vertex 0 0 3e-6
    endloop
  endfacet
endsolid test
""",
        encoding="utf-8",
    )

    result = scale_binary_stl(source, destination)
    scaled = destination.read_text(encoding="utf-8")

    assert result["format"] == "ASCII_STL"
    assert result["triangle_count"] == 1
    assert "facet normal 0 0 1" in scaled
    assert "vertex 1 0 0" in scaled
    assert "vertex 0 2 0" in scaled
    assert "vertex 0 0 3" in scaled
