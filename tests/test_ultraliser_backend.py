from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pytest
import trimesh

from utils.cfd_lumen.config import CFDLumenConfig
from utils.cfd_lumen.ultraliser_backend import (
    ULTRALISER_PACKING_ALGORITHM,
    ULTRALISER_RADIUS_SCALE,
    UltraliserCycleConflict,
    build_ultraliser_command,
    export_roi_for_ultraliser,
)
from utils.cfd_lumen.ultraliser_qc import (
    evaluate_radius_fidelity,
    evaluate_surface_topology,
    validate_source_roi,
    write_geometry_outputs,
)
from utils.sampling.sampling_types import CutPort, ROIRecord


def _roi(
    *,
    positions: np.ndarray | None = None,
    radii: np.ndarray | None = None,
    edges: np.ndarray | None = None,
    anchor_id: int = 42,
    cut_nodes: tuple[int, ...] = (0,),
) -> ROIRecord:
    positions = np.asarray(
        positions
        if positions is not None
        else ((0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (8.0, 2.0, 0.0), (8.0, -2.0, 0.0)),
        dtype=float,
    )
    radii = np.asarray(radii if radii is not None else (1.5, 1.4, 1.0, 1.1), dtype=float)
    edges = np.asarray(
        edges if edges is not None else ((0, 1), (1, 2), (1, 3)), dtype=np.int64
    ).reshape((-1, 2))
    edge_points = positions[edges]
    edge_radii = radii[edges]
    cuts = tuple(
        CutPort(
            cut_port_id=f"roi__cut_{index:03d}",
            local_node_id=node,
            global_edge_id=100
            + next(edge_index for edge_index, edge in enumerate(edges) if node in edge),
            intersection_position_um=tuple(map(float, positions[node])),
            radius_at_cut_um=float(radii[node]),
            boundary_face="test",
        )
        for index, node in enumerate(cut_nodes)
    )
    return ROIRecord(
        roi_id=f"synthetic__anchor_{anchor_id:06d}",
        source_model_id="synthetic",
        source_mouse_id="mouse",
        anchor_id=anchor_id,
        anchor_position_um=tuple(map(float, positions[0])),
        bbox_min_um=tuple(map(float, positions.min(axis=0))),
        bbox_max_um=tuple(map(float, positions.max(axis=0))),
        bbox_center_um=tuple(map(float, 0.5 * (positions.min(axis=0) + positions.max(axis=0)))),
        bbox_size_um=tuple(map(float, np.ptp(positions, axis=0))),
        global_node_ids=tuple(range(len(positions))),
        global_edge_ids=tuple(range(len(edges))),
        local_node_ids=np.arange(len(positions), dtype=np.int64),
        local_node_global_ids=np.arange(len(positions), dtype=np.int64),
        local_node_positions_um=positions,
        local_node_radius_um=radii,
        local_edges=edges,
        local_edge_ids=np.arange(len(edges), dtype=np.int64),
        local_edge_global_ids=np.arange(len(edges), dtype=np.int64),
        local_edge_points_um=edge_points,
        local_edge_radius_um=edge_radii,
        true_terminal_local_ids=(),
        true_terminal_global_ids=(),
        cut_ports=cuts,
        raw_component_count=1,
        raw_total_vessel_length_um=float(
            np.linalg.norm(np.diff(edge_points, axis=1), axis=2).sum()
        ),
        retained_component_length_um=float(
            np.linalg.norm(np.diff(edge_points, axis=1), axis=2).sum()
        ),
    )


def test_default_radius_scale_091_writes_diameter_1p82_and_recovers_0p91(
    tmp_path: Path,
) -> None:
    result = export_roi_for_ultraliser(_roi(radii=np.ones(4)), tmp_path)
    with h5py.File(result.h5_path, "r") as stream:
        diameter = np.asarray(stream["points"])[:, 3]
        assert np.allclose(diameter, 1.82)
        assert np.allclose(0.5 * diameter, 0.91)
        assert float(stream.attrs["radius_scale_for_ultraliser"]) == pytest.approx(0.91)
    mapping = result.radius_feed_mapping_path.read_text(encoding="utf-8-sig")
    assert "source_radius_um" in mapping and "feed_radius_um" in mapping


def test_radius_scale_1_preserves_unscaled_h5_reader_radius(tmp_path: Path) -> None:
    roi = _roi()
    result = export_roi_for_ultraliser(
        roi,
        tmp_path,
        radius_scale_for_ultraliser=1.0,
    )
    with h5py.File(result.h5_path, "r") as stream:
        recovered = 0.5 * np.asarray(stream["points"])[:, 3]
    expected = np.asarray(
        [roi.local_node_radius_um[node] for section in result.sections for node in section]
    )
    assert np.allclose(recovered, expected)


def test_source_swc_hash_is_independent_of_h5_radius_scale(tmp_path: Path) -> None:
    roi = _roi()
    default_directory = tmp_path / "default"
    unscaled_directory = tmp_path / "unscaled"
    default_directory.mkdir()
    unscaled_directory.mkdir()
    scaled = export_roi_for_ultraliser(roi, default_directory)
    unscaled = export_roi_for_ultraliser(
        roi,
        unscaled_directory,
        radius_scale_for_ultraliser=1.0,
    )
    scaled_hash = hashlib.sha256(scaled.swc_path.read_bytes()).hexdigest()
    unscaled_hash = hashlib.sha256(unscaled.swc_path.read_bytes()).hexdigest()
    assert scaled_hash == unscaled_hash
    metadata = json.loads(scaled.metadata_path.read_text(encoding="utf-8"))
    assert metadata["source_swc_modified"] is False
    assert metadata["source_radius_modified"] is False


def test_roi_export_is_not_bound_to_anchor_3274(tmp_path: Path) -> None:
    result = export_roi_for_ultraliser(_roi(anchor_id=8123), tmp_path)
    assert result.swc_path.is_file()
    assert result.h5_path.is_file()
    assert result.radius_feed_mapping_path.is_file()
    assert result.cut_port_mapping_path.is_file()


def test_cycle_is_rejected_without_edge_deletion(tmp_path: Path) -> None:
    roi = _roi(
        positions=np.asarray(((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 2.0, 0.0))),
        radii=np.ones(3),
        edges=np.asarray(((0, 1), (1, 2), (2, 0))),
    )
    with pytest.raises(UltraliserCycleConflict, match="cycle_rank=1"):
        export_roi_for_ultraliser(roi, tmp_path)


def test_command_contains_the_validated_official_workflow_flags(tmp_path: Path) -> None:
    command, rendered = build_ultraliser_command(
        tmp_path / "ultraVessMorpho2Mesh.exe",
        tmp_path / "roi_core.h5",
        tmp_path / "raw",
        prefix="anchor003274",
    )
    assert command[command.index("--voxels-per-micron") + 1] == "6"
    assert command[command.index("--packing-algorithm") + 1] == ULTRALISER_PACKING_ALGORITHM
    assert command[command.index("--voxelization-axis") + 1] == "xyz"
    assert command[command.index("--isosurface-technique") + 1] == "dmc"
    assert command[command.index("--optimization-iterations") + 1] == "5"
    assert command[command.index("--smooth-iterations") + 1] == "5"
    assert command[command.index("--laplacian-iterations") + 1] == "10"
    assert command[command.index("--threads") + 1] == "8"
    assert "--solid" in command and "--adaptive-optimization" in command
    assert "--export-stl-mesh" in command
    assert all(term not in rendered.lower() for term in ("polyball", "boolean", "fallback"))


def test_surface_topology_qc_accepts_clean_watertight_mesh() -> None:
    mesh = trimesh.creation.icosphere(subdivisions=2, radius=1.0)
    report = evaluate_surface_topology(mesh, CFDLumenConfig())
    assert report["status"] == "PASS"
    assert report["watertight"] is True
    assert report["component_count"] == 1
    assert report["boundary_edge_count"] == 0
    assert report["nonmanifold_edge_count"] == 0
    assert report["self_intersection_count"] == 0
    assert report["degenerate_triangle_count"] == 0


def test_radius_fidelity_qc_uses_source_branch_cross_section() -> None:
    positions = np.asarray(((0.0, 0.0, -6.0), (0.0, 0.0, 0.0), (0.0, 0.0, 6.0)))
    roi = _roi(
        positions=positions,
        radii=np.ones(3),
        edges=np.asarray(((0, 1), (1, 2))),
        cut_nodes=(0, 2),
    )
    mesh = trimesh.creation.cylinder(radius=1.0, height=12.0, sections=96)
    branches, source_report = validate_source_roi(roi)
    samples, report = evaluate_radius_fidelity(mesh, branches, roi, CFDLumenConfig())
    assert source_report["source_radius_modified"] is False
    assert len(samples) == 1
    assert report["status"] == "PASS"
    assert report["p95_absolute_relative_error"] < 0.01


def test_stl_vtp_and_metre_stl_exports_are_unit_exact(tmp_path: Path) -> None:
    mesh = trimesh.creation.icosphere(subdivisions=1, radius=3.0)
    raw = tmp_path / "raw.stl"
    mesh.export(raw)
    geometry = tmp_path / "geometry"
    written = write_geometry_outputs(raw, geometry)
    metre = trimesh.load_mesh(geometry / "lumen_surface_m.stl", process=False)
    assert (geometry / "lumen_surface_um.stl").is_file()
    assert (geometry / "lumen_surface_um.vtp").is_file()
    assert (geometry / "lumen_surface_m.stl").is_file()
    assert np.ptp(np.asarray(metre.vertices), axis=0) == pytest.approx(
        np.ptp(np.asarray(written.vertices), axis=0) * 1.0e-6,
        rel=1.0e-6,
    )


def test_formal_configuration_defaults_are_validated() -> None:
    config = CFDLumenConfig()
    config.validate()
    assert config.ultraliser.radius_scale == ULTRALISER_RADIUS_SCALE
    assert config.ultraliser.voxels_per_micron == 6.0
    assert config.ultraliser.threads == 8
