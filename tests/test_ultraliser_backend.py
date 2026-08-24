from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

from utils.cfd_lumen.ultraliser_backend import (
    ULTRALISER_PACKING_ALGORITHM,
    UltraliserCycleConflict,
    build_ultraliser_command,
    estimate_resolution,
    export_roi_for_ultraliser,
)
from utils.sampling.sampling_types import CutPort, ROIRecord


def _roi(
    *,
    positions: np.ndarray | None = None,
    radii: np.ndarray | None = None,
    edges: np.ndarray | None = None,
) -> ROIRecord:
    positions = np.asarray(
        positions
        if positions is not None
        else ((0.0, 0.0, 0.0), (4.0, 0.0, 0.0), (8.0, 2.0, 0.0), (8.0, -2.0, 0.0)),
        dtype=float,
    )
    radii = np.asarray(radii if radii is not None else (1.5, 1.4, 1.0, 1.1), dtype=float)
    edges = np.asarray(edges if edges is not None else ((0, 1), (1, 2), (1, 3)), dtype=np.int64)
    edge_points = positions[edges]
    edge_radii = radii[edges]
    cut = CutPort(
        cut_port_id="roi__cut_000",
        local_node_id=0,
        global_edge_id=100,
        intersection_position_um=tuple(map(float, positions[0])),
        radius_at_cut_um=float(radii[0]),
        boundary_face="x_min",
    )
    return ROIRecord(
        roi_id="synthetic__anchor_003274",
        source_model_id="synthetic",
        source_mouse_id="mouse",
        anchor_id=3274,
        anchor_position_um=tuple(map(float, positions[1])),
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
        true_terminal_local_ids=tuple(int(node) for node in range(len(positions)) if node != 0),
        true_terminal_global_ids=tuple(int(node) for node in range(len(positions)) if node != 0),
        cut_ports=(cut,),
        raw_component_count=1,
        raw_total_vessel_length_um=float(np.linalg.norm(np.diff(edge_points, axis=1), axis=2).sum()),
        retained_component_length_um=float(np.linalg.norm(np.diff(edge_points, axis=1), axis=2).sum()),
    )


def test_export_preserves_swc_geometry_and_h5_reader_radius(tmp_path: Path) -> None:
    roi = _roi()
    result = export_roi_for_ultraliser(roi, tmp_path)
    rows = [
        line.split()
        for line in result.swc_path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    assert len(rows) == roi.node_count
    by_swc_id = {int(row[0]): row for row in rows}
    metadata = __import__("json").loads(result.metadata_path.read_text(encoding="utf-8"))
    mapping = {int(key): int(value) for key, value in metadata["swc_node_id_by_local_node_id"].items()}
    for local_id, swc_id in mapping.items():
        row = by_swc_id[swc_id]
        assert np.allclose(np.asarray(row[2:5], dtype=float), roi.local_node_positions_um[local_id])
        assert float(row[5]) == pytest.approx(float(roi.local_node_radius_um[local_id]))

    with h5py.File(result.h5_path, "r") as stream:
        points = np.asarray(stream["points"])
        assert stream.attrs["points_fourth_column"] == "diameter_um"
        assert np.all(points[:, 3] > 0)
        recovered = 0.5 * points[:, 3]
    expanded = np.asarray(
        [roi.local_node_radius_um[node] for section in result.sections for node in section]
    )
    assert np.allclose(recovered, expanded)

    covered = {
        tuple(sorted((first, second)))
        for section in result.sections
        for first, second in zip(section[:-1], section[1:])
    }
    expected = {tuple(sorted(map(int, edge))) for edge in roi.local_edges}
    assert covered == expected


def test_cycle_is_rejected_without_edge_deletion(tmp_path: Path) -> None:
    roi = _roi(
        positions=np.asarray(((0.0, 0.0, 0.0), (2.0, 0.0, 0.0), (1.0, 2.0, 0.0))),
        radii=np.ones(3),
        edges=np.asarray(((0, 1), (1, 2), (2, 0))),
    )
    with pytest.raises(UltraliserCycleConflict, match="cycle_rank=1"):
        export_roi_for_ultraliser(roi, tmp_path)


def test_resolution_falls_back_from_12_to_10_cells() -> None:
    positions = np.asarray(
        ((0.0, 0.0, 0.0), (80.0, 0.0, 0.0), (80.0, 80.0, 0.0), (80.0, 80.0, 120.0))
    )
    roi = _roi(positions=positions, radii=np.full(4, 1.15), edges=np.asarray(((0, 1), (1, 2), (2, 3))))
    estimate = estimate_resolution(roi)
    assert estimate["requested_12_cell_case"]["predicted_voxel_count"] > 150_000_000
    assert estimate["fallback_to_10_cells_applied"] is True
    assert estimate["selected_final_case"]["predicted_voxel_count"] <= 150_000_000
    assert estimate["selected_final_case"]["voxels_per_micron"] == 5.0
    assert estimate["smoke_case"]["voxels_per_micron"] == 2.5


def test_command_uses_only_confirmed_vascular_workflow_flags(tmp_path: Path) -> None:
    executable = tmp_path / "ultraVessMorpho2Mesh.exe"
    morphology = tmp_path / "roi_core.h5"
    output = tmp_path / "raw_ultraliser"
    command, rendered = build_ultraliser_command(
        executable,
        morphology,
        output,
        prefix="roi003274_final",
        voxels_per_micron=5.0,
        threads=8,
    )
    assert command[0] == str(executable)
    for flag in (
        "--morphology",
        "--output-directory",
        "--prefix",
        "--scaled-resolution",
        "--voxels-per-micron",
        "--solid",
        "--voxelization-axis",
        "--packing-algorithm",
        "--export-stl-mesh",
    ):
        assert flag in command
    assert command[command.index("--voxelization-axis") + 1] == "xyz"
    assert command[command.index("--packing-algorithm") + 1] == ULTRALISER_PACKING_ALGORITHM
    assert "v7" not in rendered and "v8" not in rendered and "v9" not in rendered
