from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from utils.sampling.clustering import deterministic_kmeans
from utils.sampling.feature_scaling import robust_scale
from utils.sampling.radius_features import compute_radius_features
from utils.sampling.representative_selection import select_representatives
from utils.sampling.roi_extraction import (
    EdgeSpatialIndex,
    extract_candidate_rois,
    extract_connected_roi,
    generate_anchor_ids,
    global_model_from_swc,
)
from utils.sampling.roi_features import populate_roi_features
from utils.sampling.sampling_config import SamplingConfig
from utils.sampling.sampling_types import ROIRecord
from utils.sampling.structural_features import compute_structural_features


def _model(
    points: list[tuple[float, float, float]],
    parents: list[int],
    *,
    radii: list[float] | None = None,
    bounds: tuple[float, float, float, float, float, float] = (-5, 15, -5, 5, -5, 5),
):
    node_ids = np.arange(1, len(points) + 1, dtype=np.int64)
    swc = SimpleNamespace(
        node_ids=node_ids,
        points_um=np.asarray(points, dtype=float),
        radius_raw_um=np.asarray(radii or [2.0] * len(points), dtype=float),
        parent_ids=np.asarray(parents, dtype=np.int64),
    )
    return global_model_from_swc(
        swc,
        source_model_id="synthetic",
        source_mouse_id="mouse",
        model_bounds_xyz_um=bounds,
    )


def _extract(model, anchor_id: int, size=(10.0, 10.0, 10.0)) -> ROIRecord:
    roi, reason = extract_connected_roi(
        model,
        EdgeSpatialIndex(model),
        anchor_id=anchor_id,
        roi_size_um=size,
    )
    assert reason is None
    assert roi is not None
    populate_roi_features(roi)
    return roi


def test_straight_vessel_clipping_radius_and_cut_ports() -> None:
    model = _model([(-5, 0, 0), (5, 0, 0), (15, 0, 0)], [-1, 1, 2])
    roi = _extract(model, 2)
    assert roi.node_count == 3
    assert roi.edge_count == 2
    assert roi.branch_count == 1
    assert roi.true_terminal_count == 0
    assert roi.cut_port_count == 2
    assert {port.boundary_face for port in roi.cut_ports} == {"x_min", "x_max"}
    assert sorted(port.intersection_position_um[0] for port in roi.cut_ports) == [0.0, 10.0]
    assert roi.radius_features["r50"] == 2.0


def test_y_bifurcation_branch_count_and_true_terminals() -> None:
    model = _model(
        [(-5, 0, 0), (0, 0, 0), (5, 3, 0), (5, -3, 0)],
        [-1, 1, 2, 2],
        bounds=(-5, 5, -5, 5, -5, 5),
    )
    roi = _extract(model, 2)
    assert roi.branch_count == 3
    assert roi.bifurcation_count == 1
    assert roi.true_terminal_count == 3
    assert roi.cut_port_count == 0


def test_true_terminal_is_not_cut_port_after_boundary_cutting() -> None:
    model = _model(
        [(-5, 0, 0), (0, 0, 0), (5, 3, 0), (5, -3, 0)],
        [-1, 1, 2, 2],
        bounds=(-5, 5, -5, 5, -5, 5),
    )
    roi = _extract(model, 2, size=(4.0, 4.0, 4.0))
    assert roi.true_terminal_count == 0
    assert roi.cut_port_count == 3
    assert not set(roi.true_terminal_local_ids).intersection(
        port.local_node_id for port in roi.cut_ports
    )


def test_arc_length_weighting_is_resampling_invariant_for_same_vessel() -> None:
    dense = _model(
        [(float(x), 0, 0) for x in range(11)],
        [-1] + list(range(1, 11)),
        bounds=(0, 10, -5, 5, -5, 5),
    )
    sparse = _model(
        [(0, 0, 0), (5, 0, 0), (10, 0, 0)],
        [-1, 1, 2],
        bounds=(0, 10, -5, 5, -5, 5),
    )
    dense_features = compute_radius_features(_extract(dense, 6))
    sparse_features = compute_radius_features(_extract(sparse, 2))
    for name in ("r10", "r25", "r50", "r75", "r90", "radius_mean_um"):
        assert np.isclose(dense_features[name], sparse_features[name])


def _graph_roi(edges: list[tuple[int, int]], points: list[tuple[float, float, float]]) -> ROIRecord:
    positions = np.asarray(points, dtype=float)
    edge_array = np.asarray(edges, dtype=np.int64)
    return ROIRecord(
        roi_id="graph",
        source_model_id="synthetic",
        source_mouse_id="mouse",
        anchor_id=1,
        anchor_position_um=tuple(points[0]),
        bbox_min_um=(-1, -1, -1),
        bbox_max_um=(2, 2, 2),
        bbox_center_um=(0.5, 0.5, 0.5),
        bbox_size_um=(3, 3, 3),
        global_node_ids=tuple(range(len(points))),
        global_edge_ids=tuple(range(len(edges))),
        local_node_ids=np.arange(len(points), dtype=np.int64),
        local_node_global_ids=np.arange(len(points), dtype=np.int64),
        local_node_positions_um=positions,
        local_node_radius_um=np.ones(len(points)),
        local_edges=edge_array,
        local_edge_ids=np.arange(len(edges), dtype=np.int64),
        local_edge_global_ids=np.arange(len(edges), dtype=np.int64),
        local_edge_points_um=np.asarray([(positions[a], positions[b]) for a, b in edges]),
        local_edge_radius_um=np.ones((len(edges), 2)),
        true_terminal_local_ids=(),
        true_terminal_global_ids=(),
        cut_ports=(),
        raw_component_count=1,
        raw_total_vessel_length_um=float(len(edges)),
        retained_component_length_um=float(len(edges)),
    )


def test_branch_count_and_cycle_rank() -> None:
    straight = _graph_roi([(0, 1), (1, 2)], [(0, 0, 0), (1, 0, 0), (2, 0, 0)])
    y_tree = _graph_roi(
        [(0, 1), (1, 2), (1, 3)],
        [(0, 0, 0), (1, 0, 0), (2, 1, 0), (2, -1, 0)],
    )
    loop = _graph_roi(
        [(0, 1), (1, 2), (2, 3), (3, 0)],
        [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)],
    )
    straight_features = compute_structural_features(straight)
    y_features = compute_structural_features(y_tree)
    loop_features = compute_structural_features(loop)
    assert straight_features["branch_count"] == 1
    assert y_features["branch_count"] == 3
    assert y_features["bifurcation_count"] == 1
    assert straight_features["cycle_rank"] == 0
    assert loop_features["cycle_rank"] == 1


def test_sampling_is_deterministic_and_representatives_are_real(tmp_path) -> None:
    config = SamplingConfig(
        output_root=tmp_path,
        seed=17,
        roi_size_um=(10, 10, 10),
        min_anchor_distance_um=2,
        max_candidate_anchors=5,
        n_clusters=2,
        target_selected_count=2,
        max_selected_overlap=1.0,
        min_representative_distance_um=0.0,
    )
    model = _model(
        [(float(x), float(x % 3), 0) for x in range(12)],
        [-1] + list(range(1, 12)),
        bounds=(0, 11, -5, 5, -5, 5),
    )
    assert generate_anchor_ids(model, config) == generate_anchor_ids(model, config)
    first_batch = extract_candidate_rois(model, config)
    second_batch = extract_candidate_rois(model, config)
    assert first_batch.anchor_ids == second_batch.anchor_ids
    assert [roi.roi_id for roi in first_batch.candidates] == [
        roi.roi_id for roi in second_batch.candidates
    ]
    for first_roi, second_roi in zip(first_batch.candidates, second_batch.candidates):
        populate_roi_features(first_roi)
        populate_roi_features(second_roi)
        assert first_roi.radius_features == second_roi.radius_features
        assert first_roi.structural_features == second_roi.structural_features
    values = np.asarray(((0, 0), (0.1, 0.2), (4, 4), (4.2, 3.9)), dtype=float)
    scaled, _ = robust_scale(values, ("r10", "branch_count"))
    first = deterministic_kmeans(scaled, n_clusters=2, feature_names=("r10", "branch_count"), seed=17, max_iter=50)
    second = deterministic_kmeans(scaled, n_clusters=2, feature_names=("r10", "branch_count"), seed=17, max_iter=50)
    assert np.array_equal(first.assignments, second.assignments)
    rois = [_graph_roi([(0, 1)], [(index * 20.0, 0, 0), (index * 20.0 + 1, 0, 0)]) for index in range(4)]
    for index, roi in enumerate(rois):
        roi.roi_id = f"roi_{index}"
        roi.bbox_min_um = (index * 20.0, 0, 0)
        roi.bbox_max_um = (index * 20.0 + 5, 5, 5)
        roi.bbox_center_um = (index * 20.0 + 2.5, 2.5, 2.5)
    selection = select_representatives(rois, first, config)
    assert selection.selected_indices == select_representatives(rois, second, config).selected_indices
    assert all(rois[index].roi_id in {roi.roi_id for roi in rois} for index in selection.selected_indices)
