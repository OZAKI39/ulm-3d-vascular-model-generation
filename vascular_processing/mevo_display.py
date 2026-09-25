"""Adapt exact graph ROIs to the existing display data contract; no renderer edits."""
from dataclasses import dataclass
import json
from pathlib import Path

import numpy as np

from utils.rodent_vasculature.catalog import build_catalog, select_records
from utils.sampling.pipeline import load_models_from_rodent_run
from utils.sampling.roi_features import populate_roi_features
from utils.sampling.sampling_io import (create_sampling_layout, write_candidate_tables, write_global_edge_manifest,
                                       write_json, write_roi_library, write_sampling_config)
from utils.sampling.sampling_types import CutPort, ROIRecord

from .mevo_graph import graph_counts, load_exact_graph
from .mevo_pipeline import extract_file
from .mevo_roi import extract_graph
from .roi_landmarks import LandmarkError, ROI_DEFINITIONS, load_landmarks


@dataclass
class MeVODisplayRun:
    run_root: Path
    status: str
    candidates: list
    summary_path: Path
    log_path: Path


def preflight_mevo(config, landmarks, *, allow_natural_terminal=False):
    records = select_records(build_catalog(config.input_dir, config.cohort), sample_id=config.sample_id,
        parent_group_id=config.parent_group_id, split=config.split, max_samples=config.max_samples)
    if len(records) != 1:
        raise LandmarkError("Existing dual-viewport MeVO display requires one explicitly selected source; use CLI batch for multiple files")
    source = records[0].swc_path
    document, _ = load_landmarks(Path(landmarks), source)
    if document["source"].get("units", "mm") != config.swc_units:
        raise LandmarkError("Landmark units disagree with the visualization source units")
    enabled = [roi for roi in document["rois"].values() if roi["enabled"]]
    if not enabled:
        raise LandmarkError(f"尚未标注 MeVO 解剖边界：{landmarks}。请运行 tools/vascularmd_extract_mevo.py annotate 标注近端及所有远端；不会用空间采样 ROI 替代 MeVO。")
    exact = load_exact_graph(source)
    for roi in enabled:
        extract_graph(exact, roi, allow_natural_terminal=allow_natural_terminal)


def display_record(model, exact, definition, result, name, *, anatomical_status, units, rank):
    ids = [int(n) for n in model.node_ids if int(n) in result.graph]
    if set(ids) != set(result.graph):
        raise LandmarkError("MeVO ROI lies outside the saved analysis component; select its component explicitly before displaying")
    local = {n: i for i, n in enumerate(ids)}
    indexes = [model.node_index_by_id[n] for n in ids]
    points, radii = model.node_positions_um[indexes].copy(), model.node_radius_um[indexes].copy()
    factor = 1000.0 if units == "mm" else 1.0
    raw = np.asarray([exact.graph.nodes[n]["coords"] for n in ids])
    np.testing.assert_array_equal(points, raw[:, :3] * factor)
    np.testing.assert_array_equal(radii, raw[:, 3] * factor)
    selected_edges = [edge for edge in model.edges if result.graph.has_edge(edge.upstream_node_id, edge.downstream_node_id)]
    if len(selected_edges) != result.graph.number_of_edges():
        raise LandmarkError("Saved global graph edges do not match exact MeVO extraction")
    edges = np.asarray([[local[e.upstream_node_id], local[e.downstream_node_id]] for e in selected_edges], dtype=np.int64)
    global_edges = {(e.upstream_node_id, e.downstream_node_id): e.edge_id for e in model.edges}
    suffix = "" if anatomical_status == "MANUALLY_VERIFIED" else "__UNVERIFIED"
    roi_id = f"{model.source_model_id}__{name}{suffix}"
    boundaries = []
    for n in definition["proximal_nodes"]:
        parent = next(exact.graph.predecessors(n), None)
        boundaries.append((n, global_edges.get((parent, n), -1), "ANATOMICAL_PROXIMAL"))
    for n in definition["distal_nodes"]:
        outgoing = list(exact.graph.successors(n))
        if outgoing:
            boundaries.extend((n, global_edges.get((n, child), -1), "ANATOMICAL_DISTAL") for child in outgoing)
        else:
            boundaries.append((n, -1, "ANATOMICAL_DISTAL_AT_SOURCE_TERMINAL"))
    for excluded in definition["exclude_subtree_roots"]:
        parent = next(exact.graph.predecessors(excluded))
        if parent in local:
            boundaries.append((parent, global_edges.get((parent, excluded), -1), "EXCLUDED_SUBTREE_BOUNDARY"))
    ports = tuple(CutPort(f"{roi_id}__cut_{i:03d}", local[n], edge, tuple(points[local[n]]),
                         float(radii[local[n]]), "anatomical_landmark", role)
                  for i, (n, edge, role) in enumerate(boundaries))
    terminal_ids = tuple(n for n in ids if exact.graph.out_degree(n) == 0 and n not in definition["distal_nodes"])
    # Boxes are display envelopes of the exact subgraph, never anatomical classifiers/cuts.
    padding = max(float(radii.max()) * 1.2, 1.0)
    low, high = points.min(axis=0)-padding, points.max(axis=0)+padding
    length = float(sum(np.linalg.norm(points[b]-points[a]) for a, b in edges))
    record = ROIRecord(roi_id=roi_id, source_model_id=model.source_model_id, source_mouse_id=model.source_mouse_id,
        anchor_id=definition["proximal_nodes"][0], anchor_position_um=tuple(points[local[definition["proximal_nodes"][0]]]),
        bbox_min_um=tuple(low), bbox_max_um=tuple(high), bbox_center_um=tuple((low+high)/2), bbox_size_um=tuple(high-low),
        global_node_ids=tuple(ids), global_edge_ids=tuple(e.edge_id for e in selected_edges),
        local_node_ids=np.arange(len(ids), dtype=np.int64), local_node_global_ids=np.asarray(ids, dtype=np.int64),
        local_node_positions_um=points, local_node_radius_um=radii, local_edges=edges,
        local_edge_ids=np.arange(len(edges), dtype=np.int64), local_edge_global_ids=np.asarray([e.edge_id for e in selected_edges]),
        local_edge_points_um=points[edges], local_edge_radius_um=radii[edges],
        true_terminal_local_ids=tuple(local[n] for n in terminal_ids), true_terminal_global_ids=terminal_ids,
        cut_ports=ports, raw_component_count=1, raw_total_vessel_length_um=length, retained_component_length_um=length,
        cluster_id=list(ROI_DEFINITIONS).index(name), distance_to_cluster_center=0.0,
        is_representative=True, selection_rank=rank)
    populate_roi_features(record)
    # A strict anatomical root may itself be a bifurcation. The spatial ROI
    # descriptors use undirected degree and would merge its two daughter paths.
    # Keep MeVO's displayed statistics consistent with its directed manifest.
    counts = graph_counts(result.graph)
    volume = float(np.prod(high-low))
    tortuosities = []
    for start in result.graph:
        if result.graph.in_degree(start) == 1 and result.graph.out_degree(start) == 1:
            continue
        for child in result.graph.successors(start):
            path = [start, child]
            while result.graph.in_degree(path[-1]) == 1 and result.graph.out_degree(path[-1]) == 1:
                path.append(next(result.graph.successors(path[-1])))
            path_points = points[[local[n] for n in path]]
            chord = float(np.linalg.norm(path_points[-1]-path_points[0]))
            if chord > 1e-12:
                tortuosities.append(float(np.linalg.norm(np.diff(path_points, axis=0), axis=1).sum()/chord))
    record.structural_features.update(branch_count=float(counts["branch_count"]),
        bifurcation_count=float(counts["bifurcation_count"]),
        branch_density_per_um3=counts["branch_count"]/volume,
        bifurcation_density_per_um3=counts["bifurcation_count"]/volume,
        mean_branch_tortuosity=float(np.mean(tortuosities)) if tortuosities else 1.0)
    return record


def run_mevo_from_rodent_run(run_root, config, *, landmarks, allow_natural_terminal=False, verbose=False):
    sample_dirs = [p for p in (Path(run_root)/"samples").iterdir() if p.is_dir()]
    if len(sample_dirs) != 1:
        raise LandmarkError("MeVO display expects a single saved source model")
    manifest = json.loads((sample_dirs[0]/"preprocess_manifest.json").read_text())
    source = Path(manifest["record"]["swc_path"])
    layout = create_sampling_layout(config, run_label="mevo_landmarks")
    extracted = extract_file(source, Path(landmarks), layout.run_root/"mevo", allow_natural_terminal=allow_natural_terminal)
    if extracted["status"] not in {"PASS", "PASS_WITH_WARNINGS"}:
        raise LandmarkError(f"MeVO extraction {extracted['status']}; see {layout.run_root/'mevo'/extracted['manifest']}")
    doc, _ = load_landmarks(Path(landmarks), source)
    exact = load_exact_graph(source)
    models = load_models_from_rodent_run(run_root)
    model, = models
    factor = 1000.0 if extracted["units"] == "mm" else 1.0
    source_values = np.asarray([exact.graph.nodes[int(n)]["coords"] for n in model.node_ids])
    np.testing.assert_array_equal(model.node_positions_um, source_values[:, :3] * factor)
    np.testing.assert_array_equal(model.node_radius_um, source_values[:, 3] * factor)
    records = []
    for name, item in extracted["rois"].items():
        definition = doc["rois"][name]
        result = extract_graph(exact, definition, allow_natural_terminal=allow_natural_terminal)
        records.append(display_record(model, exact, definition, result, name,
                       anatomical_status=item["anatomical_status"], units=extracted["units"], rank=len(records)+1))
    write_sampling_config(layout, config)
    write_global_edge_manifest(layout, models)
    write_roi_library(layout, records)
    write_candidate_tables(layout, records)
    status = "WARNING" if extracted["needs_manual_review"] or extracted["status"] == "PASS_WITH_WARNINGS" else "PASS"
    summary = {"status": status, "mode": "anatomical_landmark_mevo", "source_sha256": extracted["source_sha256"],
               "candidate_count": len(records), "valid_candidate_count": len(records), "selected_count": len(records),
               "mevo_manifest": str(layout.run_root/"mevo"/extracted["manifest"]),
               "geometry_source": "exact original nodes in annotated source, converted once to um for existing renderer",
               "bounding_boxes": "display envelopes only; not anatomical extraction boundaries",
               "clustering_performed": False, "display_cluster_semantics": "stable anatomical ROI groups for existing C-key cycling",
               "branch_count_semantics": "directed source branches; a bifurcating proximal root remains a bifurcation",
               "spatial_sampling_performed": False, "needs_manual_review": extracted["needs_manual_review"],
               "warnings": extracted["warnings"] + [w for item in extracted["rois"].values() for w in item["warnings"]]}
    write_json(layout.summary_file, summary)
    layout.log_file.write_text(json.dumps(summary, ensure_ascii=False, indent=2)+"\n")
    return MeVODisplayRun(layout.run_root, status, records, layout.summary_file, layout.log_file)
