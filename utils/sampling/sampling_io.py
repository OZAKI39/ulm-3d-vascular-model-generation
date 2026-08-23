"""Stable sampling-run layout and machine-readable exports."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import yaml

from .sampling_config import SamplingConfig
from .sampling_types import CutPort, GlobalVascularModel, ROIRecord, SamplingExperiment


@dataclass(frozen=True, slots=True)
class SamplingOutputLayout:
    run_root: Path
    logs: Path
    config: Path
    manifests: Path
    features: Path
    clustering: Path
    figures: Path
    roi_previews: Path
    roi_library: Path
    comparison: Path
    report: Path
    log_file: Path
    summary_file: Path


def create_sampling_layout(config: SamplingConfig) -> SamplingOutputLayout:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base = config.output_root.resolve() / "sampling"
    stem = f"{timestamp}_{config.feature_mode}_k{config.n_clusters}"
    run_root = base / stem
    counter = 1
    while run_root.exists():
        run_root = base / f"{stem}_{counter:02d}"
        counter += 1
    folders = {
        name: run_root / name
        for name in (
            "logs", "config", "manifests", "features", "clustering", "figures",
            "roi_previews", "roi_library", "comparison", "report",
        )
    }
    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=False)
    return SamplingOutputLayout(
        run_root=run_root,
        logs=folders["logs"],
        config=folders["config"],
        manifests=folders["manifests"],
        features=folders["features"],
        clustering=folders["clustering"],
        figures=folders["figures"],
        roi_previews=folders["roi_previews"],
        roi_library=folders["roi_library"],
        comparison=folders["comparison"],
        report=folders["report"],
        log_file=folders["logs"] / "sampling.log",
        summary_file=folders["report"] / "sampling_summary.json",
    )


def write_json(path: Path, payload: Any) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return path


def write_rows(path: Path, rows: Iterable[dict[str, Any]], fieldnames: Iterable[str] | None = None) -> Path:
    payload = list(rows)
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(fieldnames or (payload[0].keys() if payload else ()))
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names, extrasaction="ignore")
        if names:
            writer.writeheader()
            writer.writerows(payload)
    return path


def write_sampling_config(layout: SamplingOutputLayout, config: SamplingConfig) -> Path:
    path = layout.config / "sampling_config.yaml"
    path.write_text(
        yaml.safe_dump(config.report(), sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )
    return path


def write_global_edge_manifest(
    layout: SamplingOutputLayout,
    models: list[GlobalVascularModel],
) -> Path:
    rows = []
    for model in models:
        for edge in model.edges:
            rows.append(
                {
                    "source_model_id": model.source_model_id,
                    "global_edge_id": edge.edge_id,
                    "upstream_global_node_id": edge.upstream_node_id,
                    "downstream_global_node_id": edge.downstream_node_id,
                }
            )
    return write_rows(layout.manifests / "global_edges.csv", rows)


def write_roi_library(layout: SamplingOutputLayout, rois: list[ROIRecord]) -> list[Path]:
    paths: list[Path] = []
    for roi in rois:
        path = layout.roi_library / f"{roi.roi_id}.npz"
        np.savez_compressed(
            path,
            roi_id=np.asarray(roi.roi_id),
            source_model_id=np.asarray(roi.source_model_id),
            anchor_id=np.asarray(roi.anchor_id, dtype=np.int64),
            anchor_position_um=np.asarray(roi.anchor_position_um, dtype=float),
            bbox_min_um=np.asarray(roi.bbox_min_um, dtype=float),
            bbox_max_um=np.asarray(roi.bbox_max_um, dtype=float),
            local_node_ids=roi.local_node_ids,
            local_node_global_ids=roi.local_node_global_ids,
            local_node_positions_um=roi.local_node_positions_um,
            local_node_radius_um=roi.local_node_radius_um,
            local_edges=roi.local_edges,
            local_edge_ids=roi.local_edge_ids,
            local_edge_global_ids=roi.local_edge_global_ids,
            local_edge_points_um=roi.local_edge_points_um,
            local_edge_radius_um=roi.local_edge_radius_um,
            true_terminal_local_ids=np.asarray(roi.true_terminal_local_ids, dtype=np.int64),
            true_terminal_global_ids=np.asarray(roi.true_terminal_global_ids, dtype=np.int64),
            cut_port_local_ids=np.asarray([port.local_node_id for port in roi.cut_ports], dtype=np.int64),
            cut_port_global_edge_ids=np.asarray([port.global_edge_id for port in roi.cut_ports], dtype=np.int64),
            cut_port_positions_um=np.asarray(
                [port.intersection_position_um for port in roi.cut_ports], dtype=float
            ).reshape((-1, 3)),
            cut_port_radius_um=np.asarray([port.radius_at_cut_um for port in roi.cut_ports], dtype=float),
            cut_port_boundary_faces=np.asarray([port.boundary_face for port in roi.cut_ports]),
            cut_port_boundary_roles=np.asarray([port.boundary_role for port in roi.cut_ports]),
            cluster_id=np.asarray(roi.cluster_id, dtype=np.int64),
            selection_rank=np.asarray(roi.selection_rank, dtype=np.int64),
        )
        paths.append(path)
    return paths


def write_candidate_tables(layout: SamplingOutputLayout, rois: list[ROIRecord]) -> list[Path]:
    candidates = write_rows(layout.manifests / "candidate_rois.csv", [roi.manifest_row() for roi in rois])
    selected = write_rows(
        layout.manifests / "selected_rois.csv",
        [roi.manifest_row() for roi in rois if roi.is_representative],
    )
    features = write_rows(layout.features / "roi_features.csv", [roi.feature_row() for roi in rois])
    cut_ports = write_rows(
        layout.manifests / "cut_ports.csv",
        [
            {"roi_id": roi.roi_id, "source_model_id": roi.source_model_id, **port.report()}
            for roi in rois
            for port in roi.cut_ports
        ],
        fieldnames=(
            "roi_id", "source_model_id", "cut_port_id", "local_node_id", "global_edge_id",
            "intersection_x_um", "intersection_y_um", "intersection_z_um", "radius_at_cut_um",
            "boundary_face",
        ),
    )
    return [candidates, selected, features, cut_ports]


def experiment_report(experiment: SamplingExperiment) -> dict[str, Any]:
    return {
        "feature_mode": experiment.feature_mode,
        "features": list(experiment.feature_names),
        "feature_dimension": len(experiment.feature_names),
        "scaler": experiment.scaler.report(),
        "clustering": {
            "method": experiment.clustering.method,
            "n_clusters": experiment.clustering.n_clusters,
            "inertia": experiment.clustering.inertia,
            "silhouette_score": experiment.clustering.silhouette_score,
            "cluster_sizes": list(experiment.clustering.cluster_sizes),
        },
        "selection": {
            "mode": experiment.selection.selection_mode,
            "scientific_label": experiment.selection.scientific_label,
            "requested_count": experiment.selection.requested_count,
            "selected_count": len(experiment.selection.selected_indices),
            "overlap_rejection_count": experiment.selection.overlap_rejection_count,
            "selected_candidate_indices": list(experiment.selection.selected_indices),
        },
        "validation": experiment.validation,
    }


def write_experiment(
    layout: SamplingOutputLayout,
    experiment: SamplingExperiment,
    rois: list[ROIRecord],
    *,
    subdirectory: Path | None = None,
) -> list[Path]:
    base = subdirectory or layout.run_root
    feature_dir = base / "features"
    clustering_dir = base / "clustering"
    manifest_dir = base / "manifests"
    for folder in (feature_dir, clustering_dir, manifest_dir):
        folder.mkdir(parents=True, exist_ok=True)
    scaler_path = write_json(feature_dir / "scaler.json", experiment.scaler.report())
    assignments_path = write_rows(
        clustering_dir / "cluster_assignments.csv",
        [
            {
                "roi_id": roi.roi_id,
                "cluster_id": int(experiment.clustering.assignments[index]),
                "distance_to_center": float(experiment.clustering.distances_to_center[index]),
                "selected": index in experiment.selection.selected_indices,
            }
            for index, roi in enumerate(rois)
        ],
    )
    centers_path = write_rows(
        clustering_dir / "cluster_centers.csv",
        [
            {"cluster_id": cluster_id, **{
                name: float(value)
                for name, value in zip(experiment.feature_names, center)
            }}
            for cluster_id, center in enumerate(experiment.clustering.centers)
        ],
    )
    cluster_summary = write_rows(
        manifest_dir / "cluster_summary.csv",
        [
            {
                "cluster_id": cluster_id,
                "candidate_count": experiment.clustering.cluster_sizes[cluster_id],
                "selected_count": sum(
                    int(experiment.clustering.assignments[index]) == cluster_id
                    for index in experiment.selection.selected_indices
                ),
                "representative_roi_ids": ";".join(
                    rois[index].roi_id
                    for index in experiment.selection.selected_indices
                    if int(experiment.clustering.assignments[index]) == cluster_id
                ),
            }
            for cluster_id in range(experiment.clustering.n_clusters)
        ],
    )
    report_path = write_json(base / "experiment_summary.json", experiment_report(experiment))
    return [scaler_path, assignments_path, centers_path, cluster_summary, report_path]


def load_sampling_display_rois(run_root: Path, *, selected_only: bool = False) -> list[ROIRecord]:
    """Load saved geometry for the display-only UI layer."""

    manifest_name = "selected_rois.csv" if selected_only else "candidate_rois.csv"
    manifest_path = Path(run_root) / "manifests" / manifest_name
    if not manifest_path.is_file():
        return []
    with manifest_path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    feature_path = Path(run_root) / "features" / "roi_features.csv"
    feature_by_id: dict[str, dict[str, str]] = {}
    if feature_path.is_file():
        with feature_path.open(encoding="utf-8-sig", newline="") as stream:
            feature_by_id = {row["roi_id"]: row for row in csv.DictReader(stream)}
    output: list[ROIRecord] = []
    for row in rows:
        archive_path = Path(run_root) / "roi_library" / f"{row['roi_id']}.npz"
        if not archive_path.is_file():
            continue
        with np.load(archive_path) as archive:
            cut_positions = np.asarray(archive["cut_port_positions_um"], dtype=float).reshape((-1, 3))
            cut_roles = (
                archive["cut_port_boundary_roles"]
                if "cut_port_boundary_roles" in archive.files
                else np.full(len(cut_positions), "CORE_ROI_BOUNDARY")
            )
            cut_ports = tuple(
                CutPort(
                    cut_port_id=f"{row['roi_id']}__cut_{index:03d}",
                    local_node_id=int(local_id),
                    global_edge_id=int(edge_id),
                    intersection_position_um=tuple(map(float, position)),
                    radius_at_cut_um=float(radius),
                    boundary_face=str(face),
                    boundary_role=str(role),
                )
                for index, (local_id, edge_id, position, radius, face, role) in enumerate(
                    zip(
                        archive["cut_port_local_ids"], archive["cut_port_global_edge_ids"],
                        cut_positions, archive["cut_port_radius_um"],
                        archive["cut_port_boundary_faces"], cut_roles,
                    )
                )
            )
            roi = ROIRecord(
                roi_id=row["roi_id"],
                source_model_id=row["source_model_id"],
                source_mouse_id=row.get("source_mouse_id", ""),
                anchor_id=int(row["anchor_id"]),
                anchor_position_um=tuple(map(float, archive["anchor_position_um"])),
                bbox_min_um=tuple(map(float, archive["bbox_min_um"])),
                bbox_max_um=tuple(map(float, archive["bbox_max_um"])),
                bbox_center_um=tuple((np.asarray(archive["bbox_min_um"]) + np.asarray(archive["bbox_max_um"])) * 0.5),
                bbox_size_um=tuple(np.asarray(archive["bbox_max_um"]) - np.asarray(archive["bbox_min_um"])),
                global_node_ids=tuple(int(value) for value in archive["local_node_global_ids"] if value >= 0),
                global_edge_ids=tuple(map(int, archive["local_edge_global_ids"])),
                local_node_ids=np.asarray(archive["local_node_ids"], dtype=np.int64),
                local_node_global_ids=np.asarray(archive["local_node_global_ids"], dtype=np.int64),
                local_node_positions_um=np.asarray(archive["local_node_positions_um"], dtype=float),
                local_node_radius_um=np.asarray(archive["local_node_radius_um"], dtype=float),
                local_edges=np.asarray(archive["local_edges"], dtype=np.int64),
                local_edge_ids=np.asarray(archive["local_edge_ids"], dtype=np.int64),
                local_edge_global_ids=np.asarray(archive["local_edge_global_ids"], dtype=np.int64),
                local_edge_points_um=np.asarray(archive["local_edge_points_um"], dtype=float),
                local_edge_radius_um=np.asarray(archive["local_edge_radius_um"], dtype=float),
                true_terminal_local_ids=tuple(map(int, archive["true_terminal_local_ids"])),
                true_terminal_global_ids=tuple(map(int, archive["true_terminal_global_ids"])),
                cut_ports=cut_ports,
                raw_component_count=int(row["raw_component_count"]),
                raw_total_vessel_length_um=float(row["raw_total_vessel_length_um"]),
                retained_component_length_um=float(row["retained_component_length_um"]),
                cluster_id=int(row["cluster_id"]),
                distance_to_cluster_center=float(row["distance_to_center"]),
                is_representative=str(row["selected"]).lower() == "true",
                selection_rank=int(row["selection_rank"]),
            )
            feature_row = feature_by_id.get(roi.roi_id, row)
            for key in ("r10", "r25", "r50", "r75", "r90"):
                if key in feature_row and feature_row[key]:
                    roi.radius_features[key] = float(feature_row[key])
            for key in ("branch_count", "bifurcation_count", "total_vessel_length_um", "cycle_rank"):
                if key in feature_row and feature_row[key]:
                    roi.structural_features[key] = float(feature_row[key])
            output.append(roi)
    return output
