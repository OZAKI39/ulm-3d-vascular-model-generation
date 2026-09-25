"""Readable visual acceptance artifacts for the directed Schmid graph."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection

from ..reporting.acceptance import AcceptanceResult
from .config import SchmidPKLConfig
from .model import DirectedGraphResult, VESSEL_TYPE_NAMES


ROLE_COLORS = {
    "source": "#1b9e77",
    "source_split": "#00a878",
    "sink": "#377eb8",
    "merge_sink": "#235789",
    "split": "#e6ab02",
    "merge": "#d95f02",
    "mixed_split_merge": "#984ea3",
    "unresolved_junction": "#e41a1c",
    "unresolved_terminal": "#ff5a5f",
    "connector": "#777777",
    "isolated": "#111111",
}


def _segments(result: DirectedGraphResult, smoothed: bool = True) -> list[np.ndarray]:
    return [
        item.points_smoothed_um if smoothed else item.points_raw_um for item in result.branches
    ]


def _all_points(result: DirectedGraphResult) -> np.ndarray:
    return np.concatenate(_segments(result), axis=0)


def _set_equal_3d(axis: plt.Axes, points: np.ndarray) -> None:
    minimum = points.min(axis=0)
    maximum = points.max(axis=0)
    center = (minimum + maximum) / 2
    radius = max(float(np.max(maximum - minimum)) / 2, 1.0)
    axis.set_xlim(center[0] - radius, center[0] + radius)
    axis.set_ylim(center[1] - radius, center[1] + radius)
    axis.set_zlim(center[2] - radius, center[2] + radius)
    axis.set_box_aspect((1, 1, 1))


def _project(sequence: np.ndarray, view: int) -> np.ndarray:
    axes = ((1, 2), (0, 2), (0, 1))[view]
    return sequence[:, axes]


def _render_3d(result: DirectedGraphResult, output_dir: Path) -> Path:
    path = output_dir / "directed_hierarchical_graph_3d.png"
    figure = plt.figure(figsize=(13, 10), constrained_layout=True)
    axis = figure.add_subplot(111, projection="3d")
    regular = [item.points_smoothed_um for item in result.branches if item.direction_status == "known"]
    unresolved = [item.points_smoothed_um for item in result.branches if item.direction_status != "known"]
    axis.add_collection3d(
        Line3DCollection(regular, colors="#527d99", linewidths=0.38, alpha=0.55),
        autolim=False,
    )
    if unresolved:
        axis.add_collection3d(
            Line3DCollection(unresolved, colors="#e41a1c", linewidths=1.4, alpha=0.95),
            autolim=False,
        )
    _set_equal_3d(axis, _all_points(result))
    axis.set_xlabel("Source X (um)")
    axis.set_ylabel("Source Y (um)")
    axis.set_zlabel("Source Z (um)")
    axis.set_title(
        "Directed Schmid hierarchical graph\nBlue: pressure-oriented; red: direction unresolved"
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_three_views(result: DirectedGraphResult, output_dir: Path) -> Path:
    path = output_dir / "directed_graph_three_views.png"
    figure, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    for view, axis in enumerate(axes):
        known = [_project(item.points_smoothed_um, view) for item in result.branches if item.direction_status == "known"]
        unresolved = [_project(item.points_smoothed_um, view) for item in result.branches if item.direction_status != "known"]
        axis.add_collection(LineCollection(known, colors="#537d96", linewidths=0.3, alpha=0.55))
        if unresolved:
            axis.add_collection(LineCollection(unresolved, colors="#e41a1c", linewidths=1.0, alpha=0.95))
        axis.autoscale()
        axis.set_aspect("equal", adjustable="box")
        axis.axis("off")
        axis.set_title(("X view (YZ)", "Y view (XZ)", "Z view (XY)")[view])
    figure.suptitle("Whole-network directed connectivity; unresolved branches are red")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_pressure_flow(result: DirectedGraphResult, output_dir: Path) -> Path:
    path = output_dir / "pressure_and_flow_maps.png"
    figure, axes = plt.subplots(1, 2, figsize=(16, 7), constrained_layout=True)
    sequences = [_project(item.points_smoothed_um, 2) for item in result.branches]
    pressure = result.cleanup.source.pressure_mmhg
    pressure_values = np.asarray(
        [(pressure[item.node_a] + pressure[item.node_b]) / 2 for item in result.branches]
    )
    pressure_collection = LineCollection(
        sequences,
        array=pressure_values,
        cmap="viridis",
        norm=Normalize(float(pressure_values.min()), float(pressure_values.max())),
        linewidths=0.55,
    )
    axes[0].add_collection(pressure_collection)
    axes[0].autoscale()
    figure.colorbar(pressure_collection, ax=axes[0], label="Pressure (mmHg)")
    flow_values = np.log10(
        np.maximum(
            np.asarray([item.flow_um3_per_ms for item in result.branches]),
            np.finfo(float).tiny,
        )
    )
    flow_collection = LineCollection(
        sequences,
        array=flow_values,
        cmap="plasma",
        norm=Normalize(float(flow_values.min()), float(flow_values.max())),
        linewidths=0.55,
    )
    axes[1].add_collection(flow_collection)
    axes[1].autoscale()
    figure.colorbar(flow_collection, ax=axes[1], label="log10 flow (um3/ms)")
    for axis, title in zip(axes, ("Mean endpoint pressure", "Flow magnitude"), strict=True):
        axis.set_aspect("equal", adjustable="box")
        axis.axis("off")
        axis.set_title(title)
    figure.suptitle("Z/XY projection: simulated pressure and flow")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_vessel_types(result: DirectedGraphResult, output_dir: Path) -> Path:
    path = output_dir / "vessel_type_map.png"
    palette = {
        0: "#d73027",
        1: "#4575b4",
        2: "#fc8d59",
        3: "#74add1",
        4: "#666666",
        5: "#984ea3",
    }
    figure, axis = plt.subplots(figsize=(11, 9), constrained_layout=True)
    for code in sorted(VESSEL_TYPE_NAMES):
        sequences = [
            _project(item.points_smoothed_um, 2)
            for item in result.branches
            if code in item.vessel_type_codes
        ]
        if sequences:
            axis.add_collection(
                LineCollection(
                    sequences,
                    colors=palette[code],
                    linewidths=0.45 if code == 4 else 0.85,
                    alpha=0.6 if code == 4 else 0.9,
                    label=f"{code}: {VESSEL_TYPE_NAMES[code]}",
                )
            )
    axis.autoscale()
    axis.set_aspect("equal", adjustable="box")
    axis.axis("off")
    axis.legend(loc="best", fontsize=8)
    axis.set_title("Source vessel classes (Z/XY projection)")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_node_roles(result: DirectedGraphResult, output_dir: Path) -> Path:
    path = output_dir / "directed_node_roles.png"
    figure, axis = plt.subplots(figsize=(11, 9), constrained_layout=True)
    sequences = [_project(item.points_smoothed_um, 2) for item in result.branches]
    axis.add_collection(LineCollection(sequences, colors="#c3c8cc", linewidths=0.25, alpha=0.35))
    for role, color in ROLE_COLORS.items():
        coordinates = np.asarray(
            [node.coordinates_um for node in result.nodes if node.node_role == role], dtype=float
        )
        if len(coordinates):
            axis.scatter(
                coordinates[:, 0], coordinates[:, 1], s=8, c=color, label=role, alpha=0.8
            )
    axis.autoscale()
    axis.set_aspect("equal", adjustable="box")
    axis.axis("off")
    axis.legend(loc="best", fontsize=7, ncol=2)
    axis.set_title("Sources, sinks, splits, merges, and unresolved junctions")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_cycles(result: DirectedGraphResult, output_dir: Path) -> Path:
    path = output_dir / "undirected_cycle_overview.png"
    figure, axes = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    for view, axis in enumerate(axes):
        non_cycle = [
            _project(item.points_smoothed_um, view)
            for item in result.branches
            if not item.cycle_ids
        ]
        cyclic = [
            _project(item.points_smoothed_um, view)
            for item in result.branches
            if item.cycle_ids
        ]
        if non_cycle:
            axis.add_collection(
                LineCollection(non_cycle, colors="#c7ccd0", linewidths=0.25, alpha=0.35)
            )
        if cyclic:
            axis.add_collection(
                LineCollection(cyclic, colors="#f28e2b", linewidths=0.55, alpha=0.72)
            )
        axis.autoscale()
        axis.set_aspect("equal", adjustable="box")
        axis.axis("off")
        axis.set_title(("X view (YZ)", "Y view (XZ)", "Z view (XY)")[view])
    figure.suptitle(
        f"Undirected network redundancy: cycle rank {len(result.cycles):,}; orange branches occur in the stored basis"
    )
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_direction_arrows(
    result: DirectedGraphResult, config: SchmidPKLConfig, output_dir: Path
) -> Path:
    path = output_dir / "flow_direction_arrows.png"
    candidates = [
        item
        for item in result.branches
        if item.direction_status == "known" and item.straight_distance_um > 0
    ]
    candidates.sort(key=lambda item: item.flow_um3_per_ms, reverse=True)
    if len(candidates) > config.max_direction_arrows:
        sample_indices = np.linspace(0, len(candidates) - 1, config.max_direction_arrows).astype(int)
        selected = [candidates[index] for index in sample_indices]
    else:
        selected = candidates
    figure, axis = plt.subplots(figsize=(12, 9), constrained_layout=True)
    sequences = [_project(item.points_smoothed_um, 2) for item in result.branches]
    axis.add_collection(LineCollection(sequences, colors="#bfc6cb", linewidths=0.25, alpha=0.35))
    positions: list[np.ndarray] = []
    vectors: list[np.ndarray] = []
    for branch in selected:
        points = branch.points_smoothed_um
        midpoint = len(points) // 2
        left = max(0, midpoint - 1)
        right = min(len(points) - 1, midpoint + 1)
        vector = points[right, :2] - points[left, :2]
        norm = np.linalg.norm(vector)
        if norm > np.finfo(float).eps:
            positions.append(points[midpoint, :2])
            vectors.append(vector / norm)
    if positions:
        position_array = np.asarray(positions)
        vector_array = np.asarray(vectors)
        axis.quiver(
            position_array[:, 0],
            position_array[:, 1],
            vector_array[:, 0],
            vector_array[:, 1],
            color="#d73027",
            angles="xy",
            scale_units="xy",
            scale=0.035,
            width=0.002,
            alpha=0.75,
        )
    axis.autoscale()
    axis.set_aspect("equal", adjustable="box")
    axis.axis("off")
    axis.set_title(f"Pressure-derived flow arrows (shown: {len(positions):,})")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_distributions(result: DirectedGraphResult, output_dir: Path) -> Path:
    path = output_dir / "branch_morphology_distributions.png"
    lengths = np.asarray([item.source_length_um for item in result.branches])
    radii = np.asarray([np.mean(item.radius_raw_um) for item in result.branches])
    flows = np.asarray([item.flow_um3_per_ms for item in result.branches])
    tortuosity = np.asarray(
        [item.tortuosity_raw for item in result.branches if item.tortuosity_raw is not None]
    )
    figure, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
    axes[0, 0].hist(lengths, bins=50, color="#4c78a8")
    axes[0, 0].set_title("Source branch length")
    axes[0, 0].set_xlabel("um")
    axes[0, 1].hist(radii, bins=50, color="#9c6ade")
    axes[0, 1].set_title("Mean source radius")
    axes[0, 1].set_xlabel("um")
    axes[1, 0].hist(np.log10(np.maximum(flows, np.finfo(float).tiny)), bins=50, color="#e45756")
    axes[1, 0].set_title("Flow magnitude")
    axes[1, 0].set_xlabel("log10 um3/ms")
    axes[1, 1].hist(tortuosity, bins=50, color="#54a24b")
    axes[1, 1].set_title("Centerline tortuosity")
    axes[1, 1].set_xlabel("length / straight distance")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_flow_qc(result: DirectedGraphResult, output_dir: Path) -> Path:
    path = output_dir / "flow_conservation_qc.png"
    rows = [
        row
        for row in result.flow_conservation
        if not row["is_pressure_boundary"]
        and row["unresolved_incident_edge_count"] == 0
        and row["incoming_flow_um3_per_ms"] > 0
        and row["outgoing_flow_um3_per_ms"] > 0
    ]
    relative = np.asarray([row["relative_imbalance"] for row in rows], dtype=float)
    figure, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
    positive = relative[relative > 0]
    axes[0].hist(np.log10(positive), bins=60, color="#4c78a8")
    axes[0].set_xlabel("log10 relative imbalance")
    axes[0].set_title("Internal flow-conservation error")
    counts = result.report()["node_role_counts"]
    labels = list(counts)
    values = [counts[label] for label in labels]
    axes[1].barh(labels, values, color=[ROLE_COLORS.get(label, "#777777") for label in labels])
    axes[1].set_title("Directed node roles")
    axes[1].set_xlabel("node count")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def _render_cleanup(result: DirectedGraphResult, output_dir: Path) -> Path:
    path = output_dir / "cleanup_summary.png"
    report = result.cleanup.report()
    geometry = report["geometry_status_counts"]
    actions = report["decision_action_counts"]
    figure, axes = plt.subplots(1, 2, figsize=(14, 5), constrained_layout=True)
    axes[0].bar(list(actions), list(actions.values()), color="#4c78a8")
    axes[0].tick_params(axis="x", rotation=25)
    axes[0].set_title("Step 1 decisions")
    axes[0].set_ylabel("record count")
    axes[1].bar(list(geometry), list(geometry.values()), color="#f58518")
    axes[1].tick_params(axis="x", rotation=25)
    axes[1].set_title("Source centerline geometry status")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path


def render_schmid_visualizations(
    result: DirectedGraphResult, config: SchmidPKLConfig, output_dir: Path
) -> list[Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    return [
        _render_3d(result, output_dir),
        _render_three_views(result, output_dir),
        _render_pressure_flow(result, output_dir),
        _render_vessel_types(result, output_dir),
        _render_node_roles(result, output_dir),
        _render_cycles(result, output_dir),
        _render_direction_arrows(result, config, output_dir),
        _render_distributions(result, output_dir),
        _render_flow_qc(result, output_dir),
        _render_cleanup(result, output_dir),
    ]


def render_acceptance_dashboard(
    acceptance: AcceptanceResult, output_dir: Path
) -> Path:
    path = output_dir / "directed_graph_acceptance_dashboard.png"
    counts = {
        status: sum(item.status == status for item in acceptance.checks)
        for status in ("PASS", "WARNING", "FAIL")
    }
    figure, axis = plt.subplots(figsize=(9, 5), constrained_layout=True)
    bars = axis.bar(
        list(counts), list(counts.values()), color=["#2ca02c", "#ffbf00", "#d62728"]
    )
    axis.bar_label(bars)
    axis.set_ylabel("check count")
    axis.set_title(f"Directed graph acceptance: {acceptance.overall_status}")
    figure.savefig(path, dpi=180)
    plt.close(figure)
    return path
