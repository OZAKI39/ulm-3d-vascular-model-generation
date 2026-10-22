"""Descriptive QC only: these functions never modify radius or geometry."""
from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


def arc_length(points):
    return np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(points[:, :3], axis=0), axis=1))]


def distribution(values):
    values = np.asarray(values, dtype=float)
    if not len(values):
        return {"number_of_evaluated_pairs": 0, **{key: None for key in ("median", "p90", "p95", "p99", "maximum")}}
    return {"number_of_evaluated_pairs": len(values),
            **dict(zip(("median", "p90", "p95", "p99", "maximum"), map(float, np.percentile(values, [50, 90, 95, 99, 100]))))}


def jumps(branch, points, *, raw=False):
    radius = points[:, 3]
    log_jump = abs(np.diff(np.log(radius)))
    relative = abs(np.diff(radius)) / ((radius[1:] + radius[:-1]) / 2)
    valid = np.ones(len(log_jump), dtype=bool)
    # Exclude every pair touching a shared bifurcation node, independently per branch.
    if branch.start_bif:
        valid[0] = False
    if branch.end_bif:
        valid[-1] = False
    if raw:
        for i in branch.raw_bifurcation_indices:
            if i > 0:
                valid[i - 1] = False
            if i < len(valid):
                valid[i] = False
    return log_jump, relative, valid


def graph_statistics(graph, branch_count):
    radius = np.asarray([graph.nodes[n]["coords"][3] for n in graph])
    return {"node_count": len(graph), "root_count": sum(graph.in_degree(n) == 0 for n in graph),
            "terminal_count": sum(graph.out_degree(n) == 0 for n in graph),
            "bifurcation_count": sum(graph.out_degree(n) > 1 for n in graph), "branch_count": branch_count,
            "radius_min": float(radius.min()), "radius_median": float(np.median(radius)), "radius_max": float(radius.max()),
            "total_centerline_length": float(sum(np.linalg.norm(graph.nodes[a]["coords"][:3] - graph.nodes[b]["coords"][:3]) for a, b in graph.edges))}


def compare(branches, raw_graph, smooth_graph, raw_branches=None):
    raw_branches = branches if raw_branches is None else raw_branches
    report = {"pair_exclusion": "Pairs touching a bifurcation node excluded; no threshold-based filtering",
              "sampling_caveat": "Jump magnitude depends on sampling density; original-count preserves each mapped source-path count. Native merges can change total count by duplicating short shared source trunks. Spacing mode changes density.",
              "raw": graph_statistics(raw_graph, len(raw_branches)), "smooth": graph_statistics(smooth_graph, len(branches)), "branches": []}
    report["raw_baseline"] = "Original topology branches counted once, including when native merging duplicates a source trunk across mapped paths"
    for name in ("raw", "smooth"):
        logs, relatives, gradients = [], [], []
        for branch in raw_branches if name == "raw" else branches:
            values = getattr(branch, name)
            log_jump, relative, valid = jumps(branch, values, raw=name == "raw")
            logs.extend(log_jump[valid])
            relatives.extend(relative[valid])
            distance = np.diff(arc_length(values))
            mask = valid & (distance > 0)
            gradients.extend(log_jump[mask] / distance[mask])
        report[name]["log_radius_jump"] = distribution(logs)
        report[name]["relative_radius_jump"] = distribution(relatives)
        report[name]["abs_log_radius_gradient_per_input_unit"] = distribution(gradients)
    for branch in branches:
        item = {"branch_id": f"{branch.key[0]}->{branch.key[1]}", "raw_count": len(branch.raw), "smooth_count": len(branch.smooth),
                "model_edges": [list(piece.edge) for piece in branch.pieces],
                "native_piece_lengths": [piece.length for piece in branch.pieces],
                "source_type_values": sorted(set(map(int, branch.types)))}
        for name in ("raw", "smooth"):
            log_jump, _, valid = jumps(branch, getattr(branch, name), raw=name == "raw")
            item[name + "_log_radius_jump"] = distribution(log_jump[valid])
        report["branches"].append(item)
    return report


def write_radius_csv(path: Path, branches):
    # A long table avoids fabricating raw/smooth point correspondences after resampling.
    fields = ["branch_id", "series", "sample_index", "arc_length", "raw_radius", "smooth_radius",
              "relative_jump_raw", "relative_jump_smooth", "is_near_bifurcation", "pair_evaluated", "source_node_id"]
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for branch in branches:
            for name in ("raw", "smooth"):
                values = getattr(branch, name)
                _, relative, valid = jumps(branch, values, raw=name == "raw")
                arc = arc_length(values) if name == "raw" else branch.arc
                near = np.zeros(len(values), dtype=bool) if name == "raw" else branch.near_bif.copy()
                if name == "raw":
                    for i in branch.raw_bifurcation_indices:
                        near[max(0, i - 1):i + 2] = True
                if branch.start_bif:
                    near[:2] = True
                if branch.end_bif:
                    near[-2:] = True
                for i in range(len(values)):
                    writer.writerow({"branch_id": f"{branch.key[0]}->{branch.key[1]}", "series": name, "sample_index": i,
                                     "arc_length": float(arc[i]), name + "_radius": float(values[i, 3]),
                                     "relative_jump_" + name: float(relative[i - 1]) if i > 0 else "",
                                     "is_near_bifurcation": bool(near[i]), "pair_evaluated": bool(i > 0 and valid[i - 1]),
                                     "source_node_id": branch.raw_ids[i] if name == "raw" else ""})


def plot_diagnostics(path: Path, branches, count: int = 6, *, input_units="mm", bifurcation_label="Near native bifurcation"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    def severity(branch):
        values, _, valid = jumps(branch, branch.raw, raw=True)
        return float(values[valid].max()) if valid.any() else 0.0

    selected = sorted(branches, key=severity, reverse=True)[:count]
    fig, axes = plt.subplots(len(selected), 1, figsize=(10, 3.2 * len(selected)), squeeze=False, constrained_layout=True)
    for axis, branch in zip(axes[:, 0], selected):
        axis.plot(arc_length(branch.raw), branch.raw[:, 3], "o-", color="#8d5b36", markersize=3, linewidth=1, label="Raw SWC")
        axis.plot(branch.arc, branch.smooth[:, 3], "o-", color="#087f8c", markersize=2, label="VascularMD model samples")
        near = branch.near_bif
        indexes = np.flatnonzero(near)
        if indexes.size:
            groups = np.split(indexes, np.flatnonzero(np.diff(indexes) > 1) + 1)
            for i, group in enumerate(groups):
                lo, hi = max(0, group[0] - 1), min(len(near) - 1, group[-1] + 1)
                axis.axvspan(branch.arc[lo], branch.arc[hi], alpha=0.12, color="#bb6e20", label=bifurcation_label if i == 0 else None)
        axis.set(title=f"Source path {branch.key[0]} -> {branch.key[1]}", xlabel=f"Arc length ({input_units})", ylabel=f"Radius ({input_units})")
        axis.grid(alpha=0.2)
        axis.legend(loc="best")
    fig.savefig(path, dpi=160)
    plt.close(fig)
