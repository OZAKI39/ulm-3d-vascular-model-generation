"""Fit each original branch with upstream AIC splines and fixed SWC endpoints.

This uses VascularMD's vessel approximation, not its Nfurcation reconstruction.
No merging, radius clipping, or substitute smoothing algorithm is performed.
"""
from __future__ import annotations

import logging

import numpy as np

from third_party.vascularmd.Spline import Spline
from third_party.vascularmd.utils import resample
from .swc_export import validate_tree
from .vascularmd_adapter import Piece

LOG = logging.getLogger(__name__)
JUNCTION_POLICY = (
    "Branch-wise native Spline.approximation with AIC; original shared junction, "
    "root and terminal XYZ/radius are fixed endpoint constraints. "
    "Original connections and per-branch node counts are preserved. "
    "No Nfurcation surface or nearby-bifurcation merge is constructed; "
    "cross-branch tangent continuity is not imposed."
)


def fit_branches(tree, branches, source):
    """Return native samples on the original graph and auditable branch fits."""
    topo = tree.get_topo_graph()
    edges = {(int(topo.nodes[a]["source_swc_id"]), int(topo.nodes[b]["source_swc_id"])): (a, b)
             for a, b in topo.edges}
    if set(edges) != {branch.key for branch in branches}:
        raise ValueError("Preprocessing changed original branch connectivity")
    graph = source.graph.copy()
    fits = []
    for index, branch in enumerate(branches, 1):
        tree.active_component = {"stage": "branch_aic", "source_branch": list(branch.key),
                                 "completed_branches": len(fits), "total_branches": len(branches)}
        LOG.info("AIC branch %d/%d: %s -> %s", index, len(branches), *branch.key)
        a, b = edges[branch.key]
        data = np.vstack((topo.nodes[a]["coords"], topo.edges[a, b]["coords"], topo.nodes[b]["coords"]))
        preprocessed_count = len(data)
        # Match the official __model_vessel minimum-data rule.
        if len(data) <= 6:
            data = resample(data, 6)
        np.testing.assert_allclose(data[[0, -1]], branch.raw[[0, -1]], rtol=1e-12, atol=1e-12)
        data[[0, -1]] = branch.raw[[0, -1]]
        values = np.zeros((4, 4))
        values[[0, -1]] = branch.raw[[0, -1]]
        spline = Spline()
        spline.approximation(data, [True, False, False, True], values, False,
                             radius_model=True, criterion="AIC", akaike=False, max_distance=6)
        dense = spline.get_points()
        if not np.isfinite(dense).all() or np.any(dense[:, 3] <= 0):
            raise ValueError(f"Nonpositive/nonfinite native spline on branch {branch.key}")
        piece = Piece(branch.key, spline, False)
        if not np.isfinite(piece.length) or piece.length <= 0:
            raise ValueError(f"Invalid native arc length on branch {branch.key}")
        arc = np.linspace(0, piece.length, len(branch.raw))
        sampled = piece.evaluate(arc)
        np.testing.assert_allclose(sampled[[0, -1]], branch.raw[[0, -1]], rtol=1e-10, atol=1e-10)
        # Store the exact constrained originals, avoiding endpoint round-off.
        sampled[[0, -1]] = branch.raw[[0, -1]]
        for node, point in zip(branch.raw_ids, sampled):
            graph.nodes[node]["coords"] = point.copy()
        near = np.zeros(len(sampled), dtype=bool)
        near[0], near[-1] = branch.start_bif, branch.end_bif
        branch.smooth, branch.arc, branch.near_bif = sampled, arc, near
        branch.sampled_types, branch.pieces = branch.types.copy(), [piece]
        lambdas = list(map(float, spline.get_lbd()))
        if not np.isfinite(lambdas).all():
            raise ValueError(f"Nonfinite native smoothing parameters on branch {branch.key}")
        fits.append({"source_branch": list(branch.key), "original_count": len(branch.raw),
                     "preprocessed_count": preprocessed_count, "fit_count": len(data),
                     "output_count": len(sampled), "spatial_lambda": lambdas[0], "radius_lambda": lambdas[1],
                     "control_point_count": len(spline.get_control_points()),
                     "dense_radius_min": float(dense[:, 3].min()), "arc_length": piece.length})
    validate_tree(graph)
    if set(graph) != set(source.graph) or set(graph.edges) != set(source.graph.edges):
        raise ValueError("Branch-wise fitting changed original SWC topology")
    for node in graph:
        if source.graph.in_degree(node) != 1 or source.graph.out_degree(node) != 1:
            np.testing.assert_array_equal(graph.nodes[node]["coords"], source.graph.nodes[node]["coords"])
    graph.graph["bifurcation_swc_policy"] = JUNCTION_POLICY
    return graph, fits
