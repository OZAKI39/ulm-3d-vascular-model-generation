"""Adapt the pinned upstream model to SWC without refitting its geometry."""
from __future__ import annotations

import hashlib
import logging
import inspect
from contextlib import contextmanager
from dataclasses import dataclass, field

import networkx as nx
import numpy as np

from third_party.vascularmd.ArterialTree import ArterialTree
from third_party.vascularmd.Spline import Spline
from .swc_export import Source, validate_tree

LOG = logging.getLogger(__name__)
COMMIT = "770feb8fbfb591d6d43f966db57b1465010b9a13"


class TopologyChanged(ValueError):
    def __init__(self, context):
        self.context = context
        super().__init__(f"Native model changed input topology: {context}")


class NativeApexNotFound(ValueError):
    def __init__(self, context):
        self.context = context
        super().__init__(f"Native bifurcation apex search exhausted all available outlet data: {context}")


@contextmanager
def guard_exhausted_apex_search():
    """Stop the upstream retry loop only when no additional data can be added.

    This observes native results without changing fitting or intersection math.
    Upstream retries with nb_min += 10 even after every outlet reaches its tip;
    at that point all later fits would repeat the same finite data indefinitely.
    """
    original = Spline.first_intersectionv2

    def checked(spline, other):
        result = original(spline, other)
        if 1.0 not in result[1]:
            return result
        frame = inspect.currentframe().f_back
        try:
            if frame.f_code.co_name != "__model_furcation":
                return result
            state = frame.f_locals
            tree, node = state["self"], state["n"]
            graph = tree.get_model_graph()
            counts = []
            for edge in state["e_out"]:
                count = 1 + len(graph.edges[edge]["coords"]) + 1
                tip = edge[1]
                while graph.out_degree(tip):
                    child = next(iter(graph.successors(tip)))
                    count += len(graph.edges[tip, child]["coords"]) + 1
                    tip = child
                counts.append(count)
            if state["nb_min"] >= max(counts):
                raise NativeApexNotFound({"stage": "native_apex_search", "model_node": int(node),
                    "source_swc_node": graph.nodes[node].get("source_swc_id"),
                    "requested_outlet_points": state["nb_min"], "available_outlet_points": counts})
            return result
        finally:
            del frame

    Spline.first_intersectionv2 = checked
    try:
        yield
    finally:
        Spline.first_intersectionv2 = original


@contextmanager
def cache_native_distance_points():
    """Reuse the identical native sample array within one read-only distance call.

    Upstream distance projects each point separately and each projection converts
    the whole geomdl point list to a NumPy array. Cache only for that call, then
    restore the getter before any geometry can change. Projection math is intact.
    """
    original = Spline.distance

    def cached(spline, data):
        points = spline.get_points()
        previous = spline.__dict__.get("get_points")
        spline.get_points = lambda: points
        try:
            return original(spline, data)
        finally:
            if previous is None:
                del spline.get_points
            else:
                spline.get_points = previous

    Spline.distance = cached
    try:
        yield
    finally:
        Spline.distance = original

JUNCTION_POLICY = (
    "Nfurcation trajectory splines have zero-radius placeholders. Positions use the native trajectory; "
    "radii use its corresponding native shape spline at the nearest projected centerline position. "
    "The inlet uses shape spline 0; each outlet uses its own shape spline. "
    "A shared junction has one representative radius from inlet shape 0, without averaging. "
    "SWC cannot encode the complete noncircular bifurcation; the native VTK surface is authoritative."
)


@dataclass
class Branch:
    key: tuple[int, int]
    raw_ids: list[int]
    raw: np.ndarray
    types: np.ndarray
    start_bif: bool
    end_bif: bool
    smooth: np.ndarray | None = None
    arc: np.ndarray | None = None
    near_bif: np.ndarray | None = None
    sampled_types: np.ndarray | None = None
    pieces: list = field(default_factory=list)
    raw_bifurcation_indices: list[int] = field(default_factory=list)


class TracedTree(ArterialTree):
    """Progress/error context around upstream calls; mathematical code unchanged."""
    active_component = None

    def _ArterialTree__model_furcation(self, node, *args, **kwargs):
        self.active_component = {"stage": "model_bifurcation", "model_node": int(node),
                                 "source_swc_node": self._topo_graph.nodes[node].get("source_swc_id")}
        LOG.info("Model bifurcation %s", self.active_component)
        return super()._ArterialTree__model_furcation(node, *args, **kwargs)

    def _ArterialTree__model_vessel(self, edge, *args, **kwargs):
        self.active_component = {"stage": "model_vessel", "model_edge": list(edge)}
        LOG.info("Model vessel %s", edge)
        return super()._ArterialTree__model_vessel(edge, *args, **kwargs)

    def furcation_cross_sections(self, node):
        self.active_component = {"stage": "mesh_bifurcation", "model_node": int(node)}
        LOG.info("Mesh bifurcation %s", node)
        return super().furcation_cross_sections(node)

    def vessel_cross_sections(self, edge):
        self.active_component = {"stage": "mesh_vessel", "model_edge": list(edge)}
        return super().vessel_cross_sections(edge)


def load_tree(source: Source, auto_resample: bool | None):
    # Capture the native topology/full_id mapping BEFORE official preprocessing.
    # Constructor default is True; the deferred call has the same native behavior.
    tree = TracedTree(source.path.stem, "SWC", str(source.path), automatic_resampling=False)
    if not tree.check_full_graph():
        bad = [(n, tree.get_full_graph().in_degree(n), tree.get_full_graph().out_degree(n))
               for n in tree.get_full_graph() if tree.get_full_graph().out_degree(n) > 4]
        raise ValueError(f"VascularMD check_full_graph rejected the input: {bad}")
    root = next(n for n in source.graph if source.graph.in_degree(n) == 0)
    if source.graph.out_degree(root) != 1:
        raise ValueError(f"SWC root {root}: VascularMD would treat a branching inlet as a sink; manual review required")
    topo = tree.get_topo_graph()
    for node in topo:
        topo.nodes[node]["source_swc_id"] = int(topo.nodes[node]["full_id"])
    original_topo = topo.copy()
    branches = []
    for start, end, data in original_topo.edges(data=True):
        ids = [int(original_topo.nodes[start]["full_id"]), *map(int, data["full_id"]), int(original_topo.nodes[end]["full_id"])]
        branches.append(Branch((ids[0], ids[-1]), ids, np.asarray([source.graph.nodes[n]["coords"] for n in ids]),
                               np.asarray([source.graph.nodes[n]["swc_type"] for n in ids]),
                               original_topo.out_degree(start) > 1, original_topo.out_degree(end) > 1))
    if auto_resample is not False:
        tree.automatic_resampling()
    if not tree.check_full_graph() or not nx.is_arborescence(tree.get_full_graph()):
        raise ValueError("VascularMD preprocessing produced an invalid directed tree")
    return tree, original_topo, branches


def prepare_export_topology(tree, original_topo, branches, source: Source, accept_native_merges=False):
    """Accept ONLY native contractions of original bifurcations, never lost tips/reparenting.

    VascularMD keeps source_swc_id on surviving topo nodes even when full_id is
    renumbered. Use that provenance, not coordinate proximity, to recover paths.
    """
    info = {"policy": "accept-native-merges" if accept_native_merges else "strict",
            "merged_bifurcation_swc_ids": [], "merge_groups": [], "input_topology_preserved": True}
    if not accept_native_merges:
        return original_topo, branches, info
    native_topo = tree.get_topo_graph()
    if not nx.is_arborescence(native_topo):
        raise ValueError("Native merged topology is not a connected outward tree")
    original_nodes = {int(data["source_swc_id"]) for _, data in original_topo.nodes(data=True)}
    retained = {int(data["source_swc_id"]) for _, data in native_topo.nodes(data=True)}
    if len(retained) != len(native_topo) or not retained <= original_nodes:
        raise ValueError("Native topology has duplicated or untraceable source nodes")
    raw_root = next(n for n in source.graph if source.graph.in_degree(n) == 0)
    native_root = next(n for n in native_topo if native_topo.in_degree(n) == 0)
    raw_tips = {n for n in source.graph if source.graph.out_degree(n) == 0}
    native_tips = {int(native_topo.nodes[n]["source_swc_id"]) for n in native_topo if native_topo.out_degree(n) == 0}
    if int(native_topo.nodes[native_root]["source_swc_id"]) != raw_root or native_tips != raw_tips:
        raise ValueError(f"Accepting merges never permits changing the root or terminals: missing={sorted(raw_tips-native_tips)}, added={sorted(native_tips-raw_tips)}")
    removed = sorted(original_nodes - retained)
    if any(source.graph.out_degree(node) < 2 for node in removed):
        raise ValueError(f"Native model removed a non-bifurcation node: {removed}")
    parents = {node: next(iter(source.graph.predecessors(node)), None) for node in source.graph}
    paths = []
    for a, b in native_topo.edges:
        start, end = (int(native_topo.nodes[n]["source_swc_id"]) for n in (a, b))
        reverse_path = [end]
        node = parents[end]
        while node is not None and node not in retained:
            reverse_path.append(node)
            node = parents[node]
        if node != start:
            raise ValueError(f"Native topology change is not an allowed bifurcation contraction: {start}->{end}")
        ids = [start, *reversed(reverse_path)]
        paths.append(Branch((start, end), ids,
                            np.asarray([source.graph.nodes[n]["coords"] for n in ids]),
                            np.asarray([source.graph.nodes[n]["swc_type"] for n in ids]),
                            native_topo.out_degree(a) > 1, native_topo.out_degree(b) > 1,
                            raw_bifurcation_indices=[i for i, n in enumerate(ids) if source.graph.out_degree(n) > 1]))
    groups = {}
    for node in removed:
        ancestor = parents[node]
        while ancestor not in retained:
            ancestor = parents[ancestor]
        groups.setdefault(ancestor, []).append(node)
    info.update({"input_topology_preserved": not removed, "merged_bifurcation_swc_ids": removed,
                 "merge_groups": [{"retained_source_bifurcation": node, "merged_source_bifurcations": values} for node, values in groups.items()],
                 "all_original_terminals_preserved": True,
                 "original_branch_count": original_topo.number_of_edges(), "native_branch_count": native_topo.number_of_edges(),
                 "sampling_note": "original-count follows each surviving original source path; short merged trunks may appear in more than one path, so total output node count can change."})
    return native_topo, paths, info


def signatures(graph: nx.DiGraph) -> dict:
    """Rooted topology signatures with source-labeled tips, ignoring degree-two nodes."""
    if not nx.is_arborescence(graph):
        raise ValueError("VascularMD model is not an outward tree")
    result = {}
    for node in reversed(list(nx.topological_sort(graph))):
        children = list(graph.successors(node))
        if not children:
            source_id = graph.nodes[node].get("source_swc_id")
            if source_id is None:
                raise ValueError(f"Unmapped VascularMD terminal node {node}")
            token = f"tip:{source_id}"
        elif len(children) == 1:
            token = result[children[0]]
        else:
            token = hashlib.sha256("/".join(sorted(result[c] for c in children)).encode()).hexdigest()
        if graph.in_degree(node) == 0:
            token = "root:" + token
        result[node] = token
    return result


@dataclass
class Piece:
    edge: tuple[int, int]
    spline: object
    reverse: bool
    shape: object | None = None

    @property
    def length(self):
        return float(self.spline.length())

    def evaluate(self, lengths: np.ndarray) -> np.ndarray:
        native_lengths = self.length - lengths if self.reverse else lengths
        times = self.spline.length_to_time(np.clip(native_lengths, 0, self.length).tolist())
        values = self.spline.point(times, radius=True)
        if self.shape is not None:
            # Native projection + radius evaluation, NOT another fitted radius model.
            rt = [self.shape.project_point_to_centerline(point[:3]) for point in values]
            values[:, 3] = self.shape.radius(rt)
        if not np.isfinite(values).all() or np.any(values[:, 3] <= 0):
            raise ValueError(f"Nonpositive/nonfinite VascularMD result on edge {self.edge}")
        return values


def model_piece(graph, edge) -> Piece:
    spline = graph.edges[edge].get("spline")
    if spline is None:
        raise ValueError(f"Unmodeled edge {edge}; no fallback geometry will be generated")
    shape = None
    bif_node = next((n for n in edge if graph.nodes[n]["type"] == "bif"), None)
    if bif_node is not None:
        bif = graph.nodes[bif_node]["bifurcation"]
        if bif is None:
            raise ValueError(f"Missing Nfurcation at model node {bif_node}")
        if edge[1] == bif_node:
            shape = bif.get_spl()[0]
        else:
            # Native outlet trajectory is oriented from its sep node TO the bif.
            sep = np.asarray(graph.nodes[edge[1]]["coords"])[:3]
            index = int(np.argmin([np.linalg.norm(section[0][:3] - sep) for section in bif.get_endsec()[1:]]))
            if np.linalg.norm(bif.get_endsec()[index + 1][0][:3] - sep) > 1e-6:
                raise ValueError(f"Cannot map native outlet shape for edge {edge}")
            shape = bif.get_spl()[index]
    start = np.asarray(graph.nodes[edge[0]]["coords"])[:3]
    reverse = np.linalg.norm(spline.point(1.0) - start) < np.linalg.norm(spline.point(0.0) - start)
    return Piece(edge, spline, bool(reverse), shape)


def sample_model(tree, original_topo, branches, source: Source, mode: str, spacing: float | None):
    model = tree.get_model_graph()
    if model is None:
        raise ValueError("VascularMD returned no model")
    raw_sig, model_sig = signatures(original_topo), signatures(model)
    source_by_sig = {raw_sig[n]: original_topo.nodes[n]["source_swc_id"] for n in original_topo}
    significant = [n for n in model if model.in_degree(n) != 1 or model.out_degree(n) != 1]
    if set(model_sig[n] for n in significant) != set(source_by_sig):
        missing = [source_by_sig[key] for key in source_by_sig if key not in {model_sig[n] for n in significant}]
        raise TopologyChanged({"stage": "topology_validation", "affected_source_swc_nodes_and_ancestors": missing,
                               "original_bifurcation_count": sum(original_topo.out_degree(n) > 1 for n in original_topo),
                               "model_bifurcation_count": sum(model.out_degree(n) > 1 for n in model),
                               "original_terminal_count": sum(original_topo.out_degree(n) == 0 for n in original_topo),
                               "model_terminal_count": sum(model.out_degree(n) == 0 for n in model)})
    original_id = {n: int(source_by_sig[model_sig[n]]) for n in significant}
    paths = {}
    for start in significant:
        for child in model.successors(start):
            path = [start, child]
            while path[-1] not in original_id:
                path.append(next(iter(model.successors(path[-1]))))
            paths[(original_id[start], original_id[path[-1]])] = list(zip(path[:-1], path[1:]))
    if set(paths) != {b.key for b in branches}:
        raise ValueError("Original branch connectivity was not preserved by the native model")
    graph = nx.DiGraph()
    # Use a single shared junction node. Its radius is NOT a daughter average.
    for node in significant:
        incoming = list(model.in_edges(node))
        edge = incoming[0] if incoming else next(iter(model.out_edges(node)))
        piece = model_piece(model, edge)
        value = piece.evaluate(np.asarray([piece.length if incoming else 0.0]))[0]
        if model.nodes[node]["type"] == "bif":
            value[:3] = model.nodes[node]["bifurcation"].get_X()
        source_id = original_id[node]
        graph.add_node(source_id, coords=value, swc_type=source.graph.nodes[source_id]["swc_type"])
    next_id = max(source.graph) + 1
    for branch in branches:
        pieces = [model_piece(model, edge) for edge in paths[branch.key]]
        lengths = np.asarray([piece.length for piece in pieces])
        if np.any(lengths <= 0) or not np.isfinite(lengths).all():
            raise ValueError(f"Invalid model arc length on original branch {branch.key}")
        boundaries = np.r_[0, np.cumsum(lengths)]
        count = len(branch.raw) if mode == "original-count" else max(2, int(np.ceil(boundaries[-1] / spacing)) + 1)
        arc = np.linspace(0, boundaries[-1], count)
        piece_ids = np.minimum(np.searchsorted(boundaries[1:], arc, side="right"), len(pieces) - 1)
        values, near = np.empty((count, 4)), np.zeros(count, dtype=bool)
        for i, piece in enumerate(pieces):
            mask = piece_ids == i
            if mask.any():
                values[mask] = piece.evaluate(arc[mask] - boundaries[i])
                near[mask] = piece.shape is not None
        values[0], values[-1] = graph.nodes[branch.key[0]]["coords"], graph.nodes[branch.key[1]]["coords"]
        near[0], near[-1] = branch.start_bif, branch.end_bif
        # Type inheritance by nearest normalized arc position INSIDE the original branch.
        raw_arc = np.r_[0, np.cumsum(np.linalg.norm(np.diff(branch.raw[:, :3], axis=0), axis=1))]
        fractions = raw_arc / raw_arc[-1]
        target = arc / arc[-1]
        right = np.searchsorted(fractions, target).clip(0, len(fractions) - 1)
        left = np.maximum(right - 1, 0)
        nearest = np.where(abs(fractions[left] - target) <= abs(fractions[right] - target), left, right)
        if len(branch.raw) > 2:
            nearest[1:-1] = nearest[1:-1].clip(1, len(branch.raw) - 2)
        types = branch.types[nearest]
        types[0], types[-1] = branch.types[0], branch.types[-1]
        parent = branch.key[0]
        for i in range(1, count - 1):
            graph.add_node(next_id, coords=values[i], swc_type=int(types[i]))
            graph.add_edge(parent, next_id)
            parent, next_id = next_id, next_id + 1
        graph.add_edge(parent, branch.key[1])
        branch.smooth, branch.arc, branch.near_bif, branch.sampled_types, branch.pieces = values, arc, near, types, pieces
    validate_tree(graph)
    # Rooted, source-labeled topology already matched above; counts provide an extra guard.
    if sum(source.graph.out_degree(n) == 0 for n in source.graph) != sum(graph.out_degree(n) == 0 for n in graph):
        raise ValueError("Terminal branch lost during export")
    return graph
