"""Exact source graph and explicit mappings to the installed VascularMD topology."""
from dataclasses import dataclass
from pathlib import Path

import networkx as nx
import numpy as np

from third_party.vascularmd.ArterialTree import ArterialTree
from .swc_export import Source, read_source

ROLE_CODES = {"reg": 0, "end": 1, "bif": 2, "sink": 3, "sep": 4}


@dataclass
class ExactGraph:
    source: Source
    tree: ArterialTree
    locations: dict
    edge_locations: dict
    native_valid: bool

    @property
    def graph(self):
        return self.source.graph

    @property
    def topo(self):
        return self.tree.get_topo_graph()


def load_exact_graph(path: Path) -> ExactGraph:
    source = read_source(path, allow_forest=True)
    full = source.graph.copy()
    for node in full:
        full.nodes[node]["coords"] = full.nodes[node]["coords"].copy()
        full.nodes[node]["original_swc_id"] = node
    for edge in full.edges:
        full.edges[edge]["coords"] = np.empty((0, 4))
    tree = ArterialTree(path.stem, "BraVa")
    tree.set_full_graph(full)  # No automatic_resampling or model_to_full: preserve IDs/XYZR.
    full = tree.get_full_graph()
    if set(full) != set(source.graph) or set(full.edges) != set(source.graph.edges):
        raise ValueError("Native full graph changed original nodes or directed edges")
    topo = tree.get_topo_graph()
    locations = {n: {"original_swc_id": n, "full_graph_node_id": n, "topo_node_id": -1,
                     "topo_edge_start": -1, "topo_edge_end": -1, "topo_edge_offset": -1,
                     "topology_type_code": ROLE_CODES["reg"]} for n in full}
    edge_locations = {}
    mapped = set()
    for node, data in topo.nodes(data=True):
        original = int(data["full_id"])
        np.testing.assert_array_equal(data["coords"], source.graph.nodes[original]["coords"])
        locations[original].update(topo_node_id=int(node), topology_type_code=ROLE_CODES[data["type"]])
        mapped.add(original)
    for a, b, data in topo.edges(data=True):
        ids = [int(topo.nodes[a]["full_id"]), *map(int, data["full_id"]), int(topo.nodes[b]["full_id"])]
        for i, node in enumerate(ids[1:-1]):
            if node in mapped:
                raise ValueError(f"Ambiguous original/full/topo mapping at {node}")
            np.testing.assert_array_equal(data["coords"][i], source.graph.nodes[node]["coords"])
            locations[node].update(topo_edge_start=int(a), topo_edge_end=int(b), topo_edge_offset=i)
            mapped.add(node)
        for edge in zip(ids[:-1], ids[1:]):
            edge_locations[edge] = (int(a), int(b))
    if mapped != set(source.graph) or set(edge_locations) != set(source.graph.edges):
        raise ValueError("Native topology mapping does not exactly cover source SWC")
    for i, component in enumerate(sorted(nx.weakly_connected_components(full), key=min)):
        for node in component:
            locations[node]["component_id"] = i
    return ExactGraph(source, tree, locations, edge_locations, tree.check_full_graph())


def graph_counts(graph):
    significant = [n for n in graph if graph.in_degree(n) != 1 or graph.out_degree(n) != 1]
    return {"node_count": len(graph), "edge_count": graph.number_of_edges(),
            "component_count": nx.number_weakly_connected_components(graph),
            "root_count": sum(graph.in_degree(n) == 0 for n in graph),
            "terminal_count": sum(graph.out_degree(n) == 0 for n in graph),
            "bifurcation_count": sum(graph.out_degree(n) > 1 for n in graph),
            "branch_count": sum(graph.out_degree(n) for n in significant)}
