"""Auditable sample-level cuts of frozen MeVO candidates; no semantic queries.

Coordinates, radii, types and directed edges always come from the raw SWC.
Only membership changes. See RULES for the physical-window discretization.
"""
from __future__ import annotations

from contextlib import ExitStack, contextmanager
import csv
from pathlib import Path
from unittest.mock import patch

import networkx as nx
import numpy as np

from .boundary_review_export import component_roots
from .swc_export import read_source, validate_tree, write_swc
from .topbrain_qc import sha256

MIN_DIAMETER_MM = 0.75
MIN_RADIUS_MM = 0.375
PROXIMAL_WINDOW_MM = 5.0
DISTAL_WINDOW_MM = 3.0
REBOUND_WINDOW_MM = 5.0
RULES = {
    'minimum_diameter_mm': MIN_DIAMETER_MM, 'minimum_radius_mm': MIN_RADIUS_MM,
    'radius_source': 'raw BG001.CNG.swc; D = 2*r; no interpolation or smoothing',
    'valid_prediction': 'finite probabilities/support/agreement; valid_donor_count >= 3; cached label != UNKNOWN',
    'proximal': 'pMeVO >= .60; closed [s,s+min(5,remaining)] mm window; >=3 valid original samples; >=80% pMeVO>=.60; no finite pM1>=.80 in window',
    'reversion': 'reject if any consecutive original samples with finite pM1>=.60 span >=3 mm to branch end, including cached UNKNOWN; nonfinite pM1 interrupts the run',
    'distal': 'minimal original-sample window spanning >=3 mm (ends at first sample reaching/passing 3 mm); first sample D<.75, >=2 low samples, >=80% low samples; no terminal extrapolation',
    'distal_discretization': 'No radius interpolation: actual sample span may exceed 3 mm and is recorded. 80% is sample-count based, as requested, not length-weighted.',
    'rebound': 'consecutive original samples with D>=.75 spanning >=3 mm inside closed next 5 mm; reject that candidate',
    'bifurcation': 'evaluate full root-to-tip paths; a confirmed cut removes all descendants, even other daughter subtrees; a rebound on ANY eligible continuation vetoes that shared cut',
    'branch_start': 'a sustained cut at a branch junction removes its outgoing subtree; the junction remains an endpoint; root D<.75 excludes the whole derived component',
    'meaning': 'Operational MeVO/manufacturing filter; NOT anatomical M3 endpoint or M3/M4 annotation. Tolerated dips and cut endpoints may have D<.75.',
}


@contextmanager
def frozen_upstream():
    """Allow native modeling of derived files, disable all upstream inference."""
    names = [
        'similarity_registration.register', 'semantic_ensemble.register',
        'semantic_ensemble.select_donors', 'semantic_ensemble.nearest_baseline',
        'semantic_ensemble.nearest_support', 'brava_branch_labels.aggregate_branches',
        'brava_branch_labels.select_major_branches', 'nn_support_gate.freeze_pilot',
        'nn_support_gate.calibrate',
    ]
    with ExitStack() as stack:
        mocks = [stack.enter_context(patch('vascular_processing.' + name,
                 side_effect=AssertionError('FORBIDDEN_UPSTREAM: ' + name))) for name in names]
        yield names
        if any(m.call_count for m in mocks):
            raise AssertionError('A forbidden upstream function was invoked')


def table(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields or list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def arc(nodes, coordinates):
    xyz = np.asarray([coordinates[n][:3] for n in nodes])
    return np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(xyz, axis=0), axis=1))]


def continuous_run(s, predicate, minimum=3.0):
    """A single point has zero length; no interpolation of a discrete label."""
    start = None
    for distance, value in zip(s, predicate):
        if not value:
            start = None
        else:
            if start is None:
                start = distance
            if distance - start >= minimum - 1e-9:
                return True
    return False


def proximal_start(nodes, coordinates, predictions):
    s = arc(nodes, coordinates)
    p1 = np.array([float(predictions[n]['p_M1']) for n in nodes])
    pm = np.array([float(predictions[n]['p_MeVO']) for n in nodes])
    valid = np.array([
        int(predictions[n]['valid_donor_count']) >= 3 and predictions[n]['label'] != 'UNKNOWN'
        and np.isfinite([p1[i], pm[i], float(predictions[n]['agreement']),
                         float(predictions[n]['normalized_support_distance'])]).all()
        for i, n in enumerate(nodes)])
    attempts = []
    for i, n in enumerate(nodes):
        if not valid[i] or pm[i] < .60:
            continue
        length = min(PROXIMAL_WINDOW_MM, float(s[-1] - s[i]))
        window = (s >= s[i]) & (s <= s[i] + length + 1e-9)
        count = int(np.count_nonzero(window & valid))
        fraction = float(np.mean(pm[window & valid] >= .60)) if count else 0.0
        reversal = continuous_run(s[i:], np.isfinite(p1[i:]) & (p1[i:] >= .60))
        strong = bool(np.any(p1[window] >= .80))
        status = ('PROXIMAL_REFINEMENT_INSUFFICIENT_SAMPLES' if count < 3 else
                  'M1_REVERSION' if reversal or strong else
                  'UNSTABLE_MEVO_WINDOW' if fraction < .80 else 'SUSTAINED_MEVO_START')
        result = dict(node=n, index=i, arc_mm=float(s[i]), window_length_mm=length,
                      window_valid_samples=count, window_MeVO_fraction=fraction,
                      reversion_detected=reversal or strong, status=status)
        attempts.append(result)
        if status == 'SUSTAINED_MEVO_START':
            return result, attempts
    status = ('PROXIMAL_REFINEMENT_INSUFFICIENT_SAMPLES' if attempts and all(
        a['status'] == 'PROXIMAL_REFINEMENT_INSUFFICIENT_SAMPLES' for a in attempts)
        else 'NO_SUSTAINED_MEVO_START')
    return dict(node=None, index=None, arc_mm=None, window_length_mm=None,
                window_valid_samples=0, window_MeVO_fraction=None,
                reversion_detected=any(a['reversion_detected'] for a in attempts), status=status), attempts


def branch_edges(cache):
    return {(a, b): key for key, branch in cache.branches.items()
            for a, b in zip(branch['node_ids'], branch['node_ids'][1:])}


def original_component(cache, component, source):
    keys = {key for key, c in cache.branch_component.items() if c == component}
    nodes = {n for key in keys for n in cache.branches[key]['node_ids']}
    return source.graph.subgraph(nodes).copy()


def refine_component(cache, component, source):
    """Treat each entrance independently; merge only genuinely shared kept nodes."""
    original = original_component(cache, component, source)
    roots = component_roots(cache, component)
    original_root = validate_tree(original)
    retained = set()
    rows, audit = [], []

    def search(key, boundary_branch, boundary_node, extend=False):
        nodes = cache.branches[key]['node_ids']
        cut, attempts = proximal_start(nodes, cache.coordinates, cache.points)
        for item in attempts:
            audit.append(dict(component=component, branch=key, **item))
        n = cut['node']
        p = cache.points[n] if n is not None else {}
        xyzr = cache.coordinates[n] if n is not None else [None] * 4
        # Signed along native flow: an upstream extension has negative shift.
        shift = cut['arc_mm']
        if n is not None and extend:
            shift -= float(arc(nodes, cache.coordinates)[-1])
        elif n is not None:
            shift = nx.shortest_path_length(cache.node_graph, boundary_node, n, weight='weight')
        rows.append(dict(component=component, original_boundary_node=boundary_node,
            original_boundary_branch=boundary_branch, refined_cut_branch=key,
            refined_cut_original_swc_id=n, refined_cut_x=xyzr[0], refined_cut_y=xyzr[1],
            refined_cut_z=xyzr[2], distance_from_original_boundary_mm=shift,
            p_M1=float(p['p_M1']) if p else None, p_MeVO=float(p['p_MeVO']) if p else None,
            agreement=float(p['agreement']) if p else None,
            normalized_support_distance=float(p['normalized_support_distance']) if p else None,
            window_length_mm=cut['window_length_mm'], window_MeVO_fraction=cut['window_MeVO_fraction'],
            reversion_detected=cut['reversion_detected'], status=cut['status'],
            window_valid_samples=cut['window_valid_samples'], radius_mm=xyzr[3],
            diameter_mm=2 * xyzr[3] if n is not None else None))
        if n is not None:
            if extend:
                retained.update(original)
                retained.update(nodes[cut['index']:])
            else:
                retained.update({n} | nx.descendants(original, n))
            return True
        return False

    # Only the explicitly requested part03 UNKNOWN parent is eligible for extension.
    extended = False
    if component == 3:
        if {cache.branches[k]['parent_branch'] for k in roots} != {'94_112'}:
            raise ValueError('Frozen part03 no longer has expected UNKNOWN parent 94_112')
        extended = search('94_112', ';'.join(roots), original_root, extend=True)
    if not extended:
        def descend(key, boundary_key, boundary_node):
            if search(key, boundary_key, boundary_node):
                return
            for child in cache.branch_graph.successors(key):
                if cache.branch_component.get(child) == component:
                    descend(child, boundary_key, boundary_node)
        for key in roots:
            descend(key, key, cache.branches[key]['node_ids'][0])
    forest = source.graph.subgraph(retained).copy()
    groups = sorted(nx.weakly_connected_components(forest), key=lambda ns: min(ns))
    graphs = [forest.subgraph(ns).copy() for ns in groups]
    for graph in graphs:
        validate_tree(graph)
    return original, graphs, rows, audit


def distal_candidates(path, coordinates):
    """Sample-supported >=3 mm windows with a separate exact 5 mm rebound check."""
    s = arc(path, coordinates)
    diameter = np.array([2 * coordinates[n][3] for n in path])
    low = diameter < MIN_DIAMETER_MM
    candidates = []
    for i, n in enumerate(path):
        if not low[i]:
            continue
        j = int(np.searchsorted(s, s[i] + DISTAL_WINDOW_MM - 1e-9))
        if j >= len(path):
            continue
        window = low[i:j + 1]
        fraction = float(np.mean(window))
        if np.count_nonzero(window) < 2 or fraction < .80:
            continue
        end = int(np.searchsorted(s, s[i] + REBOUND_WINDOW_MM + 1e-9, side='right'))
        rebound = continuous_run(s[i:end], ~low[i:end])
        candidates.append(dict(node=n, window_span_mm=float(s[j] - s[i]),
            window_samples=len(window), window_low_fraction=fraction,
            status='DIAMETER_REBOUND' if rebound else 'DISTAL_DIAMETER_CUT'))
    return candidates


def graph_length(graph, coordinates):
    return float(sum(np.linalg.norm(coordinates[b][:3] - coordinates[a][:3]) for a, b in graph.edges))


def prune_diameter(graph, coordinates, edge_branch, component, derived_id):
    root = validate_tree(graph)
    rows = []
    length = graph_length(graph, coordinates)
    if 2 * coordinates[root][3] < MIN_DIAMETER_MM:
        rows.append(dict(component=component, derived_component=derived_id,
            original_branch=';'.join(sorted({edge_branch[e] for e in graph.edges if root in e})),
            cut_original_swc_id=root, cut_distance_from_refined_root_mm=0.0,
            radius_mm=float(coordinates[root][3]), diameter_mm=float(2 * coordinates[root][3]),
            removed_downstream_nodes=len(graph), removed_downstream_branches=len({edge_branch[e] for e in graph.edges}),
            removed_length_mm=length, reason='COMPONENT_BELOW_DIAMETER_THRESHOLD',
            window_span_mm=0.0, window_samples=1, window_low_fraction=1.0))
        return nx.DiGraph(), rows, []
    by_node = {}
    for tip in [n for n in graph if graph.out_degree(n) == 0]:
        for item in distal_candidates(nx.shortest_path(graph, root, tip), coordinates):
            by_node.setdefault(item['node'], []).append(dict(terminal=tip, **item))
    distances = nx.single_source_dijkstra_path_length(graph, root,
        weight=lambda a, b, _: float(np.linalg.norm(coordinates[b][:3] - coordinates[a][:3])))
    kept = graph.copy()
    audit = []
    for n in sorted(by_node, key=lambda n: (distances[n], n)):
        items = by_node[n]
        rebound = any(item['status'] == 'DIAMETER_REBOUND' for item in items)
        info = max(items, key=lambda a: (a['status'] == 'DIAMETER_REBOUND', a['window_span_mm']))
        eligible = n in kept and kept.out_degree(n) > 0
        audit.append(dict(component=component, derived_component=derived_id, node=n,
            status='ALREADY_REMOVED_BY_PROXIMAL_CUT' if not eligible else
                   'DIAMETER_REBOUND' if rebound else 'DISTAL_DIAMETER_CUT', path_evidence=items))
        if not eligible:
            continue
        removed = set() if rebound else nx.descendants(kept, n)
        removed_edges = {e for e in kept.edges if e[1] in removed}
        incoming = list(graph.in_edges(n))
        branch = edge_branch[incoming[0]] if incoming else ''
        # At a topology junction the incoming branch is only the endpoint;
        # identify all outgoing affected branches in the audit row.
        outgoing_keys = {edge_branch[e] for e in graph.out_edges(n)}
        if len(outgoing_keys) > 1:
            branch = ';'.join(sorted(outgoing_keys))
        rows.append(dict(component=component, derived_component=derived_id,
            original_branch=branch, cut_original_swc_id=n,
            cut_distance_from_refined_root_mm=float(distances[n]),
            radius_mm=float(coordinates[n][3]), diameter_mm=float(2 * coordinates[n][3]),
            removed_downstream_nodes=len(removed),
            removed_downstream_branches=len({edge_branch[e] for e in removed_edges}),
            removed_length_mm=float(sum(np.linalg.norm(coordinates[b][:3] - coordinates[a][:3]) for a, b in removed_edges)),
            reason='DIAMETER_REBOUND' if rebound else 'DISTAL_DIAMETER_CUT',
            window_span_mm=info['window_span_mm'], window_samples=info['window_samples'],
            window_low_fraction=info['window_low_fraction']))
        kept.remove_nodes_from(removed)
    validate_tree(kept)
    return kept, rows, audit


def statistics(graph, coordinates, edge_branch, full_branch_lengths):
    radius = np.array([coordinates[n][3] for n in graph])
    keys = {edge_branch[e] for e in graph.edges}
    equivalent = sum(float(np.linalg.norm(coordinates[b][:3] - coordinates[a][:3])) /
                     full_branch_lengths[edge_branch[a, b]] for a, b in graph.edges)
    return dict(node_count=len(graph), length_mm=graph_length(graph, coordinates),
        branch_count=len(keys), branch_equivalent_count=float(equivalent),
        terminal_count=sum(graph.out_degree(n) == 0 for n in graph),
        min_radius_mm=float(radius.min()) if len(radius) else None,
        median_radius_mm=float(np.median(radius)) if len(radius) else None,
        max_radius_mm=float(radius.max()) if len(radius) else None,
        min_diameter_mm=float(2 * radius.min()) if len(radius) else None,
        median_diameter_mm=float(2 * np.median(radius)) if len(radius) else None,
        max_diameter_mm=float(2 * radius.max()) if len(radius) else None,
        retained_subthreshold_samples=int(np.count_nonzero(2 * radius < MIN_DIAMETER_MM)))


def export_subset(path, graph, source, edge_branch, mapping_dir, *, context_nodes=()):
    validate_tree(graph)
    for n in graph:
        np.testing.assert_array_equal(graph.nodes[n]['coords'], source.graph.nodes[n]['coords'])
        if graph.nodes[n]['swc_type'] != source.graph.nodes[n]['swc_type']:
            raise ValueError('Original TYPE changed')
    if not set(graph.edges) <= set(source.graph.edges):
        raise ValueError('Derived graph introduced non-native edges')
    path.parent.mkdir(parents=True, exist_ok=True)
    mapping = write_swc(path, graph, source.path, 'not-fitted', 'original-samples', header_lines=[
        'Frozen NN point-level MeVO refinement; original BraVa subset, not a smoothed SWC',
        f'Raw source SHA256: {sha256(source.path)}', 'XYZ/radius in mm; diameter=2*radius',
        'Raw xyz/radius/TYPE and directed edges unchanged; IDs reindexed; root parent=-1',
        'Diameter filtering is operational, NOT an anatomical M3/M4 boundary'])
    readback = read_source(path)
    if set(readback.graph.edges) != {(mapping[a], mapping[b]) for a, b in graph.edges}:
        raise ValueError('SWC round-trip edges differ')
    rows = []
    for n, new in mapping.items():
        np.testing.assert_array_equal(readback.graph.nodes[new]['coords'], source.graph.nodes[n]['coords'])
        if readback.graph.nodes[new]['swc_type'] != source.graph.nodes[n]['swc_type']:
            raise ValueError('SWC round-trip TYPE differs')
        parents = list(graph.predecessors(n))
        original_parents = list(source.graph.predecessors(n))
        rows.append(dict(roi_node_id=new, original_swc_id=n, roi_parent_id=mapping[parents[0]] if parents else -1,
            original_parent_id=original_parents[0] if original_parents else -1,
            radius_mm=float(graph.nodes[n]['coords'][3]), diameter_mm=float(2 * graph.nodes[n]['coords'][3]),
            is_context=n in context_nodes))
    table(mapping_dir / (path.stem + '_nodes.csv'), rows)
    table(mapping_dir / (path.stem + '_edges.csv'), [dict(roi_parent_id=mapping[a], roi_child_id=mapping[b],
        original_parent_id=a, original_child_id=b, original_branch_id=edge_branch[a, b],
        is_context=a in context_nodes or b in context_nodes) for a, b in graph.edges])
    return dict(path=str(path.resolve()), sha256=sha256(path), node_count=len(graph),
                root_original_swc_id=validate_tree(graph), single_inlet_stem=graph.out_degree(validate_tree(graph)) == 1,
                topology_qc=True, geometry_exact=True, radius_exact=True, type_exact=True,
                round_trip=True, context_nodes=sorted(context_nodes))
