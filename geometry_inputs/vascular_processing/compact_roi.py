"""Single frozen MeVO component -> bounded native subtree; never recompute semantics."""
from __future__ import annotations

import copy
import json
from pathlib import Path

import networkx as nx
import numpy as np

from .boundary_review_export import read_csv
from .manufacturing_roi import (geometry_stats, topology_paths, compensate_radius,
                                coordinates, export_manufacturing)
from .mevo_refinement import graph_length, branch_edges, table
from .swc_export import validate_tree

MODES = ('MINI', 'BALANCED', 'RICH')


def distances(graph, root):
    return nx.single_source_dijkstra_path_length(graph, root, weight=lambda a, b, _: float(
        np.linalg.norm(graph.nodes[a]['coords'][:3] - graph.nodes[b]['coords'][:3])))


def metrics(graph):
    stats = geometry_stats(graph)
    root = validate_tree(graph)
    bif = {n for n in graph if graph.out_degree(n) > 1}
    stats.update(branch_count=len(topology_paths(graph)), bifurcation_count=len(bif),
                 root_original_id=root, maximum_inlet_to_outlet_mm=max(distances(graph, root).values()),
                 maximum_bifurcation_depth=max(sum(n in bif for n in nx.shortest_path(graph, root, tip))
                    for tip in graph if graph.out_degree(tip) == 0))
    return stats


def load_sources(directory, raw):
    directory = Path(directory)
    summary = json.loads((directory / 'refinement_summary.json').read_text())
    sources = []
    for item in summary['exports']:
        path = Path(item['path'])
        rows = read_csv(directory / 'mappings' / (path.stem + '_nodes.csv'))
        nodes = {int(r['original_swc_id']) for r in rows}
        graph = raw.graph.subgraph(nodes).copy()
        root = validate_tree(graph)
        sources.append(dict(component=item['component'], derived_component=item['derived_component'],
            level=item['level'], path=str(path), source_sha256=item['sha256'],
            root=root, graph=graph))
    return sources


def effective_limits(config, mode):
    c = config['compact_roi']; m = config['compact_phantom'][mode.lower()]
    return dict(max_branches=min(c['max_topology_branches'], m['max_branches']),
        max_bifurcations=min(c['hard_max_bifurcations'], m['max_bifurcations']),
        max_outlets=min(c['hard_max_outlets'], m['max_outlets']),
        max_path_mm=min(c['hard_max_path_length_mm'], m['max_path_mm']),
        max_longest_dimension_mm=min(c['hard_max_longest_dimension_mm'], m['max_longest_dimension_mm']))


def branch_score(path, graph, semantic_support, config, max_diameter, max_length):
    xyzr = np.array([graph.nodes[n]['coords'] for n in path])
    lengths = np.linalg.norm(np.diff(xyzr[:, :3], axis=0), axis=1)
    # Probability comes directly from frozen branch rows, not from a new vote/aggregation.
    support = float(semantic_support.get((path[0], path[1]), 0.))
    diameter = float(np.median(2 * xyzr[:, 3])); length = float(lengths.sum())
    w = config['compact_roi']['daughter_score_weights']
    score = w['semantic'] * support + w['diameter'] * diameter / max(max_diameter, 1e-12) + w['length'] * length / max(max_length, 1e-12)
    return dict(score=float(score), frozen_semantic_support=support,
                original_median_diameter_mm=diameter, original_length_mm=length)


def compact_subtree(graph, config, mode, semantic_support):
    """Expand complete topological branches greedily, clip only at native samples.

    Bifurcation depth counts encountered anatomical junctions, even if an earlier
    junction loses a daughter. It is never used to assign M2/M3 labels.
    """
    root = validate_tree(graph); limits = effective_limits(config, mode); c = config['compact_roi']
    dist = distances(graph, root)
    paths = topology_paths(graph)
    by_start = {}
    for path in paths:
        by_start.setdefault(path[0], []).append(path)
    diam = max(float(np.median([2 * graph.nodes[n]['coords'][3] for n in p])) for p in paths)
    longest = max(dist[p[-1]] - dist[p[0]] for p in paths)
    scores = {tuple(p[:2]): branch_score(p, graph, semantic_support, config, diam, longest) for p in paths}
    retained = {root}; stops = {}; decisions = []; frontier = [root]
    expanded_bifurcations = 0

    def discard_below(node, reason):
        for a, b in graph.subgraph(nx.descendants(graph, node) | {node}).edges:
            stops.setdefault((a, b), reason)

    while frontier:
        frontier.sort(key=lambda n: (-max((scores[tuple(p[:2])]['score'] for p in by_start.get(n, [])), default=-1), n))
        start = frontier.pop(0); daughters = by_start.get(start, [])
        if not daughters:
            continue
        is_bif = len(daughters) > 1
        depth = sum(graph.out_degree(n) > 1 for n in nx.shortest_path(graph, root, start))
        if is_bif and depth > c['max_bifurcation_depth']:
            discard_below(start, 'EXCEEDS_BIFURCATION_DEPTH'); continue
        if is_bif and expanded_bifurcations >= limits['max_bifurcations']:
            discard_below(start, 'EXCEEDS_BRANCH_BUDGET'); continue
        ranked = sorted(daughters, key=lambda p: (-scores[tuple(p[:2])]['score'], p[1]))
        chosen = ranked[:c['max_daughters_per_bifurcation']] if is_bif else ranked
        for path in ranked[len(chosen):]:
            stops[path[0], path[1]] = 'LOW_MANUFACTURING_SCORE'
            discard_below(path[1], 'LOW_MANUFACTURING_SCORE')
        proposals = []
        for path in chosen:
            clipped = [n for n in path if dist[n] <= limits['max_path_mm'] + 1e-9]
            if not clipped or clipped[0] != start:
                clipped = [start]
            if len(clipped) < len(path):
                cut = clipped[-1]
                for a, b in zip(path, path[1:]):
                    if dist[b] > limits['max_path_mm'] + 1e-9:
                        stops[a, b] = 'EXCEEDS_PATH_LENGTH'
                discard_below(cut, 'EXCEEDS_PATH_LENGTH')
            if len(clipped) > 1:
                proposals.append(clipped)
        trial = graph.subgraph(retained | {n for p in proposals for n in p}).copy()
        # Keep the branch budget independent of SWC's sampling density.
        if (len(topology_paths(trial)) > limits['max_branches'] or
                sum(trial.out_degree(n) == 0 for n in trial) > limits['max_outlets'] or
                graph_length(trial, coordinates(trial)) > c['hard_max_total_centerline_mm']):
            discard_below(start, 'EXCEEDS_BRANCH_BUDGET'); continue
        retained.update(trial)
        if is_bif and len(proposals) > 1:
            expanded_bifurcations += 1
        for path in proposals:
            decisions.append(dict(first=path[0], last=path[-1], original_ids=path,
                                  **scores[tuple(path[:2])]))
            if path[-1] in by_start:
                frontier.append(path[-1])
    result = graph.subgraph(retained).copy()
    removed = []
    twig = c['terminal_twig']
    for path in topology_paths(result):
        if result.out_degree(path[-1]) or result.out_degree(path[0]) < 2:
            continue
        xyzr = np.array([result.nodes[n]['coords'] for n in path])
        length = float(np.linalg.norm(np.diff(xyzr[:, :3], axis=0), axis=1).sum())
        median = float(np.median(2 * xyzr[:, 3]))
        if length < twig['max_length_mm'] and median < twig['original_diameter_below_mm']:
            removed.append(dict(original_ids=path[1:], length_mm=length, diameter_median_mm=median,
                                reason='THIN_SHORT_TERMINAL_TWIG'))
            for edge in zip(path, path[1:]):
                stops[edge] = 'THIN_SHORT_TERMINAL_TWIG'
    for item in removed:
        result.remove_nodes_from(item['original_ids'])
    validate_tree(result)
    return result, dict(limits=limits, decisions=decisions, stops=stops, removed_twigs=removed)


def add_context(core, raw_tree, other_semantic_nodes, config):
    """Choose the actual upstream prefix nearest 10 mm, bounded by 8--15 mm."""
    root = validate_tree(core); cursor = root; total = 0.; chain = []; choices = []
    c = config['compact_phantom']['proximal_context']
    while raw_tree.in_degree(cursor):
        parent = next(raw_tree.predecessors(cursor))
        if parent in other_semantic_nodes:
            break
        edge = float(np.linalg.norm(raw_tree.nodes[parent]['coords'][:3] - raw_tree.nodes[cursor]['coords'][:3]))
        if total + edge > c['max_mm'] + 1e-9:
            break
        total += edge; chain.append(parent); cursor = parent
        choices.append((total, chain.copy()))
    feasible = [x for x in choices if x[0] >= c['min_mm']]
    selected = min(feasible or choices or [(0., [])], key=lambda x: (abs(x[0] - c['preferred_mm']), x[0]))
    result = raw_tree.subgraph(set(core) | set(selected[1])).copy()
    validate_tree(result)
    return result, dict(length_mm=selected[0], original_ids=selected[1], refined_root=root,
        final_root=validate_tree(result), min_satisfied=selected[0] >= c['min_mm'],
        purpose='Stable inlet only; native upstream samples, no added siblings, no component merging')


def gate(stats, config, mode):
    c = config['compact_roi']; lim = effective_limits(config, mode); reasons = []
    if stats['inlet_count'] != 1 or stats['connected_component_count'] != 1:
        reasons.append('SINGLE_COMPONENT_OR_INLET_REQUIRED')
    if not c['min_bifurcations'] <= stats['bifurcation_count'] <= lim['max_bifurcations']:
        reasons.append('BIFURCATION_BUDGET')
    if not 2 <= stats['outlet_count'] <= lim['max_outlets']:
        reasons.append('OUTLET_BUDGET')
    if stats['branch_count'] > lim['max_branches']:
        reasons.append('BRANCH_BUDGET')
    if stats['maximum_bifurcation_depth'] > c['max_bifurcation_depth']:
        reasons.append('BIFURCATION_DEPTH_BUDGET')
    if stats['maximum_inlet_to_outlet_mm'] > c['max_inlet_to_outlet_mm'] and not c.get('adaptive_path', {}).get('enabled'):
        reasons.append('INLET_OUTLET_LENGTH_BUDGET')
    if stats['centerline_length_mm'] > c['hard_max_total_centerline_mm']:
        reasons.append('TOTAL_LENGTH_BUDGET')
    if stats['bbox']['longest_dimension_mm'] > lim['max_longest_dimension_mm'] or stats['bbox']['diagonal_mm'] > c['hard_max_bbox_diagonal_mm']:
        reasons.append('MAXIMUM_SIZE_BUDGET')
    undersized = stats['bbox']['longest_dimension_mm'] < c['minimum_longest_dimension_mm']
    return dict(structure_passed=not reasons, undersized=undersized, reasons=reasons,
        passed=not reasons and (not undersized or c['allow_undersized_candidate']),
        size_status='BELOW_60_MM_TARGET' if undersized else 'SIZE_TARGET_MET')


def radius_metrics(original, compensated, config):
    nodes = list(original)
    r = np.array([original.nodes[n]['coords'][3] for n in nodes])
    rc = np.array([compensated.nodes[n]['coords'][3] for n in nodes]); ratio = rc / r
    length = floor_length = changed_length = 0.
    floor = config['effective']['radius_floor_mm']
    for a, b in original.edges:
        ds = float(np.linalg.norm(original.nodes[a]['coords'][:3] - original.nodes[b]['coords'][:3]))
        length += ds
        floor_length += ds * sum(float(original.nodes[n]['coords'][3] < floor - 1e-12) for n in (a, b)) / 2
        changed_length += ds * sum(float(compensated.nodes[n]['coords'][3] > original.nodes[n]['coords'][3] + 1e-12) for n in (a, b)) / 2
    return dict(modified_samples=int(np.count_nonzero(rc > r + 1e-12)), sample_fraction=float(np.mean(rc > r + 1e-12)),
        centerline_floor_fraction=floor_length / length if length else 0.,
        centerline_compensated_fraction=changed_length / length if length else 0.,
        fraction_method='Length weighted endpoint-indicator trapezoid; no anatomical interpolation',
        maximum_radius_inflation_factor=float(ratio.max()), sum_added_radius_mm=float((rc-r).sum()),
        large_compensation_ids=[int(n) for n, v in zip(nodes, ratio) if v > config['manufacturing']['max_radius_inflation_factor']],
        warning='LARGE_MANUFACTURING_COMPENSATION' if np.any(ratio > config['manufacturing']['max_radius_inflation_factor']) else None)


def fixed_candidate(source, semantic_nodes, raw, config, mode, edge_support, forbidden):
    graph, audit = compact_subtree(source['graph'], config, mode, edge_support)
    graph, context = add_context(graph, raw.graph, forbidden, config)
    comp, blend = compensate_radius(graph, config['effective']['radius_floor_mm'], config['manufacturing']['radius_blend_length_mm'])
    stats = metrics(graph); radii = radius_metrics(graph, comp, config)
    score = source['component_support'] - config['compact_roi']['candidate_compensation_penalty_weight'] * radii['centerline_floor_fraction']
    return dict(mode=mode, graph=graph, compensated=comp, stats=stats, context=context, radii=radii,
        blend=blend, audit=audit, gate=gate(stats, config, mode), score=score,
        core_nodes=set(graph) & semantic_nodes)


def generate_candidate(source, semantic_nodes, raw, config, mode, edge_support, forbidden):
    adaptive = config['compact_roi'].get('adaptive_path', {})
    if not adaptive.get('enabled'):
        return fixed_candidate(source, semantic_nodes, raw, config, mode, edge_support, forbidden)
    g = source['graph']; root = validate_tree(g); dist = distances(g, root)
    ceiling = config['compact_roi']['hard_max_total_centerline_mm']
    values = sorted(set(v for v in dist.values() if 0 < v <= ceiling))
    # Deterministic sampling of native endpoints plus every topology boundary.
    indices = np.unique(np.linspace(0, len(values)-1, min(len(values), adaptive['maximum_trials']), dtype=int))
    cuts = {values[i] for i in indices}
    cuts.update(dist[p[-1]] for p in topology_paths(g) if dist[p[-1]] <= ceiling)
    targets = config['compact_phantom'][mode.lower()]
    desired_bif = min(targets['max_bifurcations'], (targets['max_branches']-1)//2,
                      sum(g.out_degree(n)>1 for n in g))
    rows = []; best = None; best_key = None
    w = adaptive['score_weights']
    for cutoff in sorted(cuts):
        trial = copy.deepcopy(config)
        trial['compact_roi']['hard_max_path_length_mm'] = cutoff
        trial['compact_phantom'][mode.lower()]['max_path_mm'] = cutoff
        candidate = fixed_candidate(source, semantic_nodes, raw, trial, mode, edge_support, forbidden)
        s = candidate['stats']; size = s['bbox']['longest_dimension_mm']; length = s['centerline_length_mm']
        score = (w['size_deficit'] * max(0., 1-size/targets['target_longest_dimension_mm'])
            + w['structure_deficit'] * max(0, desired_bif-s['bifurcation_count']) / max(desired_bif,1)
            + w['length_excess'] * max(0., length/config['compact_roi']['preferred_total_centerline_mm']['max']-1)
            + w['compensation'] * candidate['radii']['centerline_floor_fraction'])
        row = dict(cutoff_mm=cutoff, valid=candidate['gate']['structure_passed'], cost=score,
            branches=s['branch_count'], bifurcations=s['bifurcation_count'], outlets=s['outlet_count'],
            length_mm=length, longest_dimension_mm=size, reasons=candidate['gate']['reasons'])
        rows.append(row)
        key = (not row['valid'], score, length, cutoff)
        if best_key is None or key < best_key:
            best = candidate; best_key = key
    best['adaptive_search'] = dict(enabled=True, selected_cutoff_mm=best['audit']['limits']['max_path_mm'],
        trials=rows, target_longest_dimension_mm=targets['target_longest_dimension_mm'],
        former_fixed_path_reference_mm=targets['max_path_mm'],
        former_inlet_to_outlet_reference_mm=config['compact_roi']['max_inlet_to_outlet_mm'],
        selection='Valid topology/size budgets, then size + bifurcation + total-length + compensation cost; shortest tie',
        authorization=adaptive['authorization'])
    return best


def choose_source(sources, semantic, raw, scores, cache, config):
    means = {int(r['component_id'].removeprefix('roi_part')): float(r['MeVO_mean_probability']) for r in scores}
    ebranch = branch_edges(cache)
    support = {edge: float(cache.branches[b]['p_MeVO']) for edge, b in ebranch.items()}
    c = config['compact_roi']; priority = c['component_priority']; levels = c['source_level_priority']
    trials = []
    for source in sources:
        component = source['component']; source['component_support'] = means[component]
        forbidden = set().union(*(ns for k, ns in semantic.items() if k != component)) | (semantic[component] - set(source['graph']))
        family = {mode: generate_candidate(source, semantic[component], raw, config, mode, support, forbidden) for mode in MODES}
        candidate = family['BALANCED']
        reason = candidate['gate']['reasons'].copy()
        for mode in ('MINI', 'RICH'):
            if not family[mode]['gate']['structure_passed']:
                reason.append(mode + '_MINIMUM_STRUCTURE_UNAVAILABLE')
        if means[component] < c['minimum_semantic_support']:
            reason.append('LOW_FROZEN_SEMANTIC_SUPPORT')
        valid = not reason
        item = dict(component=component, derived_component=source['derived_component'], source_level=source['level'],
            source_path=source['path'], root=source['root'], structure_eligible=valid, reasons=reason,
            size_gate=candidate['gate'], balanced_stats=candidate['stats'],
            family_structure={mode:x['gate'] for mode,x in family.items()}, candidate_score=candidate['score'])
        trials.append((source, item, forbidden))
    eligible = [t for t in trials if t[1]['structure_eligible']]
    if not eligible:
        raise ValueError('NO_STRUCTURALLY_VALID_COMPACT_SOURCE: ' + json.dumps([t[1] for t in trials]))
    # Explicit user priority is configuration, never a filename/node ID in the algorithm.
    selected = min(eligible, key=lambda t: (priority.index(t[0]['component']) if t[0]['component'] in priority else len(priority),
        -t[1]['candidate_score'], levels.index(t[0]['level']), t[0]['root']))
    source, _, forbidden = selected
    candidates = {mode: generate_candidate(source, semantic[source['component']], raw, config, mode, support, forbidden) for mode in MODES}
    assert len({source['component'] for _ in candidates}) == 1
    return source, candidates, [t[1] for t in trials]


def branch_audit(semantic_graph, selected, source, edge_branch, stops):
    rows = []
    root = source['root']
    grouped = {}
    for edge in semantic_graph.edges:
        grouped.setdefault(edge_branch[edge], []).append(edge)
    for branch, edges in sorted(grouped.items()):
        kept = [e for e in edges if selected.has_edge(*e)]; missing = [e for e in edges if e not in kept]
        reasons = set()
        for a, b in missing:
            if (a, b) in stops:
                reasons.add(stops[a, b]); continue
            if root in semantic_graph and nx.has_path(semantic_graph, root, b):
                path = nx.shortest_path(semantic_graph, root, b)
                reason = next((stops[e] for e in zip(path, path[1:]) if e in stops), 'NOT_ON_SELECTED_SUBTREE')
                reasons.add(reason)
            else:
                reasons.add('NOT_ON_SELECTED_SUBTREE')
        rows.append(dict(original_branch_id=branch, retained_edges=len(kept), discarded_edges=len(missing),
            status='RETAINED' if not missing else ('PARTIALLY_RETAINED' if kept else 'DISCARDED'),
            discard_reasons=';'.join(sorted(reasons)),
            original_node_ids=';'.join(map(str, sorted({n for e in edges for n in e}))),
            retained_length_mm=sum(float(np.linalg.norm(semantic_graph.nodes[a]['coords'][:3]-semantic_graph.nodes[b]['coords'][:3])) for a,b in kept)))
    return rows


def export_candidate(directory, candidate, source, raw, semantic_graph, edge_branch):
    directory = Path(directory); directory.mkdir(parents=True)
    exports = {}
    for key, graph in [('original_radius', candidate['graph']), ('compensated', candidate['compensated'])]:
        exports[key] = export_manufacturing(directory / (key + '.swc'), graph, raw, candidate['core_nodes'], edge_branch, key == 'compensated')
    audit = branch_audit(semantic_graph, candidate['graph'], source, edge_branch, candidate['audit']['stops'])
    table(directory / 'branch_selection.csv', audit)
    return dict(mode=candidate['mode'], status='COMPACT_PRINT_CANDIDATE', source_component=source['component'],
        source_derived_component=source['derived_component'], source_level=source['level'], source_path=source['path'],
        refined_root_original_id=source['root'], selected_root=candidate['stats']['root_original_id'],
        selected_original_ids=sorted(candidate['graph']), core_original_ids=sorted(candidate['core_nodes']),
        selected_branches=[r['original_branch_id'] for r in audit if r['retained_edges']],
        discarded_eligible_branches=[r for r in audit if r['discarded_edges']], branch_audit=audit,
        stats=candidate['stats'], compensated_stats=metrics(candidate['compensated']),
        proximal_context=candidate['context'], radius_compensation=candidate['radii'],
        effective_limits=candidate['audit']['limits'], removed_twigs=candidate['audit']['removed_twigs'],
        selected_topology_paths=candidate['audit']['decisions'], gate=candidate['gate'],
        manufacturing_score=candidate['score'], blend_windows=candidate['blend'], exports=exports,
        adaptive_path=candidate.get('adaptive_search', {'enabled': False}),
        geometry_label='MANUFACTURING_COMPENSATED_GEOMETRY', manufacturing_validated=False)
