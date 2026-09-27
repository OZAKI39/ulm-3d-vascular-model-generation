"""Three explicit layers: immutable semantic input, native topology, print radii."""
from __future__ import annotations

import copy
import csv
import json
from pathlib import Path

import networkx as nx
import numpy as np
import yaml

from .boundary_review_export import load_cache, read_csv
from .mevo_refinement import branch_edges, graph_length, table
from .swc_export import read_source, validate_tree, write_swc
from .topbrain_qc import sha256


def merge_dict(base, update):
    result = copy.deepcopy(base)
    for key, value in update.items():
        result[key] = merge_dict(result.get(key, {}), value) if isinstance(value, dict) else value
    return result


def profile(path):
    path = Path(path).resolve()
    data = yaml.safe_load(path.read_text())
    if 'inherits' in data:
        data = merge_dict(profile(path.parent / data.pop('inherits')), data)
    nozzle = float(data['printer']['nozzle_diameter_mm'])
    m = data['manufacturing']
    if nozzle <= 0 or m['global_scale'] != 1.0:
        raise ValueError('Positive nozzle and global_scale=1.0 required')
    hard = max(m['hard_min_core_diameter_mm'], m['hard_nozzle_multiplier'] * nozzle)
    preferred = max(m['preferred_min_core_diameter_mm'], m['preferred_nozzle_multiplier'] * nozzle)
    floor = max(m['manufacturing_radius_floor_mm'], preferred / 2)
    usable = np.array(data['printer']['build_volume_mm'], float) - np.array([
        2 * data['printer']['bed_margin_xy_mm'], 2 * data['printer']['bed_margin_xy_mm'],
        data['printer']['top_margin_mm']])
    if not np.all(usable > 0) or hard > preferred or floor <= 0:
        raise ValueError('Invalid manufacturing thresholds or margins')
    data['effective'] = dict(hard_min_diameter_mm=hard, preferred_min_diameter_mm=preferred,
                             radius_floor_mm=floor, usable_volume_mm=usable.tolist())
    return data


def coordinates(graph):
    return {n: graph.nodes[n]['coords'] for n in graph}


def stable_pca(points):
    points = np.asarray(points, float)
    values, vectors = np.linalg.eigh(np.cov(points, rowvar=False))
    axes = vectors[:, np.argsort(values)[::-1]].T
    for i in range(3):
        index = np.argmax(np.abs(axes[i]))
        if axes[i, index] < 0:
            axes[i] *= -1
    if np.linalg.det(axes) < 0:
        axes[-1] *= -1
    return axes


def bbox(points):
    points = np.asarray(points, float)
    low, high = points.min(axis=0), points.max(axis=0)
    extents = high - low
    rotation = stable_pca(points)
    obb = np.ptp((points - points.mean(axis=0)) @ rotation.T, axis=0)
    return dict(min_mm=low.tolist(), max_mm=high.tolist(), extents_mm=extents.tolist(),
                longest_dimension_mm=float(extents.max()), diagonal_mm=float(np.linalg.norm(extents)),
                pca_obb_extents_mm=obb.tolist(), pca_obb_rotation=rotation.tolist(),
                obb_method='PCA aligned bounding box; not an exact minimum-volume OBB')


def geometry_stats(graph):
    xyzr = np.array([d['coords'] for _, d in graph.nodes(data=True)])
    diameter = 2 * xyzr[:, 3]
    return dict(node_count=len(graph), connected_component_count=nx.number_weakly_connected_components(graph),
        inlet_count=sum(graph.in_degree(n) == 0 for n in graph),
        outlet_count=sum(graph.out_degree(n) == 0 for n in graph),
        centerline_length_mm=graph_length(graph, coordinates(graph)), bbox=bbox(xyzr[:, :3]),
        diameter_min_mm=float(diameter.min()), diameter_median_mm=float(np.median(diameter)),
        diameter_max_mm=float(diameter.max()))


def load_semantic(directory):
    directory = Path(directory).resolve()
    cache = load_cache(directory.parent)
    summary = json.loads((directory / 'refinement_summary.json').read_text())
    snapshot = dict(cache.snapshot)
    snapshot.update(json.loads((directory / 'protection_before.json').read_text()))
    manifest = json.loads((directory / 'artifact_manifest.json').read_text())
    for rel, digest in manifest['files'].items():
        path = directory / rel
        if not path.resolve().is_relative_to(directory) or sha256(path) != digest:
            raise ValueError('Frozen refinement artifact mismatch: ' + rel)
        snapshot[str(path)] = digest
    snapshot[str(directory / 'artifact_manifest.json')] = sha256(directory / 'artifact_manifest.json')
    for path, digest in snapshot.items():
        if not Path(path).is_file() or sha256(path) != digest:
            raise ValueError('Frozen protected input mismatch: ' + path)
    raw = read_source(Path(cache.manifest['source']['raw_source']))
    semantic = {}
    for item in summary['exports']:
        if item['level'] != 'semantic_refined':
            continue
        path = Path(item['path'])
        assert sha256(path) == item['sha256']
        rows = read_csv(directory / 'mappings' / (path.stem + '_nodes.csv'))
        mapping = {int(r['roi_node_id']): int(r['original_swc_id']) for r in rows}
        actual = read_source(path).graph
        for n, original in mapping.items():
            np.testing.assert_array_equal(actual.nodes[n]['coords'], raw.graph.nodes[original]['coords'])
            assert actual.nodes[n]['swc_type'] == raw.graph.nodes[original]['swc_type']
        assert all(raw.graph.has_edge(mapping[a], mapping[b]) for a, b in actual.edges)
        semantic.setdefault(item['component'], set()).update(mapping.values())
    score_path = cache.directory / 'manual_review/tables/BG001_RMCA_component_summary.csv'
    if str(score_path) not in snapshot:
        raise ValueError('Frozen component scores not bound by source hashes')
    scores = read_csv(score_path)
    return cache, raw, semantic, scores, snapshot


def select_core(semantic, frozen_scores, raw, config):
    """Rank pre-existing component scores; no branch aggregation or classification."""
    roi = config['roi']
    rows = []
    for row in frozen_scores:
        component = int(row['component_id'].removeprefix('roi_part'))
        stats = geometry_stats(raw.graph.subgraph(semantic[component]))
        rows.append(dict(component=component, mean_pmevo=float(row['MeVO_mean_probability']),
            known_fraction=float(row['known_fraction']), frozen_branch_count=int(row['branch_count']),
            semantic_geometry=stats, status='OPTIONAL_SEMANTIC_COMPONENT'))
    weights = roi['score_weights']
    for r in rows:
        factors = dict(semantic=r['mean_pmevo'], known=r['known_fraction'],
            branches=r['frozen_branch_count'] / max(x['frozen_branch_count'] for x in rows),
            length=r['semantic_geometry']['centerline_length_mm'] / max(x['semantic_geometry']['centerline_length_mm'] for x in rows),
            spatial=r['semantic_geometry']['bbox']['diagonal_mm'] / max(x['semantic_geometry']['bbox']['diagonal_mm'] for x in rows))
        r['manufacturing_score'] = sum(weights[k] * v for k, v in factors.items())
        r['eligible_core'] = r['mean_pmevo'] >= roi['core_mean_pmevo_min'] and r['known_fraction'] >= roi['core_known_fraction_min']
    ranked = sorted((r for r in rows if r['eligible_core']), key=lambda r: (-r['manufacturing_score'], r['component']))
    selected = [r['component'] for r in ranked[:roi['max_core_components']]]
    if not selected:
        raise ValueError('NO_ELIGIBLE_SEMANTIC_CORE')
    for row in rows:
        if row['component'] in selected:
            row['status'] = 'SELECTED_SEMANTIC_CORE'
        row['manufacturing_role'] = 'SEMANTIC_CORE' if row['component'] in selected else 'OPTIONAL_CONTEXT'
        row['reason'] = 'TOP_MANUFACTURING_SCORE' if row['component'] in selected else (
            'CORE_COUNT_LIMIT_OPTIONAL' if row['eligible_core'] else 'BELOW_CORE_SEMANTIC_GATE_OPTIONAL')
    return selected, rows


def minimal_subtree(tree, selected_nodes):
    root = validate_tree(tree)
    common = set(tree)
    entry_nodes = [n for n in selected_nodes if not set(tree.predecessors(n)) & selected_nodes]
    for n in entry_nodes:
        common &= nx.ancestors(tree, n) | {n}
    depths = nx.single_source_shortest_path_length(tree, root)
    lca = max(common, key=lambda n: depths[n])
    retained = set(selected_nodes)
    for n in entry_nodes:
        retained.update(nx.shortest_path(tree, lca, n))
    result = tree.subgraph(retained).copy()
    validate_tree(result)
    return result


def extend_proximal(graph, full, roi_config):
    result = graph.copy()
    initial = validate_tree(graph)
    cursor, distance = initial, 0.
    added = []
    while True:
        extent = geometry_stats(result)['bbox']['longest_dimension_mm']
        if distance >= roi_config['minimum_inlet_stem_mm'] and extent >= roi_config['preferred_longest_dimension_mm']:
            break
        parents = list(full.predecessors(cursor))
        if not parents:
            break
        parent = parents[0]
        length = float(np.linalg.norm(full.nodes[parent]['coords'][:3] - full.nodes[cursor]['coords'][:3]))
        if distance + length > roi_config['proximal_context_max_mm'] + 1e-9:
            break
        result.add_node(parent, **copy.deepcopy(full.nodes[parent]))
        result.add_edge(parent, cursor)
        added.append(parent)
        distance += length
        cursor = parent
    return result, dict(original_lca=initial, final_root=cursor, added_original_ids=added,
                        extended_length_mm=distance, maximum_mm=roi_config['proximal_context_max_mm'])


def topology_paths(graph):
    starts = [n for n in graph if graph.in_degree(n) != 1 or graph.out_degree(n) != 1]
    paths = []
    for start in sorted(starts):
        for child in sorted(graph.successors(start)):
            path = [start, child]
            while graph.out_degree(path[-1]) == 1:
                path.append(next(iter(graph.successors(path[-1]))))
            paths.append(path)
    return paths


def prune_twigs(graph, config):
    """One pass on existing terminal branches; never recursively prune connectors."""
    kept = graph.copy()
    removed = []
    for path in topology_paths(graph):
        if graph.out_degree(path[-1]) or graph.out_degree(path[0]) < 2:
            continue
        xyzr = np.array([graph.nodes[n]['coords'] for n in path])
        length = float(np.linalg.norm(np.diff(xyzr[:, :3], axis=0), axis=1).sum())
        median = float(np.median(2 * xyzr[:, 3]))
        if length < config['manufacturing']['thin_twig_length_mm'] and median < config['effective']['hard_min_diameter_mm']:
            kept.remove_nodes_from(path[1:])
            removed.append(dict(original_node_ids=path[1:], retained_branchpoint=path[0],
                                length_mm=length, original_median_diameter_mm=median,
                                reason='SHORT_THIN_TERMINAL_TWIG'))
    validate_tree(kept)
    return kept, removed


def smoothstep(t):
    t = np.clip(t, 0., 1.)
    return t * t * (3 - 2 * t)


def compensate_radius(graph, floor, blend_length):
    """Local, upward-only smoothstep shoulders; no xyz or global radius filtering.

    Each activation crossing gets a length-limited smoothstep shoulder. Taking
    the maximum with raw/floor prevents shrinking anatomy or violating the floor.
    Junction radii are shared (maximum of incident proposals), not duplicated.
    These are discrete samples; VascularMD subsequently fits continuous splines.
    """
    result = copy.deepcopy(graph)
    base = {n: max(float(graph.nodes[n]['coords'][3]), floor) for n in graph}
    proposed = dict(base)
    windows = []
    for path in topology_paths(graph):
        xyzr = np.array([graph.nodes[n]['coords'] for n in path])
        s = np.r_[0., np.cumsum(np.linalg.norm(np.diff(xyzr[:, :3], axis=0), axis=1))]
        active = xyzr[:, 3] < floor
        values = np.array([base[n] for n in path])
        for i in np.flatnonzero(active[:-1] != active[1:]):
            middle = (s[i] + s[i+1]) / 2
            a = max(0, int(np.searchsorted(s, middle - blend_length/2, side='right')) - 1)
            b = min(len(s)-1, int(np.searchsorted(s, middle + blend_length/2)))
            if b <= a or s[b] <= s[a]:
                continue
            target = values[a] + (values[b]-values[a]) * smoothstep((s[a:b+1]-s[a])/(s[b]-s[a]))
            for j, value in zip(range(a,b+1), target):
                proposed[path[j]] = max(proposed[path[j]], float(value))
            windows.append(dict(first_original_id=path[a], last_original_id=path[b],
                                actual_span_mm=float(s[b]-s[a]), requested_blend_mm=blend_length))
    for n in result:
        result.nodes[n]['coords'][3] = proposed[n]
        np.testing.assert_array_equal(result.nodes[n]['coords'][:3], graph.nodes[n]['coords'][:3])
    return result, windows


def roi_size_gate(stats, config):
    b, roi = stats['bbox'], config['roi']
    reasons = []
    if b['longest_dimension_mm'] < roi['min_longest_dimension_mm'] or b['diagonal_mm'] < roi['min_bbox_diagonal_mm']:
        reasons.append('ROI_BELOW_PROJECT_SIZE_TARGET')
    if b['longest_dimension_mm'] > roi['max_longest_dimension_mm']:
        reasons.append('ROI_ABOVE_PROJECT_SIZE_TARGET')
    if stats['inlet_count'] != 1 or stats['connected_component_count'] != 1 or stats['outlet_count'] < roi['minimum_outlets']:
        reasons.append('ROI_CONNECTIVITY_OR_PORT_GATE_FAILED')
    return dict(passed=not reasons, reasons=reasons,
        preferred_size=roi['preferred_longest_dimension_mm'] <= b['longest_dimension_mm'] <= roi['preferred_longest_dimension_upper_mm'])


def export_manufacturing(path, graph, source, core_nodes, edge_branch, compensated):
    header = ['MANUFACTURING_COMPENSATED_GEOMETRY' if compensated else 'MANUFACTURING_ROI_TOPOLOGY_WITH_ORIGINAL_RADIUS',
        'Original BraVa xyz and TYPE preserved; directed native edges only; units mm',
        'MANUFACTURING RADIUS COMPENSATION APPLIED' if compensated else 'Original BraVa radius unchanged',
        'Original BraVa radius preserved in mapping CSV', 'Not anatomical lumen ground truth',
        'Scale=1.0; no orientation transform applied to this SWC']
    ids = write_swc(path, graph, source.path, 'manufacturing-derivative', 'original-samples', header_lines=header)
    actual = read_source(path).graph
    assert set(actual.edges) == {(ids[a],ids[b]) for a,b in graph.edges}
    rows = []
    for n, new in ids.items():
        original = source.graph.nodes[n]['coords']
        value = graph.nodes[n]['coords']
        np.testing.assert_array_equal(actual.nodes[new]['coords'], value)
        np.testing.assert_array_equal(original[:3], value[:3])
        assert actual.nodes[new]['swc_type'] == source.graph.nodes[n]['swc_type']
        if not compensated:
            np.testing.assert_array_equal(original, value)
        rows.append(dict(swc_id=new, original_swc_id=n, original_radius_mm=float(original[3]),
            manufacturing_radius_mm=float(value[3]), radius_added_mm=float(value[3]-original[3]),
            diameter_original_mm=float(2*original[3]), diameter_manufacturing_mm=float(2*value[3]),
            is_anatomical_mevo=n in core_nodes, is_manufacturing_context=n not in core_nodes))
    table(path.with_name(path.stem+'_mapping.csv'), rows)
    table(path.with_name(path.stem+'_edges.csv'), [dict(parent_id=ids[a], child_id=ids[b],
        original_parent_id=a, original_child_id=b, original_branch=edge_branch[a,b],
        is_anatomical_mevo=a in core_nodes and b in core_nodes,
        is_manufacturing_context=not (a in core_nodes and b in core_nodes)) for a,b in graph.edges])
    return dict(path=str(path), sha256=sha256(path), round_trip=True, xyz_type_edges_preserved=True,
                original_radius_preserved=not compensated)
