"""Physical tetra/wall diagnostics. Distances are SI; angles are internal degrees."""
from itertools import combinations

import numpy as np
from scipy.spatial import cKDTree

from .cap_remesh import triangle_quality
from .mesh_input import FACET_NAMES
from .mesh_qc import quantiles, tetra_volumes, triangle_geometry

EDGES = np.array(list(combinations(range(4), 2)))
FACES = np.array([[1, 2, 3], [0, 2, 3], [0, 1, 3], [0, 1, 2]])


def tetra_shape(vertices):
    """Translation/scale-stable geometry, outward face normals, six inner angles."""
    x = np.asarray(vertices, dtype=float)
    lengths = np.linalg.norm(x[EDGES[:, 0]] - x[EDGES[:, 1]], axis=1)
    scale = float(lengths.max())
    if not np.isfinite(x).all() or scale <= 0:
        raise ValueError('Invalid tetra coordinates')
    y = (x - x[0]) / scale
    signed = float(np.linalg.det(y[1:]) / 6)
    if signed == 0:
        raise ValueError('Degenerate tetra')
    volume = abs(signed) * scale**3
    cross = np.cross(y[FACES[:, 1]] - y[FACES[:, 0]], y[FACES[:, 2]] - y[FACES[:, 0]])
    areas_normalized = np.linalg.norm(cross, axis=1) / 2
    if np.any(areas_normalized == 0):
        raise ValueError('Degenerate face')
    normals = cross / (2 * areas_normalized[:, None])
    inward = np.einsum('ij,ij->i', normals, y - y[FACES].mean(axis=1)) > 0
    normals[inward] *= -1
    angles = []
    for i, j in combinations(range(4), 2):
        angles.append({'edge_local_vertices': sorted(map(int, set(FACES[i]) & set(FACES[j]))),
                       'face_opposite_vertices': [i, j],
                       'degrees': float(np.degrees(np.arccos(np.clip(-normals[i] @ normals[j], -1, 1))))})
    fq, _ = triangle_quality(y, FACES)
    center = np.linalg.solve(2 * y[1:], np.einsum('ij,ij->i', y[1:], y[1:]))
    circumradius = float(np.linalg.norm(center) * scale)
    inradius = float(3 * volume / (areas_normalized.sum() * scale**2))
    return {'vertex_coordinates_m': x.tolist(), 'centroid_m': x.mean(axis=0).tolist(),
            'volume_m3': volume, 'signed_volume_m3': signed * scale**3,
            'six_edge_lengths_m': lengths.tolist(), 'edge_ratio': float(lengths.max() / lengths.min()),
            'four_face_areas_m2': (areas_normalized * scale**2).tolist(), 'four_face_quality': fq.tolist(),
            'minimum_face_quality': float(fq.min()), 'maximum_face_quality': float(fq.max()),
            'dihedral_angles': angles, 'minimum_dihedral_degrees': min(a['degrees'] for a in angles),
            'maximum_dihedral_degrees': max(a['degrees'] for a in angles),
            'circumradius_m': circumradius, 'inradius_m': inradius,
            'radius_ratio_3r_over_R': 3 * inradius / circumradius,
            'radius_ratio_definition': '3 * inradius / circumradius; regular tetra = 1'}


def face_neighbors(tetra):
    """Face-adjacent tetra; target excluded, never vertex-only adjacency."""
    faces = np.sort(np.asarray(tetra)[:, FACES].reshape(-1, 3), axis=1)
    _, inverse, count = np.unique(faces, axis=0, return_inverse=True, return_counts=True)
    if np.any(count > 2):
        raise ValueError('Nonmanifold tetra face')
    order = np.argsort(inverse, kind='stable')
    offsets = np.r_[0, np.cumsum(count)]
    repeated = np.flatnonzero(count == 2)
    left = order[offsets[repeated]] // 4
    right = order[offsets[repeated] + 1] // 4
    adjacency = [set() for _ in tetra]
    for a, b in zip(left, right):
        adjacency[a].add(int(b)); adjacency[b].add(int(a))
    return adjacency


def rings(adjacency, cell, depth=2):
    seen = {int(cell)}
    frontier = {int(cell)}
    result = []
    for _ in range(depth):
        frontier = set().union(*(adjacency[i] for i in frontier)) - seen if frontier else set()
        seen |= frontier
        result.append(sorted(seen - {int(cell)}))
    return result


def wall_adjacency(triangles, tags):
    indices = np.flatnonzero(tags == 1)
    edge_faces = {}
    vertex_faces = {}
    for i in indices:
        for v in triangles[i]:
            vertex_faces.setdefault(int(v), set()).add(int(i))
        for a, b in combinations(triangles[i], 2):
            edge_faces.setdefault(tuple(sorted((int(a), int(b)))), []).append(int(i))
    adjacency = {int(i): set() for i in indices}
    for faces in edge_faces.values():
        for i in faces:
            adjacency[i].update(set(faces) - {i})
    return adjacency, vertex_faces


def classify_reason(record, policy):
    """Conservative evidence labels, not a claim of unavoidable CFD error."""
    wall = record.get('wall_diagnostics')
    if not wall or not record.get('boundary_constraints'):
        return {'category': 'UNKNOWN', 'evidence': ['Missing direct boundary/face evidence'], 'causal_proof': False}
    cfg = policy['classification']
    bc = record['boundary_constraints']
    actual = wall['actual_wall_faces']
    locked = bc['boundary_vertex_count'] == 4 and len(actual) >= 2
    poor = any(f['q_tri'] < cfg['poor_wall_face_quality'] for f in actual)
    rapid = wall['normal_variation']['maximum_degrees'] > cfg['rapid_normal_angle_degrees']
    evidence = [f"{bc['boundary_vertex_count']}/4 vertices frozen on exterior", f'{len(actual)} actual WALL faces']
    if locked:
        evidence.append('At fixed connectivity all four vertex positions are immutable; smoothing has no degree of freedom')
        return {'category': 'WALL_TRIANGLE_CONSTRAINED', 'evidence': evidence, 'causal_proof': 'fixed-connectivity constraint only; retriangulation may still help'}
    if poor:
        evidence.append('An actual frozen WALL face is below the predeclared face-quality threshold')
    if rapid:
        evidence.append('Two edge-rings of wall faces exceed the predeclared normal-variation threshold')
    if poor and rapid:
        return {'category': 'MIXED', 'evidence': evidence, 'causal_proof': False}
    if poor:
        return {'category': 'WALL_TRIANGLE_CONSTRAINED', 'evidence': evidence, 'causal_proof': False}
    if rapid and actual:
        return {'category': 'RAPID_GEOMETRIC_VARIATION', 'evidence': evidence, 'causal_proof': False}
    if not actual and bc['boundary_vertex_count'] == 0 and record['minimum_dihedral_degrees'] < cfg['sliver_min_angle_degrees'] and record['minimum_face_quality'] >= cfg['poor_wall_face_quality']:
        return {'category': 'INTERIOR_SLIVER', 'evidence': evidence + ['No constrained exterior vertex or face; small dihedral with adequate faces'], 'causal_proof': False}
    return {'category': 'UNKNOWN', 'evidence': evidence + ['Available metrics do not isolate a cause'], 'causal_proof': False}


def diagnose(data, policy):
    points, tetra, boundary = data['points_m'], data['tetra'], data['boundary_triangles']
    tags, q = data['facet_tags'], data['min_sicn']
    low = np.flatnonzero(q < policy['quality']['low_threshold'])
    adjacency = face_neighbors(tetra)
    wall_adj, vertex_faces = wall_adjacency(boundary, tags)
    areas, centers, vectors = triangle_geometry(points, boundary)
    normals = vectors / areas[:, None]
    bq, _ = triangle_quality(points, boundary)
    wall_ids = np.flatnonzero(tags == 1)
    wall_sorted = np.sort(bq[wall_ids])
    boundary_nodes = set(map(int, np.unique(boundary)))
    bmap = {tuple(sorted(map(int, face))): i for i, face in enumerate(boundary)}
    tree, walltree = cKDTree(centers), cKDTree(centers[wall_ids])

    def face_record(i):
        coords = points[boundary[i]]
        edge = np.linalg.norm(coords[[0, 1, 2]] - coords[[1, 2, 0]], axis=1)
        return {'triangle_index': int(i), 'tag': int(tags[i]), 'vertex_ids': boundary[i].tolist(),
                'coordinates_m': coords.tolist(), 'edge_lengths_m': edge.tolist(), 'area_m2': float(areas[i]),
                'q_tri': float(bq[i]), 'aspect_ratio_longest_over_shortest_edge': float(edge.max() / edge.min()),
                'aspect_ratio_longest_over_min_altitude': float(edge.max()**2 / (2 * areas[i])),
                'wall_quality_percentile': float(100 * np.searchsorted(wall_sorted, bq[i], side='right') / len(wall_sorted))}

    records = []
    for index in low:
        r = tetra_shape(points[tetra[index]])
        distance, nearest = tree.query(r['centroid_m'])
        wall_distance, wi = walltree.query(r['centroid_m']); wi = int(wall_ids[wi])
        neighborhoods = rings(adjacency, int(index))
        actual = [bmap[tuple(sorted(map(int, face)))] for face in tetra[index][FACES]
                  if tuple(sorted(map(int, face))) in bmap]
        wall_actual = [i for i in actual if tags[i] == 1]
        wall_rings = rings(wall_adj, wi)
        patch = sorted({wi, *wall_rings[-1]})
        angle = np.degrees(np.arccos(np.clip(normals[patch] @ normals[wi], -1, 1)))
        proxies = []
        for j in patch:
            for k in wall_adj[j]:
                if k > j and k in patch:
                    sep = np.linalg.norm(centers[k] - centers[j])
                    if sep > 0:
                        proxies.append(float(np.arccos(np.clip(normals[k] @ normals[j], -1, 1)) / sep))
        local_ids = np.array([int(index)] + neighborhoods[-1])
        local_edges = np.unique(np.sort(tetra[local_ids][:, EDGES].reshape(-1, 2), axis=1), axis=0)
        local_median = float(np.median(np.linalg.norm(points[local_edges[:, 0]] - points[local_edges[:, 1]], axis=1)))
        r.update(cell_index=int(index), gmsh_element_id=int(data['gmsh_element_ids'][index]), vertex_ids=tetra[index].tolist(),
                 min_sicn=float(q[index]), nearest_boundary_triangle_id=int(nearest), nearest_boundary_tag=int(tags[nearest]),
                 nearest_boundary_name=FACET_NAMES[int(tags[nearest])], distance_to_nearest_boundary_triangle_center_m=float(distance),
                 one_ring_neighbor_tetra_ids=neighborhoods[0], two_ring_neighbor_tetra_ids=neighborhoods[1],
                 neighborhoods={f'{n + 1}_ring': {'cell_ids': ids, 'count': len(ids), 'quality': quantiles(q[ids]) if ids else None,
                                               'low_count': int(np.count_nonzero(q[ids] < policy['quality']['low_threshold']))}
                                for n, ids in enumerate(neighborhoods)},
                 local_median_tetra_edge_length_m=local_median,
                 boundary_constraints={'boundary_vertex_count': sum(int(v) in boundary_nodes for v in tetra[index]),
                                       'frozen_vertex_ids': [int(v) for v in tetra[index] if int(v) in boundary_nodes],
                                       'actual_boundary_face_ids': list(map(int, actual))},
                 wall_diagnostics={'nearest_wall_triangle': face_record(wi), 'nearest_wall_distance_m': float(wall_distance),
                                   'actual_wall_faces': [face_record(i) for i in wall_actual], 'patch_triangle_ids': patch,
                                   'neighbor_wall_quality': quantiles(bq[wall_rings[-1]]) if wall_rings[-1] else None,
                                   'patch_quality': quantiles(bq[patch]),
                                   'normal_variation': {'maximum_degrees': float(angle.max()), 'median_degrees': float(np.median(angle)),
                                                        'definition': 'Angles relative to nearest outward wall face in two edge-connected wall rings'},
                                   'local_vertex_wall_triangle_valence': {str(v): len(vertex_faces[int(v)]) for v in np.unique(boundary[patch])},
                                   'curvature_proxy': {'values_rad_per_m': proxies, 'median_rad_per_m': float(np.median(proxies)) if proxies else None,
                                                       'definition': 'Adjacent outward-normal angle / triangle-center distance',
                                                       'interpretation': 'Discrete normal-gradient diagnostic; not a fitted smooth curvature or proof of cause'}})
        r['reason'] = classify_reason(r, policy)
        records.append(r)
    return {'low_threshold': policy['quality']['low_threshold'], 'residual_count': len(records), 'cells': records,
            'wall_quality_distribution': quantiles(bq[wall_ids]), 'wall_triangle_count': len(wall_ids),
            'neighborhood_definition': 'Face adjacency; cumulative 1/2 rings exclude the target cell',
            'classification_note': 'Evidence-based diagnostic categories; no claim that all boundary-fixed topologies are impossible'}
