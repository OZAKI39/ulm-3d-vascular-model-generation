"""Compact surface orientation with inlet and distal-cantilever preferences."""
from pathlib import Path

import numpy as np
from scipy.spatial import ConvexHull
from threadpoolctl import threadpool_limits

from .manufacturing_roi import topology_paths
from .mevo_refinement import table
from .print_orientation import (branch_descriptors, orientation_rotations, score_orientation,
                                triangle_geometry, transform_points)
from .swc_export import validate_tree


def compact_preferences(graph, matrix, config, height):
    root = validate_tree(graph)
    path = next(p for p in topology_paths(graph) if p[0] == root)
    points = np.array([graph.nodes[n]['coords'][:3] for n in path])
    vector = points[-1] - points[0]
    angle = float(np.degrees(np.arcsin(np.clip(abs(vector @ matrix[2, :3]) / np.linalg.norm(vector), 0, 1))))
    limits = config['orientation']['preferred_stem_angle_deg']
    penalty = max(limits['min'] - angle, angle - limits['max'], 0.) / 90.
    inlet_height = float(transform_points(np.array([points[0]]), matrix)[0, 2])
    total = horizontal = 0.
    for p in topology_paths(graph):
        if p[0] == root:
            continue
        xyzr = np.array([graph.nodes[n]['coords'] for n in p]); vectors = np.diff(xyzr[:, :3], axis=0)
        length = np.linalg.norm(vectors, axis=1)
        if length.sum() < config['orientation']['long_branch_length_mm']:
            continue
        angles = np.degrees(np.arcsin(np.clip(np.divide(abs(vectors @ matrix[2, :3]), length, out=np.zeros_like(length), where=length>0), 0, 1)))
        total += float(length.sum())
        horizontal += float(length[angles < config['orientation']['horizontal_branch_angle_deg']].sum())
    return dict(inlet_height_mm=inlet_height, proximal_stem_angle_from_plate_deg=angle,
        inlet_height_penalty=max(0., inlet_height) / max(height, 1e-12), stem_angle_penalty=penalty,
        long_horizontal_distal_fraction=horizontal / total if total else 0.,
        stem_angle_method='Chord of the proximal topological stem, not anatomical orientation')


def optimize_compact(mesh, graph, config, output):
    output = Path(output); output.mkdir(parents=True, exist_ok=True)
    points = np.asarray(mesh.points, float); areas, normals, centers = triangle_geometry(mesh)
    hull = points[ConvexHull(points).vertices]; branches = branch_descriptors(graph)
    rotations = orientation_rotations(points, config); settings = config['orientation']
    if not settings['minimum_candidates'] <= len(rotations) <= settings['maximum_candidates']:
        raise ValueError('Compact orientation search must contain 50--150 candidates')
    results = []
    with threadpool_limits(limits=1):
        for index, (method, rotation, offsets) in enumerate(rotations):
            result = score_orientation(hull, areas, normals, centers, rotation, config, branches)
            preferences = compact_preferences(graph, np.asarray(result['transform_4x4']), config, result['height_mm'])
            w = settings['compact_score_weights']
            result['surface_score'] = result['score']
            result['score'] = (w['surface_score'] * result['surface_score'] + w['inlet_height'] * preferences['inlet_height_penalty']
                + w['stem_angle'] * preferences['stem_angle_penalty'] + w['distal_horizontal'] * preferences['long_horizontal_distal_fraction'])
            result.update(preferences, candidate_id=index, method=method, pca_offset_euler_xyz_deg=offsets)
            results.append(result)
    feasible = sorted((r for r in results if r['fits_build_volume']), key=lambda r: (r['score'], r['candidate_id']))
    top = feasible[:3]
    for rank, result in enumerate(top, 1):
        transformed = mesh.copy(); transformed.points = transform_points(points, result['transform_4x4'])
        path = output / f'orientation_{rank:02d}.stl'; transformed.save(path)
        result.update(rank=rank, stl=str(path))
        np.testing.assert_allclose(transform_points(transformed.points, np.linalg.inv(result['transform_4x4'])), points, atol=1e-10)
    baseline = mesh.copy(); baseline.points = transform_points(points, results[0]['transform_4x4'])
    baseline.save(output / 'orientation_baseline.stl')
    table(output / 'orientation_scores.csv', [dict(candidate_id=r['candidate_id'], score=r['score'], feasible=r['fits_build_volume'],
        height_mm=r['height_mm'], inlet_height_mm=r['inlet_height_mm'], stem_angle_deg=r['proximal_stem_angle_from_plate_deg'],
        downward_area_mm2=r['downward_support_area_mm2'], distal_horizontal_fraction=r['long_horizontal_distal_fraction']) for r in results])
    return dict(tested=len(results), feasible_count=len(feasible), top_candidates=top, baseline=results[0],
                all_candidates=results, orientation_scope='BALANCED only', scale=1.)
