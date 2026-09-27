#!/usr/bin/env python3
"""Validate the actual SV mesh and export compact, consistently oriented solver faces."""
import argparse
import json
import sys
from pathlib import Path
import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'mesh_generate/src'))
from sv_validation.geometry import (ROLE_IDS, polydata, triangles_of, surface_stats,
                                    enclosed_volume, closest_surface_distances, surface_samples)
from sv_validation.mesh_diagnostics import volume_validity, tetra_sicn, quality_summary
from sv_validation.validation import require, face_mapping
from sv_validation.provenance import sha256, write_json


def compact_surface(points, triangles):
    used, inverse = np.unique(triangles, return_inverse=True)
    surface = polydata(points[used], inverse.reshape(-1, 3))
    surface.point_data['GlobalNodeID'] = (used + 1).astype(np.int32)
    return surface


def one_port(triangles):
    nodes, inverse = np.unique(triangles, return_inverse=True)
    local = inverse.reshape(-1, 3)
    edges = np.sort(np.concatenate([local[:, [0, 1]], local[:, [0, 2]], local[:, [1, 2]]]), axis=1)
    unique, counts = np.unique(edges, axis=0, return_counts=True)
    graph = coo_matrix((np.ones(2*len(unique)), (np.r_[unique[:, 0], unique[:, 1]], np.r_[unique[:, 1], unique[:, 0]])),
                       shape=(len(nodes), len(nodes)))
    require(connected_components(graph, directed=False, return_labels=False) == 1, 'Split port')
    require(np.all(counts <= 2), 'Nonmanifold port')
    rim = unique[counts == 1]
    _, degrees = np.unique(rim, return_counts=True)
    require(len(rim) > 0 and np.all(degrees == 2), 'Port rim is not a closed loop')
    require(len(nodes) - len(unique) + len(triangles) == 1, 'Port is not a topological disk')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--attempt', choices=['primary', 'fallback'], default='primary')
    args = parser.parse_args()
    raw = ROOT / 'outputs/sv1/mesh_generation' / args.attempt
    generation = json.loads((raw / 'generation.json').read_text())
    require(generation['status'] == 'PASS', 'Mesh generation did not pass')
    policy_path = ROOT / 'configs/mesh_policy.json'
    policy = json.loads(policy_path.read_text())
    frozen = json.loads((ROOT / 'reports/sv1/mesh_policy_provenance.json').read_text())
    require(sha256(policy_path) == frozen['sha256'], 'Meshing policy changed')
    grid = pv.read(raw / 'volume.vtu')
    require(set(grid.cells_dict) == {pv.CellType.TETRA}, 'Only linear tetrahedra are accepted')
    points = np.asarray(grid.points, dtype=np.float64).copy()
    tetra = grid.cells_dict[pv.CellType.TETRA].astype(np.int64)
    validity = volume_validity(points, tetra)
    quality_values = tetra_sicn(points, tetra)
    quality = quality_summary(quality_values)
    require(quality['q_min'] > 0, 'Nonpositive signed element quality')

    all_faces = np.concatenate([tetra[:, ids] for ids in ([0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3])])
    owners = np.tile(np.arange(len(tetra)), 4)
    keys, first, counts = np.unique(np.sort(all_faces, axis=1), axis=0, return_index=True, return_counts=True)
    boundary = all_faces[first[counts == 1]].copy()
    adjacent = owners[first[counts == 1]]
    xyz = points[boundary]
    outward = xyz.mean(axis=1) - points[tetra[adjacent]].mean(axis=1)
    wrong = np.einsum('ij,ij->i', np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0]), outward) < 0
    boundary[wrong] = boundary[wrong][:, [0, 2, 1]]
    lookup = {tuple(sorted(row)): i for i, row in enumerate(boundary)}
    tags = np.zeros(len(boundary), np.int32)
    tree = cKDTree(points)
    maximum_match_distance = 0.
    for face_id in ROLE_IDS:
        face = pv.read(raw / ('face_' + str(face_id) + '.vtp'))
        require(face.n_cells > 0, 'Missing named face')
        distances, matched = tree.query(face.points)
        maximum_match_distance = max(maximum_match_distance, float(distances.max()))
        require(distances.max() <= 1e-6 * policy['h_ref_m'], 'Face vertices do not match actual volume vertices')
        face_triangles = matched[triangles_of(face)]
        if face_id != 1:
            one_port(face_triangles)
        for row in face_triangles:
            key = tuple(sorted(row))
            require(key in lookup, 'Labeled face is not an exterior tetrahedron facet')
            index = lookup[key]
            require(tags[index] == 0, 'Duplicate or merged exterior face label')
            tags[index] = face_id
    require(np.all(tags > 0), 'Unlabeled exterior surface')
    mapping = json.loads((ROOT / 'configs/face_map.json').read_text())
    contract = json.loads((ROOT / 'inputs/fem_reference/port_contract.json').read_text())
    face_mapping(mapping['faces'], contract['facet_names'])

    with np.load(ROOT / 'inputs/fem_reference/exterior_surface.npz') as original:
        source_points = original['points_m']
        source_triangles = original['triangles']
        source_tags = original['facet_tags']
    ports = {}
    for face_id, role in ROLE_IDS.items():
        if role == 'WALL':
            continue
        before = surface_stats(source_points, source_triangles[source_tags == face_id])
        after = surface_stats(points, boundary[tags == face_id])
        normal_dot = float(np.dot(before['normal'], after['normal']))
        require(normal_dot > 0, 'Port orientation flip')
        ports[role] = {'original': before, 'actual': after, 'normal_dot': normal_dot,
                       'area_relative_error': abs(after['area_m2']-before['area_m2'])/before['area_m2'],
                       'centroid_shift_m': float(np.linalg.norm(np.array(after['centroid_m'])-before['centroid_m'])),
                       'connected_regions': 1, 'topological_disk': True, 'orientation_flip': False}
    original_wall = compact_surface(source_points, source_triangles[source_tags == 1])
    new_wall = compact_surface(points, boundary[tags == 1])
    forward = closest_surface_distances(surface_samples(original_wall), new_wall)
    backward = closest_surface_distances(surface_samples(new_wall), original_wall)
    distances = np.r_[forward, backward]
    source_volume = abs(enclosed_volume(source_points, source_triangles))
    actual_surface_volume = enclosed_volume(points, boundary)
    require(np.isclose(actual_surface_volume, validity['volume_m3'], rtol=1e-11, atol=0), 'Surface and tetrahedral volume disagree')
    geometry = {'status': 'PASS_TOPOLOGY_GEOMETRY_METRICS_FOR_REVIEW',
                'wall': {'max_m': float(distances.max()), 'P95_m': float(np.quantile(distances, .95)),
                         'RMS_m': float(np.sqrt(np.mean(distances**2))),
                         'max_over_h_wall': float(distances.max()/policy['h_wall_m']),
                         'forward_sample_count': len(forward), 'backward_sample_count': len(backward),
                         'method': 'Every vertex and triangle centroid to the other triangle surface; sampled bidirectional distance'},
                'source_enclosed_volume_m3': source_volume, 'actual_enclosed_volume_m3': actual_surface_volume,
                'volume_relative_error': abs(actual_surface_volume-source_volume)/source_volume,
                'ports': ports, 'boundary_roles_complete': True,
                'face_vertex_max_match_distance_m': maximum_match_distance}

    output = ROOT / 'outputs/sv1/SV_MESH'
    (output / 'mesh-surfaces').mkdir(parents=True, exist_ok=True)
    # Coordinates and cells are the actual SV result; only standard identity arrays are assigned.
    grid.point_data['GlobalNodeID'] = np.arange(1, len(points)+1, dtype=np.int32)
    grid.cell_data['GlobalElementID'] = np.arange(1, len(tetra)+1, dtype=np.int32)
    grid.cell_data['minSICN'] = quality_values
    grid.save(output / 'mesh-complete.mesh.vtu')
    surface = compact_surface(points, boundary)
    surface.cell_data['ModelFaceID'] = tags
    surface.cell_data['GlobalElementID'] = (adjacent+1).astype(np.int32)
    surface.save(output / 'mesh-complete.exterior.vtp')
    for face_id, role in ROLE_IDS.items():
        selected = tags == face_id
        face = compact_surface(points, boundary[selected])
        face.cell_data['GlobalElementID'] = (adjacent[selected]+1).astype(np.int32)
        face.save(output / 'mesh-surfaces' / (role + '.vtp'))
    np.savez_compressed(output / 'mesh_arrays.npz', points_m=points, tetra=tetra,
                        boundary_triangles=boundary, facet_tags=tags, adjacent_tetra=adjacent,
                        min_sicn=quality_values)
    np.savez_compressed(output / 'wall_distances.npz', source_to_new_m=forward, new_to_source_m=backward)
    write_json(ROOT / 'reports/sv1/mesh_validity.json', dict(validity, status='PASS', generation=generation,
               boundary_roles_complete=True, raw_volume_sha256=sha256(raw/'volume.vtu'),
               solver_mesh_sha256=sha256(output/'mesh-complete.mesh.vtu')))
    write_json(ROOT / 'reports/sv1/mesh_quality.json', dict(quality, tetra=len(tetra), status='MEASURED',
               metric='minSICN: constant linear tetrahedral ideal-map signed inverse Frobenius condition number'))
    write_json(ROOT / 'reports/sv1/geometry_qc.json', geometry)
    print(json.dumps({'validity': validity, 'quality': quality, 'wall': geometry['wall'],
                      'volume_relative_error': geometry['volume_relative_error'],
                      'ports': {k: {'area_relative_error': v['area_relative_error'], 'normal_dot': v['normal_dot']} for k, v in ports.items()}}, indent=2))


if __name__ == '__main__':
    main()
