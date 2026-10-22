"""Read-only PLC seam, labels and orientation checks; no point welding."""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
import pyvista as pv
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


def surface_qa(path, radius, length):
    mesh = pv.read(path)
    points = np.asarray(mesh.points, dtype=float)
    faces = mesh.faces.reshape(-1, 4)
    assert np.all(faces[:, 0] == 3), 'Only triangular PLC supported'
    triangles = faces[:, 1:]
    tags = np.asarray(mesh.cell_data['ModelFaceID'])
    directed = np.concatenate([triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]]])
    owners = np.tile(np.arange(len(triangles)), 3)
    edges, inverse, counts = np.unique(np.sort(directed, axis=1), axis=0, return_inverse=True, return_counts=True)
    orientation_sum = np.bincount(inverse, weights=np.where(directed[:, 0] < directed[:, 1], 1, -1))
    order = np.argsort(inverse, kind='stable')
    starts = np.r_[0, np.cumsum(counts)[:-1]]
    adjacent = owners[order[np.c_[starts[counts == 2], starts[counts == 2] + 1]]]
    graph = coo_matrix((np.ones(2 * len(adjacent)),
                       (np.r_[adjacent[:, 0], adjacent[:, 1]], np.r_[adjacent[:, 1], adjacent[:, 0]])),
                      shape=(len(triangles), len(triangles))).tocsr()
    xyz = points[triangles]
    area_vec = np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]) / 2
    centers = xyz.mean(axis=1)
    area = np.linalg.norm(area_vec, axis=1)
    volume = float(np.einsum('ni,ni->n', area_vec, centers).sum() / 3)
    _, coordinate_counts = np.unique(points, axis=0, return_counts=True)
    face_components = {str(int(tag)): int(connected_components(graph[tags == tag][:, tags == tag], directed=False)[0])
                       for tag in np.unique(tags)}
    wall_vertices = np.unique(triangles[tags == 1])
    ports = {}
    for tag, name, outward, end in [(2, 'INLET', -1, 0.), (3, 'OUTLET', 1, length)]:
        port_faces = triangles[tags == tag]
        cap_edges, cap_counts = np.unique(np.sort(np.concatenate([port_faces[:, [0, 1]], port_faces[:, [1, 2]], port_faces[:, [2, 0]]]), axis=1), axis=0, return_counts=True)
        loop_edges = cap_edges[cap_counts == 1]
        loop_vertices = np.unique(loop_edges)
        _, duplicates = np.unique(points[loop_vertices], axis=0, return_counts=True)
        cap_wall = [(int(a), int(b)) for a, b in adjacent if {int(tags[a]), int(tags[b])} == {1, tag}]
        signed_area = float(area_vec[tags == tag, 2].sum())
        ports[name] = dict(face_id=tag, loop_vertex_count=len(loop_vertices), loop_edge_count=len(loop_edges),
                           duplicate_coordinate_count=int(np.sum(duplicates - 1)),
                           shared_topology_vertex_count=len(np.intersect1d(wall_vertices, loop_vertices)),
                           shared_topology_edge_count=len(cap_wall),
                           orientation='OUTWARD' if signed_area * outward > 0 and np.all(area_vec[tags == tag, 2] * outward > 0) else 'INCONSISTENT',
                           cap_area_m2=float(area[tags == tag].sum()),
                           cap_area_relative_error=abs(float(area[tags == tag].sum()) / (np.pi * radius ** 2) - 1),
                           axial_coordinate_max_error_m=float(np.max(abs(points[np.unique(port_faces), 2] - end))))
    wall_orientation = np.einsum('ni,ni->n', area_vec[tags == 1, :2], centers[tags == 1, :2])
    output = dict(source=str(path), sha256=hashlib.sha256(Path(path).read_bytes()).hexdigest(),
                  points=len(points), triangles=len(triangles), edges=len(edges),
                  free_edges=int(np.sum(counts == 1)), nonmanifold_edges=int(np.sum(counts > 2)),
                  inconsistent_oriented_edges=int(np.sum(orientation_sum != 0)),
                  global_duplicate_coordinate_count=int(np.sum(coordinate_counts - 1)),
                  zero_area_faces=int(np.sum(area == 0)), connected_surface_components=int(connected_components(graph, directed=False)[0]),
                  face_ids=list(map(int, np.unique(tags))), face_connected_components=face_components,
                  face_triangle_counts={str(int(tag)): int(np.sum(tags == tag)) for tag in np.unique(tags)},
                  ports=ports, wall_all_normals_outward=bool(np.all(wall_orientation > 0)),
                  signed_volume_m3=volume, cylinder_volume_relative_error=abs(volume / (np.pi * radius ** 2 * length) - 1),
                  wall_vertex_radius_max_error_m=float(np.max(abs(np.linalg.norm(points[wall_vertices, :2], axis=1) - radius))),
                  construction='Shared point IDs checked directly. No cleaning, tolerance welding or geometry mutation.')
    output['status'] = 'PASS' if not any(output[key] for key in ('free_edges', 'nonmanifold_edges', 'inconsistent_oriented_edges', 'global_duplicate_coordinate_count', 'zero_area_faces')) and output['connected_surface_components'] == 1 and output['wall_all_normals_outward'] and volume > 0 and all(value == 1 for value in face_components.values()) and output['face_ids'] == [1, 2, 3] and all(p['loop_vertex_count'] == p['shared_topology_vertex_count'] and p['loop_edge_count'] == p['shared_topology_edge_count'] and p['orientation'] == 'OUTWARD' for p in ports.values()) else 'FAIL'
    return output


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case_dir', type=Path)
    parser.add_argument('--surface', default='model.vtp')
    parser.add_argument('--output', default='SEAM_QA.json')
    args = parser.parse_args()
    config = json.loads((args.case_dir / 'INPUT.json').read_text())
    result = surface_qa(args.case_dir / args.surface, config['radius_m'], config['length_m'])
    (args.case_dir / args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))
