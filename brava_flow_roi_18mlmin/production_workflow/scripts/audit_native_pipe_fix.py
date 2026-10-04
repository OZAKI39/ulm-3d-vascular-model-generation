"""Check the user-defined native pipe mesh gate independently of flow accuracy."""
from pathlib import Path
import argparse
import hashlib
import json

import numpy as np
import pyvista as pv
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

from pipe_surface_qa import surface_qa

ROOT = Path(__file__).resolve().parents[1]


def main(case_name, mesh_dir):
    case = ROOT / 'straight_pipe' / case_name
    config = json.loads((case / 'INPUT.json').read_text())
    radius, length = config['radius_m'], config['length_m']
    original = surface_qa(case / 'model.vtp', radius, length)
    grid = pv.read(mesh_dir / 'volume.vtu')
    types, counts = np.unique(grid.celltypes, return_counts=True)
    assert set(types) == {10}, 'The native volume must contain linear tetrahedra'
    tetra = grid.cells_dict[10]
    native_point_dtype = grid.points.dtype
    points = np.asarray(grid.points, dtype=float)
    signed = np.linalg.det(points[tetra[:, 1:]] - points[tetra[:, :1]]) / 6
    faces = np.concatenate([tetra[:, face] for face in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]])
    owners = np.tile(np.arange(len(tetra)), 4)
    unique, first, inverse, multiplicity = np.unique(np.sort(faces, axis=1), axis=0, return_index=True, return_inverse=True, return_counts=True)
    order = np.argsort(inverse, kind='stable')
    starts = np.r_[0, np.cumsum(multiplicity)[:-1]]
    adjacent = owners[order[np.c_[starts[multiplicity == 2], starts[multiplicity == 2] + 1]]]
    graph = coo_matrix((np.ones(2 * len(adjacent)),
                       (np.r_[adjacent[:, 0], adjacent[:, 1]], np.r_[adjacent[:, 1], adjacent[:, 0]])),
                      shape=(len(tetra), len(tetra))).tocsr()
    components = int(connected_components(graph, directed=False)[0])
    boundary = faces[first[multiplicity == 1]].copy()
    boundary_owner = owners[first[multiplicity == 1]]
    triangle = points[boundary]
    normals = np.cross(triangle[:, 1] - triangle[:, 0], triangle[:, 2] - triangle[:, 0])
    flip = np.einsum('ni,ni->n', normals, triangle.mean(axis=1) - points[tetra[boundary_owner]].mean(axis=1)) < 0
    boundary[flip] = boundary[flip][:, [0, 2, 1]]
    mapping = {tuple(sorted(face)): i for i, face in enumerate(boundary)}
    tags = np.zeros(len(boundary), np.int32)
    max_error, errors = 0., []
    tree = cKDTree(points)
    for tag, name in [(1, 'WALL'), (2, 'INLET'), (3, 'OUTLET_01')]:
        face = pv.read(mesh_dir / (name + '.vtp'))
        distances, nodes = tree.query(face.points)
        max_error = max(max_error, float(distances.max()))
        for vertices in nodes[face.faces.reshape(-1, 4)[:, 1:]]:
            index = mapping.get(tuple(sorted(vertices)))
            if index is None:
                errors.append(name + ': not an exterior face')
            elif tags[index]:
                errors.append(name + ': duplicate assignment')
            else:
                tags[index] = tag
    # Reconstruct only an audit surface from tetra boundary; no change to the native grid.
    used_nodes, local_boundary = np.unique(boundary, return_inverse=True)
    exterior = pv.PolyData(points[used_nodes], np.c_[np.full(len(boundary), 3), local_boundary.reshape(-1, 3)].ravel())
    exterior.cell_data['ModelFaceID'] = tags
    exterior.save(mesh_dir / 'audited_exterior.vtp')
    surface = surface_qa(mesh_dir / 'audited_exterior.vtp', radius, length)
    wall_nodes = np.unique(boundary[tags == 1])
    allwall = np.all(np.isin(tetra, wall_nodes), axis=1)
    volume = float(signed.sum())
    relative_plc_volume = abs(volume / original['signed_volume_m3'] - 1)
    # The native exporter stores coordinates as Float32. Use two storage ULPs
    # at the largest geometric coordinate; do not apply a radius-only tolerance to z.
    tolerance_m = 2 * float(np.spacing(np.asarray(max(radius, length), dtype=native_point_dtype)))
    cap_relative = {name: abs(surface['ports'][name]['cap_area_m2'] / original['ports'][name]['cap_area_m2'] - 1)
                    for name in ('INLET', 'OUTLET')}
    checks = dict(nonempty=len(tetra) > 0, zero_volume_absent=bool(np.all(signed != 0)),
                  inverted_absent=bool(np.all(signed > 0)), connected_volume=components == 1,
                  manifold_volume=bool(np.all(multiplicity <= 2)), closed_oriented_labelled_surface=surface['status'] == 'PASS',
                  roles_preserved=not errors and bool(np.all(tags > 0)) and max_error < tolerance_m,
                  finite_volume=np.isfinite(volume).item() and volume > 0,
                  volume_matches_PLC=relative_plc_volume < 1e-6,
                  cap_areas_match_PLC=max(cap_relative.values()) < 1e-6,
                  radius_preserved=surface['wall_vertex_radius_max_error_m'] < tolerance_m,
                  ends_preserved=max(p['axial_coordinate_max_error_m'] for p in surface['ports'].values()) < tolerance_m)
    result = dict(status='PASS' if all(checks.values()) else 'FAIL', case=case_name, checks=checks,
                  points=len(points), tetra=len(tetra), native_cell_types=dict(zip(map(str, types), map(int, counts))),
                  zero_volume=int(np.sum(signed == 0)), negative_inverted=int(np.sum(signed < 0)),
                  volume_components=components, volume_m3=volume, PLC_volume_relative_error=relative_plc_volume,
                  circle_volume_relative_error=abs(volume / (np.pi * radius ** 2 * length) - 1),
                  cap_PLC_area_relative_errors=cap_relative, geometry_tolerance_m=tolerance_m,
                  native_point_dtype=str(native_point_dtype),
                  geometry_tolerance_basis='Two ULPs of native output coordinate storage at max(radius, length); no geometry modification or welding.',
                  boundary_mapping_max_error_m=max_error, boundary_mapping_errors=sorted(set(errors)),
                  all_WALL_tetra_count=int(allwall.sum()), surface=surface, original_surface=original,
                  note='This is the native topology/geometry gate only. Exact-cylinder near-wall probe coverage and CFD errors have not passed.',
                  volume_sha256=hashlib.sha256((mesh_dir / 'volume.vtu').read_bytes()).hexdigest())
    (mesh_dir / 'NATIVE_MESH_QA.json').write_text(json.dumps(result, indent=2) + '\n')
    np.savez_compressed(mesh_dir / 'BOUNDARY_TOPOLOGY.npz', points_m=points, tetra=tetra,
                        boundary_triangles=boundary, facet_tags=tags, adjacent_tetra=boundary_owner,
                        all_wall_tetra_ids=np.flatnonzero(allwall))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('case')
    parser.add_argument('mesh_dir', type=Path)
    args = parser.parse_args()
    main(args.case, args.mesh_dir)
