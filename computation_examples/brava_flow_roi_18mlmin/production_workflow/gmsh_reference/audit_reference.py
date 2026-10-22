"""Read-only scientific gates and lossless Gmsh -> existing linear-tetra interface."""
from pathlib import Path
import hashlib
import itertools
import json
import sys
import argparse

import numpy as np
import pyvista as pv
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'scripts'))
from pipe_surface_qa import surface_qa
from pipe_probe_geometry import frozen_probes, strict_cells
from audit_pipe_fidelity import boundary_profile


def stats(x):
    return dict(zip(['min', 'P5', 'median', 'P95', 'max'], map(float, np.quantile(x, [0, .05, .5, .95, 1]))))


def finite_facet_deviation(xyz, radius):
    """Exact max radial deviation on each planar triangle, not centroid sagitta."""
    xy = xyz[..., :2]
    minima = np.full(len(xy), np.inf)
    for i, j in [(0, 1), (1, 2), (2, 0)]:
        a, delta = xy[:, i], xy[:, j]-xy[:, i]
        denom = np.sum(delta**2, axis=1)
        t = np.divide(-np.sum(a*delta, axis=1), denom, out=np.zeros(len(a)), where=denom > 0)
        minima = np.minimum(minima, np.linalg.norm(a + np.clip(t, 0, 1)[:, None]*delta, axis=1))
    cross = lambda a, b: a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0]
    signs = np.stack([cross(xy[:, j]-xy[:, i], -xy[:, i]) for i,j in [(0,1),(1,2),(2,0)]], axis=1)
    inside = (np.all(signs > 0, axis=1) | np.all(signs < 0, axis=1))
    minima[inside] = 0
    return np.maximum(np.abs(radius-minima), np.max(np.abs(np.linalg.norm(xy, axis=2)-radius), axis=1))


def main(v2=False):
    plan_path = HERE/('GMSH_REFERENCE_MESH_V2_PLAN.json' if v2 else 'GMSH_REFERENCE_MESH_PLAN.json')
    plan = json.loads(plan_path.read_text())
    out = HERE/('candidate_v2' if v2 else 'candidate_0')
    identity = json.loads((out/'IDENTITY.json').read_text())
    radius, length = plan['radius_m'], plan['length_m']
    with np.load(out/'RAW_MESH.npz') as data:
        points, tetra = data['points_m'], data['tetra']
        cad, cad_role = data['cad_triangles'], data['cad_role']
        gmsh_nodes, gmsh_tets = data['gmsh_node_ids'], data['gmsh_tet_ids']
    signed = np.linalg.det(points[tetra[:, 1:]]-points[tetra[:, :1]])/6
    faces = np.concatenate([tetra[:, f] for f in [(0,1,2),(0,1,3),(0,2,3),(1,2,3)]])
    owners = np.tile(np.arange(len(tetra)), 4)
    keys, first, inv, count = np.unique(np.sort(faces, axis=1), axis=0, return_index=True, return_inverse=True, return_counts=True)
    order = np.argsort(inv, kind='stable')
    starts = np.r_[0, np.cumsum(count)[:-1]]
    shared = owners[order[np.c_[starts[count==2], starts[count==2]+1]]]
    graph = coo_matrix((np.ones(2*len(shared)),
                       (np.r_[shared[:,0],shared[:,1]], np.r_[shared[:,1],shared[:,0]])),
                      shape=(len(tetra),len(tetra))).tocsr()
    components = int(connected_components(graph, directed=False)[0])
    lookup = {tuple(k):i for i,k in enumerate(keys)}
    cad_keys = np.sort(cad, axis=1)
    cad_ids = np.array([lookup.get(tuple(k), -1) for k in cad_keys])
    matched = cad_ids >= 0
    expected = np.where(cad_role > 0, 1, 2)
    conforming = bool(np.all(matched) and np.all(count[cad_ids] == expected))
    cad_duplicate = len(cad_keys)-len(np.unique(cad_keys, axis=0))
    boundary = faces[first[count == 1]].copy()
    adjacent = owners[first[count == 1]]
    xyz = points[boundary]
    normals = np.cross(xyz[:,1]-xyz[:,0], xyz[:,2]-xyz[:,0])
    flip = np.einsum('ni,ni->n', normals, xyz.mean(1)-points[tetra[adjacent]].mean(1)) < 0
    boundary[flip] = boundary[flip][:,[0,2,1]]
    role_lookup = {tuple(k):int(v) for k,v in zip(cad_keys[cad_role > 0], cad_role[cad_role > 0])}
    tags = np.array([role_lookup.get(tuple(sorted(k)), 0) for k in boundary], np.int32)
    mesh = out/'mesh'
    mesh.mkdir(exist_ok=False)
    grid = pv.UnstructuredGrid(np.c_[np.full(len(tetra), 4), tetra].ravel(), np.full(len(tetra), 10, np.uint8), points)
    grid.point_data['GlobalNodeID'] = np.arange(1,len(points)+1, dtype=np.int32)
    grid.cell_data['GlobalElementID'] = np.arange(1,len(tetra)+1, dtype=np.int32)
    grid.point_data['GmshNodeID'] = gmsh_nodes
    grid.cell_data['GmshElementID'] = gmsh_tets
    grid.save(mesh/'volume.vtu')
    def surface(selected, path):
        node, index = np.unique(boundary[selected], return_inverse=True)
        poly = pv.PolyData(points[node], np.c_[np.full(sum(selected), 3), index.reshape(-1,3)].ravel())
        poly.point_data['GlobalNodeID'] = (node+1).astype(np.int32)
        poly.cell_data['GlobalElementID'] = (adjacent[selected]+1).astype(np.int32)
        poly.cell_data['ModelFaceID'] = tags[selected]
        poly.save(path)
        return poly
    surface(np.ones(len(boundary), dtype=bool), mesh/'audited_exterior.vtp')
    caps = {}
    for tag, name in [(1,'WALL'),(2,'INLET'),(3,'OUTLET_01')]:
        poly = surface(tags == tag, mesh/(name+'.vtp'))
        if tag != 1:
            caps[name], velocity = boundary_profile(poly, radius)
            if v2:
                center_id = int(np.argmin(np.linalg.norm(poly.points[:,:2],axis=1)))
                caps[name].update(center_radius_m=float(np.linalg.norm(poly.points[center_id,:2])),
                                 center_velocity_m_s=float(velocity[center_id]))
                caps[name]['center_pass'] = bool(caps[name]['center_radius_m']<1e-12 and abs(velocity[center_id]-.2)<1e-12)
    ext = surface_qa(mesh/'audited_exterior.vtp', radius, length)
    wall = points[boundary[tags == 1]]
    center = wall.mean(1)
    central = (center[:,2]>=.4*length) & (center[:,2]<=.6*length)
    wall_grid_comparison = None
    if v2:
        with np.load(HERE/'candidate_0/mesh/BOUNDARY_TOPOLOGY.npz') as previous:
            old_wall = previous['points_m'][previous['boundary_triangles'][previous['facet_tags']==1]]
        old_center = old_wall.mean(1)
        old_wall = old_wall[(old_center[:,2]>=.4*length)&(old_center[:,2]<=.6*length)]
        def coordinate_rows(triangles):
            # Comparison only: native mesh points are never rounded or moved.
            vertices = np.round(triangles,12)
            order = np.lexsort((vertices[:,:,2],vertices[:,:,1],vertices[:,:,0]),axis=1)
            rows = np.take_along_axis(vertices,order[:,:,None],axis=1).reshape(-1,9)
            return rows[np.lexsort(rows[:,::-1].T[::-1])]
        a,b = coordinate_rows(old_wall),coordinate_rows(wall[central])
        same = a.shape==b.shape and np.array_equal(a,b)
        wall_grid_comparison = dict(PASS=bool(same),old_facets=len(a),new_facets=len(b),
                                    comparison_roundoff_m=1e-12,
                                    definition='Coordinate triangle sets in central slab, independent of node IDs; comparison rounding only, no mesh modification')
    n = np.cross(wall[:,1]-wall[:,0],wall[:,2]-wall[:,0])
    n /= np.linalg.norm(n, axis=1)[:,None]
    altitude = np.abs(np.einsum('nvi,ni->nv',points[tetra[adjacent[tags == 1]]]-center[:,None,:],n)).max(1)
    # Include every triangle intersecting the central slab for a conservative
    # finite-facet bound; these triangles span fixed straight axial intervals.
    intersects = (wall[:,:,2].max(1)>=.4*length) & (wall[:,:,2].min(1)<=.6*length)
    deviation = finite_facet_deviation(wall[intersects], radius)
    edges = np.array(list(itertools.combinations(range(4),2)))
    distance = radius-np.linalg.norm(points[tetra].mean(1)[:,:2], axis=1)
    bulk_ids = np.flatnonzero(distance>=150e-6)
    bulk = np.linalg.norm(points[tetra[bulk_ids][:,edges[:,0]]]-points[tetra[bulk_ids][:,edges[:,1]]],axis=2).mean(1)
    z_bulk = points[tetra[bulk_ids]].mean(1)[:,2]
    central_bulk = bulk[(z_bulk>=.4*length)&(z_bulk<=.6*length)]
    central_bulk_stats = stats(central_bulk)
    height_stats, bulk_stats = stats(altitude[central]), stats(bulk)
    _, nc = np.unique(points, axis=0, return_counts=True)
    wall_nodes = np.unique(boundary[tags==1])
    allwall = np.all(np.isin(tetra, wall_nodes), axis=1)
    ratio = (central_bulk_stats if v2 else bulk_stats)['median']/plan['bulk_metric_P0_m']
    checks = dict(nonempty=len(tetra)>0, first_order_tetra_only=identity['volume_element_properties']==[['Tetrahedron 4',3,1,4]],
                  zero_volume_absent=bool(np.all(signed != 0)), inverted_absent=bool(np.all(signed > 0)),
                  duplicate_tetra_absent=len(np.unique(np.sort(tetra,axis=1),axis=0))==len(tetra),
                  connected_volume=components==1, manifold_volume=bool(np.all(count<=2)),
                  closed_oriented_labelled_surface=ext['status']=='PASS',
                  roles_preserved=set(tags)=={1,2,3}, internal_interfaces_conforming=conforming,
                  duplicate_CAD_faces_absent=cad_duplicate==0, duplicate_nodes_absent=bool(np.all(nc==1)),
                  finite_volume=bool(np.isfinite(points).all() and np.isfinite(signed).all() and signed.sum()>0),
                  ends_preserved=max(p['axial_coordinate_max_error_m'] for p in ext['ports'].values())<=1e-12,
                  finite_facet_geometry=float(deviation.max())<=.05e-6,
                  wall_altitude=20e-6<=height_stats['median']<=30e-6,
                  bulk_scale=.75<=ratio<=1.25)
    if v2:
        checks['cap_trace'] = all(c['PASS'] and c['center_pass'] and c['rim_analytic_coordinate_roundoff_m_s']<1e-12 for c in caps.values())
        checks['central_wall_grid_unchanged'] = wall_grid_comparison['PASS']
    proto = json.loads((HERE.parent/'manufactured_pipe/PROTOCOL.json').read_text())
    levels, probes, groups, planes, angles = frozen_probes(proto,radius,length)
    probe_cells, _ = strict_cells(grid,points,tetra,probes)
    if v2:
        checks['strict_frozen_probe_coverage'] = bool(np.all(probe_cells>=0))
    np.savez_compressed(mesh/'BOUNDARY_TOPOLOGY.npz',points_m=points,tetra=tetra,boundary_triangles=boundary,
                        facet_tags=tags,adjacent_tetra=adjacent,all_wall_tetra_ids=np.flatnonzero(allwall))
    np.savez_compressed(out/'QA_SAMPLES.npz', central_altitude_m=altitude[central],
                        central_wall_centers_m=center[central],bulk_mean_edge_m=bulk,
                        central_bulk_mean_edge_m=central_bulk,
                        central_facet_deviation_m=deviation, probes_m=probes,probe_cell=probe_cells)
    result = dict(status='GMSH_REFERENCE_MESH_PASS' if all(checks.values()) else 'REFERENCE_MESH_FAIL',
                  VALIDATION_ONLY=True, checks=checks, points=len(points),tetra=len(tetra),
                  zero_volume=int(sum(signed==0)), negative_inverted=int(sum(signed<0)),
                  duplicate_tetra=len(tetra)-len(np.unique(np.sort(tetra,axis=1),axis=0)),
                  volume_components=components, surface=ext, all_WALL_tetra_count=int(sum(allwall)),
                  internal_CAD_triangles=int(sum(cad_role==0)), CAD_face_multiplicity_pass=conforming,
                  duplicate_CAD_triangles=cad_duplicate,
                  wall_altitude_m=height_stats, bulk_mean_edge_m=bulk_stats, bulk_P0_ratio=ratio,
                  historical_matched_central_bulk=dict(z_over_L=[.4,.6], wall_distance_min_m=150e-6,
                      count=len(central_bulk), mean_edge_m=central_bulk_stats,
                      P0_ratio=central_bulk_stats['median']/plan['bulk_metric_P0_m'],
                      gate_pass=.75<=central_bulk_stats['median']/plan['bulk_metric_P0_m']<=1.25),
                  bulk_scope_note='Full-length statistics follow the new plan; the frozen P0 183.512362 um reference was central z/L=0.4-0.6. The separately reported central comparison matches its population. Neither can be replaced by nominal size.',
                  central_finite_facet_deviation_max_m=float(deviation.max()),
                  strict_frozen_probes=dict(total=len(probes),inside=int(sum(probe_cells>=0))),
                  cap_analytical_trace_diagnostic=caps, identity=identity,
                  volume_sha256=hashlib.sha256((mesh/'volume.vtu').read_bytes()).hexdigest(),
                  correction_permitted=bool(all(v for k,v in checks.items() if k not in ['wall_altitude','bulk_scale'])
                                            and not checks['wall_altitude'] and 15e-6<=height_stats['median']<=35e-6),
                  physics_or_BraVa_modified=False)
    if v2:
        result.update(status='GMSH_REFERENCE_MESH_V2_PASS' if all(checks.values()) else 'REFERENCE_MESH_FAIL_FINAL',
                      central_wall_grid_comparison=wall_grid_comparison,
                      correction_permitted=False,bulk_scope_note='Gate uses central z/L=0.4-0.6 and centroid d_wall>=150 um, identical to P0. Full-length bulk_mean_edge_m is diagnostic only.',
                      QA_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    result_path = HERE.parent/'GMSH_REFERENCE_MESH_V2_RESULT.json' if v2 else HERE/'GMSH_REFERENCE_MESH_RESULT.json'
    result_path.write_text(json.dumps(result,indent=2)+'\n')
    (mesh/'NATIVE_MESH_QA.json').write_text(json.dumps(dict(result,status='PASS' if all(checks.values()) else 'FAIL'),indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['identity','surface']},indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--v2',action='store_true')
    main(parser.parse_args().v2)
