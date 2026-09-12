#!/usr/bin/python3
"""Frozen STL/contract preparation. Every generated artifact goes to RUN_DIR."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import yaml

EXPECTED_SHA = '840da5e1c43ec31ac70ec781cbf75b32940bc73538b5dba83c58af1ddf5e36fb'
PORT_IDS = ['inlet', 'outlet_01', 'outlet_02', 'outlet_03']

def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n')

def read_stl(path):
    raw = Path(path).read_bytes()
    count = int.from_bytes(raw[80:84], 'little')
    assert len(raw) == 84 + count*50, 'Only exact binary STL is accepted'
    return np.frombuffer(raw, offset=84, dtype=np.dtype([
        ('normal', '<f4', 3), ('v', '<f4', (3, 3)), ('attribute', '<u2')]))['v'].astype(float)

def projected_inside(points, triangles, normal, tolerance=1e-10):
    """Barycentric test of projection onto each actual triangle; no circle/hull."""
    # Remove normal coordinate using an orthonormal plane basis.
    seed = np.eye(3)[np.argmin(abs(normal))]
    u = np.cross(normal, seed); u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    basis = np.stack([u, v], axis=1)
    reference = triangles[0, 0]
    p = (points-reference) @ basis
    tri = (triangles-reference) @ basis
    hit = np.zeros(len(points), dtype=bool)
    for a, b, c in tri:
        e, f = b-a, c-a
        den = e[0]*f[1]-e[1]*f[0]
        if abs(den) < 1e-30:
            continue
        q = p-a
        s = (q[:, 0]*f[1]-q[:, 1]*f[0])/den
        t = (e[0]*q[:, 1]-e[1]*q[:, 0])/den
        hit |= (s >= -tolerance) & (t >= -tolerance) & (s+t <= 1+tolerance)
    return hit

def edge_components(triangles):
    _, inv = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
    faces = inv.reshape(-1, 3)
    owners = {}; rows = []; cols = []
    for i, face in enumerate(faces):
        for a, b in [(0,1),(1,2),(2,0)]:
            edge = tuple(sorted([int(face[a]), int(face[b])]))
            for j in owners.get(edge, []):
                rows.extend([i,j]); cols.extend([j,i])
            owners.setdefault(edge, []).append(i)
    graph = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(faces),len(faces)))
    return connected_components(graph, directed=False)

def prepare(contract, run, config):
    geometry = json.loads((contract/'metadata/geometry_contract.json').read_text())
    boundary = json.loads((contract/'metadata/boundary_contract.json').read_text())
    units = json.loads((contract/'metadata/unit_contract.json').read_text())
    stl = contract/geometry['packaged_stl_relative_path']
    digest = hashlib.sha256(stl.read_bytes()).hexdigest()
    assert digest == EXPECTED_SHA, 'STEP1_INPUT_INTEGRITY=FAIL'
    assert geometry['coordinate_unit'] == units['FINAL_STL_UNIT'] == 'm'
    ports = boundary['inlets'] + boundary['outlets']
    assert [p['id'] for p in ports] == PORT_IDS
    tri = read_stl(stl)
    bbox = np.stack([tri.min(axis=(0,1)), tri.max(axis=(0,1))])
    extents = bbox[1]-bbox[0]
    dx = float(yaml.safe_load(config.read_text())['mesh']['dx_m'])
    assert np.isfinite(dx) and dx > 0
    refdir = int(np.argmax(extents)); refn = int(np.floor(extents[refdir]/dx+0.5))
    effective = float(extents[refdir]/refn)
    # Fixed before voxelization: four float32 ULPs cover serialized plane noise.
    plane_tol = float(4*np.max(np.spacing(abs(tri).astype(np.float32))))
    normal_cos = float(np.cos(np.deg2rad(1.0)))
    vectors = np.cross(tri[:,1]-tri[:,0], tri[:,2]-tri[:,0])
    areas = np.linalg.norm(vectors,axis=1)/2
    unit_normals = vectors/(2*areas[:,None])
    volume = float(np.einsum('ij,ij->i', tri[:,0]-bbox[0], vectors).sum()/6)
    all_ids = []; records = {}; arrays = {}
    for label, port in enumerate(ports,1):
        assert port.get('center') and port.get('normal'), 'STEP2_PORT_MAPPING=BLOCKED'
        assert port['center_unit']=='m' and port['area_unit']=='m2'
        assert port['is_cap'] and port['is_artificial_extension']
        center = np.array(port['center']); normal = np.array(port['normal'])
        assert abs(np.linalg.norm(normal)-1)<1e-8
        distance = abs((tri-center) @ normal).max(axis=1)
        aligned = abs(unit_normals @ normal) >= normal_cos
        candidates = np.flatnonzero((distance <= plane_tol) & aligned)
        assert len(candidates), f'No planar cap candidates: {port["id"]}'
        num, components = edge_components(tri[candidates])
        ranked = []
        for c in range(num):
            ids = candidates[components==c]
            contains = bool(projected_inside(center[None,:], tri[ids], normal)[0])
            ac = np.average(tri[ids].mean(axis=1), weights=areas[ids], axis=0)
            ranked.append((not contains, float(np.linalg.norm(ac-center)), c, ids))
        ranked.sort(key=lambda r:r[:3]); ids = ranked[0][3]
        assert not ranked[0][0], f'Contract center not contained in recovered cap: {port["id"]}'
        area = float(areas[ids].sum()); rel = abs(area-port['area'])/port['area']
        mean = vectors[ids].sum(axis=0); mean /= np.linalg.norm(mean)
        dot = float(mean @ normal); outward_dot = float(np.sign(volume)*dot)
        relation = 'SAME_ORIENTATION' if dot>=normal_cos else 'OPPOSITE_ORIENTATION' if dot<=-normal_cos else 'WRONG_DIRECTION'
        record = dict(port, label=label, triangle_ids=ids.tolist(),
            candidate_triangle_count=len(candidates), candidate_components=num,
            reconstructed_triangle_count=len(ids), reconstructed_cap_triangle_area_m2=area,
            cap_area_relative_error=rel, cap_area_tolerance=0.01,
            cap_area_tolerance_scope='Step 2 geometry mapping only; old distal/proximal 5% threshold is a different metric',
            reconstructed_area_centroid_m=np.average(tri[ids].mean(axis=1),weights=areas[ids],axis=0).tolist(),
            max_plane_residual_m=float(distance[ids].max()), plane_tolerance_m=plane_tol,
            raw_triangle_mean_normal=mean.tolist(), normal_dot=dot,
            normal_orientation_relation=relation, outward_normal_dot=outward_dot,
            identification_pass=bool(rel<=0.01 and outward_dot>=normal_cos))
        records[port['id']] = record; arrays[port['id']] = ids
        all_ids.extend(ids.tolist())
    disjoint = len(set(all_ids))==len(all_ids)
    settings = dict(input_stl=str(stl),input_sha256=digest,input_unit='m',
        boundary_contract=str(contract/'metadata/boundary_contract.json'),
        smoke_dx_source=str(config)+':mesh.dx_m',smoke_dx_config_sha256=hashlib.sha256(config.read_bytes()).hexdigest(),
        requested_dx_m=dx,effective_dx_m_expected=effective,ref_dir=refdir,ref_dir_n=refn,
        physical_bbox_m=bbox.tolist(),requested_effective_relative_difference=(effective-dx)/dx,
        purpose='REFERENCE_GEOMETRY_SMOKE_DX; geometry only; not production/RBC resolution',
        full_stl_signed_volume_m3=volume,cap_plane_tolerance_m=plane_tol,
        cap_alignment_max_degrees=1.0,cap_area_relative_tolerance=0.01,
        opening_rule='One 6-neighbor link from closed lumen through a recovered cap triangle; exterior voxel projection stays inside actual cap triangle union',
        center_mapping_max_distance_lu=2.0,refinement_runs_allowed_only_if_disconnected=1)
    write_json(run/'run_settings.json', settings)
    write_json(run/'diagnostics/cap_identification.json',dict(ports=records,cap_triangle_sets_disjoint=disjoint))
    np.savez(run/'diagnostics/cap_triangle_ids.npz',**arrays)
    lines=['label\tid\tcx_m\tcy_m\tcz_m\tnx\tny\tnz\tarea_m2\tsource_id']
    for i,p in enumerate(ports,1):
        lines.append('\t'.join(map(str,[i,p['id'],*p['center'],*p['normal'],p['area'],p['source_id']])))
    (run/'ports.tsv').write_text('\n'.join(lines)+'\n')
    assert disjoint and all(p['identification_pass'] for p in records.values()), 'PORT_CAP_IDENTIFICATION=FAIL'
    print(json.dumps(settings,indent=2))
    for name,p in records.items():
        print(name,'triangles',p['reconstructed_triangle_count'],'area_rel_error',p['cap_area_relative_error'],'normal_dot',p['normal_dot'])

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('contract',type=Path);ap.add_argument('run',type=Path);ap.add_argument('cfd_config',type=Path)
    a=ap.parse_args();prepare(a.contract,a.run,a.cfd_config)
