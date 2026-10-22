"""One frozen OCC/HXT validation fixture. No legacy mesh input or point snapping."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'runtime/gmsh'))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def build(plan, out, cad_only=False):
    import gmsh
    start = time.monotonic()
    cpu_start = time.process_time()
    gmsh.initialize([], readConfigFiles=False)
    try:
        gmsh.model.add('GMSH_R1_REFERENCE_MESH')
        for key, value in plan['options'].items():
            gmsh.option.setNumber(key, value)
        # OCC works in mm here; export applies only the exact unit conversion.
        R, L = plan['radius_m'] * 1000, plan['length_m'] * 1000
        radii = np.array(plan['shell_radii_m']) * 1000
        solids = []
        for quadrant in range(4):
            for radius in radii:
                tag = gmsh.model.occ.addCylinder(0, 0, 0, 0, 0, L, float(radius), angle=math.pi/2)
                gmsh.model.occ.rotate([(3, tag)], 0, 0, 0, 0, 0, 1, quadrant * math.pi/2)
                solids.append((3, tag))
        gmsh.model.occ.fragment(solids[:1], solids[1:], removeObject=True, removeTool=True)
        gmsh.model.occ.synchronize()
        volumes = gmsh.model.getEntities(3)
        assert len(volumes) == 16, volumes
        curve_counts = {}
        for _, tag in gmsh.model.getEntities(1):
            kind = gmsh.model.getType(1, tag)
            endpoints = gmsh.model.getBoundary([(1, tag)], oriented=False)
            p = np.array([gmsh.model.getValue(0, v, []) for _, v in endpoints])
            if kind == 'Circle':
                intervals = plan['circumferential_intervals'] // 4
            elif kind == 'Line' and np.ptp(p[:, 2]) > L * .9:
                intervals = plan['axial_intervals']
            else:
                distance = np.linalg.norm(p[1] - p[0])
                intervals = 1 if distance < .03 else math.ceil(distance / (plan['bulk_size_m'] * 1000))
            gmsh.model.mesh.setTransfiniteCurve(tag, intervals + 1)
            curve_counts[str(tag)] = dict(type=kind, intervals=intervals)
        roles = {'WALL': [], 'INLET': [], 'OUTLET': [], 'INTERNAL': []}
        surfaces = []
        for _, tag in gmsh.model.getEntities(2):
            up, _ = gmsh.model.getAdjacencies(2, tag)
            kind = gmsh.model.getType(2, tag)
            curves = gmsh.model.getBoundary([(2, tag)], oriented=False)
            center = np.array(gmsh.model.occ.getCenterOfMass(2, tag))
            if len(up) == 2:
                role = 'INTERNAL'
            elif kind == 'Cylinder':
                role = 'WALL'
            elif abs(center[2]) < 1e-9:
                role = 'INLET'
            elif abs(center[2] - L) < 1e-9:
                role = 'OUTLET'
            else:
                raise ValueError(('unclassified surface', tag, kind, center))
            roles[role].append(tag)
            if len(curves) == 4:
                gmsh.model.mesh.setTransfiniteSurface(tag, arrangement='AlternateLeft')
            surfaces.append(dict(tag=tag, type=kind, role=role, adjacent_volumes=list(map(int, up)),
                                 curves=[c[1] for c in curves], transfinite=len(curves) == 4))
        for key, name in enumerate(['WALL', 'INLET', 'OUTLET'], 1):
            gmsh.model.addPhysicalGroup(2, roles[name], key, name=name)
        gmsh.model.addPhysicalGroup(3, [v[1] for v in volumes], 1, name='FLUID')
        field = gmsh.model.mesh.field.add('MathEval')
        gmsh.model.mesh.field.setString(field, 'F', format(plan['bulk_size_m'] * 1000, '.17g'))
        gmsh.model.mesh.field.setAsBackgroundMesh(field)
        out.mkdir(parents=True, exist_ok=False)
        gmsh.write(str(out / 'reference.brep'))
        identity = dict(VALIDATION_ONLY=True, plan_sha256=sha(HERE / 'GMSH_REFERENCE_MESH_PLAN.json'),
                        generator_sha256=sha(__file__), gmsh_version=gmsh.__version__,
                        gmsh_binary_sha256=sha(gmsh.libpath), CAD_units='mm', export_units='m',
                        geometry='Exact OCC quarter cylinders, shared BooleanFragments; 3 radial shells and core',
                        volumes=[v[1] for v in volumes], surfaces=surfaces, roles=roles, curves=curve_counts,
                        cad_sha256=sha(out / 'reference.brep'))
        (out / 'IDENTITY.json').write_text(json.dumps(identity, indent=2) + '\n')
        if cad_only:
            print(json.dumps({'CAD_ONLY': True, 'volume_count': len(volumes), 'roles': roles}), flush=True)
            return
        gmsh.logger.start()
        mesh_start = time.monotonic()
        gmsh.model.mesh.generate(3)
        mesh_seconds = time.monotonic() - mesh_start
        gmsh.write(str(out / 'reference.msh'))
        nodes, xyz, _ = gmsh.model.mesh.getNodes()
        sort = np.argsort(nodes)
        nodes = np.asarray(nodes)[sort]
        points = np.asarray(xyz).reshape(-1, 3)[sort] / 1000
        types, element_ids, element_nodes = gmsh.model.mesh.getElements(3)
        type_info = [gmsh.model.mesh.getElementProperties(int(t))[:4] for t in types]
        assert list(types) == [4], ('volume must be first-order tet', types, type_info)
        tetra = np.searchsorted(nodes, np.asarray(element_nodes[0]).reshape(-1, 4))
        all_surfaces, entity, role_ids = [], [], []
        role_number = dict(WALL=1, INLET=2, OUTLET=3, INTERNAL=0)
        for surf in surfaces:
            stypes, _, snodes = gmsh.model.mesh.getElements(2, surf['tag'])
            assert list(stypes) == [2], ('surface must be linear triangles', surf, stypes)
            faces = np.searchsorted(nodes, np.asarray(snodes[0]).reshape(-1, 3))
            all_surfaces.append(faces)
            entity.extend([surf['tag']] * len(faces))
            role_ids.extend([role_number[surf['role']]] * len(faces))
        np.savez_compressed(out / 'RAW_MESH.npz', points_m=points, tetra=tetra,
                            gmsh_node_ids=nodes, gmsh_tet_ids=element_ids[0],
                            cad_triangles=np.concatenate(all_surfaces), cad_entity=entity, cad_role=role_ids)
        identity.update(mesh_generation_seconds=mesh_seconds, elapsed_seconds=time.monotonic()-start,
                        process_CPU_seconds=time.process_time()-cpu_start,
                        average_CPU_workers=(time.process_time()-cpu_start)/(time.monotonic()-start),
                        volume_element_properties=type_info, node_count=len(points), tetra_count=len(tetra),
                        mesh_sha256=sha(out/'reference.msh'), arrays_sha256=sha(out/'RAW_MESH.npz'))
        (out/'IDENTITY.json').write_text(json.dumps(identity, indent=2)+'\n')
        (out/'gmsh.log').write_text('\n'.join(gmsh.logger.get())+'\n')
        print(json.dumps(identity), flush=True)
    finally:
        gmsh.finalize()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cad-only', action='store_true')
    args = parser.parse_args()
    build(json.loads((HERE/'GMSH_REFERENCE_MESH_PLAN.json').read_text()),
          HERE/('cad_preflight' if args.cad_only else 'candidate_0'), args.cad_only)
