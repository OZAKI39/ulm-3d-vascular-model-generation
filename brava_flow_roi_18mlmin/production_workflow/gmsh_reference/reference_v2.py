"""One authorized cap bisection and V2; exact V1 BREP and wall shells stay frozen."""
from pathlib import Path
import argparse
import copy
import datetime
import hashlib
import json
import math
import sys
import time

import numpy as np
import pyvista as pv

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path[:0] = [str(ROOT/'runtime/gmsh'), str(ROOT/'scripts')]
from audit_pipe_fidelity import boundary_profile


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def put(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def configure(cap_size_m, cap_only=False, core_radial_intervals=3):
    import gmsh
    old = json.loads((HERE/'GMSH_REFERENCE_MESH_PLAN.json').read_text())
    identity = json.loads((HERE/'candidate_0/IDENTITY.json').read_text())
    assert sha(HERE/'candidate_0/reference.brep') == identity['cad_sha256']
    assert gmsh.__version__ == '4.15.2'
    assert sha(gmsh.libpath) == identity['gmsh_binary_sha256']
    gmsh.model.add('GMSH_R1_REFERENCE_V2')
    for key, value in old['options'].items():
        gmsh.option.setNumber(key, .200 if key == 'Mesh.MeshSizeMax' else value)
    gmsh.model.occ.importShapes(str(HERE/'candidate_0/reference.brep'))
    gmsh.model.occ.synchronize()
    length = old['length_m']*1000
    volumes = gmsh.model.getEntities(3)
    assert len(volumes) == 16
    assert core_radial_intervals in (3, 6)
    curves = {}
    core_curves = []
    for _, tag in gmsh.model.getEntities(1):
        kind = gmsh.model.getType(1,tag)
        boundary = gmsh.model.getBoundary([(1,tag)],oriented=False)
        xyz = np.array([gmsh.model.getValue(0,t,[]) for _,t in boundary])
        if kind == 'Circle':
            intervals = 68
        elif kind == 'Line' and np.ptp(xyz[:,2]) > .9*length:
            intervals = 30
        else:
            distance = np.linalg.norm(xyz[1]-xyz[0])
            if distance < .03:
                intervals = 1
            else:
                assert kind == 'Line' and abs(distance-old['shell_radii_m'][-1]*1000)<1e-9
                intervals = core_radial_intervals
                core_curves.append(tag)
        gmsh.model.mesh.setTransfiniteCurve(tag, intervals+1)
        curves[str(tag)] = dict(type=kind,intervals=intervals)
    assert len(core_curves) == 8
    roles = dict(WALL=[],INLET=[],OUTLET=[],INTERNAL=[])
    surfaces = []
    transition_faces = []
    for _,tag in gmsh.model.getEntities(2):
        up,_ = gmsh.model.getAdjacencies(2,tag)
        kind = gmsh.model.getType(2,tag)
        edges = gmsh.model.getBoundary([(2,tag)],oriented=False)
        z = gmsh.model.occ.getCenterOfMass(2,tag)[2]
        if len(up)==2: role='INTERNAL'
        elif kind=='Cylinder': role='WALL'
        elif abs(z)<1e-9: role='INLET'
        elif abs(z-length)<1e-9: role='OUTLET'
        else: raise ValueError(('unclassified boundary',tag,kind,z))
        roles[role].append(tag)
        core_transition = (core_radial_intervals == 6 and role == 'INTERNAL' and kind == 'Plane'
                           and any(t in core_curves for _,t in edges))
        if core_transition:
            # These four axial core partition faces must not extrude six radial
            # rows throughout the pipe. Their boundaries stay fixed; the native
            # 200-um Frontal-Delaunay mesh supplies a local end transition.
            transition_faces.append(tag)
        transfinite = len(edges)==4 and not core_transition
        if transfinite:
            gmsh.model.mesh.setTransfiniteSurface(tag,arrangement='AlternateLeft')
        surfaces.append(dict(tag=tag,type=kind,role=role,adjacent_volumes=list(map(int,up)),
                             curves=[x[1] for x in edges],transfinite=transfinite))
    assert len(transition_faces) == (4 if core_radial_intervals == 6 else 0)
    for key,name in enumerate(['WALL','INLET','OUTLET'],1):
        gmsh.model.addPhysicalGroup(2,roles[name],key,name=name)
    gmsh.model.addPhysicalGroup(3,[v[1] for v in volumes],1,name='FLUID')
    background = gmsh.model.mesh.field.add('MathEval')
    gmsh.model.mesh.field.setString(background,'F','.200')
    cap = gmsh.model.mesh.field.add('MathEval')
    gmsh.model.mesh.field.setString(cap,'F',format(cap_size_m*1000,'.17g'))
    restricted = gmsh.model.mesh.field.add('Restrict')
    gmsh.model.mesh.field.setNumber(restricted,'InField',cap)
    gmsh.model.mesh.field.setNumber(restricted,'IncludeBoundary',0)
    gmsh.model.mesh.field.setNumber(restricted,'IncludeEmbedded',0)
    gmsh.model.mesh.field.setNumbers(restricted,'SurfacesList',roles['INLET']+roles['OUTLET'])
    minimum = gmsh.model.mesh.field.add('Min')
    gmsh.model.mesh.field.setNumbers(minimum,'FieldsList',[background,restricted])
    gmsh.model.mesh.field.setAsBackgroundMesh(minimum)
    if cap_only:
        gmsh.model.setVisibility(gmsh.model.getEntities(),0)
        gmsh.model.setVisibility([(2,t) for t in roles['INLET']+roles['OUTLET']],1,recursive=True)
        gmsh.option.setNumber('Mesh.MeshOnlyVisible',1)
    else:
        gmsh.option.setNumber('Mesh.MeshOnlyVisible',0)
    return dict(volumes=[v[1] for v in volumes],surfaces=surfaces,roles=roles,curves=curves,
                core_radial_curves=core_curves,core_radial_intervals=core_radial_intervals,
                unstructured_core_transition_faces=transition_faces,
                reused_BREP_sha256=identity['cad_sha256'])


def nodes():
    import gmsh
    tags,xyz,_ = gmsh.model.mesh.getNodes()
    order = np.argsort(tags)
    return np.asarray(tags)[order],np.asarray(xyz).reshape(-1,3)[order]/1000


def cap_poly(tags, points, entities):
    import gmsh
    triangles = []
    for entity in entities:
        types,_,ids = gmsh.model.mesh.getElements(2,entity)
        assert list(types)==[2]
        triangles.append(np.searchsorted(tags,np.asarray(ids[0]).reshape(-1,3)))
    triangles = np.concatenate(triangles)
    used,local = np.unique(triangles,return_inverse=True)
    return pv.PolyData(points[used],np.c_[np.full(len(triangles),3),local.reshape(-1,3)].ravel())


def cap_gate(poly,radius):
    result,velocity = boundary_profile(poly,radius)
    center = int(np.argmin(np.linalg.norm(poly.points[:,:2],axis=1)))
    result.update(center_radius_m=float(np.linalg.norm(poly.points[center,:2])),
                  center_velocity_m_s=float(velocity[center]),
                  center_pass=bool(np.linalg.norm(poly.points[center,:2])<1e-12 and abs(velocity[center]-.2)<1e-12))
    result['PASS'] = bool(result['PASS'] and result['center_pass'] and result['rim_analytic_coordinate_roundoff_m_s']<1e-12)
    return result


def calibrate():
    import gmsh
    output = ROOT/'CAP_TRACE_CALIBRATION.json'
    assert not output.exists(), 'Do not repeat the bounded calibration'
    old = json.loads((HERE/'GMSH_REFERENCE_MESH_PLAN.json').read_text())
    baseline = json.loads((ROOT/'GMSH_REFERENCE_MESH_RESULT.json').read_text())
    source = HERE/'candidate_0/mesh/INLET.vtp'
    poly = pv.read(source)
    xyz = np.asarray(poly.points)
    tri = poly.faces.reshape(-1,4)[:,1:]
    edges = np.unique(np.sort(np.concatenate([tri[:,[0,1]],tri[:,[1,2]],tri[:,[2,0]]]),axis=1),axis=0)
    h = np.linalg.norm(xyz[edges[:,1]]-xyz[edges[:,0]],axis=1)
    axis = (abs(xyz[:,1])<1e-12)&(xyz[:,0]>=-1e-12)
    old_error = baseline['cap_analytical_trace_diagnostic']['INLET']
    ratio = min(.01/old_error['profile_relative_L2'],.005/old_error['Q_relative_error'])
    proposed = old['bulk_size_m']*math.sqrt(ratio)*.9
    cap_size = math.floor(proposed/10e-6)*10e-6
    result = dict(status='RUNNING',VALIDATION_ONLY=True,volume_generations=0,
                  baseline=dict(cap_sha256=sha(source),nodes=len(xyz),triangles=len(tri),
                                edge_length_quantiles_m=np.quantile(h,[0,.05,.5,.95,1]).tolist(),
                                positive_x_axis_radii_m=np.sort(xyz[axis,0]).tolist(),
                                trace=old_error),
                  sizing_rule='h_new=h_old*sqrt(min(0.01/e_L2,0.005/e_Q))*0.9, rounded down to 10 um; max 3 cap-only trials; stop on first PASS',
                  first_unrounded_size_m=proposed,cap_policy='Restrict scalar size to end-surface interiors, IncludeBoundary=0; all transfinite curves and four-edge faces unchanged',
                  trials=[],maximum_trials=3,created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
    put(output,result)
    for trial in range(3):
        start=time.monotonic()
        gmsh.initialize([],readConfigFiles=False)
        try:
            identity=configure(cap_size,cap_only=True)
            gmsh.model.mesh.generate(2)
            assert sum(len(x) for x in gmsh.model.mesh.getElements(3)[1])==0
            tags,points=nodes()
            caps={name:cap_gate(cap_poly(tags,points,identity['roles'][name]),old['radius_m']) for name in ['INLET','OUTLET']}
            row=dict(trial=trial,cap_size_m=cap_size,elapsed_s=time.monotonic()-start,caps=caps,
                     PASS=all(c['PASS'] for c in caps.values()),volume_elements=0)
            result['trials'].append(row)
            print(json.dumps(row),flush=True)
        finally:
            gmsh.finalize()
        if row['PASS']:
            result.update(status='CAP_TRACE_CALIBRATION_PASS',selected_cap_size_m=cap_size,
                          selected_trial=trial,selection='First passing sequential trial, not an exhaustive minimum-cost proof')
            put(output,result)
            return
        e_u=max(c['profile_relative_L2'] for c in caps.values())
        e_q=max(c['Q_relative_error'] for c in caps.values())
        next_size=math.floor(cap_size*math.sqrt(min(.01/e_u,.005/e_q))*.9/10e-6)*10e-6
        cap_size=max(10e-6,min(cap_size-10e-6,next_size))
        put(output,result)
    result['status']='CAP_TRACE_CALIBRATION_FAIL'
    put(output,result)


def correct_cap():
    """Exactly one new 2D trial. A durable receipt prevents accidental replay."""
    import gmsh
    output = ROOT/'CAP_TOPOLOGY_CORRECTION.json'
    assert not output.exists(), 'The single authorized topology trial is already consumed'
    old = json.loads((HERE/'GMSH_REFERENCE_MESH_PLAN.json').read_text())
    result = dict(status='RUNNING', VALIDATION_ONLY=True, trials_authorized=1,
                  cap_size_m=30e-6, core_radial_intervals_before=3,core_radial_intervals=6,
                  core_radial_spacing_m=old['shell_radii_m'][-1]/6,
                  volume_generations=0, source_BREP_sha256=sha(HERE/'candidate_0/reference.brep'),
                  generator_sha256=sha(__file__),gmsh_library_sha256=sha(gmsh.libpath),
                  cap_policy='Fixed 30 um surface-interior target; eight end-core radial curves deterministically bisected from 3 to 6 intervals; identical inlet/outlet.',
                  transition_policy='Only four internal axial core-plane faces use unstructured Algorithm 6 with 200 um bulk target. All wall/shell surfaces and axial/circumferential curve counts remain frozen; no six-row extrusion through the core.',
                  created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
    put(output,result)
    start=time.monotonic()
    gmsh.initialize([],readConfigFiles=False)
    gmsh.logger.start()
    try:
        identity=configure(30e-6,cap_only=True,core_radial_intervals=6)
        gmsh.model.mesh.generate(2)
        assert sum(len(x) for x in gmsh.model.mesh.getElements(3)[1])==0
        tags,points=nodes()
        arrays={}
        caps={}
        for name in ['INLET','OUTLET']:
            poly=cap_poly(tags,points,identity['roles'][name])
            caps[name]=cap_gate(poly,old['radius_m'])
            arrays[name+'_points_m']=poly.points
            arrays[name+'_triangles']=poly.faces.reshape(-1,4)[:,1:]
        raw=HERE/'CAP_TOPOLOGY_CORRECTION_MESH.npz'
        np.savez_compressed(raw,**arrays)
        result.update(status='CAP_TOPOLOGY_CORRECTION_PASS' if all(c['PASS'] for c in caps.values()) else 'CAP_TOPOLOGY_CORRECTION_FAIL',
                      caps=caps,identity=identity,elapsed_s=time.monotonic()-start,
                      retained_cap_mesh_sha256=sha(raw),completed_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
        put(output,result)
        print(json.dumps({k:v for k,v in result.items() if k!='identity'}),flush=True)
    except Exception as exc:
        result.update(status='CAP_TOPOLOGY_CORRECTION_FAIL',error=repr(exc),elapsed_s=time.monotonic()-start)
        put(output,result)
        raise
    finally:
        (HERE/'cap_topology_correction.log').write_text('\n'.join(gmsh.logger.get())+'\n')
        gmsh.finalize()


def freeze_plan():
    target=HERE/'GMSH_REFERENCE_MESH_V2_PLAN.json'
    assert not target.exists()
    old=json.loads((HERE/'GMSH_REFERENCE_MESH_PLAN.json').read_text())
    calibration=json.loads((ROOT/'CAP_TOPOLOGY_CORRECTION.json').read_text())
    assert calibration['status']=='CAP_TOPOLOGY_CORRECTION_PASS'
    assert calibration['generator_sha256']==sha(__file__)
    plan=copy.deepcopy(old)
    plan.update(schema='GMSH_REFERENCE_V2_PLAN',created_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                bulk_size_m=200e-6,cap_size_m=30e-6,core_radial_intervals=6,
                bulk_metric='Central z/L=[0.4,0.6], centroid d_wall>=150 um; median of arithmetic mean of six tet edges',
                cap_policy=calibration['cap_policy'],cap_calibration_sha256=sha(ROOT/'CAP_TRACE_CALIBRATION.json'),
                cap_topology_correction_sha256=sha(ROOT/'CAP_TOPOLOGY_CORRECTION.json'),
                transition_policy=calibration['transition_policy'],
                unstructured_core_transition_faces=calibration['identity']['unstructured_core_transition_faces'],
                V1_plan_sha256=sha(HERE/'GMSH_REFERENCE_MESH_PLAN.json'),
                source_BREP_sha256=sha(HERE/'candidate_0/reference.brep'),
                generator_sha256=sha(__file__),maximum_initial_candidates=1,maximum_deterministic_corrections=0,
                correction_condition='The user authorizes exactly one V2, no V3 or further correction.',
                bulk_size_reason='210*229.390453/236.894625 um, with user-frozen safety margin to 200 um.',
                physical_groups=dict(WALL=1,INLET=2,OUTLET=3,FLUID=1),
                gates=dict(first_order_tet_only=True,zero_negative_duplicate_tetra=0,fluid_components=1,
                           free_nonmanifold_exterior_edges=0,internal_face_neighbors=2,duplicate_nodes_interfaces=0,
                           central_facet_deviation_m=5e-8,cap_plane_tolerance_m=1e-12,
                           wall_altitude_median_m=[20e-6,30e-6],central_bulk_ratio=[.75,1.25],
                           cap_profile_L2=.01,cap_Q_relative=.005,cap_center_velocity_m_s=.2,rim_tolerance_m_s=1e-12,
                           strict_frozen_probes=960,central_wall_grid_unchanged=True),
                cap_scope='End-core curves and surface interiors; no cap-size volume field extension; identical central wall grid and shell radii',
                surface_meshing='V1 transfinite wall/shell surfaces unchanged. Four internal core-plane interfaces use Algorithm 6 and bulk 200 um, allowing six end-core intervals to transition locally. End core sectors use Algorithm 6, cap 30 um.',
                field_gates=dict(particle_bands_m=[.377e-6,1e-6,2e-6,5e-6,10e-6],points_per_band=96,
                                 u_L2=.05,G_Frobenius=.08,s_L2=.08,curl_L2=.08,bulk_u_L2=.05,Q_relative=.05,dpdz_relative=.08,
                                 WSS_and_0p1um='diagnostic only'),
                trajectory_gates=dict(speed_history_L2=.05,axial_displacement_relative=.05,transit_residence_relative=.05,
                                      extra_radial_drift_m=.1e-6,wall_penetration_allowed=False,
                                      comparison='Predeclared common physical times, never step-index alignment'),
                timestep_gates=dict(dt_s=[.0005,.00025,.000125],primary_pair_s=[.00025,.000125],
                                    speed_history_L2=.02,axial_displacement_relative=.02,transit_residence_relative=.02,
                                    matching_contact_terminal_and_minimum_gap_categories=True),
                mesher_decision='GMSH_VALIDATION_REFERENCE_SELECTED')
    plan['options']['Mesh.MeshSizeMax']=.200
    put(target,plan)
    (target.with_suffix('.sha256')).write_text(sha(target)+'\n')
    print('FROZEN_V2_PLAN',sha(target),flush=True)


def generate():
    import gmsh
    path=HERE/'GMSH_REFERENCE_MESH_V2_PLAN.json'
    assert sha(path)==path.with_suffix('.sha256').read_text().strip()
    plan=json.loads(path.read_text())
    assert sha(__file__)==plan['generator_sha256'], 'Generator changed after plan freeze'
    out=HERE/'candidate_v2'
    out.mkdir(exist_ok=False)
    start,cpu=time.monotonic(),time.process_time()
    gmsh.initialize([],readConfigFiles=False)
    try:
        identity=configure(plan['cap_size_m'],core_radial_intervals=plan['core_radial_intervals'])
        assert identity['unstructured_core_transition_faces']==plan['unstructured_core_transition_faces']
        # Copy exact source BREP; re-export would only change serialization/order.
        (out/'reference.brep').write_bytes((HERE/'candidate_0/reference.brep').read_bytes())
        identity.update(VALIDATION_ONLY=True,plan_sha256=sha(path),generator_sha256=sha(__file__),
                        gmsh_version=gmsh.__version__,gmsh_binary_sha256=sha(gmsh.libpath),
                        cad_sha256=sha(out/'reference.brep'),CAD_units='mm',export_units='m')
        put(out/'IDENTITY.json',identity)
        gmsh.logger.start()
        began=time.monotonic()
        gmsh.model.mesh.generate(3)
        generation_seconds=time.monotonic()-began
        gmsh.write(str(out/'reference.msh'))
        tags,points=nodes()
        types,element_ids,element_nodes=gmsh.model.mesh.getElements(3)
        assert list(types)==[4]
        tetra=np.searchsorted(tags,np.asarray(element_nodes[0]).reshape(-1,4))
        triangles,entity,role=[],[],[]
        for surface in identity['surfaces']:
            types2,_,ids=gmsh.model.mesh.getElements(2,surface['tag'])
            assert list(types2)==[2]
            faces=np.searchsorted(tags,np.asarray(ids[0]).reshape(-1,3))
            triangles.append(faces);entity.extend([surface['tag']]*len(faces))
            role.extend([dict(WALL=1,INLET=2,OUTLET=3,INTERNAL=0)[surface['role']]]*len(faces))
        np.savez_compressed(out/'RAW_MESH.npz',points_m=points,tetra=tetra,gmsh_node_ids=tags,
                            gmsh_tet_ids=element_ids[0],cad_triangles=np.concatenate(triangles),cad_entity=entity,cad_role=role)
        identity.update(mesh_generation_seconds=generation_seconds,elapsed_seconds=time.monotonic()-start,
                        process_CPU_seconds=time.process_time()-cpu,
                        average_CPU_workers=(time.process_time()-cpu)/(time.monotonic()-start),
                        volume_element_properties=[gmsh.model.mesh.getElementProperties(4)[:4]],
                        node_count=len(points),tetra_count=len(tetra),mesh_sha256=sha(out/'reference.msh'),
                        arrays_sha256=sha(out/'RAW_MESH.npz'))
        put(out/'IDENTITY.json',identity)
        (out/'gmsh.log').write_text('\n'.join(gmsh.logger.get())+'\n')
        print(json.dumps({k:identity[k] for k in ['mesh_generation_seconds','node_count','tetra_count','mesh_sha256']}),flush=True)
    finally:
        gmsh.finalize()


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('stage',choices=['correct-cap','freeze-plan','generate'])
    args=parser.parse_args()
    {'correct-cap':correct_cap,'freeze-plan':freeze_plan,'generate':generate}[args.stage]()
