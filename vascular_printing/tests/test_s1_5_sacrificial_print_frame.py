"""Persistent geometry, frozen input and local slicer safety regressions.

Run synthetic tests first with -k 'not real'. Real tests consume the independent
s1-5 integration artifacts and fail if these artifacts are absent.
"""
import copy
import json
import logging
from pathlib import Path
import zipfile

import numpy as np
import pytest
import trimesh

from vascular_processing import sacrificial_print_frame as frame
from vascular_processing import print_frame_qc as qc
from vascular_processing import support_access_qc as access
from vascular_processing import bambu_slice_adapter as bambu
from vascular_processing import sacrificial_fixture as legacy

ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi'
OUT=BASE/'print_frame_design'
EXPECTED_CORE_HASH='b9fb2d8a57b58b877f940f6c2486f923d93ac71b42793593f4092180dcb7e35e'


@pytest.fixture(scope='module')
def cfg():return frame.load_config(ROOT/'config/sacrificial_print_frame_BG001.yaml')


@pytest.fixture(scope='module')
def fixture(cfg):
    spec=json.loads((ROOT/'tests/data/simple_open_frame_fixture.json').read_text())
    nodes=np.array(spec['nodes'],dtype=float);r=spec['radius_mm']
    def tube(a,b):return trimesh.creation.cylinder(radius=r,sections=32,segment=np.array([a,b]))
    parts=[tube(nodes[a],nodes[b]) for a,b in spec['edges']]
    for p in nodes:
        sphere=trimesh.creation.icosphere(subdivisions=2,radius=r);sphere.apply_translation(p);parts.append(sphere)
    central=trimesh.boolean.union(parts,engine='manifold')
    inner=central.bounds+np.array([[-12]*3,[12]*3]);outer=inner+np.array([[-3,-3,-3],[3,3,0]])
    ports=[];full=[central]
    for p in spec['ports']:
        direction=legacy.NORMALS[p['face']];axis,side=legacy.face_axis(p['face'])
        target=nodes[p['node']].copy();target[axis]=inner[side,axis];end=target+13*direction
        full.append(tube(nodes[p['node']]-.3*direction,end))
        ports.append(dict(port_id=p['port_id'],face=p['face'],radius_mm=r,target=target,end=end,direction=direction))
    core=trimesh.boolean.union(full,engine='manifold')
    # Clean synthetic Boolean coincidences before ray tests. This fixture-only
    # operation is never applied to the accepted vascular input in production.
    core,_=legacy.union_core(core,[],cfg)
    return dict(core=core,central=central,ports=ports,old_box=dict(inner_bounds_mm=inner,outer_bounds_mm=outer))


@pytest.fixture(scope='module')
def sparse(fixture,cfg):return frame.build_frame('SPARSE_FRAME',0,fixture,cfg)


@pytest.fixture(scope='module')
def checked(fixture,sparse,cfg):return qc.frame_checks(sparse,fixture,frame.keep_zone_bounds(fixture['central'].bounds,cfg),cfg)[0]


@pytest.fixture(scope='module')
def combined(fixture,sparse,cfg):return frame.combine(fixture['core'],sparse['mesh'],cfg)[0]


def test_open_panel_generated(fixture,cfg):
    d=frame.build_frame('OPEN_PANEL',0,fixture,cfg)
    assert legacy.mesh_qc(d['mesh'],cfg)['passed']
    np.testing.assert_array_equal(d['outer_bounds'],fixture['old_box']['outer_bounds_mm'])
    assert sum(m['role']=='comparison_side_panel' for m in d['members'])==2


def test_sparse_frame_generated_without_panels(sparse):
    assert all('panel' not in m['role'] for m in sparse['members'])
    assert sum(m['role']=='external_port_anchor' for m in sparse['members'])==3


def test_opposite_faces_open(sparse,cfg):
    assert sparse['open_faces']==['-X','+X']
    for side in (0,1):assert qc.open_face_metrics(sparse,side,cfg)['passed']


def test_open_area_matches_analytic_rectangle(sparse,cfg):
    empty=copy.deepcopy(sparse);empty['members']=[]
    assert qc.open_face_metrics(empty,0,cfg)['free_area_fraction']==1
    axis=0;b=empty['inner_bounds'];lo=b[0].copy();hi=b[1].copy()
    lo[axis]=empty['outer_bounds'][0,axis];hi[axis]=b[0,axis]
    hi[1]=(lo[1]+hi[1])/2
    empty['members']=[frame.member('half','test',lo,hi)]
    assert qc.open_face_metrics(empty,0,cfg)['free_area_fraction']==pytest.approx(.5)


def test_sparse_one_component(sparse,cfg):assert legacy.mesh_qc(sparse['mesh'],cfg)['passed']
def test_union_one_component(combined,cfg):assert legacy.mesh_qc(combined,cfg)['passed']


def test_each_expected_port_overlap(checked):
    assert len(checked['port_connections'])==3
    assert all(r['status']=='PASS' and r['intersection_volume_mm3']>0 for r in checked['port_connections'])
    assert all(r['intersection_length_mm']==pytest.approx(4) for r in checked['port_connections'])


def test_no_mid_vessel_rod(checked):assert checked['unintended_mid_vessel_contact_mm3']<1e-6


def test_keep_zone_uses_central_bounds(fixture,cfg):
    keep=frame.keep_zone_bounds(fixture['central'].bounds,cfg)
    np.testing.assert_allclose(keep,fixture['central'].bounds+np.array([[-8]*3,[8]*3]))
    assert np.any(keep[1]<fixture['core'].bounds[1])


def test_exact_intrusion_volume():
    a=frame.box_mesh([[0,0,0],[4,4,4]]);b=frame.box_mesh([[2,2,2],[6,6,6]])
    assert qc.volume(qc.intersection(a,b))==pytest.approx(8)


def test_zero_intrusion_pass(checked):
    assert checked['frame_intrusion_keep_zone_mm3']==0
    assert checked['passed'],checked['failures']


def test_intrusion_rejected_without_auto_accept(fixture,sparse,cfg):
    bad=copy.deepcopy(sparse)
    intruder=frame.box_mesh([[-1,-1,-1],[1,1,1]])
    bad['mesh']=trimesh.boolean.union([bad['mesh'],intruder],engine='manifold')
    result,_,_=qc.frame_checks(bad,fixture,frame.keep_zone_bounds(fixture['central'].bounds,cfg),cfg)
    assert result['frame_intrusion_keep_zone_mm3']==pytest.approx(8)
    assert 'FRAME_INTRUDES_PDMS_KEEP_ZONE' in result['failures'] and not result['passed']


def test_visibility_ray_blocked_and_clear():
    obstacle=frame.box_mesh([[-1,-1,-1],[1,1,1]]);remote=obstacle.copy();remote.apply_translation([0,0,20])
    assert not access.line_visible(obstacle,remote,np.array([-5.,0,0]),np.array([5.,0,0]))
    assert access.line_visible(obstacle,remote,np.array([-5.,4,0]),np.array([5.,4,0]))


def test_capsule_center_and_collision():
    probe=access.capsule(np.array([0.,0,0]),np.array([10.,0,0]),2,16)
    np.testing.assert_allclose(probe.bounds[:,0],[-2,12],atol=.03)
    manager=trimesh.collision.CollisionManager();manager.add_object('frame',frame.box_mesh([[4,-1,-1],[6,1,1]]))
    assert manager.in_collision_single(probe)
    probe.apply_translation([0,6,0]);assert not manager.in_collision_single(probe)


def test_access_region_has_probe_route(fixture,sparse,cfg):
    regions=[dict(region_id='inlet',nearest_vascular_branch='inlet',center_mm=np.array([0.,-3.,0]),radius_mm=.6)]
    result,rows,paths=access.evaluate_access(sparse,fixture['core'],regions,cfg)
    assert result['accessible_region_count']==1
    assert result['probe_accessible_region_count']==1
    assert len(paths)==2 and rows[0]['tested_paths_A']>1


def test_trapped_free_space_proxy_detection(cfg):
    outer=np.array([[-10.]*3,[10.]*3]);members=[]
    for axis in range(3):
        for side in (0,1):
            lo=outer[0].copy();hi=outer[1].copy()
            if side==0:hi[axis]=-8
            else:lo[axis]=8
            members.append(frame.member(str((axis,side)),'closed_test_wall',lo,hi))
    core=trimesh.creation.icosphere(subdivisions=2,radius=.6)
    design=dict(outer_bounds=outer,inner_bounds=outer+np.array([[2]*3,[-2]*3]),members=members,open_axis=0,
                open_faces=['-X','+X'],mesh=frame.union_members(members,cfg))
    grid=access.free_space_components(design,core,cfg)
    label,exit_possible=access.point_component([4,0,0],grid)
    assert label!=0 and not exit_possible
    result,rows,_=access.evaluate_access(design,core,[dict(region_id='closed',nearest_vascular_branch='test',center_mm=np.zeros(3),radius_mm=.6)],cfg)
    assert result['trapped_proxy_region_count']==1 and rows[0]['status']=='TRAPPED_SUPPORT_PROXY'


def test_orientations_deterministic(cfg):
    a=qc.orientation_rotations(cfg);b=qc.orientation_rotations(cfg)
    assert len(a)==40
    for (name,m),(name2,m2) in zip(a,b):
        assert name==name2;np.testing.assert_array_equal(m,m2)
        np.testing.assert_allclose(m@m.T,np.eye(3),atol=1e-12);assert np.linalg.det(m)==pytest.approx(1)


def test_build_volume_gate(combined,sparse,cfg):
    proxy=dict(trapped_proxy_region_count=0,accessible_fraction=1.)
    with pytest.warns(UserWarning,match='Gimbal lock'):
        rows,top=qc.rank_orientations(combined,sparse['mesh'],[10,10,10],proxy,cfg)
    assert top==[] and not any(r['fits_build_volume'] for r in rows)


def test_bambu_unavailable_graceful(monkeypatch,cfg,tmp_path):
    monkeypatch.setattr(bambu.existing,'discover',lambda _:dict(status='BAMBU_STUDIO_NOT_FOUND'))
    d=bambu.discover_current(cfg,tmp_path)
    result=bambu.slice_top_five(d,cfg,[],tmp_path,logging.getLogger('test'))
    assert result['status']=='BAMBU_STUDIO_NOT_FOUND' and result['results']==[]


@pytest.mark.parametrize('flag',['--send','--upload','--start-print','--print','--orient','--scale'])
def test_no_printer_send_or_unapproved_transform(flag):
    with pytest.raises(ValueError):bambu.check_local_command(['bambu-studio',flag,'1'])


def test_local_slice_command_allowed():assert bambu.check_local_command(['bambu-studio','--slice','0','--export-3mf','local.3mf'])


def test_actual_gcode_parser_tracks_support_extrusion(tmp_path):
    path=tmp_path/'support.3mf'
    with zipfile.ZipFile(path,'w') as archive:
        archive.writestr('Metadata/plate_1.gcode','G90\nM83\nG1 X1 Y2 Z3\n; FEATURE: Support\nG1 X5 E1\nG1 Y6 E-1\n; FEATURE: Outer wall\nG1 X8 E1\n')
    report,segments=bambu.support_paths_from_archive(path)
    assert report['support_segment_count']==1
    np.testing.assert_array_equal(segments,[[[1,2,3],[5,2,3]]])
    assert not report['support_geometry_available']


@pytest.mark.parametrize('clockwise',[True,False])
def test_arc_direction_and_length(clockwise):
    start=np.array([1.,0,3.]);end=np.array([0.,-1 if clockwise else 1,3.])
    points=bambu.arc_polyline(start,end,dict(I=-1.,J=0.),clockwise)
    assert np.all(points[:,1]<=1e-10) if clockwise else np.all(points[:,1]>=-1e-10)
    np.testing.assert_allclose(np.linalg.norm(points[:,:2],axis=1),1,atol=1e-12)
    assert np.linalg.norm(np.diff(points,axis=0),axis=1).sum()==pytest.approx(np.pi/2,rel=.01)


def test_support_arcs_update_following_line_origin(tmp_path):
    path=tmp_path/'arc.3mf'
    with zipfile.ZipFile(path,'w') as archive:
        archive.writestr('Metadata/plate_1.gcode','G90\nM83\nG1 X1 Y0 Z3\n; FEATURE: Support\nG3 X0 Y1 I-1 J0 E1\nG1 X0 Y2 E.5\n')
    report,segments=bambu.support_paths_from_archive(path)
    assert report['support_arc_move_count']==1 and report['support_linear_move_count']==1
    assert report['support_extrusion_filament_mm']==pytest.approx(1.5)
    np.testing.assert_array_equal(segments[-1],[[0,1,3],[0,2,3]])


def test_stl_roundtrip(combined,cfg,tmp_path):
    report=legacy.export_checked_stl(combined,tmp_path/'roundtrip.stl',cfg)
    assert report['passed']


def test_source_volume_preserved(combined,fixture,cfg):
    removed=trimesh.boolean.difference([fixture['core'],combined],engine='manifold')
    assert qc.volume(removed)<=cfg['geometry']['volume_tolerance_mm3']


@pytest.fixture(scope='module')
def real():
    path=OUT/'sacrificial_print_frame_summary.json'
    assert path.is_file(),'Run s1-5_sacrificial_print_frame.py before real integration checks'
    return json.loads(path.read_text())


def test_real_source_hash_unchanged(real):
    assert legacy.sha256(Path(real['source_path']))==EXPECTED_CORE_HASH
    assert legacy.sha256(Path(real['vascular_only_file']))==EXPECTED_CORE_HASH


def test_real_interface_qc_unchanged(real):
    assert {r['port_id'] for r in real['interface_qc']}=={'I1','O1','O2','O3'}
    assert all(r['status']=='PASS' and r['before_geometry_hash']==r['after_geometry_hash'] for r in real['interface_qc'])


@pytest.mark.parametrize('name',['s1-2_swc_roi_generate_human.py','s1-3_swc_roi_generate_MeVO.py','s1-4_sacrificial_box_and_ports.py'])
def test_real_frozen_ui_unchanged(real,name):
    protection=json.loads((OUT/'protected_source_hashes.json').read_text())
    key=str(ROOT/name)
    from vascular_processing.project_paths import frozen_hash
    assert frozen_hash(protection['before_sha256'],key)==legacy.sha256(ROOT/name)==frozen_hash(protection['after_sha256'],key)


def test_real_all_protected_sources_unchanged(real):
    before=json.loads((OUT/'protected_before.json').read_text())
    from vascular_processing.project_paths import verify_migrated_snapshot
    assert verify_migrated_snapshot(before)['passed']


def test_real_mesh_and_connections(real,cfg):
    assert legacy.mesh_qc(trimesh.load_mesh(real['final_stl']),cfg)['passed']
    assert real['source_volume_loss_mm3']<=cfg['geometry']['volume_tolerance_mm3']
    assert len(real['frame']['port_connections'])==4 and all(r['status']=='PASS' for r in real['frame']['port_connections'])
    assert real['pdms_keep_zone']['frame_intrusion_mm3']==0


def test_real_visual_evidence(real):
    from vascular_processing.print_frame_review import NAMES
    for name in NAMES:assert (OUT/'QC'/(name+'.png')).stat().st_size>10000


def test_real_slicing_only_top_five_and_local(real):
    assert real['orientation']['candidate_count']==40
    assert len(real['orientation']['top_five'])==5
    records=real['bambu']['results'];assert records
    allowed={r['candidate_id'] for r in real['orientation']['top_five']}
    for record in records:
        assert record['orientation_id'] in allowed
        assert bambu.check_local_command(record['command'])
        if record['status']=='BAMBU_SLICED':assert record['geometry_qc']['rotation_scale_position_preserved']
    assert not real['printer_job_sent']


def test_real_support_preview_includes_circular_paths(real):
    record=real['bambu']['selected'];audit=record['support_archive_audit']
    assert audit['support_arc_move_count']>0 and audit['support_linear_move_count']>0
    assert len(audit['mesh_objects'])==1 and audit['model_part_types']==['normal_part']
    with zipfile.ZipFile(record['archive']) as archive:
        lines=archive.read('Metadata/plate_1.gcode').decode().splitlines()
    feature='';expected=0
    for line in lines:
        if line.startswith('; FEATURE:'):feature=line.split(':',1)[1]
        if line.startswith(('G2 ','G3 ')) and 'support' in feature.lower() and ' E' in line:expected+=1
    assert audit['support_arc_move_count']==expected
    paths=np.load(record['support_paths_file'])
    assert len(paths)==audit['support_segment_count'] and np.isfinite(paths).all()
