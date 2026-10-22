"""Permanent algorithm regressions plus independent, read-only real-run assertions."""
import copy
import json
from pathlib import Path
from types import SimpleNamespace
import zipfile
import networkx as nx
import numpy as np
import pytest
import trimesh
from vascular_processing import roi_support_metrics as m
from vascular_processing import support_aware_roi as s
from vascular_processing import support_aware_port_design as p
from vascular_processing import support_aware_evaluation as e
from vascular_processing import sacrificial_fixture as f
from vascular_processing import casting_mold_qc as q
from vascular_processing import abs_casting_mold as mold
from vascular_processing.bambu_slice_adapter import check_local_command
from vascular_processing.port_attachment_alignment import align_known_cap
from vascular_processing.surface_continuity_qc import port_continuity

FIXTURE=json.loads((Path(__file__).parent/'data/support_aware_roi_synthetic.json').read_text())


@pytest.fixture
def graph():
    g=nx.DiGraph()
    for n,*xyzr in FIXTURE['nodes']:g.add_node(n,coords=np.array(xyzr,float),swc_type=4)
    g.add_edges_from(FIXTURE['edges']);return g


@pytest.fixture
def cfg():return p.configuration()


@pytest.fixture(scope='module')
def real():
    if not (s.DEFAULT_OUTPUT/'model_results.json').exists():pytest.skip('Run s1-7 real BG001 stages first')
    return s.DEFAULT_OUTPUT


def test_hash_protection_detects_change(tmp_path):
    path=tmp_path/'source.swc';path.write_text('frozen');snapshot={str(path):s.sha256(path)}
    assert s.verify_snapshot(snapshot)['all_unchanged'];path.write_text('changed')
    assert s.verify_snapshot(snapshot)['changed']==[str(path)]


def test_source_discovery_is_frozen_d075(real):
    d=json.loads((real/'source_discovery.json').read_text())
    assert d['subject']=='BG001' and d['side']=='RMCA' and not d['semantic_recomputed']
    assert all('d075' in r['path'] for r in d['sources'])


@pytest.mark.parametrize('label,expected',[('MeVO',True),('M1',False),('UNKNOWN',False)])
def test_mevo_gate_never_guesses(graph,label,expected):
    assert s.mevo_gate(graph,{edge:'branch' for edge in graph.edges},{'branch':{'label':label}})==expected


@pytest.mark.parametrize('factor',[.9,1.000001,1.1,2.])
def test_scaling_rejected(factor):
    transform=np.eye(4);transform[0,0]=factor
    with pytest.raises(ValueError,match='SCALING'):m.rigid_transform(transform)


def test_rigid_translation_preserves_dimensions():
    t=np.eye(4);t[:3,3]=[100,20,-8];m.rigid_transform(t)
    points=np.array([[0,0,0],[60,50,10.]])
    assert m.dimensions(points)['xyz_extents_mm']==m.dimensions(points+t[:3,3])['xyz_extents_mm']


@pytest.mark.parametrize('case',FIXTURE['dimension_cases'])
def test_size_and_xy_swap(case):assert m.dimensions(case['points'])['size_target_met']==case['pass']


def test_connected_single_inlet_three_outlets(graph):
    r=m.topology_metrics(graph)
    assert r['topology_valid'] and r['inlet_count']==1 and r['outlet_count']==3 and r['connected_components']==1
    assert r['bifurcation_count']==2 and r['branch_count']==5


def test_disconnected_subtree_rejected(graph):
    graph.remove_edge(4,9);assert not m.topology_metrics(graph)['topology_valid']


def test_root_at_bifurcation_rejected(graph):
    graph.remove_node(1);assert not m.topology_metrics(graph)['topology_valid']


def test_short_terminal_filtered(graph):
    graph.nodes[5]['coords'][:3]=graph.nodes[4]['coords'][:3]+[.1,0,0]
    assert m.topology_metrics(graph)['short_terminal_count']==1 and not m.topology_metrics(graph)['topology_valid']


def test_two_outlets_allowed_but_one_not(graph):
    graph.remove_nodes_from([9,10,11]);assert m.topology_metrics(graph)['topology_valid']
    graph.remove_nodes_from([6,7,8]);assert not m.topology_metrics(graph)['topology_valid']


def test_candidate_determinism_and_native_coordinates(graph):
    data=dict(transform=np.eye(4),sources=[dict(graph=graph,id='synthetic',component=4)],
        edge_branch={edge:'b' for edge in graph.edges},branches={'b':dict(label='MeVO')},
        profile=dict(effective=dict(radius_floor_mm=.6),manufacturing=dict(radius_blend_length_mm=5)))
    first=s.generate_candidates(data,minimum=1);second=s.generate_candidates(data,minimum=1)
    assert first==second and len(first)>=20
    assert all(set(c['native_node_ids'])<=set(graph) and c['scale_factor']==1. and c['topology_valid'] for c in first)


def test_surrogate_distinguishes_horizontal_vertical(graph):
    horizontal=m.centerline_support(graph,np.eye(4));g=copy.deepcopy(graph)
    for n in g:g.nodes[n]['coords'][[0,2]]=g.nodes[n]['coords'][[2,0]]
    assert m.centerline_support(g,np.eye(4))['estimated_support_contact_mm2']<horizontal['estimated_support_contact_mm2']
    assert max(r['branch_depth'] for r in horizontal['segments'])==2


def test_thin_and_bifurcation_penalties_separate(graph):
    r=m.centerline_support(graph,np.eye(4));assert r['thin_branch_support_penalty']>0 and r['support_near_bifurcation_count']>0
    for n in graph:graph.nodes[n]['coords'][3]=1.
    assert m.centerline_support(graph,np.eye(4))['thin_branch_support_penalty']==0


def test_triangle_area_not_vertex_count(graph):
    # A 10x20 box has 200 mm² horizontal underside, independent of subdivision.
    mesh=trimesh.creation.box([10,20,2]);r,_=m.surface_support(mesh,graph,np.eye(4))
    dense=mesh.subdivide();rd,_=m.surface_support(dense,graph,np.eye(4))
    assert r['overhang_area_mm2']==pytest.approx(200) and rd['overhang_area_mm2']==pytest.approx(200)


def test_mesh_validity_and_open_mesh_fail(cfg):
    mesh=trimesh.creation.box();assert f.mesh_qc(mesh,cfg)['passed']
    mesh.update_faces(np.arange(len(mesh.faces)-1));r=f.mesh_qc(mesh,cfg)
    assert not r['watertight'] and not r['manifold'] and not r['passed']


def test_fcl_triangle_gap_not_vertex_distance():
    a=trimesh.creation.box([2,2,2]);b=trimesh.creation.box([2,2,2]);b.apply_translation([3.6,.4,.3])
    distance,_,_=q.fcl_distance(a,b);assert distance==pytest.approx(1.6,abs=1e-8)


@pytest.fixture
def cap_route(cfg):
    # Axial rings give the independent seam audit local adjacent triangles;
    # two ten-mm triangles would provide no samples in its +/-1mm window.
    theta=np.arange(24)*2*np.pi/24;xs=np.linspace(-5,5,41)
    vertices=[[x,.6*np.cos(a),.6*np.sin(a)] for x in xs for a in theta];faces=[]
    for j in range(len(xs)-1):
        for k in range(24):
            a=j*24+k;b=j*24+(k+1)%24;c=(j+1)*24+k;d=(j+1)*24+(k+1)%24
            faces.extend([[a,b,c],[b,d,c]])
    for j in (0,len(xs)-1):
        center=len(vertices);vertices.append([xs[j],0,0])
        faces.extend([[center,j*24+k,j*24+(k+1)%24] for k in range(24)])
    source=trimesh.Trimesh(vertices,faces,process=True);source.fix_normals()
    endpoint=f.Endpoint('I1','inlet',1,1,'B1',np.array([5.,0,0]),.6,np.array([1.,0,0]),5,5,.6,.6)
    align_known_cap(endpoint,source,dict(vertices=24),cfg)
    box=f.compute_box(source.bounds,cfg);target=endpoint.position.copy();target[0]=box.inner[1,0]
    route=f.make_route(endpoint,'+X',target,box,cfg,'ROUTE_A_STRAIGHT')
    return source,endpoint,box,route


def test_actual_cap_24_and_buried_collar(cap_route):
    _,endpoint,_,_=cap_route;q=endpoint.attachment_qc
    assert q['cap_vertex_count']==24 and q['collar_contained'] and not q['full_radius_backward_overlap_used']
    assert q['production_profile']=='ACTUAL_CAP_POLYGON'


def test_port_interface_and_section(cap_route,cfg):
    source,endpoint,_,route=cap_route;core,_=f.union_core(source,[route],cfg)
    checks=p.interface_qc(core,source,[endpoint],[route],cfg)
    assert checks[0]['feature_edge_20_length_mm']==0 and checks[0]['section_area_change_fraction']<.01 and checks[0]['passed']


def test_local_cap_precision_preserves_world_geometry(cap_route,cfg):
    source,endpoint,_,_=cap_route;shift=np.array([150.,120.,30.])
    source.apply_translation(shift);endpoint.position+=shift;before=source.vertices.copy()
    record=p.align_cap_local(endpoint,source,dict(vertices=24),cfg)
    assert record['collar_outside_source_volume_mm3']<=1e-6 and record['cap_vertex_count']==24
    assert not record['tolerances_changed'] and np.array_equal(source.vertices,before)
    assert endpoint.position==pytest.approx([155,120,30],abs=1e-9)


def test_new_layout_can_use_opposite_outlet_walls(cfg):
    options=[]
    for i,face in enumerate(['-Y','-X','+X']):
        mesh=trimesh.creation.box();mesh.apply_translation([10*i,0,0])
        options.append([SimpleNamespace(face=face,mesh=mesh,target=np.array([10*i,0,0.]),
            endpoint=SimpleNamespace(radius=.6),vessel_clearance=4.,length=10.,turn=20.,option_id=str(i))])
    routes,audit=p.select_new_layout(options,cfg,3.)
    assert len(routes)==3 and audit['status']=='PASS'


def test_shorter_maximum_extension_precedes_smaller_sum(cfg,monkeypatch):
    def route(name,length,face):
        return SimpleNamespace(face=face,mesh=name,target=np.array([0.,0,0]),endpoint=SimpleNamespace(radius=.6),
            vessel_clearance=4.,length=length,turn=0.,option_id=name)
    allowed={frozenset(('a10','b40')),frozenset(('a25','b30'))}
    monkeypatch.setattr(f,'surface_distance',lambda a,b:4. if frozenset((a,b)) in allowed else .5)
    rows,audit=p.select_new_layout([[route('a10',10,'-X'),route('a25',25,'-X')],
        [route('b30',30,'+Y'),route('b40',40,'+Y')]],cfg,3.)
    assert [r.option_id for r in rows]==['a25','b30']
    assert audit['maximum_extension_mm']==30 and audit['total_length_mm']==55


def test_feature_edge_gate_rejects_sharp_ring():
    continuity=dict(features=[dict(angle_deg=20,total_length_mm=1.,circumferential_ring_count=1)])
    assert not m.interface_gate(continuity,dict(passed=True),0)


def test_new_box_open_top_closed_sides_and_bottom(cfg):
    bounds=np.array([[0,0,0],[60,50,10.]])
    box,_=p.new_box(bounds,cfg);check=q.intact_box_qc(box,cfg)
    assert check['top_fully_open'] and check['four_side_walls_closed'] and check['bottom_closed']
    assert check['inner_dimensions_mm']==pytest.approx([84,74,32])
    assert check['outer_dimensions_mm']==pytest.approx([90,80,35])
    assert check['wall_thickness_mm']==3 and check['bottom_thickness_mm']==3


def test_adaptive_unconnected_faces_fit_geometry_without_moving_ports(cfg):
    central=trimesh.creation.box([60,50,10]);before=central.vertices.copy()
    seed,_=p.new_box(central.bounds,cfg)
    ports=[dict(face='+X',target=np.array([seed['inner'][1,0],0.,0.]))]
    box=p.adaptive_box(central,ports,seed,cfg);faces={r['face']:r for r in box['adaptive_faces']}
    assert np.array_equal(central.vertices,before)
    assert faces['+X']['adaptive_plane_mm']==faces['+X']['original_plane_mm']
    for name in ('-X','-Y','+Y','BOTTOM'):
        assert faces[name]['core_envelope_gap_mm']==pytest.approx(5.)
        assert faces[name]['fcl_surface_distance_mm']==pytest.approx(5.)
    assert faces['TOP_OPENING']['core_envelope_gap_mm']==pytest.approx(5.)
    qc=q.intact_box_qc(box,cfg)
    assert qc['four_side_walls_closed'] and qc['bottom_closed'] and qc['top_fully_open']
    assert qc['wall_thickness_mm']==3 and qc['bottom_thickness_mm']==3


def test_adaptive_box_preserves_port_edge_material(cfg):
    central=trimesh.creation.box([10,10,2]);seed,_=p.new_box(central.bounds,cfg)
    ports=[dict(face='+X',target=np.array([seed['inner'][1,0],7.,0.]))]
    box=p.adaptive_box(central,ports,seed,cfg);faces={r['face']:r for r in box['adaptive_faces']}
    assert faces['+Y']['reason']=='PORT_EDGE_MATERIAL_LIMIT'
    assert box['inner'][1,1]-ports[0]['target'][1]>=cfg['ports']['min_wall_edge_clearance_mm']
    assert all(r.get('absolute_clearance_pass',True) for r in faces.values())


def test_changed_box_retires_dependent_geometry_and_slices(tmp_path,cfg):
    (tmp_path/'vascular_core.stl').write_bytes(b'fixed vascular geometry')
    p.ensure_current_design(tmp_path,np.eye(4),cfg)
    (tmp_path/'manufacturing').mkdir();(tmp_path/'orientations').mkdir()
    (tmp_path/'evaluation.json').write_text('{}')
    p.ensure_current_design(tmp_path,np.eye(4),cfg)
    assert (tmp_path/'evaluation.json').exists() and not (tmp_path/'history').exists()
    changed=copy.deepcopy(cfg);changed['box']['unconnected_face_gap_mm']=4.
    p.ensure_current_design(tmp_path,np.eye(4),changed)
    assert not (tmp_path/'evaluation.json').exists()
    assert (tmp_path/'history/revision_01/evaluation.json').exists()
    assert (tmp_path/'history/revision_01/orientations').is_dir()


def test_route_feedback_is_bound_to_exact_vascular_source(tmp_path,cfg):
    source=tmp_path/'vascular_core.stl';source.write_bytes(b'fixed source')
    review=dict(vascular_sha256=s.sha256(source),allowed_faces_by_endpoint={'O2':['-X']})
    (tmp_path/'port_route_review.json').write_text(json.dumps(review))
    current=p.candidate_configuration(tmp_path,cfg)
    assert current['route_review']==review and 'route_review' not in cfg
    source.write_bytes(b'different source')
    with pytest.raises(ValueError,match='SOURCE_MISMATCH'):p.candidate_configuration(tmp_path,cfg)


def test_route_feedback_cannot_enable_z_wall(tmp_path,cfg):
    source=tmp_path/'vascular_core.stl';source.write_bytes(b'fixed source')
    (tmp_path/'port_route_review.json').write_text(json.dumps(dict(vascular_sha256=s.sha256(source),allowed_faces_by_endpoint={'O2':['+Z']})))
    with pytest.raises(ValueError,match='INVALID_REVIEW_SIDEWALL'):p.candidate_configuration(tmp_path,cfg)


def test_top_access_feedback_only_changes_affected_artificial_port():
    geometry=dict(ports=[dict(port_id='O2',face='+Y'),dict(port_id='I1',face='-Y')])
    check=dict(rows=[dict(status='TRAPPED_SUPPORT_RISK',vascular_branch='O2'),dict(status='TRAPPED_SUPPORT_RISK',vascular_branch='B1')])
    feedback=e.route_feedback(geometry,check,'surface hash')
    assert feedback['allowed_faces_by_endpoint']=={'O2':['-X','+X','-Y']}
    assert feedback['central_geometry_changed'] is False and feedback['vascular_sha256']=='surface hash'
    assert e.route_feedback(geometry,dict(rows=[check['rows'][1]]),'hash') is None


def test_top_access_feedback_stops_when_legal_walls_exhausted():
    geometry=dict(ports=[dict(port_id='O2',face='-X')],route_review=dict(allowed_faces_by_endpoint={'O2':['-X']}))
    check=dict(rows=[dict(status='TRAPPED_SUPPORT_RISK',vascular_branch='O2')])
    assert e.route_feedback(geometry,check,'hash') is None


def test_port_wall_overlap_and_one_body_no_loss(cap_route,cfg):
    source,endpoint,bounds,route=cap_route;core,_=f.union_core(source,[route],cfg)
    box,_=p.new_box(source.bounds,cfg)
    inputs=dict(core=core,ports=[dict(port_id='I1',radius_mm=.6,face='+X',target=route.target,end=route.end,direction=route.direction)])
    rows,_=q.port_overlaps(inputs,box,cfg);assert rows[0]['status']=='PASS' and rows[0]['intersection_axial_span_mm']>=2.7
    union,_=mold.combine(core,box['mesh'],cfg);audit=q.union_audit(core,box['mesh'],union,cfg)
    assert audit['connected_components']==1 and audit['passed'] and audit['core_material_loss_mm3']<=.001


@pytest.mark.parametrize('case',FIXTURE['trapped_cases'])
def test_trapped_requires_both_failures(case):assert e.trapped(case['visible'],case['probe'])==case['trapped']


def test_top_only_visibility_detects_roof(cfg):
    box,_=p.new_box(np.array([[0,0,0],[10,10,10.]]),cfg)
    core=trimesh.creation.box([1,1,1]);core.apply_translation([5,5,5]);face=int(np.argmax(core.face_normals[:,2]))
    assert e.top_visibility(core,face,box,cfg)
    roof=trimesh.creation.box([100,100,1]);roof.apply_translation([5,5,15]);obstacle=trimesh.util.concatenate([core,roof])
    assert not e.top_visibility(obstacle,face,box,cfg)


def test_reduction_and_zero_denominator():
    assert m.reduction(10,7.5)==25 and m.reduction(0,1) is None
    assert m.reduction_gate(10,7.5,20,20) and m.reduction_gate(10,9,20,14)
    assert not m.reduction_gate(10,9,20,19)


def test_bambu_unavailable_explicit_pending(tmp_path,cfg):
    r=e.slice_cached(dict(status='BAMBU_STUDIO_NOT_FOUND'),cfg,{},tmp_path,e.make_logger(tmp_path),None,None,None)
    assert not r['printer_job_sent'] and not r['results'] and 'PENDING' in r['status']


@pytest.mark.parametrize('flag',['--send','--upload','--print','--start-print','--device-id','--host'])
def test_printer_actions_forbidden(flag):
    with pytest.raises(ValueError):check_local_command(['bambu-studio',flag,'something'])


def test_local_slice_command_allowed():assert check_local_command(['bambu-studio','--slice','0','--export-3mf','result.3mf','in.stl'])


def test_warning_xml_elements_are_not_dropped(tmp_path):
    path=tmp_path/'warning.3mf'
    with zipfile.ZipFile(path,'w') as z:z.writestr('Metadata/slice_info.config','<root><warning msg="floating object"/><warning msg="bed_temperature_too_high_than_filament"/></root>')
    warnings,critical=e.archived_warnings(path);assert len(warnings)==2 and len(critical)==1


def test_real_8_vascularmd_surfaces_and_ligament(real):
    rows=json.loads((real/'model_results.json').read_text());assert len(rows)==8
    for r in rows:
        assert r['VascularMD']['model_network_called'] and r['VascularMD']['topology_preserved']
        assert r['mesh']['watertight'] and r['mesh']['manifold'] and r['mesh']['degenerate_faces']==0
        assert r['pdms']['passed'] and r['radius_qc']['passed']


def test_real_candidates_mevo_scale_and_count(real):
    rows=json.loads((real/'candidates.json').read_text());assert len(rows)>=20
    assert all(r['mevo_only'] and r['scale_factor']==1 and r['outlet_count'] in (2,3) and r['topology_valid'] for r in rows)


def test_real_baseline_resliced_same_profile(real):
    path=real/'baseline/legacy_70mm_slice_metrics.json'
    if not path.exists():pytest.skip('Baseline slicing still pending')
    r=json.loads(path.read_text());assert r['real_resliced'] and s.sha256(r['source_stl'])==r['source_sha256']
    if r['selected_slice']:
        expected=json.loads((real/'bambu/evaluation_profile/sha256.json').read_text())
        assert r['selected_slice']['evaluation_profile_sha256']==expected
        assert r['selected_slice']['support_quantity']['actual_support_segment_count']>0


def test_all_completed_new_slices_use_identical_profiles(real):
    expected=json.loads((real/'bambu/evaluation_profile/sha256.json').read_text())
    checked=0
    for path in real.glob('candidates/*/orientations/slice_*/evaluation_slice.json'):
        result=json.loads(path.read_text())
        if result.get('selected'):
            assert result['selected']['evaluation_profile_sha256']==expected;checked+=1
    if not checked:pytest.skip('Real candidate slicing pending')


def test_delivered_real_port_interfaces_no_loss(real):
    records=[json.loads(p.read_text()) for p in real.glob('candidates/*/manufacturing/geometry.json')]
    if not records:pytest.skip('Real port/box generation pending')
    for record in records:
        assert record['union']['core_material_loss_mm3']<=.001
        assert record['box']['four_side_walls_closed'] and record['box']['bottom_closed'] and record['box']['top_fully_open']
        for port in record['final_interfaces']:
            assert port['feature_edge_20_length_mm']==0 and port['new_closed_sharp_ring']==0 and port['section_area_change_fraction']<.01


def test_old_s1_6_outputs_unchanged(real):
    snapshot=json.loads((real/'protected_before.json').read_text())
    assert any('s1-6_' in p for p in snapshot) and any('final_abs_casting_mold_2mm_X1C' in p and p.endswith('.stl') for p in snapshot)
    assert s.verify_snapshot(snapshot)['all_unchanged']


def test_real_adaptive_geometry_and_face_table(real):
    paths=list(real.glob('candidates/*/manufacturing/geometry.json'))
    assert paths
    for path in paths:
        record=json.loads(path.read_text());faces=record['box']['adaptive_faces']
        assert len(faces)==6
        for face in faces:
            if not face['port_connected'] and face['face']!='TOP_OPENING':
                assert face['absolute_clearance_pass'] and face['fcl_surface_distance_mm']>=1.6
                assert face['fcl_surface_distance_mm']==pytest.approx(face['core_envelope_gap_mm'],abs=1e-4)
        table=s.read_csv(path.parent/'adaptive_box_faces.csv')
        assert len(table)==6 and 'fcl_surface_distance_mm' in table[0]
        for port in record['ports']:
            assert port['outside_stub_mm']==6
            assert port['moderate_length_target_met']==(port['length_mm']<=port['preferred_maximum_extension_mm'])
            assert port['face'] in ('-X','+X','-Y','+Y')
        previous=json.loads((real/'previous_port_box_geometry.json').read_text()).get(path.parent.parent.name)
        if previous:
            assert max(p['length_mm'] for p in record['ports'])<max(p['length_mm'] for p in previous['ports'])


def test_real_print_archive_contains_intact_model_and_gcode(real):
    evaluated=[json.loads(p.read_text()) for p in real.glob('candidates/*/evaluation.json')]
    completed=[r for r in evaluated if r.get('orientations',{}).get('selected',{}).get('slice')]
    assert len(completed)==3
    for row in completed:
        chosen=row['orientations']['selected'];record=chosen['slice']
        assert record['status']=='BAMBU_SLICED' and record['feature_preservation']['passed']
        assert not record['critical_warnings'] and chosen['orientation']['trapped_risk_count']==0
        assert chosen['orientation']['passed']
        with zipfile.ZipFile(record['archive']) as archive:
            assert archive.testzip() is None
            assert any(n.endswith('.gcode') for n in archive.namelist())
            assert any(n.endswith('.model') for n in archive.namelist())


def test_real_trapped_o2_feedback_and_corrected_access(real):
    candidate=real/'candidates/candidate_0008';matched=False
    for path in candidate.glob('history/*/manufacturing/geometry.json'):
        geometry=json.loads(path.read_text())
        if not any(p['port_id']=='O2' and p['face']=='+Y' for p in geometry['ports']):continue
        orientations=json.loads((path.parent.parent/'orientations/orientation_candidates.json').read_text())
        feedback=e.route_feedback(geometry,orientations[0]['support'],s.sha256(candidate/'vascular_core.stl'))
        assert feedback and '+Y' not in feedback['allowed_faces_by_endpoint']['O2'];matched=True
    assert matched
    corrected=json.loads((candidate/'manufacturing/upright_port_access_qc.json').read_text())
    assert corrected['assessment']['trapped_risk_count']==0
