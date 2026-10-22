"""Permanent synthetic and real BG001 casting-mold regression tests."""
import copy
import json
import logging
from pathlib import Path
import zipfile
import networkx as nx
import numpy as np
import pytest
import trimesh
from vascular_processing import abs_casting_mold as mold
from vascular_processing import casting_mold_qc as qc
from vascular_processing import support_removal_qc as access
from vascular_processing import bambu_casting_mold as bambu
from vascular_processing import sacrificial_fixture as old
from vascular_processing.sacrificial_print_frame import box_mesh

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/final_abs_casting_mold'
HASH='b9fb2d8a57b58b877f940f6c2486f923d93ac71b42793593f4092180dcb7e35e'

@pytest.fixture(scope='module')
def cfg():return mold.load_config(ROOT/'config/abs_casting_mold_BG001.yaml')

@pytest.fixture(scope='module')
def fixture(cfg):
    spec=json.loads((ROOT/'tests/data/simple_five_wall_mold.json').read_text());nodes=np.array(spec['nodes'],float);r=spec['radius_mm'];graph=nx.DiGraph();branches={}
    def tube(a,b):return trimesh.creation.cylinder(radius=r,segment=[a,b],sections=32)
    parts=[]
    for a,b in spec['edges']:parts.append(tube(nodes[a],nodes[b]));graph.add_edge(a,b);branches[a,b]=f'{a}_{b}'
    for j,p in enumerate(nodes):
        graph.nodes[j]['coords']=np.r_[p,r];ball=trimesh.creation.icosphere(subdivisions=2,radius=r);ball.apply_translation(p);parts.append(ball)
    central=trimesh.boolean.union(parts,engine='manifold');full=[central];ports=[];inner=np.array(spec['inner_bounds_mm']);outer=np.array(spec['outer_bounds_mm'])
    for p in spec['ports']:
        axis,side=old.face_axis(p['face']);direction=old.NORMALS[p['face']];target=nodes[p['node']].copy();target[axis]=inner[side,axis];end=target+10*direction
        pipe=tube(nodes[p['node']]-.3*direction,end);full.append(pipe)
        ports.append(dict(port_id=p['port_id'],radius_mm=r,face=p['face'],target=target,end=end,direction=direction,
            reference_mesh=pipe,parent_branch=p['branch'],swc_id=p['node']))
    core=trimesh.boolean.union(full,engine='manifold');core,_=old.union_core(core,[],cfg)
    return dict(core=core,central=central,ports=ports,old_box=spec,companion=dict(graph=graph,branches=branches,transform=np.eye(4)),summary={})

@pytest.fixture(scope='module')
def box(fixture,cfg):return mold.build_intact_box(fixture,cfg)
@pytest.fixture(scope='module')
def boxqc(box,cfg):return qc.intact_box_qc(box,cfg)
@pytest.fixture(scope='module')
def combined(fixture,box,cfg):return mold.combine(fixture['core'],box['mesh'],cfg)[0]
@pytest.fixture(scope='module')
def audit(fixture,box,combined,cfg):return qc.union_audit(fixture['core'],box['mesh'],combined,cfg)
@pytest.fixture(scope='module')
def partition(fixture,box,combined,cfg):return qc.partition_vascular(combined,fixture,box,cfg)


def test_four_side_walls_exist(boxqc):assert boxqc['four_side_walls_closed']
def test_bottom_exists(boxqc):assert boxqc['bottom_closed']
def test_top_open(boxqc):assert boxqc['top_fully_open'] and boxqc['blocked_top_ray_count']==0
def test_no_old_holes(boxqc):assert boxqc['no_assembly_clearance_holes'] and all(r['missing_material_mm3']==0 for r in boxqc['walls'])
def test_dimensions_reused(box,fixture):np.testing.assert_array_equal(box['outer'],fixture['old_box']['outer_bounds_mm'])
def test_top_opening_dimensions(boxqc):assert boxqc['top_opening_area_mm2']==900 and boxqc['top_opening_width_mm']==30 and boxqc['top_opening_length_mm']==30

def test_closed_roof_fails(box,cfg):
    bad=copy.copy(box);roof=box_mesh([[-18,-18,23],[18,18,24]])
    bad['mesh']=trimesh.boolean.union([box['mesh'],roof],engine='manifold')
    assert not qc.intact_box_qc(bad,cfg)['top_fully_open']

def test_positive_port_intersections(fixture,box,cfg):
    rows,_=qc.port_overlaps(fixture,box,cfg)
    assert len(rows)==3 and all(r['intersection_volume_mm3']>0 for r in rows)

def test_port_overlap_depth(fixture,box,cfg):
    rows,_=qc.port_overlaps(fixture,box,cfg)
    assert all(r['intersection_axial_span_mm']>=.9*r['wall_thickness_mm'] and r['status']=='PASS' for r in rows)

def test_insufficient_overlap_rejected(fixture,box,cfg):
    bad=dict(fixture);bad['core']=qc.intersection(fixture['core'],box_mesh([[-100,-15.5,-100],[100,100,100]]))
    row=qc.port_overlaps(bad,box,cfg)[0][0]
    assert row['intersection_axial_span_mm']==pytest.approx(.5)
    assert row['status']=='PORT_WALL_OVERLAP_INSUFFICIENT_I1'

def test_final_one_component(audit):assert audit['connected_components']==1
def test_final_watertight(audit):assert audit['watertight']
def test_final_manifold(audit):assert audit['manifold']
def test_core_material_preserved(audit):assert audit['core_material_loss_mm3']<1e-6
def test_box_material_preserved(audit):assert audit['box_material_loss_mm3']<1e-6
def test_volume_balance(audit):assert audit['volume_balance_error_mm3']<1e-5

def test_stl_roundtrip(combined,cfg,tmp_path):assert old.export_checked_stl(combined,tmp_path/'mold.stl',cfg)['passed']

def test_ligament_excludes_true_junctions(fixture,combined,partition,cfg):
    result=qc.ligament_qc(combined,fixture,partition,cfg)
    assert result['minimum_ligament_mm']>0
    assert result['same_branch_pairs_excluded']
    assert any(r['exclusion']=='LOCAL_TRUE_BIFURCATION' for r in result['all_pairs'])
    assert all(r['first']!=r['second'] for r in result['all_pairs'])

def test_wall_clearance_excludes_legal_penetration(fixture,combined,partition,box,cfg):
    rows=qc.wall_clearance(combined,partition,fixture,box,cfg)
    exempt=[r for r in rows if r['legal_port_crossing_excluded']]
    assert len(exempt)==3 and all(r['distance_mm']==pytest.approx(2) for r in exempt)
    assert all(r['distance_mm']>0 for r in rows)

def test_support_threshold_from_active_profile(cfg):
    value=access.overhang_threshold(dict(flattened=dict(process=dict(support_threshold_angle='35'))),cfg)
    assert value['angle_from_horizontal_deg']==35 and value['source']=='ACTIVE_BAMBU_PROCESS_PROFILE'

def test_support_threshold_heuristic_named(cfg):assert access.overhang_threshold({},cfg)['source']=='OVERHANG_THRESHOLD_HEURISTIC'

def test_risk_identification_and_vertical_access(box,cfg):
    sphere=trimesh.creation.icosphere(subdivisions=2,radius=1);sphere.apply_translation([0,0,9])
    region=dict(region_id='vessel',vascular_branch='vessel',face_ids=np.arange(len(sphere.faces)),center_mm=[0,0,9])
    rotation=np.diag([1,-1,-1]) # Surface facing casting top now faces print down.
    result=access.assess_regions(sphere,box,[region],rotation,dict(angle_from_horizontal_deg=30),cfg)
    assert result['candidate_regions']==1 and result['internal_support_risk_area_mm2']>0
    assert result['vertical_visible_count']==1 and result['trapped_risk_count']==0

def test_probe_to_top(box,cfg):
    sphere=trimesh.creation.icosphere(subdivisions=2,radius=1);sphere.apply_translation([0,0,9])
    manager=trimesh.collision.CollisionManager();manager.add_object('vascular',sphere);manager.add_object('box',box['mesh'])
    path=access.top_path(np.array([0.,0.,10]),np.array([0.,0.,1]),3,box,sphere,manager,cfg)
    assert path is not None and path['start'][2]>box['inner'][1,2]

def test_trapped_closed_lid_support(box,cfg):
    bad=copy.copy(box);bad['mesh']=trimesh.boolean.union([box['mesh'],box_mesh([[-18,-18,23],[18,18,24]])],engine='manifold')
    sphere=trimesh.creation.icosphere(subdivisions=2,radius=1);sphere.apply_translation([0,0,9])
    region=dict(region_id='trapped',vascular_branch='test',face_ids=np.arange(len(sphere.faces)),center_mm=[0,0,9])
    result=access.assess_regions(sphere,bad,[region],np.diag([1,-1,-1]),dict(angle_from_horizontal_deg=30),cfg)
    assert result['trapped_risk_count']==1 and result['status']=='TRAPPED_SUPPORT_RISK'

def test_orientations_deterministic(cfg):
    a=access.rotations(cfg);b=access.rotations(cfg);assert len(a)==52
    for (n,r),(nn,rr) in zip(a,b):assert n==nn;np.testing.assert_array_equal(r,rr)

def test_build_volume_gate(combined,box,cfg):
    proxy=dict(threshold=dict(angle_from_horizontal_deg=30),trapped_risk_count=0,internal_support_risk_area_mm2=1,probe_accessible_count=1,top_accessible_count=1)
    row=access.orientation_row(0,'test',np.eye(3),combined,box,[20,20,20],proxy,cfg)
    assert not row['fits_build_volume'] and not row['passed']

def test_bambu_unavailable_fallback(monkeypatch,cfg,tmp_path):
    monkeypatch.setattr(bambu.adapter.existing,'discover',lambda _:dict(status='BAMBU_STUDIO_NOT_FOUND'))
    assert bambu.discover(cfg,tmp_path)['status']=='BAMBU_STUDIO_NOT_FOUND'

@pytest.mark.parametrize('option',['--send','--upload','--start-print','--print'])
def test_no_printer_send(option):
    with pytest.raises(ValueError):bambu.adapter.check_local_command(['bambu-studio',option])

def test_actual_gcode_arc_feature_parser(tmp_path):
    path=tmp_path/'sample.3mf'
    with zipfile.ZipFile(path,'w') as z:z.writestr('Metadata/plate_1.gcode','G90\nM83\nG1 X1 Y0 Z3\n; FEATURE: Outer wall\nG3 X0 Y1 I-1 J0 E1\n; FEATURE: Support\nG1 X0 Y2 E.5\n')
    data=bambu.parse_paths(path,allow_missing_settings=True)
    assert data['arc_command_count']==1 and data['kind'][-1]==1
    assert data['extrusion_filament_mm'].sum()==pytest.approx(1.5)
    np.testing.assert_allclose(data['segments'][-1],[[0,1,3],[0,2,3]])

def test_nonzero_archived_extruder_offset_and_custom_paths(tmp_path):
    path=tmp_path/'offset.3mf'
    with zipfile.ZipFile(path,'w') as z:
        z.writestr('Metadata/project_settings.config',json.dumps(dict(extruder_offset=['0x2'])))
        z.writestr('Metadata/plate_1.gcode','G90\nM83\nG1 X10 Y18 Z3\n; FEATURE: Outer wall\nG1 X11 Y18 E1\n; FEATURE: Custom\nG1 X12 E1\n; FEATURE: Support\nG1 X13 E1\n')
    data=bambu.parse_paths(path)
    np.testing.assert_array_equal(data['segments'][0],[[10,20,3],[11,20,3]])
    assert list(data['kind'])==[0,2,1]
    assert data['coordinate_frame']['offset_added_mm']==[0,2,0]

@pytest.mark.parametrize('settings',[{},dict(extruder_offset=['0x0','1x2'])])
def test_unresolved_gcode_coordinate_mapping_rejected(tmp_path,settings):
    path=tmp_path/'unknown.3mf'
    with zipfile.ZipFile(path,'w') as z:z.writestr('Metadata/project_settings.config',json.dumps(settings))
    with pytest.raises(ValueError):bambu.parse_paths(path)

def test_support_paths_cannot_satisfy_feature_presence(cfg):
    paths=dict(segments=np.array([[[0,0,0],[1,0,0]]]),kind=np.array([1]))
    result=bambu.feature_preservation(paths,None,None,None,np.eye(4),cfg)
    assert not result['passed'] and result['reason']=='No model extrusion paths'

def test_native_warning_extraction_does_not_confuse_feature_tags():
    record=dict(native_cli_result=dict(sliced_plates=[dict(feature_type_times={'Floating vertical shell':10},warning_message='')]))
    assert bambu.native_warning_lines(record)==[]
    record['native_cli_result']['sliced_plates'][0]['warning_message']='Object has floating regions. Enable support.'
    assert bambu.native_warning_lines(record)==['Object has floating regions. Enable support.']

@pytest.fixture(scope='module')
def real():
    path=OUT/'abs_casting_mold_summary.json';assert path.is_file(),'Run real s1-6 integration first'
    return json.loads(path.read_text())

def test_real_accepted_core_untouched(real):assert old.sha256(Path(real['input']['accepted_core_path']))==HASH

def test_real_all_port_outputs_untouched(real):
    before=json.loads((OUT/'protected_before.json').read_text());subset={k:v for k,v in before.items() if 'print_fixture_design_all_ports_aligned' in k}
    assert len(subset)>=40 and old.verify_snapshot(subset)['all_unchanged']

@pytest.mark.parametrize('name',['s1-2_swc_roi_generate_human.py','s1-3_swc_roi_generate_MeVO.py','s1-4_sacrificial_box_and_ports.py','s1-5_sacrificial_print_frame.py'])
def test_real_existing_entries_untouched(real,name):
    from vascular_processing.project_paths import frozen_hash
    before=json.loads((OUT/'protected_before.json').read_text());assert old.sha256(ROOT/name)==frozen_hash(before,ROOT/name)

def test_real_o3_gap_regression(real):
    assert real['pdms']['o3_gap_mm']>=1.6 and real['pdms']['o3_deviation_mm']<=.10

def test_real_union_geometry(real):
    q=real['union'];assert q['connected_components']==1 and q['watertight'] and q['manifold'] and q['degenerate_faces']==0
    assert q['core_material_loss_mm3']<=.001 and q['box_material_loss_mm3']<=.001
    assert all(r['status']=='PASS' for r in real['ports'])

def test_real_top_opening(real):assert real['box']['top_fully_open'] and real['box']['four_side_walls_closed'] and real['box']['bottom_closed']

def test_real_no_risk_hidden(real):
    if real['support']['trapped_risk_count']:assert real['status']=='NEEDS_ADJUSTMENT'
    assert not real['support']['real_support_removal_certified']

def test_real_bambu_and_features(real):
    records=real['bambu']['results'];assert records
    assert records[0]['orientation_id']==0 and records[0]['support_mode']=='OFF'
    for row in records:
        assert bambu.adapter.check_local_command(row['command'])
        if row['status']=='BAMBU_SLICED':assert row['geometry_qc']['rotation_scale_position_preserved']
    if real['bambu'].get('selected'):assert real['bambu']['selected']['feature_preservation']['passed']
    assert len({r['orientation_id'] for r in records})<=6 # mandatory baseline plus at most five ranked orientations

def test_real_fifteen_figures(real):
    from vascular_processing.casting_mold_review import NAMES
    for name in NAMES:assert (OUT/'QC'/(name+'.png')).stat().st_size>10000

def test_real_selected_slice_and_offsets(real):
    assert real['status']=='READY_FOR_HUMAN_REVIEW' and not real['failures']
    actual=real['bambu']['selected']
    assert actual['gcode_coordinate_frame']['offset_added_mm']==[0,2,0]
    assert len(actual['feature_preservation']['rows'])==14
    assert all(r['coverage_fraction']==1 for r in actual['feature_preservation']['rows'])
    assert Path(real['final_3mf']).is_file() and Path(real['final_gcode']).is_file()

def test_real_casting_restore_and_oriented_mesh(real):
    forward=np.array(real['orientation']['selected']['transform_4x4'])
    inverse=np.array(json.loads((OUT/'casting_restore_transform.json').read_text())['transform_4x4'])
    np.testing.assert_allclose(inverse@forward,np.eye(4),atol=1e-10)
    source=trimesh.load_mesh(real['final_stl']);placed=trimesh.load_mesh(real['final_oriented_stl'])
    source.apply_transform(forward)
    np.testing.assert_allclose(source.bounds,placed.bounds,atol=2e-5)
