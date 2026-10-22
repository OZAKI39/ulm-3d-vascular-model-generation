"""Permanent geometric, integration and source-protection support regression tests."""
from pathlib import Path
import copy,json,zipfile
import numpy as np
import pytest
import trimesh
import yaml
from shapely.geometry import box as polygon_box,GeometryCollection
from vascular_processing.support_inputs import ROOT,discover,digest,verify,profile_values,swc_branches,save
from vascular_processing import support_detection as det
from vascular_processing import support_freefloating as ff
from vascular_processing import support_contact_planner as cp
from vascular_processing import support_tree_planner as tree
from vascular_processing import support_scaffold as scaffold
from vascular_processing import support_slicer_adapter as adapter
from vascular_processing.support_accessibility import component_access
from vascular_processing.support_qc import candidate_gates,compare,nozzle_gate
from vascular_processing.vascular_support_pipeline import cache_key,parameter_grid

DATA=Path(__file__).parent/'data/support'
CFG=yaml.safe_load((ROOT/'config/vascular_support_planner.yaml').read_text())
PROFILE=dict(layer_height=.1,first_layer_height=.12,line_width=.22,nozzle_diameter=.2,support_top_gap=.2,support_xy_gap=.45)

def tube(name='horizontal_tube'):return trimesh.load_mesh(DATA/(name+'.stl'))
def raw_points(length=20.,step=.25):
    return [dict(point_id=f'b:{i}',branch_id='b',sample_index=i,xyz=[float(x),0,9.4],arc_position=float(x),local_radius=.6,local_curvature=0.,
        detectors_triggered=['SWC_LONG_SPAN'],unsupported_span=length,risk_level='CRITICAL',distance_to_bifurcation=10.,surface_hit_valid=True,is_port=False) for i,x in enumerate(np.arange(0,length+step/2,step))]

def test_source_discovery_latest_2mm():
    s=discover();assert '2mm' in str(s['mold']);assert s['smoothed_swc'].is_file();assert s['branch_mapping'].is_file();assert s['core'].is_file()

def test_source_discovery_missing(tmp_path):
    with pytest.raises(FileNotFoundError,match='SOURCE_NOT_FOUND'):discover(tmp_path)

def test_sha_protection_detects_change(tmp_path):
    p=tmp_path/'source.stl';p.write_text('before');hashes={str(p):digest(p)};assert verify(hashes)['all_unchanged'];p.write_text('after')
    with pytest.raises(RuntimeError,match='UPSTREAM_GEOMETRY_MODIFIED'):verify(hashes)

def test_profile_parsing_actual_units():
    v=profile_values(dict(machine={'nozzle_diameter':['.2']},process={'layer_height':'.1','initial_layer_print_height':'.12','line_width':'110%','support_top_z_distance':'.2','support_object_xy_distance':'.45'}))
    assert v['line_width']==pytest.approx(.22);assert v['first_layer_height']==.12;assert v['support_xy_gap']==.45;assert not v['fallback_fields']

def test_profile_fallback_explicit():
    p=profile_values({});assert p['status']=='PROFILE_VALUE_FALLBACK';assert len(p['fallback_fields'])==6

def test_triangle_overhang_area_normals():
    m=trimesh.creation.box([2,3,4]);r=det.overhang(m);assert r['area'][r['required']].sum()==pytest.approx(6.);assert len(r['centroids'])==len(m.faces)

def test_layer_unsupported_overlap():
    p=polygon_box(0,0,1,1);bad,reach=det.unsupported_layer(p,p,.1,.22);assert bad.area==0;assert reach==pytest.approx(.1)
    bad,_=det.unsupported_layer(p,GeometryCollection(),.1,.22);assert bad.area==1

def test_real_triangle_slicing_horizontal_tube():
    r=det.layers(tube(),PROFILE);assert r['rows'];assert max(p['unsupported_area_mm2'] for p in r['rows'])>.1

def test_layer_field_covers_long_polygon_ends():
    # A single polygon centroid cannot stand in for a 20 mm unsupported area.
    layer=dict(levels=np.array([9.46]),unsupported_polygons=[polygon_box(0,-.2,20,.2)])
    points=np.array([[.1,0,9.4],[10,0,9.4],[19.9,0,9.4]])
    distance=det.layer_field_distance(points,np.full(3,.6),layer)
    assert np.allclose(distance,.06)

def test_freefloating_wrapper_real_rays():
    points,audit=ff.detect(tube(),.1,.7,1.);assert len(points)>0;assert audit['nodes']>0;assert audit['implementation']=='MIT_CPU_PORT_RAY_INTERVAL_GRAPH';assert np.min(points[:,2])<9.6

def test_freefloating_quantization():assert np.array_equal(ff.machine_layer([.04,.06,.14,.16],.1),[0,1,1,2])

def test_freefloating_graph_reference_reach():
    xyz=np.array([[0,0,1],[1,0,1],[2,0,1],[3,0,1]]);edges=[[(1,1.)],[(0,1.),(2,1.)],[(1,1.),(3,1.)],[(2,1.)]]
    points,ids=ff.select_graph_points(xyz,edges,1.);assert ids==[0,2];assert len(points)==2

def test_local_minimum_detection():
    m=trimesh.creation.icosphere(subdivisions=1,radius=1);m.apply_translation([0,0,10]);p=ff.local_minima(m);assert len(p);assert p[:,2].min()==pytest.approx(9.)

def write_swc(tmp_path):
    f=tmp_path/'bif.swc';f.write_text('1 3 -8 0 10 .6 -1\n2 3 0 0 10 .6 1\n3 3 10 -6 12 .6 2\n4 3 10 6 12 .6 2\n');return f

def test_swc_mapping_and_topology(tmp_path):
    p=write_swc(tmp_path);branches,t=swc_branches(p,np.eye(4));assert len(branches)==3;assert t['bifurcation_ids']==[2];assert t['inlets']==[1];assert t['outlets']==[3,4]

def test_swc_transform_keeps_radius(tmp_path):
    p=write_swc(tmp_path);m=np.eye(4);m[:3,3]=[100,0,0];b,t=swc_branches(p,m);assert b[0]['xyz'][0,0]==92;assert np.allclose(b[0]['radius'],.6)

def test_curvature_curved_branch(tmp_path):
    p=tmp_path/'curve.swc';p.write_text('1 3 0 0 10 .6 -1\n2 3 4 0 10 .6 1\n3 3 4 4 10 .6 2\n');b,_=swc_branches(p,np.eye(4));assert b[0]['curvature'].max()>.5

def test_long_span_detection():
    s=np.arange(0,20.1,.25);xyz=np.column_stack([s,np.zeros(len(s)),np.full(len(s),9.4)])
    branch=dict(branch_id='b',arc=s,tangent=np.tile([1,0,0],(len(s),1)),radius=np.full(len(s),.6));mask,angle,d,rows=det.long_spans(branch,xyz,xyz,PROFILE,CFG);assert mask.all();assert rows[0]['risk']=='CRITICAL';assert rows[0]['length_mm']==20

def test_box_support_suppresses_false_long_span():
    s=np.arange(0,10.1,.25);xyz=np.column_stack([s,np.zeros(len(s)),np.full(len(s),9.4)]);b=dict(branch_id='b',arc=s,tangent=np.tile([1,0,0],(len(s),1)),radius=np.full(len(s),.6))
    mask,_,_,rows=det.long_spans(b,xyz,xyz,PROFILE,CFG,np.ones(len(s),bool));assert not mask.any();assert rows==[]

def test_deterministic_pruning_and_over_support():
    raw=raw_points();a,c=cp.prune(raw,CFG);b,d=cp.prune(copy.deepcopy(raw),CFG);assert a==b;assert c==d;assert len(a)<len(raw)/4;assert c['percent']==100

def test_keepout_is_avoided_when_coverage_allows():
    raw=raw_points();
    for p in raw:p['distance_to_bifurcation']=abs(p['arc_position']-10)
    contacts,_=cp.prune(raw,CFG);assert all(p['distance_to_bifurcation']>=2 for p in contacts)

def test_keepout_exception_is_not_silent():
    raw=raw_points(1.)
    for p in raw:p['distance_to_bifurcation']=.1
    contacts,c=cp.prune(raw,CFG);assert c['percent']==100;assert all(p['status']=='KEEP_OUT_EXCEPTION_CRITICAL_COVERAGE' for p in contacts);assert all(p['keepout_reason'] for p in contacts)

def test_coverage_is_along_same_branch():
    a=raw_points(4.);b=[{**p,'branch_id':'other','point_id':'other'+p['point_id']} for p in a];c,_=cp.prune(a,CFG);assert cp.coverage(a+b,c,CFG)['percent']==50

def test_unsupported_arc_and_missing_coverage():
    cov=cp.coverage(raw_points(),[],CFG);assert cov['percent']==0;assert cov['max_critical_unsupported_span']>=20

def test_thin_curved_contact_influence_shorter():
    p=raw_points()[0];a=cp.influence(p,CFG);q={**p,'local_radius':.3,'local_curvature':.8};assert cp.influence(q,CFG)<a

def test_taper_actual_mesh_geometry():
    m=tree.strut([0,0,0],[0,0,5],2.5,.5);assert m.is_volume;assert np.ptp(m.vertices[m.vertices[:,2]>4.99,:2],axis=0).max()==pytest.approx(.5)

def test_load_aware_radius():assert tree.load_diameter(2.,10)>tree.load_diameter(2.,1)>=2

def test_nozzle_gate_uses_actual_line_width():
    assert nozzle_gate(.5,PROFILE)
    with pytest.raises(ValueError,match='BELOW_RELIABLE_LINE_WIDTH'):nozzle_gate(.2,PROFILE)

def test_tree_merge_angle_and_floor():
    p=tree.merge_node([0,0,10],[4,0,10],40.,0.);assert p is not None;assert p[2]<10;assert tree.merge_node([0,0,2],[40,0,2],40,0) is None

def test_colocated_contacts_do_not_create_zero_struts():
    a=np.array([0,0,10.]);m=tree.merge_node(a,a,40.,0);assert np.linalg.norm(m-a)>.1;assert tree.strut(m,a,2.,.6).is_volume

def test_pillar_generation():
    m=tree.strut([0,0,0],[0,0,10],2.5);assert m.is_volume;assert m.volume==pytest.approx(np.pi*1.25**2*10,rel=.02)

def test_collision_detection_real_mesh():
    vessel=tube();support=tree.strut([10,0,0],[10,0,10],.6);assert not tree.illegal_collision(support,vessel)['passed']

def test_support_top_gap():
    vessel=tube();support=tree.strut([10,0,0],[10,0,9.2],.6);assert tree.illegal_collision(support,vessel)['passed'];assert trimesh.proximity.closest_point(vessel,[[10,0,9.2]])[1][0]==pytest.approx(.2,abs=.01)

def test_support_box_collision():
    vessel=tube();box=trimesh.creation.box([2,2,2]);support=tree.strut([0,0,-4],[0,0,4],.6);assert not tree.illegal_collision(support,vessel,box)['passed']

def test_bridge_limit_and_generation():
    mesh,r=scaffold.bridge([0,0,5],[4,0,5],1,6);assert mesh.is_volume;assert r['status']=='EXPERIMENTAL_BRIDGE'
    with pytest.raises(ValueError,match='SPAN_EXCEEDS'):scaffold.bridge([0,0,5],[9,0,5],1,6)

def test_bridge_self_printability_not_assumed():
    mesh,_=scaffold.bridge([0,0,5],[6,0,5],1,6);r=scaffold.self_printability(mesh,PROFILE,6);assert not r['passed']

def test_scaffold_coincident_attachment_does_not_create_zero_bridge():
    vessel=tube();wall=trimesh.creation.box([1,1,1]);wall.apply_translation([50,50,0])
    contacts=[dict(contact_id=f'C{i}',xyz=[10,0,9.4],tip_diameter=.6) for i in range(2)]
    mesh,record=scaffold.build(contacts,vessel,wall,np.array([[-30,-30,0],[30,30,20]]),PROFILE)
    assert mesh.is_volume
    assert len(record['pillars'])==2
    assert record['bridges']==[]

def test_accessible_support():
    v=tube();wall=trimesh.creation.box([1,1,1]);wall.apply_translation([50,50,0]);s=tree.strut([10,0,0],[10,0,8],2);r=component_access([s],v,wall,np.array([[-10,-10,0],[30,30,20]]));assert r['trapped_support_count']==0

def test_trapped_support_closed_cavity():
    shell=tube('trapped_support_case');v=tree.strut([0,0,3],[0,0,10],1);dummy=trimesh.creation.box([1,1,1]);dummy.apply_translation([50,50,50]);r=component_access([v],dummy,shell,shell.bounds);assert r['trapped_support_count']==1

def test_enforcer_and_blocker_disjoint():
    m=tube();raw=raw_points(4.);contacts,_=cp.prune(raw,CFG);enforcer,blocker=cp.masks(m,contacts);assert enforcer.is_volume;assert blocker.is_volume
    common=trimesh.boolean.intersection([enforcer,blocker],engine='manifold');assert abs(common.volume)<1e-5

def test_project_roles_and_geometry_roundtrip(tmp_path):
    m=tube();m.export(tmp_path/'source.stl');contacts,_=cp.prune(raw_points(4.),CFG);e,b=cp.masks(m,contacts);p=tmp_path/'project.3mf';roles=adapter.write_project(p,m,e,b,dict(machine={'name':'X1C'},process={'name':'Vessel'},filament={'name':'ABS'}));assert roles==['normal_part','support_enforcer','support_blocker'];assert adapter.production_geometry(p,tmp_path/'source.stl')['passed']

def test_engine_unavailable_fallback(tmp_path):
    r=adapter.slice_local({'status':'ENGINE_UNAVAILABLE','engine':'Orca'},{},tmp_path/'missing.stl',tmp_path);assert r['status']=='ENGINE_UNAVAILABLE';assert not r['slicer_success']

def test_actual_cli_parser():assert adapter.cli_options('Use --help and --slice 0, --export-3mf file')==['--export-3mf','--help','--slice']

def test_slicer_warning_parser():
    r=adapter.warnings({'native_result':{'sliced_plates':[{'warning_message':'unsupported floating region'}],'error_string':'Success.'}});assert len(r['critical'])==1;assert 'Success.' not in r['all']

@pytest.mark.parametrize('option',['--upload','--send','--start-print','--print','--host','--device'])
def test_no_printer_upload_or_start(option):
    with pytest.raises(ValueError,match='NONLOCAL'):adapter.local_command(['engine',option,'x'])

def test_only_local_slice_allowed():assert adapter.local_command(['engine','--slice','0','--export-3mf','file'])

def test_hard_gates_prioritize_coverage():
    r=dict(source_unchanged=True,slicer_success=True,coverage={'percent':99,'max_critical_unsupported_span':6},trapped_support_count=0,no_illegal_collision=True);f=candidate_gates(r,CFG);assert 'CRITICAL_COVERAGE_INCOMPLETE' in f;assert 'CRITICAL_UNSUPPORTED_SPAN_EXCEEDED' in f

def test_comparison_never_prefers_failed_low_area():
    good=dict(contact_area_mm2=10,tip_count=3,coverage={'percent':100},support_near_bifurcation_count=0,trapped_support_count=0,top_accessible_fraction=1.,estimated_support_mass_g=1.,estimated_print_time_seconds=100,hard_gate_failures=[])
    bad={**good,'contact_area_mm2':1,'hard_gate_failures':['CRITICAL_COVERAGE_INCOMPLETE']};rows,chosen=compare({'BAMBU_NATIVE_BASELINE':good,'GUIDED_ORCA':bad});assert chosen=='BAMBU_NATIVE_BASELINE';assert rows[1]['contact_area_reduction_percent']==90

def test_cache_key_changes_with_layer_and_orientation():
    s={'hashes':{'a':'b'},'print_transform':np.eye(4).tolist()};a=cache_key(s,PROFILE,CFG);assert a!=cache_key(s,{**PROFILE,'layer_height':.12},CFG);t=copy.deepcopy(s);t['print_transform'][0][3]=1;assert a!=cache_key(t,PROFILE,CFG)

def test_grid_is_bounded_and_deterministic():
    a=parameter_grid(CFG);assert len(a)==12;assert a==parameter_grid(CFG);assert {float(p['tree_support_tip_diameter']) for p in a}=={.5,.6,.8}

@pytest.mark.parametrize('name',['horizontal_tube','curved_tube','bifurcation','long_bridge_tube','two_close_vessels','trapped_support_case','over_supported_case'])
def test_persistent_fixture_valid(name):assert tube(name).is_volume;assert (DATA/(name+'.json')).exists()

def test_upstream_s1_2_through_s1_7_and_models_unchanged():
    s=discover();p=s['root'].parent/'vascular_support_planner/protected_before.json';h=json.loads(p.read_text());assert any('s1-7_' in v for v in h);assert any('s1-2_' in v for v in h);assert str(s['mold']) in h;assert str(s['core']) in h;assert verify(h)['all_unchanged']

# Explicit production tests extend the original benchmark tests.
from vascular_processing import support_breakaway as breakaway
from vascular_processing import support_custom_graph as custom_graph
from vascular_processing import support_custom_qc as custom_qc
from vascular_processing.support_custom_detection import intervals,span_sensitive_raw
from vascular_processing.support_custom_pipeline import profiles_support_off,cache_identity
from vascular_processing.support_removal import cut_necks,removal_paths,swept_path

def test_tip_and_neck_widths_derive_from_actual_profile():
    assert breakaway.widths(PROFILE)==[.44,.66,.88]
    assert breakaway.widths({**PROFILE,'line_width':.3})==[.6,.9,1.2]

def test_coupon_all_27_cells_survive_stl_roundtrip(tmp_path):
    mesh,rows,features=breakaway.coupon(tmp_path,PROFILE)
    saved=trimesh.load_mesh(tmp_path/'support_breakaway_coupon.stl')
    assert mesh.is_volume and saved.is_volume
    assert len(rows)==27 and len(features)==54
    assert len({(r['neck_width_mm'],r['tip_width_mm'],r['top_gap_mm']) for r in rows})==27
    assert {r['top_gap_mm'] for r in rows}=={.15,.20,.25}

def test_interval_fusion_keeps_branch_and_evidence():
    raw=raw_points(4.);raw[2]['detectors_triggered']+=['FREEFLOATING','LAYER_UNSUPPORTED']
    r=intervals(raw,CFG);assert len(r)==1;assert r[0]['length_mm']==4.;assert 'FREEFLOATING' in r[0]['detectors_triggered']
    other=[{**p,'branch_id':'other','point_id':'other'+p['point_id']} for p in raw]
    assert len(intervals(raw+other,CFG))==2

def test_span_sensitivity_preserves_mandatory_evidence_without_promoting_proximity():
    raw=raw_points(4.)
    for p in raw:p['unsupported_span']=4.;p['detectors_triggered']=['FREEFLOATING'];p['freefloating_mandatory']=False
    raw[0]['freefloating_mandatory']=True
    raw[1]['detectors_triggered'].append('LOCAL_MINIMUM')
    changed=span_sensitive_raw(raw,5.)
    assert [i for i,p in enumerate(changed) if p['risk_level']=='CRITICAL']==[0,1]
    assert all(p['risk_level']=='CRITICAL' for p in span_sensitive_raw(raw,3.))

def test_local_minimum_is_not_deleted_by_arc_cover():
    raw=raw_points();raw[40]['detectors_triggered']+=['LOCAL_MINIMUM']
    selected,cov=cp.prune(raw,{**CFG,'selection_priority':'vascular'})
    assert cov['percent']==100;assert any(abs(p['arc_position']-10)<.51 for p in selected)

def test_custom_selection_deterministic_and_no_redundant_tips():
    raw=raw_points();cfg={**CFG,'selection_priority':'vascular'}
    a,v=cp.prune(raw,cfg);b,w=cp.prune(copy.deepcopy(raw),cfg);assert a==b;assert v==w
    assert all(cp.coverage(raw,[p for p in a if p!=q],cfg)['percent']<100 for q in a)

def test_cone_merge_lies_in_both_printable_cones():
    a=np.array([0.,0,12]);b=np.array([4.,0,10]);m=custom_graph.cone_intersection(a,b,35)
    for p in [a,b]:assert np.linalg.norm(p[:2]-m[:2])/(p[2]-m[2])<=np.tan(np.radians(35))+1e-8
    assert np.allclose(m,custom_graph.cone_intersection(b,a,35))

def test_unstable_edge_detected_and_widening_repairs_without_tips():
    a=[0,0,0];b=[0,0,30]
    assert custom_graph.stability(a,b,.5,1,CFG)['status']=='UNSTABLE_SUPPORT_EDGE'
    assert custom_graph.stability(a,b,3.,1,CFG)['status']=='GEOMETRIC_STABILITY_SCREEN_PASSED'

@pytest.fixture
def explicit_pair():
    # Leave enough height for the printable taper, shared cone, and breakaway
    # root; this fixture tests a feasible merge, not a short pillar fallback.
    vascular=trimesh.creation.box([20,2,1.2]);vascular.apply_translation([10,0,14])
    bottom=trimesh.creation.box([40,30,2]);bottom.apply_translation([10,0,-1]);inner=np.array([[-10,-15,0],[30,15,25]])
    contacts=[]
    for i,x in enumerate([8.,12.]):
        p=raw_points()[0].copy();p.update(contact_id=f'C{i}',point_id=f'b:{i}',xyz=[x,0,13.4],arc_position=x,distance_to_bifurcation=10)
        contacts.append(p)
    support,g,feet=custom_graph.build(contacts,vascular,bottom,inner,PROFILE,CFG,.44,.66,2.2)
    return support,g,feet,vascular,bottom,inner

def test_clever_graph_accumulates_load_and_roots(explicit_pair):
    support,g,feet,_,_,_=explicit_pair;assert support.is_volume;assert len(g['merge_nodes'])>=1
    assert all(t['root_connected'] for t in g['tips']);assert max(n['load'] for n in g['nodes'])==2
    assert g['roots'][0]['trunk_load']==2

def test_exact_triangle_top_gap(explicit_pair):
    _,g,_,vascular,_,_=explicit_pair
    assert all(t['measured_top_gap_mm']==pytest.approx(.2,abs=.003) for t in g['tips'])

def test_root_can_move_to_avoid_lower_vessel_without_new_contact():
    upper=trimesh.creation.box([20,10,1.2]);upper.apply_translation([0,0,20])
    lower=tree.strut([-3,0,10],[3,0,10],1.2);vascular=trimesh.util.concatenate([upper,lower])
    bottom=trimesh.creation.box([30,30,2]);bottom.apply_translation([0,0,-1]);inner=np.array([[-15,-15,0],[15,15,25]])
    p=raw_points()[0];p.update(contact_id='C0',xyz=[0,0,19.4])
    support,g,feet=custom_graph.build([p],vascular,bottom,inner,PROFILE,CFG,.44,.66,2.2)
    assert len([t for t in g['tips'] if t['root_connected']])==1
    assert any(r['action']=='OFFSET_PILLAR_COLLISION_REPAIR' for r in g['stability_repairs'])
    assert all(r['added_vascular_tips']==0 for r in g['stability_repairs'])
    assert custom_qc.self_printability(support,bottom,PROFILE)['passed']

def test_offset_pillar_caps_remain_horizontal():
    m=custom_graph.pillar_strut([0,0,2],[0,4,20],2.36,1.2)
    assert m.is_volume
    assert m.bounds[:,2]==pytest.approx([2,20])
    cap=m.vertices[np.isclose(m.vertices[:,2],2)]
    assert np.ptp(cap[:,:2],axis=0)==pytest.approx([2.36,2.36])

def test_body_clearance_uses_fcl_not_parameters(explicit_pair):
    _,g,_,vascular,_,_=explicit_pair
    for e in g['edges']:
        if e['kind']=='taper':continue
        m=custom_graph.edge_geometry(g,e)
        assert custom_graph.minimum_distance(m,vascular)>=.45-.015

def test_breakaway_foot_neck_and_width_are_real_geometry():
    m,foot,r=breakaway.root_geometry([0,0],2,2.2,.66,PROFILE)
    assert m.is_volume;assert r['neck_lines']==pytest.approx(3);assert r['contact_area_mm2']==pytest.approx(np.pi*r['foot_diameter_mm']**2/4)
    cut=m.section([0,0,1],r['estimated_break_location']);assert cut is not None
    assert np.ptp(cut.vertices[:,:2],axis=0).max()==pytest.approx(.66,abs=.005)

def test_root_only_allowed_on_inner_bottom():
    _,_,r=breakaway.root_geometry([0,0],2,2.2,.66,PROFILE);bounds=np.array([[-10,-10,2],[10,10,20]])
    assert breakaway.bottom_root_allowed(r,bounds)
    assert not breakaway.bottom_root_allowed({**r,'xyz':[0,0,5]},bounds)
    assert not breakaway.bottom_root_allowed({**r,'xyz':[10,0,2]},bounds)

def test_only_approved_root_union(explicit_pair):
    s,g,feet,v,b,inner=explicit_pair;r=custom_qc.root_union_audit(s,v,b,feet,g['roots'],inner)
    assert r['passed'];assert r['approved_bottom_overlap_volume_mm3']>0;assert r['vascular_and_ports_overlap_volume_mm3']<1e-6

def test_illegal_vessel_union_fails(explicit_pair):
    s,g,feet,v,b,inner=explicit_pair;bad=trimesh.boolean.union([s,tree.strut([8,0,0],[8,0,v.centroid[2]],1)],engine='manifold')
    assert not custom_qc.root_union_audit(bad,v,b,feet,g['roots'],inner)['passed']

def test_stl_roundtrip_removes_only_degenerate_faces(explicit_pair):
    import io
    s,*_=explicit_pair
    dirty=trimesh.Trimesh(vertices=s.vertices,faces=np.vstack([s.faces,[0,0,1]]),process=False)
    clean=custom_qc.canonical_stl(dirty)
    loaded=trimesh.load_mesh(io.BytesIO(clean.export(file_type='stl')),file_type='stl')
    assert loaded.is_volume
    assert loaded.volume==pytest.approx(s.volume,rel=1e-5)
    assert loaded.bounds==pytest.approx(s.bounds,abs=1e-5)

def test_illegal_bottom_area_outside_feet_fails(explicit_pair):
    s,g,feet,v,b,inner=explicit_pair;extra=tree.strut([23,0,-.2],[23,0,4],2);bad=trimesh.util.concatenate([s,extra])
    assert not custom_qc.root_union_audit(bad,v,b,feet,g['roots'],inner)['passed']

def test_neck_cut_removes_roots_but_retains_tree(explicit_pair):
    s,g,_,_,_,_=explicit_pair;parts=cut_necks(s,g['roots']);assert parts
    z=max(r['estimated_break_location'][2] for r in g['roots']);assert all(p.bounds[0,2]>=z-1e-5 for p in parts)
    assert sum(p.volume for p in parts)<s.volume

def test_swept_collision_catches_thin_obstacle_between_safe_endpoints():
    piece=trimesh.creation.box([.2,.2,.2]);wall=trimesh.creation.box([.05,2,2]);a=np.eye(4);b=np.eye(4);a[0,3]=-2;b[0,3]=2
    assert not swept_path(piece,[wall],[a,b])['passed']
    a[1,3]=b[1,3]=3;assert swept_path(piece,[wall],[a,b])['passed']

def test_tilt_can_release_support_blocked_on_pure_vertical():
    support=tree.strut([10,0,2.5],[10,0,9.],.44);vascular=tube();local,_,paths=removal_paths(support,15)
    assert not swept_path(local,[vascular],paths[0][1])['passed']
    assert any(swept_path(local,[vascular],poses)['passed'] for name,poses in paths if name.startswith('TILT'))

def test_closed_cavity_traps_whole_component():
    shell=tube('trapped_support_case');support=tree.strut([0,0,3],[0,0,10],.5);local,_,paths=removal_paths(support,shell.bounds[1,2])
    assert all(not swept_path(local,[shell],poses)['passed'] for _,poses in paths)

def test_auto_support_off_profile_preserves_nozzle():
    d={'flattened':{'machine':{'name':'X1C','nozzle_diameter':['.2']},'process':{'name':'Vessel','enable_support':'1'},'filament':{'name':'ABS'}}}
    p=profiles_support_off(d);assert p['process']['enable_support']=='0';assert p['machine']['nozzle_diameter']==['.2'];assert d['flattened']['process']['enable_support']=='1'

def test_feature_preservation_detects_dropped_tip():
    features=[dict(feature_id='tip',kind='tip',xyz=[0,0,5],width_mm=.44)]
    assert custom_qc.feature_preservation(np.array([[0,0,5.]]),features,PROFILE)['passed']
    assert custom_qc.feature_preservation(np.array([[2,0,5.]]),features,PROFILE)['status']=='SUPPORT_TIP_DROPPED'

def test_native_success_cannot_promote_custom_failure():
    r=dict(source_unchanged=True,planned_coverage={'percent':90,'max_critical_unsupported_span':6},self_printability={'passed':True},union_audit={'passed':True},gap_audit={'passed':True},removal={'trapped_support_count':0},slice={'slicer_success':True,'auto_support_off':False,'native_support_toolpath_count':1,'features':{'passed':True}},contact_area_reduction_percent=50)
    f=custom_qc.hard_gates(r,CFG);assert 'CRITICAL_INTERVAL_COVERAGE_INCOMPLETE' in f;assert 'NATIVE_AUTO_SUPPORT_NOT_OFF' in f

def test_cache_changes_with_gap_or_full_profile_hash():
    inputs={'hashes':{'mesh':'abc'},'print_transform':np.eye(4).tolist()};a=cache_identity(inputs,PROFILE,CFG)
    assert a!=cache_identity(inputs,{**PROFILE,'support_top_gap':.3},CFG)

@pytest.mark.parametrize('name',['long_curved_tube','S_curve','single_bifurcation','double_bifurcation','thin_vessel','close_parallel_vessels','unstable_tree','pillar_bridge','breakaway_root','trapped_tree'])
def test_custom_persistent_fixtures(name):
    assert tube(name).is_volume;assert json.loads((DATA/(name+'.json')).read_text())['synthetic']
