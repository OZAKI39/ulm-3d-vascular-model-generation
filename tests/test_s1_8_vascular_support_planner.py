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

# FDM_NATIVE_VASCULAR_SUPPORT_V2: V1 regressions above remain untouched.
from vascular_processing import support_v2_orientation as vo
from vascular_processing import support_v2_primitives as vp
from vascular_processing import support_v2_layers as vl
from vascular_processing import support_v2_planner as vg
from vascular_processing import support_v2_pipeline as vpipe
V2PROFILE={**PROFILE,'bridge_limit_mm':2.,'bridge_calibration_status':'BRIDGE_LIMIT_UNCALIBRATED'}


def test_v2_deterministic_full_3d_orientation_space():
    a=np.array(vo.rotations());assert a.shape==(128,3,3);assert np.array_equal(a,vo.rotations())
    assert len(np.unique(np.round(a[:,2],7),axis=0))==32
    for r in a:
        m=np.eye(4);m[:3,:3]=r;assert vo.rigid(m)


def test_v2_casting_restore_and_machine_volume():
    m=trimesh.creation.box([60,50,25]);bounds=np.array([[0,0,0],[256,256,250.]])
    for r in vo.rotations():
        tf=vo.bed_transform(m.vertices,r,bounds);inverse=np.linalg.inv(tf)
        q=trimesh.transform_points(m.vertices,tf)
        assert np.allclose(trimesh.transform_points(q,inverse),m.vertices)
        assert vo.fits(m.vertices,tf,bounds)
    assert not vo.fits(m.vertices,tf,[[0,0,0],[10,10,10]])
    scaled=tf.copy();scaled[:3,:3]*=1.001;assert not vo.rigid(scaled)


@pytest.mark.parametrize('kind',vp.KINDS)
def test_v2_primitive_self_buildability(kind):
    m=vp.primitive_gallery(V2PROFILE)[kind];assert m.is_volume
    m.apply_translation([0,0,-m.bounds[0,2]])
    assert vg.self_buildability(m,None,V2PROFILE)['passed']


def test_v2_short_bridge_requires_two_anchors():
    span=polygon_box(0,0,2,.8)
    both=polygon_box(-1,0,0,.8).union(polygon_box(2,0,3,.8))
    assert vl.bridge_evidence(span,both,2.,.22)['bridgeable']
    assert not vl.bridge_evidence(span,polygon_box(-1,0,0,.8),2.,.22)['bridgeable']
    assert not vl.bridge_evidence(span,both,1.,.22)['bridgeable']


def test_v2_full_polygon_coverage_not_centroid_or_arc():
    poly=polygon_box(0,-.2,20,.2)
    field=dict(rows=[dict(region_id='r',classification='SUPPORT_REQUIRED',z_mm=10.)],polygons=[poly])
    roof=dict(roof_id='tiny',top_z_mm=9.8,gap_mm=.2,_footprint=polygon_box(9.5,-.2,10.5,.2))
    c=vl.coverage(field,[roof],V2PROFILE);assert not c['all_required_covered'];assert c['percent']<10
    roof['_footprint']=poly
    assert vl.coverage(field,[roof],V2PROFILE)['all_required_covered']
    assert not vl.coverage(field,[roof],V2PROFILE,actual_ids=set())['all_required_covered']


def test_v2_quantization_and_layer_prism_rejection():
    assert vp.quantize(.71,.22)/.22==pytest.approx(4)
    assert vg.prismatic_layer_check(vp.branch([0,0,0],[1,0,8],1.32),V2PROFILE)
    assert not vg.prismatic_layer_check(vp.branch([0,0,0],[20,0,2],1.32),V2PROFILE)


def test_v2_redundant_roof_and_bifurcation_soft_penalty():
    poly=polygon_box(0,0,1,1);field=dict(rows=[dict(region_id='r',classification='SUPPORT_REQUIRED',z_mm=1.)],polygons=[poly])
    base=dict(valid=True,_footprint=poly,top_z_mm=.8,gap_mm=.2,distance_to_bifurcation=0.,curvature=0.)
    selected,_=vg.select_roofs([{**base,'roof_id':'a'},{**base,'roof_id':'b'}],[field],V2PROFILE)
    assert len(selected)==1;assert selected[0]['distance_to_bifurcation']==0.


@pytest.mark.parametrize('domain',['BUILD_PLATE','BOX_BOTTOM','BOX_SIDE_WALL','SACRIFICIAL_BOX_SURFACE'])
def test_v2_legal_root_domains_have_breakaway_necks(domain):
    mesh,foot,r=vp.breakaway_foot([0,0],0,1.76,V2PROFILE,domain,'root')
    assert mesh.is_volume and foot.is_volume
    assert r['neck_width_mm']==pytest.approx(3*.22)
    assert r['neck_height_mm']==pytest.approx(3*.1)
    assert r['approved_box_contact']==(domain!='BUILD_PLATE')


@pytest.mark.parametrize('name',['straight_horizontal_tube','long_curved_tube','S_tube','Y_bifurcation','high_arch','low_arch','external_port','tube_over_wall','rail_supported_tube','bridge_supported_tube','trapped_scaffold'])
def test_v2_permanent_fixture_meshes(name):
    m=trimesh.load_mesh(DATA/'v2'/(name+'.stl'));assert m.is_volume
    assert json.loads((DATA/'v2'/(name+'.json')).read_text())['synthetic']


def test_v2_trapped_fixture_fails_swept_opening_path():
    from vascular_processing.support_removal import swept_path
    closed=trimesh.load_mesh(DATA/'v2/trapped_scaffold.stl')
    part=trimesh.creation.box([1,1,2]);initial=np.eye(4);end=np.eye(4);end[2,3]=12.
    assert not swept_path(part,[closed],[initial,end])['passed']


def test_v2_actual_coverage_cannot_be_replaced_by_planned():
    c=dict(vascular_coverage={'all_required_covered':True},box_coverage={'all_required_covered':True},
        self_buildability={'passed':True},union_audit={'passed':True},build_volume_fit=True,removal={'trapped_pieces':0},
        slice=dict(slicer_success=True,auto_support_off=True,all_features_preserved=True,all_roots_have_start_paths=True,vascular_toolpaths_preserved=True,
                   actual_vascular_coverage={'all_required_covered':False},actual_box_coverage={'all_required_covered':True},native_support_toolpaths=0,support_only_slicer_success=True,bridge_toolpaths={'passed':True}))
    assert 'ACTUAL_ROOF_COVERAGE_INCOMPLETE' in vpipe.gates(c)


def test_v2_fdm_branch_feature_samples_follow_slant():
    from vascular_processing.support_v2_validation import graph_features
    p=dict(primitive_id='p',kind='SELF_SUPPORTING_BRANCH',bbox_mm=[[-1,-1,0],[5,1,10]],start=[0,0,0],end=[4,0,10],width_mm=2.)
    rows=graph_features({'primitives':[p]},V2PROFILE)
    assert np.allclose([r['xyz'] for r in rows],[[1,0,2.5],[2,0,5],[3,0,7.5]])


def test_v2_v1_cannot_enter_production():
    from vascular_processing.support_custom_pipeline import run
    with pytest.raises(RuntimeError,match='LEGACY_EXPLICIT_SUPPORT_V1_RETIRED'):run(stage='geometry')


def test_v2_bridge_gcode_filter_preserves_absolute_position_and_offsets(tmp_path):
    from vascular_processing.support_v2_validation import bridge_paths
    archive=tmp_path/'bridges.3mf'
    with zipfile.ZipFile(archive,'w') as z:
        z.writestr('Metadata/project_settings.config',json.dumps({'extruder_offset':['-1x0']}))
        z.writestr('Metadata/plate_1.gcode','M83\nG90\n; FEATURE: Outer wall\nG1 X2 Y3 Z1 E1\n; FEATURE: Bridge\nG1 X5 Y3 E1\n; FEATURE: Inner wall\nG1 X8 Y3 E1\n')
    segments=bridge_paths(archive);assert segments.shape==(1,2,3)
    assert np.allclose(segments[0],[[1,3,1],[4,3,1]])


def test_v2_removal_restores_casting_and_removes_upper_obstacle_first(tmp_path):
    from vascular_processing.support_v2_removal import run_removal
    from vascular_processing.support_removal import swept_path
    from scipy.spatial.transform import Rotation
    low=vp.pillar([0,0,1.2],[0,0,3.],1.0);high=vp.pillar([0,0,3.2],[0,0,5.],1.0)
    initial=np.eye(4);end=np.eye(4);end[2,3]=10
    assert not swept_path(low,[high],[initial,end])['passed']
    base=vp.loft([0,0,0],[0,0,1],[18,18]);walls=[]
    for x,y,size in [(0,-8.5,[18,1]),(0,8.5,[18,1]),(-8.5,0,[1,18]),(8.5,0,[1,18])]:walls.append(vp.loft([x,y,1],[x,y,8],size))
    box=trimesh.boolean.union([base,*walls],engine='manifold');box.export(tmp_path/'box.stl')
    vessel=trimesh.creation.box([1,1,1]);vessel.apply_translation([25,0,3]);vessel.export(tmp_path/'core.stl')
    transform=np.eye(4);transform[:3,:3]=Rotation.from_euler('x',90,degrees=True).as_matrix();transform[:3,3]=[20,20,20]
    support=trimesh.util.concatenate([low,high]);support.apply_transform(transform);support.export(tmp_path/'custom_support.stl')
    (tmp_path/'support_graph_v2.json').write_text(json.dumps(dict(segment_necks=[],roots=[])));(tmp_path/'roofs.json').write_text('[]')
    r=run_removal(dict(candidate_id='SYNTHETIC',folder=str(tmp_path),orientation={'casting_restore_transform':np.linalg.inv(transform).tolist()}),dict(core=tmp_path/'core.stl',box=tmp_path/'box.stl',inner_bounds=[[-8,-8,1],[8,8,8]]))
    assert r['trapped_pieces']==0 and len(r['sequence'])==2 and r['other_fragments_are_obstacles']
    rows={s['segment_id']:s for s in r['segments']};assert rows[r['sequence'][0]]['bounds_mm'][0][2]>rows[r['sequence'][1]]['bounds_mm'][0][2]
    assert all(rows[i]['path']['coordinate_frame']=='FROZEN_CASTING' for i in r['sequence'])


def test_v2_two_roofs_can_bridge_but_actual_requires_bridge_extrusions():
    field=dict(rows=[dict(region_id='r',classification='SUPPORT_REQUIRED',z_mm=1.)],polygons=[polygon_box(0,0,4,.8)])
    roofs=[dict(roof_id='left',top_z_mm=.8,gap_mm=.2,_footprint=polygon_box(0,0,1,.8)),dict(roof_id='right',top_z_mm=.8,gap_mm=.2,_footprint=polygon_box(3,0,4,.8))]
    assert vl.coverage(field,roofs,V2PROFILE)['all_required_covered']
    assert not vl.coverage(field,roofs,V2PROFILE,{'left','right'})['all_required_covered']
    bridges=np.array([[[1,y,1.],[3,y,1.]] for y in [.1,.3,.5,.7]])
    assert vl.coverage(field,roofs,V2PROFILE,{'left','right'},bridges)['all_required_covered']


def test_v2_bridge_requires_anchors_across_whole_region_not_only_centroid():
    missing=polygon_box(0,0,2,2)
    anchors=polygon_box(-1,.9,0,1.1).union(polygon_box(2,.9,3,1.1))
    assert not vl.bridge_evidence(missing,anchors,2.1,.22)['bridgeable']


def test_v2_segmented_trunk_is_printable_and_cuts_create_multiple_fragments():
    from vascular_processing.support_v2_refinement import segmented_pillar
    from vascular_processing.support_v2_removal import fragment_support
    shapes,records,cuts=segmented_pillar([0,0,0],[0,0,24],1.76,V2PROFILE)
    mesh=trimesh.boolean.union(shapes,engine='manifold')
    assert cuts and all(c['role']=='TRUNK_BREAK' for c in cuts)
    assert vg.self_buildability(mesh,None,V2PROFILE)['passed']
    assert len(fragment_support(mesh,cuts))==len(cuts)+1
    assert mesh.volume<1.76**2*24


def test_v2_perimeter_bridge_requires_actual_previous_layer_end_anchors():
    from vascular_processing.support_v2_validation import bridge_preservation
    graph={'primitives':[dict(primitive_id='b',kind='SHORT_BRIDGE',span_mm=2.,bridge_direction_xy=[1.,0.],bbox_mm=[[0,-.5,1.],[4,.5,1.4]])]}
    paths=np.array([[[0,0,1.],[1,0,1.]],[[3,0,1.],[4,0,1.]],[[.5,0,1.1],[3.5,0,1.1]]])
    r=bridge_preservation(np.empty((0,2,3)),graph,paths,V2PROFILE)
    assert r['passed'] and r['rows'][0]['actual_two_end_anchored_span_extrusions']>0
    assert not bridge_preservation(np.empty((0,2,3)),graph,paths[[0,2]],V2PROFILE)['passed']


def test_v2_foot_box_fast_rejection_matches_unapproved_geometry():
    from vascular_processing.support_custom_qc import intersection_volume
    mesh,approved,_=vp.breakaway_foot([0,0],1.,1.76,V2PROFILE,'BOX_BOTTOM','r')
    upper=vp.foot_above_approved_pad([0,0],1.,1.76,V2PROFILE,'BOX_BOTTOM')
    shallow=vp.loft([0,0,0],[0,0,1.2],5.)
    intruding=vp.loft([0,0,0],[0,0,1.6],5.)
    assert intersection_volume(upper,shallow)==0.
    assert intersection_volume(upper,intruding)>0.
    assert intersection_volume(mesh,shallow)==pytest.approx(intersection_volume(approved,shallow),abs=1e-6)


def test_v2_virtual_fracture_does_not_certify_an_unweakened_neighbour():
    from vascular_processing.support_v2_refinement import segmented_pillar
    from vascular_processing.support_v2_removal import break_design_audit
    shapes,records,cuts=segmented_pillar([0,0,0],[0,0,12],1.76,V2PROFILE)
    mesh=trimesh.boolean.union(shapes,engine='manifold');necks=[m for m,r in zip(shapes,records) if r['kind']=='BREAKAWAY_SEGMENT_NECK']
    assert break_design_audit(mesh,cuts,necks)['passed']
    filled=trimesh.boolean.union([mesh,vp.pillar([.4,0,0],[.4,0,12],1.76)],engine='manifold')
    assert not break_design_audit(filled,cuts,necks)['passed']


def test_v2_fragment_export_failure_remains_a_failed_gate(monkeypatch):
    from vascular_processing import support_v2_removal as removal
    from vascular_processing.support_v2_pipeline import gates
    def failed_export(mesh):
        raise ValueError('SUPPORT_STL_PRECISION_REPAIR_FAILED')
    monkeypatch.setattr(removal,'support_stl',failed_export)
    fragments=removal.fragment_support(vp.pillar([0,0,0],[0,0,2],1.),[])
    assert len(fragments)==1 and fragments[0].is_volume
    assert fragments[0].metadata['fragment_stl_precision_failure']=='SUPPORT_STL_PRECISION_REPAIR_FAILED'
    assert 'FRAGMENT_STL_PRECISION_FAILURE' in gates({'removal':{'precision_failed_fragments':1}})


def test_v2_internal_fragment_must_completely_exit_casting_top(tmp_path):
    from vascular_processing.support_v2_removal import opening_exit_audit,cavity_overlap
    from vascular_processing.support_v2_pipeline import gates
    folder=tmp_path/'removal_fragments';folder.mkdir()
    mesh=vp.pillar([0,0,0],[0,0,2],1.);mesh.export(folder/'fragment_0000.stl')
    initial=np.eye(4);final=np.eye(4);final[2,3]=7.
    inner=[[-4,-4,-1],[4,4,10]];box_bounds=[[-5,-5,-2],[5,5,10]]
    r={'segments':[dict(segment_id='S0000',status='REMOVAL_PATH_FOUND',bounds_mm=mesh.bounds.tolist(),path={'transforms':[initial.tolist(),final.tolist()]})]}
    check=opening_exit_audit(tmp_path,r,inner,box_bounds)
    assert not check['passed']
    assert 'REMOVAL_PATH_DOES_NOT_FULLY_EXIT_CASTING_OPENING' in gates({'removal':{'opening_exit_audit':check}})
    final[2,3]=11.;r['segments'][0]['path']['transforms'][-1]=final.tolist()
    assert opening_exit_audit(tmp_path,r,inner,box_bounds)['passed']
    assert cavity_overlap([[3,-1,1],[8,1,2]],inner)  # centroid outside, but straddles cavity


def test_v2_coupon_calibration_requires_matching_speed_and_temperature():
    from vascular_processing.support_v2_inputs import calibration_signature,calibrated_span
    profile=dict(printer='X1 Carbon',process_name='Vessel',nozzle_diameter=.2,layer_height=.1,line_width=.22,material='ABS',bridge_settings={'bridge_flow':1},speeds={'bridge_speed':20},temperatures={'nozzle_temperature':265})
    calibration=dict(physical_test_passed=True,profile=calibration_signature(profile),maximum_bridge_span_mm=4.)
    assert calibrated_span(calibration,profile)==4.
    assert calibrated_span(calibration,{**profile,'speeds':{'bridge_speed':50}}) is None
    assert calibrated_span(calibration,{**profile,'temperatures':{'nozzle_temperature':270}}) is None
    assert calibrated_span({**calibration,'physical_test_passed':False},profile) is None
    assert calibrated_span({**calibration,'maximum_bridge_span_mm':float('inf')},profile) is None


def test_v2_break_necks_do_not_fake_shorter_unbraced_columns():
    from vascular_processing.support_v2_sensitivity import unbraced_trunk_metrics
    from vascular_processing.support_v2_refinement import segmented_pillar
    before=[dict(root_id='r',kind='VERTICAL_PILLAR',width_mm=1.76,bbox_mm=[[-.88,-.88,0],[.88,.88,24]])]
    shapes,records,_=segmented_pillar([0,0,0],[0,0,24],1.76,V2PROFILE)
    after=[{**r,'root_id':'r','bbox_mm':s.bounds.tolist()} for r,s in zip(records,shapes)]
    a=unbraced_trunk_metrics(before);b=unbraced_trunk_metrics(after)
    assert b['max_unbraced_height_mm']==a['max_unbraced_height_mm']
    assert b['max_slenderness']>=a['max_slenderness']

# V3 regressions are appended; V1/V2 cases above remain permanent.
from vascular_processing.support_v3_fixtures import NAMES as V3_FIXTURES, fixture as v3_fixture
from vascular_processing.support_v3_geometry import curtain as v3_curtain, cut_section as v3_cut_section
V3PROFILE={**V2PROFILE,'roof_gap_mm':.20,'bridge_limit_mm':0.}

@pytest.mark.parametrize('name',V3_FIXTURES)
def test_v3_fixture_closed_quantized_and_segmented(name):
    rib,_=v3_fixture(name,V3PROFILE)
    assert rib['mesh'].is_volume
    assert rib['width_mm']/V3PROFILE['line_width']==pytest.approx(round(rib['width_mm']/V3PROFILE['line_width']))
    assert len(rib['windows'])>=1
    assert all(w['clear_width_mm']>=3. for w in rib['windows'])
    assert len(rib['cuts'])>len(rib['windows'])

@pytest.mark.parametrize('name',['straight_tube_long_rib','curved_tube_segmented_rib','S_tube_rib','rib_self_support','one_line_rib','two_line_rib','rib_roof_transition'])
def test_v3_vertical_curved_ribs_have_no_recursive_support(name):
    from vascular_processing.support_v2_layers import analyse
    rib,_=v3_fixture(name,V3PROFILE)
    field=analyse(rib['mesh'],V3PROFILE,'support',bridge_limit=0.)
    assert not [r for r in field['rows'] if r['classification']=='SUPPORT_REQUIRED']


def test_v3_sparse_pad_actual_coverage_fails_continuous_extrusion_passes():
    from shapely.geometry import box,LineString
    from shapely.ops import unary_union
    from vascular_processing.support_v3_validation import coverage
    polygon=box(0,-.10,20,.10)
    field={'rows':[dict(region_id='tube',classification='SUPPORT_REQUIRED',z_mm=10.,area_mm2=polygon.area)],'polygons':[polygon]}
    sparse={9.8:unary_union([LineString([[x,0],[x+1,0]]).buffer(.11) for x in [0,5,10,15]])}
    continuous={9.8:unary_union([LineString([[x,0],[min(20,x+5.2),0]]).buffer(.11) for x in [0,5,10,15]])}
    old=coverage(field,sparse,V3PROFILE)[0];new=coverage(field,continuous,V3PROFILE)[0]
    assert old['percent']<50 and not old['all_required_covered']
    assert new['percent']==pytest.approx(100.) and new['all_required_covered']


def test_v3_inclined_member_requires_support_but_rib_does_not():
    from vascular_processing.support_v2_layers import analyse
    branch=vp.loft([0,0,0],[12,0,4],.44)
    assert any(r['classification']=='SUPPORT_REQUIRED' for r in analyse(branch,V3PROFILE,'support',bridge_limit=0)['rows'])
    rib,_=v3_fixture('rib_self_support',V3PROFILE)
    assert not analyse(rib['mesh'],V3PROFILE,'support',bridge_limit=0)['rows']


def test_v3_cut_plane_catches_neighbouring_unweakened_member():
    from shapely.geometry import box
    from vascular_processing.support_v3_geometry import union
    intended=vp.loft([0,0,0],[0,0,2],[2,.22]);other=vp.loft([0,1,0],[0,1,2],[2,.66])
    aperture=box(-2,-2,2,2)
    assert v3_cut_section(intended,intended,[0,0,1],[0,0,1],aperture)['passed']
    bad=v3_cut_section(union([intended,other]),intended,[0,0,1],[0,0,1],aperture)
    assert not bad['passed'] and bad['status']=='CUT_SECTION_CROSSES_UNWEAKENED_SUPPORT'
    assert bad['other_support_cross_section_area_mm2']>.5


def test_v3_window_is_void_but_roof_crown_is_continuous():
    rib,_=v3_fixture('rib_window',V3PROFILE);window=rib['windows'][0]
    x,y,z=window['xyz'];q=np.array([[x,y,window['clear_shoulder_z_mm']-.2],[x,y,z+.15]])
    assert rib['mesh'].contains(q).tolist()==[False,True]


def test_v3_top_tool_access_is_separate_from_fragment_removal():
    from vascular_processing.support_v3_removal import probe_access
    box=vp.loft([0,0,-1],[0,0,0],[20,20]);obstacle=vp.loft([0,0,0],[0,0,5],[.44,5])
    access=probe_access([0,0,1],[1,0,0],.22,box,[box,obstacle])
    assert access['passed'] and access['physical_tool_validated'] is False
    assert access['diameter_mm']==3


def test_v3_trapped_fragment_regression():
    from vascular_processing.support_v2_removal import paths_for
    from vascular_processing.support_v3_geometry import union
    from vascular_processing.support_removal import swept_path
    # Enclosed fragment cannot cross the top plate; regression must never
    # accept a clear endpoint after teleporting through an obstacle.
    fragment=vp.loft([0,0,0],[0,0,4],[6,.44])
    lid=vp.loft([0,0,5],[0,0,6],[40,40]);floor=vp.loft([0,0,-2],[0,0,-1],[40,40])
    local,paths=paths_for(fragment,union([floor,lid]),False)
    checks=[swept_path(local,[lid,floor],poses) for _,poses in paths]
    status='REMOVAL_PATH_FOUND' if any(r['passed'] for r in checks) else 'TRAPPED_FRAGMENT'
    assert status=='TRAPPED_FRAGMENT'


def test_v3_gap_uncalibrated_does_not_itself_fail_software_gate():
    from vascular_processing.support_v3_pipeline import gates
    record=dict(build_volume_fit=True,scale_factor=1.,vascular_intersection_mm3=0.,cuts_valid=True,source_unchanged=True,
        self_buildability={'unsupported_region_count':0},removal={'tool_access_passed':True,'trapped_pieces':0,'all_internal_exit_top':True},
        slice={'actual_vascular_coverage':{'all_required_covered':True},'actual_box_coverage':{'all_required_covered':True},'slicer_success':True,
            'support_only_slicer_success':True,'auto_support_off':True,'native_support_toolpaths':0,'all_features_preserved':True,'unresolved_warning_count':0,'actual_self_buildability':{'unsupported_region_count':0}},
        physical_gap_status='PHYSICAL_GAP_UNCALIBRATED',physical_breakaway='UNVALIDATED')
    assert gates(record)==[]
    record['slice']['actual_vascular_coverage']={'all_required_covered':False,'percent':88.}
    assert 'ACTUAL_SUPPORT_COVERAGE_INCOMPLETE' in gates(record)


def test_v3_actual_model_paths_outside_rib_cannot_count_as_support():
    from vascular_processing.support_v3_validation import extrusion_sections
    rib=vp.loft([0,0,0],[0,0,2],[10,.44])
    segments=np.array([[[-4,0,1.02],[4,0,1.02]],[[-4,5,1.02],[4,5,1.02]]])
    footprints,_=extrusion_sections(segments,rib,V3PROFILE)
    from shapely.geometry import Point
    assert footprints[1.02].intersects(Point(0,0))
    assert not footprints[1.02].intersects(Point(0,5))


def test_v3_actual_extrusion_layer_growth_detects_new_floating_material():
    from shapely.geometry import box
    from vascular_processing.support_v3_validation import actual_self_support
    plate=box(0,0,10,.66)
    stable={.12:plate,.22:plate,.32:plate}
    assert actual_self_support(stable,{},V3PROFILE)['passed']
    floating={**stable,.42:box(0,0,10,1.2)}
    result=actual_self_support(floating,{},V3PROFILE)
    assert not result['passed'] and result['unsupported_region_count']>0


def test_v3_invalid_diagnostic_assembly_is_reported_without_aborting(tmp_path):
    import pickle
    from vascular_processing.support_v3_verify import geometry_audit
    x=np.linspace(0,6,31)
    rib=v3_curtain(np.c_[x,np.zeros_like(x)],np.full(len(x),8.),0,V3PROFILE,2,1.32)
    rib.update(rib_id='test',is_port=False)
    with (tmp_path/'rib_geometry_cache.pkl').open('wb') as f:pickle.dump([rib],f)
    model=vp.loft([3,0,0],[3,0,1],[8,6])
    broken=model.copy();broken.update_faces(np.arange(len(broken.faces)-1))
    for name,mesh in [('custom_support_v3.stl',rib['mesh']),('box.stl',model),('mold_without_support.stl',model),('mold_with_custom_support_v3.stl',broken),
                      ('support_ribs.stl',rib['ribs']),('support_roofs.stl',rib['roofs']),('support_cut_necks.stl',rib['necks']),('support_exterior.stl',rib['mesh'])]:mesh.export(tmp_path/name)
    record={'folder':str(tmp_path),'orientation':{'print_transform':np.eye(4).tolist(),'casting_restore_transform':np.eye(4).tolist()}}
    inputs={'inner_bounds':[[0,-2,1],[6,2,8]],'outer_bounds':[[-1,-3,0],[7,3,8]]}
    audit=geometry_audit(record,inputs)
    assert not audit['passed'] and audit['frozen_model_removed_mm3'] is None
    assert any(r['file']=='mold_with_custom_support_v3.stl' and not r['closed_positive_volume'] for r in audit['exports'])

# V3.1 additions: all V1/V2/V3 tests above remain unchanged.
from vascular_processing import support_v3_1_orientation as v31o
from vascular_processing import support_v3_1_residual as v31r
from vascular_processing import support_v3_1_fixtures as v31f
from vascular_processing.support_v3_1_release import plate_release as v31release


def test_v31_exact_inversion_matrix():
    r=v31o.rotation_for('INV_0');assert np.allclose(r.T@r,np.eye(3));assert np.linalg.det(r)==pytest.approx(1)
    assert np.allclose(r@[0,0,1],[0,0,-1])


def test_v31_open_normal_down():
    m=np.eye(4);m[:3,:3]=v31o.rotation_for('INV_0');a=v31o.orientation_audit('INV_0',m)
    assert a['label_matrix_matches'] and a['angle_to_print_down_deg']==pytest.approx(0)
    assert a['bottom_normal_print'][2]==pytest.approx(1)


@pytest.mark.parametrize('axis',['X','Y'])
@pytest.mark.parametrize('tilt',[15,30,45])
def test_v31_casting_local_tilt(axis,tilt):
    from scipy.spatial.transform import Rotation
    name=f'INV_{axis}P{tilt}';m=np.eye(4);m[:3,:3]=v31o.rotation_for(name)
    assert np.allclose(m[:3,:3],v31o.rotation_for('INV_0')@Rotation.from_euler(axis.lower(),tilt,degrees=True).as_matrix())
    a=v31o.orientation_audit(name,m);assert a['angle_to_print_down_deg']==pytest.approx(tilt);assert a['label_matrix_matches']


@pytest.mark.parametrize('axis',['X','Y'])
def test_v31_plus_minus_axis(axis):
    p=v31o.rotation_for(f'INV_{axis}P30')@[0,0,1];m=v31o.rotation_for(f'INV_{axis}M30')@[0,0,1]
    assert np.allclose(p[:2],-m[:2]);assert p[2]==pytest.approx(m[2])


def test_v31_orientation_label_matrix_mismatch():
    m=np.eye(4);m[:3,:3]=v31o.rotation_for('INV_XM30')
    assert v31o.orientation_audit('INV_XP30',m)['orientation_status']=='ORIENTATION_LABEL_MATRIX_MISMATCH'


def test_v31_cavity_visibility_requires_real_opening():
    box=v31f.inverted_box();m=np.diag([1.,-1.,-1.,1.]);m[2,3]=6.
    # Casting cavity z=1..6 -> inverted print cavity z=0..5.
    result=v31o.cavity_visibility(box,[[-5,-5,1],[5,5,6]],m,samples=5)
    assert result['fraction']==1
    upright=v31o.cavity_visibility(box,[[-5,-5,0],[5,5,5]],np.eye(4),samples=5)
    assert upright['fraction']==0  # a free ray that does not cross the mouth is insufficient


def test_v31_direct_buildplate_root_ray():
    box=v31f.inverted_box();assert v31o.vertical_clear(box,[[0,0,4]])[0]
    obstacle=v31f.block([2,2,1],[0,0,2]);assert not v31o.vertical_clear(obstacle,[[0,0,4]])[0]


def test_v31_box_ceiling_detection():
    mesh=v31f.inverted_box();down=mesh.face_normals[:,2]<-.71
    inside=(mesh.triangles_center[:,2]>4.9)&(mesh.triangles_center[:,2]<5.1)
    assert mesh.area_faces[down&inside].sum()==pytest.approx(100)


def test_v31_plate_release_translation():
    result=v31release(v31f.inverted_box(),v31f.fixture('plate_release_pass'))
    assert result['passed'];assert result['level']==1;assert result['path']['full_mold_exit']


def test_v31_plate_release_tilt():
    mold,support=v31f.tilt_fixture();result=v31release(mold,support)
    assert not result['attempts'][0]['passed'];assert result['passed'];assert result['level']==2


def test_v31_plate_release_trapped():
    from vascular_processing.support_v3_geometry import union
    mold=union([v31f.inverted_box(),v31f.block([6,6,.5],[0,0,2])])
    assert not v31release(mold,v31f.fixture('plate_release_fail'))['passed']


def test_v31_fixture_files_closed():
    for p in (DATA/'v3_1').glob('*.stl'):assert trimesh.load_mesh(p).is_volume
    assert len(list((DATA/'v3_1').glob('*.stl')))==6


def test_v31_polygon_coverage_difference_hole_preserved():
    required=polygon_box(0,0,3,3);actual=required.difference(polygon_box(1,1,2,2))
    assert v31r.difference(required,actual).area==pytest.approx(1)
    assert v31r.difference(required,required).is_empty


def test_v31_residual_grain_clustering_branches_and_adjacency():
    branches=[dict(branch_id='a',xyz=np.array([[0,0,1],[0,1,1]]),arc=np.array([0.,1.]),tangent=np.array([[0,1,0],[0,1,0]])),
              dict(branch_id='b',xyz=np.array([[4,0,1],[4,1,1]]),arc=np.array([0.,1.]),tangent=np.array([[0,1,0],[0,1,0]]))]
    rows=[dict(component='VASCULAR',z_mm=z,uncovered_wkt=polygon_box(x,0,x+.1,.1).wkt) for x,z in [(0,1),(0,1.1),(4,1),(0,3)]]
    grains=v31r.cluster(rows,branches,dict(PROFILE,roof_gap_mm=.2))
    assert len(grains)==3;assert sorted(len(g['source_rows']) for g in grains)==[1,1,2]
    assert {'a','b'}=={g['branch'] for g in grains}


def test_v31_rib_width_quantization():
    assert v31r.quantized_width(.45,.22)==pytest.approx(.66)
    assert v31r.quantized_width(.66,.22)==pytest.approx(.66)


def test_v31_one_line_local_roof_extension():
    from vascular_processing.support_v3_geometry import curtain
    from vascular_processing.support_v3_1_repair import local_variant
    profile=dict(PROFILE,roof_gap_mm=.2,cut_neck_lines=2)
    rib=curtain(np.c_[np.linspace(0,10,21),np.zeros(21)],np.full(21,6.),0.,profile,3,1.54,5,3)
    rib.update(rib_id='x',interval_id='x',branch_id='x',is_port=False)
    g=dict(footprint_wkt=polygon_box(1,.6,2,1).wkt,bbox=[[1,.6,6.1],[2,1,6.2]],component='VASCULAR')
    t=local_variant(rib,g,'EXTEND_ROOF',profile,{})
    extra=trimesh.boolean.difference([t['mesh'],rib['mesh']],engine='manifold')
    assert extra.volume>0;assert extra.bounds[0,0]>=1-2*.22-.001;assert extra.bounds[1,0]<=2+2*.22+.001
    assert extra.bounds[1,1]<=1.54/2+.22/2+.001


def test_v31_top_down_propagation_stops_at_actual_material():
    p=polygon_box(0,0,1,1);profile=dict(PROFILE,roof_gap_mm=.2)
    result=v31r.top_down_layers(p,1.,{.6:p},profile)
    assert result['connected_to']=='ACTUAL_PREVIOUS_MATERIAL';assert result['base_z_mm']==pytest.approx(.6)
    no_material=v31r.top_down_layers(p,1.,{},profile)
    assert no_material['base_z_mm']==0;assert no_material['connected_to']=='LEGAL_ROOT_REQUIRED'


def test_v31_support_self_closure_detects_lateral_island():
    from vascular_processing.support_v3_validation import actual_self_support
    p=polygon_box(0,0,1,1);q=p.union(polygon_box(3,0,4,1))
    profile=dict(PROFILE,roof_gap_mm=.2)
    assert actual_self_support({.12:p,.22:q},{},profile)['unsupported_region_count']==1
    assert actual_self_support({.12:q,.22:q},{},profile)['unsupported_region_count']==0


def test_v31_repair_iteration_convergence():
    assert v31r.iteration_stop([1,.5],True)=='ALL_COVERAGE_PASSED'
    assert v31r.iteration_stop([1,.5,.2]) is None
    assert v31r.iteration_stop([10,8,6,4,3,2])=='REPAIR_ITERATION_LIMIT'


def test_v31_repair_stagnation():
    assert v31r.iteration_stop([1,.995,.994])=='COVERAGE_REPAIR_STAGNATED'
    assert v31r.iteration_stop([1,1.1,1.2])=='COVERAGE_REPAIR_STAGNATED'
    assert v31r.iteration_stop([1,.5,.499]) is None


def test_v31_removal_aware_repair_rejects_vessel_intersection():
    from vascular_processing.support_v3_geometry import curtain
    from vascular_processing.support_v3_1_repair import candidate_preflight
    p=dict(PROFILE,roof_gap_mm=.2,cut_neck_lines=2)
    rib=curtain([[0,0],[3,0]],[5,5],0.,p,3,1.54)
    v=v31f.block([2,2,2],[1,0,3]);box=v31f.block([20,20,1],[0,0,-1])
    r=candidate_preflight(rib,trimesh.Trimesh(),v,box,np.eye(4),[[-10,-10,0],[10,10,6]],p)
    assert not r['passed'];assert r['reason']=='FROZEN_VASCULAR_INTERSECTION'


def test_v31_cut_point_relocation_changes_real_mesh():
    from vascular_processing.support_v3_geometry import curtain
    from vascular_processing.support_v3_1_repair import move_cut_up
    p=dict(PROFILE,roof_gap_mm=.2,cut_neck_lines=2)
    rib=curtain([[0,0],[3,0]],[6,6],0.,p,3,1.54)
    old=np.array(rib['cuts'][0]['origin']);trial=move_cut_up(rib,0,.88,p)
    assert trial['cuts'][0]['origin'][2]==pytest.approx(old[2]+.88)
    assert trial['mesh'].volume<rib['mesh'].volume;assert trial['cuts'][0]['intended'].volume>0
    assert trial['approved_retained_stubs'].volume==pytest.approx(rib['approved_retained_stubs'].volume)


def test_v31_tolerance_derived_and_never_relaxed():
    assert v31r.area_tolerance(.22,150)<=1e-5
    assert v31r.area_tolerance(.22,1)<v31r.area_tolerance(.22,150)


def test_v31_void_shell_is_not_free_material_component():
    from vascular_processing.support_v3_1_validation import solid_components
    outer=trimesh.creation.box([4,4,4]);cavity=trimesh.creation.box([1,1,1]);cavity.invert()
    solid=trimesh.util.concatenate([outer,cavity]);parts=solid_components(solid)
    assert len(parts)==1;assert parts[0].volume==pytest.approx(63)
    assert len(parts[0].faces)==len(solid.faces)


def test_v31_actual_brim_role_filter_keeps_coordinates(tmp_path):
    from vascular_processing.support_v3_1_validation import actual_brim
    p=tmp_path/'brim.3mf'
    with zipfile.ZipFile(p,'w') as z:
        z.writestr('Metadata/project_settings.config',json.dumps({'extruder_offset':['0x0']}))
        z.writestr('Metadata/plate_1.gcode','M83\nG90\nG1 X0 Y0 Z0.12\n; FEATURE: Outer wall\nG1 X10 E1\n; FEATURE: Brim\nG1 Y10 E1\n; FEATURE: Outer wall\nG1 X0 E1\n')
    brim=actual_brim(p,PROFILE,.25)
    assert brim.area==pytest.approx(2.5);assert brim.bounds[0]==pytest.approx(9.875)


def test_v31_tilt_release_cannot_cross_buildplate():
    from vascular_processing.support_v3_1_release import release_poses,bed_clearance
    mold,support=v31f.tilt_fixture()
    for axis in 'xy':
        for angle in [-15,15]:
            local,poses=release_poses(mold,support,axis,angle)
            assert bed_clearance(local,poses)['passed']


def test_v31_bed_check_rejects_rotation_without_lift():
    from vascular_processing.support_v3_1_release import bed_clearance
    from scipy.spatial.transform import Rotation
    local=v31f.inverted_box();a=np.eye(4);b=np.eye(4);b[:3,:3]=Rotation.from_euler('x',15,degrees=True).as_matrix()
    assert not bed_clearance(local,[a,b])['passed']


def test_v31_self_repair_ends_in_previous_layer():
    from vascular_processing.support_v3_geometry import curtain
    from vascular_processing.support_v3_1_repair import local_variant
    p=dict(PROFILE,roof_gap_mm=.2,cut_neck_lines=2)
    rib=curtain([[0,0],[3,0]],[6,6],0.,p,3,1.54)
    rib.update(rib_id='a',interval_id='a',branch_id='a',is_port=True)
    grain=dict(footprint_wkt=polygon_box(1,.54,2,.65).wkt,bbox=[[1,.54,.22],[2,.65,.22]],component='SUPPORT_SELF')
    trial=local_variant(rib,grain,'TOP_DOWN_THIN_WEB',p,{})
    added=trimesh.boolean.difference([trial['mesh'],rib['mesh']],engine='manifold')
    assert added.volume>0;assert added.bounds[1,2]==pytest.approx(.12)

# RFUS: retained historical tests above; no sparse selection in production.
def rfus_profile():
    return dict(line_width=.22,layer_height=.1,first_layer_height=.12,roof_gap_mm=.2)


def rfus_config():
    import yaml
    from vascular_processing.support_inputs import ROOT
    return yaml.safe_load((ROOT/'config/robust_underbody.yaml').read_text())


@pytest.mark.parametrize('name',['full_support_straight_tube','full_support_curved_tube','full_support_S_tube','full_support_Y_branch','full_support_external_port','vertical_overlap_vessels','curtain_lateral_offset','twin_curtain','segmented_curtain','staggered_seam','removal_handle','handle_collision','base_neck_isolation','curtain_removal','trapped_curtain'])
def test_rfus_synthetic_fixture_closed(name):
    from vascular_processing.support_rfus_fixtures import fixture
    f=fixture(name,rfus_profile())
    assert f['support']['mesh'].is_volume and f['vascular'].is_volume


@pytest.mark.parametrize('name',['full_support_straight_tube','full_support_curved_tube','full_support_S_tube','full_support_external_port'])
def test_rfus_every_arc_sample_is_target(name):
    from vascular_processing.support_rfus_fixtures import fixture
    from vascular_processing.support_rfus_geometry import corridor
    f=fixture(name,rfus_profile());b=f['branches'][0]
    c=corridor(b,f['vascular'],np.array([[80,80,0],[130,130,20]]),rfus_profile(),rfus_config())
    assert c['arc'][-1]==pytest.approx(b['arc'][-1])
    assert np.diff(c['arc']).max()<=rfus_profile()['line_width']
    assert len(c['arc'])>200


def test_rfus_lower_envelope_uses_first_surface_not_lower_crossing():
    from vascular_processing.support_rfus_fixtures import tube
    from vascular_processing.support_rfus_geometry import lower_points
    from vascular_processing.support_v3_geometry import union
    m=union([tube([[0,0,8],[10,0,8]]),tube([[5,-5,4],[5,5,4]])])
    p=lower_points(m,np.array([[5.,0,8.]]))
    assert p[0,2]==pytest.approx(7.4,abs=.01)


def test_rfus_vertical_overlap_requires_offset_or_twin():
    from vascular_processing.support_rfus_fixtures import fixture
    from vascular_processing.support_rfus_geometry import route_segment
    from vascular_processing.support_v3_geometry import union
    p=rfus_profile();f=fixture('vertical_overlap_vessels',p);s=f['support'];v=union([f['vascular']]+f['obstacles'])
    box=trimesh.creation.box([5,5,1]);box.apply_translation([70,70,-1])
    routed,attempts=route_segment(s['path_xy'],s['top'],s['base'],s['lateral'],p,rfus_config(),v,box,2,2,1)
    assert not attempts[0]['passed']
    assert routed is not None and routed['route'] in ('LATERAL_OFFSET','TWIN_CURTAINS','LOCAL_PILLAR_TRANSITION','TWIN_CURTAINS_AROUND_CROSSING')
    assert routed['mesh'].is_volume
    assert trimesh.boolean.intersection([routed['mesh'],v],engine='manifold').volume<1e-5


def test_rfus_twin_curtain_retains_central_roof():
    from vascular_processing.support_rfus_fixtures import fixture
    f=fixture('twin_curtain',rfus_profile());s=f['support']
    assert s['twin'] and s['roof_width_mm']==pytest.approx(.44)
    assert s['roofs'].bounds[0,1]<100<s['roofs'].bounds[1,1]
    assert len(s['mesh'].split())==1


def test_rfus_seams_split_whole_height():
    from vascular_processing.support_rfus_fixtures import fixture
    for name in ['segmented_curtain','staggered_seam']:
        s=fixture(name,rfus_profile())['support']['mesh']
        assert len(s.split())==2


def test_rfus_handle_collision_detected():
    from vascular_processing.support_rfus_fixtures import fixture
    from vascular_processing.support_custom_qc import intersection_volume
    f=fixture('handle_collision',rfus_profile())
    assert intersection_volume(f['support']['handles'],f['obstacles'][0])>0


def test_rfus_base_neck_bypass_rejected():
    from vascular_processing.support_rfus_fixtures import fixture
    from vascular_processing.support_rfus_geometry import neck_audit
    from vascular_processing.support_v3_geometry import union
    f=fixture('base_neck_isolation',rfus_profile());s=f['support']
    assert neck_audit(s)['passed']
    block=trimesh.creation.box([2,1,2]);block.apply_translation([105,100,1]);s['mesh']=union([s['mesh'],block])
    assert not neck_audit(s)['passed']


def test_rfus_curtain_removal_sideways_then_up():
    from vascular_processing.support_rfus_fixtures import fixture
    from vascular_processing.support_rfus_geometry import disconnect_base
    from vascular_processing.support_rfus_removal import paths_for
    from vascular_processing.support_removal import swept_path
    f=fixture('curtain_removal',rfus_profile());s=f['support'];p=disconnect_base(s);local,paths=paths_for(s,p,15)
    assert not swept_path(local,[f['vascular']],paths[0][1])['passed']
    assert any(swept_path(local,[f['vascular']],poses)['passed'] for _,poses in paths[1:])


def test_rfus_trapped_curtain_has_no_vertical_release():
    from vascular_processing.support_rfus_fixtures import fixture
    from vascular_processing.support_rfus_geometry import disconnect_base
    from vascular_processing.support_rfus_removal import paths_for
    from vascular_processing.support_removal import swept_path
    f=fixture('trapped_curtain',rfus_profile());local,paths=paths_for(f['support'],disconnect_base(f['support']),15)
    assert not swept_path(local,f['obstacles'],paths[0][1])['passed']


def test_rfus_no_nominal_stl_credit_and_nonzero_missing_arc():
    from vascular_processing.support_rfus_fixtures import fixture
    from vascular_processing.support_rfus_geometry import corridor
    from vascular_processing.support_rfus_validation import underbody_coverage
    f=fixture('full_support_straight_tube',rfus_profile());b=corridor(f['branches'][0],f['vascular'],np.array([[80,80,0],[130,130,20]]),rfus_profile(),rfus_config())
    r=underbody_coverage([b],{}, {},rfus_profile())
    assert r['percent']==0 and r['max_uncovered_underbody_arc_mm']==pytest.approx(30)


def test_rfus_actual_width_and_travel_e_modes(tmp_path):
    from vascular_processing.support_rfus_travel import events
    import zipfile,json
    p=tmp_path/'paths.3mf'
    text='G90\nM82\nG92 E0\n; FEATURE: Outer wall\n; LINE_WIDTH: 0.31\nG1 X0 Y0 Z1\nG1 X2 E1\nG1 X4 E1\nG92 E0\nM83\nG1 X6 E-0.3\n; LINE_WIDTH: 0.23\nG1 X8 E0.4\n'
    with zipfile.ZipFile(p,'w') as z:
        z.writestr('Metadata/project_settings.config',json.dumps({'extruder_offset':['1x2'],'line_width':.22,'layer_height':.1,'nozzle_diameter':[.2]}));z.writestr('plate.gcode',text)
    r=events(p)
    assert len(r['extrusion'])==2 and len(r['travel'])==2
    assert r['extrusion_widths'].tolist()==pytest.approx([.31,.23])
    assert r['extrusion'][0,0].tolist()==pytest.approx([1,2,1])
    assert r['travel_ids'][0]>r['extrusion_ids'][0]


def test_rfus_variable_width_changes_actual_footprint():
    from vascular_processing.support_rfus_travel import actual_sections
    s=trimesh.creation.box([10,1,1]);s.apply_translation([5,0,.5]);p=rfus_profile()
    data={'features':['Outer wall'],'extrusion':np.array([[[1,0,.52],[9,0,.52]]]),'extrusion_widths':np.array([.4])}
    a=actual_sections(data,s,p)[.52]
    data['extrusion_widths']=np.array([.2]);b=actual_sections(data,s,p)[.52]
    assert a.area==pytest.approx(2*b.area)


def test_rfus_rejects_single_line_tall_curtain():
    from vascular_processing.support_rfus_geometry import make_segment
    with pytest.raises(ValueError,match='TWO_LINE'):
        make_segment(np.array([[0,0],[10,0]]),np.array([10,10]),0,np.array([[0,1],[0,1]]),rfus_profile(),1)


def test_rfus_y_parent_and_both_daughters_remain_targets():
    from vascular_processing.support_rfus_fixtures import fixture
    f=fixture('full_support_Y_branch',rfus_profile())
    assert {b['branch_id'] for b in f['branches']}=={'Y_PARENT','Y_DAUGHTER_1','Y_DAUGHTER_2'}
    assert sum(b['arc'][-1] for b in f['branches'])>48
    assert len(f['support']['mesh'].split())>=2


def test_rfus_double_precision_boolean_accepts_readonly_and_empty():
    from vascular_processing.support_rfus_geometry import boolean
    m=trimesh.creation.box([1,1,1]);m.vertices.flags.writeable=False
    n=boolean([m,m.copy()])
    assert n.is_volume and n.volume==pytest.approx(1)
    assert not len(boolean([n,trimesh.Trimesh()],'intersection').faces)
    assert boolean([n,trimesh.Trimesh()],'difference').volume==pytest.approx(1)


def test_rfus_diagnostic_mesh_failure_cannot_pass_production():
    from vascular_processing.support_rfus_validation import gates
    r=dict(source_unchanged=True,scale_factor=1.,orientation={'orientation_id':'O000','print_transform':np.eye(4).tolist()},support_export_error='NONMANIFOLD_EXPORT')
    assert 'SUPPORT_EXPORT_FAILED' in gates(r)


def test_rfus_cap_plane_rays_keep_both_arc_endpoints():
    from vascular_processing.support_rfus_fixtures import tube
    from vascular_processing.support_rfus_geometry import lower_points
    xyz=np.array([[90.,100,10],[120.,100,10]])
    lower=lower_points(tube(xyz),xyz)
    assert np.isfinite(lower).all()
    assert np.max(np.linalg.norm(lower[:,:2]-xyz[:,:2],axis=1))<=np.spacing(np.float32(120.))*1.01
    assert lower[:,2].tolist()==pytest.approx([9.4,9.4],abs=1e-5)


def test_rfus_native_comparison_uses_actual_support_layer_height():
    from shapely.geometry import box
    from vascular_processing.support_rfus_review import native_self_support
    parsed={'extrusion':np.array([[[0,0,.32],[10,0,.32]]]),'extrusion_widths':np.array([.22]),'extrusion_heights':np.array([.2])}
    r=native_self_support(parsed,np.array([True]),{.12:box(-1,-1,11,1)}, {},rfus_profile())
    assert r['unsupported_region_count']==0


def test_rfus_detached_sliver_is_not_one_removable_curtain():
    from vascular_processing.support_rfus_fixtures import fixture
    from vascular_processing.support_rfus_removal import piece_topology
    from vascular_processing.support_rfus_geometry import union
    s=fixture('curtain_removal',rfus_profile())['support']
    assert piece_topology(s)['removable_as_one_component']
    sliver=trimesh.creation.box([.1,.1,1]);sliver.apply_translation([105,105,5])
    s['mesh']=union([s['mesh'],sliver])
    r=piece_topology(s)
    assert r['positive_components']==2 and not r['removable_as_one_component']


def test_rfus_multiple_curtains_require_ordered_removal():
    from vascular_processing.support_rfus_fixtures import fixture
    from vascular_processing.support_rfus_geometry import disconnect_base
    from vascular_processing.support_rfus_removal import paths_for
    from vascular_processing.support_removal import swept_path
    import copy
    f=fixture('curtain_removal',rfus_profile());b=f['support'];a=copy.deepcopy(b)
    for key in ['mesh','curtains','roofs','feet','necks']:
        a[key].apply_translation([0,1.5,0])
    a['path_xy']+=np.array([0,1.5]);va=f['vascular'].copy();va.apply_translation([0,1.5,0])
    la,pa=paths_for(a,disconnect_base(a),15);lb,pb=paths_for(b,disconnect_base(b),15)
    route=lambda paths:next(p for name,p in paths if name=='LATERAL_+3_THEN_TOP')
    assert not swept_path(lb,[f['vascular'],va,a['mesh']],route(pb))['passed']
    assert swept_path(la,[f['vascular'],va,b['mesh']],route(pa))['passed']
    assert swept_path(lb,[f['vascular'],va],route(pb))['passed']


def test_rfus_deep_curtain_receives_accessible_top_handle():
    from vascular_processing.support_rfus_fixtures import fixture
    from vascular_processing.support_rfus_geometry import make_segment,add_handle
    from vascular_processing.support_v3_removal import probe_access
    p=rfus_profile();f=fixture('removal_handle',p);original=f['support']
    s=make_segment(original['path_xy'],original['top'],2,original['lateral'],p)
    outer=trimesh.creation.box([60,30,25]);outer.apply_translation([105,100,12.5])
    cavity=trimesh.creation.box([56,26,27]);cavity.apply_translation([105,100,15.5])
    box=trimesh.boolean.difference([outer,cavity],engine='manifold')
    cfg={**rfus_config(),'handle_lines':3,'handle_access_fraction':.45}
    assert add_handle(s,f['vascular'],box,[],np.array([[77,87,2],[133,113,25]]),p,cfg)
    assert s['handle_target'][2]>=s['handle_access_zone_min_z_mm']
    access=probe_access(s['handle_target'],np.array([0.,1,0]),s['handle_width_mm']/2,box,[box,f['vascular'],s['mesh']],3)
    assert access['passed']
