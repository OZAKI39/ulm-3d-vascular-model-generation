"""Permanent BCF geometry, negative-path and delivered BG001 artifact regressions."""
import copy,json,zipfile
from pathlib import Path
import numpy as np
import pytest
import trimesh
import yaml
from vascular_processing.support_inputs import ROOT,digest,discover
from vascular_processing.casting_fixture_geometry import rigid,source_layout,build_frame,bounds_box,shape_mesh,mesh_qc
from vascular_processing.casting_fixture_synthetic import CASES,layout_for
from vascular_processing.casting_fixture_audit import body_constraints,perturbation_fit,sweep_translation
from vascular_processing.casting_fixture_pipeline import default_output,dish_check,summary_gates
from vascular_processing.support_slicer_adapter import local_command

CFG=yaml.safe_load((ROOT/'config/vascular_casting_fixture.yaml').read_text())

@pytest.fixture(scope='module')
def synthetic():
    layout=layout_for('four_side_ports');geometry=build_frame(layout,CFG,.22);return layout,geometry

@pytest.fixture(scope='module')
def actual():
    return json.loads((default_output()/'vascular_casting_fixture_summary.json').read_text())

@pytest.mark.parametrize('case',CASES)
def test_synthetic_fixtures_are_real_closed_meshes(case):
    p=ROOT/'tests/data/casting_fixture'/case;assert json.loads((p/'fixture.json').read_text())['synthetic'];assert mesh_qc(trimesh.load_mesh(p/'core.stl'))['passed']

@pytest.mark.parametrize('scale',[.9,1.001,2.])
def test_no_scaling(scale):
    m=np.eye(4);m[0,0]=scale
    with pytest.raises(ValueError,match='SCALING'):rigid(m)

def test_source_discovery_uses_production_not_support_copy():
    x=discover();assert x['selected_candidate']=='candidate_0008';assert 'support_aware_roi_60x50' in str(x['core']);assert 'vascular_support_planner' not in str(x['core'])

def test_core_copy_byte_identical_and_original_sha(actual):
    assert digest(Path(actual['vascular']['source']))==digest(Path(actual['vascular']['independent_core_copy']))==actual['vascular']['source_hash']=='84f1989392fb8d210432750a8022a830943b9bff890b044aa9d310822b619ce0'

def test_real_placement_rigid(actual):
    t=json.loads((default_output()/'assembly/casting_core_transform.json').read_text());rigid(t['transform_4x4']);assert t['scale_factor']==1

@pytest.mark.parametrize('axis,side',[(0,-1),(0,1),(1,-1),(1,1)])
def test_four_walls_present(synthetic,axis,side):
    _,g=synthetic;m=shape_mesh(g['frame']);point=np.array([0.,0.,3.]);point[axis]=side*11;assert m.contains([point])[0]

def test_exact_cavity_empty_certifies_no_floor_or_crossbar(synthetic):
    l,g=synthetic;cavity=bounds_box(l['low']+[1e-5]*3,l['high']-[1e-5]*3)
    assert g['frame'].intersect(cavity).Volume()<1e-9
    for s in g['closures'].values():assert s.intersect(cavity).Volume()<1e-9

def test_closed_bottom_negative_case_is_detected(synthetic):
    l,g=synthetic;bad=g['frame'].fuse(bounds_box([-10,-10,0],[10,10,1]));assert bad.intersect(bounds_box([-9,-9,.1],[9,9,2])).Volume()>200

def test_no_bottom_export_rays_and_top_open(synthetic):
    _,g=synthetic;m=shape_mesh(g['frame']);xy=np.array([[0,0],[-8,8],[8,-8],[9,9]],float)
    for sign in (-1,1):assert not m.ray.intersects_any(np.c_[xy,np.full(4,4.)],np.tile([0,0,sign],(4,1))).any()

def test_outward_flange_and_closed_seal(synthetic):
    l,g=synthetic;assert g['flange_bounds'][0,0]<l['low'][0];assert np.allclose(g['seal_loop'][0],g['seal_loop'][-1]);assert g['seal_shape'].intersect(bounds_box(l['low'],l['high'])).Volume()==pytest.approx(0,abs=1e-8)

def test_minimum_dish_diameter_has_matching_lower_bound(synthetic):
    _,g=synthetic;diagonal=np.linalg.norm(g['flange_bounds'][1,:2]-g['flange_bounds'][0,:2]);assert g['minimum_dish_diameter']==pytest.approx(diagonal);radius=diagonal/2;assert np.linalg.norm(g['flange_bounds'][1,:2])==pytest.approx(radius)

def test_missing_dish_does_not_invent_diameter(tmp_path):
    r=dish_check(tmp_path,CFG,130);assert not r['actual_dimensions_known'] and r['status']=='DISH_FIT_NOT_VERIFIED'

def test_too_small_measured_dish_rejected(tmp_path):
    p=tmp_path/'dish.json';p.write_text(json.dumps(dict(dish_inner_flat_diameter=80,dish_wall_height=20,dish_floor_thickness=2)));assert not dish_check(tmp_path,{**CFG,'dish_measurement_file':str(p)},130)['fit_verified']

def test_detect_any_number_of_ports(synthetic):
    for case,count in [('one_inlet_two_outlets',3),('four_side_ports',4)]:assert len(layout_for(case)['ports'])==count

def test_real_straight_length_safe_distances(actual):
    for p in actual['ports']['designs']:assert p['straight_length_mm']>6 and p['inner_distance_from_transition_mm']>=-1e-5 and p['end_reserve_mm']>=-1e-5

def test_too_short_straight_source_is_rejected(tmp_path):
    inputs=discover();meta=json.loads(inputs['geometry'].read_text());p=meta['ports'][0];p['end']=(np.array(p['target'])+np.array(p['direction'])*1).tolist();f=tmp_path/'geometry.json';f.write_text(json.dumps(meta));inputs['geometry']=f
    with pytest.raises(ValueError,match='PORT_EXTENSION_TOO_SHORT'):source_layout(inputs,CFG,.22)

def test_slots_open_saddles_and_upper_closures(synthetic):
    l,g=synthetic;m=shape_mesh(g['frame'])
    for p in l['ports']:
        c=np.array(p['capture_center']);assert not m.ray.intersects_any([c],[[0,0,1]])[0];assert m.ray.intersects_any([c],[[0,0,-1]])[0]
        closure=shape_mesh(g['closures'][p['port_id']]);assert closure.ray.intersects_any([c],[[0,0,1]])[0]

def test_master_has_one_stop_secondary_floats(actual):
    p=actual['ports'];assert p['master_port']=='O2';assert sum(r['role']=='MASTER_DATUM' for r in p['designs'])==1;assert all(p['axial_float'][k]>=.5 for k in p['secondary_ports'])

def test_worst_case_float_improves_over_all_fixed(synthetic):
    l,_=synthetic;errors=[np.array(p['axis'])*.8 if p['role']!='MASTER_DATUM' else np.zeros(3) for p in l['ports']]
    assert perturbation_fit(l['ports'],errors,True);assert not perturbation_fit(l['ports'],errors,False)

def test_excess_radial_error_not_hidden_by_axial_float(synthetic):
    l,_=synthetic;errors=[np.array([0,0,2 if i in (1,2) else -2]) for i,p in enumerate(l['ports'])];assert not perturbation_fit(l['ports'],errors,True)

def test_real_top_and_all_closure_paths(actual):
    a=actual['assembly'];assert a['top_loading_path']['passed'] and a['closure_insertion_pass'];r=json.loads((default_output()/'tables/assembly_audit.json').read_text());assert all(v['passed'] for v in r['closure_insertions']+r['pin_insertions'])

def test_thin_obstacle_between_endpoints_rejected():
    moving=trimesh.creation.box([1,1,1]);obstacle=trimesh.creation.box([2,2,.03]);obstacle.apply_translation([0,0,2.]);assert not sweep_translation(moving,[obstacle],[0,0,4])['passed']

def test_central_contact_and_cavity_intrusion_absent(actual):
    assert actual['assembly']['central_vascular_contact_count']==0;assert actual['assembly']['illegal_fixture_intrusion_volume']<1e-8;assert actual['assembly']['central_min_clearance_mm']>5

def test_contact_zones_are_straight_only(actual):
    for p in actual['ports']['designs']:
        n=np.array(p['axis']);c=np.array(p['capture_center']);a=np.array(p['straight_start']);b=np.array(p['end']);assert 0<(c-a)@n<(b-a)@n;assert np.linalg.norm((c-a)-((c-a)@n)*n)<1e-5

def test_bottom_wall_clearance(actual):assert actual['assembly']['core_to_dish_clearance_mm']==pytest.approx(10) and actual['assembly']['central_min_clearance_mm']>0

def test_anti_float_uses_no_secondary_axial_lock(synthetic,actual):
    l,_=synthetic;r,_,_=body_constraints(l['ports']);assert r['matrix_rank']==6;assert r['geometric_Z_lock'] and r['pitch_lock'] and r['roll_lock'];assert not r['physical_force_validated'];assert actual['anti_float']['closures_positive_Z_locked'];assert all(r['collision'] for r in actual['anti_float']['finite_pose_checks'])

@pytest.mark.parametrize('group',['frame','closures','pins','vascular'])
def test_delivered_individual_mesh_qc(actual,group):
    q=actual['mesh_qc'][group];assert q['passed'] if group in ('frame','vascular') else all(v['passed'] for v in q.values())

def test_individual_stl_round_trip():
    for p in [default_output()/'frame/bottomless_PLA_casting_frame.stl',*(default_output()/'closures').glob('port_*.stl')]:assert mesh_qc(trimesh.load_mesh(p))['passed']

@pytest.mark.parametrize('part',['frame_slice','closure_slice'])
def test_actual_pla_auto_off_slice_archive(actual,part):
    r=actual['bambu'][part];assert r['slicer_success'] and r['support_free_slice_passed'] and r['native_support_moves']==0;assert r['actual_profile']['filament_type']==['PLA'];assert r['actual_profile']['nozzle_diameter']==['0.2'];assert not r['printer_start'] and not r['printer_upload']
    with zipfile.ZipFile(r['archive']) as z:assert z.testzip() is None

def test_core_only_actual_abs_reference(actual):
    r=actual['bambu']['core_reference_slice'];assert r['slicer_success'];assert r['actual_profile']['filament_type']==['ABS']

@pytest.mark.parametrize('flag',['--send','--upload','--start-print','--print'])
def test_no_printer_actions_allowed(flag):
    with pytest.raises(ValueError,match='NONLOCAL'):local_command(['bambu',flag,'file.3mf'])

def test_cleanup_protected_input_hashes(actual):
    c=json.loads((ROOT/'cleanup_result.json').read_text());assert c['deleted_file_count']>0 and c['protected_hash_mismatches']==0;assert actual['source_protection']['all_unchanged']

def test_status_cannot_pass_missing_closure_validation(actual):
    s=copy.deepcopy(actual);s['assembly']['closure_insertion_pass']=False;assert summary_gates(s)[0]=='NEEDS_ADJUSTMENT'

def test_qc_and_step_all_exist():
    from vascular_processing.casting_fixture_review import NAMES
    for n in NAMES:assert (default_output()/'QC'/(n+'.png')).stat().st_size>20000
    assert (default_output()/'assembly/full_fixture_assembled.step').stat().st_size>1000

# V2 additions. V1 tests above remain intact and refer to the archived V1 run.
from vascular_processing.casting_fixture_v2_pipeline import default_output as v2_output,gates as v2_gates
from vascular_processing.casting_fixture_v2_geometry import build as v2_build,straight_port,port_transform
from vascular_processing.casting_fixture_v2_audit import manifold,mesh as manifold_mesh,choose_strategy,contact_zone
from vascular_processing.casting_fixture_v2_synthetic import CASES as V2_CASES,load as load_v2
from vascular_processing.casting_mold_qc import fcl_distance
V2_CFG=yaml.safe_load((ROOT/'config/internal_wall_supported_fixture_v2.yaml').read_text())

@pytest.fixture(scope='module')
def v2_actual():return json.loads((v2_output()/'internal_wall_supported_fixture_summary.json').read_text())

@pytest.fixture(scope='module')
def v2_meshes():
    o=v2_output();return dict(core=trimesh.load_mesh(o/'vascular_core/assembly_placed_ABS_core.stl'),frame=trimesh.load_mesh(o/'frame/closed_wall_bottomless_frame_with_support_blades.stl'),before=trimesh.load_mesh(o/'diagnostic/PDMS_before_trim.stl'),after=trimesh.load_mesh(o/'diagnostic/PDMS_after_trim.stl'),retainers=[trimesh.load_mesh(p) for p in sorted((o/'retainers').glob('port_*.stl'))])

@pytest.mark.parametrize('case',V2_CASES)
def test_v2_permanent_synthetic_geometry(case):
    l,raw=load_v2(case);assert l['core'].is_watertight;assert l['central'].is_watertight;assert raw

@pytest.mark.parametrize('case,error',[('curved_extension_rejected','CURVED'),('short_extension_rejected','TOO_SHORT')])
def test_v2_reject_unapproved_extension(case,error):
    _,raw=load_v2(case)
    with pytest.raises(ValueError,match=error):straight_port(raw[0],V2_CFG)

@pytest.mark.parametrize('case,expected',[('two_ports_same_wall','SIDE_TRIM_STRIP'),('side_strip_trim','SIDE_TRIM_STRIP'),('local_notch_trim','LOCAL_PORT_NOTCH_CUT')])
def test_v2_same_wall_strategy_has_explicit_geometric_volume_reason(case,expected):
    from shapely.geometry import box
    from shapely.ops import unary_union
    l,_=load_v2(case);half=V2_CFG['trim_notch_half_width_mm'];depth=6;local=unary_union([box(0,p['trim_origin'][1]-half,depth,p['trim_origin'][1]+half) for p in l['ports']]);strip=box(0,l['low'][1],depth,l['high'][1]);local=local.intersection(strip)
    name,reason=choose_strategy(local.area,strip.area,2,True);assert name==expected;assert '10%' in reason

def test_v2_invalid_shared_strip_is_never_selected():assert choose_strategy(100,100,2,False)[0]=='LOCAL_PORT_NOTCH_CUT'

def test_v2_blade_collision_negative():
    l,_=load_v2('support_blade_collision');g=v2_build(l,V2_CFG,.22);assert fcl_distance(l['core'],shape_mesh(g['frame']))[0]<=0

def test_v2_central_near_wall_is_reported_without_anatomy_change():
    l,_=load_v2('central_tree_near_wall');g=v2_build(l,V2_CFG,.22);assert fcl_distance(l['central'],shape_mesh(g['frame']))[0]<=0

def test_v2_wrong_branch_cut_is_detected():
    l,_=load_v2('trim_plane_hits_wrong_branch');p=l['ports'][0];pt=p['trim_origin'];cut=trimesh.creation.box([.1,8.4,27]);cut.apply_translation([pt[0],pt[1],13.5]);assert (manifold(cut)^manifold(l['central'])).volume()>.01

def test_v2_full_core_inside_walls_independently(v2_meshes,v2_actual):
    l=json.loads((v2_output()/'assembly/layout.json').read_text());c=shape_mesh(bounds_box(l['low'],l['high']));assert (manifold(v2_meshes['core'])-manifold(c)).volume()<1e-7;assert v2_actual['placement']['full_core_inside']

def test_v2_all_four_solid_wall_slabs(v2_meshes):
    l=json.loads((v2_output()/'assembly/layout.json').read_text());lo=np.array(l['low']);hi=np.array(l['high']);shell=bounds_box(lo-[2,2,0],hi+[2,2,0]).cut(bounds_box(lo-[0,0,1],hi+[0,0,1]));assert (manifold(shape_mesh(shell))-manifold(v2_meshes['frame'])).volume()<1e-7

def test_v2_bottom_and_top_open_in_retained_cavity(v2_meshes):
    f=v2_meshes['frame'];assert not f.ray.intersects_any([[0,0,5]],[[0,0,-1]])[0];assert not f.ray.intersects_any([[0,0,5]],[[0,0,1]])[0]

def test_v2_continuous_top_insertion_real_mesh(v2_meshes,v2_actual):
    assert sweep_translation(v2_meshes['core'],[v2_meshes['frame']],[0,0,47])['passed'];assert v2_actual['assembly']['top_insertion']['min_clearance_mm']>.2

def test_v2_every_contact_approved_and_central_clear(v2_actual):
    a=v2_actual['assembly'];assert a['all_design_contacts_on_extensions'] and a['central_contact_count']==0
    for r in a['contacts']:
        assert r['central_clearance_mm']>0
        if r['design_contact']:assert r['approved_extension_zone']

def test_v2_bad_contact_is_not_silently_accepted(v2_actual):
    p=v2_actual['ports'][0];assert not contact_zone(np.array(p['straight_start'])-np.array(p['axis'])*3,p)

def test_v2_cradle_positive_clearance_and_secondary_float(v2_actual):
    assert sum(p['role']=='MASTER' for p in v2_actual['ports'])==1
    for p in v2_actual['ports']:
        assert p['cradle_clearance_mm']>0 and p['anti_float_gap_mm']>0
        if p['role']=='SECONDARY':assert p['axial_float_mm']==pytest.approx(1.)

def test_v2_retainer_top_insertion_then_internal_slide(v2_actual):
    for r in v2_actual['assembly']['retainer_insertions']:assert r['passed'] and r['min_clearance_mm']>.1 and r['sequence']==['TOP_DOWN','SHORT_WALLWARD_SLIDE']

def test_v2_retainer_and_core_anti_float_without_force_claim(v2_actual):
    a=v2_actual['anti_float'];assert a['geometric_Z_constraint'];assert not a['physical_force_validated'];assert all(x['blocked'] for x in a['finite_pose_checks']);assert all(x['up_0p8mm_blocked_by_inner_guides'] for x in a['retainer_lock_checks'])

def test_v2_vertical_perpendicular_trim_avoids_bifurcation(v2_actual):
    for p in v2_actual['ports']:assert p['trim_normal'][2]==0;assert np.dot(p['trim_normal'],p['axis'])==pytest.approx(1);assert p['distance_to_bifurcation_mm']>20

def test_v2_usable_straight_length(v2_actual):
    for p in v2_actual['ports']:assert p['usable_straight_length_after_trim_mm']==pytest.approx(p['straight_length_mm']/2);assert p['distance_to_curvature_mm']==pytest.approx(4.5)

def test_v2_rib_same_thickness_as_blades():
    import csv
    rows=list(csv.DictReader((v2_output()/'tables/wall_support_design.csv').open()))
    for r in rows:assert float(r['rib_thickness_mm'])==float(r['support_blade_thickness_mm']);assert r['rib_inner_face_on_trim_plane']=='True'

def test_v2_ribs_extend_vertically_and_stay_sacrificial():
    r=trimesh.load_mesh(v2_output()/'trim/trim_alignment_ribs.stl');l=json.loads((v2_output()/'assembly/layout.json').read_text());assert r.bounds[1,2]>l['high'][2]-.6
    final=trimesh.load_mesh(v2_output()/'diagnostic/FINAL_PDMS_SIDE.stl');assert (manifold(r)^manifold(final)).volume()<1e-7

def test_v2_final_pdms_has_no_fixture_ghost_channels(v2_meshes,v2_actual):
    final=manifold(trimesh.load_mesh(v2_output()/'diagnostic/FINAL_PDMS_SIDE.stl'));ideal=final-manifold(v2_meshes['core']);actual=manifold(v2_meshes['after']);assert (ideal-actual).volume()<1e-6;assert (actual-ideal).volume()<1e-6;assert v2_actual['trim']['fixture_imprint_remaining_after_trim_mm3']<1e-7

def test_v2_trim_exposes_every_port_and_keeps_central_tree(v2_actual):
    t=v2_actual['trim'];assert t['all_ports_exposed'];assert t['central_vascular_removed_mm3']<1e-7;assert t['minimum_nonport_vascular_trim_clearance_mm']>4;assert t['final_pdms_components']==1

def test_v2_trim_volume_reduces_pdms_and_is_separated(v2_actual):
    t=v2_actual['trim'];assert t['before_pdms_volume_mm3']>t['after_pdms_volume_mm3']>0;assert t['peripheral_to_final_boolean_intersection_mm3']<1e-7

@pytest.mark.parametrize('name',['frame_slice','retainer_slice'])
def test_v2_actual_pla_support_off_archives(v2_actual,name):
    b=v2_actual['bambu'][name];assert b['support_free_slice_passed'];assert b['native_support_moves']==0;assert b['actual_profile']['filament_type']==['PLA'];assert b['actual_profile']['nozzle_diameter']==['0.2'];assert not b['warning_audit']['critical']
    with zipfile.ZipFile(b['archive']) as z:assert z.testzip() is None

def test_v2_two_and_three_lines_real_toolpath_comparison(v2_actual):
    comparison=json.loads((v2_output()/'bambu/blade_line_comparison.json').read_text());assert [x['lines'] for x in comparison]==[2,3];assert v2_actual['supports']['blade_line_count']==2
    for r in comparison:assert r['slice']['slicer_success'];assert r['blade_preservation']['passed'];assert all(x['missing_layers']==0 for x in r['blade_preservation']['rows'])

@pytest.mark.parametrize('key',['SOURCE_UNCHANGED','FRAME_REAL_SLICE','TRIM_IMPRINT_REMOVED','ANTI_FLOAT'])
def test_v2_ready_requires_hard_gates(v2_actual,key):
    s=copy.deepcopy(v2_actual)
    if key=='SOURCE_UNCHANGED':s['source_protection']['all_unchanged']=False
    if key=='FRAME_REAL_SLICE':s['bambu']['frame_slice']['support_free_slice_passed']=False
    if key=='TRIM_IMPRINT_REMOVED':s['trim']['fixture_imprint_remaining_after_trim_mm3']=.01
    if key=='ANTI_FLOAT':s['assembly']['anti_float']['geometric_Z_constraint']=False
    status,_,failed=v2_gates(s);assert status=='NEEDS_ADJUSTMENT' and key in failed

def test_v2_original_frozen_hash_and_ready(v2_actual):
    assert digest(Path(v2_actual['vascular']['source']))==digest(Path(v2_actual['vascular']['independent_frozen_copy']))=='84f1989392fb8d210432750a8022a830943b9bff890b044aa9d310822b619ce0';assert v2_actual['status']=='READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_SEAL_PENDING';assert not v2_actual['failures']

def test_v2_all_twenty_one_qc_and_separate_step():
    from vascular_processing.casting_fixture_v2_review import NAMES
    for n in NAMES:assert (v2_output()/'QC'/(n+'.png')).stat().st_size>20000
    assert (v2_output()/'assembly/assembled_casting_fixture.step').stat().st_size>1000

def test_v2_actual_finite_cut_planes_exclude_corner_exterior(v2_actual):
    t=v2_actual['trim'];assert t['trim_plane_count']==11
    counts={p:sum(r['port']==p for r in t['actual_cut_planes']) for p in ['I1','O1','O2','O3']}
    assert counts=={'I1':3,'O1':3,'O2':2,'O3':3};assert t['maximum_unrelated_port_cut_volume_mm3']<1e-7

def test_v2_no_vertical_undercut_in_retained_pdms(v2_actual):
    t=v2_actual['trim'];assert t['no_obvious_vertical_interlock'];assert t['vertical_release_prism_PLA_intrusion_mm3']<1e-7

def test_v2_real_alternative_clearance_cad():
    import csv
    p=v2_output()/'config/clearance_study';r=list(csv.DictReader((p/'comparison.csv').open()));assert len(r)==12
    assert {round(float(x['clearance_mm']),2) for x in r}=={.11,.22,.33}
    assert all(x['collision_free']=='True' for x in r)
    assert len(list(p.glob('*/*.step')))==12

def test_v2_pdms_before_trim_has_closed_negative_void_shells(v2_meshes,v2_actual):
    m=v2_meshes['before'];assert m.is_watertight
    shells=m.split(only_watertight=False);assert sum(q.volume>0 for q in shells)==1
    assert sum(q.volume<0 for q in shells)==5  # ABS core plus four separate retainers.
    assert v2_actual['trim']['diagnostic_mesh_qc_passed']

def test_v2_master_stop_exception_is_only_distal_face(v2_actual):
    p=next(x for x in v2_actual['ports'] if x['role']=='MASTER')
    assert contact_zone(p['end'],p,distal_stop=True)
    assert not contact_zone(p['straight_start'],p,distal_stop=True)
