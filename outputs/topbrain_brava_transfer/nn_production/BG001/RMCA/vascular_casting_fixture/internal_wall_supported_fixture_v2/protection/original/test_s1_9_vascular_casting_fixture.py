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
