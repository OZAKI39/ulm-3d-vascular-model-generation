"""Permanent all-ABS regression checks, including actual BG001 and bad cases."""
import json,zipfile
from pathlib import Path
import numpy as np
import pytest,trimesh
from vascular_processing.support_inputs import ROOT,discover,digest,verify
from vascular_processing.all_abs_pipeline import default_output,gates
from vascular_processing.all_abs_synthetic import run_case,vapor_scope,NAMES
from vascular_processing.all_abs_audit import network_access,graph_status
from vascular_processing.casting_fixture_geometry import rigid

@pytest.fixture(scope='module')
def out():return default_output()
@pytest.fixture(scope='module')
def s(out):return json.loads((out/'all_abs_sacrificial_casting_summary.json').read_text())
@pytest.fixture(scope='module')
def sources():return discover()


def test_source_discovery(sources):
    assert sources['candidate'].name=='candidate_0008'
    assert all(sources[x].is_file() for x in ['core','compensated_swc','smoothed_swc','geometry'])

def test_frozen_copy_bytes(out,sources):assert digest(out/'vascular_core/BG001_RMCA_frozen_ABS_core.stl')==digest(sources['core'])
def test_source_hash_freeze(out):
    x=json.loads((out/'source_protection.json').read_text());assert x['mismatch_count']==0;assert verify(x['start_hashes'])['all_unchanged'];assert x['start_hashes']==x['end_hashes']

def test_vascular_only_rigid(out,sources):
    m=np.asarray(json.loads((out/'assembly/casting_core_transform.json').read_text())['transform_4x4']);rigid(m);a=trimesh.load_mesh(sources['core']);a.apply_transform(m);b=trimesh.load_mesh(out/'vascular_core/assembly_placed_ABS_core.stl');assert abs(a.volume-b.volume)<.0001;assert a.bounds==pytest.approx(b.bounds,abs=1e-5)

def test_scale_is_one(s):assert s['vascular']['scale_factor']==1
@pytest.mark.parametrize('key',['walls_closed','bottom_closed','top_open'])
def test_closed_mold_open_top(s,key):assert s['outer_mold'][key]
def test_no_wall_openings(s):assert s['assembly']['wall_missing_volume_mm3']==0

def test_complete_core_inside(s):assert abs(s['trim']['core_outside_cavity_mm3'])<1e-8

def test_top_insertion_is_continuous(s):
    a=s['assembly']['top_insertion'];assert a['passed'];assert a['positive_clearance_without_lands']['distance_queries']>10;assert a['full_mold_to_0p02mm_standoff']['passed'];assert a['downward_face_condition'];assert a['final_interference_mm3']<1e-6

def test_contact_extension_only(s):assert s['assembly']['all_contacts_on_approved_extensions'];assert s['assembly']['fixture_contacts']==4

def test_central_no_auxiliary_contact(s):assert s['assembly']['central_tree_contacts']==0;assert s['assembly']['minimum_central_auxiliary_clearance_mm']>5

def test_retainer_continuous_insertion(s):assert all(x['passed'] and x['min_clearance_mm']>.1 for x in s['assembly']['retainer_insertions'])
def test_geometric_antifloat(s):assert s['assembly']['anti_float_geometric_pass'];assert len(s['assembly']['anti_float_poses'])==5;assert all(x['blocked'] for x in s['assembly']['retainer_locks'])

def test_mold_manifold_actual_file(out):
    m=trimesh.load_mesh(out/'outer_mold/ABS_sacrificial_outer_mold.stl');assert m.is_watertight and m.is_winding_consistent;assert len(m.split())==1
@pytest.mark.parametrize('k',['I1','O1','O2','O3'])
def test_retainer_manifold_actual_file(out,k):
    m=trimesh.load_mesh(out/f'retainers/port_{k}_ABS_retainer.stl');assert m.is_watertight and m.is_winding_consistent;assert len(m.split())==1
@pytest.mark.parametrize('key',['mold_auto_support_off_slice','retainer_auto_support_off_slice'])
def test_actual_abs_support_off(s,key):
    r=s['bambu'][key];assert r['slice_passed'] and r['auto_support_off'];assert r['native_support_extrusion_segments']==0;assert r['actual_profile']['filament_type']==['ABS'];assert r['actual_profile']['nozzle_diameter']==['0.2'];assert r['archive_geometry']['passed']
    with zipfile.ZipFile(r['archive']) as z:assert z.testzip() is None

def test_actual_core_reference(s):assert s['bambu']['vascular_reference_slice']['slice_passed'];assert s['bambu']['vascular_reference_slice']['actual_profile']['filament_type']==['ABS']
def test_vapor_scope_excludes_core(s):assert s['smoothing']['vascular_core_forbidden'];assert not s['smoothing']['dimensionally_calibrated'];assert not s['smoothing']['software_compensation_applied']

def test_access_graph_actual_edges(out):
    import networkx as nx
    r=json.loads((out/'solvent_access/acetone_access_graph.json').read_text());g=nx.DiGraph();g.add_edges_from((x['source'],x['target']) for x in r['edges']);assert graph_status(g,['port_'+k for k in ['I1','O1','O2','O3']])['status']=='SOLVENT_ACCESS_GEOMETRIC_PASS';assert len(list(nx.node_disjoint_paths(g,'external_bath','vascular_component_0')))==4

def test_access_edges_have_real_contact(s):
    for r in s['solvent_access']['port_paths']:
        assert r['contact_area_mm2']>.01 and r['contact_surface_error_mm']<1e-5 and r['interference_mm3']<1e-6
        assert r['wall_fixture_attachment_mm3']>0 and r['fixture_land_attachment_mm3']>0 and r['retainer_exposed_volume_mm3']>0

def test_all_ports_and_retainers_reachable(s):assert s['solvent_access']['all_ports_accessible'] and s['solvent_access']['all_retainers_accessible']
def test_all_real_swc_branches_reachable(s):
    for key in ['network','compensated_network']:
        n=s['solvent_access'][key];assert n['all_vascular_branches_reachable'];assert n['branch_count']==5;assert all(len(r['reachable_ports'])==4 for r in n['rows'])

def test_aux_void_calculation(s):assert s['voids']['auxiliary_void_volume_before_trim']>0;assert abs(s['voids']['auxiliary_void_volume_after_trim'])<1e-6

def test_trim_hits_extensions(out):
    import csv
    rows=list(csv.DictReader((out/'trim/port_trim_sections.csv').open()));assert len(rows)==4;assert all(float(r['section_area_mm2'])>1. and float(r['axis_plane_angle_degrees'])==90. for r in rows)

def test_trim_spares_central_and_unrelated(s):assert s['trim']['central_vascular_removed_mm3']==0;assert s['trim']['maximum_unrelated_port_cut_mm3']==0;assert s['trim']['minimum_nonport_clearance']>0

def test_aux_void_removed_without_ghosts(s):assert abs(s['trim']['auxiliary_void_volume_after_trim'])<1e-6;assert abs(s['trim']['ghost_missing_volume_mm3'])<1e-6;assert abs(s['trim']['ghost_extra_volume_mm3'])<1e-6

def test_all_ports_exposed(s):assert s['trim']['all_ports_exposed']
def test_usable_port_lengths(s):assert s['trim']['usable_port_lengths']==dict(I1=4.5,O1=4.5,O2=4.5,O3=4.5)
@pytest.mark.parametrize('name',['PDMS_as_cast','PDMS_after_all_ABS_removed','PDMS_after_trim'])
def test_pdms_exports_watertight(out,name):
    m=trimesh.load_mesh(out/f'diagnostic_PDMS/{name}.stl');assert m.is_watertight and m.is_winding_consistent

def test_exports_step_and_independent_assembly(out):
    for p in [out/'outer_mold/ABS_sacrificial_outer_mold.step',out/'assembly/all_ABS_casting_system.step',*[out/f'retainers/port_{k}_ABS_retainer.step' for k in ['I1','O1','O2','O3']]]:
        assert p.is_file() and p.stat().st_size>1000;assert 'ISO-10303-21' in p.open().read(40)
    names=set()
    with (out/'assembly/all_ABS_casting_system.step').open() as f:
        for line in f:
            if 'PRODUCT(' in line:
                for name in ['FROZEN_ABS_CORE','ABS_SACRIFICIAL_MOLD','ABS_RETAINER_I1','ABS_RETAINER_O1','ABS_RETAINER_O2','ABS_RETAINER_O3']:
                    if name in line:names.add(name)
    assert len(names)==6

def test_no_printer_delivery(s):
    assert s['no_printer_delivery']
    for k in ['mold_auto_support_off_slice','retainer_auto_support_off_slice','vascular_reference_slice']:
        r=s['bambu'][k];assert not r['printer_upload'] and not r['printer_start'];assert '--slice' in r['command'];assert not any(x in ' '.join(r['command']).lower() for x in ['--send','--upload','--start-print'])

def test_s1_9_files_protected(out):
    h=json.loads((out/'source_protection_start.json').read_text());assert any('s1-9_' in x for x in h);assert any('internal_wall_supported_fixture_report.md' in x for x in h);assert verify(h)['all_unchanged']

def test_thin_features_in_actual_toolpaths(s):
    a=s['bambu']['feature_preservation'];assert a['passed'];assert {r['feature'] for r in a['rows']}>=set(['solvent_web','conformal_access_land','top_exposed_retainer_stem','transferred_trim_guide'])

def test_24_review_figures(out):
    from vascular_processing.all_abs_review import NAMES
    assert len(NAMES)==24
    for n in NAMES:assert (out/'QC'/(n+'.png')).stat().st_size>10000

def test_state_requires_every_gate(s):
    status,checks,failed=gates(s);assert status=='READY_FOR_HUMAN_REVIEW_WITH_PROCESS_PARAMETERS_UNCALIBRATED';assert not failed;assert not any(s['physical_validation'].values())

def test_missing_route_blocks_readiness(s):
    import copy
    bad=copy.deepcopy(s);bad['solvent_access']['all_ports_accessible']=False;assert gates(bad)[0]=='NEEDS_ADJUSTMENT'

def test_nonrigid_transform_rejected():
    m=np.eye(4);m[0,0]=1.01
    with pytest.raises(ValueError):rigid(m)
@pytest.mark.parametrize('name',NAMES)
def test_synthetic_geometric_cases(name,tmp_path):
    r=run_case(name,tmp_path/name)
    if name=='blocked_port_access':assert r['solvent_status']=='SOLVENT_ACCESS_BLOCKED';assert r['core_gap_mm']>.49
    elif name=='fixture_touches_tree_fail':assert r['central_tree_contact'];assert r['solvent_status']=='SOLVENT_ACCESS_BLOCKED'
    elif name=='trim_hits_branch_fail':assert r['central_cut_volume_mm3']>0
    else:assert r['solvent_status']=='SOLVENT_ACCESS_GEOMETRIC_PASS';assert not r['central_tree_contact'];assert r['central_cut_volume_mm3']==0
    assert r['disconnected_network_rejected'];assert r['core_smoothing_rejected']

def test_fully_encapsulated_core_no_access(tmp_path):
    r=run_case('blocked_port_access',tmp_path);assert r['solvent_status']=='SOLVENT_ACCESS_BLOCKED';assert not r['branch_network']['all_vascular_branches_reachable']

def test_smoothing_policy():assert not vapor_scope(['INNER_BOTTOM','VASCULAR_CORE'])['passed'];assert vapor_scope(['INNER_SIDE_WALLS'])['passed']

def test_reproducible_synthetic_inputs():
    for name in NAMES:
        folder=ROOT/'tests/data/all_abs_casting'/name;assert (folder/'input.json').is_file();assert (folder/'core.stl').is_file();assert (folder/'result.json').is_file()

def test_stale_contact_metadata_cannot_bridge_a_real_gap():
    from vascular_processing.all_abs_audit import contact_surface_evidence
    from vascular_processing.casting_fixture_geometry import bounds_box,shape_mesh
    core=shape_mesh(bounds_box([-1,-1,1],[1,1,2]));land=shape_mesh(bounds_box([-.5,-.5,0],[.5,.5,1]));w=[[0,0,1]];cfg=dict(contact_area_tolerance_mm2=.01,contact_distance_tolerance_mm=1e-5)
    assert contact_surface_evidence(core,land,w,1.,cfg)['passed']
    land.apply_translation([0,0,-.2]);assert not contact_surface_evidence(core,land,w,1.,cfg)['passed']

@pytest.mark.parametrize('name',['all_ABS_voids_inside_cavity','auxiliary_voids_before_trim'])
def test_void_exports_have_no_zero_thickness_wall_sheets(out,name):
    m=trimesh.load_mesh(out/f'diagnostic_PDMS/{name}.stl')
    for c in m.split(only_watertight=False):assert c.is_watertight;assert np.min(c.extents)>1e-7

def test_diagnostic_regularization_does_not_modify_production(s):
    t=s['trim'];assert not t['production_geometry_changed_for_sewing'];assert t['diagnostic_contact_added_outside_occupied_domain_mm3']<1e-6;assert t['diagnostic_sewing_PDMS_symmetric_difference_mm3']<t['diagnostic_PDMS_difference_limit_mm3'];assert not t['frozen_core_snapped']

# COMPACT_C_RETAINER_V2 regression extension. All V1 tests above stay permanent.
@pytest.fixture(scope='module')
def compact_out():
    from vascular_processing.compact_c_pipeline import default_output
    return default_output()

@pytest.fixture(scope='module')
def cs(compact_out):return json.loads((compact_out/'compact_c_retainer_summary.json').read_text())

@pytest.fixture(scope='module')
def compact_geometry(cs):
    import yaml
    from vascular_processing.compact_c_geometry import load_frozen,candidate
    b,s,i,c,l,g,line=load_frozen();cfg=yaml.safe_load((ROOT/'config/compact_c_retainer_v2.yaml').read_text());return l,g,candidate(l,g,cfg,line,cs['selected_wrap_angle_deg'])


def test_c_retainer_not_closed_ring(cs):
    from vascular_processing.compact_c_geometry import validate_parameters
    assert all(0<p['c_retainer_wrap_angle_deg']<=180 and p['retainer_opening_angle_deg']>=180 for p in cs['ports'])
    for angle in [0,181,359,360]:
        with pytest.raises(ValueError,match='CLOSED_OR_EXCESSIVE_WRAP'):validate_parameters(angle,.22,.44,1.2)


def test_c_retainer_positive_clearance(cs,compact_geometry):
    from vascular_processing.casting_mold_qc import fcl_distance
    from vascular_processing.casting_fixture_geometry import shape_mesh
    l,old,new=compact_geometry
    for p in cs['ports']:assert p['upper_clearance_mm']==.22 and p['guide_clearance_mm']==.12
    assert all(fcl_distance(l['core'],shape_mesh(m))[0]>.21 for m in new['retainers'].values())


def test_c_retainer_top_insertion(cs):
    for p in cs['insertion']['rows']:
        assert p['passed'] and p['positive_clearance_path']['distance_queries']>10
        assert p['final_seating_path_without_walls']['passed'] and p['wall_seating_halfspace_proof']
        assert p['mold_interference_mm3']<1e-6 and p['slide_mm']==1.6


def test_c_retainer_no_snap_required(cs):assert cs['insertion']['no_elastic_deformation'];assert all(p['no_elastic_snap_required'] for p in cs['ports'])

def test_c_retainer_extension_only_contact(cs):
    assert cs['contact']['extension_only_contacts'];assert all(r['approved_extension_zone'] and r['distance_to_transition_point_mm']>0 and r['distance_to_nearest_bifurcation_mm']>0 for r in cs['contact']['rows'])


def test_c_retainer_no_tree_contact(cs,compact_geometry):
    from vascular_processing.casting_mold_qc import fcl_distance
    from vascular_processing.casting_fixture_geometry import shape_mesh
    l,old,new=compact_geometry
    assert cs['contact']['central_tree_contacts']==0
    assert all(fcl_distance(l['central'],shape_mesh(v))[0]>10 for v in new['retainers'].values())


def test_c_retainer_master_constraint(cs):
    master=[r for r in cs['axial_float']['rows'] if r['role']=='MASTER'];assert len(master)==1;assert master[0]['port_id']=='O2';assert master[0]['required_axial_float_mm']==0
    assert 0<master[0]['new_max_mm']<.23;assert master[0]['new_max_mm']==master[0]['old_max_mm']


def test_c_retainer_secondary_axial_float(cs):
    rows=[r for r in cs['axial_float']['rows'] if r['role']=='SECONDARY'];assert len(rows)==3
    for r in rows:
        assert r['new_min_mm']<=-1 and r['new_max_mm']>=1
        assert r['new_min_mm']==r['old_min_mm'] and r['new_max_mm']==r['old_max_mm']
        assert r['required_float_continuous_paths']['new']['passed'] and r['land_stays_on_straight_at_required_float']
        assert not any(r['new_search_censored'])


def test_c_retainer_Z_stop(cs):
    assert cs['anti_float']['poses'][0]['blocked'];assert all(p['port_up_1mm_blocked'] and p['retainer_up_0p8_blocked'] for p in cs['anti_float']['ports'])


def test_c_retainer_pitch_roll(cs):
    assert {p['mode'] for p in cs['anti_float']['poses'][1:]}=={'pitch+','pitch-','roll+','roll-'}
    assert all(p['blocked'] for p in cs['anti_float']['poses'][1:]);assert cs['anti_float_not_worse']


def test_c_retainer_height_reduction(cs):assert cs['compact']['max_height']<3;assert cs['compact']['height_reduction_percent']>90

def test_c_retainer_volume_reduction(cs):assert cs['compact']['volume_reduction_percent']>80;assert all(r['new_volume']<r['old_volume'] for r in cs['comparison'])


def test_wall_solvent_route_without_retainer_stem(compact_out,cs):
    import networkx as nx
    d=json.loads((compact_out/'solvent_access/acetone_access_graph.json').read_text());G=nx.DiGraph();G.add_edges_from((e['source'],e['target']) for e in d['edges']);G.remove_nodes_from([n for n in G if n.startswith('retainer_')])
    for k in ['I1','O1','O2','O3']:
        path=nx.shortest_path(G,'external_bath','port_'+k);assert path==['external_bath','wall_'+k,'fixture_'+k,'access_land_'+k,'port_'+k]
    assert len(list(nx.node_disjoint_paths(G,'external_bath','vascular_component_0')))==4;assert not cs['solvent']['retainer_stem_required']


def test_all_ports_wall_solvent_access(cs):assert all(r['wall_based_solvent_access']=='PASS' and not r['retainer_route_required'] for r in cs['solvent']['wall_based_routes'])


def test_SWC_reachable_without_retainer_stem(cs):
    for key in ['network','compensated_network']:
        n=cs['solvent'][key];assert n['branch_count']==5 and n['all_vascular_branches_reachable'];assert all(len(r['reachable_ports'])==4 for r in n['rows'])


def test_retainer_void_not_increased(cs):assert all(r['new_retainer_auxiliary_void_mm3']<r['old_retainer_auxiliary_void_mm3'] for r in cs['retainer_voids']['rows']);assert cs['PDMS']['auxiliary_void_new']<cs['PDMS']['auxiliary_void_old']

def test_retainer_void_trimmed_out(cs):assert abs(cs['PDMS']['auxiliary_void_after_trim'])<1e-6;assert all(abs(r['new_retainer_void_outside_existing_trim_mm3'])<1e-6 for r in cs['retainer_voids']['rows'])


def test_retainer_auto_support_off(cs):
    r=cs['bambu']['selected_slice'];assert r['slice_passed'] and r['auto_support_off'] and r['native_support_extrusion_segments']==0;assert r['actual_profile']['filament_type']==['ABS'];assert r['actual_profile']['nozzle_diameter']==['0.2'];assert r['archive_geometry']['passed']
    with zipfile.ZipFile(r['archive']) as z:assert z.testzip() is None


def test_retainer_feature_preservation(cs):
    from vascular_processing.bambu_casting_mold import parse_paths
    f=cs['bambu']['features'];assert f['passed'];assert len(f['rows'])==12;assert all(r['query_count']>0 and r['missing_query_count']==0 for r in f['rows']);assert all(r['passed'] and r['estimated_bead_envelope_radial_clearance_mm']>.17 for r in f['axial_float_slots'])
    parsed=parse_paths(cs['bambu']['selected_slice']['archive']);assert (parsed['kind']==0).sum()==f['model_segment_count'];assert (parsed['kind']==1).sum()==0


def test_buried_C_retainer_self_removal_has_finite_wall_contact(cs):
    assert all(p['below_nominal_fill'] for p in cs['ports']);assert cs['solvent']['all_buried_retainers_accessible']
    assert all(p['wall_seated_contact'] and p['finite_contact_area_mm2']>1. and p['no_stem'] for p in cs['solvent']['retainer_removal_routes'])


def test_compact_protection_includes_legacy_and_frozen_mold(compact_out):
    d=json.loads((compact_out/'source_protection.json').read_text());assert d['mismatch_count']==0 and d['start_hashes']==d['end_hashes'];assert verify(d['start_hashes'])['all_unchanged'];assert any('legacy_retainer_v1' in k for k in d['start_hashes']);assert any(k.endswith('outer_mold/ABS_sacrificial_outer_mold.stl') for k in d['start_hashes'])


def test_compact_trim_planes_and_lengths_are_exactly_unchanged(cs,s):
    assert cs['PDMS']['trim']['actual_cut_planes']==s['trim']['actual_cut_planes'];assert cs['PDMS']['trim']['usable_port_lengths']==s['trim']['usable_port_lengths'];assert cs['PDMS']['trim']['all_ports_exposed']


@pytest.mark.parametrize('name',['PDMS_after_C_retainer_dissolution','PDMS_after_C_retainer_trim'])
def test_compact_PDMS_export_actual_mesh(compact_out,name):
    m=trimesh.load_mesh(compact_out/f'diagnostic_PDMS/{name}.stl');assert m.is_watertight and m.is_winding_consistent;assert m.volume>0


@pytest.mark.parametrize('k',['I1','O1','O2','O3'])
def test_compact_retainer_actual_STL_STEP(compact_out,k):
    m=trimesh.load_mesh(compact_out/f'retainers/port_{k}_C_retainer.stl');assert m.is_watertight and m.is_winding_consistent and len(m.split())==1;assert (compact_out/f'retainers/port_{k}_C_retainer.step').stat().st_size>1000


def test_compact_independent_assembly_products(compact_out):
    path=compact_out/'assembly/all_ABS_casting_system_C_retainer.step';names=set()
    with path.open() as f:
        for line in f:
            if 'PRODUCT(' in line:
                for n in ['FROZEN_ABS_CORE','ABS_SACRIFICIAL_MOLD','ABS_RETAINER_I1','ABS_RETAINER_O1','ABS_RETAINER_O2','ABS_RETAINER_O3']:
                    if n in line:names.add(n)
    assert len(names)==6


def test_C_retainer_lexicographic_selection_ignores_CAD_roundoff(cs):
    from vascular_processing.compact_c_pipeline import selection_key
    winner=min(cs['candidates'],key=selection_key);assert winner['wrap_angle_deg']==cs['selected_wrap_angle_deg']==90
    assert len(cs['candidates'])==3 and all(c['hard_gates_passed'] for c in cs['candidates'])
    assert selection_key(dict(hard_gates_passed=False,max_height_above_port_mm=0,total_volume_mm3=0))>selection_key(winner)


def test_C_readiness_rejects_missing_evidence(cs):
    import copy
    from vascular_processing.compact_c_pipeline import gates
    assert gates(cs)[0]=='READY_FOR_HUMAN_REVIEW_WITH_PROCESS_PARAMETERS_UNCALIBRATED'
    for section,key in [('solvent','all_ports_wall_access'),('solvent','all_buried_retainers_accessible'),('bambu','critical_features_preserved'),('axial_float','passed')]:
        bad=copy.deepcopy(cs);bad[section][key]=False;assert gates(bad)[0]=='NEEDS_ADJUSTMENT'


def test_C_retainer_review_18_figures(compact_out):
    from vascular_processing.compact_c_review import NAMES
    assert len(NAMES)==18
    for n in NAMES:assert (compact_out/'QC'/(n+'.png')).stat().st_size>10000
    cam=json.loads((compact_out/'QC/legacy_compact_shared_camera.json').read_text());assert cam['parallel_scale']>0


def test_C_no_vapor_no_physical_claim_no_printer_delivery(cs):
    assert cs['smoothing']['retainers_forbidden'] and cs['smoothing']['vascular_core_forbidden'];assert not cs['smoothing']['dimensionally_calibrated'];assert not any(cs['physical_validation'].values());assert cs['no_printer_delivery'];r=cs['bambu']['selected_slice'];assert not r['printer_upload'] and not r['printer_start']


from vascular_processing.compact_c_synthetic import NAMES as C_SYNTHETICS
@pytest.mark.parametrize('name',C_SYNTHETICS)
def test_c_retainer_synthetic_geometric_controls(name,tmp_path):
    from vascular_processing.compact_c_synthetic import run_case
    r=run_case(name,tmp_path/name)
    assert r['parameter_rejected']==(name in ['too_tight_C','closed_ring_invalid'])
    assert r['central_tree_contact']==(name=='C_hits_tree')
    assert r['insertion']['passed']==(name!='C_insertion_blocked')
    assert r['wall_route_pass']==(name!='wall_solvent_route_blocked')
    assert r['master_positive_stop']==(name=='master_port')
    if name in ['single_port','floating_secondary']:assert r['secondary_float_pass'] and r['port_up_1mm_blocked'] and r['positive_core_gap_mm']>.21
    if name=='wall_solvent_route_blocked':assert not r['contact_evidence']['passed']
    folder=ROOT/'tests/data/all_abs_casting/c_retainer'/name
    assert (folder/'input.json').is_file() and (folder/'result.json').is_file() and (folder/'C_retainer.stl').is_file()
