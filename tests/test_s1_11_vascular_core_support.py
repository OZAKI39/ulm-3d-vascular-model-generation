"""Permanent s1-11 regression: real data, synthetic collisions, native G-code."""
import copy,json,pickle,zipfile
from pathlib import Path
import numpy as np
import pytest,trimesh
from scipy.spatial import cKDTree
from shapely.geometry import box
from vascular_processing.core_support_inputs import default_output,discover_source,rotations,ROOT
from vascular_processing.core_support_geometry import envelope,route,overlap,intervals,ruled
from vascular_processing.core_support_fixtures import fixture,NAMES
from vascular_processing.core_support_slicing import actual_underbody,actual_self,write_project
from vascular_processing.core_support_removal import wrap_audit,notch_access
from vascular_processing.core_support_pipeline import hard_gates
from vascular_processing.support_inputs import digest,verify
from vascular_processing.support_rfus_geometry import lower_points
from vascular_processing.support_rfus_fixtures import tube
from vascular_processing.support_removal import swept_path


@pytest.fixture(scope='module')
def context():return pickle.load((default_output()/'source/context.pkl').open('rb'))


@pytest.fixture(scope='module')
def real():
    p=default_output()/'candidates/O001_R25_C2'
    return p,pickle.load((p/'geometry.pkl').open('rb'))


def test_source_discovery_reads_current_s1_10_only():
    source=discover_source()
    assert 'all_abs_sacrificial_casting_system' in str(source['core'])
    assert 'vascular_support_planner' not in str(source['core'])
    assert not any(k in source for k in ['box','mold','retainers','pdms'])
    assert source['compensated_swc'].is_file() and source['smoothed_swc'].is_file()


def test_all_s1_10_hashes_frozen():
    h=json.loads((default_output()/'validation/protected_before.json').read_text())
    assert verify(h)['all_unchanged']
    assert any('outer_mold' in p for p in h) and any('retainers' in p for p in h)
    assert any('compensated.swc' in p for p in h)


def test_current_saved_profiles_are_not_modified():
    r=json.loads((default_output()/'bambu/profiles/profiles.json').read_text())
    assert r['values']['source']=='CURRENT_SAVED_VESSEL_PROFILE'
    assert r['actual']['process']['enable_support']=='0'
    assert r['actual']['machine']['nozzle_diameter']==['0.2']
    for source in r['values']['profile_sources']:
        assert digest(source['path'])==source['sha256']


def test_hash_mutation_detected(tmp_path):
    p=tmp_path/'source.stl';p.write_bytes(b'frozen');h={str(p):digest(p)};p.write_bytes(b'changed')
    with pytest.raises(RuntimeError,match='UPSTREAM_GEOMETRY_MODIFIED'):verify(h)


def test_orientation_deterministic_rigid_and_three_dimensional():
    a=rotations();b=rotations();assert len(a)==48
    for x,y in zip(a,b):
        assert np.array_equal(x,y)
        assert x.T@x==pytest.approx(np.eye(3),abs=1e-10)
        assert np.linalg.det(x)==pytest.approx(1.)
    assert len({tuple(np.round(r[2],5)) for r in a})>=20


def test_frozen_vascular_not_remeshed_or_scaled(context,real):
    p,g=real;M=np.asarray(json.loads((p/'vascular_print_transform.json').read_text())['transform_4x4']);inv=np.asarray(json.loads((p/'inverse_transform.json').read_text())['transform_4x4']);assert M@inv==pytest.approx(np.eye(4),abs=1e-10)
    source=trimesh.load_mesh(context['inputs']['core']);source.apply_transform(M);actual=trimesh.load_mesh(p/'vascular_core.stl')
    assert len(source.faces)==len(actual.faces)
    assert cKDTree(source.vertices).query(actual.vertices)[0].max()<2e-4
    assert cKDTree(actual.vertices).query(source.vertices)[0].max()<2e-4


@pytest.mark.parametrize('name',NAMES)
def test_synthetic_frozen_stl_ray_envelope(name,context):
    mesh,bs=fixture(name);env=envelope(mesh,bs,context['profile'],.2)
    assert mesh.is_volume
    for b in env:
        assert b['valid'].all()
        assert np.diff(b['arc']).max()<=context['profile']['line_width']
        assert np.max(b['lower'][:,2]-b['xyz'][:,2])<1e-4
    assert (ROOT/'tests/data/vascular_core_support'/(name+'.stl')).exists()


def test_horizontal_ray_is_real_mesh_surface_not_radius_guess(context):
    mesh,bs=fixture('horizontal_tube');q=lower_points(mesh,bs[0]['xyz'])
    assert q[:,2]==pytest.approx(9.4,abs=2e-5)
    altered=copy.deepcopy(bs);altered[0]['radius']*=4
    other=envelope(mesh,altered,context['profile'],.2)[0]['lower']
    assert other[:,2]==pytest.approx(9.4,abs=2e-5)


@pytest.mark.parametrize('name',['horizontal_tube','curved_tube','S_tube','long_horizontal_port'])
def test_continuous_roof_and_vertical_curtain(name,context):
    mesh,bs=fixture(name);b=envelope(mesh,bs,context['profile'],.2)[0]
    s=route(b['xyz'][:,:2],b['roof_top'],mesh,context['profile'],2,2)
    assert s['route']!='UNRESOLVED' and s['mesh'].is_volume
    assert overlap(s['mesh'],mesh)<1e-7
    roof=lower_points(s['roofs'],b['lower']+[0,0,.01])
    assert np.isfinite(roof).all()
    assert np.max(np.abs(roof[:,2]-b['roof_top']))<.002
    assert s['mesh'].bounds[0,2]==pytest.approx(0.)


@pytest.mark.parametrize('name',['crossing_tubes','stacked_tubes'])
def test_projection_collision_repaired_without_deleting_upper_roof(name,context):
    mesh,bs=fixture(name);b=envelope(mesh,bs,context['profile'],.2)[0];p={**context['profile'],'roof_gap_mm':.2}
    direct=ruled(b['xyz'][:,:2],b['roof_top'],.44,.44,.1)
    assert overlap(direct['mesh'],mesh)>0
    s=route(b['xyz'][:,:2],b['roof_top'],mesh,p,2,2)
    assert s['route']!='UNRESOLVED'
    assert overlap(s['mesh'],mesh)<1e-7
    hit=lower_points(s['roofs'],b['lower']+[0,0,.01]);assert np.isfinite(hit).all()


def test_single_lateral_offset_is_tried_before_twins(context):
    from vascular_processing.core_support_geometry import union
    upper,bs=fixture('horizontal_tube');lower=tube([[100,110,5],[130,110,5]],radius=.1);mesh=union([upper,lower]);b=envelope(mesh,bs,context['profile'],.2)[0]
    s=route(b['xyz'][:,:2],b['roof_top'],mesh,context['profile'],2,2)
    assert s['route']=='LATERAL_OFFSET' and abs(s['offset_mm'])<=2*context['profile']['support_line_width']
    assert overlap(s['mesh'],mesh)<1e-7


def test_adaptive_segmentation_and_gap(real):
    _,g=real
    for b in g['branches']:
        intervals_=intervals(b)
        assert intervals_[0][0]==0 and intervals_[-1][1]==len(b['arc'])-1
        assert all(x[1]==y[0] for x,y in zip(intervals_[:-1],intervals_[1:]))
        assert np.max(b['lower'][:,2]-b['roof_top'])==pytest.approx(.25+.05+.00001)


def test_actual_30mm_horizontal_underbody_is_100_percent():
    p=default_output()/'validation/horizontal_actual_result.json';r=json.loads(p.read_text())
    assert r['slice']['slicer_success'] and r['auto_support_off']
    assert r['coverage']['centerline_underbody_percent']==100
    assert r['coverage']['max_uncovered_arc_mm']==0
    assert r['support_self']['self_unsupported_regions']==0
    assert not r['layer_witness']['credited_as_underbody']


def test_gcode_missing_middle_arc_is_never_credited(context):
    p=default_output()/'validation/actual_horizontal_tube';g=pickle.load((p/'geometry.pkl').open('rb'));fp=pickle.load((p/'actual_footprints.pkl').open('rb'))['support'];sparse={z:poly.difference(box(105,100,125,120)) for z,poly in fp.items()}
    result=actual_underbody(g['branches'],sparse,g['profile'])
    assert not result['passed'] and result['max_uncovered_arc_mm']>19


def test_lower_wall_layer_cannot_hide_a_roof_that_is_too_close(context):
    mesh,bs=fixture('horizontal_tube');env=envelope(mesh,bs,context['profile'],.2)
    footprint=box(99,109,131,111)
    good=actual_underbody(env,{9.12:footprint},context['profile'])
    bad=actual_underbody(env,{9.12:footprint,9.32:footprint},context['profile'])
    assert good['passed'] and not bad['passed']
    assert all(r['too_close_support_layer_tops']==[9.32] for r in bad['rows'])


def test_support_self_cannot_borrow_vascular_foundation(context):
    result=actual_self({2.:box(0,0,1,1)},context['profile'])
    assert not result['passed'] and result['self_unsupported_regions']>0


def test_bed_roots_distinguish_cavity_from_floating_material():
    from vascular_processing.core_support_geometry import bed_roots
    body=trimesh.creation.box([4,4,4]);body.apply_translation([0,0,2])
    cavity=trimesh.creation.box([1,1,1]);cavity.apply_translation([0,0,2]);cavity.invert()
    rows=bed_roots(trimesh.util.concatenate([body,cavity]))
    assert len(rows)==1 and rows[0]['bed_rooted']
    island=trimesh.creation.box([1,1,1]);island.apply_translation([10,0,2])
    rows=bed_roots(trimesh.util.concatenate([body,cavity,island]))
    assert len(rows)==2 and sum(not r['bed_rooted'] for r in rows)==1


def test_wrap_trap_rejected(context):
    mesh,bs=fixture('horizontal_tube');env=envelope(mesh,bs,context['profile'],.2)
    ring=trimesh.creation.annulus(r_min=.8,r_max=1.,height=.6,sections=64)
    M=trimesh.geometry.align_vectors([0,0,1],[1,0,0]);M[:3,3]=[115,110,10];ring.apply_transform(M)
    r=wrap_audit(mesh,ring,env,context['profile']);assert not r['passed'] and r['maximum_wrap_angle_deg']>120


def test_cut_lip_accessible_from_open_side(context):
    mesh,bs=fixture('horizontal_tube');b=envelope(mesh,bs,context['profile'],.2)[0];s=ruled(b['xyz'][:,:2],b['roof_top'],.44,.44,.1)
    n=dict(notch_id='TEST',point=[115,110,5],lateral=[0,1,0],tangent=[1,0,0],cut_width_mm=.44)
    assert notch_access(n,mesh,s['mesh'])['passed']


def test_simple_downward_segment_removal(context):
    mesh,bs=fixture('horizontal_tube');b=envelope(mesh,bs,context['profile'],.2)[0];s=ruled(b['xyz'][20:40,:2],b['roof_top'][20:40],.44,.44,.1)
    start=np.eye(4);end=np.eye(4);end[2,3]=-15
    assert swept_path(s['mesh'],[mesh],[start,end])['passed']


def test_native_3mf_parts_are_separate_auto_support_off(real):
    p,_=real
    with zipfile.ZipFile(p/'bambu/sliced.3mf') as z:
        import xml.etree.ElementTree as ET
        settings=json.loads(z.read('Metadata/project_settings.config'));cfg=ET.fromstring(z.read('Metadata/model_settings.config'))
        assert z.testzip() is None and settings['enable_support']=='0'
        assert len(cfg.findall('.//part'))==2
        assert all(x.attrib['subtype']=='normal_part' for x in cfg.findall('.//part'))


def test_real_coverage_and_floating_warning_not_hidden(real):
    p,_=real;e=json.loads((p/'evaluation.json').read_text())
    coverage=json.loads((p/'coverage.json').read_text())
    assert any(r['too_close_support_layer_tops'] for r in coverage['rows'])
    assert not coverage['passed'] and e['coverage']['centerline_underbody_percent']<100
    assert e['feature_preservation']['passed']
    assert not e['passed'] and e['unresolved_warnings']


def test_travel_audit_uses_prior_real_extrusion_only(real):
    p,_=real;r=json.loads((p/'travel_near_miss.json').read_text());assert r['chronological']
    assert r['travel_segments']>0
    assert all(x['nearest_extrusion_event']<x['event_id'] for x in r['rows'])


def test_99_percent_cannot_unlock_production(real):
    p,_=real;r=json.loads((p/'geometry.json').read_text());e=json.loads((p/'evaluation.json').read_text());rem=json.loads((p/'removal.json').read_text());bad=copy.deepcopy(e);bad['coverage'].update(passed=False,centerline_underbody_percent=99.)
    gates,failures=hard_gates(r,bad,rem,{'all_unchanged':True});assert not gates['CENTERLINE_100'] and 'CENTERLINE_100' in failures


def test_native_slicing_never_uploads_or_starts_printer(real):
    p,_=real;s=json.loads((p/'bambu/slice_result.json').read_text())
    assert not s['printer_upload'] and not s['printer_start']
    assert '--export-3mf' in s['command']
    assert not any('upload' in x.lower() or 'start-print' in x.lower() for x in s['command'])

# Closure / safe-removal regressions. Previous s1-11 tests remain unchanged.
def test_coverage_local_closure_uses_actual_extrusion(context):
    from vascular_processing.core_support_inputs import default_output
    folder=default_output()/'validation/actual_horizontal_tube'
    g=pickle.load((folder/'geometry.pkl').open('rb'));fp=pickle.load((folder/'actual_footprints.pkl').open('rb'))['support']
    broken={z:p.difference(box(112,109,113,111)) for z,p in fp.items()}
    assert not actual_underbody(g['branches'],broken,g['profile'])['passed']
    assert actual_underbody(g['branches'],fp,g['profile'])['centerline_underbody_percent']==100.


def test_coverage_never_sacrificed_for_removal():
    from vascular_processing.core_support_safe_geometry import removal_preserves_coverage
    base=dict(coverage=dict(passed=True,centerline_underbody_percent=100.,port_underbody_percent=100.,max_uncovered_arc_mm=0.),layer_coverage=dict(passed=True,percent=100.),support_self=dict(self_unsupported_regions=0))
    bad=copy.deepcopy(base);bad['coverage']['centerline_underbody_percent']=99.999
    assert removal_preserves_coverage(base,base)
    assert not removal_preserves_coverage(base,bad)


def test_support_self_closure(context):
    p=context['profile'];h=p['layer_height'];z=p['first_layer_height']
    assert not actual_self({round(z+h,5):box(0,0,1,1)},p)['passed']
    assert actual_self({z:box(0,0,1,1),round(z+h,5):box(0,0,1,1)},p)['passed']


def test_planned_segment_actually_separates():
    from vascular_processing.core_support_safe_fixtures import fixture
    m=fixture('continuous_roof_two_segments');assert len(m.split())==1
    cutter=trimesh.creation.box([.22,2,8]);cutter.apply_translation([0,0,3])
    after=trimesh.boolean.difference([m,cutter],engine='manifold')
    components=[p for p in after.split() if p.volume>1e-6]
    assert len(components)==2 and all(p.extents[0]<6 for p in components)


def test_roof_seam_print_support_continuity():
    result=json.loads((default_output()/'closure_v2/validation/actual_thin_web_C3/thin_web_result.json').read_text())
    assert result['evaluation']['slice']['slicer_success']
    assert result['evaluation']['coverage']['centerline_underbody_percent']==100
    assert result['evaluation']['auto_support_off']


def test_weak_link_preserved_by_real_slicer():
    result=json.loads((default_output()/'closure_v2/validation/actual_thin_web_C3/thin_web_result.json').read_text())
    assert result['actual_preservation']['passed']
    assert result['actual_preservation']['rows']


def test_weak_link_cross_section_measured_from_mesh():
    from vascular_processing.core_support_safe_removal import cross_section_area,weak_ratio
    m=trimesh.creation.box([4,2,3])
    area=cross_section_area(m,[0,0,0],[0,0,1])
    assert area==pytest.approx(8.)
    assert weak_ratio(.44,.3,area)==pytest.approx(.0165)
    assert weak_ratio(.44,.3,None) is None


def test_local_vascular_cross_section_excludes_other_branches():
    from vascular_processing.core_support_safe_removal import cross_section_area
    a=trimesh.creation.box([2,2,4]);b=trimesh.creation.box([3,3,4]);b.apply_translation([8,0,0])
    mesh=trimesh.util.concatenate([a,b])
    assert cross_section_area(mesh,[0,0,0],[0,0,1])==pytest.approx(4.)
    assert cross_section_area(mesh,[8,0,0],[0,0,1])==pytest.approx(9.)
    assert cross_section_area(mesh,[4,0,0],[0,0,1]) is None


def test_sequential_tool_access_rechecks_actual_obstacles():
    from vascular_processing.core_support_safe_removal import tool_access
    vessel=trimesh.creation.box([1,1,1]);vessel.apply_translation([0,20,0])
    b=trimesh.creation.box([1,1,1]);a=b.copy();a.apply_translation([10,0,0])
    assert not tool_access([.5,0,0],vessel,{'A':a,'B':b},allowed=[[1,0,0]])['passed']
    assert tool_access([.5,0,0],vessel,{'B':b},allowed=[[1,0,0]])['passed']


def test_sequential_component_removal_updates_obstacles():
    from vascular_processing.core_support_safe_removal import sequential_plan
    parts={'B':'inner','A':'outer'};contacts={k:{'planned_cut_points':[]} for k in parts}
    def access(name,remaining):return dict(passed=name=='A' or 'A' not in remaining)
    seen=[]
    def release(name,obstacles):seen.append((name,sorted(obstacles)));return dict(passed=True)
    result=sequential_plan(parts,None,contacts,access,release)
    assert result['complete'] and [x['segment_id'] for x in result['steps']]==['A','B']
    assert seen==[('A',['B']),('B',[])]
    assert not result['all_notches_accessible_at_start_required']


def test_lever_arm_measurement_is_perpendicular_distance():
    from vascular_processing.core_support_safe_removal import lever_arm
    assert lever_arm([3,4,10],[0,0,0],[0,0,-1])==pytest.approx(5.)
    assert lever_arm([0,0,20],[0,0,0],[0,0,-1])==0.


def test_last_contact_piece_detection_distinguishes_long_lever():
    from vascular_processing.core_support_safe_fixtures import fixture
    assert fixture('short_contact_piece').extents[0]<3
    assert fixture('long_lever_support').extents[0]>20


def test_bifurcation_safe_direction():
    from vascular_processing.core_support_safe_removal import direction_check
    assert direction_check([0,0,-1],[0,0,-1],[1,0,0],True)['passed']
    assert not direction_check([1,0,0],[0,0,-1],[1,0,0],True)['passed']


def test_wrap_over_150_rejected():
    from vascular_processing.core_support_safe_removal import wrap_status
    assert wrap_status(90)=='ACCEPTABLE'
    assert wrap_status(120)=='ACCEPTABLE'
    assert wrap_status(140)=='HIGH_REMOVAL_RISK'
    assert wrap_status(160)=='WRAP_OVER_150_REJECTED'


def test_closed_cage_rejected_by_actual_section_topology():
    from vascular_processing.core_support_safe_fixtures import fixture
    from vascular_processing.core_support_safe_removal import closed_cage
    assert closed_cage(fixture('closed_ring_support'),[0,0,0],[1,0,0])
    assert not closed_cage(fixture('wrap_160'),[0,0,0],[1,0,0])


def test_removal_requires_vessel_deflection_rejected():
    from vascular_processing.core_support_safe_removal import release_path
    vessel=trimesh.creation.box([2,2,2]);support=trimesh.creation.box([1,1,1])
    result=release_path(support,vessel,{})
    assert not result['passed'] and not result['vessel_deflection_permitted']


def test_gap_025_not_selected_if_coverage_lower():
    from vascular_processing.core_support_closure import choose_production
    good=dict(coverage=dict(passed=True,centerline_underbody_percent=100.,port_underbody_percent=100.,max_uncovered_arc_mm=0.),layer_coverage=dict(passed=True,percent=100.),support_self=dict(self_unsupported_regions=0))
    bad=copy.deepcopy(good);bad['coverage']['centerline_underbody_percent']=99.9
    a=(dict(gap_mm=.20,support_volume_mm3=100),good,[])
    b=(dict(gap_mm=.25,support_volume_mm3=50),bad,[])
    assert choose_production([b,a])[0]['gap_mm']==.20


def test_closed_routes_compare_safety_before_gap_or_material():
    from vascular_processing.core_support_safe_removal import choose_safe_state
    def state(gap,volume,lever,coverage=True,sequence=True):
        return dict(record=dict(gap_mm=gap,support_volume_mm3=volume),gates=dict(COVERAGE_100=coverage,COMPLETE_SEQUENTIAL_REMOVAL=sequence),removal=dict(max_wrap_angle=80.,max_lever_arm_mm=lever,max_removal_component_arc_span_mm=2.,cut_count=35))
    a=state(.2,200.,1.);b=state(.25,10.,.5,coverage=False)
    assert choose_safe_state([b,a]) is a
    b=state(.25,10.,.5,sequence=False)
    assert choose_safe_state([b,a]) is a
    b=state(.25,10.,2.)
    assert choose_safe_state([b,a]) is a
    b=state(.25,250.,.5)
    assert choose_safe_state([a,b]) is b


def test_two_real_round_stagnation_not_one():
    from vascular_processing.core_support_closure import stagnated
    assert not stagnated([(22,1.3,45),(22,1.3,44)])
    assert not stagnated([(22,1.3,45),(22,1.3,44),(22,1.2,44)])
    assert stagnated([(22,1.3,45),(22,1.3,44),(23,1.4,40)])


def test_review_replays_later_real_result_after_runtime_failure(tmp_path,monkeypatch):
    import vascular_processing.core_support_closure as closure
    import vascular_processing.core_support_safe_removal as removal
    from vascular_processing.support_inputs import save
    base=tmp_path/'base';out=base/'closure_v2';(out/'reports/baseline').mkdir(parents=True)
    job=base/'candidates/BASE';job.mkdir(parents=True)
    e=dict(coverage=dict(unsupported_samples=10,centerline_underbody_percent=90.,port_underbody_percent=90.,max_uncovered_arc_mm=1.,passed=False),layer_coverage=dict(uncovered_area_sum_mm2=1.,percent=90.,passed=False),support_self=dict(self_unsupported_regions=1),slice=dict(slicer_success=True))
    save(job/'evaluation.json',e);save(job/'geometry.json',dict(candidate_id='BASE',folder=str(job)))
    save(out/'candidates/BASE_I01/runtime_error.json',dict(error='roundtrip failed'))
    trial=out/'candidates/BASE_I02';later=copy.deepcopy(e);later['coverage'].update(unsupported_samples=9,centerline_underbody_percent=91.);later['layer_coverage'].update(uncovered_area_sum_mm2=.9,percent=91.)
    save(trial/'evaluation.json',later);save(trial/'geometry.json',dict(candidate_id=trial.name,folder=str(trial),closure_parent=str(job)))
    monkeypatch.setattr(closure,'discover_source',lambda:{})
    monkeypatch.setattr(closure,'freeze',lambda *args:{})
    monkeypatch.setattr(closure,'current_profiles',lambda *args:({}, {}, {}))
    monkeypatch.setattr(closure,'baseline_jobs',lambda *args:[job])
    def finish(base,out,record,e,iterations,stop,*args,**kwargs):return dict(record=record,iterations=iterations,stop=stop)
    monkeypatch.setattr(removal,'finish',finish)
    result=closure.run(base,max_iterations=3,review_only=True)
    assert result['record']['candidate_id']=='BASE_I02'
    assert [r['iteration'] for r in result['iterations']]==[2]
    assert result['stop'][0]['iteration']==1
    assert 'RUNTIME_FAILURE' in result['stop'][0]['reason']


def test_removal_geometry_locked_before_coverage_closure(tmp_path):
    from vascular_processing.core_support_safe_geometry import optimize_after_closure
    with pytest.raises(ValueError,match='LOCKED_UNTIL_100'):
        optimize_after_closure({},dict(coverage={}),tmp_path,{}, {})


@pytest.mark.parametrize('name',['continuous_roof_two_segments','staggered_roof_seam','thin_web_seam','long_lever_support','short_contact_piece','bifurcation_pull_bad','bifurcation_downward_remove_good','wrap_90','wrap_120','wrap_160','closed_ring_support'])
def test_safe_removal_synthetic_fixtures_are_closed(name):
    from vascular_processing.core_support_safe_fixtures import fixture
    assert fixture(name).is_volume


def test_weak_link_dropped_by_slicer_is_rejected():
    result=json.loads((default_output()/'closure_v2/validation/actual_thin_web/thin_web_result.json').read_text())
    assert not result['actual_preservation']['passed']
    assert any(not r['preserved'] for r in result['actual_preservation']['rows'])


def test_weak_link_missing_intermediate_layer_is_rejected(context):
    from vascular_processing.core_support_safe_removal import preserved_feature
    feature=dict(witness_points=[[0,0,.12],[0,0,.22],[0,0,.32]])
    check=preserved_feature(feature,{.12:box(-1,-1,1,1),.32:box(-1,-1,1,1)},context['profile'])
    assert not check['passed']
    assert [r['preserved'] for r in check['rows']]==[True,False,True]


def test_ascii_support_fallback_does_not_expand_geometry():
    from vascular_processing.core_support_closure import canonical_support
    mesh=trimesh.creation.box([1e-6,.44,.3]);mesh.apply_translation([128,128,10])
    repaired,encoding=canonical_support(mesh)
    assert encoding=='stl_ascii' and repaired.is_volume
    assert np.max(np.abs(repaired.bounds-mesh.bounds))<1e-12


def test_removal_surface_normal_comes_from_frozen_facet_not_vertical_ray():
    from vascular_processing.core_support_safe_removal import actual_surface_normal
    from scipy.spatial.transform import Rotation
    mesh=trimesh.creation.box([2,2,.4]);M=np.eye(4);M[:3,:3]=Rotation.from_euler('y',30,degrees=True).as_matrix();mesh.apply_transform(M)
    i=int(np.argmin(mesh.face_normals[:,2]));normal,face=actual_surface_normal(mesh,mesh.triangles_center[i])
    assert normal==pytest.approx(mesh.face_normals[i])
    assert abs(normal[0])>.4 and normal[2]<-.8


def test_closure_delivery_keeps_actual_failures_and_separate_native_parts():
    from vascular_processing.core_support_safe_review import NAMES as QC_NAMES
    from PIL import Image
    import xml.etree.ElementTree as ET
    out=default_output()/'closure_v2'
    s=json.loads((out/'safe_removal_summary.json').read_text())
    required=['support_closure_iterations.csv','safe_removal_segments.csv','safe_removal_sequence.json','removal_load_paths.csv','weak_links.csv','final_contact_pieces.csv','bifurcation_removal_audit.csv','coverage_final.json','support_self_final.json','safe_removal_summary.json','SAFE_SUPPORT_REMOVAL_GUIDE.md','vascular_core_support_safe_removal_report.md']
    assert all((out/name).stat().st_size>0 for name in required)
    c=json.loads((out/'coverage_final.json').read_text())
    assert c['unsupported_samples']==sum(not r['supported'] for r in c['rows'])
    assert len(s['failures'])==sum(not value for value in s['gates'].values())
    sequence=json.loads((out/'safe_removal_sequence.json').read_text())
    if s['failures']:
        assert s['status']=='NEEDS_ADJUSTMENT' and s['production_selected'] is None
        assert sequence['production_plan'] is None
    assert s['source_protection']['all_unchanged'] and not s['source_protection']['changed']
    for name in ['vascular_core.stl','custom_support.stl']:
        assert trimesh.load_mesh(out/'selected'/name).is_volume
    for name in QC_NAMES:
        picture=np.asarray(Image.open(out/'QC_SAFE'/(name+'.png')).convert('RGB'))
        assert picture.shape[0]>100 and picture.std()>1
    with zipfile.ZipFile(out/'selected/X1C_project.3mf') as archive:
        assert archive.testzip() is None
        settings=json.loads(archive.read('Metadata/project_settings.config'))
        config=ET.fromstring(archive.read('Metadata/model_settings.config'))
        assert settings['enable_support']=='0'
        assert float(settings['layer_height'])==.1
        assert settings['nozzle_diameter']==['0.2']
        assert len(config.findall('.//part'))==2
        assert all(x.attrib['subtype']=='normal_part' for x in config.findall('.//part'))


def test_real_closure_stops_only_at_limit_or_measured_stagnation():
    from vascular_processing.core_support_closure import stagnated,baseline_jobs
    s=json.loads((default_output()/'closure_v2/safe_removal_summary.json').read_text())
    terminal=[r for r in s['stop_conditions'] if r['reason'] in ['ITERATION_LIMIT','TWO_CONSECUTIVE_REAL_ROUNDS_WITHOUT_SAMPLE_OR_AREA_IMPROVEMENT']]
    # These real routes did not close. Every route must have an explicit,
    # measured stop, instead of an early runtime error pretending to converge.
    if s['production_selected'] is None and not s['gates']['COVERAGE_100']:
        assert {r['baseline'] for r in terminal}=={p.name for p in baseline_jobs(default_output())}
    for row in terminal:
        if row['reason']=='ITERATION_LIMIT':assert row['iteration']==8
        else:assert stagnated(row['history'])


def test_final_load_path_points_are_on_actual_material_surfaces():
    out=default_output()/'closure_v2'
    diagnostic=json.loads((out/'safe_removal_diagnostic.json').read_text())
    vascular=trimesh.load_mesh(out/'selected/vascular_core.stl')
    for row in diagnostic['rows']:
        part=trimesh.load_mesh(out/'virtual_components'/(row['segment_id']+'.stl'))
        assert part.nearest.on_surface([row['force_point']])[1][0]<1e-4
        assert vascular.nearest.on_surface([row['contact_point']])[1][0]<1e-4
        assert row['closed_cage_checked_samples']>0


# Final practical closure regressions. Historical strict tests above remain
# unchanged: they describe the earlier algorithm, not the new READY policy.
def practical_profile():
    return dict(roof_gap_mm=.25,layer_height=.1,first_layer_height=.12,
                support_line_width=.22,line_width=.22,gcode_xy_half_cell_mm=.0005)


def practical_branch(points):
    q=np.asarray(points,dtype=float)
    return dict(branch_id='TEST',lower=q,xyz=q+[0,0,1],is_port=False,
                arc=np.r_[0,np.linalg.norm(np.diff(q,axis=0),axis=1).cumsum()])


def test_practical_true_missing_vs_too_close_classification():
    from vascular_processing.core_support_practical_audit import structural_samples
    b=practical_branch([[0,0,.4],[1,0,.4],[2,0,.4]])
    r=structural_samples([b],{.22:box(-.1,-.1,.1,.1)},practical_profile())
    assert r['rows'][0]['new_class']=='SUPPORT_TOO_CLOSE'
    assert not r['rows'][0]['true_missing_support']
    assert r['true_missing_support_count']==2


def test_practical_too_far_never_borrows_low_wall():
    from vascular_processing.core_support_practical_audit import structural_samples
    r=structural_samples([practical_branch([[0,0,1.]])],{.12:box(-1,-1,1,1)},practical_profile())
    assert r['rows'][0]['new_class']=='SUPPORT_TOO_FAR'
    assert r['true_missing_support_count']==1
    assert r['rows'][0]['farther_support_witness']['gap_mm']==pytest.approx(.88)


def test_practical_numerical_boundary_exception():
    from vascular_processing.core_support_practical_audit import structural_samples
    b=practical_branch([[-.1,-.01,.42],[0,.001,.42],[.1,-.01,.42]])
    r=structural_samples([b],{.12:box(-1,-1,1,0)},practical_profile())
    assert r['numerical_exceptions']==1 and r['true_missing_support_count']==0
    assert r['rows'][1]['exception_evidence']['same_actual_connected_footprint']
    b['lower'][1,1]=.01
    r=structural_samples([b],{.12:box(-1,-1,1,0)},practical_profile())
    assert r['numerical_exceptions']==0 and r['true_missing_support_count']==1


def test_practical_numerical_two_point_gap_not_exempted():
    from vascular_processing.core_support_practical_audit import structural_samples
    b=practical_branch([[-.2,-.01,.42],[-.1,.001,.42],[.1,.001,.42],[.2,-.01,.42]])
    r=structural_samples([b],{.12:box(-1,-1,1,0)},practical_profile())
    assert r['numerical_exceptions']==0 and r['true_missing_support_count']==2


def test_practical_coverage_repair_only_true_missing():
    out=default_output()/'final_practical_v3'
    for i in [1,2,3]:
        folder=out/'candidates'/f'PRACTICAL_REPAIR_R{i}'
        g=json.loads((folder/'geometry.json').read_text());parent=Path(g['practical_parent'])
        audit_path=parent/'practical/final_true_missing_support.json'
        if not audit_path.exists():audit_path=out/'baseline_audits'/parent.name/'final_true_missing_support.json'
        old=json.loads(audit_path.read_text())
        missing={r['sample_id'] for r in old['rows'] if r['true_missing_support']}
        repair=json.loads((folder/'true_missing_repairs.json').read_text())
        assert all(r['sample_id'] in missing for r in repair['changes'])
        assert not repair['too_close_targets_modified'] and not repair['global_roof_regeneration']
        assert digest(folder/'vascular_core.stl')==digest(parent/'vascular_core.stl')


def test_practical_micro_self_deficit_not_macro_fail():
    from vascular_processing.core_support_practical_audit import macro_printability
    p=practical_profile();fp={.12:box(0,0,1,1),.22:box(0,0,1.12,1)}
    strict=actual_self(fp,p)
    r=macro_printability(fp,p,strict_rows=strict['rows'])
    assert r['micro_deficit_count']>0 and r['macro_fail_count']==0


def test_practical_floating_support_island_macro_fail():
    from vascular_processing.core_support_practical_audit import macro_printability
    p=practical_profile();fp={.12:box(0,0,1,1),.22:box(3,0,4,1)}
    r=macro_printability(fp,p,strict_rows=actual_self(fp,p)['rows'])
    assert r['floating_support_islands']>0 and r['macro_fail_count']>0


def test_practical_missing_layer_is_floating_not_micro():
    from vascular_processing.core_support_practical_audit import macro_printability
    p=practical_profile();r=macro_printability({.12:box(0,0,1,1),.42:box(0,0,1,1)},p)
    assert r['floating_support_islands']==1


def test_practical_vascular_start_requires_real_support():
    from vascular_processing.core_support_practical_audit import vascular_layer_starts
    p=practical_profile();v={.52:box(0,0,1,1)}
    assert vascular_layer_starts(v,{},p)['floating_vascular_region_count']==1
    assert vascular_layer_starts(v,{.22:box(0,0,1,1)},p)['floating_vascular_region_count']==0


def practical_bar():
    mesh=trimesh.creation.box([8,.8,1]);mesh.apply_translation([0,0,.5]);b=practical_branch(np.c_[np.linspace(-4,4,81),np.zeros(81),np.full(81,1.25)])
    return mesh,dict(branches=[b],profile=practical_profile())


def test_practical_strategic_cut_increases_components():
    from vascular_processing.core_support_practical_geometry import material_components
    from vascular_processing.core_support_practical_topology import cut_effect
    from vascular_processing.core_support_geometry import boolean
    mesh,g=practical_bar();cutter=trimesh.creation.box([.22,2,3]);pieces=material_components(boolean([mesh,cutter],'difference'))
    result=cut_effect([mesh],pieces,g)
    assert len(pieces)==2 and result['effective']
    assert result['new_max_contact_arc_mm']<result['old_max_contact_arc_mm']


def test_practical_ineffective_cut_rejected():
    from vascular_processing.core_support_practical_topology import cut_effect
    mesh,g=practical_bar();assert not cut_effect([mesh],[mesh.copy()],g)['effective']


def test_practical_thin_web_C3_preserved():
    p=default_output()/'final_practical_v3/candidates/PRACTICAL_TOPOLOGY_R1/practical_removal/removal_audit.json'
    r=json.loads(p.read_text());assert r['thin_web_checks']
    assert all(x['fixture']=='actual_thin_web_C3' and x['passed'] for x in r['thin_web_checks'])


def test_practical_major_component_split_keeps_internal_cavity():
    from vascular_processing.core_support_practical_geometry import material_components
    from vascular_processing.core_support_geometry import boolean
    outer=trimesh.creation.box([3,3,3]);inner=trimesh.creation.box([1,1,1]);hollow=boolean([outer,inner],'difference')
    parts=material_components(hollow)
    assert len(parts)==1 and parts[0].volume==pytest.approx(26)
    assert not parts[0].contains([[0,0,0]])[0]


def test_practical_last_contact_piece_bifurcation_relative_priority():
    from vascular_processing.core_support_practical_removal import last_piece_check
    assert last_piece_check(5,8)['passed']
    assert not last_piece_check(5,8,True)['passed']
    assert not last_piece_check(40,8)['passed']


def ready_practical_input():
    a=dict(structural=dict(true_missing_support_count=0,max_true_unsupported_arc_mm=0,structural_support_percent=99.9,gap_quality_counts=dict(TOO_CLOSE=10)),macro=dict(macro_fail_count=0,micro_deficit_count=20),floating_vascular_regions_unresolved=0,native_slice_pass=True,auto_support_off=True,vascular_features_preserved=True,warnings=dict(critical_unresolved=0))
    r=dict(closed_cage_count=0,major_components_successfully_separated=True,sequence=dict(complete=True),short_final_contact_pieces=True,thin_web_preserved=True)
    return a,r


def test_practical_closed_cage_hard_fail():
    from vascular_processing.core_support_practical_removal import practical_gates,wrap_risk
    a,r=ready_practical_input();r['closed_cage_count']=1
    assert not all(practical_gates(a,r).values()) and wrap_risk(100,True)=='CLOSED_CAGE'
    assert wrap_risk(140)=='MODERATE_HIGH' and wrap_risk(160)=='HIGH_RISK'


def test_practical_sequential_access_after_removal():
    from vascular_processing.core_support_safe_removal import sequential_plan
    parts={'A':object(),'B':object()};seen=[]
    def access(name,remaining):
        seen.append((name,set(remaining)));return dict(passed=name=='A' or 'A' not in remaining)
    r=sequential_plan(parts,None,{n:dict(planned_cut_points=[]) for n in parts},access,lambda n,o:dict(passed=True),dict(A=True,B=True))
    assert r['complete'] and ('B',{'B'}) in seen


def test_practical_bifurcation_no_forced_lateral_pull():
    from vascular_processing.core_support_safe_removal import direction_check
    assert direction_check([0,0,-1],[0,0,-1],[1,0,0],True)['passed']
    assert not direction_check([1,0,-.1],[0,0,-1],[0,1,0],True)['passed']


def test_practical_ready_without_mathematical_100pct():
    from vascular_processing.core_support_practical_removal import practical_gates
    a,r=ready_practical_input();assert all(practical_gates(a,r).values())
    a['structural']['true_missing_support_count']=1
    assert not all(practical_gates(a,r).values())


def test_practical_warning_not_critical_for_verified_gap_only():
    from vascular_processing.core_support_practical_audit import warning_classification
    assert warning_classification(['floating region'],dict(true_missing_support_count=0),dict(macro_fail_count=0),0)['critical_unresolved']==0
    assert warning_classification(['floating region'],dict(true_missing_support_count=1),dict(macro_fail_count=0),0)['critical_unresolved']==1


def test_practical_round_limit_rejects_fourth():
    from vascular_processing.core_support_practical_pipeline import run
    with pytest.raises(ValueError,match='THREE_REPAIR'):run(max_rounds=4)


def test_practical_contact_measure_uses_surface_with_sparse_triangles():
    from vascular_processing.core_support_practical_topology import contact_measure
    mesh=trimesh.creation.box([8,2,1]);mesh.apply_translation([0,0,.5])
    b=practical_branch(np.c_[np.linspace(-4,4,81),np.zeros(81),np.full(81,1.25)])
    r=contact_measure([mesh],dict(branches=[b],profile=practical_profile()),exact=True)[0]
    assert r['contact_arc_mm']==pytest.approx(8.)


def test_practical_bulk_cut_rejects_other_vascular_branch():
    from vascular_processing.core_support_practical_removal import bulk_cut_proposal
    part=trimesh.creation.box([8,1,4]);part.apply_translation([0,0,2])
    other=trimesh.creation.box([1,1,1]);other.apply_translation([0,0,3.25])
    b=practical_branch([[-3,0,5],[3,0,5]])
    r=bulk_cut_proposal(part,dict(component_id='M00',sample_indices=[0,1]),dict(branches=[b],profile=practical_profile(),vascular=other),dict(M00=part))
    assert r['status']=='FINITE_BULK_BLADE_INTERSECTS_ANOTHER_VASCULAR_BRANCH'
    assert not r['executed'] and r['finite_blade_vascular_intersection_mm3']>0


def test_practical_incomplete_bulk_sequence_cannot_be_ready():
    from vascular_processing.core_support_practical_removal import practical_gates
    a,r=ready_practical_input();r['sequence']['complete']=False
    assert not practical_gates(a,r)['COMPLETE_SAFE_REMOVAL']


def test_practical_preserved_parent_curtain_is_not_isolated_C3():
    p=default_output()/'final_practical_v3/candidates/PRACTICAL_TOPOLOGY_R1/practical_removal/removal_audit.json'
    r=json.loads(p.read_text())
    assert r['actual_extrusion_thin_link_witnesses_preserved']
    assert any(not x['isolated_C3_geometry'] for x in r['thin_web_checks'])
    assert not r['thin_web_preserved']


def test_practical_prepared_extrusion_clipping_is_exact():
    import shapely
    from vascular_processing.core_support_slicing import clip_extrusion_lines
    region=box(0,0,4,4).difference(box(1,1,3,3))
    lines=shapely.linestrings(np.array([[[.2,.2],[.8,.8]],[[-1,2],[5,2]],[[2,1.5],[2,2.5]],[[-1,-1],[0,0]],[[0,0],[4,0]]]))
    expected=shapely.intersection(lines,region)
    actual=clip_extrusion_lines(lines,region)
    assert all(shapely.equals(a,b) for a,b in zip(actual,expected))


def test_practical_archive_recovery_never_calls_native_again(monkeypatch):
    import vascular_processing.core_support_slicing as slicing
    out=default_output()/'final_practical_v3';folder=out/'candidates/PRACTICAL_TOPOLOGY_R3';profiles=json.loads((out/'bambu/profiles/profiles.json').read_text())
    original=digest(folder/'bambu/input.3mf')
    def forbidden(*args,**kwargs):pytest.fail('Archive recovery must not write or slice a new project')
    class StopAfterArchiveVerification(Exception):pass
    monkeypatch.setattr(slicing,'write_project',forbidden);monkeypatch.setattr(slicing,'slice_local',forbidden)
    def parsed(*args):raise StopAfterArchiveVerification()
    monkeypatch.setattr(slicing,'events',parsed)
    with pytest.raises(StopAfterArchiveVerification):
        slicing.evaluate(json.loads((folder/'geometry.json').read_text()),profiles['actual'],profiles['engine'],reuse_native=True)
    assert digest(folder/'bambu/input.3mf')==original


def test_practical_final_roof_residue_excludes_bed_curtain_shells():
    out=default_output()/'final_practical_v3';summary=json.loads((out/'final_practical_summary.json').read_text());folder=Path(summary['candidate_folder']);audit=json.loads((folder/'practical_removal/removal_audit.json').read_text());g=pickle.load((folder/'geometry.pkl').open('rb'))
    assert audit['audit_version'].endswith('construction-roof-residue')
    for row in audit['final_contact_pieces']:
        # Match the established ASCII precision policy: default 1e-8 welding
        # collapses tiny Boolean seams and changes these diagnostic surfaces.
        mesh=trimesh.load_mesh(folder/'practical_removal/contact_residues'/(row['piece_id']+'.stl'),process=False)
        mesh.merge_vertices(digits_vertex=14)
        assert mesh.is_watertight and mesh.volume>0
        assert mesh.bounds[0,2]>g['profile']['first_layer_height']


def test_practical_delivery_integrity_and_no_false_production_freeze():
    out=default_output()/'final_practical_v3';summary=json.loads((out/'final_practical_summary.json').read_text());manifest=json.loads((out/'production_candidate_freeze.json').read_text())
    assert summary['repair_round_count']<=3 and summary['topology_round_count']<=3
    assert summary['actual_native_slice_attempt_count']<=6
    assert manifest['vascular_source_sha256']==digest(discover_source()['core'])
    assert manifest['production_frozen']==all(summary['gates'].values())
    if summary['failures']:assert manifest['production_candidate'] is None
    assert verify({str(out/p):h for p,h in manifest['files'].items()})['all_unchanged']
    for name in ['final_project.3mf','final_sliced.3mf']:
        with zipfile.ZipFile(out/'selected'/name) as z:assert z.testzip() is None
    assert len(list((out/'QC_FINAL').glob('*.png')))==18
