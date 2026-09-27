"""Permanent synthetic and delivered-file tests for the four-port review stage."""
import copy
import json
from pathlib import Path
import numpy as np
import pytest
import trimesh

from vascular_processing import sacrificial_fixture as f
from vascular_processing import port_attachment_alignment as alignment
from vascular_processing import surface_continuity_qc as qc
from vascular_processing.all_port_attachment_run import PORTS, WALLS, saved_routes

BASE=f.ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi'
OUT=BASE/'print_fixture_design_all_ports_aligned'


@pytest.fixture(scope='module')
def cfg():
    return f.load_config(f.ROOT/'config/sacrificial_box_BG001.yaml')


@pytest.fixture(scope='module')
def curved_four_terminal_source():
    """One curved solid tree, four explicit 24-fan caps, actual last-5mm fits."""
    arms=[];endpoints=[]
    for j,degrees in enumerate((0,20,40,60)):
        radial=np.array([np.cos(j*np.pi/2),np.sin(j*np.pi/2),0.])
        vertical=np.array([0.,0.,1.]);azimuth=np.cross(vertical,radial)
        s=np.r_[np.linspace(0,8,17)[:-1],np.linspace(8,13,101)]
        theta=2*np.radians(degrees)
        centers=[];tangents=[]
        for arc in s:
            if arc <= 8 or theta == 0:
                center=arc*radial;tangent=radial
            else:
                angle=(arc-8)*theta/5;radius=5/theta
                center=8*radial+radius*(np.sin(angle)*radial+(1-np.cos(angle))*vertical)
                tangent=np.cos(angle)*radial+np.sin(angle)*vertical
            centers.append(center);tangents.append(tangent)
        centers=np.array(centers);tangents=np.array(tangents);vertices=[]
        for p,t in zip(centers,tangents):
            second=np.cross(t,azimuth)
            vertices.extend(p+.6*(np.cos(a)*azimuth+np.sin(a)*second) for a in np.linspace(0,2*np.pi,24,endpoint=False))
        vertices.extend([centers[0],centers[-1]]);faces=[]
        for k in range(len(s)-1):
            for a in range(24):
                b=(a+1)%24;lo=k*24;hi=(k+1)*24
                faces.extend([[lo+a,lo+b,hi+b],[lo+a,hi+b,hi+a]])
        for a in range(24):
            b=(a+1)%24;last=(len(s)-1)*24
            faces.extend([[len(vertices)-2,b,a],[len(vertices)-1,last+a,last+b]])
        arm=trimesh.Trimesh(vertices,faces,process=True);arm.fix_normals();assert arm.is_volume
        arms.append(arm)
        fitted=np.linalg.lstsq(np.column_stack([s[-101:]-8,np.ones(101)]),centers[-101:],rcond=None)[0][0]
        mean=f.unit(fitted)
        e=f.Endpoint(f'TERMINAL_{j}', 'outlet',j+1,j+1,f'arm{j}',centers[-1],.6,mean,5.,101,.6,.6,[j+1])
        endpoints.append(e)
    source=trimesh.boolean.union([trimesh.creation.icosphere(subdivisions=2,radius=1.)]+arms,engine='manifold')
    assert source.is_volume and len(source.split())==1
    return source,endpoints


@pytest.mark.parametrize('index,expected',enumerate((0,20,40,60)))
def test_generic_curved_cap_alignment_non_o3(cfg,curved_four_terminal_source,index,expected):
    source,eps=curved_four_terminal_source;e=copy.deepcopy(eps[index]);before=source.vertices.copy()
    result=f.align_known_cap(e,source,{'vertices':24},cfg)
    assert result['before_axis_error_deg']==pytest.approx(expected,abs=.002)
    assert result['provenance']['unique'] and result['cap_vertex_count']==24
    assert result['collar_contained'] and result['collar_outside_source_volume_mm3']<=1e-6
    assert result['production_profile']=='ACTUAL_CAP_POLYGON'
    assert result['full_radius_backward_overlap_used'] is False
    assert e.radius==.6 and np.linalg.norm(e.position-eps[index].position)<1e-4
    assert result['maximum_plane_projection_mm']<1e-4
    np.testing.assert_array_equal(source.vertices,before)


@pytest.mark.parametrize('volumes,chosen', [([0.],.2),([2e-6,0.],.15),([2e-6,2e-6,0.],.1)])
def test_collar_fallback_is_descending_and_fixed_fraction(cfg,curved_four_terminal_source,monkeypatch,volumes,chosen):
    source,eps=curved_four_terminal_source;calls=[];v=iter(volumes)
    def measured(collar,mesh,settings):
        calls.append(settings['attachment_alignment']['collar_inner_radius_fraction']);return next(v)
    monkeypatch.setattr(alignment,'collar_outside_volume',measured)
    result=f.align_known_cap(copy.deepcopy(eps[0]),source,{'vertices':24},cfg)
    assert result['collar_length_mm']==chosen
    assert [r['length_mm'] for r in result['collar_attempts']]==[.2,.15,.1][:len(volumes)]
    assert calls==[.65]*len(volumes)


def test_failed_shortest_collar_stops_without_old_overlap(cfg,curved_four_terminal_source,monkeypatch):
    source,eps=curved_four_terminal_source;e=copy.deepcopy(eps[0]);original=e.position.copy()
    monkeypatch.setattr(alignment,'collar_outside_volume',lambda *a:2e-6)
    with pytest.raises(ValueError,match='CAP_COLLAR_NOT_FEASIBLE'):
        f.align_known_cap(e,source,{'vertices':24},cfg)
    assert e.attachment_collar is None and e.cap_ring is None
    np.testing.assert_array_equal(e.position,original)


def test_ambiguous_center_fails_provenance(cfg,curved_four_terminal_source):
    source,eps=curved_four_terminal_source;source=source.copy()
    source.vertices=np.vstack([source.vertices,eps[0].position+[1e-5,0,0]])
    with pytest.raises(ValueError,match='CAP_PROVENANCE_UNRESOLVED_TERMINAL_0'):
        f.align_known_cap(copy.deepcopy(eps[0]),source,{'vertices':24},cfg)


def test_actual_ring_sweep_and_radius_after_stl_roundtrip(cfg,curved_four_terminal_source,tmp_path):
    source,eps=curved_four_terminal_source;e=copy.deepcopy(eps[0]);f.align_known_cap(e,source,{'vertices':24},cfg)
    target=e.position+20*e.tangent;box=f.BoxBounds(np.array([[-50,-50,-50],[33,50,50.]]),np.array([[-53,-53,-53],[36,53,50.]]))
    route=f.make_route(e,'+X',target,box,cfg,'ROUTE_A_STRAIGHT')
    path=tmp_path/'port.stl';f.export_checked_stl(route.mesh,path,cfg);route.mesh=trimesh.load_mesh(path)
    result=f.attachment_section_qc(route,cfg)
    assert result['passed'] and result['contour_hausdorff_mm']<1e-4 and result['center_offset_mm']<1e-4


def step_mesh():
    # A real annular shoulder has a sharp closed circle, unlike faceted sides.
    return trimesh.creation.revolve(np.array([[0,-2],[.6,-2],[.6,0],[1.,0],[1.,2],[0,2]]),sections=48)


def test_closed_sharp_ring_and_circumference_fraction():
    mesh=step_mesh();features,poly=qc.feature_edges(mesh,np.zeros(3),[0,0,1],1.2,.5,20)
    assert features['circumferential_ring_count']==2
    fraction=qc.circumferential_fraction(poly,np.zeros(3),[0,0,1],.6)
    assert fraction['sharp_feature_circumference_fraction']>2


def test_longitudinal_facets_are_not_circumferential_seams():
    mesh=trimesh.creation.cylinder(radius=.6,height=4,sections=24)
    features,poly=qc.feature_edges(mesh,np.zeros(3),[0,0,1],1.,1.,10)
    # Long edges are intentionally clipped by the local window; use explicit lines.
    import pyvista as pv
    poly=pv.PolyData(np.array([[.6,0,-.2],[.6,0,.2]]),lines=np.array([2,0,1]))
    assert qc.circumferential_fraction(poly,np.zeros(3),[0,0,1],.6)['sharp_feature_circumference_fraction']==0


def test_normal_jump_and_cross_section_measure_real_geometry():
    mesh=step_mesh();metrics,_=qc.normal_jumps(mesh,np.zeros(3),[0,0,1],1.2,2.)
    assert metrics['maximum_deg']==pytest.approx(90) and metrics['p95_deg']>0
    rows=[qc.cross_section(mesh,[0,0,s],[0,0,1],.1) for s in [-.2,-.1,.1,.2]]
    assert qc.area_jump(rows)>1 and rows[0]['equivalent_radius_mm']==pytest.approx(.6,abs=.002)


def test_gate_does_not_hide_a_ring_with_low_p95():
    before=dict(features=[{},dict(circumferential_ring_count=0,total_length_mm=0)],normal_jumps=dict(p95_deg=15))
    after=copy.deepcopy(before);after['maximum_area_jump_fraction']=.11;after['features'][1]['circumferential_ring_count']=1
    result=qc.all_port_gate('O1',before,after,4)
    assert not result['passed'] and 'CROSS_SECTION_VARIATION_WARNING' in result['warnings']


@pytest.fixture(scope='module')
def real_inputs(cfg):
    return f.load_inputs(cfg)


@pytest.fixture(scope='module')
def real_summary():
    assert (OUT/'all_port_attachment_summary.json').exists(), 'Run s1-4 --all-ports-aligned before real artifact tests'
    return json.loads((OUT/'all_port_attachment_summary.json').read_text())


@pytest.fixture(scope='module')
def real_core():
    return trimesh.load_mesh(OUT/'core/BG001_RMCA_BALANCED_core_with_ports_all_aligned.stl')


@pytest.mark.parametrize('port',PORTS)
def test_real_all_four_caps_have_unique_provenance_and_radii(real_inputs,port):
    e=next(e for e in real_inputs['endpoints'] if e.endpoint_id==port);q=e.attachment_qc
    assert q['provenance']['unique'] and q['cap_vertex_count']==24
    assert q['collar_contained'] and q['collar_outside_source_volume_mm3']<=1e-6
    assert q['maximum_plane_projection_mm']<1e-4 and q['center_adjustment_mm']<1e-4
    assert not q['radius_changed'] and not q['frozen_source_modified']


def test_real_i1_projection_uses_known_center_not_looser_tolerance(real_inputs,cfg):
    e=real_inputs['endpoints'][0];q=e.attachment_qc
    assert q['plane_anchor']=='KNOWN_CAP_CENTER_VERTEX_AXIAL_PLANE'
    assert q['maximum_plane_projection_mm']<7e-6
    assert cfg['attachment_alignment']['collar_outside_volume_tolerance_mm3']==1e-6
    assert q['collar_length_mm']==.2 and q['collar_outside_source_volume_mm3']<1e-6


def test_real_o3_geometry_is_identical_to_prior_alignment(real_inputs):
    old=json.loads((BASE/'print_fixture_design_aligned/attachment_alignment_qc.json').read_text())['O3']
    e=next(e for e in real_inputs['endpoints'] if e.endpoint_id=='O3');new=e.attachment_qc
    for key in ('aligned_ring_print_mm','cap_normal','collar_length_mm','collar_inner_radius_fraction','collar_outside_source_volume_mm3','radius_mm'):
        np.testing.assert_array_equal(new[key],old[key])


@pytest.mark.parametrize('port',PORTS)
def test_real_fixed_targets_radii_and_clearances(real_summary,real_inputs,port):
    d=real_summary['per_port'][port];row=saved_routes(BASE/'print_fixture_design_aligned')[port]
    assert row['face']==WALLS[port]
    assert d['wall_target_mm']==[float(row['target_'+k]) for k in 'xyz']
    assert not d['radius_changed'] and d['clearance_mm'] >= (1.6 if port=='O3' else 3)
    assert d['minimum_bend_radius_mm']>=6.12 and d['section_matches_actual_ring']['passed']
    e=next(e for e in real_inputs['endpoints'] if e.endpoint_id==port)
    assert d['radius_mm']==e.radius


@pytest.mark.parametrize('port',PORTS)
def test_real_continuity_improvement_and_o3_regression(real_summary,port):
    d=real_summary['per_port'][port];old,new=d['baseline_continuity'],d['new_continuity']
    assert d['gate']['passed']
    assert new['features'][1]['circumferential_ring_count']==0
    if port=='O3':
        assert new['normal_jumps']['p95_deg']<=old['normal_jumps']['p95_deg']+1
        assert new['features'][1]['total_length_mm']==0
    else:
        assert new['normal_jumps']['p95_deg']<=20 or new['normal_jumps']['p95_deg']<=old['normal_jumps']['p95_deg']*.5
        assert new['features'][1]['sharp_feature_circumference_fraction']<old['features'][1]['sharp_feature_circumference_fraction']


def test_real_baseline_is_actual_byte_identical_old_file(real_summary):
    b=real_summary['baseline'];assert b['sha256']==b['copy_sha256']==f.sha256(Path(b['path']))
    assert b['rebuilt'] is False


def test_real_baseline_normal_and_feature_metrics_reproducible(real_summary,real_inputs):
    mesh=trimesh.load_mesh(real_summary['baseline']['path'])
    for e in real_inputs['endpoints']:
        old=real_summary['per_port'][e.endpoint_id]['baseline_continuity']
        normals,_=qc.normal_jumps(mesh,e.position,e.tangent,e.radius*1.7,1.)
        edges,_=qc.feature_edges(mesh,e.position,e.tangent,e.radius*1.7,1.,20)
        assert normals['p95_deg']==pytest.approx(old['normal_jumps']['p95_deg'],abs=1e-8)
        assert edges['total_length_mm']==pytest.approx(old['features'][1]['total_length_mm'],abs=1e-8)


@pytest.mark.parametrize('port',PORTS)
def test_real_final_metrics_and_cross_sections_from_delivered_file(real_summary,real_core,real_inputs,port):
    e=next(e for e in real_inputs['endpoints'] if e.endpoint_id==port)
    new=real_summary['per_port'][port]['new_continuity'];normals,_=qc.normal_jumps(real_core,e.position,e.tangent,e.radius*1.7,1.)
    assert normals['p95_deg']==pytest.approx(new['normal_jumps']['p95_deg'],abs=1e-8)
    rows=f.read_csv(OUT/'tables'/(port+'_cross_section_profile.csv'))
    for version in ('baseline','aligned'):
        group=[r for r in rows if r['version']==version]
        assert len(group)==16 and float(group[0]['arc_mm'])==-.5 and float(group[-1]['arc_mm'])==1.
        assert all(float(r['area_mm2'])>0 for r in group)
    assert new['maximum_area_jump_fraction']<=.1 or 'CROSS_SECTION_VARIATION_WARNING' in real_summary['per_port'][port]['gate']['warnings']


def test_real_final_core_watertight_manifold_one_component_and_source_preserved(real_core,real_inputs,cfg,real_summary):
    result=f.mesh_qc(real_core,cfg)
    assert result['passed'] and result['connected_components']==1 and result['degenerate_faces']==0
    missing=trimesh.boolean.difference([real_inputs['mesh'],real_core],engine='manifold')
    volume=abs(float(missing.volume)) if len(missing.faces) else 0.
    assert volume <= .001 and volume==pytest.approx(real_summary['source_volume_loss_mm3'],abs=1e-8)


def test_real_staged_union_and_pair_clearance(real_summary):
    assert [q['stage'] for q in real_summary['staged_union']]==['ADD_'+p for p in PORTS]
    assert all(q['passed'] and q['connected_components']==1 for q in real_summary['staged_union'])
    assert len(real_summary['pairwise_clearance'])==6
    assert min(q['clearance_mm'] for q in real_summary['pairwise_clearance'])>=3


def test_real_protected_hashes_and_no_upstream_or_final_printing(real_summary):
    evidence=json.loads((OUT/'protected_source_hashes.json').read_text())
    assert evidence['all_unchanged'] and real_summary['protected_sources_unchanged']
    assert real_summary['protected_file_count']>1000
    assert not real_summary['box_union_performed'] and not real_summary['bambu_slicing_performed']
    assert not real_summary['production_vmtk_ramp_used'] and not real_summary['upstream_recomputed']
    assert real_summary['source_sha256']==f.sha256(Path(real_summary['source_stl']))


def test_real_review_package_has_all_images_same_camera_and_chinese_report(real_summary):
    from PIL import Image
    files=real_summary['visual_review']['png_files'];assert len(files)==14
    for name in files:
        with Image.open(OUT/'QC'/name) as image:
            assert min(image.size)>=600
    cameras=json.loads((OUT/'camera_settings.json').read_text())
    assert cameras['before_after_identical'] and set(cameras['per_port'])==set(PORTS)
    assert real_summary['status']=='READY_FOR_HUMAN_REVIEW'
    report=(OUT/'all_port_attachment_report.md').read_text()
    assert all(f'## {n}.' in report for n in range(1,9))


def test_local_numerical_weld_removes_float_slivers_without_remeshing(cfg):
    mesh=trimesh.creation.revolve(np.array([[0,-3],[.6,-3],[.6,0],[.600002,0],[.600002,.1],[.6,.1],[.6,3],[0,3]]),sections=24)
    source=trimesh.creation.cylinder(radius=.6,height=6,sections=24)
    e=f.Endpoint('SYNTHETIC','outlet',1,1,'b',np.zeros(3),.6,np.array([0.,0.,1.]),5.,10,.6,.6,[1])
    before=mesh.vertices.copy();frozen=source.vertices.copy()
    cleaned,audit=alignment.clean_boolean_attachment_slivers(mesh,source,[e],cfg)
    assert audit['merged_vertices']>0 and audit['maximum_displacement_mm']<=1e-5
    assert audit['mesh_qc']['passed'] and not audit['whole_model_remeshed']
    np.testing.assert_array_equal(mesh.vertices,before);np.testing.assert_array_equal(source.vertices,frozen)
    outside=before[np.abs(before[:,2])>1]
    assert all(np.any(np.all(cleaned.vertices==p,axis=1)) for p in outside)
    feature,_=qc.feature_edges(cleaned,np.zeros(3),e.tangent,1.,1.,20)
    assert feature['total_length_mm']==0


def test_real_boolean_weld_stays_inside_existing_numerical_limit(real_summary,cfg):
    audit=real_summary['numerical_seam_cleanup']
    assert audit['maximum_displacement_mm']<=cfg['geometry']['boolean_cleanup_tolerance_mm']
    assert not audit['source_modified'] and not audit['outside_seam_vertices_modified']
    assert not audit['whole_model_remeshed']
    assert all(d['new_continuity']['features'][1]['total_length_mm']==0 for d in real_summary['per_port'].values())


def test_cap_ring_center_must_match_identified_swc_center(cfg):
    source=trimesh.creation.cylinder(radius=.6,height=4,sections=24)
    v=source.vertices.copy();rim=(v[:,2]>1.99)&(np.linalg.norm(v[:,:2],axis=1)>.5)
    v[rim,0]+=.01;source.vertices=v
    e=f.Endpoint('SHIFTED','outlet',1,1,'b',np.array([0.,0.,2.]),.6,np.array([0.,0.,1.]),5.,10,.6,.6,[1])
    with pytest.raises(ValueError,match='KNOWN_CAP_RING_CENTER_MISMATCH'):
        f.align_known_cap(e,source,{'vertices':24},cfg)
