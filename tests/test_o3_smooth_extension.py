"""Permanent O3 preflight regressions; production ramp tests require that runtime.

No passing test below claims that the clinical O3 seam has been corrected.
"""
from pathlib import Path
import json
import sys

import numpy as np
import pytest
import trimesh
from scipy.spatial.transform import Rotation

from vascular_processing import surface_continuity_qc as qc
from vascular_processing.o3_vmtk_bridge import invoke, require_ramp
from vascular_processing.o3_smooth_preflight import bridge_smoke
from vascular_processing.topbrain_qc import sha256

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi/print_fixture_design_o3_smooth'
PYTHON=Path('/mnt/d/anaconda3/envs/vmtk-env/python.exe')


def stepped_tube():
    # Two closed cylindrical segments with an exposed annular shoulder at z=0.
    return trimesh.creation.revolve(np.array([[0,-2],[.6,-2],[.6,0],[.8,0],[.8,2],[0,2]]),sections=48)


def uniform_tube():
    return trimesh.creation.cylinder(radius=.6,height=4,sections=48)


def test_watertight_shoulder_is_not_mistaken_for_smoothness():
    mesh=stepped_tube()
    assert mesh.is_watertight and mesh.is_volume
    feature,_=qc.feature_edges(mesh,[0,0,0],[0,0,1],1.1,.7,20)
    assert feature['circumferential_ring_count']==2
    assert feature['total_length_mm']>8


def test_uniform_tube_has_no_twenty_degree_ring_near_middle():
    feature,_=qc.feature_edges(uniform_tube(),[0,0,0],[0,0,1],1.,.8,20)
    assert feature['circumferential_ring_count']==0
    assert feature['edge_count']==0


def test_shoulder_normal_jump_exceeds_uniform_reference():
    first,_=qc.normal_jumps(stepped_tube(),[0,0,0],[0,0,1],1.1,2.1)
    second,_=qc.normal_jumps(uniform_tube(),[0,0,0],[0,0,1],1.1,1.1)
    assert first['maximum_deg']>=89
    assert second['maximum_deg']<8


def test_feature_measurements_invariant_under_rigid_transform():
    mesh=stepped_tube();a,_=qc.feature_edges(mesh,[0,0,0],[0,0,1],1.1,.7,20)
    t=np.eye(4);t[:3,:3]=Rotation.from_euler('xyz',[23,51,11],degrees=True).as_matrix();t[:3,3]=[100,20,-8]
    moved=mesh.copy();moved.apply_transform(t)
    b,_=qc.feature_edges(moved,t[:3,3],t[:3,2],1.1,.7,20)
    assert a['circumferential_ring_count']==b['circumferential_ring_count']
    assert b['total_length_mm']==pytest.approx(a['total_length_mm'],rel=1e-4)


def test_neighboring_branch_not_counted_as_local_feature_ring():
    neighbor=stepped_tube();neighbor.apply_translation([5,0,0])
    mesh=trimesh.util.concatenate([uniform_tube(),neighbor])
    features,_=qc.feature_edges(mesh,[0,0,0],[0,0,1],1.1,.8,20)
    assert features['edge_count']==0


def test_sections_detect_annular_area_discontinuity():
    mesh=stepped_tube()
    rows=[qc.cross_section(mesh,[0,0,z],[0,0,1],1.) for z in [-.25,.25]]
    assert rows[0]['equivalent_radius_mm']==pytest.approx(.6,rel=.01)
    assert rows[1]['equivalent_radius_mm']==pytest.approx(.8,rel=.01)
    assert qc.area_jump(rows)>.7


def test_section_chooses_known_local_vessel_not_another_contour():
    neighbor=trimesh.creation.cylinder(radius=2,height=4);neighbor.apply_translation([6,0,0])
    mesh=trimesh.util.concatenate([uniform_tube(),neighbor])
    value=qc.cross_section(mesh,[0,0,.1],[0,0,1],maximum_center_distance=.5)
    assert value['equivalent_radius_mm']==pytest.approx(.6,rel=.01)
    assert abs(value['centroid_x'])<1e-5


def test_normal_improvement_alone_cannot_override_sharp_ring():
    assert not qc.continuity_pass({'circumferential_ring_count':1},80,5,.01)


def test_normal_improvement_alone_cannot_override_area_jump():
    assert not qc.continuity_pass({'circumferential_ring_count':0},80,5,.2)


def test_ramp_setter_is_required_not_an_integer_enum_guess():
    with pytest.raises(RuntimeError,match='VMTK_RAMP_UNAVAILABLE'):
        require_ramp({'ramp_available':False,'preserve_shape_available':True,'capper_available':True})


def test_preserve_shape_api_is_required_for_requested_comparison():
    with pytest.raises(RuntimeError,match='VMTK_PRESERVE_SHAPE_API_UNAVAILABLE'):
        require_ramp({'ramp_available':True,'preserve_shape_available':False,'capper_available':True})


@pytest.mark.skipif(not PYTHON.exists(),reason='Existing Windows VMTK environment unavailable')
def test_real_vtp_roundtrip_preserves_geometry_and_data_arrays(tmp_path):
    before=set(sys.modules)
    result=bridge_smoke(PYTHON,tmp_path/'bridge')
    assert result['passed']
    assert not any(n.startswith('vmtk') for n in set(sys.modules)-before)


@pytest.mark.skipif(not PYTHON.exists(),reason='Existing Windows VMTK environment unavailable')
def test_real_worker_cannot_silently_substitute_tps_for_ramp(tmp_path):
    result=invoke(PYTHON,tmp_path/'blocked',operation='extension')
    assert result['status']=='FAIL' and result['returncode']==2
    assert result['error'] in ['VMTK_RAMP_UNAVAILABLE','O3_EXTENSION_PIPELINE_NOT_IMPLEMENTED']
    assert not (tmp_path/'blocked/roundtrip_output.vtp').exists()


@pytest.fixture
def real_audit():
    path=OUTPUT/'o3_smoothing_summary.json'
    if not path.exists():pytest.skip('Run s1-4 --o3-preflight to generate the real BG001 audit')
    return json.loads(path.read_text())


def test_real_audit_does_not_publish_unrepaired_stl_as_success(real_audit):
    assert real_audit['status']=='NEEDS_ADJUSTMENT'
    assert real_audit['failure_code']=='VMTK_RAMP_UNAVAILABLE'
    assert real_audit['final_stl'] is None
    assert real_audit['new_geometry_generated'] is False
    assert real_audit['boolean_union_performed'] is False
    assert real_audit['bambu_slicing_performed'] is False


def test_source_balanced_and_both_visualizers_are_unchanged(real_audit):
    assert real_audit['protection']['all_unchanged']
    assert real_audit['source_BALANCED_sha256']=='163f4bae4216d1f0c6a8e3da48e4d99cb8822d6d17e3105c01ed7153c93887bb'
    assert real_audit['ui_sha256']=={
        's1-2_swc_roi_generate_human.py':'3361cea44a97e8bfa907ff1504b0f977e329d35eb59cb9f57cc1b3a2a55cd729',
        's1-3_swc_roi_generate_MeVO.py':'348f4e9ae709f6d9aa934bb00d56a5881269dca3648bb13508f6f997347fccfa'}


@pytest.mark.parametrize('port',['I1','O1','O2'])
def test_frozen_port_audit_files_are_unchanged(real_audit,port):
    item=real_audit['frozen_routes'][port]
    assert item['before_sha256']==item['after_sha256']==sha256(OUTPUT/'baseline'/('frozen_'+port+'.stl'))


def test_clearance_floor_is_not_lowered(real_audit):
    assert real_audit['required_clearance_floor_mm']>=1.6
    assert real_audit['current_baseline']['prior_vessel_clearance_mm']>=1.6
    assert real_audit['production_clipping_performed'] is False
