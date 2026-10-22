"""Real-core orientation search, finite clearance and no automatic slicing."""
import json
from pathlib import Path
import numpy as np
import pytest
import yaml

from vascular_processing import straight_core_orientation as o
from vascular_processing import support_slicer_adapter as slicer


@pytest.fixture(scope='module')
def study(tmp_path_factory):
    dest=tmp_path_factory.mktemp('orientations')
    config=yaml.safe_load((o.q.ROOT/'config/straight_two_top_core_orientation.yaml').read_text())
    config['output']=str(dest/'result')
    path=dest/'study.yaml';path.write_text(yaml.safe_dump(config))
    with pytest.MonkeyPatch.context() as mp:
        def forbidden(*a,**k):raise AssertionError('NATIVE_SLICING_NOT_AUTHORIZED')
        mp.setattr(slicer,'discover_engine',forbidden);mp.setattr(slicer,'slice_local',forbidden)
        mp.setattr(o,'package_study',lambda out:dest/'not_published.zip')
        result=o.run_study(path)
    return Path(config['output']),result


def test_search_is_diverse_unique_and_contains_exact_baselines():
    poses=o.candidates(.5,1000)
    assert len(poses)==1721
    assert len({tuple(np.round(r[2],10)) for _,r in poses})==len(poses)
    assert {'tilt45','vertical_previous','horizontal_clearance','horizontal_original'}<={n for n,_ in poses}
    for _,r in poses:
        assert np.allclose(r.T@r,np.eye(3),atol=1e-8)
        assert np.linalg.det(r)==pytest.approx(1.)


def test_selected_pose_reduces_geometric_burden_without_claiming_sliced_volume(study):
    _,r=study
    assert r['tested_candidates']==1721
    assert r['selected']['main_angle_to_bed_deg']==pytest.approx(45,abs=.5)
    assert r['selected']['main_clearance_mm']>.5
    assert r['selected']['first_layer']['area_mm2']>=1
    a=r['selected']['metrics_by_angle']['45'];b=r['previous']['metrics_by_angle']['45']
    assert a['overhang_surface_mm2']<b['overhang_surface_mm2']*.2
    assert a['height_weighted_projected_proxy_mm3']<b['height_weighted_projected_proxy_mm3']*.1
    assert not r['global_optimality_proven'] and not r['slicer_executed']
    assert r['generated_support_meshes']==0 and not r['physical_print_tested']


def test_horizontal_geometry_is_supported_by_original_collars_not_new_mesh(study):
    _,r=study;record=r['deliveries'][1];g=record['geometry']
    assert g['main_straight_angle_to_bed_deg']==pytest.approx(0,abs=1e-8)
    assert g['main_straight_clearance_mm']==pytest.approx(1.2,abs=.001)
    assert g['print_height_mm']==pytest.approx(4,abs=.001)
    assert record['metrics']['first_layer']['components']==2
    assert g['same_triangles'] and not g['extra_feet_or_support_meshes_added']


def test_both_real_projects_preserve_volume_and_editable_settings(study):
    out,r=study
    for record in r['deliveries']:
        folder=out/record['folder'];p=folder/record['project']['file']
        objs,s,plates,crc=o.q.read_objects(p)
        assert crc and plates==1 and len(objs)==1
        expected=json.loads((folder/'effective_settings.json').read_text())['settings']
        assert all(s[k]==v for k,v in expected.items())
        assert objs[0]['metadata']['raft_layers']=='0'
        assert objs[0]['metadata']['support_type']==expected['support_type']
        assert s['support_threshold_angle']=='45' and s['layer_height']=='0.08'
        assert s['nozzle_diameter']==['0.4'] and s['filament_type']==['ABS']
        assert record['geometry']['maximum_inverse_error_mm']<1e-7
        assert record['geometry']['source_volume_mm3']==pytest.approx(record['geometry']['print_volume_mm3'])
        assert abs(objs[0]['mesh'].bounds[0,2])<1e-7


def test_pose_audit_rejects_scaled_mesh(study):
    out,r=study;record=r['deliveries'][0]
    p=out/record['folder']/record['project']['file'];after=o.q.read_objects(p)[0][0]['mesh']
    T=np.array(record['geometry']['casting_to_print_transform']);before=after.copy();before.apply_transform(np.linalg.inv(T))
    after.apply_scale(1.01)
    config=o.q.load_config(o.q.ROOT/'config/straight_channel_two_top_ports_1p6_p1s_abs.yaml')
    with pytest.raises(ValueError,match='RIGID_CORE_GEOMETRY_CHECK_FAILED'):o.pose_audit(before,after,T,config)
