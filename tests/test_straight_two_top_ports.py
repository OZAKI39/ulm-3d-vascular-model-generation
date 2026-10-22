"""Real CAD, mesh, labelling and native project checks for the added variant."""
import copy
from pathlib import Path

import numpy as np
import pyvista as pv
import pytest
import trimesh
import yaml

from vascular_processing import straight_channel_casting as s
from vascular_processing import straight_two_top_ports as top
from vascular_processing import bottomless_channel_fixture as f
from vascular_processing import rhino_channel_export as rhino
from vascular_processing.all_abs_print_plate import read_objects


@pytest.fixture(scope='module')
def real(tmp_path_factory):
    cfg=s.load_config(s.ROOT/'config/history/before_magnet_slots_20261012/straight_channel_two_top_ports_1p6_p1s_abs.yaml')
    out=tmp_path_factory.mktemp('two_top_ports')
    parts,g=s.geometry(cfg,out)
    channel=trimesh.load_mesh(out/'assembly/Straight_channel_without_collars_reference.stl')
    return cfg,out,parts,g,channel


def test_exactly_two_top_openings_and_no_horizontal_stub(real):
    cfg,out,parts,g,channel=real
    assert [p['name'] for p in g['ports']]==['INLET','OUTLET1']
    assert all(p['outward_axis']==[0,0,1] and p['minimum_cap_normal_dot']>.99999 for p in g['ports'])
    assert g['closed_chamber_qc']['top_opening_count']==2
    assert channel.bounds[1,0]==pytest.approx(12.8,abs=2e-5)
    assert not channel.contains([[15,0,2.3],[12,0,2.3]]).any()
    assert channel.contains([[-7,0,2.3],[0,0,2.3],[12,0,8]]).all()
    assert all(m.is_volume and len(m.split())==1 for m in parts.values())
    assert len(parts)==6 and all(g['gates'].values())
    assert 'horizontal_extension_inside_pdms' not in g['gates']
    assert g['design']['horizontal_end'] is None
    assert not g['closed_chamber_qc']['manual_horizontal_cut_required']


def test_vtp_has_no_outlet2_and_cap_labels_are_on_top(real):
    _,out,_,g,_=real
    m=pv.read(out/'labelled_straight_lumen_mm.vtp')
    assert set(m.cell_data['ModelFaceID'])=={1,2,3}
    assert list(m.field_data['ModelFaceNames'])==['WALL','INLET','OUTLET1']
    assert float(m.field_data['coordinate_scale_to_m'][0])==.001
    for label in [2,3]:
        port=m.extract_cells(np.where(m.cell_data['ModelFaceID']==label)[0])
        assert np.allclose(port.points[:,2],18.4,atol=1e-5)
    assert g['closed_chamber_qc']['volume_relative_error']<.003


def test_smooth_elbows_straight_run_and_risers_are_1p6(real):
    _,_,_,g,channel=real
    assert len(g['section_checks'])==23
    assert g['design']['straight_centerline_length_mm']==pytest.approx(24.2)
    for row in g['section_checks']:
        assert row['equivalent_diameter_mm']==pytest.approx(1.6,abs=.003)
    path=np.array(top.centerline_samples(g['design']))
    assert path[:,0].max()<=12 and channel.contains(path).all()
    assert g['ports'][0]['top_center_mm']==pytest.approx([-17,0,18.4])
    assert g['ports'][1]['top_center_mm']==pytest.approx([12,0,18.4])


def test_old_three_port_channel_is_rejected_as_two_ports(real):
    _,_,_,g,_=real
    old=s.load_config(s.ROOT/'config/straight_channel_1p6_p1s_abs.yaml')
    channel=trimesh.load_mesh(old['output']/'assembly/Straight_channel_without_collars_reference.stl')
    assert channel.contains([[15,0,2.3]]).all()
    # Even if only two names are supplied, the actual third opening/stub fails.
    result=top.channel_qc(g['design'],channel,g['ports'])
    assert not result['passed']


def test_existing_inlet_and_low_straight_body_are_preserved(real):
    _,_,parts,g,_=real
    old=s.load_config(s.ROOT/'config/straight_channel_1p6_p1s_abs.yaml')
    old_core=trimesh.load_mesh(old['output']/'assembly/Straight_lumen_core.stl')
    clip=s.box([[-35,-4,0],[8.5,4,21]])
    a=trimesh.boolean.intersection([old_core,clip],engine='manifold')
    b=trimesh.boolean.intersection([parts['Straight_lumen_core'],clip],engine='manifold')
    assert a.volume+b.volume-2*s.volume_intersection(a,b)<.005
    assert parts['Straight_lumen_core'].bounds[0,2]==pytest.approx(1.5,abs=1e-5)


def test_frame_shrinks_and_raised_retainer_connections_still_fit(real):
    _,_,parts,g,_=real
    assert g['cavity_xy_size_mm']==pytest.approx([36.6,9.6],abs=1e-5)
    assert parts['Straight_bottomless_box_wall2mm'].extents==pytest.approx([49.8,21.8,10.4],abs=2e-5)
    assert g['pdms_fill_height_mm']==pytest.approx(9.6)
    assert g['wall_mm']==2 and g['base_mm']==0
    assert all(g['fixture_qc']['gates'].values())
    gap=g['fixture_qc']['retainer_rim_clearance_qc']
    assert gap['passed'] and all(r['measured_wall_to_arm_gaps_mm']==pytest.approx([3,3,3],abs=1e-5) for r in gap['parts'])
    assert g['manufacturing_height_qc']['retainer_top_z_mm']==pytest.approx(17.8)
    assert g['manufacturing_height_qc']['core_top_z_mm']==pytest.approx(18.4)


def test_rhino_and_native_project_include_six_parts(real):
    cfg,out,parts,g,_=real
    meshes,records=s.arrange(parts,out)
    rhino.export_models(cfg,parts,g,out)
    assert rhino.verify_exports(out)['passed']
    dest=out/'test.3mf'
    s.write_project(dest,meshes,records,{k:{'name':k} for k in ['machine','process','filament']},cfg)
    objects,_,plates,crc=read_objects(dest)
    assert crc and plates==1 and len(objects)==6
    assert {o['metadata']['name'] for o in objects}==set(parts)
    for obj in objects:
        meta=obj['metadata'];name=meta['name']
        assert meta['enable_support']==('1' if 'core' in name or 'retainer' in name else '0')
        assert meta['raft_layers']==('2' if 'core' in name else '0')


def test_variant_is_additional_and_rejects_horizontal_outlet_config(real,tmp_path):
    cfg,_,_,_,_=real
    original=s.load_config(s.ROOT/'config/straight_channel_1p6_p1s_abs.yaml')
    assert not s.two_top_ports(original) and original['output']!=cfg['output']
    bad=copy.deepcopy(cfg);bad['horizontal_extension_after_outlet1_mm']=11.15
    file=tmp_path/'mixed.yaml'
    file.write_text(yaml.safe_dump({k:str(v) if isinstance(v,Path) else v for k,v in bad.items()}))
    with pytest.raises(ValueError,match='TWO_TOP_PORTS_MUST_NOT_DECLARE_HORIZONTAL_OUTLET'):
        s.load_config(file)
