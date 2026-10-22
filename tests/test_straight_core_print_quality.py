"""Actual frozen 1.6 mm core: rigid pose change, usable GUI metadata, no slicing."""
import json
from pathlib import Path
import zipfile

import numpy as np
import pytest
import yaml

from vascular_processing import straight_core_print_quality as q
from vascular_processing import support_slicer_adapter as slicer
from vascular_processing.all_abs_print_plate import read_objects


@pytest.fixture(scope='module')
def revision(tmp_path_factory):
    dest=tmp_path_factory.mktemp('quality')
    config=yaml.safe_load((q.ROOT/'config/straight_two_top_core_print_quality.yaml').read_text())
    config['output']=str(dest/'result')
    path=dest/'settings.yaml';path.write_text(yaml.safe_dump(config))
    # The production exporter must not discover or invoke a slicer at all.
    with pytest.MonkeyPatch.context() as mp:
        def forbidden(*a,**k):raise AssertionError('SLICING_IS_FORBIDDEN')
        mp.setattr(slicer,'discover_engine',forbidden)
        mp.setattr(slicer,'slice_local',forbidden)
        mp.setattr(q,'package',lambda out:dest/'not_published.zip')
        report=q.run(path)
    return Path(config['output']),report,config


def test_critical_tube_clear_of_bed_and_rigid_geometry(revision):
    _,r,_=revision;g=r['geometry']
    assert g['passed'] and r['original_files_unchanged']
    assert g['main_straight_bed_clearance_mm']==pytest.approx(4.4,abs=1e-5)
    assert g['main_straight_axis_in_print']==pytest.approx([0,0,1],abs=1e-9)
    assert g['print_height_mm']==pytest.approx(33,abs=1e-5)
    assert g['maximum_inverse_vertex_error_mm']<1e-7
    assert g['source_volume_mm3']==pytest.approx(g['print_volume_mm3'],rel=1e-8)
    assert g['diameter_mm']==1.6 and not g['main_segment_contacts_bed']
    assert not r['slicer_executed'] and not r['generated_support_meshes']


def test_real_archive_settings_and_object_override_precedence(revision):
    out,r,_=revision;p=out/r['project']['file']
    objects,settings,plates,crc=read_objects(p)
    assert crc and plates==1 and len(objects)==1
    with zipfile.ZipFile(p) as z:
        assert set(z.namelist())==q.PROJECT_MEMBERS
    assert settings['layer_height']=='0.08'
    assert settings['nozzle_diameter']==['0.4'] and settings['filament_type']==['ABS']
    for key,value in [('raft_layers','0'),('support_type','tree(auto)'),('support_style','tree_slim'),('enable_support','1')]:
        assert settings[key]==objects[0]['metadata'][key]==value
    assert settings['support_on_build_plate_only']=='1'
    assert float(settings['support_top_z_distance'])==pytest.approx(3*float(settings['layer_height']))
    assert settings['filament_colour'] and settings['filament_self_index']==['1','1']
    assert settings['outer_wall_speed']==['20','20']
    assert settings['outer_wall_acceleration']==['500','500']
    assert r['project']['gui_filament_qc']['passed']
    # Bambu restores unlisted fields from the system preset during GUI load.
    # A full snapshot must explicitly pin the fields in the corresponding group.
    for group in settings['different_settings_to_system']:
        assert set(revision[2]['settings']).issubset(set(group.split(';')))


def test_reject_mutated_source_before_export(tmp_path):
    cfg=yaml.safe_load((q.ROOT/'config/straight_two_top_core_print_quality.yaml').read_text())
    cfg['source_project_sha256']='0'*64;cfg['output']=str(tmp_path/'never_written')
    path=tmp_path/'wrong.yaml';path.write_text(yaml.safe_dump(cfg))
    with pytest.raises(ValueError,match='FROZEN_SOURCE_HASH_MISMATCH'):q.run(path)
    assert not Path(cfg['output']).exists()


def test_reject_horizontal_pose(revision):
    out,r,cfg=revision;source_cfg=q.load_config(q.ROOT/cfg['source_config'])
    old=q.ROOT/r['source_project'];before=read_objects(old)[0][0]['mesh']
    T=np.array(json.loads((old.parent/'print_layout.json').read_text())['records'][0]['casting_to_plate_transform'])
    with pytest.raises(ValueError,match='STRAIGHT_AXIS_MUST_BE_VERTICAL'):
        q.geometry_audit(before,before,T,T,source_cfg)


def test_gui_effective_check_rejects_silent_preset_reset(revision,tmp_path):
    out,r,cfg=revision;source=out/r['project']['file']
    with zipfile.ZipFile(source) as z:members={n:z.read(n) for n in z.namelist()}
    settings=json.loads(members['Metadata/project_settings.config'])
    settings['outer_wall_speed']=['60','60']
    members['Metadata/project_settings.config']=json.dumps(settings).encode()
    wrong=tmp_path/'reset.3mf'
    with zipfile.ZipFile(wrong,'w',zipfile.ZIP_DEFLATED) as z:
        for name,data in members.items():z.writestr(name,data)
    with pytest.raises(ValueError,match='GUI_SETTINGS_RESET'):
        q.verify_gui_roundtrip(source,wrong,cfg)


def test_gui_checks_effective_object_support_over_global_defaults(revision,tmp_path):
    out,r,cfg=revision;source=out/r['project']['file']
    with zipfile.ZipFile(source) as z:members={n:z.read(n) for n in z.namelist()}
    settings=json.loads(members['Metadata/project_settings.config'])
    settings.update(support_type='normal(auto)',support_style='default')
    members['Metadata/project_settings.config']=json.dumps(settings).encode()
    saved=tmp_path/'global_normalized.3mf'
    with zipfile.ZipFile(saved,'w',zipfile.ZIP_DEFLATED) as z:
        for name,data in members.items():z.writestr(name,data)
    result=q.verify_gui_roundtrip(source,saved,cfg)
    assert result['settings_verified']['support_type']=='tree(auto)'
    assert result['global_support_settings']['support_type']=='normal(auto)'
