"""ABS-only cores, PLA fixtures, and a hard boundary against automatic slicing."""
import copy
import json
from pathlib import Path
import re
import shutil
import zipfile
import xml.etree.ElementTree as ET

import numpy as np
import pytest
import rhino3dm
import trimesh
import yaml

from vascular_processing import y_channel_casting as y,straight_channel_casting as s
from vascular_processing import channel_manual_delivery as manual
from vascular_processing.all_abs_print_plate import read_objects
from vascular_processing import support_slicer_adapter

CONFIGS=['y_channel_p1s_abs.yaml','straight_channel_1p6_p1s_abs.yaml','straight_channel_two_top_ports_1p6_p1s_abs.yaml']


@pytest.fixture(params=CONFIGS)
def current(request):
    module=y if request.param.startswith('y_') else s
    cfg=module.load_config(y.ROOT/'config'/request.param);out=cfg['output']
    summary=json.loads((out/'manual_slicing_summary.json').read_text())
    return module,cfg,out,summary


def test_delivered_projects_have_correct_real_material_and_no_toolpaths(current):
    _,cfg,out,summary=current
    assert summary['status']=='READY_FOR_USER_SLICING' and not summary['slicer_executed']
    assert summary['actual_print_time_seconds'] is None
    assert not list(out.glob('*sliced.3mf')) and not (out/'bambu').exists()
    seen=[]
    for row in summary['projects']:
        path=out/row['file']
        with zipfile.ZipFile(path) as z:
            assert set(z.namelist())==manual.PROJECT_MEMBERS and z.testzip() is None
            name_qc=manual.validate_bambu_names(path.name,z.read('3D/3dmodel.model'),z.read('Metadata/model_settings.config'))
        assert name_qc['passed'] and len(name_qc['plate_names'])==1
        assert row['material'] in name_qc['plate_names'][0]
        objects,settings,count,crc=read_objects(path)
        assert manual.validate_gui_filaments(settings)['passed']
        assert count==1 and crc and settings['nozzle_diameter']==['0.4']
        assert settings['filament_type']==[row['material']]
        for obj in objects:
            name=obj['metadata']['name'];seen.append(name)
            assert row['material']==('ABS' if 'core' in name else 'PLA')
            expected=trimesh.load_mesh(out/'print_parts'/(name+'.stl'))
            # STL round-trip may reorder vertices, so compare actual surfaces.
            assert obj['mesh'].bounds==pytest.approx(expected.bounds,abs=2e-5)
            assert trimesh.proximity.closest_point(expected,obj['mesh'].vertices)[1].max()<2e-5
        if row['material']=='PLA':
            assert 'PLA' in settings['filament_settings_id'][0]
            assert settings['nozzle_temperature']==['220','220']
    assert len(seen)==len(set(seen))==len(summary['part_materials'])


def test_rhino_material_attributes_and_geometry_preserved(current):
    _,cfg,out,summary=current
    source=y.ROOT/'logs/manual_slice_ABS_PLA_20261012/baseline.json'
    baseline=next(r for r in json.loads(source.read_text()).values() if r['new_output']==str(out.relative_to(y.ROOT)))
    for relative,h in baseline['assembly_sha256'].items():
        if relative.endswith('.step'):
            previous=y.ROOT/baseline['old_output']/relative
            if not previous.exists():
                previous=y.ROOT/'history/before_manual_slice_ABS_PLA_20261012'/baseline['old_output']/relative
            assert y.digest(previous)==h
            # Open CASCADE updates the creation timestamp on every STEP export.
            # Compare every other byte, including the complete BRep entity data.
            def without_creation_time(p):
                return re.sub(rb"(FILE_NAME\('[^']*',')[^']*(')",rb"\1<creation time>\2",p.read_bytes(),count=1)
            assert without_creation_time(out/relative)==without_creation_time(previous)
        else:
            assert y.digest(out/relative)==h
    for path in out.glob('*.3dm'):
        doc=rhino3dm.File3dm.Read(str(path));assert doc.Settings.ModelUnitSystem==rhino3dm.UnitSystem.Millimeters
        for obj in doc.Objects:
            assert obj.Attributes.GetUserString('print_material')==('ABS' if 'core' in obj.Attributes.Name else 'PLA')
    assert 'abs_dimensional_allowance' not in cfg
    for axis in ['width','height']:
        q=y.dimension_budget(cfg,axis)
        assert q['required_cad_minimum_mm']==9.6 and q['material']=='PLA'
        assert q['assumed_material_shrink_fraction'] is None


def test_project_stage_cannot_invoke_engine_or_slice(current,tmp_path,monkeypatch):
    module,cfg,out,_=current;cfg=copy.deepcopy(cfg)
    dest=tmp_path/'output';shutil.copytree(out,dest)
    cfg['output']=dest
    y.save(dest/'effective_config.json',cfg)
    config=tmp_path/'config.yaml';config.write_text(yaml.safe_dump({k:str(v) if isinstance(v,Path) else v for k,v in cfg.items()}))
    def forbidden(*a,**kw):raise AssertionError('NO_SLICER_MAY_BE_CALLED')
    monkeypatch.setattr(y,'slice_plate',forbidden);monkeypatch.setattr(s,'slice_plate',forbidden)
    monkeypatch.setattr(support_slicer_adapter,'discover_engine',forbidden)
    monkeypatch.setattr(support_slicer_adapter,'slice_local',forbidden)
    # Current real delivery ZIP is checked separately; keep test ZIPs out of
    # the user's deliverables directory while exercising the real export.
    monkeypatch.setattr(manual,'package',lambda *args:None)
    result=module.run(config,stage='project')
    assert result['status']=='READY_FOR_USER_SLICING' and not result['slicer_executed']


def test_slice_entry_is_blocked_and_stale_toolpath_members_are_rejected(current,tmp_path):
    module,cfg,out,summary=current
    with pytest.raises(ValueError,match='AUTOMATIC_SLICING_DISABLED'):module.run(Path('not_needed.yaml'),stage='slice')
    with pytest.raises(ValueError,match='AUTOMATIC_SLICING_DISABLED'):module.slice_plate(cfg,{},None,tmp_path)
    path=tmp_path/'bad.3mf';shutil.copy2(out/summary['projects'][0]['file'],path)
    with zipfile.ZipFile(path,'a') as z:z.writestr('Metadata/plate_1.gcode','G1 X1 Y1 E1')
    with pytest.raises(ValueError,match='UNSLICED_PROJECT_HAS_UNEXPECTED_MEMBERS'):
        manual.validate_project(path,{},'ABS',cfg)


def test_package_contains_models_and_only_unsliced_projects(current):
    _,_,out,_=current
    receipt=json.loads((out/'delivery_receipt.json').read_text())
    archive=y.ROOT/'deliverables'/Path(receipt['zip']).name
    assert y.digest(archive)==receipt['zip_sha256'] and not receipt['sliced_artifacts_included']
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        assert not any(n.endswith(('.gcode','.bgcode')) or n.endswith('sliced.3mf') or '/bambu/' in n for n in z.namelist())
        assert len([n for n in z.namelist() if n.endswith('_project.3mf')])==2


@pytest.mark.parametrize('bad',list('<>:"/\\|?*')+['\n','\t'])
@pytest.mark.parametrize('field',['plater_name','name'])
def test_bambu_forbidden_names_rejected_after_xml_unescaping(bad,field):
    config=ET.Element('config');plate=ET.SubElement(config,'plate')
    ET.SubElement(plate,'metadata',key=field,value='Example'+bad+'ABS')
    with pytest.raises(ValueError,match='BAMBU_INVALID_NAME'):
        manual.validate_bambu_names('good_project.3mf',b'<model/>',ET.tostring(config))


def test_legacy_pipe_name_is_rejected_in_real_project(current,tmp_path):
    _,cfg,out,summary=current;src=out/summary['projects'][0]['file']
    with zipfile.ZipFile(src) as z:members={n:z.read(n) for n in z.namelist()}
    xml=ET.fromstring(members['Metadata/model_settings.config'])
    for node in xml.findall('.//plate/metadata'):
        if node.get('key')=='plater_name':node.set('value','Y | core ABS | UNSLICED')
    members['Metadata/model_settings.config']=ET.tostring(xml)
    bad=tmp_path/'legacy_name.3mf'
    with zipfile.ZipFile(bad,'w',zipfile.ZIP_DEFLATED) as z:
        for n,b in members.items():z.writestr(n,b)
    with pytest.raises(ValueError,match='BAMBU_INVALID_NAME: plater_name'):
        manual.validate_project(bad,{},'ABS',cfg)


@pytest.mark.parametrize('field',['filament_colour','filament_self_index','filament_ids'])
def test_incomplete_project_spool_metadata_is_rejected(current,tmp_path,field):
    _,cfg,out,summary=current;src=out/summary['projects'][0]['file']
    with zipfile.ZipFile(src) as z:members={n:z.read(n) for n in z.namelist()}
    settings=json.loads(members['Metadata/project_settings.config']);settings.pop(field)
    members['Metadata/project_settings.config']=json.dumps(settings).encode()
    bad=tmp_path/'missing_gui_metadata.3mf'
    with zipfile.ZipFile(bad,'w',zipfile.ZIP_DEFLATED) as z:
        for n,b in members.items():z.writestr(n,b)
    with pytest.raises(ValueError,match='BAMBU_GUI_'):
        manual.validate_project(bad,{},'ABS',cfg)


def test_two_flow_variants_belong_to_one_spool(current):
    _,_,out,summary=current
    for row in summary['projects']:
        _,settings,_,_=read_objects(out/row['file'])
        assert len(settings['filament_extruder_variant'])==2
        assert settings['filament_self_index']==['1','1']
        assert len(settings['filament_colour'])==len(settings['filament_settings_id'])==1
        bad=dict(settings,filament_self_index=['1','2'])
        with pytest.raises(ValueError,match='BAMBU_GUI_FILAMENT_VARIANT_MAPPING_INVALID'):
            manual.validate_gui_filaments(bad)
