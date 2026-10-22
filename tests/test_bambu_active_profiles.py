"""Current user profiles must not silently fall back to historical hardware."""
import json
import pytest
from vascular_processing import bambu_manufacturing as b


@pytest.fixture
def presets(tmp_path):
    def write(kind,name,values,user=False):
        p=tmp_path/('user/123' if user else 'system/BBL')/kind/(name+'.json')
        p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(dict(name=name,**values,**({'from':'User'} if user else {}))))
        return p
    write('machine','Start template',dict(machine_start_gcode='G28 ; 0.2 mm template'))
    write('machine','X1 0.4',dict(printer_model='Bambu Lab X1 Carbon',nozzle_diameter=['0.4'],
        printable_area=['0x0','256x0','256x256','0x256'],printable_height='250',machine_start_gcode='old'))
    write('machine','X1 0.2',dict(inherits='X1 0.4',nozzle_diameter=['0.2'],include=['Start template']))
    write('process','Standard',dict(layer_height='0.1',line_width='0.22',outer_wall_speed=['60','70'],compatible_printers=['X1 0.2']))
    write('process','Vessel',dict(inherits='Standard',outer_wall_speed=['30','nil']),True)
    write('filament','ABS',dict(filament_type=['ABS'],nozzle_temperature=['260','270'],compatible_printers=['X1 0.2']))
    write('filament','ABS Vessel',dict(inherits='ABS',nozzle_temperature=['265','nil']),True)
    return tmp_path,dict(machine='X1 0.2',process='Vessel',filaments=['ABS Vessel']),write


def test_current_user_profiles_and_include_templates(presets):
    root,active,_=presets;r=b.active_profiles(root,active,'ABS')
    assert r['nozzle_diameter_mm']==.2 and r['model']=='Bambu Lab X1 Carbon'
    p=r['flattened'];assert p['machine']['machine_start_gcode']=='G28 ; 0.2 mm template'
    assert p['process']['layer_height']=='0.1' and p['process']['outer_wall_speed']==['30','70']
    assert p['filament']['nozzle_temperature']==['265','270']
    assert len(r['profile_sources'])==7
    assert all('inherits' not in v and 'include' not in v for v in p.values())


def test_wrong_active_material_rejected(presets):
    root,active,write=presets;write('filament','ABS Vessel',dict(filament_type=['PLA']),True)
    with pytest.raises(ValueError,match='MATERIAL_MISMATCH'):b.active_profiles(root,active,'ABS')


def test_04_process_cannot_match_02_via_machine_inheritance(presets):
    root,active,write=presets;write('process','Vessel',dict(inherits='Standard',compatible_printers=['X1 0.4']),True)
    with pytest.raises(ValueError,match='PRINTER_MISMATCH'):b.active_profiles(root,active,'ABS')


def test_missing_active_profile_does_not_fall_back(presets):
    root,active,_=presets;active['process']='Missing'
    with pytest.raises(FileNotFoundError):b.active_profiles(root,active,'ABS')


def test_conflicting_profile_names_rejected(presets):
    root,active,write=presets;write('process','Vessel',dict(layer_height='0.2'))
    with pytest.raises(ValueError,match='AMBIGUOUS'):b.active_profiles(root,active,'ABS')


def test_include_cycle_rejected(presets):
    _,_,write=presets;p=write('process','Cycle',dict(include=['Cycle']))
    with pytest.raises(ValueError,match='Cyclic'):b.flatten(p,resolve_includes=True)


def test_local_settings_override_included_template(presets):
    _,_,write=presets
    p=write('machine','Override',dict(inherits='X1 0.4',include=['Start template'],machine_start_gcode='explicit'))
    assert b.flatten(p,resolve_includes=True)[0]['machine_start_gcode']=='explicit'
