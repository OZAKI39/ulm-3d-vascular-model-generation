"""P1S printer identity, separate-output protection and real regeneration evidence."""
import copy
import json
import zipfile
from pathlib import Path
import pytest
from vascular_processing.core_support_p1s import default_output, validate_profiles, evaluated, PRESETS
from vascular_processing.support_inputs import digest

@pytest.fixture
def profiles():
    return json.loads((default_output()/'bambu/profiles/profiles.json').read_text())['actual']

def test_shared_process_name_does_not_change_machine(profiles):
    assert 'X1C' in profiles['process']['name']
    assert PRESETS['machine'] in profiles['process']['compatible_printers']
    assert validate_profiles(profiles)
    assert profiles['machine']['printer_model']=='Bambu Lab P1S'

@pytest.mark.parametrize('kind,key,value,message',[
    ('machine','printer_model','Bambu Lab X1 Carbon','P1S_MACHINE'),
    ('machine','nozzle_diameter',['0.2'],'P1S_04_NOZZLE'),
    ('filament','filament_type',['PLA'],'ABS_REQUIRED'),
    ('process','compatible_printers',['Bambu Lab X1 Carbon 0.4 nozzle'],'PROCESS_NOT_COMPATIBLE'),
])
def test_wrong_hardware_or_material_is_rejected(profiles,kind,key,value,message):
    candidate=copy.deepcopy(profiles);candidate[kind][key]=value
    with pytest.raises(ValueError,match=message):validate_profiles(candidate)

def test_support_dimensions_follow_new_extrusion_width():
    out=default_output();record=json.loads((out/'candidates/P1S_04_BASE/geometry.json').read_text())
    assert record['curtain_width_mm']==pytest.approx(.84)
    assert record['roof_width_mm']==pytest.approx(.84)
    assert record['gap_mm']==pytest.approx(.2)

def test_review_cannot_start_a_missing_slice(tmp_path,monkeypatch):
    def forbidden(*args,**kwargs):
        pytest.fail('review must not launch the native slicer')
    monkeypatch.setattr('vascular_processing.core_support_p1s.evaluate',forbidden)
    with pytest.raises(FileNotFoundError,match='REVIEW_REQUIRES_COMPLETED_EVALUATION'):
        evaluated(tmp_path,{},readonly=True)

def test_archived_profile_review_does_not_probe_old_platform(tmp_path,monkeypatch):
    from vascular_processing.core_support_p1s import profiles as resolve_profiles
    import shutil
    for p in (default_output()/'bambu/profiles').glob('*.json'):
        shutil.copy2(p,tmp_path/p.name)
    def forbidden(*args,**kwargs):
        pytest.fail('historical archive review must not resolve or launch a new engine')
    monkeypatch.setattr('vascular_processing.core_support_p1s.native.active_profiles',forbidden)
    monkeypatch.setattr('vascular_processing.support_slicer_adapter.discover_engine',forbidden)
    monkeypatch.setenv('VASCULAR_BAMBU_BIN','/missing/original/windows/slicer.exe')
    assert resolve_profiles(tmp_path,review_only=True)['actual']['machine']['nozzle_diameter']==['0.4']

def test_migrated_review_does_not_refreeze_or_rewrite_sources(monkeypatch):
    from vascular_processing.core_support_p1s import run
    def forbidden(*args,**kwargs):
        pytest.fail('review must preserve archived source paths and bytes')
    monkeypatch.setattr('vascular_processing.core_support_p1s.freeze',forbidden)
    monkeypatch.setattr('vascular_processing.core_support_p1s.deliver',forbidden)
    path=default_output()/'source/discovery.json';before=digest(path)
    assert run(stage='review')['status']=='NEEDS_ADJUSTMENT'
    assert digest(path)==before

def test_real_p1s_delivery_and_frozen_anatomy():
    out=default_output();summary=json.loads((out/'regeneration_summary.json').read_text())
    assert summary['scope']=='VASCULAR_CORE_AND_CUSTOM_SUPPORT_ONLY'
    assert summary['source_protection']['all_unchanged']
    assert summary['historical_protection']['all_unchanged']
    assert summary['saved_profile_protection']['all_unchanged']
    assert not summary['box_or_retainers_generated'] and not summary['anatomy_scaled']
    assert summary['archive_validation']['passed']
    assert summary['production_frozen']==all(summary['gates'].values())
    assert digest(out/'selected/frozen_vascular_core.stl')=='84f1989392fb8d210432750a8022a830943b9bff890b044aa9d310822b619ce0'
    for name,expected in summary['files'].items():assert digest(out/name)==expected
    with zipfile.ZipFile(out/'selected/BG001_RMCA_core_support_P1S_04_ABS_sliced.3mf') as z:
        assert z.testzip() is None
        cfg=json.loads(z.read('Metadata/project_settings.config'))
        assert cfg['printer_model']=='Bambu Lab P1S'
        assert cfg['nozzle_diameter']==['0.4'] and cfg['filament_type']==['ABS']
        assert cfg['enable_support']=='0'
        assert len([n for n in z.namelist() if n.endswith('.gcode')])==1
